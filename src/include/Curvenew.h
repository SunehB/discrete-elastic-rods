#include "glm/exponential.hpp"
#include "glm/fwd.hpp"
#include "glm/geometric.hpp"
#include "glm/matrix.hpp"
#include "polyscope/curve_network.h"
#include "polyscope/utilities.h"
#include <cmath>
#include <cstddef>
// #include <iostream>
#include <mach/task_info.h>
#include <numeric>
// #include <ostream>
#include <vector>

using namespace polyscope;
class Curven {
    using point = glm::vec3;

  public:
    Curven(){};
    ~Curven(){};

    void resize_vectors() {
        edge_length.resize(num_edges);
        arc_length.resize(num_edges);
        turn_angles.resize(num_controlpoints);
        tangent_on_edges.resize(num_edges);
        normal_on_edges.resize(num_edges);
        binormal_on_edges.resize(num_edges);
        curvature.resize(num_controlpoints);
        darboux.resize(num_controlpoints);
        vertex_weight.resize(num_controlpoints);

        twist_thetas.resize(num_edges);
        material_u.resize(num_edges);
        material_v.resize(num_edges);

        velocity.resize(num_controlpoints);
        acceleration.resize(num_controlpoints);
        bending_force.resize(num_controlpoints);
        twisting_force.resize(num_controlpoints);
    }


    void initcurve() {
        size_t nSamples = 10;
        double dt = 1.0 / nSamples + 1;

        for (size_t i = 0; i < nSamples; ++i) {
            double t = i * dt;
            point p;
            p.x = cos(2 * M_PI * t);
            p.y = sin(2 * M_PI * t);
            p.z = 0.3 * sin(4 * M_PI * t);
            controlpoints.push_back(p);
        }

        num_controlpoints = controlpoints.size();
        num_edges = nSamples;
        resize_vectors();

        for (size_t i = 0; i < nSamples; i++) {
            curveEdgesId.push_back({i, (i + 1) % nSamples});
            curveEdges.push_back(controlpoints[(i + 1) % nSamples] - controlpoints[i]);
            edge_length[i] = glm::length(curveEdges[i]);
        }

        curve = registerCurveNetwork("Spiral", controlpoints, curveEdgesId);

        cal_tangent();
        normal_on_edges[0] = glm::normalize(glm::cross(reference, tangent_on_edges[0]));
        binormal_on_edges[0] = glm::normalize(glm::cross(tangent_on_edges[0], normal_on_edges[0]));

        update_bishop();
        update_material_frame();
        cal_attrs();
        calforce();
    }

    void cal_tangent() {
        for (size_t i = 0; i < num_edges; i++) {
            tangent_on_edges[i] = glm::normalize(curveEdges[i]);
        }
        curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
    }

    void update_bishop() {
        for (size_t i = 1; i < num_edges; i++) {
            point t_i = tangent_on_edges[i];
            point t_i_1 = tangent_on_edges[(i - 1) % num_edges];

            point axis = glm::cross(t_i_1, t_i);
            float theta = std::atan2(glm::length(axis), glm::dot(t_i_1, t_i));
            glm::mat3 rotationMatrix = glm::rotate(glm::mat4(1.0f), theta, axis);

            if (theta > 1e-6) {
                normal_on_edges[i] = rotationMatrix * normal_on_edges[(i - 1) % num_edges];
            } else {
                normal_on_edges[i] = normal_on_edges[(i - 1) % num_edges];
            }
        }
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 0; i < num_edges; i++) {
            binormal_on_edges[i] = glm::normalize(glm::cross(tangent_on_edges[i], normal_on_edges[i]));
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);
    }

    void update_material_frame() {
        for (size_t i = 0; i < num_edges; i++) {
            glm::mat3 rotationMatrix = glm::rotate(glm::mat4(1.0f), totaltwist, tangent_on_edges[i]);
            material_u[i] = rotationMatrix * normal_on_edges[i];
            material_v[i] = rotationMatrix * binormal_on_edges[i];
        }
        curve->addEdgeVectorQuantity("material_u", material_u);
        curve->addEdgeVectorQuantity("material_v", material_v);
    }

    void cal_attrs() {
        // arc_length[0] = 0;
        // for (size_t i = 1; i < num_edges; i++) {
        //     arc_length[i] = arc_length[i - 1] + glm::length(curveEdges[i]) / totallength;
        // }
        // curve->addEdgeScalarQuantity("arc_length", arc_length);

        for (size_t i = 0; i < num_controlpoints; i++) {
            point Ledge = curveEdges[(i - 1 + num_controlpoints) % num_controlpoints];
            point Redge = curveEdges[(i) % num_controlpoints];

            vertex_weight[i] = 0.5 * (glm::length(Ledge) + glm::length(Redge));
            // turn_angles[i] = PI - acos(glm::dot(glm::normalize(Ledge), glm::normalize(Redge)));

            darboux[i] =
                2.f * glm::cross(Ledge, Redge) / (glm::length(Ledge) * glm::length(Redge) + glm::dot(Ledge, Redge));
            curvature[i] = glm::dot(darboux[i], darboux[i]);
        }
        curve->addNodeScalarQuantity("vertex_weight", vertex_weight);
        // curve->addNodeScalarQuantity("turning angles", turn_angles);
        curve->addNodeScalarQuantity("curvature", curvature);
        curve->addNodeVectorQuantity("darboux", darboux);
    }

    void cal_bending_force() {
        bending_force.clear();
        bending_force.resize(num_controlpoints);

        for (size_t i = 0; i < num_controlpoints; i++) {
            auto cur_edge = curveEdges[i];
            auto prev_edge = curveEdges[(i - 1 + num_controlpoints) % num_controlpoints];
            float co = 2 / (vertex_weight[i] *
                            (glm::length(cur_edge) * glm::length(prev_edge) + glm::dot(cur_edge, prev_edge)));

            bending_force[i] += 2.f * glm::cross(cur_edge, darboux[i]) + glm::dot(darboux[i], cur_edge) * darboux[i];
            bending_force[i] += 2.f * glm::cross(prev_edge, darboux[i]) - glm::dot(darboux[i], prev_edge) * darboux[i];
            bending_force[i] *= co;
        }
        curve->addNodeVectorQuantity("bending force", bending_force);
    }

    void cal_twisting_force() {
        twisting_force.clear();
        twisting_force.resize(num_controlpoints);

        for (size_t i = 1; i < num_controlpoints; i++) {
            auto prevd = 0.5f * darboux[i] / edge_length[(i - 1 + num_controlpoints) % num_controlpoints];
            auto nextd = -0.5f * darboux[i] / edge_length[i];
            float L = 0.5 * std::accumulate(vertex_weight.begin(), vertex_weight.begin() + i, 0.0f);
            twisting_force[i] = -1 * totaltwist * (nextd + prevd) / (L + eps);
        }
        twisting_force[0] = -1 * totaltwist * 0.5f * darboux[0] / vertex_weight[0];
        curve->addNodeVectorQuantity("twisting force", twisting_force);
    }

    void calforce() {
        cal_bending_force();
        cal_twisting_force();
    }

    void cal_velocity() {
        acceleration.clear();
        acceleration.resize(num_controlpoints);
        velocity.clear();
        velocity.resize(num_controlpoints);

        float dt = 0.0001;
        for (size_t i = 0; i < num_controlpoints; i++) {
            point F_total = -bending_force[i] + twisting_force[i];
            acceleration[i] = F_total / vertex_weight[i];
            point new_v = velocity[i] + acceleration[i] * dt;
            velocity[i] = new_v;
        }
        curve->addNodeVectorQuantity("velocity", velocity);
        curve->addNodeVectorQuantity("acceleration", acceleration);

        std::vector<point> newpoints;
        newpoints.resize(num_controlpoints);

        for (size_t i = 0; i < num_controlpoints; i++) {
            point new_x = controlpoints[i] + velocity[i] * dt;
            newpoints[i] = new_x;
        }
        controlpoints = newpoints;

        std::vector<point> newedges;
        std::vector<float> newedge_length;
        newedges.resize(num_edges);
        newedge_length.resize(num_edges);

        for (size_t i = 0; i < num_controlpoints - 1; i++) {
            newedges[i] = newpoints[i + 1] - newpoints[i];
            newedge_length[i] = glm::length(newedges[i]);
        }

        curveEdges = newedges;
        edge_length = newedge_length;
        curve->updateNodePositions(controlpoints);

        // std::vector<point> newtangent_on_edges;
        // newtangent_on_edges.resize(num_edges);
        // for (size_t i = 0; i < num_edges; i++) {
        //     newtangent_on_edges[i] = glm::normalize(newpoints[curveEdgesId[i][1]] - newpoints[curveEdgesId[i][0]]);
        // }
    }
    void loop() {
        cal_velocity();

        cal_tangent();
        update_bishop();
        update_material_frame();
        cal_attrs();
        calforce();
        // curve->updateNodePositions(controlpoints);
    }


  private:
    size_t num_controlpoints;
    size_t num_edges;

    float eps = 1e-6;

    float alpha = 1;
    float beta = 1;

    const point reference = point(0, 0, 1);
    float totaltwist = 10 * PI;
    float totallength;

    std::vector<point> controlpoints;
    std::vector<std::vector<size_t>> curveEdgesId;
    std::vector<point> curveEdges;
    std::vector<float> edge_length;
    CurveNetwork* curve;

    std::vector<float> arc_length;

    std::vector<point> tangent_on_edges;
    std::vector<point> normal_on_edges;
    std::vector<point> binormal_on_edges;


    std::vector<glm::mat3> parallel_transport;
    std::vector<float> twist_thetas;
    std::vector<point> material_v;
    std::vector<point> material_u;


    std::vector<float> turn_angles;
    std::vector<float> curvature;

    std::vector<point> velocity;

    std::vector<point> darboux;
    std::vector<float> vertex_weight;
    std::vector<point> bending_force;
    std::vector<point> twisting_force;

    std::vector<point> acceleration;
    std::vector<float> old_edge_length;
};