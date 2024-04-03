#include "polyscope/curve_network.h"
#include "polyscope/utilities.h"
#include <cmath>
#include <cstddef>
#include <iostream>
#include <mach/task_info.h>
#include <numeric>
#include <ostream>
#include <vector>

using namespace polyscope;
class Curve {
    using point = glm::vec3;

  public:
    Curve(){};
    ~Curve(){};

    void resize_vectors() {

        edge_length.resize(num_edges);
        arc_length.resize(num_edges);
        turn_angles.resize(num_controlpoints);
        tangent_on_edges.resize(num_edges);
        normal_on_edges.resize(num_edges);
        binormal_on_edges.resize(num_edges);
        curvature.resize(num_edges);
        darboux.resize(num_edges);

        twist_angles.resize(num_edges);
        material_u.resize(num_edges);
        material_v.resize(num_edges);

        velocity.resize(num_controlpoints);
        acceleration.resize(num_controlpoints);
        bending_force.resize(num_controlpoints);
        twisting_force.resize(num_controlpoints);
    }

    void init_straight_line() {
        for (int i = 0; i < 10; i++) {
            controlpoints.push_back({i, 1, 0});
        }
        num_controlpoints = 10;
        num_edges = 9;
        for (size_t i = 0; i < num_controlpoints - 1; i++) {
            curveEdgesId.push_back({i, i + 1});
            curveEdges.push_back(controlpoints[i + 1] - controlpoints[i]);
        }
        curve = registerCurveNetwork("Straight Line", controlpoints, curveEdgesId);

        resize_vectors();
    }


    void init_sprial() {
        int count = 0;
        for (double theta = 0; theta < 4 * M_PI; theta += M_PI / 20) {
            controlpoints.push_back({cos(theta), sin(theta), theta / 10});
            count += 1;
        }
        num_controlpoints = count;
        num_edges = count - 1;

        for (size_t i = 0; i < num_controlpoints - 1; i++) {
            curveEdgesId.push_back({i, i + 1});
            curveEdges.push_back(controlpoints[i + 1] - controlpoints[i]);
        }
        curve = registerCurveNetwork("Sprial", controlpoints, curveEdgesId);

        resize_vectors();
    }


    void animate() {
        static float increment = 0.01f;
        for (auto& pt : controlpoints) {
            pt.z += increment;
            if (pt.z < 0 || pt.z > 3) increment = -increment; // Reverse direction when reaching limits
        }
        curve->updateNodePositions(controlpoints);
    }

    void arc_length_parameterization() {
        arc_length[0] = 0;

        for (size_t i = 1; i < num_edges; i++) {
            totallength += glm::length(curveEdges[i]);
            edge_length[i] = glm::length(curveEdges[i]);
        }
        curve->addEdgeScalarQuantity("edge length", edge_length);

        for (size_t i = 1; i < num_edges; i++) {
            arc_length[i] = arc_length[i - 1] + glm::length(curveEdges[i]) / totallength;
        }
        curve->addEdgeScalarQuantity("arc length", arc_length);
    }


    void update_bishop() {

        for (size_t iE = 0; iE < num_edges; iE++) {
            size_t i0 = curveEdgesId[iE][0];
            size_t i1 = curveEdgesId[iE][1];
            tangent_on_edges[iE] = glm::normalize(controlpoints[i1] - controlpoints[i0]);
        }
        curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);

        normal_on_edges[0] = glm::normalize(glm::cross(reference, tangent_on_edges[0]));
        for (size_t i = 1; i < num_edges; i++) {
            point d_tangent = tangent_on_edges[i] - tangent_on_edges[i - 1];
            if (glm::length(d_tangent) > 1e-6) { // 非零向量
                normal_on_edges[i] = glm::normalize(d_tangent);
            } else { // 零向量，沿用前一条边的法线
                normal_on_edges[i] = normal_on_edges[i - 1];
            }
        }
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 0; i < num_edges; i++) {
            binormal_on_edges[i] = glm::normalize(glm::cross(tangent_on_edges[i], normal_on_edges[i]));
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);
    }

    void calc_turn_angle() {
        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            auto cur_edge = curveEdges[i];
            auto prev_edge = curveEdges[i - 1];
            auto angle = acos(glm::dot(glm::normalize(cur_edge), glm::normalize(prev_edge)));
            turn_angles[i] = angle;
        }
    }

    // void gen_parallel_transport() {

    //     for (size_t i = 1; i < num_controlpoints; i++) {
    //         auto prev_tan = tangent_on_edges[i - 1];
    //         point next_tan = tangent_on_edges[(i) % tangent_on_edges.size()];

    //         glm::vec3 axis = glm::cross(prev_tan, next_tan);
    //         if (glm::length(axis) > 0.0001) {
    //             axis = glm::normalize(axis);
    //             float angle = acos(glm::dot(glm::normalize(prev_tan), glm::normalize(next_tan)));
    //             turn_angles[i] = angle;
    //             glm::mat3 rotationMatrix = glm::rotate(glm::mat4(1.0f), angle, axis);
    //             parallel_transport[i] = rotationMatrix;
    //         } else {
    //             parallel_transport[i] = glm::mat3(1.0f);
    //         }
    //     }
    //     curve->addNodeScalarQuantity("turning angles", turn_angles);
    // }


    void test_run() {
        update_bishop();

        arc_length_parameterization();
        calc_turn_angle();
        gen_vertex_weight();
        cal_curvature();
        cal_darboux();

        init_twist();
        update_material_frame();

        cal_bending_force();
        cal_twisting_force();
        // symEuler();
    }


    void cal_curvature() {

        for (size_t i = 0; i < num_edges; i++) {
            curvature[i] = 2 * glm::tan(turn_angles[i] / 2);
        }
        curve->addEdgeScalarQuantity("curvature", curvature);
    }

    void gen_vertex_weight() {
        vertex_weight.resize(num_controlpoints);
        vertex_weight[0] = 0.5 * (glm::length(curveEdges[0]));
        vertex_weight[num_controlpoints - 1] = 0.5 * (glm::length(curveEdges[curveEdges.size() - 1]));
        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            vertex_weight[i] = 0.5 * (edge_length[i - 1] + edge_length[i]);
        }
        curve->addNodeScalarQuantity("vertex_weight", vertex_weight);
    }

    void cal_darboux() {

        for (size_t i = 1; i < num_edges; i++) {
            darboux[i] = curvature[i] * binormal_on_edges[i] + eps;
        }
        curve->addEdgeVectorQuantity("darboux", darboux);
    }


    void init_twist() {
        for (size_t i = 1; i < num_edges; i++) {
            twist_angles[i] = i * (totaltwist / num_edges);
        }
        curve->addEdgeScalarQuantity("twisting angle", twist_angles);
    }


    void update_material_frame() {
        for (size_t i = 1; i < num_edges; i++) {
            glm::mat3 matrix = glm::rotate(glm::mat4(1.0f), twist_angles[i], tangent_on_edges[i]);
            material_u[i] = matrix * binormal_on_edges[i];
            material_v[i] = matrix * normal_on_edges[i];
        }
        curve->addEdgeVectorQuantity("material_v", material_v);
        curve->addEdgeVectorQuantity("material_u", material_u);
    }

    void cal_bending_force() {
       
        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            auto cur_edge = curveEdges[i];
            auto prev_edge = curveEdges[i - 1];
            float co = 2 / (vertex_weight[i] *
                            (glm::length(cur_edge) * glm::length(prev_edge) + glm::dot(cur_edge, prev_edge)));
            if (i != 1) {
                bending_force[i] +=
                    2.f * glm::cross(-1.f * cur_edge, darboux[i]) + glm::dot(darboux[i], cur_edge) * darboux[i];
            }
            if (i != num_controlpoints - 1) {
                bending_force[i] +=
                    2.f * glm::cross(-1.f * prev_edge, darboux[i]) - glm::dot(darboux[i], prev_edge) * darboux[i];
            }
            bending_force[i] *= co;
        }
        curve->addNodeVectorQuantity("bending force", bending_force);
    }


    void cal_twisting_force() {
        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            auto prevd = 0.5f * darboux[i] / glm::length(curveEdges[i - 1]);
            auto nextd = -0.5f * darboux[i] / glm::length(curveEdges[i]);
            float L = 0.5 * std::accumulate(vertex_weight.begin(), vertex_weight.begin() + i, 0.0f);
            twisting_force[i] = -1 * totaltwist * (nextd + prevd) / (L + eps);
        }
        curve->addNodeVectorQuantity("twisting force", twisting_force);
    }


    void adjusttwist(float newtwist) {
        float totaltwist = newtwist;
        // cal_twisting_force();
    }

    void symEuler() {
        float dt = 0.0001;
        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            point F_total = bending_force[i] + twisting_force[i];
            acceleration[i] = F_total / vertex_weight[i];
            point new_v = velocity[i] + acceleration[i] * dt;
            velocity[i] = new_v;
        }
        curve->addNodeVectorQuantity("velocity", velocity);
        curve->addNodeVectorQuantity("acceleration", acceleration);
        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            point new_x = controlpoints[i] + velocity[i] * dt;
            controlpoints[i] = new_x;
        }
        curve->updateNodePositions(controlpoints);
    }

    void tst_gd() {
        float lr = 0.01;
        int iterations = 10;

        for (int iter = 0; iter < iterations; iter++) {
            for (size_t i = 0; i < twist_angles.size(); i++) {
                float gradient = 2 * twist_angles[i];
                twist_angles[i] -= lr * gradient;
            }
            updatematerial_frame();

            curve->addEdgeScalarQuantity("twisting angle", twist_angles);
            float loss = std::accumulate(twist_angles.begin(), twist_angles.end(), 0.0f,
                                         [](float acc, float angle) { return acc + angle * angle; });
            std::cout << "Iteration " << iter << ", Loss: " << loss << std::endl;
        }
    }

    void tst_gd_bend() {

        float v0 = 0.0;
        float lr = 0.01;
        float dt = 0.01;

        for (size_t iter = 0; iter < 1; iter++) {
            for (size_t i = 1; i < num_controlpoints - 1; i++) {
                acceleration[i] = lr * bending_force[i] / vertex_weight[i];
                point new_x = controlpoints[i] + v0 * dt + 0.5f * acceleration[i] * dt * dt;
                controlpoints[i] = new_x;
            }
            curve->updateNodePositions(controlpoints);
            update_bishop();
            gen_vertex_weight();
            update_material_frame();
            cal_curvature();
            cal_darboux();
            cal_bending_force();
        }
    }

    // void updatupdate_all() {
    //     for (size_t i = 1; i < num_edges; i++) {
    //         glm::mat3 rot = glm::rotate(glm::mat4(1.0f), twist_angles[i], tangent_on_edges[i]);
    //         material_u[i] = rot * binormal_on_edges[i];
    //         material_v[i] = rot * normal_on_edges[i];
    //     }
    //     curve->addEdgeVectorQuantity("material_v", material_v);
    //     curve->addEdgeVectorQuantity("material_u", material_u);
    // }

    void updatematerial_frame() {
        for (size_t i = 1; i < num_edges; i++) {
            glm::mat3 rot = glm::rotate(glm::mat4(1.0f), twist_angles[i], tangent_on_edges[i]);
            material_u[i] = rot * binormal_on_edges[i];
            material_v[i] = rot * normal_on_edges[i];
        }
        curve->addEdgeVectorQuantity("material_v", material_v);
        curve->addEdgeVectorQuantity("material_u", material_u);
    }


  private:
    size_t num_controlpoints;
    size_t num_edges;

    float eps = 1e-6;

    float alpha = 0.2;
    float beta = 100;

    const point reference = point(0, 0, 1);
    float totaltwist = 100 * PI;
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
    std::vector<float> twist_angles;
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
};