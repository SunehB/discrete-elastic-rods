#include "glm/exponential.hpp"
#include "glm/fwd.hpp"
#include "glm/geometric.hpp"
#include "glm/matrix.hpp"
#include "polyscope/curve_network.h"
#include "polyscope/utilities.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
// #include <iostream>
#include <iostream>
#include <mach/task_info.h>
#include <numeric>
// #include <ostream>
#include <Eigen/Sparse>
#include <vector>

// #include <iostream>
#include <thread>
#include <chrono>

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
        size_t nSamples = 100;
        double dt = 1.0 / nSamples;

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
        totallength = std::accumulate(edge_length.begin(), edge_length.end(), 0.0f);

        curve = registerCurveNetwork("Spiral", controlpoints, curveEdgesId);

        cal_tangent();
        normal_on_edges[0] = glm::normalize(glm::cross(reference, tangent_on_edges[0]));
        binormal_on_edges[0] = glm::normalize(glm::cross(tangent_on_edges[0], normal_on_edges[0]));
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

        for (size_t i = 1; i < num_edges; i++) {
            binormal_on_edges[i] = glm::normalize(glm::cross(tangent_on_edges[i], normal_on_edges[i]));
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);


        // update_bishop();
        update_material_frame();
        cal_attrs();
        calforce();
    }

    void cal_tangent() {
        for (size_t i = 0; i < num_edges; i++) {
            tangent_on_edges[i] = glm::normalize(curveEdges[i]);
        }
        auto u = curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
        // std::cout << u << std::endl;
    }

    void update_bishop() {
        for (size_t i = 0; i < num_edges; i++) {
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
        for (size_t i = 0; i < num_controlpoints; i++) {
            point Ledge = curveEdges[(i - 1 + num_controlpoints) % num_controlpoints];
            point Redge = curveEdges[(i) % num_controlpoints];

            vertex_weight[i] = 0.5 * (glm::length(Ledge) + glm::length(Redge));

            darboux[i] =
                2.f * glm::cross(Ledge, Redge) / (glm::length(Ledge) * glm::length(Redge) + glm::dot(Ledge, Redge));
            curvature[i] = glm::sqrt(glm::dot(darboux[i], darboux[i]));
        }
        curve->addNodeScalarQuantity("vertex_weight", vertex_weight);
        curve->addNodeScalarQuantity("curvature", curvature);
        curve->addNodeVectorQuantity("darboux", darboux);
    }

    void cal_bending_force() {
        bending_force.clear();
        bending_force.resize(num_controlpoints);

        for (size_t i = 0; i < num_controlpoints; i++) {
            for (size_t j = 0; j < num_controlpoints; j++) {
                point edgeR = curveEdges[j];
                point edgeL = curveEdges[(j - 1 + num_controlpoints) % num_controlpoints];
                auto w = vertex_weight[j];
                auto coeff = 0.25f / (w * (glm::length(edgeR) * glm::length(edgeL) + glm::dot(edgeR, edgeL)));
                // std::cout << "coeff " << coeff << std::endl;
                auto kbj = darboux[j];

                if (i == j - 1) {
                    bending_force[i] += coeff * (2.f * glm::cross(-edgeR, kbj) + glm::dot(kbj, edgeR) * kbj);
                } else if (i == j) {
                    bending_force[i] += coeff * (2.f * glm::cross(-edgeR, kbj) - glm::dot(kbj, edgeR) * kbj +
                                                 2.f * glm::cross(-edgeL, kbj) + glm::dot(kbj, edgeL) * kbj);
                } else if (i == j + 1) {
                    bending_force[i] += coeff * (2.f * glm::cross(-edgeL, kbj) + glm::dot(kbj, edgeL) * kbj);
                }
            }
        }
        // std::cout << "bending force " << bending_force[0] << std::endl;
        curve->addNodeVectorQuantity("bending force", bending_force);
    }

    void cal_twisting_force() {
        twisting_force.clear();
        twisting_force.resize(num_controlpoints);

        for (size_t i = 0; i < num_controlpoints; i++) {
            auto prevd = 0.5f * darboux[i] / edge_length[(i - 1 + num_controlpoints) % num_controlpoints];
            auto nextd = -0.5f * darboux[i] / edge_length[i];
            twisting_force[i] = -1 * totaltwist * (nextd + prevd + eps) / totallength;
        }
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

        for (size_t i = 0; i < num_controlpoints; i++) {
            point F_total = -bending_force[i] + twisting_force[i];
            acceleration[i] = 0.1f * F_total / vertex_weight[i];
            point new_v = velocity[i] + acceleration[i] * dt;
            velocity[i] = new_v;
        }
        curve->addNodeVectorQuantity("velocity", velocity);
        curve->addNodeVectorQuantity("acceleration", acceleration);

        // std::cout << "velocity " << velocity[0]*dt << std::endl;

        newpoints.resize(num_controlpoints);

        for (size_t i = 0; i < num_controlpoints; i++) {
            newpoints[i] = controlpoints[i] + velocity[i] * dt;
        }
        controlpoints = newpoints;

        newedges.resize(num_edges);
        newedge_length.resize(num_edges);

        for (size_t i = 0; i < num_controlpoints; i++) {
            newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
            newedge_length[i] = glm::length(newedges[i]);
        }

        curveEdges = newedges;
        edge_length = newedge_length;
        
    }

    void fastprojection() {
        // 100 * 1 
        constraint.clear();
        constraint.resize(num_edges);

        for (size_t i = 0; i < num_edges; i++) {
            constraint[i] = newedge_length[i] - edge_length[i];
        }

        std::vector<float> abs_con(constraint.size());
        std::transform(constraint.begin(), constraint.end(), abs_con.begin(),
                       [](int value) { return std::abs(value); });
        auto max = *std::max_element(abs_con.begin(), abs_con.end());

        // 100 * 1 
        Eigen::MatrixXf constraintMatrix(num_edges, 1);
        for (size_t i = 0; i < num_edges; i++) {
            constraintMatrix(i, 0) = constraint[i];
        }

        Eigen::SparseMatrix<float> Mass(3 * 100, 3 * 100);
        std::vector<Eigen::Triplet<float>> M_tripletList;
        for (size_t i = 0; i < num_edges; i++) {
            auto w = vertex_weight[i];
            M_tripletList.push_back(Eigen::Triplet<float>(i * 3, i * 3, w));
            M_tripletList.push_back(Eigen::Triplet<float>(i * 3 + 1, i * 3 + 1, w));
            M_tripletList.push_back(Eigen::Triplet<float>(i * 3 + 2, i * 3 + 2, w));
        }
        Mass.setFromTriplets(M_tripletList.begin(), M_tripletList.end());

        while (max > 1e-10) {
            // 300 * 100
            Eigen::SparseMatrix<float> constraintGrad(3 * 100, 100);
            std::vector<Eigen::Triplet<float>> tripletList;
            for (size_t i = 0; i < num_edges; i++) {
                auto edge = newedges[i];
                edge *= 2;
                tripletList.push_back(Eigen::Triplet<float>(i, i * 3, -edge.x));
                tripletList.push_back(Eigen::Triplet<float>(i, i * 3 + 1, -edge.y));
                tripletList.push_back(Eigen::Triplet<float>(i, i * 3 + 2, -edge.z));
                tripletList.push_back(Eigen::Triplet<float>(i, (i * 3 + 3) % (3 * num_edges), edge.x));
                tripletList.push_back(Eigen::Triplet<float>(i, (i * 3 + 4) % (3 * num_edges), edge.y));
                tripletList.push_back(Eigen::Triplet<float>(i, (i * 3 + 5) % (3 * num_edges), edge.z));
            }

            // 100 * 300
            Eigen::SparseMatrix<float> constraintGradT = constraintGrad.transpose();

            Eigen::SparseLU<Eigen::SparseMatrix<float>> MinvDCsolver;
            MinvDCsolver.compute(Mass);
            if (MinvDCsolver.info() != Eigen::Success) {
                return;
            }

            // 300 * 100 
            Eigen::SparseMatrix<float> MinvDC = MinvDCsolver.solve(constraintGradT);
            if (MinvDCsolver.info() != Eigen::Success) {
                return;
            }
            // 100 * 100 
            Eigen::SparseMatrix<float> DCMinvDC = constraintGrad * MinvDC;
            Eigen::SparseLU<Eigen::SparseMatrix<float>> dLambdasolver;
            dLambdasolver.compute(DCMinvDC);
            if (dLambdasolver.info() != Eigen::Success) {
                return;
            }
            // 100 *1
            Eigen::MatrixXf dLambda = dLambdasolver.solve(constraintMatrix);
            if (MinvDCsolver.info() != Eigen::Success) {
                return;
            }

            // 300 * 1
            auto dxmatrix = -MinvDC * dLambda;

            // reshape
            std::vector<point> dx;
            for (size_t i = 0; i < num_controlpoints; i++) {
                dx.push_back(point(dxmatrix(i * 3), dxmatrix(i * 3 + 1), dxmatrix(i * 3 + 2)));
            }
            for (size_t i = 0; i < num_controlpoints; i++) {
                newpoints[i] = newpoints[i] + dx[i];
            }
            for (size_t i = 0; i < num_controlpoints; i++) {
                newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
                newedge_length[i] = glm::length(newedges[i]);
            }
            for (size_t i = 0; i < num_controlpoints; i++) {
                constraint[i] = newedge_length[i] - edge_length[i];
            }
            std::transform(constraint.begin(), constraint.end(), abs_con.begin(),
                           [](float value) { return std::abs(value); });
            max = *std::max_element(abs_con.begin(), abs_con.end());
        }

      
        controlpoints = newpoints;
        curveEdges = newedges;
        edge_length = newedge_length;

        for (size_t i = 0; i < num_edges; i++) {
            velocity[i] = (newpoints[i] - controlpoints[i])/dt;
        }
    }

    void loop() {

        for (size_t i = 0; i < 1000; i++) {
            cal_velocity();
            fastprojection();

            cal_tangent();
            update_bishop();
            update_material_frame();
            cal_attrs();
            calforce();
            // std::cout << "iter " << i << "done" <<std::endl;
            // std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }
        std::cout << "done" << std::endl;
          curve->updateNodePositions(controlpoints);
    }


  private:
    size_t num_controlpoints;
    size_t num_edges;

    float eps = 1e-6;

    float alpha = 1;
    float beta = 1;
    float dt = 0.001;

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

    std::vector<point> newpoints;
    std::vector<point> newedges;
    std::vector<float> newedge_length;

    std::vector<float> constraint;
};