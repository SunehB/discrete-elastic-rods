#include "glm/glm.hpp"
#include "polyscope/curve_network.h"
#include <Eigen/Sparse>
#include <algorithm>
#include <iostream>
#include <vector>

using namespace polyscope;
using namespace glm;

std::ostream& operator<<(std::ostream& os, const glm::dvec3& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

class Curven {
    using point = dvec3;

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
        dt = 1.0 / nSamples;

        for (size_t i = 0; i < nSamples; ++i) {
            double t = i * dt;
            point p;
            p.x = cos(2 * M_PI * t);
            p.y = sin(2 * M_PI * t);
            p.z = 0.3 * sin(4 * M_PI * t);
            controlpoints.push_back(p);
            // std::cout << "controlpoints " << i << ": " << p << std::endl;
        }

        num_controlpoints = controlpoints.size();
        num_edges = nSamples;
        resize_vectors();

        for (size_t i = 0; i < nSamples; i++) {
            curveEdgesId.push_back({i, (i + 1) % nSamples});
            curveEdges.push_back(controlpoints[(i + 1) % nSamples] - controlpoints[i]);
            edge_length[i] = length(curveEdges[i]);
            // std::cout << "edge_length " << i << ": " << edge_length[i] << std::endl;
        }
        totallength = std::accumulate(edge_length.begin(), edge_length.end(), 0.0f);

        curve = registerCurveNetwork("Spiral", controlpoints, curveEdgesId);

        cal_tangent();
        normal_on_edges[0] = normalize(cross(tangent_on_edges[0], reference));
        binormal_on_edges[0] = normalize(cross(tangent_on_edges[0], normal_on_edges[0]));
        // std::cout << "normal on edges 0: " << normal_on_edges[0] << std::endl;
        // std::cout << "binormal on edges 0: " << binormal_on_edges[0] << std::endl;
        for (size_t i = 1; i < num_edges; i++) {
            point t_i = tangent_on_edges[i];
            point t_i_1 = tangent_on_edges[(i - 1) % num_edges];

            point axis = cross(t_i_1, t_i);
            // std::cout << "axis " << i << ": " << axis << std::endl;
            double theta = std::atan2(length(axis), dot(t_i_1, t_i));
            dmat3 rotationMatrix = rotate(dmat4(1), theta, normalize(axis));
            // for (size_t j = 0; j < 3; j++) {
            //     for (size_t k = 0; k < 3; k++) {
            //         std::cout << rotationMatrix[j][k] << " ";
            //     }
            // }
            // std::cout << "rotationMatrix "  << std::endl;
            // std::cout << "theta " << i << ": " << theta << std::endl;
            if (theta > 1e-10) {
                normal_on_edges[i] = normalize(rotationMatrix * normal_on_edges[i - 1]);

            } else {
                normal_on_edges[i] = normal_on_edges[i - 1];
            }
            // std::cout << "normal on edges " << i << ": " << normal_on_edges[i] << std::endl;
        }
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 1; i < num_edges; i++) {
            binormal_on_edges[i] = normalize(cross(tangent_on_edges[i], normal_on_edges[i]));
            // std::cout << "binormal on edges " << i << ": " << binormal_on_edges[i] << std::endl;
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);


        // update_bishop();
        update_material_frame();
        cal_attrs();
        calforce();
        // cal_velocity();
    }

    void cal_tangent() {
        for (size_t i = 0; i < num_edges; i++) {
            tangent_on_edges[i] = normalize(curveEdges[i]);
            // std::cout << "tangent on edges " << i << ": " << tangent_on_edges[i] << std::endl;
        }
        auto u = curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
        // std::cout << u << std::endl;
    }

    void update_bishop() {
        for (size_t i = 0; i < num_edges; i++) {
            point t_i = tangent_on_edges[i];
            point t_i_1 = tangent_on_edges[(i - 1) % num_edges];

            point axis = cross(t_i_1, t_i);
            double theta = std::atan2(length(axis), dot(t_i_1, t_i));
            dmat3 rotationMatrix = rotate(dmat4(1), theta, axis);

            if (theta > 1e-10) {
                normal_on_edges[i] = rotationMatrix * normal_on_edges[(i - 1) % num_edges];
            } else {
                normal_on_edges[i] = normal_on_edges[(i - 1) % num_edges];
            }
        }
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 0; i < num_edges; i++) {
            binormal_on_edges[i] = normalize(cross(tangent_on_edges[i], normal_on_edges[i]));
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);
    }

    void update_material_frame() {
        for (size_t i = 0; i < num_edges; i++) {
            dmat3 rotationMatrix = rotate(dmat4(1.0f), totaltwist, tangent_on_edges[i]);
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

            vertex_weight[i] = 0.5 * (length(Ledge) + length(Redge));

            darboux[i] = 2.0 * cross(Ledge, Redge) / (length(Ledge) * length(Redge) + dot(Ledge, Redge));
            // std::cout << "darboux at " << i << ": " << darboux[i] << std::endl;
            curvature[i] = sqrt(dot(darboux[i], darboux[i]));
            // std::cout << "curvature at " << i << ": " << curvature[i] << std::endl;
        }
        curve->addNodeScalarQuantity("vertex_weight", vertex_weight);
        curve->addNodeScalarQuantity("curvature", curvature);
        curve->addNodeVectorQuantity("darboux", darboux);
    }

    dmat3 outerProduct(const vec3& A, const vec3& B) {
        dmat3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result[i][j] = A[i] * B[j];
            }
        }
        return result;
    }
    point dmat3vec3(const dmat3& A, const vec3& B) {
        point result;
        for (int i = 0; i < 3; ++i) {
            result[i] = 0;
            for (int j = 0; j < 3; ++j) {
                result[i] += A[i][j] * B[j];
            }
        }
        return result;
    }

    void cal_bending_force() {
        bending_force.clear();
        bending_force.resize(num_controlpoints);

        for (size_t i = 0; i < num_controlpoints; i++) {
            for (size_t j = 0; j < num_controlpoints; j++) {
                point edgeR = curveEdges[j];
                point edgeL = curveEdges[(j - 1 + num_controlpoints) % num_controlpoints];

                dmat3 Rmatrix = transpoe(edgeR);
                dmat3 Lmatrix = transpoe(edgeL);

                dmat3 RT = transpose(Rmatrix);
                dmat3 LT = transpose(Lmatrix);

                double w = 0.5f * vertex_weight[j];
                auto coeff = -2 / (w * (length(edgeR) * length(edgeL) + dot(edgeR, edgeL)));
                // std::cout << "coeff " << coeff << std::endl;
                auto kbj = darboux[j];
                // dmat3 kbjT = transpose(transpoe(kbj));
                dmat3 kbjTR = outerProduct(kbj, edgeR);
                dmat3 kbjTRT = transpose(kbjTR);

                dmat3 kbjTL = outerProduct(kbj, edgeL);
                dmat3 kbjTLT = transpose(kbjTL);
                // for (int i = 0; i < 3; ++i) {
                //     for (int j = 0; j < 3; ++j) {
                //         std::cout << kbjTRT[i][j] << " ";
                //     }
                // }
                // std::cout << std::endl;
                point kbjTRTkbj = dmat3vec3(kbjTRT, kbj);
                point kbjTLTkbj = dmat3vec3(kbjTLT, kbj);
                // std::cout << " kbj " << kbjTRTkbj << std::endl;
                // std::cout << " kbj " << kbjTLTkbj << std::endl;
                point p;
                if (i == j - 1) {
                    p = 2.0 * cross(edgeR, kbj) + kbjTRTkbj;
                } else if (i == j) {
                    p = -1.0 * (2.0 * cross(edgeR, kbj) + kbjTRTkbj + 2.0 * cross(edgeL, kbj) - kbjTLTkbj);
                } else if (i == j + 1) {
                    p = 2.0 * cross(edgeL, kbj) - kbjTLTkbj;
                }
                bending_force[i] += p * coeff;
                // std::cout << "bending force " << p << std::endl;
            }
            bending_force[i] *= -1;
            // std::cout << "bending force " << bending_force[i] << std::endl;
        }
        // std::cout << "bending force " << bending_force[0] << std::endl;
        curve->addNodeVectorQuantity("bending force", bending_force);

        // auto maxlist = std::transform(bending_force.begin(), bending_force.end(), [](point a, point b) {
        //     return length(a) < length(b);
        // });
        // auto max  = *std::max_element(bending_force.begin(), bending_force.end(), [](point a, point b) {
        //     return length(a) > length(b);
        // });
        // std::cout << "max bending force: " << max << std::endl;
        // std::cout << "bending force done" << std::endl;
    }

    dmat3 transpoe(point v) {
        return dmat3(v.x, 0, 0, v.y, 0, 0, v.z, 0, 0);
    }

    void cal_twisting_force() {
        twisting_force.clear();
        twisting_force.resize(num_controlpoints);

        for (size_t i = 0; i < num_controlpoints; i++) {
            auto prevd = 0.5 * darboux[i] / edge_length[(i - 1 + num_controlpoints) % num_controlpoints];
            auto nextd = -0.5 * darboux[i] / edge_length[i];
            twisting_force[i] = -1 * totaltwist * (nextd + prevd + eps) / (0.5f * totallength);
            // std::cout << "twisting force " << twisting_force[i] << std::endl;
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

        // for (size_t i = 0; i < num_controlpoints; i++) {
        //     // acceleration[i] =  (bending_force[i] + twisting_force[i]) / vertex_weight[i];
        //     // velocity[i] += acceleration[i] * dt;
        //     // std::cout << "velocity " << i << ": " << velocity[i] << std::endl;
        //     // std::cout << "acceleration " << i << ": " << acceleration[i] << std::endl;
        // }


        for (size_t i = 0; i < num_controlpoints; i++) {
            // std::cout << "velocity " << i << ": " << acceleration[i] * dt << std::endl;
            acceleration[i] = (bending_force[i] + twisting_force[i]) / vertex_weight[i];
            // velocity[i] = acceleration[i] / 1000;
            // std::cout << acceleration[i].x /1000  << std::endl;
            // std::cout << acceleration[i].y * dt  << std::endl;
            // std::cout << acceleration[i].z * dt  << std::endl;
            point v = point(acceleration[i].x / 1000, acceleration[i].y / 1000, acceleration[i].z / 1000);
            velocity[i] += v;
            // std::cout << "velocity " << i << ": " << acceleration[i] * dt << std::endl;
            // std::cout << "velocity2 " << i << ": " << velocity[i] << std::endl;
            // std::cout << "acceleration " << i << ": " << acceleration[i] << std::endl;
        }
        curve->addNodeVectorQuantity("velocity", velocity);
        curve->addNodeVectorQuantity("acceleration", acceleration);

        // std::cout << "velocity " << velocity[0]*dt << std::endl;

        newpoints.resize(num_controlpoints);

        for (size_t i = 0; i < num_controlpoints; i++) {
            // newpoints[i] = controlpoints[i] + velocity[i] * dt;
            point p = point(velocity[i].x / 1000, velocity[i].y / 1000, velocity[i].z / 1000);
            newpoints[i] = controlpoints[i] + p;
            // std::cout << "newpoints at " << i << ": " << newpoints[i] << std::endl;
        }
        // controlpoints = newpoints;

        newedges.resize(num_edges);
        newedge_length.resize(num_edges);

        for (size_t i = 0; i < num_controlpoints; i++) {
            newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
            // std::cout << "newedges at " << i << ": " << newedges[i] << std::endl;
            newedge_length[i] = length(newedges[i]);
        }

        // curveEdges = newedges;
        // edge_length = newedge_length;
    }

    void fastprojection() {
        // 100 * 1
        constraint.clear();
        constraint.resize(num_edges);

        for (size_t i = 0; i < num_edges; i++) {
            constraint[i] = pow(newedge_length[i], 2) - pow(edge_length[i], 2);
            // std::cout << "constraint at " << i << ": " << constraint[i] << std::endl;
        }

        std::vector<double> abs_con(constraint.size());
        std::transform(constraint.begin(), constraint.end(), abs_con.begin(),
                       [](double value) { return std::abs(value); });
        auto max = *std::max_element(abs_con.begin(), abs_con.end());
        // std::cout << "max:" << max << std::endl;

        // 100 * 1
        Eigen::MatrixXd constraintMatrix(num_edges, 1);
        for (size_t i = 0; i < num_edges; i++) {
            constraintMatrix(i, 0) = constraint[i];
        }

        Eigen::SparseMatrix<double> Mass(3 * nSamples, 3 * nSamples);
        std::vector<Eigen::Triplet<double>> M_tripletList;
        for (size_t i = 0; i < num_edges; i++) {
            auto w = vertex_weight[i];
            M_tripletList.push_back(Eigen::Triplet<double>(i * 3, i * 3, w));
            M_tripletList.push_back(Eigen::Triplet<double>(i * 3 + 1, i * 3 + 1, w));
            M_tripletList.push_back(Eigen::Triplet<double>(i * 3 + 2, i * 3 + 2, w));
        }
        Mass.setFromTriplets(M_tripletList.begin(), M_tripletList.end());
        // int count = 0;
        while (max > 1e-10) {
            // 100 * 300
            Eigen::SparseMatrix<double> constraintGrad(nSamples, 3 * nSamples);
            std::vector<Eigen::Triplet<double>> tripletList;
            for (size_t i = 0; i < num_edges; i++) {
                auto edge = newedges[i];
                // std::cout << "edge" << edge << std::endl;
                edge *= 2;
                tripletList.push_back(Eigen::Triplet<double>(i, i * 3, -edge.x));
                tripletList.push_back(Eigen::Triplet<double>(i, i * 3 + 1, -edge.y));
                tripletList.push_back(Eigen::Triplet<double>(i, i * 3 + 2, -edge.z));
                tripletList.push_back(Eigen::Triplet<double>(i, (i * 3 + 3) % (3 * num_edges), edge.x));
                tripletList.push_back(Eigen::Triplet<double>(i, (i * 3 + 4) % (3 * num_edges), edge.y));
                tripletList.push_back(Eigen::Triplet<double>(i, (i * 3 + 5) % (3 * num_edges), edge.z));
            }
            constraintGrad.setFromTriplets(tripletList.begin(), tripletList.end());
            // 100 * 300
            Eigen::SparseMatrix<double> constraintGradT = constraintGrad.transpose();

            Eigen::SparseLU<Eigen::SparseMatrix<double>> MinvDCsolver;
            MinvDCsolver.compute(Mass);
            if (MinvDCsolver.info() != Eigen::Success) {
                std::cout << "MinvDCsolver decomposition failed" << std::endl;
            }

            // 300 * 100
            Eigen::SparseMatrix<double> MinvDC = MinvDCsolver.solve(constraintGradT);
            if (MinvDCsolver.info() != Eigen::Success) {
                std::cout << "MinvDCsolver solve failed" << std::endl;
            }
            // 100 * 100
            Eigen::SparseMatrix<double> DCMinvDC = constraintGrad * MinvDC;
            Eigen::SparseLU<Eigen::SparseMatrix<double>> dLambdasolver;
            dLambdasolver.compute(DCMinvDC);
            if (dLambdasolver.info() != Eigen::Success) {
                std::cout << "dLambdasolver decomposition failed" << std::endl;
            }

            for (size_t i = 0; i < num_edges; i++) {
                constraintMatrix(i, 0) = constraint[i];
            }
            // 100 *1
            Eigen::MatrixXd dLambda = dLambdasolver.solve(constraintMatrix);
            if (MinvDCsolver.info() != Eigen::Success) {
                std::cout << "dLambdasolver solve failed" << std::endl;
            }

            // 300 * 1
            Eigen::MatrixXd dxmatrix = -1.0 * MinvDC * dLambda;

            // reshape
            // std::vector<point> dx;
            for (size_t i = 0; i < num_controlpoints; i++) {
                point p = point(dxmatrix(i * 3), dxmatrix(i * 3 + 1), dxmatrix(i * 3 + 2));
                // dx.push_back(p);
                newpoints[i] += p;
                // std::cout << "dx"<< p << std::endl;

                // std::cout << "dx at "<< i  << ": "<< p << std::endl;
            }

            for (size_t i = 0; i < num_controlpoints; i++) {
                newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
                // std::cout << "newedges at " << i << ": " << newedges[i] << std::endl;
                newedge_length[i] = length(newedges[i]);
            }
            for (size_t i = 0; i < num_controlpoints; i++) {
                constraint[i] = pow(newedge_length[i], 2) - pow(edge_length[i], 2);
                // std::cout << "constraint at " << i << ": " << constraint[i] << std::endl;
            }
            std::transform(constraint.begin(), constraint.end(), abs_con.begin(),
                           [](double value) { return std::abs(value); });
            max = *std::max_element(abs_con.begin(), abs_con.end());
            // double sum = std::accumulate(abs_con.begin(), abs_con.end(), 0.0f);
            // std::cout << "sum:" << sum << std::endl;
            // std::cout << "max:" << max << std::endl;
            // break;
            // count += 1;
            // if (count > 1) {
            //     break;
            // }
        }


        controlpoints = newpoints;
        curveEdges = newedges;
        edge_length = newedge_length;

        for (size_t i = 0; i < num_edges; i++) {
            velocity[i] = (newpoints[i] - controlpoints[i]) / dt;
        }
    }

    void loop() {

        // for (size_t i = 0; i < 4; i++) {
        cal_velocity();
        fastprojection();
        curve->updateNodePositions(controlpoints);
        // std::this_thread::sleep_for(std::chrono::milliseconds(1000));

        cal_tangent();
        update_bishop();
        update_material_frame();
        cal_attrs();
        calforce();
        // std::cout << "iter " << i << "done" <<std::endl;

        // }
        // std::cout << "done" << std::endl;
    }

    void set_alpha(double a) { alpha = a; }
    void set_beta(double b) { beta = b; }
    void set_dt(double d) { dt = d; }
    void set_nSamples(size_t n) { nSamples = n; }
    void set_totaltwist(double t) { totaltwist = t; }

    void reset() {
        controlpoints.clear();
        curveEdgesId.clear();
        curveEdges.clear();
        edge_length.clear();
        arc_length.clear();
        tangent_on_edges.clear();
        normal_on_edges.clear();
        binormal_on_edges.clear();
        parallel_transport.clear();
        twist_thetas.clear();
        material_v.clear();
        material_u.clear();
        turn_angles.clear();
        curvature.clear();
        velocity.clear();
        darboux.clear();
        vertex_weight.clear();
        bending_force.clear();
        twisting_force.clear();
        acceleration.clear();
        old_edge_length.clear();
        newpoints.clear();
        newedges.clear();
        newedge_length.clear();
        constraint.clear();

        initcurve();
    }


  private:
    size_t num_controlpoints;
    size_t num_edges;

    double eps = 1e-6;

    double alpha = 1;
    double beta = 1;
    double dt = 1e-3;
    size_t nSamples = 100;

    const point reference = point(0, 0, 1);
    double totaltwist = PI;
    double totallength;

    std::vector<point> controlpoints;
    std::vector<std::vector<size_t>> curveEdgesId;
    std::vector<point> curveEdges;
    std::vector<double> edge_length;
    CurveNetwork* curve;

    std::vector<double> arc_length;

    std::vector<point> tangent_on_edges;
    std::vector<point> normal_on_edges;
    std::vector<point> binormal_on_edges;


    std::vector<dmat3> parallel_transport;
    std::vector<double> twist_thetas;
    std::vector<point> material_v;
    std::vector<point> material_u;


    std::vector<double> turn_angles;
    std::vector<double> curvature;
    std::vector<point> velocity;


    std::vector<point> darboux;
    std::vector<double> vertex_weight;
    std::vector<point> bending_force;
    std::vector<point> twisting_force;

    std::vector<point> acceleration;
    std::vector<double> old_edge_length;

    std::vector<point> newpoints;
    std::vector<point> newedges;
    std::vector<double> newedge_length;

    std::vector<double> constraint;
};