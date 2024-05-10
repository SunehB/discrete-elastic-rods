#include "Curve.h"
#include "polyscope/curve_network.h"
#include "polyscope/curve_network.ipp"
#include <Eigen/Sparse>
#include <cmath> // For std::cosh and std::sqrt
#include <iostream>
#include <vector>


std::ostream& operator<<(std::ostream& os, const glm::dvec3& v) {
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

using namespace glm;
using namespace polyscope;

class Catenary {
  public:
    Catenary(){};
    ~Catenary(){};

    void init() {
        point p1 = {0, 0, 6};
        point p2 = {20, 0, 6};
        double d = sqrt(dot(p2 - p1, p2 - p1)); // Correct distance calculation
        double a = 5;
        size_t nsample = 12;

        for (size_t i = 0; i < nsample; i++) {
            double t = i / (nsample - 1.0); // Remove +1e-6 to simplify
            double x = p1.x + t * (p2.x - p1.x);
            double s = (x - (p1.x + p2.x) / 2) / a; // Adjust cosh parameter
            double z = p1.z + a * cosh(s);
            point p = point(x, 0, z); // Adjust y and z coordinates correctly
            // std::cout << "initial controlpoints " << i << ": " << p << std::endl;
            controlpoints.push_back(p);
        }

        for (size_t i = 0; i < nsample - 1; i++) {
            curveEdgesId.push_back({i, i + 1});
        }
        curve = registerCurveNetwork("catenary", controlpoints, curveEdgesId);
        num_edges = nsample - 1;

        initcurve();
    }
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
        bool curveloop = false;
        num_controlpoints = controlpoints.size();
        if (curveloop) {
            nSamples = num_controlpoints;
        } else {
            nSamples = num_controlpoints - 1;
        }
        resize_vectors();

        for (size_t i = 0; i < nSamples; i++) {
            if (curveloop) {
                u_edgesid.push_back({i, (i + 1) % nSamples});
                v_edgesid.push_back({i, (i + 1) % nSamples});
                curveEdgesId.push_back({i, (i + 1) % nSamples});
                curveEdges.push_back(controlpoints[(i + 1) % nSamples] - controlpoints[i]);
            } else {
                u_edgesid.push_back({i, i + 1});
                v_edgesid.push_back({i, i + 1});
                curveEdgesId.push_back({i, i + 1});
                curveEdges.push_back(controlpoints[i + 1] - controlpoints[i]);
            }
        }
        for (size_t i = 0; i < num_edges; i++) {
            edge_length[i] = length(curveEdges[i]);
        }

        totallength = std::accumulate(edge_length.begin(), edge_length.end(), 0.0f);
        for (size_t i = 0; i < num_edges; i++) {
            arc_length[i] = std::accumulate(edge_length.begin(), edge_length.begin() + i, 0.0f);
        }

        cal_tangent();
        normal_on_edges[0] = normalize(cross(tangent_on_edges[0], reference));
        binormal_on_edges[0] = normalize(cross(tangent_on_edges[0], normal_on_edges[0]));
        for (size_t i = 1; i < num_edges; i++) {
            point t_i = tangent_on_edges[i];
            point t_i_1;
            if (curveloop) {
                t_i_1 = tangent_on_edges[(i - 1) % num_edges];
            } else {
                t_i_1 = tangent_on_edges[i - 1];
            }

            point axis = cross(t_i_1, t_i);
            double theta = std::atan2(length(axis), dot(t_i_1, t_i));
            dmat3 rotationMatrix = rotate(dmat4(1), theta, normalize(axis));
            if (theta > 1e-8) {
                normal_on_edges[i] = normalize(rotationMatrix * normal_on_edges[i - 1]);
            } else {
                normal_on_edges[i] = normal_on_edges[i - 1];
            }
            rdy = true;
        }
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 1; i < num_edges; i++) {
            binormal_on_edges[i] = normalize(cross(tangent_on_edges[i], normal_on_edges[i]));
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);

        for (size_t i = 0; i < num_edges; i++) {
            double curtwist = totaltwist * arc_length[i] / totallength;
            dmat3 rotationMatrix = rotate(dmat4(1.0f), curtwist, tangent_on_edges[i]);
            material_u[i] = normalize(rotationMatrix * normal_on_edges[i]);
            material_v[i] = normalize(rotationMatrix * binormal_on_edges[i]);
        }
        curve->addEdgeVectorQuantity("material_u", material_u);
        curve->addEdgeVectorQuantity("material_v", material_v);

        material_u.resize(num_edges);
        material_v.resize(num_edges);

      
        cal_attrss();
        cal_bending_force();
        cal_twisting_force();

        cal_velocity();
    }
    dmat3 transpoe(point v) {
        return dmat3(v.x, 0, 0, v.y, 0, 0, v.z, 0, 0);
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
        bending_force.resize(num_controlpoints, point(0, 0, 0));
        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            for (size_t j = 1; j < num_edges; j++) {
                if (abs(i - j) > 1) {
                    continue;
                }
                point edgeR = curveEdges[j];
                point edgeL = curveEdges[(j - 1) % num_edges];
                double w = 0.5f * vertex_weight[j];
                auto coeff = -2 / (w * (length(edgeR) * length(edgeL) + dot(edgeR, edgeL)));
                auto kbj = darboux[j];
                point kbjTRTkbj = dmat3vec3(transpose(outerProduct(kbj, edgeR)), kbj);
                point kbjTLTkbj = dmat3vec3(transpose(outerProduct(kbj, edgeL)), kbj);

                point p;
                if (i == j - 1) {
                    p = 2.0 * cross(edgeR, kbj) + kbjTRTkbj;
                } else if (i == j) {
                    p = -1.0 * (2.0 * cross(edgeR, kbj) + kbjTRTkbj + 2.0 * cross(edgeL, kbj) - kbjTLTkbj);
                } else if (i == j + 1) {
                    p = 2.0 * cross(edgeL, kbj) - kbjTLTkbj;
                }
                bending_force[i] += p * coeff;
            }
        }
        curve->addNodeVectorQuantity("bending force", bending_force);
    }

    void update_bishop() {
        for (size_t i = 1; i < num_edges; i++) {
            point t_i = tangent_on_edges[i];
            point t_i_1 = tangent_on_edges[(i - 1)];

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

    void loop() {
        cal_velocity();
        fastprojection();
        curve->updateNodePositions(controlpoints);
        cal_tangent();
        update_bishop();
        // // update_material_frame();
        cal_attrss();
        cal_bending_force();
        cal_twisting_force();
    }


    void cal_twisting_force() {
        twisting_force.clear();
        twisting_force.resize(num_controlpoints);

        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            auto prevd = 0.5 * darboux[i] / edge_length[i-1];
            auto nextd = -0.5 * darboux[i] / edge_length[i];
            twisting_force[i] = -1 * totaltwist * (nextd + prevd) / (0.5f * totallength);
            // std::cout << "twisting force " << twisting_force[i] << std::endl;
        }
        curve->addNodeVectorQuantity("twisting force", twisting_force);
    }

    void cal_tangent() {
        tangent_on_edges.resize(num_edges);
        for (size_t i = 0; i < num_edges; i++) {
            tangent_on_edges[i] = normalize(curveEdges[i]);
        }
        auto u = curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
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
        std::cout << "0max:" << max << std::endl;

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
        // int count = 0
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
            // std::cout << "constraintGrad good "<< std::endl;
            // 100 * 300
            Eigen::SparseMatrix<double> constraintGradT = constraintGrad.transpose();

            Eigen::SparseLU<Eigen::SparseMatrix<double>> MinvDCsolver;
            MinvDCsolver.compute(Mass);
            if (MinvDCsolver.info() != Eigen::Success) {
                std::cout << "MinvDCsolver decomposition failed" << std::endl;
            }
            // std::cout << "MinvDCsolver good "<< std::endl;

            // 300 * 100
            Eigen::SparseMatrix<double> MinvDC = MinvDCsolver.solve(constraintGradT);
            if (MinvDCsolver.info() != Eigen::Success) {
                std::cout << "MinvDCsolver solve failed" << std::endl;
            }

            // std::cout << "MinvDC good "<< std::endl;
            // 100 * 100
            Eigen::SparseMatrix<double> DCMinvDC = constraintGrad * MinvDC;
            Eigen::SparseLU<Eigen::SparseMatrix<double>> dLambdasolver;
            dLambdasolver.compute(DCMinvDC);
            if (dLambdasolver.info() != Eigen::Success) {
                std::cout << "dLambdasolver decomposition failed" << std::endl;
            }

            // std::cout << "dLambdasolver good "<< std::endl;

            for (size_t i = 0; i < num_edges; i++) {
                constraintMatrix(i, 0) = constraint[i];
            }
            // 100 *1
            Eigen::MatrixXd dLambda = dLambdasolver.solve(constraintMatrix);
            if (MinvDCsolver.info() != Eigen::Success) {
                std::cout << "dLambdasolver solve failed" << std::endl;
            }
            // std::cout << "dLambda good "<< std::endl;

            // 300 * 1
            Eigen::MatrixXd dxmatrix = -1.0 * MinvDC * dLambda;

            // reshape
            for (size_t i = 1; i < num_controlpoints - 1; i++) {
                point p = point(dxmatrix(i * 3), dxmatrix(i * 3 + 1), dxmatrix(i * 3 + 2));
                newpoints[i] += p;
                // std::cout << "dx"<< p.x << " " << p.y << " " << p.z << std::endl;
            }

            for (size_t i = 1; i < num_edges; i++) {
                newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
                // std::cout << "newedges at " << i << ": " << newedges[i] << std::endl;
                newedge_length[i] = length(newedges[i]);
            }
            for (size_t i = 0; i < num_edges; i++) {
                constraint[i] = pow(newedge_length[i], 2) - pow(edge_length[i], 2);
                // std::cout << "constraint at " << i << ": " << constraint[i] << std::endl;
            }
            std::transform(constraint.begin(), constraint.end(), abs_con.begin(),
                           [](double value) { return std::abs(value); });
            max = *std::max_element(abs_con.begin(), abs_con.end());
            std::cout << "max:" << max << std::endl;
            break;
        }

        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            controlpoints[i] = newpoints[i];
        }

        curveEdges = newedges;
        edge_length = newedge_length;

        totallength = std::accumulate(edge_length.begin(), edge_length.end(), 0.0f);
        for (size_t i = 0; i < nSamples; i++) {
            arc_length[i] = std::accumulate(edge_length.begin(), edge_length.begin() + i, 0.0f);
            // std::cout << "arc_length " << i << ": " << arc_length[i] << std::endl;
        }

        for (size_t i = 0; i < num_edges; i++) {
            velocity[i] = (newpoints[i] - controlpoints[i]) / dt;
        }
    }

    void cal_velocity() {
        acceleration.clear();
        acceleration.resize(num_controlpoints);

        newpoints.resize(num_controlpoints);
        newpoints[0] = controlpoints[0];
        newpoints[num_controlpoints - 1] = controlpoints[num_controlpoints - 1];

        for (size_t i = 1; i < num_controlpoints - 1; i++) {
            acceleration[i] = (bending_force[i] + twisting_force[i]) / vertex_weight[i];
            velocity[i] += acceleration[i] * dt;
            newpoints[i] = controlpoints[i] + velocity[i] * dt;
        }
        curve->addNodeVectorQuantity("velocity", velocity);
        curve->addNodeVectorQuantity("acceleration", acceleration);
        newedges.resize(num_edges);
        newedge_length.resize(num_edges);

        for (size_t i = 0; i < num_controlpoints - 1; i++) {
            newedges[i] = newpoints[i + 1] - newpoints[i];
            newedge_length[i] = length(newedges[i]);
        }
    }

    void cal_attrss() {
        for (size_t i = 0; i < num_controlpoints; i++) {
            point Ledge;
            point Redge;
            if (i == 0) {
                Ledge = point(0, 0, 0);
            } else {
                Ledge = curveEdges[i - 1];
            }
            if (i == num_controlpoints - 1) {
                Redge = point(0, 0, 0);

            } else {
                Redge = curveEdges[i];
            }
            if (i == 0) {
                vertex_weight[i] = length(Redge);
                darboux[i] = 2.0 * cross(Redge, Redge) / (length(Redge) * length(Redge) + dot(Redge, Redge));
            } else if (i == num_controlpoints - 1) {
                vertex_weight[i] = length(Ledge);
                darboux[i] = 2.0 * cross(Ledge, Ledge) / (length(Ledge) * length(Ledge) + dot(Ledge, Ledge));
            } else {
                vertex_weight[i] = 0.5 * (length(Ledge) + length(Redge));
                darboux[i] = 2.0 * cross(Ledge, Redge) / (length(Ledge) * length(Redge) + dot(Ledge, Redge));
            }

            // darboux[i] = 2.0 * cross(Ledge, Redge) / (length(Ledge) * length(Redge) + dot(Ledge, Redge));
            // std::cout << "darboux at " << i << ": " << darboux[i] << std::endl;
            curvature[i] = sqrt(dot(darboux[i], darboux[i]));
            // std::cout << "curvature at " << i << ": " << curvature[i] << std::endl;
        }
        curve->addNodeScalarQuantity("vertex_weight", vertex_weight);
        curve->addNodeScalarQuantity("curvature", curvature);
        curve->addNodeVectorQuantity("darboux", darboux);
    }

  private:
    bool rdy = false;
    size_t num_controlpoints;
    size_t num_edges;

    double eps = 1e-6;

    double alpha = 1;
    double beta = 1;
    double dt = 1e-3;
    size_t nSamples = 100;

    const point reference = point(0, 0, 1);
    double totaltwist = 0;
    double totallength;

    std::vector<point> controlpoints;
    std::vector<std::vector<size_t>> curveEdgesId;
    std::vector<point> curveEdges;
    std::vector<double> edge_length;
    CurveNetwork* curve;
    CurveNetwork* ucurve;
    CurveNetwork* vcurve;

    std::vector<double> arc_length;

    std::vector<point> tangent_on_edges;
    std::vector<point> normal_on_edges;
    std::vector<point> binormal_on_edges;


    std::vector<dmat3> parallel_transport;
    std::vector<double> twist_thetas;
    std::vector<point> material_v;
    std::vector<point> material_u;

    std::vector<point> u_points;
    std::vector<point> v_points;
    std::vector<std::vector<size_t>> u_edgesid;
    std::vector<std::vector<size_t>> v_edgesid;


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
    std::vector<bool> isfixed;
};
