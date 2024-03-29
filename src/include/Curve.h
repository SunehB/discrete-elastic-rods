#include "glm/geometric.hpp"
#include "glm/glm.hpp"
#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include <cstddef>
#include <ostream>
#include <vector>
#include <Eigen/Sparse>

using namespace polyscope;
class Curve {
    using point = glm::vec3;

  public:
    Curve(){};
    ~Curve(){};

    void init() {
        for (double theta = 0; theta < 4 * M_PI; theta += M_PI / 20) {
            controlpoints.push_back({cos(theta), sin(theta), theta / 10});
        }
        for (size_t i = 0; i < controlpoints.size() - 1; i++) {
            curveEdges.push_back({i, i + 1});
        }

        curve = registerCurveNetwork("Sprial", controlpoints, curveEdges);

        gen_bishop();   
        init_twist();
        gen_frenet();

        cal_trun_angles();
        cal_curvature();

        cal_darboux();
        cal_voroni_region_length();
        cal_bending_energy();
        cal_twisting_energy();
    }

    void animate() {
        static float increment = 0.01f;
        for (auto& pt : controlpoints) {
            pt.z += increment;
            if (pt.z < 0 || pt.z > 3) increment = -increment; // Reverse direction when reaching limits
        }
        curve->updateNodePositions(controlpoints);

        
    }

    void gen_bishop() {
        tangent_on_edges.resize(curveEdges.size());
        normal_on_edges.resize(curveEdges.size());
        binormal_on_edges.resize(curveEdges.size());

        for (size_t iE = 0; iE < curveEdges.size(); iE++) {
            size_t i0 = curveEdges[iE][0];
            size_t i1 = curveEdges[iE][1];
            tangent_on_edges[iE] = glm::normalize(controlpoints[i1] - controlpoints[i0]);
        }
        curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
        for (size_t i = 1; i < curveEdges.size(); i++) {
            point d_tangent = tangent_on_edges[i] - tangent_on_edges[i - 1];
            normal_on_edges[i] = glm::normalize(d_tangent);
        }
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 0; i < curveEdges.size(); i++) {
            binormal_on_edges[i] = glm::cross(tangent_on_edges[i], normal_on_edges[i]);
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);
    }

    void init_twist() {
        twist_angles.clear();
        double step = (4 * M_PI) / curveEdges.size();
        for (size_t i = 0; i < curveEdges.size(); ++i) {
            twist_angles.push_back(i * step);
        }
        curve->addEdgeScalarQuantity("twist angles", twist_angles);
    }

    void gen_frenet() {
        rotated_normal_on_edges.resize(curveEdges.size());
        rotated_binormal_on_edges.resize(curveEdges.size());
        for (size_t i = 0; i < curveEdges.size(); i++) {
            float theta = twist_angles[i];
            rotated_normal_on_edges[i] =
                glm::normalize(normal_on_edges[i] * cos(theta) + binormal_on_edges[i] * sin(theta));
            rotated_binormal_on_edges[i] =
                glm::normalize(-normal_on_edges[i] * sin(theta) + binormal_on_edges[i] * cos(theta));
        }
        curve->addEdgeVectorQuantity("rotated normal on edges", rotated_normal_on_edges);
        curve->addEdgeVectorQuantity("rotated binormal on edges", rotated_binormal_on_edges);
    }

    void cal_trun_angles() {
        turn_angles.resize(curveEdges.size());
        for (size_t i = 1; i < curveEdges.size(); i++) {
            turn_angles[i] = glm::acos(glm::dot(tangent_on_edges[i], tangent_on_edges[i - 1]));
        }
        curve->addEdgeScalarQuantity("turn angles", turn_angles);
    }

    void cal_curvature() {
        curvature.resize(curveEdges.size());
        for (size_t i = 0; i < curveEdges.size(); i++) {
            curvature[i] = 2 * glm::tan(turn_angles[i] / 2);
        }
        curve->addEdgeScalarQuantity("curvature", curvature);
    }

    void cal_voroni_region_length() {
        voroni_region_length.resize(curveEdges.size());
        voroni_region_length[0] = 0;
        for (size_t i = 1; i < curveEdges.size()-1; i++) {
            voroni_region_length[i] =
                0.5 * (glm::length(controlpoints[curveEdges[i][0]] - controlpoints[curveEdges[i][1]]) +
                       glm::length(controlpoints[curveEdges[i - 1][0]] - controlpoints[curveEdges[i - 1][1]]));
        }
        curve->addEdgeScalarQuantity("voroni region length", voroni_region_length);
    }

    void cal_darboux() {
        darboux.resize(curveEdges.size());
        for (size_t i = 1; i < curveEdges.size(); i++) {
            auto edge = controlpoints[curveEdges[i][1]] - controlpoints[curveEdges[i - 1][0]];
            auto egde_prev = controlpoints[curveEdges[i - 1][1]] - controlpoints[curveEdges[i - 1][0]];
            auto deo = 2.f * glm::cross(edge, egde_prev);
            float base = glm::dot(edge, egde_prev) + glm::length(edge) * glm::length(egde_prev);
            darboux[i] = deo / base;
        }
        curve->addEdgeVectorQuantity("darboux", darboux);
    }


    void cal_bending_energy() {
        bending_energy.resize(curveEdges.size());

        for (size_t i = 1; i < curveEdges.size(); i++) {
            bending_energy[i] = bending_energy[i - 1] + 0.5 * glm::length(darboux[i]) / voroni_region_length[i];
        }
        curve->addEdgeScalarQuantity("bending energy", bending_energy);
    }

    void cal_twisting_energy() {
        twisting_energy.resize(curveEdges.size());
        for (size_t i = 1; i < curveEdges.size(); i++) {
            auto diff_twist = twist_angles[i] - twist_angles[i - 1];
            twisting_energy[i] = twisting_energy[i - 1] + glm::pow(diff_twist,2) / voroni_region_length[i];
        }
        curve->addEdgeScalarQuantity("twisting energy", twisting_energy);
    }

    void init_M() {
        size_t n = 3 * curveEdges.size() +2 ;
        M.resize(n, n);

        std::vector<Eigen::Triplet<double>> tripletList_M ;
        for (size_t i = 0; i < curveEdges.size(); i++) {
            tripletList_M.push_back(Eigen::Triplet<double>(3 * i, 3 * i, 1));
            tripletList_M.push_back(Eigen::Triplet<double>(3 * i, 3 * i + 1, 0));
            tripletList_M.push_back(Eigen::Triplet<double>(3 * i, 3 * i + 2, 0));

            tripletList_M.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i, 0));
            tripletList_M.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 1, 1));
            tripletList_M.push_back(Eigen::Triplet<double>(3 * i + 1, 3 * i + 2, 0));

            tripletList_M.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i, 0));
            tripletList_M.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 1, 0));
            tripletList_M.push_back(Eigen::Triplet<double>(3 * i + 2, 3 * i + 2, 1));
        }
        tripletList_M.push_back(Eigen::Triplet<double>(3 * curveEdges.size(), 3 * curveEdges.size(), 1));
    }

    void solve_x () {
        Eigen::VectorXd x;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(M);
        x = solver.solve(Eigen::VectorXd::Zero(M.rows()));
    }

    void update() {
        curve->updateNodePositions(controlpoints);
    }

    void cal_nabla_psi() {
        nabla_psi.resize(curveEdges.size());
        for (size_t i = 1; i < curveEdges.size(); i++) {
            auto prev_edge_length = glm::length(controlpoints[curveEdges[i - 1][0]] - controlpoints[curveEdges[i - 1][1]]);
            auto next_egde_length = glm::length(controlpoints[curveEdges[i][0]] - controlpoints[curveEdges[i][1]]);
            auto prev_nabla_psi = 0.5f *darboux[i - 1] / prev_edge_length;
            auto next_nabla_psi = -0.5f *darboux[i] / next_egde_length;

            auto cur_nabla_psi = -(prev_nabla_psi + next_nabla_psi);
        }
        curve->addEdgeVectorQuantity("nabla psi", nabla_psi);
    }

    void cal_nabla_Psi(){
        nabla_Psi.resize(curveEdges.size());
        nabla_Psi[0] = nabla_psi[0];
        for (size_t i = 1; i < curveEdges.size(); i++) {
            nabla_Psi[i] += nabla_Psi[i - 1] + nabla_psi[i];
        }
        curve->addEdgeVectorQuantity("nabla Psi", nabla_Psi);
    }

  private:
    std::vector<point> controlpoints;
    std::vector<std::vector<size_t>> curveEdges;
    CurveNetwork* curve;

    std::vector<point> tangent_on_edges;
    std::vector<point> normal_on_edges;
    std::vector<point> binormal_on_edges;

    std::vector<double> twist_angles;
    std::vector<point> rotated_normal_on_edges;
    std::vector<point> rotated_binormal_on_edges;


    std::vector<float> turn_angles;
    std::vector<double> curvature;

    std::vector<point> darboux;
    std::vector<double> voroni_region_length;
    std::vector<double> bending_energy;
    std::vector<double> twisting_energy;



    std::vector<point> nabla_psi;
    std::vector<point> nabla_Psi;



    Eigen::SparseMatrix<double> M;
};