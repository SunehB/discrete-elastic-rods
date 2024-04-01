#include "polyscope/curve_network.h"
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

    void init() {
        for (double theta = 0; theta < 4 * M_PI; theta += M_PI / 20) {
            controlpoints.push_back({cos(theta), sin(theta), theta / 10});
        }
        for (size_t i = 0; i < controlpoints.size() - 1; i++) {
            curveEdgesId.push_back({i, i + 1});
            curveEdges.push_back(controlpoints[i + 1] - controlpoints[i]);
        }

        tangent_on_edges.resize(curveEdgesId.size());
        normal_on_edges.resize(curveEdgesId.size());
        binormal_on_edges.resize(curveEdgesId.size());

        parallel_transport.resize(controlpoints.size());
        turn_angles.resize(controlpoints.size());
        arc_length.resize(curveEdgesId.size());

        material_u.resize(curveEdgesId.size());
        material_v.resize(curveEdgesId.size());

        acceleration.resize(controlpoints.size());

        curve = registerCurveNetwork("Sprial", controlpoints, curveEdgesId);

        arc_length_parameterization();
        gen_bishop_0();
        gen_vertex_weight();

        gen_twist_angle();

        cal_curvature();
        cal_darboux();
        cal_bending_force();

        cal_twisting_force();

        twisting_force.resize(curveEdgesId.size());
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
        auto sum = 0.0;
        for (size_t i = 1; i < curveEdgesId.size(); i++) {
            sum += glm::length(curveEdges[i]);
        }
        for (size_t i = 1; i < curveEdgesId.size(); i++) {
            arc_length[i] = arc_length[i - 1] + glm::length(curveEdges[i]) / sum;
        }
        curve->addEdgeScalarQuantity("arc length", arc_length);
    }  


    void gen_bishop_0() {


        for (size_t iE = 0; iE < curveEdgesId.size(); iE++) {
            size_t i0 = curveEdgesId[iE][0];
            size_t i1 = curveEdgesId[iE][1];
            tangent_on_edges[iE] = glm::normalize(controlpoints[i1] - controlpoints[i0]);
        }
        curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
        for (size_t i = 1; i < curveEdgesId.size(); i++) {
            point d_tangent = tangent_on_edges[i] - tangent_on_edges[i - 1];
            normal_on_edges[i] = glm::normalize(d_tangent);
        }
        normal_on_edges[0] = glm::normalize(glm::cross(reference, tangent_on_edges[0]));
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 0; i < curveEdgesId.size(); i++) {
            binormal_on_edges[i] = glm::normalize(glm::cross(tangent_on_edges[i], normal_on_edges[i]));
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);
    }

    void update_bishop(){

        for (size_t iE = 0; iE < curveEdgesId.size(); iE++) {
            size_t i0 = curveEdgesId[iE][0];
            size_t i1 = curveEdgesId[iE][1];
            tangent_on_edges[iE] = glm::normalize(controlpoints[i1] - controlpoints[i0]);
        }
        curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
        for (size_t i = 1; i < curveEdgesId.size(); i++) {
            point d_tangent = tangent_on_edges[i] - tangent_on_edges[i - 1];
            normal_on_edges[i] = glm::normalize(d_tangent);
        }
        normal_on_edges[0] = glm::normalize(glm::cross(reference, tangent_on_edges[0]));
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 0; i < curveEdgesId.size(); i++) {
            binormal_on_edges[i] = glm::normalize(glm::cross(tangent_on_edges[i], normal_on_edges[i]));
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);
    
    }

    void gen_parallel_transport() {
       

        for (size_t i = 1; i < controlpoints.size(); i++) {
            auto prev_tan = tangent_on_edges[i - 1];
            point next_tan = tangent_on_edges[(i) % tangent_on_edges.size()];

            glm::vec3 axis = glm::cross(prev_tan, next_tan);
            if (glm::length(axis) > 0.0001) {
                axis = glm::normalize(axis);
                float angle = acos(glm::dot(glm::normalize(prev_tan), glm::normalize(next_tan)));
                turn_angles[i] = angle;
                glm::mat3 rotationMatrix = glm::rotate(glm::mat4(1.0f), angle, axis);
                parallel_transport[i] = rotationMatrix;
            } else {
                parallel_transport[i] = glm::mat3(1.0f);
            }
        }
        curve->addNodeScalarQuantity("turning angles", turn_angles);
    }

    void gen_twist_angle() {
        gen_parallel_transport();

        twist_angles.resize(curveEdgesId.size());
        twist_angles[0] = PI;
        glm::mat3 matrix = glm::rotate(glm::mat4(1.0f), twist_angles[0], tangent_on_edges[0]);
        material_u[0] = glm::normalize(matrix * binormal_on_edges[0]);
        material_v[0] = glm::normalize(glm::cross(tangent_on_edges[0], material_u[0]));

        for (size_t i = 1; i < curveEdgesId.size(); i++) {
            material_u[i] = glm::normalize(parallel_transport[i] * material_u[i - 1]);
            twist_angles[i] = glm::acos(glm::dot(material_u[i], binormal_on_edges[i]));
            material_v[i] = glm::normalize(glm::cross(tangent_on_edges[i], material_u[i]));
        }
        curve->addEdgeScalarQuantity("twisting angle", twist_angles);
        curve->addEdgeVectorQuantity("material_v", material_v);
        curve->addEdgeVectorQuantity("material_u", material_u);
    }


    void cal_curvature() {
        curvature.resize(curveEdgesId.size());
        for (size_t i = 0; i < curveEdgesId.size(); i++) {
            curvature[i] = 2 * glm::tan(turn_angles[i] / 2);
        }
        curve->addEdgeScalarQuantity("curvature", curvature);
    }

    void gen_vertex_weight() {
        vertex_weight.resize(controlpoints.size());
        vertex_weight[0] = 0.5 * (glm::length(curveEdges[0]));
        vertex_weight[controlpoints.size() - 1] = 0.5 * (glm::length(curveEdges[curveEdges.size() - 1]));
        for (size_t i = 1; i < controlpoints.size() - 1; i++) {
            vertex_weight[i] = 0.5 * (glm::length(curveEdges[i - 1]) + glm::length(curveEdges[i]));
        }

        curve->addNodeScalarQuantity("vertex_weight", vertex_weight);
    }

    void cal_darboux() {
        darboux.resize(curveEdgesId.size());
        for (size_t i = 1; i < curveEdgesId.size(); i++) {
            darboux[i] = curvature[i] * binormal_on_edges[i];
        }
        curve->addEdgeVectorQuantity("darboux", darboux);
    }

    void cal_bending_force() {
        bending_force.resize(curveEdgesId.size());
        for (size_t i = 1; i < curveEdgesId.size() - 1; i++) {
            auto cur_edge = curveEdges[i];
            auto next_edge = curveEdges[i + 1];
            float co =
                2 / (vertex_weight[i] * glm::length(cur_edge) * glm::length(next_edge) + glm::dot(cur_edge, next_edge));
            for (size_t j = 0; j < curveEdgesId.size() - 1; j++) {
                if (j == i + 1) {
                    bending_force[i] +=
                        glm::cross(-1.f * next_edge, darboux[j]) + glm::dot(darboux[j], next_edge) * darboux[j];
                } else if (j == i) {
                    bending_force[i] +=
                        glm::cross(-1.f * cur_edge, darboux[j]) + glm::dot(darboux[j], cur_edge) * darboux[j] +
                        glm::cross(-1.f * next_edge, darboux[j]) + glm::dot(darboux[j], next_edge) * darboux[j];
                } else if (j == i - 1) {
                    bending_force[i] +=
                        glm::cross(-1.f * cur_edge, darboux[j]) + glm::dot(darboux[j], cur_edge) * darboux[j];
                } else {
                    continue;
                }
            }
            bending_force[i] *= co;
            bending_force[i] *= -2;
        }
        curve->addEdgeVectorQuantity("bending force", bending_force);
    }

    // void cal_twisting_force() {
    //     twisting_force.resize(curveEdgesId.size());
    //     for (size_t i = 1; i < curveEdgesId.size(); i++) {
    //         float L = 0.5 * std::accumulate(vertex_weight.begin(), vertex_weight.begin() + i, 0.0f);
    //         auto d_twist = twist_angles[i] - twist_angles[0];
    //         auto delta_prev_Psi = 0.5f * darboux[i] / glm::length(curveEdges[i - 1]);
    //         auto delta_next_Psi = -0.5f * darboux[i] / glm::length(curveEdges[i]);
    //         auto delta_Psi = -(delta_prev_Psi + delta_next_Psi);
    //         twisting_force[i] = d_twist / L * delta_Psi;
    //     }
    //     curve->addEdgeVectorQuantity("twisting force", twisting_force);
    // }

    void cal_twisting_force() {
        
        for (size_t i = 1; i < curveEdgesId.size(); i++) {
            float L = 0.5 * std::accumulate(vertex_weight.begin(), vertex_weight.begin() + i, 0.0f);
            auto d_twist = twist_angles[i] - twist_angles[0];
            twisting_force[i] = 2 * d_twist / L;
        }
        curve->addEdgeScalarQuantity("twisting force", twisting_force);
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

        for(size_t iter = 0; iter < 1; iter++){
            for (size_t i = 1; i < controlpoints.size()-1; i++) {
                acceleration[i] = lr * bending_force[i]/vertex_weight[i]; 
                point new_x =  controlpoints[i] + v0 * dt + 0.5f * acceleration[i] * dt * dt;
                controlpoints[i] = new_x;
            }
            curve->updateNodePositions(controlpoints);
            update_bishop();
            gen_vertex_weight();
            gen_twist_angle();
            cal_curvature();
            cal_darboux();
            cal_bending_force();

        } 
    }

    void updatupdate_all() {
        for (size_t i = 1; i < curveEdgesId.size(); i++) {
            glm::mat3 rot = glm::rotate(glm::mat4(1.0f), twist_angles[i], tangent_on_edges[i]);
            material_u[i] = rot * binormal_on_edges[i];
            material_v[i] = rot * normal_on_edges[i];
        }
        curve->addEdgeVectorQuantity("material_v", material_v);
        curve->addEdgeVectorQuantity("material_u", material_u);
    }

    void updatematerial_frame() {
        for (size_t i = 1; i < curveEdgesId.size(); i++) {
            glm::mat3 rot = glm::rotate(glm::mat4(1.0f), twist_angles[i], tangent_on_edges[i]);
            material_u[i] = rot * binormal_on_edges[i];
            material_v[i] = rot * normal_on_edges[i];
        }
        curve->addEdgeVectorQuantity("material_v", material_v);
        curve->addEdgeVectorQuantity("material_u", material_u);
    }


  private:

    const point reference = point(0, 0, 1);


    std::vector<point> controlpoints;
    std::vector<std::vector<size_t>> curveEdgesId;
    std::vector<point> curveEdges;
    CurveNetwork* curve;

    std::vector<double> arc_length;

    std::vector<point> tangent_on_edges;
    std::vector<point> normal_on_edges;
    std::vector<point> binormal_on_edges;


    std::vector<glm::mat3> parallel_transport;
    std::vector<float> twist_angles;
    std::vector<point> material_v;
    std::vector<point> material_u;


    std::vector<float> turn_angles;
    std::vector<float> curvature;

    std::vector<point> darboux;
    std::vector<float> vertex_weight;
    std::vector<point> bending_force;
    std::vector<float> twisting_force;

    std::vector<point> acceleration;
};