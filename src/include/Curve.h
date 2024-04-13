#include "polyscope/curve_network.h"
#include <Eigen/Sparse>

#include <vector>

using namespace polyscope;
using namespace glm;

// std::ostream& operator<<(std::ostream& os, const glm::dvec3& v) {
//     os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
//     return os;
// }

class Curve {
    using point = dvec3;

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


    void initcurve(int inputn) {
        nSamples = inputn;
        double dx = 1.0 / nSamples;

        for (size_t i = 0; i < nSamples; ++i) {
            double t = i * dx;
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
        for (size_t i = 0; i < nSamples; i++) {
            arc_length[i] = std::accumulate(edge_length.begin(), edge_length.begin() + i, 0.0f);
            // std::cout << "arc_length " << i << ": " << arc_length[i] << std::endl;
        }
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
            rdy = true;
        }
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 1; i < num_edges; i++) {
            binormal_on_edges[i] = normalize(cross(tangent_on_edges[i], normal_on_edges[i]));
            // std::cout << "binormal on edges " << i << ": " << binormal_on_edges[i] << std::endl;
        }
        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);

        update_material_frame();
        cal_attrs();
        cal_bending_force();
        cal_twisting_force();
    }

    void cal_tangent();

    void update_bishop();

    void update_material_frame();

    void cal_attrs();

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

    void cal_bending_force();

    dmat3 transpoe(point v) {
        return dmat3(v.x, 0, 0, v.y, 0, 0, v.z, 0, 0);
    }

    void cal_twisting_force();


    void cal_velocity();

    void fastprojection();
    void fastsuiteprojection();

    void loop() {
        cal_velocity();
        fastprojection();
        curve->updateNodePositions(controlpoints);
        cal_tangent();
        update_bishop();
        update_material_frame();
        cal_attrs();
        cal_bending_force();
        cal_twisting_force();
    }

    void set_alpha(double a) {
        alpha = a;
    }
    void set_beta(double b) {
        beta = b;
    }
    void set_dt(double d) {
        dt = d;
    }
    void set_nSamples(size_t n) {
        nSamples = n;
    }
    void set_totaltwist(double t) {
        totaltwist = t;
    }

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

        initcurve(nSamples);
    }

    bool status() {
        return rdy;
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