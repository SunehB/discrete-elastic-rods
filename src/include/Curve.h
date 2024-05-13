#include "polyscope/curve_network.h"
#include "polyscope/utilities.h"
#include <Eigen/Sparse>

#include <vector>

using namespace polyscope;
using namespace glm;

class Curve {
    using point = dvec3;

  public:
    Curve(){};
    ~Curve(){};


    size_t nSamples = 100;
    bool closed = true;


    void initcurve() {
        double dx = 1.0 / nSamples;

        for (size_t i = 0; i < nSamples; ++i) {
            double t = i * dx;
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
            edge_length[i] = length(curveEdges[i]);
        }
        totallength = std::accumulate(edge_length.begin(), edge_length.end(), 0.0f);
        for (size_t i = 0; i < nSamples; i++) {
            arc_length[i] = std::accumulate(edge_length.begin(), edge_length.begin() + i, 0.0f);
        }
        curve = registerCurveNetwork("Spiral", controlpoints, curveEdgesId);

        cal_tangent();
        normal_on_edges[0] = normalize(cross(tangent_on_edges[0], reference));
        binormal_on_edges[0] = normalize(cross(tangent_on_edges[0], normal_on_edges[0]));

        for (size_t i = 1; i < num_edges; i++) {
            point t_i = tangent_on_edges[i];
            point t_i_1 = tangent_on_edges[(i - 1) % num_edges];
            point axis = cross(t_i_1, t_i);
            double theta = std::atan2(length(axis), dot(t_i_1, t_i));
            dmat3 rotationMatrix = rotate(dmat4(1), theta, normalize(axis));
            if (theta > 1e-10)
                normal_on_edges[i] = normalize(rotationMatrix * normal_on_edges[i - 1]);
            else
                normal_on_edges[i] = normal_on_edges[i - 1];
            rdy = true;
        }
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 1; i < num_edges; i++)
            binormal_on_edges[i] = normalize(cross(tangent_on_edges[i], normal_on_edges[i]));

        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);

        update_material_frame();
        cal_attrs();
        cal_bending_force();
        cal_twisting_force();
    }

    void init_catenary() {
        double a = 10.0;                                // 控制悬链线弯曲程度的参数
        double x_start = -15.0;                         // x的起始值
        double x_end = 15.0;                            // x的结束值
        double dx = (x_end - x_start) / (nSamples - 1); // x增量


        for (size_t i = 0; i < nSamples; ++i) {
            double x = x_start + i * dx;
            point p;
            p.x = x;
            p.y = 0;
            p.z = a * std::cosh(x / a);
            controlpoints.push_back(p);
        }

        num_controlpoints = controlpoints.size();
        num_edges = nSamples - 1;
        resize_vectors();

        for (size_t i = 0; i < num_edges; i++) {
            curveEdgesId.push_back({i, i + 1});
            curveEdges.push_back(controlpoints[i + 1] - controlpoints[i]);
            edge_length[i] = length(curveEdges[i]);
        }
        totallength = std::accumulate(edge_length.begin(), edge_length.end(), 0.0f);
        for (size_t i = 0; i < nSamples - 1; i++) {
            arc_length[i] = std::accumulate(edge_length.begin(), edge_length.begin() + i, 0.0f);
        }
        curve = registerCurveNetwork("Catenary", controlpoints, curveEdgesId);

        cal_tangent();
        normal_on_edges[0] = normalize(cross(tangent_on_edges[0], reference));
        binormal_on_edges[0] = normalize(cross(tangent_on_edges[0], normal_on_edges[0]));

        for (size_t i = 1; i < num_edges; i++) {
            point t_i = tangent_on_edges[i];
            point t_i_1 = tangent_on_edges[(i - 1) % num_edges];
            point axis = cross(t_i_1, t_i);
            double theta = std::atan2(length(axis), dot(t_i_1, t_i));
            dmat3 rotationMatrix = rotate(dmat4(1), theta, normalize(axis));
            if (theta > 1e-10)
                normal_on_edges[i] = normalize(rotationMatrix * normal_on_edges[i - 1]);
            else
                normal_on_edges[i] = normal_on_edges[i - 1];
            rdy = true;
        }
        curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 1; i < num_edges; i++)
            binormal_on_edges[i] = normalize(cross(tangent_on_edges[i], normal_on_edges[i]));

        curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);

        update_material_frame();
        cal_attrs();
        cal_bending_force();
        cal_twisting_force();
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
        darbouxdenom.resize(num_controlpoints);
        vertex_weight.resize(num_controlpoints);

        twist_thetas.resize(num_edges);
        material_u.resize(num_edges);
        material_v.resize(num_edges);

        velocity.resize(num_controlpoints);
        acceleration.resize(num_controlpoints);
        bending_force.resize(num_controlpoints);
        twisting_force.resize(num_controlpoints);

        newpoints.resize(num_controlpoints);
        newedges.resize(num_edges);
        newedge_length.resize(num_edges);
    }
    void cal_tangent();
    void update_bishop();
    void update_material_frame();
    void cal_attrs();
    void cal_bending_force();
    void cal_twisting_force();
    void cal_velocity();
    void fastprojection();
    void fastsuiteprojection();

    dmat3 skew(point v) {
        return dmat3(0, -v.z, v.y, v.z, 0, -v.x, -v.y, v.x, 0);
    }


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

        initcurve();
    }

    bool status() {
        return rdy;
    };

  private:
    bool rdy = false;
    size_t num_controlpoints;
    size_t num_edges;

    double eps = 1e-6;

    double alpha = 1;
    double beta = 1;
    double dt = 0.1;

    point last_tangent;

    const point reference = point(0, 0, 1);
    double totaltwist = 2 * PI;
    double totallength;

    std::vector<point> controlpoints;
    std::vector<std::vector<size_t>> curveEdgesId;
    std::vector<point> curveEdges;
    std::vector<double> edge_length;
    std::vector<double> dual_length;
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
    std::vector<double> darbouxdenom;
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