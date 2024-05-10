#include "glm/fwd.hpp"
#include "glm/geometric.hpp"
#include "polyscope/curve_network.h"
#include <Eigen/Sparse>

#include <algorithm>
#include <vector>

using namespace polyscope;
using namespace glm;


using point = dvec3;

struct FrameEdge;

struct Vertex {
    point pos;
    double vertex_weight;
    point darboux;
    double curvature;
    double turning_angle;

    FrameEdge* Ledge;
    FrameEdge* Redge;

    point acceleration;
    point velocity;

    point bending_force;
    point twisting_force;

    Vertex* Vnext;
    Vertex* Vprev; 
};

std::ostream& operator<<(std::ostream& os, const Vertex& v) {
    os << "Vertex Position: (" << v.pos.x << ", " << v.pos.y << ", " << v.pos.z << ")\n"
       << "Weight: " << v.vertex_weight << "\n"
       << "Curvature: " << v.curvature << "\n"
       << "Turning Angle: " << v.turning_angle << "\n"
       << "Acceleration: (" << v.acceleration.x << ", " << v.acceleration.y << ", " << v.acceleration.z << ")\n"
       << "Velocity: (" << v.velocity.x << ", " << v.velocity.y << ", " << v.velocity.z << ")\n"
       << "Bending Force: (" << v.bending_force.x << ", " << v.bending_force.y << ", " << v.bending_force.z << ")\n"
       << "Twisting Force: (" << v.twisting_force.x << ", " << v.twisting_force.y << ", " << v.twisting_force.z << ")";
    return os;
}

struct FrameEdge {
    point edge;
    double length;
    double arc_length;

    double twist_theta;

    Vertex* Lvertex;
    Vertex* Rvertex;

    point tangent;
    point bishop_normal;
    point bishop_binormal;

    point material_u;
    point material_v;

    point constraint;
};

std::ostream& operator<<(std::ostream& os, const FrameEdge& e) {
    os << "Edge: (" << e.edge.x << ", " << e.edge.y << ", " << e.edge.z << ")\n"
       << "Length: " << e.length << "\n"
       << "Arc Length: " << e.arc_length << "\n"
       << "Twist Theta: " << e.twist_theta << "\n"
       << "Tangent: (" << e.tangent.x << ", " << e.tangent.y << ", " << e.tangent.z << ")\n"
       << "Bishop Normal: (" << e.bishop_normal.x << ", " << e.bishop_normal.y << ", " << e.bishop_normal.z << ")\n"
       << "Bishop Binormal: (" << e.bishop_binormal.x << ", " << e.bishop_binormal.y << ", " << e.bishop_binormal.z << ")\n"
       << "Material U: (" << e.material_u.x << ", " << e.material_u.y << ", " << e.material_u.z << ")\n"
       << "Material V: (" << e.material_v.x << ", " << e.material_v.y << ", " << e.material_v.z << ")\n"
       << "Constraint: (" << e.constraint.x << ", " << e.constraint.y << ", " << e.constraint.z << ")";
    return os;
}


class Curveoop {

  public:
    Curveoop(){};
    ~Curveoop(){};


    void init(int inputn) {
        nSamples = inputn;
        double dx = 1.0 / nSamples;

        for (size_t i = 0; i < nSamples; ++i) {
            double t = i * dx;
            point p;
            p.x = cos(2 * M_PI * t);
            p.y = sin(2 * M_PI * t);
            p.z = 0.3 * sin(4 * M_PI * t);
            Vertex v;
            v.pos = p;
            Vertices.push_back(v);
        }

        for (size_t i = 0; i < nSamples; ++i) {

            FrameEdge e;
            e.Lvertex = &Vertices[i];
            e.Rvertex = &Vertices[(i + 1) % nSamples];
            e.edge = e.Rvertex->pos - e.Lvertex->pos;
            e.length = length(e.edge);
            e.tangent = normalize(e.edge);
            Edges.push_back(e);
            totallength += e.length;

            Vertices[i].Ledge = &Edges[(i - 1 + nSamples) % nSamples];
            Vertices[i].Vprev = &Vertices[(i - 1 + nSamples) % nSamples];
            Vertices[i].Redge = &Edges[i];
            Vertices[i].Vnext = &Vertices[(i + 1) % nSamples];
        }

        cal_attrs();

        std::transform(Vertices.begin(), Vertices.end(), controlpoints.begin(), [](const Vertex& v) { return v.pos; });
        std::transform(Edges.begin(), Edges.end(), edges.begin(), [](const FrameEdge& e) { return e.edge; });

        curve = registerCurveNetwork("Michelle", controlpoints, edges);

        FrameEdge initframe;
        initframe.bishop_normal = cross(Edges[0].tangent, reference);
        initframe.bishop_binormal = cross(Edges[0].tangent, initframe.bishop_normal);

        update_bishop(initframe);
        update_material_frame();


        cal_bending_force();
        cal_twisting_force();
    }
    void twist_holonomy();


    void update_bishop(FrameEdge initframe) {

        Edges[0].bishop_normal = initframe.bishop_normal;
        Edges[0].bishop_binormal = initframe.bishop_binormal;

        for (size_t i = 1; i < nSamples; ++i) {
            point axis = cross(Edges[i].Lvertex->Ledge->tangent, Edges[i].tangent);
            double angle = acos(dot(Edges[i].Lvertex->Ledge->tangent, Edges[i].tangent));

            dmat3 R = rotate(dmat4(1), angle, normalize(axis));
            Edges[i].bishop_normal = R * Edges[i].Lvertex->Ledge->bishop_normal;
            Edges[i].bishop_binormal = R * Edges[i].Lvertex->Ledge->bishop_binormal;
        }
    }

    void update_material_frame() {
        for (size_t i = 0; i < nSamples; i++) {
            double curtwist = totaltwist * Edges[i].arc_length / totallength;
            dmat3 rot = rotate(dmat4(1.0f), curtwist, Edges[i].tangent);
            Edges[i].material_u = rot * Edges[i].bishop_normal;
            Edges[i].material_v = rot * Edges[i].bishop_binormal;
        }

        std::transform(Edges.begin(), Edges.end(), material_u.begin(), [](const FrameEdge& e) { return e.material_u; });
        std::transform(Edges.begin(), Edges.end(), material_v.begin(), [](const FrameEdge& e) { return e.material_v; });
        curve->addEdgeVectorQuantity("material_u", material_u);
        curve->addEdgeVectorQuantity("material_v", material_v);
    }

    void cal_attrs() {
        // totallength = 0;
        // for (size_t i = 0; i < nSamples; ++i) {
        //     // Edges[i].length = length(Edges[i].edge);
        //     totallength += Edges[i].length;
        // }

        for (size_t i = 0; i < nSamples; ++i) {
            if (i == 0) {
                Edges[i].arc_length = 0;
            } else {
                Edges[i].arc_length = Edges[i - 1].arc_length + Edges[i].length / totallength;
            }
            FrameEdge* edgeL = Vertices[i].Ledge;
            FrameEdge* edgeR = Vertices[i].Redge;
            std::cout << "edgeL " << &edgeL->length << std::endl;
            std::cout << "edgeR " << &edgeR->length << std::endl;
            Vertices[i].vertex_weight = 0.5 * (edgeL->length + edgeR->length);

            Vertices[i].turning_angle = acos(dot(edgeL->tangent, edgeR->tangent));
            Vertices[i].darboux = cross(edgeL->tangent, edgeR->tangent) /
                                  (edgeL->length * edgeR->length + dot(edgeL->edge, edgeR->edge));
        }
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


    dmat3 transpoe(point v) {
        return dmat3(v.x, 0, 0, v.y, 0, 0, v.z, 0, 0);
    }


    void cal_bending_force(){
        for (size_t i = 0; i < nSamples; i++) {
        for (size_t j = 0; j < nSamples; j++) {
            FrameEdge* edgeL = Vertices[i].Ledge;
            FrameEdge* edgeR = Vertices[i].Redge;

            dmat3 Rmatrix = transpoe(edgeR->edge);
            dmat3 Lmatrix = transpoe(edgeL->edge);

            dmat3 RT = transpose(Rmatrix);
            dmat3 LT = transpose(Lmatrix);

            double w = 0.5 * Vertices[i].vertex_weight;
            auto coeff = -2 / (w * (edgeR->length* edgeL->length + dot(edgeR->edge, edgeL->edge)));
            // std::cout << "coeff " << coeff << std::endl;
            auto kbj =Vertices[i].darboux;
            // dmat3 kbjT = transpose(transpoe(kbj));
            dmat3 kbjTR = outerProduct(kbj, edgeR->edge);
            dmat3 kbjTRT = transpose(kbjTR);

            dmat3 kbjTL = outerProduct(kbj, edgeL->edge);
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
                p = 2.0 * cross(edgeR->edge, kbj) + kbjTRTkbj;
            } else if (i == j) {
                p = -1.0 * (2.0 * cross(edgeR->edge, kbj) + kbjTRTkbj + 2.0 * cross(edgeL->edge, kbj) - kbjTLTkbj);
            } else if (i == j + 1) {
                p = 2.0 * cross(edgeL->edge, kbj) - kbjTLTkbj;
            }
            Vertices[i].bending_force += coeff * p;
            // std::cout << "bending force " << p << std::endl;
        }
        Vertices[i].bending_force *= -1;
        // std::cout << "bending force " << bending_force[i] << std::endl;
    }
    // std::cout << "bending force " << bending_force[0] << std::endl;
    
    std::transform(Vertices.begin(), Vertices.end(), bending_force.begin(), [](const Vertex& v) { return v.bending_force; });
    curve->addNodeVectorQuantity("bending force", bending_force);


    }

    void cal_twisting_force() {
        // twisting_force.clear();
        // twisting_force.resize(nSamples);

        for (size_t i = 0; i < nSamples; i++) {
            auto Vprevd = 0.5 * Vertices[i].Vprev->darboux / Vertices[i].Ledge->length;
            auto Vnextd = -0.5 * Vertices[i].Vprev->darboux / Vertices[i].Redge->length;
            Vertices[i].twisting_force = -1 * totaltwist * (Vnextd + Vprevd + eps) / (0.5f * totallength);
            // std::cout << "twisting force " << twisting_force[i] << std::endl;
        }
        
        std::transform(Vertices.begin(), Vertices.end(), twisting_force.begin(), [](const Vertex& v) { return v.twisting_force; });
        curve->addNodeVectorQuantity("twisting force", twisting_force);
    }


    void cal_velocity(){
        for (size_t i = 0; i < nSamples; i++) {
            Vertices[i].acceleration = (alpha * Vertices[i].bending_force + beta* Vertices[i].twisting_force)/Vertices[i].vertex_weight;
            Vertices[i].velocity = Vertices[i].velocity + dt * Vertices[i].acceleration;
        }

    }

    // void fastprojection(){


    // };
    void fastsuiteprojection();

    void run_loop() {
        cal_velocity();
        // fastprojection();
        curve->updateNodePositions(controlpoints);
        FrameEdge initframe;
        initframe.bishop_normal = Edges[0].bishop_normal;
        initframe.bishop_binormal = Edges[0].bishop_binormal;

        update_bishop(initframe);
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


    bool status() {
        return rdy;
    }

  private:
    bool rdy = false;
    bool loop = false;

    double eps = 1e-8;
    double alpha = 1;
    double beta = 1;
    double dt = 1e-3;
    size_t nSamples = 100;

    const point reference = point(0, 0, 1);
    double totaltwist = PI;
    double totallength;

    CurveNetwork* curve;

    std::vector<Vertex> Vertices;
    std::vector<FrameEdge> Edges;

    std::vector<point> controlpoints;
    std::vector<point> edges;
    std::vector<point> material_u;
    std::vector<point> material_v;

    std::vector<point> bending_force;
    std::vector<point> twisting_force;

    std::vector<double> constraint;
};