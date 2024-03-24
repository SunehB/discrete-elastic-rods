#include "polyscope/curve_network.h"
#include <vector>

// Declaration
class BezierPoint;
class Beziermanager;


    using point = glm::vec3;

class BezierPoint{
public:

    point default_shift = 2.f * point(1,1,1);
    BezierPoint(point pos):position(pos){
        L_handle = pos - default_shift;
        R_handle = pos + default_shift;
        L_enabled = false;
        R_enabled = false;
    }
    ~BezierPoint(){};

    point position;
    point L_handle;
    point R_handle;
    bool L_enabled;
    bool R_enabled;
private:

};



class Beziermanager{


public:

    Beziermanager(){
    };

    void init(){
        BezierPoint p1(point(1,1,1));
        BezierPoint p2(point(1,0,1));
        BezierPoint p3(point(4,2,4));
        BezierPoint p4(point(3,5,3));
        controlPoints.push_back(p1);
        controlPoints.push_back(p2);
        // controlPoints.push_back(p3);
        // controlPoints.push_back(p4);
    }
    ~Beziermanager(){};

    void generateCurve() {
        if (controlPoints.size() < 2) return; 

        for (size_t i = 0; i < controlPoints.size() - 1; ++i) {
            auto& p0 = controlPoints[i].position;
            auto& p1 = controlPoints[i].R_handle;
            auto& p2 = controlPoints[i + 1].L_handle;
            auto& p3 = controlPoints[i + 1].position;

            for (float t = 0; t <= 1; t += dt) {
                point p = calculateBezierPoint(t, p0, p1, p2, p3);
                curvePoints.push_back(p);
            }
        }   

        accmu_arcLengths.push_back(0);
        for (size_t i = 1; i < curvePoints.size(); ++i) {
            double distance = glm::distance(curvePoints[i], curvePoints[i-1]);
            accmu_arcLengths.push_back(accmu_arcLengths.back() + distance); // Add distance to the last value
        }

        double totalLength = accmu_arcLengths.back();
        
        int numNewPoints = curvePoints.size();
        std::vector<point> newPoints;
        
        for (size_t i = 0; i < curvePoints.size() - 1; ++i) {
            curveEdges.push_back({i, i + 1});
        }

         random_curve = polyscope::registerCurveNetwork("random_curve", curvePoints,curveEdges);
        return;
    }

    void gen_randomvec(){
        std::vector<point> randomvector_nodes(curvePoints.size());
        for (size_t iN = 0; iN < curvePoints.size(); iN++) {
            randomvector_nodes[iN] = glm::vec3{polyscope::randomUnit() - .5, polyscope::randomUnit() - .5, polyscope::randomUnit() - .5};
        }
        random_curve->addNodeVectorQuantity("random vector on nodes", randomvector_nodes);
    }

    void gen_binormal(){
        tangent_on_edges.resize(curveEdges.size());
        normal_on_edges.resize(curveEdges.size());
        binormal_on_edges.resize(curveEdges.size());

        for (size_t iE = 0; iE < curveEdges.size(); iE++) {
            size_t i0 = curveEdges[iE][0];
            size_t i1 = curveEdges[iE][1];
            tangent_on_edges[iE] = glm::normalize(curvePoints[i1] - curvePoints[i0]);
        }
        random_curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
       
        for (size_t i = 0; i < curveEdges.size(); i++) {
            std::cout << tangent_on_edges[i].x << " " << tangent_on_edges[i].y << " " << tangent_on_edges[i].z << std::endl;
        }
        for (size_t i = 1; i < curveEdges.size(); i++) {
            point d_tangent = tangent_on_edges[i] - tangent_on_edges[i - 1];
            point edgeNormal = glm::normalize(d_tangent);
            normal_on_edges[i] = edgeNormal;
            if (glm::dot(normal_on_edges[i], normal_on_edges[i-1]) < 0) {
                normal_on_edges[i] = -normal_on_edges[i]; // 反转方向
            }
        }
        normal_on_edges[0] = glm::normalize(tangent_on_edges[0]);
        random_curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

        for (size_t i = 0; i < curveEdges.size(); i++) {
            binormal_on_edges[i] = glm::cross(tangent_on_edges[i], normal_on_edges[i]);
        }
        random_curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);
    }
    
    // void gen_truningAngles(){
    //     truningAngles_on_nodes.resize(curvePoints.size());
    //     for (size_t iN = 1; iN < curvePoints.size()-1; iN++) {
    //         point edg0 = curvePoints[iN] - curvePoints[iN - 1];
    //         point edg1 = curvePoints[iN + 1] - curvePoints[iN];
    //         truningAngles_on_nodes[iN] = glm::acos(glm::dot(edg0, edg1) / (glm::length(edg0) * glm::length(edg1)));
    //     }
    //     random_curve->addNodeScalarQuantity("truningAngles on nodes", truningAngles_on_nodes);
    // }

    void get_bishopframe(){ 
        rotated_binormal_on_egdes.resize(curveEdges.size());
        rotated_normal_on_egdes.resize(curveEdges.size());

        int n = curvePoints.size();
        twistAngles_on_nodes = std::vector<double>(n, 0);
        double angleIncrement = 2 * M_PI / (n - 1);
        for (int i = 0; i < n; ++i) {
            twistAngles_on_nodes[i] = i * angleIncrement;
        }

        for (size_t iN = 1; iN < curveEdges.size()-1; iN++) {
            float cur_theta = twistAngles_on_nodes[iN];
            float cos_theta = glm::cos(cur_theta);
            float sin_theta = glm::sin(cur_theta);
            point m_1 = normal_on_edges[iN] * cos_theta + binormal_on_edges[iN] * sin_theta;
            point m_2 = -normal_on_edges[iN] * sin_theta + binormal_on_edges[iN] * cos_theta;
            rotated_normal_on_egdes[iN] = m_1;
            rotated_binormal_on_egdes[iN] = m_2;
        }

        random_curve->addEdgeVectorQuantity("rotated_binormal", rotated_binormal_on_egdes);
        random_curve->addEdgeVectorQuantity("rotated_normal", rotated_normal_on_egdes);
    }

private:

    double dt = 0.02 ;
    std::vector<double>accmu_arcLengths;

    polyscope::CurveNetwork* random_curve;
    std::vector<BezierPoint> controlPoints;
    std::vector<point> curvePoints;
    std::vector<std::array<size_t, 2>> curveEdges;
    std::vector<point> tangent_on_edges;
    std::vector<point> normal_on_edges;
    std::vector<point> binormal_on_edges;

    // std::vector<double> truningAngles_on_nodes;
    std::vector<double> twistAngles_on_nodes;

    std::vector<point> rotated_normal_on_egdes;
    std::vector<point> rotated_binormal_on_egdes;;


    


    point calculateBezierPoint(float t, const point& p0, const point& p1, const point& p2, const point& p3) {
         float u = 1 - t;
        float tt = t * t;
        float uu = u * u;
        float uuu = uu * u;
        float ttt = tt * t;

        point p = uuu * p0; // 第一个项
        p += 3 * uu * t * p1; // 第二个项
        p += 3 * u * tt * p2; // 第三个项
        p += static_cast<float>(glm::pow(u, 3)) * p3; // 第四个项

        return p;
    }
};