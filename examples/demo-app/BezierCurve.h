#include "polyscope/curve_network.h"
#include <vector>

// Declaration
class BezierPoint;
class Beziermanager;


    using point = glm::vec3;

class BezierPoint{
public:
    BezierPoint(point pos):position(pos){
        L_handle = pos;
        R_handle = pos;
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
        BezierPoint p1(point(0,0,0));
        BezierPoint p2(point(1,1,1));
        controlPoints.push_back(p1);
        controlPoints.push_back(p2);
    }
    ~Beziermanager(){};

    void generateCurve() {
        if (controlPoints.size() < 2) return; // 如果控制点不足2个，直接返回空集合

        // 遍历每对相邻的控制点
        for (size_t i = 0; i < controlPoints.size() - 1; ++i) {
            auto& p0 = controlPoints[i].position;
            auto& p1 = controlPoints[i].R_handle;
            auto& p2 = controlPoints[i + 1].L_handle;
            auto& p3 = controlPoints[i + 1].position;

            for (float t = 0; t <= 1; t += dt) {
                curvePoints.push_back(calculateBezierPoint(t, p0, p1, p2, p3));
            }
        }   
        for(size_t i = 0; i < curvePoints.size() - 1; ++i){
            curveEdges.push_back(std::make_pair(curvePoints[i], curvePoints[i+1]));
        }

        return;
    }

    void registercurve(){
        random_curve = polyscope::registerCurveNetwork("random_curve", curvePoints,curveEdges);

    }   
    double dt = 0.1 ;

private:
    polyscope::CurveNetwork* random_curve;
    std::vector<BezierPoint> controlPoints;
    std::vector<point> curvePoints;
    std::vector<std::pair<point, point>> curveEdges;

    point calculateBezierPoint(float t, const point& p0, const point& p1, const point& p2, const point& p3) {
         float u = 1 - t;
        float tt = t * t;
        float uu = u * u;
        float uuu = uu * u;
        float ttt = tt * t;

        point p = uuu * p0; // 第一个项
        p += 3 * uu * t * p1; // 第二个项
        p += 3 * u * tt * p2; // 第三个项
        p += ttt * p3; // 第四个项

        return p;
    }
};


// Point2D PointOnCubicBezier( Point2D* cp, float t )
// {
//     float   ax, bx, cx;
//     float   ay, by, cy;
//     float   tSquared, tCubed;
//     Point2D result;

//     /*計算多項式係數*/

//     cx = 3.0 * (cp[1].x - cp[0].x);
//     bx = 3.0 * (cp[2].x - cp[1].x) - cx;
//     ax = cp[3].x - cp[0].x - cx - bx;

//     cy = 3.0 * (cp[1].y - cp[0].y);
//     by = 3.0 * (cp[2].y - cp[1].y) - cy;
//     ay = cp[3].y - cp[0].y - cy - by;

//     /*計算位於參數值t的曲線點*/

//     tSquared = t * t;
//     tCubed = tSquared * t;

//     result.x = (ax * tCubed) + (bx * tSquared) + (cx * t) + cp[0].x;
//     result.y = (ay * tCubed) + (by * tSquared) + (cy * t) + cp[0].y;

//     return result;
// }


