#include "glm/ext/vector_float3.hpp"
#include "imgui.h"
#include "include/Bez.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/view.h"
#include <vector>

static bool isRunning = false;
static int currentIteration = 0;
static float alpha = 1.f;
static float beta = 1.f;
static float dt = 0.001f;
static float totaltwist = 0.0f;
static int nSample = 100;
static int nIteration = 100;
static bool edit = false;

using namespace polyscope;
using point = glm::vec3;

PointCloud* pc = nullptr;
CurveNetwork* curve = nullptr;
std::vector<point> points;
std::vector<std::vector<size_t>> curveEdgesId;

void mouse_unfc(ImGuiIO& io) {
    std::cout << "Mouse unfc" << std::endl;
    if (io.MouseDoubleClicked[0]) { // if the left mouse button was clicked
        std::cout << "Mouse unfc0" << std::endl;
        // gather values
        glm::vec2 screenCoords{io.MousePos.x, io.MousePos.y};
        glm::vec3 worldRay = polyscope::view::screenCoordsToWorldRay(screenCoords);
        glm::vec3 worldPos = polyscope::view::screenCoordsToWorldPosition(screenCoords);
        // std::pair<polyscope::Structure*, size_t> pickPair =
        //     polyscope::pick::evaluatePickQuery(screenCoords.x, screenCoords.y);

        // print some values
        std::cout << "    io.MousePos.x: " << io.MousePos.x << " io.MousePos.y: " << io.MousePos.y << std::endl;
        std::cout << "    screenCoords.x: " << screenCoords.x << " screenCoords.y: " << screenCoords.y << std::endl;
        std::cout << "    worldRay: ";
        polyscope::operator<<(std::cout, worldRay) << std::endl;
        std::cout << "    worldPos: ";
        polyscope::operator<<(std::cout, worldPos) << std::endl;
        // if (pickPair.first == nullptr) {
        //     std::cout << "    structure: "
        //               << "none" << std::endl;
        // } else {
        //     std::cout << "    structure: " << pickPair.first << " element id: " << pickPair.second << std::endl;
        // }
    }
}

void callback() {
    ImGui::Begin("Curve Parameters");
    ImGuiIO& io = ImGui::GetIO();
    ImGui::SliderFloat("Alpha", &alpha, 0.0f, 10.0f);
    if (ImGui::Checkbox("Edit mode", &edit)) {
        std::cout << "Begin Select " << edit << std::endl;
        view::setNavigateStyle(NavigateStyle::None);
        mouse_unfc(io);
    };
    ImGui::Text("G for move, E for extrude, R for rotate");
    if (ImGui::IsKeyPressed(ImGuiKey_Tab)) {
        edit = !edit;
        std::cout << "Edit Mode" << std::endl;
        view::setNavigateStyle(NavigateStyle::None);
        mouse_unfc(io);
    }
    if (ImGui::IsKeyPressed(ImGuiKey_G)) {
        std::cout << "Begin Move " << std::endl;
    }
    if (ImGui::IsKeyPressed(ImGuiKey_E)) {
        std::cout << "Extrude Once" << std::endl;
    }
    if (ImGui::IsKeyPressed(ImGuiKey_R)) {
        std::cout << "Rotate" << std::endl;
    }

    ImGui::End();
}

// void myKeyboardFunction(unsigned int key, int action, int mods) {
//     // 检查键盘事件类型，GLFW_PRESS, GLFW_RELEASE, GLFW_REPEAT
//     if (action == GLFW_PRESS) {
//         switch (key) {
//         case GLFW_KEY_A:
//             std::cout << "A key pressed" << std::endl;
//             break;
//         case GLFW_KEY_B:
//             std::cout << "B key pressed" << std::endl;
//             break;
//         default:
//             break;
//         }
//     }
// }

int main() {
    polyscope::init();
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
    polyscope::state::userCallback = [&]() { callback(); };
    // polyscope::state::us = myKeyboardFunction;

    point p0 = {0.0f, 0.0f, 0.0f};
    point p1 = {1.0f, 0.0f, 0.0f};
    points = {p0, p1};

    polyscope::PointCloud* psCloud = polyscope::registerPointCloud("really great points", points);
    psCloud->setPointRadius(0.05);

    curveEdgesId.push_back({0, 1});
    CurveNetwork* psCurve = registerCurveNetwork("really great curve", points, curveEdgesId);

    polyscope::show();
    return 0;
}
