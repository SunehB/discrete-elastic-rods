#include "polyscope/polyscope.h"
#include "include/Curve.h"
#include "imgui.h"


void DERcallback(){
    ImGui::Begin("DERcallback");
    ImGui::Text("Hello, world!");
    ImGui::End();
}

int main() {

    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
   
    polyscope::init();
    Curve curve;
    curve.init();

    // polyscope::state::userCallback =  [&]() {
    //     curve.animate();
    // };

    polyscope::show(); 
    return 0;
}

