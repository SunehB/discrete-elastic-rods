#include "polyscope/polyscope.h"
#include "include/Curve.h"
#include "imgui.h"



Curve curve;


void DERcallback(){
    ImGui::Begin("DERcallback");
    ImGui::Text("Hello, world!");
    if(ImGui::Button("Button")){
        curve.tst_gd();
    }
    ImGui::End();
}

int main() {

    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
   
    polyscope::init();
    
    curve.init();

    polyscope::state::userCallback = DERcallback;
    // [&]() {
    //     curve.animate();
    // };

    polyscope::show(); 
    return 0;
}

