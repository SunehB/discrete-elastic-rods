#include "polyscope/polyscope.h"
#include "include/Curve.h"
#include "imgui.h"


Curve sprial;
Curve straight_line;

void DERcallback(){
    ImGui::Begin("DERcallback");
    ImGui::Text("Hello, wsssssorld!");
    if(ImGui::Button("test_run")){
        sprial.test_run();
    }
    static float newtwist = 0.0;
    if(ImGui::SliderFloat("twistonline", &newtwist, 0, 4*PI)){
        sprial.adjusttwist(newtwist);
    }

    ImGui::End();
}


int main() {

    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
   
    polyscope::init();

    Curve sprial;
    sprial.init_sprial();

    Curve straight_line;
    straight_line.init_straight_line();

    sprial.test_run();
    straight_line.test_run();

    polyscope::state::userCallback = DERcallback;
    // [&]() {
    //     curve.animate();
    // };

    polyscope::show(); 
    return 0;
}



