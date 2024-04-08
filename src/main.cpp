#include "imgui.h"
#include "include/Curvenew.h"
#include "polyscope/polyscope.h"


// Curve sprial;
// Curve straight_line;

void DERcallback(Curven& cur) {
    ImGui::Begin("DERcallback");
    ImGui::Text("Hello, wsssssorld!");
    if (ImGui::Button("test_run")) {
        for (int i = 0; i < 1; i++) {
            // cur.symEuler();
            // // cur.manifoldProjection();
            // cur.test_run();
            cur.loop();
        }
    }

    // static float newtwist = 0.0;
    // if(ImGui::SliderFloat("twistonline", &newtwist, 0, 4*PI)){
    //     cur.adjusttwist(newtwist);
    // }

    ImGui::End();
}


int main() {

    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);


    polyscope::init();

    // Curve sprial;
    // sprial.init_sprial();
    // sprial.test_run();

    // Curve straight_line;
    // straight_line.init_straight_line();
    // straight_line.test_run();


    Curven curve;
    curve.initcurve();
    // curve.test_run();

    polyscope::state::userCallback = [&]() { DERcallback(curve); };

    polyscope::show();
    return 0;
}
