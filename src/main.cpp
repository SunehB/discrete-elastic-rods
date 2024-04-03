#include "imgui.h"
#include "include/Curve.h"
#include "polyscope/polyscope.h"


Curve sprial;
Curve straight_line;

void DERcallback(Curve &cur ){
    ImGui::Begin("DERcallback");
    ImGui::Text("Hello, wsssssorld!");
    if(ImGui::Button("test_run")){
        cur.symEuler();
        cur.test_run();
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

    Curve sprial;
    sprial.init_sprial();
    sprial.test_run();

    // Curve straight_line;
    // straight_line.init_straight_line();
    // straight_line.test_run();

    // 使用lambda表达式设置回调，以捕获并使用上面的实例
    polyscope::state::userCallback = [&]() {
        DERcallback(sprial);
    };

    polyscope::show();
    return 0;
}
