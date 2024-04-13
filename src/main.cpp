#include "imgui.h"
#include "include/Curve.h"
#include "polyscope/polyscope.h"

static bool isRunning = false;
static int currentIteration = 0;
static float alpha = 1.f;
static float beta = 1.f;
static float dt = 0.001f;
static float totaltwist = 0.0f;
static int nSample = 100;
static int nIteration = 100;
static Curve cur;

void DERcallback() {
    ImGui::Begin("11Michelle Buckling Curve D");
    ImGui::InputInt("Set Sample Number", &nSample);
    if (ImGui::Button("Init Curve")) {
        if (cur.status()){
            polyscope::warning("Dont init curve twice!");
            return;
        }
        cur.initcurve(nSample);
    }
    ImGui::SameLine();
    if (ImGui::Button("Reset Curve")) {
        if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        cur.reset();
    }
    if (ImGui::SliderFloat("Set bending modulus (alpha)", &alpha, 0.0f, 1.0f)) {
         if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        cur.set_alpha(alpha);
    }
    if (ImGui::SliderFloat("Set twisting modulus (beta)", &beta, 0.0f, 1.0f)) {
        if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        cur.set_beta(beta);
    }
    if (ImGui::InputFloat("Set twist angle", &totaltwist, 0.0f, 1.0f)) {
         if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        cur.set_totaltwist(totaltwist);
    }
    if (ImGui::SliderFloat("Set time stamp (dt)", &dt, 0.0f, 1.0f)) {
         if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        cur.set_dt(dt);
    }
    ImGui::InputInt("Set Iteration number", &nIteration);
    if (ImGui::Button("Run Loop")) {
        // std::cout <<cur.status()<<std::endl;
        if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        isRunning = true;
        currentIteration = 0;
    }

    if (isRunning && currentIteration < nIteration) {
        cur.loop();
        // std::cout << "loop " << currentIteration << " done" << std::endl;
        currentIteration++;
        if (currentIteration >= nIteration) {
            isRunning = false;
            // std::cout << "All loops done." << std::endl;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Pause/Resume")) {
        isRunning = !isRunning;
    };
    ImGui::End();
}

int main() {
    polyscope::init();
    polyscope::view::setUpDir(polyscope::UpDir::ZUp);
    polyscope::view::setFrontDir(polyscope::FrontDir::NegYFront);
    polyscope::state::userCallback = [&]() { DERcallback(); };
    polyscope::show();
    return 0;
}
