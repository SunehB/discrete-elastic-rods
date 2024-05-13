#include "imgui.h"
#include "include/Curve.h"
#include "polyscope/polyscope.h"

static bool isRunning = false;
static int currentIteration = 0;
static float alpha = 1.f;
static float beta = 1.f;
static float dt = 1e-3f;
static float totaltwist = 2 * PI;
static int nSample = 100;
static int nIteration = 100;
static const char * timestep = "1e-3\01e-4\01e-5\01e-6\0";
static int cur_dt = 0;
static Curve cur;

void DERcallback() {
    ImGui::Begin("Mitchelle Buckling Curve");
    if(ImGui::InputInt("Set Sample Number", &nSample)){
        if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        cur.nSamples = nSample;
    
    };
    if (ImGui::Button("Init closed Curve")) {
        if (cur.status()){
            polyscope::warning("Dont init curve twice!");
            return;
        }
        cur.initcurve();
    }
     if (ImGui::Button("Init open Curve")) {
        if (cur.status()){
            polyscope::warning("Dont init curve twice!");
            return;
        }   
        cur.closed = false;
        cur.init_catenary();
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
    if (ImGui::SliderFloat("Set twisting modulus (beta)", &beta, 1.0f, 100.0f)) {
        if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        cur.set_beta(beta);
    }
    if (ImGui::InputFloat("Set twist angle", &totaltwist, 1.0f, 100.0f)) {
         if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        cur.set_totaltwist(totaltwist);
    }
    if(ImGui::Combo("Time Setp", &cur_dt, timestep, 4)){
        if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        switch (cur_dt){
            case 0:
                dt = 1e-3f;
                break;
            case 1:
                dt = 1e-4f;
                break;
            case 2:
                dt = 1e-5f;
                break;
            case 3:
                dt = 1e-6f;
                break;
        }
        cur.set_dt(dt);
    }
    ImGui::InputInt("Set Iteration number", &nIteration);
    if (ImGui::Button("Run Loop")) {
        if (!cur.status()){
            polyscope::warning("Init curve first!");
            return;
        }
        isRunning = true;
        currentIteration = 0;
    }

    if (isRunning && currentIteration < nIteration) {
        cur.loop();
        currentIteration++;
        if (currentIteration >= nIteration) {
            isRunning = false;
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
