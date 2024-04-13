#include "Curve.h"
#include <algorithm>
#include <iostream>
// #ifdef USE_SUITESPARSE
#include <cholmod.h>
#include <cstddef>
// #endif
// using namespace cholmod;

using point = dvec3;

void Curve::cal_tangent() {
    for (size_t i = 0; i < num_edges; i++) {
        tangent_on_edges[i] = normalize(curveEdges[i]);
        // std::cout << "tangent on edges " << i << ": " << tangent_on_edges[i] << std::endl;
    }
    auto u = curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
    // std::cout << u << std::endl;
}

void Curve::update_bishop() {
    for (size_t i = 0; i < num_edges; i++) {
        point t_i = tangent_on_edges[i];
        point t_i_1 = tangent_on_edges[(i - 1) % num_edges];

        point axis = cross(t_i_1, t_i);
        double theta = std::atan2(length(axis), dot(t_i_1, t_i));
        dmat3 rotationMatrix = rotate(dmat4(1), theta, axis);

        if (theta > 1e-10) {
            normal_on_edges[i] = rotationMatrix * normal_on_edges[(i - 1) % num_edges];
        } else {
            normal_on_edges[i] = normal_on_edges[(i - 1) % num_edges];
        }
    }
    curve->addEdgeVectorQuantity("normal on edges", normal_on_edges);

    for (size_t i = 0; i < num_edges; i++) {
        binormal_on_edges[i] = normalize(cross(tangent_on_edges[i], normal_on_edges[i]));
    }
    curve->addEdgeVectorQuantity("binormal on edges", binormal_on_edges);
}

void Curve::update_material_frame() {
    for (size_t i = 0; i < num_edges; i++) {
        double curtwist = totaltwist * arc_length[i] / totallength;
        dmat3 rotationMatrix = rotate(dmat4(1.0f), curtwist, tangent_on_edges[i]);
        material_u[i] = rotationMatrix * normal_on_edges[i];
        material_v[i] = rotationMatrix * binormal_on_edges[i];
    }
    curve->addEdgeVectorQuantity("material_u", material_u);
    curve->addEdgeVectorQuantity("material_v", material_v);
}

void Curve::cal_attrs() {
    for (size_t i = 0; i < num_controlpoints; i++) {
        point Ledge = curveEdges[(i - 1 + num_controlpoints) % num_controlpoints];
        point Redge = curveEdges[(i) % num_controlpoints];

        vertex_weight[i] = 0.5 * (length(Ledge) + length(Redge));

        darboux[i] = 2.0 * cross(Ledge, Redge) / (length(Ledge) * length(Redge) + dot(Ledge, Redge));
        // std::cout << "darboux at " << i << ": " << darboux[i] << std::endl;
        curvature[i] = sqrt(dot(darboux[i], darboux[i]));
        // std::cout << "curvature at " << i << ": " << curvature[i] << std::endl;
    }
    curve->addNodeScalarQuantity("vertex_weight", vertex_weight);
    curve->addNodeScalarQuantity("curvature", curvature);
    curve->addNodeVectorQuantity("darboux", darboux);
}


void Curve::cal_bending_force() {
    bending_force.clear();
    bending_force.resize(num_controlpoints);

    for (size_t i = 0; i < num_controlpoints; i++) {
        for (size_t j = 0; j < num_controlpoints; j++) {
            point edgeR = curveEdges[j];
            point edgeL = curveEdges[(j - 1 + num_controlpoints) % num_controlpoints];

            dmat3 Rmatrix = transpoe(edgeR);
            dmat3 Lmatrix = transpoe(edgeL);

            dmat3 RT = transpose(Rmatrix);
            dmat3 LT = transpose(Lmatrix);

            double w = 0.5f * vertex_weight[j];
            auto coeff = -2 / (w * (length(edgeR) * length(edgeL) + dot(edgeR, edgeL)));
            // std::cout << "coeff " << coeff << std::endl;
            auto kbj = darboux[j];
            // dmat3 kbjT = transpose(transpoe(kbj));
            dmat3 kbjTR = outerProduct(kbj, edgeR);
            dmat3 kbjTRT = transpose(kbjTR);

            dmat3 kbjTL = outerProduct(kbj, edgeL);
            dmat3 kbjTLT = transpose(kbjTL);
            // for (int i = 0; i < 3; ++i) {
            //     for (int j = 0; j < 3; ++j) {
            //         std::cout << kbjTRT[i][j] << " ";
            //     }
            // }
            // std::cout << std::endl;
            point kbjTRTkbj = dmat3vec3(kbjTRT, kbj);
            point kbjTLTkbj = dmat3vec3(kbjTLT, kbj);
            // std::cout << " kbj " << kbjTRTkbj << std::endl;
            // std::cout << " kbj " << kbjTLTkbj << std::endl;
            point p;
            if (i == j - 1) {
                p = 2.0 * cross(edgeR, kbj) + kbjTRTkbj;
            } else if (i == j) {
                p = -1.0 * (2.0 * cross(edgeR, kbj) + kbjTRTkbj + 2.0 * cross(edgeL, kbj) - kbjTLTkbj);
            } else if (i == j + 1) {
                p = 2.0 * cross(edgeL, kbj) - kbjTLTkbj;
            }
            bending_force[i] += p * coeff;
            // std::cout << "bending force " << p << std::endl;
        }
        bending_force[i] *= -1;
        // std::cout << "bending force " << bending_force[i] << std::endl;
    }
    // std::cout << "bending force " << bending_force[0] << std::endl;
    curve->addNodeVectorQuantity("bending force", bending_force);

    // auto maxlist = std::transform(bending_force.begin(), bending_force.end(), [](point a, point b) {
    //     return length(a) < length(b);
    // });
    // auto max  = *std::max_element(bending_force.begin(), bending_force.end(), [](point a, point b) {
    //     return length(a) > length(b);
    // });
    // std::cout << "max bending force: " << max << std::endl;
    // std::cout << "bending force done" << std::endl;
}

dmat3 transpoe(point v) {
    return dmat3(v.x, 0, 0, v.y, 0, 0, v.z, 0, 0);
}

void Curve::cal_twisting_force() {
    twisting_force.clear();
    twisting_force.resize(num_controlpoints);

    for (size_t i = 0; i < num_controlpoints; i++) {
        auto prevd = 0.5 * darboux[i] / edge_length[(i - 1 + num_controlpoints) % num_controlpoints];
        auto nextd = -0.5 * darboux[i] / edge_length[i];
        twisting_force[i] = -1 * totaltwist * (nextd + prevd + eps) / (0.5f * totallength);
        // std::cout << "twisting force " << twisting_force[i] << std::endl;
    }
    curve->addNodeVectorQuantity("twisting force", twisting_force);
}


void Curve::cal_velocity() {
    acceleration.clear();
    acceleration.resize(num_controlpoints);

    // for (size_t i = 0; i < num_controlpoints; i++) {
    //     // acceleration[i] =  (bending_force[i] + twisting_force[i]) / vertex_weight[i];
    //     // velocity[i] += acceleration[i] * dt;
    //     // std::cout << "velocity " << i << ": " << velocity[i] << std::endl;
    //     // std::cout << "acceleration " << i << ": " << acceleration[i] << std::endl;
    // }


    for (size_t i = 0; i < num_controlpoints; i++) {
        acceleration[i] = (bending_force[i] + twisting_force[i]) / vertex_weight[i];
        velocity[i] += acceleration[i] * dt;
        // std::cout << "velocity " << i << ": " << acceleration[i] * dt << std::endl;
        // std::cout << "velocity " << i << ": " << velocity[i].x << " " << velocity[i].y << " " << velocity[i].z <<
        // std::endl; std::cout << "acceleration " << i << ": " << acceleration[i].x << " " << acceleration[i].y << " "
        // << acceleration[i].z << std::endl;
    }
    curve->addNodeVectorQuantity("velocity", velocity);
    curve->addNodeVectorQuantity("acceleration", acceleration);

    // std::cout << "velocity " << velocity[0]*dt << std::endl;

    newpoints.resize(num_controlpoints);

    for (size_t i = 0; i < num_controlpoints; i++) {
        newpoints[i] = controlpoints[i] + velocity[i] * dt;
        // std::cout << "newpoints at " << i << ": " << newpoints[i].x << " " << newpoints[i].y << " " << newpoints[i].z
        // << std::endl;
    }
    // controlpoints = newpoints;

    newedges.resize(num_edges);
    newedge_length.resize(num_edges);

    for (size_t i = 0; i < num_controlpoints; i++) {
        newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
        // std::cout << "newedges at " << i << ": " << newedges[i] << std::endl;
        newedge_length[i] = length(newedges[i]);
    }
    // totallength = std::accumulate(newedge_length.begin(), newedge_length.end(), 0.0f);
    // for (size_t i = 0; i < num_edges; i++) {
    //     arc_length[i] = std::accumulate(newedge_length.begin(), newedge_length.begin() + i, 0.0f);
    //     // std::cout << "arc_length " << i << ": " << arc_length[i] << std::endl;
    // }

    // curveEdges = newedges;
    // edge_length = newedge_length;
}
// #ifdef USE_SUITESPARSE
void Curve::fastsuiteprojection() {
    std::cout << "suite" << std::endl;

    cholmod_common c;
    cholmod_start(&c);
    // 100 * 1
    constraint.clear();
    constraint.resize(num_edges);

    for (size_t i = 0; i < num_edges; i++) {
        constraint[i] = pow(newedge_length[i], 2) - pow(edge_length[i], 2);
        // std::cout << "constraint at " << i << ": " << constraint[i] << std::endl;
    }

    std::vector<double> abs_con(constraint.size());
    std::transform(constraint.begin(), constraint.end(), abs_con.begin(), [](double value) { return std::abs(value); });
    auto max = *std::max_element(abs_con.begin(), abs_con.end());
    // std::cout << "max:" << max << std::endl;

    // 100 * 1
    cholmod_dense* con = cholmod_allocate_dense(num_edges, 1, num_edges, CHOLMOD_REAL, &c);
    double* matrixData = (double*)con->x;
    for (size_t i = 0; i < num_edges; ++i) {
        matrixData[i] = constraint[i];
    }
    std::cout << "good con " << std::endl;
    cholmod_triplet* M_tripletList =
        cholmod_allocate_triplet(3 * num_edges, 3 * num_edges, 3 * num_edges, 0, CHOLMOD_REAL, &c);
    for (size_t k = 0; k < 3 * num_edges; k++) {
        int idx = round(k / 3);
        ((int*)M_tripletList->i)[M_tripletList->nnz] = idx;
        ((int*)M_tripletList->j)[M_tripletList->nnz] = idx;
        ((double*)M_tripletList->x)[M_tripletList->nnz] = vertex_weight[idx];
        M_tripletList->nnz++; // Increment the number of non-zero entries
    }
    cholmod_sparse* Mass = cholmod_triplet_to_sparse(M_tripletList, 3 * num_edges, &c);

    std::cout << "good mass " << std::endl;
    while (max > 1e-10) {
        // 100 * 300
        // cholmod_triplet* triplet = cholmod_l_allocate_triplet(num_edges, 3 * nSamples, 10 * num_edges, 0, CHOLMOD_REAL, &c);
        cholmod_triplet* triplet = cholmod_allocate_triplet( num_edges, 3 * num_edges, 6 * num_edges, 0, CHOLMOD_REAL, &c);
        std::cout << " init good triplet " << std::endl;
        for (size_t k = 0; k < num_edges; k++) {
            point egde = newedges[k];
            double edge_x = egde.x;
            double edge_y = egde.y;
            double edge_z = egde.z;
            std::cout << edge_x << " " << edge_y << " " << edge_z << std::endl;

            int base_idx = k * 3;
            ((int*)triplet->i)[triplet->nnz] = k;
            ((int*)triplet->j)[triplet->nnz] = base_idx;
            ((double*)triplet->x)[triplet->nnz++] = -edge_x;
            ((long int*)triplet->i)[triplet->nnz] = k;
            ((long int*)triplet->j)[triplet->nnz] = base_idx + 1;
            ((double*)triplet->x)[triplet->nnz++] = -edge_y;
            ((long int*)triplet->i)[triplet->nnz] = k;
            ((long int*)triplet->j)[triplet->nnz] = base_idx + 2;
            ((double*)triplet->x)[triplet->nnz++] = -edge_z;

            ((long int*)triplet->i)[triplet->nnz] = k;
            ((long int*)triplet->j)[triplet->nnz] = (base_idx + 3) % (3 * nSamples);
            ((double*)triplet->x)[triplet->nnz++] = edge_x;
            ((long int*)triplet->i)[triplet->nnz] = k;
            ((long int*)triplet->j)[triplet->nnz] = (base_idx + 4) % (3 * nSamples);
            ((double*)triplet->x)[triplet->nnz++] = edge_y;
            ((long int*)triplet->i)[triplet->nnz] = k;
            ((long int*)triplet->j)[triplet->nnz] = (base_idx + 5) % (3 * nSamples);
            ((double*)triplet->x)[triplet->nnz++] = edge_z;  ((long int*)triplet->i)[triplet->nnz] = k;
        }

        cholmod_sparse* constraintGrad = cholmod_triplet_to_sparse(triplet, triplet->nnz, &c);
        std::cout << "good constraintGrad " << std::endl;

        // 300 * 100
        cholmod_sparse* conGradT = cholmod_transpose(constraintGrad, 1, &c);
        cholmod_dense* b = cholmod_sparse_to_dense(conGradT, &c);
        cholmod_factor* factor = cholmod_analyze(Mass, &c);
        cholmod_factorize(Mass, factor, &c);

        // 300 * 100
        cholmod_dense* MinvDC = cholmod_solve(CHOLMOD_A, factor, b, &c);
        cholmod_sparse* MinvDCSp = cholmod_dense_to_sparse(MinvDC, 1, &c);
        std::cout << "good MinvDC " << std::endl;

        // 100 * 100
        cholmod_sparse* DCMinvDC = cholmod_ssmult(constraintGrad, MinvDCSp, 0, 1, 1, &c);
        cholmod_factor* dLfactor = cholmod_analyze(DCMinvDC, &c);
        cholmod_factorize(DCMinvDC, dLfactor, &c);
        std::cout << "good DCMinvDC " << std::endl;

        // 100 * 1
        cholmod_dense* dLambda = cholmod_solve(CHOLMOD_A, dLfactor, con, &c);
        std::cout << "good dLambda " << std::endl;

        // 300 * 1
        // Eigen::MatrixXd dxmatrix = -1.0 * MinvDC * dLambda;
        cholmod_dense* dxmatrix = cholmod_allocate_dense(MinvDC->nrow, dLambda->ncol, MinvDC->nrow, CHOLMOD_REAL, &c);
        double alpha[2] = {1.0, 0.0}, beta[2] = {0.0, 0.0};
        cholmod_sdmult(MinvDCSp, 0, alpha, beta, dLambda, dxmatrix, &c);
        std::cout << "good dxmatrix " << std::endl;

        double* matrixData = (double*)dxmatrix->x;
        for (size_t i = 0; i < num_controlpoints; i++) {
            point p = point(matrixData[i * 3], matrixData[i * 3 + 1], matrixData[i * 3 + 2]);
            newpoints[i] += p;
        }

        for (size_t i = 0; i < num_controlpoints; i++) {
            newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
            newedge_length[i] = length(newedges[i]);
        }
        for (size_t i = 0; i < num_controlpoints; i++) {
            constraint[i] = pow(newedge_length[i], 2) - pow(edge_length[i], 2);
            matrixData[i] = constraint[i];
        }

        std::transform(constraint.begin(), constraint.end(), abs_con.begin(),
                       [](double value) { return std::abs(value); });
        max = *std::max_element(abs_con.begin(), abs_con.end());
    }


    controlpoints = newpoints;
    curveEdges = newedges;
    edge_length = newedge_length;

    totallength = std::accumulate(edge_length.begin(), edge_length.end(), 0.0f);
    for (size_t i = 0; i < nSamples; i++) {
        arc_length[i] = std::accumulate(edge_length.begin(), edge_length.begin() + i, 0.0f);
    }

    for (size_t i = 0; i < num_edges; i++) {
        velocity[i] = (newpoints[i] - controlpoints[i]) / dt;
    }
}
// #else
void Curve::fastprojection() {
     std::cout << "eigen" << std::endl;
    // 100 * 1
    constraint.clear();
    constraint.resize(num_edges);

    for (size_t i = 0; i < num_edges; i++) {
        constraint[i] = pow(newedge_length[i], 2) - pow(edge_length[i], 2);
        // std::cout << "constraint at " << i << ": " << constraint[i] << std::endl;
    }

    std::vector<double> abs_con(constraint.size());
    std::transform(constraint.begin(), constraint.end(), abs_con.begin(), [](double value) { return std::abs(value);
    }); auto max = *std::max_element(abs_con.begin(), abs_con.end()); std::cout << "0max:" << max << std::endl;

    // 100 * 1
    Eigen::MatrixXd constraintMatrix(num_edges, 1);
    for (size_t i = 0; i < num_edges; i++) {
        constraintMatrix(i, 0) = constraint[i];
    }

    Eigen::SparseMatrix<double> Mass(3 * nSamples, 3 * nSamples);
    std::vector<Eigen::Triplet<double>> M_tripletList;
    for (size_t i = 0; i < num_edges; i++) {
        auto w = vertex_weight[i];
        M_tripletList.push_back(Eigen::Triplet<double>(i * 3, i * 3, w));
        M_tripletList.push_back(Eigen::Triplet<double>(i * 3 + 1, i * 3 + 1, w));
        M_tripletList.push_back(Eigen::Triplet<double>(i * 3 + 2, i * 3 + 2, w));
    }
    Mass.setFromTriplets(M_tripletList.begin(), M_tripletList.end());
    // int count = 0;
    while (max > 1e-10) {
        // 100 * 300
        Eigen::SparseMatrix<double> constraintGrad(nSamples, 3 * nSamples);
        std::vector<Eigen::Triplet<double>> tripletList;
        for (size_t i = 0; i < num_edges; i++) {
            auto edge = newedges[i];
            // std::cout << "edge" << edge << std::endl;
            edge *= 2;
            tripletList.push_back(Eigen::Triplet<double>(i, i * 3, -edge.x));
            tripletList.push_back(Eigen::Triplet<double>(i, i * 3 + 1, -edge.y));
            tripletList.push_back(Eigen::Triplet<double>(i, i * 3 + 2, -edge.z));
            tripletList.push_back(Eigen::Triplet<double>(i, (i * 3 + 3) % (3 * num_edges), edge.x));
            tripletList.push_back(Eigen::Triplet<double>(i, (i * 3 + 4) % (3 * num_edges), edge.y));
            tripletList.push_back(Eigen::Triplet<double>(i, (i * 3 + 5) % (3 * num_edges), edge.z));
        }
        constraintGrad.setFromTriplets(tripletList.begin(), tripletList.end());
        // 100 * 300
        Eigen::SparseMatrix<double> constraintGradT = constraintGrad.transpose();

        Eigen::SparseLU<Eigen::SparseMatrix<double>> MinvDCsolver;
        MinvDCsolver.compute(Mass);
        if (MinvDCsolver.info() != Eigen::Success) {
            std::cout << "MinvDCsolver decomposition failed" << std::endl;
        }

        // 300 * 100
        Eigen::SparseMatrix<double> MinvDC = MinvDCsolver.solve(constraintGradT);
        if (MinvDCsolver.info() != Eigen::Success) {
            std::cout << "MinvDCsolver solve failed" << std::endl;
        }
        // 100 * 100
        Eigen::SparseMatrix<double> DCMinvDC = constraintGrad * MinvDC;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> dLambdasolver;
        dLambdasolver.compute(DCMinvDC);
        if (dLambdasolver.info() != Eigen::Success) {
            std::cout << "dLambdasolver decomposition failed" << std::endl;
        }

        for (size_t i = 0; i < num_edges; i++) {
            constraintMatrix(i, 0) = constraint[i];
        }
        // 100 *1
        Eigen::MatrixXd dLambda = dLambdasolver.solve(constraintMatrix);
        if (MinvDCsolver.info() != Eigen::Success) {
            std::cout << "dLambdasolver solve failed" << std::endl;
        }

        // 300 * 1
        Eigen::MatrixXd dxmatrix = -1.0 * MinvDC * dLambda;

        // reshape
        for (size_t i = 0; i < num_controlpoints; i++) {
            point p = point(dxmatrix(i * 3), dxmatrix(i * 3 + 1), dxmatrix(i * 3 + 2));
            newpoints[i] += p;
            // std::cout << "dx"<< p.x << " " << p.y << " " << p.z << std::endl;
        }

        for (size_t i = 0; i < num_controlpoints; i++) {
            newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
            // std::cout << "newedges at " << i << ": " << newedges[i] << std::endl;
            newedge_length[i] = length(newedges[i]);
        }
        for (size_t i = 0; i < num_controlpoints; i++) {
            constraint[i] = pow(newedge_length[i], 2) - pow(edge_length[i], 2);
            // std::cout << "constraint at " << i << ": " << constraint[i] << std::endl;
        }
        std::transform(constraint.begin(), constraint.end(), abs_con.begin(),
                       [](double value) { return std::abs(value); });
        max = *std::max_element(abs_con.begin(), abs_con.end());
        // std::cout << "max:" << max << std::endl;
    }


    controlpoints = newpoints;
    curveEdges = newedges;
    edge_length = newedge_length;

    totallength = std::accumulate(edge_length.begin(), edge_length.end(), 0.0f);
    for (size_t i = 0; i < nSamples; i++) {
        arc_length[i] = std::accumulate(edge_length.begin(), edge_length.begin() + i, 0.0f);
        // std::cout << "arc_length " << i << ": " << arc_length[i] << std::endl;
    }

    for (size_t i = 0; i < num_edges; i++) {
        velocity[i] = (newpoints[i] - controlpoints[i]) / dt;
    }
}
// #endif