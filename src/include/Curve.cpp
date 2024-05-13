#include "Curve.h"
// #ifdef USE_SUITESPARSE
#include <cholmod.h>
// #endif


void Curve::cal_tangent() {
    last_tangent  = tangent_on_edges[0];
    for (size_t i = 0; i < num_edges; i++) {
        tangent_on_edges[i] = normalize(curveEdges[i]);
    }
    auto u = curve->addEdgeVectorQuantity("tangent on edges", tangent_on_edges);
}

void Curve::update_bishop() {

    point timePtAxis = cross(tangent_on_edges[0], last_tangent);
    double timePtTheta = std::atan2(length(timePtAxis), dot(tangent_on_edges[0], tangent_on_edges[1]));
    dmat3 time_pt_rotationMatrix = rotate(dmat4(1), timePtTheta, normalize(timePtAxis));
    normal_on_edges[0] = time_pt_rotationMatrix * normalize(cross(tangent_on_edges[0], reference));
    binormal_on_edges[0] = time_pt_rotationMatrix * normalize(cross(tangent_on_edges[0], normal_on_edges[0]));

    for (size_t i = 1; i < num_edges; i++) {
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
    // for (size_t i = 0; i < num_controlpoints; i++) {
    //     // if (closed && i == 0) continue;
    //     // if (closed && i == num_controlpoints - 1) continue;
    //     point Ledge = curveEdges[(i - 1 + num_controlpoints) % num_controlpoints];
    //     point Redge = curveEdges[(i) % num_controlpoints];

    //     vertex_weight[i] = 0.5 * (length(Ledge) + length(Redge));
    //     darbouxdenom[i] = length(Ledge) * length(Redge) + dot(Ledge, Redge);
    //     darboux[i] = 2.0 * cross(Ledge, Redge) / darbouxdenom[i];
    //     curvature[i] = sqrt(dot(darboux[i], darboux[i]));
    // }
    for (size_t i = 1; i < num_controlpoints - 1; i++) {
        // if (closed && i == 0) continue;
        // if (closed && i == num_controlpoints - 1) continue;
        point Ledge, Redge;

        Ledge = curveEdges[i - 1];
        Redge = curveEdges[i];

        vertex_weight[i] = 0.5 * (length(Ledge) + length(Redge));
        darbouxdenom[i] = length(Ledge) * length(Redge) + dot(Ledge, Redge);
        darboux[i] = 2.0 * cross(Ledge, Redge) / darbouxdenom[i];
        curvature[i] = sqrt(dot(darboux[i], darboux[i]));
    }
    vertex_weight[0] = 0.5 * (length(curveEdges[1]));

    curve->addNodeScalarQuantity("vertex_weight", vertex_weight);
    curve->addNodeScalarQuantity("curvature", curvature);
    curve->addNodeVectorQuantity("darboux", darboux);
}


void Curve::cal_bending_force() {
    bending_force.clear();
    bending_force.resize(num_controlpoints);
    for (size_t j = 1; j < num_controlpoints - 1; j++) {
        // size_t cur_idx = j;
        // size_t prev_idx = (j - 1 + num_controlpoints) % num_controlpoints;
        // size_t next_idx = (j + 1) % num_controlpoints;

        size_t cur_idx = j;
        size_t prev_idx = j - 1;
        size_t next_idx = j + 1;

        auto cur = curveEdges[cur_idx];
        auto prev = curveEdges[prev_idx];
        auto next = curveEdges[next_idx];

        double denom = darbouxdenom[cur_idx];

        double l = 2 * vertex_weight[prev_idx];
        point kb = darboux[prev_idx];
        dmat3 grad_prev = (2.0 * skew(cur) + outerProduct(cur, kb)) / denom;
        point bend_prev = -2 / l * (transpose(grad_prev) * kb);
        if (j == 1) bend_prev = point(0, 0, 0);

        l = 2 * vertex_weight[next_idx];
        kb = darboux[next_idx];
        auto grad_next = (2.0 * skew(prev) - outerProduct(prev, kb)) / denom;
        auto bend_next = -2 / l * (transpose(grad_next) * kb);
        if (j == num_controlpoints - 2) bend_next = point(0, 0, 0);

        l = 2 * vertex_weight[cur_idx];
        kb = darboux[cur_idx];
        auto grad_cur = -(grad_prev + grad_next);
        auto bend_cur = -2 / l * (transpose(grad_cur) * kb);

        bending_force[cur_idx] =  (bend_prev + bend_cur + bend_next);
        // std::cout << "bending force at " << j << ": " << bending_force[cur_idx][0] << " " << bending_force[cur_idx][1]
        //           << " " << bending_force[cur_idx][2] << std::endl;
    }
    curve->addNodeVectorQuantity("bending force", bending_force);
}


void Curve::cal_twisting_force() {
    twisting_force.clear();
    twisting_force.resize(num_controlpoints);

    // for (size_t i = 1; i < num_controlpoints-1; i++) {
    //     auto curr_idx = i;
    //     auto kb = darboux[curr_idx];

    //     auto grad_prev = 0.5 * kb / edge_length[curr_idx];
    //     auto grad_next = -0.5 * kb / edge_length[i];

    //     auto grad_cur = -(grad_prev + grad_next);

    //     twisting_force[curr_idx] = (grad_prev  + grad_next) * totaltwist / (0.5f * totallength);
    // }
    // curve->addNodeVectorQuantity("twisting force", twisting_force);
}


void Curve::cal_velocity() {
    acceleration.clear();
    acceleration.resize(num_controlpoints);

    for (size_t i = 1; i < num_controlpoints - 1; i++) {
        acceleration[i] = (bending_force[i] + twisting_force[i]) / vertex_weight[i];
        velocity[i] += acceleration[i] * dt;
    }
    curve->addNodeVectorQuantity("velocity", velocity);
    curve->addNodeVectorQuantity("acceleration", acceleration);

    for (size_t i = 0; i < num_controlpoints; i++) {
        newpoints[i] = controlpoints[i] + velocity[i] * dt;
    }
    for (size_t i = 0; i < num_edges; i++) {
        newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
        newedge_length[i] = length(newedges[i]);
    }
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
        // cholmod_triplet* triplet = cholmod_l_allocate_triplet(num_edges, 3 * nSamples, 10 * num_edges, 0,
        // CHOLMOD_REAL, &c);
        cholmod_triplet* triplet =
            cholmod_allocate_triplet(num_edges, 3 * num_edges, 6 * num_edges, 0, CHOLMOD_REAL, &c);
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
            ((double*)triplet->x)[triplet->nnz++] = edge_z;
            ((long int*)triplet->i)[triplet->nnz] = k;
        }

        cholmod_sparse* constraintGrad = cholmod_triplet_to_sparse(triplet, triplet->nnz, &c);
        // std::cout << "good constraintGrad " << std::endl;

        // 300 * 100
        cholmod_sparse* conGradT = cholmod_transpose(constraintGrad, 1, &c);
        cholmod_dense* b = cholmod_sparse_to_dense(conGradT, &c);
        cholmod_factor* factor = cholmod_analyze(Mass, &c);
        cholmod_factorize(Mass, factor, &c);

        // 300 * 100
        cholmod_dense* MinvDC = cholmod_solve(CHOLMOD_A, factor, b, &c);
        cholmod_sparse* MinvDCSp = cholmod_dense_to_sparse(MinvDC, 1, &c);
        // std::cout << "good MinvDC " << std::endl;

        // 100 * 100
        cholmod_sparse* DCMinvDC = cholmod_ssmult(constraintGrad, MinvDCSp, 0, 1, 1, &c);
        cholmod_factor* dLfactor = cholmod_analyze(DCMinvDC, &c);
        cholmod_factorize(DCMinvDC, dLfactor, &c);
        // std::cout << "good DCMinvDC " << std::endl;

        // 100 * 1
        cholmod_dense* dLambda = cholmod_solve(CHOLMOD_A, dLfactor, con, &c);
        // std::cout << "good dLambda " << std::endl;

        // 300 * 1
        // Eigen::MatrixXd dxmatrix = -1.0 * MinvDC * dLambda;
        cholmod_dense* dxmatrix = cholmod_allocate_dense(MinvDC->nrow, dLambda->ncol, MinvDC->nrow, CHOLMOD_REAL, &c);
        double alpha[2] = {1.0, 0.0}, beta[2] = {0.0, 0.0};
        cholmod_sdmult(MinvDCSp, 0, alpha, beta, dLambda, dxmatrix, &c);
        // std::cout << "good dxmatrix " << std::endl;

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
    std::cout << "0max:" << max << std::endl;

    // 100 * 1
    Eigen::MatrixXd constraintMatrix(num_edges, 1);
    for (size_t i = 0; i < num_edges; i++) {
        constraintMatrix(i, 0) = constraint[i];
    }

    Eigen::SparseMatrix<double> Mass(3 * (num_edges), 3 * (num_edges));
    std::vector<Eigen::Triplet<double>> M_tripletList;
    for (size_t i = 0; i < num_edges; i++) {
        auto w = vertex_weight[i];
        // std::cout << "w: " << w << std::endl;
        M_tripletList.push_back(Eigen::Triplet<double>(i * 3, i * 3, w));
        M_tripletList.push_back(Eigen::Triplet<double>(i * 3 + 1, i * 3 + 1, w));
        M_tripletList.push_back(Eigen::Triplet<double>(i * 3 + 2, i * 3 + 2, w));
    }
    Mass.setFromTriplets(M_tripletList.begin(), M_tripletList.end());

    while (max > 1e-10) {
        // 100 * 300
        Eigen::SparseMatrix<double> constraintGrad(num_edges, 3 * num_edges);
        std::vector<Eigen::Triplet<double>> tripletList;
        for (size_t i = 0; i < num_edges; i++) {
            auto edge = newedges[i];
            edge *= 2;
            tripletList.push_back(Eigen::Triplet<double>(i, i * 3, -edge.x));
            tripletList.push_back(Eigen::Triplet<double>(i, i * 3 + 1, -edge.y));
            tripletList.push_back(Eigen::Triplet<double>(i, i * 3 + 2, -edge.z));
            tripletList.push_back(Eigen::Triplet<double>(i, (i * 3 + 3) % (3 * num_edges), edge.x));
            tripletList.push_back(Eigen::Triplet<double>(i, (i * 3 + 4) % (3 * num_edges), edge.y));
            tripletList.push_back(Eigen::Triplet<double>(i, (i * 3 + 5) % (3 * num_edges), edge.z));
        }
        constraintGrad.setFromTriplets(tripletList.begin(), tripletList.end());
        // std::cout << constraintGrad << std::endl;
        // 100 * 300
        Eigen::SparseMatrix<double> constraintGradT = constraintGrad.transpose();

        Eigen::SparseLU<Eigen::SparseMatrix<double>> MinvDCsolver;
        MinvDCsolver.analyzePattern(Mass);
        MinvDCsolver.factorize(Mass);
        MinvDCsolver.compute(Mass);
        if (MinvDCsolver.info() != Eigen::Success) {
            std::cout << "Factorization failed: " << MinvDCsolver.lastErrorMessage() << std::endl;
        }

        // 300 * 100
        Eigen::SparseMatrix<double> MinvDC = MinvDCsolver.solve(constraintGradT);
        if (MinvDCsolver.info() != Eigen::Success) {
            std::cout << "MinvDCsolver solve failed" << std::endl;
        }
        // 100 * 100
        Eigen::SparseMatrix<double> DCMinvDC = constraintGrad * MinvDC;
        // std::cout << DCMinvDC << std::endl;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> dLambdasolver;
        dLambdasolver.analyzePattern(DCMinvDC);
        dLambdasolver.factorize(DCMinvDC);
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
        // std::cout << "dLambda: " << dLambda.size() << std::endl;
        // 300 * 1
        Eigen::MatrixXd dxmatrix = -1.0 * MinvDC * dLambda;
        // std::cout << "dxmatrix: " << dxmatrix.size() << std::endl;

        // reshape
        for (size_t i = 0; i < num_controlpoints - 1; i++) {
            point p = point(dxmatrix(i * 3), dxmatrix(i * 3 + 1), dxmatrix(i * 3 + 2));
            newpoints[i] += p;
        }

        for (size_t i = 0; i < num_edges; i++) {
            newedges[i] = newpoints[(i + 1) % num_controlpoints] - newpoints[i];
            newedge_length[i] = length(newedges[i]);
        }
        for (size_t i = 0; i < num_edges; i++) {
            constraint[i] = pow(newedge_length[i], 2) - pow(edge_length[i], 2);
        }
        std::transform(constraint.begin(), constraint.end(), abs_con.begin(),
                       [](double value) { return std::abs(value); });
        max = *std::max_element(abs_con.begin(), abs_con.end());
    }

    controlpoints = newpoints;
    curveEdges = newedges;
    edge_length = newedge_length;

    totallength = std::accumulate(edge_length.begin(), edge_length.end(), 0.0f);

    for (size_t i = 0; i < num_edges; i++) velocity[i] = (newpoints[i] - controlpoints[i]) / dt;
}
// #endif