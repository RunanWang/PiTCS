//
// Created by 王润安 on 2024/1/11.
//

#include "testHeader.h"
#include "../Algorithm/PivotalTruss.h"
#include "../Algorithm/Algo.h"
#include "../Utils/logging.h"
#include "../Utils/Timer.h"
#include "../constant.h"
#include "stdlib.h"
#include "queue"
#include "algorithm"

// PiTCS using PiTIndex
double evalGroundTruthCommunityByF1(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, NODE_TYPE query,
                                    double &bestden, vector<NODE_TYPE> &gt) {
    double bestf1 = 0;
    bestden = 0;
    EDGE_TYPE best_k = 0;
    LAYER_TYPE bestComboId = 0;
    vector<vector<NODE_TYPE> > qus;
    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        EDGE_TYPE kmax = getMaxTrussness(mlg, index[i], query);
        for (EDGE_TYPE k = 0; k <= kmax; k++) {
            queryNodeCommunityAllowZero(mlg, index[i], query, k, qus);
            for (const auto &qu: qus) {
                double f1 = f1Score(qu, gt);
                // double den = density(qu, mlg, DEN_BETA);
                // double den = relativeDensity(qu, mlg, DEN_BETA);
                double den = jaccard(qu, gt);
                if (f1 > bestf1) {
                    bestComboId = i;
                    bestf1 = f1;
                    bestden = den;
                    best_k = k;
                }
            }
        }
    }
    if (ci->isPair) {
        LOG(INFO) << fmt::format("Best community of vertex {} is pivot-lamb-k {}-{}-{} (rank={}) ({:.2f})", query,
                                 ci->pivotalList[bestComboId].pivotLayer, ci->order[bestComboId][best_k].lamb,
                                 ci->order[bestComboId][best_k].k, best_k, bestf1);
    } else {
        LOG(INFO) << fmt::format("Best community of vertex {} is {}-{} ({:.2f})", query,
                                 ci->pivotalList[bestComboId].toString(),
                                 best_k, bestf1);
    }
    return bestf1;
}

void connect(MLGWithSchema *mlg, vector<EDGE_TYPE> &subgraph, vector<vector<EDGE_TYPE> > &res, NODE_TYPE query) {
    NODE_TYPE u, v, w;
    EDGE_TYPE e, e1, e2;
    auto neighbors = new NODE_TYPE[mlg->sumGraph->maxDegree];
    NODE_TYPE neighborSize = 0;
    set<EDGE_TYPE> sg(subgraph.begin(), subgraph.end());
    auto queryNodeNeighborStart = mlg->sumGraph->start(query);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(query);
    set<EDGE_TYPE> edgesAttachToQueryNode;
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        set<EDGE_TYPE> currentEdgesInC;
        // 从一个连接查询点的边作为出发
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(query, neighbor);
        if (sg.find(edge) == sg.end() or edgesAttachToQueryNode.find(edge) != edgesAttachToQueryNode.end()) {
            continue;
        }
        edgesAttachToQueryNode.insert(edge);
        queue<EDGE_TYPE> toInsertEdges;
        toInsertEdges.push(edge);
        while (not toInsertEdges.empty()) {
            e = toInsertEdges.front();
            toInsertEdges.pop();
            if (currentEdgesInC.find(e) != currentEdgesInC.end()) {
                continue;
            }
            edgesTriangleConnectedByGivenEdgeBasic(mlg, e, neighbors, neighborSize);
            mlg->schema->idToEdge(e, &u, &v);
            for (NODE_TYPE j = 0; j < neighborSize; j++) {
                w = neighbors[j];
                e1 = mlg->schema->edgeToId(u, w);
                e2 = mlg->schema->edgeToId(v, w);
                if (sg.find(e1) != sg.end() and sg.find(e2) != sg.end()) {
                    toInsertEdges.push(e1);
                    toInsertEdges.push(e2);
                }
                if (w == query) {
                    edgesAttachToQueryNode.insert(e1);
                    edgesAttachToQueryNode.insert(e2);
                }
            }
            currentEdgesInC.insert(e);
        }
        res.resize(res.size() + 1);
        res[res.size() - 1].assign(currentEdgesInC.begin(), currentEdgesInC.end());
    }
}

// RPiTCS using RPiTIndex
double evalGroundTruthCommunityByFuncInd(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, NODE_TYPE query,
                                         double &bestden, vector<NODE_TYPE> &gt) {
    double bestf1 = 0;
    bestden = 0;
    NODE_TYPE u, v;
    EDGE_TYPE best_k = 0;
    LAYER_TYPE bestComboId = 0;
    vector<vector<EDGE_TYPE> > qus;
    PivotalCombo c = new PivotalCombo();
    c.hasPivot = true;

    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        EDGE_TYPE rank = getMaxTrussness(mlg, index[i], query);
        EDGE_TYPE ki = ci->order[i][rank].k;
        c.lamb = ci->order[i][rank].lamb;
        c.pivotLayer = ci->pivotalList[i].pivotLayer;

        queryNodeCommunityEdge(mlg, index[i], query, rank, qus);
        for (const auto &qu: qus) {
            vector<vector<EDGE_TYPE> > ccs;
            if (ki != 0) {
                vector<EDGE_TYPE> result;
                peelRk(mlg, qu, c, ki, result);
                connect(mlg, result, ccs, query);
            } else {
                ccs.emplace_back(qu);
            }
            for (auto cc: ccs) {
                set<NODE_TYPE> nodeInCommunity;
                nodeInCommunity.insert(query);
                for (auto edge: cc) {
                    mlg->schema->idToEdge(edge, &u, &v);
                    nodeInCommunity.insert(u);
                    nodeInCommunity.insert(v);
                }
                vector<NODE_TYPE> com;
                com.assign(nodeInCommunity.begin(), nodeInCommunity.end());
                double f1 = f1Score(com, gt);
                // double den = density(qu, mlg, DEN_BETA);
                // double den = relativeDensity(qu, mlg, DEN_BETA);
                double den = jaccard(com, gt);
                if (f1 > bestf1) {
                    bestComboId = i;
                    bestf1 = f1;
                    bestden = den;
                    best_k = rank;
                }
                LOG(INFO) << fmt::format("pivot-lamb-k {}-{}-{} (rank={}) ({:.2f}) comm={}", c.pivotLayer, c.lamb, ki,
                                         rank,
                                         f1, fmt::join(com.begin(), com.end(), ", "));
            }
        }
    }
    if (ci->isPair) {
        LOG(INFO) << fmt::format("Best community of vertex {} is pivot-lamb-k {}-{}-{} (rank={}) ({:.2f})", query,
                                 ci->pivotalList[bestComboId].pivotLayer, ci->order[bestComboId][best_k].lamb,
                                 ci->order[bestComboId][best_k].k, best_k, bestf1);
    } else {
        LOG(INFO) << fmt::format("Best community of vertex {} is {}-{} ({:.2f})", query,
                                 ci->pivotalList[bestComboId].toString(),
                                 best_k, bestf1);
    }
    return bestf1;
}

double evalGroundTruthCommunity(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, NODE_TYPE query,
                                double &bestden, vector<NODE_TYPE> &gt) {
     return evalGroundTruthCommunityByF1(mlg, index, ci, query, bestden, gt);
}

double evalGroundTruthCommunityRank(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, NODE_TYPE query,
                                    double &bestden, vector<NODE_TYPE> &gt) {
    return evalGroundTruthCommunityByFuncInd(mlg, index, ci, query, bestden, gt);
}


void evalCommunityByDensity(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, NODE_TYPE query,
                            double &bestden, double &bestdiam, Timer *timer) {
    bestden = 0;
    bestdiam = 0;
    EDGE_TYPE best_k = 0;
    LAYER_TYPE bestComboId = 1;
    vector<vector<NODE_TYPE> > qus;
    vector<NODE_TYPE> bestqu;
    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        // LOG(INFO) << fmt::format("Evaluating query = {}, combo = {}", query, ci->idToComboName[i]);

        auto trussnessBound = getMaxTrussness(mlg, index[i], query);
        EDGE_TYPE begin = 1;
        for (EDGE_TYPE k = begin; k <= trussnessBound; k++) {
            timer->startTimer();
            queryNodeCommunity(mlg, index[i], query, k, qus);
            timer->endTimer();
            for (const auto &qu: qus) {
                double den = density(qu, mlg, DEN_BETA);
                // LOG(INFO) << fmt::format("Combo = {}-{}, den = {}, community = [{}]", ci->pivotalList[i].toString(), k,
                //                          den, fmt::join(qu.begin(), qu.end(), ", "));
                if (den > bestden) {
                    bestComboId = i;
                    bestden = den;
                    best_k = k;
                    bestqu.assign(qu.begin(), qu.end());
                }
            }
        }
    }
    if (ci->isPair) {
        LOG(INFO) << fmt::format("Best community of vertex {} is pivot-lamb-k {}-{}-{} (rank={}) ({:.2f})", query,
                                 ci->pivotalList[bestComboId].pivotLayer, ci->order[bestComboId][best_k].lamb,
                                 ci->order[bestComboId][best_k].k, best_k, bestden);
    } else {
        LOG(INFO) << fmt::format("Best community of vertex {} is {}-{} ({:.2f})", query,
                                 ci->pivotalList[bestComboId].toString(),
                                 best_k, bestden);
    }
    // bestdiam = double(diameterMultilayer(bestqu, mlg));
    //    LOG(INFO) << fmt::format("Diameter = {:.4f}/{}", bestdiam, bestqu.size());
    bestdiam = 1;
}


void evalCommunityByUtilityAndDenAcc(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, NODE_TYPE query,
                                     double &bestden, double &bestdiam, Timer *timer) {
    bestden = 0;
    bestdiam = 0;
    EDGE_TYPE best_k = 0;
    LAYER_TYPE bestComboId = 1;
    vector<vector<NODE_TYPE> > qus;
    vector<NODE_TYPE> bestqu;
    PivotalCombo c = new PivotalCombo();
    c.hasPivot = true;
    NODE_TYPE u, v;

    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        EDGE_TYPE rank = getMaxTrussness(mlg, index[i], query);
        EDGE_TYPE ki = ci->order[i][rank].k;
        c.lamb = ci->order[i][rank].lamb;
        c.pivotLayer = ci->pivotalList[i].pivotLayer;
        // LOG(INFO) << fmt::format("Evaluating query = {}, combo = {}", query, ci->idToComboName[i]);
        timer->startTimer();
        queryNodeCommunityEdge(mlg, index[i], query, rank, qus);
        //        queryNodeCommunity(mlg, index[i], query, 1, qu);
        timer->endTimer();
        for (const auto &qu: qus) {
            vector<vector<EDGE_TYPE> > ccs;
            if (ki != 0) {
                vector<EDGE_TYPE> result;
                timer->startTimer();
                peelRk(mlg, qu, c, ki, result);
                connect(mlg, result, ccs, query);
                timer->endTimer();
            } else {
                ccs.emplace_back(qu);
            }
            for (auto cc: ccs) {
                set<NODE_TYPE> nodeInCommunity;
                nodeInCommunity.insert(query);
                for (auto edge: cc) {
                    mlg->schema->idToEdge(edge, &u, &v);
                    nodeInCommunity.insert(u);
                    nodeInCommunity.insert(v);
                }
                vector<NODE_TYPE> com;
                com.assign(nodeInCommunity.begin(), nodeInCommunity.end());
                double den = density(com, mlg, DEN_BETA);
                if (den > bestden) {
                    bestden = den;
                    bestqu.assign(qu.begin(), qu.end());
                    bestComboId = i;
                    best_k = rank;
                }
            }
        }
    }
    if (ci->isPair) {
        LOG(INFO) << fmt::format("Best community of vertex {} is pivot-lamb-k {}-{}-{} (rank={}) ({:.2f})", query,
                                 ci->pivotalList[bestComboId].pivotLayer, ci->order[bestComboId][best_k].lamb,
                                 ci->order[bestComboId][best_k].k, best_k, bestden);
    } else {
        LOG(INFO) << fmt::format("Best community of vertex {} is {}-{} ({:.2f})", query,
                                 ci->pivotalList[bestComboId].toString(),
                                 best_k, bestden);
    }
    // bestdiam = double(diameterMultilayer(bestqu, mlg));
    //    LOG(INFO) << fmt::format("Diameter = {:.4f}/{}", bestdiam, bestqu.size());
    bestdiam = 1;
}


void testPiCommunityByDensityUsingGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                               const vector<NODE_TYPE> &queries) {
    auto path = filesystem::path(input_dir_path);
    auto graph_path = path.append(fmt::format("{}.txt", dataset_name));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(graph_path);
    LAYER_TYPE maxComboSize = ALLOW_NO_PIVOT
                                  ? mlg->layersNum * (mlg->layersNum - 1) + mlg->layersNum
                                  : mlg->layersNum * (mlg->layersNum - 1);
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new PivotalComboIndex(mlg->layersNum);

    PivotTrussDecomposition(mlg, trussness, ci, ALLOW_NO_PIVOT);
    auto **index = new EqualTree *[ci->comboSize + 1];
    PivotIndexNoZero(mlg, index, ci, trussness);
    double sumDen = 0.0;
    double sumDiam = 0.0;
    NODE_TYPE havingDiam = 0;
    double den, diam;
    auto timer = new Timer();
    LOG(INFO) << fmt::format("Finding Communities for q in Queries = [{}]", fmt::join(queries, ", "));
    for (auto q: queries) {
        LOG(INFO) << fmt::format("Finding community for query vertex {}...", q);
        evalCommunityByDensity(mlg, index, ci, q, den, diam, timer);
        sumDen += den;
        sumDiam += diam;
        if (diam != 0) {
            havingDiam++;
        }
    }
    double avgDen = sumDen / double(queries.size());
    double avgDiam = sumDiam / double(havingDiam);
    double avgQueryTime = timer->getTimerSecond() / double(queries.size());
    LOG(INFO) << fmt::format("Avg Density of {} is {:.4f}, Diameter is {}, cost time = {:.4f}s.", dataset_name, avgDen,
                             avgDiam, avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        delete index[i];
    }
    delete ci;
}

void testRankCommunityByDensityUsingGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                                 const vector<NODE_TYPE> &queries) {
    auto path = filesystem::path(input_dir_path);
    auto graph_path = path.append(fmt::format("{}.txt", dataset_name));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(graph_path);
    LAYER_TYPE maxComboSize = ALLOW_NO_PIVOT
                                  ? mlg->layersNum * (mlg->layersNum - 1) + mlg->layersNum
                                  : mlg->layersNum * (mlg->layersNum - 1);
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *oci = new PivotalComboIndex(mlg->layersNum);

    PivotTrussDecomposition(mlg, trussness, oci);
    auto **rank = new EDGE_TYPE *[mlg->layersNum];
    auto *ci = new PivotalComboIndex(mlg->layersNum);
    generateRank(mlg, trussness, oci, rank, ci);
    auto **index = new EqualTree *[ci->comboSize + 1];
    PivotIndexNoZero(mlg, index, ci, rank);
    double sumDen = 0.0;
    double sumDiam = 0.0;
    NODE_TYPE havingDiam = 0;
    double den, diam;
    auto timer = new Timer();
    LOG(INFO) << fmt::format("Finding Communities for q in Queries = [{}]", fmt::join(queries, ", "));
    for (auto q: queries) {
        LOG(INFO) << fmt::format("Finding community for query vertex {}...", q);
        evalCommunityByUtilityAndDenAcc(mlg, index, ci, q, den, diam, timer);
        sumDen += den;
        sumDiam += diam;
        if (diam != 0) {
            havingDiam++;
        }
    }
    double avgDen = sumDen / double(queries.size());
    double avgDiam = sumDiam / double(havingDiam);
    double avgQueryTime = timer->getTimerSecond() / double(queries.size());
    LOG(INFO) << fmt::format("Avg Density of {} is {:.4f}, Diameter is {}, cost time = {:.4f}s.", dataset_name, avgDen,
                             avgDiam, avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        delete index[i];
    }
    delete ci;
}


void testPiRMCommunity(const filesystem::path &input_dir_path) {
    auto path = filesystem::path(input_dir_path);
    auto RM_path = path.append(fmt::format("RM.txt"));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(RM_path);
    LAYER_TYPE maxComboSize = ALLOW_NO_PIVOT
                                  ? mlg->layersNum * (mlg->layersNum - 1) + mlg->layersNum
                                  : mlg->layersNum * (mlg->layersNum - 1);
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new PivotalComboIndex(mlg->layersNum);

    PivotTrussDecomposition(mlg, trussness, ci, ALLOW_NO_PIVOT);
    auto **index = new EqualTree *[ci->comboSize + 1];
    PivotIndexAll(mlg, index, ci, trussness);
    double sumF1 = 0.0;
    double sumDen = 0.0;
    double den;
    auto rmgt1 = new vector<NODE_TYPE>{
        12, 36, 38, 39, 40, 41, 43, 45, 46, 47, 48, 51, 54, 58, 60, 62, 64, 65, 68,
        69, 72, 73, 76, 80, 82, 84, 85
    };
    auto rmgt2 = new vector<NODE_TYPE>{
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
        23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 37, 42, 44, 49, 50, 52,
        53, 55, 56, 57, 59, 61, 63, 66, 67, 70, 71, 74, 75, 77, 78, 79, 81, 83, 86,
        87, 88, 89, 90, 91
    };
    auto timer = new Timer();
    timer->startTimer();
    LOG(INFO) << fmt::format("Finding community for each query vertex...");
    for (auto q: *rmgt1) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt1);
        sumDen += den;
    }
    for (auto q: *rmgt2) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt2);
        sumDen += den;
    }
    timer->endTimer();
    double avgF1 = sumF1 / 91.0;
    double avgDen = sumDen / 91.0;
    double avgQueryTime = timer->getTimerSecond() / 91.0;
    LOG(INFO) << fmt::format("Avg F1 of RM is {:.4f}, avg Den is {:.4f}, cost time = {:.4f}s.", avgF1, avgDen,
                             avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        delete index[i];
    }
    delete ci;
}


void testPiWRMCommunity(const filesystem::path &input_dir_path) {
    auto path = filesystem::path(input_dir_path);
    auto RM_path = path.append(fmt::format("RM.txt"));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(RM_path);
    LAYER_TYPE maxComboSize = ALLOW_NO_PIVOT
                                  ? mlg->layersNum * (mlg->layersNum - 1) + mlg->layersNum
                                  : mlg->layersNum * (mlg->layersNum - 1);
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *oci = new PivotalComboIndex(mlg->layersNum);

    PivotTrussDecomposition(mlg, trussness, oci);
    auto **rank = new EDGE_TYPE *[mlg->layersNum];
    auto *ci = new PivotalComboIndex(mlg->layersNum);
    generateRank(mlg, trussness, oci, rank, ci);
    auto **index = new EqualTree *[ci->comboSize + 1];
    PivotIndexAll(mlg, index, ci, rank);
    double sumF1 = 0.0;
    double sumDen = 0.0;
    double den;
    auto rmgt1 = new vector<NODE_TYPE>{
        12, 36, 38, 39, 40, 41, 43, 45, 46, 47, 48, 51, 54, 58, 60, 62, 64, 65, 68,
        69, 72, 73, 76, 80, 82, 84, 85
    };
    auto rmgt2 = new vector<NODE_TYPE>{
        1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22,
        23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 37, 42, 44, 49, 50, 52,
        53, 55, 56, 57, 59, 61, 63, 66, 67, 70, 71, 74, 75, 77, 78, 79, 81, 83, 86,
        87, 88, 89, 90, 91
    };
    auto timer = new Timer();
    timer->startTimer();
    LOG(INFO) << fmt::format("Finding community for each query vertex...");
    for (auto q: *rmgt1) {
        sumF1 += evalGroundTruthCommunityRank(mlg, index, ci, q, den, *rmgt1);
        sumDen += den;
    }
    for (auto q: *rmgt2) {
        sumF1 += evalGroundTruthCommunityRank(mlg, index, ci, q, den, *rmgt2);
        sumDen += den;
    }
    timer->endTimer();
    double avgF1 = sumF1 / 91.0;
    double avgDen = sumDen / 91.0;
    double avgQueryTime = timer->getTimerSecond() / 91.0;
    LOG(INFO) << fmt::format("Avg F1 of RM is {:.4f}, avg Den is {:.4f}, cost time = {:.4f}s.", avgF1, avgDen,
                             avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        delete index[i];
    }
    delete ci;
}


void testPiTerroristCommunity(const filesystem::path &input_dir_path) {
    auto path = filesystem::path(input_dir_path);
    auto RM_path = path.append(fmt::format("terrorist.txt"));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(RM_path);
    LAYER_TYPE maxComboSize = ALLOW_NO_PIVOT
                                  ? mlg->layersNum * (mlg->layersNum - 1) + mlg->layersNum
                                  : mlg->layersNum * (mlg->layersNum - 1);
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new PivotalComboIndex(mlg->layersNum);

    PivotTrussDecomposition(mlg, trussness, ci, ALLOW_NO_PIVOT);
    auto **index = new EqualTree *[ci->comboSize + 1];
    PivotIndexAll(mlg, index, ci, trussness);
    double sumF1 = 0.0;
    double sumDen = 0.0;
    double den;
    auto rmgt1 = new vector<NODE_TYPE>{3, 18, 24, 36, 38, 39, 44, 48, 50, 55, 56, 57, 63, 64};
    auto rmgt2 = new vector<NODE_TYPE>{1, 2, 9, 12, 14, 15, 16, 19, 29, 30, 35, 62, 66};
    auto rmgt3 = new vector<NODE_TYPE>{4, 20, 21, 27, 42, 47, 53, 60, 65};
    auto rmgt4 = new vector<NODE_TYPE>{
        5, 6, 7, 8, 10, 11, 13, 17, 22, 23, 25, 26, 28, 31, 33, 34, 40, 46, 49, 51, 52,
        54, 58, 59, 61, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79
    };
    auto rmgt5 = new vector<NODE_TYPE>{32, 37, 41, 43, 45};
    auto timer = new Timer();
    timer->startTimer();
    LOG(INFO) << fmt::format("Finding community for each query vertex...");
    for (auto q: *rmgt1) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt1);
        sumDen += den;
    }
    for (auto q: *rmgt2) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt2);
        sumDen += den;
    }
    for (auto q: *rmgt3) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt3);
        sumDen += den;
    }
    for (auto q: *rmgt4) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt4);
        sumDen += den;
    }
    for (auto q: *rmgt5) {
        sumF1 += evalGroundTruthCommunity(mlg, index, ci, q, den, *rmgt5);
        sumDen += den;
    }
    timer->endTimer();
    double avgF1 = sumF1 / 79.0;
    double avgDen = sumDen / 79.0;
    double avgQueryTime = timer->getTimerSecond() / 79.0;
    LOG(INFO)
            << fmt::format("Avg F1 of terrorist is {:.4f}, avg density is {:.4f}, cost time = {:.4f}s.", avgF1, avgDen,
                           avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        delete index[i];
    }
    delete ci;
}


void testPiWTerroristCommunity(const filesystem::path &input_dir_path) {
    auto path = filesystem::path(input_dir_path);
    auto RM_path = path.append(fmt::format("terrorist.txt"));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(RM_path);
    LAYER_TYPE maxComboSize = ALLOW_NO_PIVOT
                                  ? mlg->layersNum * (mlg->layersNum - 1) + mlg->layersNum
                                  : mlg->layersNum * (mlg->layersNum - 1);
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *oci = new PivotalComboIndex(mlg->layersNum);

    PivotTrussDecomposition(mlg, trussness, oci);
    auto **rank = new EDGE_TYPE *[mlg->layersNum];
    auto *ci = new PivotalComboIndex(mlg->layersNum);
    generateRank(mlg, trussness, oci, rank, ci);
    auto **index = new EqualTree *[ci->comboSize + 1];
    PivotIndexAll(mlg, index, ci, rank);
    double sumF1 = 0.0;
    double sumDen = 0.0;
    double den;
    auto rmgt1 = new vector<NODE_TYPE>{3, 18, 24, 36, 38, 39, 44, 48, 50, 55, 56, 57, 63, 64};
    auto rmgt2 = new vector<NODE_TYPE>{1, 2, 9, 12, 14, 15, 16, 19, 29, 30, 35, 62, 66};
    auto rmgt3 = new vector<NODE_TYPE>{4, 20, 21, 27, 42, 47, 53, 60, 65};
    auto rmgt4 = new vector<NODE_TYPE>{
        5, 6, 7, 8, 10, 11, 13, 17, 22, 23, 25, 26, 28, 31, 33, 34, 40, 46, 49, 51, 52,
        54, 58, 59, 61, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79
    };
    auto rmgt5 = new vector<NODE_TYPE>{32, 37, 41, 43, 45};
    auto timer = new Timer();
    timer->startTimer();
    LOG(INFO) << fmt::format("Finding community for each query vertex...");
    for (auto q: *rmgt1) {
        sumF1 += evalGroundTruthCommunityRank(mlg, index, ci, q, den, *rmgt1);
        sumDen += den;
    }
    for (auto q: *rmgt2) {
        sumF1 += evalGroundTruthCommunityRank(mlg, index, ci, q, den, *rmgt2);
        sumDen += den;
    }
    for (auto q: *rmgt3) {
        sumF1 += evalGroundTruthCommunityRank(mlg, index, ci, q, den, *rmgt3);
        sumDen += den;
    }
    for (auto q: *rmgt4) {
        sumF1 += evalGroundTruthCommunityRank(mlg, index, ci, q, den, *rmgt4);
        sumDen += den;
    }
    for (auto q: *rmgt5) {
        sumF1 += evalGroundTruthCommunityRank(mlg, index, ci, q, den, *rmgt5);
        sumDen += den;
    }
    timer->endTimer();
    double avgF1 = sumF1 / 79.0;
    double avgDen = sumDen / 79.0;
    double avgQueryTime = timer->getTimerSecond() / 79.0;
    LOG(INFO)
            << fmt::format("Avg F1 of terrorist is {:.4f}, avg density is {:.4f}, cost time = {:.4f}s.", avgF1, avgDen,
                           avgQueryTime);

    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        delete index[i];
    }
    delete ci;
}


void getCommunityForGivenPiCombo(NODE_TYPE queryNode, MLGWithSchema *mlg, EDGE_TYPE *trussness,
                                 vector<vector<NODE_TYPE> > &ans) {
    NODE_TYPE u, v, w;
    EDGE_TYPE e, e1, e2;
    auto neighbors = new NODE_TYPE[mlg->sumGraph->maxDegree];
    NODE_TYPE neighborSize = 0;
    set<NODE_TYPE> trussNodes;
    vector<EDGE_TYPE> trussEdges;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    ans.clear();
    vector<NODE_TYPE> tempCommunity;
    // 先找到对应的最大的k
    EDGE_TYPE trussnessBound = 0;
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        trussnessBound = max(trussness[edge], trussnessBound);
    }
    if (trussnessBound == 0) return;


    // 然后找到所有满足k的社区
    EDGE_TYPE begin = 0;
    for (EDGE_TYPE k = begin; k <= trussnessBound; k++) {
        set<EDGE_TYPE> edgesAttachToQueryNode;
        for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
            set<EDGE_TYPE> currentEdgesInC;
            // 从一个连接查询点的边作为出发
            auto neighbor = mlg->sumGraph->edges[neighbor_i];
            auto edge = mlg->schema->edgeToId(queryNode, neighbor);
            if (trussness[edge] < k or edgesAttachToQueryNode.find(edge) != edgesAttachToQueryNode.end()) {
                continue;
            }
            edgesAttachToQueryNode.insert(edge);
            queue<EDGE_TYPE> toInsertEdges;
            toInsertEdges.push(edge);
            while (not toInsertEdges.empty()) {
                e = toInsertEdges.front();
                toInsertEdges.pop();
                if (currentEdgesInC.find(e) != currentEdgesInC.end()) {
                    continue;
                }
                edgesTriangleConnectedByGivenEdgeBasic(mlg, e, neighbors, neighborSize);
                mlg->schema->idToEdge(e, &u, &v);
                for (NODE_TYPE j = 0; j < neighborSize; j++) {
                    w = neighbors[j];
                    e1 = mlg->schema->edgeToId(u, w);
                    e2 = mlg->schema->edgeToId(v, w);
                    if (trussness[e1] >= k and trussness[e2] >= k) {
                        toInsertEdges.push(e1);
                        toInsertEdges.push(e2);
                    }
                    if (w == queryNode) {
                        edgesAttachToQueryNode.insert(e1);
                        edgesAttachToQueryNode.insert(e2);
                    }
                }
                currentEdgesInC.insert(e);
            }
            ans.resize(ans.size() + 1);
            // edges to nodes
            set<NODE_TYPE> nodeInCommunity;
            nodeInCommunity.insert(queryNode);
            for (auto ee: currentEdgesInC) {
                mlg->schema->idToEdge(ee, &u, &v);
                nodeInCommunity.insert(u);
                nodeInCommunity.insert(v);
            }
            ans[ans.size() - 1].assign(nodeInCommunity.begin(), nodeInCommunity.end());
            // double den = density(ans[ans.size() - 1], mlg, DEN_BETA);
            // LOG(INFO) << fmt::format("k={}, den={:.2f}, comm ={}", k, den, fmt::join(ans[ans.size() - 1].begin(), ans[ans.size() - 1].end(), ", "));
        }
    }
}


void getPiCommunity(NODE_TYPE queryNode, MLGWithSchema *mlg, EDGE_TYPE **trussness, PivotalComboIndex *ci,
                    double &bestden, Timer *timer) {
    bestden = 0;
    vector<NODE_TYPE> community;
    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        // 循环全部可能的combo，找到其中density最大的一个
        vector<vector<NODE_TYPE> > qus;
        timer->startTimer();
        getCommunityForGivenPiCombo(queryNode, mlg, trussness[i], qus);
        timer->endTimer();
        if (not qus.empty()) {
            LOG(INFO) << fmt::format("Found all trusses for {} in combo {}.", queryNode, ci->pivotalList[i].toString());
        }
        for (const auto &qu: qus) {
            double den = density(qu, mlg, DEN_BETA);
            //            LOG(INFO) << fmt::format("Combo = {}-{}, den = {}, community = [{}]", ci->idToComboName[i], trussnessBound,
            //                                     den, fmt::join(qu.begin(), qu.end(), ", "));
            if (den > bestden) {
                bestden = den;
            }
        }
    }
}

void searchPiCommByTrussnessForGivenQueries(const filesystem::path &input_dir_path, const string &dataset_name,
                                            const vector<NODE_TYPE> &queries) {
    auto path = filesystem::path(input_dir_path);
    auto graph_path = path.append(fmt::format("{}.txt", dataset_name));
    auto mlg = new MLGWithSchema();
    mlg->LoadFromFile(graph_path);
    auto maxComboSize = LAYER_TYPE(pow(2, mlg->layersNum)) + mlg->layersNum * LAYER_TYPE(pow(2, mlg->layersNum - 1));
    auto trussness = new EDGE_TYPE *[maxComboSize];
    auto *ci = new PivotalComboIndex(mlg->layersNum);
    PivotTrussDecomposition(mlg, trussness, ci);
    double sumDen = 0.0;
    double sumDiam = 0.0;
    double den, diam;
    auto timer = new Timer();
    LOG(INFO) << fmt::format("Finding Communities for q in Queries = [{}]", fmt::join(queries, ", "));
    for (auto q: queries) {
        getPiCommunity(q, mlg, trussness, ci, den, timer);
        // LOG(INFO) << fmt::format("Finding community for query vertex {}... den = {:.2f}", q, den);
        sumDen += den;
    }
    double avgDen = sumDen / double(queries.size());
    double avgDiam = sumDiam / double(queries.size());
    double avgQueryTime = timer->getTimerSecond() / double(queries.size());
    LOG(INFO) << fmt::format("Avg Density of {} is {:.4f}, Diameter is {}, cost time = {:.4f}s.", dataset_name, avgDen,
                             avgDiam, avgQueryTime);
    delete ci;
}
