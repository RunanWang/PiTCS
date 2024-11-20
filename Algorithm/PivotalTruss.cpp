//
// Created by 王润安 on 2024/10/17.
//

#include "PivotalTruss.h"

#include <algorithm>

#include "Algo.h"
#include "../DataStructure/Bucket.h"
#include "set"
#include "../Utils/logging.h"
#include "../Utils/Timer.h"
#include "vector"
#include "queue"
#include "../constant.h"
#include "iostream"
#include "math.h"
#include <fmt/ranges.h>  // 支持打印容器

/*
 * Decomposition
 */
inline EDGE_TYPE calInd(EDGE_TYPE edge, LAYER_TYPE layer, LAYER_TYPE layerNum) {
    return edge * layerNum + layer;
}

inline double weightLambdaK(LAYER_TYPE lambda, EDGE_TYPE k, LAYER_TYPE layerNum, EDGE_TYPE maxK) {
    return WEIGHT_ALPHA * double(lambda) / double(layerNum) + double(k) / double(maxK);
}


struct lkPairCompare {
    double weight(LAYER_TYPE lambda, EDGE_TYPE k) const {
        return pow(double(lambda), WEIGHT_ALPHA) * double(k);
        // return pow(double(lambda), WEIGHT_ALPHA) * double(k + 0.1);
    }

    bool operator()(const lkPair &lhs, const lkPair &rhs) const {
        return weight(lhs.lamb, lhs.k) < weight(rhs.lamb, rhs.k);
    }
};


void generateLKPairOrderByValueFunction(const set<lkPair> &allPairs, vector<lkPair> &order) {
    set<lkPair, lkPairCompare> allPossiblePair;
    for (auto p: allPairs) {
        allPossiblePair.insert(p);
    }
    order.clear();
    order.assign(allPossiblePair.begin(), allPossiblePair.end());
}


EDGE_TYPE trussDecompositionForGivenCombo(MLGWithSchema *mlg, EDGE_TYPE *trussness, EDGE_TYPE *fatherRes,
                                          const PivotalCombo &combo) {
    // Init
    LAYER_TYPE lamb = combo.lamb;
    LAYER_TYPE layer = 0;
    NODE_TYPE u, v, w;
    EDGE_TYPE edge, edge1, edge2, newK;
    EDGE_TYPE maxK = 0;
    EDGE_TYPE maxTrussness = 0;

    LOG(INFO) << fmt::format("Decomposition for combo {}.", combo.toString());
    auto triangles = new EDGE_TYPE[mlg->schema->edgesNum * mlg->layersNum]; // 存每条边目前的三角形数目
    auto subGraph = new bool[mlg->schema->edgesNum]; // 存每条边目前是否在子图中
    auto toCheck = set<EDGE_TYPE>(); // 需要top-k检查的临时存储边set
    auto commonNeighbor = new NODE_TYPE[mlg->sumGraph->maxDegree * 2];
    bool *layerIsFocus = new bool[mlg->layersNum]; // 一个focus-layer的map,接下来对它初始化
    for (layer = 0; layer < mlg->layersNum; layer++) {
        layerIsFocus[layer] = false;
    }
    if (combo.hasPivot) {
        layerIsFocus[combo.pivotLayer] = true;
    }

    // triangles和trussness初始化
    mlg->genTriangle(triangles, fatherRes);
    for (edge = 0; edge < mlg->schema->edgesNum; edge++) {
        trussness[edge] = findTopK(&triangles[calInd(edge, 0, mlg->layersNum)], mlg->layersNum, lamb);
        trussness[edge] = min(trussness[edge], fatherRes[edge]);
        if (combo.hasPivot) {
            trussness[edge] = min(trussness[edge], triangles[calInd(edge, combo.pivotLayer, mlg->layersNum)]);
        }
        maxK = max(maxK, trussness[edge]);
        if (fatherRes[edge] == 0) { subGraph[edge] = false; } else { subGraph[edge] = true; }
    }
    auto bucket = new Bucket(mlg->schema->edgesNum, maxK);
    bucket->construct(trussness);

    for (NODE_TYPE k = 0; k < maxK + 1; k++) {
        auto toRemoveEdges = bucket->getBucket(k);
        if (bucket->length[k] != 0) {
            maxTrussness = k;
        }
        for (NODE_TYPE i = 0; i < bucket->length[k]; i++) {
            edge = toRemoveEdges[i];
            if (not subGraph[edge]) continue; // 这个edge被father已经移除了，它的trussness=0，并且triangle里没有计数，不用再peel
            subGraph[edge] = false;
            //            trussness[edge] = k;
            for (layer = 0; layer < mlg->layersNum; layer++) {
                if (not mlg->schema->edgeInSchema[edge][layer]) continue;
                mlg->schema->idToEdge(edge, &u, &v);
                auto u_start = mlg->graph_layers[layer].start(u);
                auto u_end = mlg->graph_layers[layer].end(u);
                auto uList = &mlg->graph_layers[layer].edges[u_start];
                auto uLen = u_end - u_start;
                auto v_start = mlg->graph_layers[layer].start(v);
                auto v_end = mlg->graph_layers[layer].end(v);
                auto vList = &mlg->graph_layers[layer].edges[v_start];
                auto vLen = v_end - v_start;
                NODE_TYPE commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
                for (NODE_TYPE neighbor_i = 0; neighbor_i < commonNeighborLen; neighbor_i++) {
                    w = commonNeighbor[neighbor_i];
                    edge1 = mlg->schema->edgeToId(u, w);
                    edge2 = mlg->schema->edgeToId(v, w);
                    if (subGraph[edge1] and subGraph[edge2]) {
                        // edge1
                        triangles[calInd(edge1, layer, mlg->layersNum)] -= 1;
                        newK = triangles[calInd(edge1, layer, mlg->layersNum)];
                        newK = min(newK, fatherRes[edge1]);
                        if (layerIsFocus[layer] and newK < trussness[edge1] and newK > k) {
                            bucket->changeItemBucket(edge1, trussness[edge1], newK);
                            trussness[edge1] = newK;
                        }
                        if (newK + 1 == trussness[edge1])
                            toCheck.insert(edge1);

                        // edge2
                        triangles[calInd(edge2, layer, mlg->layersNum)] -= 1;
                        newK = triangles[calInd(edge2, layer, mlg->layersNum)];
                        newK = min(newK, fatherRes[edge2]);
                        if (layerIsFocus[layer] and newK < trussness[edge2] and newK > k) {
                            bucket->changeItemBucket(edge2, trussness[edge2], newK);
                            trussness[edge2] = newK;
                        }
                        if (newK + 1 == trussness[edge2])
                            toCheck.insert(edge2);
                    }
                }
            }
            for (auto checkEdge: toCheck) {
                auto newTopK = findTopK(&triangles[calInd(checkEdge, 0, mlg->layersNum)], mlg->layersNum, lamb);
                newTopK = min(newTopK, fatherRes[checkEdge]);
                if (newTopK != trussness[checkEdge] and newTopK <= k) {
                    bucket->changeItemBucket(checkEdge, trussness[checkEdge], k);
                    trussness[checkEdge] = k;
                } else if (newTopK < trussness[checkEdge]) {
                    bucket->changeItemBucket(checkEdge, trussness[checkEdge], newTopK);
                    trussness[checkEdge] = newTopK;
                } else {
                    continue;
                }
            }
            toCheck.clear();
        }
    }
    LOG(INFO) << fmt::format("max trussness = {}", maxTrussness + 2);
    delete[] triangles;
    delete[] subGraph;
    delete[] layerIsFocus;
    delete bucket;
    delete[] commonNeighbor;
    return maxTrussness;
}

void PivotTrussDecomposition(MLGWithSchema *mlg, EDGE_TYPE **result, PivotalComboIndex *ci, bool allowNoPivot) {
    auto timer = new Timer();
    timer->startTimer();
    PivotalCombo *c;
    EDGE_TYPE maxTrussness = 0; // 用于判断分解的truss是否为空
    EDGE_TYPE schemaNum = mlg->schema->edgesNum; // Edge-schema的总数
    LAYER_TYPE comboSize = 0;

    // Init
    EDGE_TYPE *trussness; // 用于暂存分解结果
    auto fatherRes = new EDGE_TYPE[schemaNum]; // 父亲combo的trussness min值
    if (allowNoPivot) {
        for (LAYER_TYPE lamb = 1; lamb < mlg->layersNum; lamb++) {
            c = new PivotalCombo(lamb, false, 0);

            if (lamb != 1) {
                memcpy(fatherRes, result[lamb - 2], mlg->schema->edgesNum * sizeof(EDGE_TYPE));
            } else {
                for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                    fatherRes[edge] = mlg->maxNodesNum;
                }
            }
            trussness = new EDGE_TYPE[schemaNum];
            maxTrussness = trussDecompositionForGivenCombo(mlg, trussness, fatherRes, c);
            if (maxTrussness <= MAX_PRUNE_TRUSSNESS) {
                delete[] trussness;
                delete c;
                break;
            } else {
                result[comboSize] = trussness;
                comboSize++;
                ci->addCombo(c);
                // LOG(INFO) << fmt::format("Check: Id {} vs. Size {}.", ci->comboToInd(c), ci->comboSize);
            }
        }
    }


    for (LAYER_TYPE pLayer = 0; pLayer < mlg->layersNum; pLayer++) {
        for (LAYER_TYPE lamb = 1; lamb < mlg->layersNum; lamb++) {
            c = new PivotalCombo(lamb, true, pLayer);
            for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                fatherRes[edge] = mlg->maxNodesNum;
            }
            if (lamb != 1) {
                memcpy(fatherRes, result[ci->eachLayerFirstComboPos[c->pivotLayer] + lamb - 2],
                       schemaNum * sizeof(EDGE_TYPE));
            }
            if (allowNoPivot) {
                for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                    fatherRes[edge] = min(fatherRes[edge], result[lamb - 1][edge]);
                }
            }
            trussness = new EDGE_TYPE[schemaNum];
            maxTrussness = trussDecompositionForGivenCombo(mlg, trussness, fatherRes, c);
            if (maxTrussness <= MAX_PRUNE_TRUSSNESS) {
                delete[] trussness;
                delete c;
                break;
            } else {
                result[comboSize] = trussness;
                comboSize++;
                ci->addCombo(c);
                // LOG(INFO) << fmt::format("Check: Id {} vs. Size {}.", ci->comboToInd(c), ci->comboSize);
            }
        }
    }

    timer->endTimer();
    LOG(INFO) << fmt::format("Decomposition cost: {:.4f}s.", timer->getTimerSecond());

    auto result_size = double(comboSize * mlg->schema->edgesNum * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;

    LOG(INFO) << fmt::format("Non-empty meaningful Trussness size: {:.4f}MB.", result_size);
    LOG(INFO) << fmt::format("Calculated combo num: {}.", comboSize);
    LAYER_TYPE prunedSons = allowNoPivot
                                ? mlg->layersNum * (mlg->layersNum - 1) + mlg->layersNum
                                : mlg->layersNum * (mlg->layersNum - 1);
    prunedSons = prunedSons - ci->comboSize;
    LOG(INFO) << fmt::format("Pruned sons = {},  Useful Combos = {}.", prunedSons,
                             ci->comboSize);

    // statIntersect(result, ci, schemaNum);

    // De-allocate
    delete[] fatherRes;
    delete timer;
}


void generateRank(MLGWithSchema *mlg, EDGE_TYPE **trussness, PivotalComboIndex *ci, EDGE_TYPE **rank,
                  PivotalComboIndex *nci) {
    EDGE_TYPE schemaNum = mlg->schema->edgesNum;
    PivotalCombo combo;
    combo.hasPivot = true;
    LAYER_TYPE cind = 0;
    nci->isPair = true;

    for (LAYER_TYPE pLayer = 0; pLayer < mlg->layersNum; pLayer++) {
        combo.pivotLayer = pLayer;
        LAYER_TYPE maxLambdaOfThisPivotLayer = ci->getMaxLambda(pLayer, false);
        if (maxLambdaOfThisPivotLayer > 0) {
            // 先定序
            rank[pLayer] = new EDGE_TYPE[schemaNum];
            set<lkPair> allPairs;
            vector<lkPair> pairOrder;
            for (LAYER_TYPE lambda = 1; lambda <= maxLambdaOfThisPivotLayer; lambda++) {
                combo.lamb = lambda;
                cind = ci->comboToInd(combo);
                for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                    auto p = lkPair(lambda, trussness[cind][edge]);
                    allPairs.insert(p);
                }
            }
            // 下面的函数将allPairs中的全部pair返回一个顺序
            generateLKPairOrderByValueFunction(allPairs, pairOrder);
            LOG(INFO) << fmt::format("Order of pivot-layer {} as follow (lamb-k)...", pLayer);
            LOG(INFO) << fmt::format("{}", pairOrder);
            unordered_map<lkPair, EDGE_TYPE, lkPairHash> invertedIndex;
            for (EDGE_TYPE i = 0; i < pairOrder.size(); ++i) {
                invertedIndex[pairOrder[i]] = i;
            }

            // 然后把点洗入序中
            for (EDGE_TYPE edge = 0; edge < schemaNum; edge++) {
                vector<EDGE_TYPE> allPossibleRankOfThisEdge;
                for (LAYER_TYPE lambda = 1; lambda <= maxLambdaOfThisPivotLayer; lambda++) {
                    combo.lamb = lambda;
                    cind = ci->comboToInd(combo);
                    auto p = lkPair(lambda, trussness[cind][edge]);
                    allPossibleRankOfThisEdge.emplace_back(invertedIndex[p]);
                }
                // 这里可以改：min，或者其他洗入方法
                EDGE_TYPE maxRank = *max_element(allPossibleRankOfThisEdge.begin(), allPossibleRankOfThisEdge.end());
                // EDGE_TYPE maxRank = *min_element(allPossibleRankOfThisEdge.begin(), allPossibleRankOfThisEdge.end());
                rank[pLayer][edge] = maxRank;
            }
            auto c = new PivotalCombo(0, true, pLayer);
            nci->addCombo(c);
            nci->order.emplace_back(pairOrder);
        }
    }
}

/*
 * Index
 */
void PivotIndexAll(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, EDGE_TYPE **trussness) {
    LOG(INFO) << fmt::format("Build PivotTruss Index For all non-empty PivotTruss...");
    double totalIndexSize = 0;
    auto timer = new Timer();
    timer->startTimer();
    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        index[i] = new EqualTree(mlg->schema->edgesNum);
        TriangleConnectedTrussAllowZero(mlg, trussness[i], index[i], mlg->schema);
        delete[] trussness[i];
        totalIndexSize += index[i]->getMemUsage();
    }
    timer->endTimer();

    LOG(INFO) << fmt::format("Index Construction cost: {:.4f}s.", timer->getTimerSecond());
    LOG(INFO) << fmt::format("Index size: {:.4f}MB, index/graph = {:.2f}.", totalIndexSize,
                             (totalIndexSize / mlg->getMemUsage()));
    delete[] trussness;
    delete timer;
}


void PivotIndexNoZero(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, EDGE_TYPE **trussness) {
    LOG(INFO) << fmt::format("Build PivotTruss Index For all non-empty PivotTruss...");
    double totalIndexSize = 0;
    auto timer = new Timer();
    timer->startTimer();
    for (LAYER_TYPE i = 0; i < ci->comboSize; i++) {
        index[i] = new EqualTree(mlg->schema->edgesNum);
        TriangleConnectedTruss(mlg, trussness[i], index[i], mlg->schema);
        delete[] trussness[i];
        totalIndexSize += index[i]->getMemUsage();
    }
    timer->endTimer();

    LOG(INFO) << fmt::format("Index Construction cost: {:.4f}s.", timer->getTimerSecond());
    LOG(INFO) << fmt::format("Index size: {:.4f}MB, index/graph = {:.2f}.", totalIndexSize,
                             (totalIndexSize / mlg->getMemUsage()));
    delete[] trussness;
    delete timer;
}


void peelRk(MLGWithSchema *mlg, const vector<EDGE_TYPE> &edgesRemain, const PivotalCombo &combo, EDGE_TYPE rk,
            vector<EDGE_TYPE> &result) {
    // Init
    LAYER_TYPE lamb = combo.lamb;
    LAYER_TYPE layer = 0;
    NODE_TYPE u, v, w;
    EDGE_TYPE edge, edge1, edge2, newK;
    auto triangles = new EDGE_TYPE[mlg->schema->edgesNum * mlg->layersNum]; // 存每条边目前的三角形数目
    bool *subGraph = new bool[mlg->schema->edgesNum];
    for (edge = 0; edge < mlg->schema->edgesNum; edge++) {
        subGraph[edge] = false;
    }
    // bool subGraph[mlg->schema->edgesNum] = {};
    auto toCheck = set<EDGE_TYPE>(); // 需要top-k检查的临时存储边set
    auto commonNeighbor = new NODE_TYPE[mlg->sumGraph->maxDegree * 2];
    bool *layerIsFocus = new bool[mlg->layersNum]; // 一个focus-layer的map,接下来对它初始化
    for (layer = 0; layer < mlg->layersNum; layer++) {
        layerIsFocus[layer] = false;
    }
    layerIsFocus[combo.pivotLayer] = true;
    EDGE_TYPE* fatherRes = new EDGE_TYPE[mlg->schema->edgesNum];
    for (edge = 0; edge < mlg->schema->edgesNum; edge++) {
        fatherRes[edge] = 0;
    }
    vector<EDGE_TYPE> trussness(mlg->schema->edgesNum, 0);

    for (auto e: edgesRemain) {
        fatherRes[e] = rk;
    }
    mlg->genTriangle(triangles, fatherRes);

    queue<EDGE_TYPE> toPeel;
    for (auto e: edgesRemain) {
        trussness[e] = findTopK(&triangles[calInd(e, 0, mlg->layersNum)], mlg->layersNum, lamb);
        trussness[e] = min(trussness[e], triangles[calInd(e, combo.pivotLayer, mlg->layersNum)]);
        subGraph[e] = true;
        if (trussness[e] < rk) toPeel.emplace(e);
    }


    while (not toPeel.empty()) {
        edge = toPeel.front();
        toPeel.pop();
        if (not subGraph[edge]) continue; // 这个edge被father已经移除了，它的trussness=0，并且triangle里没有计数，不用再peel
        subGraph[edge] = false;
        for (layer = 0; layer < mlg->layersNum; layer++) {
            if (not mlg->schema->edgeInSchema[edge][layer]) continue;
            mlg->schema->idToEdge(edge, &u, &v);
            auto u_start = mlg->graph_layers[layer].start(u);
            auto u_end = mlg->graph_layers[layer].end(u);
            auto uList = &mlg->graph_layers[layer].edges[u_start];
            auto uLen = u_end - u_start;
            auto v_start = mlg->graph_layers[layer].start(v);
            auto v_end = mlg->graph_layers[layer].end(v);
            auto vList = &mlg->graph_layers[layer].edges[v_start];
            auto vLen = v_end - v_start;
            NODE_TYPE commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
            for (NODE_TYPE neighbor_i = 0; neighbor_i < commonNeighborLen; neighbor_i++) {
                w = commonNeighbor[neighbor_i];
                edge1 = mlg->schema->edgeToId(u, w);
                edge2 = mlg->schema->edgeToId(v, w);
                if (subGraph[edge1] and subGraph[edge2]) {
                    // edge1
                    triangles[calInd(edge1, layer, mlg->layersNum)] -= 1;
                    newK = triangles[calInd(edge1, layer, mlg->layersNum)];
                    if (layerIsFocus[layer] and newK < trussness[edge1]) {
                        trussness[edge1] = newK;
                        if (newK < rk) {
                            toPeel.push(edge1);
                        }
                    }
                    if (newK + 1 == trussness[edge1])
                        toCheck.insert(edge1);

                    // edge2
                    triangles[calInd(edge2, layer, mlg->layersNum)] -= 1;
                    newK = triangles[calInd(edge2, layer, mlg->layersNum)];
                    if (layerIsFocus[layer] and newK < trussness[edge2]) {
                        trussness[edge2] = newK;
                        if (newK < rk) {
                            toPeel.push(edge2);
                        }
                    }
                    if (newK + 1 == trussness[edge2])
                        toCheck.insert(edge2);
                }
            }
        }
        for (auto checkEdge: toCheck) {
            auto newTopK = findTopK(&triangles[calInd(checkEdge, 0, mlg->layersNum)], mlg->layersNum, lamb);
            if (newTopK != trussness[checkEdge] and newTopK < rk) {
                toPeel.push(checkEdge);
                trussness[checkEdge] = newTopK;
            } else if (newTopK < trussness[checkEdge]) {
                trussness[checkEdge] = newTopK;
            } else {
                continue;
            }
        }
        toCheck.clear();
    }
    for (auto e: edgesRemain) {
        if (trussness[e] >= rk) {
            result.emplace_back(e);
        }
    }
    delete[] triangles;
    delete[] layerIsFocus;
    delete[] commonNeighbor;
    delete[] fatherRes;
    delete[] subGraph;
}
