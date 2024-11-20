//
// Created by Ryanw on 24-10-17.
//

#include "TriangleConnTruss.h"
#include <climits> // 包含 UINT_MAX


void edgesTriangleConnectedByGivenEdgeBasic(MLGWithSchema *mlg, EDGE_TYPE edge, NODE_TYPE *ret, NODE_TYPE &size) {
    // 最基础的连通
    NODE_TYPE u, v, w;
    size = 0;
    mlg->schema->idToEdge(edge, &u, &v);
    auto u_start = mlg->sumGraph->start(u);
    auto u_end = mlg->sumGraph->end(u);
    auto uList = &mlg->sumGraph->edges[u_start];
    auto uLen = u_end - u_start;
    auto v_start = mlg->sumGraph->start(v);
    auto v_end = mlg->sumGraph->end(v);
    auto vList = &mlg->sumGraph->edges[v_start];
    auto vLen = v_end - v_start;
    auto commonNeighbor = new NODE_TYPE[max(uLen, vLen)];
    NODE_TYPE commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
    for (NODE_TYPE neighbor_i = 0; neighbor_i < commonNeighborLen; neighbor_i++) {
        w = commonNeighbor[neighbor_i];
        ret[size] = w;
        size++;
    }
    delete[] commonNeighbor;
}


void edgesTriangleConnectedByGivenEdgeShare(MLGWithSchema *mlg, EDGE_TYPE edge, NODE_TYPE *ret, NODE_TYPE &size) {
    NODE_TYPE u, v, w;
    size = 0;
    set<NODE_TYPE> commonNeighborSet;
    mlg->schema->idToEdge(edge, &u, &v);
    for (LAYER_TYPE layer = 0; layer < mlg->layersNum; layer++) {
        auto u_start = mlg->graph_layers[layer].start(u);
        auto u_end = mlg->graph_layers[layer].end(u);
        auto uList = &mlg->graph_layers[layer].edges[u_start];
        auto uLen = u_end - u_start;
        auto v_start = mlg->graph_layers[layer].start(v);
        auto v_end = mlg->graph_layers[layer].end(v);
        auto vList = &mlg->graph_layers[layer].edges[v_start];
        auto vLen = v_end - v_start;
        auto commonNeighbor = new NODE_TYPE[max(uLen, vLen)];
        NODE_TYPE commonNeighborLen = intersect(uList, vList, uLen, vLen, commonNeighbor);
        for (NODE_TYPE neighbor_i = 0; neighbor_i < commonNeighborLen; neighbor_i++) {
            w = commonNeighbor[neighbor_i];
            commonNeighborSet.insert(w);
        }
        delete[] commonNeighbor;
    }
    for (auto node: commonNeighborSet) {
        ret[size] = node;
        size++;
    }
}


void TriangleConnectedTruss(MLGWithSchema *mlg, EDGE_TYPE *trussness, EqualTree *eqTree, Schema *usingSchema) {
    //    auto eqTree = new EqualTree();
    NODE_TYPE u, v, w;
    EDGE_TYPE edge, e1, e2, pvt, compactId, superNodeId, mt;
    EDGE_TYPE maxTrussness = 0;
    EDGE_TYPE schemaNum = usingSchema->edgesNum;
    auto edgeIdToCompactId = new EDGE_TYPE[schemaNum + 1];
    auto compactIdToEdgeId = new EDGE_TYPE[schemaNum + 1];
    for (EDGE_TYPE i = 0; i < schemaNum + 1; i++) {
        edgeIdToCompactId[i] = 0;
    }
    for (EDGE_TYPE i = 0; i < schemaNum + 1; i++) {
        compactIdToEdgeId[i] = 0;
    }
    set<EDGE_TYPE> parentPivots;
    auto neighbors = new NODE_TYPE[mlg->sumGraph->maxDegree];
    NODE_TYPE neighborSize = 0;

    // 先生成有序的边集
    for (edge = 0; edge < schemaNum; edge++) {
        maxTrussness = max(maxTrussness, trussness[edge]);
    }
    // auto trussnessToNum = *new vector<EDGE_TYPE>(maxTrussness + 1, 0);
    vector<EDGE_TYPE> trussnessToNum;
    trussnessToNum.assign(maxTrussness + 1, 0);
    for (edge = 0; edge < schemaNum; edge++) {
        trussnessToNum[trussness[edge]] += 1;
    }
    trussnessToNum[0] = 0;
    for (auto i = 2; i < trussnessToNum.size(); i++) {
        trussnessToNum[i] += trussnessToNum[i - 1];
    }
    for (auto i = trussnessToNum.size() - 1; i > 0; i--) {
        trussnessToNum[i] = trussnessToNum[i - 1] + 1;
    }
    for (edge = 0; edge < schemaNum; edge++) {
        if (trussness[edge] == 0) continue;
        edgeIdToCompactId[edge] = trussnessToNum[trussness[edge]];
        compactIdToEdgeId[trussnessToNum[trussness[edge]]] = edge;
        trussnessToNum[trussness[edge]] += 1;
    }
    trussnessToNum[0] = 1;

    auto uf = new UnionFind(trussnessToNum[maxTrussness] + 1);
    auto compactIdToSuperNodeId = new EDGE_TYPE[trussnessToNum[maxTrussness] + 1];
    for (EDGE_TYPE i = 0; i < trussnessToNum[maxTrussness] + 1; i++) {
        compactIdToSuperNodeId[i] = 0;
    }

    for (EDGE_TYPE k = maxTrussness; k > 0; k--) {
        EDGE_TYPE compactIdStart = trussnessToNum[k - 1];
        EDGE_TYPE compactIdEnd = trussnessToNum[k];
        //        LOG(INFO) << fmt::format("generating equal tree for k={}.", k);
        // 寻找parent
        //        LOG(INFO) << "find parent ...";
        for (compactId = compactIdStart; compactId < compactIdEnd; compactId++) {
            edge = compactIdToEdgeId[compactId];
            //            if (trussness[edge] != k) LOG(ERROR) << fmt::format("trussness={} in bucket={}.", trussness[edge], k);
            if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(basic)) {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            } else if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(shareEdge)) {
                edgesTriangleConnectedByGivenEdgeShare(mlg, edge, neighbors, neighborSize);
            } else {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            }
            usingSchema->idToEdge(edge, &u, &v);
            for (NODE_TYPE j = 0; j < neighborSize; j++) {
                w = neighbors[j];
                e1 = usingSchema->edgeToId(u, w);
                e2 = usingSchema->edgeToId(v, w);
                mt = min(trussness[e1], trussness[e2]);
                mt = min(mt, trussness[edge]);
                if (trussness[edge] == mt) {
                    pvt = uf->getPivot(edgeIdToCompactId[e1]);
                    parentPivots.insert(pvt);
                    pvt = uf->getPivot(edgeIdToCompactId[e2]);
                    parentPivots.insert(pvt);
                }
            }
        }

        // 连通
        //        LOG(INFO) << "union find ...";
        for (compactId = compactIdStart; compactId < compactIdEnd; compactId++) {
            edge = compactIdToEdgeId[compactId];
            if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(basic)) {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            } else if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(shareEdge)) {
                edgesTriangleConnectedByGivenEdgeShare(mlg, edge, neighbors, neighborSize);
            } else {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            }
            usingSchema->idToEdge(edge, &u, &v);
            for (NODE_TYPE j = 0; j < neighborSize; j++) {
                w = neighbors[j];
                e1 = usingSchema->edgeToId(u, w);
                e2 = usingSchema->edgeToId(v, w);
                if (trussness[e1] >= k and trussness[e2] >= k) {
                    uf->toUnion(compactId, edgeIdToCompactId[e1]);
                    uf->toUnion(compactId, edgeIdToCompactId[e2]);
                }
            }
        }

        // 建树
        //        LOG(INFO) << "construct tree...";
        for (compactId = compactIdStart; compactId < compactIdEnd; compactId++) {
            pvt = uf->getPivot(compactId);
            if (compactIdToSuperNodeId[pvt] == 0) {
                superNodeId = eqTree->superNodeNum;
                eqTree->addSuperNode(k);
                compactIdToSuperNodeId[pvt] = superNodeId;
            }
            superNodeId = compactIdToSuperNodeId[pvt];
            compactIdToSuperNodeId[compactId] = superNodeId;
            EDGE_TYPE edgeId = compactIdToEdgeId[compactId];
            eqTree->addEdgeToSuperNode(superNodeId, edgeId);
        }

        // 树的父子关系
        //        LOG(INFO) << "setup son-parent ...";
        for (auto e: parentPivots) {
            pvt = uf->getPivot(e);
            auto child = compactIdToSuperNodeId[e];
            auto parent = compactIdToSuperNodeId[pvt];
            if (child != parent) {
                eqTree->addParent(child, parent);
                eqTree->addChild(parent, child);
            }
        }

        parentPivots.clear();
    }

    delete uf;
    delete[] compactIdToEdgeId;
    delete[] edgeIdToCompactId;
    delete[] compactIdToSuperNodeId;
    delete[] neighbors;
    eqTree->compactTree();
    trussnessToNum.clear();
    parentPivots.clear();
}


void TriangleConnectedTrussAllowZero(MLGWithSchema *mlg, EDGE_TYPE *trussness, EqualTree *eqTree, Schema *usingSchema) {
    //    auto eqTree = new EqualTree();
    NODE_TYPE u, v, w;
    EDGE_TYPE edge, e1, e2, pvt, compactId, superNodeId, mt;
    EDGE_TYPE maxTrussness = 0;
    EDGE_TYPE schemaNum = usingSchema->edgesNum;
    auto edgeIdToCompactId = new EDGE_TYPE[schemaNum + 1];
    auto compactIdToEdgeId = new EDGE_TYPE[schemaNum + 1];
    for (EDGE_TYPE i = 0; i < schemaNum + 1; i++) {
        edgeIdToCompactId[i] = 0;
    }
    for (EDGE_TYPE i = 0; i < schemaNum + 1; i++) {
        compactIdToEdgeId[i] = 0;
    }
    set<EDGE_TYPE> parentPivots;
    auto neighbors = new NODE_TYPE[mlg->sumGraph->maxDegree];
    NODE_TYPE neighborSize = 0;

    // 先生成有序的边集
    for (edge = 0; edge < schemaNum; edge++) {
        maxTrussness = max(maxTrussness, trussness[edge]);
    }
    // auto trussnessToNum = *new vector<EDGE_TYPE>(maxTrussness + 1, 0);
    vector<EDGE_TYPE> trussnessToNum;
    trussnessToNum.assign(maxTrussness + 2, 0);
    for (edge = 0; edge < schemaNum; edge++) {
        trussnessToNum[trussness[edge]] += 1;
    }
    // trussnessToNum[0] = 0;
    for (auto i = 1; i < trussnessToNum.size(); i++) {
        trussnessToNum[i] += trussnessToNum[i - 1];
    }
    for (auto i = trussnessToNum.size() - 1; i > 0; i--) {
        trussnessToNum[i] = trussnessToNum[i - 1] + 1;
    }
    trussnessToNum[0] = 1;
    for (edge = 0; edge < schemaNum; edge++) {
        // if (trussness[edge] == 0) continue;
        edgeIdToCompactId[edge] = trussnessToNum[trussness[edge]];
        compactIdToEdgeId[trussnessToNum[trussness[edge]]] = edge;
        trussnessToNum[trussness[edge]] += 1;
    }
    // trussnessToNum[0] = 1;

    auto uf = new UnionFind(trussnessToNum[maxTrussness] + 1);
    auto compactIdToSuperNodeId = new EDGE_TYPE[trussnessToNum[maxTrussness] + 1];
    for (EDGE_TYPE i = 0; i < trussnessToNum[maxTrussness] + 1; i++) {
        compactIdToSuperNodeId[i] = 0;
    }

    for (EDGE_TYPE k = maxTrussness; k != UINT_MAX; k--) {
        EDGE_TYPE compactIdStart;
        if (k == 0) {
            compactIdStart = 0;
        } else {
            compactIdStart = trussnessToNum[k - 1];
        }
        EDGE_TYPE compactIdEnd = trussnessToNum[k];
        //        LOG(INFO) << fmt::format("generating equal tree for k={}.", k);
        // 寻找parent
        //        LOG(INFO) << "find parent ...";
        for (compactId = compactIdStart; compactId < compactIdEnd; compactId++) {
            edge = compactIdToEdgeId[compactId];
            //            if (trussness[edge] != k) LOG(ERROR) << fmt::format("trussness={} in bucket={}.", trussness[edge], k);
            if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(basic)) {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            } else if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(shareEdge)) {
                edgesTriangleConnectedByGivenEdgeShare(mlg, edge, neighbors, neighborSize);
            } else {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            }
            usingSchema->idToEdge(edge, &u, &v);
            for (NODE_TYPE j = 0; j < neighborSize; j++) {
                w = neighbors[j];
                e1 = usingSchema->edgeToId(u, w);
                e2 = usingSchema->edgeToId(v, w);
                mt = min(trussness[e1], trussness[e2]);
                mt = min(mt, trussness[edge]);
                if (trussness[edge] == mt) {
                    pvt = uf->getPivot(edgeIdToCompactId[e1]);
                    parentPivots.insert(pvt);
                    pvt = uf->getPivot(edgeIdToCompactId[e2]);
                    parentPivots.insert(pvt);
                }
            }
        }

        // 连通
        //        LOG(INFO) << "union find ...";
        for (compactId = compactIdStart; compactId < compactIdEnd; compactId++) {
            edge = compactIdToEdgeId[compactId];
            if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(basic)) {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            } else if (DEFAULT_TRIANGLE_CONNECT_LEVEL == triangleConnectLevelEnum(shareEdge)) {
                edgesTriangleConnectedByGivenEdgeShare(mlg, edge, neighbors, neighborSize);
            } else {
                edgesTriangleConnectedByGivenEdgeBasic(mlg, edge, neighbors, neighborSize);
            }
            usingSchema->idToEdge(edge, &u, &v);
            for (NODE_TYPE j = 0; j < neighborSize; j++) {
                w = neighbors[j];
                e1 = usingSchema->edgeToId(u, w);
                e2 = usingSchema->edgeToId(v, w);
                if (trussness[e1] >= k and trussness[e2] >= k) {
                    uf->toUnion(compactId, edgeIdToCompactId[e1]);
                    uf->toUnion(compactId, edgeIdToCompactId[e2]);
                }
            }
        }

        // 建树
        //        LOG(INFO) << "construct tree...";
        for (compactId = compactIdStart; compactId < compactIdEnd; compactId++) {
            pvt = uf->getPivot(compactId);
            if (compactIdToSuperNodeId[pvt] == 0) {
                superNodeId = eqTree->superNodeNum;
                eqTree->addSuperNode(k);
                compactIdToSuperNodeId[pvt] = superNodeId;
            }
            superNodeId = compactIdToSuperNodeId[pvt];
            compactIdToSuperNodeId[compactId] = superNodeId;
            EDGE_TYPE edgeId = compactIdToEdgeId[compactId];
            eqTree->addEdgeToSuperNode(superNodeId, edgeId);
        }

        // 树的父子关系
        //        LOG(INFO) << "setup son-parent ...";
        for (auto e: parentPivots) {
            pvt = uf->getPivot(e);
            auto child = compactIdToSuperNodeId[e];
            auto parent = compactIdToSuperNodeId[pvt];
            if (child != parent) {
                eqTree->addParent(child, parent);
                eqTree->addChild(parent, child);
            }
        }

        parentPivots.clear();
    }

    delete uf;
    delete[] compactIdToEdgeId;
    delete[] edgeIdToCompactId;
    delete[] compactIdToSuperNodeId;
    delete[] neighbors;
    eqTree->compactTree();
    trussnessToNum.clear();
    parentPivots.clear();
}


void queryNodeCommunity(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode, EDGE_TYPE trussnessBound,
                        vector<vector<NODE_TYPE> > &result) {
    result.clear();
    if (trussnessBound == 0) return;
    NODE_TYPE u, v;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    vector<EDGE_TYPE> edgesAttachedToQueryNode;
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        if (eqTree->getEdgeTrussness(edge) >= trussnessBound) {
            edgesAttachedToQueryNode.emplace_back(edge);
        }
    }
    auto superNodes = eqTree->edgesToSuperNodes(edgesAttachedToQueryNode, trussnessBound);
    result.resize(superNodes.size());
    EDGE_TYPE resultId = 0;
    for (auto superNode: superNodes) {
        // auto nodeInCommunity = new set<NODE_TYPE>();
        set<NODE_TYPE> nodeInCommunity;
        nodeInCommunity.insert(queryNode);
        auto temp = eqTree->getSubgraphBySuperNodeID(superNode);
        for (auto edge: temp) {
            mlg->schema->idToEdge(edge, &u, &v);
            nodeInCommunity.insert(u);
            nodeInCommunity.insert(v);
        }
        result[resultId].assign(nodeInCommunity.begin(), nodeInCommunity.end());
        resultId++;
    }
    superNodes.clear();
    superNodes.shrink_to_fit();
}


void queryNodeCommunityAllowZero(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode, EDGE_TYPE trussnessBound,
                                 vector<vector<NODE_TYPE> > &result) {
    result.clear();
    // if (trussnessBound == 0) return;
    NODE_TYPE u, v;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    vector<EDGE_TYPE> edgesAttachedToQueryNode;
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        if (eqTree->getEdgeTrussness(edge) >= trussnessBound) {
            edgesAttachedToQueryNode.emplace_back(edge);
        }
    }
    auto superNodes = eqTree->edgesToSuperNodes(edgesAttachedToQueryNode, trussnessBound);
    result.resize(superNodes.size());
    EDGE_TYPE resultId = 0;
    for (auto superNode: superNodes) {
        // auto nodeInCommunity = new set<NODE_TYPE>();
        set<NODE_TYPE> nodeInCommunity;
        nodeInCommunity.insert(queryNode);
        auto temp = eqTree->getSubgraphBySuperNodeID(superNode);
        for (auto edge: temp) {
            mlg->schema->idToEdge(edge, &u, &v);
            nodeInCommunity.insert(u);
            nodeInCommunity.insert(v);
        }
        result[resultId].assign(nodeInCommunity.begin(), nodeInCommunity.end());
        resultId++;
    }
    superNodes.clear();
    superNodes.shrink_to_fit();
}


void queryNodeCommunityEdge(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode, EDGE_TYPE trussnessBound,
                            vector<vector<EDGE_TYPE> > &result) {
    result.clear();
    // if (trussnessBound == 0) return;
    NODE_TYPE u, v;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    vector<EDGE_TYPE> edgesAttachedToQueryNode;
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        if (eqTree->getEdgeTrussness(edge) >= trussnessBound) {
            edgesAttachedToQueryNode.emplace_back(edge);
        }
    }
    auto superNodes = eqTree->edgesToSuperNodes(edgesAttachedToQueryNode, trussnessBound);
    result.resize(superNodes.size());
    EDGE_TYPE resultId = 0;
    for (auto superNode: superNodes) {
        set<EDGE_TYPE> edgeInCommunity;
        auto temp = eqTree->getSubgraphBySuperNodeID(superNode);
        for (auto edge: temp) {
            mlg->schema->idToEdge(edge, &u, &v);
            edgeInCommunity.insert(edge);
        }
        result[resultId].assign(edgeInCommunity.begin(), edgeInCommunity.end());
        resultId++;
    }
    superNodes.clear();
    superNodes.shrink_to_fit();
}


EDGE_TYPE getNodeMaxTrussness(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode) {
    NODE_TYPE u, v;
    EDGE_TYPE trussness = 0;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    vector<EDGE_TYPE> edgesAttachedToQueryNode;
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        edgesAttachedToQueryNode.emplace_back(edge);
    }
    for (auto edge: edgesAttachedToQueryNode) {
        auto t = eqTree->getEdgeTrussness(edge);
        trussness = max(trussness, t);
    }
    return trussness;
}


EDGE_TYPE getMaxTrussness(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode) {
    NODE_TYPE u, v;
    EDGE_TYPE maxTrussness = 0;
    auto queryNodeNeighborStart = mlg->sumGraph->start(queryNode);
    auto queryNodeNeighborEnd = mlg->sumGraph->end(queryNode);
    // auto edgesAttachedToQueryNode = new vector<EDGE_TYPE>();
    // auto result = new vector<NODE_TYPE>();
    // auto nodeInCommunity = new set<NODE_TYPE>();
    for (EDGE_TYPE neighbor_i = queryNodeNeighborStart; neighbor_i < queryNodeNeighborEnd; neighbor_i++) {
        auto neighbor = mlg->sumGraph->edges[neighbor_i];
        auto edge = mlg->schema->edgeToId(queryNode, neighbor);
        maxTrussness = max(eqTree->getEdgeTrussness(edge), maxTrussness);
    }
    return maxTrussness;
}
