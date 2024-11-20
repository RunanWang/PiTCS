//
// Created by 王润安 on 2023/12/28.
//

#include "EqualTree.h"
#include "queue"
#include "set"
#include "algorithm"

EqualTree::EqualTree(EDGE_TYPE edgeSize) {
    edgeIdToSuperNodeId.assign(edgeSize, 0);
    // 下面把superId=0空出来
    edgesInSuperNodes.resize(1);
    superNodeChild.resize(1);
    superNodeParent.assign(1, 0);
    superNodeIdToTrussness.assign(1, 0);
};

void EqualTree::addSuperNode(EDGE_TYPE trussness) {
    edgesInSuperNodes.resize(edgesInSuperNodes.size() + 1);
    superNodeChild.resize(superNodeChild.size() + 1);
    superNodeParent.emplace_back(0);
    superNodeIdToTrussness.emplace_back(trussness);
    superNodeNum++;
}

void EqualTree::addEdgeToSuperNode(EDGE_TYPE superNodeId, EDGE_TYPE toAddEdgeId) {
    edgesInSuperNodes[superNodeId].emplace_back(toAddEdgeId);
    edgeIdToSuperNodeId[toAddEdgeId] = superNodeId;
}

void EqualTree::addChild(EDGE_TYPE thisSuperNodeId, EDGE_TYPE childSuperNodeId) {
    superNodeChild[thisSuperNodeId].emplace_back(childSuperNodeId);
}

void EqualTree::addParent(EDGE_TYPE thisSuperNodeId, EDGE_TYPE parentSuperNodeId) {
    superNodeParent[thisSuperNodeId] = parentSuperNodeId;
}

void EqualTree::compactTree() {
    superNodeIdToTrussness.shrink_to_fit();
    superNodeParent.shrink_to_fit();
    superNodeChild.shrink_to_fit();
    for (auto item: superNodeChild) {
        item.shrink_to_fit();
    }
    edgesInSuperNodes.shrink_to_fit();
    for (auto item: edgesInSuperNodes) {
        item.shrink_to_fit();
    }
}

EDGE_TYPE EqualTree::getSuperNodeID(EDGE_TYPE edgeId) {
    return edgeIdToSuperNodeId[edgeId];
}

EDGE_TYPE EqualTree::getEdgeTrussness(EDGE_TYPE edgeId) {
    EDGE_TYPE snid = edgeIdToSuperNodeId[edgeId];
    return superNodeIdToTrussness[snid];
}

vector<EDGE_TYPE> EqualTree::getSubgraphBySuperNodeID(EDGE_TYPE superNodeID) {
    EDGE_TYPE snid;
    auto ans = vector<EDGE_TYPE>();
    queue<EDGE_TYPE> q;
    q.emplace(superNodeID);
    while (not q.empty()) {
        snid = q.front();
        q.pop();
        for (auto child: superNodeChild[snid]) {
            q.emplace(child);
        }
        ans.insert(ans.end(), edgesInSuperNodes[snid].begin(), edgesInSuperNodes[snid].end());
    }
    return ans;
}

vector<EDGE_TYPE> EqualTree::edgesToSuperNodes(const vector<EDGE_TYPE> &edges) {
    // 返回superNode集合，这些集合中的superNode是frontier的，彼此不是对方的子树
    auto ans = vector<EDGE_TYPE>();
    auto snSet = set<EDGE_TYPE>();
    // 先把edgeID换成压缩后的superNodeId，并去重
    for (auto edgeId: edges) {
        snSet.insert(edgeIdToSuperNodeId[edgeId]);
    }
    vector<EDGE_TYPE> tempAns;
    for (auto snId: snSet) {
        tempAns.emplace_back(snId);
    }
    std::sort(tempAns.begin(), tempAns.end());
    // 接下来找这些superNode组成的树，对于每个边snIdi，我们检查有没有边snIdj，使j是i的父亲节点
    for (EDGE_TYPE i = 0; i < tempAns.size(); i++) {
        auto snIdi = tempAns[i];
        bool iIsRoot = true;
        for (EDGE_TYPE j = i + 1; j < tempAns.size(); j++) {
            auto snIdj = tempAns[j];
            if (superNodeIdToTrussness[snIdi] <= superNodeIdToTrussness[snIdj]) continue;
            auto checkSnId = snIdi;
            while (checkSnId != 0) {
                if (checkSnId == snIdj) {
                    iIsRoot = false;
                    break;
                }
                checkSnId = superNodeParent[checkSnId];
            }
            if (not iIsRoot) break;
        }
        if (iIsRoot) {
            ans.emplace_back(snIdi);
        }
    }
    return ans;
}


vector<EDGE_TYPE> EqualTree::edgesToSuperNodes(const vector<EDGE_TYPE> &edges, const EDGE_TYPE trussbound) {
    // 返回superNode集合，这些集合中的superNode是frontier的，彼此不是对方的子树
    auto ans = vector<EDGE_TYPE>();
    auto snSet = set<EDGE_TYPE>();
    // 先把edgeID换成压缩后的superNodeId，并去重
    for (auto edgeId: edges) {
        if(getEdgeTrussness(edgeId) != trussbound) {
            EDGE_TYPE snid = edgeIdToSuperNodeId[edgeId];
            while(superNodeIdToTrussness[snid] > trussbound) {
                snid = superNodeParent[snid];
            }
            snSet.insert(snid);
        }
        else {
            snSet.insert(edgeIdToSuperNodeId[edgeId]);
        }
    }
    vector<EDGE_TYPE> tempAns;
    for (auto snId: snSet) {
        tempAns.emplace_back(snId);
    }
    std::sort(tempAns.begin(), tempAns.end());
    // 接下来找这些superNode组成的树，对于每个边snIdi，我们检查有没有边snIdj，使j是i的父亲节点
    for (EDGE_TYPE i = 0; i < tempAns.size(); i++) {
        auto snIdi = tempAns[i];
        bool iIsRoot = true;
        for (EDGE_TYPE j = i + 1; j < tempAns.size(); j++) {
            auto snIdj = tempAns[j];
            if (superNodeIdToTrussness[snIdi] <= superNodeIdToTrussness[snIdj]) continue;
            auto checkSnId = snIdi;
            while (checkSnId != 0) {
                if (checkSnId == snIdj) {
                    iIsRoot = false;
                    break;
                }
                checkSnId = superNodeParent[checkSnId];
            }
            if (not iIsRoot) break;
        }
        if (iIsRoot) {
            ans.emplace_back(snIdi);
        }
    }
    return ans;
}


double EqualTree::getMemUsage() {
    double memInMB = 0;
    for (const auto &item: edgesInSuperNodes) {
        memInMB += double(item.size() * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;
    }
    memInMB += double(superNodeParent.size() * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;
    for (const auto &item: superNodeChild) {
        memInMB += double(item.size() * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;
    }
    memInMB += double(edgeIdToSuperNodeId.size() * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;
    memInMB += double(superNodeIdToTrussness.size() * sizeof(EDGE_TYPE)) / 1000.0 / 1000.0;
    return memInMB;
}


//FoTreeIndex::FoTreeIndex() {
//    idToCombo = *new vector<Combo *>();
//    comboNameToId = *new unordered_map<string, LAYER_TYPE>();
//}
//
//FoTreeIndex::~FoTreeIndex() {
//    idToCombo.clear();
//    idToCombo.shrink_to_fit();
//    comboNameToId.clear();
//    index.clear();
//}
//
//vector<LAYER_TYPE> FoTreeIndex::getOneHopFatherIds(Combo *childCombo) {
//    auto fatherCombo = new Combo(childCombo->layerNum);
//    string fatherName;
//    auto fathers = new vector<LAYER_TYPE>();
//    if (childCombo->lamb > 1 and childCombo->lamb != childCombo->focusNum) {
//        fatherCombo->copyCombo(childCombo);
//        fatherCombo->lamb -= 1;
//        fatherName = fatherCombo->toString();
//        if (comboNameToId.find(fatherName) != comboNameToId.end()) {
//            fathers->emplace_back(comboNameToId[fatherName]);
//        }
//    }
//    if (childCombo->focusNum != 0) {
//        for (LAYER_TYPE i = 0; i < childCombo->focusNum; i++) {
//            fatherCombo->copyCombo(childCombo);
//            fatherCombo->focusLayers[i] = fatherCombo->layerNum + 1;
//            fatherName = fatherCombo->toString();
////            LOG(INFO) << fmt::format("FatherName of {}: {}", idToCombo[childId].toString(), fatherName);
//            if (comboNameToId.find(fatherName) != comboNameToId.end()) {
//                fathers->emplace_back(comboNameToId[fatherName]);
//            }
//        }
//    }
//    delete fatherCombo;
//    return *fathers;
//}
//
//bool FoTreeIndex::allFatherInIndex(Combo *childCombo) {
//    auto fatherCombo = new Combo(childCombo->layerNum);
//    string fatherName;
//    if (childCombo->lamb > 1 and childCombo->lamb != childCombo->focusNum) {
//        fatherCombo->copyCombo(childCombo);
//        fatherCombo->lamb -= 1;
//        fatherName = fatherCombo->toString();
//        if (comboNameToId.find(fatherName) == comboNameToId.end()) {
//            return false;
//        }
//    }
//    if (childCombo->focusNum > 0) {
//        for (LAYER_TYPE i = 0; i < childCombo->focusNum; i++) {
//            fatherCombo->copyCombo(childCombo);
//            fatherCombo->focusLayers[i] = fatherCombo->layerNum + 1;
//            fatherName = fatherCombo->toString();
//            if (comboNameToId.find(fatherName) == comboNameToId.end()) {
//                return false;
//            }
//        }
//    }
//    delete fatherCombo;
//    return true;
//}
//
//void FoTreeIndex::insert(Combo *childCombo, EqualTree *eqTree) {
//    comboNameToId[childCombo->toString()] = idToCombo.size();
//    Combo ins = new Combo(childCombo);
//    idToCombo.emplace_back(&ins);
//    index.emplace_back(eqTree);
//}
//
//void FoTreeIndex::dump(const filesystem::path &path) {
//
//}
//
//void FoTreeIndex::load(const filesystem::path &path) {
//
//}