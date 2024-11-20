//
// Created by 王润安 on 2024/10/17.
//

#include "PivotalCombo.h"


PivotalCombo::PivotalCombo() {
    lamb = 0;
    pivotLayer = 0;
}

PivotalCombo::PivotalCombo(PivotalCombo *c) {
    hasPivot = c->hasPivot;
    lamb = c->lamb;
    pivotLayer = c->pivotLayer;
}

PivotalCombo::PivotalCombo(LAYER_TYPE l, bool has, LAYER_TYPE pl) {
    hasPivot = has;
    lamb = l;
    pivotLayer = pl;
}

string PivotalCombo::toString() const {
    if (hasPivot) return fmt::format("Pivot[{}]-Background{}", pivotLayer, lamb);
    else return fmt::format("Pivot[]-Background{}", lamb);
}

PivotalComboIndex::PivotalComboIndex(LAYER_TYPE numLayer) {
    eachLayerFirstComboPos.assign(numLayer, numLayer * (numLayer + 2));
    maxLambda.assign(numLayer + 1, 0);
}

LAYER_TYPE PivotalComboIndex::comboToInd(const PivotalCombo &combo) {
    // 检查是否出界，出界返回maxId+1
    if (combo.hasPivot) {
        if (combo.lamb > maxLambda[combo.pivotLayer]) return comboSize + 10;
        return eachLayerFirstComboPos[combo.pivotLayer] + combo.lamb - 1;
    } else {
        if (combo.lamb > maxLambda[maxLambda.size() - 1]) return comboSize + 10;
        return combo.lamb - 1;
    }
}

void PivotalComboIndex::addCombo(PivotalCombo *combo) {
    pivotalList.emplace_back(combo);
    if (combo->hasPivot) {
        LAYER_TYPE thisId = pivotalList.size() - 1;
        if (thisId < eachLayerFirstComboPos[combo->pivotLayer]) {
            eachLayerFirstComboPos[combo->pivotLayer] = thisId;
        }
        maxLambda[combo->pivotLayer] = max(maxLambda[combo->pivotLayer], combo->lamb);
    } else {
        maxLambda[maxLambda.size() - 1] = max(maxLambda[maxLambda.size() - 1], combo->lamb);
    }
    this->comboSize++;
}

LAYER_TYPE PivotalComboIndex::getMaxLambda(LAYER_TYPE pivotLayer, bool hasPivot) {
    if (hasPivot) {
        return maxLambda[maxLambda.size() - 1];
    } else {
        return maxLambda[pivotLayer];
    }
}

lkPair::lkPair(LAYER_TYPE il, EDGE_TYPE ik) {
    lamb = il;
    k = ik;
}

bool lkPair::operator==(const lkPair &other) const {
    return (lamb == other.lamb) && (k == other.k);
}

EDGE_TYPE lkPair::combine() const {
    return (lamb * 100) + k;
}

bool lkPair::operator<(const lkPair &other) const {
    return this->combine() < other.combine();
}

