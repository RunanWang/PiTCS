//
// Created by 王润安 on 2024/10/17.
//

#ifndef FOTRUSS_PIVOTALTRUSS_H
#define FOTRUSS_PIVOTALTRUSS_H

#include "../DataStructure/MLGWithSchema.h"
#include "../DataStructure/PivotalCombo.h"
#include "../DataStructure/EqualTree.h"
#include "TriangleConnTruss.h"

// Decomposition Algorithms
void PivotTrussDecomposition(MLGWithSchema *mlg, EDGE_TYPE **result, PivotalComboIndex *ci, bool allowNoPivot = false);

void generateRank(MLGWithSchema *mlg, EDGE_TYPE **trussness, PivotalComboIndex *ci, EDGE_TYPE **rank, PivotalComboIndex *nci);

// Index
void PivotIndexAll(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, EDGE_TYPE **trussness);

void PivotIndexNoZero(MLGWithSchema *mlg, EqualTree **index, PivotalComboIndex *ci, EDGE_TYPE **trussness);

void peelRk(MLGWithSchema *mlg, const vector<EDGE_TYPE> &edgesRemain, const PivotalCombo &combo, EDGE_TYPE rk,
            vector<EDGE_TYPE> &result);

#endif //FOTRUSS_PIVOTALTRUSS_H
