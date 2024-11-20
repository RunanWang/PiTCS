//
// Created by Ryanw on 24-10-17.
//

#ifndef TRIANGLECONNTRUSS_H
#define TRIANGLECONNTRUSS_H

#include "../DataStructure/MLGWithSchema.h"
#include "../DataStructure/EqualTree.h"
#include "../constant.h"
#include "Algo.h"
#include "../DataStructure/UnionFind.h"

void edgesTriangleConnectedByGivenEdgeBasic(MLGWithSchema *mlg, EDGE_TYPE edge, NODE_TYPE *ret, NODE_TYPE &size);

void edgesTriangleConnectedByGivenEdgeShare(MLGWithSchema *mlg, EDGE_TYPE edge, NODE_TYPE *ret, NODE_TYPE &size);

void TriangleConnectedTruss(MLGWithSchema *mlg, EDGE_TYPE *trussness, EqualTree *eqTree, Schema *usingSchema);

void queryNodeCommunity(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode, EDGE_TYPE trussnessBound,
                        vector<vector<NODE_TYPE> > &result);

void TriangleConnectedTrussAllowZero(MLGWithSchema *mlg, EDGE_TYPE *trussness, EqualTree *eqTree, Schema *usingSchema);

void queryNodeCommunityAllowZero(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode, EDGE_TYPE trussnessBound,
                                 vector<vector<NODE_TYPE> > &result);

void queryNodeCommunityEdge(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode, EDGE_TYPE trussnessBound,
                            vector<vector<EDGE_TYPE> > &result);

EDGE_TYPE getNodeMaxTrussness(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode);

EDGE_TYPE getMaxTrussness(MLGWithSchema *mlg, EqualTree *eqTree, NODE_TYPE queryNode);

#endif //TRIANGLECONNTRUSS_H
