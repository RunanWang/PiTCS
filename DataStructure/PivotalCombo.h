//
// Created by 王润安 on 2024/10/17.
//

#ifndef FOTRUSS_PIVOTALCOMBO_H
#define FOTRUSS_PIVOTALCOMBO_H

#include "../constant.h"
#include "string"
#include <fmt/format.h>
#include "vector"
#include <fmt/core.h>


using namespace std;

struct lkPair {
    LAYER_TYPE lamb;
    EDGE_TYPE k;

    lkPair(LAYER_TYPE il, EDGE_TYPE ik);

    EDGE_TYPE combine() const;

    bool operator==(const lkPair& other) const;
    bool operator<(const lkPair& other) const;
};

struct lkPairHash {
    std::size_t operator()(const lkPair& p) const {
        // 使用 std::hash 组合 a 和 b 的哈希值
        return hash<LAYER_TYPE>()(p.lamb) ^ (hash<EDGE_TYPE>()(p.k) << 1);
    }
};

template <>
struct fmt::formatter<lkPair> {
    constexpr auto parse(fmt::format_parse_context& ctx) { return ctx.begin(); }

    template <typename FormatContext>
    auto format(const lkPair& s, FormatContext& ctx) const {
        return fmt::format_to(ctx.out(), "{{{}-{}}}", s.lamb, s.k);
    }
};

class PivotalCombo {
public:
    LAYER_TYPE lamb{0};
    bool hasPivot{false};
    LAYER_TYPE pivotLayer{0};

    PivotalCombo();

    PivotalCombo(PivotalCombo *c);

    PivotalCombo(LAYER_TYPE l, bool has, LAYER_TYPE pl);

    string toString() const;
};

class PivotalComboIndex {
public:
    LAYER_TYPE comboSize{0};
    vector<PivotalCombo> pivotalList;
    vector<LAYER_TYPE> eachLayerFirstComboPos;
    vector<LAYER_TYPE> maxLambda;

    bool isPair = false;
    vector<vector<lkPair>> order;

    PivotalComboIndex(LAYER_TYPE numLayer);

    LAYER_TYPE comboToInd(const PivotalCombo &combo);

    void addCombo(PivotalCombo *combo);

    LAYER_TYPE getMaxLambda(LAYER_TYPE pivotLayer, bool hasPivot);
};

#endif //FOTRUSS_PIVOTALCOMBO_H
