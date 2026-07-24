#pragma once

#include <unordered_set>

#include "linear_code.hh"
#include "utils.hh"

auto remove_equivalent_codes(
    const unordered_set<pair<vector<int>, int>> &candidates,
    Params &params) -> vector<LCode>;
