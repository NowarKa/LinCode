#pragma once

#include "linear_code.hh"
#include "utils.hh"

auto remove_equivalent_codes(const std::vector<LCode> &candidates,
                             Params &params) -> std::vector<LCode>;
