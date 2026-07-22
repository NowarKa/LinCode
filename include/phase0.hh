#pragma once

#include "utils.hh"
#include <memory>

auto check_feasibility(const Equations &system, int r,
                       shared_ptr<long long int> time_scip) -> bool;
