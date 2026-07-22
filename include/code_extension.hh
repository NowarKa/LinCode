#pragma once

#include "linear_code.hh"
#include "utils.hh"

#include <memory>
#include <vector>

auto extend_code(LCode &code, Params &params, int thread_id) -> void;
