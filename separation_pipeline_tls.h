// separation_pipeline_tls.h

#pragma once
#include <vector>
#include "cycle_key.h"
#include "reheat_pool.h"
#include <unordered_set>

extern thread_local int g_sp_cycles_accepted;

// Reheat pool: stash non-violated cycle keys to retry in future fractional rounds
extern thread_local ReheatPool g_reheat_pool;

extern thread_local std::unordered_set<fmkey::CycleKey, fmkey::CycleKeyHash, fmkey::CycleKeyEq> g_reheat_inflight;
