// separation_pipeline_tls.h

#pragma once
#include <vector>
#include "cycle_key.h"
#include <unordered_set>

extern thread_local int g_sp_cycles_accepted;

extern thread_local std::unordered_set<fmkey::CycleKey, fmkey::CycleKeyHash, fmkey::CycleKeyEq> g_reheat_inflight;
