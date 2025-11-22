// File: fm_debug.h
#pragma once

// ── Master switch: build for experiments (no probes) ──────────────────────────
// Turn on with: -DFM_EXPERIMENT_MODE=1
#ifndef FM_EXPERIMENT_MODE
#  define FM_EXPERIMENT_MODE 0
#endif

// Default probe state depends on master switch
#if FM_EXPERIMENT_MODE
#  define FM_PROBE_DEFAULT 0
#else
#  define FM_PROBE_DEFAULT 1
#endif

// ── Fine-grained toggles (can be overridden via compiler flags) ──────────────
#ifndef FM_PROBE_PIPELINE     // [PIPELINE], high-level round info
#  define FM_PROBE_PIPELINE   FM_PROBE_DEFAULT
#endif
#ifndef FM_PROBE_SWITCH       // [SWITCH-PROBE] degree lists, sign hash
#  define FM_PROBE_SWITCH     FM_PROBE_DEFAULT
#endif
#ifndef FM_PROBE_PHASE        // [PHASE-PROBE], [OMEGA'-SAMPLE], [OMEGA'-STATS]
#  define FM_PROBE_PHASE      FM_PROBE_DEFAULT
#endif
#ifndef FM_PROBE_TBB          // [TBB-FILTER], [TBB-STATS], [TBB-SELECT], [ANNEAL]
#  define FM_PROBE_TBB        FM_PROBE_DEFAULT
#endif
#ifndef FM_PROBE_SP           // [SP-GATE], [TRI_CYC-PROFILE], [SP-BATCH], [SP-PROBE]
#  define FM_PROBE_SP         FM_PROBE_DEFAULT
#endif
#ifndef FM_PROBE_CALLBACK     // [CALLBACK] tallies from CPLEX cut callback
#  define FM_PROBE_CALLBACK   FM_PROBE_DEFAULT
#endif
#ifndef FM_PROBE_GREEDY       // [GREEDY-B], [ALIGN] heuristic prints
#  define FM_PROBE_GREEDY     FM_PROBE_DEFAULT
#endif
#ifndef FM_PROBE_ROUND        // [ROUND-POLICY]
#  define FM_PROBE_ROUND      FM_PROBE_DEFAULT
#endif

// ── Helpers to compile-out blocks & declarations completely ───────────────────
#if FM_PROBE_PIPELINE
#  define DBG_PIPE(...) do { __VA_ARGS__ } while(0)
#else
#  define DBG_PIPE(...) do {} while(0)
#endif

#if FM_PROBE_SWITCH
#  define DBG_SWITCH(...) do { __VA_ARGS__ } while(0)
#else
#  define DBG_SWITCH(...) do {} while(0)
#endif

#if FM_PROBE_PHASE
#  define DBG_PHASE(...) do { __VA_ARGS__ } while(0)
#  define DBG_PHASE_DECL(...) __VA_ARGS__
#else
#  define DBG_PHASE(...) do {} while(0)
#  define DBG_PHASE_DECL(...)
#endif

#if FM_PROBE_TBB
#  define DBG_TBB(...) do { __VA_ARGS__ } while(0)
#  define DBG_TBB_DECL(...) __VA_ARGS__
#else
#  define DBG_TBB(...) do {} while(0)
#  define DBG_TBB_DECL(...)
#endif

#if FM_PROBE_SP
#  define DBG_SP(...) do { __VA_ARGS__ } while(0)
#  define DBG_SP_DECL(...) __VA_ARGS__
#else
#  define DBG_SP(...) do {} while(0)
#  define DBG_SP_DECL(...)
#endif

#if FM_PROBE_CALLBACK
#  define DBG_CB(...) do { __VA_ARGS__ } while(0)
#  define DBG_CB_DECL(...) __VA_ARGS__
#else
#  define DBG_CB(...) do {} while(0)
#  define DBG_CB_DECL(...)
#endif

#if FM_PROBE_GREEDY
#  define DBG_GREEDY(...) do { __VA_ARGS__ } while(0)
#  define DBG_GREEDY_DECL(...) __VA_ARGS__
#else
#  define DBG_GREEDY(...) do {} while(0)
#  define DBG_GREEDY_DECL(...)
#endif

#if FM_PROBE_ROUND
#  define DBG_ROUND(...) do { __VA_ARGS__ } while(0)
#else
#  define DBG_ROUND(...) do {} while(0)
#endif
