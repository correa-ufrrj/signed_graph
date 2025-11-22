cmake_minimum_required(VERSION 3.10)

project(article_borneinf VERSION 0.1 DESCRIPTION "Frustration index computation with CPLEX and igraph")

# Headers from build dir (for config.h)
configure_file(config.h.in config.h)

# Linker search paths for CPLEX and CONCERT
link_directories(
    /opt/ibm/ILOG/CPLEX_Studio2211/cplex/lib/x86-64_linux/static_pic
    /opt/ibm/ILOG/CPLEX_Studio2211/concert/lib/x86-64_linux/static_pic
)

# Set include directories
include_directories(
    "${PROJECT_BINARY_DIR}"
    "${PROJECT_SOURCE_DIR}"
    /usr/include/igraph
    /usr/include/x86_64-linux-gnu/igraph
    /opt/ibm/ILOG/CPLEX_Studio2211/cplex/include
    /opt/ibm/ILOG/CPLEX_Studio2211/concert/include
)

# Gather source files
set(SOURCES
    main.cpp
    frustration_model.cpp
    frustration_model_xy.cpp
    frustration_model_x.cpp
    modularity_model.cpp
    triangle_bucket_batch.cpp
    negative_cycle_batch.cpp
    anneal_tri_budget.cpp
#    reheat_pool.cpp
    signed_graph.cpp
    signed_graph_mip.cpp
    separation_pipeline.cpp
)

# ---------------- Debug/Probe control (matches fm_debug.h) ----------------
#set_target_properties(frustration_index PROPERTIES
#  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")
#message(STATUS "FM_EXPERIMENT_MODE=${FM_EXPERIMENT_MODE} "
#               "→ PIPELINE=${FM_PROBE_PIPELINE} SWITCH=${FM_PROBE_SWITCH} "
#               "PHASE=${FM_PROBE_PHASE} TBB=${FM_PROBE_TBB} "
#               "SP=${FM_PROBE_SP} CALLBACK=${FM_PROBE_CALLBACK} "
#               "GREEDY=${FM_PROBE_GREEDY} ROUND=${FM_PROBE_ROUND}")

# Master switch: build for experiments (no probes)
option(FM_EXPERIMENT_MODE "Disable all probes for reproducible experiments" OFF)

# Per-probe switches (override master). If not set by the user, default to:
#   OFF when FM_EXPERIMENT_MODE=ON; otherwise ON.
macro(_default_probe VAR)
  if(NOT DEFINED ${VAR})
    if(FM_EXPERIMENT_MODE)
      set(${VAR} OFF CACHE BOOL "Enable ${VAR}" FORCE)
    else()
      set(${VAR} ON  CACHE BOOL "Enable ${VAR}" FORCE)
    endif()
  endif()
endmacro()

_default_probe(FM_PROBE_PIPELINE)   # [PIPELINE]
_default_probe(FM_PROBE_SWITCH)     # [SWITCH-PROBE]
_default_probe(FM_PROBE_PHASE)      # [PHASE-PROBE], [OMEGA'-*]
_default_probe(FM_PROBE_TBB)        # [TBB-*]
_default_probe(FM_PROBE_SP)         # [SP-*], [TRI_CYC-PROFILE]
_default_probe(FM_PROBE_CALLBACK)   # [CALLBACK]
_default_probe(FM_PROBE_GREEDY)     # [GREEDY-*], [ALIGN]
_default_probe(FM_PROBE_ROUND)      # [ROUND-POLICY]

# Convert ON/OFF to 1/0 for preprocessor
function(_bool_to_num in out)
  if(${in})
    set(${out} 1 PARENT_SCOPE)
  else()
    set(${out} 0 PARENT_SCOPE)
  endif()
endfunction()

_bool_to_num(FM_EXPERIMENT_MODE FM_EXPERIMENT_MODE_NUM)
_bool_to_num(FM_PROBE_PIPELINE   FM_PROBE_PIPELINE_NUM)
_bool_to_num(FM_PROBE_SWITCH     FM_PROBE_SWITCH_NUM)
_bool_to_num(FM_PROBE_PHASE      FM_PROBE_PHASE_NUM)
_bool_to_num(FM_PROBE_TBB        FM_PROBE_TBB_NUM)
_bool_to_num(FM_PROBE_SP         FM_PROBE_SP_NUM)
_bool_to_num(FM_PROBE_CALLBACK   FM_PROBE_CALLBACK_NUM)
_bool_to_num(FM_PROBE_GREEDY     FM_PROBE_GREEDY_NUM)
_bool_to_num(FM_PROBE_ROUND      FM_PROBE_ROUND_NUM)

# --- existing sources & target ---
add_executable(frustration_index ${SOURCES})

# Export the defines to the compiler for this target
target_compile_definitions(frustration_index PRIVATE
  FM_EXPERIMENT_MODE=${FM_EXPERIMENT_MODE_NUM}
  FM_PROBE_PIPELINE=${FM_PROBE_PIPELINE_NUM}
  FM_PROBE_SWITCH=${FM_PROBE_SWITCH_NUM}
  FM_PROBE_PHASE=${FM_PROBE_PHASE_NUM}
  FM_PROBE_TBB=${FM_PROBE_TBB_NUM}
  FM_PROBE_SP=${FM_PROBE_SP_NUM}
  FM_PROBE_CALLBACK=${FM_PROBE_CALLBACK_NUM}
  FM_PROBE_GREEDY=${FM_PROBE_GREEDY_NUM}
  FM_PROBE_ROUND=${FM_PROBE_ROUND_NUM}
)

# Optionally: different defaults by build type (uncomment if you like)
# if(CMAKE_BUILD_TYPE STREQUAL "Release")
#   target_compile_definitions(frustration_index PRIVATE FM_EXPERIMENT_MODE=1)
# endif()

# --- existing linking lines below ---
target_link_libraries(frustration_index ilocplex cplex concert igraph m pthread)
