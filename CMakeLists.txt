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
#    frustration_model.cpp
#    frustration_model_xy.cpp
#    frustration_model_x.cpp
#    modularity_model.cpp
    triangle_bucket_batch.cpp
    negative_cycle_batch.cpp
    anneal_tri_budget.cpp
    reheat_pool.cpp
    signed_graph.cpp
    signed_graph_mip.cpp
    separation_pipeline.cpp
)

# Main executable
add_executable(frustration_index ${SOURCES})

target_link_libraries(frustration_index ilocplex cplex concert igraph m pthread)
