
# file(GLOB SOURCE_V0 *.cpp)
# file(GLOB SOURCE_V0 graph.cpp)
file(GLOB SOURCE_V0 graph_opt.cpp)
# file(GLOB SOURCE_V0 graph_opt_seq.cpp)

project(graph LANGUAGES CXX)

add_executable(graph ${SOURCE_V0})

target_compile_options(graph PRIVATE /Qopt-report:max /QxCORE-AVX2 /ZI /DEBUG /Qipo- /Qopenmp /O3)

message(STATUS ${CMAKE_CURRENT_SOURCE_DIR})
