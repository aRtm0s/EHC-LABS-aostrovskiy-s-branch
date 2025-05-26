
file(GLOB SOURCE *.cpp)
add_executable(nbody-prj ${SOURCE})
message(STATUS ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(nbody-prj PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}\\..\\include)
target_compile_options(nbody-prj PRIVATE /ZI /DEBUG /Qopt-report:max)
