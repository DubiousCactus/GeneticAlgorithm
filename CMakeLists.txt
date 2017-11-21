cmake_minimum_required(VERSION 3.8)
project(GeneticAlgorithm)

set(CMAKE_C_STANDARD 99)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(SOURCE_FILES main.c)

add_executable(GeneticAlgorithm ${SOURCE_FILES})

target_link_libraries(GeneticAlgorithm m ncurses)
