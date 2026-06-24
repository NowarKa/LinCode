include(FetchContent)
cmake_minimum_required(VERSION 3.23...3.31.6)

project(LinCode LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

FetchContent_Declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG v2.4.2
)

FetchContent_MakeAvailable(CLI11)

option(ENABLE_ASAN "Enable AddressSanitizer" ON)

if(ENABLE_ASAN)
  message(STATUS "AddressSanitizer and UndefinedBehaviorSanitizer enabled")

  add_compile_options(
    -fsanitize=address,undefined
    -fno-omit-frame-pointer
    -g
  )

  add_link_options(
    -fsanitize=address,undefined
  )

  set(ENV{ASAN_OPTIONS} "detect_leaks=1:halt_on_error=1")
endif()

#################################
###### FiniteField Library ######
#################################

add_subdirectory(finite_field)

#################################
######## LinCode Library ########
#################################

add_library(LinCodeLibrary

  src/constraints_phase1.cc
  src/constructed_codes_table.cc
  src/linear_code.cc
  src/phase0.cc
)

target_include_directories(LinCodeLibrary
  PUBLIC
  ${CMAKE_SOURCE_DIR}/include
)

target_link_directories(LinCodeLibrary
  PRIVATE
  /usr/local/lib
)

target_link_libraries(LinCodeLibrary
  PUBLIC
  FiniteField
  scip
  CLI11::CLI11
)

#################################
####### LinCode Executable ######
#################################

add_executable(LinCode
  src/main.cc
)

target_link_libraries(LinCode
  PRIVATE
  LinCodeLibrary
)

#################################
######## Main project tests #####
#################################

# add_executable(...)
# target_link_libraries(...)

#################################
######## Utility targets ########
#################################

add_custom_target(output_clean
  COMMAND rm -f solutions problem_phase1.txt
)

add_custom_target(data_clean
  COMMAND rm -rf data/
)
