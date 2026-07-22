include(FetchContent)
cmake_minimum_required(VERSION 3.23...3.31.6)

project(LinCode LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#################################
######## Format and Tidy ########
#################################

find_program(CLANG_FORMAT_EXE NAMES clang-format)
find_program(CLANG_TIDY_EXE NAMES clang-tidy)

file(GLOB_RECURSE PROJECT_SOURCES
  src/*.cc
  src/*.hh
  include/*.hh
)

if(CLANG_FORMAT_EXE)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXE} -i ${PROJECT_SOURCES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running clang-format"
  )

  add_custom_target(format-check
    COMMAND ${CLANG_FORMAT_EXE}
      --dry-run
      --Werror
      ${PROJECT_SOURCES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Checking clang-format"
  )
else()
  message(WARNING "clang-format not found")
endif()

if(CLANG_TIDY_EXE)
  add_custom_target(tidy
    COMMAND ${CLANG_TIDY_EXE}
      -p=${CMAKE_BINARY_DIR}
      ${PROJECT_SOURCES}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    COMMENT "Running clang-tidy"
  )
else()
  message(WARNING "clang-tidy not found")
endif()

#################################
############## END ##############
#################################

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
  src/code_extension.cc
  src/code_equivalence.cc
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
######## Utility targets ########
#################################

add_custom_target(output_clean
  COMMAND rm -f solutions problem_phase1.txt
)

add_custom_target(data_clean
  COMMAND rm -rf data/
)
