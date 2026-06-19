include(FetchContent)
cmake_minimum_required(VERSION 3.23...3.31.6)

project(LinCode LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
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

add_compile_options(
  -L/usr/local/lib
  -lscip
  )
#################################
###### GaloisField Library ######
#################################

#file(GLOB GALOIS_SOURCES CONFIGURE_DEPENDS
#  finite_field/*.cc
#)

#add_library(GaloisField STATIC ${GALOIS_SOURCES})

#target_include_directories(GaloisField PUBLIC
#  ${CMAKE_SOURCE_DIR}/finite_field
#)

#################################
######## LinCode Library ########
#################################

add_library(LinCodeLibrary
  src/field.cc
  src/modular_arith.cc
  src/polynomial_arith.cc
  src/field_element.cc
  src/constraints_phase1.cc  
  src/field_vector.cc
  src/hyperplane.cc
  src/linear_code.cc
  src/projective_space.cc
  src/projective_point.cc
  src/constructed_codes_table.cc
  src/phase0.cc
)

target_include_directories(LinCodeLibrary PUBLIC
  ${CMAKE_SOURCE_DIR}/include
)

target_link_directories(LinCodeLibrary PRIVATE
  /usr/local/lib
)

target_link_libraries(LinCodeLibrary PUBLIC
  scip
  CLI11::CLI11
)

#################################
####### LinCode Executable ######
#################################

add_executable(LinCode src/main.cc)

target_link_libraries(LinCode PRIVATE 
  LinCodeLibrary
  #GaloisField
)

#################################
####### Finte Field tests #######
#################################

add_executable(test_field
    tests/test_field.cc
)

target_link_libraries(test_field PRIVATE 
  LinCodeLibrary
  #GaloisField
)

target_compile_features(test_field PRIVATE cxx_std_17)

#################################
##### Projectif Space tests #####
#################################

add_executable(test_projectif_space
    tests/test_projectif_space.cc
)

target_link_libraries(test_projectif_space PRIVATE 
  LinCodeLibrary
)

target_compile_features(test_projectif_space PRIVATE cxx_std_17)
