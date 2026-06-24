cmake_minimum_required(VERSION 3.23)

project(FiniteField LANGUAGES CXX)

file(GLOB FINITE_FIELD_SOURCES CONFIGURE_DEPENDS
  src/*.cc
)

add_library(FiniteField STATIC
  ${FINITE_FIELD_SOURCES}
)

target_include_directories(FiniteField
  PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(FiniteField
  PUBLIC
  cxx_std_20
)

#################################
########### Tests ###############
#################################

add_executable(test_field
  tests/test_field.cc
)

target_link_libraries(test_field PRIVATE
  FiniteField
)

add_executable(test_projective_space
  tests/test_projectif_space.cc
)

target_link_libraries(test_projective_space PRIVATE
  FiniteField
)
