add_subdirectory(gtest EXCLUDE_FROM_ALL)

set(CMAKE_CXX_CLANG_TIDY "")

if(${CURVE} STREQUAL "BN128")
  include_directories(ate-pairing/include)
  include_directories(xbyak)
  add_library(
    zm
    STATIC

    ate-pairing/src/zm.cpp
    ate-pairing/src/zm2.cpp
  )
endif()
