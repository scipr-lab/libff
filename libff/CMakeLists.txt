set(FF_EXTRASRCS)
set(FF_EXTRALIBS)
if(${CURVE} STREQUAL "BN128")
  set(
    FF_EXTRASRCS

    ${FF_EXTRASRCS}
    algebra/curves/bn128/bn128_fields.cpp
    algebra/curves/bn128/bn128_g1.cpp
    algebra/curves/bn128/bn128_g2.cpp
    algebra/curves/bn128/bn128_gt.cpp
    algebra/curves/bn128/bn128_init.cpp
    algebra/curves/bn128/bn128_pairing.cpp
    algebra/curves/bn128/bn128_pp.cpp
  )
  set(
    FF_EXTRALIBS

    ${FF_EXTRALIBS}
    zm
  )
endif()

add_library(
  ff
  STATIC

  algebra/curves/alt_bn128/alt_bn128_fields.cpp
  algebra/curves/alt_bn128/alt_bn128_g1.cpp
  algebra/curves/alt_bn128/alt_bn128_g2.cpp
  algebra/curves/alt_bn128/alt_bn128_init.cpp
  algebra/curves/alt_bn128/alt_bn128_pairing.cpp
  algebra/curves/alt_bn128/alt_bn128_pp.cpp
  algebra/curves/edwards/edwards_fields.cpp
  algebra/curves/edwards/edwards_g1.cpp
  algebra/curves/edwards/edwards_g2.cpp
  algebra/curves/edwards/edwards_init.cpp
  algebra/curves/edwards/edwards_pairing.cpp
  algebra/curves/edwards/edwards_pp.cpp
  algebra/curves/mnt/mnt4/mnt4_fields.cpp
  algebra/curves/mnt/mnt4/mnt4_g1.cpp
  algebra/curves/mnt/mnt4/mnt4_g2.cpp
  algebra/curves/mnt/mnt4/mnt4_init.cpp
  algebra/curves/mnt/mnt4/mnt4_pairing.cpp
  algebra/curves/mnt/mnt4/mnt4_pp.cpp
  algebra/curves/mnt/mnt46_common.cpp
  algebra/curves/mnt/mnt6/mnt6_fields.cpp
  algebra/curves/mnt/mnt6/mnt6_g1.cpp
  algebra/curves/mnt/mnt6/mnt6_g2.cpp
  algebra/curves/mnt/mnt6/mnt6_init.cpp
  algebra/curves/mnt/mnt6/mnt6_pairing.cpp
  algebra/curves/mnt/mnt6/mnt6_pp.cpp

  common/double.cpp
  common/profiling.cpp
  common/utils.cpp

  ${FF_EXTRASRCS}
)
target_link_libraries(
  ff

  GMP::gmp
  ${PROCPS_LIBRARIES}
  ${FF_EXTRALIBS}
)
target_include_directories(
  ff
  PUBLIC ..
)

# sodium for binary fields
find_path(
  SODIUM_INCLUDE_DIR
  sodium.h
  HINTS ${LIBSODIUM_INCLUDE_DIRS}
  /usr/local/include /opt/local/include /opt/include
)
find_library(
  SODIUM_LIBRARY
  NAMES sodium
  HINTS ${LIBSODIUM_LIBRARY_DIRS}
  /usr/local/lib /opt/local/lib /opt/lib
)

target_include_directories(
  ff
  PUBLIC
  ${SODIUM_INCLUDE_DIR}
)
target_compile_options(
  ff
  PUBLIC
  ${LIBSODIUM_CFLAGS} ${LIBSODIUM_CFLAGS_OTHER}
)
target_link_libraries(
  ff
  ${SODIUM_LIBRARY} ${LIBSODIUM_LDFLAGS_OTHER}
)

install(
  DIRECTORY "" DESTINATION "include/libff"
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.tcc"
  PATTERN "tests" EXCLUDE
  PATTERN "examples" EXCLUDE
)

install(
  TARGETS ff DESTINATION lib
)

# Tests
if ("${IS_LIBFF_PARENT}")
  # Tests for prime fields
  add_executable(
    algebra_bilinearity_test
    EXCLUDE_FROM_ALL

    algebra/curves/tests/test_bilinearity.cpp
  )
  target_link_libraries(
    algebra_bilinearity_test

    ff
    gtest_main
  )

  add_executable(
    algebra_groups_test
    EXCLUDE_FROM_ALL

    algebra/curves/tests/test_groups.cpp
  )
  target_link_libraries(
    algebra_groups_test

    ff
    gtest_main
  )

  add_executable(
    algebra_prime_fields_test
    EXCLUDE_FROM_ALL

    algebra/fields/tests/test_prime_fields.cpp
  )
  target_link_libraries(
    algebra_prime_fields_test

    ff
    gtest_main
  )

  add_executable(
    algebra_all_fields_test
    EXCLUDE_FROM_ALL

    algebra/fields/tests/test_all_fields.cpp
  )
  target_link_libraries(
    algebra_all_fields_test

    ff
    gtest_main
  )

  include(CTest)
  add_test(
    NAME algebra_bilinearity_test
    COMMAND algebra_bilinearity_test
  )
  add_test(
    NAME algebra_groups_test
    COMMAND algebra_groups_test
  )
  add_test(
    NAME algebra_prime_fields_test
    COMMAND algebra_prime_fields_test
  )
  add_test(
    NAME algebra_all_fields_test
    COMMAND algebra_all_fields_test
  )

  add_dependencies(check algebra_bilinearity_test)
  add_dependencies(check algebra_groups_test)
  add_dependencies(check algebra_prime_fields_test)
  add_dependencies(check algebra_all_fields_test)

  # Tests for binary fields
  add_executable(
    algebra_binary_fields_test
    EXCLUDE_FROM_ALL

    algebra/fields/tests/test_binary_fields.cpp
  )
  target_link_libraries(
    algebra_binary_fields_test

    ff
    gtest_main
  )

  add_test(
    NAME algebra_binary_fields_test
    COMMAND algebra_binary_fields_test
  )

  add_dependencies(check algebra_binary_fields_test)

  # Profiling for prime fields
  add_executable(
    multiexp_profile
    EXCLUDE_FROM_ALL

    algebra/scalar_multiplication/multiexp_profile.cpp
  )
  target_link_libraries(
    multiexp_profile

    ${OPENSSL_LIBRARIES}
    ff
  )

  add_dependencies(profile multiexp_profile)
endif()
