## Pending

### Features
- #58 Add a defined API for every field type, and ensure all fields implement it (Thanks @alexander-zw)
- #71 Add BLS12-381 (Thanks @yelhousni)
- #79 Separate field initialization from curves (Thanks @alexander-zw)
- #80 Add clang-tidy checks to library and CI
- #82 Convert tests to use Google test (Thanks @alexander-zw)
- #85 Add more unit tests for fields (Thanks @alexander-zw)
- #86 Add binary fields from [libiop](https://github.com/scipr-lab/libiop) (Thanks @alexander-zw)

### Bug fixes
- #75 Get rid of warning for unused constant PI, in complex field
- #78 Reduce prints when inhibit_profiling_info is set
- #79 Use std::size_t for all curves, fix bugs introduced by #58

## v1.1.0

_Special thanks to all downstream projects upstreaming their patches!_

### Breaking Changes
- File structure changed: All field utils are now in `libff/algebra/field_utils/`, `Fp_model` is
  now in `libff/algebra/fields/prime_base/`, and all other F_p^n fields in
  `libff/algebra/fields/prime_extension/`.
- The function `base_field_char()` of all fields and curves have been renamed to `field_char()`.
- The provided fields used in curves have been moved to separate files so that they can be imported
  separately from `[field name]_fields.hpp`. However they are still accessible from the init file.

### Features
- #20 Improve operator+ speed for alt_bn, correct the corresponding docs, and reduce code duplication.
- #50 Add mul_by_cofactor to elliptic curve groups
- #50 Add sage scripts for altbn and mnt curves, to verify cofactors and generators
- #52 Change default procps build flags to work with Mac OS

### Bug fixes
- #19 Fix is_little_endian always returning true
- #20 Fix operator+ not contributing to alt_bn_128 profiling opcount 
- #26 Remove unused warnings in release build
- #39 Update Travis Config for newer Ubuntu mirror defaults
- #50 Fix incorrect mnt4 g2 generator
- #54 Fix is_power_of_two for n > 2^32
- #55 Throw informative error for division by zero in div_ceil
