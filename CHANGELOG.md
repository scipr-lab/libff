## v1.1.0

_Special thanks to all downstream projects upstreaming their patches!_

### Breaking Changes
- None!

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