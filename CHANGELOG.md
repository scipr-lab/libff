## PENDING

_Special thanks to all downstream projects upstreaming their patches!_

### Breaking Changes
* None!

### Features
- #20 Improve operator+ speed for alt_bn, correct the corresponding docs, and reduce code duplication.
- #52 Change default build flags to work with Mac OS, and add instructions for finding macOS openSSL installation

### Bug fixes
- #54 Fix is_power_of_two for n > 2^32
- #19 Fix is_little_endian always returning true
- #20 Fix operator+ not contributing to alt_bn_128 profiling opcount 
- #26 Remove unused warnings in release build
- #39 Update Travis Config for newer Ubuntu mirror defaults