<h1 align="center">libff</h1>
<h4 align="center">C++ library for Finite Fields and Elliptic Curves</h4>

___libff___ is a C++ library for finite fields and elliptic curves. The library is developed by [SCIPR Lab] and contributors (see [AUTHORS] file) and is released under the MIT License (see [LICENSE] file).

## Table of contents

- [Directory structure](#directory-structure)
- [Elliptic curve choices](#elliptic-curve-choices)
- [Build guide](#build-guide)

## Directory structure

The directory structure is as follows:

* [__libff__](libff): C++ source code, containing the following modules:
  * [__algebra__](libff/algebra): fields and elliptic curve groups
  * [__common__](libff/common): miscellaneous utilities
* [__depends__](depends): dependency libraries

## Elliptic curve choices

The libsnark library currently provides three options:

* `edwards`:
   an instantiation based on an Edwards curve, providing 80 bits of security.

* `bn128`:
   an instantiation based on a Barreto-Naehrig curve, providing 128
   bits of security. The underlying curve implementation is
   \[ate-pairing], which has incorporated our patch that changes the
   BN curve to one suitable for SNARK applications.

    *   This implementation uses dynamically-generated machine code for the curve
        arithmetic. Some modern systems disallow execution of code on the heap, and
        will thus block this implementation.

        For example, on Fedora 20 at its default settings, you will get the error
        `zmInit ERR:can't protect` when running this code. To solve this,
        run `sudo setsebool -P allow_execheap 1` to allow execution,
        or use `make CURVE=ALT_BN128` instead.

* `alt_bn128`:
   an alternative to `bn128`, somewhat slower but avoids dynamic code generation.

Note that `bn128` requires an x86-64 CPU while the other curve choices
should be architecture-independent.

## Build guide

The library has the following dependencies:

* [Boost](http://www.boost.org/)
* [CMake](http://cmake.org/)
* [GMP](http://gmplib.org/)
* [libsodium](https://libsodium.gitbook.io/doc/)
* [libprocps](http://packages.ubuntu.com/trusty/libprocps-dev) (turned off by default)

The library has been tested on Linux, but it is compatible with Windows and MacOS.

### Installation

On Ubuntu 14.04 LTS:

```
sudo apt-get install build-essential git libboost-all-dev cmake libgmp3-dev libssl-dev libprocps3-dev pkg-config libsodium-dev
```


On MacOS, all of the libraries from the previous section can be installed with `brew`, except for `libprocps`, which is turned off by default.

Fetch dependencies from their GitHub repos:

```
git submodule init && git submodule update
```

### Compilation

To compile, starting at the project root directory, create the build directory and Makefile:

```
mkdir build && cd build
cmake ..
```

If you are on macOS, change the cmake command to be

```
cmake .. -DOPENSSL_ROOT_DIR=$(brew --prefix openssl)
```

Other build flags include:
| Flag | Value | Description |
| ---- | ----- | ----------- |
| MAKE_INSTALL_PREFIX | (your path) | Specifies the desired install location. |
| CMAKE_BUILD_TYPE | Debug | Enables asserts. Note that tests now use gtest instead of asserts. |
| WITH_PROCPS | ON | Enables `libprocps`, which is by default turned off since it is not supported on some systems such as MacOS. |

Then, to compile and install the library, run this within the build directory:
```
make
make install
```

This will install `libff.a` into `/install/path/lib`; so your application should be linked using `-L/install/path/lib -lff`. It also installs the requisite headers into `/install/path/include`; so your application should be compiled using `-I/install/path/include`.

## Testing

To build and execute the tests for this library, run:
```
make check
```

## Code formatting and linting

To run clang-tidy on this library, specify the variable `USE_CLANG_TIDY` (eg. `cmake .. -D USE_CLANG_TIDY=ON`).
Then, run:
```
make clang-tidy
```

One can specify which clang-tidy checks to run and which files to run clang-tidy on using the `.clang-tidy` file in the root directory of the project.

## Profile

To compile the multi-exponentiation profiler in this library, run:
```
make profile
```
The resulting profiler is named `multiexp_profile` and can be found in the `libff` folder under the build directory.

[SCIPR Lab]: http://www.scipr-lab.org/ (Succinct Computational Integrity and Privacy Research Lab)

[LICENSE]: LICENSE (LICENSE file in top directory of libff distribution)

[AUTHORS]: AUTHORS (AUTHORS file in top directory of libff distribution)
