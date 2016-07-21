<h1 align="center">libff</h1>
<h4 align="center">C++ library for Finite Fields and Elliptic Curves</h4>

___libff___ is a C++ library for finite fields and elliptic curves. The library is developed by [SCIPR Lab] and contributors (see [AUTHORS] file) and is released under the MIT License (see [LICENSE] file).

## Table of contents

- [Directory Structure](#directory-structure)
- [Elliptic curve choices](#elliptic-curve-choices)
- [Build guide](#build-guide)

## Directory structure

The directory structure is as follows:

* [__src__](src): C++ source code, containing the following modules:
  * [__algebra__](src/algebra): fields and elliptic curve groups
  * [__common__](src/common): miscellaneous utilities
* [__third\_party__](third_party): third party libraries

## Elliptic curve choices

The libsnark library currently provides three options:

* "edwards":
   an instantiation based on an Edwards curve, providing 80 bits of security.

* "bn128":
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

* "alt_bn128":
   an alternative to "bn128", somewhat slower but avoids dynamic code generation.

Note that bn128 requires an x86-64 CPU while the other curve choices
should be architecture-independent.

## Build guide

The library has the following dependencies:

* [Boost](http://www.boost.org/)
* [CMake](http://cmake.org/)
* [GMP](http://gmplib.org/)
* libcrypto
* [libprocps](http://packages.ubuntu.com/trusty/libprocps-dev)

The library has been tested on Linux, but it is compatible with Windows and Mac OS X.

### Installation

On Ubuntu 14.04 LTS:

```
sudo apt-get install build-essential git libboost-all-dev cmake libgmp3-dev libssl-dev libprocps3-dev
```

Fetch dependencies from their GitHub repos:

    $ git submodule init && git submodule update

Create the Makefile:

    $ mkdir build && cd build && cmake .. 

Then, to compile the library, tests, and profiling harness, run this within the `build directory:

    $ make

### Using libff as a library

To build and install the libff library:

	$ DESTDIR=/install/path make install

This will install libff.a into /install/path/lib; so your application should be linked using -L/install/path/lib -lff. It also installs the requisite headers into /install/path/include; so your application should be compiled using -I/install/path/include.

[SCIPR Lab]: http://www.scipr-lab.org/ (Succinct Computational Integrity and Privacy Research Lab)

[LICENSE]: LICENSE (LICENSE file in top directory of libff distribution)

[AUTHORS]: AUTHORS (AUTHORS file in top directory of libff distribution)