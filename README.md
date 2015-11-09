BGEN reference implementation
========

This repository contains a reference implementation of the [BGEN format](http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format_v1.2.html) in C++, 
originally sourced from the [qctool](https://bitbucket.org/gavinband/bgen) implementation.  A utility program, `bgen_to_vcf`, is also provided as an example; as the name suggests it converts a BGEN file to VCF.

Branches
========

This repo follows the standard branch naming practice in which `master` represents the most up-to-date code considered in a 'releasable' state.  If you are interested in using bgen code in your own project, we therefore recommend cloning the `master` branch.  Code development takes place in the `default` branch and/or in feature branches branched from the `default` branch.

Code can be cloned using mercurial like this:

`hg clone https://gavinband@bitbucket.org/gavinband/bgen`

Compilation
=====

To compile the code, either use make:
```sh
make
```

Or use the supplied waf build tool:
```sh
./waf-1.8.13 configure
./waf-1.8.13
```
Results will appear under the `build/` directory.

It's also possible to supply an installation prefix
```sh
./waf-1.8.13 configure --prefix /my/installation/dir
./waf-1.8.13 install
```

This will install the example and test applications into the specified directory.

Testing
=====

BGEN's tests can be run by typing 
```sh
./build/test/test_bgen
```
If all goes well a message like `All tests passed` should be printed.

The example program provided reads a bgen file (v1.1 or v1.2) and outputs it as a VCF file to stdout.  You can try running it
by typing
```sh
./build/example/bgen_to_vcf file.bgen
```
<<<<<<< local
which should output vcf-formatted data to stdout.  We've provided example bgen files in the `example/` subdirectory.

More information
=====
See the [changelog](https://bitbucket.org/gavinband/bgen/src/master/CHANGELOG.md),
the [source code](https://bitbucket.org/gavinband/bgen/src) or
the [Wiki](https://bitbucket.org/gavinband/bgen/wiki/Home).
