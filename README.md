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

We've supplied example files named `example.v11.bgen` (A BGEN v1.1-formatted file) and `example.<n>bits.bgen` (A BGEN v1.2-formatted file with <n> bits per probability)  in the `example/` subdirectory.  If all goes well the command

```sh
./build/bgen_to_vcf example/example.v11.bgen
```
and
```sh
./build/bgen_to_vcf example/example.<n>bits.bgen
```
should each run and output vcf-formatted data to stdout.

More information
=====
See the [changelog](https://bitbucket.org/gavinband/bgen/src/master/CHANGELOG.md), the [source code](https://bitbucket.org/gavinband/bgen/src) or the [Wiki](https://bitbucket.org/gavinband/bgen/wiki/Home).