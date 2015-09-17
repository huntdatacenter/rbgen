BGEN reference implementation
========

This repository contains a reference implementation of the [BGEN format](http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format_v1.2.html) in C++, 
originally sourced from the [qctool](https://bitbucket.org/gavinband/bgen) implementation.  A utility program, `bgen_to_vcf`, is also provided as an example using the API; as the name suggests it converts a BGEN file to VCF.

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

Testing
=====

If all goes well the command

```sh
./build/bgen_to_vcf example/example.bgen
./build/bgen_to_vcf example/example.bgen_v12
```
should run and output vcf-formatted data to stdout.