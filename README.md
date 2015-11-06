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

If all goes well the command

```sh
./build/bgen_to_vcf example/example.bgen
./build/bgen_to_vcf example/example.bgen_v12
```
should run and output vcf-formatted data to stdout.

More information
=====
See the [source code](https://bitbucket.org/gavinband/bgen/src) or the [Wiki](https://bitbucket.org/gavinband/bgen/wiki/Home).

History
====
6 Nov 2015
Major changes in revision 218f28a0cda6:

1. I’ve changed the behaviour of BGEN v1.2 with respect to samples with missing data: they are now stored with dummy zero probabilities.  The spec is now in 'beta' which means I don’t have any other planned changes to make; unless major issues are uncovered this will be the final version of the format.

2. I’ve revamped the setter api of parse_probability_data somewhat.  It is documented in the code and here [on the wiki](https://bitbucket.org/gavinband/bgen/wiki/The_Setter_API).  The main breaking changes are:
- Renamed operator() to set_value(), and given it an index argument; I think these make the API more consistent.
- Added an initial ploidy argument to set_number_of_entries() as requested.  (The type of data - phased or unphased - is already reported in the order_type argument so I don’t think another argument is needed).
- Added two new method calls, which are optional: set_min_max_ploidy() (useful for setting storage) and finalise().  See the docs for info.

3. I’ve also got rid of the max_id_size option to write_snp_identifying_data().  (This is now not needed because writing BGEN v1.0 files is no longer supported.)

4. I’ve also added some test code (using the [catch framework](https://github.com/philsquared/catch), which seems pretty good).  Tests are not exhaustive but hopefully a start.