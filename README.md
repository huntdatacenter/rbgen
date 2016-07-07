BGEN reference implementation
========

This repository contains a reference implementation of the [BGEN format](http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format_v1.2.html) in C++, 
originally sourced from the [qctool](https://bitbucket.org/gavinband/bgen) implementation.  In addition two applications - `bgenix`, which indexes BGEN files, and `cat-bgen` which efficiently concatenates BGEN files - are provided.  See the [wiki](https://bitbucket.org/gavinband/bgen/wiki/Home) for documentation on these programs.

A utility program, `bgen_to_vcf`, is also provided; as the name suggests it converts a BGEN file to VCF.  Currently this is mainly intended as an example program that shows how to use the BGEN file reading API.

License
========
This BGEN implementation is released under the Boost Software License v1.0.  This is a relatively permissive open-source license that is compatible with many other open-source licenses.  See [this page](http://www.boost.org/users/license.html) and the file [LICENSE_1_0.txt](https://bitbucket.org/gavinband/bgen/src/tip/LICENSE_1_0.txt) for full details.

This repository also contains code from  the [sqlite](www.sqlite.org) and [boost](www.boost.org) libraries.  The former is [available in the public domain](http://www.sqlite.org/copyright.html) and the latter under the boost software license.  These libraries are not used in the core BGEN implementation, but may be used in the example programs provided.

Apps
=====

The following programs are currently built with the BGEN repository.

* [bgenix](https://bitbucket.org/gavinband/bgen/wiki/bgenix) - a tool to index and efficiently retrieve subsets of a BGEN file. 
* [cat-bgen](https://bitbucket.org/gavinband/bgen/wiki/cat-bgen) - a tool to efficiently concatenate BGEN files.

Click on the names to see the wiki pages.

Download
========

You can download a tarball of the latest master branch by visiting this link:

Use mercurial to download the master branch as follows:
```sh
hg clone https://gavinband@bitbucket.org/gavinband/bgen -u master
```


Compilation
=====

To compile the code, use the supplied waf build tool:
```sh
./waf-1.8.13 configure
./waf-1.8.13
```
Results will appear under the `build/` directory.  For more detail on compilation see [the wiki](https://bitbucket.org/gavinband/bgen/wiki/Troubleshooting_compilation).

Testing
=====

BGEN's tests can be run by typing 
```sh
./build/test/test_bgen
```
If all goes well a message like `All tests passed` should be printed.

The example program `bgen_to_vcf` reads a bgen file (v1.1 or v1.2) and outputs it as a VCF file to stdout.  You can try running it
by typing
```sh
./build/example/bgen_to_vcf example/example.8bits.bgen
```
which should output vcf-formatted data to stdout.  We've provided further example bgen files in the `example/` subdirectory.

Branches
========

This repo follows the branch naming practice in which `master` represents the most up-to-date code considered in a 'releasable' state.  If you are interested in using bgen code in your own project, we therefore recommend cloning the `master` branch.  Code development takes place in the `default` branch and/or in feature branches branched from the `default` branch.  The command given above downloads the master branch, which is what most people will want.

More information
=====
See the [changelog](https://bitbucket.org/gavinband/bgen/src/master/CHANGELOG.md),
the [source code](https://bitbucket.org/gavinband/bgen/src) or
the [Wiki](https://bitbucket.org/gavinband/bgen/wiki/Home).