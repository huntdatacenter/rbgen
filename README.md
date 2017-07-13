BGEN reference implementation
========

This repository contains a reference implementation of the [BGEN format](http://www.well.ox.ac.uk/~gav/bgen_format/bgen_format_v1.2.html) in C++, 
originally sourced from the [qctool](https://bitbucket.org/gavinband/bgen) implementation.  In addition three utilities - `bgenix`, which indexes BGEN files, `cat-bgen` which efficiently concatenates BGEN files, and `edit-bgen` which is used to manipulate BGEN metadata - are provided.
See the [wiki](https://bitbucket.org/gavinband/bgen/wiki/Home) for documentation on these programs.

An example program, `bgen_to_vcf`, is also provided; as the name suggests it converts a BGEN file to VCF.  This is intended as an example program that shows how to use the BGEN file reading API.

**!! Important note on the UK Biobank data**

The UK Biobank has released imputed data for the full release in BGEN format, with accompanying `bgenix` index files.  However, *these index files are not named in the way `bgenix` expects by default*.  Options for fixing this are:

1. Use `bgenix` to recreate the index files - e.g. by running `bgenix -g filename.bgen -index`.  This typically takes a few minutes per file, and is the recommended option because it includes additional metadata in the index file.)
2. Rename or copy each index files to the expected name, e.g. rename `ukb_imp_chr[N]_v2.bgi` to  `ukb_imp_chr[N]_v2.bgen.bgi`.
3. We will shortly add an option to `bgenix` to allow specifying the index files. 

See the [BGEN in the UK Biobank](https://bitbucket.org/gavinband/bgen/wiki/BGEN_in_the_UK_Biobank) page for further technical information on the UK Biobank data.

License
========
This BGEN implementation is released under the Boost Software License v1.0.  This is a relatively permissive open-source license that is compatible with many other open-source licenses.  See [this page](http://www.boost.org/users/license.html) and the file [LICENSE_1_0.txt](https://bitbucket.org/gavinband/bgen/src/tip/LICENSE_1_0.txt) for full details.

This repository also contains code from  the [sqlite](www.sqlite.org) and [boost](www.boost.org) libraries.  The former is [available in the public domain](http://www.sqlite.org/copyright.html) and the latter under the boost software license.  These libraries are not used in the core BGEN implementation, but may be used in the example programs provided.

Apps
=====

The following programs are currently built with the BGEN repository.

* [bgenix](https://bitbucket.org/gavinband/bgen/wiki/bgenix) - a tool to index and efficiently retrieve subsets of a BGEN file. 
* [cat-bgen](https://bitbucket.org/gavinband/bgen/wiki/cat-bgen) - a tool to efficiently concatenate BGEN files.
* [edit-bgen](https://bitbucket.org/gavinband/bgen/wiki/edit-bgen) - a tool to edit BGEN file metadata.

R package
========

An experimental R package called [rbgen](https://bitbucket.org/gavinband/bgen/wiki/rbgen) is also constructed in the build directory.  See the [rbgen wiki page](https://bitbucket.org/gavinband/bgen/wiki/rbgen) for more information on using this package.

Download
========

A tarball of the latest master branch is available here: http://bitbucket.org/gavinband/bgen/get/master.tar.gz.

Alternatively, use mercurial to download the master branch as follows:
```sh
hg clone https://gavinband@bitbucket.org/gavinband/bgen -u master
```
(This command can take a while.)

Additionally, pre-built version of the bgen utilities may be available from [this page](http://www.well.ox.ac.uk/~gav/resources/).  **Note**: the recommended use is to download and compile bgenix for your platform; these binaries are provided for convenience in getting started quickly.

Compilation
=====

To compile the code, use the supplied waf build tool:
```sh
./waf-1.8.13 configure
./waf-1.8.13
```
Results will appear under the `build/` directory.  

Note: a full build requires a compiler that supports C++-11, e.g. gcc v4.7 or above.  To specify the compiler used, set the `CXX` environment variable during the configure step.  For example (if your shell is `bash`):
```
CXX=/path/to/g++ ./waf-1.8.13 configure
./waf-1.8.13
```

We have tested compilation on gcc 4.9.3 and 5.4.0.

If you don't have access to a compiler with C++-11 support, you can still build the core bgen implementation, but won't be able to build the applications or example programs.  See [the wiki](https://bitbucket.org/gavinband/bgen/wiki/Troubleshooting_compilation) for more information.

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

Installation
========

The command
```sh
./waf-1.8.13 install
```
will install the applications listed above into a specified system or user directory.  By default this is `/usr/local`.  To change it, specify the prefix at the configure step:
```sh
./waf-1.8.13 configure --prefix=/path/to/installation/directory
./waf-1.8.13 install
```
The programs listed above will be installed into a folder called `bin/` under the prefix dir, e.g. `bgenix` will be installed as `/path/to/installation/directory/bin/bgenix` etc.

Note that in many cases there's no need for installation; the executables are self-contained.  The install step simply copies them into the destination directory.

(The installation prefix need not be a system-wide directory.  For example, I typically specify an installation directory within my home dir, e.g. `~gav/projects/software/`.

Branches
========

This repo follows the branch naming practice in which `master` represents the most up-to-date code considered in a 'releasable' state.  If you are interested in using bgen code in your own project, we therefore recommend cloning the `master` branch.  Code development takes place in the `default` branch and/or in feature branches branched from the `default` branch.  The command given above downloads the master branch, which is what most people will want.

More information
=====
See the [changelog](https://bitbucket.org/gavinband/bgen/src/master/CHANGELOG.md),
the [source code](https://bitbucket.org/gavinband/bgen/src) or
the [Wiki](https://bitbucket.org/gavinband/bgen/wiki/Home).