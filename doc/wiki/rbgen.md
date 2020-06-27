This page documents the `rbgen` R package.  For more on the bgen library in general, go [here](https://bitbucket.org/gavinband/bgen/wiki/Home).

# Introduction #

The `rbgen` package enables loading data from BGEN format files into [R](https://www.r-project.org/).  For this to work `rbgen` currently requires the BGEN file to have been indexed (by running [bgenix -index](http://bitbucket.org/gavinband/bgen/wiki/bgenix)).

**Note:** This package is currently considered experimental - if you use it, be prepared for crashes, breakages, loss of data, loss of time, loss of hair, etc.  (On the other hand, to the best of my knowledge it works fine.  Please report any issues to the [issue tracker](https://bitbucket.org/gavinband/bgen/issues).)

## How it works ##

`rbgen` uses the [Rcpp](http://www.rcpp.org/) package to glue R to the C++ code in the bgen repo. If you're interested in how this fits together, I wrote a [blog post about it](https://gavinband.github.io/biobank/bgen/2017/05/16/Getting_biobank_data_into_R.html).

# Obtaining and installing rbgen #

There are two ways to obtain the rbgen source package.  Firstly, you might find a downloadable version on [my website](http://www.well.ox.ac.uk/~gav/resources/).  A quick way to install it would be to do this in R:
```R
install.packages( "http://www.well.ox.ac.uk/~gav/resources/rbgen_v1.1.5.tgz", repos = NULL, type = "source" )
```

where you should adjust the version number to pick the most recent available version.  If all the planets have aligned, R will download, compile, and install the package, ending with a message like "`* DONE(rbgen)`".  You can now start R and load the library with `library(rbgen)`.

# Installing from a local copy #
Alternatively, assuming you have already downloaded the [bgen reference library](http://bitbucket.org/gavinband/bgen), you can choose to create the rbgen R package yourself, as follows:

* Configure and (optionally) build the bgen library, as described on the [overview page](https://bitbucket.org/gavinband/bgen/overview) and the [wiki page on troubleshooting compilation](http://bitbucket.org/gavinband/bgen/wiki/Troubleshooting_compilation).  Generally this involves running `./waf configure` followed by `./waf`.
* Run `./waf build_rbgen` to create the rbgen source package.

The last step creates an R package tarball in a temporary folder, and prints out its location.  You can then install it like this from the shell:

```sh
R CMD INSTALL /path/to/rbgen_<version>.tgz
```
or like this from within R:
```
install.packages( "/path/to/rbgen_<version>.tgz", repos = NULL )
```

# Package documentation #

The full docs are [here](http://www.well.ox.ac.uk/~gav/resources/bgen.load.html).

# A quick tutorial #

rbgen is a simple package containing a single function, `bgen.load()`. This function gets data for specific genomic ranges from an indexed bgen file.  Let's look at an example of use:
```R
library( rbgen )
ranges = data.frame(
  chromosome = "01",
  start = 1000,
  end = 2000
)
data = bgen.load( "example/example.16bits.bgen", ranges )
```
This command returns (hopefully) all the data you could want - look let's have a look at the result:
```R
> str(data)
List of 5
 $ variants:'data.frame':	2 obs. of  6 variables:
  ..$ chromosome       : Factor w/ 1 level "01": 1 1
  ..$ position         : int [1:2] 1001 2000
  ..$ rsid             : Factor w/ 2 levels "RSID_101","RSID_2": 1 2
  ..$ number_of_alleles: int [1:2] 2 2
  ..$ allele0          : Factor w/ 1 level "A": 1 1
  ..$ allele1          : Factor w/ 1 level "G": 1 1
 $ samples : chr [1:500] "sample_001" "sample_002" "sample_003" "sample_004" ...
 $ ploidy  : int [1:2, 1:500] 2 2 2 2 2 2 2 2 2 2 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:2] "RSID_101" "RSID_2"
  .. ..$ : chr [1:500] "sample_001" "sample_002" "sample_003" "sample_004" ...
 $ phased  : logi [1:2] FALSE FALSE
 $ data    : num [1:2, 1:500, 1:3] 0 NA 0.00784 0.02745 0.99608 ...
  ..- attr(*, "dimnames")=List of 3
  .. ..$ : chr [1:2] "RSID_101" "RSID_2"
  .. ..$ : chr [1:500] "sample_001" "sample_002" "sample_003" "sample_004" ...
  .. ..$ : chr [1:3] "g=0" "g=1" "g=2"
```
I.e. you get back the list of variants, the list of sample IDs, the ploidy for each sample at each variant, an indicator of whether data is phased or unphased at each variants, and the genotype probabilities.  The latter is stored as a 3d array with the first dimension being the variants, the second the samples, and the third the genotype. (Genotypes come in the order defined in the [bgen file spec](http://www.well.ox.ac.uk/~gav/bgen_format/)).

For example, here are the genotype probabilities for the first ten samples at the first variant:
```R
> data$data['RSID_101',1:10,]
                   g=0         g=1         g=2
sample_001 0.000000000 0.003921569 0.996078431
sample_002 0.007843137 0.000000000 0.992156863
sample_003 0.996078431 0.003921569 0.000000000
sample_004 0.000000000 0.996078431 0.003921569
sample_005 0.003921569 0.996078431 0.000000000
sample_006 0.003921569 0.000000000 0.996078431
sample_007 0.011764706 0.980392157 0.007843137
sample_008 0.007843137 0.992156863 0.000000000
sample_009 0.003921569 0.992156863 0.003921569
sample_010 0.003921569 0.976470588 0.019607843
```

## Handling multiallelic markers or non-diploid samples ##
In the example above all samples were diploid and all variants biallelic, so three probabilities per sample was enough to store the data.  To avoid wasting memory this is the default setting in `rbgen`.  In more complex examples you will need to explicitly ask for more storage:
```
data = bgen.load(
  "example/complex.bgen",
  ranges = data.frame( chromosome = '01', start = 0, end = 10000 ),
  max_entries = 36
)
```

In this file there are several multiallelic markers and some samples have non-diploid ploidy:
```R
> data$ploidy
    sample_0 sample_1 sample_2 sample_3
V1         1        2        2        2
V2         1        1        1        1
V3         1        2        2        2
M4         1        2        2        2
M5         1        3        3        2
M6         1        1        1        1
M7         1        1        1        1
M8         1        1        1        1
M9         1        1        1        2
M10        4        4        4        4
```

Unused genotype probability fields are filled with `NA`s, e.g.:
```R
> data$data['V1',,]
         g=0 g=1 g=2 g=3 g=4 g=5 g=6 g=7 g=8 g=9 g=10 g=11 g=12 g=13 g=14 g=15
sample_0   1   0  NA  NA  NA  NA  NA  NA  NA  NA   NA   NA   NA   NA   NA   NA
sample_1   1   0   0  NA  NA  NA  NA  NA  NA  NA   NA   NA   NA   NA   NA   NA
sample_2   1   0   0  NA  NA  NA  NA  NA  NA  NA   NA   NA   NA   NA   NA   NA
sample_3   0   1   0  NA  NA  NA  NA  NA  NA  NA   NA   NA   NA   NA   NA   NA
         g=16 g=17 g=18 g=19 g=20 g=21 g=22 g=23 g=24 g=25 g=26 g=27 g=28 g=29
sample_0   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA
sample_1   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA
sample_2   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA
sample_3   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA   NA
         g=30 g=31 g=32 g=33 g=34 g=35
sample_0   NA   NA   NA   NA   NA   NA
sample_1   NA   NA   NA   NA   NA   NA
sample_2   NA   NA   NA   NA   NA   NA
sample_3   NA   NA   NA   NA   NA   NA
```

Note that support for these more complex settings is still being worked on and is not yet feature complete.

