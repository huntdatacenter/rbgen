#!/bin/bash
R="$1"
${R} CMD INSTALL build/R/rbgen
${R} --vanilla << HERE_DOC
library( rbgen )
D = bgen.load( "example/example.16bits.bgen", data.frame( chromosome = '01', start = 0, end = 100000 ))
str( D )
head( D\$variants )
D\$data[1,1:10,1:3]

D = bgen.load(  "example/complex.bgen", data.frame( chromosome = '01', start = 0, end = 1000000 ), max_entries_per_sample = 50 )
str(D)
HERE_DOC
