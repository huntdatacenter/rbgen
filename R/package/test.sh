#!/bin/bash
R CMD INSTALL build/R/rbgen
R --vanilla << HERE_DOC
library( rbgen )
D = bgen.load( "example/example.16bits.bgen", data.frame( chromosome = '01', start = 0, end = 100000 ))
str( D )
head( D\$variants )
D\$data[1,1:10,1:3]

HERE_DOC
