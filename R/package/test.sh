#!/bin/bash
R CMD INSTALL ./
R --vanilla << HERE_DOC
library( rbgen )
D = giveMeMyData( "example/example.16bits.bgen", chromosome = '01', start = 0, end = 100000 )
str( D )
head( D\$variants )
D\$data[1:2,1:10,1]
HERE_DOC
