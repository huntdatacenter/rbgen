[bgenix](https://code.enkre.net/bgen.org/gavinband/bgen/src/default/apps/bgenix.cpp) is a tool to create an index of variants in a bgen file and to use that index for efficient retrieval of data for specific variants or regions.  

Cheat Sheet
====

Here's a quick list of common bgenix command lines and what they do:

| :Command line or option | :What it does |
| ----- | ------------ |
| :`bgenix -help` 											| :Print help on the various options bgenix supports |
| :`bgenix -g file.bgen -index` 							| :Don't output data; instead create an index file for the given bgen file.  It will be named `file.bgen.bgi`
| :`bgenix -g file.bgen` 									| :Output data from file.bgen (which must be already indexed) in the same format as file.bgen but ordered as in the index. |
| :`bgenix -g file.bgen -incl-rsids rs1` 					| :Restrict the output to variants with the specified rsid |
| :`bgenix -g file.bgen -incl-range 11:3500000-6500000` 	| :Restrict to variants in the given genomic range. |
| :`bgenix -g file.bgen -vcf` 								| :Transcode data to VCF format. |
| :`bgenix -g file.bgen -v11` 								| :Transcode data to BGEN v1.1 format. |
| :`bgenix -g file.bgen -list`								| :Don't output genotype data, just list the variants in the index. |

Query options and output format options can of course be combined, e.g. the command
```
bgenix -g file.bgen -list -incl-range 11:3500000-6500000
```
will list all variants in the given range, while
```
bgenix -g file.bgen -vcf -incl-range 11:3500000-6500000
```
will output a VCF file for that region.

**Note:** bgenix writes its output to stdout.  A full command line will therefore often redirect the output to a file, as in:
```
bgenix -g file.bgen -incl-range 11:3500000-6500000 > output.bgen
```
or pipe it to another command, as in:
```
bgenix -g file.bgen -incl-range 11:3500000-6500000 | qctool -g - -filetype bgen -snp-stats -osnp stats.txt
```

Basic usage
====

bgenix can be used to construct an index file like this:
```sh
bgenix -g file.bgen -index
```
This creates an index file called `file.bgen.bgi` containing the index.  Subsequently, data can be retrieved by further calls, e.g.
```sh
bgenix -g file.bgen -incl-range 11:0-1000000 > region.bgen
```
As the command suggests this outputs a bgen file containing only data in the specified region.

See `bgenix -help` for a full list of supported options.  See [[The bgenix index file format]] for a fuller description of the index file format.

Listing variants
====

If the `-list` option is given, bgenix will list variants instead of outputting a bgen file.  For example, using the file `complex.bgen` included in the `example/` folder in the bgen repository, the command:

`bgenix example/complex.bgen -incl-range 01:0- -list`

produces this output:
```
\# bgenix: started 2016-07-06 09:01:15
alternate_ids	rsid	chromosome	position	number_of_alleles	first_allele
.	V1	01	1	2	A	G
V2.1	V2	01	2	2	A	G
.	V3	01	3	2	A	G
.	M4	01	4	3	A	G,T
.	M5	01	5	2	A	G
.	M6	01	7	4	A	G,GT,GTT
.	M7	01	7	6	A	G,GT,GTT,GTTT,GTTTT
.	M8	01	8	7	A	G,GT,GTT,GTTT,GTTTT,GTTTTT
.	M9	01	9	8	A	G,GT,GTT,GTTT,GTTTT,GTTTTT,GTTTTTT
.	M10	01	10	2	A	G
\# bgenix: success, total 10 variants.
```

As described below, another way to list variants is to query the index directly using `sqlite3`, e.g.:

`sqlite3 -header -csv example/complex.bgen.bgi "SELECT * FROM Variant"`

*Note*: you need sqlite3 version 3.8.2 or above for this to work out of the box - see below for more information.

Queries
====

Bgenix can pull out data based on chromosome and position, or by variant identifier.   
In general, a variant will be output if it satisfies at least one *inclusion*  (`-incl-*`) condition, and it does not satisfy any *exclusion* (`-excl-*`) condition.  The relevant options are:

| Query type |     Syntax      | Notes | Example(s) |
| -----------| --------------- | ----- | ------- |
| Include specific identifiers | `-incl-rsids <a> [<b>...]` | For each argument, bgenix looks to see if the argument is the name of a readable file.  If so it opens the file and reads whitespace-separated identifiers from it.  Otherwise it assumes the argument itself is an identifier. | `bgenix x.bgen -incl-rsids rs8176719 myids.txt` |
| Exclude specific identifiers | `-excl-rsids <a> [<b>...]` | " | `bgenix x.bgen -excl-rsids rs8176719 myids.txt` |
| Include range | `-incl-range [chr]:[pos]-[pos]` or `-incl-range <filename>` | One or both positions can be missing in which case a half-open interval is assumed.  As for -incl-rsids, bgenix checks to see if the argument is the name of a readable file and if so reads ranges from it; otherwise it treats the argument itself as a range. | `bgenix x.bgen -incl-range 11:0-10000000` `bgenix x.bgen -incl-range 11:-10000000` `bgenix x.bgen -incl-range 11:10000000-` |
| Exclude range | `-excl-range [chr]:[pos]-[pos]` | " | `bgenix x.bgen -excl-range 11:10000001-` |

Transcoding to other formats
====

By default Bgenix outputs data in the same format it is stored.  This is the most efficient way of using bgenix since it does not involve any decompression or recompression.  However, optionally you can ask bgenix to transcode data to other formats instead.  Two formats are currently supported:

1. the [VCF](Link URL) format (specified using the `-vcf` option).
2. BGEN v1.1 format (specified using the `-v11` option).  Currently this is only supported when the input data is in a specific format, namely BGEN with 'layout=2' blocks, 8-bit probability encoding, and all samples are diploid.

Details of the index file format
====

`bgenix` stores its index in a plain [sqlite3](http://www.sqlite.org) file.  This has a number of advantages.
For example, you can inspect the index file using the sqlite3 command-line program.  E.g. to get a list of variants:
```sh
sqlite3 file.bgen.bgi "SELECT * FROM Variant"
```
You can use sqlite3 to do whatever you want with the file - for example, you could create a new table containing only a specific subset of variants for use in your project (the `-table` option can be used to tell `bgenix` to use this new table).

**Note**: For performance reasons bgenix uses ["WITHOUT ROWID" tables](https://www.sqlite.org/withoutrowid.html) to implement the index.  This means you need sqlite3 version 3.8.2 or greater to inspect the file - otherwise you'll get a message like *"Error: malformed database schema"*".

Alternatively, you can use the `-with-rowid` option when building the index:
```sh
bgenix -g myfile.bgen -index -with-rowid
```
This builds the index in a table with a `rowid` column, useable with earlier versions of sqlite.

See [[The bgenix index file format]] for full details.

Acknowledgements
====
`bgenix` is motivated by and in some respects designed to mimic [tabix](http://www.htslib.org/doc/tabix.html), the htslib tool for indexing tab-delimited files.  The key functionality of `bgenix` is all implemented using the [sqlite3](http://www.sqlite.org) library.  Thank you, sqlite authors!

