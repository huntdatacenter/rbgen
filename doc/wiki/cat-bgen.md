`cat-bgen` is a tool to concatenate bgen files.

Usage:

```
cat-bgen -g file1.bgen [file2.bgen...] -og concatenated.bgen
```

cat-bgen assumes the bgen files all have the same header information (i.e. they are encoded with the same version of the BGEN format, they have the same number of samples, etc.).  It works by copying header information from the first file to stdout and then copying the data content of each file to stdout ignoring all the other headers.  (Some checking is done to check that headers are compatible.)

`cat-bgen` has a couple of other options:

| Option | Description|
|--------|------------|
| `-clobber` | overwrite an existing output file, even if it already exists. |
| `-omit-sample-identifier-block` | do not write the sample identifier block, even if there is one in the first file specified to `-g` |
| `-set-free-data <string>` | set the BGEN free data area (as defined in [the spec](https://bgenformat.org)) to the specified string. |


