pnicorr
=======

High performance auto-correlation of fMRI data

- Intial version supports .1D | .1D.gz format as exported from [afni](http://afni.nimh.nih.gov)'s SUMA program.
- Uses multi-threaded BLAS to do the correlation
- Allows one to provide an upper RAM limit

From the help:

    usage: pnicorr file.1D[.gz] [options]
     -norm|-nonorm:         normalize rows (or don't)

     -mem=MB:               memory (MB)
                            Smaller means more file activity; computing is done
                            in stages. Default is the same as -mem=4000  (4G)

     -iotype=1D|1Dgz|mat:   1D: same as input (SUMA ascii)
                          1Dgz: same as input, gzipped
                           mat:  matlab .mat file

Really the input can be any matrix of values, where each row is a time-series to be correlated with all other rows. The format is:

     [# comments]
     [int int int ] fl.oat fl.oat ...

     ([] means optional.)

Note that only the metadata columns must be ints and must have no decimal. Data must have a decimal, even if it is just `.0`. This is the quick and dirty way the file is parsed. This is the format saved by the `3dVol2Surf` program from [afni](http://afni.nimh.nih.gov). Other input formats will be supported asap.

Each row is a timeseries with its values separated by whitespace. The rows can be any timeseries data. No spatial information is used (nor available.)



