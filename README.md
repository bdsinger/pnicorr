pnicorr
=======

High performance auto-correlation of fMRI data

- Intial version supports .1D | .1D.gz format as exported from [afni](http://afni.nimh.nih.gov)'s SUMA program.
- Uses multi-threaded BLAS to do the correlation
- Allows one to provide an upper RAM limit

From the help:

    usage: pnicorr file.1D[.gz] [options]
     -nonorm:	do not normalize rows

     -savetext:	save text rather than binary files

     -gzout:	compress the output
		    -gzout1 .. -gzout9 sets level
		    -gzout alone is same as -gzout4

     -mem:	memory (MB)
		    -mem1 .. -mem999999 asks for 1M through nearly 1T.
		    Smaller means more file activity; computing is done 
		    in stages. Default is the same as -mem4000  (4G)

Really the input can be any matrix of values, where each row is a time-series to be correlated with all other rows. The format is:

     [# comments]
     [int int int ] fl.oat fl.oat ...

     ([] means optional.)

Note that metadata is typically restricted to the first few columns as output by the afni/suma program `3dVol2Surf`.

Each row is a timeseries with its values separated by whitespace. `3dVol2Surf` interpolates voxel timeseries data onto each node (vertex) on a MGH FreeSurfer surface. But the rows could be any timeseries data; no spatial information is
used (nor is it supplied by any recognized input arguments).



