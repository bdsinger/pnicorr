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
