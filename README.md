REASON2
=======

Matlab codes for REASON2, paper can be found at http://arxiv.org/abs/1402.5131


This repository provides Matlab codes for running REASON 2. This paper will be presented in NIPS 2014 http://nips.cc/Conferences/2014/Program/event.php?ID=4676
The long version can be accessed at
http://arxiv.org/abs/1402.5131

The code uses the PROPACK package from IALM to make a fair comparison with IALM. To access the package use the following link
http://perception.csl.illinois.edu/matrix-rank/sample_code.html
and download IALM codes.

Our code consists of two files
1.prepare.m 
you make the matrix M here. 
2.REASON2.m 
This file accesses the sample matrix, runs REASON2 and shows individual reconstruction errors at every iteration of the REASON2. The code can be used to run the algorithm for a specific time span or for a specific number of iterations.

Code is written by Hanie Sedghi. 


