////////////////////////////////////////////////////////
//How to use and compile Miryam's Gaussian Blur Filter//
///////////////////////////////////////////////////////
////////////////
// To compile //
////////////////
1) on EPFL computer room GRB001 : 
Load the nvidia modules as below to make the nvcc compiler and cuda libraries available.
Change path to the right directory, and then use make

$bash
$export PATH=$PATH:/usr/local/cuda/bin/
$export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
$cd Project_chaabouni/
$make 

2) on your own machine : 
Set the right paths for cuda libraries in the following files headers : 
Makefile, helpers.h 

Then type make in the terminal 

////////////////
//To execute //
///////////////
In the terminal use the following command : 

$./main <<inputFileName>> <<outputFileName>>

You can use any image as an input, but it must be in PNG format. 
You have to indicate an output file name. I
f it does not exist it will be created.
If it already exist, it will be overwritten
Example : 

$./main Lena.png LenaBlur.png


CUDA : no __global__ in class

too compare : diff 
diff <(./fluid | gnuplot) <(./fluid-cuda | gnuplot)

diff <(./fluid) <(./fluid-cuda)

Pas same entry --> c'est a cause des flags compilateurs --> meme flags
