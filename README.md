# GRTOOL 0.1 2015/9/18

A Research-grade tool to study performance of threaded, GPU accelerated rigid rotation of 3D coordinates

This is a very early alpha stage of GRTOOL, a utility to calculate the RMSD of a MD trajectory.

# To make:

$ make install

which will produce the qrtool binary in your current directory. Note that you need the OpenCL library, and a static netCDF library for GRTOOL to compile correctly.


# To run: 

Determine the settings needed in src/defines.h. These are the "parameters" in the thesis, and control the behavior of RMSD calculation. If these are changed, GRTOOL will need to be recompiled (make install). 

Note that the trajectory to process is also in the defines.h file. Unfortunately, the means that GRTOOL needs to be recompiled each time a different trajectory is processed. I currently work around this by making a link to the desired trajectory and calling the link wt.nc:

ln -sf /home/user/md_trajectory.nc wt.nc

The output is currently stored in qalign.out, and is hard coded.



# Some known issues:

1. Does not work well (maybe at all) if the blocking factory is not a clean multiple of the trajectory size. Weird memory stuff happens.

2. Most user friendliness is removed, such as hardcoding output editing defines.h and recompiling. This will be changed in future versions.

3. I have experienced weird hangups very randomly (seems 1/1000 runs), perhaps from issue #1, but I think it has to do more with a memory management bug.


