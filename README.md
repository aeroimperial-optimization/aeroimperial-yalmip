aeroimperial-yalmip
====================

A MATLAB toolbox for optimization modeling, forked from [yalmip](https://yalmip.github.io). This fork has extra functionality, modules, and support for the multiple precision solver SDPA-GMP. See below for documented updates, and explore the code (or contact the developers) for undocumented ones.

Please bear in mind that this is research code in active development.

Using SDPA-GMP
--------------
If you want to use the multiple precision solver SDPA-GMP, please update the path to the your `sdpa_gmp` executable file in the function [`path2sdpagmp`](https://github.com/aeroimperial-optimization/aeroimperial-yalmip/blob/develop/extras/path2sdpagmp.m).
To ensure that your `sdpa_gmp` executable can indeed be found, the recommended way of doing this is to use [`setpath2sdpagmp`](https://github.com/aeroimperial-optimization/aeroimperial-yalmip/blob/develop/extras/setpath2sdpagmp.m).

**Important:** To avoid sharing the details of your local configuration with the online world next time you push changes to this repository, it is *highly recommended* that you tell `git` to ignore the changes you made to `path2sdpagmp.m`. To do so, run:
````
$ cd <your-local-path>/aeroimperial-yalmip
$ git update-index --assume-unchanged ./extras/path2sdpagmp.m
````
