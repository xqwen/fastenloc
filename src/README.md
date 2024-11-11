
## Compile Source Code

To compile the source code in this directory, simply run 
```
make
``` 
The ``Makefile`` is adaptable to your operating system. The resulting executable will be named ``fastenloc``.

By default, the compiled version does **not** include multi-threading support via the OpenMP library. To enable multi-thread parallel processing, use:
```
make USE_OPENMP=1
```
Only the multi-thread version of the executable supports command line option ``-thread``.

## Required Library 

The following libraries are assumed for compiling 

+ GNU scientific library (GSL)
+ Boost library 
+ OpenMP library (required only if multi-thread parallel processing is specified) 

To set up these libraries on MacOS, refer to the instructions [here](MacOS_library_setup.md).
