
## Compile the source code

To compile the source code in this directory, simply run 
```
make
``` 
The ``Makefile`` is adaptable to your operating system. The resulting executable will be named ``fastenloc``.

By default, the compiled version does **not** include multi-threading support via the OpenMP library. To enable multi-threaded parallel processing, use:
```
make USE_OPENMP=1
```
Only the multi-thread version of the executable supports command line option ``-thread``.

## Required library 

``fastENLOC`` requires GNU scientific library (GSL), Boost library. 
The OpenMP library is also required if multi-threaded parallel processing is desired. 

To set up these libraries on MacOS, refer to [this file](MacOS_library_setup.md).
