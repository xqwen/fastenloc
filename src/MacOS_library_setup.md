## Compile the source code in MacOS

1. Install homebrew 
2. Install GNU Scientific Library
```
brew install gsl
```
3. Install Boost library
```
brew install boost
```
4. Install OpenMP library
```
brew install libomp
```
5. You may need reconfigure your library paths to ensure a smooth compiling. A simple way is to link installed library files to ```/usr/local/include``` and ```/usr/local/lib/```.
6. Replace ``Makefile`` in ``src`` with ``Makefile.MacOS`` provided in this directory (rename it to ``Makefile``) and run 
```
make
```


You will probably receive a few warning messages but the source code will compile. 
