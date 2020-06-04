## Compile the source code in MacOS

1. Install Xcode 
```
xcode-select --install
```
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
brew install llvm
brew install libomp
```
5. Replace ``Makefile`` in ``src`` with ``Makefile.MacOS`` provided in this directory (rename it to ``Makefile``) and run 
```
make
```


You will probably receive a few warning messages but the source code will compile. 

