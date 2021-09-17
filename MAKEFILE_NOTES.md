Installed on Mac OS:
* brew install gsl
* brew install boost
* enable openMP with -Xpreprocessor -fopenmp on clang - How can we make this work for both GCC and Clang?
* Add homebrew's header location to includes and homebrews library location to LDFLAGS

Should i be adding include to CXXFLAGS? Maybe.
OH! CPPFLAGS is c pre-processor flags