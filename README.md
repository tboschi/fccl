# Feldman Cousins confidence belt constructor

This project provides header libraries to build multidimensional confidence level regions using the
(Feldman Cousins approach ([arXiv:physics/9711021](https://arxiv.org/abs/physics/9711021)).

Each dimension corresponds to a discovery channel, as multiple conditions may be required to be true
in order to determine special events.
The confidence belt is implemented as a recursive templated object.
Thanks to metaprogramming, N nested vectors are created at compile time, where N is the required
number of channels.

## Requirements

* c++11 enabled compiler
* cmake 3.5


## Compiling examples

As this is a header only library, there is no need to compile any software.
There are some examples in the `src` folder that illustrates how to use the library.
To compile those, the typical `cmake` procedure must be followed
```
mkdir build
cd build/
cmake ..
make
make install
```
The last command copies the built binaries in the build folder, under `build/bin`.
