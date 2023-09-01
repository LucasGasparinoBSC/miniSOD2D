# miniSOD

Small driver for running the main kernels of SOD2D. Generates a mesh with a specified number of elements of chosen order. Runs and times the kernels multiple times, then outputs a file containing the timings of each run for both kernels, as well as the average time per kernel. In debug mode, the program will also output the results of the kernels to a file, as well as other relevant information.

## Requirements

* [Make](https://www.gnu.org/software/make/)
* [CMake](https://cmake.org/) (>= 3.10)
* [GCC](https://gcc.gnu.org/) (>= 7.5.0)
* [MPI](https://www.open-mpi.org/) (>= 3.1.3)

`Make` and `CMake` are used to build the project. `GCC` is used as the compiler. `MPI` is required to use the `MPI_Wtime()` function for timing the kernels. `Make` can be substituted by [Ninja](https://ninja-build.org/) if desired. `GCC` can be substituted by the Intel [OneAPI](https://software.intel.com/content/www/us/en/develop/tools/oneapi.html) compiler if desired, which includes an MPI implementation.

## Building

To build the project, run the following commands:

```bash
mkdir build && cd build
cmake ..
make
```

Assuming that `GCC` and `MPI` are properly set up, this should build the project. The executable will be located in the `build/src/app` directory. Debug mode can be activated as follows:

```bash
cmake -DCMAKE_BUILD_TYPE=Debug ..
```

Once again, this generates a lot of debug information, so it is not recommended to use this mode for large meshes.

If you want to use `Ninja` instead of `Make`, run the following commands:

```bash
cmake -GNinja ..
ninja
```

To use the Intel `OneAPI` compiler, make sure that the environment for it is well set. We recommend using [Lmod](https://lmod.readthedocs.io/en/latest/) to manage the environment, as well as the script provided with OneAPI to generate proper modulefiles. The compilation commands are otherwise identical. Compiler flags are set in `cmake/compilerOps.cmake`.

The app and its libraries can also be installed on the system after the build is complete. To do so, run the following command:

```bash
make install
```

Or, if using `Ninja`:

```bash
ninja install
```

we recommend using `ccmake` to configure the install directory, as well as other options. The default install directory is `/usr/local/`, which is NOT a recommended install location.

## Usage

In `app/sources/main.f90`, just change the variable `nelem` to the desired number of elements. The element order is set inside `utils/sources/mod_constants.f90`, under the variable `porder`. Any positive integer is accepted. The program will generate a mesh with `nelem` elements of order `porder`, and run the kernels on it. The number of times the kernels are run is manually set to `100` in the driver, and can be changed if execution takes too long. Just remember to adjust the averaging of the timers to reflect the change.

Once the program is built, run:

```bash
mpirun -np 1 <path to executable>/miniSOD
```

Trying to run with more than 1 processor will crash the app. If the app has been installed, remember to set its environment on both PATH and LD_LIBRARY_PATH. we recommend simply running the app from the build directory, since no environment setting is required.

## Adding libraries

Through the `CMake` infrastructure, it is possible to add libraries to the project. For examples on how to do so, check the `CMakeLists.txt` files in the `src` and `utils` directories. The `CMake` documentation is also a good resource for this.