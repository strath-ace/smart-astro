# SMART-ASTRO [![Build Status](https://travis-ci.com/strath-ace-labs/smart-astro.svg?token=BGZD8zC9w6stvzyuFQhA&branch=master)](https://travis-ci.com/strath-ace-labs/smart-astro)



(Strathclyde Mechanical and Aerospace Research Toolbox for Astrodynamics)

`SMART` is a collection of toolboxes developed and maintained since 2015 by the department of Mechanical and Aerospace Engineer of Strathclyde University. `SMART-ASTRO` in particular, is a collection of astrodynamics and mathematical core routines.

`SMART` is aimed at making code development inside the team more efficient, reusable and easy maintainabile.

`SMART` is CMake-based C++ project. It can be used as the basis for new projects. 


Features
------

  - General directory structure common to C++ projects
  - Example of CMake module
  - Testing framework ([Catch](https://www.github.com/philsquared/Catch "Catch Github repository"))
  - Install script (`make install`)
  - CPack script for packaging (`make package`)
  - Automatic API documentation ([Doxygen](http://www.doxygen.org "Doxygen homepage"))
  - Continuous Integration ([Travis CI](https://travis-ci.org/ "Travis CI homepage"))  (PENDING)
  - Example of how to include external dependencies (using `ExternalProject` module)

Requirements
------

To install this project, please ensure that you have installed the following (install guides are provided on the respective websites):

  - [Git](http://git-scm.com)
  - A C++ compiler, e.g., [GCC](https://gcc.gnu.org/), [clang](http://clang.llvm.org/), [MinGW](http://www.mingw.org/)
  - [CMake](http://www.cmake.org)
  - [Doxygen](http://www.doxygen.org "Doxygen homepage") (optional)

`SMART-ASTRO` depends on the following libraries:

  - [CATCH](https://www.github.com/philsquared/Catch) (unit testing library necessary for `BUILD_TESTS` build option, optional)
  - [BOOST](http://www.boost.org/) odeint library (mandatory)

This dependency will be downloaded and configured automagically if not already present locally (requires an internet connection).

Installation
------

Run the following commands to download, build, and install this project.

```shell
    git clone https://github.com/space-art/smart-astro.git
    (Introcude your username and password)
    cd smart-astro
```

Run the `before_install.sh` script to install prerequisites to compile the project components.

```shell
    mkdir build && cd build
    cmake .. && make
```

To install the header files, run the following from within the `build` directory:

```shell
    make install
```

Note that dependencies are installed by fetching them online, in case they cannot be detected on your local system. If the build process fails, check the error log given. Typically, building fails due to timeout. Simply run the `cmake --build .` command once more.

Build options
-------------

You can pass the following, general command-line options when running CMake:

  - `-DCMAKE_INSTALL_PREFIX[=$install_dir]`: set path prefix for install script (`make install`); if not set, defaults to usual locations
  - `-DBUILD_SHARED_LIBS=[on|off (default)]`: build shared libraries instead of static
  - `-DBUILD_MAIN[=on|off (default)]`: build the main-function
  - `-DBUILD_DOCS[=ON|OFF (default)]`: build the [Doxygen](http://www.doxygen.org "Doxygen homepage") documentation ([LaTeX](http://www.latex-project.org/) must be installed with `amsmath` package)
  - `-DBUILD_TESTS[=ON|OFF (default)]`: build tests (execute tests from build-directory using `ctest -V`)
  - `-DBUILD_DEPENDENCIES[=ON|OFF (default)]`: force local build of dependencies, instead of first searching system-wide using `find_package()`

The following command is conditional and can only be set if `BUILD_TESTS = ON`:


requires [GCC](https://gcc.gnu.org/) compiler; execute coverage analysis from build-directory using `make coverage`)

Pass these options either directly to the `cmake ..` build command or run `ccmake ..` instead to bring up the interface that can be used to toggle options.

Contributing
------------

Once you've made your great commits:

  1. [Fork](https://github.com/space-art/smart-astro//fork) `smart-astro`
  2. Create a topic branch - `git checkout -b my_branch`
  3. Push to your branch - `git push origin my_branch`
  4. Create a [Pull Request](http://help.github.com/pull-requests/) from your branch
  5. That's it!

Disclaimer
------

The copyright holders are not liable for any damage(s) incurred due to improper use of `smart-astro`.
