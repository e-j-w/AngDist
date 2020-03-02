# Angular Distribution Calculator

Maintainer: Jonathan Williams

Contributors: K. Starosta, J. Williams


## Description

Calculates the directional distribution for gamma radiation emitted from an axially symmetric oriented source.  In part based on `g77` code by K. Starosta for gamma-gamma angular correlations.

The `ang_dist` code prompts the user to input parameters, while the `ang_dist_cmd` code takes command line aguments as parameters (running the code without arguments displays a list of arguments needed).

## Installation

Use `make` to compile.  Tested and seems to work with gfortran on Ubuntu 14.04/16.04 and Scientific Linux 6.  Requires `libmathlib` (in the `cernlib` package on Ubuntu 14.04/16.04 and SL6).
If on SL6, you may need to replace the LIB line in the `Makefile` with:

```
LIB = /usr/lib64/cernlib/2006/lib/libmathlib.so.2_gfortran
```

On CentOS 7, one must manually install the cernlib package (available at http://download-ib01.fedoraproject.org/pub/epel/6/x86_64/Packages/c/cernlib-2006-35.el6.x86_64.rpm) since it is not available in the official repos.

### Notes from porting the code from `g77` to `gfortran`:

* used -std=legacy flag to supress warnings
* used -ffixed-line-length-none to avoid errors caused by lines that go over the default 72 character limit
