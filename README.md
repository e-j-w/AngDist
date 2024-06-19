# Angular Distribution Calculator

Maintainer: Jonathan Williams

## Description

Calculates the directional distribution for gamma radiation emitted from an axially symmetric oriented source.  Based on code written for the `g77` compiler by K. Starosta to calculate gamma-gamma angular correlations.

## Build instructions

Use `make` to compile.  To run the program from anywhere, move the resulting `ang_dist` executable to any directory under your `$PATH` environment variable.

This shouldn't depend on any external libraries, just the `gfortran` compiler.  Tested and seems to work on Ubuntu 14.04/16.04, Scientific Linux 6, and CentOS 7.

## Usage

The program can be run using command line argments or using an interactive prompt.  to use the interactive prompt, run the program without any arguments:

```
./ang_dist
```

Running the program in this way also prints out the syntax used to run the program with command line arguments, which is:

```
 ./ang_dist I_final I_init L delta sigmaj q2 q4 q6
```

Only the first 3 argments (I_final, I_init, L) are required - omitting later argments will cause default values to be used instead.  The argument list is described in the following table:

|**Argument**|**Description**|
|:---:|:---:|
| I_final | Final spin of the electromagnetic transition. |
| I_init  | Initial spin of the electromagnetic transition. |
| L       | Multipolarity (1=dipole, 2=quadrupole, etc.) of the electromagnetic transition. |
| delta   | The mixing ratio with the L+1 multipole (L+1/L ratio, 0 for no mixing). Default value: 0 |
| sigmaj  | De-orientation parameter specifying the width of the initial distribution of magnetic sub-states.  Default value: 0 |
| q2      | Attenuation factor (multiplicative coefficient) for the 2nd order Legendre polynomial term.  Default value: 1 |
| q4 | Attenuation factor (multiplicative coefficient) for the 4th order Legendre polynomial term.  Default value: 1 |
| q6 | Attenuation factor (multiplicative coefficient) for the 6th order Legendre polynomial term.  Default value: 1 |

If all arguments are specified and an extra argument (of any value) is added at the end, this causes the program to run in a mode where only the a0, a2, a4, etc. coefficents are reported, which is useful for interfacing with scripts.


## Contributors 

* K. Starosta (original codebase)
* J. Williams (more features, port to `gfortran`)
* This code uses components of the [CERNLIB](https://cernlib.web.cern.ch/) math library, which are provided in the cernlib [directory](cernlib/) of this repo under their original license.

## References

* W.D Hamilton et al., 'The Electromagnetic	Interaction in Nuclear Spectroscopy', particularly eq. 12.197

### Notes from porting the code from `g77` to `gfortran`:

* used -std=legacy flag to supress warnings
* used -ffixed-line-length-none to avoid errors caused by lines that go over the default 72 character limit