# Mistle

About Mistle: w.i.p.

Mistle is a fast spectral search engine. It uses a fragment-indexing technique and SIMD intrinsics match experimental MS2 spectra to large spectral libraries at a high performance.
## Requirements

C++20
Cmake (tested for version 3.19.3)
gcc (10.2.0)

## Build

For building the project, please create (mkdir) a separate directory. From there run:
    cmake /path/to/mistle/
In order to make use SIMD instruction AVX2 or AVX512 build with -DAVX_2=ON or -DAVX_512=ON compiler flag. Check if CPU supports these. If necessary adjust CMakeList.txt according to the CPU.

## Usage

### Mistle build

Build mistle's fragment ion index for spectral matching.

    mistle-build [OPTION...] [optional args]

### Mistle search

Search experimental mass spectra in mistle fragment ion index.


    mistle-search [OPTION...] [optional args]

Use *-h* flag to print the help message. 

## Important for usage on Linux
Input files coming from Windows distribution may have a line ending with \r\n (carriage return). Linux requires \n as line end only.
Remove \r character (char 13) using the following commad line
* *tr -d '\r' < FILE.mgf > FILE_FIXED.mgf*
