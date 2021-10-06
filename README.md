# Mistle

About Mistle: w.i.p.

Mistle is a fast spectral search engine. It uses a fragment-indexing technique and SIMD intrinsics match experimental MS2 spectra to large spectral libraries at a high performance.
## Requirements

C++20

## Usage

### Build

Build mistle's fragment ion index for spectral matching.

    mistle-build [OPTION...] [optional args]

### Search

Search experimental mass spectra in mistle fragment ion index.


    mistle-search [OPTION...] [optional args]

Use *-h* flag to print the help message. 

## Important for usage on Linux
Input files coming from Windows distribution may have a line ending with \r\n (carriage return). Linux requires \n as line end only.
Remove \r character (char 13) using the following commad line
* *tr -d '\r' < FILE.mgf > FILE_FIXED.mgf*
