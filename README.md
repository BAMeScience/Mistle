# Mistle

About Mistle: w.i.p.

Development of a search engine for spectral libraries consisting of predicted precursor_mass spectra with a focus on large meta-proteomic libraries.

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
