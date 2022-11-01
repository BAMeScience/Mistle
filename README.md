# Mistle

Mistle is a fast spectral search engine. It uses a fragment-indexing technique and SIMD intrinsics to match experimental MS2 spectra to large spectral libraries at a high performance.

## Requirements
Tested only on linux (debian) for the specified versions:

* C++20
* Cmake (version 3.19.3)
* gcc (10.2.0)

## Build

For building the project, please create (mkdir) a separate build directory. Change into the build directory and run:

    cmake /path/to/mistle/
    cmake --build .
    
In order to make use of SIMD instruction AVX2 or AVX512 build with -DAVX_2=ON or -DAVX_512=ON compiler flag. Check if your CPU supports these. If necessary adjust CMakeList.txt according to the preferences of your CPU.

Optionally, export the directory where *mistle* was build to executable PATH in *~/.bashrc* file. Add line
    
    export PATH="/home/$USER/path/to/mistle/build:$PATH"


## Usage

### Mistle build

Build Mistle's fragment ion index from spectral library.

    mistle-build [OPTION...] [optional args]

### Mistle search

Search experimental mass spectra in Mistle's fragment ion index.


    mistle-search [OPTION...] [optional args]

Use *-h* flag to print the help message. Refer to the example [README](example/README.md) and directory to see an example usage and to test the program.

## Output format

Peptide spectrum matches (PSMs) are provided in tab separated format. 
First line is a commend tagged by # name the exact command and parameters used to produce the output. 

Next line is the header listing all tracked attributes (tab separated).

    id	spectrum	charge	hit_rank	match	peptide	isomers	similarity	bias    [...]


Below, all matched experimental spetra are listed and indexed by their scan name and the rank of the matched library spectrum. (Rank R is append with /R to the name). See example [output](example/index/example_results_control.csv).








## Important for usage on Linux


## Known issues

Input files coming from Windows distribution may have a line ending with \r\n (carriage return). Linux requires \n as line end only.
Remove \r character (char 13) using the following commad line
* *tr -d '\r' < FILE.mgf > FILE_FIXED.mgf*

Similarity is the preferred score baseline. For the dot product certain properties (e.g. bias, lgamma score) are not defined or ill-defined
