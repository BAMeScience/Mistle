# Mistle

Mistle is a fast spectral search engine. It uses a fragment-indexing technique and SIMD intrinsics to match experimental MS2 spectra to large spectral libraries at a high performance.

## Requirements
Tested only on linux (debian) for the specified versions:

* C++20
* Cmake (version 3.19.3)
* g++ (10.2.1)

## Build

For building the project, please create (mkdir) a separate build directory. Change into the build directory and run:

    cmake /path/to/mistle/
    cmake --build .
    
In order to make use of SIMD instruction AVX2 or AVX512 build with -DAVX_2=ON or -DAVX_512=ON compiler flag. Check if your CPU supports these. If necessary adjust CMakeList.txt according to the preferences of your CPU.

Optionally, export the directory where *mistle* was built as an executable PATH in the *~/.bashrc* file. Add the following line:
    
    export PATH="/home/$USER/path/to/mistle/build:$PATH"


## Usage

### Mistle build

Build Mistle's fragment ion index from spectral library.

    mistle-build -i /path/to/library/ -o /path/to/index/ [optional args]

Required arguments are the input directory, which must contain spectral library files (.msp format), and the output directory for the fragment index. 

### Mistle search

Search experimental mass spectra in Mistle's fragment ion index.


    mistle-search -s /path/to/search_file.mgf -i /path/to/index/ [optional args]

Required arguments are the search file (.mgf or .msp format) and the path to the fragment index. Additionally, output directory and formats can be specified as well as various search parameters. Use *-h* flag to print the help message for more information. Also, refer to the [EXAMPLE README](example/README.md) and the example directory to test the program.

## Output format

Peptide spectrum matches (PSMs) are provided in tab separated format. 
First line (comment tagged by #) names the exact shell command and parameters used to produce the output. 

The next line is the header listing all tracked attributes (tab separated).

    id	spectrum	charge	hit_rank	match	peptide	isomers	similarity	bias    [...]

A large number of scores and statistics are appended as additional columns (marked [...]). The most important scores are ...  

Below, all matched experimental spetra are listed and indexed by their scan name and the rank of the matched library spectrum. (Rank R is appended with /R to the scan name). See example [output](example/index/example_results_control.csv).






## Known issues

### On linux

Input files coming from Windows distributions may have a line ending with \r\n (carriage return). Linux and Mistle require \n as the exclusive line ending.
Remove \r character (char 13) using the following commad line
* *tr -d '\r' < FILE.mgf > FILE_FIXED.mgf*

### Scores 

Similarity is the preferred baseline score. For the dot product certain properties (e.g. bias, lgamma score) are not defined or defined based on the similarity.

