# Mistle

Mistle is a fast spectral search engine. It uses a fragment-indexing technique and SIMD intrinsics to match experimental MS2 spectra to large spectral libraries at a high performance. Find out more about Mistle in our publication:

>**Mistle: bringing spectral library predictions to metaproteomics with an efficient search index**  
> Yannek Nowatzky, Philipp Benner, Knut Reinert, Thilo Muth  
> Bioinformatics, Volume 39, Issue 6, June 2023, btad376, https://doi.org/10.1093/bioinformatics/btad376 

Please use the above citation, if you are using Mistle.

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

Required arguments are the input directory, which must contain spectral library files (.msp or .mgf format), and the output directory for the fragment index. 

### Mistle search

Search experimental mass spectra in Mistle's fragment ion index.


    mistle-search -s /path/to/search_file.mgf -i /path/to/index/ [optional args]

Required arguments are the search file (.mgf or .msp format) and the path to the fragment index. Additionally, output directory and formats can be specified as well as various search parameters. Use *-h* flag to print the help message for more information. Also, refer to the [EXAMPLE README](example/README.md) and the example directory to test the program.

## Output format

Peptide spectrum matches (PSMs) are provided in tab separated format. 
First line (comment tagged by #) names the exact shell command and parameters used to produce the output. 

The next line is the header listing all tracked attributes (tab separated).

    id	spectrum	charge	hit_rank	match	peptide	isomers	similarity	bias    [...]

A large number of scores and statistics are appended as additional columns (marked [...]). A detailed explanation of the scores can be found in the next section.

Below the header, all matched experimental spetra are listed and indexed by their scan name and the rank of the matched library spectrum. (Rank R is appended with /R to the scan name). See example [output](example/example_results_control.csv).

Alternatively, a pin-tab format that is readable by Percolator (Käll *et al.*, 2007) can be produced, listing the same scores as features. To obtain this output format, the user needs to specify the output path (*-o*) during mistle-search with the file extension *.pin*. Note that the library label needs to be set correctly at index construction (1: target, -1: decoy libary) and the *results.pin* files of target and decoy search need to be concatenated or merged before using Percolator. It's recommended to use the this python [script](scripts/merge_pin_output.py) to merge the query results and correctly update delta scores.

## Scores 

*Similarity* is the preferred baseline score, which is a refined version of the normalized dot product based on square root transformed peak intensities. A *bias* measurement highlights how biased the *similarity* is on a few matching peaks, and a *delta_similarity* score describes the *similarity* difference between the top hit and second-best hit. Additionally, an *annotation_similarity* version of these scores exists, which accounts only for peak intensities matching reference peaks. This is useful when the library consists of fewer annotated or predicted peaks and is less noisy than the query spectra. 

As a high-quality discriminant scoring function we suggest the *avg_bias_adjusted_similarity*, which is composed equally of the *similarity* and *annotation similarity* metrics. Specifically, a *bias-adjusted similarity* (*sim2*) is calculated by the product of *similarity* and *(1-bias)* and is averaged between standard and annotation version. This scoring function provides excellent discrimination between target and decoy matches.




## Known issues

### On linux

Input files coming from Windows distributions may have a line ending with \r\n (carriage return). Linux and Mistle require \n as the exclusive line ending.
Remove \r character (char 13) using the following commad line
* *tr -d '\r' < FILE.mgf > FILE_FIXED.mgf*

