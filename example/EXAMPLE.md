# Mistle Example Usage

A toy example is provided to test the program. A library of simulated mass spectra (predicted by Prosit) of 1000 yeast peptides is used for reference. 3 experimental spectra matching the species are retrieved from 9MM FASP dataset (Study )

## Commands

Change into build directory. Construct fragment index from spectral library:

    mistle-build <path to mistle/examples/index/> -n 4 ...

Perform example searches.
    
    mistle-search ...

Compare output file to xxx