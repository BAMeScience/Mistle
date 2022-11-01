# Mistle Example Usage

A toy example is provided to test the program. A library of simulated mass spectra (predicted by Prosit) of 1000 yeast peptides is used for reference. 3 experimental spectra matching the species are retrieved from 9MM FASP dataset (Study )

## Commands

Open terminal or change into this example directory. To construct the fragment index from the example spectral library, run

    mistle-build -i yeast_1000.msp -o index/ -n 4 -t 1

Perform example searches.
    
    mistle-search -s yeast_exp.mgf -i index/ -o example_results.csv -p 10 -b 0.2

Compare output file to xxx