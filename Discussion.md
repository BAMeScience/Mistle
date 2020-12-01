# Discussion

This document is meant to discuss ideas and progress for the project.

## Outline

### First

Find/Research suitable (rather simplistic) algorithm for spectrum scoring (experimental to predicted)
* dot product (compare spectraST)
* cross correlation (compare Sequest)
* do more research

### Second

Implement datastructures to support precursor_mass spectrum input
* read raw msp. file input into datastructure
* index data and save it compressed
* find smart way to retrieve data pieces from indexed/compressed database

### Third

Run the search, and benchmark matches against spectraST and DB search. On single species and 9MM data set.

* matches (FDR, consensus)
* runtime
* RAM and disk space usage

### Fourth

Pull out the big boys
* Think about validation
* taxonomic/functional profiling
* Let's see how far we can get




