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

## Open Questions
### Reproduce SpectraST results
* why are my results different
* how do they check mz-difference in search and why
* weird bin numbers, where does the + 1 come from
* using sqrt(x) to deemphasize dominant peaks -> does spectrast do this with Prosit output too, if so, does that make sense?
* spectraST dot-bias = 0.0 for all psms.
* BIG Q: peak neighbor-spanning for both query target spectra, or which one??
* try run spectrast -sR (for caching) 
*
* Puzzling: 786.912536621094_3658.89820000002_2 match to AMANLLSNILNENR/2
* Why does spectrast have 0.4+ while I have 0.2
*
* peaks in precusor region might be "sequence non-specific neutral losses of the precursor"

### Scoring
* dot-product seems to be the go to method for spectral search
* dot-product oversensitive to dominant peaks (2*peak size => 4*dot)
    * square root intensities
    * dot bias
*
* Later: Tackle mods and chimeric spectra (MS Fragger method?)
* MS Fragger: Complementary ions (?)
* Narrow window searches 
* ProteoStorm = MS-Fragger with 2-step. First tryptic. Then semi-tryptic (just more peaks in the same spectrum)


### Improvements so far
* capture best hit, don't sort all hits (twice)
* don't divide / maginitudes in dot-calc, do it before hand

### The benchmark (pyrofur)
spectrast:
* normal mode: 10min (suspect, that ST expects mz-sorted query file)
* RAM mode: 105 seconds (sp4 dot only, unbinned, req. ca. 1GB RAM)

* my standard search 50 seconds

* my indexed search: 12 seconds, req. ca. 600 MB Ram (approx)


### Post-Christmas TODOs
* Include OpenMS (look at handout on their git) !! Has everything I want
    * challenging BUILD Windows/Linux
    * First try to get external project example running in \OpenMS\share\OpenMS\examples\external_code\
    * Then copy/migrate code and integrate my specific application
    * Compare results, reproducing once again, what I have done and what spectraST does
* Install Tools like Proteostorm, MS-Fragger
* See where we could implement spectral matching into the source code
* Probably, make own fragment ion index with MS-Fragger as role model
* Compare to Proteostorm/MSFragger/Spectrast

### February TODO's
* Get the OpenMS library working
* Test performance on 9MM (BAM server)
* Look into paralization, splitting up ion-mz index by parent masses.
  * look at *buffer watcher scheduling*
  * idea: group queries (feasible number) by parent masses.
  * load ion-mz index specified for the parentmass
  * run all these queries parallel
  * look at SeqAN Dream-Yara
* Parallel computing in registers SIMD
  * computing the dot product on bit vectors in one operation for multiple data points (smaller bit vectors inside)
  * look at knuts VL
  * figure out what scoring suites best for performance
* Mass shift claibration on marker peptide
  * some implementation in OpenMS
* Intensity scaling/scoring
  * check OpenMS here too
  
### Problems over problems
How to make the sub-indexes for 300GB MS2
* load in batches, make bins and index for the mz windows
  * issue: decide on bins beforehand 
  * issue: update FIindexes on the fly
* go through input multiple times and regard only window of peptides and build full subindex in 1 go
  * issue: 2000 bins makes 2000 times going over 300GB
  * issue: decide on bins beforhand
* sort file beforehand, then only bit by bit 
  * issue: ugly + doubles data


## Implementation

### Data fromats (w.i.p.)

precursor_idx.csv
* header: Num: <number precursors> (needed for parsing)
* then precursors (spectrum bookmarks) line by line
  * encoded: ID, RANK, MZ, CHARGE, PEPTIDE
* every info that should be outputted has to be preserved here some how
* For reconstruction, the ranking vector is reconstructed by precursor rank field while reading

frag_idx_<X>.csv
* TODO


### Active todos
* Improve file reading (not stringstream? getline)
* preset borders for precursor idx based on sub idx limits
* make limits floats
* figure out what to do with library class (delete it *-.-*)
* have a static search class (maybe)
* move load() and save() functions to file_writer classes (config, precursor idx)
* try to figure out what is going on with out of bound spectra (limits)