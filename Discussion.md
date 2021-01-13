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


### Implementation
* TODO: consider charges
* TODO: binary searches
