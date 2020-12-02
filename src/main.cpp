#include <iostream>
#include "spectrum.h"
#include "msp_reader.h"
#include "library.h"

using namespace std;

int main() {
    cout << "Welcome, welcome" << endl;
    string msp_file = R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\PyroFur_Complete_simulatedSpectra\pyrofur_2019.msp)";
    //vector<spectrum*> library = msp_reader::read_file(R"(C:\Users\ynowatzk\Desktop\data\9MM\simulated_spectra\Brevibacillus+laterosporus.msp)");

    library *lib = new library(msp_file);


    spectrum *first = lib->spectrum_list[0];
    cout << first->peptide << endl;
    cout << "peak_positions: ";
    for (float p : first->peak_positions) {
        cout << p << " ";
    }
    cout << endl << "inten: ";
    for (float i : first->intensities) {
        cout << i << " ";
    }
    cout << endl << "bins: ";
    for (float i : first->bins) {
        cout << i << " ";
    }
    cout << endl;
    delete lib;
    return 0;
}

//spectrum::spectrum *spectrum = new spectrum();

