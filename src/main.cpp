#include <iostream>
#include "spectrum.h"
#include "msp_reader.h"

using namespace std;

int main() {
    cout << "Welcome, welcome" << endl;
    vector<spectrum*> library = msp_reader::read_file(R"(C:\Users\ynowatzk\Desktop\data\pyrococcus_furiosus\PyroFur_Complete_simulatedSpectra\pyrofur_2019.msp)");

    spectrum *first = library[0];
    cout << first->peptide << endl;
    cout << "peaks: ";
    for (float p : first->peaks) {
        cout << p << " ";
    }
    cout << endl << "inten: ";
    for (float i : first->intensities) {
        cout << i << " ";
    }
    cout << endl;

    return 0;
}

//spectrum::spectrum *spectrum = new spectrum();

