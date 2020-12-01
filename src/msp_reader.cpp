#include "msp_reader.h"
#include <iostream>
#include <fstream>


msp_reader::msp_reader() {

}

vector<spectrum *> msp_reader::read_file(string path, msp_read_mode read_mode) {

    vector<spectrum*> mass_spectra;
    fstream infile;

    infile.open(path, ios::in);
    if (!infile) {
        cerr << "Could not open file at " << path << endl;
    }

    string tag, value;
    spectrum* c_spectrum = nullptr;
    while(tag != "Num peaks") { // what if no colon -> colon_pos == string::npos
        string line;
        getline(infile, line);

        size_t colon_pos = line.find(':');

        tag = line.substr(0, colon_pos);
        value = line.substr(colon_pos + 2, string::npos);
        // parse information
        if (tag == "Name") {
            c_spectrum = new spectrum();
            c_spectrum->peptide = value.substr(0, value.find('/'));
        } else if (tag == "MW") {
            c_spectrum->precursor_mass = stof(value);
        } else if (tag == "Comment") {

        }
    }
    //else case: tag = Num peaks
    int num_peaks = stoi(value);
    for (int i = 0; i < num_peaks; ++i) {
        string line;
        getline(infile, line);

        // parse peaks and intensities

        std::size_t tab_pos = line.find('\t');
        float peak = stof(line.substr(0, tab_pos));
        float intensity = stof(line.substr(tab_pos + 1, line.find('\t')));

        c_spectrum->peaks.push_back(peak);
        c_spectrum->intensities.push_back(intensity);
    }


    mass_spectra.push_back(c_spectrum);

    infile.close();
    return mass_spectra;
}
