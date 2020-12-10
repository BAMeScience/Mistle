#include "msp_reader.h"
#include <iostream>
#include <fstream>
#include <sstream>


msp_reader::msp_reader() {

}

vector<spectrum *> msp_reader::read_file(string &path, msp_read_mode read_mode) {

    vector<spectrum *> spectrum_list;
    fstream infile;

    infile.open(path, ios::in);
    if (!infile) {
        cerr << "Could not open file at " << path << endl;
    }

    string line;
    while (!infile.eof()) {

        string tag, value;
        spectrum *c_spectrum = nullptr;
        while (tag != "Num peaks") { // what if no colon -> colon_pos == string::npos
            if (infile.eof())
                exit(12);
            getline(infile, line);

            // split up line to identify comment tags
            size_t colon_pos = line.find(':');

            tag = line.substr(0, colon_pos);
            value = line.substr(colon_pos + 2, string::npos);

            // parse information
            if (tag == "Name") {
                c_spectrum = new spectrum();
                c_spectrum->name = value;
                c_spectrum->peptide = value.substr(0, value.find('/'));
                c_spectrum->charge = stoi(value.substr(value.rfind('/') + 1, string::npos));
            } else if (tag == "MW") {
                c_spectrum->precursor_mass = stof(value);
            } else if (tag == "Comment") {

            }
        }
        //else case: tag = Num peak_positions
        int num_peaks = stoi(value);
        for (int i = 0; i < num_peaks; ++i) {
            getline(infile, value, '\t');
            float peak = stof(value);
            getline(infile, value, '\t');
            float intensity = stof(value);

            getline(infile, value);




            // parse peak_positions and intensities

            /*std::size_t tab_pos = line.find('\t');
            float peak = stof(line.substr(0, tab_pos));
            float intensity = stof(line.substr(tab_pos + 1, line.find('\t')));
            */
            if (c_spectrum->name == "AMANLLSNILNENR/2") {
                cout << peak << " " << intensity << endl;
            }


            c_spectrum->peak_positions.push_back(peak);
            c_spectrum->intensities.push_back(intensity);
        }
        //c_spectrum->intensity_bin_spanning_factor = -1.f; //TODO figure out if neighbor_spanning here
        c_spectrum->bin_peaks(true,true);
        spectrum_list.push_back(c_spectrum);

    }


    infile.close();
    return spectrum_list;
}
