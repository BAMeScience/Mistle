#include "mgf_reader.h"
#include <fstream>
#include <iostream>

mgf_reader::mgf_reader() {

}

vector<spectrum *> mgf_reader::read_file(string path) {
    //todo w.i.p.
    vector<spectrum *> spectrum_list;
    fstream infile;

    infile.open(path, ios::in);
    if (!infile) {
        cerr << "Could not open file at " << path << endl;
    }

    spectrum *c_spectrum = nullptr;
    while (!infile.eof()) {
        string line;
        getline(infile, line);

        if (line == "END IONS") {
            // Post-process and save the current spectrum
            //c_spectrum->intensity_bin_spanning_factor = -1.f; //TODO figure out if neighbor_spanning here
            c_spectrum->bin_peaks(true,true);
            spectrum_list.push_back(c_spectrum);
            continue;
        }

        if (line == "BEGIN IONS") {
            c_spectrum = new spectrum();
            continue;
        }

        // split up line to identify comment tags
        string tag, value;
        size_t separator_pos = line.find('=');

        if (separator_pos != string::npos) {
            tag = line.substr(0, separator_pos);
            value = line.substr(separator_pos + 1, string::npos);

            //parse information
            if (tag == "TITLE") {
                c_spectrum->name = value;
            } else if (tag == "PEPMASS") {
                c_spectrum->precursor_mass = stof(value); //todo check if that's actually true
            } else if (tag == "RTINSECONDS") {

            } else if (tag == "CHARGE") {
                c_spectrum->charge = stoi(value);
            }
        }
        else {
                // No separator: Assume peak information is noted down in the line
                if (line.empty())
                    continue;
                std::size_t space_pos = line.find(' ');
                if (space_pos == string::npos)
                    continue;
                float pos = stof(line.substr(0, space_pos));
                c_spectrum->peak_positions.push_back(pos);
                float intensity = stof(line.substr(space_pos, string::npos));
                c_spectrum->intensities.push_back(intensity);
            }
        }
    infile.close();
    return spectrum_list;
}
