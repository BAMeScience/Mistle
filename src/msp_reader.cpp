#include "msp_reader.h"
#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <utility>
#include <memory>

using namespace std;

std::fstream msp_reader::infile;


bool msp_reader::read_file(string &path, vector<spectrum *> &output_spectra, msp_read_mode read_mode) {

    infile.open(path, ios::in);
    if (!infile) {
        cerr << "Could not open file at " << path << endl;
        return false;
    }

    string line;
    while (!infile.eof()) {
        string tag, value;
        spectrum *c_spectrum = nullptr;
        while (tag != "Num peaks") { // what if no colon -> colon_pos == string::npos
            if (infile.eof())
                return false;
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


            c_spectrum->peak_positions.push_back(peak);
            c_spectrum->intensities.push_back(intensity);
        }
        //c_spectrum->intensity_bin_spanning_factor = -1.f; //TODO figure out if neighbor_spanning here
        //c_spectrum->bin_peaks(true,true);
        c_spectrum->bin_peaks_sparse(true, true);
        output_spectra.push_back(c_spectrum);

    }


    infile.close();
    return true;
}

bool msp_reader::read_spectra_from_positions(string &path, vector<precursor *> &precursor_list, vector<spectrum *> &output_spectra) {

    output_spectra.reserve(precursor_list.size());
    infile.open(path);
    for (int i = 0; i < precursor_list.size(); ++i) {

        unsigned long start = precursor_list[i]->offset_begin;
        unsigned long end = precursor_list[i]->offset_end;

        infile.seekg(start);
        std::string s;
        //cout << precursor_list[i]->name << endl;
        //cout << start << " " <<  end << " " << end - start << endl;
        //if (end > 1844674407)
        //    continue;
        s.resize(end - start);
        infile.read(&s[0], end - start);

        output_spectra.push_back(read_spectrum_from_buffer(s).get()); //TODO .get() is messing with pointers again!!! REFACTOR

        delete output_spectra.back();
        output_spectra.pop_back();
    }

    return false;
}

shared_ptr<spectrum> msp_reader::read_spectrum_from_buffer(const string& buffer) {

    string line, tag, value;
    shared_ptr<spectrum> c_spectrum = std::make_shared<spectrum>();

    stringstream ss(buffer);
    while (tag != "Num peaks") { // what if no colon -> colon_pos == string::npos
        if (ss.eof())
            return nullptr;
        getline(ss, line);

        // split up line to identify comment tags
        size_t colon_pos = line.find(':');

        tag = line.substr(0, colon_pos);
        value = line.substr(colon_pos + 2, string::npos);

        // parse information
        if (tag == "Name") {
            c_spectrum = std::make_shared<spectrum>();
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
        getline(ss, value, '\t');
        float peak = stof(value);
        getline(ss, value, '\t');
        float intensity = stof(value);

        getline(ss, value);




        // parse peak_positions and intensities

        /*std::size_t tab_pos = line.find('\t');
        float peak = stof(line.substr(0, tab_pos));
        float intensity = stof(line.substr(tab_pos + 1, line.find('\t')));
        */


        c_spectrum->peak_positions.push_back(peak);
        c_spectrum->intensities.push_back(intensity);
    }
    //c_spectrum->intensity_bin_spanning_factor = -1.f; //TODO figure out if neighbor_spanning here
    //c_spectrum->bin_peaks(true,true);
    c_spectrum->bin_peaks_sparse(true, true);

    return c_spectrum;
}


bool msp_reader::read_file_precursors(string &path, vector<precursor *> &precursor_list) {

    infile.open(path, ios::in);
    if (!infile) {
        cerr << "Could not open file at " << path << endl;
        return false;
    }

    string line;
    while (!infile.eof()) {
        string tag, value;
        precursor *parent;
        while (tag != "Num peaks") { // what if no colon -> colon_pos == string::npos
            if (infile.eof())
                return false;
            getline(infile, line);

            // split up line to identify comment tags
            size_t colon_pos = line.find(':');

            tag = line.substr(0, colon_pos);
            value = line.substr(colon_pos + 2, string::npos);

            // parse information
            if (tag == "Name") {
                parent = new precursor();
                parent->offset_begin = infile.tellg();
                parent->offset_begin -= (line.length() + 1);
                parent->name = value;
                //parent->peptide = value.substr(0, value.find('/'));
                parent->charge = stoi(value.substr(value.rfind('/') + 1, string::npos));
            } else if (tag == "MW") {
                parent->mass = stof(value);
            } else if (tag == "Comment") {

            }
        }
        //else case: tag = Num peak_positions
        int num_peaks = stoi(value);
        for (int i = 0; i < num_peaks; ++i) {
            //Skip lines
            getline(infile, value);
        }
        parent->offset_end = infile.tellg();
        precursor_list.push_back(parent);
    }
    infile.clear();
    infile.seekg(0, ios::end);
    precursor_list.back()->offset_end = infile.tellg();


    infile.close();
    return true;
}

bool msp_reader::read_file_precursors_efficient(string &path, vector<precursor *> &precursor_list) {



    infile.open(path, ios::in);
    if (!infile) {
        cerr << "Could not open file at " << path << endl;
        return false;
    }


    precursor *parent = new precursor();
    string line;
    unsigned int id = 0;
    while (!infile.eof()) {
        getline(infile, line);

        if (line.rfind("Name:", 0) == 0) { // rightmost match, but starting at pos 0 (or earlier), i.e. prefix

            unsigned long start = infile.tellg();
            start -= (line.length() + 1);

            //Push back previous element
            parent->offset_end = start - 1;
            precursor_list.push_back(parent);

            //Init new element
            parent = new precursor();
            parent->id = id;
            parent->offset_begin = start;
            parent->charge = stoi(line.substr(line.rfind('/') + 1, string::npos));
            ++id;
        }

        if (line.rfind("MW:", 0) == 0) {
            parent->mass = stof(line.substr(4, string::npos));
        }

    }

    // Adding last element
    infile.clear();
    infile.seekg(0, ios::end);
    parent->offset_end = infile.tellg();
    precursor_list.push_back(parent);

    // Deleting first element (dummy)
    precursor_list.erase(precursor_list.begin()); //inefficient (but called once)

    infile.close();

    return true;
}

bool msp_reader::read_next_entry_into_buffer(ifstream &f, string &buffer) {
    /*
     * Requires open filestream and reads until end a new entry
     */

    buffer.clear();

    string line;

    if(!getline(f, line)) {
        return false;
    }

    if (line.rfind("Name:", 0) != 0) {
        cout << line << endl;
        cerr << "entry does not start with Name:" << endl;
        return false;
    }
    buffer.append(line + "\n");

    while (getline(f, line)) {
        if (line.rfind("Name:", 0) == 0) { // rightmost match, but starting at pos 0 (or earlier), i.e. prefix

            //Jump back to line beginning and return buffer
            f.seekg(-(line.length() + 1), ios::cur);
            return true;

        }

        buffer.append(line + "\n");

    }


    return true;
}

