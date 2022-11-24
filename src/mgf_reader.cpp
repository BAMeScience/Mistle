#include "mgf_reader.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include "settings.h"

using namespace std;

mgf_reader::mgf_reader() {

}

bool mgf_reader::read_file(string path, vector<std::shared_ptr<spectrum>> &output_spectra) {

    fstream infile;

    infile.open(path, ios::in);
    if (!infile) {
        return false;
    }

    std::shared_ptr<spectrum> c_spectrum = std::make_shared<spectrum>();
    while (!infile.eof()) {
        string line;
        getline(infile, line);

        if (line == "END IONS") {
            // Post-process and save the current spectrum
            //c_spectrum->intensity_bin_spanning_factor = -1.f; //TODO figure out if neighbor_spanning here
            //c_spectrum->bin_peaks(true,true); //TODO comment out
            if (settings::apply_topX_in_window_denoising)
                c_spectrum->denoise_mz_window(settings::peaks_per_window, settings::window_size); //TODO this exists only for .mgf search file spectra
            c_spectrum->bin_peaks_sparse(true, true);
            c_spectrum->root_scale_intensities();
            c_spectrum->normalize_intensities(); //TODO put into one somehow
            output_spectra.push_back(c_spectrum);
            continue;
        }

        if (line == "BEGIN IONS") {
            c_spectrum = std::make_shared<spectrum>();
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
        else {  // TODO make sure this works with long numbers. Otherwise implement like in msp reader
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


                /*
                 * alternatively run (but it is slower)
                istringstream ss(line);
                float pos, intensity;
                 ss >> pos >> intensity;
                */
            }
        }
    infile.close();
    return true;
}

bool mgf_reader::read_next_entry_into_buffer(ifstream &f, string &buffer) { //TODO
    /*
     * Requires open filestream and reads until end of a new entry
     */

    buffer.clear();

    string line;

    if(!getline(f, line)) {
        return false;
    }

    if (line == "BEGIN IOONS") {
        cout << line << endl;
        cerr << "entry does not start with BEGIN IONS" << endl;
        return false;
    }
    buffer.append(line + "\n");

    while (getline(f, line)) {
        
        buffer.append(line + "\n");
        if (line.rfind("END IONS", 0) == 0) {
            return true;
        }
    }


    return true;
}


shared_ptr<spectrum> mgf_reader::read_spectrum_from_buffer(const string& buffer) { //TODO

    std::string line, tag, value;
    shared_ptr<spectrum> c_spectrum = std::make_shared<spectrum>();

    stringstream ss(buffer);
    while (getline(ss, line)) { // what if no colon -> colon_pos == string::npos
        if (ss.eof())
            return nullptr;

        // parse information
        if (line == "END IONS") {
            //c_spectrum->bin_peaks_sparse(true, true); // TODO YES NO?
            c_spectrum->root_scale_intensities();
            c_spectrum->normalize_intensities(); 
            return c_spectrum;
        }
        if (line == "BEGIN IONS") {
            continue;
        }

        // split up line to identify comment tags
        size_t separator_pos = line.find('=');

        if (separator_pos != string::npos) {
            tag = line.substr(0, separator_pos);
            value = line.substr(separator_pos + 1, string::npos);

            //parse information
            if (tag == "TITLE") {
                c_spectrum->name = value;
                c_spectrum->peptide = value.substr(0, value.find('/')); //TODO Does this work?
            } else if (tag == "PEPMASS") {
                c_spectrum->precursor_mass = stof(value); //todo check if that's actually true
            } else if (tag == "RTINSECONDS") {

            } else if (tag == "CHARGE") {
                c_spectrum->charge = stoi(value);
            }
        }
        else {  // TODO make sure this works with long numbers. Otherwise implement like in msp reader
                // No separator: Assume peak information is noted down in the line


            if (line.empty())
                continue;

            std::size_t space_pos = line.find_first_of(" \t"); //Finds space or tab
            if (space_pos == string::npos)
                continue;


            float pos = stof(line.substr(0, space_pos));
            c_spectrum->peak_positions.push_back(pos);


            float intensity = stof(line.substr(space_pos, string::npos));
            c_spectrum->intensities.push_back(intensity);
        
        }
    }
    return nullptr;
}



bool mgf_reader::read_file_batch(fstream &infile, vector<std::shared_ptr<spectrum>> &output_spectra, int batch_size) {

    if (!infile) {
        cerr << "Error reading file" << endl;
        exit(1);
    }

    int count = 0;
    std::shared_ptr<spectrum> c_spectrum = std::make_shared<spectrum>();
    while (!infile.eof()) {
        string line;
        getline(infile, line);

        if (line == "END IONS") {
            // Post-process and save the current spectrum
            //c_spectrum->intensity_bin_spanning_factor = -1.f; //TODO figure out if neighbor_spanning here
            //c_spectrum->bin_peaks(true,true); //TODO comment out
            c_spectrum->bin_peaks_sparse(true, true);
            c_spectrum->root_scale_intensities();
            c_spectrum->normalize_intensities(); //TODO put into one somehow
            output_spectra.push_back(c_spectrum);
            ++count;
            if (count==batch_size) {
                return false; //Indicating more batches to come
            }
            continue;
        }

        if (line == "BEGIN IONS") {
            c_spectrum = std::make_shared<spectrum>();
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
        else {  // TODO make sure this works with long numbers. Otherwise implement like in msp reader
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


            /*
             * alternatively run (but it is slower)
            istringstream ss(line);
            float pos, intensity;
             ss >> pos >> intensity;
            */
        }
    }
    infile.close();
    return true; //Indicating last batch was reached

}
