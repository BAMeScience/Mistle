#include <iostream>
#include <cxxopts.hpp>
#include "indexing_manager.h"

using namespace std;


cxxopts::ParseResult parseArgs(int argc, const char* argv[], std::vector<std::string> &input_directories, const std::shared_ptr<configuration>& config) {
    try {
        for (int i = 0; i < argc; ++i) {
            config->build_command += argv[i];
            config->build_command += " ";
        }
        config->build_command.pop_back();

        cxxopts::Options options("mistle-build", "Build mistle's fragment ion index for spectral matching");

        options.positional_help("[optional args]").show_positional_help();

        options.add_options()
                ("h, help", "Print this help message")
                ("i,input", "input directory containing mass spectra (.msp format)", cxxopts::value<std::string>(), "PATH")
                ("o,output", "output directory where indices will be generated", cxxopts::value<std::string>(), "PATH")
                ("n,num_indices", "number of buckets the fragment ion index will be split in", cxxopts::value<unsigned int>()->default_value("64"), "NUM")
                ("min_pep_length", "Minimum peptide length for the reference spectrum to be loaded into the index", cxxopts::value<unsigned int>()->default_value("7"), "NUM")
                ("t,threads", "number of threads (experimental)\n - 1 thread for reading, other threads for processing. Has increased RAM costs (try using more threads or GLIBC_TUNABLES=glibc.malloc.tcache_count=0 for compensation)", cxxopts::value<int>()->default_value("1"), "NUM");

        options.parse_positional({"input", "output"});

        auto result = options.parse(argc,argv);


        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            exit(0);
        }
        if (result.count("threads")) {
            config->num_build_threads = result["threads"].as<int>();
        }
        if (result.count("num_indices")) {
            config->num_indices = result["num_indices"].as<unsigned int>();
        }
        if (result.count("input")) {
            // Parse list of input directories (separated by black space)
            std::string dir_list = result["input"].as<std::string>();
            std::string::size_type start_pos = 0;
            for (auto end_pos = 0; (end_pos = dir_list.find(' ', end_pos)) != std::string::npos; ++end_pos)
            {
                input_directories.push_back(dir_list.substr(start_pos, end_pos - start_pos));
                start_pos = end_pos + 1;
            }

            input_directories.push_back(dir_list.substr(start_pos));
        }
        if (result.count("output")) {
            config->idx_path = result["output"].as<std::string>();
        }
        config->minimum_peptide_length = result["min_pep_length"].as<unsigned int>();


        return result;

    }
    catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}


int main(int argc, const char* argv[]) {
    cout << "+++ Mistle Build +++" << endl;

    /*
     * Args
     */
    std::vector<std::string> input_directories;
    std::shared_ptr<configuration> config = std::make_shared<configuration>();

    parseArgs(argc, argv, input_directories, config);

    /*
     * Build indices
     */


    auto start = chrono::high_resolution_clock::now();
    indexing_manager im(input_directories, config);
    im.build_indices();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " seconds" << endl;

    return 0;

}