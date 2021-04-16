#include <iostream>
#include <cxxopts.hpp>
#include "indexing_manager.h"

using namespace std;


cxxopts::ParseResult parseArgs(int argc, const char* argv[], std::string &input_directory, const std::shared_ptr<configuration>& config) {
    try {
        cxxopts::Options options("mistle-build", "Build mistle's fragment ion index for spectral matching");

        options.positional_help("[optional args]").show_positional_help();

        options.add_options()
                ("h, help", "Print this help message")
                ("i,input", "input directory containing mass spectra (.msp format)", cxxopts::value<std::string>(), "PATH")
                ("o,output", "output directory where indices will be generated", cxxopts::value<std::string>(), "PATH")
                ("n,num_indices", "number of buckets the fragment ion index will be split in", cxxopts::value<unsigned int>()->default_value("64"), "NUM")
                ("t,threads", "number of threads (NOT IN USE)", cxxopts::value<int>()->default_value("1"), "NUM");

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
            input_directory = result["input"].as<std::string>();
        }
        if (result.count("output")) {
            config->idx_path = result["output"].as<std::string>();
        }


        return result;

    }
    catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}


int main(int argc, const char* argv[]) {
    cout << "Hello World Builder" << endl;


    /*
     * Args
     */
    std::string input_directory = "/home/ynowatzk/data/9MM/msp";
    std::shared_ptr<configuration> config = std::make_shared<configuration>();

    parseArgs(argc, argv, input_directory, config);


    /*
     * Build indices
     */


    indexing_manager im(input_directory, config);
    auto start = chrono::high_resolution_clock::now();
    im.build_indices();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " seconds" << endl;

    return 0;

}