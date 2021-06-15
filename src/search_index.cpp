#include <iostream>
#include <chrono>
#include <utility>
#include <cxxopts.hpp>
#include "settings.h"
#include "search_manager.h"
#include "thread_pool.h"

using namespace std;

cxxopts::ParseResult parseArgs(int argc, const char* argv[], std::string &search_path, std::string &index_path) {
    try {
        cxxopts::Options options("mistle-search", "Search experimental mass spectra in mistle fragment ion index");

        options.positional_help("[optional args]").show_positional_help();

        options.add_options()
                ("h, help", "Print this help message")
                ("s,search", "search file or directory ", cxxopts::value<std::string>(), "PATH")
                ("i,index", "index directory (must contain config.txt and binary index files)", cxxopts::value<std::string>(), "PATH")
                ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"), "NUM")
                ("p,ppm_tolerance", "precursor mz tolerance given in ppm", cxxopts::value<float>()->default_value("10"), "NUM")
                ("m,mz_tolerance", "precursor mz tolerance (absolut value in Da)", cxxopts::value<float>(), "NUM")
                ("b,bin_size", "bin size for fragment ion binning (in Da)", cxxopts::value<float>()->default_value("1"), "NUM")
                ("neighbors", "number of neighboring bins intensity is carried over (on search spectrum peaks)", cxxopts::value<int>()->default_value("0"), "NUM")
                ("neighbors_intensity_factor", "fraction [0, 1] of intensity carried over to neighboring bin(s)", cxxopts::value<float>()->default_value("0.5"), "NUM");

        options.parse_positional({"search", "index"});

        auto result = options.parse(argc,argv);


        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }
        if (result.count("search")) {
            search_path = result["search"].as<std::string>();
        }
        if (result.count("index")) {
            index_path = result["index"].as<std::string>();
        }
        if (result.count("threads")) {
            settings::num_threads = result["threads"].as<int>();
            settings::parallel = (settings::num_threads > 1);
        }
        if (result.count("mz_tolerance")) {
            settings::mz_tolerance = result["mz_tolerance"].as<float>();
            settings::use_ppm_tolerance = false;
        }
        if (result.count("ppm_tolerance")) {
            settings::use_ppm_tolerance = true;
            settings::ppm_tolerance = result["ppm_tolerance"].as<float>();
            settings::ppm_factor = settings::ppm_tolerance / 1000000.0f;
            if (result.count("mz_tolerance")) {
                cerr << "Precursor mz tolerance given in ppm and Dalton. Please choose one." << endl;
                exit(1);
            }
        }
        if (result.count("bin_size")) {
            settings::bin_size = result["bin_size"].as<float>();
        }
        if (result.count("neighbors")) {
            settings::neighbors = result["neighbors"].as<int>();
            settings::neighbors_intensity_factor = result["neighbors_intensity_factor"].as<float>();
        }



        return result;

    }
    catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, const char* argv[]) {

    cout << "+++ Mistle Search +++" << endl;


    /*
     * Args
     */

#if USE_AVX_2
    std::cout << "USING AVX2" << endl;
#endif
#if USE_AVX_512
    std::cout << "USING AVX512" << endl;
#endif
    std::string search_file = "/home/ynowatzk/data/9MM/mgf/9MM_FASP.mgf";
    std::string index_dir = "./test/";
    parseArgs(argc, argv, search_file, index_dir);

    auto start = chrono::high_resolution_clock::now();


    /*
     * Preparation and Search
     */

    search_manager sm(search_file, index_dir);

    std::cout << "Preparing libraries and indices ..." << std::endl;
    sm.prepare_search_library();
    std::cout << "Loading precursor index" << std::endl;
    auto check_point = chrono::high_resolution_clock::now();
    sm.prepare_precursor_index();
    std::cout << "Loading time (index): " << duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - check_point).count() << " seconds" << std::endl;

    std::cout << "Searching fragment-ion-indices in batches" << endl;
    sm.perform_searches();
    std::cout << "Merging overlapping results" << std::endl;
    sm.merge_matches();
    //std::cout << "Writing results to file" << std::endl;
    sm.save_search_results_to_file(index_dir + "results.csv");


    cout << "Inner search time elapsed: " << sm.get_time_spent_in_inner_search() << " seconds" << endl;


    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " seconds" << endl;

    return 0;

}