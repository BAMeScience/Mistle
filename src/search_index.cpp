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
                ("avx2", "Use avx2 instructions for arithmetic operations of 8 floats simultaneously")
                ("avx512", "Use avx512 instructions for arithmetic operations of 16 floats simultaneously")
                ("m,mz_tolerance", "mz tolerance for candidate spectra", cxxopts::value<float>()->default_value("3.0"), "NUM");

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
        }
        settings::avx2 = result.count("avx2");
        settings::avx512 = result.count("avx512");



        return result;

    }
    catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, const char* argv[]) {

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
    cout << "Hello World Explorer" << endl;


    auto start = chrono::high_resolution_clock::now();




    /*
     * Preparation and Search
     */

    search_manager sm(search_file, index_dir);

    std::cout << "Preparing libraries and indices" << std::endl;
    sm.prepare_search_library();
    std::cout << "Loading precursor index" << std::endl;
    sm.prepare_precursor_index();


    std::cout << "Searching fragment-ion-indices in batches" << endl;
    sm.perform_searches();
    std::cout << "Merging overlapping results" << std::endl;
    sm.merge_matches();
    std::cout << "Writing results to file" << std::endl;
    sm.save_search_results_to_file(index_dir + "results.csv");




    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Total time elapsed: " <<  duration.count() << " seconds" << endl;


    /*indexing_manager im(directory);
    auto start = chrono::high_resolution_clock::now();
    im.build_indices();
    auto stop = chrono::high_resolution_clock::now();
    auto duration = duration_cast<chrono::seconds>(stop - start);
    cout << "Loading Time: " <<  duration.count() << " seconds" << endl;
    */



}