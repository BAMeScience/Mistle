#include <iostream>
#include <chrono>
#include <utility>
#include <cxxopts.hpp>
#include "settings.h"
#include "search_manager.h"
#include "thread_pool.h"

using namespace std;

cxxopts::ParseResult parseArgs(int argc, const char* argv[]) {
    try {
        cxxopts::Options options("mistle-search", "Search experimental mass spectra in mistle fragment ion index");

        options.positional_help("[optional args]").show_positional_help();

        options.add_options()
                ("h, help", "Print this help message")
                ("s,search", "search file or directory ", cxxopts::value<std::string>())
                ("i,index", "index directory (must contain config.txt and index files )", cxxopts::value<std::string>())
                ("t,threads", "number of threads", cxxopts::value<int>()->default_value("1"));

        options.parse_positional({"search", "input"});

        auto result = options.parse(argc,argv);


        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            exit(0);
        }
        if (result.count("threads")) {
            settings::num_threads = result["threads"].as<int>();
        }

        return result;

    }
    catch (const cxxopts::OptionException& e) {
        std::cout << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}

int main(int argc, const char* argv[]) {

    parseArgs(argc, argv);
    cout << "Hello World Explorer" << endl;


    auto start = chrono::high_resolution_clock::now();

    /*
     * Args
     */

    std::string search_file = "/home/ynowatzk/data/9MM/mgf/9MM_FASP.mgf";
    std::string index_dir = "./test/";


    /*
     * Preparation and Search
     */

    search_manager sm(search_file, index_dir);

    std::cout << "Preparing libraries and indices" << std::endl;
    sm.prepare_search_library();
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