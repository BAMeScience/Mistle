#ifndef SIMPLE_EXAMPLE_SETTINGS_H
#define SIMPLE_EXAMPLE_SETTINGS_H



class settings {

public:
    /*
     * Run parameters
     */

    static bool parallel;
    static int num_threads;
    static int batch_size;

    /*
     * Search properties
     */
    static float mz_tolerance;
    static bool use_ppm_tolerance;
    static float ppm_tolerance;
    static float ppm_factor;
    //todo add parameters like neighbors bin scoring. sqrt normalization etc.
    static float bin_size;
    static int neighbors;
    static float neighbors_intensity_factor;



};


#endif //SIMPLE_EXAMPLE_SETTINGS_H
