#include "settings.h"

bool settings::parallel = false;
int settings::num_threads = 4;
int settings::batch_size = 100000;
bool settings::load_batches = false;


float settings::mz_tolerance = 3.0f; // By default unused
bool settings::use_ppm_tolerance = true;
float settings::ppm_tolerance = 10.0f;
float settings::ppm_factor = 10.0f / 1000000.0f;
float settings::bin_size = 1.0f;
int settings::num_hit_ranks = 2;

int settings::neighbors = 0;
float settings::neighbors_intensity_factor = 0.5f;

std::string settings::search_path;
std::string settings::index_path;
std::string settings::output_name;


/*
 * Debugging and testing
 */

std::string settings::search_command;
bool settings::save_search_command = true;
bool settings::turn_off_fragment_intensities = false; //when loading fragment, set intensity to 1
bool settings::apply_topX_in_window_denoising = false;

int settings::peaks_per_window = 5;
float settings::window_size = 100.f;
