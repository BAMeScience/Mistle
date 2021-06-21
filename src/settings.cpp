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

int settings::neighbors = 0;
float settings::neighbors_intensity_factor = 0.5f;