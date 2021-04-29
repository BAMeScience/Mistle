#include "settings.h"

bool settings::parallel = false;
int settings::num_threads = 4;
int settings::batch_size = 100000;


float settings::mz_tolerance = 3.0f;
float settings::ppm_tolerance = 10.0f;
float settings::bin_size = 1.0f;