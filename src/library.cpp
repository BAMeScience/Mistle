#include "library.h"

#include <utility>
#include "msp_reader.h"

library::library() {

}

library::library(string path) {
    load_library_from_file(path);
}

bool library::load_library_from_file(string path) {

    //todo multiple file extensions, if endswith .msp do ...
    spectrum_list = msp_reader::read_file(std::move(path));

    return true;
}
