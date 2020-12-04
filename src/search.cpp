#include <iostream>
#include "search.h"

search::search() {

}

search::search(library *search_lib) : search_lib(search_lib) {

}

bool search::search_target_library(library *target_lib) {
    this->target_lib = target_lib;
    cout << "Begin searching target library" << endl;
    
    return true;
}

