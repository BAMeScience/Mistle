#include <iostream>
#include "spectrum.h"
#include "msp_reader.h"

using namespace std;

int main() {
    cout << "Welcome, welcome" << endl;
    vector<spectrum*> library = msp_reader::read_file("test");

    return 0;
}

//spectrum::spectrum *spectrum = new spectrum();

