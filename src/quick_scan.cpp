#include <iostream>
#include "scanner.h"

using namespace std;

int main() {

    cout << "Hello Quick Scan" << endl;

    string directory = "/home/ynowatzk/data/9MM/msp/";

    scanner sc;

    sc.scan_directory(directory);
    return 0;
}