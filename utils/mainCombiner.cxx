// This executable serves as a wrapper to populate the Combiner class
#include "HZZWorkspace/Combiner.h"
// Populate the Combiner class using iostream file parsing functions
#include <iostream>

///////////////////////////////////
// This is the main function of the workspace framework, used to create a workspace based on your config file.
// For documentation on how to build your config file, please see
// https://gitlab.cern.ch/HZZ/HZZSoftware/HZZWorkspace/wikis/Config/Introduction
//
// Usage:
// > mainCombiner config.ini wsname.root
// > mainCombiner help

// Use the standard library namespace
using namespace std;

int main(int argc, char* argv[]){
    // Print a statement to remind the user of the execution syntax
    if (argc > 1 && string(argv[1]) == "help"){
        cout << argv[0] << " config_file combined.root " << endl;
        return 0;
    }

    // Default configuration file name is "test.ini", in case none is provided
    string configname = "test.ini";

    // Default workspace output file name is "combined.root"
    string out_name = "combined.root";

    // Update configuration file name
    if (argc > 1){ configname = string(argv[1]); }

    // Update output workspace file name
    if(argc > 2) out_name = string(argv[2]);

    // Print the configuration file ID to the command line
    cout << Form( "Configuration file: %s", configname.c_str() ) << endl;

    // Create an instance of Combiner to build the workspace
    Combiner* combiner = new Combiner(out_name.c_str(), configname.c_str());

    // Read the configuration file and populate the Combiner class
    combiner->combine();

    // Clean up
    delete combiner;

    // Exit
    return 0;
}
