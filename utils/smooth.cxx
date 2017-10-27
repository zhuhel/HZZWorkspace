#include "Hzzws/SmoothMan.h"

int main(int argc, char **argv) 
{
    if (argc > 1 && string(argv[1]) == "help"){
        cout << argv[0] << " config_file " << endl;
        return 0;
    }
    string config_name("testSmooth.ini");
    if (argc > 1){
        config_name = string(argv[1]);
    }
    
    SmoothMan *sm = new SmoothMan(config_name.c_str());
    sm->process();
    delete sm;
    return 0;
}
