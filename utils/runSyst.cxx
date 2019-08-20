#include "HZZWorkspace/SysProd.h"
#include <iostream>


int main(int argc, char **argv) 
{
    if (argc > 1 && string(argv[1]) == "help"){
        cout << argv[0] << " config_file " << endl;
        return 0;
    }
    string config_name("config.ini");
    if (argc > 1){
        config_name = string(argv[1]);
    }

    SysProd *s = new SysProd(config_name.c_str());
    s->process();
  
    delete s;
    return 0;
}

