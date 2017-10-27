#include <stdlib.h>
#include <string>
#include <map>
#include <fstream>
#include <iostream>
#include <iomanip> 

#include "Hzzws/Helper.h"

using namespace std;
int main(int argc, char* argv[]){
    string config_file = "scale.txt";
    string input_name = "yields_8TeV_200.txt";
    string output_name = "yields_13TeV_200.txt";
    if (argc > 1 && string(argv[1]) == "help"){
        cout << argv[0] << " scale.txt input.txt out.txt" << endl;
        return 0;
    }
    if (argc > 3) {
        config_file = string(argv[1]);
        input_name = string(argv[2]);
        output_name = string(argv[3]);
    }
    
    map<string, double> scale_dic;
    Helper::readScaleFile(config_file.c_str(), scale_dic);
    map<string, map<string, double> > all_norm_dic;
    Helper::readNormTable(input_name.c_str(), all_norm_dic);
    Helper::printDic<double>( all_norm_dic );
    // get category names
    fstream input_file(input_name.c_str(), fstream::in);
    string line;
    getline(input_file, line);
    vector<string> category_names;
    Helper::tokenizeString(line, '&', category_names);
    
    fstream out_file(output_name.c_str(), fstream::out);
    out_file << line << endl;
    for (const auto& key : all_norm_dic){
        const string& sample_name = key.first;
        const map<string, double>& sample_dic = key.second;
        out_file << sample_name;  
        double scale_value = 1.0;
        try {
            scale_value = scale_dic.at(sample_name);
        } catch (const out_of_range& oor) {
            cout <<"WARNNING: " << sample_name <<" does not have scale factor!!" << endl;
        }
        for (const auto& cat_name : category_names){
            out_file << " & " << setprecision(4) << sample_dic.at(cat_name) * scale_value ;
        }
        out_file << endl;
    }
    out_file.close();
}
