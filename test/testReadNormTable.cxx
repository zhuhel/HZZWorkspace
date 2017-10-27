#include <stdlib.h>
#include <string>
#include <iostream>
#include <map>

#include "Hzzws/RooStatsHelper.h"
#include "Hzzws/Helper.h"

using namespace std;
int main(int argc, char** argv){
    map<string, map<string, double> > all_norm_dic;
    Helper::readNormTable("examples/yields_8TeV_200.txt", all_norm_dic);
    cout << all_norm_dic.size() << endl;
    Helper::printDic<double>(all_norm_dic);
    return 1;
}
