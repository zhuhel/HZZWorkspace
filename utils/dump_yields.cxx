#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"

using namespace std;
int main(int argc, char* argv[] ) {
    if (argc < 2 || string(argv[1]).compare("help")==0 ) {
        cout << argv[0] << " input_file " << endl;
        return 1;
    }
    TFile* file = TFile::Open(argv[1], "read");
    vector<string> tree_names;
    tree_names.push_back("tree_ggF_4mu");
    tree_names.push_back("tree_ggF_2mu2e");
    tree_names.push_back("tree_ggF_2e2mu");
    tree_names.push_back("tree_ggF_4e");
    tree_names.push_back("tree_VBF");
    tree_names.push_back("tree_VH_Lep");
    tree_names.push_back("tree_VH_Had");
    for (const auto& tree_name : tree_names){
        TTree* tree = (TTree*) file->Get(tree_name.c_str());
        TH1F* h1 = new TH1F("h1", "h1", 100, 0, 1000);
        tree->Draw("m4l_constrained >> h1", "weight*(m4l");
        double n_events = 0, errors = 0;
        n_events = h1->IntegralAndError(1, 100, errors);
        printf("%s %.2f +/- %.2f\n", tree_name.c_str(), n_events, errors);
        delete h1;
    }
    file->Close();
    return 0;
}
