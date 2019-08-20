// 
//    Description:  Overall control of all systematics.
//    It tells Sample which systematics to add
// 
#ifndef __HZZWS_SYSTEMATICSMANAGER_H__
#define __HZZWS_SYSTEMATICSMANAGER_H__

#include <vector>
#include "SampleBase.h"
#include "Coefficient.h"
#include <RooWorkspace.h>
#include <TString.h>

using namespace std;

class SystematicsManager{
    private:
        vector<TString>* all_nps;

    public:
        SystematicsManager();
        SystematicsManager(const char* fileName);
        virtual ~SystematicsManager();
        void readNPs(const char* fileName);
        vector<TString>* add_sys(SampleBase*);
        inline int totalNP(){return all_nps->size();}
        inline vector<TString>* getNP() { return all_nps; }
};
#endif
