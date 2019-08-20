/* 
 * Bascially the same as SampleBase, but added MC stat
 * */
#ifndef __HZZWS_SAMPLECOUNT_H__
#define __HZZWS_SAMPLECOUNT_H__
#include <string>
#include <fstream>
#include <map>

#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <RooAbsPdf.h>
#include <RooHistPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooProduct.h>
#include <RooAbsReal.h>
#include <RooFitExtensions/RooStarMomentMorph.h>
#include <RooFitExtensions/RooMCHistConstraint.h>
#include "RooBinning.h"

#include "HZZWorkspace/SampleBase.h"

using namespace std;
class SampleCount : public SampleBase {

    public:
        typedef map<TString, vector<TH1*> > ShapeDic;

        SampleCount(const char* name // used to construct PDF
                );
        virtual ~SampleCount();

        void setMCCThreshold(float thresh); // if thresh < 0, not use mc constraint

        /* return MC constraint terms */
        virtual RooAbsPdf* get_mc_constraint();

    private:


    protected:
        bool use_mcc_ ; // use MC constraint if true
        float thresh_ ; // threshold value for MC constraint


};
#endif
