/* ***************************
 * Used to model the signal+interference.
 * */
#ifndef __HZZWS_SAMPLEHISTPARAM_H__
#define __HZZWS_SAMPLEHISTPARAM_H__
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
#include "HZZWorkspace/RooHistParamPdf.h"

using namespace std;
class SampleHistParam : public SampleBase {

    public:
        typedef map<TString, vector<TH1*> > ShapeDic;

        SampleHistParam(const char* name, // used to construct PDF
                const char* input        // root file contains smoothed histograms
                );
        virtual ~SampleHistParam();

        virtual bool setChannel(const RooArgSet&, const char* channelName, bool with_sys);

        virtual RooAbsPdf* getPDF();

private:
        RooHistParamPdf* makeNominalPDF();
        bool getHistograms();



    protected:
        TFile* sig_file_;
        TFile* bkg_file_;
        RooRealVar* mu_;

        //////////////////////////////////////// 
        // Following variables are dependent on category
        //////////////////////////////////////// 
        //
        string obsname;

        TH1* h_sig_;
        TH1* h_hb_;
        TH1* h_hH_;
        TH1* h_bonly_;
        
};
#endif
