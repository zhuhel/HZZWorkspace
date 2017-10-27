/* Shapes are obtained from histograms, implemented as RooHistPdf
 * and shape variations implemented as RooStarMomentMorph
 * */
#ifndef __HZZWS_SAMPLEHIST_H__
#define __HZZWS_SAMPLEHIST_H__
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
#include <RooStarMomentMorph.h>
#include <RooMCHistConstraint.h>
#include "RooBinning.h"

#include "Hzzws/SampleBase.h"

using namespace std;
class SampleHist : public SampleBase {

    public:
        typedef map<TString, vector<TH1*> > ShapeDic;

        SampleHist(const char* name, // used to construct PDF
                const char* input,        // root file contains smoothed histograms
                const char* shape_sys    // root file contains shape variations 
                );
        virtual ~SampleHist();

        virtual bool setChannel(const RooArgSet&, const char* channelName, bool with_sys);
        void setMCCThreshold(float thresh); // if thresh < 0, not use mc constraint
        void setInterpOrder(int o);

        virtual bool addShapeSys(const TString& npName);

        int getStats() const {return nom_hist->GetEntries();};

        virtual RooAbsPdf* getPDF();

    private:
        RooAbsPdf* makeHistPdf(TH1*, const char* base_name, bool is_norm = false);
        void getShapeSys();
        RooStarMomentMorph* createRooStarMomentMorph(const string& outputName);
        //virtual void getExpectedValue();
        float cut_sys(float var);
        bool getNominalHist();
        bool hasNominalHist(const char* hist_name);


    protected:
        bool use_mcc_ ; // use MC constraint if true
        float thresh_ ; // threshold value for MC constraint
        TFile* hist_files_;
        TFile* shape_files_;
        int m_interp; //interpolation order

        //////////////////////////////////////// 
        // Following variables are dependent on category
        //////////////////////////////////////// 
        // RooArgList obs_list_ ;
        string obsname;
        
        TH1* nom_hist; // nominal histogram
        TH1* raw_hist;
        //////////////////////////////////////// 
        // PDF systematics
        //////////////////////////////////////// 
        ShapeDic shapes_dic;
        RooAbsPdf* nom_pdf;
        vector<pair<RooAbsPdf*, RooAbsPdf*> > sysPdfs;
        vector<string> paramNames;

};
#endif
