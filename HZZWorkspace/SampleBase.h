////////////////////////////////////////////////////////////
// SampleBase defines general interface for a sample
// It knows how to implement the systematics and acquire PDFs
////////////////////////////////////////////////////////////
#ifndef __HZZWS_SAMPLEBASE_H__
#define __HZZWS_SAMPLEBASE_H__
#include <string>
#include <fstream>
#include <map>

#include <TString.h>
#include <RooAbsPdf.h>
#include <RooArgSet.h>
#include <RooArgList.h>
#include <RooAbsReal.h>
#include <RooFitExtensions/RooMCHistConstraint.h>

#include "HZZWorkspace/Coefficient.h"
#include "HZZWorkspace/SysText.h"

using namespace std;

class SampleBase{

    public:
        typedef map< TString, vector< float > >  NormDic;
        virtual ~SampleBase();

        /* tell sample to provide information in _channelName_
         * return false, if no such channel.
         * */
        virtual bool setChannel( const RooArgSet& observable, const char* channelName, bool with_sys );

        /* add shape and normalization systematic*/
        virtual bool addShapeSys( const TString& npName ) { log_warn( "Careful! You are calling addShapeSys (\"%s\") for a sample \"%s\" which uses the SampleBase::addShapeSys implementation - will do NOTHING!", npName.Data(), pdf_name_.c_str() );
            return false;
        }  // MG: Is this supposed to do nothing??
        virtual bool addNormSys(const TString& npName);

        /* return final PDF for this Sample in _channelName_ */
        virtual RooAbsPdf* getPDF();

        /* return MC constraint terms */
        virtual RooAbsPdf* get_mc_constraint();

        virtual RooAbsReal* getCoefficient(const char* customname = "" );
        void addCoefficient( Coefficient* );

        // Use adaptive binning
        void useAdaptiveBinning();

        // Mass point manipulations
        double get_mass() const { return mass_; }
        void set_mass( double mass ) { mass_ = mass; }

        const string& get_pdf_name() { return pdf_name_; }

        const string& get_channel() { return category_name_; }

        // Get the associated obs/nps
        const RooArgList& getGlobal() { return global_obs_list_; }
        const RooArgList& getNuisance() { return nuisance_obs_list_; }

        // To get the stats that is used to make the histrogram
        virtual int getStats() const { return -1; };

        /* name: used to construct PDF
         * nickname: used to name the signal strength, if it's signal
         * */
        SampleBase(const char* name);
    protected:

        string pdf_name_; // used in pdf (e.g. ATLAS_Signal_ggH)
        string nickname_; // used for some variables to keep name short
        bool use_adpt_bin_; // use adaptive binning if true
        double mass_;

        RooArgList global_obs_list_;    // associated global obs
        RooArgList nuisance_obs_list_;  // associated nps

        ////////////////////////////////////////
        // Following variables depend on category
        ////////////////////////////////////////
        string category_name_;
        RooArgList obs_list_ ;
        TString base_name_ ; // name_categoryName
        RooAbsPdf* nom_pdf;
        Coefficient* coef_;

        ////////////////////////////////////////
        // Constraint terms for each binning systematics
        ////////////////////////////////////////
        RooMCHistConstraint* mc_constraint;
};

#endif
