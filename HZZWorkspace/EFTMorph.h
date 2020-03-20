/* Based on RooLagrangianMorphFunc
 * */
#ifndef __HZZWS_EFTMORPH_H__
#define __HZZWS_EFTMORPH_H__

#include <string>
#include <vector>

#include "HZZWorkspace/SampleBase.h"
#include "HZZWorkspace/Coefficient.h"

#include "RooAbsPdf.h"
#include "RooArgList.h"
#include <RooLagrangianMorphing/RooLagrangianMorphing.h>
// #include "RooLagrangianMorphFunc.h"

using namespace std;
class EFTMorph : public SampleBase {

    public:

        EFTMorph(const char* name, // used to construct PDF
                const char* configfile, bool _shape_BSM_only=false);
        virtual ~EFTMorph();

        virtual bool setChannel(const RooArgSet&, const char* channelName, bool with_sys);

        virtual RooAbsPdf* getPDF();

        virtual RooAbsReal* getCoefficient(const char* customname="");
        virtual bool addShapeSys(const TString& npName);

    private:

        RooArgSet createHCMorphParaSet(std::string parlist);
        RooArgSet createGeneralMorphParaSet(std::string parlist);
        RooArgSet couplingsDatabase;
        RooArgSet formulaDatabase;
        RooAbsReal* getOverallNormalization();


    protected:

        std::map<std::string, std::map<std::string, std::string> > morph_dic;
        std::map<std::string, Coefficient*> m_sampleCoefMap;
        RooLagrangianMorphFunc* m_eftfunc;

        RooArgList *m_morphcoefs, *m_morphfuncs;

// default is false: norm is BSM sensitive. with this bool to true remove norm dependence on BSM couplings -> only shape is BSM sensitive
        bool onlyShapeBSMsensitive;

};
#endif
