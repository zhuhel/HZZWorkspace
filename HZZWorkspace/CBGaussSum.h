
#ifndef __HZZWS_CBGAUSSSUM_H__
#define __HZZWS_CBGAUSSSUM_H__

#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooAbsPdf.h>

#include "HZZWorkspace/SampleBase.h"
#include "HZZWorkspace/Coefficient.h"


class CBGaussSum : public SampleBase {

  public:

    explicit CBGaussSum(
            const char* name, // used to construct PDF
            const char* input        // input text file contains parameters
        );
    virtual ~CBGaussSum();

    virtual RooAbsPdf* getPDF();

    virtual bool setChannel(const RooArgSet& _obs, const char* _ch_name, bool with_sys);

    virtual bool addShapeSys(const TString& npName);

  protected:
    std::map<std::string, SampleBase*> signalContainer_;
};
#endif
