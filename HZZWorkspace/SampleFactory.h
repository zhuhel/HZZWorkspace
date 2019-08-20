#ifndef __HZZWS_SAMPLEFACTORY_H__
#define __HZZWS_SAMPLEFACTORY_H__


#include "HZZWorkspace/SampleBase.h"
#include "HZZWorkspace/SampleCount.h"
#include "HZZWorkspace/SampleHist.h"
#include "HZZWorkspace/SampleHistParam.h"
#include "HZZWorkspace/SampleKeys.h"
#include "HZZWorkspace/CBGauss.h"
#include "HZZWorkspace/CBGaussSum.h"
#include "HZZWorkspace/ParametrizedSample.h"
#include "HZZWorkspace/ExpLandau.h"
#include "HZZWorkspace/AnalyticHMBkg.h"
#include "HZZWorkspace/EFTMorph.h"
#include "HZZWorkspace/SimpleMorph.h"
#include "HZZWorkspace/Helper.h"
#include <map>

typedef SampleBase* (*SBConstructor)(strvec&);

namespace SampleFactory{

  // If you wish to add a new sample type, follow the pattern to add it here, and in Root/SampleFactory.cxx
  // Any questions? Please contact graham.cree@cern.ch, xiangyang.Ju@cern.ch


  SampleBase* FactorySampleCount(strvec& args);
  SampleBase* FactorySampleHist(strvec& args);
  SampleBase* FactorySampleHistParam(strvec& args);
  SampleBase* FactorySampleKeys(strvec& args);
  SampleBase* FactoryCBGauss(strvec& args);
  SampleBase* FactoryCBGaussSum(strvec& args);
  SampleBase* FactoryParametrizedSample(strvec& args);
  SampleBase* FactoryExpLandau(strvec& args);
  SampleBase* FactoryAnalyticHMBkg(strvec& args);
  SampleBase* FactoryEFTMorph(strvec& args);
  SampleBase* FactorySimpleMorph(strvec& args);
  // SampleBase* FactoryYourSampleType(strvec& args);

  SampleBase* CreateSample(const std::string& type, strvec& args);



  bool AddKeysSample(SampleBase*, strvec& args);

  strvec& Categories(strvec* in=NULL);
}

#endif
