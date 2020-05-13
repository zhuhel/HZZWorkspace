#ifndef __HZZWS_LINKDEF_H__
#define __HZZWS_LINKDEF_H__

#include <HZZWorkspace/RooStatsHelper.h>

#ifdef __ROOTCLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;

#pragma link C++ namespace RooStatsHelper;
#pragma link C++ function RooStatsHelper::minimize;


#pragma link C++ class RelativisticBW+;
#pragma link C++ class RelativisticBWInt+;
#pragma link C++ class RooGravitonRBWPdf+;
#pragma link C++ class AnalyticHMBkg+;
#pragma link C++ class RooHistParamPdf+;
// #pragma link C++ class RooLagrangianMorphFunc+;
// #pragma link C++ class RooLagrangianMorphOptimizer+;
// #pragma link C++ class RooHCggfWWMorphFunc+;
// #pragma link C++ class RooHCvbfWWMorphFunc+;
// #pragma link C++ class RooHCggfZZMorphFunc+;
// #pragma link C++ class RooHCvbfZZMorphFunc+;
// #pragma link C++ class RooHCvbfMuMuMorphFunc+;

#endif

#endif
