// A set of helper functions for doing the intergration,
// used in modeling interference in Highmass Search
//
// Example of defining a 2D function in RooFit
// x-axis: mass, y-axis is width
// run_example()
// Helper functions for interference studies.

#ifndef _HZZWS_INTEGRALHELPER_H
#define _HZZWS_INTEGRALHELPER_H

#include <TF1.h>

namespace IntegralHelper{
    double L_gg (double x);
    Double_t getHiggsXSTF1(Double_t *x, Double_t *par);
    Double_t getSignalIntegralTF1(Double_t *x, Double_t *par);
    Double_t getTotalIntegralTF2(Double_t *x, Double_t *par);

    Double_t getHiggsXS(Double_t mH, Double_t width);
    Double_t getSignalIntegral(Double_t mH, Double_t width, Double_t acc_0, Double_t acc_1, Double_t acc_2, Double_t acc_3);
    Double_t getTotalIntegral(Double_t mH, Double_t width, Double_t kappa, Double_t acc_0, Double_t acc_1, Double_t acc_2, Double_t acc_3, Double_t a_0, Double_t a_1, Double_t a_2, Double_t a_3, Double_t a_4, Double_t b_0, Double_t b_1, Double_t b_2, Double_t b_3, Double_t b_4);
}

#endif
