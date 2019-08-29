#ifndef __ATLASUTILS_C_
#define __ATLASUTILS_C_

#include <iostream>
#include <cmath>

/// #include "AtlasUtils.h"

#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TPave.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"

/// Put these methods, which used to be in a .c macro,
/// into a class as static methods, to make them accessible 
/// from within pyROOT. 
class AtlasUtils{
public: 
static void ATLAS_LABEL(Double_t x,Double_t y,Color_t color) ;

static TGraphErrors* myTGraphErrorsDivide(TGraphErrors* g1,TGraphErrors* g2);


static TGraphAsymmErrors* myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2);

static TGraphAsymmErrors* myMakeBand(TGraphErrors* g0, TGraphErrors* g1,TGraphErrors* g2);

static void myAddtoBand(TGraphErrors* g1, TGraphAsymmErrors* g2);

static TGraphErrors* TH1TOTGraph(TH1 *h1);

static void myText(Double_t x,Double_t y,Color_t color, const char *text, bool use_ndc = true);

static void myBoxText(Double_t x, Double_t y,Double_t textsize,Int_t mcolor,const char *text) ;

static void myBoxText(Double_t x, Double_t y, const TAttFill& h1 ,const char *text, Double_t textsize) ;

static void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle, const char *text,Float_t msize) ;

static void myLineText(Double_t x, Double_t y, Int_t color, Int_t mstyle, const char* text, Float_t msize);

static void myLineText(Double_t x, Double_t y, const TAttLine& line, const char* text, float textsize = 0.04);

static void myLineMarkerText(Double_t x, Double_t y, const TH1& h1, const char* text);

static TLegend* myLegend(double x1, double y1, double x2, double y2);

static void SetXTitle(TH1* h1, const char* title);

static void AddLine(TH1* h1, double y_value, int color = 11, int style = 2);

// some clever person decided to do #define ATLAS 1, so we can not call our method ATLAS :-( 
static void drawATLAS(double x, double y, int opt = 0);

static void DATA(double x, double y, int cms, double lumi);

};

#endif
