#include "TColor.h"

#ifndef PLOTHELP_H
#define PLOTHELP_H

void ATLASLabel(Double_t x,Double_t y,const char* text,Color_t color);
void myText(Double_t x,Double_t y,Color_t color, const char *text, float size); 
void myLineText(Double_t x,Double_t y,Int_t color,Int_t lstyle, float lsize, const char *text,Float_t tsize);
void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle, float msize, const char *text,Float_t tsize, bool error=false);
void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor, Int_t boxstyle, Int_t lcolor, const char *text, float tsize);

#endif
