#include "HZZWorkspace/PlotHelp.h"
#include "TLatex.h"
#include "TPave.h"
#include "TLine.h"
#include "TMarker.h"
#include "TPad.h"

///********************************************************************************************
///********************************************************************************************

void ATLASLabel(Double_t x,Double_t y,const char* text,Color_t color) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);

  double delx = 0.115*696*gPad->GetWh()/(472*gPad->GetWw());

  l.DrawLatex(x,y,"ATLAS");
  if (text) {
    TLatex p; 
    p.SetNDC();
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+delx,y,text);
  }
}
void myText(Double_t x,Double_t y,Color_t color, const char *text, float size) {

  TLatex l; 
  //l.SetTextAlign(11);
  l.SetTextSize(size); 
  l.SetTextFont(42);
  l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
 

void myBoxText(Double_t x, Double_t y,Double_t boxsize,Int_t mcolor, Int_t boxstyle, Int_t lcolor, const char *text, float tsize) 
{

  TLatex l; l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetTextFont(42);
  l.SetNDC();
  l.DrawLatex(x,y,text);

  Double_t y1=y-0.50*tsize;
  Double_t y2=y+0.30*tsize;
  Double_t x2=x-0.3*tsize;
  Double_t x1=x2-boxsize;

  //printf("x1= %f x2= %f y1= %f y2= %f \n",x1,x2,y1,y2);

  TPave *mbox= new TPave(x1,y1,x2,y2,1,"NDC");

  mbox->SetLineColor(lcolor);
  mbox->SetFillColor(mcolor);
  mbox->SetFillStyle(boxstyle);
  mbox->Draw();

  //TLine mline;
  //mline.SetLineWidth(4);
  //mline.SetLineColor(1);
  //mline.SetLineStyle(1);
  //Double_t y_new=(y1+y2)/2.;
  //mline.DrawLineNDC(x1,y_new,x2,y_new);

}


void myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle, float msize,  const char *text,Float_t tsize, bool error) 
{
  TMarker *marker = new TMarker(x-1.0*tsize, y, mstyle);
  marker->SetMarkerColor(color);  
  marker->SetNDC(true);
  marker->SetMarkerSize(msize);
  marker->Draw();

  TLine* line = new TLine(x-1.0*tsize, y+0.02, x-1.0*tsize, y-0.02);
  line->SetNDC(true);
  line->SetLineStyle(kSolid);
  line->SetLineColor(color);
  line->SetLineWidth(1);
  if (error) line->Draw();

  TLatex l; l.SetTextAlign(12); 
  l.SetTextSize(tsize); 
  l.SetTextFont(42);
  l.SetNDC();
  l.DrawLatex(x,y,text);
}

void myLineText(Double_t x,Double_t y,Int_t color,Int_t lstyle, float lsize,  const char *text,Float_t tsize) 
{
  TLine* line = new TLine(x-1.1*tsize,y,x-0.5*tsize,y);
  line->SetLineStyle(lstyle);
  line->SetLineColor(color);
  line->SetNDC(true);
  line->SetLineWidth(lsize);
  line->Draw();

  TLatex l; l.SetTextAlign(12); 
  l.SetTextSize(tsize); 
  l.SetTextFont(42);
  l.SetNDC();
  l.DrawLatex(x,y,text);
}

