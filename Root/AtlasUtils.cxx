#include "HZZWorkspace/AtlasUtils.h"

void AtlasUtils::ATLAS_LABEL(Double_t x,Double_t y,Color_t color) 
{
  TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize); 
  l.SetNDC();
  l.SetTextFont(72);
  l.SetTextColor(color);
  l.DrawLatex(x,y,"ATLAS");
}

TGraphErrors* AtlasUtils::myTGraphErrorsDivide(TGraphErrors* g1,TGraphErrors* g2) {
 
  const Int_t debug=0; 

  if (!g1) printf("**myTGraphErrorsDivide: g1 does not exist !  \n"); 
  if (!g2) printf("**myTGraphErrorsDivide: g2 does not exist !  \n"); 


  Int_t n1=g1->GetN();
  Int_t n2=g2->GetN();

  if (n1!=n2) {
   printf("**myTGraphErrorsDivide: vector do not have same number of entries !  \n"); 
  }

  TGraphErrors* g3= new TGraphErrors();

  Double_t  x1=0., y1=0., x2=0., y2=0.;
  Double_t dx1=0.,dy1=0.,       dy2=0.;

  Int_t iv=0;
  for (Int_t i1=0; i1<n1; i1++) {
   for (Int_t i2=0; i2<n2; i2++) {
     //if (debug) printf("**myTGraphErrorsDivide: %d  %d !  \n",i1,i2);

    g1->GetPoint(i1,x1,y1);
    g2->GetPoint(i2,x2,y2);
    if (x1!=x2) {
      //printf("**myTGraphErrorsDivide: %d x1!=x2  %f %f  !  \n",iv,x1,x2);
    }else{
      //if (debug) printf("**myTGraphErrorsDivide: %d x1=x2  %f %f  !  \n",iv,x1,x2);
     dx1  = g1->GetErrorX(i1);
     if (y1!=0) dy1  = g1->GetErrorY(i1)/y1;
     if (y2!=0) dy2  = g2->GetErrorY(i2)/y2;
   
     if (debug)
      printf("**myTGraphErrorsDivide: %d x1=%f x2=%f y1=%f y2=%f  \n",iv,x1,x2,y1,y2);

     if (y2!=0.) g3->SetPoint(iv, x1,y1/y2);
     else        g3->SetPoint(iv, x1,y2);
   
     Double_t e=0.;
     if (y1!=0 && y2!=0) e=std::sqrt(dy1*dy1+dy2*dy2)*(y1/y2); 
     g3->SetPointError(iv,dx1,e);


     if (debug) {
       //Double_t g3y, g3x,g3e;
       //g3->GetPoint(iv, g3y,g3x);
       //g3e=g3->GetErrorY(iv);
       //printf("%d g3y= %f g3e=%f  \n",iv,g3y,g3e);
     }
     iv++;
    }
    //    printf("**myTGraphErrorsDivide: ...next  \n");
   }
  }  
  return g3;

}


TGraphAsymmErrors* AtlasUtils::myTGraphErrorsDivide(TGraphAsymmErrors* g1,TGraphAsymmErrors* g2) {

  const Int_t debug=0; 

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();
  Int_t n1=g1->GetN();
  Int_t n2=g2->GetN();

  if (n1!=n2) {
    printf(" vectors do not have same number of entries !  \n");
   return g3;
  }

  Double_t   x1=0.,   y1=0., x2=0., y2=0.;
  Double_t dx1h=0., dx1l=0.;
  Double_t dy1h=0., dy1l=0.;
  Double_t dy2h=0., dy2l=0.;

  // Double_t* X1 = g1->GetX();
  // Double_t* Y1 = g1->GetY();
  Double_t* EXhigh1 = g1->GetEXhigh();
  Double_t* EXlow1 =  g1->GetEXlow();
  Double_t* EYhigh1 = g1->GetEYhigh();
  Double_t* EYlow1 =  g1->GetEYlow();

  // Double_t* X2 = g2->GetX();
  // Double_t* Y2 = g2->GetY();
  // Double_t* EXhigh2 = g2->GetEXhigh();
  // Double_t* EXlow2 =  g2->GetEXlow();
  Double_t* EYhigh2 = g2->GetEYhigh();
  Double_t* EYlow2 =  g2->GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i,x1,y1);
    g2->GetPoint(i,x2,y2);
    dx1h  = EXhigh1[i];
    dx1l  = EXlow1[i];
    if (y1!=0.) dy1h  = EYhigh1[i]/y1;
    else        dy1h  = 0.;
    if (y2!=0.) dy2h  = EYhigh2[i]/y2;
    else        dy2h  = 0.;
    if (y1!=0.) dy1l  = EYlow1 [i]/y1;
    else        dy1l  = 0.;
    if (y2!=0.) dy2l  = EYlow2 [i]/y2;
    else        dy2l  = 0.;
   
    //if (debug)
    //printf("%d x1=%f x2=%f y1=%f y2=%f  \n",i,x1,x2,y1,y2);
    if (debug)
      printf("%d dy1=%f %f dy2=%f %f sqrt= %f %f \n",i,dy1l,dy1h,dy2l,dy2h,
         std::sqrt(dy1l*dy1l+dy2l*dy2l), std::sqrt(dy1h*dy1h+dy2h*dy2h));

    if (y2!=0.) g3->SetPoint(i, x1,y1/y2);
    else       g3->SetPoint(i, x1,y2);
    Double_t el=0.; Double_t eh=0.;

    if (y1!=0. && y2!=0.) el=std::sqrt(dy1l*dy1l+dy2l*dy2l)*(y1/y2);
    if (y1!=0. && y2!=0.) eh=std::sqrt(dy1h*dy1h+dy2h*dy2h)*(y1/y2);

    if (debug) printf("dx1h=%f  dx1l=%f  el=%f  eh=%f \n",dx1h,dx1l,el,eh);
    g3->SetPointError(i,dx1h,dx1l,el,eh);

  }  
  return g3;

}

TGraphAsymmErrors* AtlasUtils::myMakeBand(TGraphErrors* g0, TGraphErrors* g1,TGraphErrors* g2) {
  // default is g0
    //const Int_t debug=0;

  TGraphAsymmErrors* g3= new TGraphAsymmErrors();

  Double_t  x1=0., y1=0., x2=0., y2=0., y0=0, x3=0.;
  //Double_t dx1=0.;
  Double_t dum;
  for (Int_t i=0; i<g1->GetN(); i++) {
    g0->GetPoint(i, x1,y0);
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);

    // if (y1==0) y1=1;
    //if (y2==0) y2=1;

    if (i==g1->GetN()-1) x2=x1;
    else                 g2->GetPoint(i+1,x2,dum);

    if (i==0)            x3=x1;
    else                 g2->GetPoint(i-1,x3,dum);

    Double_t tmp=y2;
    if (y1<y2) {y2=y1; y1=tmp;}
    //Double_t y3=1.;
    Double_t y3=y0;
    g3->SetPoint(i,x1,y3);

    Double_t binwl=(x1-x3)/2.;
    Double_t binwh=(x2-x1)/2.;
    if (binwl==0.)  binwl= binwh;
    if (binwh==0.)  binwh= binwl;
    g3->SetPointError(i,binwl,binwh,(y3-y2),(y1-y3));

  }
  return g3;

}

void AtlasUtils::myAddtoBand(TGraphErrors* g1, TGraphAsymmErrors* g2) {

  Double_t  x1=0., y1=0.,  y2=0., y0=0;
  //Double_t dx1=0.;
  //Double_t dum;

  if (g1->GetN()!=g2->GetN())
    std::cout << " graphs have not the same # of elements " << std::endl;

  Double_t* EYhigh = g2-> GetEYhigh();
  Double_t* EYlow  = g2-> GetEYlow();

  for (Int_t i=0; i<g1->GetN(); i++) {
    g1->GetPoint(i, x1,y1);
    g2->GetPoint(i, x1,y2);
    
    if ( y1==0 || y2==0 ) { 
      std::cerr << "check these points very carefully : myAddtoBand() : point " << i << std::endl;  
    }
    //    if (y1==0) y1=1;
    //    if (y2==0) y2=1;

    //    if (i==g1->GetN()-1) x2=x1;
    //    else                 g2->GetPoint(i+1,x2,dum);
    //    if (i==0)            x3=x1;
    //    else                 g2->GetPoint(i-1,x3,dum);

    Double_t eyh=0., eyl=0.;
    //if (y1<y2) {y2=y1; y1=tmp;}
    //Double_t y3=1.;

    //printf("%d: y1=%f y2=%f Eyhigh= %f Eylow= %f \n",i,y1,y2,EYhigh[i],EYlow[i]);

    y0=y1-y2;
    if (y0!=0) {
     if (y0>0){
      eyh=EYhigh[i];
      eyh=std::sqrt(eyh*eyh+y0*y0);
      //printf("high: %d: y0=%f eyh=%f  \n",i,y0,eyh);
      g2->SetPointEYhigh(i,eyh);
     } else {
      eyl=EYlow[i];
      eyl=std::sqrt(eyl*eyl+y0*y0);
      // printf("low: %d: y0=%f eyl=%f  \n",i,y0,eyl);
      g2->SetPointEYlow (i,eyl);
     }
    }
  }
  return;

}

TGraphErrors* AtlasUtils::TH1TOTGraph(TH1 *h1){


  if (!h1) std::cout << "TH1TOTGraph: histogram not found !" << std::endl;

 TGraphErrors* g1= new TGraphErrors();

 Double_t x, y, ex, ey;
 for (Int_t i=1 ; i<=h1->GetNbinsX(); i++) {
   y=h1->GetBinContent(i);
   ey=h1->GetBinError(i);
   x=h1->GetBinCenter(i);
   ex=h1->GetBinWidth(i);
   
  //   cout << " x,y = " << x << " " << y << " ex,ey = " << ex << " " << ey << endl;

   g1->SetPoint(i-1,x,y);
   g1->SetPointError(i-1,ex,ey);

 }

 //g1->Print();

 return g1;
}

void AtlasUtils::myText(Double_t x,Double_t y,Color_t color, const char *text, bool use_ndc) {

  Double_t tsize=0.04;
  TLatex l; //l.SetTextAlign(12); 
  l.SetTextSize(tsize); 
  l.SetTextFont(42);
  if(use_ndc) l.SetNDC();
  l.SetTextColor(color);
  l.DrawLatex(x,y,text);
}
 

void AtlasUtils::myBoxText(Double_t x, Double_t y,Double_t textsize,Int_t mcolor,const char *text) 
{

  Double_t tsize=0.06;
  Double_t boxsize = 0.07;

  TLatex l; 
  l.SetTextAlign(12); 
  l.SetTextFont(42);
  l.SetTextSize(textsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);

  Double_t y1=y-0.25*tsize;
  Double_t y2=y+0.25*tsize;
  Double_t x2=x-0.3*tsize;
  Double_t x1=x2-boxsize;

  // printf("x1= %f x2= %f y1= %f y2= %f \n",x1,x2,y1,y2);

  TPave *mbox = new TPave(x1,y1,x2,y2,0,"NDC");

  mbox->SetFillColor(mcolor);
  mbox->SetFillStyle(1001);
  mbox->Draw();

  TLine mline;
  mline.SetLineWidth(4);
  mline.SetLineColor(1);
  mline.SetLineStyle(1);
  Double_t y_new=(y1+y2)/2.;
  mline.DrawLineNDC(x1,y_new,x2,y_new);

}

void AtlasUtils::myBoxText(Double_t x, Double_t y, const TAttFill& h1 ,const char *text, Double_t textsize) 
{
    myBoxText(x, y, textsize, h1.GetFillColor(), text);
}


void AtlasUtils::myMarkerText(Double_t x,Double_t y,Int_t color,Int_t mstyle, const char *text,Float_t msize) 
{
  Double_t tsize=0.06;
  TMarker marker(x-(0.4*tsize),y,8);
  marker.SetMarkerColor(color);  
  marker.SetNDC();
  marker.SetMarkerStyle(mstyle);
  marker.SetMarkerSize(msize);
  marker.Draw();

  TLatex l; 
  l.SetTextAlign(12); 
  // l.SetTextSize(tsize); 
  l.SetNDC();
  l.DrawLatex(x,y,text);
}

void AtlasUtils::myLineText(Double_t x, Double_t y, Int_t color, Int_t mstyle, const char* text, Float_t msize){
    Double_t tsize = 0.06*0.1;
    Double_t x1 = x - msize - tsize;
    TLine line;
    line.SetLineColor(color);
    line.SetLineWidth(2);
    line.SetLineStyle(mstyle);
    line.DrawLineNDC(x1, y, x-tsize, y);

    TLatex l; 
    l.SetTextAlign(12); 
    l.SetTextSize(0.04); 
    l.SetNDC();
    l.DrawLatex(x,y,text);
}

void AtlasUtils::myLineText(Double_t x, Double_t y, const TAttLine& line, const char* text, float /* textsize */)
{
    myLineText(x, y, line.GetLineColor(), line.GetLineStyle(), text, 0.04);
}

void AtlasUtils::myLineMarkerText(Double_t x, Double_t y, const TH1& h1, const char* text)
{
    double msize = 0.04;
    Double_t x1 = x - msize - 0.06*0.1;
    Double_t m_x = x1 + msize/2.0;
    TMarker* marker = new TMarker(m_x,y,8);
    marker->SetMarkerColor(h1.GetMarkerColor());  
    marker->SetNDC();
    marker->SetMarkerStyle(h1.GetMarkerStyle());
    marker->SetMarkerSize(msize-0.01);
    marker->Draw();
    myLineText(x, y, h1, text);
}

TLegend* AtlasUtils::myLegend(double x1, double y1, double x2, double y2)
{
/* https://twiki.cern.ch/twiki/bin/view/AtlasProtected/PubComPlotStyle#Legend */
    TLegend* legend = new TLegend(x1, y1, x2, y2);
    legend->SetBorderSize(0);
    legend->SetFillColor(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    return legend;
}

void AtlasUtils::SetXTitle(TH1* h1, const char* title)
{
    h1->SetXTitle(title);
    h1->GetXaxis()->SetTitleOffset(1.4);
    h1->GetXaxis()->SetTitleFont(42);
    h1->GetXaxis()->SetTitleSize(0.05);
}

void AtlasUtils::AddLine(TH1* h1, double y_value, int color, int style)
{
    double x_low = h1->GetBinLowEdge(1);
    double x_hi = h1->GetBinLowEdge(h1->GetNbinsX()+1);
    TLine equal_line;
    equal_line.SetLineColor(color);
    equal_line.SetLineStyle(style);
    equal_line.SetLineWidth(2);
    equal_line.DrawLine(x_low, y_value, x_hi, y_value);
}

void AtlasUtils::drawATLAS(double x, double y, int opt)
{
    switch(opt){
        case 1:
            myText(x, y, 1, "#bf{#it{ATLAS}} Preliminary");
            break;
        case 2:
            myText(x, y, 1, "#bf{#it{ATLAS}}");
            break;
        default:
            myText(x, y, 1, "#bf{#it{ATLAS}} Internal");
            break;
    }
}

void AtlasUtils::DATA(double x, double y, int cms, double lumi){
    myText(x, y, 1, Form("#sqrt{s} = 13 TeV: #scale[0.55]{#int}Ldt = %.1f fb^{-1}", lumi));
    myText(x, y, 1, Form("#sqrt{s} = %d TeV, %.1f fb^{-1}", cms, lumi));
}
