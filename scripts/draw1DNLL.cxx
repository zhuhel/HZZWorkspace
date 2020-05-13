/*
 * =====================================================================================
 *
 *       Filename:  draw1DNLL.cxx
 *
 *    Description:  get the 1D NLL plot (TGraph) from File 
 *
 *        Version:  1.0
 *        Created:  12/15/2013 10:03:50 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Xiangyang Ju (), xiangyang.ju@gmail.com
 *   Organization:  
 *
 * =====================================================================================
 */
#ifndef draw1DNLL_cxx
#define draw1DNLL_cxx

#include <stdlib.h>
#include "TChain.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TString.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TAxis.h"
#include "TColor.h"
#include "TROOT.h"
#include "TLine.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TTree.h"


#include <string>
#include <iostream>
#include <fstream>
#include "TMath.h"

TGraph* getGraphFromFile(const char* filename, const char* treename, 
        const char* nllName, const char* variableName, 
        double& minMass, double& masslow, double& masshi, double yvalue = 1.0)
{

  TFile* fin = TFile::Open(filename,"read");
  if(!fin){
    cout<< filename<<" doesn't exist"<<endl;
    return NULL;
  }
  TTree* physics = (TTree*) fin->Get(treename);
  if(!physics){
    fin->Close();
    cout<<treename<<" doesn't exist in "<< filename<<endl;
    return NULL;
  }

  vector<double> massvalue;
  vector<double> nllvalue;
  
  double _nll, _mass;
  physics ->SetBranchAddress(nllName,&_nll);
  physics ->SetBranchAddress(variableName,&_mass);
  masslow = -999, masshi = -999; 
  int nentries = physics ->GetEntries();

  double lowMassError, hiMassError;
  double lownll, hinll;
  double cvLow = 99999, cvHi = 99999;
  double minNLL = -999;
  cout<<treename<<" has "<< nentries<<" entries "<<endl;

  for(int ientry = 0; ientry < nentries; ientry++){
    physics ->GetEntry(ientry);
    if(ientry == 0) {
      minNLL = _nll;
      minMass = _mass;
    }else{
      if(TMath::IsNaN(_nll)){
          cout<<_mass<<" NLL nan"<<endl;
          continue;
      }
      double value = 2* (_nll - minNLL);
      value = value<0?0:value;
      massvalue.push_back(_mass);
      nllvalue.push_back(value);

      if(fabs(value- yvalue) < cvLow && _mass < minMass){
        cvLow = fabs(value- yvalue);
        lowMassError = _mass;
        lownll = value;
      }
      if(fabs(value- yvalue) < cvHi && _mass > minMass){
          cvHi = fabs(value- yvalue); 
          hiMassError = _mass;
          hinll = value;
      }
    }
  }
  TGraph* hmassNLL =new TGraph(massvalue.size(),&massvalue[0],&nllvalue[0]);
  cout<<"Low: "<< lowMassError<<" "<<lownll<<endl;
  cout<<"Hi:  "<< hiMassError <<" "<< hinll <<endl;
  //find the errors
  int ntrys = 0;
  double oldnll = lownll;
  double epsilon = 0.01, steps = 0.01;
  while(fabs(hmassNLL->Eval(lowMassError) -  yvalue) > epsilon &&
          (hmassNLL->Eval(lowMassError) -  yvalue)*(oldnll -  yvalue) > 0 &&
          ntrys++ < 100)
  {
    if(lownll >  yvalue) lowMassError += steps;
    else lowMassError -= steps;
  }

  ntrys = 0;
  oldnll = hinll;
  while(fabs(hmassNLL->Eval(hiMassError) -  yvalue) > epsilon
          &&(hmassNLL->Eval(hiMassError) -  yvalue)*(oldnll -  yvalue) > 0 
          && ntrys++ < 100)
  {
    if(hinll >  yvalue) hiMassError -= steps;
    else hiMassError += steps;
  }
  masshi = hiMassError;
  masslow = lowMassError;
  cout<<"Best "<<variableName<<": "<< minMass<<" +"<<masshi-minMass<<" -"<<minMass - masslow<<" "<<endl;
  delete physics;
  fin->Close();

  return hmassNLL;
}
#endif
