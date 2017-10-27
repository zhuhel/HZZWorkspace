#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <exception>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TVectorD.h"
#include "TBranch.h"

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"

#include "Hzzws/Helper.h"
#include "Hzzws/SysProd.h"

using namespace std;
using namespace RooFit;

SysProd::SysProd(const string& configFile) {
readConfig(configFile);
}

SysProd::~SysProd() {
}

void SysProd::readConfig(const string& configFile) {
   Helper::readConfig(configFile.c_str(), '=', p_dic);
   Helper::printDic(p_dic);
}

bool SysProd::checkConfig() {
   bool check(false);
   if(p_dic["main"].count("categories") && p_dic["main"].count("samples") && p_dic["main"].count("path") && p_dic["main"].count("treename")  && p_dic["main"].count("NPlist"))
     check = true;
   // check samples
   vector<string> samples; 
   Helper::tokenizeString(p_dic["main"]["samples"].c_str(), ',', samples);
   for(int i=0; i<(int)samples.size(); ++i){
      if(!p_dic["samples"].count(samples[i].c_str()))
      {log_err("ERROR: look at samples!");  check = false;}
   }
   // check categories
   vector<string> categories;
   Helper::tokenizeString(p_dic["main"]["categories"].c_str(), ',', categories);
   for(int i=0; i<(int)categories.size(); ++i){
      if(!p_dic.count(categories[i].c_str()))
      {log_err("ERROR: look at categories!");  check = false;}
   } 
   if (p_dic["main"].count("outdir")) m_outDir = p_dic["main"]["outdir"];
   else m_outDir="./";
   log_info("configured outDir to be :%s",m_outDir.c_str());

   if (p_dic["main"].count("sysDir")==0) p_dic["main"]["sysDir"]="Systematics";

   return check;
}

void SysProd::process() {

   log_info("Processing...");

   bool check = checkConfig();
   if(!check){
      log_err("ERROR: Check config file!");
      return;
   }

   // read NP list
   Helper::readConfig(p_dic["main"]["NPlist"].c_str(), '=', NP_dic);
   Helper::printDic(NP_dic);
   // dump NP list and dirs into vectors
   vector<string> systDirs; 
   vector<string> outSystHistName;
   vector<string> outNormSystWeightName;
   vector<string> outNormSystHistName;
   typedef map<string, string>::iterator it_type;
   for(it_type iterator = NP_dic["shapeLike"].begin(); iterator != NP_dic["shapeLike"].end(); iterator++) {
      if(iterator->first.compare("JET_JER_SINGLE_NP__1down")==0){
        systDirs.push_back("JET_JER_SINGLE_NP__1up");
        outSystHistName.push_back(iterator->second);
      }
      else{
        systDirs.push_back(iterator->first);
        outSystHistName.push_back(iterator->second);
      }
   } 
   for(it_type iterator = NP_dic["normLike"].begin(); iterator != NP_dic["normLike"].end(); iterator++) {
      outNormSystWeightName.push_back(iterator->first);
      outNormSystHistName.push_back(iterator->second);
   }

    // weight to use
   string weightVarName = p_dic["main"]["weightName"].c_str();
   log_info("INFO: Using weight var: %s",weightVarName.c_str());


   // check shape like NP
   if(!checkNP(systDirs))              { log_err("ERROR: please check shape-like NP in the config file (up/down, etc.)!"); return;} 
   // check norm like NP
   if(!checkNP(outNormSystWeightName)) { log_err("ERROR: please check norm-like NP in the config file (up/down, etc.)!"); return;} 
   
   // get sample list
   vector<string> samplesList, samplesNames; 
   Helper::tokenizeString(p_dic["main"]["samples"].c_str(), ',', samplesList);
   for(int i=0; i<(int)samplesList.size(); ++i){
      samplesNames.push_back(samplesList[i].c_str()); 
      samplesList[i]  = p_dic["samples"][samplesList[i].c_str()].c_str();
      }

   // loop over samples
   for(int i=0; i<(int)samplesList.size(); ++i){
      log_info("Working on %s",(samplesList[i]).c_str());

      // loop over categories
      string obsname=""; 
      vector<string> catList; 
      Helper::tokenizeString(p_dic["main"]["categories"].c_str(), ',', catList);

      int nSystDirs = (int)systDirs.size();
      int nNormSyst = (int)outNormSystWeightName.size();

      vector<TH1*> systHistVec(nSystDirs*(int)catList.size());
      vector<TH1*> systHistVecNorm(nNormSyst*(int)catList.size());

      // vectors to store original shapes for comparison
      vector<TH1*> systHistVecOut(nSystDirs*(int)catList.size());
      vector<TH1*> systHistVecNom;

      vector<double> normalizationValues(nSystDirs*(int)catList.size());
      vector<double> normalizationValuesW(nNormSyst*(int)catList.size());

      // vectors to store NP effect for mean and sigma (high mass)
      vector<double> normalizationValues_mean(nSystDirs*(int)catList.size());
      vector<double> normalizationValuesW_mean(nNormSyst*(int)catList.size());
      vector<double> normalizationValues_sigma(nSystDirs*(int)catList.size());
      vector<double> normalizationValuesW_sigma(nNormSyst*(int)catList.size());


      int nomHistCounter(0);
      int systHistCounter(0);
      int systNormHistCounter(0);


       for(int j=0; j<(int)catList.size(); ++j){
          cout << "cat.:" << catList[j] <<endl;

             // check if smooth
             string smooth="";
             try{
                smooth = p_dic.at(catList[j].c_str()).at("smooth");
             }catch(const out_of_range& orr){log_info("Using binned histogram (not smoothing)!");}

             // get nominal histogram
             string fname  = p_dic["main"]["path"] + "/Nominal/" + samplesList[i].c_str();
 
             vector<string> fNameList;
             vector<string> flist;
             Helper::tokenizeString(samplesList[i].c_str(), ',', flist);

             fNameList.clear();
             for(size_t fI = 0; fI < flist.size(); fI++)
             {
                 fNameList.push_back(p_dic["main"]["path"] + "/Nominal/" + flist[fI]);
                 log_info("Using File: %s",fNameList[fI].c_str());
             }


             string hName = obsname+"-"+samplesNames[i]+"-"+catList[j];      
             systHistVecNom.push_back( getHist(hName, samplesNames[i], catList[j].c_str(), fNameList, obsname, p_dic[catList[j]]["cuts"].c_str(), weightVarName, smooth) );
             hName = obsname+"-"+samplesNames[i]+"-"+catList[j]; //unique name, obsname not empty
             systHistVecNom[nomHistCounter]->SetNameTitle(hName.c_str(), hName.c_str());
             TH1* histNom_norm  = getHist(hName+"_norm", samplesNames[i], catList[j].c_str(), fNameList, obsname, p_dic[catList[j]]["cuts"].c_str(), weightVarName, "");

             // get shape systematics histograms
             for(int s=0; s<nSystDirs; s++){
               fname  = p_dic["main"]["path"] + "/" + p_dic["main"]["sysDir"] + "/" + systDirs[s] + "/" + samplesList[i].c_str();
               fNameList.clear();
               for(size_t fI = 0; fI < flist.size(); fI++)
               {
                   fNameList.push_back(p_dic["main"]["path"] + "/" + p_dic["main"]["sysDir"] +  "/" + systDirs[s] + "/" +flist[fI]);
                   cout<<"Using File for sys: "<<fNameList[fI]<<endl;
               }

               string systHistName  = "";
               if(s%2==0){systHistName = obsname+"-"+outSystHistName[s]+"-"+catList[j]+"-down";}
               else      {systHistName = obsname+"-"+outSystHistName[s]+"-"+catList[j]+"-up";} 
               // get histogram
               systHistVec[systHistCounter] = getHist(systHistName, samplesNames[i], catList[j].c_str(), fNameList, obsname, p_dic[catList[j]]["cuts"].c_str(), weightVarName, smooth);       
               systHistVec[systHistCounter]->SetNameTitle(systHistName.c_str(), systHistName.c_str());
               systHistVecOut[systHistCounter] = (TH1*)systHistVec[systHistCounter]->Clone(); // make a copy for later use

               // get normalization
               if(s!=0 && s%2==1 && outSystHistName[s].find("ATLAS_JER")!=std::string::npos){
                  // symmetriez error for JER
                  double symErr = TMath::Abs(1.-systHistVec[systHistCounter]->Integral()/systHistVecNom[nomHistCounter]->Integral());
                  normalizationValues[systHistCounter-1] = 1. - symErr; // down
                  normalizationValues[systHistCounter]   = 1. + symErr; // up
                  // mean and sigma variations for high mass
                  double symErr_mean  = TMath::Abs(1.-systHistVec[systHistCounter]->GetMean()/systHistVecNom[nomHistCounter]->GetMean());
                  double symErr_sigma = TMath::Abs(1.-systHistVec[systHistCounter]->GetRMS()/systHistVecNom[nomHistCounter]->GetRMS());
                  normalizationValues_mean[systHistCounter-1]  = 1. - symErr_mean;  // down
                  normalizationValues_mean[systHistCounter]    = 1. + symErr_mean;  // up
                  normalizationValues_sigma[systHistCounter-1] = 1. - symErr_sigma; // down
                  normalizationValues_sigma[systHistCounter]   = 1. + symErr_sigma; // up
               }
               else{
                  normalizationValues[systHistCounter] = systHistVec[systHistCounter]->Integral()/systHistVecNom[nomHistCounter]->Integral();
                  // mean and sigma variations for high mass
                  normalizationValues_mean[systHistCounter]  = systHistVec[systHistCounter]->GetMean()/systHistVecNom[nomHistCounter]->GetMean();
                  normalizationValues_sigma[systHistCounter] = systHistVec[systHistCounter]->GetRMS()/systHistVecNom[nomHistCounter]->GetRMS();
               } 

               // scale to nominal hist
               systHistVec[systHistCounter]->Scale(systHistVecNom[nomHistCounter]->Integral()/systHistVec[systHistCounter]->Integral());
               systHistVecOut[systHistCounter]->Scale(systHistVecNom[nomHistCounter]->Integral()/systHistVecOut[systHistCounter]->Integral());

               bool divOk = systHistVec[systHistCounter]->Divide(systHistVecNom[nomHistCounter]);
               if(!divOk) { cout << "ERROR: dividing histograms!" << endl; return;}
               fillEmptyBins(catList[j], systHistVec[systHistCounter]); // 1st argument used to define dimentionality of the hist
              
               systHistCounter++;     
               }

             // get normalization systematics histograms integrals
             // check norm like NP
             fname  = p_dic["main"]["path"] + "/" + p_dic["main"]["sysDir"] + "/NormSystematic/" + samplesList[i].c_str();
             fNameList.clear();
             for(size_t fI = 0; fI < flist.size(); fI++)
             {
                 fNameList.push_back(p_dic["main"]["path"] + "/" + p_dic["main"]["sysDir"] + "/NormSystematic/" +flist[fI]);
                 cout<<"Using File for norm sys: "<<fNameList[fI]<<endl;
             }

             //	if(!checkNormNP(outNormSystWeightName, fname, p_dic["main"]["treename"].c_str())) { cout << "ERROR: please check norm-like NP in the config file!" <<endl; return;} 

             for(int s=0; s<nNormSyst; s++){
                 string systHistName  = "";
                 if(s%2==0){systHistName = obsname+"-"+outNormSystHistName[s]+"-"+catList[j]+"-down";}
                 else      {systHistName = obsname+"-"+outNormSystHistName[s]+"-"+catList[j]+"-up";}

                 string normWeightVar = weightVarName;
                 normWeightVar += "*";
                 normWeightVar += outNormSystWeightName[s].c_str();

                 systHistVecNorm[systNormHistCounter] = getHist(systHistName, samplesNames[i], catList[j].c_str(), fNameList, obsname, p_dic[catList[j]]["cuts"].c_str(), normWeightVar.c_str(), "");

                 normalizationValuesW[systNormHistCounter] = systHistVecNorm[systNormHistCounter]->Integral()/histNom_norm->Integral(); 
                 if (normalizationValuesW[systNormHistCounter]==0 || std::isnan(normalizationValuesW[systNormHistCounter])) normalizationValuesW[systNormHistCounter]=1.0;

                 std::cout<<"set value of "<<systHistName<<" to "<<normalizationValuesW[systNormHistCounter]<<std::endl;

                 // mean and sigma variations for high mass
                 normalizationValuesW_mean[systNormHistCounter]  = systHistVecNorm[systNormHistCounter]->GetMean()/histNom_norm->GetMean(); 
                 normalizationValuesW_sigma[systNormHistCounter] = systHistVecNorm[systNormHistCounter]->GetRMS()/histNom_norm->GetRMS(); 

                 if (normalizationValuesW_mean[systNormHistCounter]==0 || std::isnan(normalizationValuesW_mean[systNormHistCounter])) normalizationValuesW_mean[systNormHistCounter]=1.0;
                 if (normalizationValuesW_sigma[systNormHistCounter]==0 || std::isnan(normalizationValuesW_sigma[systNormHistCounter])) normalizationValuesW_sigma[systNormHistCounter]=1.0;


                 systNormHistCounter++;
             }

      nomHistCounter++;
  
      } // end cat.


      // dump shape systematics histograms
      string sName = m_outDir+"/"+"syst_"+samplesNames[i]+".root";
      TFile* fsyst = new TFile(sName.c_str(), "RECREATE");
      for(int t=0; t<(int)systHistVec.size(); ++t) { systHistVec[t]->Write(); }
      fsyst->Close();
      delete fsyst;

      // dump nominal and down/up shape original histograms
      string fName = m_outDir+"/"+"outputs_"+samplesNames[i]+".root";
      TFile* fout = new TFile(fName.c_str(), "RECREATE");

      for(int t=0; t<(int)systHistVecNom.size(); ++t) { systHistVecNom[t]->Write(); }  
      for(int t=0; t<(int)systHistVecOut.size(); ++t) { systHistVecOut[t]->Write(); }
      fout->Close();
      delete fout;
  
      // write normalization file
      string normFileName = m_outDir+"/"+"norm_"+samplesNames[i]+".txt";
      fillNormFile(normFileName, catList, outSystHistName, outNormSystHistName,
                   normalizationValues, normalizationValuesW);

      if(p_dic["main"].count("doMeanSigma") && p_dic["main"]["doMeanSigma"].compare("true")==0){
         string meanFileName = m_outDir+"/"+"mean_"+samplesNames[i]+".txt";
         fillNormFile(meanFileName, catList, outSystHistName, outNormSystHistName,
                   normalizationValues_mean, normalizationValuesW_mean);

         string sigmaFileName = m_outDir+"/"+"sigma_"+samplesNames[i]+".txt";
         fillNormFile(sigmaFileName, catList, outSystHistName, outNormSystHistName,
                   normalizationValues_sigma, normalizationValuesW_sigma);
      }

 
      // delete
      for(int t=0; t<(int)systHistVecNom.size(); ++t)  {delete systHistVecNom[t];}
      for(int t=0; t<(int)systHistVecOut.size(); ++t)  {delete systHistVecOut[t];}
      for(int t=0; t<(int)systHistVec.size(); ++t)     {delete systHistVec[t];}
      for(int t=0; t<(int)systHistVecNorm.size(); ++t) {delete systHistVecNorm[t];}
      systHistVecNom.clear();
      systHistVecOut.clear();
      systHistVec.clear();
      systHistVecNorm.clear();
      normalizationValues.clear();
      normalizationValuesW.clear();
      normalizationValues_mean.clear();
      normalizationValuesW_mean.clear();
      normalizationValues_sigma.clear();
      normalizationValuesW_sigma.clear();

   } // end samples

    
}

TH1* SysProd::getHist(string hName, string sample, string cat, vector<string> fname, string& obsname, string cuts, string weightLabel, string smooth) {

    TChain* tcut = Helper::loader(fname, p_dic["main"]["treename"].c_str());

    // hard-coded to add cuts
    RooArgSet treeobs;
    getObs(cat, obsname, treeobs); 
    RooArgSet obsAndCut(treeobs);


    // get histogram
    RooArgList obsList(treeobs);
    TH1* hist=nullptr;

    //if we know we won't smooth we can shortcut very easily
    if (smooth.empty() || smooth.find(sample)==string::npos){
        RooRealVar *x = (RooRealVar*)obsList.at(0); 

        std::cout<<"going to draw "<<x->GetName()<<" with weight "<< Form("%s*(%s)",weightLabel.c_str(),cuts.c_str())<<std::endl;


        if (obsList.getSize()==1){
            hist = new TH1F(Form("%s", hName.c_str()),Form("%s", hName.c_str()), x->getBinning().numBins(), x->getMin(), x->getMax());
            tcut->Draw(Form("%s>>%s",x->GetName(),hName.c_str()),Form("%s*(%s)",weightLabel.c_str(),cuts.c_str()));
        }
        else if (obsList.getSize()==2){
            RooRealVar *y = (RooRealVar*)obsList.at(1); 
            hist = new TH2F(Form("%s", hName.c_str()),Form("%s", hName.c_str()), x->getBinning().numBins(), x->getMin(), x->getMax(), y->getBinning().numBins(), y->getMin(), y->getMax());
            tcut->Draw(Form("%s:%s>>%s",x->GetName(),y->GetName(),hName.c_str()),Form("%s*(%s)",weightLabel.c_str(),cuts.c_str()));
        }
    }
    //we are going to smooth, so we need to make fancy datasets
    else {

        //Add variables to cut on
        strvec name_list;
        Helper::getListOfNames(cuts, name_list);
        name_list.push_back(weightLabel.c_str());
        for(auto name : name_list){
            if (!obsAndCut.find(name.c_str())){
                obsAndCut.add(*(new RooRealVar(name.c_str(), name.c_str(), -1E6, 1E6)));
            }
        }

        std::cout<<"variables to be read:"<<std::endl;
        obsAndCut.Print("v");

        RooDataSet *ds = new RooDataSet(Form("%s_RooDataSet", obsname.c_str()), "dataset", 
                obsAndCut, RooFit::Import(*tcut),
                RooFit::Cut(cuts.c_str()),
                RooFit::WeightVar(weightLabel.c_str()));

        RooRealVar *x = (RooRealVar*)obsList.at(0); 

        if(obsList.getSize()==1){     // 1D
            double rho(1.);
            getRho(sample, smooth, rho); 
            RooKeysPdf *keyspdf = new RooKeysPdf(Form("%s_RooKeysPdf", obsname.c_str()), "keyspdf", *x, *ds, RooKeysPdf::NoMirror, rho);
            hist = (TH1*)keyspdf->createHistogram(Form("%s", hName.c_str()), *x, RooFit::Binning(x->getBinning()));
            delete keyspdf;

        }
        else if(obsList.getSize()==2){ // 2D
            RooRealVar *y = (RooRealVar*)obsList.at(1); 
            RooCmdArg arg = RooFit::YVar(*y, RooFit::Binning(y->getBinning()));

            double rho_x(1.), rho_y(1.);                 
            getRhoVec(sample, smooth, rho_x, rho_y); 
            TVectorD rhovec(2);
            rhovec(0) = rho_x;
            rhovec(1) = rho_y;
            RooNDKeysPdf *keyspdf = new RooNDKeysPdf(Form("%s_RooNDKeysPdf", obsname.c_str()), "keyspdf", treeobs, *ds, rhovec, "m");
            hist =(TH1*)keyspdf->createHistogram(Form("%s", hName.c_str()), *x, RooFit::Binning(x->getBinning()), arg); 
            delete keyspdf;
            rhovec.Delete();
        }
        delete ds;
    }
    delete tcut;

    hist->SetName("tmp_hist_name");

    cout<<"cuts :"<<cuts<<endl;
    cout<<"Entries :"<<hist->GetEntries()<<" integral: "<<hist->Integral()<<endl;
    // check if apply variable binning (works for 1D)
    if(p_dic[cat].count("bins")){
        vector<string> tmpBins;
        Helper::tokenizeString(p_dic[cat]["bins"].c_str(), '/', tmpBins);
        double xbins[tmpBins.size()];
        for(int i=0; i<(int)tmpBins.size(); ++i){ xbins[i] = (double)::atof(tmpBins[i].c_str()); }
        return hist->Rebin(tmpBins.size()-1, hName.c_str(), xbins);
    }
    else {
        hist->SetName(hName.c_str());
        return hist;
    }
}

void SysProd::getRho(string sample, string smooth, double &rho){
   rho = 1.;
   vector<string> rho_vec;  
   Helper::tokenizeString(smooth.c_str(), ';', rho_vec);
   for(int i=0; i<(int)rho_vec.size(); ++i){
      if(rho_vec[i].find(sample) != string::npos){
      vector<string> rho_val; 
      Helper::tokenizeString(rho_vec[i].c_str(), ':', rho_val);
      rho = (double)::atof(rho_val.at(1).c_str()); 
      }
   }
}

void SysProd::getRhoVec(string sample, string  smooth, double &rho_x, double &rho_y){
   rho_x = 1.; 
   rho_y = 1.;
   vector<string> rho_vec; 
   Helper::tokenizeString(smooth.c_str(), ';', rho_vec);
   for(int i=0; i<(int)rho_vec.size(); ++i){
      if(rho_vec[i].find(sample) != string::npos){
      vector<string> rho_tmp, rho_val; 
      Helper::tokenizeString(rho_vec[i].c_str(), ':', rho_tmp);
      Helper::tokenizeString(rho_tmp.at(1).c_str(), ',', rho_val);
      rho_x = (double)::atof(rho_val.at(0).c_str()); 
      rho_y = (double)::atof(rho_val.at(1).c_str());  
      }
   }
}

void SysProd::getObs(string cat, string &oname, RooArgSet &treeobs) 
{
    oname = "";
    vector<string> branch;  
    Helper::tokenizeString(p_dic[cat]["observables"], ';', branch);

    for (int i=0; i<(int)branch.size(); ++i){
        vector<string> tmp_obs;
        string tmp_branch;
        Helper::readObservable(branch[i], tmp_obs, tmp_branch);

        if (tmp_obs.size() != 4) {
            cout << "Wrong number of parameters for observable, skipping." << endl;
            return;
        }

        RooRealVar *v = new RooRealVar(tmp_branch.c_str(), tmp_branch.c_str(), atof(tmp_obs.at(2).c_str()), atof(tmp_obs.at(3).c_str()));
        v->setBins( atoi(tmp_obs.at(1).c_str()) ); 
        treeobs.add(*v);
        if((int)branch.size()==1)
          oname = tmp_obs.at(0);
        else if(i==((int)branch.size()-1))  
          oname += tmp_obs.at(0);
        else
          oname += tmp_obs.at(0) + "_";
    }
}

void SysProd::fillEmptyBins(string cat, TH1 *hist){

   RooArgSet treeobs;
   string obsname = "";
   getObs(cat, obsname, treeobs); 
   RooArgList obsList(treeobs);

   if(obsList.getSize()==1){      // 1D
     for(int ibinx = 1; ibinx <= hist->GetNbinsX(); ibinx++) {
        if(hist->GetBinContent(ibinx) == 0.0) {
           hist->SetBinContent(ibinx, 1.); 
           hist->SetBinError(ibinx, sqrt(1.)); 
        }
      }
   }
   else if(obsList.getSize()==2){ // 2D
     for(int ibinx = 1; ibinx <= hist->GetNbinsX(); ibinx++) {
        for(int ibiny = 1; ibiny <= hist->GetNbinsY(); ibiny++) {
           if(hist->GetBinContent(ibinx, ibiny) == 0.0) {
             hist->SetBinContent(ibinx, ibiny, 1.); 
             hist->SetBinError(ibinx, ibiny, 0.5);
           }
        }
      }
   }
}


void SysProd::fillNormFile(const string& normFileName, vector<string> catList, vector<string> outSystHistName, vector<string> outNormSystHistName,
     vector<double> normalizationValues, vector<double> normalizationValuesW){

   std::remove(normFileName.c_str());
   ofstream sysOut(normFileName.c_str());
   if(!sysOut){
      cout<<"ERROR: Cannot create "<< normFileName <<" file"<<endl;
      return;
   }
   else{
      cout<<"Created file: "<< normFileName <<endl;
   }

   int systHistCounter=0;
   int systNormHistCounter=0;

      for(int j=0; j<(int)catList.size(); ++j){
         sysOut << "[" << catList[j] << "]" <<endl;

         for(int s=0; s<(int)outSystHistName.size(); s++){  
            if(s%2==0){sysOut << outSystHistName[s] << " = " << Form("%.6f",normalizationValues[systHistCounter]); }
            else      {sysOut << " " << Form("%.6f",normalizationValues[systHistCounter]) <<endl;}
            systHistCounter++;
         }
       
         for(int s=0; s<(int)outNormSystHistName.size(); s++){
            if(s%2==0){sysOut << outNormSystHistName[s] << " = " << Form("%.6f",normalizationValuesW[systNormHistCounter]); }
            else      {sysOut << " " << Form("%.6f",normalizationValuesW[systNormHistCounter]) <<endl;}
            systNormHistCounter++;
         }
        
      }
      sysOut.close();
}

bool SysProd::checkNP(vector<string> systDirs){
   if(systDirs.size() == 1) { cout << "ERROR: You must have up/down dirs for NP!" <<endl; return false;}

   int countUpDown=0; 
   for(int i=0; i<(int)systDirs.size(); i++){

      if(systDirs[i].compare("JET_JER_SINGLE_NP__1up")==0) continue; // skip JER

      if(systDirs[i].find("1up")!=std::string::npos) ++countUpDown;
      else if(systDirs[i].find("1down")!=std::string::npos) --countUpDown;
   }

   if(countUpDown!=0) {cout << "ERROR: You must have up/down dirs (branches) for each NP!" <<endl; return false;}
   else               {return true;}
}

bool SysProd::checkNormNP(vector<string> outNormSystWeightName, string fName, string tName){

   TFile *f = new TFile(fName.c_str(), "read");
   TTree *input = (TTree*)f->Get(tName.c_str());

   for(int i=0; i<(int)outNormSystWeightName.size(); i++){
      TBranch* branch = (TBranch*)input->GetListOfBranches()->FindObject(outNormSystWeightName[i].c_str());
      if(!branch) { 
         cout << "ERROR: Branch " << outNormSystWeightName[i] << " does not exist in the tree!" <<endl; 
         delete branch;
         delete input;    
         f->Close();
         delete f;
         return false;}
      else {delete branch;} 
   }

   delete input;    
   f->Close();
   delete f;

   return true;
}


