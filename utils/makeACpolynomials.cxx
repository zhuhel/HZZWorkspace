#include "TSystem.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TAxis.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "Hzzws/PlotHelp.h"
#include "Rtypes.h"
#include <map>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMath.h"
#include "TROOT.h"

#include <iostream>
#include <fstream>
#include <stdio.h>

int main(){

  std::string minitreeDir = "/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/MiniTrees/Prod_v10/mc/Nominal/";
  //std::string minitreeDir = "/afs/cern.ch/work/d/ddenysiu/public/forGraham/Prod_v05/";
  std::cout<<"looking for files in "<<minitreeDir<<std::endl;

  std::string outTFileName = "ACplots.root";
  std::string outFileName = "polyNorm.txt";
  std::string outFileSIName = "polyNormSI.txt";
  std::string outFileSysName = "polySys.txt";

  float minMH = 120;  //Range to user MC mass points into fit //FIXME
  float maxMH = 130;
  float minPlot = 122; //Range to plot fits //FIXME MELA uses 300, LWA uses 400,  NWA uses 200
  float maxPlot = 128;
  float minPlotY=0.05; //set to -1 to use auto
  float maxPlotY=0.2; //set to -1 to use auto

  // float lumi=1;
  // bool storedLumi=false;

  std::vector<std::string> widths;
  //widths.push_back("NW");
  widths.push_back("ZZ4lep");
  //widths.push_back("w5");
  //widths.push_back("w10");
  //widths.push_back("w15");

  std::vector<std::string> prod;
  prod.push_back("ggH");
  //prod.push_back("HZZ4l");
  //prod.push_back("VBFH");

  std::vector<std::string> prodlabel;
  prodlabel.push_back("ggF");
  //prodlabel.push_back("ggF (15\%)");
  //prodlabel.push_back("VBF");

  std::vector<std::string> categories;
  std::vector<std::string> cuts;
  std::vector<std::string> catlabels;
  std::vector<int> color;
  std::vector<int> marker;


  //Mass Msmt
 
  categories.push_back("ggF_4mu"); color.push_back(kViolet); marker.push_back(24); catlabels.push_back("4#mu");
  categories.push_back("ggF_4e"); color.push_back(kOrange); marker.push_back(25); catlabels.push_back("4e");
  categories.push_back("ggF_2mu2e"); color.push_back(kGreen+2); marker.push_back(26); catlabels.push_back("2#mu2e");
  categories.push_back("ggF_2e2mu"); color.push_back(kCyan); marker.push_back(27); catlabels.push_back("2e2#mu");
  cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_constrained&&m4l_constrained<135&&event_type==0)");
  cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_constrained&&m4l_constrained<135&&event_type==1)");
  cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_constrained&&m4l_constrained<135&&event_type==2)");
  cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_constrained&&m4l_constrained<135&&event_type==3)");
  //cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_fsr&&m4l_fsr<135&&event_type==0)");
  //cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_fsr&&m4l_fsr<135&&event_type==1)");
  //cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_fsr&&m4l_fsr<135&&event_type==2)");
  //cuts.push_back("(pass_vtx4lCut==1 && 110<m4l_fsr&&m4l_fsr<135&&event_type==3)");

  
  //NWA
  /*
  categories.push_back("ggF_4mu"); color.push_back(kViolet); marker.push_back(24); catlabels.push_back("4#mu ggF-like");
  categories.push_back("ggF_4e"); color.push_back(kOrange); marker.push_back(25); catlabels.push_back("4e ggF-like");
  categories.push_back("ggF_2mu2e"); color.push_back(kGreen+2); marker.push_back(26); catlabels.push_back("2#mu2e ggF-like");
  categories.push_back("VBF"); color.push_back(kRed); marker.push_back(32); catlabels.push_back("VBF-like");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<3000&&event_type==0 && !(dijet_invmass>400 && dijet_deltaeta>3.3))");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<3000&&event_type==1 && !(dijet_invmass>400 && dijet_deltaeta>3.3))");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<3000&&(event_type==3||event_type==2) && !(dijet_invmass>400 && dijet_deltaeta>3.3))");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<3000 && (dijet_invmass>400 && dijet_deltaeta>3.3))");
  */
  
  
  
  //LWA
  /*
  categories.push_back("ggF_4mu"); color.push_back(kViolet); marker.push_back(24); catlabels.push_back("4#mu ggF-like");
  categories.push_back("ggF_4e"); color.push_back(kOrange); marker.push_back(25); catlabels.push_back("4e ggF-like");
  categories.push_back("ggF_2mu2e"); color.push_back(kGreen+2); marker.push_back(26); catlabels.push_back("2#mu2e ggF-like");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&event_type==0)");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&event_type==1)");
  cuts.push_back("(pass_vtx4lCut==1 && 140<m4l_constrained_HM&&m4l_constrained_HM<1300&&(event_type==3||event_type==2))");
  */

  //MELA
/* 
  categories.push_back("ggF_4mu_MELA1"); color.push_back(kMagenta+4); marker.push_back(20); catlabels.push_back("4#mu MELA 0.0-0.2");
  categories.push_back("ggF_4mu_MELA2"); color.push_back(kMagenta+3); marker.push_back(21); catlabels.push_back("4#mu MELA 0.2-0.4");
  categories.push_back("ggF_4mu_MELA3"); color.push_back(kMagenta+2); marker.push_back(22); catlabels.push_back("4#mu MELA 0.4-0.6");
  categories.push_back("ggF_4mu_MELA4"); color.push_back(kMagenta+1); marker.push_back(23); catlabels.push_back("4#mu MELA 0.6-0.8");
  categories.push_back("ggF_4mu_MELA5"); color.push_back(kMagenta); marker.push_back(24); catlabels.push_back("4#mu MELA 0.8-1.0");
  categories.push_back("ggF_4e_MELA1"); color.push_back(kCyan+4); marker.push_back(25); catlabels.push_back("4e MELA 0.0-0.2");
  categories.push_back("ggF_4e_MELA2"); color.push_back(kCyan+3); marker.push_back(26); catlabels.push_back("4e MELA 0.2-0.4");
  categories.push_back("ggF_4e_MELA3"); color.push_back(kCyan+2); marker.push_back(27); catlabels.push_back("4e MELA 0.4-0.6");
  categories.push_back("ggF_4e_MELA4"); color.push_back(kCyan+1); marker.push_back(28); catlabels.push_back("4e MELA 0.6-0.8");
  categories.push_back("ggF_4e_MELA5"); color.push_back(kCyan); marker.push_back(29); catlabels.push_back("4e MELA 0.8-1.0");
  categories.push_back("ggF_2mu2e_MELA1"); color.push_back(kRed+4); marker.push_back(30); catlabels.push_back("2mu2e MELA 0.0-0.2");
  categories.push_back("ggF_2mu2e_MELA2"); color.push_back(kRed+3); marker.push_back(31); catlabels.push_back("2mu2e MELA 0.2-0.4");
  categories.push_back("ggF_2mu2e_MELA3"); color.push_back(kRed+2); marker.push_back(32); catlabels.push_back("2mu2e MELA 0.4-0.6");
  categories.push_back("ggF_2mu2e_MELA4"); color.push_back(kRed+1); marker.push_back(33); catlabels.push_back("2mu2e MELA 0.6-0.8");
  categories.push_back("ggF_2mu2e_MELA5"); color.push_back(kRed); marker.push_back(34); catlabels.push_back("2mu2e MELA 0.8-1.0");
  //const char* MELA = "(1/(1 + c_factor(m4l_constrained_HM,event_type)*exp(-1*KD_discriminant)))";
  const char* MELA = "MELA";
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.0 && %s<0.2)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.2 && %s<0.4)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.4 && %s<0.6)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.6 && %s<0.8)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==0) && %s>=0.8 && %s<=1.0)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.0 && %s<0.2)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.2 && %s<0.4)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.4 && %s<0.6)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.6 && %s<0.8)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&event_type==1) && %s>=0.8 && %s<=1.0)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.0 && %s<0.2)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.2 && %s<0.4)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.4 && %s<0.6)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.6 && %s<0.8)",MELA,MELA));
  cuts.push_back(Form("((pass_vtx4lCut==1 && 200<m4l_constrained_HM&&m4l_constrained_HM<1200&&(event_type==2||event_type==3)) && %s>=0.8 && %s<=1.0)",MELA,MELA));
  */
  
  
  TFile* file = new TFile(Form("plots/%s",outTFileName.c_str()),"RECREATE");
  std::ofstream outfile;
  std::ofstream outfileSI;
  std::ofstream outfileSys;
  outfile.open(outFileName.c_str(), std::ios::out);
  outfileSI.open(outFileSIName.c_str(), std::ios::out);
  outfileSys.open(outFileSysName.c_str(), std::ios::out);

  //
  // LOOP OVER PRODUCTION MODES
  //
  TCanvas can("can","",600,600);
  can.SetGrid(2,2);
  can.SetLeftMargin(0.15);

  for (unsigned int p(0);p<prod.size();++p){

    for (unsigned int w(0);w<widths.size();++w){
      std::cout<<"\nprod: "<<prod[p]<<" width: "<<widths[w]<<std::endl;

      bool hasTau=true;
      if (widths[w]=="NW") hasTau= true;
      if (widths[w]=="ZZ4lep") hasTau=false;

      //
      // BUILD MAP OF (MASS POINT VS FILE)
      //
      std::map< float , std::string > fileMap;
      std::string str;
      const char* entry;
      void* dirp = gSystem->OpenDirectory(minitreeDir.c_str());


      while ((entry=(char*)gSystem->GetDirEntry(dirp))) {
        str = entry;
        size_t pos = str.find(prod[p].c_str());
        size_t pos2 = str.find(widths[w].c_str());

        if (pos!=std::string::npos && pos2!=std::string::npos){
          pos+=prod[p].size()-1;
          std::string massString="";
          while (!isdigit(str[++pos])){}
          while (isdigit(str[pos])){ massString+=str[pos++]; }
          if (str.substr(pos,2)=="p5") massString+=".5";
          if (str.substr(pos,3)=="p25") massString+=".25";
          float mass = atof(massString.c_str());
          if (mass<minMH || mass>maxMH) continue;
          std::cout<<"from file:"<<str<<" extracted mass="<<mass<<std::endl;

          if (hasTau&&str.find("noTau")!=std::string::npos) continue;
          if (!hasTau && str.find("noTau")==std::string::npos) continue;

          fileMap[mass] = minitreeDir+str;
        }
        //else
          //std::cout<<"on file "<<str<<" found pos="<<pos<<" and pos2="<<pos2<<std::endl;
      }

      std::cout<<"built a map of file names vs mass. size="<<fileMap.size()<<std::endl;
      if (fileMap.size()==0) continue;

      //BUILD A MAP OF (MASS POINT VS ACC)
      std::map<std::string, std::vector<float> > norm, norm_error;
      std::vector<float> zero, masses;

      for (std::map<float,std::string>::iterator it=fileMap.begin(); it!= fileMap.end(); ++it){

        std::cout<<"for mass="<<it->first<<", using file "<<it->second<<std::endl;

        TFile* file = new TFile(it->second.c_str(),"READ");
        TTree* tree = (TTree*)file->Get("tree_incl_all");

        masses.push_back((*it).first);
        zero.push_back(0.);

        for (unsigned int c(0);c<categories.size();++c){
          std::string cat = categories[c];
          TH1F* h = new TH1F("h","h",1,-9999,9999);
          tree->Draw("m4l_constrained_HM>>h",Form("%s*(9./4.)*(weight/(w_lumi*w_xs*w_br))",cuts[c].c_str()));
          norm[cat].push_back(h->GetBinContent(1));
          norm_error[cat].push_back(h->GetBinError(1));
          delete h;
        }

      }

      std::cout<<"built a map of acceptance vs mass. size="<<masses.size()<<std::endl;

      std::cout<<"MC acceptances for "<<prod[p]<<":"<<std::endl;
      for (unsigned int c(0);c<categories.size();++c){
        std::cout<<"\t"<<categories[c]<<std::endl;
        for (unsigned int m(0);m<masses.size();++m){
          std::cout<<"\t\t"<<masses[m]<<":\t"<<norm[categories[c]][m]<<"\t+-\t"<<norm_error[categories[c]][m]<<std::endl;
        }
      }

      const char* pol="pol0";
      if (masses.size()>1) pol="pol1";
      if (masses.size()>2) pol="pol2";
      //if (masses.size()>3) pol="pol3";//FIXME
      //if (masses.size()>4) pol="pol4";
      //if (masses.size()>5) pol="pol5";
      //if (masses.size()>5) pol="pol6";

      //BUILD A GRAPH OF MASS POINT VS NORMS
      TMultiGraph* mg = new TMultiGraph();
      can.cd();
      for (unsigned int c(0);c<categories.size();++c){
        std::string cat = categories[c];

        TGraphErrors* graph = new TGraphErrors(masses.size(), &(masses[0]), &(norm[cat][0]), &(zero[0]), &(norm_error[cat][0]));
        TFitResultPtr fit = graph->Fit(pol,"qs");
        mg->Add(graph);
        file->cd();
        std::string graphName = Form("AC_%s_%s_%s",prod[p].c_str(),widths[w].c_str(),cat.c_str());
        //std::cout<<"going to save graph as "<<graphName<<std::endl;
        graph->Clone(graphName.c_str())->Write();


        //ESTIMATE SYSTEMATIC
        TF1* polfit = graph->GetFunction(pol);
        std::vector<float> dev;
        for (unsigned int m(0);m<masses.size();++m){
          dev.push_back(polfit->Eval(masses[m]) - norm[cat][m]);
          //std::cout<<"actual dev "<<masses[m]<<"\t = "<<norm[cat][m] / polfit->Eval(masses[m])<<std::endl;
        }
        float rms = TMath::RMS(dev.size(),&(dev[0]))/sqrt(masses.size()-1);
        dev.clear();
        for (unsigned int m(0);m<masses.size();++m){
          float nom = polfit->Eval(masses[m]);
          outfileSys<<prod[p]<<"  "<<cat<<"  "<<masses[m]<<"  "<<Form("%.8f",(nom-rms)/nom)<<"   "<<Form("%.8f",(nom+rms)/nom)<<std::endl;
        }
        
        //Printout for conf note table
        std::cout<<"prod "<<prod[p]<<" in category "<<cat<<" acceptance @125GeV = "<<100*polfit->Eval(125)<<"%,    @126GeV = "<<100*polfit->Eval(126)<<"%"<<std::endl;

        can.cd();
        graph->SetMinimum(0);
        graph->GetFunction(pol)->SetLineColor(color[c]);
        //graph->SetTitle(Form("%s %s",prod[p].c_str(), cat.c_str()));
        graph->SetTitle("");
        graph->GetXaxis()->SetTitle("m_{H} [GeV]");
        graph->GetXaxis()->SetRangeUser(masses[0],masses[masses.size()-1]);
        graph->GetXaxis()->CenterTitle();
        graph->GetXaxis()->SetTitleFont(42);
        graph->GetXaxis()->SetTitleSize(0.04);
        graph->GetXaxis()->SetLabelSize(0.035);
        graph->GetXaxis()->SetRangeUser(minPlot,maxPlot);
        graph->GetYaxis()->SetTitle(Form("%s  acceptance",prodlabel[p].c_str()));
        graph->GetYaxis()->SetTitleFont(42);
        graph->GetYaxis()->SetTitleOffset(1.4);
        graph->GetYaxis()->SetTitleSize(0.04);
        graph->GetYaxis()->SetLabelSize(0.035);
        if (minPlotY!=-1) graph->SetMinimum(minPlotY);
        if (maxPlotY!=-1) graph->SetMaximum(maxPlotY);
        graph->SetMarkerStyle(marker[c]);
        graph->SetMarkerSize(0.75);
        graph->SetMarkerColor(color[c]);
        graph->SetLineColor(color[c]);
        graph->Draw("ap");

        TLine* line = new TLine(0.65,0.75,0.68,0.75);
        line->SetNDC(true);
        line->SetLineColor(color[c]);
        line->SetLineWidth(3);
        line->Draw("same");
        myText(0.57,0.74,kBlack,catlabels[c].c_str(),0.04);
        //ATLASLabel(0.18,0.84,"Internal",kBlack);

        gSystem->Exec(Form("mkdir -p $PWD/plots/acceptance/%s/%s/",prod[p].c_str(),widths[w].c_str()));
        //can.Print(Form("$PWD/plots/acceptance/%s/%s/acceptance_%s_%s_%s_%s.png",prod[p].c_str(),widths[w].c_str(),prod[p].c_str(),widths[w].c_str(),cat.c_str(),pol)); 
        can.Print(Form("$PWD/plots/acceptance/%s/%s/acceptance_%s_%s_%s_%s.eps",prod[p].c_str(),widths[w].c_str(),prod[p].c_str(),widths[w].c_str(),cat.c_str(),pol)); 

        //STORE POLY PARAMS
        std::vector<double> polyParameters = fit->Parameters();
        outfile << Form("[%s %s %s]\n",prod[p].c_str(),widths[w].c_str(),cat.c_str());
        outfileSI << Form("[%s %s %s]\n",prod[p].c_str(),widths[w].c_str(),cat.c_str());
        std::string outline = "";
        std::string outlineSI = "& ";
        for (unsigned int r(0);r<polyParameters.size();++r) {
          outline=outline+" "+Form("%.20f",polyParameters[r]);
          outlineSI=outlineSI+Form("%.3e",polyParameters[r])+" & ";
        }
        outlineSI+="\n";
        outline+="\n";
        outfile << outline;
        outfileSI<<outlineSI;
      }

      mg->Draw("ap");
      if (minPlotY!=-1) mg->SetMinimum(minPlotY);
      if (maxPlotY!=-1) mg->SetMaximum(maxPlotY);
      mg->GetXaxis()->SetRangeUser(masses[0],masses[masses.size()-1]);
      mg->GetXaxis()->SetTitle("m_{H} [GeV]"); //m_S for high mass //FIXME
      mg->GetXaxis()->CenterTitle();
      mg->GetXaxis()->SetTitleFont(42);
      mg->GetXaxis()->SetTitleSize(0.04);
      mg->GetXaxis()->SetLabelSize(0.035);
      mg->GetXaxis()->SetRangeUser(minPlot,maxPlot);
      mg->GetYaxis()->SetTitle(Form("%s  acceptance",prodlabel[p].c_str()));
      mg->GetYaxis()->SetTitleFont(42);
      mg->GetYaxis()->SetTitleOffset(1.4);
      mg->GetYaxis()->SetTitleSize(0.04);
      mg->GetYaxis()->SetLabelSize(0.035);

      float x(0.25),y(0.85);
      ATLASLabel(x,y,"Internal",kBlack);
      for (unsigned int c(0);c<categories.size();++c)
          myLineText(x,y-0.05-0.04*c,color[c],kSolid,2,catlabels[c].c_str(),0.035);

      can.Update();
      can.Print(Form("$PWD/plots/acceptance/%s/%s/acceptance_%s_%s_%s.eps",prod[p].c_str(),widths[w].c_str(),prod[p].c_str(),widths[w].c_str(),pol)); 

    } //end loop over widths

  } //end loop over prods

  outfile.close();
  outfileSI.close();
  outfileSys.close();
}

