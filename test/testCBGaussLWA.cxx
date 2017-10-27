#include "Hzzws/CBGauss.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "RooAbsPdf.h"
#include "TFile.h"
#include "RooDataSet.h"
#include "RooWorkspace.h"

EColor cols[]={kGreen,kSpring,kYellow,kOrange,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,
               kGreen,kSpring,kYellow,kOrange,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,
               kGreen,kSpring,kYellow,kOrange,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,
               kGreen,kSpring,kYellow,kOrange,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal};

using namespace RooFit;

int main(){

  strmap args;
  args["global"]="N(10)";

  auto c = new Coefficient(args);

  Helper::getInputPath("/afs/cern.ch/atlas/groups/HSG2/H4l/run2/2016/Workspaces/HighMass/Prod_v04/2016_07_11/input04_LWA/");
  auto p = new CBGauss("pdf", "analyticalParams", "", true, false);
  p->addCoefficient(c);


  RooRealVar x("m4l_constrained","m4l_constrained",140,1500);
  RooRealVar event_type("event_type","event_type",0,3);
  RooRealVar weight("weight","weight",0,999);
  RooArgSet obs(x);


  RooWorkspace* wsp = new RooWorkspace("wsp","");


  const char* cnames[3]={"ggF_2mu2e_13TeV","ggF_4mu_13TeV","ggF_4e_13TeV"};

  for (float gamma_rel_pct(0.01);gamma_rel_pct<10.1;gamma_rel_pct*=10.){
    std::map<std::string, RooPlot*> frame;
    frame["ggF_4mu_13TeV"] = x.frame(Title("4#mu"),Range(140,1500));
    frame["ggF_4e_13TeV"] = x.frame(Title("4e"),Range(140,1500));
    frame["ggF_2mu2e_13TeV"] = x.frame(Title("2#mu2e"),Range(140,1500));

    for (int c(0);c<3;++c){
      p->setChannel(obs, cnames[c], false);
      wsp->import(*p->getPDF(),RecycleConflictNodes());
      int k(0);
      for (float mH(350);mH<1001;mH+=50.){
        ++k;
        float gamma = 0.01*gamma_rel_pct*mH;
        wsp->var("mH")->setVal(mH);
        wsp->var("gamma")->setVal(gamma);

        auto pdf = wsp->pdf(Form("pdf_%s_conv",cnames[c]));
        std::cout<<"tried to get pdf with name: "<<Form("pdf_%s_conv",cnames[c])<<std::endl;
        std::cout<<"pdf:"<<pdf; pdf->Print();
        pdf->plotOn(frame[cnames[c]],LineColor(cols[k]));
        frame[cnames[c]]->SetMinimum(1e-6);
      }
    }
    TCanvas* can = new TCanvas("c","c",900,1200);
    can->Divide(1,3);
    can->cd(1)->SetLogy(); frame["ggF_2mu2e_13TeV"]->Draw();
    can->cd(2)->SetLogy(); frame["ggF_4mu_13TeV"]->Draw();
    can->cd(3)->SetLogy(); frame["ggF_4e_13TeV"]->Draw();
    can->Print(Form("plots/CBGaussLWA_Gamma%2.2fpct.png",gamma_rel_pct));
  }


}
