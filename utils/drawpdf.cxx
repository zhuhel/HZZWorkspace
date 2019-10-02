#include <stdlib.h>
#include <iostream>
#include <string>
#include <utility>
#include <TStyle.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <RooWorkspace.h>
#include <RooRealVar.h>
#include <TFile.h>
#include <TTree.h>
#include <RooPlot.h>
#include <RooSimultaneous.h>

#include "HZZWorkspace/RooStatsHelper.h"
#include "HZZWorkspace/Helper.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;

//------------------------------------------------------------------------------
// Draws PDFs in each category within a workspace on top of data
// --> Prefit data-MC comparison
//
// Takes 6 required arguments:
// - workspace ROOT file name
// - ModelConfig name (usually "ModelConfig")
// - Dataname ("asimovData" "or "obsData")
// - wsName (default "combined" )
// - Plotting options: "log" - plot in logarithm scale 
// *** Currently out of  service! Check back soon OR better yet! Get involved! ***
//------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    if (argc < 7){
        log_err(" Please provide the following 6 args: fileName - mcName - dataName - obsName - wsName - option ");
        return EXIT_FAILURE;

    }
    string fileName(argv[1]);
    string mcName(argv[2]);
    string dataName(argv[3]);
    string obsName(argv[4]);
    string wsName(argv[5]);
    string option(argv[6]);
    vector<string> obsVec;
    obsVec.clear();
    Helper::tokenizeString(obsName,',',obsVec);

    TFile* f = new TFile(fileName.data());
    RooWorkspace* w = (RooWorkspace*)(f->Get(wsName.data()));
    RooStats::ModelConfig* mc = (RooStats::ModelConfig*)(w->obj(mcName.data()));

    //RooDataSet* data = AsymptoticCalculator::GenerateAsimovData( *mc->GetPdf(), *mc->GetObservables() );
    RooAbsData* data = w->data(dataName.data());

    RooSimultaneous* pdf = (RooSimultaneous*)mc->GetPdf();

    RooAbsCategoryLValue*  cat = (RooAbsCategoryLValue*)&pdf->indexCat();
    int nCat = cat->numBins(0);
    TList* dataList = data->split( *cat, true );

    for ( int i= 0; i < nCat; i++ ) {

        cat->setBin(i);
        RooAbsPdf* pdfi = pdf->getPdf(cat->getLabel());
        RooDataSet* datai = (RooDataSet*)( dataList->At(i) );

        RooPlot* xframe = w->var(obsVec.at(i).data())->frame();

        datai->plotOn(xframe, Name("Data"));

        pdfi->plotOn(xframe, LineColor(kBlack),Name("nominal"));
        RooArgSet* check = pdfi->getComponents();
        check->Print();
        TIterator* iter = check->createIterator();
        TObject* arg = NULL;
        vector<string> nComp;
        nComp.clear();
        int icolor = 1;
        while((arg=(TObject*)iter->Next())){
           string str(arg->GetName());
           if((str.find("ATLAS") != std::string::npos) && (str.find("fiv") == std::string::npos && str.find("alpha") == std::string::npos && str.find("nTot") == std::string::npos)){
           icolor++;
           cout <<arg->GetName()<<arg->ClassName()<<endl;
           nComp.push_back(str);
           pdfi->plotOn(xframe,Components(str.data()),LineColor(icolor),Name(str.data()));
         }
        }

        TCanvas *c1 = new TCanvas("c1","c1",800,600);
        if(option.find("log") != std::string::npos) c1->SetLogy();
        string saveName = "Pdf_"+string(pdfi->GetName());
        xframe->SetMaximum(xframe->GetMaximum()*2);
        xframe->GetXaxis()->SetTitle(obsVec.at(i).data());
        xframe->GetXaxis()->SetTitleSize(0.06);
        xframe->GetXaxis()->SetLabelSize(0.06);
        xframe->GetXaxis()->SetTitleOffset(1.2);
        xframe->GetYaxis()->SetTitle("Events");
        xframe->GetYaxis()->SetTitleSize(0.06);
        xframe->GetYaxis()->SetLabelSize(0.06);
        xframe->GetYaxis()->SetTitleOffset(1.2);

        TLegend *lg = new TLegend(0.20, 0.60, 0.75, 0.85);
        lg->SetTextSize(0.033);
        lg->SetBorderSize(0);
        lg->SetFillColor(10);
        lg->AddEntry(xframe->findObject("Data"),"Input Data","PE");
        lg->AddEntry(xframe->findObject("nominal"),"sigal+bkg nominal pdf", "L");
        for(size_t icomp=0; icomp<nComp.size(); icomp++)
          lg->AddEntry(xframe->findObject(nComp.at(icomp).data()),nComp.at(icomp).data(), "L");
        xframe->Draw();
        lg->Draw();

        c1->SaveAs(Form("%s.png",saveName.data()));
        c1->SaveAs(Form("%s.pdf",saveName.data()));

        delete c1;
 }
}
