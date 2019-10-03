#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include <memory>

#include "RooFitResult.h"
#include "RooMinimizer.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooCmdArg.h"
#include "RooRealVar.h"
#include "RooCategory.h"
#include "RooRandom.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/ModelConfig.h"
#include "TH2.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "HZZWorkspace/RooStatsHelper.h"
#include "HZZWorkspace/Helper.h"

//-----------------------------------------------------------------------------
// Statistics helper class to
// * facilitate minimization
// * create NLL
// * create Asimov data
// * estimate significance
// * generate toys
// * scan POI
// * get observed number of events
//-----------------------------------------------------------------------------

void RooStatsHelper::setDefaultMinimize(){
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(1);
  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);
}


RooFitResult* RooStatsHelper::minimize(RooNLLVar* nll, RooWorkspace* ){
    log_warn("Careful! You called minimize (RooNLLVar*, RooWorkspace* ) - be aware that this signature is deprecated in favour of minimize(RooNLLVar* nll, \
        bool save, const RooArgSet* minosSet), since the WS is not actually used anywhere. Please update your method call to make this warning go away.");
    return minimize(nll);
}

RooFitResult* RooStatsHelper::minimize(RooNLLVar* nll,
        bool save, const RooArgSet* minosSet)
{
    nll->enableOffsetting(true);
    int printLevel = ROOT::Math::MinimizerOptions::DefaultPrintLevel();
    // RooFit::MsgLevel msglevel = RooMsgService::instance().globalKillBelow();

    if (printLevel < 0) RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);

    int strat = ROOT::Math::MinimizerOptions::DefaultStrategy();

    RooMinimizer minim(*nll);
    minim.optimizeConst(2);
    minim.setStrategy(strat);
    minim.setPrintLevel(printLevel);
    // minim.setProfile();  // print running time
    // minim.setEps(0.001);

    // minim.setErrorLevel(1E-3);


    int status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());


    if (status != 0 && status != 1 && strat < 2)
    {
        strat++;
        cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
        minim.setStrategy(strat);
        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }

    if (status != 0 && status != 1 && strat < 2)
    {
        strat++;
        cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
        minim.setStrategy(strat);
        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
    }
    // //switch minuit version and try again
    if (status != 0 && status != 1)
    {
        string minType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
        string newMinType;
        if (minType == "Minuit2") newMinType = "Minuit";
        else newMinType = "Minuit2";

        cout << "Switching minuit type from " << minType << " to " << newMinType << endl;

        ROOT::Math::MinimizerOptions::SetDefaultMinimizer(newMinType.c_str());
        strat = ROOT::Math::MinimizerOptions::DefaultStrategy();
        minim.setStrategy(strat);

        status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());

        if (status != 0 && status != 1 && strat < 2)
        {
            strat++;
            cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        }

        if (status != 0 && status != 1 && strat < 2)
        {
            strat++;
            cout << "Fit failed with status " << status << ". Retrying with strategy " << strat << endl;
            minim.setStrategy(strat);
            status = minim.minimize(ROOT::Math::MinimizerOptions::DefaultMinimizerType().c_str(), ROOT::Math::MinimizerOptions::DefaultMinimizerAlgo().c_str());
        }

        ROOT::Math::MinimizerOptions::SetDefaultMinimizer(minType.c_str());
    }
    if (minosSet != NULL) {
        minim.minos(*minosSet);
    }
    // minim.minos();
    if (save && status == 0) return minim.save();
    else return NULL;
}

void RooStatsHelper::setVarfixed(RooWorkspace* combined, const char* varName, double value)
{
    RooRealVar* _var = (RooRealVar*) combined->var(varName);
    if (_var) {
        cout<<varName<<" fixed to " << value<<endl;
        _var->setVal(value);
        _var->setConstant(1);
    }
    else {
        cout << "Error: cannot find " << varName <<endl;
    }
}

void RooStatsHelper::setVarFree(RooWorkspace* combined, const char* varName)
{
    RooRealVar* _var = (RooRealVar*) combined->var(varName);
    if (_var) {
        _var->setConstant(kFALSE);
    }
    else {
        cout << "Error: RooStatsHelper cannot find " << varName << endl;
    }
}

pair<double,double> RooStatsHelper::getVarVal(const RooWorkspace& combined, const char* varName)
{
    if(varName == NULL) return make_pair(-9999., -9999.);
    RooRealVar* mhiggs = dynamic_cast<RooRealVar*>(combined.var(varName));
    if (mhiggs) {
        return make_pair(mhiggs ->getVal(),mhiggs->getError());
    } else {
        return make_pair(-9999.,-9999.);
    }
}

RooNLLVar* RooStatsHelper::createNLL(RooAbsData* _data, RooStats::ModelConfig* _mc)
{
    const RooArgSet& nuis = *_mc->GetNuisanceParameters();
    RooNLLVar* nll = (RooNLLVar*)_mc->GetPdf()->createNLL(*_data,
            RooFit::Constrain(nuis),
            RooFit::GlobalObservables(*_mc->GetGlobalObservables())
            );
    // RooCmdArg agg = condVarSet.getSize() > 0?RooFit::ConditionalObservables(condVarSet):RooCmdArg::none(); // for conditional RooFit
    // RooCmdArg agg = RooCmdArg::none();

    return nll;
}

double RooStatsHelper::getPvalue(RooWorkspace* combined,
        RooStats::ModelConfig* mc,
        RooAbsData* data,
        const char* muName,
        bool isRatioLogLikelihood)
{
    if(!combined || !mc || !data){
        log_err("check inputs");
        return -9999;
    }
    TString dataname(data->GetName());

    cout<<"Getting pvalue for "<< data->GetName()<<endl;
    auto combPdf = mc->GetPdf();
    if(!combPdf) {
        log_err("overall pdf does not exist!");
        return -9999;
    } else {
        combPdf->Print();
    }

    if (data->numEntries() <= 0) {
        log_err("total number of events is less than 0: %d",
                data->numEntries());
        return 1.0;
    }
    else {
        cout<<"total events: "<< data->numEntries() <<" sumEntries = "<< data ->sumEntries()<<endl;
    }

    RooRealVar* mu = (RooRealVar*) combined ->var(muName);
    if (!mu){
        log_err("%s does not exist", mu->GetName());
        return -9999;
    }

    if(!isRatioLogLikelihood){
        mu ->setConstant(0);
        mu->Print();
        cout<<"Fitting with "<<mu->GetName()<<" free"<<endl;
    }else{
        mu ->setVal(1);
        mu ->setConstant(1);
        cout<<"Fitting with "<<mu->GetName()<<"=1"<<endl;
    }


    PrintExpEvts(combPdf, mu, mc->GetObservables(), data);
    RooNLLVar* nll = createNLL(data, mc);

    minimize(nll);
    double obs_nll_min = nll ->getVal();
    cout << "mu_hat for " << data->GetName() << " " << mu->getVal() <<
        " " << mu->getError() << " " << obs_nll_min << endl;
    delete nll;
    PrintExpEvts(combPdf, mu, mc->GetObservables(), data);
    bool reverse = (mu ->getVal() < 0);

    // combined->saveSnapshot("unconditionalFit", combined->allVars());
    // combined->writeToFile("UnCond_XS_ggF_combined_HZZ_1200GeV_llqq_vvqq_afterPara.root");

    cout<<"Fitting background only hypothesis "<< mu->GetName()<<endl;

    mu ->setVal(1.0e-200);
    mu ->setConstant(1);
    RooNLLVar* nllCond = createNLL(data, mc);
    minimize(nllCond);
    double obs_nll_min_bkg = nllCond ->getVal();
    cout<<"NLL for background only: "<< obs_nll_min_bkg<<endl;
    delete nllCond;

    double obs_q0 = 2*(obs_nll_min_bkg - obs_nll_min);
    if(reverse) obs_q0 = -obs_q0;
    double sign = int(obs_q0 == 0 ? 0 : obs_q0 / fabs(obs_q0));
    double obs_sig = sign*sqrt(fabs(obs_q0));
    cout<<"NLL min: "<< obs_nll_min <<" "<< obs_nll_min_bkg <<" "<< obs_sig<<endl;

    return RooStats::SignificanceToPValue(obs_sig);
}

double RooStatsHelper::calculateSignificance(double obs_nll_min, double obs_nll_bkg){
    double obs_q0 = 2*(obs_nll_bkg - obs_nll_min);
    double sign = int(obs_q0 == 0 ? 0 : obs_q0 / fabs(obs_q0));
    double obs_sig = sign*sqrt(fabs(obs_q0));
    return obs_sig;
}

RooDataSet* RooStatsHelper::makeAsimovData(RooWorkspace* combined,
        double muval,
        double profileMu,
        const char* muName,
        const char* mcname,
        const char* dataname,
        bool doprofile)
{
    RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);
    cout<<"In making asimov-Data"<<endl;
    stringstream muStr("_");
    muStr << muval;
    if(muval==1){
        if(profileMu == 0) muStr<<"_paz";
        if(profileMu > 1) muStr <<"_pamh";
    }

    RooStats::ModelConfig* mcInWs = (RooStats::ModelConfig*) combined->obj(mcname);
    RooRealVar* mu = (RooRealVar*) combined->var(muName);
    mu->setVal(muval);

    const RooArgSet& mc_globs = *mcInWs->GetGlobalObservables();
    const RooArgSet& mc_nuis  = *mcInWs->GetNuisanceParameters();

    RooDataSet* combData = NULL;
    try {
        combData = (RooDataSet*) combined->data(dataname);
    } catch (...) {}

    //save the snapshots of nominal parameters, but only if they're not already saved
    if (!combined->loadSnapshot("nominalGlobs"))
    {
        cout << "nominalGlobs doesn't exist. Saving snapshot." << endl;
        // set global observable to 0 first
        TIter glob_iter (mc_globs.createIterator());
        RooRealVar* glob;
        while ( (glob = (RooRealVar*) glob_iter()) ){
            if(TString(glob->GetName()).Contains("gamma_stat")){
                continue;
            } else{
                glob->setVal(0.);
            }
        }
        combined->saveSnapshot("nominalGlobs",*mcInWs->GetGlobalObservables());
    }
    if (!combined->loadSnapshot("nominalNuis"))
    {
        cout << "nominalNuis doesn't exist. Saving snapshot." << endl;
        TIter nuis_iter (mc_nuis.createIterator());
        RooRealVar* nuis;
        while ( (nuis= (RooRealVar*) nuis_iter()) ){
            TString nuis_ts(nuis->GetName());
            if( nuis_ts.Contains("gamma_stat") ){
                continue;
            } else if( nuis_ts.Contains("_norm_") ) {
                nuis->setVal(1.);
            }else{
                nuis->setVal(0.);
            }
        }
        combined->saveSnapshot("nominalNuis",*mcInWs->GetNuisanceParameters());
    }

    mu ->setVal(profileMu);
    if (profileMu > 1)
    { // would profile to the \hat_mu,
        mu ->setConstant(kFALSE);
        // mu ->setRange(-40,40); // !!!!Do that outside the function!!!
        RooNLLVar* conditioning_nll = createNLL(combData, mcInWs);
        minimize(conditioning_nll);//find the \hat_mu
        delete conditioning_nll;
    }
    mu ->setConstant(kTRUE); // Fix mu at the profileMu
    // mu ->setRange(0,40); // !!!!Do that outside the function!!!

    if (doprofile)
    { // profile nuisance parameters at mu = profileMu
        RooNLLVar* conditioning_nll = createNLL(combData, mcInWs);
        minimize(conditioning_nll);
    }

    // loop over the nui/glob list, grab the corresponding variable from the tmp ws,
    // and set the glob to the value of the nui
    TIter glob_iter (mc_globs.createIterator());
    TIter nuis_iter (mc_nuis.createIterator());
    RooRealVar* nuis;
    RooRealVar* glob;
    while ((
                glob = (RooRealVar*) glob_iter(),
                nuis = (RooRealVar*) nuis_iter()
           ))
    {
        if (!nuis || !glob) continue;
        // cout<<" nuisance Name: " << nuis->GetName() << endl;
        // cout<<" global Name: " << glob->GetName() << endl;
        glob->setVal(nuis->getVal());
    }
    combined->saveSnapshot(("conditionalGlobs"+muStr.str()).c_str(), mc_globs);
    combined->saveSnapshot(("conditionalNuis" +muStr.str()).c_str(), mc_nuis);


    mu->setVal(muval);
    cout<<"Making asimov data, "<< muval <<endl;

    auto adata =  makeUnconditionalAsimov(combined, mcInWs, Form("asimovData_%s",muStr.str().c_str()));

    combined->loadSnapshot("nominalGlobs");

    return adata;
}

RooDataSet* RooStatsHelper::makeUnconditionalAsimov(RooWorkspace* combined, RooStats::ModelConfig* mcInWs, const char * dataname)
{


  const char* weightName="weightVar";
  RooArgSet obsAndWeight;
  obsAndWeight.add(*mcInWs->GetObservables());
  RooRealVar* weightVar = NULL;
  if (!(weightVar = combined->var(weightName)))
  {
    combined->import(*(new RooRealVar(weightName, weightName, 1,0,100000000)));
    weightVar = combined->var(weightName);
  }
  obsAndWeight.add(*combined->var(weightName));

  RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mcInWs->GetPdf());
  if(!simPdf) {
    log_err("cannot cast %s to RooSimultaneous", mcInWs->GetPdf()->GetName());
    return NULL;
  }
  else {
    log_info("using %s to generate asimov:", simPdf->GetName());
    simPdf->Print();
  }
  map<string, RooDataSet*> asimovDataMap;

  RooCategory* channelCat = (RooCategory*)&simPdf->indexCat();
  TIterator* iter = channelCat->typeIterator() ;
  RooCatType* tt = NULL;
  int nrIndices = 0;
  while((tt=(RooCatType*) iter->Next())) {
    nrIndices++;
  }

  int iFrame=0;
  for (int i=0;i<nrIndices;i++){
    channelCat->setIndex(i);
    iFrame++;
    // Get pdf associated with state from simpdf
    RooAbsPdf* pdftmp = simPdf->getPdf(channelCat->getLabel()) ;
    if(!pdftmp){
      cout<<"ERRORs no pdf associated with "<< channelCat ->getLabel()<<endl;
    }
    // Generate observables defined by the pdf associated with this state
    RooArgSet* obstmp = pdftmp->getObservables(*mcInWs->GetObservables()) ;

    RooDataSet* obsDataUnbinned = new RooDataSet(
        Form("combAsimovData%d",iFrame),
        Form("combAsimovData%d",iFrame),
        RooArgSet(obsAndWeight,*channelCat),
        RooFit::WeightVar(*weightVar)
        );
    double expectedEvents = pdftmp->expectedEvents(*obstmp);

    int dim = obstmp->getSize();
    if (dim == 0) { //counting
      std::cout<<"building asimov for counting"<<std::endl;
      obsDataUnbinned->add(*obstmp, expectedEvents);
    }
    else if (dim <= 3) { //1D,2D,3D pdfs

      std::cout<<"building asimov for "<<dim<<"D data"<<std::endl;

      RooArgList obs(*obstmp);
      RooRealVar* X = (RooRealVar*)obs.at(0);
      RooRealVar* Y = (RooRealVar*)obs.at(1);
      RooRealVar* Z = (RooRealVar*)obs.at(2);

      RooCmdArg argY = ( Y ? RooFit::YVar( *Y ) : RooCmdArg::none() );
      RooCmdArg argZ = ( Z ? RooFit::ZVar( *Z ) : RooCmdArg::none() );

      // TH1* hist = pdftmp->createHistogram("htemp", *X, RooFit::IntrinsicBinning(false), argY, argZ);
      unique_ptr<TH1> hist( pdftmp->createHistogram("htemp", *X, RooFit::IntrinsicBinning(false), argY, argZ));

      hist->Scale(expectedEvents/hist->Integral()); //scale histogram to expectation
      for (int ix(1);ix<=hist->GetNbinsX(); ++ix){
        for (int iy(1);iy<=hist->GetNbinsY(); ++iy){
          for (int iz(1);iz<=hist->GetNbinsZ(); ++iz){
            if (X) X->setVal(hist->GetXaxis()->GetBinCenter(ix));
            if (Y) Y->setVal(hist->GetYaxis()->GetBinCenter(iy));
            if (Z) Z->setVal(hist->GetZaxis()->GetBinCenter(iz));
            double cellWidth = (X ? hist->GetXaxis()->GetBinWidth(ix) : 1) * 
                               (Y ? hist->GetYaxis()->GetBinWidth(iy) : 1) * 
                               (Z ? hist->GetZaxis()->GetBinWidth(iz) : 1);

            obsDataUnbinned->add(*obstmp, hist->GetBinContent(ix,iy,iz) * cellWidth);
          }
        }
      }
    }
    else { //4D or more :|
      std::cout<<"trying to build asimov for "<<obstmp->getSize()<<"D data: not supported!"<<std::endl;
      exit(1);
    }

    if(TMath::IsNaN(obsDataUnbinned->sumEntries())){
      cout << "sum entries is nan"<<endl;
      exit(1);
    }

    std::cout<<"in category "<<channelCat->getLabel()<<" using asimov dataset:"<<std::endl;
    obsDataUnbinned->Print();

    asimovDataMap[string(channelCat->getLabel())] = obsDataUnbinned;//tempData;
  }

  RooDataSet* asimovData = new RooDataSet(
      dataname,
      dataname,
      RooArgSet(obsAndWeight,*channelCat),
      RooFit::Index(*channelCat),
      RooFit::Import(asimovDataMap),
      RooFit::WeightVar(*weightVar)
      );
  combined->import(*asimovData);
  cout<<"AsimovData is created: "<< asimovData->GetName() << endl;
  //asimovData ->Print("v");
  //cout<< asimovData ->sumEntries()<<endl;

  return asimovData;
}


double RooStatsHelper::getRoughSig(double s, double b)
{
    return TMath::Sqrt(2*((s+b)*TMath::Log(1+s/b)-s));
}

void RooStatsHelper::randomizeSet(RooAbsPdf* pdf, RooArgSet* globs, int seed)
{
    RooRandom::randomGenerator() -> SetSeed(seed) ; // This step is necessary
    RooDataSet *one= pdf->generate(*globs, 1);
    const RooArgSet *values= one->get(0);
    RooArgSet *allVars=pdf->getVariables();
    *allVars=*values;
    delete one;
    delete allVars;
}

void RooStatsHelper::SetRooArgSetConst(RooArgSet& argset, bool flag)
{
    argset.setAttribAll("Constant",flag);
}
void RooStatsHelper::fitData(RooWorkspace* w, const char* mcName,
        const char* dataName,
        const char* muName, double poival, map<string,double>& result)
{
    RooAbsData* data = w->data(dataName);
    fitData(w, mcName, data, muName, poival, result);
}

void RooStatsHelper::fitData(RooWorkspace* w, const char* mcName,
        RooAbsData* data,
        const char* muName, double poival, map<string,double>& result)
{

  ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(1);

  cout<<"in "<<__func__<<std::endl;

  RooStats::ModelConfig *mc=(RooStats::ModelConfig*)w->obj(mcName);
  RooArgSet* nuisanceParameters = (RooArgSet*) mc->GetNuisanceParameters();
  // auto pois = const_cast<RooArgSet*>(mc->GetParametersOfInterest());

  //perform S+B fit (with mu=poival fixed)
  cout<<"S+B fixed"<<endl;
  RooRealVar* muVar =(RooRealVar*) w ->var(muName);
  muVar->setVal(poival);
  muVar->setConstant(true);
  RooNLLVar* nll_SBfixed = createNLL(data, mc);
  int status = (minimize(nll_SBfixed) == NULL)?0:1;

  result["nll_SB_fixed"] = nll_SBfixed->getVal();
  result["poi_SB_fixed"] = muVar->getVal();
  result["status_SB_fixed"] = status*1.0;
  delete nll_SBfixed;


  //perform S+B fit
  cout<<"S+B free"<<endl;
  muVar->setVal(0.);
  muVar->setConstant(false);

  RooNLLVar* nll_SBfree = createNLL(data, mc);
  status = (minimize(nll_SBfree) == NULL)?0:1;
  result["nll_SB_free"] = nll_SBfree->getVal();
  result["poi_SB_free"] = muVar->getVal();
  result["status_SB_free"] = status*1.0;
  delete nll_SBfree;
  w->saveSnapshot("bestNP", *nuisanceParameters);

  std::cout<<"done fitting data"<<std::endl;
}



void RooStatsHelper::generateToy(RooWorkspace* w ,
        const char* poi_name, double poi_value, int seed,
        map<string, double>& result)
{
    /**
     * before call this function,
     * perform conditional fit first
     * and saved snapshot of NP and GO(global observables)
     * w->saveSnapshot("condNP", *(mc->GetNuisanceParameters()));
     * w->saveSnapshot("condGO", *(mc->GetGlobalObservables()));
     ***/

    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);
    cout<<"in generate toys"<<endl;

    RooStats::ModelConfig *mc = (RooStats::ModelConfig*)w->obj("ModelConfig");
    auto pois = const_cast<RooArgSet*>(mc->GetParametersOfInterest());

    RooRealVar* muVar =(RooRealVar*) w ->var(poi_name);
    double mu_old = muVar->getVal();

    // set other pois to zero
    // and fix, when generating toys,
    // but let them fly afterwards
    TIter itr_poi(pois->createIterator());
    RooRealVar* poi_;
    while ((poi_ = (RooRealVar*) itr_poi())) {
        if(TString(poi_->GetName()) == poi_name) continue;
        poi_->setVal(0.);
        poi_->setConstant(true);
    }
    muVar->setVal(poi_value);
    std::cout<<"before generating toy:"<<std::endl;pois->Print("v");
    RooDataSet* toyData = dynamic_cast<RooDataSet*>(generatePseudoData(w, poi_name, seed));

    //perform S+B fit (everything free) // unconditional
    std::cout<<"S+B free (everything free)"<<std::endl;
    pois->setAttribAll("Constant", false);
    RooNLLVar* nll_SBfree = createNLL(toyData, mc);
    int status = (minimize(nll_SBfree) == NULL)?0:1;
    result["nll_SB_free"] = nll_SBfree->getVal();
    result["poi_SB_free"] = muVar->getVal();
    result["status_SB_free"] = status*1.0;
    delete nll_SBfree;

    //perform S+B fit (with mu fixed) // conditional
    std::cout<<"S+B fixed (everything free except "<<poi_name<<" fixed at "<<poi_value<<")"<<std::endl;
    muVar->setVal(poi_value);
    muVar->setConstant(true);
    RooNLLVar* nll_SBfixed = createNLL(toyData, mc);
    status = (minimize(nll_SBfixed) == NULL)?0:1;
    result["nll_SB_fixed"] = nll_SBfixed->getVal();
    result["poi_SB_fixed"] = muVar->getVal();
    result["status_SB_fixed"] = status*1.0;
    delete nll_SBfixed;

    //perform B-only fit (with mu fixed=0)
    std::cout<<"B-only fixed (everything free except "<<poi_name<<" fixed at "<<0.<<")"<<std::endl;

    itr_poi.Reset();
    while ((poi_ = (RooRealVar*) itr_poi())) {
        poi_->setVal(0.);
        poi_->setConstant(true);
    }

    RooNLLVar* nll_Bonly = createNLL(toyData, mc);
    status = (minimize(nll_Bonly) == NULL)?0:1;
    result["nll_B_fixed"] = nll_Bonly->getVal();
    result["poi_B_fixed"] = muVar->getVal();
    result["status_B_fixed"] = status*1.0;
    delete nll_Bonly;

    /* reset */
    pois->setAttribAll("Constraint", false);
    muVar ->setVal(mu_old);
    delete  toyData;
    cout << "out of generating toys: " << endl;
}

bool RooStatsHelper::ScanPOI(RooWorkspace* ws,
        const string& data_name,
        const string& poi_name,
        int total, double low, double hi,
        TTree* tree)
{
    if(!ws || !tree) {
        log_err("ws or tree is null");
        return false;
    }
    RooStats::ModelConfig*  mc_config = (RooStats::ModelConfig*) ws->obj("ModelConfig");
    RooDataSet* dataset = (RooDataSet*) ws->data(data_name.c_str());
    if(!dataset) {
        log_err("dataset: %s does not exist", data_name.c_str());
        return false;
    }
    RooRealVar* poi = (RooRealVar*) ws->var(poi_name.c_str());
    if(!poi) {
        log_err("POI: %s does not exist", poi_name.c_str());
        return false;
    }
    RooSimultaneous* simPdf = dynamic_cast<RooSimultaneous*>(mc_config->GetPdf());

    // Is it really necessary?
    // RooRealVar* m4l = (RooRealVar*) ws->var("m4l");
    // m4l->setRange("signal", 118, 129);
    //
    bool is_fixed_poi = poi->isConstant();
    poi->setConstant(false);

    double val_nll, val_poi, val_status, obs_sig;
    tree->Branch("NLL", &val_nll, "NLL/D");
    tree->Branch("Status", &val_status, "NLL/D");
    tree->Branch(poi_name.c_str(), &val_poi, Form("%s/D",poi_name.c_str()));
    tree->Branch("obs_sig", &obs_sig, "obs_sig/D");

    TStopwatch timer;
    timer.Start();
    //get best fit
    RooNLLVar* nll = createNLL(dataset, mc_config);

    int status = (minimize(nll) == NULL)?0:1 ;
    timer.Stop();
    cout<<"One fit takes: "<< endl; Helper::printStopwatch(timer);
    timer.Reset(); timer.Start();

    val_nll = nll ->getVal();
    val_poi = poi->getVal();
    val_status = status;
    bool do_subrange = false;
    obs_sig = GetObsNevtsOfSignal(simPdf, poi, mc_config->GetObservables(), do_subrange);
    tree->Fill();

    double step = (hi-low)/total;
    for (int i = 0; i < total+1; i++) {
        double poi_value = low + i*step;
        poi->setVal(poi_value);
        poi->setConstant(true);
        val_status = (minimize(nll)==NULL)?0:1;
        val_nll = nll ->getVal();
        val_poi = poi->getVal();
        obs_sig = GetObsNevtsOfSignal(simPdf, poi, mc_config->GetObservables(), do_subrange);
        tree->Fill();
    }
    if(!is_fixed_poi) poi->setConstant(false);
    delete nll;
    return true;
}

void RooStatsHelper::unfoldConstraints(
        RooArgSet& initial, RooArgSet& final_,
        RooArgSet& obs, RooArgSet& nuis, int& counter
        )
{
  if (counter > 50)
  {
    cout << "ERROR::Couldn't unfold constraints!" << endl;
    cout << "Initial: " << endl;
    initial.Print("v");
    cout << endl;
    cout << "Final: " << endl;
    final_.Print("v");
    exit(1);
  }
  TIterator* itr = initial.createIterator();
  RooAbsPdf* pdf;
  while ((pdf = (RooAbsPdf*)itr->Next()))
  {
    RooArgSet nuis_tmp = nuis;
    RooArgSet constraint_set(*pdf->getAllConstraints(obs, nuis_tmp, false));
    string className(pdf->ClassName());
    if (className != "RooGaussian" && className != "RooLognormal" && className != "RooGamma" && className != "RooPoisson" && className != "RooBifurGauss")
    {
      counter++;
      unfoldConstraints(constraint_set, final_, obs, nuis, counter);
    }
    else
    {
      final_.add(*pdf);
    }
  }
  delete itr;
}

void RooStatsHelper::PrintExpEvts(RooAbsPdf* inputPdf,
        RooRealVar* mu, const RooArgSet* obs, RooAbsData* data)
{
    auto simPdf = dynamic_cast<RooSimultaneous*>(inputPdf);
    if (!simPdf) return;

    const RooCategory& category = dynamic_cast<const RooCategory&>(simPdf->indexCat());
    TIter cat_iter(category.typeIterator());
    RooCatType* obj;
    double old_mu = mu->getVal();
    TList* data_lists = NULL;
    if(data) {
        data_lists = data->split(category, true);
    }
    printf("category allExpected ExpectedSignal ExpectedBackground\n");
    while( (obj= (RooCatType*)cat_iter()) ){
        const char* label_name = obj->GetName();
        RooAbsPdf* pdf = simPdf->getPdf(label_name);
        mu->setVal(0.0);
        double bkg_evts = pdf->expectedEvents(*obs);
        mu->setVal(old_mu);
        double all_evts = pdf->expectedEvents(*obs);
        double data_ch = (data_lists != NULL)? (dynamic_cast<RooDataSet*>(data_lists->At(obj->getVal())))->sumEntries() : 0;

        printf("%s %.3f %.3f %.3f %.3f\n", label_name, all_evts, all_evts-bkg_evts, bkg_evts, data_ch);
    }
    mu->setVal(old_mu);
}

double RooStatsHelper::GetObsNevtsOfSignal(RooSimultaneous* simPdf, RooRealVar* mu, const RooArgSet* obs, bool subrange)
{
    if(!simPdf || !obs || !mu){
        return -999.;
    }
    double result = 0;
    const RooCategory& category = *dynamic_cast<const RooCategory*>(&simPdf->indexCat());
    TIter cat_iter(category.typeIterator());
    RooCatType* obj;
    double old_mu = mu->getVal();
    while( (obj= (RooCatType*)cat_iter()) ){
        const char* label_name = obj->GetName();
        RooAbsPdf* pdf = simPdf->getPdf(label_name);
        mu->setVal(0.0);
        double bkg_evts = pdf->expectedEvents(*obs);
        double fraction_bkg = 1.0;
        if (subrange) fraction_bkg = pdf->createIntegral(*obs, RooFit::Range("signal"))->getVal();
        mu->setVal(old_mu);
        double all_evts = pdf->expectedEvents(*obs);
        double fraction_all = 1.0;
        if(subrange) fraction_all = pdf->createIntegral(*obs, RooFit::Range("signal"))->getVal();
        result += (all_evts*fraction_all - bkg_evts*fraction_bkg);
    }
    mu->setVal(old_mu);
    return result;
}


RooAbsData* RooStatsHelper::generatePseudoData(RooWorkspace* w, const char* poi_name, int seed)
{
    cout<<"in generate pseudo-data with seed: "<< seed << ", poi name: " << poi_name << endl;
    /**
     * before call this function,
     * perform conditional fit first
     * and saved snapshot of NP and GO(global observables)
     * w->saveSnapshot("condNP", *(mc->GetNuisanceParameters()));
     * w->saveSnapshot("condGO", *(mc->GetGlobalObservables()));
     ***/
    auto* poi = (RooRealVar*) w->var(poi_name);
    if(!poi) {
        cout <<"POI is not defined" << endl;
        return NULL;
    }
    double old_poi_val = poi->getVal();
    try {
        w->loadSnapshot("condNP");
        w->loadSnapshot("condGO");
        w->loadSnapshot("condPOI");
    } catch (...) {
        return NULL;
    }

    RooStats::ModelConfig *mc=(RooStats::ModelConfig*)w->obj("ModelConfig");
    RooSimultaneous* combPdf =  (RooSimultaneous*) mc->GetPdf() ;
    RooArgSet* nuisanceParameters = (RooArgSet*) mc->GetNuisanceParameters();
    RooArgSet* globalObservables = (RooArgSet*) mc->GetGlobalObservables();

    /* Fix nuisance parameters,
     * set global observables to np values,
     * need to match the global to the nuisance!
     * make asimov data **/
    SetRooArgSetConst(*nuisanceParameters);
    RooRealVar* nuis;
    RooRealVar* glob;
    TIter nuis_iter (nuisanceParameters->createIterator());

    nuisanceParameters->Print("v");
    globalObservables->Print("v");
    // sort global names according to the length of their names
    vector<RooRealVar*> global_vec;
    TIter glob_iter(globalObservables->createIterator());
    while((glob = (RooRealVar*) glob_iter())){
        global_vec.push_back(glob);
    }
    sort(global_vec.begin(), global_vec.end(), compare_TObject_byName);
    while ((nuis = (RooRealVar*) nuis_iter()))
    {
        if(!nuis || TString(nuis->GetName()).Contains("gamma") ) continue;
        for(auto glob_tmp : global_vec){
            if(TString(glob_tmp->GetName()).Contains(TString(nuis->GetName()).ReplaceAll("alpha_","")))
            {
                glob_tmp->setVal(nuis->getVal());
                break;
            }
        }
    }
    nuisanceParameters->Print("v");
    globalObservables->Print("v");
    mc->GetParametersOfInterest()->Print("v");

    // randomize the global values
    randomizeSet(combPdf, globalObservables, seed);
    SetRooArgSetConst(*globalObservables); //set global observables to the randomized values, and constant
    SetRooArgSetConst(*nuisanceParameters, false);

    const char* weightName="weightVar";
    RooArgSet obsAndWeight;
    obsAndWeight.add(*mc->GetObservables());
    RooRealVar* weightVar = NULL;
    if (!(weightVar = w->var(weightName)))
    {
        w->import(*(new RooRealVar(weightName, weightName, 1.,0,1e3)));
        weightVar = w->var(weightName);
    }
    obsAndWeight.add(*weightVar);

    map<string, RooDataSet*> toyDataMap;
    RooCategory* channelCat = (RooCategory*)&combPdf->indexCat();
    TIterator* iter = channelCat->typeIterator() ;
    RooCatType* tt = NULL;
    int nrIndices = 0;
    while((tt=(RooCatType*) iter->Next())) {
        nrIndices++;
    }
    delete iter;

    poi->setVal(old_poi_val);
    int iFrame=0;
    for (int i=0; i < nrIndices; i++)
    {
        channelCat->setIndex(i);
        iFrame ++;
        // Get pdf associated with state from simpdf
        RooAbsPdf* pdftmp = combPdf->getPdf(channelCat->getLabel()) ;
        if(!pdftmp){
            cout<<"ERRORs no pdf associated with "<< channelCat ->getLabel()<<endl;
        }
        // Generate observables defined by the pdf associated with this state

        RooArgSet* obstmp = pdftmp->getObservables(*mc->GetObservables()) ;
        RooDataSet* obsDataUnbinned = pdftmp->generate(*obstmp, RooFit::Extended());

        int nentries = obsDataUnbinned->numEntries();
        if(TMath::IsNaN(obsDataUnbinned->sumEntries())){
            cout << "sum entries is nan"<<endl;
            exit(1);
        }
        cout <<channelCat->getLabel() << ": nentries: " << nentries << " " << obsDataUnbinned->isWeighted() << endl;

        //Create a new dataset with the weightVar inside.
        RooDataSet* data_category = new RooDataSet(
                Form("combtoyData%d",iFrame),
                Form("combtoyData%d",iFrame),
                RooArgSet(obsAndWeight,*channelCat),
                RooFit::WeightVar(*weightVar)
                );
        for(int i = 0; i < nentries; i ++){
            auto args = obsDataUnbinned->get(i);
            // cout <<"value: "<< ((RooRealVar*)args->first())->getVal()<< "weight: " <<  obsDataUnbinned->weight() << endl; // only for debugs
            data_category->add(*args, obsDataUnbinned->weight());
        }
        toyDataMap[string(channelCat->getLabel())] = data_category;//tempData;
    }
    RooDataSet* toyData = new RooDataSet(
            Form("toyData_%d",seed),
            Form("toyData_%d",seed),
            RooArgSet(obsAndWeight,*channelCat),
            RooFit::Index(*channelCat),
            RooFit::Import(toyDataMap),
            RooFit::WeightVar(*weightVar)
            );

    try {
        w->loadSnapshot("nominalNP");
        w->loadSnapshot("nominalGO");
    } catch(...) {
        log_err("save a snapshot of nominal values");
    }

    // clean the data map
    for(auto item : toyDataMap) {
        if(item.second) delete item.second;
    }
    toyDataMap.clear();

    // Set Observable to 0.
    // SetRooArgSetValue(*globalObservables, 0.);
    // w->loadSnapshot("condGO");
    toyData->Print("v");
    return toyData;
}

bool RooStatsHelper::fixTermsWithPattern(RooStats::ModelConfig* mc, const char* pat)
{
    RooArgSet nuis(*mc->GetNuisanceParameters());
    TIter iter(nuis.createIterator());
    RooRealVar* par = NULL;
    while( (par = (RooRealVar*) iter()) ) {
        if(string(par->GetName()).find(pat) != string::npos) {
            par->setConstant();
        }
    }
    return true;
}

void RooStatsHelper::fixVariables(RooWorkspace* workspace, const string& options, RooStats::ModelConfig* mc)
{
    // options can be like: "mG:750,GkM:0.02,noSys"
    // or "mH=750,width=200,noSys"
    if (!workspace || options == "") return;
    vector<string> tokens;
    Helper::tokenizeString(options, ',', tokens);
    if (tokens.size() < 1) return ;

    for(auto iter = tokens.begin(); iter != tokens.end(); iter++)
    {
        string token(*iter);
        size_t delim_pos = token.find(':');
        if(delim_pos == string::npos){
            delim_pos = token.find('=');
        }
        string var_name = token;
        double var_val = nan("NaN");
        if(delim_pos != string::npos){
            var_name = token.substr(0, delim_pos);
            var_val = atof( token.substr(delim_pos+1, token.size()).c_str());
        } else {
            if(token.find("noSys") != string::npos && mc) {
                RooArgSet nuis(*mc->GetNuisanceParameters());
                TIter iter(nuis.createIterator());
                RooRealVar* par = NULL;
                while( (par = (RooRealVar*) iter()) ) {
                    par->setConstant();
                    log_info("%s fixed to %.2f", par->GetName(), par->getVal());
                }
                continue;
            }
        }

        if(token.find("gamma_stat") != string::npos && mc) {
            const char* pat = "gamma_stat";
            RooArgSet nuis(*mc->GetNuisanceParameters());
            TIter iter(nuis.createIterator());
            RooRealVar* par = NULL;
            while( (par = (RooRealVar*) iter()) ) {
                if(string(par->GetName()).find(pat) != string::npos) {
                    if(!std::isnan(var_val)) { // use std::isnan to avoid ambigulity
                        par->setVal(var_val);
                    }
                    par->setConstant();
                    log_info("%s fixed to %.2f", par->GetName(), par->getVal());
                }
            }
            continue;
        }

        auto par = (RooRealVar*) workspace->var(var_name.c_str());
        if(!par) {
            log_warn("%s does not exist", var_name.c_str());
            continue;
        } else {
            if(!std::isnan(var_val)) {
                double low_val = var_val - 1, hi_val = var_val + 1;
                par->setRange(low_val, hi_val);
                par->setVal(var_val);
            }
            par->setConstant();
            log_info("%s fixed to %.2f", var_name.c_str(), par->getVal());
        }
    }
}

bool RooStatsHelper::CheckNuisPdfConstraint(const RooArgSet* nuis, const RooArgSet* pdfConstraint)
{
    if(!nuis || nuis->getSize() < 1) {
        cout << "no nusiance parameter is found" << endl;
        return true;
    }
    if(nuis->getSize() != pdfConstraint->getSize())
    {
        log_err("number of nuisance %d; pdfConstraint size: %d", nuis->getSize(), pdfConstraint->getSize());

        TIterator* iterNuis = nuis->createIterator();
        TIterator* iterPdfCons = pdfConstraint->createIterator();
        RooRealVar* nuisVal;
        RooAbsPdf* pdfCon;
        while((nuisVal = (RooRealVar*)iterNuis->Next())){
            TString constrainName(Form("%sConstraint",nuisVal->GetName()));
            if(pdfConstraint->find(constrainName.Data()) == NULL){
                cout <<nuisVal->GetName()<<" not Constraint"<<endl;
            }
        }
        while((pdfCon = (RooAbsPdf*)iterPdfCons->Next())){
            TString nuisName(pdfCon ->GetName());
            nuisName.ReplaceAll("Constraint","");
            if(nuis->find(nuisName.Data()) == NULL){
                cout <<nuisName<<" is omitted"<<endl;
            }
        }
        return false;
    } else {
        cout <<"total number of nuisance parameters: "<<nuis->getSize()<<endl;
    }
    return true;
}

void RooStatsHelper::setOtherPOIs(
        RooArgSet* pois, const string& veto,
        float other_poi_value, bool is_other_poi_const)
{
    TIterator* itr_poi = pois->createIterator();
    RooRealVar* poi_;
    while ((poi_ = (RooRealVar*) itr_poi->Next())) {
        if (TString(veto.c_str()).Contains(poi_->GetName())) continue;
        poi_->setVal(other_poi_value);
        poi_->setConstant(is_other_poi_const);
    }
    delete itr_poi;
}

void RooStatsHelper::SetRooArgSetValue(RooArgSet& argset, float value){
    TIter iter(argset.createIterator());
    RooRealVar* arg = NULL;
    while ( (arg == (RooRealVar*) iter()) ){
        arg->setVal(value);
    }
}

bool RooStatsHelper::compare_TObject_byName(TObject* a, TObject* b){
    return strlen(a->GetName()) < strlen(b->GetName());
}
