#include "RooStats/BayesianCalculator.h"
#include "RooUniform.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
#include "RooAbsArg.h"
#include "RooProdPdf.h"
#include "RooAddPdf.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TFile.h"
#include "TROOT.h"
#include "RooPlot.h"
#include "TCanvas.h"

using namespace RooStats;

void setParsRanges(RooArgSet *pars, int ns){
   TIterator *it = pars->createIterator();
   RooAbsArg *aarg;
   while (aarg = (RooAbsArg*)it->Next()) { 
	RooRealVar *arg = dynamic_cast<RooRealVar*>(aarg);
	arg->setMax(arg->getVal()+ns*arg->getError());
	arg->setMin(arg->getVal()-ns*arg->getError());
   }

}

void quickTest(){

  //gROOT->SetBatch(1);
  gROOT->ProcessLine(".x paperStyle.C");

  TFile *fi = TFile::Open("envelopews_wsignal_toy1_110to150.root");
  RooWorkspace *multipdf = (RooWorkspace*)fi->Get("multipdf");
  RooRealVar *x    = multipdf->var("CMS_hgg_mass");

  // The data we will use (a Toy dataset);
  RooDataHist *datatoy = (RooDataHist *)multipdf->data("roohist_data_mass_cat1_toy1_cutrange__CMS_hgg_mass");
  
  // Build signal model
  RooRealVar mean("mean","mean",125); mean.setConstant();
  RooRealVar sigma("sigma","sigma",1.19); sigma.setConstant();
  RooRealVar nsignal_const("nsignal_const","nsignal_const",50.8); nsignal_const.setConstant();

  RooGaussian sig_pdf("gaus","gaus",*x,mean,sigma);
  
  // Should add this to the RooWorkspace
  RooRealVar nbkg("nbkg","nbkg",datatoy->sumEntries(),0,10000);

  RooRealVar mu("r","r",1,0,5);
  RooFormulaVar nsig("nsig","nsig","@0*@1",RooArgList(mu,nsignal_const));
  //mu.setConstant();
  RooPlot *plot = mu.frame();

  std::string bkgfunctions[2] = {"env_pdf_1_8TeV_pow1","env_pdf_1_8TeV_exp1"};
  int color[2] = {4,2};

  for (int id=0;id<2;id++){

  // start wih a model !
   RooAbsPdf *bkg_pdf = multipdf->pdf(bkgfunctions[id].c_str());
   RooAddPdf spdf_a("Asplusb","Asplusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));

   // start someplace reasonable
   spdf_a.fitTo(*datatoy);
   //RooPlot *pl1 = x->frame();
   //datatoy->plotOn(pl1);
   //spdf_a.plotOn(pl1);
   //pl1->Draw();

   RooWorkspace *wsp = new RooWorkspace("NewWorkspace");

   wsp->import(spdf_a,RooFit::RecycleConflictNodes());
   // Now start messing around with the parameters
   int nsigma = 4;
   setParsRanges(bkg_pdf->getParameters(*datatoy),4);
   //wsp->var("env_pdf_1_8TeV_pow1_p1")->setMin(wsp->var("env_pdf_1_8TeV_pow1_p1")->getVal()-4*wsp->var("env_pdf_1_8TeV_pow1_p1")->getError());
   //wsp->var("env_pdf_1_8TeV_pow1_p1")->setMax(wsp->var("env_pdf_1_8TeV_pow1_p1")->getVal()+4*wsp->var("env_pdf_1_8TeV_pow1_p1")->getError());
   wsp->var("nbkg")->setMin(wsp->var("nbkg")->getVal()-nsigma*wsp->var("nbkg")->getError());
   wsp->var("nbkg")->setMax(wsp->var("nbkg")->getVal()+nsigma*wsp->var("nbkg")->getError());

   // Prior for signal strength + nuisances 
   RooUniform uniform_mu("p_mu","p_nbkg",RooArgSet(*(wsp->var(mu.GetName()))));

   // Nuisance parameters will be background parameters
   RooArgSet nuisances(*(bkg_pdf->getParameters(*datatoy)));
   nuisances.add(*(wsp->var(nbkg.GetName())));
   RooUniform prior_pdf("prior_nuis","nuisance prior",nuisances);

   // PDF is usual s+b pdf * constraint (priors)
   RooProdPdf spdf("splusb","splusb",RooArgList(spdf_a,prior_pdf));
   wsp->import(spdf,RooFit::RecycleConflictNodes());

   // Model Config 
   RooArgSet poi(*(wsp->var(mu.GetName())));

   RooStats::ModelConfig mc("mc","mc",wsp);
   mc.SetObservables(RooArgSet(*x));
   mc.SetPdf(*(wsp->pdf(spdf.GetName())));
   mc.SetPriorPdf(uniform_mu);
   mc.SetNuisanceParameters(nuisances);
   mc.SetParametersOfInterest(poi);
   mc.SetGlobalObservables(RooArgSet());

   mc.Print("v");

   RooStats::BayesianCalculator bc(*datatoy,mc);
   bc.SetConfidenceLevel(0.683);
   bc.SetScanOfPosterior(50);
   bc.SetIntegrationType("MISER");
   RooAbsReal *pdfpost =  bc.GetPosteriorFunction();
   pdfpost->Print("vT");
   pdfpost->plotOn(plot,RooFit::LineColor(color[id]));
//  spdf.Print("v");
  }

  //RooPlot * plot = bc.GetPosteriorPlot();
  TCanvas *can = new TCanvas();
  plot->Draw();
}
