#include "RooStats/BayesianCalculator.h"
#include "RooUniform.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooAbsPdf.h"
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
  RooRealVar nbkg("nbkg","nbkg",datatoy->sumEntries(),0,1000);

  RooRealVar mu("r","r",1,0,5);
  RooFormulaVar nsig("nsig","nsig","@0*@1",RooArgList(mu,nsignal_const));
  //mu.setConstant();

 // start wih a model !
  RooAbsPdf *bkg_pdf = multipdf->pdf("env_pdf_1_8TeV_pow1");
  RooAddPdf spdf_a("Asplusb","Asplusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));
  //RooAddPdf spdf("Asplusb","Asplusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));


  // make priors
 // multipdf->factory("Uniform::prior_mu(mu,0,10)");
 // multipdf->factory("Uniform::prior_nbkg(nbkg,0,10000)");
 // multipdf->factory("Uniform::prior_par(env_pdf_1_8TeV_pow1_p1,-10,1)");
//  multipdf->factory("Uniform::prior_mu");
//  multipdf->factory("Uniform::prior_nbkg");
//  multipdf->factory("Uniform::prior_par");
  RooUniform uniform_nbkg("p_nbkg","p_nbkg",nbkg);
  RooUniform uniform_par("p_par","p_nbkg",*(multipdf->var("env_pdf_1_8TeV_pow1_p1")));
  RooUniform uniform_mu("p_mu","p_nbkg",mu);
  multipdf->var("env_pdf_1_8TeV_pow1_p1")->setMin(-6);
  multipdf->var("env_pdf_1_8TeV_pow1_p1")->setMax(-2);
  multipdf->Print();

  //multipdf->factory("PROD::priorpdf(prior_nbkg,prior_par)");
  RooProdPdf prior_pdf("prior_nuis","nuisance prior",RooArgList(uniform_nbkg,uniform_par));
  //RooProdPdf prior_pdf("prior_nuis","nuisance prior",RooArgList(uniform_nbkg));

  RooProdPdf spdf("splusb","splusb",RooArgList(spdf_a,prior_pdf));
  //multipdf->import(spdf,RooFit::RecycleConflictNodes());
  spdf_a.fitTo(*datatoy);
  RooPlot *pl1 = x->frame();
  datatoy->plotOn(pl1);
  spdf_a.plotOn(pl1);
   pl1->Draw();

  TCanvas *can = new TCanvas();
  // Nuisance parameters will be background parameters
  RooArgSet nuisances(*(bkg_pdf->getParameters(*datatoy)));
  nuisances.add(nbkg);
  nuisances.Print("V");

//  multipdf->pdf("prior_mu")->Print();
 datatoy->Print();
spdf.Print();
mu.Print();
uniform_mu.Print();

  RooArgSet poi(mu);

  RooStats::ModelConfig mc("mc","mc",multipdf);
  mc.SetPdf(spdf);
  mc.SetPriorPdf(uniform_mu);

  //RooStats::BayesianCalculator bc(*datatoy,spdf,poi,uniform_mu,(&nuisances));
  //bc.Print("v");
  //bc.SetConfidenceLevel(0.683);
  //RooAbsPdf * posterior = bc.GetPosteriorPdf();
  //RooPlot *pl =  bc.GetPosteriorPlot() ;//mu.frame();
//  posterior->plotOn(pl);
  RooStats::MCMCCalculator MCMC(*data,mc);
  MCMC.

 pl->Draw();
}
