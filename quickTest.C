#include "RooStats/BayesianCalculator.h"
#include "RooUniform.h"
using namespace RooStats;
void quickTest(){

  //gROOT->SetBatch(1);
  gROOT->ProcessLine(".x paperStyle.C");

  TFile *fi = TFile::Open("envelopews_wsignal_toy1_110to150.root");
  RooWorkspace *multipdf = fi->Get("multipdf");
  RooRealVar *x    = multipdf->var("CMS_hgg_mass");

  // The data we will use (a Toy dataset);
  RooDataHist *datatoy = multipdf->data("roohist_data_mass_cat1_toy1__CMS_hgg_mass");

  // Build signal model
  RooRealVar mean("mean","mean",125); mean.setConstant();
  RooRealVar sigma("sigma","sigma",1.19); sigma.setConstant();
  RooRealVar nsignal_const("nsignal_const","nsignal_const",50.8); nsignal_const.setConstant();

  RooGaussian sig_pdf("gaus","gaus",*x,mean,sigma);
  
  // Should add this to the RooWorkspace
  RooRealVar nbkg("nbkg","nbkg",datatoy->sumEntries(),0,10E6);

  RooRealVar mu("r","r",1,-10,10);
  RooFormulaVar nsig("nsig","nsig","@0*@1",RooArgList(mu,nsignal_const));
  mu.setConstant();

 // start wih a model !
  RooAbsPdf *bkg_pdf = multipdf->pdf("env_pdf_1_8TeV_pow1");
  RooAddPdf spdf_a("Asplusb","Asplusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));

  multipdf->import(spdf,RooFit::RecycleConflictNodes());

  // make priors
  multipdf->factory("Uniform::prior_mu(mu,0,10)");
  multipdf->factory("Uniform::prior_nbkg(nbkg,0,10000)");
  multipdf->factory("Uniform::prior_par(env_pdf_1_8TeV_pow1_p1,-10,1)");
  multipdf->factory("PROD::priorpdf(prior_nbkg,prior_par,prior_mu)");

  RooProdPdf spdf("splusb","splusb",RooArgList(spdf_a,*multipdf->pdf("priorpdf")));

  // Nuisance parameters will be background parameters
  RooArgSet nuisances(*(bkg_pdf->getParameters(*datatoy)));
  nuisances.add(nbkg);
  nuisances.Print("V");

  multipdf->pdf("priorpdf")->Print();
  RooStats::BayesianCalculator bc(*datatoy,spdf,RooArgSet(mu),(*multipdf->pdf("priorpdf")),RooArgSet(nuisances));
  RooAbsReal * posterior = bc.GetPosteriorFunction();
  RooPlot *pl = x->frame();
  posterior->plotOn(pl);
  pl->Draw();
}
