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
#include "TMath.h"
#include <stdio.h>
#include "getChisq.C"

using namespace RooStats;


double gstep = 0.02;

void setParsRanges(RooArgSet *pars, int ns){
   TIterator *it = pars->createIterator();
   RooAbsArg *aarg;
   while (aarg = (RooAbsArg*)it->Next()) { 
	RooRealVar *arg = dynamic_cast<RooRealVar*>(aarg);
	arg->setMax(arg->getVal()+ns*arg->getError());
	arg->setMin(arg->getVal()-ns*arg->getError());
   }

}


void makeGraph(TGraph *gr, RooRealVar &poi, RooDataHist &dat, RooRealVar *var, RooAbsPdf *pdf, RooAbsPdf &pi, double step = gstep){
   
  int p = 0;
  double range = poi.getMax() - poi.getMin();
  for (double x = poi.getMin(); x<=poi.getMax();x+=step){
	poi.setVal(x);
	double nllvalue = getChisq(dat,*pdf,*var,0);
	 gr->SetPoint(p,x,10E45*(1./range)*TMath::Exp(-0.5*nllvalue));
	p++;
  }

  gr->SetLineWidth(2);
  gr->GetXaxis()->SetTitle(poi.GetName());
}

void test_cartoon(){

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

  RooRealVar mu("r","r",1,-1,3);
  RooFormulaVar nsig("nsig","nsig","@0*@1",RooArgList(mu,nsignal_const));
  //mu.setConstant();
  RooPlot *plot = mu.frame();

  //std::string bkgfunctions[2] = {"env_pdf_1_8TeV_pow1","env_pdf_1_8TeV_exp1"};
  int color[5]     = {2,2,4,2,2};
  //int bkgvalues[5] = {-4.9-10*0.1775,-4.9-5*0.1775,-4.9,-4.9+5*0.1775,-4.9+10*0.1775};  //change the slope. 
  
  RooAbsPdf *bkg_pdf = multipdf->pdf("env_pdf_1_8TeV_pow1");
  RooAddPdf spdf_a("Asplusb","Asplusb",RooArgList(sig_pdf,*bkg_pdf),RooArgList(nsig,nbkg));
  spdf_a.fitTo(*datatoy);

  double cen = nbkg.getVal();
  double err = nbkg.getError();
  double bestfit = mu.getVal();
  int bkgvalues[5] = {cen-2*err,cen-1*err,cen,cen+1*err,cen+2*err};  //change the slope. 

   RooWorkspace *wsp = new RooWorkspace("NewWorkspace");
   wsp->import(spdf_a);
   RooArgSet nuisances(*(wsp->var("nbkg")));
   RooUniform prior_pdf("prior_nuis","nuisance prior",nuisances);
   RooProdPdf spdf("splusb","splusb",RooArgList(*(wsp->pdf(spdf_a.GetName())),prior_pdf));
   wsp->import(spdf,RooFit::RecycleConflictNodes());

   RooUniform uniform_mu("p_mu","p_nbkg",RooArgSet(*(wsp->var(mu.GetName()))));
   wsp->import(uniform_mu);
   // Model Config 
   RooArgSet poi(*(wsp->var(mu.GetName())));

   RooStats::ModelConfig mc("mc","mc",wsp);
   mc.SetObservables(RooArgSet(*(wsp->var(x->GetName()))));
   mc.SetPdf(*(wsp->pdf(spdf.GetName())));
   mc.SetPriorPdf(*(wsp->pdf(uniform_mu.GetName())));
   mc.SetNuisanceParameters(nuisances);
   mc.SetParametersOfInterest(poi);
   mc.SetGlobalObservables(RooArgSet());

   // Nuisance parameters will be background parameters
   //nuisances.add(*(wsp->var(nbkg.GetName())));
   //RooArgSet nuisances(*(wsp->var(nbkg.GetName())));

   // PDF is usual s+b pdf * constraint (priors)


  // lets make our own guy?
  RooAbsReal *nll = wsp->pdf(spdf_a.GetName())->createNLL(*datatoy);

  //TGraph *grs[5] = {0,0,0,0,0};
  std::cout<< "READY TO GO " << std::endl; 
  //TGraph *tot = new TGraph();
  int ip=0;
  for (double xx = mu.getMin(); xx<=mu.getMax();xx+=gstep){
	ip+=1;
  }
  TH1D *tot = new TH1D("sum","",ip,mu.getMin()-gstep/2,mu.getMax()-gstep/2);
  TH1D *cenh = new TH1D("cen","",ip,mu.getMin()-gstep/2,mu.getMax()-gstep/2);
  RooPlot *pl1 = x->frame();
 
  TCanvas *can = new TCanvas();
  wsp->var("nbkg")->setConstant(true);
  std::cout<< "READY TO GO ECHECK POITNS" << std::endl;

 
  for (int id=0;id<5;id++){
   // start someplace reasonable
   datatoy->plotOn(pl1,RooFit::Binning(80));
   int nsigma = 3;
   setParsRanges(bkg_pdf->getParameters(*datatoy),3);

  std::cout<< "READY TO GO MAKING STUFF" << std::endl; 
   wsp->var("env_pdf_1_8TeV_pow1_p1")->setConstant(true);
   wsp->var("nbkg")->setConstant(true);
   wsp->var("nbkg")->setVal(bkgvalues[id]);
 
   TGraph *gr = new TGraph();
  std::cout<< "READY TO GO PLOTTING A GRAPH " << std::endl; 
   makeGraph(gr,*(wsp->var(mu.GetName())),*datatoy,wsp->var(x->GetName()),wsp->pdf(spdf_a.GetName()),*(wsp->pdf(uniform_mu.GetName())));
  
  std::cout<< "READY TO GO PLOTTING A GRAPH (DONE) " << std::endl; 
   double xx,yy;
   for (int i = 1; i<=tot->GetNbinsX();i++){
	//tot->GetPoint(i,xx,yy);
	//tot->SetPoint(i,xx,yy+gr->Eval(xx)/5);
	std::cout<< tot->GetBinCenter(i)<< " " <<gr->Eval(tot->GetBinCenter(i))<< std::endl;
	tot->SetBinContent(i,tot->GetBinContent(i)+gr->Eval(tot->GetBinCenter(i))/5);
	if (id==2)cenh->SetBinContent(i,gr->Eval(cenh->GetBinCenter(i)));
   }
   
   //grs[id] = gr;

   if (id!=2) { 
	gr->SetLineStyle(2);
//	gr->SetLineColor(2);
   }
   gr->SetLineColor(color[id]);
   if (id==0) gr->Draw("AL");
   else gr->Draw("L");

   //mc.Print("v");
   mu.setVal(bestfit);
   std::cout << "nBKG = " << (wsp->var("nbkg"))->getVal() << std::endl;
   if (id==2) wsp->pdf(spdf_a.GetName())->plotOn(pl1,RooFit::LineColor(color[id]),RooFit::Normalization((wsp->var("nbkg"))->getVal()+nsig.getVal(),RooAbsReal::NumEvent));
   else wsp->pdf(spdf_a.GetName())->plotOn(pl1,RooFit::LineColor(color[id]),RooFit::LineStyle(2),RooFit::Normalization((wsp->var("nbkg"))->getVal()+nsig.getVal(),RooAbsReal::NumEvent));
  }

  TCanvas *can3 = new TCanvas();
  // Finaly bayescalc this guy
  //RooPlot *plot = mu.frame();
  wsp->var("nbkg")->setConstant(false);
  // fit again 
  wsp->pdf(spdf_a.GetName())->fitTo(*datatoy);

  setParsRanges(&RooArgSet(*(wsp->var("nbkg"))),4);

  RooStats::BayesianCalculator bc(*datatoy,mc);
  bc.SetConfidenceLevel(0.683);
  bc.SetScanOfPosterior(25);
  bc.SetIntegrationType("MISER");
  RooAbsReal *pdfpost =  bc.GetPosteriorFunction();
  //mu.setRange("RNGE",mu.getMin(),mu.getMax());
  //pdfpost->plotOn(plot,RooFit::LineColor(kMagenta),RooFit::Normalization(1.,RooAbsReal::NumEvent));
  //plot->Draw();
  // Make the curve ourselves because RooFit SUCKS?!
 
  double wd = tot->GetBinWidth(1);
  TH1D *tot_B = (TH1D*) tot->Clone();

  for (int b=1;b<=tot->GetNbinsX();b++){
	wsp->var(mu.GetName())->setVal(tot_B->GetBinCenter(b));
	tot_B->SetBinContent(b,pdfpost->getVal());
  }
  tot->SetMaximum(0.12);
  tot_B->Scale(1./(wd*tot_B->Integral()));
  tot->Scale(1./(wd*tot->Integral()));
  cenh->Scale(1./(wd*cenh->Integral()));

  tot->GetXaxis()->SetTitle("#mu");
  tot->GetYaxis()->SetTitle("Posterior probability density");

  //tot->SetLineStyle(2);
  tot->SetLineWidth(3);
  tot_B->SetLineWidth(3);
  cenh->SetLineWidth(3);
  tot_B->SetLineColor(kMagenta);
  cenh->SetLineColor(kBlue);

  tot->Draw("hist");
  ///tot->Draw("histsame");
  tot_B->Draw("histsameL");
  cenh->Draw("histsameL");

  TCanvas *c = new TCanvas();
  pl1->Draw();
}
