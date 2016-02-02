#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooHist.h"
using namespace RooFit ;


int main()
{

  
  RooRealVar x("x","x",-10,10) ;

 
  RooRealVar sigma("sigma","sigma",3,0.1,10) ;
  RooRealVar mean("mean","mean",0,-10,10) ;
  RooGaussian gauss("gauss","gauss",x,RooConst(0),sigma);

  RooDataSet* data = gauss.generate(x,100);

  sigma=3.15 ;


  
  RooPlot* frame1 = x.frame(Title("Data with distorted Gaussian pdf"),Bins(10)) ;
  data->plotOn(frame1, DataError(RooAbsData::SumW2)) ; // DataError(RooAbsData::SumW2) legt den Typ der Verteilung 
                                                       //   Fehler fest!
  gauss.plotOn(frame1) ;


  
  std::cout << "chi^2 = " << frame1->chiSquare() << std::endl ;


  
  RooHist* hresid = frame1->residHist() ;

  
  RooHist* hpull = frame1->pullHist() ;

  
  RooPlot* frame2 = x.frame(Title("Residual Distribution")) ;
  frame2->addPlotable(hresid,"P") ;

  
  RooPlot* frame3 = x.frame(Title("Pull Distribution")) ;
  frame3->addPlotable(hpull,"P") ;



  TCanvas* c = new TCanvas("rf109_chi2residpull","rf109_chi2residpull",900,300) ;
  c->Divide(3) ;
  c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame1->GetYaxis()->SetTitleOffset(1.6) ; frame1->Draw() ;
  c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
  c->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.6) ; frame3->Draw() ;
  c->SaveAs("/home/chasenberg/bachelor-template/pull.pdf");
  return 0;
}
