#include "string"
#include "vector"

#include "RooBMixDecay.h"
#include "RooCategory.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGenericPdf.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooSimPdfBuilder.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooSimWSTool.h"

#include "doofit/builder/EasyPdf/EasyPdf.h"
#include "doofit/config/CommonConfig.h"
#include "doofit/plotting/fitresult/FitResultPrinter.h"
#include "doofit/fitter/easyfit/EasyFit.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStd.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStdConfig.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStd.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStdConfig.h"




int main(int argc, char *argv[]) {

  doofit::builder::EasyPdf *epdf = new doofit::builder::EasyPdf();

  //Observable 
  epdf->Var("obsMass").setVal(5280);
  epdf->Var("obsMass").setRange(5200.,5500.);
  epdf->Var("obsTime").setVal(1.3);
  epdf->Var("obsTime").setRange(0.3, 18.3);
  
  //Zusammenfassen der Parameter in einem RooArgSet
  RooArgSet Observables;
  Observables.add(RooArgSet(epdf->Var("obsMass"), epdf->Var("obsTime")));

  //Anlegen der Parameter
  /*epdf->Var("pdf_sig_mass_mean").setRange(5230, 5330);
  epdf->Var("pdf_sig_mass_mean").setVal(5280);
  epdf->Var("pdf_sig_mass_width").setRange(1.0, 15);
  epdf->Var("pdf_sig_mass_width").setVal(10);
  epdf->Var("pdf_bkg_mass_lambda").setRange(-0.0074, 0.0);
  epdf->Var("pdf_bkg_mass_lambda").setVal(-0.001);
  epdf->Var("sig_Yield").setRange(0.0, 10000);
  epdf->Var("sig_Yield").setVal(5000);
  epdf->Var("bkg_Yield").setRange(0.0, 10000); 
  epdf->Var("bkg_Yield").setVal(5000);*/
  
  ///////////////////Generiere PDF's/////////////////////
  
  //Signal PDF's//
  epdf->Gaussian("pdf_sig_mass_gauss", epdf->Var("obsMass"), epdf->Var("pdf_sig_mass_mean"),epdf->Var("pdf_sig_mass_width"));
  //Untergrund PDF
  epdf->Exponential("pdf_bkg_mass_expo", epdf->Var("obsMass"),epdf->Var("pdf_bkg_mass_lambda"));
  

  //Zerfalls PDF's//
  //Resolution Model
  epdf->GaussModel("pdf_resolution", epdf->Var("obsTime"), epdf->Var("pdf_res_mean"), epdf->Var("obsTimeErr"));
  //Zerfalls PDF des Untergrunds
  epdf->Decay("pdf_bkg_time_decay", epdf->Var("obsTime"),epdf->Var("pdf_bkg_time_tau"), epdf->Model("pdf_resolution"));
  //Zerfall PDF des Signals 
  epdf->Decay("pdf_sig_time_decay", epdf->Var("obsTime"),epdf->Var("pdf_sig_time_tau"), epdf->Model("pdf_resolution"));


  //Multipliziere Signal und Untergrund PDF mit ihrer jeweiligen Zerfalls PDF//
  //Untergrund * Zerfall
  epdf->Product("pdf_bkg_mass*pdf_bkg_time_decay", RooArgSet(epdf->Pdf("pdf_bkg_mass_expo"), epdf->Pdf("pdf_bkg_time_decay")));
  //Signal * Zerfall
  epdf->Product("pdf_sig_mass_gauss*pdf_sig_time_decay", RooArgSet(epdf->Pdf("pdf_sig_mass_gauss"),epdf->Pdf("pdf_sig_time_decay")));




  //Addiere PDF's
  epdf->Add("pdf_total", RooArgSet(epdf->Pdf("pdf_sig_mass_gauss*pdf_sig_time_decay"), epdf->Pdf("pdf_bkg_mass*pdf_bkg_time_decay")), RooArgSet(epdf->Var("bkg_Yield"),epdf->Var("sig_Yield")));
  
 


  RooWorkspace ws;
  ws.import(epdf->Pdf("pdf_total"));
  ws.defineSet("Observables",Observables, true);

  ws.Print();

   doofit::config::CommonConfig cfg_com("common");
   cfg_com.InitializeOptions(argc, argv);
   doofit::toy::ToyFactoryStdConfig cfg_tfac("toyfac");
   cfg_tfac.InitializeOptions(cfg_com);
   doofit::toy::ToyStudyStdConfig cfg_tstudy("toystudy");
   cfg_tstudy.InitializeOptions(cfg_tfac);

   // set a previously defined workspace to get PDF from (not mandatory, but convenient)
   cfg_tfac.set_workspace(&ws);

   // Check for a set --help flag and if so, print help and exit gracefully
   // (recommended).
   cfg_com.CheckHelpFlagAndPrintHelp();

   // More custom code, e.g. to set options internally.
   // Not required as configuration via command line/config file is enough.
   

   // Print overview of all options (optional)
   //7cfg_com.PrintAll();
   
   // Initialize the toy factory module with the config objects and start
   // generating toy samples.
   doofit::toy::ToyFactoryStd tfac(cfg_com, cfg_tfac);
   doofit::toy::ToyStudyStd tstudy(cfg_com, cfg_tstudy);
   

   for(int i=0;i<10  ;i++)
   {
   RooDataSet* data = tfac.Generate();
   epdf->Pdf("pdf_total").getParameters(data)->readFromFile("/home/chasenberg/Repository/bachelor-template/parameters_generate.txt");
   epdf->Pdf("pdf_total").getParameters(data)->writeToFile("/home/chasenberg/Repository/bachelor-template/parameters_generate.txt.new");
   
   // fitting on the generated dataset is your responsibility
   RooFitResult* fit_result = epdf->Pdf("pdf_total").fitTo(*data, RooFit::Save(true));

   //Speichern der Ergebnisse
  tstudy.StoreFitResult(fit_result);
   }
  tstudy.FinishFitResultSaving();
  tstudy.ReadFitResults();
  tstudy.EvaluateFitResults();
  tstudy.PlotEvaluatedParameters();
  
  
  
 }
