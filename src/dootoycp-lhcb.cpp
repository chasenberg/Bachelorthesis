#include "string"
#include "vector"
#include "iostream"
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
#include "RooExtendPdf.h"
#include "doofit/builder/EasyPdf/EasyPdf.h"
#include "doofit/config/CommonConfig.h"
#include "doofit/plotting/fitresult/FitResultPrinter.h"
#include "doofit/fitter/easyfit/EasyFit.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStd.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStdConfig.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStd.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStdConfig.h"

#include "doofit/plotting/Plot/Plot.h"
#include "doofit/plotting/Plot/PlotConfig.h"

#include "RooAbsPdf.h"
#include "Urania/DecRateCoeff.h"

int main(int argc, char* argv[]) {

     doofit::builder::EasyPdf *epdf = new doofit::builder::EasyPdf();
     
     //Mass 
     //epdf->Var("obsMass").setVal(5280);
     //epdf->Var("obsMass").setRange(5200.,5500.);
     //epdf->Var("obsMass").SetTitle("Masse");


	//decay time
	epdf->Var("obsTime");
	epdf->Var("obsTime").SetTitle("t_{#kern[-0.3]{B}^{#kern[-0.3]{0}}}");
	epdf->Var("obsTime").setUnit("ps");
	epdf->Var("obsTime").setRange(0.,16.);

	// tag, respectively the initial state of the produced B meson
	epdf->Cat("obsTag");
	epdf->Cat("obsTag").defineType("B",1);
	epdf->Cat("obsTag").defineType("Bbar",-1);

	//finalstate
	epdf->Cat("catFinalState");
	epdf->Cat("catFinalState").defineType("f",1);
	epdf->Cat("catFinalState").defineType("fbar",-1);

     epdf->Var("sig_Yield");
     epdf->Var("sig_Yield").setVal(5000);

	//Zusammenfassen der Parameter in einem RooArgSet
     //RooArgSet Observables;
     //Observables.add(RooArgSet(epdf->Var("obsMass"), epdf->Var("obsTime"), epdf->Cat("catFinalState"), epdf->Cat("obsTag") ));


     //Koeffizienten
	DecRateCoeff *coeff_c = new DecRateCoeff("coef_cos","coef_cos",DecRateCoeff::CPOdd,epdf->Cat("catFinalState"),epdf->Cat("obsTag"),epdf->Var("C_f"),epdf->Var("C_fbar"),epdf->Var("tageff"),epdf->Var("eta"),epdf->Var("etabar"),epdf->Var("asym_prod"),epdf->Var("asym_det"),epdf->Var("asym_tageff"));
	DecRateCoeff *coeff_s = new DecRateCoeff("coef_sin","coef_sin",DecRateCoeff::CPOdd,epdf->Cat("catFinalState"),epdf->Cat("obsTag"),epdf->Var("S_f"),epdf->Var("S_fbar"),epdf->Var("tageff"),epdf->Var("eta"),epdf->Var("etabar"),epdf->Var("asym_prod"),epdf->Var("asym_det"),epdf->Var("asym_tageff"));
	DecRateCoeff *coeff_sh = new DecRateCoeff("coef_sinh","coef_sinh",DecRateCoeff::CPEven,epdf->Cat("catFinalState"),epdf->Cat("obsTag"),epdf->Var("f1_f"),epdf->Var("f1_fbar"),epdf->Var("tageff"),epdf->Var("eta"),epdf->Var("etabar"),epdf->Var("asym_prod"),epdf->Var("asym_det"),epdf->Var("asym_tageff"));
	DecRateCoeff *coeff_ch = new DecRateCoeff("coef_cosh","coef_cosh",DecRateCoeff::CPEven,epdf->Cat("catFinalState"),epdf->Cat("obsTag"),epdf->Var("f0_f"),epdf->Var("f0_fbar"),epdf->Var("tageff"),epdf->Var("eta"),epdf->Var("etabar"),epdf->Var("asym_prod"),epdf->Var("asym_det"),epdf->Var("asym_tageff"));

	epdf->AddRealToStore(coeff_ch);
	epdf->AddRealToStore(coeff_sh);
	epdf->AddRealToStore(coeff_c);
	epdf->AddRealToStore(coeff_s);


	///////////////////Generiere PDF's/////////////////////
	//Zeit
	epdf->GaussModel("resTimeGauss",epdf->Var("obsTime"),epdf->Var("allTimeResMean"),epdf->Var("allTimeReso"));
	epdf->BDecay("pdfSigTime",epdf->Var("obsTime"),epdf->Var("tau"),epdf->Var("dgamma"),epdf->Real("coef_cosh"),epdf->Real("coef_sinh"),epdf->Real("coef_cos"),epdf->Real("coef_sin"),epdf->Var("deltaM"),epdf->Model("resTimeGauss"));
	//Zeit Untergrund
	//epdf->Exponential("pdf_bkg_mass_time", epdf->Var("obsTime"), epdf->Var("tau_bkg"));
     epdf->Extend("pdfExtend", epdf->Pdf("pdfSigTime"), epdf->Real("sig_Yield"));
	
	//PDF für SIgnal und UNtergrund der Masse
	/*epdf->Gaussian("pdf_sig_mass_gauss", epdf->Var("obsMass"), epdf->Var("pdf_sig_mass_mean"),epdf->Var"pdf_sig_mass_width"));
	epdf->Exponential("pdf_bkg_mass_expo", epdf->Var("obsMass"),epdf->Var("pdf_bkg_mass_lambda"));*/




     
     
     //Zusammenfassen der Parameter in einem RooArgSet
     RooArgSet Observables;
     Observables.add(RooArgSet( epdf->Var("obsTime"),epdf->Cat("catFinalState"),epdf->Cat("obsTag")));

     
     
  
     

     //Multipliziere Signal und Untergrund PDF mit ihrer jeweiligen Zerfalls PDF//
     //Untergrund * Zerfall
     /*epdf->Product("pdf_bkg", RooArgSet(epdf->Pdf("pdf_bkg_mass_expo"), epdf->Pdf("pdf_bkg_mass_time")));
     //Signal * Zerfall
     epdf->Product("pdf_sig", RooArgSet(epdf->Pdf("pdf_sig_mass_gauss"),epdf->Pdf("pdfSigTime")));
	//Addiere PDF's
     epdf->Add("pdf_total", RooArgSet(epdf->Pdf("pdf_sig_mass_gauss*pdf_sig_time_decay"), epdf->Pdf("pdf_bkg_mass*pdf_bkg_time_decay")), RooArgSet(epdf->Var("bkg_Yield"),epdf->Var("sig_Yield")));*/


    	    


     
     
     RooWorkspace ws;
     ws.import(epdf->Pdf("pdfExtend"));
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
     // cfg_com.PrintAll();

     // Initialize the toy factory module with the config objects and start
     // generating toy samples.
     doofit::toy::ToyFactoryStd tfac(cfg_com, cfg_tfac);
     doofit::toy::ToyStudyStd tstudy(cfg_com, cfg_tstudy);


     //for(int i=0;i<100;i++)
     //{
     RooDataSet* data = tfac.Generate();
     epdf->Pdf("pdfExtend").getParameters(data)->readFromFile("/home/chasenberg/Repository/bachelor-template/ToyStudy/dootoycp-parameter.txt");
     epdf->Pdf("pdfExtend").getParameters(data)->writeToFile("/home/chasenberg/Repository/bachelor-template/ToyStudy/dootoycp-parameter.txt.new");

     
    
     // fitting on the generated dataset is your responsibility
     RooFitResult* fit_result = epdf->Pdf("pdfExtend").fitTo(*data, RooFit::Save(true));
  

	using namespace doofit::plotting;

     PlotConfig cfg_plot("cfg_plot");
     cfg_plot.InitializeOptions();

     // plot PDF and directly specify component
     Plot myplot(cfg_plot, epdf->Var("obsTime"), *data, RooArgList(epdf->Pdf("pdfExtend")));
     myplot.PlotIt();









     //Speichern der Ergebnisse
     tstudy.StoreFitResult(fit_result);
     //}
     /*tstudy.FinishFitResultSaving();
     tstudy.ReadFitResults();
     tstudy.EvaluateFitResults();
     tstudy.PlotEvaluatedParameters();*/
    


 }

