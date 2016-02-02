#include "doofit/config/CommonConfig.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStd.h"
#include "doofit/toy/ToyFactoryStd/ToyFactoryStdConfig.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStd.h"
#include "doofit/toy/ToyStudyStd/ToyStudyStdConfig.h"

int main(int argc, char *argv[]) {

	doofit::config::CommonConfig cfg_com("common");
   cfg_com.InitializeOptions(argc, argv);
   doofit::toy::ToyFactoryStdConfig cfg_tfac("toyfac");
   cfg_tfac.InitializeOptions(cfg_com);
   doofit::toy::ToyStudyStdConfig cfg_tstudy("toystudy");
   cfg_tstudy.InitializeOptions(cfg_tfac);

   cfg_com.CheckHelpFlagAndPrintHelp();
   cfg_com.PrintAll();

   doofit::toy::ToyStudyStd tstudy(cfg_com, cfg_tstudy);

   // do the complete automated analysis of toy fits
   tstudy.ReadFitResults();
   tstudy.EvaluateFitResults();
   tstudy.PlotEvaluatedParameters();
}