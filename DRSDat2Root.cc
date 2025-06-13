// C++ includes
#include <string>
#include <assert.h>

// ROOT includes
#include <TROOT.h>

//LOCAL INCLUDES
#include "DRSAnalyzer.hh"

int main(int argc, char **argv) 
{
  gROOT->SetBatch();

  DRSAnalyzer analyzer;
  analyzer.GetCommandLineArgs(argc, argv);
  analyzer.InitOutput();
  analyzer.InitLoop();
  analyzer.RunEventsLoop();
  
  return 0;
}
