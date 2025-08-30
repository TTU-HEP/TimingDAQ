// C++ includes
#include <string>
#include <assert.h>

// ROOT includes
#include <TROOT.h>

//LOCAL INCLUDES
#include "DRSAnalyzer.hh"

TString ParseCommandLine( int argc, char* argv[], TString opt )
{
  TString out = "";
  for (int i = 1; i < argc && out==""; i++ ) 
  {
    TString tmp( argv[i] );
    if ( tmp.Contains("--"+opt) ) 
    {
      out = tmp.Contains("=") ? tmp(tmp.First("=") + 1, tmp.Length()) : TString("true");
    }            
  }
  return out;
}

void GetCommandLineArgs(int argc, char **argv, DRSAnalyzer& a) 
{
  std::cout << "-------------- TimingDAQ (DRS2Root) --------------" << std::endl;

  TString aux;
  aux = ParseCommandLine( argc, argv, "verbose" );
  aux.ToLower();
  if(aux == "true") a.verbose = true;

  a.input_file_path = ParseCommandLine( argc, argv, "input_file" );
  a.output_file_path = ParseCommandLine( argc, argv, "output_file" );

  aux = ParseCommandLine( argc, argv, "config" );
  a.config = new Configuration(aux.Data(), a.verbose);

  aux = ParseCommandLine( argc, argv, "N_evts" );
  long int N_tmp = aux.Atoi();
  a.N_evts = N_tmp;

  aux = ParseCommandLine( argc, argv, "start_evt" );
  if(aux != "") 
  {
    a.start_evt = aux.Atoi();
    std::cout << "[INFO] Starting from event: " << a.start_evt << std::endl;
  }

  aux = ParseCommandLine( argc, argv, "N_evt_expected" );
  a.N_evt_expected = aux.Atoi();

  aux = ParseCommandLine( argc, argv, "correctForTimeOffsets" );
  aux.ToLower();
  if(aux == "true") a.correctForTimeOffsets = true;

  aux = ParseCommandLine( argc, argv, "save_raw" );
  aux.ToLower();
  if(aux == "true")
  {
    a.save_raw = true;
    std::cout << "save raw" << std::endl;
  }

  aux = ParseCommandLine( argc, argv, "draw_debug_pulses" );
  aux.ToLower();
  if(aux != "false" && aux != "") 
  {
    if(aux != "true") a.img_format = aux;
    std::cout << "[INFO]: Saving debug pulses in " << a.img_format.Data() << std::endl;
    a.draw_debug_pulses =  true;
  }
}

int main(int argc, char **argv) 
{
  gROOT->SetBatch();
  DRSAnalyzer a;
  GetCommandLineArgs(argc, argv, a);

  a.file_in = new TFile(a.input_file_path,"READ");
  a.file = new TFile(a.output_file_path.Data(), "RECREATE");
  a.tree_in = (TTree*)a.file_in->Get("EventTree");
  a.tree = new TTree("EventTree", "Digitized waveforms");
  a.branches = a.tree_in->GetListOfBranches();
  TPRegexp channelRegex("^DRS_Board[0-9]+_Group[0-9]+_Channel[0-9]+$");

  for (unsigned int i = 0; i < a.branches->GetEntries(); ++i) 
  {
    auto* br = (TBranch*)a.branches->At(i);
    const auto& name = br->GetName();
    auto* leaf = br->GetLeaf(name);
    const TString& typeName = leaf->GetTypeName();

    
    if (channelRegex.Match(name))
    {
      auto* vec = new std::vector<float>;
      a.tree_in->SetBranchAddress(name, &vec);
      a.channelMap[name] = vec;
    }
    else if(TString(name).BeginsWith("FERS") && typeName.Contains("vector"))
    {
      if(typeName.Contains("unsigned short"))
      {
	      auto* vec = new std::vector<unsigned short>;
	      a.tree_in->SetBranchAddress(name, &vec);
	      a.channelFERSUSMap[name] = vec;
      }
    }
    else
    {
      a.branch_names.emplace_back(name);
    }
  }

  a.InitLoop();
  a.tree_in->GetEntry(0);  
  unsigned int evt_progress_print_rate = a.verbose ? 1 : 1000;
  unsigned int N_written_evts = 0;
  int n_evt_tree = a.tree_in->GetEntries();

  for(int i_aux = a.start_evt; i_aux < n_evt_tree && (a.N_evts==0 || i_aux<a.N_evts); i_aux++)
  {
    if (i_aux % 500 == 0) std::cout << "Processing Event " << i_aux << std::endl;
    a.GetChannelsMeasurement(i_aux);
    a.Analyze();
    a.tree->Fill();
    N_written_evts++;
    a.event_n++;
  }

  a.file->Write();

  return 0;
}
