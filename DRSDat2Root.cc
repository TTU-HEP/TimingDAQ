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

  aux = ParseCommandLine( argc, argv, "save_meas" );
  aux.ToLower();
  if(aux == "true") a.save_meas = true;

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

  a.file = new TFile(a.output_file_path.Data(), "RECREATE");
  a.tree = new TTree("EventTree", "Digitized waveforms");
  a.file_in = new TFile(a.input_file_path,"READ");
  a.tree_in = (TTree*)a.file_in->Get("EventTree");
  
  a.GetDim(a.tree_in, "channel", a.NUM_CHANNELS, a.NUM_SAMPLES);
  a.GetDim(a.tree_in, "time", a.NUM_TIMES, a.NUM_SAMPLES);
  for(unsigned int i = 0; i < a.NUM_CHANNELS; i++){a.active_ch.emplace_back(i);}

  a.InitLoop();
  a.tree_in->SetBranchAddress("event_n", &a.event_n);
  a.tree_in->SetBranchAddress("channel", &(a.channel[0][0]));
  a.tree_in->SetBranchAddress("time", &(a.time[0][0]));
  a.tree_in->SetBranchAddress("timeoffsets", &(a.timeOffset[0]));
  
  a.tree->Branch("event_n", &a.event_n, "event_n/i");
  a.tree->Branch("timeoffsets", &(a.timeOffset[0]), Form("timeoffsets[%d]/F", a.NUM_CHANNELS));

  // // Get the list of branches
  // TObjArray* branchList = tree_in->GetListOfBranches();

  // std::vector<std::string> DRS_Channel_Names;

  // std::cout << "Branches in tree:" << std::endl;
  // for (int i = 0; i < branchList->GetEntries(); ++i) 
  // {
  //   TBranch* branch = (TBranch*)branchList->At(i);
  //   std::string brName = std::string(branch->GetName());
  //   if (brName.find("DRS") != std::string::npos && brName.find("Channel") != std::string::npos) 
  //   {
  //     DRS_Channel_Names.push_back(brName);
  //   }
  // }

  // for(const auto& n : DRS_Channel_Names)
  // {
  //   std::cout<<n<<std::endl;
  // }

  unsigned int evt_progress_print_rate = a.verbose ? 100 : 1000;
  unsigned int N_written_evts = 0;
  int n_evt_tree = a.tree_in->GetEntries();

  for(int i_aux = a.start_evt; i_aux < n_evt_tree && (a.N_evts==0 || i_aux<a.N_evts); i_aux++)
  {
    if (i_aux % 500 == 0) std::cout << "Processing Event " << i_aux << "\n" << std::endl;
    a.GetChannelsMeasurement(i_aux);
    a.Analyze();
    a.tree->Fill();
    N_written_evts++;
    a.event_n++;
  }
  a.file->Write();
  
  return 0;
}
