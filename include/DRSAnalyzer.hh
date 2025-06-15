#ifndef DRSAnalyzer_HH
#define DRSAnalyzer_HH

// STD INCLUDES
#include <iostream>
#include <string>
#include <ctime>
#include <unistd.h>

// SYS includes
#include <sys/types.h>
#include <sys/stat.h>
#include <typeinfo>

// ROOT INCLUDES
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TGraph.h"
#include "TAxis.h"
#include "TText.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TVectorF.h"
#include "TMatrixF.h"
#include <TComplex.h>
#include "TDecompChol.h"
#include "TF1.h"
#include "TLeaf.h"
#include <TRegexp.h>
#include <TPRegexp.h>

// LOCAL INCLUDES
#include "Configuration.hh"

#define F_LOW 0//low frequency in GHz
#define F_HIGH 5//High frequency in GHz

// This is the base class for .dat --> .root converters.

class DRSAnalyzer {
    public:
        DRSAnalyzer(int numChannels=999, int numTimes=999, int numSamples=999, int res=1, float scale=1.0) :
        NUM_CHANNELS(numChannels), NUM_TIMES(numTimes), NUM_SAMPLES(numSamples), DAC_RESOLUTION(res), DAC_SCALE(scale),
        file(0), tree(0) 
        {
        }
        ~DRSAnalyzer(){ if (file) file->Close(); }

        void Analyze();
        void InitLoop();
        void RunEventsLoop();
        void ResetVar();
        void ResetAnalysisVariables();

        float GetPulseIntegral(std::vector<float> a, std::vector<float> t, unsigned int i_st, unsigned int i_stop);
        unsigned int GetIdxClosest(float value, std::vector<float> v, unsigned int i_st, int direction=+1);
        unsigned int GetIdxFirstCross(float value, std::vector<float> v, unsigned int i_st, int direction=+1);
        void AnalyticalPolinomialSolver(unsigned int Np, float* in_x, float* in_y, unsigned int deg, float* &out_coeff);
        float PolyEval(float x, float* coeff, unsigned int deg);
        float WSInterp(float t, int N, std::vector<float> tn, std::vector<float> cn);
        std::string split(const std::string& half, const std::string& s, const std::string& h) const;
        void GetDim(TTree* const tree, const std::string& var, unsigned int& f, unsigned int& s);
        int GetChannelsMeasurement(int i_aux);
        unsigned int GetTimeIndex(unsigned int n_ch) { return 0; }
        inline bool exists_test2(const std::string& name) { return ( access( name.c_str(), F_OK ) != -1 ); }

    // protected:

        unsigned int NUM_CHANNELS;
        unsigned int NUM_TIMES;
        unsigned int NUM_SAMPLES;
        const unsigned int DAC_RESOLUTION; // DAC resolution (2^[bit])
        const float DAC_SCALE; // [V] total scale of the DAC
        float scale_minimum = -500; // [mV] Voltage value corresponding to 0 DAC counts
        unsigned int N_warnings = 0;
        unsigned int N_warnings_to_print = 15;

        // Set by command line arguments or default
        Configuration* config = nullptr;

        TString input_file_path;
        TString output_file_path;
        unsigned long int N_evts = 0;
        unsigned int start_evt = 0;

        bool verbose = false;
        bool save_raw = false;
        bool save_meas = false;
        bool draw_debug_pulses = false;
        bool correctForTimeOffsets = false;
        TString img_format = ".png";
        long int N_evt_expected = -1;

        // Analysis variables
        std::vector<float> time;
        std::map<TString, std::vector<float>*> channelMap;
        std::vector<float> timeOffset;

        // Output tree vars
        unsigned int event_n = 0;

        // Output root file
        TFile* file;
        TTree* tree;
        //Input root files (optional, now use by ETL simulation)
        TFile* file_in;
        TTree* tree_in;
        TObjArray* branches;

        std::map<TString, float> var;
        std::vector<TString> var_names = {
          "baseline",
          "baseline_RMS",
          "noise",
          "amp",
          "t_peak",
          "integral",
          "intfull",
          "risetime",
          "decaytime"
        };
        std::vector<TString> DRS_Channel_Names;
        std::map<TString, void*> branch_buffers;
        std::vector<TString> branch_names_OnlyCopy;
        std::vector<TString> branch_names_DRS;

        // vector<int> active_ch;
};

#endif
