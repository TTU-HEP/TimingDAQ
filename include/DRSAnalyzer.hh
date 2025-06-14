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

// LOCAL INCLUDES
#include "Configuration.hh"

#define F_LOW 0//low frequency in GHz
#define F_HIGH 5//High frequency in GHz

// This is the base class for .dat --> .root converters.

class DRSAnalyzer {
    public:
        DRSAnalyzer(int numChannels=999, int numTimes=999, int numSamples=999, int res=1, float scale=1.0, int numFsamples=0) :
        NUM_CHANNELS(numChannels), NUM_TIMES(numTimes), NUM_SAMPLES(numSamples), DAC_RESOLUTION(res), DAC_SCALE(scale), NUM_F_SAMPLES(numFsamples),
        file(0), tree(0) 
        {
        }

        ~DRSAnalyzer()
        {
          cout << "In DRSAnalyzer destructor" << endl;
          if (file)
          {
            std::cout<<"Close File"<<std::endl;
            file->Close();
          }
          std::cout<<"Done with DRSAnalyzer destructor"<<std::endl;
        }

        int getNumChannels() { return NUM_CHANNELS; }
        int getNumTimes() { return NUM_TIMES; }
        int getNumSamples() { return NUM_SAMPLES; }

        void setNumChannels(unsigned int numChannels) {NUM_CHANNELS = numChannels;}
        void setNumTimes(unsigned int numTimes) {NUM_TIMES = numTimes;}
        void setNumSamples(unsigned int numSamples) {NUM_SAMPLES = numSamples;}

        void Analyze();
        void InitLoop();
        void RunEventsLoop();
        void ResetVar(unsigned int n_ch);
        void ResetAnalysisVariables();

        float GetPulseIntegral(float *a, float *t, unsigned int i_st, unsigned int i_stop); //returns charge in pC asssuming 50 Ohm termination
        unsigned int GetIdxClosest(float value, float* v, unsigned int i_st, int direction=+1);
        unsigned int GetIdxFirstCross(float value, float* v, unsigned int i_st, int direction=+1);
        void AnalyticalPolinomialSolver(unsigned int Np, float* in_x, float* in_y, unsigned int deg, float* &out_coeff, float* err = 0);
        float PolyEval(float x, float* coeff, unsigned int deg);
        float WSInterp(float t, int N, float* tn, float* cn);
        float FrequencySpectrum(double freq, double tMin, double tMax, int ich, int t_index);
        float FrequencySpectrum(double freq, double tMin, double tMax, unsigned int n_samples, float* my_channel, float* my_time);
        std::string split(const std::string& half, const std::string& s, const std::string& h) const;
        void GetDim(TTree* const tree, const std::string& var, unsigned int& f, unsigned int& s);
        int GetChannelsMeasurement(int i_aux);
        unsigned int GetTimeIndex(unsigned int n_ch) { return 0; }

        inline bool exists_test2(const std::string& name) 
        {
          return ( access( name.c_str(), F_OK ) != -1 );
        }

    // protected:

        unsigned int NUM_CHANNELS;
        unsigned int NUM_TIMES;
        unsigned int NUM_SAMPLES;
        const unsigned int NUM_F_SAMPLES;//Fourier samples
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
        float* AUX_time;
        float* AUX_channel;
        float* AUX_channel_spectrum;

        float** time;
        float** channel;
        float** channel_spectrum;
        float* frequency;
        float* timeOffset;

        // Output tree vars
        int event_n = 0;

        // Output root file
        TFile* file;
        TTree* tree;
        //Input root files (optional, now use by ETL simulation)
        TFile* file_in;
        TTree* tree_in;

        std::map<TString, float*> var;
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

        vector<int> active_ch;
};

#endif
