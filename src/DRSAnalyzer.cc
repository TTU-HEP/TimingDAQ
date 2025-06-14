#include "DRSAnalyzer.hh"

using namespace std;

void DRSAnalyzer::Analyze()
{
  /*************************************
  LOOP OVER CHANNELS
  **************************************/

  for(unsigned int i=0; i<NUM_CHANNELS; i++) 
  {
    ResetVar(i);
    if ( !config->hasChannel(i) ) 
    {
      continue;
    }
    TString name = Form("pulse_event%d_ch%d", event_n, i);
    
    //-------------------------------------------------------------------------------------
    // Get the attenuation/amplification scale factor and convert ADC counts to mV
    //IMPORTANT: polarity is taking into account here! if negative it will switch the pulse
    //-------------------------------------------------------------------------------------
    float scale_factor = (1.0 * DAC_SCALE / (float)DAC_RESOLUTION) * config->getChannelMultiplicationFactor(i);
    
    // ------- Get baseline ------
    unsigned int bl_st_idx = static_cast<unsigned int>((config->channels[i].baseline_time[0])*NUM_SAMPLES);
    unsigned int bl_en_idx = static_cast<unsigned int>((config->channels[i].baseline_time[1])*NUM_SAMPLES);
    unsigned int bl_length = bl_en_idx - bl_st_idx;

    float baseline = 0;
    for(unsigned int j = bl_st_idx; j < bl_en_idx; j++) 
    {
        baseline += channel[i][j];
    }
    
    if(bl_length <=1) cout << "WARNING: Baseline window is trivially short, probably configured incorrectly"<<endl;
    baseline /= (float) bl_length;
    TF1* f = nullptr;
    if(config->channels[i].algorithm.Contains("HNR")) 
    {
      // Perform a sin fit for the baseline
      auto gr_bl = TGraph(bl_length, &(time[GetTimeIndex(i)][bl_st_idx]), &(channel[i][bl_st_idx]));
      f = new TF1("f_bl", "[0]+[1]*sin([2]+[3]*x)");
      f->SetParameter(0, baseline);
      f->SetParameter(1, 40./scale_factor);
      f->SetParameter(2, 1.);
      f->SetParameter(3, 2*TMath::Pi()/75.);
      f->SetParNames("const", "A", "phi0", "omega");
      auto r = gr_bl.Fit(f, "SQN");
      if (draw_debug_pulses) 
      {
        cout << "Harmonic Noise Removal:\n";
        cout << Form("const = %.2f\n", f->GetParameter(0)*scale_factor);
        cout << Form("A = %.2f\n", f->GetParameter(1)*scale_factor);
        cout << Form("phi_0 = %.2f\n", f->GetParameter(2));
        cout << Form("T = %.2f\n", 2*TMath::Pi()/f->GetParameter(3));
      }

      baseline = f->GetParameter(0);
    }

    if (config->channels[i].v_baseline.size() < 200) 
    {
      config->channels[i].v_baseline.push_back(baseline);
    }
    else 
    {
      float mean = TMath::Mean(config->channels[i].v_baseline.size(), &(config->channels[i].v_baseline[0]));
      float rms = TMath::RMS(config->channels[i].v_baseline.size(), &(config->channels[i].v_baseline[0]));

      if (config->channels[i].v_baseline.size() < 800)
      {
        config->channels[i].v_baseline.push_back(baseline);
      }
    }
    var["baseline"][i] = scale_factor * baseline;
    // ------------- Get minimum position, max amplitude and scale the signal
    unsigned int idx_min = 0;
    float amp = 0;
    for(unsigned int j=0; j<NUM_SAMPLES; j++) {
      if(config->channels[i].algorithm.Contains("HNR")) 
      {
        channel[i][j] = scale_factor * (channel[i][j] - f->Eval(time[GetTimeIndex(i)][j]));//baseline subtraction
      }
      else channel[i][j] = scale_factor * (channel[i][j] - baseline);//baseline subtraction
      bool range_check = j>bl_st_idx+bl_length && j<(int)(0.9*NUM_SAMPLES);
      bool max_check = true;
      if ( config->channels[i].counter_auto_pol_switch > 0 ) 
      {
        max_check = fabs(channel[i][j]) > fabs(amp);
      }
      else 
      {
        max_check = channel[i][j] < amp;
      }

      if(( range_check && max_check) || j == bl_st_idx+bl_length) 
      {
      // if(( j>bl_st_idx+bl_length && j<(int)(0.9*NUM_SAMPLES) && fabs(channel[i][j]) > fabs(amp)) || j == bl_st_idx+bl_length) {
        idx_min = j;
        amp = channel[i][j];
      }
    }

    //std::cout << "amp: " << amp << std::endl;
    if(config->channels[i].algorithm.Contains("HNR")) {delete f;}

    //************************************************************************************
    //If the minimum point is the first sample, then the channel is bad, and we skip it.
    //************************************************************************************
    if (idx_min == 0) continue;

    var["t_peak"][i] = time[GetTimeIndex(i)][idx_min];
    var["amp"][i] = -amp;

    float baseline_RMS = 0;
    var["noise"][i] = channel[i][bl_st_idx+5];
    for(unsigned int j=bl_st_idx; j<=(bl_st_idx+bl_length); j++) 
    {
      baseline_RMS += channel[i][j]*channel[i][j];
    }
    baseline_RMS = sqrt(baseline_RMS/bl_length);
    var["baseline_RMS"][i] = baseline_RMS;

    // --------------- Define pulse graph
    float * yerr = new float[NUM_SAMPLES];
    for(unsigned j = 0; j < NUM_SAMPLES; j++) yerr[j] = 0 * var["baseline_RMS"][i];
    TGraphErrors* pulse = new TGraphErrors(NUM_SAMPLES, time[GetTimeIndex(i)], channel[i], 0, yerr);
    pulse->SetNameTitle("g_"+name, "g_"+name);

    // Variables used both by analysis and pulse drawer
    unsigned int j_90_pre = 0, j_10_pre = 0;
    unsigned int j_90_post = 0, j_10_post = 0;
    unsigned int j_area_pre = 0, j_area_post = 0;
    vector<float*> coeff_poly_fit;
    vector<pair<int, int>> poly_bounds;
    float Re_b, Re_slope;
    bool fittable = true;
    fittable *= idx_min < (int)(NUM_SAMPLES*0.8);
    fittable *= fabs(amp) > 3 * baseline_RMS;
    fittable *= fabs(channel[i][idx_min+1]) > 2*baseline_RMS;
    fittable *= fabs(channel[i][idx_min-1]) > 2*baseline_RMS;
    fittable *= fabs(channel[i][idx_min+2]) > 2*baseline_RMS;
    fittable *= fabs(channel[i][idx_min-2]) > 2*baseline_RMS;
    // fittable *= fabs(channel[i][idx_min+3]) > 2*baseline_RMS;
    // fittable *= fabs(channel[i][idx_min-3]) > 2*baseline_RMS;
    
    //-------------------------------------------------------------
    //If pulse is still positive change it automatically
    //No sure this is a good idea (CP)
    //-------------------------------------------------------------
    if( fittable  && !config->channels[i].algorithm.Contains("None")) 
    {
      if( var["amp"][i] < 0 && config->channels[i].counter_auto_pol_switch > 0 ) 
      {
        config->channels[i].polarity *= -1;
        amp = -amp;
        var["amp"][i] = -var["amp"][i];
        scale_factor = -scale_factor;
        var["baseline"][i] = -var["baseline"][i];
        for(unsigned int j=0; j<NUM_SAMPLES; j++) 
        {
          channel[i][j] = -channel[i][j];
        }
        delete pulse;
        pulse = new TGraphErrors(NUM_SAMPLES, time[GetTimeIndex(i)], channel[i], 0, yerr);
        pulse->SetNameTitle("g_"+name, "g_"+name);

        if ( config->channels[i].counter_auto_pol_switch == 10 ) 
        {
          cout << "[WARNING] Channel " << i << ": automatic polarity switched more than 10 times" << endl;
          cout << "[WARNING] Channel " << i << ": gonna keep inverting it for you. Better check your pulse polarity!!" << endl;
        }
        config->channels[i].counter_auto_pol_switch ++;
      }
    
      // if( fittable  && !config->channels[i].algorithm.Contains("None")) {
      /************************************
       //Get 10% of the amplitude crossings
       ************************************
       */
      j_10_pre = GetIdxFirstCross(amp*0.1, channel[i], idx_min, -1);
      j_10_post = GetIdxFirstCross(amp*0.1, channel[i], idx_min, +1);

      // -------------- Integrate the pulse
      j_area_pre = GetIdxFirstCross(amp*0.05, channel[i], idx_min, -1);
      j_area_post = GetIdxFirstCross(0. , channel[i], idx_min, +1);
      var["integral"][i] = GetPulseIntegral(channel[i], time[GetTimeIndex(i)], j_area_pre, j_area_post);
      var["intfull"][i] = GetPulseIntegral(channel[i], time[GetTimeIndex(i)], 5, NUM_SAMPLES-5);

      // -------------- Compute rise and falling time
      float* coeff;

      j_90_pre = GetIdxFirstCross(amp*0.9, channel[i], j_10_pre, +1);
      AnalyticalPolinomialSolver(j_90_pre-j_10_pre+1, &(time[GetTimeIndex(i)][j_10_pre]), &(channel[i][j_10_pre]), 1, coeff);
      var["risetime"][i] = abs(0.8*var["amp"][i] / coeff[1]);
      delete [] coeff;

      j_90_post = GetIdxFirstCross(amp*0.9, channel[i], j_10_post, -1);
      AnalyticalPolinomialSolver(j_10_post-j_90_post+1, &(time[GetTimeIndex(i)][j_90_post]), &(channel[i][j_90_post]), 1, coeff);
      var["decaytime"][i] = coeff[1];
      delete [] coeff;

      /************************************
      // -------------- Global Time-offsets
      *************************************/
      double myTimeOffset = 0;
      if (correctForTimeOffsets) 
      {
        myTimeOffset = timeOffset[i];
      }

      /************************************
      // -------------- Do the gaussian fit
      *************************************/
      if( config->channels[i].algorithm.Contains("G") ) 
      {
        float frac = config->channels[i].gaus_fraction;
        unsigned int j_down = GetIdxFirstCross(amp*frac, channel[i], idx_min, -1);
        unsigned int j_up = GetIdxFirstCross(amp*frac, channel[i], idx_min, +1);
        if( j_up - j_down < 4 ) 
        {
          j_up = idx_min + 1;
          j_down = idx_min - 1;
        }

        TF1* fpeak = new TF1("fpeak"+name, "gaus", time[GetTimeIndex(i)][j_down], time[GetTimeIndex(i)][j_up]);

        float ext_sigma = time[GetTimeIndex(i)][j_up] - time[GetTimeIndex(i)][j_down];
        if (amp*frac < -baseline_RMS) ext_sigma *= 0.25;
        fpeak->SetParameter(0, amp * sqrt(2*3.14) * ext_sigma );
        fpeak->SetParameter(1, time[GetTimeIndex(i)][idx_min]);
        fpeak->SetParameter(2, ext_sigma);
        fpeak->SetLineColor(kBlue);

        TString opt = "R";
        if ( draw_debug_pulses ) opt += "+";
        else opt += "QN0";
        pulse->Fit("fpeak"+name, opt);

        var["gaus_mean"][i] = fpeak->GetParameter(1) + myTimeOffset;
        var["gaus_sigma"][i] = fpeak->GetParameter(2);
        var["gaus_chi2"][i] = pulse->Chisquare(fpeak, "R");

        delete fpeak;
      }
      /*******************************
      // -------------- Do  linear fit
      ********************************/
      if(config->channels[i].algorithm.Contains("Re") ) 
      {
        unsigned int i_min = GetIdxFirstCross(config->channels[i].re_bounds[0]*amp, channel[i], idx_min, -1); // GetIdxClosest(min(channel[i]), channel[i], idx_min, -1);
        unsigned int i_max = GetIdxFirstCross(config->channels[i].re_bounds[1]*amp, channel[i], i_min  , +1);
        float t_min = time[GetTimeIndex(i)][i_min];
        float t_max = time[GetTimeIndex(i)][i_max];

        TF1* flinear = new TF1("flinear"+name, "[0]*x+[1]", t_min, t_max);
        flinear->SetLineColor(2);

        TString opt = "R";
        if ( draw_debug_pulses ) opt += "+";
        else opt += "QN0";
        pulse->Fit("flinear"+name, opt);
        Re_slope = flinear->GetParameter(0);
        Re_b     = flinear->GetParameter(1);

        for ( auto f : config->constant_fraction ) 
        {
          var[Form("linear_RE_%d", (int)(100*f))][i] = (f*amp-Re_b)/Re_slope + myTimeOffset;
        }

        for ( auto thr : config->constant_threshold ) 
        {
            var[Form("linear_RE__%dmV", (int)(fabs(thr)))][i] = (thr-Re_b)/Re_slope + myTimeOffset;
        }

        delete flinear;
      }
      /*************************************
      // -------------- Local polinomial fit
      **************************************/
      if ( config->constant_fraction.size() ) 
      {
        float start_level =  - 3 * baseline_RMS;
        unsigned int j_start =  GetIdxFirstCross( start_level, channel[i], idx_min, -1);

        for(auto f : config->constant_fraction) 
        {
          unsigned int j_st = j_start;
          if ( amp*f > start_level ) {
            if ( amp*f > -baseline_RMS && verbose) 
            {
              if(N_warnings< N_warnings_to_print) 
              {
                N_warnings++;
                cout << Form("[WARNING] ev:%d ch:%d - fraction %.2f below noise RMS", event_n, i, f) << endl;
              }
              else if (N_warnings_to_print == N_warnings) 
              {
                N_warnings++;
                cout << "[WARNING] Max number of warnings passed. No more warnings will be printed." << endl;;
              }
            }
            j_st =  GetIdxFirstCross( amp*f, channel[i], idx_min, -1);
          }

          unsigned int j_close = GetIdxFirstCross(amp*f, channel[i], j_st, +1);
          if ( fabs(channel[i][j_close-1] - f*amp) < fabs(channel[i][j_close] - f*amp) ) j_close--;

          for(auto n : config->channels[i].PL_deg) 
          {
            unsigned int span_j = (int) (min( j_90_pre-j_close , j_close-j_st)/1.5);

            if (j_90_pre - j_10_pre <= 3*n) 
            {
              span_j = max((unsigned int)(n*0.5), span_j);
              span_j = max((unsigned int)1, span_j);
            }
            else 
            {
              span_j = max((unsigned int) n, span_j);
            }

            if( j_close < span_j || j_close + span_j >= NUM_SAMPLES ) 
            {
              cout << Form("[WARNING] evt %d ch %d:  Short span around the closest point. Analytical fit not performed.", event_n, i) << endl;
              continue;
            }

            float* coeff;
            int N_add = 1;
            if (span_j + N_add + j_close < j_90_pre) 
            {
              N_add++;
            }
            AnalyticalPolinomialSolver( 2*span_j + N_add , &(channel[i][j_close - span_j]), &(time[GetTimeIndex(i)][j_close - span_j]), n, coeff);

            var[Form("LP%d_%d", n, (int)(100*f))][i] = PolyEval(f*amp, coeff, n) + myTimeOffset;

            if(draw_debug_pulses) 
            {
              coeff_poly_fit.push_back(coeff);
              poly_bounds.push_back(pair<int,int>(j_close-span_j, j_close+span_j+N_add-1));
            }
            else delete [] coeff;
          }
        }
      }
      if ( config->constant_threshold.size() ) 
      {
        float start_level =  - 3 * baseline_RMS;
        unsigned int j_start =  GetIdxFirstCross( start_level, channel[i], idx_min, -1);

        for(auto thr : config->constant_threshold) 
        {
          if (thr < amp ) continue;
          unsigned int j_st = j_start;
          if ( thr > start_level ) 
          {
            if (thr > -baseline_RMS && verbose) 
            {
              if(N_warnings< N_warnings_to_print) 
              {
                N_warnings++;
                cout << Form("[WARNING] ev:%d ch:%d - thr %.2f mV below noise RMS", event_n, i, thr) << endl;
              }
              else if (N_warnings_to_print == N_warnings) 
              {
                N_warnings++;
                cout << "[WARNING] Max number of warnings passed. No more warnings will be printed." << endl;;
              }
            }
            j_st =  GetIdxFirstCross( thr, channel[i], idx_min, -1);
          }

          unsigned int j_close = GetIdxFirstCross(thr, channel[i], j_st, +1);
          if ( fabs(channel[i][j_close-1] - thr) < fabs(channel[i][j_close] - thr) ) j_close--;

          for(auto n : config->channels[i].PL_deg) 
          {
            unsigned int span_j = (int) (min( j_90_pre-j_close , j_close-j_st)/1.5);

            if (j_90_pre - j_10_pre <= 3*n) 
            {
              span_j = max((unsigned int)(n*0.5), span_j);
              span_j = max((unsigned int)1, span_j);
            }
            else 
            {
              span_j = max((unsigned int) n, span_j);
            }

            if( j_close < span_j || j_close + span_j >= NUM_SAMPLES ) 
            {
              if (verbose) 
              {
                cout << Form("[WARNING] evt %d ch %d:  Short span around the closest point. Analytical fit not performed.", event_n, i) << endl;
                cout << j_close << "  " << span_j << "  " << thr << endl;
              }
              continue;
            }

            float* coeff;
            int N_add = 1;
            if (span_j + N_add + j_close < j_90_pre) 
            {
              N_add++;
            }
            AnalyticalPolinomialSolver( 2*span_j + N_add , &(channel[i][j_close - span_j]), &(time[GetTimeIndex(i)][j_close - span_j]), n, coeff);

            var[Form("LP%d_%dmV", n, (int)(fabs(thr)))][i] = PolyEval(thr, coeff, n) + myTimeOffset;

            if(draw_debug_pulses) 
            {
              coeff_poly_fit.push_back(coeff);
              poly_bounds.push_back(pair<int,int>(j_close-span_j, j_close+span_j+N_add-1));
            }
            else delete [] coeff;
          }
        }
      }
    }

    /*********************************************
    // ===================  Draw plot of the pulse
    **********************************************/
    if(draw_debug_pulses) 
    {
      cout << "========= Event: " << event_n << " - ch: " << i << endl;

      TCanvas* c =  new TCanvas("c_"+name, "c_"+name, 1600, 600);
      c->Divide(2);

      TLine* line = new TLine();

      // ---------- All range plot
      c->cd(1);
      c->SetGrid();
      // Draw pulse
      pulse->SetMarkerStyle(4);
      pulse->SetMarkerSize(0.5);
      pulse->GetYaxis()->SetTitle("Amplitude [mV]");
      pulse->GetXaxis()->SetTitle("Time [ns]");
      pulse->Draw("APE1");
      // Draw baseline
      line->SetLineWidth(1);
      line->SetLineColor(46);
      line->SetLineStyle(7);
      line->DrawLine(time[GetTimeIndex(i)][0], 0, time[GetTimeIndex(i)][NUM_SAMPLES-1], 0);
      line->SetLineStyle(1);
      line->DrawLine(time[GetTimeIndex(i)][bl_st_idx], 0, time[GetTimeIndex(i)][bl_st_idx+bl_length], 0);
      line->SetLineColor(47);
      line->DrawLine(time[GetTimeIndex(i)][0], var["baseline_RMS"][i], time[GetTimeIndex(i)][NUM_SAMPLES-1], var["baseline_RMS"][i]);
      line->DrawLine(time[GetTimeIndex(i)][0], -var["baseline_RMS"][i], time[GetTimeIndex(i)][NUM_SAMPLES-1], -var["baseline_RMS"][i]);

      // Draw peak
      line->SetLineColor(8);
      line->SetLineStyle(4);
      line->DrawLine(time[GetTimeIndex(i)][0], amp, var["t_peak"][i], amp);
      line->DrawLine(var["t_peak"][i], 0, var["t_peak"][i], amp);

      // Draw 10% and 90% lines;
      TLine* line_lvs = new TLine();
      line_lvs->SetLineWidth(1);
      line_lvs->SetLineColor(4);
      line_lvs->DrawLine(time[GetTimeIndex(i)][0], 0.1*amp, time[GetTimeIndex(i)][NUM_SAMPLES-1], 0.1*amp);
      line_lvs->DrawLine(time[GetTimeIndex(i)][0], 0.9*amp, time[GetTimeIndex(i)][NUM_SAMPLES-1], 0.9*amp);
      // Draw constant fractions lines
      line_lvs->SetLineColor(38);
      line_lvs->SetLineStyle(10);
      for(auto f : config->constant_fraction) 
      {
        line_lvs->DrawLine(time[GetTimeIndex(i)][0], f*amp, time[GetTimeIndex(i)][NUM_SAMPLES-1], f*amp);
      }
      // Draw constant threshold lines
      line_lvs->SetLineColor(28);
      for(auto thr : config->constant_threshold) 
      {
        line_lvs->DrawLine(time[GetTimeIndex(i)][0], thr, time[GetTimeIndex(i)][NUM_SAMPLES-1], thr);
      }


      // Draw integral area
      int N_tot_integral = j_area_post-j_area_pre;
      if( N_tot_integral > 0 ) 
      {
        vector<float> aux_time = {time[GetTimeIndex(i)][j_area_pre]};
        vector<float> aux_volt = {0};
        for(unsigned int j = j_area_pre; j <= j_area_post; j++) 
        {
          aux_time.push_back(time[GetTimeIndex(i)][j]);
          aux_volt.push_back(channel[i][j]);
        }
        aux_time.push_back(time[GetTimeIndex(i)][j_area_post]);
        aux_volt.push_back(0);
        TGraph * integral_pulse = new TGraph(aux_time.size(), &(aux_time[0]), &(aux_volt[0]));
        integral_pulse->SetFillColor(40);
        integral_pulse->SetFillStyle(3144);
        integral_pulse->Draw("FC");
        TText* t_int = new TText(var["t_peak"][i]+3, amp, Form("Integral Pulse = %1.2f,  Integral Full =%1.2f ", -var["integral"][i], -var["intfull"][i]));
        t_int->SetTextAlign(kHAlignLeft+kVAlignBottom);
        t_int->Draw();

        TText* aNt = new TText(var["t_peak"][i]+3, amp+6, Form("LP2_50 = %1.2f, t peak = %1.2f,  amp =%1.2f ", var["LP2_50"][i], var["t_peak"][i], amp));
        aNt->SetTextAlign(kHAlignLeft+kVAlignBottom);
        aNt->Draw();

        // Draw 90% and 10% pre and post points
        TGraph* gr_pre_post = new TGraph(4);
        gr_pre_post->SetPoint(0, time[GetTimeIndex(i)][j_10_pre], channel[i][j_10_pre]);
        gr_pre_post->SetPoint(1, time[GetTimeIndex(i)][j_90_pre], channel[i][j_90_pre]);
        gr_pre_post->SetPoint(2, time[GetTimeIndex(i)][j_10_post], channel[i][j_10_post]);
        gr_pre_post->SetPoint(3, time[GetTimeIndex(i)][j_90_post], channel[i][j_90_post]);
        gr_pre_post->SetMarkerColor(4);
        gr_pre_post->Draw("P*");


        // ---------- Rising edge only inverted!! -----
        c->cd(2);
        c->SetGrid();

        if( config->channels[i].algorithm.Contains("Re") ) 
        {
          unsigned int i_min = GetIdxFirstCross(config->channels[i].re_bounds[0]*amp, channel[i], idx_min, -1);
          unsigned int i_max = GetIdxFirstCross(config->channels[i].re_bounds[1]*amp, channel[i], i_min, +1);
          float y[2], x[2];
          x[0] = channel[i][i_min];
          x[1] = channel[i][i_max];
          y[0] = (channel[i][i_min] - Re_b)/Re_slope;
          y[1] = (channel[i][i_max] - Re_b)/Re_slope;

          TGraph* gr_Re = new TGraph(2, x, y);

          gr_Re->SetLineColor(46);
          gr_Re->SetLineWidth(1);
          gr_Re->SetLineStyle(7);
          gr_Re->Draw("CP");
        }

        unsigned int j_begin = j_10_pre - 6;
        unsigned int j_span = j_90_pre - j_10_pre + 10;
        //if(!(j_begin >= 0 && j_begin < sizeof(time) / sizeof(time[0]))) continue;
        TGraphErrors* inv_pulse = new TGraphErrors(j_span, &(channel[i][j_begin]), &(time[GetTimeIndex(i)][j_begin]), yerr);
        inv_pulse->SetNameTitle("g_inv"+name, "g_inv"+name);
        inv_pulse->SetMarkerStyle(5);
        inv_pulse->GetXaxis()->SetTitle("Amplitude [mV]");
        inv_pulse->GetYaxis()->SetTitle("Time [ns]");
        inv_pulse->Draw("APE1");

        vector<float> t_WS;
        vector<float> c_WS;
        float overstep = 6;
        for(unsigned int jj = j_begin; jj < j_begin + j_span; jj++) 
        {
          for(unsigned int kk = 0; kk < overstep; kk++) 
          {
            float tt = time[GetTimeIndex(i)][jj] + (time[GetTimeIndex(i)][jj+1] - time[GetTimeIndex(i)][jj]) * kk/overstep;
            float cc = WSInterp(tt, NUM_SAMPLES, time[GetTimeIndex(i)], channel[i]);

            t_WS.push_back(tt);
            c_WS.push_back(cc);
          }
        }
        TGraph* inv_pulse_WS = new TGraph(t_WS.size(), &(c_WS[0]), &(t_WS[0]));
        inv_pulse_WS->SetMarkerStyle(7);
        inv_pulse_WS->SetMarkerColor(2);
        // inv_pulse_WS->Draw("P");

        TGraph* gr_inv_pre_post = new TGraph(2);
        gr_inv_pre_post->SetPoint(0, channel[i][j_10_pre], time[GetTimeIndex(i)][j_10_pre]);
        gr_inv_pre_post->SetPoint(1, channel[i][j_90_pre], time[GetTimeIndex(i)][j_90_pre]);
        gr_inv_pre_post->SetMarkerColor(4);
        gr_inv_pre_post->Draw("P*");

        // -------------- If exist, draw local polinomial fit
        unsigned int count = 0;
        vector<int> frac_colors = {2, 6, 8, 5, 40, 46, 4, 9, 12};
        while(frac_colors.size() < config->constant_fraction.size() + config->constant_threshold.size()) 
        {
          frac_colors.push_back(2);
        }

        for( unsigned int kk = 0; kk < config->constant_fraction.size(); kk++) 
        {
          float f = config->constant_fraction[kk];
          line_lvs->SetLineColor(frac_colors[kk]);
          line_lvs->DrawLine(amp*f, time[GetTimeIndex(i)][j_begin], amp*f, time[GetTimeIndex(i)][j_90_pre + 3]);
          for(auto n : config->channels[i].PL_deg) 
          {
            vector<float> polyval;
            for(unsigned int j = poly_bounds[count].first; j <= poly_bounds[count].second; j++) 
            {
              polyval.push_back(PolyEval(channel[i][j], coeff_poly_fit[count], n));
            }

            TGraph* g_poly = new TGraph(polyval.size(), &(channel[i][poly_bounds[count].first]), &(polyval[0]) );
            g_poly->SetLineColor(frac_colors[kk]);
            g_poly->SetLineWidth(2);
            g_poly->SetLineStyle(7);
            g_poly->Draw("C");

            count++;
          }
        }

        for( unsigned int kk = 0; kk < config->constant_threshold.size(); kk++) 
        {
          float thr = config->constant_threshold[kk];
          if (thr < amp ) continue;
          line_lvs->SetLineColor(frac_colors[kk + config->constant_fraction.size()]);
          line_lvs->DrawLine(thr, time[GetTimeIndex(i)][j_begin], thr, time[GetTimeIndex(i)][j_90_pre + 3]);
          for(auto n : config->channels[i].PL_deg) 
          {
            vector<float> polyval;
            for(unsigned int j = poly_bounds[count].first; j <= poly_bounds[count].second; j++) 
            {
              polyval.push_back(PolyEval(channel[i][j], coeff_poly_fit[count], n));
            }

            TGraph* g_poly = new TGraph(polyval.size(), &(channel[i][poly_bounds[count].first]), &(polyval[0]) );
            g_poly->SetLineColor(frac_colors[kk + config->constant_fraction.size()]);
            g_poly->SetLineWidth(2);
            g_poly->SetLineStyle(7);
            g_poly->Draw("C");

            count++;
          }
        }
      }

      c->SetGrid();
      c->SaveAs("./pulses_imgs/"+name+img_format);
      delete c;
    }

    delete [] yerr;
    delete pulse;
  }
}

void DRSAnalyzer::InitLoop() 
{
    /*
    ************************
    Define PULSHAPES
    ************************
    */
    std::cout << "Define pulse shapes" << std::endl;
    AUX_time = new float[NUM_TIMES*NUM_SAMPLES];
    AUX_channel = new float[NUM_CHANNELS*NUM_SAMPLES];

    time    = new float*[NUM_TIMES];
    channel = new float*[NUM_CHANNELS];
    timeOffset = new float[NUM_CHANNELS];

    if ( NUM_F_SAMPLES > 0 )
    {
      AUX_channel_spectrum = new float[NUM_CHANNELS*NUM_F_SAMPLES];
      channel_spectrum     = new float*[NUM_F_SAMPLES];
      frequency = new float[NUM_F_SAMPLES];
    }

    for(unsigned int i=0; i<NUM_CHANNELS; i++)
    {
      channel[i] = &(AUX_channel[i*NUM_SAMPLES]);
      if(i<NUM_TIMES) time[i] = &(AUX_time[i*NUM_SAMPLES]);
      if ( NUM_F_SAMPLES ) channel_spectrum[i] = &(AUX_channel_spectrum[i*NUM_F_SAMPLES]);
    }

    if ( NUM_F_SAMPLES > 0 )
    {
      float f_step = (F_HIGH-F_LOW)/float(NUM_F_SAMPLES);
      for (unsigned int i = 0; i < NUM_F_SAMPLES; i++)
      {
        frequency[i] = F_LOW + f_step*float(i);
      }
    }

    if ( verbose ) 
    {
      cout << "NUM_CHANNELS: " << NUM_CHANNELS << endl;
      cout << "NUM_TIMES: " << NUM_TIMES << endl;
      cout << "NUM_SAMPLES: " << NUM_SAMPLES << endl;
      cout << "NUM_F_SAMPLES: " << NUM_F_SAMPLES << endl;
      cout << "DAC_SCALE : " << DAC_SCALE << "\n";
      cout << "DAC_RESOLUTION : " << DAC_RESOLUTION << "\n";
    }

    /*
    ************************
    SAVE PULSHAPES
    ************************
    */
    if(save_meas)
    {
      tree->Branch("channel", &(channel[0][0]), Form("channel[%d][%d]/F", NUM_CHANNELS, NUM_SAMPLES));
      tree->Branch("time", &(time[0][0]), Form("time[%d][%d]/F", NUM_TIMES, NUM_SAMPLES));
      if ( NUM_F_SAMPLES > 0 )
      {
        tree->Branch("channel_spectrum", &(channel_spectrum[0][0]), Form("channel_spectrum[%d][%d]/F", NUM_CHANNELS, NUM_F_SAMPLES));
        tree->Branch("frequency", &(frequency[0]), Form("frequency[%d]/F", NUM_F_SAMPLES));
      }
    }

    /*
    **********************************
    Obtain Algorithms from config file
    **********************************
    */
    bool at_least_1_gaus_fit = false;
    bool at_least_1_rising_edge = false;
    int at_least_1_LP[3] = {false};
    for(auto c : config->channels) 
    {
      if( c.second.algorithm.Contains("G")) at_least_1_gaus_fit = true;
      if( c.second.algorithm.Contains("Re")) at_least_1_rising_edge = true;
      if( c.second.algorithm.Contains("LP1")) at_least_1_LP[0] = true;
      if( c.second.algorithm.Contains("LP2")) at_least_1_LP[1] = true;
      if( c.second.algorithm.Contains("LP3")) at_least_1_LP[2] = true;
    }

    for(unsigned int i = 0; i< 3; i++) 
    {
      if(at_least_1_LP[i]) 
      {
        for (auto f : config->constant_fraction) 
        {
          var_names.push_back(Form("LP%d_%d", i+1, (int)(100*f)));
        }
        for (auto thr : config->constant_threshold) 
        {
          var_names.push_back(Form("LP%d_%dmV", i+1, (int)(fabs(thr))));
        }
      }
    }
    /*
    ************
    Gaussian Fit
    ************
    */
    if( at_least_1_gaus_fit ) 
    {
      var_names.push_back("gaus_mean");
      var_names.push_back("gaus_sigma");
      var_names.push_back("gaus_chi2");
    }
    /*
    **********************
    Linear Rising Edge Fit(good old stuff!)
    **********************
    */
    if( at_least_1_rising_edge ) 
    {
      for (auto f : config->constant_fraction) {
        var_names.push_back(Form("linear_RE_%d", (int)(100*f)));
      }

      for (auto thr : config->constant_threshold) {
        var_names.push_back(Form("linear_RE__%dmV", (int)(fabs(thr))));
      }
    }

    /*
    *******************************************************
    // Create the tree branches an the associated variables
    *******************************************************
    */
    if ( verbose ) { cout << "Initializing all tree variables" << endl; }
    if (save_meas && verbose) 
    {
      cout << "   channel\n   time" << endl;
    }
    for(TString n : var_names)
    {
      var[n] = new float[NUM_CHANNELS];
      tree->Branch(n, &(var[n][0]), n+Form("[%d]/F", NUM_CHANNELS));
      if( verbose ) { cout << "   " << n.Data() << endl; }
    }

    for(unsigned int i = 0; i < NUM_CHANNELS; i++) ResetVar(i);
};

int DRSAnalyzer::GetChannelsMeasurement(int i_aux) 
{
  ResetAnalysisVariables();
  tree_in->GetEntry(i_aux);
  return 0;
}

void DRSAnalyzer::ResetVar(unsigned int n_ch) 
{
  for(auto n: var_names) 
  {
    var[n][n_ch] = 0;
  }
}

void DRSAnalyzer::ResetAnalysisVariables() 
{
  for(unsigned int i=0; i<NUM_CHANNELS; i++) 
  {
    for(unsigned int j=0; j<NUM_SAMPLES; j++) 
    {
      channel[i][j] = 0;
      if(i < NUM_TIMES) time[i][j] = 0;
    }
  }
}

float DRSAnalyzer::GetPulseIntegral(float *a, float *t, unsigned int i_st, unsigned int i_stop) //returns charge in pC asssuming 50 Ohm termination
{
  //Simpson's Rule for equaled space with Cartwright correction for unequaled space
  float integral = 0.;
  for (unsigned int i=i_st; i < i_stop-2 ; i+=2) 
  {
    float aux = ( 2-(t[i+2]-t[i+1])/(t[i+1]-t[i]) ) * a[i];
    aux += (t[i+2]-t[i])*(t[i+2]-t[i])/((t[i+2]-t[i+1])*(t[i+1]-t[i])) * a[i+1];
    aux += ( 2-(t[i+1]-t[i])/(t[i+2]-t[i+1]) ) * a[i+2];

    integral += ( (t[i+2]-t[i]) / 6.0 ) * aux;
  }

  integral *= -1.0;//1e-9 * 1e-3 * (1.0/50.0) * 1e12; //in units of pC, for 50 [Ohms] termination
  return integral;
}

unsigned int DRSAnalyzer::GetIdxClosest(float value, float* v, unsigned int i_st, int direction) 
{
  unsigned int idx_end = direction>0 ? NUM_SAMPLES-1 : 0;
  unsigned int i = i_st;
  unsigned int i_min = i_st;
  float min_distance = fabs(value - v[i]);

  while( i != idx_end ) 
  {
    i += direction;
    float d = fabs(value - v[i]);
    if( d < min_distance ) 
    {
      min_distance = d;
      i_min = i;
    }
  }

  return i_min;
}

unsigned int DRSAnalyzer::GetIdxFirstCross(float value, float* v, unsigned int i_st, int direction) 
{
  unsigned int idx_end = direction>0 ? NUM_SAMPLES-1 : 0;
  bool rising = value > v[i_st]? true : false; // Check if the given max value is greater than the given waveform amplitude at a given start time (i_st).

  // if it is: the value of the waveform needs to rise and otherwise it doesn't need so the function will return the i_st.
  unsigned int i = i_st;
  while( (i != idx_end)) 
  { // loop over the time bins
    if(rising && v[i] > value) break; // check if the rising variable is true and if the waveform value is going higher than the given amplitude value. And stops the loop.
    else if( !rising && v[i] < value) break; // check if the rising variable is false and if the waveform value is still lower than the given amplitude value. And stops the loop.
    i += direction; // Otherwise it moves to the next time bin.
  }

  return i;
}

void DRSAnalyzer::AnalyticalPolinomialSolver(unsigned int Np, float* in_x, float* in_y, unsigned int deg, float* &out_coeff, float* err) 
{
  if(deg <= 0 || deg>3) { cout << "[ERROR]: You don't need AnalyticalPolinomialSolver for this" << endl; exit(0);}
  if(Np < deg+1) return;

  TVectorF x, x2, x3;
  x.Use(Np, in_x);

  TVectorF y;
  y.Use(Np, in_y);

  TMatrixF A(Np, deg+1);

  TMatrixFColumn(A, 0) = 1.;
  TMatrixFColumn(A, 1) = x;

  float *in_x2 = new float[Np];
  float *in_x3 = new float[Np];
  if( deg >= 2 ) 
  {
    for(unsigned int i = 0; i < Np; i++) in_x2[i] = in_x[i]*in_x[i];
    x2.Use(Np, in_x2);
    TMatrixFColumn(A, 2) = x2;
  }
  if( deg >= 3 ) 
  {
    for(unsigned int i = 0; i < Np; i++) in_x3[i] = in_x2[i] * in_x[i];
    x3.Use(Np, in_x3);
    TMatrixFColumn(A, 3) = x3;
  }

  const TVectorD c_norm = NormalEqn(A,y);

  out_coeff = new float[deg+1];
  for(unsigned int i = 0; i<= deg; i++) {
    out_coeff[i] = c_norm[i];
  }

  delete [] in_x2;
  delete [] in_x3;
  return;
}

float DRSAnalyzer::PolyEval(float x, float* coeff, unsigned int deg) 
{
  float out = coeff[0] + x*coeff[1];
  for(unsigned int i=2; i<=deg; i++) 
  {
    out += coeff[i]*pow(x, i);
  }
  return out;
}

float DRSAnalyzer::WSInterp(float t, int N, float* tn, float* cn) 
{
  float out = 0;
  float dt = (tn[0] - tn[N-1])/N;
  for(unsigned i = 0; i < N; i++) 
  {
    float x = (t - tn[i])/dt;
    out += cn[i] * sin(3.14159265358 * x) / (3.14159265358 * x);
  }
  return out;
};

float DRSAnalyzer::FrequencySpectrum(double freq, double tMin, double tMax, int ich, int t_index)
{
	const int range = 0; // extension of samples to be used beyond [tMin, tMax]
	double deltaT = (time[t_index][NUM_SAMPLES - 1] - time[t_index][0])/(double)(NUM_SAMPLES - 1); // sampling time interval
	double fCut = 0.5/deltaT; // cut frequency = 0.5 * sampling frequency from WST
	int n_min = floor(tMin/deltaT) - range; // first sample to use
	int n_max = ceil(tMax/deltaT) + range; // last sample to use
	n_min = std::max(n_min,0); // check low limit
	n_max = std::min(n_max, (int)NUM_SAMPLES - 1); // check high limit
	int n_0 = (n_min + n_max)/2;

	TComplex s(0.,0.); // Fourier transform at freq
	TComplex I(0.,1.); // i

	for(int n = n_min; n <= n_max; n++)
	{
		s += deltaT*(double)channel[ich][n]*TComplex::Exp(-I*(2.*TMath::Pi()*freq*(n-n_0)*deltaT));//maybe don't need n_0 here, I think it will just add a phase to the fourier transform
	}
  return s.Rho();
};

float DRSAnalyzer::FrequencySpectrum(double freq, double tMin, double tMax, unsigned int n_samples, float* my_channel, float* my_time)
{
  const int range = 0;
	double deltaT = (my_time[n_samples - 1] - my_time[0])/(double)(n_samples - 1);
	double fCut = 0.5/deltaT;
	int n_min = floor(tMin/deltaT) - range;
	int n_max = ceil(tMax/deltaT) + range;
	n_min = std::max(n_min,0);
	n_max = std::min(n_max, (int)(n_samples - 1));
	int n_0 = (n_min + n_max)/2;

	TComplex s(0.,0.);
	TComplex I(0.,1.);

	for(int n = n_min; n <= n_max; n++)
	{
    s += deltaT*(double)my_channel[n]*TComplex::Exp(-I*(2.*TMath::Pi()*freq*(n-n_0)*deltaT));
	}
  return s.Rho();
};

std::string DRSAnalyzer::split(const std::string& half, const std::string& s, const std::string& h) const
{
  if(s.find(h) != std::string::npos)
  {
    std::string token;
    if      ("first"==half) token = s.substr(0, s.find(h));
    else if ("last" ==half) token = s.substr(s.find(h) + h.length(), std::string::npos);
    return token;
  }
  else
  {
    return s;
  }
}

void DRSAnalyzer::GetDim(TTree* const tree, const std::string& var, unsigned int& f, unsigned int& s)
{
  TBranch* branch = tree->GetBranch(var.c_str());
  TObjArray *lol = branch->GetListOfLeaves();
  TLeaf *leaf = (TLeaf*)lol->UncheckedAt(0);
  std::string title = leaf->GetTitle();
  std::string firstdim  = split("last", split("first", title, "]"), "[");
  std::string seconddim = split("first", split("last", title, "]["), "]");
  f = static_cast<unsigned int>(std::atoi(firstdim.c_str()));
  s = static_cast<unsigned int>(std::atoi(seconddim.c_str()));
}
