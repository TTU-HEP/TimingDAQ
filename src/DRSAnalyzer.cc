#include "DRSAnalyzer.hh"

void DRSAnalyzer::Analyze()
{
  ResetVar();

  /*************************************
  LOOP OVER CHANNELS
  **************************************/
  unsigned int i = 0;
  for(const auto& [n, chVec] : channelMap) 
  {
    const auto& chName = n + "_";
    auto channel = std::vector<float>(chVec->begin(), chVec->begin() + std::min((size_t)NUM_SAMPLES, chVec->size()));

    if ( !config->hasChannel(i) ) continue;
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
        baseline += channel[j];
    }
    
    if(bl_length <=1) std::cout << "WARNING: Baseline window is trivially short, probably configured incorrectly"<<std::endl;
    baseline /= (float) bl_length;
    TF1* f = nullptr;
    if(config->channels[i].algorithm.Contains("HNR")) 
    {
      // Perform a sin fit for the baseline
      auto gr_bl = TGraph(bl_length, &(time[bl_st_idx]), &(channel[bl_st_idx]));
      f = new TF1("f_bl", "[0]+[1]*sin([2]+[3]*x)");
      f->SetParameter(0, baseline);
      f->SetParameter(1, 40./scale_factor);
      f->SetParameter(2, 1.);
      f->SetParameter(3, 2*TMath::Pi()/75.);
      f->SetParNames("const", "A", "phi0", "omega");
      auto r = gr_bl.Fit(f, "SQN");
      if (draw_debug_pulses) 
      {
        std::cout << "Harmonic Noise Removal:\n";
        std::cout << Form("const = %.2f\n", f->GetParameter(0)*scale_factor);
        std::cout << Form("A = %.2f\n", f->GetParameter(1)*scale_factor);
        std::cout << Form("phi_0 = %.2f\n", f->GetParameter(2));
        std::cout << Form("T = %.2f\n", 2*TMath::Pi()/f->GetParameter(3));
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
    var[chName+"baseline"] = scale_factor * baseline;
    // ------------- Get minimum position, max amplitude and scale the signal
    unsigned int idx_min = 0;
    float amp = 0;
    for(unsigned int j=0; j<NUM_SAMPLES; j++) 
    {
      if(config->channels[i].algorithm.Contains("HNR")) 
      {
        channel[j] = scale_factor * (channel[j] - f->Eval(time[j]));//baseline subtraction
      }
      else channel[j] = scale_factor * (channel[j] - baseline);//baseline subtraction
      bool range_check = j>bl_st_idx+bl_length && j<(int)(NUM_SAMPLES);
      bool max_check = true;
      if ( config->channels[i].counter_auto_pol_switch > 0 ) 
      {
        max_check = fabs(channel[j]) > fabs(amp);
      }
      else 
      {
        max_check = channel[j] < amp;
      }

      if(( range_check && max_check) || j == bl_st_idx+bl_length) 
      {
        idx_min = j;
        amp = channel[j];
      }
    }

    if(config->channels[i].algorithm.Contains("HNR")) {delete f;}

    //************************************************************************************
    //If the minimum point is the first sample, then the channel is bad, and we skip it.
    //************************************************************************************
    // if (idx_min == 0) continue;

    var[chName+"t_peak"] = time[idx_min];
    var[chName+"amp"] = -amp;

    float baseline_RMS = 0;
    var[chName+"noise"] = channel[bl_st_idx+5];
    for(unsigned int j=bl_st_idx; j<=(bl_st_idx+bl_length); j++) 
    {
      baseline_RMS += channel[j]*channel[j];
    }
    baseline_RMS = sqrt(baseline_RMS/bl_length);
    var[chName+"baseline_RMS"] = baseline_RMS;

    // --------------- Define pulse graph
    float * yerr = new float[NUM_SAMPLES];
    for(unsigned j = 0; j < NUM_SAMPLES; j++) yerr[j] = 0 * var[chName+"baseline_RMS"];
    TGraphErrors* pulse = new TGraphErrors(NUM_SAMPLES, time.data(), channel.data(), 0, yerr);
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
    fittable *= fabs(amp) > 5 * baseline_RMS;
    fittable *= fabs(channel[idx_min+1]) > 2*baseline_RMS;
    fittable *= fabs(channel[idx_min-1]) > 2*baseline_RMS;
    fittable *= fabs(channel[idx_min+2]) > 1*baseline_RMS;
    fittable *= fabs(channel[idx_min-2]) > 1*baseline_RMS;
    // fittable *= fabs(channel[idx_min+3]) > 2*baseline_RMS;
    // fittable *= fabs(channel[idx_min-3]) > 2*baseline_RMS;
    
    //-------------------------------------------------------------
    //If pulse is still positive change it automatically
    //-------------------------------------------------------------
    if( fittable  && !config->channels[i].algorithm.Contains("None")) 
    {
      if( var[chName+"amp"] < 0 && config->channels[i].counter_auto_pol_switch > 0 ) 
      {
        config->channels[i].polarity *= -1;
        amp = -amp;
        var[chName+"amp"] = -var[chName+"amp"];
        scale_factor = -scale_factor;
        var[chName+"baseline"] = -var[chName+"baseline"];
        for(unsigned int j=0; j<NUM_SAMPLES; j++) 
        {
          channel[j] = -channel[j];
        }
        delete pulse;
        pulse = new TGraphErrors(NUM_SAMPLES, time.data(), channel.data(), 0, yerr);
        pulse->SetNameTitle("g_"+name, "g_"+name);

        if ( config->channels[i].counter_auto_pol_switch == 10 ) 
        {
          std::cout << "[WARNING] Channel " << i << ": automatic polarity switched more than 10 times" << std::endl;
          std::cout << "[WARNING] Channel " << i << ": gonna keep inverting it for you. Better check your pulse polarity!!" << std::endl;
        }
        config->channels[i].counter_auto_pol_switch ++;
      }

      /************************************
      //Get 10% of the amplitude crossings
      ************************************
      */
      j_10_pre = GetIdxFirstCross(amp*0.1, channel, idx_min, -1);
      j_10_post = GetIdxFirstCross(amp*0.1, channel, idx_min, +1);

      // -------------- Integrate the pulse
      j_area_pre = GetIdxFirstCross(amp*0.05, channel, idx_min, -1);
      j_area_post = GetIdxFirstCross(0. , channel, idx_min, +1);
      var[chName+"integral"] = GetPulseIntegral(channel, time, j_area_pre, j_area_post);
      var[chName+"intfull"] = GetPulseIntegral(channel, time, 5, NUM_SAMPLES-5);

      // -------------- Compute rise and falling time
      float* coeff;

      j_90_pre = GetIdxFirstCross(amp*0.9, channel, j_10_pre, +1);
      AnalyticalPolinomialSolver(j_90_pre-j_10_pre+1, &(time[j_10_pre]), &(channel[j_10_pre]), 1, coeff);
      var[chName+"risetime"] = abs(0.8*var[chName+"amp"] / coeff[1]);
      delete [] coeff;

      j_90_post = GetIdxFirstCross(amp*0.9, channel, j_10_post, -1);
      AnalyticalPolinomialSolver(j_10_post-j_90_post+1, &(time[j_90_post]), &(channel[j_90_post]), 1, coeff);
      var[chName+"decaytime"] = coeff[1];
      delete [] coeff;

      /************************************
      // -------------- Global Time-offsets
      *************************************/
      double myTimeOffset = 0;
      // We need to make the code smarter to handle time offsets
      // if (correctForTimeOffsets) 
      // {
      //   myTimeOffset = timeOffset[i];
      // }

      /************************************
      // -------------- Do the gaussian fit
      *************************************/
      if( config->channels[i].algorithm.Contains("G") ) 
      {
        float frac = config->channels[i].gaus_fraction;
        unsigned int j_down = GetIdxFirstCross(amp*frac, channel, idx_min, -1);
        unsigned int j_up = GetIdxFirstCross(amp*frac, channel, idx_min, +1);
        if( j_up - j_down < 4 ) 
        {
          j_up = idx_min + 1;
          j_down = idx_min - 1;
        }

        TF1* fpeak = new TF1("fpeak"+name, "gaus", time[j_down], time[j_up]);

        float ext_sigma = time[j_up] - time[j_down];
        if (amp*frac < -baseline_RMS) ext_sigma *= 0.25;
        fpeak->SetParameter(0, amp * sqrt(2*3.14) * ext_sigma );
        fpeak->SetParameter(1, time[idx_min]);
        fpeak->SetParameter(2, ext_sigma);
        fpeak->SetLineColor(kBlue);

        TString opt = "R";
        if ( draw_debug_pulses ) opt += "+";
        else opt += "QN0";
        pulse->Fit("fpeak"+name, opt);

        var[chName+"gaus_mean"] = fpeak->GetParameter(1) + myTimeOffset;
        var[chName+"gaus_sigma"] = fpeak->GetParameter(2);
        var[chName+"gaus_chi2"] = pulse->Chisquare(fpeak, "R");

        delete fpeak;
      }
      /*******************************
      // -------------- Do  linear fit
      ********************************/
      if(config->channels[i].algorithm.Contains("Re") ) 
      {
        unsigned int i_min = GetIdxFirstCross(config->channels[i].re_bounds[0]*amp, channel, idx_min, -1);
        unsigned int i_max = GetIdxFirstCross(config->channels[i].re_bounds[1]*amp, channel, i_min  , +1);
        float t_min = time[i_min];
        float t_max = time[i_max];

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
          var[chName+Form("linear_RE_%d", (int)(100*f))] = (f*amp-Re_b)/Re_slope + myTimeOffset;
        }

        for ( auto thr : config->constant_threshold ) 
        {
            var[chName+Form("linear_RE__%dmV", (int)(fabs(thr)))] = (thr-Re_b)/Re_slope + myTimeOffset;
        }

        delete flinear;
      }
      /*************************************
      // -------------- Local polinomial fit
      **************************************/
      if ( config->constant_fraction.size() ) 
      {
        float start_level =  - 3 * baseline_RMS;
        unsigned int j_start =  GetIdxFirstCross( start_level, channel, idx_min, -1);

        for(auto f : config->constant_fraction) 
        {
          unsigned int j_st = j_start;
          if ( amp*f > start_level ) 
          {
            if ( amp*f > -baseline_RMS && verbose) 
            {
              if(N_warnings< N_warnings_to_print) 
              {
                N_warnings++;
                std::cout << Form("[WARNING] ev:%d ch:%d - fraction %.2f below noise RMS", event_n, i, f) << std::endl;
              }
              else if (N_warnings_to_print == N_warnings) 
              {
                N_warnings++;
                std::cout << "[WARNING] Max number of warnings passed. No more warnings will be printed." << std::endl;;
              }
            }
            j_st =  GetIdxFirstCross( amp*f, channel, idx_min, -1);
          }

          unsigned int j_close = GetIdxFirstCross(amp*f, channel, j_st, +1);
          if ( fabs(channel[j_close-1] - f*amp) < fabs(channel[j_close] - f*amp) ) j_close--;

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
              std::cout << Form("[WARNING] evt %d ch %d:  Short span around the closest point. Analytical fit not performed.", event_n, i) << std::endl;
              continue;
            }

            float* coeff;
            int N_add = 1;
            if (span_j + N_add + j_close < j_90_pre) 
            {
              N_add++;
            }
            AnalyticalPolinomialSolver( 2*span_j + N_add , &(channel[j_close - span_j]), &(time[j_close - span_j]), n, coeff);

            var[chName+Form("LP%d_%d", n, (int)(100*f))] = PolyEval(f*amp, coeff, n) + myTimeOffset;

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
        unsigned int j_start =  GetIdxFirstCross( start_level, channel, idx_min, -1);

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
                std::cout << Form("[WARNING] ev:%d ch:%d - thr %.2f mV below noise RMS", event_n, i, thr) << std::endl;
              }
              else if (N_warnings_to_print == N_warnings) 
              {
                N_warnings++;
                std::cout << "[WARNING] Max number of warnings passed. No more warnings will be printed." << std::endl;;
              }
            }
            j_st =  GetIdxFirstCross( thr, channel, idx_min, -1);
          }

          unsigned int j_close = GetIdxFirstCross(thr, channel, j_st, +1);
          if ( fabs(channel[j_close-1] - thr) < fabs(channel[j_close] - thr) ) j_close--;

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
                std::cout << Form("[WARNING] evt %d ch %d:  Short span around the closest point. Analytical fit not performed.", event_n, i) << std::endl;
                std::cout << j_close << "  " << span_j << "  " << thr << std::endl;
              }
              continue;
            }

            float* coeff;
            int N_add = 1;
            if (span_j + N_add + j_close < j_90_pre) 
            {
              N_add++;
            }
            AnalyticalPolinomialSolver( 2*span_j + N_add , &(channel[j_close - span_j]), &(time[j_close - span_j]), n, coeff);

            var[chName+Form("LP%d_%dmV", n, (int)(fabs(thr)))] = PolyEval(thr, coeff, n) + myTimeOffset;

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
      std::cout << "========= Event: " << event_n << " - ch: " << i << std::endl;

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
      line->DrawLine(time[0], 0, time[NUM_SAMPLES-1], 0);
      line->SetLineStyle(1);
      line->DrawLine(time[bl_st_idx], 0, time[bl_st_idx+bl_length], 0);
      line->SetLineColor(47);
      line->DrawLine(time[0], var[chName+"baseline_RMS"], time[NUM_SAMPLES-1], var[chName+"baseline_RMS"]);
      line->DrawLine(time[0], -var[chName+"baseline_RMS"], time[NUM_SAMPLES-1], -var[chName+"baseline_RMS"]);

      // Draw peak
      line->SetLineColor(8);
      line->SetLineStyle(4);
      line->DrawLine(time[0], amp, var[chName+"t_peak"], amp);
      line->DrawLine(var[chName+"t_peak"], 0, var[chName+"t_peak"], amp);

      // Draw 10% and 90% lines;
      TLine* line_lvs = new TLine();
      line_lvs->SetLineWidth(1);
      line_lvs->SetLineColor(4);
      line_lvs->DrawLine(time[0], 0.1*amp, time[NUM_SAMPLES-1], 0.1*amp);
      line_lvs->DrawLine(time[0], 0.9*amp, time[NUM_SAMPLES-1], 0.9*amp);
      // Draw constant fractions lines
      line_lvs->SetLineColor(38);
      line_lvs->SetLineStyle(10);
      for(auto f : config->constant_fraction) 
      {
        line_lvs->DrawLine(time[0], f*amp, time[NUM_SAMPLES-1], f*amp);
      }
      // Draw constant threshold lines
      line_lvs->SetLineColor(28);
      for(auto thr : config->constant_threshold) 
      {
        line_lvs->DrawLine(time[0], thr, time[NUM_SAMPLES-1], thr);
      }

      // Draw integral area
      int N_tot_integral = j_area_post-j_area_pre;
      if( N_tot_integral > 0 ) 
      {
        vector<float> aux_time = {time[j_area_pre]};
        vector<float> aux_volt = {0};
        for(unsigned int j = j_area_pre; j <= j_area_post; j++) 
        {
          aux_time.push_back(time[j]);
          aux_volt.push_back(channel[j]);
        }
        aux_time.push_back(time[j_area_post]);
        aux_volt.push_back(0);
        TGraph * integral_pulse = new TGraph(aux_time.size(), &(aux_time[0]), &(aux_volt[0]));
        integral_pulse->SetFillColor(40);
        integral_pulse->SetFillStyle(3144);
        integral_pulse->Draw("FC");
        TText* t_int = new TText(var[chName+"t_peak"]+3, amp, Form("Integral Pulse = %1.2f,  Integral Full =%1.2f ", -var[chName+"integral"], -var[chName+"intfull"]));
        t_int->SetTextAlign(kHAlignLeft+kVAlignBottom);
        t_int->Draw();

        TText* aNt = new TText(var[chName+"t_peak"]+3, amp+6, Form("LP2_50 = %1.2f, t peak = %1.2f,  amp =%1.2f ", var[chName+"LP2_50"], var[chName+"t_peak"], amp));
        aNt->SetTextAlign(kHAlignLeft+kVAlignBottom);
        aNt->Draw();

        // Draw 90% and 10% pre and post points
        TGraph* gr_pre_post = new TGraph(4);
        gr_pre_post->SetPoint(0, time[j_10_pre], channel[j_10_pre]);
        gr_pre_post->SetPoint(1, time[j_90_pre], channel[j_90_pre]);
        gr_pre_post->SetPoint(2, time[j_10_post], channel[j_10_post]);
        gr_pre_post->SetPoint(3, time[j_90_post], channel[j_90_post]);
        gr_pre_post->SetMarkerColor(4);
        gr_pre_post->Draw("P*");


        // ---------- Rising edge only inverted!! -----
        c->cd(2);
        c->SetGrid();

        if( config->channels[i].algorithm.Contains("Re") ) 
        {
          unsigned int i_min = GetIdxFirstCross(config->channels[i].re_bounds[0]*amp, channel, idx_min, -1);
          unsigned int i_max = GetIdxFirstCross(config->channels[i].re_bounds[1]*amp, channel, i_min, +1);
          float y[2], x[2];
          x[0] = channel[i_min];
          x[1] = channel[i_max];
          y[0] = (channel[i_min] - Re_b)/Re_slope;
          y[1] = (channel[i_max] - Re_b)/Re_slope;

          TGraph* gr_Re = new TGraph(2, x, y);

          gr_Re->SetLineColor(46);
          gr_Re->SetLineWidth(1);
          gr_Re->SetLineStyle(7);
          gr_Re->Draw("CP");
        }

        unsigned int j_begin = j_10_pre - 6;
        unsigned int j_span = j_90_pre - j_10_pre + 10;
        TGraphErrors* inv_pulse = new TGraphErrors(j_span, &(channel[j_begin]), &(time[j_begin]), yerr);
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
            float tt = time[jj] + (time[jj+1] - time[jj]) * kk/overstep;
            float cc = WSInterp(tt, NUM_SAMPLES, time, channel);

            t_WS.push_back(tt);
            c_WS.push_back(cc);
          }
        }
        TGraph* inv_pulse_WS = new TGraph(t_WS.size(), &(c_WS[0]), &(t_WS[0]));
        inv_pulse_WS->SetMarkerStyle(7);
        inv_pulse_WS->SetMarkerColor(2);
        // inv_pulse_WS->Draw("P");

        TGraph* gr_inv_pre_post = new TGraph(2);
        gr_inv_pre_post->SetPoint(0, channel[j_10_pre], time[j_10_pre]);
        gr_inv_pre_post->SetPoint(1, channel[j_90_pre], time[j_90_pre]);
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
          line_lvs->DrawLine(amp*f, time[j_begin], amp*f, time[j_90_pre + 3]);
          for(auto n : config->channels[i].PL_deg) 
          {
            vector<float> polyval;
            for(unsigned int j = poly_bounds[count].first; j <= poly_bounds[count].second; j++) 
            {
              polyval.push_back(PolyEval(channel[j], coeff_poly_fit[count], n));
            }

            TGraph* g_poly = new TGraph(polyval.size(), &(channel[poly_bounds[count].first]), &(polyval[0]) );
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
          line_lvs->DrawLine(thr, time[j_begin], thr, time[j_90_pre + 3]);
          for(auto n : config->channels[i].PL_deg) 
          {
            vector<float> polyval;
            for(unsigned int j = poly_bounds[count].first; j <= poly_bounds[count].second; j++) 
            {
              polyval.push_back(PolyEval(channel[j], coeff_poly_fit[count], n));
            }

            TGraph* g_poly = new TGraph(polyval.size(), &(channel[poly_bounds[count].first]), &(polyval[0]) );
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
    i++;
  }
}

void DRSAnalyzer::InitLoop() 
{
  for (const auto& name : branch_names) 
  {
      auto* br = tree_in->GetBranch(name);
      auto* leaf = br->GetLeaf(name);
      const TString& typeName = leaf->GetTypeName();
      void* buffer = nullptr;

      if (typeName == "Int_t") {
          buffer = new Int_t;
          tree_in->SetBranchAddress(name, (Int_t*)buffer);
          tree->Branch(name, (Int_t*)buffer, name + "/I")->SetBasketSize(64 * 1024 * 1024); // Increase buffer size
      } else if (typeName == "UInt_t") {
          buffer = new UInt_t;
          tree_in->SetBranchAddress(name, (UInt_t*)buffer);
          tree->Branch(name, (UInt_t*)buffer, name + "/i")->SetBasketSize(64 * 1024 * 1024); // Increase buffer size
      } else if (typeName == "ULong64_t") {
          buffer = new ULong64_t;
          tree_in->SetBranchAddress(name, (ULong64_t*)buffer);
          tree->Branch(name, (ULong64_t*)buffer, name + "/l")->SetBasketSize(64 * 1024 * 1024); // Increase buffer size
      } else if (typeName == "Float_t") {
          buffer = new Float_t;
          tree_in->SetBranchAddress(name, (Float_t*)buffer);
          tree->Branch(name, (Float_t*)buffer, name + "/F")->SetBasketSize(64 * 1024 * 1024); // Increase buffer size
      } else if (typeName == "Double_t") {
          buffer = new Double_t;
          tree_in->SetBranchAddress(name, (Double_t*)buffer);
          tree->Branch(name, (Double_t*)buffer, name + "/D")->SetBasketSize(64 * 1024 * 1024); // Increase buffer size
      } else if (typeName.BeginsWith("vector")) {
          auto* vec = new std::vector<float>;
          buffer = vec;
          tree_in->SetBranchAddress(name, &vec);
          tree->Branch(name, &vec)->SetBasketSize(64 * 1024 * 1024); // Increase buffer size
      } else {
          std::cerr << "Unsupported type: " << typeName << " for branch " << name << std::endl;
          continue;
      }

      branch_buffers[name] = buffer;
  }

  for(auto& [name, channel] : channelMap) 
  {
    tree_in->SetBranchAddress(name, &channel);
    tree->Branch(name, &channel)->SetBasketSize(64 * 1024 * 1024); // Increase buffer size
  }
  tree_in->GetEntry(0);

  NUM_CHANNELS = channelMap.size();
  std::cout<<"Number of Channels: "<<NUM_CHANNELS<<std::endl;
  NUM_TIMES = 0;
  NUM_SAMPLES = 900;
  // NUM_SAMPLES = channelMap.begin()->second->size();
  time = std::vector<float>(NUM_SAMPLES);
  for (unsigned int i = 0; i < NUM_SAMPLES; i++){ time[i] = (200.0 / 1000.0) * i; }

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

  if( at_least_1_gaus_fit ) 
  {
    var_names.push_back("gaus_mean");
    var_names.push_back("gaus_sigma");
    var_names.push_back("gaus_chi2");
  }

  if( at_least_1_rising_edge ) 
  {
    for (auto f : config->constant_fraction) 
    {
      var_names.push_back(Form("linear_RE_%d", (int)(100*f)));
    }

    for (auto thr : config->constant_threshold) 
    {
      var_names.push_back(Form("linear_RE__%dmV", (int)(fabs(thr))));
    }
  }

  for(const auto& ch : channelMap) 
  {
    const auto& chName = ch.first + "_";

    for(const auto& vn : var_names)
    {
      TString n = chName + vn;
      var[n] = float(-999.9);
      tree->Branch(n, &var[n], n+"/F")->SetBasketSize(64 * 1024 * 1024); // Increase buffer size
    }
  }
}

int DRSAnalyzer::GetChannelsMeasurement(int i_aux) 
{
  ResetAnalysisVariables();
  tree_in->GetEntry(i_aux);
  return 0;
}

void DRSAnalyzer::ResetVar() 
{
  for(const auto& ch : channelMap) 
  {
    const auto& chName = ch.first + "_";
    for(const auto& vn : var_names) 
    {
      TString n = chName + vn;
      var[n] = 0;
    }
  }
}

void DRSAnalyzer::ResetAnalysisVariables() 
{
  for (auto& ch : channelMap) 
  {
    ch.second->clear();
  }
}

float DRSAnalyzer::GetPulseIntegral(std::vector<float> a, std::vector<float> t, unsigned int i_st, unsigned int i_stop) //returns charge in pC asssuming 50 Ohm termination
{
  //Simpson's Rule for equaled space with Cartwright correction for unequaled space
  float integral = 0.;
  for (unsigned int i=i_st; i < i_stop-2 ; i+=2) 
  {
    float aux = ( 2-(t[i+2]-t[i+1])/(t[i+1]-t[i]) ) * a.at(i);
    aux += (t[i+2]-t[i])*(t[i+2]-t[i])/((t[i+2]-t[i+1])*(t[i+1]-t[i])) * a.at(i+1);
    aux += ( 2-(t[i+1]-t[i])/(t[i+2]-t[i+1]) ) * a.at(i+2);

    integral += ( (t[i+2]-t[i]) / 6.0 ) * aux;
  }

  integral *= -1.0;//1e-9 * 1e-3 * (1.0/50.0) * 1e12; //in units of pC, for 50 [Ohms] termination
  return integral;
}

unsigned int DRSAnalyzer::GetIdxClosest(float value, std::vector<float> v, unsigned int i_st, int direction) 
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

unsigned int DRSAnalyzer::GetIdxFirstCross(float value, std::vector<float> v, unsigned int i_st, int direction) 
{
  unsigned int idx_end = direction>0 ? NUM_SAMPLES-1 : 0;
  bool rising = value > v.at(i_st)? true : false; // Check if the given max value is greater than the given waveform amplitude at a given start time (i_st).

  // if it is: the value of the waveform needs to rise and otherwise it doesn't need so the function will return the i_st.
  unsigned int i = i_st;
  while( (i != idx_end)) 
  { // loop over the time bins
    if(rising && v.at(i) > value) break; // check if the rising variable is true and if the waveform value is going higher than the given amplitude value. And stops the loop.
    else if( !rising && v.at(i) < value) break; // check if the rising variable is false and if the waveform value is still lower than the given amplitude value. And stops the loop.
    i += direction; // Otherwise it moves to the next time bin.
  }

  return i;
}

void DRSAnalyzer::AnalyticalPolinomialSolver(unsigned int Np, float* in_x, float* in_y, unsigned int deg, float*& out_coeff)
{
    if (deg <= 0 || deg > 3) {
        std::cerr << "[ERROR]: You don't need AnalyticalPolinomialSolver for this\n";
        out_coeff = nullptr;
        return;
    }
    if (Np < deg + 1) {
        std::cerr << "[WARNING]: Not enough points for requested polynomial degree\n";
        out_coeff = nullptr;
        return;
    }

    TMatrixD A(Np, deg + 1);
    for (unsigned int i = 0; i < Np; ++i) {
        double x = in_x[i];
        A(i, 0) = 1.0;
        if (deg >= 1) A(i, 1) = x;
        if (deg >= 2) A(i, 2) = x * x;
        if (deg >= 3) A(i, 3) = x * x * x;
    }

    TVectorD y(Np);
    for (unsigned int i = 0; i < Np; ++i)
        y[i] = in_y[i];

    TMatrixD At = TMatrixD(TMatrixD::kTransposed, A);
    TMatrixD AtA = At * A;
    TVectorD Aty = At * y;

    // Suppress ROOT error messages
    int oldErrorLevel = gErrorIgnoreLevel;
    gErrorIgnoreLevel = kBreak;  // suppress all ROOT errors

    TDecompChol chol(AtA);
    Bool_t ok;
    TVectorD coeffs = chol.Solve(Aty, ok);

    // Restore original error level
    gErrorIgnoreLevel = oldErrorLevel;

    out_coeff = new float[deg + 1];

    if (!ok) {
        // std::cerr << "[WARNING]: Matrix not positive definite. Filling zeros.\n";
        for (unsigned int i = 0; i <= deg; ++i)
            out_coeff[i] = 0.0;
    } else {
        for (unsigned int i = 0; i <= deg; ++i)
            out_coeff[i] = coeffs[i];
    }
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

float DRSAnalyzer::WSInterp(float t, int N, std::vector<float> tn, std::vector<float> cn) 
{
  float out = 0;
  float dt = (tn[0] - tn[N-1])/N;
  for(unsigned i = 0; i < N; i++) 
  {
    float x = (t - tn[i])/dt;
    out += cn.at(i) * sin(3.14159265358 * x) / (3.14159265358 * x);
  }
  return out;
};
