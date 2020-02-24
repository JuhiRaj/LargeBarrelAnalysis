/**
 *  @copyright Copyright 2018 The J-PET Framework Authors. All rights reserved.
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may find a copy of the License in the LICENCE file.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *  @file EventCategorizer.cpp
 */

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include <JPetMCHit/JPetMCHit.h>
#include "EventCategorizerTools.h"
#include "EventCategorizer.h"
#include <iostream>
#include "TMath.h"
#include "math.h"
using namespace jpet_options_tools;
using namespace std;

EventCategorizer::EventCategorizer(const char* name): JPetUserTask(name) {}


EventCategorizer::~EventCategorizer() {}


bool EventCategorizer::init()
{
  INFO("Event analysis started.");
  
  // Input events type
  fOutputEvents = new JPetTimeWindow("JPetEvent");

	//Juhi created histograms//


    getStatistics().createHistogram(new TH1F("Multiplicity", "Multiplicity",
      11, -0.5, 10.5));
    getStatistics().getHisto1D("Multiplicity")->GetXaxis()->SetTitle("Multiplicity");
    getStatistics().getHisto1D("Multiplicity")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Multiplicity_4", "Multiplicity_4",
      11, -0.5, 10.5));
    getStatistics().getHisto1D("Multiplicity_4")->GetXaxis()->SetTitle("Multiplicity");
    getStatistics().getHisto1D("Multiplicity_4")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Multiplicity_Gen_Mul_Cut", "Multiplicity_Gen_Mul_Cut",
      11, -0.5, 10.5));
    getStatistics().getHisto1D("Multiplicity_Gen_Mul_Cut")->GetXaxis()->SetTitle("Multiplicity_Gen_Mul_Cut");
    getStatistics().getHisto1D("Multiplicity_Gen_Mul_Cut")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("True_MC_Multiplicity", "True_MC_Multiplicity",
      401, -10.5, 410.5));
    getStatistics().getHisto1D("True_MC_Multiplicity")->GetXaxis()->SetTitle("True_MC_Multiplicity");
    getStatistics().getHisto1D("True_MC_Multiplicity")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("True_MC_Multiplicity_Cut", "True_MC_Multiplicity_Cut",
      401, -10.5, 410.5));
    getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->GetXaxis()->SetTitle("True_MC_Multiplicity_Cut");
    getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Energy_Primary_hit", "Energy_Primary",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_Primary_hit")->GetXaxis()->SetTitle("Energy_{i} [keV]");
    getStatistics().getHisto1D("Energy_Primary_hit")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Energy_Secondary_hit", "Energy_Secondary",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_Secondary_hit")->GetXaxis()->SetTitle("Energy'_{i} [keV]");
    getStatistics().getHisto1D("Energy_Secondary_hit")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Energy_Primary_mc_hit", "Energy_Primary",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_Primary_mc_hit")->GetXaxis()->SetTitle("Energy_{i} [keV]");
    getStatistics().getHisto1D("Energy_Primary_mc_hit")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Energy_Secondary_mc_hit", "Energy_Secondary",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_Secondary_mc_hit")->GetXaxis()->SetTitle("Energy'_{i} [keV]");
    getStatistics().getHisto1D("Energy_Secondary_mc_hit")->GetYaxis()->SetTitle("Counts");

    
    getStatistics().createHistogram(new TH1F("ZPosition_Primary_hit", "ZPosition_Primary_hit",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("ZPosition_Primary_hit")->GetXaxis()->SetTitle("Z_Pos [cm]");
    getStatistics().getHisto1D("ZPosition_Primary_hit")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("ZPosition_Secondary_hit", "ZPosition_Secondary_hit",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("ZPosition_Secondary_hit")->GetXaxis()->SetTitle("Z_Pos' [cm]");
    getStatistics().getHisto1D("ZPosition_Secondary_hit")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("ZPosition_Primary_mc_hit", "ZPosition_Primary_mc_hit",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("ZPosition_Primary_mc_hit")->GetXaxis()->SetTitle("Z_Pos [cm]");
    getStatistics().getHisto1D("ZPosition_Primary_mc_hit")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("ZPosition_Secondary_mc_hit", "ZPosition_Secondary_mc_hit",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("ZPosition_Secondary_mc_hit")->GetXaxis()->SetTitle("Z_Pos [cm]");
    getStatistics().getHisto1D("ZPosition_Secondary_mc_hit")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("Z_resolution", "Z_resolution",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("Z_resolution")->GetXaxis()->SetTitle("Z_Pos [cm]");
    getStatistics().getHisto1D("Z_resolution")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("Y_resolution", "Y_resolution",
      1100, -5.5, 5.5));
    getStatistics().getHisto1D("Y_resolution")->GetXaxis()->SetTitle("Y_Pos [cm]");
    getStatistics().getHisto1D("Y_resolution")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("X_resolution", "X_resolution",
      1100, -5.5, 5.5));
    getStatistics().getHisto1D("X_resolution")->GetXaxis()->SetTitle("X_Pos [cm]");
    getStatistics().getHisto1D("X_resolution")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Energy_resolution", "Energy_resolution",
      2010, -100.5, 100.5));
    getStatistics().getHisto1D("Energy_resolution")->GetXaxis()->SetTitle("Energy [keV]");
    getStatistics().getHisto1D("Energy_resolution")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("EV_Resolution", "EV_Resolution",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_Resolution")->GetXaxis()->SetTitle("Expectation Value");
    getStatistics().getHisto1D("EV_Resolution")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Time_resolution", "Time_resolution",
      10010, -500.5, 500.5));
    getStatistics().getHisto1D("Time_resolution")->GetXaxis()->SetTitle("Time [ns]");
    getStatistics().getHisto1D("Time_resolution")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Scintillator_ID", "Scintillator_ID",
      193, -0.5, 192.5));
    getStatistics().getHisto1D("Scintillator_ID")->GetXaxis()->SetTitle("Scintillator_ID");
    getStatistics().getHisto1D("Scintillator_ID")->GetYaxis()->SetTitle("Counts");



    getStatistics().createHistogram(new TH1F("Energy_all_primary", "Energy_all_primary",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_all_primary")->GetXaxis()->SetTitle("Energy_{i} [keV]");
    getStatistics().getHisto1D("Energy_all_primary")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("Energy_secondary", "Energy_secondary",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_secondary")->GetXaxis()->SetTitle("Energy_{i} [keV]");
    getStatistics().getHisto1D("Energy_secondary")->GetYaxis()->SetTitle("Counts");	

    getStatistics().createHistogram(new TH2F("Angle_3D_True","Angle_3D_True",191, -0.5, 190., 191, -0.5, 190.));
    getStatistics().getHisto2D("Angle_3D_True")->GetXaxis()->SetTitle("#theta_{i} [deg]");
    getStatistics().getHisto2D("Angle_3D_True")->GetYaxis()->SetTitle("#theta_{j} [deg]");

    getStatistics().createHistogram(new TH2F("Angle_3D_Hit","Angle_3D_Hit",191, -0.5, 190., 191, -0.5, 190.));
    getStatistics().getHisto2D("Angle_3D_Hit")->GetXaxis()->SetTitle("#theta_{i} [deg]");
    getStatistics().getHisto2D("Angle_3D_Hit")->GetYaxis()->SetTitle("#theta_{j} [deg]");


    getStatistics().createHistogram(new TH2F("Angle_3D_True_Cut","Angle_3D_True_Cut",191, -0.5, 190., 191, -0.5, 190.));
    getStatistics().getHisto2D("Angle_3D_True_Cut")->GetXaxis()->SetTitle("#theta_{i} [deg]");
    getStatistics().getHisto2D("Angle_3D_True_Cut")->GetYaxis()->SetTitle("#theta_{j} [deg]");

    getStatistics().createHistogram(new TH2F("Angle_3D_Hit_Cut","Angle_3D_Hit_Cut",191, -0.5, 190., 191, -0.5, 190.));
    getStatistics().getHisto2D("Angle_3D_Hit_Cut")->GetXaxis()->SetTitle("#theta_{i} [deg]");
    getStatistics().getHisto2D("Angle_3D_Hit_Cut")->GetYaxis()->SetTitle("#theta_{j} [deg]");





    getStatistics().createHistogram(new TH1F("Scatter_Hit_Lowest", "Scatter_Hit_Lowest",
      1201, -600.5, 600.5));
    getStatistics().getHisto1D("Scatter_Hit_Lowest")->GetXaxis()->SetTitle("Scatter_Hit_Lowest [ns]");
    getStatistics().getHisto1D("Scatter_Hit_Lowest")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Scatter_MC_Lowest", "Scatter_MC_Lowest",
      1201, -600.5, 600.5));
    getStatistics().getHisto1D("Scatter_MC_Lowest")->GetXaxis()->SetTitle("Scatter_MC_Lowest [ns]");
    getStatistics().getHisto1D("Scatter_MC_Lowest")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH2F("Scatter_Hit_21","Scatter_Hit_21",1201, -600.5, 600.5, 1201, -600.5, 600.5));
    getStatistics().getHisto2D("Scatter_Hit_21")->GetXaxis()->SetTitle("Scatter_Test_Value_{2} [ns]");
    getStatistics().getHisto2D("Scatter_Hit_21")->GetYaxis()->SetTitle("Scatter_Test_Value_{1} [ns]");

    getStatistics().createHistogram(new TH2F("Scatter_MC_21","Scatter_MC_21",1201, -600.5, 600.5, 1201, -600.5, 600.5));
    getStatistics().getHisto2D("Scatter_MC_21")->GetXaxis()->SetTitle("Scatter_Test_Value_{2} [ns]");
    getStatistics().getHisto2D("Scatter_MC_21")->GetYaxis()->SetTitle("Scatter_Test_Value_{1} [ns]");

    getStatistics().createHistogram(new TH1F("Scatter_Hit_Cut", "Scatter_Hit_Cut",
      1201, -600.5, 600.5));
    getStatistics().getHisto1D("Scatter_Hit_Cut")->GetXaxis()->SetTitle("Scatter_Hit_Cut [ns]");
    getStatistics().getHisto1D("Scatter_Hit_Cut")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Scatter_MC_Cut", "Scatter_MC_Cut",
      1201, -600.5, 600.5));
    getStatistics().getHisto1D("Scatter_MC_Cut")->GetXaxis()->SetTitle("Scatter_MC_Cut [ns]");
    getStatistics().getHisto1D("Scatter_MC_Cut")->GetYaxis()->SetTitle("Counts");

	//Cosine alpha - absolutely plain//
  
   getStatistics().createHistogram(new TH1F("EV_Hit", "EV_Hit",
      25, -1.0, 1.0));
    getStatistics().getHisto1D("EV_Hit")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_Hit")->GetYaxis()->SetTitle("Counts");

   getStatistics().createHistogram(new TH1F("EV_MC", "EV_MC",
      25, -1.0, 1.0));
    getStatistics().getHisto1D("EV_MC")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_MC")->GetYaxis()->SetTitle("Counts");

//21//
   getStatistics().createHistogram(new TH1F("E2K1_Hit", "E2K1_Hit",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E2K1_Hit")->GetXaxis()->SetTitle("Cosine #alpha_{21}");
    getStatistics().getHisto1D("E2K1_Hit")->GetYaxis()->SetTitle("Counts");

  getStatistics().createHistogram(new TH1F("E2K1_MC", "E2K1_MC",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E2K1_MC")->GetXaxis()->SetTitle("Cosine #alpha_{21}");
    getStatistics().getHisto1D("E2K1_MC")->GetYaxis()->SetTitle("Counts");

//20//

 getStatistics().createHistogram(new TH1F("E2K0_Hit", "E2K0_Hit",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E2K0_Hit")->GetXaxis()->SetTitle("Cosine #alpha_{20}");
    getStatistics().getHisto1D("E2K0_Hit")->GetYaxis()->SetTitle("Counts");

  getStatistics().createHistogram(new TH1F("E2K0_MC", "E2K0_MC",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E2K0_MC")->GetXaxis()->SetTitle("Cosine #alpha_{20}");
    getStatistics().getHisto1D("E2K0_MC")->GetYaxis()->SetTitle("Counts");

//10//

    getStatistics().createHistogram(new TH1F("E1K0_Hit", "E1K0_Hit",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E1K0_Hit")->GetXaxis()->SetTitle("Cosine #alpha_{10}");
    getStatistics().getHisto1D("E1K0_Hit")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("E1K0_MC", "E1K0_MC",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E1K0_MC")->GetXaxis()->SetTitle("Cosine #alpha_{10}");
    getStatistics().getHisto1D("E1K0_MC")->GetYaxis()->SetTitle("Counts");

//12//

    getStatistics().createHistogram(new TH1F("E1K2_Hit", "E1K2_Hit",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E1K2_Hit")->GetXaxis()->SetTitle("Cosine #alpha_{12}");
    getStatistics().getHisto1D("E1K2_Hit")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("E1K2_MC", "E1K2_MC",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E1K2_MC")->GetXaxis()->SetTitle("Cosine #alpha_{12}");
    getStatistics().getHisto1D("E1K2_MC")->GetYaxis()->SetTitle("Counts");


//02//

    getStatistics().createHistogram(new TH1F("E0K2_Hit", "E0K2_Hit",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E0K2_Hit")->GetXaxis()->SetTitle("Cosine #alpha_{02}");
    getStatistics().getHisto1D("E0K2_Hit")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("E0K2_MC", "E0K2_MC",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E0K2_MC")->GetXaxis()->SetTitle("Cosine #alpha_{02}");
    getStatistics().getHisto1D("E0K2_MC")->GetYaxis()->SetTitle("Counts");


//01//

    getStatistics().createHistogram(new TH1F("E0K1_Hit", "E0K1_Hit",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E0K1_Hit")->GetXaxis()->SetTitle("Cosine #alpha_{01}");
    getStatistics().getHisto1D("E0K1_Hit")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("E0K1_MC", "E0K1_MC",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("E0K1_MC")->GetXaxis()->SetTitle("Cosine #alpha_{01}");
    getStatistics().getHisto1D("E0K1_MC")->GetYaxis()->SetTitle("Counts");





     // First Cut - Z-Position//
    getStatistics().createHistogram(new TH1F("Z_position_before", "Z_position_before",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("Z_position_before")->GetXaxis()->SetTitle("Z_Pos [cm]");
    getStatistics().getHisto1D("Z_position_before")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Z_position_after", "Z_position_after",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("Z_position_after")->GetXaxis()->SetTitle("Z_Pos [cm]");
    getStatistics().getHisto1D("Z_position_after")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("EV_Hit_Z_Pos", "EV_Hit_Z_Pos",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_Hit_Z_Pos")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_Hit_Z_Pos")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("EV_MC_Z_Pos", "EV_MC_Z_Pos",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_MC_Z_Pos")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_MC_Z_Pos")->GetYaxis()->SetTitle("Counts");


    // Second Cut - Z-Position//

    getStatistics().createHistogram(new TH1F("Energy_before", "Energy_before",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_before")->GetXaxis()->SetTitle("Energy_{i} [keV]");
    getStatistics().getHisto1D("Energy_before")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Energy_after", "Energy_after",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_after")->GetXaxis()->SetTitle("Energy_{i} [keV]");
    getStatistics().getHisto1D("Energy_after")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("EV_Hit_Energy", "EV_Hit_Energy",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_Hit_Energy")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_Hit_Energy")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("EV_MC_Energy", "EV_MC_Energy",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_MC_Energy")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_MC_Energy")->GetYaxis()->SetTitle("Counts");


    // Third Cut - TOF //
 
    getStatistics().createHistogram(new TH1F("TOF_before", "TOF_before",
      2010, -100.5, 100.5));
    getStatistics().getHisto1D("TOF_before")->GetXaxis()->SetTitle("TOF_before [ns]");
    getStatistics().getHisto1D("TOF_before")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("TOF_after", "TOF_after",
      2010, -100.5, 100.5));
    getStatistics().getHisto1D("TOF_after")->GetXaxis()->SetTitle("TOF_after [ns]");
    getStatistics().getHisto1D("TOF_after")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("EV_Hit_TOF", "EV_Hit_TOF",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_Hit_TOF")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_Hit_TOF")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("EV_MC_TOF", "EV_MC_TOF",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_MC_TOF")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_MC_TOF")->GetYaxis()->SetTitle("Counts");

     getStatistics().createHistogram(new TH1F("TOF_Resolution", "TOF_Resolution",
      2010, -100.5, 100.5));
    getStatistics().getHisto1D("TOF_Resolution")->GetXaxis()->SetTitle("TOF_Resolution [ns]");
    getStatistics().getHisto1D("TOF_Resolution")->GetYaxis()->SetTitle("Counts");



   // Fourth Cut - Distance Surface //

    getStatistics().createHistogram(new TH1F("Distance_Surface_before", "Distance_Surface_before",
      510, -0.5, 50.5));
    getStatistics().getHisto1D("Distance_Surface_before")->GetXaxis()->SetTitle("Distance_Surface_before [cm]");
    getStatistics().getHisto1D("Distance_Surface_before")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Distance_Surface_after", "Distance_Surface_after",
      510, -0.5, 50.5));
    getStatistics().getHisto1D("Distance_Surface_after")->GetXaxis()->SetTitle("Distance_Surface_after [cm]");
    getStatistics().getHisto1D("Distance_Surface_after")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("EV_Hit_Distance", "EV_Hit_Distance",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_Hit_Distance")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_Hit_Distance")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("EV_MC_Distance", "EV_MC_Distance",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_MC_Distance")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_MC_Distance")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Distance_Resolution", "Distance_Resolution",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("Distance_Resolution")->GetXaxis()->SetTitle("Distance_Resolution [cm]");
    getStatistics().getHisto1D("Distance_Resolution")->GetYaxis()->SetTitle("Counts");


    // Fifth Cut - Angle 3D //

    getStatistics().createHistogram(new TH2F("Angle_Sum_before","Angle_Sum_before",251, -0.5, 250., 251, -0.5, 250.));
    getStatistics().getHisto2D("Angle_Sum_before")->GetXaxis()->SetTitle("#theta_{i+j} [deg]");
    getStatistics().getHisto2D("Angle_Sum_before")->GetYaxis()->SetTitle("#theta_{i-j} [deg]");

    getStatistics().createHistogram(new TH2F("Angle_Sum_after","Angle_Sum_after",251, -0.5, 250., 251, -0.5, 250.));
    getStatistics().getHisto2D("Angle_Sum_after")->GetXaxis()->SetTitle("#theta_{i+j} [deg]");
    getStatistics().getHisto2D("Angle_Sum_after")->GetYaxis()->SetTitle("#theta_{i-j} [deg]");

    getStatistics().createHistogram(new TH1F("EV_Hit_Angle", "EV_Hit_Angle",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_Hit_Angle")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_Hit_Angle")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("EV_MC_Angle", "EV_MC_Angle",
      3000, -1.5, 1.5));
    getStatistics().getHisto1D("EV_MC_Angle")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_MC_Angle")->GetYaxis()->SetTitle("Counts");

     getStatistics().createHistogram(new TH1F("Angle_Resolution", "Angle_Resolution",
      2100, -10.5, 10.5));
    getStatistics().getHisto1D("Angle_Resolution")->GetXaxis()->SetTitle("Angle_Resolution [deg]");
    getStatistics().getHisto1D("Angle_Resolution")->GetYaxis()->SetTitle("Counts");




    getStatistics().createHistogram(new TH1F("pol_true", "pol_true",
      5000, -0.1, 0.02));
    getStatistics().getHisto1D("pol_true")->GetXaxis()->SetTitle("pol_true");
    getStatistics().getHisto1D("pol_true")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("pol_reco", "pol_reco",
     300, -1.5, 1.5));
    getStatistics().getHisto1D("pol_reco")->GetXaxis()->SetTitle("pol_reco");
    getStatistics().getHisto1D("pol_reco")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH2F("Polarization_Check","Polarization_Check",300, -1.5, 1.5, 300, -1.5, 1.5));
    getStatistics().getHisto2D("Polarization_Check")->GetXaxis()->SetTitle("Reconstructed (kxk')");
    getStatistics().getHisto2D("Polarization_Check")->GetYaxis()->SetTitle("Monte Carlo True Polarization");

    getStatistics().createHistogram(new TH2F("Polarization_Angle_Check","Polarization_Angle_Check",181, -0.5, 180.5, 181, -0.5, 180.5));
    getStatistics().getHisto2D("Polarization_Angle_Check")->GetXaxis()->SetTitle("Reconstructed (kxk')");
    getStatistics().getHisto2D("Polarization_Angle_Check")->GetYaxis()->SetTitle("Monte Carlo True Polarization");






  return true;


}


Bool_t comparison(const pair < double, JPetHit > & a, const pair < double, JPetHit > & b) {
  
	return a.first < b.first;

}

Bool_t comparison_MC(const pair < double, JPetMCHit > & a, const pair < double, JPetMCHit > & b) {
  
	return a.first < b.first;

}



Bool_t comparison3(const pair < double, pair < JPetHit, JPetHit >> & a, const pair < double, pair < JPetHit, JPetHit >> & b) {

  return abs(a.first) < abs(b.first);

}


Bool_t comparison_scatter_Hit(const pair < double, pair < JPetHit, int >> & a, const pair < double, pair < JPetHit, int >> & b) {

  return abs(a.first) < abs(b.first);

}


Bool_t comparison_scatter_MC(const pair < double, pair < JPetMCHit, int >> & a, const pair < double, pair < JPetMCHit, int >> & b) {

  return abs(a.first) < abs(b.first);

}



double CalExpecValueMC(const JPetMCHit Hit1, const JPetMCHit Hit2, const JPetMCHit Hit3)
{

 
  	TVector3 k1 = 	Hit1.getMomentum();
	TVector3 k1p = 	Hit2.getMomentum();
	TVector3 k2 = 	Hit3.getMomentum();

       TVector3 epsilon = k1.Cross(k1p);

       double ExpecValue = ((epsilon.Dot(k2)) / ((epsilon.Mag()) * (k2.Mag())));

       return ExpecValue;

}

	TVector3 Center(0.0, 0.0, 0.0);  //Center of the Geometry

double CalExpecValueHit(const JPetHit Hit1, const JPetHit Hit2, const JPetHit Hit3)
{

  TVector3 k1(Hit1.getPosX() - Center.X() , Hit1.getPosY() - Center.Y() , Hit1.getPosZ() - Center.Z());
  TVector3 k1p(Hit2.getPosX() - Hit1.getPosX(), Hit2.getPosY() - Hit1.getPosY(), Hit2.getPosZ() - Hit1.getPosZ());
  TVector3 k2(Hit3.getPosX() - Center.X() , Hit3.getPosY() - Center.Y(), Hit3.getPosZ() - Center.Z());

       TVector3 epsilon = k1.Cross(k1p);

       double ExpecValue = ((epsilon.Dot(k2)) / ((epsilon.Mag()) * (k2.Mag())));

       return ExpecValue;

}


double CalScatterTestHit(const JPetHit & Hit1,
  const JPetHit & Hit2)

{

  double timeDiff = ((Hit2.getTime() / 1000.0) - (Hit1.getTime() / 1000.0));
  float Dist = sqrt(pow(Hit2.getPosX() - Hit1.getPosX(), 2) +
    pow(Hit2.getPosY() - Hit1.getPosY(), 2) +
    pow(Hit2.getPosZ() - Hit1.getPosZ(), 2));
  Dist = Dist / 29.979246;
  double Scat = (timeDiff - Dist);
  return Scat;

}

double CalScatterTestTrueMC(const JPetMCHit & Hit1,
  const JPetMCHit & Hit2)

{

  double timeDiff = ((Hit2.getTime() / 1000.0) - (Hit1.getTime() / 1000.0));
  float Dist = sqrt(pow(Hit2.getPosX() - Hit1.getPosX(), 2) +
    pow(Hit2.getPosY() - Hit1.getPosY(), 2) +
    pow(Hit2.getPosZ() - Hit1.getPosZ(), 2));
  Dist = Dist / 29.979246;
  double Scat = (timeDiff - Dist);
  return Scat;

}



double CalcTOFHit(const JPetHit & Hit, TVector3 Center) {

  double Dist = sqrt(pow(Hit.getPosX() - Center.X(), 2) + pow(Hit.getPosY() - Center.Y(), 2) + pow(Hit.getPosZ() - Center.Z(), 2));

  double CalTime = Dist / 29.979246; //velocity of light in cm/s, returns time in ns

  double HitTime = Hit.getTime() / 1000.0;

  double TOF = (HitTime - CalTime);

  return TOF;

}



double CalcTOFMC(const JPetMCHit & Hit, TVector3 Center) {

  double Dist = sqrt(pow(Hit.getPosX() - Center.X(), 2) + pow(Hit.getPosY() - Center.Y(), 2) + pow(Hit.getPosZ() - Center.Z(), 2));

  double CalTime = Dist / 29.979246; //velocity of light in cm/s, returns time in ns

  double HitTime = Hit.getTime() / 1000.0;

  double TOF = (HitTime - CalTime);

  return TOF;

}




double CalDistofSurfaceHit(const JPetHit Hit1, const JPetHit Hit2, const JPetHit Hit3, TVector3 Center)

{

  TVector3 vec1(Hit2.getPosX() - Hit1.getPosX(), Hit2.getPosY() - Hit1.getPosY(), Hit2.getPosZ() - Hit1.getPosZ());
  TVector3 vec2(Hit3.getPosX() - Hit2.getPosX(), Hit3.getPosY() - Hit2.getPosY(), Hit3.getPosZ() - Hit2.getPosZ());
  TVector3 crossProd = vec1.Cross(vec2);

  double Dcoeef = -crossProd(0) * Hit2.getPosX() - crossProd(1) * Hit2.getPosY() - crossProd(2) * Hit2.getPosZ();

  double distancefromZero = fabs(crossProd(0) * Center.X() + crossProd(1) * Center.Y() + crossProd(2) * Center.Z() + Dcoeef) / crossProd.Mag();
  return distancefromZero;

}


double CalDistofSurfaceMC(const JPetMCHit Hit1, const JPetMCHit Hit2, const JPetMCHit Hit3, TVector3 Center)

{

  TVector3 vec1(Hit2.getPosX() - Hit1.getPosX(), Hit2.getPosY() - Hit1.getPosY(), Hit2.getPosZ() - Hit1.getPosZ());
  TVector3 vec2(Hit3.getPosX() - Hit2.getPosX(), Hit3.getPosY() - Hit2.getPosY(), Hit3.getPosZ() - Hit2.getPosZ());
  TVector3 crossProd = vec1.Cross(vec2);

  double Dcoeef = -crossProd(0) * Hit2.getPosX() - crossProd(1) * Hit2.getPosY() - crossProd(2) * Hit2.getPosZ();

  double distancefromZero = fabs(crossProd(0) * Center.X() + crossProd(1) * Center.Y() + crossProd(2) * Center.Z() + Dcoeef) / crossProd.Mag();
  return distancefromZero;

}






double Calc3DAngleHit(const JPetHit & Hit1,const JPetHit & Hit2 , TVector3 Center) {

  double scalarProd = (Hit1.getPosX() - Center.X()) * (Hit2.getPosX() - Center.X()) + (Hit1.getPosY() - Center.Y()) * (Hit2.getPosY() - Center.Y()) + (Hit1.getPosZ() - Center.Z()) * (Hit2.getPosZ() - Center.Z());

  double magProd = sqrt((pow(Hit1.getPosX() - Center.X(), 2) +
      pow(Hit1.getPosY() - Center.Y(), 2) +
      pow(Hit1.getPosZ() - Center.Z(), 2)) *
    (pow(Hit2.getPosX() - Center.X(), 2) +
      pow(Hit2.getPosY() - Center.Y(), 2) +
      pow(Hit2.getPosZ() - Center.Z(), 2)));

  double Angle = acos(scalarProd / magProd) * 180 / 3.14159;

  return Angle;
}


double Calc3DAngleMC(const JPetMCHit & Hit1,const JPetMCHit & Hit2 , TVector3 Center) {

  double scalarProd = (Hit1.getPosX() - Center.X()) * (Hit2.getPosX() - Center.X()) + (Hit1.getPosY() - Center.Y()) * (Hit2.getPosY() - Center.Y()) + (Hit1.getPosZ() - Center.Z()) * (Hit2.getPosZ() - Center.Z());

  double magProd = sqrt((pow(Hit1.getPosX() - Center.X(), 2) +
      pow(Hit1.getPosY() - Center.Y(), 2) +
      pow(Hit1.getPosZ() - Center.Z(), 2)) *
    (pow(Hit2.getPosX() - Center.X(), 2) +
      pow(Hit2.getPosY() - Center.Y(), 2) +
      pow(Hit2.getPosZ() - Center.Z(), 2)));

  double Angle = acos(scalarProd / magProd) * 180 / 3.14159;

  return Angle;
}






bool EventCategorizer::exec()
{
  //if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
  const JPetTimeWindowMC* time_window_mc = nullptr;
  if (time_window_mc = dynamic_cast<const JPetTimeWindowMC*>(fEvent)) {
    fIsMC = true;    
  }

  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    
    for (uint i = 0; i < timeWindow->getNumberOfEvents(); i++) {
     const auto& event = dynamic_cast<const JPetEvent&>(timeWindow->operator[](i));

	//Juhi altered analysis//

	
	if(fIsMC){

	int Multiplicity = event.getHits().size();
	getStatistics().getHisto1D("Multiplicity")->Fill(Multiplicity);				//Hit Multiplicity

	

	if( event.getHits().size() == 4 ){
        getStatistics().getHisto1D("Multiplicity_4")->Fill(Multiplicity);			//Select events with only 4-Hits



	auto mc_hit_1 = time_window_mc->getMCHit<JPetMCHit>(event.getHits().at(0).getMCindex());
	auto mc_hit_2 = time_window_mc->getMCHit<JPetMCHit>(event.getHits().at(1).getMCindex());
	auto mc_hit_3 = time_window_mc->getMCHit<JPetMCHit>(event.getHits().at(2).getMCindex());
	auto mc_hit_4 = time_window_mc->getMCHit<JPetMCHit>(event.getHits().at(3).getMCindex());


	int Multiplicity_1 = mc_hit_1.getGenGammaMultiplicity();
	int Multiplicity_2 = mc_hit_2.getGenGammaMultiplicity();
	int Multiplicity_3 = mc_hit_3.getGenGammaMultiplicity();
	int Multiplicity_4 = mc_hit_4.getGenGammaMultiplicity();
		

	
	getStatistics().getHisto1D("True_MC_Multiplicity")->Fill(Multiplicity_1);		//Check the MC true multiplicity flag
	getStatistics().getHisto1D("True_MC_Multiplicity")->Fill(Multiplicity_2);		//Check the MC true multiplicity flag
	getStatistics().getHisto1D("True_MC_Multiplicity")->Fill(Multiplicity_3);		//Check the MC true multiplicity flag
	getStatistics().getHisto1D("True_MC_Multiplicity")->Fill(Multiplicity_4);		//Check the MC true multiplicity flag



	std::vector<int> Multi;
	Multi.push_back(Multiplicity_1);		
	Multi.push_back(Multiplicity_2);
	Multi.push_back(Multiplicity_3);
	Multi.push_back(Multiplicity_4);

	//Sort and seletect only events with two '2' flag candidates and two with '102' flag candidate//
 	
	int count_p = 0;
	int count_s = 0;

	for(int m=0; m<Multi.size(); m++)
		{


		if(Multi.at(m)==3)
		{
			count_p++;

				}

		else if(Multi.at(m)==103)
		{

			count_s++;
		
			}


		}
	
	

	if(count_p==3 && count_s==1)
	{

	getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->Fill(Multiplicity_1);		//Cut on the MC true multiplicity
	getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->Fill(Multiplicity_2);		//Cut on the MC true multiplicity
	getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->Fill(Multiplicity_3);		//Cut on the MC true multiplicity
	getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->Fill(Multiplicity_4);		//Cut on the MC true multiplicity
 	getStatistics().getHisto1D("Multiplicity_Gen_Mul_Cut")->Fill(Multiplicity);


	//std::vector <pair<int, JPetMCHit>> Hit_vector;

	std::vector <pair<double,pair<JPetHit,JPetMCHit>>> Hit_vector;

	Hit_vector.push_back({Multiplicity_1,{event.getHits().at(0), mc_hit_1}});
	Hit_vector.push_back({Multiplicity_2,{event.getHits().at(1), mc_hit_2}});
	Hit_vector.push_back({Multiplicity_3,{event.getHits().at(2), mc_hit_3}});
	Hit_vector.push_back({Multiplicity_4,{event.getHits().at(3), mc_hit_4}});


	std::sort(Hit_vector.begin(), Hit_vector.end(), comparison3);

	//Control Spectra//
	
	//Smeared Energy
	getStatistics().getHisto1D("Energy_Primary_hit")->Fill(Hit_vector.at(0).second.first.getEnergy());		//Energy of the Primary_1
	getStatistics().getHisto1D("Energy_Primary_hit")->Fill(Hit_vector.at(1).second.first.getEnergy());		//Energy of the Primary_2
	getStatistics().getHisto1D("Energy_Primary_hit")->Fill(Hit_vector.at(2).second.first.getEnergy());		//Energy of the Primary_3
	getStatistics().getHisto1D("Energy_Secondary_hit")->Fill(Hit_vector.at(3).second.first.getEnergy());		//Energy of the Secondary
	
	//True Energy
	getStatistics().getHisto1D("Energy_Primary_mc_hit")->Fill(Hit_vector.at(0).second.second.getEnergy());		//Energy of the Primary_1
	getStatistics().getHisto1D("Energy_Primary_mc_hit")->Fill(Hit_vector.at(1).second.second.getEnergy());		//Energy of the Primary_2
	getStatistics().getHisto1D("Energy_Primary_mc_hit")->Fill(Hit_vector.at(2).second.second.getEnergy());		//Energy of the Primary_3
	getStatistics().getHisto1D("Energy_Secondary_mc_hit")->Fill(Hit_vector.at(3).second.second.getEnergy());	//Energy of the Secondary


	//Smeared Z_Position
	getStatistics().getHisto1D("ZPosition_Primary_hit")->Fill(Hit_vector.at(0).second.first.getPosZ());		//Energy of the Primary_1
	getStatistics().getHisto1D("ZPosition_Primary_hit")->Fill(Hit_vector.at(1).second.first.getPosZ());		//Energy of the Primary_2
	getStatistics().getHisto1D("ZPosition_Primary_hit")->Fill(Hit_vector.at(2).second.first.getPosZ());		//Energy of the Primary_3
	getStatistics().getHisto1D("ZPosition_Secondary_hit")->Fill(Hit_vector.at(3).second.first.getPosZ());		//Energy of the Secondary



	//True Z_Position
	getStatistics().getHisto1D("ZPosition_Primary_mc_hit")->Fill(Hit_vector.at(0).second.second.getPosZ());		//Energy of the Primary_1
	getStatistics().getHisto1D("ZPosition_Primary_mc_hit")->Fill(Hit_vector.at(1).second.second.getPosZ());		//Energy of the Primary_2
	getStatistics().getHisto1D("ZPosition_Primary_mc_hit")->Fill(Hit_vector.at(2).second.second.getPosZ());		//Energy of the Primary_3
	getStatistics().getHisto1D("ZPosition_Secondary_mc_hit")->Fill(Hit_vector.at(3).second.second.getPosZ());	//Energy of the Secondary


	//Z Resolution Check
	getStatistics().getHisto1D("Z_resolution")->Fill(Hit_vector.at(0).second.first.getPosZ() - Hit_vector.at(0).second.second.getPosZ());		//Resolution of the Primary_1
	getStatistics().getHisto1D("Z_resolution")->Fill(Hit_vector.at(1).second.first.getPosZ() - Hit_vector.at(1).second.second.getPosZ());		//Resolution of the Primary_2
	getStatistics().getHisto1D("Z_resolution")->Fill(Hit_vector.at(2).second.first.getPosZ() - Hit_vector.at(2).second.second.getPosZ());		//Resolution of the Primary_3
	getStatistics().getHisto1D("Z_resolution")->Fill(Hit_vector.at(3).second.first.getPosZ() - Hit_vector.at(3).second.second.getPosZ());		//Resolution of the Secondary


	//True MC Expectation Value//
	
	//Momentum Directions

	TVector3 k_0 = Hit_vector.at(0).second.second.getMomentum();
	TVector3 k_1 = Hit_vector.at(1).second.second.getMomentum();
	TVector3 k_2 = Hit_vector.at(2).second.second.getMomentum();
	TVector3 s_1 = Hit_vector.at(3).second.second.getMomentum();

	double Energy_k_0 = k_0.Mag();
	double Energy_k_1 = k_1.Mag();
	double Energy_k_2 = k_2.Mag();
	double Scatter_1 = s_1.Mag();
	
	double Energy_k_0_dep = Hit_vector.at(0).second.second.getEnergy();
	double Energy_k_1_dep = Hit_vector.at(1).second.second.getEnergy();
	double Energy_k_2_dep = Hit_vector.at(2).second.second.getEnergy();
	double Scatter_1_dep = Hit_vector.at(3).second.second.getEnergy();

	
	//Angle between the primary photons - true info//
	double Angle01 = k_0.Angle(k_1)*TMath::RadToDeg();
	double Angle12 = k_1.Angle(k_2)*TMath::RadToDeg();
	double Angle20 = k_2.Angle(k_0)*TMath::RadToDeg();

	getStatistics().getHisto2D("Angle_3D_True")->Fill(Angle01, Angle12);		//Angle 3D
	getStatistics().getHisto2D("Angle_3D_True")->Fill(Angle12, Angle20);		//Angle 3D
	getStatistics().getHisto2D("Angle_3D_True")->Fill(Angle20, Angle01);		//Angle 3D

	//Angle between the primary photons - hit smeared//
	double Angle_01_hit = Calc3DAngleHit(Hit_vector.at(0).second.first, Hit_vector.at(1).second.first, Center); 
	double Angle_12_hit = Calc3DAngleHit(Hit_vector.at(1).second.first, Hit_vector.at(2).second.first, Center); 
	double Angle_20_hit = Calc3DAngleHit(Hit_vector.at(2).second.first, Hit_vector.at(0).second.first, Center); 

	getStatistics().getHisto2D("Angle_3D_Hit")->Fill(Angle_01_hit, Angle_12_hit);		//Angle 3D
	getStatistics().getHisto2D("Angle_3D_Hit")->Fill(Angle_12_hit, Angle_20_hit);		//Angle 3D
	getStatistics().getHisto2D("Angle_3D_Hit")->Fill(Angle_20_hit, Angle_01_hit);		//Angle 3D

	
	getStatistics().getHisto1D("Energy_all_primary")->Fill(Energy_k_0);		//Energy all primary
	getStatistics().getHisto1D("Energy_all_primary")->Fill(Energy_k_1);		//Energy all primary
	getStatistics().getHisto1D("Energy_all_primary")->Fill(Energy_k_2);		//Energy all primary
	getStatistics().getHisto1D("Energy_secondary")->Fill(Scatter_1);		//Energy all primary


	//Sorting the primary photons based on energy//

	//Energy of true hit vector
	std::vector <pair<double,pair<JPetHit,JPetMCHit>>> Energy_vector_mc_hit;
	Energy_vector_mc_hit.push_back({Energy_k_0,{Hit_vector.at(0).second.first, Hit_vector.at(0).second.second}});
	Energy_vector_mc_hit.push_back({Energy_k_1,{Hit_vector.at(1).second.first, Hit_vector.at(1).second.second}});
	Energy_vector_mc_hit.push_back({Energy_k_2,{Hit_vector.at(2).second.first, Hit_vector.at(2).second.second}});
	std::sort(Energy_vector_mc_hit.begin(), Energy_vector_mc_hit.end(), comparison3);	

	

	//Before doing the scatter test, we need the time, energy and X(), Y(), resolutions//

	//X() Resolution//
	getStatistics().getHisto1D("X_resolution")->Fill(Hit_vector.at(0).second.first.getPosX() - Hit_vector.at(0).second.second.getPosX());		//Resolution of the Primary_1
	getStatistics().getHisto1D("X_resolution")->Fill(Hit_vector.at(1).second.first.getPosX() - Hit_vector.at(1).second.second.getPosX());		//Resolution of the Primary_2
	getStatistics().getHisto1D("X_resolution")->Fill(Hit_vector.at(2).second.first.getPosX() - Hit_vector.at(2).second.second.getPosX());		//Resolution of the Primary_3
	getStatistics().getHisto1D("X_resolution")->Fill(Hit_vector.at(3).second.first.getPosX() - Hit_vector.at(3).second.second.getPosX());		//Resolution of the Secondary


	
	getStatistics().getHisto1D("Scintillator_ID")->Fill(Hit_vector.at(0).second.second.getScintillator().getID());
	getStatistics().getHisto1D("Scintillator_ID")->Fill(Hit_vector.at(1).second.second.getScintillator().getID());
	getStatistics().getHisto1D("Scintillator_ID")->Fill(Hit_vector.at(2).second.second.getScintillator().getID());
	getStatistics().getHisto1D("Scintillator_ID")->Fill(Hit_vector.at(3).second.second.getScintillator().getID());



	//Y() Resolution//
	getStatistics().getHisto1D("Y_resolution")->Fill(Hit_vector.at(0).second.first.getPosY() - Hit_vector.at(0).second.second.getPosY());		//Resolution of the Primary_1
	getStatistics().getHisto1D("Y_resolution")->Fill(Hit_vector.at(1).second.first.getPosY() - Hit_vector.at(1).second.second.getPosY());		//Resolution of the Primary_2
	getStatistics().getHisto1D("Y_resolution")->Fill(Hit_vector.at(2).second.first.getPosY() - Hit_vector.at(2).second.second.getPosY());		//Resolution of the Primary_3
	getStatistics().getHisto1D("Y_resolution")->Fill(Hit_vector.at(3).second.first.getPosY() - Hit_vector.at(3).second.second.getPosY());		//Resolution of the Secondary


	//Energy Resolution//
	getStatistics().getHisto1D("Energy_resolution")->Fill(Hit_vector.at(0).second.first.getEnergy() - Hit_vector.at(0).second.second.getEnergy());		//Resolution of the Primary_1
	getStatistics().getHisto1D("Energy_resolution")->Fill(Hit_vector.at(1).second.first.getEnergy() - Hit_vector.at(1).second.second.getEnergy());		//Resolution of the Primary_2
	getStatistics().getHisto1D("Energy_resolution")->Fill(Hit_vector.at(2).second.first.getEnergy() - Hit_vector.at(2).second.second.getEnergy());		//Resolution of the Primary_3
	getStatistics().getHisto1D("Energy_resolution")->Fill(Hit_vector.at(3).second.first.getEnergy() - Hit_vector.at(3).second.second.getEnergy());		//Resolution of the Secondary
	

       //Time Resolution//
	getStatistics().getHisto1D("Time_resolution")->Fill(((Hit_vector.at(0).second.first.getTime()*-1./100.) - (Hit_vector.at(0).second.second.getTime())));		//Resolution of the Primary_1
	getStatistics().getHisto1D("Time_resolution")->Fill(((Hit_vector.at(1).second.first.getTime()*-1./100.) - (Hit_vector.at(1).second.second.getTime())));		//Resolution of the Primary_2
	getStatistics().getHisto1D("Time_resolution")->Fill(((Hit_vector.at(2).second.first.getTime()*-1./100.) - (Hit_vector.at(2).second.second.getTime())));		//Resolution of the Primary_3
	getStatistics().getHisto1D("Time_resolution")->Fill(((Hit_vector.at(3).second.first.getTime()*-1./100.) - (Hit_vector.at(3).second.second.getTime())));		//Resolution of the Secondary
	



//	cout<<(Hit_vector.at(0).second.second.getTime())<<"	"<<(Hit_vector.at(0).second.first.getTime()*-1./100.)<<endl;
	//cout<<Energy_vector_mc_hit.at(0).first<<"	"<<Energy_vector_mc_hit.at(1).first<<"	"<<Energy_vector_mc_hit.at(2).first<<endl;
	//at(2) is the highest energy photon and at(0) is the lowest energy photon

	//TOF correcting times to see if scatter test is valid//

	double TOF_hit0 = CalcTOFHit(Hit_vector.at(0).second.first, Center);
	double TOF_hit1 = CalcTOFHit(Hit_vector.at(1).second.first, Center);
	double TOF_hit2 = CalcTOFHit(Hit_vector.at(2).second.first, Center);
	double TOF_hit3 = CalcTOFHit(Hit_vector.at(3).second.first, Center);







	//Scatter Test - With Hit smeared//
	
	double ScatterTest_2_Hit = CalScatterTestHit(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first);
	double ScatterTest_1_Hit = CalScatterTestHit(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first);
	double ScatterTest_0_Hit = CalScatterTestHit(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first);

	std::vector <pair<double,pair<JPetHit,int>>> Scatter_Vector_Hit;
	Scatter_Vector_Hit.push_back({ScatterTest_2_Hit,{Energy_vector_mc_hit.at(0).second.first, 2}});
	Scatter_Vector_Hit.push_back({ScatterTest_1_Hit,{Energy_vector_mc_hit.at(1).second.first, 1}});
	Scatter_Vector_Hit.push_back({ScatterTest_0_Hit,{Energy_vector_mc_hit.at(2).second.first, 0}});

	std::sort(Scatter_Vector_Hit.begin(), Scatter_Vector_Hit.end(), comparison_scatter_Hit);	
	//cout<<Scatter_Vector.at(0).first<<"	"<<Scatter_Vector.at(1).first<<"	"<<Scatter_Vector.at(2).first<<endl;
	//at(0) is the lowest modulus scatter value, at(2) is the largest modulus scatter value// 

	getStatistics().getHisto1D("Scatter_Hit_Lowest")->Fill(Scatter_Vector_Hit.at(0).first);
	getStatistics().getHisto2D("Scatter_Hit_21")->Fill(ScatterTest_2_Hit, ScatterTest_1_Hit);
	
	
	//Scatter Test - With True MC Hit//
	
	double ScatterTest_2_MC = CalScatterTestTrueMC(Energy_vector_mc_hit.at(2).second.second, Hit_vector.at(3).second.second);
	double ScatterTest_1_MC = CalScatterTestTrueMC(Energy_vector_mc_hit.at(1).second.second, Hit_vector.at(3).second.second);
	double ScatterTest_0_MC = CalScatterTestTrueMC(Energy_vector_mc_hit.at(0).second.second, Hit_vector.at(3).second.second);

	std::vector <pair<double,pair<JPetMCHit,int>>> Scatter_Vector_MC;
	Scatter_Vector_MC.push_back({ScatterTest_2_MC,{Energy_vector_mc_hit.at(0).second.second, 0}});
	Scatter_Vector_MC.push_back({ScatterTest_1_MC,{Energy_vector_mc_hit.at(1).second.second, 1}});
	Scatter_Vector_MC.push_back({ScatterTest_0_MC,{Energy_vector_mc_hit.at(2).second.second, 2}});

	std::sort(Scatter_Vector_MC.begin(), Scatter_Vector_MC.end(), comparison_scatter_MC);	
	//cout<<Scatter_Vector.at(0).first<<"	"<<Scatter_Vector.at(1).first<<"	"<<Scatter_Vector.at(2).first<<endl;
	//at(0) is the lowest modulus scatter value, at(2) is the largest modulus scatter value// 

	getStatistics().getHisto1D("Scatter_MC_Lowest")->Fill(Scatter_Vector_MC.at(0).first);
	getStatistics().getHisto2D("Scatter_MC_21")->Fill(ScatterTest_2_MC, ScatterTest_1_MC);


	if(Scatter_Vector_MC.at(0).first <= 0.5 && Scatter_Vector_MC.at(0).first >= -0.5)
	{

	getStatistics().getHisto1D("Scatter_MC_Cut")->Fill(Scatter_Vector_MC.at(0).first);
	getStatistics().getHisto1D("Scatter_Hit_Cut")->Fill(Scatter_Vector_Hit.at(0).first);
	

	//-----Final set of hits with their corresponding scatter hit----//
	
	//Energy_vector_mc_hit.at(0).second.first ----- hit of least energetic primary
	//Energy_vector_mc_hit.at(0).second.second -----mchit of least energetic primary
	
	//Energy_vector_mc_hit.at(1).second.first ----- hit of second most energetic primary
	//Energy_vector_mc_hit.at(1).second.second -----mchit of second most energetic primary

	//Energy_vector_mc_hit.at(2).second.first ----- hit of most energetic primary
	//Energy_vector_mc_hit.at(2).second.second -----mchit of most energetic primary

	//Hit_vector.at(3).second.first ----- hit of scatter photon
	//Hit_vector.at(3).second.second ----- mchit of scatter photon

	//Scatter_Vector_Hit.at(0).second.first ----- the primary it belong's to
	//Scatter_Vector_Hit.at(0).second.second ----- the flag (0 is least and 2 is most)

	//---------------------------------------------------------------//


	//Create a vector to fill//



	//------Expectation Value -- loop bit-----//

	

	if(Scatter_Vector_Hit.at(0).second.second == 2)
	{


	double E2K1_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(1).second.first);
	double E2K1_MC = CalExpecValueMC(Energy_vector_mc_hit.at(2).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(1).second.second);
		getStatistics().getHisto1D("EV_Hit")->Fill(E2K1_Hit);
		getStatistics().getHisto1D("EV_MC")->Fill(E2K1_MC);


		getStatistics().getHisto1D("E2K1_Hit")->Fill(E2K1_Hit);
		getStatistics().getHisto1D("E2K1_MC")->Fill(E2K1_MC);

	//Expectation Value Resolution//
		getStatistics().getHisto1D("EV_Resolution")->Fill(E2K1_Hit - E2K1_MC);



	double E2K0_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(0).second.first);
	double E2K0_MC = CalExpecValueMC(Energy_vector_mc_hit.at(2).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(0).second.second);

		getStatistics().getHisto1D("E2K0_Hit")->Fill(E2K0_Hit);
		getStatistics().getHisto1D("E2K0_MC")->Fill(E2K0_MC);


	//Polarization Magnitude Check//

	double polarization_mc =  Energy_vector_mc_hit.at(2).second.second.getPolarization().Mag();
	double pr_mc = Energy_vector_mc_hit.at(2).second.second.getPos().Mag();
	double sr_mc = Hit_vector.at(3).second.second.getPos().Mag();
	double pol_true_mc =  (polarization_mc)/(pr_mc*sr_mc);
	getStatistics().getHisto1D("pol_true")->Fill(pol_true_mc);
	double polarization_reco = ((Energy_vector_mc_hit.at(2).second.second.getPos()).Cross(Hit_vector.at(3).second.second.getPos())).Mag();
	double pol_reco_mc =  (polarization_reco)/(pr_mc*sr_mc);
	getStatistics().getHisto1D("pol_reco")->Fill(pol_reco_mc);
	getStatistics().getHisto2D("Polarization_Check")->Fill(pol_reco_mc, pol_true_mc);



	//Polarization Vector Angle Check//

	TVector3 poltrue = Energy_vector_mc_hit.at(2).second.second.getPolarization();
	TVector3 momentum_pr = Energy_vector_mc_hit.at(2).second.second.getMomentum();
	TVector3 pol_reco = (Energy_vector_mc_hit.at(2).second.second.getMomentum()).Cross(Hit_vector.at(3).second.second.getMomentum());

	double True_Angle = poltrue.Angle(momentum_pr)*TMath::RadToDeg();
	double Reco_Angle = pol_reco.Angle(momentum_pr)*TMath::RadToDeg();

	getStatistics().getHisto2D("Polarization_Angle_Check")->Fill(Reco_Angle, True_Angle);
	//cout<<True_Angle<<"	"<<Reco_Angle<<endl;

	}

	else if(Scatter_Vector_Hit.at(0).second.second == 1)
	{


	double E1K0_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(0).second.first);
	double E1K0_MC = CalExpecValueMC(Energy_vector_mc_hit.at(1).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(0).second.second);

	getStatistics().getHisto1D("EV_Hit")->Fill(E1K0_Hit);
	getStatistics().getHisto1D("EV_MC")->Fill(E1K0_MC);

	getStatistics().getHisto1D("E1K0_Hit")->Fill(E1K0_Hit);
	getStatistics().getHisto1D("E1K0_MC")->Fill(E1K0_MC);


		//Expectation Value Resolution//
		getStatistics().getHisto1D("EV_Resolution")->Fill(E1K0_Hit - E1K0_MC);


	double E1K2_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(2).second.first);
	double E1K2_MC = CalExpecValueMC(Energy_vector_mc_hit.at(1).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(2).second.second);

	getStatistics().getHisto1D("E1K2_Hit")->Fill(E1K2_Hit);
	getStatistics().getHisto1D("E1K2_MC")->Fill(E1K2_MC);



	}

	else if(Scatter_Vector_Hit.at(0).second.second == 0)
	{


	double E0K2_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(2).second.first);
	double E0K2_MC = CalExpecValueMC(Energy_vector_mc_hit.at(0).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(2).second.second);

	
	getStatistics().getHisto1D("EV_Hit")->Fill(E0K2_Hit);
	getStatistics().getHisto1D("EV_MC")->Fill(E0K2_MC);

	getStatistics().getHisto1D("E0K2_Hit")->Fill(E0K2_Hit);
	getStatistics().getHisto1D("E0K2_MC")->Fill(E0K2_MC);

		//Expectation Value Resolution//
		getStatistics().getHisto1D("EV_Resolution")->Fill(E0K2_Hit - E0K2_MC);



        double E0K1_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(1).second.first);
	double E0K1_MC = CalExpecValueMC(Energy_vector_mc_hit.at(0).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(1).second.second);

	getStatistics().getHisto1D("E0K1_Hit")->Fill(E0K1_Hit);
	getStatistics().getHisto1D("E0K1_MC")->Fill(E0K1_MC);



	}

	//---------------------------------------//



	
	//First Cut - Z_Position//
	
	getStatistics().getHisto1D("Z_position_before")->Fill(Energy_vector_mc_hit.at(0).second.first.getPosZ());
	getStatistics().getHisto1D("Z_position_before")->Fill(Energy_vector_mc_hit.at(1).second.first.getPosZ());
	getStatistics().getHisto1D("Z_position_before")->Fill(Energy_vector_mc_hit.at(2).second.first.getPosZ());
	getStatistics().getHisto1D("Z_position_before")->Fill(Hit_vector.at(3).second.first.getPosZ());


	if(abs(Energy_vector_mc_hit.at(0).second.first.getPosZ()) <= 22. && abs(Energy_vector_mc_hit.at(1).second.first.getPosZ()) <= 22. && abs(Energy_vector_mc_hit.at(2).second.first.getPosZ()) <= 22. && abs(Hit_vector.at(3).second.first.getPosZ()) <= 22.)
	{

	getStatistics().getHisto1D("Z_position_after")->Fill(Energy_vector_mc_hit.at(0).second.first.getPosZ());
	getStatistics().getHisto1D("Z_position_after")->Fill(Energy_vector_mc_hit.at(1).second.first.getPosZ());
	getStatistics().getHisto1D("Z_position_after")->Fill(Energy_vector_mc_hit.at(2).second.first.getPosZ());
	getStatistics().getHisto1D("Z_position_after")->Fill(Hit_vector.at(3).second.first.getPosZ());



	
	//------Expectation Value -- loop bit after Z_Position Cut-----//

	if(Scatter_Vector_Hit.at(0).second.second == 2)
	{


	double E2K1_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(1).second.first);
	double E2K1_MC = CalExpecValueMC(Energy_vector_mc_hit.at(2).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(1).second.second);
		getStatistics().getHisto1D("EV_Hit_Z_Pos")->Fill(E2K1_Hit);
		getStatistics().getHisto1D("EV_MC_Z_Pos")->Fill(E2K1_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 1)
	{


	double E1K0_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(0).second.first);
	double E1K0_MC = CalExpecValueMC(Energy_vector_mc_hit.at(1).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(0).second.second);

	getStatistics().getHisto1D("EV_Hit_Z_Pos")->Fill(E1K0_Hit);
	getStatistics().getHisto1D("EV_MC_Z_Pos")->Fill(E1K0_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 0)
	{


	double E0K2_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(2).second.first);
	double E0K2_MC = CalExpecValueMC(Energy_vector_mc_hit.at(0).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(2).second.second);

	
	getStatistics().getHisto1D("EV_Hit_Z_Pos")->Fill(E0K2_Hit);
	getStatistics().getHisto1D("EV_MC_Z_Pos")->Fill(E0K2_MC);



	}

	//---------------------------------------//



	//Second Cut - Energy Cut//

	getStatistics().getHisto1D("Energy_before")->Fill(Energy_vector_mc_hit.at(0).second.first.getEnergy());
	getStatistics().getHisto1D("Energy_before")->Fill(Energy_vector_mc_hit.at(1).second.first.getEnergy());
	getStatistics().getHisto1D("Energy_before")->Fill(Energy_vector_mc_hit.at(2).second.first.getEnergy());
	getStatistics().getHisto1D("Energy_before")->Fill(Hit_vector.at(3).second.first.getEnergy());

	if(abs(Energy_vector_mc_hit.at(0).second.first.getEnergy()) <= 340. && abs(Energy_vector_mc_hit.at(1).second.first.getEnergy()) <= 340. && abs(Energy_vector_mc_hit.at(2).second.first.getEnergy()) <= 340. && abs(Hit_vector.at(3).second.first.getEnergy()) <= 340.)
	{


	getStatistics().getHisto1D("Energy_after")->Fill(Energy_vector_mc_hit.at(0).second.first.getEnergy());
	getStatistics().getHisto1D("Energy_after")->Fill(Energy_vector_mc_hit.at(1).second.first.getEnergy());
	getStatistics().getHisto1D("Energy_after")->Fill(Energy_vector_mc_hit.at(2).second.first.getEnergy());
	getStatistics().getHisto1D("Energy_after")->Fill(Hit_vector.at(3).second.first.getEnergy());



	//------Expectation Value -- loop bit after Energy Cut-----//

	if(Scatter_Vector_Hit.at(0).second.second == 2)
	{


	double E2K1_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(1).second.first);
	double E2K1_MC = CalExpecValueMC(Energy_vector_mc_hit.at(2).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(1).second.second);
		getStatistics().getHisto1D("EV_Hit_Energy")->Fill(E2K1_Hit);
		getStatistics().getHisto1D("EV_MC_Energy")->Fill(E2K1_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 1)
	{


	double E1K0_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(0).second.first);
	double E1K0_MC = CalExpecValueMC(Energy_vector_mc_hit.at(1).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(0).second.second);

	getStatistics().getHisto1D("EV_Hit_Energy")->Fill(E1K0_Hit);
	getStatistics().getHisto1D("EV_MC_Energy")->Fill(E1K0_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 0)
	{


	double E0K2_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(2).second.first);
	double E0K2_MC = CalExpecValueMC(Energy_vector_mc_hit.at(0).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(2).second.second);

	
	getStatistics().getHisto1D("EV_Hit_Energy")->Fill(E0K2_Hit);
	getStatistics().getHisto1D("EV_MC_Energy")->Fill(E0K2_MC);



	}

	//---------------------------------------//


	//Third Cut - Emission Time Cut//

	double TOF_2 = CalcTOFHit(Energy_vector_mc_hit.at(2).second.first, Center);
	double TOF_1 = CalcTOFHit(Energy_vector_mc_hit.at(1).second.first, Center);
	double TOF_0 = CalcTOFHit(Energy_vector_mc_hit.at(0).second.first, Center);

	double TOF_2_MC = CalcTOFMC(Energy_vector_mc_hit.at(2).second.second, Center);
	double TOF_1_MC = CalcTOFMC(Energy_vector_mc_hit.at(1).second.second, Center);
	double TOF_0_MC = CalcTOFMC(Energy_vector_mc_hit.at(0).second.second, Center);

		vector <pair< double, JPetHit >> TOF_Vector;
		vector <pair< double, JPetMCHit >> TOF_Vector_MC;

	TOF_Vector.push_back({TOF_2, Energy_vector_mc_hit.at(2).second.first});	
	TOF_Vector.push_back({TOF_1, Energy_vector_mc_hit.at(1).second.first});	
	TOF_Vector.push_back({TOF_0, Energy_vector_mc_hit.at(0).second.first});	

	TOF_Vector_MC.push_back({TOF_2_MC, Energy_vector_mc_hit.at(2).second.second});	
	TOF_Vector_MC.push_back({TOF_1_MC, Energy_vector_mc_hit.at(1).second.second});	
	TOF_Vector_MC.push_back({TOF_0_MC, Energy_vector_mc_hit.at(0).second.second});	


		std::sort(TOF_Vector.begin(), TOF_Vector.end(), comparison);
		std::sort(TOF_Vector_MC.begin(), TOF_Vector_MC.end(), comparison_MC);

	double TOF = (TOF_Vector.at(2).first - TOF_Vector.at(0).first);
	double TOF_MC = (TOF_Vector_MC.at(2).first - TOF_Vector_MC.at(0).first);

	double TOF_resolution = (TOF - TOF_MC);

	getStatistics().getHisto1D("TOF_Resolution")->Fill(TOF_resolution);

	getStatistics().getHisto1D("TOF_before")->Fill(TOF);

	if((TOF_Vector.at(2).first - TOF_Vector.at(0).first) <= 2.)
	{

	getStatistics().getHisto1D("TOF_after")->Fill(TOF);


       
	//------Expectation Value -- loop bit after TOF Cut-----//

	if(Scatter_Vector_Hit.at(0).second.second == 2)
	{


	double E2K1_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(1).second.first);
	double E2K1_MC = CalExpecValueMC(Energy_vector_mc_hit.at(2).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(1).second.second);
		getStatistics().getHisto1D("EV_Hit_TOF")->Fill(E2K1_Hit);
		getStatistics().getHisto1D("EV_MC_TOF")->Fill(E2K1_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 1)
	{


	double E1K0_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(0).second.first);
	double E1K0_MC = CalExpecValueMC(Energy_vector_mc_hit.at(1).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(0).second.second);

	getStatistics().getHisto1D("EV_Hit_TOF")->Fill(E1K0_Hit);
	getStatistics().getHisto1D("EV_MC_TOF")->Fill(E1K0_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 0)
	{


	double E0K2_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(2).second.first);
	double E0K2_MC = CalExpecValueMC(Energy_vector_mc_hit.at(0).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(2).second.second);

	
	getStatistics().getHisto1D("EV_Hit_TOF")->Fill(E0K2_Hit);
	getStatistics().getHisto1D("EV_MC_TOF")->Fill(E0K2_MC);



	}

	//---------------------------------------//


	//Fourth Cut - Distance Plane Cut//


	double Distance = CalDistofSurfaceHit(Energy_vector_mc_hit.at(0).second.first, Energy_vector_mc_hit.at(1).second.first, Energy_vector_mc_hit.at(2).second.first, Center);


	double Distance_MC = CalDistofSurfaceMC(Energy_vector_mc_hit.at(0).second.second, Energy_vector_mc_hit.at(1).second.second, Energy_vector_mc_hit.at(2).second.second, Center);

	getStatistics().getHisto1D("Distance_Resolution")->Fill(Distance - Distance_MC);	
	getStatistics().getHisto1D("Distance_Surface_before")->Fill(Distance);

	if(Distance <= 3.)
	{
	
	getStatistics().getHisto1D("Distance_Surface_after")->Fill(Distance);

	//------Expectation Value -- loop bit after Distance Cut-----//

	if(Scatter_Vector_Hit.at(0).second.second == 2)
	{


	double E2K1_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(1).second.first);
	double E2K1_MC = CalExpecValueMC(Energy_vector_mc_hit.at(2).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(1).second.second);
		getStatistics().getHisto1D("EV_Hit_Distance")->Fill(E2K1_Hit);
		getStatistics().getHisto1D("EV_MC_Distance")->Fill(E2K1_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 1)
	{


	double E1K0_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(0).second.first);
	double E1K0_MC = CalExpecValueMC(Energy_vector_mc_hit.at(1).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(0).second.second);

	getStatistics().getHisto1D("EV_Hit_Distance")->Fill(E1K0_Hit);
	getStatistics().getHisto1D("EV_MC_Distance")->Fill(E1K0_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 0)
	{


	double E0K2_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(2).second.first);
	double E0K2_MC = CalExpecValueMC(Energy_vector_mc_hit.at(0).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(2).second.second);

	
	getStatistics().getHisto1D("EV_Hit_Distance")->Fill(E0K2_Hit);
	getStatistics().getHisto1D("EV_MC_Distance")->Fill(E0K2_MC);



	}

	//---------------------------------------//



	//Fifth Cut - Sum Angle 3D Cut//	

	double Angle_01 = Calc3DAngleHit(Energy_vector_mc_hit.at(0).second.first, Energy_vector_mc_hit.at(1).second.first, Center);
	double Angle_12 = Calc3DAngleHit(Energy_vector_mc_hit.at(1).second.first, Energy_vector_mc_hit.at(2).second.first, Center);
	double Angle_20 = Calc3DAngleHit(Energy_vector_mc_hit.at(2).second.first, Energy_vector_mc_hit.at(0).second.first, Center);

	double Angle_01_MC = Calc3DAngleMC(Energy_vector_mc_hit.at(0).second.second, Energy_vector_mc_hit.at(1).second.second, Center);
	double Angle_12_MC = Calc3DAngleMC(Energy_vector_mc_hit.at(1).second.second, Energy_vector_mc_hit.at(2).second.second, Center);
	double Angle_20_MC = Calc3DAngleMC(Energy_vector_mc_hit.at(2).second.second, Energy_vector_mc_hit.at(0).second.second, Center);

	getStatistics().getHisto1D("Angle_Resolution")->Fill(Angle_01-Angle_01_MC);
	getStatistics().getHisto1D("Angle_Resolution")->Fill(Angle_12-Angle_12_MC);
	getStatistics().getHisto1D("Angle_Resolution")->Fill(Angle_20-Angle_20_MC);

	vector <double> Angle_Vector;

	Angle_Vector.push_back(Angle_01);
	Angle_Vector.push_back(Angle_12);
	Angle_Vector.push_back(Angle_20);	
	
	sort( Angle_Vector.begin(), Angle_Vector.end() );


	//cout<<Angle_Vector.at(0)<<"	"<<Angle_Vector.at(1)<<"	"<<Angle_Vector.at(2)<<endl;
	//at(0) is smallest Angle and at(2) is the largest angle//
		
	getStatistics().getHisto2D("Angle_Sum_before")->Fill((Angle_Vector.at(1)+Angle_Vector.at(0)), (Angle_Vector.at(1)-Angle_Vector.at(0)));		
	

	if((Angle_Vector.at(1)+Angle_Vector.at(0)) >= 190.)
	{ 
	
	getStatistics().getHisto2D("Angle_Sum_after")->Fill((Angle_Vector.at(1)+Angle_Vector.at(0)), (Angle_Vector.at(1)-Angle_Vector.at(0)));		
         getStatistics().getHisto2D("Angle_3D_True_Cut")->Fill(Angle01, Angle12);            //Angle 3D
         getStatistics().getHisto2D("Angle_3D_True_Cut")->Fill(Angle12, Angle20);            //Angle 3D
         getStatistics().getHisto2D("Angle_3D_True_Cut")->Fill(Angle20, Angle01);            //Angle 3D

	 getStatistics().getHisto2D("Angle_3D_Hit_Cut")->Fill(Angle_01_hit, Angle_12_hit);           //Angle 3D
         getStatistics().getHisto2D("Angle_3D_Hit_Cut")->Fill(Angle_12_hit, Angle_20_hit);           //Angle 3D
         getStatistics().getHisto2D("Angle_3D_Hit_Cut")->Fill(Angle_20_hit, Angle_01_hit);           //Angle 3D





	
	//------Expectation Value -- loop bit after Angle Cut-----//

	if(Scatter_Vector_Hit.at(0).second.second == 2)
	{


	double E2K1_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(1).second.first);
	double E2K1_MC = CalExpecValueMC(Energy_vector_mc_hit.at(2).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(1).second.second);
		getStatistics().getHisto1D("EV_Hit_Angle")->Fill(E2K1_Hit);
		getStatistics().getHisto1D("EV_MC_Angle")->Fill(E2K1_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 1)
	{


	double E1K0_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(0).second.first);
	double E1K0_MC = CalExpecValueMC(Energy_vector_mc_hit.at(1).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(0).second.second);

	getStatistics().getHisto1D("EV_Hit_Angle")->Fill(E1K0_Hit);
	getStatistics().getHisto1D("EV_MC_Angle")->Fill(E1K0_MC);


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 0)
	{


	double E0K2_Hit = CalExpecValueHit(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(2).second.first);
	double E0K2_MC = CalExpecValueMC(Energy_vector_mc_hit.at(0).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(2).second.second);

	
	getStatistics().getHisto1D("EV_Hit_Angle")->Fill(E0K2_Hit);
	getStatistics().getHisto1D("EV_MC_Angle")->Fill(E0K2_MC);



	}

	//---------------------------------------//










												}//Loop of Sum Angle



											}//Loop of the Distance



										}//Loop for the Emission Time Cut 

										
										
									}//Loop for the Energy Cut




								    }//Loop for the Z_Position Cut
						
									
									


								} //Loop for Scatter Test

								



							} //True Multiplicity Cut


						}


					}
      





    }


    
  } else { return false; }
  return true;



}

bool EventCategorizer::terminate()
{
  INFO("Event analysis completed.");
  return true;
}

void EventCategorizer::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events) { fOutputEvents->add<JPetEvent>(event); }
}

void EventCategorizer::initialiseHistograms(){
}




