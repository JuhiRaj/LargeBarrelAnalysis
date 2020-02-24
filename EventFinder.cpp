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
 *  @file EventFinder.cpp
 */

using namespace std;

#include <JPetOptionsTools/JPetOptionsTools.h>
#include <JPetWriter/JPetWriter.h>
#include <JPetMCHit/JPetMCHit.h>
#include "EventFinder.h"
#include <iostream>
#include "TMath.h"
#include "math.h"
#include "TRandom3.h"

using namespace jpet_options_tools;

EventFinder::EventFinder(const char* name): JPetUserTask(name) {}

EventFinder::~EventFinder() {}

bool EventFinder::init()
{
  INFO("Event finding started.");

  fOutputEvents = new JPetTimeWindow("JPetEvent");

  // Reading values from the user options if available
  // Getting bool for using corrupted hits
  if (isOptionSet(fParams.getOptions(), kUseCorruptedHitsParamKey)) {
    fUseCorruptedHits = getOptionAsBool(fParams.getOptions(), kUseCorruptedHitsParamKey);
    if(fUseCorruptedHits){
      WARNING("Event Finder is using Corrupted Hits, as set by the user");
    } else{
      WARNING("Event Finder is NOT using Corrupted Hits, as set by the user");
    }
  } else {
    WARNING("Event Finder is not using Corrupted Hits (default option)");
  }
  // Event time window
  if (isOptionSet(fParams.getOptions(), kEventTimeParamKey)) {
    fEventTimeWindow = getOptionAsFloat(fParams.getOptions(), kEventTimeParamKey);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %lf.",
      kEventTimeParamKey.c_str(), fEventTimeWindow
    ));
  }
  // Minimum number of hits in an event to save an event
  if (isOptionSet(fParams.getOptions(), kEventMinMultiplicity)) {
    fMinMultiplicity = getOptionAsInt(fParams.getOptions(), kEventMinMultiplicity);
  } else {
    WARNING(Form(
      "No value of the %s parameter provided by the user. Using default value of %d.",
      kEventMinMultiplicity.c_str(), fMinMultiplicity
    ));
  }
  // Getting bool for saving histograms
  if (isOptionSet(fParams.getOptions(), kSaveControlHistosParamKey)){
    fSaveControlHistos = getOptionAsBool(fParams.getOptions(), kSaveControlHistosParamKey);
  }

  // Initialize histograms
  if (fSaveControlHistos) { initialiseHistograms(); }
  return true;
}

//Functions for the check by Juhi//




bool EventFinder::exec()
{
  if (auto timeWindow = dynamic_cast<const JPetTimeWindow* const>(fEvent)) {
    saveEvents(buildEvents(*timeWindow));
  } else { return false; }
  return true;
}

bool EventFinder::terminate()
{
  INFO("Event fiding ended.");
  return true;
}

void EventFinder::saveEvents(const vector<JPetEvent>& events)
{
  for (const auto& event : events){
    fOutputEvents->add<JPetEvent>(event);
  }
}

/**
 * Main method of building Events - Hit in the Time slot are groupped
 * within time parameter, that can be set by the user
 */

	TVector3 Center_EF(0.0, 0.0, 0.0);  //Center of the Geometry


Bool_t comparison_gen_mul(const pair < double, pair < JPetHit, JPetHit >> & a, const pair < double, pair < JPetHit, JPetHit >> & b) {

  return abs(a.first) < abs(b.first);

}




double Calc3DAngleHit_EF(const JPetHit & Hit1,const JPetHit & Hit2 , TVector3 Center_EF) {

  double scalarProd = (Hit1.getPosX() - Center_EF.X()) * (Hit2.getPosX() - Center_EF.X()) + (Hit1.getPosY() - Center_EF.Y()) * (Hit2.getPosY() - Center_EF.Y()) + (Hit1.getPosZ() - Center_EF.Z()) * (Hit2.getPosZ() - Center_EF.Z());

  double magProd = sqrt((pow(Hit1.getPosX() - Center_EF.X(), 2) +
      pow(Hit1.getPosY() - Center_EF.Y(), 2) +
      pow(Hit1.getPosZ() - Center_EF.Z(), 2)) *
    (pow(Hit2.getPosX() - Center_EF.X(), 2) +
      pow(Hit2.getPosY() - Center_EF.Y(), 2) +
      pow(Hit2.getPosZ() - Center_EF.Z(), 2)));

  double Angle = acos(scalarProd / magProd) * 180 / 3.14159;

  return Angle;
}



double CalScatterTestHit_EF(const JPetHit & Hit1,
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

double CalScatterTestTrueMC_EF(const JPetMCHit & Hit1,
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


double CalExpecValueMC_EF(const JPetMCHit Hit1, const JPetMCHit Hit2, const JPetMCHit Hit3)
{

 
  	TVector3 k1 = 	Hit1.getMomentum();
	TVector3 k1p = 	Hit2.getMomentum();
	TVector3 k2 = 	Hit3.getMomentum();

       TVector3 epsilon = k1.Cross(k1p);

       double ExpecValue = ((epsilon.Dot(k2)) / ((epsilon.Mag()) * (k2.Mag())));

       return ExpecValue;

}



double CalExpecValueHit_EF(const JPetHit Hit1, const JPetHit Hit2, const JPetHit Hit3)
{

  TVector3 k1(Hit1.getPosX() - Center_EF.X() , Hit1.getPosY() - Center_EF.Y() , Hit1.getPosZ() - Center_EF.Z());
  TVector3 k1p(Hit2.getPosX() - Hit1.getPosX(), Hit2.getPosY() - Hit1.getPosY(), Hit2.getPosZ() - Hit1.getPosZ());
  TVector3 k2(Hit3.getPosX() - Center_EF.X() , Hit3.getPosY() - Center_EF.Y(), Hit3.getPosZ() - Center_EF.Z());

       TVector3 epsilon = k1.Cross(k1p);

       double ExpecValue = ((epsilon.Dot(k2)) / ((epsilon.Mag()) * (k2.Mag())));

       return ExpecValue;

}





Bool_t comparison_scatter_Hit_EF(const pair < double, pair < JPetHit, int >> & a, const pair < double, pair < JPetHit, int >> & b) {

  return abs(a.first) < abs(b.first);

}


Bool_t comparison_scatter_MC_EF(const pair < double, pair < JPetMCHit, int >> & a, const pair < double, pair < JPetMCHit, int >> & b) {

  return abs(a.first) < abs(b.first);

}

 vector<JPetEvent> eventVec;

vector<JPetEvent> EventFinder::buildEvents(const JPetTimeWindow& timeWindow)
{
 


  const unsigned int nHits = timeWindow.getNumberOfEvents();
  unsigned int count = 0;
  while(count<nHits){
    auto hit = dynamic_cast<const JPetHit&>(timeWindow.operator[](count));
    if(!fUseCorruptedHits && hit.getRecoFlag()==JPetHit::Corrupted){
      count++;
      continue;
    }
    // Creating new event with the first hit
    JPetEvent event;
    event.setEventType(JPetEventType::kUnknown);
    event.addHit(hit);
    if(hit.getRecoFlag() == JPetHit::Good) {
      event.setRecoFlag(JPetEvent::Good);
    } else if(hit.getRecoFlag() == JPetHit::Corrupted){
      event.setRecoFlag(JPetEvent::Corrupted);
    }
    // Checking, if following hits fulfill time window condition,
    // then moving interator 
    unsigned int nextCount = 1;
    while(count+nextCount < nHits){
      auto nextHit = dynamic_cast<const JPetHit&>(timeWindow.operator[](count+nextCount));
      if (fabs(nextHit.getTime() - hit.getTime()) < fEventTimeWindow) {
        if(nextHit.getRecoFlag() == JPetHit::Corrupted) {
          event.setRecoFlag(JPetEvent::Corrupted);
        }
        event.addHit(nextHit);
        nextCount++;
      } else { break; }
    }
    count+=nextCount;
    if(fSaveControlHistos) {
      getStatistics().getHisto1D("hits_per_event_all")->Fill(event.getHits().size());
      if(event.getRecoFlag()==JPetEvent::Good){
        getStatistics().getHisto1D("good_vs_bad_events")->Fill(1);
      } else if(event.getRecoFlag()==JPetEvent::Corrupted){
        getStatistics().getHisto1D("good_vs_bad_events")->Fill(2);
      } else {
        getStatistics().getHisto1D("good_vs_bad_events")->Fill(3);
      }
    }

 if(event.getHits().size() == 4){
     int Multiplicity = event.getHits().size();
     getStatistics().getHisto1D("Multiplicity")->Fill(Multiplicity);
     

  const JPetTimeWindowMC* time_window_mc = nullptr;
  if (time_window_mc = dynamic_cast<const JPetTimeWindowMC*>(fEvent)) {
    fIsMC = true;    
  	}
     	
	if(fIsMC){
	auto mc_hit_0 = time_window_mc->getMCHit<JPetMCHit>(event.getHits().at(0).getMCindex());
	auto mc_hit_1 = time_window_mc->getMCHit<JPetMCHit>(event.getHits().at(1).getMCindex());
	auto mc_hit_2 = time_window_mc->getMCHit<JPetMCHit>(event.getHits().at(2).getMCindex());
	auto mc_hit_3 = time_window_mc->getMCHit<JPetMCHit>(event.getHits().at(3).getMCindex());
	
	int Multiplicity_0 = mc_hit_0.getGenGammaMultiplicity();
	int Multiplicity_1 = mc_hit_1.getGenGammaMultiplicity();
	int Multiplicity_2 = mc_hit_2.getGenGammaMultiplicity();
	int Multiplicity_3 = mc_hit_3.getGenGammaMultiplicity();
	
	getStatistics().getHisto1D("True_MC_Multiplicity")->Fill(Multiplicity_0);		//Check the MC true multiplicity flag
	getStatistics().getHisto1D("True_MC_Multiplicity")->Fill(Multiplicity_1);		//Check the MC true multiplicity flag
	getStatistics().getHisto1D("True_MC_Multiplicity")->Fill(Multiplicity_2);		//Check the MC true multiplicity flag
	getStatistics().getHisto1D("True_MC_Multiplicity")->Fill(Multiplicity_3);

	std::vector<int> Multi;
	Multi.push_back(Multiplicity_0);		
	Multi.push_back(Multiplicity_1);
	Multi.push_back(Multiplicity_2);
	Multi.push_back(Multiplicity_3);

	//Sort and seletect only events with three '3' flag candidates and one '103' flag candidate//
 	
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


	getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->Fill(Multiplicity_0);		//Cut on the MC true multiplicity
	getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->Fill(Multiplicity_1);		//Cut on the MC true multiplicity
	getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->Fill(Multiplicity_2);		//Cut on the MC true multiplicity
	getStatistics().getHisto1D("True_MC_Multiplicity_Cut")->Fill(Multiplicity_3);		//Cut on the MC true multiplicity
 	
	
	std::vector <pair<double,pair<JPetHit,JPetMCHit>>> Hit_vector;

	Hit_vector.push_back({Multiplicity_0,{event.getHits().at(0), mc_hit_0}});
	Hit_vector.push_back({Multiplicity_1,{event.getHits().at(1), mc_hit_1}});
	Hit_vector.push_back({Multiplicity_2,{event.getHits().at(2), mc_hit_2}});
	Hit_vector.push_back({Multiplicity_3,{event.getHits().at(3), mc_hit_3}});

	std::sort(Hit_vector.begin(), Hit_vector.end(), comparison_gen_mul);

	
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
	double Angle_01_hit = Calc3DAngleHit_EF(Hit_vector.at(0).second.first, Hit_vector.at(1).second.first, Center_EF); 
	double Angle_12_hit = Calc3DAngleHit_EF(Hit_vector.at(1).second.first, Hit_vector.at(2).second.first, Center_EF); 
	double Angle_20_hit = Calc3DAngleHit_EF(Hit_vector.at(2).second.first, Hit_vector.at(0).second.first, Center_EF); 

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
	std::sort(Energy_vector_mc_hit.begin(), Energy_vector_mc_hit.end(), comparison_gen_mul);	



	//Before doing the scatter test, we need the time, energy and X(), Y(), resolutions//

	//X() Resolution//
	getStatistics().getHisto1D("X_resolution")->Fill(Hit_vector.at(0).second.first.getPosX() - Hit_vector.at(0).second.second.getPosX());		//Resolution of the Primary_1
	getStatistics().getHisto1D("X_resolution")->Fill(Hit_vector.at(1).second.first.getPosX() - Hit_vector.at(1).second.second.getPosX());		//Resolution of the Primary_2
	getStatistics().getHisto1D("X_resolution")->Fill(Hit_vector.at(2).second.first.getPosX() - Hit_vector.at(2).second.second.getPosX());		//Resolution of the Primary_3
	getStatistics().getHisto1D("X_resolution")->Fill(Hit_vector.at(3).second.first.getPosX() - Hit_vector.at(3).second.second.getPosX());		//Resolution of the Secondary


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
	getStatistics().getHisto1D("Time_resolution")->Fill(((Hit_vector.at(0).second.first.getTime())/1000. - (Hit_vector.at(0).second.second.getTime()))/1000.);		//Resolution of the Primary_1
	getStatistics().getHisto1D("Time_resolution")->Fill(((Hit_vector.at(1).second.first.getTime())/1000. - (Hit_vector.at(1).second.second.getTime()))/1000.);		//Resolution of the Primary_2
	getStatistics().getHisto1D("Time_resolution")->Fill(((Hit_vector.at(2).second.first.getTime())/1000. - (Hit_vector.at(2).second.second.getTime()))/1000.);		//Resolution of the Primary_3
	getStatistics().getHisto1D("Time_resolution")->Fill(((Hit_vector.at(3).second.first.getTime())/1000. - (Hit_vector.at(3).second.second.getTime()))/1000.);		//Resolution of the Secondary
	


	//Scatter Test - With Hit smeared//
	
	double ScatterTest_2_Hit = CalScatterTestHit_EF(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first);
	double ScatterTest_1_Hit = CalScatterTestHit_EF(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first);
	double ScatterTest_0_Hit = CalScatterTestHit_EF(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first);

	std::vector <pair<double,pair<JPetHit,int>>> Scatter_Vector_Hit;
	Scatter_Vector_Hit.push_back({ScatterTest_2_Hit,{Energy_vector_mc_hit.at(0).second.first, 2}});
	Scatter_Vector_Hit.push_back({ScatterTest_1_Hit,{Energy_vector_mc_hit.at(1).second.first, 1}});
	Scatter_Vector_Hit.push_back({ScatterTest_0_Hit,{Energy_vector_mc_hit.at(2).second.first, 0}});

	std::sort(Scatter_Vector_Hit.begin(), Scatter_Vector_Hit.end(), comparison_scatter_Hit_EF);	
	
	getStatistics().getHisto1D("Scatter_Hit_Lowest")->Fill(Scatter_Vector_Hit.at(0).first);
	getStatistics().getHisto2D("Scatter_Hit_21")->Fill(ScatterTest_2_Hit, ScatterTest_1_Hit);
	


	//Scatter Test - With True MC Hit//
	
	double ScatterTest_2_MC = CalScatterTestTrueMC_EF(Energy_vector_mc_hit.at(2).second.second, Hit_vector.at(3).second.second);
	double ScatterTest_1_MC = CalScatterTestTrueMC_EF(Energy_vector_mc_hit.at(1).second.second, Hit_vector.at(3).second.second);
	double ScatterTest_0_MC = CalScatterTestTrueMC_EF(Energy_vector_mc_hit.at(0).second.second, Hit_vector.at(3).second.second);

	std::vector <pair<double,pair<JPetMCHit,int>>> Scatter_Vector_MC;
	Scatter_Vector_MC.push_back({ScatterTest_2_MC,{Energy_vector_mc_hit.at(0).second.second, 0}});
	Scatter_Vector_MC.push_back({ScatterTest_1_MC,{Energy_vector_mc_hit.at(1).second.second, 1}});
	Scatter_Vector_MC.push_back({ScatterTest_0_MC,{Energy_vector_mc_hit.at(2).second.second, 2}});

	std::sort(Scatter_Vector_MC.begin(), Scatter_Vector_MC.end(), comparison_scatter_MC_EF);	

getStatistics().getHisto1D("Scatter_MC_Lowest")->Fill(Scatter_Vector_MC.at(0).first);
	getStatistics().getHisto2D("Scatter_MC_21")->Fill(ScatterTest_2_MC, ScatterTest_1_MC);


	//Cut on the Scatter Test//

	if(Scatter_Vector_MC.at(0).first <= 0.5 && Scatter_Vector_MC.at(0).first >= -0.5)
	{

	getStatistics().getHisto1D("Scatter_MC_Cut")->Fill(Scatter_Vector_MC.at(0).first);
	getStatistics().getHisto1D("Scatter_Hit_Cut")->Fill(Scatter_Vector_Hit.at(0).first);
	

			
	//Generate a random seed Uniformly//

	double rdm = r.Uniform(0., 1.);

	getStatistics().getHisto1D("random_generator")->Fill(rdm);
	
		

		

	if(Scatter_Vector_Hit.at(0).second.second == 2)
	{



	double E2K1_Hit = CalExpecValueHit_EF(Energy_vector_mc_hit.at(2).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(1).second.first);
	double E2K1_MC = CalExpecValueMC_EF(Energy_vector_mc_hit.at(2).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(1).second.second);
	getStatistics().getHisto1D("EV_Hit_EF")->Fill(E2K1_Hit);
	getStatistics().getHisto1D("EV_MC_EF")->Fill(E2K1_MC);



	//--------------Loop for inducing asymmetry-----------------//


	if(E2K1_MC < 0.)
	{

	//-----------Most important Vector Push Back--------------//

		eventVec.push_back(event);
	//--------------------------------------------------------//
	getStatistics().getHisto1D("EV_MC_EF_Induced")->Fill(E2K1_MC);
	getStatistics().getHisto1D("random_generator_negative")->Fill(rdm);

			}

	if(E2K1_MC > 0.)
	{


	if(rdm >= 0.2)
	{


	getStatistics().getHisto1D("random_generator_positive")->Fill(rdm);
		getStatistics().getHisto1D("EV_MC_EF_Induced")->Fill(E2K1_MC);
	//-----------Most important Vector Push Back--------------//

		eventVec.push_back(event);
	//--------------------------------------------------------//
	
				}

		
			}

	if(E2K1_MC == 0.)
	{

	
		if(rdm >= 0.5)
	{

	getStatistics().getHisto1D("random_generator_zero")->Fill(rdm);
	getStatistics().getHisto1D("EV_MC_EF_Induced")->Fill(E2K1_MC);
	//-----------Most important Vector Push Back--------------//

		eventVec.push_back(event);
	//--------------------------------------------------------//
	
				}

	




			}


				
				

	}

	else if(Scatter_Vector_Hit.at(0).second.second == 1)
	{




	double E1K0_Hit = CalExpecValueHit_EF(Energy_vector_mc_hit.at(1).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(0).second.first);
	double E1K0_MC = CalExpecValueMC_EF(Energy_vector_mc_hit.at(1).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(0).second.second);

	getStatistics().getHisto1D("EV_Hit_EF")->Fill(E1K0_Hit);
	getStatistics().getHisto1D("EV_MC_EF")->Fill(E1K0_MC);


	
	//--------------Loop for inducing asymmetry-----------------//


	if(E1K0_MC < 0.)
	{

	getStatistics().getHisto1D("random_generator_negative")->Fill(rdm);
	getStatistics().getHisto1D("EV_MC_EF_Induced")->Fill(E1K0_MC);
	//-----------Most important Vector Push Back--------------//

		eventVec.push_back(event);
	//--------------------------------------------------------//
	


			}

	if(E1K0_MC > 0.)
	{


	if(rdm >= 0.2)
	{


	getStatistics().getHisto1D("random_generator_positive")->Fill(rdm);
	getStatistics().getHisto1D("EV_MC_EF_Induced")->Fill(E1K0_MC);
	//-----------Most important Vector Push Back--------------//

		eventVec.push_back(event);
	//--------------------------------------------------------//
	
				}

		
			}

	if(E1K0_MC == 0.)
	{

	
		if(rdm >= 0.5)
	{

	getStatistics().getHisto1D("random_generator_zero")->Fill(rdm);
	getStatistics().getHisto1D("EV_MC_EF_Induced")->Fill(E1K0_MC);
	//-----------Most important Vector Push Back--------------//

		eventVec.push_back(event);
	//--------------------------------------------------------//
	
				}

	




			}

			
		


	}

	else if(Scatter_Vector_Hit.at(0).second.second == 0)
	{



	double E0K2_Hit = CalExpecValueHit_EF(Energy_vector_mc_hit.at(0).second.first, Hit_vector.at(3).second.first, Energy_vector_mc_hit.at(2).second.first);
	double E0K2_MC = CalExpecValueMC_EF(Energy_vector_mc_hit.at(0).second.second,Hit_vector.at(3).second.second,Energy_vector_mc_hit.at(2).second.second);

	
	getStatistics().getHisto1D("EV_Hit_EF")->Fill(E0K2_Hit);
	getStatistics().getHisto1D("EV_MC_EF")->Fill(E0K2_MC);




	
		//--------------Loop for inducing asymmetry-----------------//


	if(E0K2_MC < 0.)
	{

	getStatistics().getHisto1D("EV_MC_EF_Induced")->Fill(E0K2_MC);
	getStatistics().getHisto1D("random_generator_negative")->Fill(rdm);
	//-----------Most important Vector Push Back--------------//

		eventVec.push_back(event);
	//--------------------------------------------------------//
	


			}

	if(E0K2_MC > 0.)
	{


	if(rdm >= 0.2)
	{

	getStatistics().getHisto1D("EV_MC_EF_Induced")->Fill(E0K2_MC);
	getStatistics().getHisto1D("random_generator_positive")->Fill(rdm);
	//-----------Most important Vector Push Back--------------//

		eventVec.push_back(event);
	//--------------------------------------------------------//
	
				}

		
			}
	
	if(E0K2_MC == 0.)
	{

	
		if(rdm >= 0.5)
	{

	getStatistics().getHisto1D("random_generator_zero")->Fill(rdm);
	getStatistics().getHisto1D("EV_MC_EF_Induced")->Fill(E0K2_MC);
	//-----------Most important Vector Push Back--------------//

		eventVec.push_back(event);
	//--------------------------------------------------------//
	
				}

	




			}

		
		


	}

	//---------------------------------------//
	
			}//Scatter Test Cut


		}//Loop for counting gen multiplicity

	}//Loop for hits if MC

		
      if(fSaveControlHistos) {
        getStatistics().getHisto1D("hits_per_event_selected")->Fill
(event.getHits().size());


      }

    }

  }


  return eventVec;

}






void EventFinder::initialiseHistograms(){
  getStatistics().createHistogram(
    new TH1F("hits_per_event_all", "Number of Hits in an all Events", 20, 0.5, 20.5)
  );
  getStatistics().getHisto1D("hits_per_event_all")->GetXaxis()->SetTitle("Hits in Event");
  getStatistics().getHisto1D("hits_per_event_all")->GetYaxis()->SetTitle("Number of Hits");

  getStatistics().createHistogram(
    new TH1F("hits_per_event_selected", "Number of Hits in selected Events (min. multiplicity)",
    20, fMinMultiplicity-0.5, fMinMultiplicity+19.5)
  );
  getStatistics().getHisto1D("hits_per_event_selected")->GetXaxis()->SetTitle("Hits in Event");
  getStatistics().getHisto1D("hits_per_event_selected")->GetYaxis()->SetTitle("Number of Hits");

  getStatistics().createHistogram(
    new TH1F("good_vs_bad_events", "Number of good and corrupted Events created", 3, 0.5, 3.5)
  );
  getStatistics().getHisto1D("good_vs_bad_events")->GetXaxis()->SetBinLabel(1,"GOOD");
  getStatistics().getHisto1D("good_vs_bad_events")->GetXaxis()->SetBinLabel(2,"CORRUPTED");
  getStatistics().getHisto1D("good_vs_bad_events")->GetXaxis()->SetBinLabel(3,"UNKNOWN");
  getStatistics().getHisto1D("good_vs_bad_events")->GetYaxis()->SetTitle("Number of Events");

//Histograms for Juhi//

getStatistics().createHistogram(new TH1F("Multiplicity", "Multiplicity",
      11, -0.5, 10.5));
    getStatistics().getHisto1D("Multiplicity")->GetXaxis()->SetTitle("Multiplicity");
    getStatistics().getHisto1D("Multiplicity")->GetYaxis()->SetTitle("Counts");

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


    getStatistics().createHistogram(new TH2F("Angle_3D_True","Angle_3D_True",191, -0.5, 190., 191, -0.5, 190.));
    getStatistics().getHisto2D("Angle_3D_True")->GetXaxis()->SetTitle("#theta_{i} [deg]");
    getStatistics().getHisto2D("Angle_3D_True")->GetYaxis()->SetTitle("#theta_{j} [deg]");


    getStatistics().createHistogram(new TH2F("Angle_3D_Hit","Angle_3D_Hit",191, -0.5, 190., 191, -0.5, 190.));
    getStatistics().getHisto2D("Angle_3D_Hit")->GetXaxis()->SetTitle("#theta_{i} [deg]");
    getStatistics().getHisto2D("Angle_3D_Hit")->GetYaxis()->SetTitle("#theta_{j} [deg]");



    getStatistics().createHistogram(new TH1F("Energy_all_primary", "Energy_all_primary",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_all_primary")->GetXaxis()->SetTitle("Energy_{i} [keV]");
    getStatistics().getHisto1D("Energy_all_primary")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("Energy_secondary", "Energy_secondary",
      601, -0.5, 600.5));
    getStatistics().getHisto1D("Energy_secondary")->GetXaxis()->SetTitle("Energy_{i} [keV]");
    getStatistics().getHisto1D("Energy_secondary")->GetYaxis()->SetTitle("Counts");	



    getStatistics().createHistogram(new TH1F("Y_resolution", "Y_resolution",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("Y_resolution")->GetXaxis()->SetTitle("Y_Pos [cm]");
    getStatistics().getHisto1D("Y_resolution")->GetYaxis()->SetTitle("Counts");


    getStatistics().createHistogram(new TH1F("X_resolution", "X_resolution",
      1010, -50.5, 50.5));
    getStatistics().getHisto1D("X_resolution")->GetXaxis()->SetTitle("X_Pos [cm]");
    getStatistics().getHisto1D("X_resolution")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Energy_resolution", "Energy_resolution",
      2010, -100.5, 100.5));
    getStatistics().getHisto1D("Energy_resolution")->GetXaxis()->SetTitle("Energy [keV]");
    getStatistics().getHisto1D("Energy_resolution")->GetYaxis()->SetTitle("Counts");

    getStatistics().createHistogram(new TH1F("Time_resolution", "Time_resolution",
      10010, -500.5, 500.5));
    getStatistics().getHisto1D("Time_resolution")->GetXaxis()->SetTitle("Time [ns]");
    getStatistics().getHisto1D("Time_resolution")->GetYaxis()->SetTitle("Counts");

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
  
   getStatistics().createHistogram(new TH1F("EV_Hit_EF", "EV_Hit_EF",
      25, -1.0, 1.0));
    getStatistics().getHisto1D("EV_Hit_EF")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_Hit_EF")->GetYaxis()->SetTitle("Counts");

   getStatistics().createHistogram(new TH1F("EV_MC_EF", "EV_MC_EF",
      25, -1.0, 1.0));
    getStatistics().getHisto1D("EV_MC_EF")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_MC_EF")->GetYaxis()->SetTitle("Counts");
	
    getStatistics().createHistogram(new TH1F("EV_MC_EF_Induced", "EV_MC_EF_Induced",
      25, -1.0, 1.0));
    getStatistics().getHisto1D("EV_MC_EF_Induced")->GetXaxis()->SetTitle("Cosine #alpha_{Three Combinations}");
    getStatistics().getHisto1D("EV_MC_EF_Induced")->GetYaxis()->SetTitle("Counts");




//Random Generator Checks//



    getStatistics().createHistogram(new TH1F("random_generator", "random_generator",
      500, -2.0, 2.0));
    getStatistics().getHisto1D("random_generator")->GetXaxis()->SetTitle("random_generator");
    getStatistics().getHisto1D("random_generator")->GetYaxis()->SetTitle("Counts");



    getStatistics().createHistogram(new TH1F("random_generator_positive", "random_generator_positive",
      500, -2.0, 2.0));
    getStatistics().getHisto1D("random_generator_positive")->GetXaxis()->SetTitle("random_generator_positive");
    getStatistics().getHisto1D("random_generator_positive")->GetYaxis()->SetTitle("Counts");



    getStatistics().createHistogram(new TH1F("random_generator_zero", "random_generator_zero",
      500, -2.0, 2.0));
    getStatistics().getHisto1D("random_generator_zero")->GetXaxis()->SetTitle("random_generator_zero");
    getStatistics().getHisto1D("random_generator_zero")->GetYaxis()->SetTitle("Counts");

	
    getStatistics().createHistogram(new TH1F("random_generator_negative", "random_generator_negative",
      500, -2.0, 2.0));
    getStatistics().getHisto1D("random_generator_negative")->GetXaxis()->SetTitle("random_generator_negative");
    getStatistics().getHisto1D("random_generator_negative")->GetYaxis()->SetTitle("Counts");




}
