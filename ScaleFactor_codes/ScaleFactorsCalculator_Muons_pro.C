#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <math.h>

using std::cout;
using std::endl;

vector<double> calculate(TString file_name, TString outfile_name, int eta_bin, int pt_bin);


void ScaleFactorsCalculator_Muons_pro()
{

  TString outfile_post = "ScaleFactors_postVFP_Muon_Trigger.txt";
  TString outfile_pre = "ScaleFactors_preVFP_Muon_Trigger.txt";
  ofstream outfile("EffectiveScaleFactors_2016UL_Muon_Trigger.txt");
  
  TString file_pre = "/home/bowerbird/Yash/Work/ZZ_Analysis/Scale_Factors/2016UL/Medium_pT_Muon/Trigger/Efficiencies_muon_generalTracks_Z_Run2016_UL_preVFP_SingleMuonTriggers.root";
  TString file_post = "/home/bowerbird/Yash/Work/ZZ_Analysis/Scale_Factors/2016UL/Medium_pT_Muon/Trigger/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers_postVFP.root";

  //Id and Iso
  int eta_bin=4;
  int pt_bin=7;

  vector<double> SF_pre = calculate(file_pre, outfile_pre, eta_bin, pt_bin);
  vector<double> SF_post = calculate(file_post, outfile_post, eta_bin, pt_bin);

  for(int k =0; k<eta_bin*pt_bin; k++){
    double effScaleFactor = ( SF_pre.at(k)*20 + SF_post.at(k)*16 )/(36);
    outfile<<setprecision(18);
    outfile<<k+1<<"\t"<<effScaleFactor;    
    if(k != eta_bin*pt_bin-1)
      outfile<<endl;
  }
  

  
 
}

vector<double> calculate(TString file_name, TString outfile_name, int eta_bin, int pt_bin){

  TString file1 = file_name;
  TString file2 = file1;
  ofstream file(outfile_name);

  //ID
  //TString plotname1 = "NUM_MediumID_DEN_TrackerMuons_abseta_pt_efficiencyData";
  //TString plotname2 = "NUM_MediumID_DEN_TrackerMuons_abseta_pt_efficiencyMC";

  //Isolation
  //TString plotname1 = "NUM_TightRelIso_DEN_MediumID_abseta_pt_efficiencyData";
  //TString plotname2 = "NUM_TightRelIso_DEN_MediumID_abseta_pt_efficiencyMC";

  //Trigger
  TString plotname1 = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt_efficiencyData";
  TString plotname2 = "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_abseta_pt_efficiencyMC";

  //Now let us open the files
  TFile *file_1 = new TFile(file1);
  TFile *file_2 = new TFile(file2);
  
  //Now open the respective histograms from the files
  TH2F *h1 = (TH2F*)file_1->Get(plotname1);
  TH2F *h2 = (TH2F*)file_2->Get(plotname2);
  
  int binno = 0;
  vector<double> SF_array;
  SF_array.clear();
  
  for(int i=1; i<eta_bin+1 ;i++){          //Reads the bins of x axis:Eta
    for(int j=1; j<pt_bin+1 ;j++){         //Reads the bins of y axis:Pt
      double Num = h1->GetBinContent(i,j);//Efficiency of the Data corresponding to bin(i,j)
      double Den = h2->GetBinContent(i,j);//Efficiency of the MC corresponding to bin(i,j)
      double SF = Num/Den;                //Scale Factor= Data eff/ MC eff
      binno++;                            // Calculates the total bin numbers
      
      cout<<setprecision(18);
      cout<<"Bin: "<<binno<<endl;
      cout<<"Numerator is : "<<Num<<endl;
      cout<<"Denominator is: "<<Den<<endl;
      cout<<"ScaleFactor is: "<<SF<<endl;
      
      file<<setprecision(18);
      file<<binno<<"\t";
      file<<Num<<"\t";
      file<<Den<<"\t";
      file<<SF;     
      
      if(binno!=eta_bin*pt_bin)
	file<<endl;

      SF_array.push_back(SF);
      
    }
  }
  
  file.close();
  return(SF_array);

}

