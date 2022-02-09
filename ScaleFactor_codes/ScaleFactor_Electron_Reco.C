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


void ScaleFactor_Electron_Reco()
{

  /*TString outfile_pre  = "SF_preVFP_Electron_Reco_above20.txt";
  TString outfile_post = "SF_postVFP_Electron_Reco_above20.txt";
  ofstream outfile("EffectiveSF_2016UL_Electron_Reco_Above20.txt");*/
  TString outfile_pre  = "SF_preVFP_Electron_Reco_below20.txt";
  TString outfile_post = "SF_postVFP_Electron_Reco_below20.txt";
  ofstream outfile("EffectiveSF_2016UL_Electron_Reco_below20.txt");
  
  //TString file_pre = "/home/bowerbird/Yash/Work/ZZ_Analysis/Scale_Factors/2016UL/Electrons/Reco/EGM2D_ptAbove20_UL2016preVFP.root";
  TString file_pre = "/home/bowerbird/Yash/Work/ZZ_Analysis/Scale_Factors/2016UL/Electrons/Reco/EGM2D_ptBelow20_UL2016preVFP.root";
  //TString file_post = "/home/bowerbird/Yash/Work/ZZ_Analysis/Scale_Factors/2016UL/Electrons/Reco/EGM2D_ptAbove20_UL2016postVFP.root";
  TString file_post = "/home/bowerbird/Yash/Work/ZZ_Analysis/Scale_Factors/2016UL/Electrons/Reco/EGM2D_ptBelow20_UL2016postVFP.root";

  //Above 20
  //int eta_bin=12;
  //int pt_bin=4;

  //Below 20
  int eta_bin=10;
  int pt_bin=1;

  vector<double> SF_pre = calculate(file_pre, outfile_pre, eta_bin, pt_bin);
  vector<double> SF_post = calculate(file_post, outfile_post, eta_bin, pt_bin);

  for(int k =0; k<eta_bin*pt_bin; k++){
    double ScaleFactor = ( SF_pre.at(k)*19.893068320 + SF_post.at(k)*16.392692052 )/(36.285760373);
    outfile<<setprecision(18);
    outfile<<k+1<<"\t"<<ScaleFactor;    
    if(k != eta_bin*pt_bin-1)
      outfile<<endl;
  }
  

  
 
}

vector<double> calculate(TString file_name, TString outfile_name, int eta_bin, int pt_bin){

  TString file1 = file_name;
  ofstream file(outfile_name);

  //Trigger
  TString plotname1 = "EGamma_SF2D";
  
  //Now let us open the files
  TFile *file_1 = new TFile(file1);
    
  //Now open the respective histograms from the files
  TH2F *h1 = (TH2F*)file_1->Get(plotname1);
  
  int binno = 0;
  vector<double> sf_array;
  sf_array.clear();

  for(int j=1; j<pt_bin+1 ;j++){         //Reads the bins of y axis:Pt
    for(int i=1; i<eta_bin+1 ;i++){          //Reads the bins of x axis:Eta
      double Number = h1->GetBinContent(i,j);//Efficiency of the Data corresponding to bin(i,j)
      binno++;                            // Calculates the total bin numbers
      
      cout<<setprecision(18);
      cout<<"Bin: "<<binno<<endl;
      cout<<"Scale factor is : "<<Number<<endl;
      
      file<<setprecision(18);
      file<<binno<<"\t";
      file<<Number;
      
      if(binno!=eta_bin*pt_bin)
	file<<endl;

      sf_array.push_back(Number);
      
    }
  }
  
  file.close();
  return(sf_array);

}
