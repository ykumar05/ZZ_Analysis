#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
//#include <format>
//#include <iomanip>
using namespace std;
//using std::cout;
//using std::endl;

/*****************************************************************************************************************************
 *                                            User Defined Functions                                                         *
 *****************************************************************************************************************************/

//Define the function decorate for histograms
void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) {

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);
  h->SetLineColor(linecolor); 
  h->SetLineWidth(linewidth);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  if(tofill==1) h->SetFillColor(markercolor);

  h->SetMarkerSize(1);
  h->SetTitle(title);
}

//Define the function decorate for legends
void decorate(TLegend *g, float textSize, TString legendheader)
{
  g->SetTextSize(textSize);
  g->SetFillStyle(1);
  g->SetFillColor(10);
  g->SetBorderSize(1);
  g->SetLineColor(10);
  //Usually legends should not have headers
  //g->SetHeader(legendheader);
}


// Here are a couple of other utility functions

// For a given histogram hst, return the number of entries between bin_lo and bin_hi
float get_nevents(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents += hst->GetBinContent(i);

  return nevents;
}
// Partner function for above, returning the error for the above nevents
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi)
{
  int bin_width = hst->GetBinWidth(1);
  int ibin_begin = 1;
  float nevents = 0, nevents_err = 0;
  while ( hst->GetBinCenter(ibin_begin) < bin_lo )
    ibin_begin++;
  int ibin_end = ibin_begin;
  while ( hst->GetBinCenter(ibin_end) < bin_hi )
    ibin_end++;
  for(int i=ibin_begin; i<ibin_end; i++)
    nevents_err += pow(hst->GetBinError(i),2);
  nevents_err = sqrt(nevents_err);

  return nevents;
}

float luminosity(float num, float xsec){
  //put sigma in pb^{-1}, which is smaller than fb^{-1}
  float numerator = num;
  float denominator = xsec;
  float lumi = numerator/denominator;
  return lumi;
}

float round_int(float var) //rounding up the integral on the legend
{
    char str[40];
    sprintf(str, "%.2f", var);
    sscanf(str, "%f", &var);
    return var;
}

/***************************************************************************************************************
 *                                            Main Function                                                    *
 ***************************************************************************************************************/

void ovlStack()
{

  //First declare the file names:SET1
  /*TString zzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set1/hst_ZZ.root";
  TString ttzfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set1/hst_TTZ.root";
  TString ttwfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set1/hst_TTW.root";
  TString wzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set1/hst_WZ.root";
  TString datafile = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set1/hst_data.root";*/

  //First declare the file names:SET2
  /*TString zzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set2/hst_ZZ.root";
  TString ttzfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set2/hst_TTZ.root";
  TString ttwfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set2/hst_TTW.root";
  TString wzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set2/hst_WZ.root";
  TString datafile = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set2/hst_data.root";*/

  //First declare the file names:SET3
  /*  TString zzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set3/hst_ZZ.root";
  TString ttzfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set3/hst_TTZ.root";
  TString ttwfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set3/hst_TTW.root";
  TString wzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set3/hst_WZ.root";
  TString datafile = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set3/hst_data.root";*/

  //First declare the file names:SET4
  /*TString zzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set4/hst_ZZ.root";
  TString ttzfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set4/hst_TTZ.root";
  TString ttwfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set4/hst_TTW.root";
  TString wzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set4/hst_WZ.root";
  TString datafile = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set4/hst_data.root";*/

  //First declare the file names:SET5
  TString zzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set5/hst_ZZ.root";
  TString ttzfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set5/hst_TTZ.root";
  TString ttwfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set5/hst_TTW.root";
  TString wzfile   = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set5/hst_WZ.root";
  TString datafile = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set5/hst_data.root";

  //Luminosity(num, xsex) = Number of events / (production cross-section *1000)
  float zzlumi  = luminosity(17733356,13.74)/1000;
  float ttzlumi = luminosity(2962856,0.2439)/1000;
  float ttwlumi = luminosity(1800823,0.2161)/1000;
  float wzlumi  = luminosity(6890010,5.052)/1000;  
  float dtlumi  = 30.497663487;

  cout<< "Luminosity of the MC samples are :" << endl;
  cout<< "Lumi_ZZ= "  << zzlumi << endl;
  cout<< "Lumi_TTZ= " << ttzlumi<< endl;
  cout<< "Lumi_TTW= " << ttwlumi<< endl;
  cout<< "Lumi_WZ= "  << wzlumi << endl;
  cout<< "Lumi_data= "<< dtlumi << endl;

  cout<< "scale Factors :" << endl;
  cout<< "Scale_ZZ= "  << (dtlumi/zzlumi) << endl;
  cout<< "Scale_TTZ= " << (dtlumi/ttzlumi)<< endl;
  cout<< "Scale_TTW= " << (dtlumi/ttwlumi)<< endl;
  cout<< "Scale_WZ= "  << (dtlumi/wzlumi) << endl;
  cout<< "Scale_data= "<< (dtlumi/dtlumi) << endl;
  
  
  //Now let us open the files
  TFile *filezz   = new TFile(zzfile);
  TFile *filettz  = new TFile(ttzfile);
  TFile *filettw  = new TFile(ttwfile);
  TFile *filewz   = new TFile(wzfile);
  TFile *filedata = new TFile(datafile);

  //Histograms you want to plot

  /*****************************************
   *              4 Muons                  *
   *****************************************/
  TString chplots[9] = {"rA_mu_pT1","rA_mu_pT2","rA_mu_pT3","rA_mu_pT4","rA_mu_NbJets","rA_mu_MET_pT","rA_mu_PuppiMET_pT","rA_mu_DiMuonInvariantMass_ChosenPairing","rA_mu_QuadMuon_Mass"};
  //TString chplots[9] = {"rB_mu_pT1","rB_mu_pT2","rB_mu_pT3","rB_mu_pT4","rB_mu_NbJets","rB_mu_MET_pT","rB_mu_PuppiMET_pT","rB_mu_DiMuonInvariantMass_ChosenPairing","rB_mu_QuadMuon_Mass"};
  //TString chplots[9] = {"rC_mu_pT1","rC_mu_pT2","rC_mu_pT3","rC_mu_pT4","rC_mu_NbJets","rC_mu_MET_pT","rC_mu_PuppiMET_pT","rC_mu_DiMuonInvariantMass_ChosenPairing","rC_mu_QuadMuon_Mass"};

  /*****************************************
   *       2 Electrons and 2 Muons         *
   *****************************************/
  //TString chplots[9] = {"rA_emu_pT1_Muon","rA_emu_pT2_Muon","rA_emu_pT1_Electron","rA_emu_pT2_Electron","rA_emu_NbJets","rA_emu_MET_pT","rA_emu_PuppiMET_pT","rA_emu_DileptonInvariantMass","rA_emu_Quadlepton_Mass"};

  /*****************************************
   *              Combined 4L              *
   *****************************************/
  //TString chplots[4] = {"NbJets_4L","MET_4L","PuppiMET_4L","DiMuonInvMass_4L"};
  
  TString ytitle = "Entries"; // Or "Events"
  TString xtitle[9] = {"pT1","pT2","pT3","pT4","NbJets","MET_pT","PuppiMET_pT","DiMuonInvariantMass_ChosenPairing","QuadMuon_Mass"};
  //TString xtitle[9] = {"pT1_Muon","pT2_Muon","pT1_Electron","pT2_Electron","NbJets","MET_pT","PuppiMET_pT","DileptonInvariantMass","QuadMuonMass_Mass"};
  //TString xtitle[4] = {"NbJets_4L","MET_4L","PuppiMET_4L","DileptonInvariantMass_4L"};

  TString plotname1;
  TH1F *h0, *h1, *h2, *h3, *h4, *g0;
  TH1F *h0new, *h1new, *h2new, *h3new, *h4new, *g0new;
  int Nbins=4;
  
  for(int i=0; i<9; i++){
    plotname1 = chplots[i];
    
    //Now open the respective histograms from the files
    h0 = (TH1F*)filezz   ->Get(plotname1);
    h1 = (TH1F*)filettz  ->Get(plotname1);
    h2 = (TH1F*)filettw  ->Get(plotname1);
    h3 = (TH1F*)filewz   ->Get(plotname1);
    g0 = (TH1F*)filedata ->Get(plotname1);
    
    decorate(h0,xtitle[i],ytitle,"",kGreen+3,2,kGreen+3,20,0);//ZZ
    decorate(h1,xtitle[i],ytitle,"",kRed+1,2,kRed+1,20,0);//TTZ
    decorate(h2,xtitle[i],ytitle,"",kViolet,2,kViolet,20,0);//TTW
    decorate(h3,xtitle[i],ytitle,"",kBlue,2,kBlue,20,0);//WZ
    decorate(g0,xtitle[i],ytitle,"",kBlack,2,kBlack,20,0);//Data

    //Scaling:
    h0->Scale(dtlumi/zzlumi);//ZZ
    h1->Scale(dtlumi/ttzlumi);//TTZ
    h2->Scale(dtlumi/ttwlumi);//TTW
    h3->Scale(dtlumi/wzlumi);//WZ

    //Rebinning:
    
    h0new = (TH1F*)h0->Rebin(Nbins);//ZZ
    h1new = (TH1F*)h1->Rebin(Nbins);//TTZ
    h2new = (TH1F*)h2->Rebin(Nbins);//TTW
    h3new = (TH1F*)h3->Rebin(Nbins);//WZ
    g0new = (TH1F*)g0->Rebin(Nbins);//Data   


    //set last bin as overflow
    h0new->SetBinContent(h0new->GetNbinsX(),h0new->GetBinContent(h0new->GetNbinsX()+1)+h0new->GetBinContent(h0new->GetNbinsX()));//ZZ
    h1new->SetBinContent(h1new->GetNbinsX(),h1new->GetBinContent(h1new->GetNbinsX()+1)+h1new->GetBinContent(h1new->GetNbinsX()));//TTZ
    h2new->SetBinContent(h2new->GetNbinsX(),h2new->GetBinContent(h2new->GetNbinsX()+1)+h2new->GetBinContent(h2new->GetNbinsX()));//TTW
    h3new->SetBinContent(h3new->GetNbinsX(),h3new->GetBinContent(h3new->GetNbinsX()+1)+h3new->GetBinContent(h3new->GetNbinsX()));//WZ
    g0new->SetBinContent(g0new->GetNbinsX(),g0new->GetBinContent(g0new->GetNbinsX()+1)+g0new->GetBinContent(g0new->GetNbinsX()));//Data
    
    THStack *stack = new THStack("Stacked","");
    h0new->SetFillColor(kAzure+7);
    h1new->SetFillColor(kGreen+3);
    h2new->SetFillColor(kRed+1);
    h3new->SetFillColor(kViolet);
    stack->Add(h0new);
    stack->Add(h1new);
    stack->Add(h2new);
    stack->Add(h3new);
    
    //Now let us declare a canvas 
    TCanvas *c1 = new TCanvas(plotname1,plotname1,800,600);
    //Draw stack:
    stack->Draw("hist");
    //Draw overlay:
    g0new->Draw("ep same");
    //c1->SetLogy();
    c1->SetTickx(1);
    c1->SetTicky(1);

    //x-axis:
    stack->GetXaxis()->SetTitle(xtitle[i]);
    
    //y-axis:
    stack->GetYaxis()->SetTitle(ytitle);
    stack->GetYaxis()->SetTitleOffset(1.75);
    stack->GetYaxis()->SetTitleSize(0.03);
    stack->GetYaxis()->SetLabelSize(0.03);
    
    //stack->SetMaximum(100000);
    //stack->SetMinimum(0.01);

    TLegend *lg1 = new TLegend(0.6,0.65,0.85,0.87,NULL,"NDC");
    decorate(lg1,0.03,"");
    TString name;

    //Putting numbers on the legend:
    
    name = "Data [" + to_string((int)g0new->Integral()) + "]"; 
    lg1->AddEntry(g0new,name,"f");
    
    name = "ZZ [" + to_string((int)h0new->Integral()) + "]";
    lg1->AddEntry(h0new,name,"f");

    name = "TTZ [" + to_string((int)h1new->Integral()) + "]";
    lg1->AddEntry(h1new,name,"f");

    name = "TTW [" + to_string((int)h2new->Integral()) + "]";
    lg1->AddEntry(h2new,name,"f");

    name = "WZ [" + to_string((int)h3new->Integral()) + "]";
    lg1->AddEntry(h3new,name,"f");

    lg1->Draw();
  }
}
