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

using std::cout;
using std::endl;

void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill);
void decorate(TLegend *g, float textSize, TString legendheader);
float get_nevents(TH1F *hst, float bin_lo, float bin_hi);
float get_nevents_err(TH1F *hst, float bin_lo, float bin_hi);

void overlay4()
{
  //Set the input files :
  TString zzfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set1/hst_ZZ.root";
  TString ttzfile = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set1/hst_TTZ.root";
  TString ttwfile = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set1/hst_TTW.root";
  TString wzfile  = "/home/bowerbird/Yash/Work/ZZ_Analysis/hst_files/Set1/hst_WZ.root";

  //Set the output plot name :
  TString output = "overlay_PuppiMET"; //(for example: metpT.png)
  //DON'T add space and the file extension.

  //Select the plots :
  TString plotname1 = "DiMuonInvMass_Chosen";
  TString plotname2 = plotname1;
  TString plotname3 = plotname1;
  TString plotname4 = plotname1;
  TString xtitle = "DiMuonInvMass_Chosen";
  TString ytitle = "Events";

  //Set Legend Entries :
  TString legendheader = "DiMuonInvMass_Chosen";
  TString entry1 = "ZZ";
  TString entry2 = "TTZ";
  TString entry3 = "TTW";
  TString entry4 = "WZ";

  //Choose which plot is from which file :
  TFile *file1 = new TFile(zzfile); //for plot1
  TFile *file2 = new TFile(ttzfile); //for plot2
  TFile *file3 = new TFile(ttwfile); //for plot3
  TFile *file4 = new TFile(wzfile); //for plot3

  TH1F *h1 = (TH1F*)file1->Get(plotname1);
  TH1F *h2 = (TH1F*)file2->Get(plotname2);
  TH1F *h3 = (TH1F*)file3->Get(plotname3);
  TH1F *h4 = (TH1F*)file4->Get(plotname4);

  //Set your colour scheme, design etc here :
  decorate(h1,xtitle,ytitle,"",kBlack,2,kBlack,20,0);
  decorate(h2,xtitle,ytitle,"",kRed,2,kRed,20,0);
  decorate(h3,xtitle,ytitle,"",kViolet, 2, kViolet, 20, 0);
  decorate(h4,xtitle,ytitle,"",kBlue, 2, kBlue, 20, 0);

  //Normalisation : (uncomment this section, if you want to compare only the shapes)
  /*h1->Scale(1.0/h1->Integral());
  h2->Scale(1.0/h2->Integral());
  h3->Scale(1.0/h3->Integral());
  h4->Scale(1.0/h4->Integral());*/

  TCanvas *c1 = new TCanvas(output,legendheader,800,600);

  TLegend *lg1 = new TLegend(0.55,0.50,0.85,0.75,NULL,"NDC");
  decorate(lg1,0.05,legendheader); // Decorated the legend using function below.
  lg1->AddEntry(h1, entry1, "lf"); // Added the two entries for the two histograms
  lg1->AddEntry(h2, entry2, "lf");           // we shall be drawing.
  lg1->AddEntry(h3, entry3, "lf");
  lg1->AddEntry(h4, entry4, "lf");
  
  //We set the stat box for the first one to be invisible first
  h1->SetStats(0);

  // Now to draw the histograms and the legend
  // NOTE : Draw the bigger one first (do hit and trial).
  h1->Draw("hist");
  h2->Draw("hist same");            
  h3->Draw("hist same");
  h4->Draw("hist same");

  lg1->Draw();

  /* To draw with stat boxes for each histogram
     -- Dont use SetStats(0)
     -- Then draw them, first one with option Draw(), next ones with option Draw("sames")
     -- The s at the end is for stats box
  */
}

void decorate(TH1*h,const char* xtitle, const char* ytitle, const char* title,
	      int linecolor, int linewidth, int markercolor, int markerstyle, int tofill) {

  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);

  h->SetLineColor(linecolor); 
  h->SetLineWidth(linewidth);
  h->SetMarkerColor(markercolor);
  h->SetMarkerStyle(markerstyle);
  if(tofill==1) h->SetFillColor(markercolor);

  h->SetMarkerSize(1.0);
  h->SetTitle(title);
  //h->SetMinimum(0);
  //h->SetMaximum(5500);
}

void decorate(TLegend *g, float textSize, TString legendheader)
{
  g->SetTextSize(textSize);
  g->SetFillStyle(1);
  g->SetFillColor(10);
  g->SetBorderSize(1);
  g->SetLineColor(10);
  //Usually legends should not have headers
  //but if you want one, uncomment the next line.
  g->SetHeader(legendheader);
}
