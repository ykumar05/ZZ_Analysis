#define ZZAna_cxx
// The class definition in ZZAna.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("ZZAna.C")
// root> T->Process("ZZAna.C","some options")
// root> T->Process("ZZAna.C+")
//


#include "ZZAna.h"
#include <TH2.h>
#include <TStyle.h>

void ZZAna::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
}

void ZZAna::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   nEvtTotal = 0;
   nEvtRan = 0;
   nCount4L=0;
   nCount2L=0;
   nCount3L=0;
   nCount2L2L=0;
   nCount_LHEweight_pos=0;
   nCount_LHEweight_neg=0;
   
   //Create the histogram file
   _HstFile = new TFile(_HstFileName,"recreate");
   //Call the function to book the histograms we declared in Hists.
   BookHistograms();

}

void ZZAna::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

  _HstFile->Write();
  _HstFile->Close();
  int n_eff = nCount_LHEweight_pos-nCount_LHEweight_neg;
  //Output to screen.
  cout<<"Total events ran = "<<nEvtRan<<endl;
  cout<<"Total events = "<<nEvtTotal<<endl;
  cout<<"Total events with atleast 4 muons(trigger applied) are: "<<nCount4L<<endl;
  cout<<"Total events with two muons and two electrons(triggers applied) are: "<<nCount2L2L<<endl;
  cout<<"Total 4L events are "<<nCount4L+nCount2L2L<<endl;
  cout<<"Number of events for which there is positive weight: "<<nCount_LHEweight_pos<<endl;
  cout<<"Number of events for which there is negative weight: "<<nCount_LHEweight_neg<<endl;
  cout<<"Effective number of events are: "<<n_eff<<endl;
  //Open the text output file
  ofstream fout(_SumFileName);
  //Put text output in the summary file.
  fout<<"Total events ran = "<<nEvtRan<<endl;
  fout<<"Total events  = "<<nEvtTotal<<endl;
}

void ZZAna::Terminate(){
  
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}

Bool_t ZZAna::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

  fReader.SetLocalEntry(entry);
  if(_data == 0)
    fReader_MC  .SetLocalEntry(entry);
  if(_data == 1)
    fReader_Data.SetLocalEntry(entry);
  float evtweight = 1;

   //Output processing information to screen based on verbosity level.
   if(_verbosity==0 && nEvtTotal%100000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;      
   else if(_verbosity>0 && nEvtTotal%100000==0)cout<<"Processed "<<nEvtTotal<<" event..."<<endl;
   
   nEvtTotal++;
   nEvtRan++;

   if(_data == 0){
     h.weights[0]->Fill(*LHEWeight_originalXWGTUP);
     if(*LHEWeight_originalXWGTUP>=0)
       nCount_LHEweight_pos++;
     if(*LHEWeight_originalXWGTUP<0)
       nCount_LHEweight_neg++;
   }
     

   GenEle.clear();
   GenMu.clear();
   llep.clear();
   
   /****************************************
    *               Gen Muons               *           
    *****************************************/
   if(_data==0){
     for(UInt_t i=0; i<*nGenPart; i++){
       Lepton temp;
       temp.v.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
       temp.ind=i; temp.id=GenPart_pdgId[i];
       bool PassCuts_Mu = temp.v.Pt()>5 && fabs(temp.v.Eta())<2.4;
       bool truemuon = false;
       if(fabs(GenPart_pdgId[i])==13 && GenPart_status[i]==1)
	 truemuon = true;
       
       if(PassCuts_Mu && truemuon){
	 GenMu.push_back(temp);
       }
     }
   }
   h.ngen[0]->Fill((int)GenMu.size());
   for(int i=0; i<(int)GenMu.size(); i++){
     h.genPart[0]->Fill(GenMu.at(i).v.Pt());
     h.genPart[1]->Fill(GenMu.at(i).v.Eta());
     h.genPart[2]->Fill(GenMu.at(i).v.Phi());
   }
   
   
   /****************************************
    *             Gen Electrons             *           
    *****************************************/
   if(_data==0){
     for(UInt_t i=0; i<*nGenPart; i++){
       Lepton temp;
       temp.v.SetPtEtaPhiM(GenPart_pt[i], GenPart_eta[i], GenPart_phi[i], GenPart_mass[i]);
       temp.ind=i; temp.id=GenPart_pdgId[i];
       bool PassCuts_Ele = temp.v.Pt()>5 && temp.v.Eta()<2.5;
       bool true_electron = false;
       if(fabs(GenPart_pdgId[i])==11 && GenPart_status[i]==1)
	 true_electron = true;
       
       if(true_electron && PassCuts_Ele){
	 GenEle.push_back(temp);
       }	
     }
   }
   h.ngen[1]->Fill((int)GenEle.size());
   for(int i=0; i<(int)GenEle.size(); i++){
     h.genPart[3]->Fill(GenEle.at(i).v.Pt());
     h.genPart[4]->Fill(GenEle.at(i).v.Eta());
     h.genPart[5]->Fill(GenEle.at(i).v.Phi());
   }
   
   
   /**********************************************************************
    *                          RECO PARTICLES                            *
    **********************************************************************/
   
   /****************************************
    *                Muons                 *           
    *****************************************/
   GoodMu.clear();
   GoodMu_pos.clear();
   GoodMu_neg.clear();
   for(unsigned int i=0; i<(*nMuon); i++){
     Lepton temp; temp.v.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i],Muon_phi[i],0.105);
     temp.id = -13*Muon_charge[i]; temp.ind = i;  temp.charge = Muon_charge[i];

     bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Muon_mediumId[i];
     passCuts = passCuts && Muon_pfRelIso04_all[i]<0.15;
     passCuts = passCuts && fabs(Muon_dxy[i])<0.05 && fabs(Muon_dz[i])<0.1;

     //dR matching with GenMuons
     if(_data==0){
       float dR_Mu_GenMu = 999.0;
       for(int i=0; i<(int)GenMu.size(); i++){
	 float dR_value= temp.v.DeltaR(GenMu.at(i).v);
	 if(dR_value<dR_Mu_GenMu){
	 dR_Mu_GenMu=dR_value;
	 }
       }
       h.dR[0]->Fill(dR_Mu_GenMu);
       passCuts = passCuts && dR_Mu_GenMu<0.1;
     }
     if(passCuts){
       GoodMu.push_back(temp);
       llep.push_back(temp);
     }
     if(passCuts && temp.charge == 1){
       GoodMu_pos.push_back(temp);//Creating mu+ array
     }
     if(passCuts && temp.charge == -1){
       GoodMu_neg.push_back(temp);//Creating mu- array
     }
   }
   SortPt(0);
   SortPt(1);
   SortPt(2);

   //Muon Array
   h.muons[0]->Fill((int)GoodMu.size());
   for(int i=0; i<(int)GoodMu.size(); i++){
     h.muons[1]->Fill(GoodMu.at(i).v.Pt());
     h.muons[2]->Fill(GoodMu.at(i).v.Eta());
     h.muons[3]->Fill(GoodMu.at(i).v.Phi());
   }
   
   
   /****************************************
    *               Electrons              *           
    *****************************************/
   GoodEle.clear();
   for(unsigned int i=0; i<(*nElectron); i++){
     Lepton temp; temp.v.SetPtEtaPhiM(Electron_pt[i],Electron_eta[i],Electron_phi[i],0.000511); 
     temp.id = -11*Electron_charge[i]; temp.ind = i; temp.charge = Electron_charge[i];
     
     bool isprompt = false;
     if(fabs(temp.v.Eta())<=1.479)
       if(fabs(Electron_dxy[i])<0.05 && fabs(Electron_dz[i])<0.1)
	 isprompt = true;      
     if(fabs(temp.v.Eta())>1.479)
       if(fabs(Electron_dxy[i])<0.1 && fabs(Electron_dz[i])<0.2)
	 isprompt = true;

     bool passCuts = temp.v.Pt()>10 && fabs(temp.v.Eta())<2.4 && Electron_cutBased[i]>2;
     passCuts = passCuts && isprompt;

     //dR matching with GenElectrons
     if(_data==0){
       float dR_Ele_GenEle = 999.0;
       for(int j=0; j<(int)GenEle.size(); j++){
	 float dR_value= temp.v.DeltaR(GenEle.at(j).v);
	 if(dR_value<dR_Ele_GenEle){
	   dR_Ele_GenEle=dR_value;
	 }
       }
       h.dR[1]->Fill(dR_Ele_GenEle);
       passCuts = passCuts && dR_Ele_GenEle<0.1;
     }
     if(passCuts){
       GoodEle.push_back(temp);
       llep.push_back(temp);
     }
   }
   SortPt(3);
   SortPt(4);
   
   h.electrons[0]->Fill((int)GoodEle.size());
   for(int i=0; i<(int)GoodEle.size(); i++){
     h.electrons[1]->Fill(GoodEle.at(i).v.Pt());
     h.electrons[2]->Fill(GoodEle.at(i).v.Eta());
     h.electrons[3]->Fill(GoodEle.at(i).v.Phi());
   }

   //Lightleptons
   h.llep[0]->Fill((int)llep.size());

   
   /****************************************
    *                  Taus                *           
    ****************************************/
   GoodTau.clear();
    for(unsigned int i=0; i< (*nTau); i++){
      if(Tau_idDecayModeOldDMs[i]&&(Tau_decayMode[i]<3||Tau_decayMode[i]>9)){
	//Tau energy scale correction
	float tlv_corr = 1.;
	if(_year==2016){
	  if(Tau_decayMode[i]==0) tlv_corr = 0.994;
	  if(Tau_decayMode[i]==1) tlv_corr = 0.995;
	  if(Tau_decayMode[i]>9)  tlv_corr = 1;
	}     
	if(_year==2017){
	  if(Tau_decayMode[i]==0) tlv_corr = 1.007;
	  if(Tau_decayMode[i]==1) tlv_corr = 0.998;
	  if(Tau_decayMode[i]==10) tlv_corr = 1.001;
	  if(Tau_decayMode[i]==11) tlv_corr = 0.999;
	}
	if(_year==2018){
	  if(Tau_decayMode[i]==0) tlv_corr = 0.987;
	  if(Tau_decayMode[i]==1) tlv_corr = 0.995;
	  if(Tau_decayMode[i]==10) tlv_corr = 0.998;
	  if(Tau_decayMode[i]==11) tlv_corr = 1;
	}
	Lepton temp; temp.v.SetPtEtaPhiM(Tau_pt[i],Tau_eta[i],Tau_phi[i],1.77);
	temp.v *= tlv_corr; //energy correction
	temp.id = -15*Tau_charge[i]; temp.ind = i; temp.charge = Tau_charge[i];
	
	bool passCuts = temp.v.Pt()>20 && fabs(temp.v.Eta()<2.3);
	passCuts = passCuts && (Tau_decayMode[i]<3 || Tau_decayMode[i]>9);
	passCuts = passCuts && fabs(Tau_dz[i])<0.2;
	passCuts = passCuts && Tau_idDeepTau2017v2p1VSe[i] >= 15 && Tau_idDeepTau2017v2p1VSmu[i] >= 3;//medium and above WP
	passCuts = passCuts && Tau_idDeepTau2017v2p1VSjet[i] >= 31; //medium WP

	//dR cleaning for tau : Muons
	bool ismuonclose = false;
	for(int j=0; j<(int)GoodMu.size(); j++){
	  float dR_mutau =temp.v.DeltaR(GoodMu.at(j).v);
	  if(dR_mutau<0.4)
	    ismuonclose = true;
	}
	//dR cleaning for tau : Electrons
	bool iselectronclose = false;
	for(int j=0; j<(int)GoodEle.size(); j++){
	  float dR_etau = temp.v.DeltaR(GoodEle.at(j).v);
	  if(dR_etau<0.4)
	    iselectronclose = true;
	}
	
	if(passCuts && !ismuonclose && !iselectronclose){
	  GoodTau.push_back(temp);
	}
      }
    }
    


   /****************************************
    *                  Jets                *           
    ****************************************/
   GoodJet.clear();
   bJet.clear();
   for(unsigned int i=0; i<(*nJet); i++ ){
     Lepton temp; temp.v.SetPtEtaPhiM(Jet_pt[i],Jet_eta[i],Jet_phi[i],Jet_mass[i]); 
     temp.ind = i; 
     bool passCuts = temp.v.Pt()>30 && fabs(temp.v.Eta())<2.4;

     //dR cleaning of the i-th jet : Muons
     float dR_mu_jet=999.0;
     for(int j=0; j<(int)GoodMu.size(); j++){
       float dR_value = temp.v.DeltaR(GoodMu.at(j).v);
       if(dR_value<dR_mu_jet){
	 dR_mu_jet=dR_value;
       }
     }
     h.dR[2]->Fill(dR_mu_jet);
     
     //dR cleaning of the i-th jet : Electrons
     float dR_ele_jet=999.0;
     for(int j=0; j<(int)GoodEle.size(); j++){
       float dR_value =temp.v.DeltaR(GoodEle.at(j).v);
       if(dR_value<dR_ele_jet){
	 dR_ele_jet=dR_value;
       }
     }
     h.dR[3]->Fill(dR_ele_jet);

     //dR cleaning of the i-th jet : Taus
     float dR_tau_jet=999.0;
     for(int j=0; j<(int)GoodTau.size(); j++){
       float dR_value =temp.v.DeltaR(GoodTau.at(j).v);
       if(dR_value<dR_tau_jet){
	 dR_tau_jet=dR_value;
       }
     }
     h.dR[4]->Fill(dR_tau_jet);
     
     passCuts = passCuts && dR_mu_jet>0.4 && dR_ele_jet>0.4 && dR_tau_jet>0.4;

     if(passCuts){
       GoodJet.push_back(temp);
       //b Tagged Jets
       h.jets[4]->Fill(Jet_btagDeepB[i]);
     }

     if(passCuts && Jet_btagDeepB[i]>=0.5847){
       bJet.push_back(temp);
     }
     
   }
   
   h.jets[0]->Fill((int)GoodJet.size());
   for(int i=0; i<(int)GoodJet.size(); i++){
     h.jets[1]->Fill(GoodJet.at(i).v.Pt());
     h.jets[2]->Fill(GoodJet.at(i).v.Eta());
     h.jets[3]->Fill(GoodJet.at(i).v.Phi());
   }

   h.jets[5]->Fill((int)bJet.size());
   for(int i=0; i<(int)bJet.size(); i++){
     h.jets[6]->Fill(bJet.at(i).v.Pt());
     h.jets[7]->Fill(bJet.at(i).v.Eta());
     h.jets[8]->Fill(bJet.at(i).v.Phi());
   }

   /****************************************
    *                3 Muons               *           
    ****************************************/
   
   
   if((int)GoodMu.size()>2){

     if(_data == 0){
       float generator_weight  = fabs(*Generator_weight)/ *Generator_weight;
       
       float evtweight_Reco_mu    = ( (getScaleFactors_Muons_Reco(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt())) );

       
       float evtweight_ID_mu      = ( (getScaleFactors_Muons_ID(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt())) );
       
       float evtweight_Iso_mu     = ( (getScaleFactors_Muons_Iso(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt())) );
       
       float evtweight_Trigger_mu = 1-( (1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt())) );
       
       
       evtweight = generator_weight*evtweight_Reco_mu*evtweight_ID_mu*evtweight_Iso_mu*evtweight_Trigger_mu;
       

       h.analysis[17]->Fill((GoodMu.at(1).v+GoodMu.at(2).v).M(),evtweight);
     }
   }
   
   
   /****************************************
    *             4L Analysis              *           
    ****************************************/

   MassArray_chosen.clear();   //MassArray_chosen tells us the masses of the chosen pairs(ZZ->llll)
   MassArray_unchosen.clear();
   MassArray_2e2mu.clear();   //MassArray_2e2mu tells us the masses of the Z pairs(ZZ->llLL)
   MassArray_4L.clear();

   //Clearing the arrays of pairings
   pairing1.clear();//[0123]
   pairing2.clear();//[0213]
   pairing3.clear();//[0312]

   //Initializing different masses for 4 muons state
   float mass_01_12=0.0;//01:GoodMu.at(0),GoodMu.at(1); 12:Pairing1 and Pairing2
   float mass_23_12=0.0;
   float mass_02_12=0.0;
   float mass_13_12=0.0;
   float mass_01_13=0.0;
   float mass_23_13=0.0;
   float mass_03_13=0.0;
   float mass_12_13=0.0;
   float mass_02_23=0.0;
   float mass_13_23=0.0;
   float mass_03_23=0.0;
   float mass_12_23=0.0;

   //Initializing different masses for 2 muons and 2 electrons state
   float mass_2e2mu_2mu=0.0;
   float mass_2e2mu_2e =0.0;
   							  
   //Defining booleans for 4 muons state
   bool mu4       = false;
   bool chosen1   = false;
   bool chosen2   = false;
   bool chosen3   = false;
   bool chosen4   = false;
   bool chosen5   = false;
   bool chosen6   = false;
   bool pairing12 = false;
   bool pairing13 = false;
   bool pairing23 = false;

   //Defining booleans for 2e2mu state
   bool e2mu2         = false;
   bool chosen1_2e2mu = false;
   bool chosen2_2e2mu = false;
   
   /****************************************
    *                4 Muons               *           
    ****************************************/

   if((int)GoodMu.size()>3 && GoodMu.at(0).v.Pt()>26){
     mu4 = true;
     nCount4L++;
     if(GoodMu.at(0).id*GoodMu.at(1).id == -169){
       pairing1.push_back(GoodMu.at(0));
       pairing1.push_back(GoodMu.at(1));
       pairing1.push_back(GoodMu.at(2));
       pairing1.push_back(GoodMu.at(3));
     }
     if(GoodMu.at(0).id*GoodMu.at(2).id == -169){
       pairing2.push_back(GoodMu.at(0));
       pairing2.push_back(GoodMu.at(2));
       pairing2.push_back(GoodMu.at(1));
       pairing2.push_back(GoodMu.at(3));
     }
     if(GoodMu.at(0).id*GoodMu.at(3).id == -169){
       pairing3.push_back(GoodMu.at(0));
       pairing3.push_back(GoodMu.at(3));
       pairing3.push_back(GoodMu.at(1));
       pairing3.push_back(GoodMu.at(2));
     }

     if(pairing1.size()>3 && pairing2.size()>3){
       pairing12= true;
       mass_01_12 = (GoodMu.at(0).v+GoodMu.at(1).v).M();
       mass_23_12 = (GoodMu.at(2).v+GoodMu.at(3).v).M();
       mass_02_12 = (GoodMu.at(0).v+GoodMu.at(2).v).M();
       mass_13_12 = (GoodMu.at(1).v+GoodMu.at(3).v).M();
       float massdiff_pairing1 = fabs(mass_01_12-91.2)+fabs(mass_23_12-91.2);
       float massdiff_pairing2 = fabs(mass_02_12-91.2)+fabs(mass_13_12-91.2);
       if(massdiff_pairing1<massdiff_pairing2){
	 chosen1 = true;
	 h.analysis[0]->Fill(mass_01_12);
	 h.analysis[1]->Fill(mass_23_12);
       }
       else{
	 chosen2 = true;
	 h.analysis[2]->Fill(mass_02_12);
	 h.analysis[3]->Fill(mass_13_12);
       }
     }
     
     if(pairing1.size()>3 && pairing3.size()>3){
       pairing13 = true;
       mass_01_13 = (GoodMu.at(0).v+GoodMu.at(1).v).M();
       mass_23_13 = (GoodMu.at(2).v+GoodMu.at(3).v).M();
       mass_03_13 = (GoodMu.at(0).v+GoodMu.at(3).v).M();
       mass_12_13 = (GoodMu.at(1).v+GoodMu.at(2).v).M();
       float massdiff_pairing1 = fabs(mass_01_13-91.2)+fabs(mass_23_13-91.2);
       float massdiff_pairing3 = fabs(mass_03_13-91.2)+fabs(mass_12_13-91.2);
       if(massdiff_pairing1<massdiff_pairing3){
	 chosen3 = true;
	 h.analysis[4]->Fill(mass_01_13);
	 h.analysis[5]->Fill(mass_23_13);
       }
       else{
	 chosen4 = true;
	 h.analysis[6]->Fill(mass_03_13);
	 h.analysis[7]->Fill(mass_12_13);
       }
     }
     if(pairing2.size()>3 && pairing3.size()>3){
       pairing23 = true;
       mass_02_23 = (GoodMu.at(0).v+GoodMu.at(2).v).M();
       mass_13_23 = (GoodMu.at(1).v+GoodMu.at(3).v).M();
       mass_03_23 = (GoodMu.at(0).v+GoodMu.at(3).v).M();
       mass_12_23 = (GoodMu.at(1).v+GoodMu.at(2).v).M();
       float massdiff_pairing2 = fabs(mass_02_23-91.2)+fabs(mass_13_23-91.2);
       float massdiff_pairing3 = fabs(mass_03_23-91.2)+fabs(mass_12_23-91.2);
       if(massdiff_pairing2<massdiff_pairing3){
	 chosen5 = true;
	 h.analysis[8]->Fill(mass_02_23);
	 h.analysis[9]->Fill(mass_13_23);
       }
       else{
	 chosen6 = true;
	 h.analysis[10]->Fill(mass_03_23);
	 h.analysis[11]->Fill(mass_12_23);
       }
     }
     
     if(pairing12 == true){
       if(chosen1 == true){
	 MassArray_chosen.push_back(mass_01_12);//DiMuon mass of 0th and 1st muon if the event satisfies pairing 1 and pairing 2
	 MassArray_chosen.push_back(mass_23_12);
	 MassArray_4L.push_back(mass_01_12);
	 MassArray_4L.push_back(mass_23_12);
       }
       if(chosen2 == true){
	 MassArray_chosen.push_back(mass_02_12);
	 MassArray_chosen.push_back(mass_13_12);
	 MassArray_4L.push_back(mass_02_12);
	 MassArray_4L.push_back(mass_13_12);
       }
       if(chosen1 == false){
	 MassArray_unchosen.push_back(mass_01_12);
	 MassArray_unchosen.push_back(mass_23_12);
       }
       if(chosen2 == false){
	 MassArray_unchosen.push_back(mass_02_12);
	 MassArray_unchosen.push_back(mass_13_12);
       }
     }
     if(pairing13 == true){
       if(chosen3 == true){
	 MassArray_chosen.push_back(mass_01_13);
	 MassArray_chosen.push_back(mass_23_13);
	 MassArray_4L.push_back(mass_01_13);
	 MassArray_4L.push_back(mass_23_13);
       }
       if(chosen4 == true){
	 MassArray_chosen.push_back(mass_03_13);
	 MassArray_chosen.push_back(mass_12_13);
	 MassArray_4L.push_back(mass_03_13);
	 MassArray_4L.push_back(mass_12_13);
       }
       if(chosen3 == false){
	 MassArray_unchosen.push_back(mass_01_13);
	 MassArray_unchosen.push_back(mass_23_13);
       }
       if(chosen4 == false){
	 MassArray_unchosen.push_back(mass_03_13);
	 MassArray_unchosen.push_back(mass_12_13);
       }
     }
     if(pairing23 == true){
       if(chosen5 == true){
	 MassArray_chosen.push_back(mass_02_23);
	 MassArray_chosen.push_back(mass_13_23);
	 MassArray_4L.push_back(mass_02_23);
	 MassArray_4L.push_back(mass_13_23);
       }
       if(chosen6 == true){
	 MassArray_chosen.push_back(mass_03_23);
	 MassArray_chosen.push_back(mass_12_23);
	 MassArray_4L.push_back(mass_03_23);
	 MassArray_4L.push_back(mass_12_23);
       }
       if(chosen5 == false){
	 MassArray_unchosen.push_back(mass_02_23);
	 MassArray_unchosen.push_back(mass_13_23);
       }
       if(chosen6 == false){
	 MassArray_unchosen.push_back(mass_03_23);
	 MassArray_unchosen.push_back(mass_12_23);
       }
     }
   }
   
   if(mu4){

     if(_data == 0){
       float generator_weight  = fabs(*Generator_weight)/ *Generator_weight;
       
       float evtweight_Reco_mu    = ( (getScaleFactors_Muons_Reco(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       float evtweight_ID_mu      = ( (getScaleFactors_Muons_ID(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       float evtweight_Iso_mu     = ( (getScaleFactors_Muons_Iso(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       float evtweight_Trigger_mu = 1-( (1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       
       evtweight = generator_weight*evtweight_Reco_mu*evtweight_ID_mu*evtweight_Iso_mu*evtweight_Trigger_mu;
       
     }
     
     //cout<<"The evtweight is: "<<evtweight<<endl;
     
     h.ZZ_mu_h[0][0]->Fill(GoodMu.at(0).v.Pt(),evtweight);
     h.ZZ_mu_h[0][1]->Fill(GoodMu.at(1).v.Pt(),evtweight);
     h.ZZ_mu_h[0][2]->Fill(GoodMu.at(2).v.Pt(),evtweight);
     h.ZZ_mu_h[0][3]->Fill(GoodMu.at(3).v.Pt(),evtweight);
     h.ZZ_mu_h[0][4]->Fill((int)bJet.size(),evtweight);
     h.ZZ_mu_h[0][5]->Fill(*MET_pt,evtweight);
     h.ZZ_mu_h[0][6]->Fill(*PuppiMET_pt,evtweight);
     for(int i=0; i<(int)MassArray_chosen.size(); i++){
       if(MassArray_chosen.at(i)>15){
	 h.ZZ_mu_h[0][7]->Fill(MassArray_chosen.at(i),evtweight);
       }
     }
     /*for(int i=0; i<(int)MassArray_unchosen.size(); i++){
       h.analysis[13]->Fill(MassArray_unchosen.at(i)); 
       }*/
     h.ZZ_mu_h[0][8]->Fill((GoodMu.at(0).v+GoodMu.at(1).v+GoodMu.at(2).v+GoodMu.at(3).v).M(),evtweight);
   }
   
   if(mu4 && (int)bJet.size()==0){

     if(_data == 0){
       float generator_weight  = fabs(*Generator_weight)/ *Generator_weight;
       
       float evtweight_Reco_mu    = ( (getScaleFactors_Muons_Reco(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       float evtweight_ID_mu      = ( (getScaleFactors_Muons_ID(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       float evtweight_Iso_mu     = ( (getScaleFactors_Muons_Iso(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       float evtweight_Trigger_mu = 1-( (1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       
       evtweight = generator_weight*evtweight_Reco_mu*evtweight_ID_mu*evtweight_Iso_mu*evtweight_Trigger_mu;
       
     }
       
     h.ZZ_mu_h[1][0]->Fill(GoodMu.at(0).v.Pt(),evtweight);
     h.ZZ_mu_h[1][1]->Fill(GoodMu.at(1).v.Pt(),evtweight);
     h.ZZ_mu_h[1][2]->Fill(GoodMu.at(2).v.Pt(),evtweight);
     h.ZZ_mu_h[1][3]->Fill(GoodMu.at(3).v.Pt(),evtweight);
     h.ZZ_mu_h[1][4]->Fill((int)bJet.size(),evtweight);
     h.ZZ_mu_h[1][5]->Fill(*MET_pt,evtweight);
     h.ZZ_mu_h[1][6]->Fill(*PuppiMET_pt,evtweight);
     for(int i=0; i<(int)MassArray_chosen.size(); i++){
       if(MassArray_chosen.at(i)>15){
	 h.ZZ_mu_h[1][7]->Fill(MassArray_chosen.at(i),evtweight);
       }
     }
     h.ZZ_mu_h[1][8]->Fill((GoodMu.at(0).v+GoodMu.at(1).v+GoodMu.at(2).v+GoodMu.at(3).v).M(),evtweight);
   }

   if(mu4 && (int)bJet.size()==0 && *MET_pt<50){
    
    if(_data == 0){
       float generator_weight  = fabs(*Generator_weight)/ *Generator_weight;
       
       float evtweight_Reco_mu    = ( (getScaleFactors_Muons_Reco(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
				      *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       float evtweight_ID_mu      = ( (getScaleFactors_Muons_ID(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       float evtweight_Iso_mu     = ( (getScaleFactors_Muons_Iso(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       float evtweight_Trigger_mu = 1-( (1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(2).v.Eta()),GoodMu.at(2).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(3).v.Eta()),GoodMu.at(3).v.Pt())) );
       
       
       evtweight = generator_weight*evtweight_Reco_mu*evtweight_ID_mu*evtweight_Iso_mu*evtweight_Trigger_mu;
       
     }

     h.ZZ_mu_h[2][0]->Fill(GoodMu.at(0).v.Pt(),evtweight);
     h.ZZ_mu_h[2][1]->Fill(GoodMu.at(1).v.Pt(),evtweight);
     h.ZZ_mu_h[2][2]->Fill(GoodMu.at(2).v.Pt(),evtweight);
     h.ZZ_mu_h[2][3]->Fill(GoodMu.at(3).v.Pt(),evtweight);
     h.ZZ_mu_h[2][4]->Fill((int)bJet.size(),evtweight);
     h.ZZ_mu_h[2][5]->Fill(*MET_pt,evtweight);
     h.ZZ_mu_h[2][6]->Fill(*PuppiMET_pt,evtweight);
     for(int i=0; i<(int)MassArray_chosen.size(); i++){
       if(MassArray_chosen.at(i)>15){
	 h.ZZ_mu_h[2][7]->Fill(MassArray_chosen.at(i),evtweight);
       }
     }
     h.ZZ_mu_h[2][8]->Fill((GoodMu.at(0).v+GoodMu.at(1).v+GoodMu.at(2).v+GoodMu.at(3).v).M(),evtweight);
   }

   /****************************************
    *             2 Ele & 2 Muons          *           
    ****************************************/

   if((int)GoodMu.size()>1 && (int)GoodEle.size()>1 && GoodMu.at(0).v.Pt()>26 && GoodEle.at(0).v.Pt()>30){
     if(GoodMu.at(0).id*GoodMu.at(1).id == -169){
       chosen1_2e2mu = true;
       mass_2e2mu_2mu = (GoodMu.at(0).v+GoodMu.at(1).v).M();       
     }
     if(GoodEle.at(0).id*GoodEle.at(1).id == -121){
       chosen2_2e2mu = true;
       mass_2e2mu_2e = (GoodEle.at(0).v+GoodEle.at(1).v).M();
     }
     if(chosen1_2e2mu && chosen2_2e2mu){
       e2mu2 = true;
       nCount2L2L++;
       MassArray_2e2mu.push_back(mass_2e2mu_2mu);
       MassArray_2e2mu.push_back(mass_2e2mu_2e);
       MassArray_4L.push_back(mass_2e2mu_2mu);
       MassArray_4L.push_back(mass_2e2mu_2e);
     }
   }
   
   if(e2mu2){
     if(_data == 0){
       float generator_weight  = fabs(*Generator_weight)/ *Generator_weight;
       
       float evtweight_Reco_emu    = ( (getScaleFactors_Muons_Reco(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				       *(getScaleFactors_Muons_Reco(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt()))
				       *(getScaleFactors_Electrons_Reco(GoodEle.at(0).v.Eta(),GoodEle.at(0).v.Pt()))
				       *(getScaleFactors_Electrons_Reco(GoodEle.at(1).v.Eta(),GoodEle.at(1).v.Pt())) );
       
       float evtweight_ID_mu      = ( (getScaleFactors_Muons_ID(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_ID(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt())) );
       
       
       float evtweight_Iso_mu     = ( (getScaleFactors_Muons_Iso(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
				      *(getScaleFactors_Muons_Iso(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt())) );
       
       float evtweight_IDIso_e    = ( (getScaleFactors_Electrons_IDIso(GoodEle.at(0).v.Eta(),GoodEle.at(0).v.Pt()))
				      *(getScaleFactors_Electrons_IDIso(GoodEle.at(1).v.Eta(),GoodEle.at(1).v.Pt())) );
       
       
       float evtweight_Trigger_mu = 1-( (1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(0).v.Eta()),GoodMu.at(0).v.Pt()))
					*(1-getEfficiency_Muons_Trigger(fabs(GoodMu.at(1).v.Eta()),GoodMu.at(1).v.Pt())) );
       
       evtweight = generator_weight*evtweight_Reco_emu*evtweight_ID_mu*evtweight_IDIso_e*evtweight_Iso_mu*evtweight_Trigger_mu;
     }

     
     h.ZZ_emu_h[0][0]->Fill(GoodMu.at(0).v.Pt(),evtweight);
     h.ZZ_emu_h[0][1]->Fill(GoodMu.at(1).v.Pt(),evtweight);
     h.ZZ_emu_h[0][2]->Fill(GoodEle.at(0).v.Pt(),evtweight);
     h.ZZ_emu_h[0][3]->Fill(GoodEle.at(1).v.Pt(),evtweight);
     h.ZZ_emu_h[0][4]->Fill((int)bJet.size(),evtweight);
     h.ZZ_emu_h[0][5]->Fill(*MET_pt,evtweight);
     h.ZZ_emu_h[0][6]->Fill(*PuppiMET_pt,evtweight);
     for(int i=0; i<(int)MassArray_2e2mu.size(); i++){
       h.ZZ_emu_h[0][7]->Fill(MassArray_2e2mu.at(i),evtweight);
     }
     h.ZZ_emu_h[0][8]->Fill((GoodMu.at(0).v+GoodMu.at(1).v+GoodEle.at(0).v+GoodEle.at(1).v).M(),evtweight);
   }

   
   /****************************************
    *     Combined Final State Analysis    *        
    ****************************************/
   h.analysis[12]->Fill((int)MassArray_4L.size());
   for(int k=0; k<(int)MassArray_4L.size(); k++){
     h.analysis[13]->Fill(MassArray_4L.at(k));
   }
   if(mu4||e2mu2){
     h.analysis[14]->Fill((int)bJet.size());
     h.analysis[15]->Fill(*MET_pt);
     h.analysis[16]->Fill(*PuppiMET_pt);
   }
   
   return kTRUE;
}

//user Defined Functions



   
void ZZAna::BookHistograms(){
  //Gen Muons
  h.ngen[0]        = new TH1F("ngenMu","ngenMu",10,0,10);
  h.genPart[0]     = new TH1F("GenMuPt","GenMuon_pT",200,0,200);
  h.genPart[1]     = new TH1F("GenMuEta","GenMuon_Eta",50,-4,4);
  h.genPart[2]     = new TH1F("GenMuPhi","GenMuon_Phi",50,-4,4);
  
  //Gen Electrons
  h.ngen[1]        = new TH1F("ngenEle", "ngenEle", 10, 0, 10);
  h.genPart[3]     = new TH1F("GenElePt", "GenElectron_pT", 200, 0, 200);
  h.genPart[4]     = new TH1F("GenEleEta", "GenElectron_Eta",50, -4, 4);
  h.genPart[5]     = new TH1F("GenElePhi", "GenElectron_Phi",50,-4,4);
  
  //LightLeptons
  h.llep[0]        = new TH1F("Nllep","N_llep",10,0,10);

  //Muons
  h.muons[0]       = new TH1F("Nmuons","NMuons",10,0,10);
  h.muons[1]       = new TH1F("MuonPt","Muon_pT", 200, 0, 200);
  h.muons[2]       = new TH1F("MuonEta","Muon_Eta",50,-4,4);
  h.muons[3]       = new TH1F("MuonPhi", "Muon_Phi",50,-4,4);
  //h.muons[4]       = new TH1F("DiMuon_InvariantMass", "DiMuon_InvariantMass",200,0,200);
  
  //Electrons
  h.electrons[0]   = new TH1F("Nelectrons","NElectrons",10,0,10);
  h.electrons[1]   = new TH1F("ElectronPt","Electron_pT", 200, 0, 200);
  h.electrons[2]   = new TH1F("ElectronEta","Electron_Eta",50,-4,4);
  h.electrons[3]   = new TH1F("ElectronPhi", "Electron_Phi",50,-4,4);
  //h.electrons[4]   = new TH1F("DiElectron_InvariantMass", "DiElectron_InvariantMass",200,0,200);

  //Jets
  h.jets[0]        = new TH1F("Njets","NJets",10,0,10);
  h.jets[1]        = new TH1F("JetPt","Jet_pT", 200, 0, 200);
  h.jets[2]        = new TH1F("JetEta","Jet_Eta",50,-4,4);
  h.jets[3]        = new TH1F("JetPhi","Jet_Phi",50,-4,4);
  h.jets[4]        = new TH1F("Jet_btagDeepB","Jet_btagDeepB",400,-2,2);

  //bJets
  h.jets[5]        = new TH1F("Nbjets","NbJets",10,0,10);
  h.jets[6]        = new TH1F("bJetPt","bJet_pT", 200, 0, 200);
  h.jets[7]        = new TH1F("bJetEta","bJet_Eta",50,-4,4);
  h.jets[8]        = new TH1F("bJetPhi","bJet_Phi",50,-4,4);
    
  //dR Plots
  h.dR[0]          = new TH1F("dR_Mu_GenMu","dR_Mu_GenMu",20,0,1);
  h.dR[1]          = new TH1F("dR_Ele_GenEle","dR_Ele_GenEle",20,0,1);
  h.dR[2]          = new TH1F("dR_mu_jet","dR_mu_jet",100,0,1);
  h.dR[3]          = new TH1F("dR_ele_jet","dR_ele_jet",100,0,1);
  h.dR[4]          = new TH1F("dR_tau_jet","dR_tau_jet",100,0,1);

  //Analysis Plots
  h.analysis[0]    = new TH1F("DiMuonInvMass01_pairing12","DiMuonInvMass01_pairing12",200,0,200);
  h.analysis[1]    = new TH1F("DiMuonInvMass23_pairing12","DiMuonInvMass23_pairing12",200,0,200);
  h.analysis[2]    = new TH1F("DiMuonInvMass02_pairing12","DiMuonInvMass02_pairing12",200,0,200);
  h.analysis[3]    = new TH1F("DiMuonInvMass13_pairing12","DiMuonInvMass13_pairing12",200,0,200);
  h.analysis[4]    = new TH1F("DiMuonInvMass01_pairing13","DiMuonInvMass01_pairing13",200,0,200);
  h.analysis[5]    = new TH1F("DiMuonInvMass23_pairing13","DiMuonInvMass23_pairing13",200,0,200);
  h.analysis[6]    = new TH1F("DiMuonInvMass03_pairing13","DiMuonInvMass12_pairing13",200,0,200);
  h.analysis[7]    = new TH1F("DiMuonInvMass12_pairing13","DiMuonInvMass12_pairing13",200,0,200);
  h.analysis[8]    = new TH1F("DiMuonInvMass02_pairing23","DiMuonInvMass02_pairing23",200,0,200);
  h.analysis[9]    = new TH1F("DiMuonInvMass13_pairing23","DiMuonInvMass02_pairing23",200,0,200);
  h.analysis[10]   = new TH1F("DiMuonInvMass03_pairing23","DiMuonInvMass02_pairing03",200,0,200);
  h.analysis[11]   = new TH1F("DiMuonInvMass12_pairing23","DiMuonInvMass02_pairing12",200,0,200);
  //Combined Final State Analysis Plots
  h.analysis[12]   = new TH1F("MassArray_4L_size","MassArray_4L_size",10,0,10);
  h.analysis[13]   = new TH1F("DiMuonInvMass_4L","DimuonInvMass_4L",200,0,200);
  h.analysis[14]   = new TH1F("NbJets_4L","NbJets_4L",10,0,10);
  h.analysis[15]   = new TH1F("MET_4L","MET_4L",200,0,200);
  h.analysis[16]   = new TH1F("PuppiMET_4L","PuppiMET_4L",200,0,200);
  //3L Final state
  h.analysis[17]   = new TH1F("DiMuonInvMass_3L","DimuonInvMass_3L",300,0,15);

  
  //Four muon final state
  TString crname_mu[4] = {"rA_mu_","rB_mu_","rC_mu_","rD_mu_"};
  TString plotname_mu[9] = {"pT1","pT2","pT3","pT4","NbJets","MET_pT","PuppiMET_pT","DiMuonInvariantMass_ChosenPairing","QuadMuon_Mass"};
  int nbins_mu[9] = {200,200,200,200,5,100,100,200,600};
  float blo_mu[9] = {0,0,0,0,0,0,0,0,0};
  float bhi_mu[9] = {200,200,200,200,5,100,100,200,600};
  for(int icr_mu=0; icr_mu<3; icr_mu++){
    for(int iplot_mu=0; iplot_mu<9; iplot_mu++){
      TString name_mu = crname_mu[icr_mu] + plotname_mu[iplot_mu];
      h.ZZ_mu_h[icr_mu][iplot_mu] = new TH1F(name_mu,name_mu,nbins_mu[iplot_mu],blo_mu[iplot_mu],bhi_mu[iplot_mu]);
      h.ZZ_mu_h[icr_mu][iplot_mu]->Sumw2();
    }
  }

  // Two electron and two muon final state
  TString crname_emu[4] = {"rA_emu_","rB_emu_","rC_emu_","rD_emu_"};
  TString plotname_emu[9] = {"pT1_Muon","pT2_Muon","pT1_Electron","pT2_Electron","NbJets","MET_pT","PuppiMET_pT","DileptonInvariantMass","Quadlepton_Mass"};
  int nbins_emu[9] = {200,200,200,200,5,100,100,200,600};
  float blo_emu[9] = {0,0,0,0,0,0,0,0,0};
  float bhi_emu[9] = {200,200,200,200,5,100,100,200,600};
  for(int icr_emu=0; icr_emu<1; icr_emu++){
    for(int iplot_emu=0; iplot_emu<9; iplot_emu++){
      TString name_emu = crname_emu[icr_emu] + plotname_emu[iplot_emu];
      h.ZZ_emu_h[icr_emu][iplot_emu] = new TH1F(name_emu,name_emu,nbins_emu[iplot_emu],blo_emu[iplot_emu],bhi_emu[iplot_emu]);
      h.ZZ_emu_h[icr_emu][iplot_emu]->Sumw2();
    }
  }
  
  //Weight Plots
  h.weights[0]     = new TH1F("LHEweight","LHEweight",40,-20,20);
}


void ZZAna::SortPt(int opt){
  if(opt==0){
    for(int i=0; i<(int)GoodMu.size()-1; i++){
      for(int j=i+1; j<(int)GoodMu.size(); j++){
	if( GoodMu[i].v.Pt() < GoodMu[j].v.Pt() ) swap(GoodMu.at(i),GoodMu.at(j));
      }
    }
  }
  if(opt==1){
    for(int i=0; i<(int)GoodMu_pos.size()-1; i++){
      for(int j=i+1; j<(int)GoodMu_pos.size(); j++){
	if( GoodMu_pos[i].v.Pt() < GoodMu_pos[j].v.Pt() ) swap(GoodMu_pos.at(i),GoodMu_pos.at(j));
      }
    }
  }
  if(opt==2){
    for(int i=0; i<(int)GoodMu_neg.size()-1; i++){
      for(int j=i+1; j<(int)GoodMu_neg.size(); j++){
	if( GoodMu_neg[i].v.Pt() < GoodMu_neg[j].v.Pt() ) swap(GoodMu_neg.at(i),GoodMu_neg.at(j));
      }
    }
  }
  if(opt==3){
    for(int i=0; i<(int)llep.size()-1; i++){
      for(int j=i+1; j<(int)llep.size(); j++){
	if( llep[i].v.Pt() < llep[j].v.Pt() ) swap(llep.at(i),llep.at(j));
      }
    }
  }
  if(opt==4){
    for(int i=0; i<(int)GoodEle.size()-1; i++){
      for(int j=i+1; j<(int)GoodEle.size(); j++){
	if(GoodEle[i].v.Pt() < GoodEle[j].v.Pt() ) swap(GoodEle.at(i),GoodEle.at(j));
      }
    }
  }
  
}

int ZZAna::GenMother(int ind, int mom_ind)
{
  int p_id = GenPart_pdgId[ind];
  int m_id = GenPart_pdgId[mom_ind];
  while(p_id==m_id){
    ind = mom_ind;
    mom_ind = GenPart_genPartIdxMother[ind];
    p_id = GenPart_pdgId[ind];
    m_id = GenPart_pdgId[mom_ind];
  }
  return m_id;
}

//RecoScaleFactors(Muons): 2016UL preVFP
double ZZAna::getScaleFactors_Muons_Reco_preVFP(float eta, float pt){
  double scale_factor=1.0;

  if(fabs(eta)<0.9){
    if( 2.0<pt && pt<3.25 )
      scale_factor = 1.0 ;
    if( 3.25<pt && pt<3.50 )
      scale_factor = 0.9723833248354335 ;
    if( 3.50<pt && pt<3.75 )
      scale_factor = 0.9627713730034183 ;
    if( 3.75<pt && pt<4.0 )
      scale_factor = 0.993506549937225 ;
    if( 4.0<pt && pt<4.5 )
      scale_factor = 1.0145450385374175 ;
    if( 4.5<pt && pt<5.0 )
      scale_factor = 0.9946165315675735 ;
    if( 5.0<pt && pt<6.0 )
      scale_factor = 1.0000394613204677 ;
    if( 6.0<pt && pt<8.0 )
      scale_factor = 0.9989333128307578 ;
    if( 8.0<pt && pt<10.0 )
      scale_factor = 1.0015481462708145 ;
    if( 10.0<pt && pt<15.0 )
      scale_factor = 1.0043436388511215 ;
    if( 15.0<pt && pt<20.0 )
      scale_factor = 0.9964178299208555 ;
    if( 20.0<pt && pt<30.0 )
      scale_factor = 1.0056647298108703 ;
    if( 30.0<pt && pt<200.0 )
      scale_factor = 0.9716164523623753 ;
    
  }
  
  else if(0.9<fabs(eta) && fabs(eta)<1.2){
    if( 2.0<pt && pt<3.25 )
      scale_factor = 1.0 ;
    if( 3.25<pt && pt<3.50 )
      scale_factor = 0.9432095069740233 ;
    if( 3.50<pt && pt<3.75 )
      scale_factor = 0.9305925272487823 ;
    if( 3.75<pt && pt<4.0 )
      scale_factor = 0.9908085986835078 ;
    if( 4.0<pt && pt<4.5 )
      scale_factor = 0.9424219731249323 ;
    if( 4.5<pt && pt<5.0 )
      scale_factor = 1.0154736097598496 ;
    if( 5.0<pt && pt<6.0 )
      scale_factor = 1.0106650703062428 ;
    if( 6.0<pt && pt<8.0 )
      scale_factor = 1.0011631885473775 ;
    if( 8.0<pt && pt<10.0 )
      scale_factor = 0.9919371050807125 ;
    if( 10.0<pt && pt<15.0 )
      scale_factor = 1.0080520133076307 ;
    if( 15.0<pt && pt<20.0 )
      scale_factor = 0.984825408074564 ;
    if( 20.0<pt && pt<30.0 )
      scale_factor = 0.9947357807826922 ;
    if( 30.0<pt && pt<200.0 )
      scale_factor = 1.0 ;
  }

  
  else if(1.2<fabs(eta) && fabs(eta)<2.1){
    if( 2.0<pt && pt<3.25 )
      scale_factor = 1.0 ;
    if( 3.25<pt && pt<3.50 )
      scale_factor = 1.004737916749597 ;
    if( 3.50<pt && pt<3.75 )
      scale_factor = 1.004760337502596 ;
    if( 3.75<pt && pt<4.0 )
      scale_factor = 1.0040852467808588 ;
    if( 4.0<pt && pt<4.5 )
      scale_factor = 0.9914968447284168 ;
    if( 4.5<pt && pt<5.0 )
      scale_factor = 0.9937390545144157 ;
    if( 5.0<pt && pt<6.0 )
      scale_factor = 1.0029312074501653 ;
    if( 6.0<pt && pt<8.0 )
      scale_factor = 1.002298685118436 ;
    if( 8.0<pt && pt<10.0 )
      scale_factor = 0.9647145163097162 ;
    if( 10.0<pt && pt<15.0 )
      scale_factor = 0.9771290781301517 ;
    if( 15.0<pt && pt<20.0 )
      scale_factor = 0.9569388437429848 ;
    if( 20.0<pt && pt<30.0 )
      scale_factor = 1.003006053615245 ;
    if( 30.0<pt && pt<200.0 )
      scale_factor = 1.0 ;
  }
  
  else if(2.1<fabs(eta) && fabs(eta)<2.4){
    if( 2.0<pt && pt<2.5 )
      scale_factor = 1.002066802835317 ;
    if( 2.5<pt && pt<2.75 )
      scale_factor = 0.8764773980483044 ;
    if( 2.75<pt && pt<3.0 )
      scale_factor =0.8606295076084071  ;
    if( 3.0<pt && pt<3.25 )
      scale_factor = 0.9937382769055424 ;
    if( 3.25<pt && pt<3.50 )
      scale_factor = 0.92034046525871 ;
    if( 3.50<pt && pt<3.75 )
      scale_factor = 0.9999997933903989 ;
    if( 3.75<pt && pt<4.0 )
      scale_factor = 0.9212975269110639 ;
    if( 4.0<pt && pt<4.5 )
      scale_factor = 0.9817523842593994 ;
    if( 4.5<pt && pt<5.0 )
      scale_factor = 0.9551181212050973 ;
    if( 5.0<pt && pt<6.0 )
      scale_factor = 0.9999999999020464 ;
    if( 6.0<pt && pt<8.0 )
      scale_factor = 0.9457822923500252 ;
    if( 8.0<pt && pt<10.0 )
      scale_factor = 0.7911488099525276 ;
    if( 10.0<pt && pt<15.0 )
      scale_factor = 1.0003541598832884 ;
    if( 15.0<pt && pt<20.0 )
      scale_factor = 1.0 ;
    if( 20.0<pt && pt<30.0 )
      scale_factor = 1.0 ;
    if( 30.0<pt && pt<200.0 )
      scale_factor = 1.0 ;
  }
  return scale_factor;
}

//RecoScaleFactors(Muons): 2016UL postVFP
double ZZAna::getScaleFactors_Muons_Reco_postVFP(float eta, float pt){
  double scale_factor=1.0;

  if(fabs(eta)<0.9){
    if( 2.0<pt && pt<3.25 )
      scale_factor = 1.0 ;
    if( 3.25<pt && pt<3.50 )
      scale_factor = 1.2735928307942586 ;
    if( 3.50<pt && pt<3.75 )
      scale_factor = 1.038536665892633 ;
    if( 3.75<pt && pt<4.0 )
      scale_factor = 1.0179272663626346 ;
    if( 4.0<pt && pt<4.5 )
      scale_factor = 1.0064396337144235 ;
    if( 4.5<pt && pt<5.0 )
      scale_factor = 0.9952610380311436 ;
    if( 5.0<pt && pt<6.0 )
      scale_factor = 0.9886617311126424 ;
    if( 6.0<pt && pt<8.0 )
      scale_factor = 1.007158929005355 ;
    if( 8.0<pt && pt<10.0 )
      scale_factor = 1.006647891756716 ;
    if( 10.0<pt && pt<15.0 )
      scale_factor = 1.0066631686152465 ;
    if( 15.0<pt && pt<20.0 )
      scale_factor = 1.0058808727351949 ;
    if( 20.0<pt && pt<30.0 )
      scale_factor = 1.0058399541473997 ;
    if( 30.0<pt && pt<200.0 )
      scale_factor = 1.0047417021410996 ;   
  }

    else if(0.9<fabs(eta) && fabs(eta)<1.2){
    if( 2.0<pt && pt<3.50 )
      scale_factor = 1.0 ;
    if( 3.50<pt && pt<3.75 )
      scale_factor = 0.6185253969882053 ;
    if( 3.75<pt && pt<4.0 )
      scale_factor = 0.9469166415007371 ;
    if( 4.0<pt && pt<4.5 )
      scale_factor = 1.0161847682127974 ;
    if( 4.5<pt && pt<5.0 )
      scale_factor = 0.9310484767223539 ;
    if( 5.0<pt && pt<6.0 )
      scale_factor = 1.009897339326928 ;
    if( 6.0<pt && pt<8.0 )
      scale_factor = 1.007027508403878 ;
    if( 8.0<pt && pt<10.0 )
      scale_factor = 1.0062041003034161 ;
    if( 10.0<pt && pt<15.0 )
      scale_factor = 1.0099097794293355 ;
    if( 15.0<pt && pt<20.0 )
      scale_factor = 1.0065290007296392 ;
    if( 20.0<pt && pt<30.0 )
      scale_factor = 1.0047119930493842 ;
    if( 30.0<pt && pt<200.0 )
      scale_factor = 1.0 ;
  }

  
  else if(1.2<fabs(eta) && fabs(eta)<2.1){
    if( 2.0<pt && pt<2.5 )
      scale_factor = 1.0 ;
    if( 2.5<pt && pt<2.75 )
      scale_factor = 1.0 ;
    if( 2.75<pt && pt<3.0 )
      scale_factor = 1.0 ;
    if( 3.0<pt && pt<3.25 )
      scale_factor = 1.0 ;
    if( 3.25<pt && pt<3.50 )
      scale_factor = 0.9455725240692203 ;
    if( 3.50<pt && pt<3.75 )
      scale_factor = 1.0046655328938918 ;
    if( 3.75<pt && pt<4.0 )
      scale_factor = 0.9912556590431282 ;
    if( 4.0<pt && pt<4.5 )
      scale_factor = 0.961200145010141 ;
    if( 4.5<pt && pt<5.0 )
      scale_factor = 1.0035208543548175 ;
    if( 5.0<pt && pt<6.0 )
      scale_factor = 1.0028236667446566 ;
    if( 6.0<pt && pt<8.0 )
      scale_factor = 1.002243448904213 ;
    if( 8.0<pt && pt<10.0 )
      scale_factor = 0.9842686913959845 ;
    if( 10.0<pt && pt<15.0 )
      scale_factor = 1.0045166286019926 ;
    if( 15.0<pt && pt<20.0 )
      scale_factor = 1.003137641796382  ;
    if( 20.0<pt && pt<30.0 )
      scale_factor = 1.0023685314073512 ;
    if( 30.0<pt && pt<200.0 )
      scale_factor = 1.0 ;
  }
  
  else if(2.1<fabs(eta) && fabs(eta)<2.4){
    if( 2.0<pt && pt<2.5 )
      scale_factor = 1.003101641969965 ;
    if( 2.5<pt && pt<2.75 )
      scale_factor = 1.0015080593564991 ;
    if( 2.75<pt && pt<3.0 )
      scale_factor = 1.0020953201750005 ;
    if( 3.0<pt && pt<3.25 )
      scale_factor = 1.0010844193841386 ;
    if( 3.25<pt && pt<3.50 )
      scale_factor = 1.0014976086683698 ;
    if( 3.50<pt && pt<3.75 )
      scale_factor = 1.0013997659539167 ;
    if( 3.75<pt && pt<4.0 )
      scale_factor = 1.0015739076237173 ;
    if( 4.0<pt && pt<4.5 )
      scale_factor = 1.000494083914997 ;
    if( 4.5<pt && pt<5.0 )
      scale_factor = 1.0011590620754274 ;
    if( 5.0<pt && pt<6.0 )
      scale_factor = 1.000262108470239 ;
    if( 6.0<pt && pt<8.0 )
      scale_factor = 1.0004848395299055 ;
    if( 8.0<pt && pt<10.0 )
      scale_factor = 1.000224534964638 ;
    if( 10.0<pt && pt<200.0 )
      scale_factor = 1.0 ;
  }
  return scale_factor;
}

//Calculating effective reconstruction scale factors(Muons)
double ZZAna::getScaleFactors_Muons_Reco(float eta, float pt){
  double scale_factor=1.0;

  scale_factor = (getScaleFactors_Muons_Reco_preVFP(float (eta), float (pt))*20 + getScaleFactors_Muons_Reco_postVFP(float (eta), float (pt))*16)/36;
  return scale_factor;
}

//Muon_Mediumid Scale factors: 2016UL(preVFP and postVFP combined)
double ZZAna::getScaleFactors_Muons_ID(float eta, float pt){
  double scale_factor=1.0;

  if(fabs(eta)<0.9){
    if( 15<pt && pt<20 )
      scale_factor = 0.998179904934898721 ;
    if( 20<pt && pt<25 )
      scale_factor = 0.99830329517085814 ;
    if( 25<pt && pt<30 )
      scale_factor = 0.998278154988515354 ;
    if( 30<pt && pt<40 )
      scale_factor = 0.997988832677700888 ;
    if( 40<pt && pt<50 )
      scale_factor = 0.997956039111710713 ;
    if( 50<pt && pt<60 )
      scale_factor = 0.997480004570111656 ;
    if( 60<pt && pt<500 )
      scale_factor = 0.997534393125491969 ;
  }
  
  else if(0.9<fabs(eta) && fabs(eta)<1.2){
    if( 15<pt && pt<20 )
      scale_factor = 0.993889945950189868 ;
    if( 20<pt && pt<25 )
      scale_factor = 0.994401793012379187 ;
    if( 25<pt && pt<30 )
      scale_factor = 0.995124811857250213 ;
    if( 30<pt && pt<40 )
      scale_factor = 0.994979724623258344 ;
    if( 40<pt && pt<50 )
      scale_factor = 0.995367143033980217 ;
    if( 50<pt && pt<60 )
      scale_factor = 0.995259384451119389 ;
    if( 60<pt && pt<500 )
      scale_factor = 0.994569540695449783 ;
  }

  
  else if(1.2<fabs(eta) && fabs(eta)<2.1){
    if( 15<pt && pt<20 )
      scale_factor = 0.994883065784978982 ;
    if( 20<pt && pt<25 )
      scale_factor = 0.994933568843053617 ;
    if( 25<pt && pt<30 )
      scale_factor = 0.995294446998933546 ;
    if( 30<pt && pt<40 )
      scale_factor = 0.99563985911437114 ;
    if( 40<pt && pt<50 )
      scale_factor = 0.996378944864988569 ;
    if( 50<pt && pt<60 )
      scale_factor = 0.996014481509890359 ;
    if( 60<pt && pt<500 )
      scale_factor = 0.996157333877533802 ;
  }
  
  else if(2.1<fabs(eta) && fabs(eta)<2.4){
    if( 15<pt && pt<20 )
      scale_factor = 0.976747405936561175 ;
    if( 20<pt && pt<25 )
      scale_factor = 0.975688840459463691 ;
    if( 25<pt && pt<30 )
      scale_factor = 0.975987906281726714 ;
    if( 30<pt && pt<40 )
      scale_factor = 0.976947217047071437 ;
    if( 40<pt && pt<50 )
      scale_factor = 0.976347433680206378 ;
    if( 50<pt && pt<60 )
      scale_factor = 0.976876997398022251 ;
    if( 60<pt && pt<500 )
      scale_factor = 0.975152286467743545 ;
  }
  return scale_factor;
}

//MuonTightIso_Scale factors: 2016UL(preVFP and postVFP combined)
double ZZAna::getScaleFactors_Muons_Iso(float eta, float pt){
  double scale_factor=1.0;

  if(fabs(eta)<0.9){
    if( 15<pt && pt<20 )
      scale_factor = 0.987778547702043985 ;
    if( 20<pt && pt<25 )
      scale_factor = 0.991494879152131303 ;
    if( 25<pt && pt<30 )
      scale_factor = 0.988953382957504123 ;
    if( 30<pt && pt<40 )
      scale_factor = 0.994468703426339307 ;
    if( 40<pt && pt<50 )
      scale_factor = 0.996726918748909729 ;
    if( 50<pt && pt<60 )
      scale_factor = 0.997054906988794509 ;
    if( 60<pt && pt<500 )
      scale_factor = 0.998868311321585267 ;
  }
  
  else if(0.9<fabs(eta) && fabs(eta)<1.2){
    if( 15<pt && pt<20 )
      scale_factor = 0.992284497975155499 ;
    if( 20<pt && pt<25 )
      scale_factor = 0.994445238716516133 ;
    if( 25<pt && pt<30 )
      scale_factor = 1.00015108444872491 ;
    if( 30<pt && pt<40 )
      scale_factor = 0.998685892116171758 ;
    if( 40<pt && pt<50 )
      scale_factor = 0.997541391524236265 ;
    if( 50<pt && pt<60 )
      scale_factor = 0.998522557263133659 ;
    if( 60<pt && pt<500 )
      scale_factor = 0.999511515676385964 ;
  }

  
  else if(1.2<fabs(eta) && fabs(eta)<2.1){
    if( 15<pt && pt<20 )
      scale_factor = 1.00060948451767384 ;
    if( 20<pt && pt<25 )
      scale_factor = 1.00247062995998437 ;
    if( 25<pt && pt<30 )
      scale_factor = 1.00337960717140229 ;
    if( 30<pt && pt<40 )
      scale_factor = 1.00037891036521942 ;
    if( 40<pt && pt<50 )
      scale_factor = 0.999359722227062264 ;
    if( 50<pt && pt<60 )
      scale_factor = 0.998656508751751715 ;
    if( 60<pt && pt<500 )
      scale_factor = 0.999743845471649295 ;
  }
  
  else if(2.1<fabs(eta) && fabs(eta)<2.4){
    if( 15<pt && pt<20 )
      scale_factor = 1.01090285270219526 ;
    if( 20<pt && pt<25 )
      scale_factor = 1.00762322825897188 ;
    if( 25<pt && pt<30 )
      scale_factor = 1.00462853050708389 ;
    if( 30<pt && pt<40 )
      scale_factor = 1.00190767245697021 ;
    if( 40<pt && pt<50 )
      scale_factor = 1.00102989223772743 ;
    if( 50<pt && pt<60 )
      scale_factor = 1.00050965081981347 ;
    if( 60<pt && pt<500 )
      scale_factor = 1.00069328135112667 ;
  }
  return scale_factor;
}

//MuonTriggerEfficiency
double ZZAna::getEfficiency_Muons_Trigger(float eta, float pt){
  double eff_factor = 1.0;

  if(fabs(eta)<0.9){
    if( 26<pt && pt<30 )
      eff_factor = 0.891430358091990116 ;
    if( 30<pt && pt<40 )
      eff_factor = 0.916997505558861614 ;
    if( 40<pt && pt<50 )
      eff_factor = 0.929323958026038288 ;
    if( 50<pt && pt<60 )
      eff_factor = 0.932949966854519364 ;
    if( 60<pt && pt<120 )
      eff_factor = 0.932453400558895584 ;
    if( 120<pt && pt<200 )
      eff_factor = 0.924933989842732784 ;
  }
  
  else if(0.9<fabs(eta) && fabs(eta)<1.2){
    if( 26<pt && pt<30 )
      eff_factor = 0.881071481439802406 ;
    if( 30<pt && pt<40 )
      eff_factor = 0.91953304078843856 ;
    if( 40<pt && pt<50 )
      eff_factor = 0.935559259520636677 ;
    if( 50<pt && pt<60 )
      eff_factor = 0.939758843845791336 ;
    if( 60<pt && pt<120 )
      eff_factor = 0.939229230086008671 ;
    if( 120<pt && pt<200 )
      eff_factor = 0.926802284187740799 ;
  }

  
  else if(1.2<fabs(eta) && fabs(eta)<2.1){
    if( 26<pt && pt<30 )
      eff_factor = 0.819473167260487911 ;
    if( 30<pt && pt<40 )
      eff_factor = 0.860544012652503132 ;
    if( 40<pt && pt<50 )
      eff_factor = 0.881642282009124756 ;
    if( 50<pt && pt<60 )
      eff_factor = 0.887126167615254757 ;
    if( 60<pt && pt<120 )
      eff_factor = 0.88759789864222205 ;
    if( 120<pt && pt<200 )
      eff_factor = 0.888963898022969601 ;
  }
  
  else if(2.1<fabs(eta) && fabs(eta)<2.4){
    if( 26<pt && pt<30 )
      eff_factor = 0.711451828479766846 ;
    if( 30<pt && pt<40 )
      eff_factor = 0.774474415514204262 ;
    if( 40<pt && pt<50 )
      eff_factor = 0.804473486211564781 ;
    if( 50<pt && pt<60 )
      eff_factor = 0.81278272469838464 ;
    if( 60<pt && pt<120 )
      eff_factor = 0.818142698870764851 ;
    if( 120<pt && pt<200 )
      eff_factor = 0.808433486355675579 ;
  }
  return eff_factor;
}

double ZZAna::getScaleFactors_Electrons_Reco(float eta, float pt){
  double scale_factor = 1.0;
  
  if(10<pt && pt<20){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 1.01448949178059888;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 0.990329477522108315;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 1.17230590184529615;
    if(-1.444<eta && eta<-1.0) 
      scale_factor = 0.991350465350680832;
    if(-1.0<eta && eta<0.0) 
      scale_factor = 1.04162266519334579;
    if(0.0<eta && eta<1.0) 
      scale_factor = 1.04162266519334579;
    if(1.0<eta && eta<1.444) 
      scale_factor = 0.991350465350680832;
    if(1.444<eta && eta<1.566) 
      scale_factor = 1.17230590184529615;
    if(1.566<eta && eta<2.0) 
      scale_factor = 0.990329477522108315;
    if(2.0<eta && eta<2.4) 
      scale_factor = 1.01448949178059888;
  }
  if(20<pt && pt<45){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 1.01393263207541562;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 0.99270998107062447;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 0.971759047773149254;
    if(-1.444<eta && eta<-1.0) 
      scale_factor = 0.986505568027496338;
    if(-1.0<eta && eta<0.5) 
      scale_factor = 0.982211887836456299;
    if(-0.5<eta && eta<0.0) 
      scale_factor = 0.978960871696472168;
    if(0.0<eta && eta<0.5) 
      scale_factor = 0.985260983308156368;
    if(0.5<eta && eta<1.0) 
      scale_factor = 0.986813724040985107;
    if(1.0<eta && eta<1.444) 
      scale_factor = 0.987445645862155441;
    if(1.444<eta && eta<1.566) 
      scale_factor = 0.973707669311099533;
    if(1.566<eta && eta<2.0) 
      scale_factor = 0.991490264733632443;
    if(2.0<eta && eta<2.4) 
      scale_factor = 0.99618503782484269;
    }
    
  if(45<pt && pt<75){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 1.00442620780732894;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 0.991372115082211014;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 0.961046099662780762;
    if(-1.444<eta && eta<-1.0) 
      scale_factor = 0.986727171474032883;
    if(-1.0<eta && eta<0.5) 
      scale_factor = 0.983468506071302651;
    if(-0.5<eta && eta<0.0) 
      scale_factor = 0.98229103618197966;
    if(0.0<eta && eta<0.5) 
      scale_factor = 0.986387530962626102;
    if(0.5<eta && eta<1.0) 
      scale_factor = 0.987345370981428383;
    if(1.0<eta && eta<1.444) 
      scale_factor = 0.987421433130900028;
    if(1.444<eta && eta<1.566) 
      scale_factor = 0.967947648631201862;
    if(1.566<eta && eta<2.0) 
      scale_factor = 0.989836785528394936;
    if(2.0<eta && eta<2.4) 
      scale_factor = 0.992115318775177002;
    }

  if(75<pt && pt<100){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 0.993540048599243164;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 0.979585521750979904;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 1.00285185707940006;
    if(-1.444<eta && eta<-1.0) 
      scale_factor = 0.982317573494381424;
    if(-1.0<eta && eta<0.5) 
      scale_factor = 0.976644403404659744;
    if(-0.5<eta && eta<0.0) 
      scale_factor = 0.982822967900170208;
    if(0.0<eta && eta<0.5) 
      scale_factor = 0.982822967900170208;
    if(0.5<eta && eta<1.0) 
      scale_factor = 0.976644403404659744;
    if(1.0<eta && eta<1.444) 
      scale_factor = 0.982317573494381424;
    if(1.444<eta && eta<1.566) 
      scale_factor = 1.00285185707940006;
    if(1.566<eta && eta<2.0) 
      scale_factor = 0.979585521750979904;
    if(2.0<eta && eta<2.4) 
      scale_factor = 0.993540048599243164;
    }
  if(100<pt && pt<1000){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 1.01000316937764478;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 1.00773619280921078;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 1.00749973456064867;
    if(-1.444<eta && eta<-1.0) 
      scale_factor = 0.985203742980957031;
    if(-1.0<eta && eta<0.5) 
      scale_factor = 0.98755870925055611;
    if(-0.5<eta && eta<0.0) 
      scale_factor = 0.988535126050313351;
    if(0.0<eta && eta<0.5) 
      scale_factor = 0.988535126050313351;
    if(0.5<eta && eta<1.0) 
      scale_factor = 0.98755870925055611;
    if(1.0<eta && eta<1.444) 
      scale_factor = 0.985203742980957031;
    if(1.444<eta && eta<1.566) 
      scale_factor = 1.00749973456064867;
    if(1.566<eta && eta<2.0) 
      scale_factor = 1.00773619280921078;
    if(2.0<eta && eta<2.4) 
      scale_factor = 1.01000316937764478;
    }
  return scale_factor;
}

double ZZAna::getScaleFactors_Electrons_IDIso(float eta, float pt){
  double scale_factor = 1.0;
  
  if(10<pt && pt<20){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 1.06316983699798584;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 1.045610162946913;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 1;
    if(-1.444<eta && eta<-0.8) 
      scale_factor = 1.02309005790286589;
    if(-0.8<eta && eta<0.0) 
      scale_factor = 0.972585002581278446;
    if(0.0<eta && eta<0.8) 
      scale_factor = 0.97979246907764006;
    if(0.8<eta && eta<1.444) 
      scale_factor = 1.03856579462687182;
    if(1.444<eta && eta<1.566) 
      scale_factor = 1;
    if(1.566<eta && eta<2.0) 
      scale_factor = 1.03034557236565494;
    if(2.0<eta && eta<2.4) 
      scale_factor = 1.05637967586517334;
  }

  if(20<pt && pt<35){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 1.02417140536838103;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 1.00863300429450131;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 1;
    if(-1.444<eta && eta<-0.8) 
      scale_factor = 0.990166293250189899;
    if(-0.8<eta && eta<0.0) 
      scale_factor = 0.968974007500542522;
    if(0.0<eta && eta<0.8) 
      scale_factor = 0.983436869250403523;
    if(0.8<eta && eta<1.444) 
      scale_factor = 0.995318664444817425;
    if(1.444<eta && eta<1.566) 
      scale_factor = 1;
    if(1.566<eta && eta<2.0) 
      scale_factor = 0.999914500448438881;
    if(2.0<eta && eta<2.4) 
      scale_factor = 1.01638425721062564;
  }

  if(35<pt && pt<50){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 1.01238840156131316;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 1.0072176456451416;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 1;
    if(-1.444<eta && eta<-0.8) 
      scale_factor = 0.986762649483150955;
    if(-0.8<eta && eta<0.0) 
      scale_factor = 0.970998260709974526;
    if(0.0<eta && eta<0.8) 
      scale_factor = 0.982610225677490234;
    if(0.8<eta && eta<1.444) 
      scale_factor = 0.990672462516360763;
    if(1.444<eta && eta<1.566) 
      scale_factor = 1;
    if(1.566<eta && eta<2.0) 
      scale_factor = 1.00322821405198837;
    if(2.0<eta && eta<2.4) 
      scale_factor = 1.00469050142500138;
  }

  if(50<pt && pt<100){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 1.00294347604115797;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 1.00411701202392578;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 1;
    if(-1.444<eta && eta<-0.8) 
      scale_factor = 0.98661394913991296;
    if(-0.8<eta && eta<0.0) 
      scale_factor = 0.970484183894263386;
    if(0.0<eta && eta<0.8) 
      scale_factor = 0.979886399375067829;
    if(0.8<eta && eta<1.444) 
      scale_factor = 0.98747599787182283;
    if(1.444<eta && eta<1.566) 
      scale_factor = 1;
    if(1.566<eta && eta<2.0) 
      scale_factor = 1.0015653106901381;
    if(2.0<eta && eta<2.4) 
      scale_factor = 0.997389793395996094;
  }

  if(100<pt && pt<1000){
    if(-2.4<eta && eta<-2.1) 
      scale_factor = 1.02094837029774976;
    if(-2.1<eta && eta<-1.566) 
      scale_factor = 1.02457796202765561;
    if(-1.566<eta && eta<-1.444) 
      scale_factor = 1;
    if(-1.444<eta && eta<-0.8) 
      scale_factor = 1.00173119703928637;
    if(-0.8<eta && eta<0.0) 
      scale_factor = 0.995304425557454464;
    if(0.0<eta && eta<0.8) 
      scale_factor = 1.0050323406855266;
    if(0.8<eta && eta<1.444) 
      scale_factor = 1.0106508731842041;
    if(1.444<eta && eta<1.566) 
      scale_factor = 1;
    if(1.566<eta && eta<2.0) 
      scale_factor = 0.992624911997053383;
    if(2.0<eta && eta<2.4) 
      scale_factor = 1.02222976419660783;
  }
  return scale_factor;
    
}
