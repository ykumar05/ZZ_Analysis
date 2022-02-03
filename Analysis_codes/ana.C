#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
/*
This is a driver script.
It decides which code to run over which sample, the names
of output files and so on.
*/


void ana(int sample=1)
{
  const char *hstfilename, *sumfilename;
  //Declare a chain for input files.
  TChain *chain = new TChain("Events");
  //Declare an instance of our code class
  ZZAna m_selec;
  
  if(sample==0){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/phazarik/nanoAOD_UL2016_Data/SingleMuon_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v1/*.root"); //B
    chain->Add("/home/work/phazarik/nanoAOD_UL2016_Data/SingleMuon_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v1/*.root"); //C
    chain->Add("/home/work/phazarik/nanoAOD_UL2016_Data/SingleMuon_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v1/*.root");//D
    chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016E_UL2016_Data/*.root");//E
    chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016F_UL2016_Data/*.root");//F
    chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016G_UL2016_Data/*.root");//G
    chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016H_UL2016_Data/*.root");//H
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_data.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_data.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  
  if(sample==1){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/ykumar/Samples/2016/ZZTo4L_20UL16_nanoAODv9_MC/*");
    //chain->Add("/home/work/ykumar/Samples/2016/ZZTo4L_20UL16_nanoAODv9_MC/1F555B22-8D51-134A-9A22-56EBBB62A9E6.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_ZZ.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_ZZ.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  if(sample==2){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/ykumar/Samples/2016/TTZToLLNuNu_20UL16_nanoAODv9_MC/*.root");
    //chain->Add("/home/work/ykumar/Samples/2016/TTZToLLNuNu_20UL16_nanoAODv9_MC/11438748-9954-664F-A675-4E305B6A5A91.root");
    
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_TTZ.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_TTZ.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  if(sample==3){
    //Add one file to chain. This is the input file.
    //chain->Add("/home/work/phazarik/nanoAOD_UL2016_MC/WZTo3LNU_UL2016_MC/EE8014A0-3EA4-B148-BC2F-83B6FD1E233A.root");
    chain->Add("/home/work/phazarik/nanoAOD_UL2016_MC/WZTo3LNU_UL2016_MC/*.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_WZ.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_WZ.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  if(sample==4){
    //Add one file to chain. This is the input file.
    //chain->Add("/home/work/phazarik/nanoAOD_UL2016_MC/TTWJetsToLNu_UL2016_MC/0CFB3EF2-F10F-694E-AAF1-678D040FF732.root");
    chain->Add("/home/work/phazarik/nanoAOD_UL2016_MC/TTWJetsToLNu_UL2016_MC/*.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_TTW.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_TTW.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(0); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }

  if(sample==-1){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/phazarik/nanoAOD_UL2016_Data/SingleMuon_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v1/*.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_DataC.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_DataC.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  if(sample==-2){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/phazarik/nanoAOD_UL2016_Data/SingleMuon_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v1/*.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_DataD.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_DataD.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  if(sample==-3){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016E_UL2016_Data/*.root");
    //chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016F_UL2016_Data/EC12C4E4-4F1F-4441-86B8-5F15678EDC79.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_DataE.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_DataE.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  if(sample==-4){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016F_UL2016_Data/*.root");
    //chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016F_UL2016_Data/EC12C4E4-4F1F-4441-86B8-5F15678EDC79.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_DataF.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_DataF.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  if(sample==-5){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016G_UL2016_Data/*.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_DataG.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_DataG.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  if(sample==-6){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/ykumar/Samples/2016/SingleMuon_Run2016H_UL2016_Data/*.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_DataH.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_DataH.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  if(sample==-7){
    //Add one file to chain. This is the input file.
    chain->Add("/home/work/phazarik/nanoAOD_UL2016_Data/SingleMuon_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v1/*.root");
    //Set names of output files.
    hstfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/hst_DataB.root";
    sumfilename = "/home/work/ykumar/ZZ_Analysis/hst_files/sum_DataB.txt";
    //Set some options.. data is 0 since this input file is simulation.
    m_selec.SetData(1); //0 - running over MC, 1 - running over Data
    m_selec.SetYear(2016);
  }
  
  std::cout<<"Output files are "<<hstfilename<<" and "<<sumfilename<<std::endl;
  // Set some more options.. set the output file names.
  m_selec.SetHstFileName(hstfilename);
  m_selec.SetSumFileName(sumfilename);
  m_selec.SetVerbose(10);//set verbosity level for output.
  // Call the process function which runs the code.
  chain->Process(&m_selec);
}
