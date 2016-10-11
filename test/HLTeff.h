#include "TTree.h"

//HltInfo
int           HLT_Run;
ULong64_t     HLT_Event;
int           HLT_LumiBlock;
int          HLT_PAFull_HM120;
int          HLT_PAFull_HM150;
int          HLT_PAFull_HM185;
int          HLT_PAFull_HM220;
int          HLT_PAFull_HM250;
int         HLT_PAFull_HM280;
int         HLT_PAFull_HM110_pt8;
int         HLT_PAFull_HM110_pt16;
int         HLT_PAFull_HFSumEt005_pt8;
int         HLT_PAFull_HFSumEt005_pt16;
int         HLT_L1SingleJet8;
int         HLT_L1SingleJet12;
int         HLT_L1SingleJet16;
void setHltTreeBranch(TTree* hltroot)
{
    hltroot->SetBranchAddress("Run",&HLT_Run);
    hltroot->SetBranchAddress("Event",&HLT_Event);
    hltroot->SetBranchAddress("LumiBlock",&HLT_LumiBlock);
    hltroot->SetBranchAddress("HLT_PAFullTracks_Multiplicity120_v1",&HLT_PAFull_HM120);
    hltroot->SetBranchAddress("HLT_PAFullTracks_Multiplicity150_v1",&HLT_PAFull_HM150);
    hltroot->SetBranchAddress("HLT_PAFullTracks_Multiplicity185_v1",&HLT_PAFull_HM185);
    hltroot->SetBranchAddress("HLT_PAFullTracks_Multiplicity220_v1",&HLT_PAFull_HM220);
    hltroot->SetBranchAddress("HLT_PAFullTracks_Multiplicity250_v1",&HLT_PAFull_HM250);
    hltroot->SetBranchAddress("HLT_PAFullTracks_Multiplicity280_v1",&HLT_PAFull_HM280);
    hltroot->SetBranchAddress("HLT_PAFullTracks_Multiplicity110_HighPt8_v1",&HLT_PAFull_HM110_pt8);
    hltroot->SetBranchAddress("HLT_PAFullTracks_Multiplicity110_HighPt16_v1",&HLT_PAFull_HM110_pt16);
    hltroot->SetBranchAddress("HLT_PAFullTracks_HFSumEt005_HighPt8_v1",&HLT_PAFull_HFSumEt005_pt8);
    hltroot->SetBranchAddress("HLT_PAFullTracks_HFSumEt005_HighPt16_v1",&HLT_PAFull_HFSumEt005_pt16);
    hltroot->SetBranchAddress("HLT_L1SingleJet8_v1",&HLT_L1SingleJet8);
    hltroot->SetBranchAddress("HLT_L1SingleJet12_v1",&HLT_L1SingleJet12);
    hltroot->SetBranchAddress("HLT_L1SingleJet16_v1",&HLT_L1SingleJet16);
}

void SetHlttreestatus(TTree * hltroot)
{
    hltroot->SetBranchStatus("*",0);
    hltroot->SetBranchStatus("Event",1);
    hltroot->SetBranchStatus("LumiBlock",1);
    hltroot->SetBranchStatus("Run",1);
    
    hltroot->SetBranchStatus("HLT_PAFullTracks_*",1);
    hltroot->SetBranchStatus("HLT_L1SingleJet*",1);
}

