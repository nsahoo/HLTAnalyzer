//#include "makeMultiPanelCanvas.C"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "HLTeff.h"

#include <vector>

void HLTeff_HighPt()
{
    TFile* file0 = TFile::Open("HLTTree_newEPOS_MB_noJEC.root");
    
    //Trigger threshold
    int PixelThreshold = 75;
    int MultiplicityThreshold = 110;
    float HFSumThreshold = 55;
    float L1JetThreshold = 12;
    float PtThreshold[2] = {8,16};
    int NPt = 2;
    
    //initialize eff plots
    int Nplot = 2;
    TH1D* h = new TH1D("","",60,0,60);
    TH1D* hFull = (TH1D*)h.Clone();
    TH1D* hFullHFsum = (TH1D*)h.Clone();
    TH1D* hHMThres1[100]; //HM+HighPt
    TH1D* hHMThres2[100]; //HFsum+HighPt
    TH1D* hHMFullThres1[100]; //HM+HighPt
    TH1D* hHMFullThres2[100]; //HFsum+HighPt
    TH1D* hHMJet; //L1Jet for HM+HighPt
    TH1D* hHFJet; //L1Jet for HFsum+HighPt
    hHMJet = (TH1D*)h.Clone();
    hHFJet = (TH1D*)h.Clone();
    for(int i=0;i<Nplot;i++)
    {
        hHMThres1[i] = (TH1D*)h.Clone();
        hHMThres2[i] = (TH1D*)h.Clone();
        hHMFullThres1[i] = (TH1D*)h.Clone();
        hHMFullThres2[i] = (TH1D*)h.Clone();
    }
    
    TGraphAsymmErrors* grHMThres1[100];
    TGraphAsymmErrors* grHMThres2[100];
    TGraphAsymmErrors* grHMFullThres1[100];
    TGraphAsymmErrors* grHMFullThres2[100];
    TGraphAsymmErrors* grHMJet = new TGraphAsymmErrors();
    TGraphAsymmErrors* grHFJet = new TGraphAsymmErrors();
    for(int i=0;i<Nplot;i++)
    {
        grHMThres1[i] = new TGraphAsymmErrors();
        grHMThres2[i] = new TGraphAsymmErrors();
        grHMFullThres1[i] = new TGraphAsymmErrors();
        grHMFullThres2[i] = new TGraphAsymmErrors();
    }
    
    //Get Tree branch
    TNtuple* Tree = (TNtuple*)file0->Get("HLTinfo/HLTeff");
    
    int           Run;
    ULong64_t     Event;
    int           LumiBlock;
    int Ntrkoffline;
    int NtrkFull;
    int NtrkPixel;
    double OfflineVtxZ;
    double HLTVtxZ;
    double OfflineLeadingPt;
    double HLTLeadingPt;
    double HFsumET;
    
    Tree->SetBranchAddress("Run",&Run);
    Tree->SetBranchAddress("Event",&Event);
    Tree->SetBranchAddress("LumiBlock",&LumiBlock);
    Tree->SetBranchAddress("Ntrkoffline",&Ntrkoffline);
    Tree->SetBranchAddress("NtrkPixel",&NtrkPixel);
    Tree->SetBranchAddress("NtrkFull",&NtrkFull);
    Tree->SetBranchAddress("OfflineVtxZ",&OfflineVtxZ);
    Tree->SetBranchAddress("HLTVtxZ",&HLTVtxZ);
    Tree->SetBranchAddress("OfflineLeadingPt",&OfflineLeadingPt);
    Tree->SetBranchAddress("HLTLeadingPt",&HLTLeadingPt);
    Tree->SetBranchAddress("HFsumET",&HFsumET);
    
    //Get hlt decision branch
    TTree* hltroot = (TTree*)file0->Get("hltbitanalysis/HltTree");
    SetHlttreestatus(hltroot);
    setHltTreeBranch(hltroot);
    
    //loop for histograms
    int nentry = Tree->GetEntries();
    for(int it=0;it<nentry;it++)
    {
        Tree->GetEntry(it);
        hltroot->GetEntry(it);

        if(it%10000==0) cout<<"processing "<<it<<"/"<<nentry<<" event"<<endl;
        
        if(fabs(OfflineVtxZ)>15) continue;
        
        //fill LeadingPt distribution without any cut for HM+HighPt
        if(Ntrkoffline>=120) hFull->Fill(OfflineLeadingPt);
        //fill LeadingPt distribution with cuts for HM+HighPt
        if(NtrkFull>=MultiplicityThreshold && NtrkPixel>=PixelThreshold)
        {
            for(int j=0;j<NPt;j++)
            {
                if(fabs(HLTVtxZ)>15) continue;
                if(HLTLeadingPt<PtThreshold[j]) continue;
                if(Ntrkoffline<120) continue;
                hHMThres1[j]->Fill(OfflineLeadingPt);
            }
        }
        //fill LeadingPt distribution with L1 jet for HM+HighPt
        if(Ntrkoffline>=120 && Run==HLT_Run && Event==HLT_Event && LumiBlock==HLT_LumiBlock && HLT_L1SingleJet8==1 && fabs(HLTVtxZ)<=15) hHMJet->Fill(OfflineLeadingPt);
        //fill leadingpt distribution for full eff
        if(Ntrkoffline>=120 && Run==HLT_Run && Event==HLT_Event && LumiBlock==HLT_LumiBlock && fabs(HLTVtxZ)<=15 && HLTLeadingPt>PtThreshold[0]) hHMFullThres1[0]->Fill(OfflineLeadingPt);
        if(Ntrkoffline>=120 && Run==HLT_Run && Event==HLT_Event && LumiBlock==HLT_LumiBlock && fabs(HLTVtxZ)<=15 && HLTLeadingPt>PtThreshold[1] && HLT_L1SingleJet8==1) hHMFullThres1[1]->Fill(OfflineLeadingPt);
        
        //fill LeadingPt distribution without any cut for HFsum+HighPt
        if(HFsumET>HFSumThreshold) hFullHFsum->Fill(OfflineLeadingPt);
        //fill LeadingPt distribution with cuts for HFsum+HighPt
        if(HFsumET>HFSumThreshold)
        {
            for(int j=0;j<NPt;j++)
            {
                if(fabs(HLTVtxZ)>15) continue;
                if(HLTLeadingPt<PtThreshold[j]) continue;
                hHMThres2[j]->Fill(OfflineLeadingPt);
            }
        }
        //fill LeadingPt distribution with L1 jet for HFsum+HighPt
        if(HFsumET>HFSumThreshold && Run==HLT_Run && Event==HLT_Event && LumiBlock==HLT_LumiBlock && HLT_L1SingleJet8==1 && fabs(HLTVtxZ)<=15) hHFJet->Fill(OfflineLeadingPt);
        //fill leadingpt distribution for full eff
        if(HFsumET>HFSumThreshold && Run==HLT_Run && Event==HLT_Event && LumiBlock==HLT_LumiBlock && fabs(HLTVtxZ)<=15 && HLTLeadingPt>PtThreshold[0]) hHMFullThres2[0]->Fill(OfflineLeadingPt);
        if(HFsumET>HFSumThreshold && Run==HLT_Run && Event==HLT_Event && LumiBlock==HLT_LumiBlock && HLT_L1SingleJet8==1 && fabs(HLTVtxZ)<=15 && HLTLeadingPt>PtThreshold[1]) hHMFullThres2[1]->Fill(OfflineLeadingPt);
        
    }
    
    //make eff turn-ons
    for(int i=0;i<NPt;i++)
    {
        grHMThres1[i]->Divide(hHMThres1[i],hFull);
        grHMThres2[i]->Divide(hHMThres2[i],hFullHFsum);
        grHMFullThres1[i]->Divide(hHMFullThres1[i],hFull);
        grHMFullThres2[i]->Divide(hHMFullThres2[i],hFullHFsum);
    }
    grHMJet->Divide(hHMJet,hFull);
    grHFJet->Divide(hHFJet,hFullHFsum);
    
    //plot eff turn-ons
    TH1D* HLTeff = new TH1D("HLTeff","HLTeff",60,0,60);
    HLTeff->GetXaxis()->SetTitle("pT_{offline}^{Leading}");
    HLTeff->GetYaxis()->SetTitle("HLT eff");
    HLTeff->SetTitle("");
    
    TH1D* L1eff = new TH1D("L1eff","L1eff",60,0,60);
    L1eff->GetXaxis()->SetTitle("pT_{offline}^{Leading}");
    L1eff->GetYaxis()->SetTitle("L1 eff");
    L1eff->SetTitle("");
   
    TH1D* Fulleff = new TH1D("Fulleff","Fulleff",60,0,60);
    Fulleff->GetXaxis()->SetTitle("pT_{offline}^{Leading}");
    Fulleff->GetYaxis()->SetTitle("L1+HLT eff");
    Fulleff->SetTitle("");
    
    //add lines at where the trigger is expected to be efficient
    TLine* l = new TLine(0,1,100,1);
    l->SetLineStyle(3);
    
    TLine* l1 = new TLine(8,0,8,1.05);
    TLine* l2 = new TLine(16,0,16,1.05);
    l1->SetLineStyle(7);
    l2->SetLineStyle(7);
    
    //create canvas
    TCanvas* c1 = new TCanvas("HM","HM",600,600);
    TCanvas* c2 = new TCanvas("HFsum","HFsum",600,600);
    TCanvas* c3 = new TCanvas("HMJet","HMJet",600,600);
    TCanvas* c4 = new TCanvas("HFJet","HFJet",600,600);
    TCanvas* c5 = new TCanvas("HMFull","HMFull",600,600);
    TCanvas* c6 = new TCanvas("HFsumFull","HFsumFull",600,600);
    
    //HM+HighPt HLT eff
    c1->cd();
    HLTeff->Draw();
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(grHMThres1[i],Form("Ntrkonline>=%d",MultiplicityThreshold),"");
    leg->SetFillStyle(0);
    for(int i=0;i<NPt;i++)
    {
        grHMThres1[i]->SetMarkerColor(i+1);
        grHMThres1[i]->SetMarkerSize(0.8);
        grHMThres1[i]->Draw("PSAME");
        leg->AddEntry(grHMThres1[i],Form("HLT_Pt%.1f",PtThreshold[i]),"p");
    }
    leg->Draw();
    l->Draw("LSAME");
    l1->Draw("LSAME");
    l2->Draw("LSAME");
    
    //HFsum+HighPt HLT eff
    c2->cd();
    HLTeff->Draw();
    TLegend* leg1 = new TLegend(0.7,0.7,0.9,0.9);
    leg1->AddEntry(grHMThres2[i],Form("HFsumET>%.1f",HFSumThreshold),"");
    leg1->SetFillStyle(0);
    for(int i=0;i<NPt;i++)
    {
        grHMThres2[i]->SetMarkerColor(i+1);
        grHMThres2[i]->SetMarkerSize(0.8);
        grHMThres2[i]->Draw("PSAME");
        leg1->AddEntry(grHMThres2[i],Form("HLT_Pt%.1f",PtThreshold[i]),"p");
    }
    leg1->Draw();
    l->Draw("LSAME");
    l1->Draw("LSAME");
    l2->Draw("LSAME");
    
    //HM L1 jet eff
    c3->cd();
    L1eff->Draw();
    TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
    leg2->AddEntry(grHMJet,Form("Ntrkoffline>=%d",120),"");
    leg2->SetFillStyle(0);
    leg2->AddEntry(grHMJet,Form("L1_SingleJet%.0f",L1JetThreshold),"p");
    leg2->Draw();
    grHMJet->Draw("PSAME");
    l2->Draw("LSAME");
    
    //HF L1 jet eff
    c4->cd();
    L1eff->Draw();
    TLegend* leg3 = new TLegend(0.7,0.7,0.9,0.9);
    leg3->AddEntry(grHFJet,Form("HFsumET>%.1f",HFSumThreshold),"");
    leg3->SetFillStyle(0);
    leg3->AddEntry(grHFJet,Form("L1_SingleJet%.0f",L1JetThreshold),"p");
    leg3->Draw();
    grHFJet->Draw("PSAME");
    l2->Draw("LSAME");
    
    //HM+HighPt HLT eff
    c5->cd();
    Fulleff->Draw();
    for(int i=0;i<NPt;i++)
    {
        grHMFullThres1[i]->SetMarkerColor(i+1);
        grHMFullThres1[i]->SetMarkerSize(0.8);
        grHMFullThres1[i]->Draw("PSAME");
    }
    leg->Draw();
    l->Draw("LSAME");
    l1->Draw("LSAME");
    l2->Draw("LSAME");
    
    //HFsum+HighPt HLT eff
    c6->cd();
    HLTeff->Draw();
    for(int i=0;i<NPt;i++)
    {
        grHMFullThres2[i]->SetMarkerColor(i+1);
        grHMFullThres2[i]->SetMarkerSize(0.8);
        grHMFullThres2[i]->Draw("PSAME");
    }
    leg1->Draw();
    l->Draw("LSAME");
    l1->Draw("LSAME");
    l2->Draw("LSAME");

}

