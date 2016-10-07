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

#include <vector>

void HLTeff_HighPt()
{
    TFile* file0 = TFile::Open("HLTeff.root");
    
    //Trigger threshold
    int PixelThreshold = 40;
    float PtThreshold[2] = {8,16};
    int MultiplicityThreshold = 110;
    float HFSumThreshold = 55;
    int NPt = 2;
    
    //initialize eff plots
    int Nplot = 2;
    TH1D* h = new TH1D("","",60,0,60);
    TH1D* hFull = (TH1D*)h.Clone();
    TH1D* hFullHFsum = (TH1D*)h.Clone();
    TH1D* hHMThres1[100]; //HM+HighPt
    TH1D* hHMThres2[100]; //HFsum+HighPt
    for(int i=0;i<Nplot;i++)
    {
        hHMThres1[i] = (TH1D*)h.Clone();
        hHMThres2[i] = (TH1D*)h.Clone();
    }
    
    TGraphAsymmErrors* grHMThres1[100];
    TGraphAsymmErrors* grHMThres2[100];
    for(int i=0;i<Nplot;i++)
    {
        grHMThres1[i] = new TGraphAsymmErrors();
        grHMThres2[i] = new TGraphAsymmErrors();
    }
    
    //Get Tree branch
    TNtuple* Tree = (TNtuple*)file0->Get("HLTinfo/HLTeff");
    
    int Ntrkoffline;
    int NtrkFull;
    int NtrkPixel;
    double OfflineVtxZ;
    double HLTVtxZ;
    double OfflineLeadingPt;
    double HLTLeadingPt;
    double HFsumET;
    
    Tree->SetBranchAddress("Ntrkoffline",&Ntrkoffline);
    Tree->SetBranchAddress("NtrkPixel",&NtrkPixel);
    Tree->SetBranchAddress("NtrkFull",&NtrkFull);
    Tree->SetBranchAddress("OfflineVtxZ",&OfflineVtxZ);
    Tree->SetBranchAddress("HLTVtxZ",&HLTVtxZ);
    Tree->SetBranchAddress("OfflineLeadingPt",&OfflineLeadingPt);
    Tree->SetBranchAddress("HLTLeadingPt",&HLTLeadingPt);
    Tree->SetBranchAddress("HFsumET",&HFsumET);
    
    //loop for histograms
    int nentry = Tree->GetEntries();
    for(int it=0;it<nentry;it++)
    {
        Tree->GetEntry(it);
        
        if(it%5000==0) cout<<"processing "<<it<<"th event"<<endl;
        
        if(fabs(OfflineVtxZ)>15) continue;
        
        //fill Ntrkoffline distribution without any cut for HM+HighPt
        if(Ntrkoffline>=120) hFull->Fill(OfflineLeadingPt);
        //fill Ntrkoffline distribution with cuts for HM+HighPt
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
        
        //fill Ntrkoffline distribution without any cut for HFsum+HighPt
        if(HFsumET>HFSumThreshold) hFullHFsum->Fill(OfflineLeadingPt);
        //fill Ntrkoffline distribution with cuts for HFsum+HighPt
        if(HFsumET>HFSumThreshold)
        {
            for(int j=0;j<NPt;j++)
            {
                if(fabs(HLTVtxZ)>15) continue;
                if(HLTLeadingPt<PtThreshold[j]) continue;
                hHMThres2[j]->Fill(OfflineLeadingPt);
            }
        }
    }
    
    //make eff turn-ons
    for(int i=0;i<NPt;i++)
    {
        grHMThres1[i]->Divide(hHMThres1[i],hFull);
        grHMThres2[i]->Divide(hHMThres2[i],hFullHFsum);
    }
    
    //plot eff turn-ons
    TH1D* HLTeff = new TH1D("HLTeff","HLTeff",60,0,60);
    HLTeff->GetXaxis()->SetTitle("pT_{offline}^{Leading}");
    HLTeff->GetYaxis()->SetTitle("HLT eff");
    HLTeff->SetTitle("");
   
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
    
    //HM+HighPt eff
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
    
    //HFsum+HighPt eff
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
}

