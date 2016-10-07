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

void HLTeff_HighMultiplicity()
{
    TFile* file0 = TFile::Open("HLTTree_newEPOS_MB.root");
    
    //Trigger threshold
    int PixelThreshold = 40;
    int TowerCountThreshold[3] = {0,62,74};
    int MultiplicityThreshold[5] = {120,150,185,220,250};
    int NTC = 3;
    int NMult = 5;
    //Triggre combination, index corresponds to above array
    int indexMult[5] = {0,1,2,3,4};
    int indexTowerCount[5] = {0,0,1,2,2};
    int Nindex = 5;
    
    //initialize eff plots
    int Nplot = 2+5*2;
    TH1D* h = new TH1D("","",80,0,400);
    TH1D* hFull = h.Clone();
    TH1D* hMult[100];
    TH1D* hTowerCount[100];
    TH1D* hAll[100];
    for(int i=0;i<Nplot;i++)
    {
        hMult[i] = (TH1D*)h.Clone();
        hTowerCount[i] = (TH1D*)h.Clone();
        hAll[i] = (TH1D*)h.Clone();
    }
    
    TGraphAsymmErrors* grMult[100];
    TGraphAsymmErrors* grTowerCount[100];
    TGraphAsymmErrors* grAll[100];
    for(int i=0;i<Nplot;i++)
    {
        grMult[i] = new TGraphAsymmErrors();
        grTowerCount[i] = new TGraphAsymmErrors();
        grAll[i] = new TGraphAsymmErrors();
    }
    
    //Get Tree branch
    TNtuple* Tree = (TNtuple*)file0->Get("HLTinfo/HLTeff");
    
    int Ntrkoffline;
    int NtrkFull;
    int NtrkPixel;
    double TowerCount;
    double OfflineVtxZ;
    double HLTVtxZ;
    
    Tree->SetBranchAddress("Ntrkoffline",&Ntrkoffline);
    Tree->SetBranchAddress("NtrkPixel",&NtrkPixel);
    Tree->SetBranchAddress("NtrkFull",&NtrkFull);
    Tree->SetBranchAddress("TowerCount",&TowerCount);
    Tree->SetBranchAddress("OfflineVtxZ",&OfflineVtxZ);
    Tree->SetBranchAddress("HLTVtxZ",&HLTVtxZ);
    
    //loop for histograms
    int nentry = Tree->GetEntries();
    for(int it=0;it<nentry;it++)
    {
        Tree->GetEntry(it);
        
        if(it%5000==0) cout<<"processing "<<it<<"th event"<<endl;
        
        if(fabs(OfflineVtxZ)>15) continue;
        
        //fill Ntrkoffline distribution without any cut
        hFull->Fill(Ntrkoffline);
        //fill Ntrkoffline distribution with cuts
        for(int i=0;i<NTC;i++)
        {
            if(TowerCount<=TowerCountThreshold[i]) continue;
            hTowerCount[i]->Fill(Ntrkoffline);
        }
        for(int j=0;j<NMult;j++)
        {
            if(fabs(HLTVtxZ)>15) continue;
            if(NtrkPixel<PixelThreshold) continue;
            if(NtrkFull<MultiplicityThreshold[j]) continue;
            hMult[j]->Fill(Ntrkoffline);
        }
        for(int k=0;k<Nindex;k++)
        {
            if(fabs(HLTVtxZ)>15) continue;
            if(TowerCount<=TowerCountThreshold[indexTowerCount[k]]) continue;
            if(NtrkPixel<PixelThreshold) continue;
            if(NtrkFull<MultiplicityThreshold[indexMult[k]]) continue;
            hAll[k]->Fill(Ntrkoffline);
        }
    }
    
    //make eff turn-ons
    for(int i=0;i<NTC;i++)
    {
        grTowerCount[i]->Divide(hTowerCount[i],hFull);
    }
    for(int j=0;j<NMult;j++)
    {
        grMult[j]->Divide(hMult[j],hFull);
    }
    for(int k=0;k<Nindex;k++)
    {
        grAll[k]->Divide(hAll[k],hFull);
    }
    
    //plot eff turn-ons
    TH1D* L1eff = new TH1D("L1eff","L1eff",350,0,350);
    L1eff->GetXaxis()->SetTitle("Ntrkoffline");
	L1eff->GetYaxis()->SetTitle("L1 eff");
    L1eff->SetTitle("");
	
    TH1D* HLTeff = new TH1D("HLTeff","HLTeff",350,0,350);
    HLTeff->GetXaxis()->SetTitle("Ntrkoffline");
    HLTeff->GetYaxis()->SetTitle("HLT eff");
    HLTeff->SetTitle("");

    TH1D* Fulleff = new TH1D("Fulleff","Fulleff",350,0,350);
    Fulleff->GetXaxis()->SetTitle("Ntrkoffline");
    Fulleff->GetYaxis()->SetTitle("L1+HLT eff");
    Fulleff->SetTitle("");
   
    //add lines at where the trigger is expected to be efficient
    TLine* l = new TLine(0,1,300,1);
    l->SetLineStyle(3);
    
    TLine* l1 = new TLine(120,0,120,1.05);
    TLine* l2 = new TLine(150,0,150,1.05);
    TLine* l3 = new TLine(185,0,185,1.05);
    TLine* l4 = new TLine(220,0,220,1.05);
    TLine* l5 = new TLine(250,0,250,1.05);
    l1->SetLineStyle(7);
    l2->SetLineStyle(7);
    l3->SetLineStyle(7);
    l4->SetLineStyle(7);
    l5->SetLineStyle(7);
    
    //create canvas
    TCanvas* c1 = new TCanvas("c1","c1",600,600);
    TCanvas* c2 = new TCanvas("c2","c2",600,600);
    TCanvas* c3 = new TCanvas("c3","c3",600,600);
    
    //L1 eff
    c1->cd();
    L1eff->Draw();
    TLegend* leg = new TLegend(0.2,0.7,0.4,0.9);
    leg->SetFillStyle(0);
    for(int i=0;i<NTC;i++)
    {
        grTowerCount[i]->SetMarkerColor(i+1);
        grTowerCount[i]->SetMarkerSize(0.8);
        grTowerCount[i]->Draw("PSAME");
        leg->AddEntry(grTowerCount[i],Form("TowerCount > %d",TowerCountThreshold[i]),"p");
    }
    leg->Draw();
    l->Draw("LSAME");
    l1->Draw("LSAME");
    l2->Draw("LSAME");
    l3->Draw("LSAME");
    l4->Draw("LSAME");
    l5->Draw("LSAME");
    
    //HLT eff
    c2->cd();
    HLTeff->Draw();
    TLegend* leg1 = new TLegend(0.2,0.7,0.4,0.9);
    leg1->SetFillStyle(0);
    for(int i=0;i<NMult;i++)
    {
        grMult[i]->SetMarkerColor(i+1);
        grMult[i]->SetMarkerSize(0.8);
        grMult[i]->Draw("PSAME");
        leg1->AddEntry(grMult[i],Form("Ntrkonlie #geq %d",MultiplicityThreshold[i]),"p");
    }
    leg1->Draw();
    l->Draw("LSAME");
    l1->Draw("LSAME");
    l2->Draw("LSAME");
    l3->Draw("LSAME");
    l4->Draw("LSAME");
    l5->Draw("LSAME");
    
    //L1+HLT eff
    c3->cd();
    Fulleff->Draw();
    TLegend* leg2 = new TLegend(0.2,0.7,0.4,0.9);
    leg2->SetFillStyle(0);
    for(int i=0;i<NMult;i++)
    {
        grAll[i]->SetMarkerColor(i+1);
        grAll[i]->SetMarkerSize(0.8);
        grAll[i]->Draw("PSAME");
        leg2->AddEntry(grAll[i],Form("HLT_HM%d",MultiplicityThreshold[i]),"p");
    }
    leg2->Draw();
    l->Draw("LSAME");
    l1->Draw("LSAME");
    l2->Draw("LSAME");
    l3->Draw("LSAME");
    l4->Draw("LSAME");
    l5->Draw("LSAME");

}

