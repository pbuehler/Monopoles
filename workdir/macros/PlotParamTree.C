#include "PlotParamTree.h"

// usage:
//  aliroot -n -b -l -q  "PlotParamTree.C(\"/data0/physics/projects/IMBA/data/Monopoles/workdir/results/m0100/g1/0010/paramTree_0100_1_0010.root\", 230.)"
//
//  root> PlotParamTree("/data0/physics/projects/IMBA/data/Monopoles/paramTree_0100_1_0010.root", 230.)
//  root> PlotParamTree("/data0/physics/projects/IMBA/data/Monopoles/paramTree_0100_2_0011.root", 230.)
//  root> PlotParamTree("/data0/physics/projects/IMBA/data/Monopoles/workdir/results/m0100/g1/0010/paramTree_0100_1_0010.root", 230.)
//  root> PlotParamTree("/data0/physics/projects/IMBA/data/Monopoles/workdir/results/m0100/g3/0017/paramTree_0100_3_0017.root", 230.)

// ----------------------------------------------------------------------------
void PlotParamTree (char* fname, Double_t rLim)
{

  Int_t nClusterLim = 90;
  
  // get mass, gch, and label
  TString sfn = Form("%s", fname);
  auto toks = sfn.Tokenize("/");
  auto fnbase = ((TObjString*)toks->At(toks->GetEntries()-1))->GetString();
  toks = fnbase.Tokenize("_");
  Double_t mass = ((TObjString*)toks->At(1))->GetString().Atof();
  Int_t gch = ((TObjString*)toks->At(2))->GetString().Atoi();
  auto rest = ((TObjString*)toks->At(toks->GetEntries()-1))->GetString();
  toks = rest.Tokenize(".");
  TString label =  ((TObjString*)toks->At(0))->GetString();
  
  // get access to data
  TFile *f1 = TFile::Open(fname, "READ");
  if (!f1) {
    return;
  }

  TTree *tr = (TTree*)f1->Get("ParamTree");
  if (!tr) {
    printf("Can not retrieve ParamTree in %s\n", fname);
    return;
  }

  Int_t cnt;
  Double_t rapidity, eta, ekin, pt;
  Int_t numCluster;
  Double_t rMin, rMax;
  Double_t dMin, dMax;
  Double_t length;
  tr->SetBranchAddress("Counter", &cnt);
  tr->SetBranchAddress("Rapidity", &rapidity);
  tr->SetBranchAddress("Eta", &eta);
  tr->SetBranchAddress("Pt", &pt);
  tr->SetBranchAddress("Ekin", &ekin);
  tr->SetBranchAddress("nCluster", &numCluster);
  tr->SetBranchAddress("rMin", &rMin);
  tr->SetBranchAddress("rMax", &rMax);
  tr->SetBranchAddress("dMin", &dMin);
  tr->SetBranchAddress("dMax", &dMax);
  tr->SetBranchAddress("Length", &length);

  // prepare histograms
  TH1F *h0a = new TH1F("radius", "radius", 60., 0., 300.);
  h0a->GetYaxis()->SetTitle("relative");
  h0a->GetYaxis()->SetRangeUser(0.0, 1.0);
  h0a->GetXaxis()->SetTitle("r_{max} [cm]");
  h0a->SetLineColor(4);

  TH2F *h0b = new TH2F("r_vs_n", "r_vs_n", 60, 0., 300., 100, 0., 500.);
  h0b->GetXaxis()->SetTitle("r_{max} [cm]");
  h0b->GetYaxis()->SetTitle("n_{Cluster}");
  h0b->SetLineColor(4);

  TH2F *h0c = new TH2F("l_vs_n","l_vs_n", 60, 0., 300., 100, 0., 500.);
  h0c->GetXaxis()->SetTitle("length [cm]");
  h0c->GetYaxis()->SetTitle("n_{Cluster}");
  h0c->SetLineColor(4);

  TH2F *h0d = new TH2F("l_vs_r","l_vs_r", 60, 0., 300., 60, 0., 300.);
  h0d->GetXaxis()->SetTitle("length [cm]");
  h0d->GetYaxis()->SetTitle("r_{max} [cm]");
  h0d->SetLineColor(4);

  TH2F *h1 = new TH2F("ekin_vs_rapidity_g","ekin_vs_rapidity_g",50, 0., 500., 150, -1.5, 1.5);
  h1->GetXaxis()->SetTitle("E_{kin} [GeV/c^{2}]");
  h1->GetYaxis()->SetTitle("eta [1]");
  h1->GetZaxis()->SetTitle("relative");
  TH2F *h2 = new TH2F("ekin_vs_rapidity_b","ekin_vs_rapidity_b",50, 0., 500., 150, -1.5, 1.5);
  
  // distribution of cluster amplitudes
  TH1F *clusterAmplitudes = NULL;

  // tracks
  TVector3* trP = NULL;
  TGraph2D *gr = NULL;
  TH2D* hgr = NULL;

  // Now loop over all events
  Long64_t nev = tr->GetEntries();
  printf("Number of events: %lli\n", nev);
  
  // prepare plots
  TString ttxt;
  TObject* obj = NULL;
  ALICEStyle(0);
  Double_t small = 1.E-5;
  fnbase = Form("ALICEAcceptance_%4.4i_%2.2i_%s.ps", (Int_t) mass, gch, label.Data());
  Bool_t islog = kTRUE;
  
  TLatex latex;
  TCanvas *c1 = new TCanvas();
  c1->SaveAs(Form("%s(", fnbase.Data()));
  
  Double_t etamean = 0.0;
  Double_t d2u;
  Double_t weight = 1./nev;
  for (auto ii=0; ii<nev; ii++) {
    tr->GetEntry(ii);
    
    // update histogram
    h0a->Fill(rMax, weight);
    
    d2u = 0.0;
    if (dMax > 0 && dMin < 999.) {
      d2u = dMax-dMin;
    }
    h0b->Fill(rMax, numCluster, 1.0);
    h0c->Fill(length, numCluster, 1.0);
    h0d->Fill(length, rMax, 1.0);
    if ((rMax > rLim) || (numCluster > nClusterLim)) {
      h1->Fill(ekin, rapidity, 1.0);
    } else {
      h2->Fill(ekin, rapidity, 1.0);
    }
    
    // update etamean
    etamean += eta;
    
    // update clusterAmplitude histogram
    label = Form("clusterAmplitude%4.4i", cnt);
    obj = f1->Get(label.Data());
    if (obj) {
      if (!clusterAmplitudes) {
        clusterAmplitudes = (TH1F*)obj->Clone("clusterAmplitudes");
      } else {
        clusterAmplitudes->Add((TH1F*)obj);
      }
    }
    
    // make track plot
    // get trackPoints
    label = Form("trackPoints%4.4i", cnt);
    obj = f1->Get(label.Data());
    if (obj) {
      auto trackPoints = (TList*)obj;

      // fill TGraph
      c1->Clear();
      c1->Divide(1,1);

      // fill TGraph
      auto nPoints = trackPoints->GetEntries();
      if (nPoints > 0) {
        gr = new TGraph2D(nPoints);
        hgr = new TH2D("track", "track", 500, -250., 250., 500, -250., 250.);
        gr->SetHistogram(hgr);
        for (auto jj=0; jj<nPoints; jj++) {
          trP = (TVector3*)trackPoints->At(jj);
          gr->SetPoint(jj, trP->Z(), trP->Y(), trP->X());
        }
        
        gr->SetMarkerStyle(8);
        gr->Draw("P");
        gPad->Update();
        gr->GetXaxis()->SetTitle("z [cm]");
        gr->GetXaxis()->SetLabelSize(0.03);
        gr->GetXaxis()->SetLabelOffset(5.0E-3);
        gr->GetXaxis()->SetTitleSize(0.04);
        gr->GetXaxis()->SetTitleOffset(1.5);

        gr->GetYaxis()->SetTitle("y [cm]");
        gr->GetYaxis()->SetLabelSize(0.03);
        gr->GetYaxis()->SetLabelOffset(2.0E-3);
        gr->GetYaxis()->SetTitleSize(0.04);
        gr->GetYaxis()->SetTitleOffset(1.5);

        gr->GetZaxis()->SetTitle("x [cm]");
        gr->GetZaxis()->SetLabelSize(0.03);
        gr->GetZaxis()->SetLabelOffset(5.0E-3);
        gr->GetZaxis()->SetTitleSize(0.04);
        gr->GetZaxis()->SetTitleOffset(1.2);
        gr->GetZaxis()->SetLimits(-250., 250.);
      }

      latex.SetTextSize(0.03);
      ttxt = Form("| mass = %6.1f [GeV/c^{2}] | gch = %2i | E_{kin} = %6.1f [GeV] | p_{t} = %6.1f [GeV/c] |",  mass, gch, ekin, pt);
      latex.DrawLatexNDC(0.1,0.96,ttxt.Data());
      ttxt = Form("|Ok = %i | r_{max} = %6.1f | n_{cluster} = %i |", ((rMax > rLim) || (numCluster > nClusterLim) ? 1 : 0), rMax, numCluster);
      latex.DrawLatexNDC(0.1,0.92,ttxt.Data());

      c1->Update();
      c1->SaveAs(Form("%s", fnbase.Data()));
      
      // clean up
      if (nPoints > 0) {
        delete gr;
        delete hgr;
      }
    }
  }
  etamean /= nev;
  
  // compute efficiency
  TH2F *h3 = (TH2F*)h1->Clone("total");
  h3->Add(h2, 1.);
  h1->Divide(h3);  
  
  // plot histogram
  c1->Clear();
  c1->Divide(1,1);
  
  h1->Draw("COLZ");
  latex.SetTextAlign(22);
  ttxt = Form("|mass = %6.1f [GeV/c^{2}] | gch = %2i |", mass, gch);
  latex.DrawLatexNDC(0.5,0.85,ttxt.Data());

  c1->Update();
  c1->SaveAs(Form("%s(", fnbase.Data()));
  
  c1->Clear();
  c1->Divide(2, 2, small, small);
  
  c1->cd(1);
  h0a->Draw("E");
  ttxt = Form("|mass = %6.1f [GeV/c^{2}] | gch = %2i | eta ~ %4.1f |", mass, gch, abs(etamean));
  latex.SetTextSize(0.05);
  latex.DrawLatexNDC(0.5,0.85,ttxt.Data());

  c1->cd(2);
  h0b->Draw("BOX");

  c1->cd(3);
  h0c->Draw("BOX");

  c1->cd(4);
  h0d->Draw("BOX");

  c1->Update();
  c1->SaveAs(Form("%s", fnbase.Data()));

  c1->Clear();
  c1->Divide(2, 2, small, small);
  
  c1->cd(1);
  h1->Draw("COLZ");
  latex.SetTextAlign(22);
  ttxt = Form("|mass = %6.1f [GeV/c^{2}] | gch = %2i |", mass, gch);
  latex.DrawLatexNDC(0.5,0.85,ttxt.Data());

  c1->cd(3);
  h0a->Draw("E");

  if (clusterAmplitudes) {
    c1->cd(4);
    c1->cd(4)->SetLogy(kTRUE);
    clusterAmplitudes->GetXaxis()->SetTitle("cluster amplitude");
    clusterAmplitudes->Draw("E");
  }

  // close ps file
  c1->Update();
  c1->SaveAs(Form("%s)", fnbase.Data()));

  delete c1;
  delete h0a;
  delete h0b;
  delete h1;
  delete h2;
  delete h3;
  
}
  
// ----------------------------------------------------------------------------
void doemall ()
{  
  
  Double_t rLim = 230;
  TString ddir = Form("/data0/physics/projects/IMBA/data/Monopoles/workdir/results/");
  
  PlotParamTree(Form("%s/m0005/g1/0053/paramTree_0005_1_0053.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0010/g1/0052/paramTree_0010_1_0052.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0020/g1/0051/paramTree_0020_1_0051.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0050/g1/0050/paramTree_0050_1_0050.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0010/paramTree_0100_1_0010.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0250/g1/0021/paramTree_0250_1_0021.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m1000/g1/0022/paramTree_1000_1_0022.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m1500/g1/0023/paramTree_1500_1_0023.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m2000/g1/0024/paramTree_2000_1_0024.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m3000/g1/0025/paramTree_3000_1_0025.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m4000/g1/0026/paramTree_4000_1_0026.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m5000/g1/0027/paramTree_5000_1_0027.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m6000/g1/0028/paramTree_6000_1_0028.root", ddir.Data()), rLim);

  PlotParamTree(Form("%s/m0100/g1/0031/paramTree_0100_1_0031.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0032/paramTree_0100_1_0032.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0033/paramTree_0100_1_0033.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0034/paramTree_0100_1_0034.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0035/paramTree_0100_1_0035.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0036/paramTree_0100_1_0036.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0037/paramTree_0100_1_0037.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0038/paramTree_0100_1_0038.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0039/paramTree_0100_1_0039.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0040/paramTree_0100_1_0040.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0041/paramTree_0100_1_0041.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g1/0042/paramTree_0100_1_0042.root", ddir.Data()), rLim);

  PlotParamTree(Form("%s/m0100/g2/0011/paramTree_0100_2_0011.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g3/0012/paramTree_0100_3_0012.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g4/0013/paramTree_0100_4_0013.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g5/0014/paramTree_0100_5_0014.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g6/0015/paramTree_0100_6_0015.root", ddir.Data()), rLim);
  PlotParamTree(Form("%s/m0100/g7/0016/paramTree_0100_7_0016.root", ddir.Data()), rLim);

}

// ----------------------------------------------------------------------------
