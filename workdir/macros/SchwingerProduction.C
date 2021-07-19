#include "SchwingerProduction.h"

// ----------------------------------------------------------------------------
Double_t ekinmin(Double_t mass, Int_t gch, Double_t eta)
{
    
  if (abs(eta) <= 0.4) {
    return 30.0;
  } else if (abs(eta) <= 0.7)  {
    return 50.;
  } else if (abs(eta) <= 0.8)  {
    return 80.;
  } else {
    return 1.E6;
  }
  
}

// ----------------------------------------------------------------------------
// according to Acharya et al., arXiv:2106.11933v1
//              equation (2)
Double_t FPAsigmaSP (Double_t mass, Int_t gch)
{
  Double_t w = 73.0;        // GeV
  Double_t ccentral = 0.0;  // 1
  Double_t B = 7.6;         // GeV^2
  Double_t RPb = 6.62;      // fm
  Double_t gD = 68.5*gch;   // 1/(2*alpha) = 68.5
  Double_t cf = 2.E7/8.99;
  
  Double_t nSigma = (w/mass + ccentral);
  nSigma *= 2.*TMath::Power(gD*B*RPb/mass,4);
  nSigma *= TMath::Exp(-4.*mass/w);
  nSigma /= 9.*TMath::Power(TMath::Pi()*w,2);   // nb
  nSigma *= cf;
    
  printf("mass: %f, gch: %i, nSigma: %e\n", mass, gch, nSigma);
  
  return nSigma;
  
}

// ----------------------------------------------------------------------------
void SPGenerate(Int_t nev, Double_t mass)
{

  // some parameters
  Int_t gch = 1;
  Double_t lumi = 0.235;  // 1/nb
  Double_t sigmaSP = FPAsigmaSP(mass, gch);

  // random generator
  auto rnd = new TRandom();

  // Schwinger-pair momentum distribution
  // w is a parameter of the momentum distribution
  // in Gould et al., arXiv:2103.14454   : w ~ 80
  // in Acharya et a., arXiv:2106.11933v1: w = 73 +- 3
  Double_t w = 73.0;

  // according to equ. (31) in 	Gould et al., arXiv:2103.14454
  Double_t pmax;
  if (mass < 1000) {
    pmax = 500.;
  } else {
    pmax = mass/2.;
  }
  auto f1 = new TF1("SPP","exp(-4/[0]*(sqrt([1]*[1]+x*x)-[1]))",0., pmax);
  f1->SetParName(0,"w");
  f1->SetParName(0,"mass");
  f1->SetParameter(0,w);
  f1->SetParameter(1,mass);

  // prepare histograms
  TString ttxt = Form("mass = %6.1f [GeV/c^{2}]", mass);
  TH1F h1("ekin","ekin",60, 0.0, 300.);
  h1.GetXaxis()->SetTitle("E_{kin} [GeV/c^{2}]");
  h1.GetYaxis()->SetTitle("relative");
  h1.GetYaxis()->SetRangeUser(1.E-6, 1.E0);
  h1.SetLineColor(4);

  TH1F h2("eta","eta",50, -5.0, 5.0);
  h2.GetXaxis()->SetTitle("eta [1]");
  h2.GetYaxis()->SetTitle("relative");
  h2.GetYaxis()->SetRangeUser(1.E-6, 1.E0);
  h2.SetLineColor(4);

  TH1F h3("pt","pt",100, 0.0, pmax);
  h3.GetXaxis()->SetTitle("p_{t} [GeV/c]");
  h3.GetYaxis()->SetTitle("relative");
  h3.GetYaxis()->SetRangeUser(1.E-6, 1.E0);
  h3.SetLineColor(4);

  TH1F h4("beta","beta",50, 0.0, 1.0);
  h4.GetYaxis()->SetTitle("relative");
  h4.GetXaxis()->SetTitle("beta [1]");
  h4.GetYaxis()->SetRangeUser(1.E-6, 1.E0);
  h4.SetLineColor(4);
  
  TH1F h5("p_o_m","p_o_m",50, 0.0, 1.0);
  h5.GetYaxis()->SetTitle("relative");
  h5.GetXaxis()->SetTitle("p/m");
  h5.SetLineColor(4);
  
  // efficiency map
  TH2F *heff1 = new TH2F("SPefficieny_g", "SPefficieny_g", 30, 0., 300., 50., -1.5, 1.5);
  heff1->GetXaxis()->SetTitle("E_{kin} [GeV/c^{2}]");
  heff1->GetYaxis()->SetTitle("eta [1]");
  heff1->GetZaxis()->SetTitle("relative");
  TH2F *heff2 = new TH2F("SPefficieny_b", "SPefficieny_b", 30, 0., 300., 50., -1.5, 1.5);
  
  // generate nev monopoles
  Double_t p, px, py, pz, ekin;
  Double_t eff = 0.0;
  
  Double_t weight = 1./nev;

  for (auto ii=0; ii<nev; ii++) {
    p = f1->GetRandom();
    rnd->Sphere(px, py, pz, p);
    
    // update histograms
    auto ene = sqrt(mass*mass+p*p);
    auto ekin = ene-mass;
    auto lv = TLorentzVector(px, py, pz, ene);
    h1.Fill(ekin,weight);
    h2.Fill(lv.Eta(),weight);
    h3.Fill(lv.Perp(), weight);
    h4.Fill(lv.Beta(), weight);
    h5.Fill(lv.P()/mass, 1.);
    
    if (ekin > ekinmin(mass, gch, lv.Eta())) {
      heff1->Fill(ekin, lv.Eta(), 1.0);
      eff += 1.0;
    } else {
      heff2->Fill(ekin, lv.Eta(), 1.0);
    }
  }
  eff /= nev;
  
  // normalize h5
  auto norm = 1./h5.GetBinContent(h5.GetMaximumBin());
  printf("Norm %f\n", norm);
  h5.Scale(norm);
  
  Double_t nevExpected = lumi*sigmaSP*eff;

  // compute efficiency
  TH2F *heff3 = (TH2F*)heff1->Clone("total");
  heff3->Add(heff2, 1.);
  heff1->Divide(heff3);

  // plot distributions
  ALICEStyle(0);
  Bool_t islog = kTRUE;
  
  Float_t small = 5e-3;
  auto pad = new TPad();
  TLatex latex;
  latex.SetTextSize(0.05);
  TCanvas *c1 = new TCanvas();
  c1->Divide(2,1,small,small);
  
  pad = (TPad*) c1->cd(1);
  pad->SetLogy(islog);
  h1.Draw("E");
  latex.SetTextAlign(22);
  latex.DrawLatexNDC(0.5,0.85,"Schwinger-Pair production");

  pad = (TPad*) c1->cd(2);
  pad->SetLogy(islog);
  h2.Draw("E");

  c1->Update();
  TString fname = Form("SchwingerProduction_%4.4i.ps(", (Int_t) mass);
  c1->SaveAs(fname.Data());
  
  c1->Clear();
  pad = (TPad*) c1->cd(1);
  pad->SetLogy(islog);
  latex.SetTextAlign(12);
  latex.DrawLatexNDC(0.2,0.25,ttxt.Data());
  h3.Draw("E");

  pad = (TPad*) c1->cd(2);
  pad->SetLogy(islog);
  latex.DrawLatexNDC(0.2,0.25,ttxt.Data());
  h4.Draw("E");
  
  c1->Update();
  fname = Form("SchwingerProduction_%4.4i.ps", (Int_t) mass);
  c1->SaveAs(fname.Data());
  
  c1->Clear();
  c1->Divide(1,1);
  
  pad = (TPad*) c1->cd(2);
  pad->SetLogy(kFALSE);
  h5.GetYaxis()->SetRangeUser(0.0, 1.0);
  h5.Draw("E");
  c1->SaveAs(fname.Data());
  
  c1->Clear();
  
  heff1->Draw("COLZ");
  latex.SetTextAlign(12);
  latex.DrawLatexNDC(0.2,0.85, "Schwinger-Pair production");
  latex.DrawLatexNDC(0.2,0.80, Form("|mass = %6.1f [GeV/c^{2}] | gch = %2i |", mass, gch));
  latex.DrawLatexNDC(0.2,0.75, Form("efficiency = %6.4f", eff));
  latex.DrawLatexNDC(0.2,0.70, Form("L = %4.2f [nb^{-1}], sigma_{SP} = %10.2e [nb]", lumi, sigmaSP));
  latex.DrawLatexNDC(0.2,0.65, Form("Expected number = %10.2e", nevExpected));
  
  c1->Update();
  fname = Form("SchwingerProduction_%4.4i.ps)", (Int_t) mass);
  c1->SaveAs(fname.Data());
  
  // clean up
  delete c1;

}

// ----------------------------------------------------------------------------
void genemall ()
{
  Int_t nev = 10000000;
  
  SPGenerate(nev,    5.);
  SPGenerate(nev,   10.);
  SPGenerate(nev,   20.);
  SPGenerate(nev,   50.);
  SPGenerate(nev,  100.);
  SPGenerate(nev,  250.);
  SPGenerate(nev,  500.);
  SPGenerate(nev, 1000.);
  SPGenerate(nev, 1500.);
  SPGenerate(nev, 2000.);
  SPGenerate(nev, 3000.);
  SPGenerate(nev, 4000.);
  SPGenerate(nev, 5000.);
  SPGenerate(nev, 6000.);
  
  TString cmd = Form("~/tools/.meps2pdf SchwingerProduction_*.eps");
  gSystem->Exec(cmd.Data());
  cmd = Form("pdfunite SchwingerProduction_*.pdf SchwingerProduction.pdf");
  gSystem->Exec(cmd.Data());
  
}

// ----------------------------------------------------------------------------
