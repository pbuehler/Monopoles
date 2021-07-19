#include <ALICEStyle.h>

#include "TSystem.h"
#include "TRandom.h"
#include "TF1.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLatex.h"

Double_t ekinmin(Double_t mass, Int_t gch, Double_t eta);
Double_t FPAsigmaSP (Double_t mass, Int_t gch);
void SPGenerate(Int_t nev, Double_t mass);
void doemall();
