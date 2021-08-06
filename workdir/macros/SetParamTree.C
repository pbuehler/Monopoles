#include "stdio.h"
#include "iomanip"
#include "stdlib.h"
#include "TSystem.h" 
#include "TSystemDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TList.h"

//
// usage:
//  aliroot -n -b -l -q  "${MonopolesDir}/macros/SetParamTree.C( 100,1,10,\"/workdir/results\", kTRUE, kTRUE)"
//  aliroot -n -b -l -q  "${MonopolesDir}/macros/SetParamTree.C( 100,2,11,\"/workdir/results\", kTRUE, kTRUE)"

void SetParamTree (Double_t massIn, Int_t chargeIn, Int_t cc, char* dpath, Bool_t withHistos = kFALSE, Bool_t withTrackPoints = kFALSE)
{

  // some parameters
  Double_t radiusLimit = 230.;
  Double_t clusterAmplitudeLimit = 900.;
  
  Double_t mass = massIn;
  Int_t charge = chargeIn;

  Int_t cnt = 0;
  Double_t eta, rapidity, mom, pt, Ekin;
  Double_t rMin, rMax, rtmp;
  Double_t dMin, dMax, dtmp, length;
  Int_t numCluster, criteria, ID;
  TVector3 pMin, pMax;

  // prepare result tree
  TTree *t= new TTree("ParamTree","Parameters of acceptance studies");
  t->Branch("Counter", &cnt, "Counter/I");
  t->Branch("Mass", &mass, "Mass/D");
  t->Branch("Charge", &charge, "Charge/I");
  t->Branch("Rapidity", &rapidity, "Rapidity/D");
  t->Branch("Eta", &eta, "Eta/D");
  t->Branch("Mom", &mom, "Mom/D");
  t->Branch("Pt", &pt, "Pt/D");
  t->Branch("Ekin", &Ekin, "Ekin/D");
  t->Branch("nCluster", &numCluster, "nCluster/I");
  t->Branch("rMin", &rMin, "rMin/D");
  t->Branch("rMax", &rMax, "rMax/D");
  t->Branch("dMin", &dMin, "dMin/D");
  t->Branch("dMax", &dMax, "dMax/D");
  t->Branch("Length", &length, "length/D");
  t->Branch("Criteria", &criteria, "Criteria/I");
  
  // check for subdirectories of ddir
  TString label = Form("%4.4d", cc);
  TString ddir = Form("%s/m%4.4i/g%i/%s",dpath, (Int_t)massIn, chargeIn, label.Data());
  Long_t id,size,flags,mt;
  gSystem->GetPathInfo(ddir.Data(),&id,&size,&flags,&mt);
  
  TList *files = NULL;
  if (size <= 0) {
    printf("%s is empty!\n", ddir.Data());
    return;
  } else {
    TSystemDirectory dir(ddir.Data(), ddir.Data());
    files = dir.GetListOfFiles();
    if (!files) {
      printf("%s is empty!\n", ddir.Data());
      return;
    }
  }

  // open output file
  TString ofname = Form("paramTree_%4.4i_%i_%s.root", (Int_t)massIn, chargeIn, label.Data());
  TFile *ofile = TFile::Open(ofname.Data(),"RECREATE");

  // loop over subdirectories of ddir
  TSystemFile *file;
  TString fname;
  TIter next(files);
  while ((file=(TSystemFile*)next()))
  {
    fname = file->GetName();
    if (!file->IsDirectory()) {
      continue;
    }
    if (strcmp(fname.Data(),".")==0 || strcmp(fname.Data(),"..")==0) {
      continue;
    }
    auto sdir = fname;

    // open input files
    // TPC.RecPoints.root
    fname = ddir+"/"+sdir+"/"+TString("TPC.RecPoints.root");
    TFile *f1 = TFile::Open(fname.Data());
    printf("File %s opened.\n", fname.Data());
    
    if(!f1) {
      printf ("%s could not be opened!\n", fname.Data());
      continue;
    }
    if(f1->IsZombie()){
      printf("Trying to recover %s\n", fname.Data());
      f1->Recover();
      if(f1->TestBit(TFile::kRecovered)) {
        printf("  The file has been recovered!\n");
      } else {
        printf("  The file cold no be recovered!\n");
        continue;
      }
    }
    
    // Kinematics.root
    fname = ddir+"/"+sdir+"/"+TString("Kinematics.root");
    TFile *f2 = TFile::Open(fname.Data());
    if (!f2) {
      printf ("%s could not be opened!\n", fname.Data());
      continue;
    }
    if (f2->IsZombie()){
      printf("Trying to recover %s\n", fname.Data());
      f2->Recover();
      if (f2->TestBit(TFile::kRecovered)) {
        printf("  The file has been recovered!\n");
      } else {
        printf("  The file could no be recovered!\n");
        continue;
      }
    }
    printf("Processing files in %s ...",  (ddir+"/"+sdir).Data());
    
    // loop over events
    TH1F *hnc = NULL;
    TList *trackPoints = NULL;
    for(auto ne=0; ne<f1->GetNkeys(); ne++)
    {
      cnt = 10*ne + sdir.Atoi();
      
      if (withHistos) {
        auto hname = Form("clusterAmplitude%4.4i", cnt);
        hnc = new TH1F(hname, hname, 1200, 0, 1200);
      }
      
      if (withTrackPoints) {
        auto lname = Form("trackPoints%4.4i", cnt);
        trackPoints = new TList();
        trackPoints->SetName(lname);
      }
      
      // read TPC.RecPoints
      auto d1 = (TDirectory*)f1->Get(Form("Event%d",ne));
      auto treeRecP = (TTree*)d1->Get("TreeR");
      auto row = new AliTPCClustersRow();
      treeRecP->GetBranch("Segment")->SetAddress(&row);
      
      AliTPCclusterMI* cluster = NULL;
      
      // reset r and d
      numCluster = 0;
      rMin=999.;
      rMax=0.;
      dMin=999.;
      dMax=0.;
      pMin = TVector3(0., 0., 0.);
      pMax = TVector3(0., 0., 0.);
      for(auto it=0; it<treeRecP->GetEntries(); it++) {
        treeRecP->GetEntry(it);
        for(auto nCl = 0; nCl < row->GetArray()->GetEntriesFast(); nCl++){
	        
          // only clusters with a minimum amplitude are considered
          cluster = static_cast<AliTPCclusterMI*> (row->GetArray()->UncheckedAt(nCl));
          if (withHistos) {
            hnc->Fill(cluster->GetMax());
          }
          printf("cluster amplitude: %f\n", cluster->GetMax());
          if (cluster->GetMax() <= clusterAmplitudeLimit) {
            continue;
          }
          numCluster++;
          
          // save trackPoint
          if (withTrackPoints) {
            auto trackPoint = new TVector3(cluster->GetX(),cluster->GetY(),cluster->GetZ());
            trackPoints->Add((TObject*)trackPoint);
          }
          
	        // compute distance from origin
          rtmp = (cluster->GetX())*(cluster->GetX())+(cluster->GetY())*(cluster->GetY());
	        dtmp = TMath::Sqrt(rtmp+(cluster->GetZ())*(cluster->GetZ()));
	        rtmp = TMath::Sqrt(rtmp);
	        
          // save minimum/maximum values
          if (rtmp > rMax) {
            rMax = rtmp;
          }
          if (rtmp < rMin) {
            rMin = rtmp;
          }
	        if (dtmp > dMax) {
            dMax = dtmp;
            pMax = TVector3(cluster->GetX(), cluster->GetY(), cluster->GetZ());
          }
	        if (dtmp < dMin) {
            dMin = dtmp;
            pMin = TVector3(cluster->GetX(), cluster->GetY(), cluster->GetZ());
          }
        }
      }
      length = (pMax-pMin).Mag();
      
      // track has to reach a radius of radiusLimit
      if (rMax >= radiusLimit) {
        criteria = 1;
      } else {
        criteria = -1;
      }
      
      // read  Kinematics
      TParticle *particle = NULL;
      auto d2 = (TDirectory*)f2->Get(Form("Event%d",ne));
      auto treeKine = (TTree*)d2->Get("TreeK");
      treeKine->SetBranchAddress("Particles",&particle);
      for(auto i=0; i<treeKine->GetEntries(); i++){
        treeKine->GetEntry(i);
        ID = particle->GetPdgCode();
        if (abs(ID) == 60000000) {
          rapidity  = particle->Y();
          eta  = particle->Eta();
          mom = particle->P();
          pt = particle->Pt();
          Ekin = TMath::Sqrt(mass*mass+mom*mom)-mass;
          break;
        }
      }

      // print summary line
      printf("event [%i] nCluster: %i r[min / max]: %f / %f d[min / max]: %f / %f length: %f\n", ne, numCluster, rMin, rMax, dMin, dMax, length);
      
      // fill result tree
      t->Fill();
      
      // save histogram
      if (withHistos) {
        ofile->cd();
        hnc->Write();
        delete hnc;
      }
      
      // save trackPoints
      if (withTrackPoints) {
        ofile->cd();
        trackPoints->Write(trackPoints->GetName(), 1);
        trackPoints->SetOwner(kTRUE);
        delete trackPoints;
      }
    }
  
    // close input files
    f1->cd();
    f1->Close();
    f2->cd();
    f2->Close();
    printf("done!\n");
  
  }
  
  // save result tree
  ofile->cd();
  t->Write();  
  ofile->Close();
  
}
  
