#include "headers.h"
#include "Utilities/Ntuple/VertexCompositeTree.h"

const double mMuon = 0.10566;

const std::vector<std::string>& generateFilenames(const std::string& prefix, int num = 1000) {
    static std::vector<std::string> filenames;

    filenames.clear();

    for (int i = 1; i <= num; ++i) {
        std::stringstream ss;
        ss << prefix << "000" << (i / 1000) << "/" << "dimuana_mc_" << i << ".root";
        filenames.push_back(ss.str());
    }

    return filenames;
}

void plotMuons()
{
  const std::string DirName = "/eos/cms/store/group/phys_heavyions/huangxi/PhaseII/SkimAndTree/STARlight_CohJpsi2MuMu_PbPb5TeV_SkimAndTree_v2/CRAB_UserFiles/STARlight_CohJpsi2MuMu_PbPb5TeV_SkimAndTree_v2/240619_142810/";
  const std::vector<std::string> &FileNames = generateFilenames(DirName);
  const auto &treeDir = "dimucontana_mc"; // For MC use dimucontana_mc

  TH2F *hMC_mu = new TH2F("hMC_mu", "MC Muons;#eta;p_{T}/GeV", 50, -3.2, 3.2, 100, 0, 1.8);
  TH2F *hRC_mu = new TH2F("hRC_mu", "Reconstructed Muons;#eta;p_{T}/GeV", 50, -3.2, 3.2, 100, 0, 1.8);
  TH2F *hEff_mu = new TH2F("hEff_mu", "Single Muon Efficiency;#eta;p_{T}/GeV", 50, -3.2, 3.2, 100, 0, 1.8);
  hMC_mu->SetStats(kFALSE);
  hEff_mu->SetStats(kFALSE);

  // Extract the tree
  VertexCompositeTree tree;
  if (!tree.GetTree(FileNames, treeDir)) { std::cout << "Invalid tree!" << std::endl; return; }

  // Loop over the tree
  for(Long64_t jentry=0; jentry<tree.GetEntries(); jentry++){// Get the entry
    if (tree.GetEntry(jentry)<0) { std::cout << "Invalid entry!" << std::endl; return; }

    Double_t posPt_gen  = tree.chargeD1_gen()[0] > 0 ? tree.pTD1_gen()[0] : tree.pTD2_gen()[0];
    Double_t posEta_gen = tree.chargeD1_gen()[0] > 0 ? tree.EtaD1_gen()[0] : tree.EtaD2_gen()[0];
    Double_t posPhi_gen = tree.chargeD1_gen()[0] > 0 ? tree.PhiD1_gen()[0] : tree.PhiD2_gen()[0]; 
    Double_t negPt_gen  = tree.chargeD1_gen()[0] < 0 ? tree.pTD1_gen()[0] : tree.pTD2_gen()[0]; 
    Double_t negEta_gen = tree.chargeD1_gen()[0] < 0 ? tree.EtaD1_gen()[0] : tree.EtaD2_gen()[0]; 
    Double_t negPhi_gen = tree.chargeD1_gen()[0] < 0 ? tree.PhiD1_gen()[0] : tree.PhiD2_gen()[0];

    hMC_mu->Fill(posEta_gen, posPt_gen);
    hMC_mu->Fill(negEta_gen, negPt_gen);

    TLorentzVector posFourMom_gen, negFourMom_gen;
    posFourMom_gen.SetPtEtaPhiM(posPt_gen, posEta_gen, posPhi_gen, mMuon);
    negFourMom_gen.SetPtEtaPhiM(negPt_gen, negEta_gen, negPhi_gen, mMuon);

    if(tree.candSize_mu() > 2) continue;
    if(tree.candSize_mu() == 2 && tree.charge_mu()[0] == tree.charge_mu()[1]) continue;
    
    for (uint imu = 0; imu < tree.candSize_mu(); imu++)
    {
      if(tree.charge_mu()[imu] > 0) hRC_mu->Fill(posEta_gen, posPt_gen);
      else hRC_mu->Fill(negEta_gen, negPt_gen);
    }
  }

  hEff_mu->Divide(hRC_mu, hMC_mu);

  TCanvas* c1 = new TCanvas("c1","",0,0,800,600);
  hMC_mu->Draw("colz");
  c1->SaveAs("MCMuon.png");
  TCanvas *c2 = new TCanvas("c2", "", 0, 0, 800, 600);
  hRC_mu->Draw("colz");
  c2->SaveAs("RCMuon.png");
  TCanvas *c3 = new TCanvas("c3", "", 0, 0, 800, 600);
  hEff_mu->Draw("colz");
  c3->SaveAs("SingleMuonEff.png");
}
