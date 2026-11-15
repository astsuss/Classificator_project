#if defined(__CLING__)
R__LOAD_LIBRARY(/home/piotrek/onnxruntime-linux-x64-1.31.1/lib/libonnxruntime.so)
#endif
    
    #include "TFile.h"
    #include "TTree.h"
    #include "TLorentzVector.h"
    #include "TF1.h"
    #include "TCanvas.h"
    #include "TFitResultPtr.h"
    #include "TMath.h"

    #include <vector>
    #include <iostream>
    #include <cmath>
    #include "/home/piotrek/projekt_d/klasyfikator/RFclassifier.cxx"

    using namespace std;

    void simulation()
    {
        TFile *starlight = new TFile("/home/piotrek/inzynierka/pliki_root/mc23_5p36TeV.860181.StarlightEG_gammaA2JpsiallCoherent_CTFilter.merge.AOD.e8587_e8586_s4357_s4258_r15716_r15516_500k.root");
        TTree *t = (TTree*)starlight->Get("G2TauTree");

        vector<float>* pion_pt = nullptr;
        t->SetBranchAddress("truth_pion_pt", &pion_pt);

        //__Zmienne track__
        vector<float>* pt = nullptr;
        vector<float>* eta = nullptr;
        vector<float>* phi = nullptr;
        vector<int>* charge = nullptr;
        vector<float>* d0 = nullptr;
        Bool_t passed;
        vector<bool>* islooseprimary = nullptr;

        t->SetBranchAddress("track_eta", &eta);
        t->SetBranchAddress("track_pt", &pt);
        t->SetBranchAddress("track_phi", &phi);
        t->SetBranchAddress("track_charge", &charge);
        t->SetBranchAddress("track_d0", &d0);
        t->SetBranchAddress("passed_HLT_mb_excl_1trk5_pt1_hi_FgapAC5_L1TRT_VTE20", &passed);
        t->SetBranchAddress("track_is_LoosePrimary", &islooseprimary);

        Double_t muon_mass=0.1056583745; // GeV/c^2

        TLorentzVector track1, track2, TRACK;

        //__Zmienne topo__
        vector<float>* topo_pt = nullptr;
        vector<float>* topo_eta = nullptr;
        vector<float>* topo_phi = nullptr;
        vector<float>* topo_lambda = nullptr;
        vector<float>* topo_lambda2 = nullptr;
        vector<float>* topo_r2 = nullptr;
        vector<float>* topo_emProb = nullptr;
        vector<bool>* pass_sig = nullptr;

        t->SetBranchAddress("topo_cluster_eta", &topo_eta);
        t->SetBranchAddress("topo_cluster_pt", &topo_pt);
        t->SetBranchAddress("topo_cluster_phi", &topo_phi);
        t->SetBranchAddress("topo_cluster_lambda", &topo_lambda);
        t->SetBranchAddress("topo_cluster_lambda2", &topo_lambda2);
        t->SetBranchAddress("topo_cluster_r2", &topo_r2);
        t->SetBranchAddress("topo_cluster_EM_prob", &topo_emProb);
        t->SetBranchAddress("topo_cluster_pass_sig_cut", &pass_sig);

        //__Hists__
        TH1F *h_class_pion_1 = new TH1F("h_class_pion_1", "Pion classification probability for track 1", 100, 0, 1);
        TH1F *h_class_pion_2 = new TH1F("h_class_pion_2", "Pion classification probability for track 2", 100, 0, 1);
        TH1F *h_class_pion_avg = new TH1F("h_class_pion_avg", "Average pion classification probability", 100, 0, 1);

        RFClassifier clf("/home/piotrek/projekt_d/klasyfikator/RFMuonElectron.onnx");

        TH1F* histograms[3] = {h_class_pion_1, h_class_pion_2, h_class_pion_avg}; 
        int nentries = t->GetEntries();
        
        for (int i=0; i<3; ++i) {
            //------Markery------
            histograms[i]->SetMarkerColor(9);
            histograms[i]->SetMarkerSize(1.1);
            histograms[i]->SetMarkerStyle(20);
            //------Podpisanie osi------
            histograms[i]->GetYaxis()->SetTitle("Entries");
        }

        for (TH1F *histogram : {h_class_pion_1, h_class_pion_2, h_class_pion_avg}) {
            histogram->SetStats(false);
            histogram->SetFillColor(kAzure - 9);
            histogram->SetLineColor(kAzure + 2);
            histogram->SetLineWidth(2);
            histogram->GetXaxis()->SetTitle("Muon probability");
        }

        
        int cluster_candidates = 0;
        int cluster_invalid_data = 0;
        int track1_cluster_matches = 0;
        int track2_cluster_matches = 0;
        int count = 0;

        //__Pętla__
        for(int i=0;i<nentries;i++){
            //__Pobieranie danych__
            t->GetEntry(i);

            //__Sprawdzanie pionów__
            if (!pion_pt || pion_pt->size() < 2) continue;
            if ((*pion_pt)[0] <= 0.1 || (*pion_pt)[1] <= 0.1) continue;
            if ((*pion_pt)[0] >= 1.6 || (*pion_pt)[1] >= 1.6) continue;
            count ++;

            //__Podstawowe breaki__
            if (i<0) break;
            if (passed == 0) continue;
            if (!pt || !eta || !phi || !charge || !d0 || !islooseprimary) continue;
            if (eta->size() != pt->size() || phi->size() != pt->size() ||
                charge->size() != pt->size() || d0->size() != pt->size() ||
                islooseprimary->size() != pt->size()) continue;

            vector<int> selected_tracks;
            for (int j = 0; j < pt->size(); ++j) {
                if (fabs((*d0)[j]) >= 1) continue;
                if ((*pt)[j] <= 0.1) continue;
                if (fabs((*eta)[j]) >= 2.5) continue;
                if (!(*islooseprimary)[j]) continue;
                selected_tracks.push_back(j);
            }
            if (selected_tracks.size() != 2) continue;
            int track_index1 = selected_tracks[0];
            int track_index2 = selected_tracks[1];

            //__Ustawianie Lorentz Vectorwo__
            track1.SetPtEtaPhiM((*pt)[track_index1], (*eta)[track_index1], (*phi)[track_index1], muon_mass);
            track2.SetPtEtaPhiM((*pt)[track_index2], (*eta)[track_index2], (*phi)[track_index2], muon_mass);

            TRACK = track1 + track2;

            //__Masa niezmiennicza__
            double inv = TRACK.M();
            double PT = TRACK.Pt();
            double ETA = fabs(TRACK.Eta());
            double PHI = TRACK.Phi();
            //__Warunki początkowe__
            if ((*charge)[track_index1] * (*charge)[track_index2] >= 0 ) continue;

            //__Druga czesc cutow__
            if ((*pt)[track_index1] <= 1 || (*pt)[track_index2] <= 1) continue;
            if (PT >= 0.2) continue;
            if (inv <= 2.5 || inv >= 3.5) continue;

            cluster_candidates++;
            if (!topo_pt || !topo_eta || !topo_phi || !topo_lambda ||
                !topo_lambda2 || !topo_r2 || !topo_emProb || !pass_sig ||
                topo_eta->size() != topo_pt->size() ||
                topo_phi->size() != topo_pt->size() ||
                topo_lambda->size() != topo_pt->size() ||
                topo_lambda2->size() != topo_pt->size() ||
                topo_r2->size() != topo_pt->size() ||
                topo_emProb->size() != topo_pt->size() ||
                pass_sig->size() != topo_pt->size()) {
                cluster_invalid_data++;
                continue;
            }

            double max_pt = 0, max_lam = 0, max_lam2 = 0, max_r2 = 0, max_emProb = 0;
            double max_pt_2 = 0, max_lam_2 = 0, max_lam2_2 = 0, max_r2_2 = 0, max_emProb_2 = 0;
            bool foundCluster1 = false;
            bool foundCluster2 = false;

            //__Cluster__
            for (int j=0;j<topo_pt->size();j++){
                if ((*pass_sig)[j] == 0) continue;
                if (fabs((*topo_eta)[j]) >= 4.9) continue;
                if ((*topo_pt)[j] <= 0.1) continue;

                float dr1 = sqrt(TMath::ACos(cos((*phi)[track_index1] - (*topo_phi)[j]))*TMath::ACos(cos((*phi)[track_index1] - (*topo_phi)[j])) + fabs((*eta)[track_index1] - (*topo_eta)[j])*fabs((*eta)[track_index1] - (*topo_eta)[j]));
                float dr2 = sqrt(TMath::ACos(cos((*phi)[track_index2] - (*topo_phi)[j]))*TMath::ACos(cos((*phi)[track_index2] - (*topo_phi)[j])) + fabs((*eta)[track_index2] - (*topo_eta)[j])*fabs((*eta)[track_index2] - (*topo_eta)[j]));

                if (dr1 < 0.8 && (!foundCluster1 || (*topo_pt)[j] > max_pt)) {
                    foundCluster1 = true;
                    max_pt = (*topo_pt)[j];
                    max_lam = (*topo_lambda)[j];
                    max_lam2 = (*topo_lambda2)[j];
                    max_r2 = (*topo_r2)[j];
                    max_emProb = (*topo_emProb)[j];
                }

                if (dr2 < 0.8 && (!foundCluster2 || (*topo_pt)[j] > max_pt_2)) {
                    foundCluster2 = true;
                    max_pt_2 = (*topo_pt)[j];
                    max_lam_2 = (*topo_lambda)[j];
                    max_lam2_2 = (*topo_lambda2)[j];
                    max_r2_2 = (*topo_r2)[j];
                    max_emProb_2 = (*topo_emProb)[j];
                }
            }

            if (!foundCluster1 || !foundCluster2) continue;

            float pMuon_1 = clf.Predict(max_pt/(*pt)[track_index1], max_emProb, max_lam, max_lam2, max_r2);
            float pMuon_2 = clf.Predict(max_pt_2/(*pt)[track_index2], max_emProb_2, max_lam_2, max_lam2_2, max_r2_2);
            
            h_class_pion_1->Fill(pMuon_1);
            h_class_pion_2->Fill(pMuon_2);

            float avg_prob = (pMuon_1 + pMuon_2) / 2.0;
            h_class_pion_avg->Fill(avg_prob);
        }

        //__Rysowanie__
        TCanvas *c4 = new TCanvas("c4","Classified particle",1200,800);
        c4->SetLogy();
        c4->Clear();
        h_class_pion_1->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticlePionStarlight.pdf(");
        c4->Clear();
        h_class_pion_2->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticlePionStarlight.pdf");
        c4->Clear();
        h_class_pion_avg->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticlePionStarlight.pdf)");

        cout << "-------------------Pion count: " << count << "-------------------" << endl;

    }