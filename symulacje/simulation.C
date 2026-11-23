#if defined(__CLING__)
R__LOAD_LIBRARY(/home/piotrek/onnxruntime-linux-x64-1.31.1/lib/libonnxruntime.so)
#endif
    
    #include "TFile.h"
    #include "TTree.h"
    #include "TF1.h"
    #include "TCanvas.h"
    #include "TFitResultPtr.h"
    #include "TMath.h"

    #include <algorithm>
    #include <random>
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
        vector<float>* d0 = nullptr;
        Bool_t passed;
        vector<bool>* islooseprimary = nullptr;

        t->SetBranchAddress("track_eta", &eta);
        t->SetBranchAddress("track_pt", &pt);
        t->SetBranchAddress("track_phi", &phi);
        t->SetBranchAddress("track_d0", &d0);
        t->SetBranchAddress("passed_HLT_mb_excl_1trk5_pt1_hi_FgapAC5_L1TRT_VTE20", &passed);
        t->SetBranchAddress("track_is_LoosePrimary", &islooseprimary);

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
        TH1F *h_class_pion_pt_gt1 = new TH1F("h_class_pion_pt_gt1", "Muon probability for tracks in entries with truth pion pT > 1 GeV", 100, 0, 1);
        TH1F *h_class_pion_pair_avg = new TH1F("h_class_pion_pair_avg", "Average muon probability for random pairs of tracks", 100, 0, 1);
        vector<float> pion_track_probabilities;

        RFClassifier clf("/home/piotrek/projekt_d/klasyfikator/RFMuonElectron.onnx");

        int nentries = t->GetEntries();

        h_class_pion_pt_gt1->SetStats(false);
        h_class_pion_pt_gt1->SetFillColor(kAzure - 9);
        h_class_pion_pt_gt1->SetLineColor(kAzure + 2);
        h_class_pion_pt_gt1->SetLineWidth(2);
        h_class_pion_pt_gt1->GetXaxis()->SetTitle("Muon probability");
        h_class_pion_pt_gt1->GetYaxis()->SetTitle("Entries");

        h_class_pion_pair_avg->SetStats(false);
        h_class_pion_pair_avg->SetFillColor(kAzure - 9);
        h_class_pion_pair_avg->SetLineColor(kAzure + 2);
        h_class_pion_pair_avg->SetLineWidth(2);
        h_class_pion_pair_avg->GetXaxis()->SetTitle("Average muon probability");
        h_class_pion_pair_avg->GetYaxis()->SetTitle("Entries");

        
        int cluster_candidates = 0;
        int cluster_invalid_data = 0;
        int track1_cluster_matches = 0;
        int track2_cluster_matches = 0;
        int entries_with_pion_above1GeV = 0;

        int t_count = 0;
        int b_count = 0;

        //__Pętla__
        for(int i=0;i<nentries;i++){
            //__Pobieranie danych__
            t->GetEntry(i);

            //__Sprawdzanie pionów__
            if (!pion_pt || pion_pt->empty()) continue;
            //if ((*pion_pt)[0] >= 1.6 || (*pion_pt)[1] >= 1.6) continue;
            bool hasPionAbove1GeV = false;
            for (float pionPt : *pion_pt) {
                if (pionPt > 1.0f) {
                    hasPionAbove1GeV = true;
                    break;
                }
            }
            if (!hasPionAbove1GeV) continue;
            entries_with_pion_above1GeV++;

            //__Podstawowe breaki__
            if (i<0) break;
            if (passed == 0) continue;
            if (!pt || !eta || !phi || !d0 || !islooseprimary) continue;
            if (eta->size() != pt->size() || phi->size() != pt->size() ||
                d0->size() != pt->size() ||
                islooseprimary->size() != pt->size()) continue;

            vector<int> selected_tracks;
            for (int j = 0; j < pt->size(); ++j) {
                if (fabs((*d0)[j]) >= 1) continue;
                if ((*pt)[j] <= 1.0) continue;
                if (fabs((*eta)[j]) >= 2.5) continue;
                if (!(*islooseprimary)[j]) continue;
                selected_tracks.push_back(j);
            }
            if (selected_tracks.empty()) continue;

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

            vector<float> track_probabilities(selected_tracks.size(), 0.0f);
            vector<bool> track_has_cluster(selected_tracks.size(), false);
            bool foundAnyCluster = false;
            for (size_t track_position = 0; track_position < selected_tracks.size(); ++track_position) {
                int track_index = selected_tracks[track_position];
                double max_pt = 0, max_lam = 0, max_lam2 = 0, max_r2 = 0, max_emProb = 0;
                bool foundCluster = false;

                //__Cluster dla pojedynczego tracka__
                for (int j = 0; j < topo_pt->size(); ++j) {
                    if ((*pass_sig)[j] == 0) continue;
                    if (fabs((*topo_eta)[j]) >= 4.9) continue;
                    if ((*topo_pt)[j] <= 0.1) continue;

                    double dphi = TMath::ACos(cos((*phi)[track_index] - (*topo_phi)[j]));
                    double deta = (*eta)[track_index] - (*topo_eta)[j];
                    double dr = sqrt(dphi * dphi + deta * deta);
                    if (dr < 0.8 && (!foundCluster || (*topo_pt)[j] > max_pt)) {
                        foundCluster = true;
                        max_pt = (*topo_pt)[j];
                        max_lam = (*topo_lambda)[j];
                        max_lam2 = (*topo_lambda2)[j];
                        max_r2 = (*topo_r2)[j];
                        max_emProb = (*topo_emProb)[j];
                    }
                }

                if (!foundCluster) continue;

                track_probabilities[track_position] =
                    clf.Predict(max_pt / (*pt)[track_index], max_emProb, max_lam, max_lam2, max_r2);
                track_has_cluster[track_position] = true;
                foundAnyCluster = true;
            }
            if (!foundAnyCluster) continue;

            for (size_t track_position = 0; track_position < track_probabilities.size(); ++track_position) {
                if (track_has_cluster[track_position]) {
                    h_class_pion_pt_gt1->Fill(track_probabilities[track_position]);
                    pion_track_probabilities.push_back(track_probabilities[track_position]);
                }
            }
        }

        random_device random_seed;
        mt19937 random_generator(random_seed());
        shuffle(pion_track_probabilities.begin(), pion_track_probabilities.end(), random_generator);
        for (size_t i = 0; i + 1 < pion_track_probabilities.size(); i += 2) {
            double pair_average =
                (pion_track_probabilities[i] + pion_track_probabilities[i + 1]) / 2.0;
            h_class_pion_pair_avg->Fill(pair_average);
        }

        //__Rysowanie__
        TCanvas *c4 = new TCanvas("c4","Classified particle",1200,800);
        c4->SetLogy();
        h_class_pion_pt_gt1->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticlePionStarlight.pdf");

        TCanvas *c5 = new TCanvas("c5", "Average classified particle", 1200, 800);
        c5->SetLogy();
        h_class_pion_pair_avg->Draw("HIST");
        c5->Update();
        c5->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticlePionStarlightAverage.pdf");

        cout << "---------Entries with pion pT > 1 GeV: "
             << entries_with_pion_above1GeV << "---------" << endl;
        cout << "---------Random pairs averaged: "
             << h_class_pion_pair_avg->GetEntries() << "---------" << endl;


    }