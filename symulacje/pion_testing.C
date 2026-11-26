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
    #include <limits>
    #include <random>
    #include <vector>
    #include <iostream>
    #include <cmath>
    #include "/home/piotrek/projekt_d/klasyfikator/RFclassifier.cxx"

    using namespace std;

    void pion_testing()
    {
        //__Pliki__
        TFile *starlight = new TFile("/home/piotrek/inzynierka/pliki_root/mc23_5p36TeV.860181.StarlightEG_gammaA2JpsiallCoherent_CTFilter.merge.AOD.e8587_e8586_s4357_s4258_r15716_r15516_500k.root");
        TTree *t = (TTree*)starlight->Get("G2TauTree");

        //__Zapisywanie__
        const string& out_file = "testpion_output.root";
        TFile* output = new TFile(out_file.data(), "recreate");

        //__Zmienne truth__
        vector<float>* pion_pt = nullptr;
        vector<float>* pion_eta = nullptr;
        vector<float>* pion_phi = nullptr;
        vector<float>* pion_pdg = nullptr;
        t->SetBranchAddress("truth_pion_pt", &pion_pt);
        t->SetBranchAddress("truth_pion_eta", &pion_eta);
        t->SetBranchAddress("truth_pion_phi", &pion_phi);
        t->SetBranchAddress("truth_pion_pdg", &pion_pdg);

        //__Zmienne track__
        vector<float>* pt = nullptr;
        vector<float>* eta = nullptr;
        vector<float>* phi = nullptr;
        vector<float>* d0 = nullptr;
        vector<float>* charge = nullptr;
        Bool_t passed;
        vector<bool>* islooseprimary = nullptr;

        t->SetBranchAddress("track_eta", &eta);
        t->SetBranchAddress("track_pt", &pt);
        t->SetBranchAddress("track_phi", &phi);
        t->SetBranchAddress("track_d0", &d0);
        t->SetBranchAddress("passed_HLT_mb_excl_1trk5_pt1_hi_FgapAC5_L1TRT_VTE20", &passed);
        t->SetBranchAddress("track_is_LoosePrimary", &islooseprimary);
        t->SetBranchAddress("track_charge", &charge);

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

        const double cuts[] = {0.02, 0.05, 0.075, 0.10, 0.15, 0.20};

        //__Hists__
        TH1F *R_values = new TH1F("R_values", "R values for pion candidates", 50, 0, 0.2);
        TH1F *matched_after_cuts = new TH1F("matched_after_cuts", "Matched tracks/pions after dR cuts", 6, 0.5, 6.5);
        const char* cut_labels[] = {"0.020", "0.050", "0.075", "0.100", "0.150", "0.200"};
        for (int c = 0; c < 6; c++) {
            matched_after_cuts->GetXaxis()->SetBinLabel(c + 1, cut_labels[c]);
        }
        TH1F *pt_05 = new TH1F("pt_05", "p_T track / p_T truth for dR < 0.05", 50, 0, 5);
        TH1F *pt_1 = new TH1F("pt_1", "p_T track / p_T truth for dR < 0.1", 50, 0, 5);
        TH1F *class_pion_01 = new TH1F("class_pion_01", "Classified particle pion dR < 0.1", 50, 0, 1);
        TH1F *class_pion_05 = new TH1F("class_pion_05", "Classified particle pion dR < 0.05", 50, 0, 1);
        TH1F *class_pion_all = new TH1F("class_pion_all", "Classified particle pion dR < all", 50, 0, 1);


        RFClassifier clf("/home/piotrek/projekt_d/klasyfikator/RFMuonElectron.onnx");

        int nentries = t->GetEntries();
        
        //__Cutflow__
        int all_tracks = 0;
        int tracks_after_d0 = 0;
        int tracks_after_pt = 0;
        int trakcs_after_eta = 0;
        int tracks_after_primary = 0;

        int all_pions = 0;
        int pions_after_pt = 0;
        int pions_after_pdg = 0;

        long long tracks_with_pion_candidate = 0;
        long long tracks_without_pion_candidate = 0;

        int charge_ok_05 = 0;
        int charge_nok_05 = 0;
        int charge_ok_1 = 0;
        int charge_nok_1 = 0;

        long long duplicate_pions_005 = 0;
        long long extra_track_assignments_005 = 0;
        long long duplicate_pions_010 = 0;
        long long extra_track_assignments_010 = 0;

        int cluster_candidates = 0;
        int cluster_invalid_data = 0;

        //__Pętla__
        for(int i=0;i<nentries;i++){
            //__Pobieranie danych__
            t->GetEntry(i);

            //__Podstawowe breaki__
            if (i<0) break;
            if (passed == 0) continue;
            if (!pt || !eta || !phi || !d0 || !islooseprimary || !charge) continue;
            if (eta->size() != pt->size() || phi->size() != pt->size() ||
                d0->size() != pt->size() ||
                islooseprimary->size() != pt->size() || pt->size() != charge->size()) continue;
            if (!pion_pt || !pion_eta || !pion_phi || !pion_pdg) continue;
            if (pion_eta->size() != pion_pt->size() ||
                pion_phi->size() != pion_pt->size() ||
                pion_pdg->size() != pion_pt->size()) continue;
            if (eta->size() != pt->size() || phi->size() != pt->size()) continue;

            vector<int> pion_matches_005(pion_pt->size(), 0);
            vector<int> pion_matches_010(pion_pt->size(), 0);
            all_tracks = all_tracks + pt->size();
            vector<double> best_dr_by_track(pt->size(), std::numeric_limits<double>::infinity());

            vector<int> selected_tracks;
            for (int j = 0; j < pt->size(); ++j) {
                if (fabs((*d0)[j]) >= 1) continue;
                tracks_after_d0++;
                if ((*pt)[j] <= 1.0) continue;
                tracks_after_pt++;
                if (fabs((*eta)[j]) >= 2.5) continue;
                trakcs_after_eta++;
                if (!(*islooseprimary)[j]) continue;
                tracks_after_primary++;
                selected_tracks.push_back(j);
            }
            if (selected_tracks.empty()) continue;

            //__Matchowanie trackow i pionow__
            for (size_t track_position = 0; track_position < selected_tracks.size(); ++track_position){
                int track_index = selected_tracks[track_position];
                double best_dr =  numeric_limits<double>::infinity();
                int pion_idx = -1;
                for(int pion = 0;pion<pion_pt->size();pion++){
                    all_pions++;
                    if ((*pion_pt)[pion] <= 1.0) continue;
                    pions_after_pt++;
                    if (abs((*pion_pdg)[pion]) != 211) continue;
                    pions_after_pdg++;

                    float deta = ((*eta)[track_index] - (*pion_eta)[pion]);
                    float dphi_raw = ((*phi)[track_index] - (*pion_phi)[pion]);
                    float dphi = atan2(sin(dphi_raw), cos(dphi_raw));
                    double dr = sqrt(deta*deta + dphi*dphi);

                    if (dr < best_dr){
                        best_dr = dr;
                        pion_idx = pion;
                    }
                }

                if (!isfinite(best_dr)) {
                    ++tracks_without_pion_candidate;
                    continue;
                }

                ++tracks_with_pion_candidate;
                R_values->Fill(best_dr);

                for(int c = 0; c < 6; c++){
                    if (best_dr < cuts[c]) matched_after_cuts->Fill(c + 1);
                }

                if (pion_idx >= 0 && best_dr < 0.05) {
                    if((*pion_pdg)[pion_idx] * (*charge)[track_index] > 0 ) charge_ok_05 ++;
                    else charge_nok_05 ++;

                    float comp05 = (*pt)[track_index] / (*pion_pt)[pion_idx];
                    pt_05->Fill(comp05);

                    ++pion_matches_005[pion_idx];
                }
                if (pion_idx >= 0 && best_dr < 0.1) {
                    if((*pion_pdg)[pion_idx] * (*charge)[track_index] > 0 ) charge_ok_1 ++;
                    else charge_nok_1 ++;     

                    float comp1 = (*pt)[track_index] / (*pion_pt)[pion_idx];
                    pt_1->Fill(comp1);  

                    ++pion_matches_010[pion_idx];
                }

                best_dr_by_track[track_index] = best_dr;

            }

            for (size_t pion = 0; pion < pion_matches_005.size(); ++pion) {
                if (pion_matches_005[pion] > 1) {
                    ++duplicate_pions_005;
                    extra_track_assignments_005 += pion_matches_005[pion] - 1;
                }
                if (pion_matches_010[pion] > 1) {
                    ++duplicate_pions_010;
                    extra_track_assignments_010 += pion_matches_010[pion] - 1;
                }
            }
            

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

                const float score = clf.Predict(
                    max_pt / (*pt)[track_index],
                    max_emProb,
                    max_lam,
                    max_lam2,
                    max_r2
                );

                const double pion_dr = best_dr_by_track[track_index];

                // All selected tracks with a valid cluster; not a pion-only sample.
                class_pion_all->Fill(score);

                // Truth-pion-matched classifier samples.
                if (isfinite(pion_dr) && pion_dr < 0.05) {
                    class_pion_05->Fill(score);
                }
                if (isfinite(pion_dr) && pion_dr < 0.10) {
                    class_pion_01->Fill(score);
                }

            }

        }

        //__Rysowanie__
        TCanvas *c_r = new TCanvas("c_r","R values for each track",1200,800);
        R_values->Draw("HIST");
        c_r->Update();
        c_r->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/Rvalues_lowerthreshold.pdf");
        c_r->Clear();
        matched_after_cuts->Draw("HIST");
        c_r->Update();
        c_r->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/Rvalues_thresholds.pdf");

        TCanvas *c_pt = new TCanvas("c_pt", "Differences in p_T", 1200, 800);
        pt_05->Draw("Hist");
        c_pt->Update();
        c_pt->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/pTcomp.pdf(");
        c_pt->Clear();
        pt_1->Draw("HIST");
        c_pt->Update();
        c_pt->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/pTcomp.pdf)"); 

        TCanvas *c_class = new TCanvas("c_class", "Classified particle pion", 1200, 800);
        class_pion_all->Draw("Hist");
        c_class->Update();
        c_class->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticle_TestPion.pdf(");
        c_class->Clear();
        class_pion_05->Draw("HIST");
        c_class->Update();
        c_class->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticle_TestPion.pdf"); 
        c_class->Clear();
        class_pion_01->Draw("HIST");
        c_class->Update();
        c_class->SaveAs("/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticle_TestPion.pdf)"); 

        //__Dane__
        cout << "All tracks: " << all_tracks << endl;
        cout << "Tracks after d0: " << tracks_after_d0 << endl;
        cout << "Tracks after pt: " << tracks_after_pt << endl;
        cout << "Tracks after eta: " << trakcs_after_eta << endl;
        cout << "Tracks after primary: " << tracks_after_primary<< endl;
        cout << "Pion count: " << all_pions << endl;
        cout << "Pion count after pt: " << pions_after_pt << endl;
        cout << "Pion count after pdg: " << pions_after_pdg << endl;
        cout << "Tracks with pion candidate: " << tracks_with_pion_candidate << endl;
        cout << "Tracks without pion candidate: " << tracks_without_pion_candidate << endl;
        cout << "dR cut 0.05 matching charges: " << charge_ok_05 << endl;
        cout << "dR cut 0.05 not matching charges: " << charge_nok_05 << endl;
        cout << "dR cut 0.1 matching charges: " << charge_ok_1 << endl;
        cout << "dR cut 0.1 not matching charges: " << charge_nok_1 << endl;
        cout << "Duplicate pions dR < 0.05: " << duplicate_pions_005 << endl;
        cout << "Duplicate pions dR < 0.1: " << duplicate_pions_010 << endl;
        cout << "dR < 0.05 entries: " << class_pion_05->GetEntries() << endl;
        cout << "dR < 0.1 entries: " << class_pion_01->GetEntries() << endl;
        cout << "dR < 0.all entries: " << class_pion_all->GetEntries() << endl;
        

    output->cd();
    output->Write();
    output->Close();


    }