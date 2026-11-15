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
    #include "TLatex.h"

    #include <vector>
    #include <iostream>
    #include <cmath>
    #include "/home/piotrek/projekt_d/klasyfikator/RFclassifier.cxx"

    using namespace std;

    void rozgrzewka()
    {
        //__Wczytywanie plikow__
        TFile *real_data = new TFile("/home/piotrek/inzynierka/pliki_root/data23_2trk_moreTCvars.root");
        TFile *mc_ee = new TFile("/home/piotrek/inzynierka/pliki_root/mc_jpsi_ee_860171.AOD.e8587_e8586_s4441_r15716_sig_sampling.root");
        TFile *mc_mumu = new TFile("/home/piotrek/inzynierka/pliki_root/mc_jpsi_mumu_860173.AOD.e8587_e8586_s4357_s4258_r15716_r15516_sig_sampling.root");
        TFile *yyee = new TFile("/home/piotrek/inzynierka/pliki_root/yyee_starlight_sig_sampling.root");
        TFile *yymumu = new TFile("/home/piotrek/inzynierka/pliki_root/yymumu_starlight_sig_sampling.root");

        //---------------------------------------------------------------------------------------------------
        //-------------------------------------------Dane rzeczywiste----------------------------------------
        //---------------------------------------------------------------------------------------------------
        //__Drzewka__
        TTree *t = (TTree*)real_data->Get("G2TauTree");

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

        //------Inicjowanie histogramów------
        
        TH1F *h_lam_1 = new TH1F("h_lam","Topo cluster lambda for first particle",100,0,4);
        TH1F *h_pt_1 = new TH1F("h_pt", "Topo cluster PT for first particle", 100, 0, 3.5);
        TH1F *h_lam2_1 = new TH1F("h_lam2", "Topo cluster lambda 2 for first particle", 100, 0, 600);
        TH1F *h_lam_2 = new TH1F("h_lam_2","Topo cluster lambda for second particle",100,0,4);
        TH1F *h_pt_2 = new TH1F("h_pt2", "Topo cluster PT for second particle", 100, 0, 3.5);
        TH1F *h_lam2_2 = new TH1F("h_lam2_2", "Topo cluster lambda 2 for second particle", 100, 0, 600);
        
        TH1F *h_classified_f = new TH1F("h_classified_f", "Classified particle (first particle)", 100, -0.001, 1.001);
        TH1F *h_classified_f_2 = new TH1F("h_classified_f_2", "Classified particle (second particle)", 100, -0.001, 1.001);
        TH1F *h_classfied_avg = new TH1F("h_classified_avg", "Classified particle (average)", 100, -0.001, 1.001);

        
        TH1F *h_class_pion_1 = new TH1F("h_class_pion_1", "Classified particle pion (first particle)", 100, -0.001, 1.001);
        TH1F *h_class_pion_2 = new TH1F("h_class_pion_2", "Classified particle pion (second particle)", 100, -0.001, 1.001);
        TH1F *h_class_pion_avg = new TH1F("h_class_pion_avg", "Classified particle pion (average)", 100, -0.001, 1.001);
        
        /*
        TH1F *h_classified_f = new TH1F("h_classified_f", "Classified particle (first particle) MC ee simulation", 100, -0.001, 1.001);
        TH1F *h_classified_f_2 = new TH1F("h_classified_f_2", "Classified particle (second particle) MC ee simulation" , 100, -0.001, 1.001);
        TH1F *h_classfied_avg = new TH1F("h_classified_avg", "Classified particle (average) MC ee simulation", 100, -0.001, 1.001);
*/
        
        TH1F* histograms[12] = {h_lam_1, h_pt_1, h_lam2_1, h_lam_2, h_pt_2, h_lam2_2, h_classified_f, h_classified_f_2, h_classfied_avg, h_class_pion_1, h_class_pion_2, h_class_pion_avg}; 
        int nentries = t->GetEntries();
        
        for (int i=0; i<12; ++i) {
            //------Markery------
            histograms[i]->SetMarkerColor(9);
            histograms[i]->SetMarkerSize(1.1);
            histograms[i]->SetMarkerStyle(20);
            //------Podpisanie osi------
            histograms[i]->GetYaxis()->SetTitle("Entries");
        }

        for (TH1F *histogram : {h_classified_f, h_classified_f_2, h_classfied_avg, h_class_pion_1, h_class_pion_2, h_class_pion_avg}) {
            histogram->SetStats(false);
            histogram->SetFillColor(kAzure - 9);
            histogram->SetLineColor(kAzure + 2);
            histogram->SetLineWidth(2);
            histogram->GetXaxis()->SetTitle("Muon probability");
        }

        //__Klasyfikator__

        RFClassifier clf("/home/piotrek/projekt_d/klasyfikator/RFMuonElectron.onnx");

        int cluster_candidates = 0;
        int cluster_invalid_data = 0;
        int track1_cluster_matches = 0;
        int track2_cluster_matches = 0;
        int muon_count = 0;
        int electron_count = 0;
        long long cf_trigger = 0;
        long long cf_track_data = 0;
        long long cf_track_sizes = 0;
        long long cf_two_tracks = 0;
        long long cf_charge = 0;
        long long cf_track_pt = 0;
        long long cf_pair_pt = 0;
        long long cf_topo_data = 0;
        long long cf_both_clusters = 0;
        long long cf_classifier = 0;

        //__Pętla__
        for(int i=0;i<nentries;i++){
            //__Pobieranie danych__
            t->GetEntry(i);

            //__Podstawowe breaki__
            if (i<0) break;
            if (passed == 0) continue;
            cf_trigger++;
            if (!pt || !eta || !phi || !charge || !d0 || !islooseprimary) continue;
            cf_track_data++;
            if (eta->size() != pt->size() || phi->size() != pt->size() ||
                charge->size() != pt->size() || d0->size() != pt->size() ||
                islooseprimary->size() != pt->size()) continue;
            cf_track_sizes++;

            vector<int> selected_tracks;
            for (int j = 0; j < pt->size(); ++j) {
                if (fabs((*d0)[j]) >= 1) continue;
                if ((*pt)[j] <= 0.1) continue;
                if (fabs((*eta)[j]) >= 2.5) continue;
                if (!(*islooseprimary)[j]) continue;
                selected_tracks.push_back(j);
            }
            if (selected_tracks.size() != 2) continue;
            cf_two_tracks++;
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
            cf_charge++;

            //__Druga czesc cutow__
            if ((*pt)[track_index1] <= 1 || (*pt)[track_index2] <= 1) continue;
            cf_track_pt++;
            if (PT >= 0.2) continue;
            cf_pair_pt++;
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
            cf_topo_data++;

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

            if (foundCluster1) {
                track1_cluster_matches++;
                h_lam_1->Fill(max_lam/(10*10*10));
                h_lam2_1->Fill(max_lam2/(10*10*10));
                h_pt_1->Fill(max_pt);
            }

            if (foundCluster2) {
                track2_cluster_matches++;
                h_lam_2->Fill(max_lam_2/(10*10*10));
                h_lam2_2->Fill(max_lam2_2/(10*10*10));
                h_pt_2->Fill(max_pt_2);
            }

            if (!foundCluster1 || !foundCluster2) continue;
            cf_both_clusters++;

            float pMuon_1 = clf.Predict(max_pt/(*pt)[track_index1], max_emProb, max_lam, max_lam2, max_r2);
            float pMuon_2 = clf.Predict(max_pt_2/(*pt)[track_index2], max_emProb_2, max_lam_2, max_lam2_2, max_r2_2);

            
            h_classified_f->Fill(pMuon_1);
            h_classified_f_2->Fill(pMuon_2);

            float avg_prob = (pMuon_1 + pMuon_2) / 2.0;
            h_classfied_avg->Fill(avg_prob);
            cf_classifier++;

            if (avg_prob > 0.5) {
                muon_count++;
            } else {
                electron_count++;
            }
            
            /*
            h_class_pion_1->Fill(pMuon_1);
            h_class_pion_2->Fill(pMuon_2);

            float avg_prob = (pMuon_1 + pMuon_2) / 2.0;
            h_class_pion_avg->Fill(avg_prob);
            */
            
        }
/*
        h_class_pion_1->GetYaxis()->SetTitle("Entries");
        h_class_pion_2->GetYaxis()->SetTitle("Entries");
        h_class_pion_avg->GetYaxis()->SetTitle("Entries");
*/
        /*
        TCanvas *c1 = new TCanvas("c1","Topo cluster lambda",1200,800);
        c1->Clear();
        h_lam_1->Draw();
        c1->Update();
        c1->SaveAs("TopoLambda.pdf(");
        c1->Clear();
        h_lam_2->Draw();
        c1->Update();
        c1->SaveAs("TopoLambda.pdf)");
                
        TCanvas *c2 = new TCanvas("c2","Topo cluster lambda2",1200,800);
        c2->Clear();
        h_lam2_1->Draw();
        c2->Update();
        c2->SaveAs("TopoLambda2.pdf(");
        c2->Clear();
        h_lam2_2->Draw();
        c2->Update();
        c2->SaveAs("TopoLambda2.pdf)");

        TCanvas *c3 = new TCanvas("c3","Topo cluster PT",1200,800);
        c3->Clear();
        h_pt_1->Draw();
        c3->Update();
        c3->SaveAs("TopoPT.pdf(");
        c3->Clear();
        h_pt_2->Draw();
        c3->Update();
        c3->SaveAs("TopoPT.pdf)");
        */

        /*
        TCanvas *c4 = new TCanvas("c4","Classified particle MC ee simulation",1200,800);
        c4->SetLogy();
        c4->Clear();
        h_classified_f->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classified/ClassifiedParticleMCEE.pdf(");
        c4->Clear();
        h_classified_f_2->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classified/ClassifiedParticleMCEE.pdf");
        c4->Clear();
        h_classfied_avg->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classified/ClassifiedParticleMCEE.pdf)");
        */

        
        TCanvas *c4 = new TCanvas("c4","Classified particle",1200,800);
        c4->SetLogy();
        c4->Clear();
        h_classified_f->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classified/ClassifiedParticle.pdf(");
        c4->Clear();
        h_classified_f_2->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classified/ClassifiedParticle.pdf");
        c4->Clear();
        h_classfied_avg->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classified/ClassifiedParticle.pdf)");
        
/*
        TCanvas *c4 = new TCanvas("c4","Classified particle",1200,800);
        c4->SetLogy();
        c4->Clear();
        h_class_pion_1->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classified/ClassifiedParticlePion.pdf(");
        c4->Clear();
        h_class_pion_2->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classified/ClassifiedParticlePion.pdf");
        c4->Clear();
        h_class_pion_avg->Draw("HIST");
        c4->Update();
        c4->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classified/ClassifiedParticlePion.pdf)");
*/
        cout << "-------------------Muon count: " << muon_count << "-------------------" << endl;
        cout << "-------------------Electron count: " << electron_count << "-------------------" << endl;
        cout << "Cluster matches: track 1 = " << track1_cluster_matches
             << ", track 2 = " << track2_cluster_matches
             << ", invalid cluster entries = " << cluster_invalid_data << endl;

        auto report_cut = [](const char *name, long long before, long long after) {
            cout << name << ": " << after << " passed, "
                 << before - after << " rejected" << endl;
        };
        cout << "Event cutflow (cumulative):" << endl;
        report_cut("All events -> trigger", nentries, cf_trigger);
        report_cut("Track branches present", cf_trigger, cf_track_data);
        report_cut("Track vector sizes match", cf_track_data, cf_track_sizes);
        report_cut("Exactly two selected tracks", cf_track_sizes, cf_two_tracks);
        report_cut("Opposite-sign charge", cf_two_tracks, cf_charge);
        report_cut("Both track pT > 1 GeV", cf_charge, cf_track_pt);
        report_cut("Pair pT < 0.2 GeV", cf_track_pt, cf_pair_pt);
        report_cut("Mass in (2.5, 3.5) GeV", cf_pair_pt, cluster_candidates);
        report_cut("Topology branches valid", cluster_candidates, cf_topo_data);
        report_cut("Clusters matched to both tracks", cf_topo_data, cf_both_clusters);
        report_cut("Classifier histograms filled", cf_both_clusters, cf_classifier);

        TH1F *h_cutflow = new TH1F("h_cutflow_warmup",
                                   ";Selection stage;Events remaining (log scale)",
                                   12, 0.5, 12.5);
        const char *cut_labels[12] = {
            "All events", "Trigger", "Branches valid", "Sizes valid",
            "Exactly 2 tracks", "Opposite sign", "Both pT > 1 GeV",
            "Pair pT < 0.2 GeV", "Mass 2.5-3.5 GeV",
            "Topology valid", "Both clusters", "Classifier"
        };
        const long long cut_counts[12] = {
            nentries, cf_trigger, cf_track_data, cf_track_sizes,
            cf_two_tracks, cf_charge, cf_track_pt, cf_pair_pt,
            cluster_candidates, cf_topo_data, cf_both_clusters, cf_classifier
        };
        for (int bin = 1; bin <= 12; ++bin) {
            h_cutflow->GetXaxis()->SetBinLabel(bin, cut_labels[bin - 1]);
            h_cutflow->SetBinContent(bin, cut_counts[bin - 1]);
        }
        h_cutflow->SetStats(false);
        h_cutflow->SetFillColor(kAzure - 9);
        h_cutflow->SetLineColor(kAzure + 2);
        h_cutflow->SetLineWidth(2);
        h_cutflow->SetBarWidth(0.72);
        h_cutflow->SetBarOffset(0.14);
        const double cutflow_min = cut_counts[11] > 0 ? cut_counts[11] * 0.5 : 1;
        h_cutflow->SetMinimum(cutflow_min);
        h_cutflow->SetMaximum(nentries * 2.0);
        h_cutflow->GetXaxis()->SetLabelFont(42);
        h_cutflow->GetYaxis()->SetLabelFont(42);
        h_cutflow->GetXaxis()->SetTitleFont(42);
        h_cutflow->GetYaxis()->SetTitleFont(42);
        h_cutflow->GetXaxis()->SetLabelSize(0.03);
        h_cutflow->GetYaxis()->SetLabelSize(0.035);
        h_cutflow->GetXaxis()->SetTitleSize(0.04);
        h_cutflow->GetYaxis()->SetTitleSize(0.04);
        h_cutflow->GetXaxis()->SetTitleOffset(1.7);
        h_cutflow->GetYaxis()->SetTitleOffset(1.35);

        TCanvas *c5 = new TCanvas("c5", "Warmup event selection cutflow", 1800, 1000);
        c5->SetLeftMargin(0.12);
        c5->SetRightMargin(0.04);
        c5->SetBottomMargin(0.28);
        c5->SetTopMargin(0.18);
        c5->SetLogy();
        c5->SetGridy();
        c5->SetTicks(1, 1);
        h_cutflow->Draw("BAR");

        TLatex cutflow_header;
        cutflow_header.SetNDC();
        cutflow_header.SetTextFont(42);
        cutflow_header.SetTextAlign(22);
        cutflow_header.SetTextSize(0.032);
        cutflow_header.DrawLatex(0.54, 0.96, "Warmup event selection cutflow");
        cutflow_header.SetTextSize(0.021);
        cutflow_header.DrawLatex(0.54, 0.925,
            "Track preselection: |d_{0}| < 1 mm, p_{T} > 0.1 GeV, |#eta| < 2.5, LoosePrimary");

        cutflow_header.SetNDC(false);
        cutflow_header.SetTextSize(0.018);
        for (int bin = 1; bin <= 12; ++bin) {
            const double label_y = cut_counts[bin - 1] > 0
                ? cut_counts[bin - 1] * 1.25
                : cutflow_min * 1.25;
            cutflow_header.DrawLatex(bin, label_y, Form("%lld", cut_counts[bin - 1]));
        }
        c5->Update();
        c5->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classfied/Cutflow.pdf");
        c5->SaveAs("/home/piotrek/projekt_d/rozgrzewka/classfied/Cutflow.png");

    }