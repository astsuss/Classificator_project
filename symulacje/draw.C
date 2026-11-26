#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPaletteAxis.h"
#include "TEfficiency.h"
#include "TStyle.h"

#include <algorithm>
#include <cstdio>
#include <initializer_list>
#include <iostream>

void draw(const char *input_path = "testpion_output.root",
            const char *output_path = "/home/piotrek/projekt_d/symulacje/ploty/Rvalues_lowerthreshold.pdf",
            const char *output_dr_path = "/home/piotrek/projekt_d/symulacje/ploty/Rvalues_thresholds.pdf",
            const char *output_class_05 = "/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticle_TestPion05.pdf",
            const char *output_class_01 = "/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticle_TestPion01.pdf",
            const char *output_class_all = "/home/piotrek/projekt_d/symulacje/ploty/ClassifiedParticle_TestPionALL.pdf")
{
    TFile input(input_path, "READ");
    if (input.IsZombie()) {
        std::cerr << "Cannot open ROOT input file: " << input_path << std::endl;
        return;
    }

    TH1F *R_values = nullptr;
    TH1F *treshholds = nullptr;
    TH1F *class_pion = nullptr;
    TH1F *class_pion_01 = nullptr;
    TH1F *class_pion_05 = nullptr;
    input.GetObject("R_values", R_values);
    input.GetObject("matched_after_cuts", treshholds);
    input.GetObject("class_pion_all", class_pion);
    input.GetObject("class_pion_01", class_pion_01);
    input.GetObject("class_pion_05", class_pion_05);
    if (!R_values || !class_pion || !class_pion_01 || !class_pion_05) {
        std::cerr << "Missing required histograms in "
                  << input_path << std::endl;
        return;
    }

    gStyle->SetOptStat(0);

    R_values->SetTitle("dR values for truth pion variables");
    R_values->GetXaxis()->SetTitle("dR value");
    R_values->GetYaxis()->SetTitle("Tracks");
    R_values->GetXaxis()->CenterTitle(true);
    R_values->GetYaxis()->CenterTitle(true);
        
    treshholds->SetTitle("dR cuts thresholds");
    treshholds->GetXaxis()->SetTitle("dR cuts");
    treshholds->GetYaxis()->SetTitle("Number of matched tracks");
    treshholds->GetXaxis()->CenterTitle(true);
    treshholds->GetYaxis()->CenterTitle(true);

    class_pion->SetTitle("Classified particle Pion dR < all");
    class_pion->GetXaxis()->SetTitle("Muon probability");
    class_pion->GetYaxis()->SetTitle("Events");
    class_pion->GetXaxis()->CenterTitle(true);
    class_pion->GetYaxis()->CenterTitle(true);

    class_pion_01->SetTitle("Classified particle Pion dR < 0.1");
    class_pion_01->GetXaxis()->SetTitle("Muon probability");
    class_pion_01->GetYaxis()->SetTitle("Events");
    class_pion_01->GetXaxis()->CenterTitle(true);
    class_pion_01->GetYaxis()->CenterTitle(true);

    class_pion_05->SetTitle("Classified particle Pion dR < 0.05");
    class_pion_05->GetXaxis()->SetTitle("Muon probability");
    class_pion_05->GetYaxis()->SetTitle("Events");
    class_pion_05->GetXaxis()->CenterTitle(true);
    class_pion_05->GetYaxis()->CenterTitle(true);

    for (TH1F *histogram : {R_values, treshholds, class_pion, class_pion_01, class_pion_05}) {
        histogram->SetStats(false);
        histogram->GetXaxis()->SetLabelFont(42);
        histogram->GetYaxis()->SetLabelFont(42);
        histogram->GetXaxis()->SetTitleFont(42);
        histogram->GetYaxis()->SetTitleFont(42);
        histogram->GetXaxis()->SetLabelSize(0.03);
        histogram->GetYaxis()->SetLabelSize(0.03);
        histogram->GetXaxis()->SetTitleSize(0.035);
        histogram->GetYaxis()->SetTitleSize(0.035);
        histogram->GetXaxis()->SetTitleOffset(1.1);
        histogram->GetYaxis()->SetTitleOffset(1.45);
    }

    R_values->SetLineColor(kAzure + 2);
    R_values->SetMarkerColor(kAzure + 2);
    R_values->SetFillColorAlpha(kAzure + 1, 0.30);
    R_values->SetFillStyle(1001);
    R_values->SetLineWidth(2);
    R_values->SetMarkerStyle(20);
    R_values->SetMarkerSize(0.65);

    treshholds->SetLineColor(kAzure + 2);
    treshholds->SetMarkerColor(kAzure + 2);
    treshholds->SetFillColorAlpha(kAzure + 1, 0.30);
    treshholds->SetFillStyle(1001);
    treshholds->SetLineWidth(2);
    treshholds->SetMarkerStyle(20);
    treshholds->SetMarkerSize(0.65);
    
    class_pion->SetLineColor(kAzure + 2);
    class_pion->SetMarkerColor(kAzure + 2);
    class_pion->SetFillColorAlpha(kAzure + 1, 0.30);
    class_pion->SetFillStyle(1001);
    class_pion->SetLineWidth(2);
    class_pion->SetMarkerStyle(20);
    class_pion->SetMarkerSize(0.65);

    for (TH1F *histogram : {class_pion_01, class_pion_05}) {
        histogram->SetLineColor(kAzure + 2);
        histogram->SetMarkerColor(kAzure + 2);
        histogram->SetFillColorAlpha(kAzure + 1, 0.30);
        histogram->SetFillStyle(1001);
        histogram->SetLineWidth(2);
        histogram->SetMarkerStyle(20);
        histogram->SetMarkerSize(0.65);
    }

    TCanvas canvas("R_values", "R values for truth pions", 1800, 1000);
    canvas.SetLeftMargin(0.12);
    canvas.SetBottomMargin(0.14);
    canvas.SetTopMargin(0.09);
    canvas.SetRightMargin(0.12);
    canvas.SetTicks(1, 1);
    canvas.SetLogy();
    R_values->Draw("HIST");

    TCanvas canvas3("Class_pion", "Classified particle pion dR < all", 1800, 1000);
    canvas3.SetLeftMargin(0.12);
    canvas3.SetBottomMargin(0.14);
    canvas3.SetTopMargin(0.09);
    canvas3.SetRightMargin(0.12);
    canvas3.SetTicks(1, 1);
    canvas3.SetLogy();
    class_pion->Draw("HIST");

    canvas3.Update();
    canvas3.SaveAs(output_class_all);

    TCanvas canvas4("Class_pion_01", "Classified particle pion dR < 0.1", 1800, 1000);
    canvas4.SetLeftMargin(0.12);
    canvas4.SetBottomMargin(0.14);
    canvas4.SetTopMargin(0.09);
    canvas4.SetRightMargin(0.12);
    canvas4.SetTicks(1, 1);
    canvas4.SetLogy();
    class_pion_01->Draw("HIST");
    canvas4.Update();
    canvas4.SaveAs(output_class_01);

    TCanvas canvas5("Class_pion_05", "Classified particle pion dR < 0.05", 1800, 1000);
    canvas5.SetLeftMargin(0.12);
    canvas5.SetBottomMargin(0.14);
    canvas5.SetTopMargin(0.09);
    canvas5.SetRightMargin(0.12);
    canvas5.SetTicks(1, 1);
    canvas5.SetLogy();
    class_pion_05->Draw("HIST");
    canvas5.Update();
    canvas5.SaveAs(output_class_05);
/*
    TLine l02(0.2, 20, 0.2, 3e5);
    TLine l015(0.15, 20, 0.15, 3e5);
    TLine l01(0.1, 20, 0.1, 3e5);
    TLine l005(0.05, 20, 0.05, 3e5);

    l02.Draw();
    l015.Draw();
    l01.Draw();
    l005.Draw();
*/
    TLatex cut_label;
    cut_label.SetNDC();
    cut_label.SetTextFont(42);
    cut_label.SetTextSize(0.029);
    cut_label.DrawLatex(0.65, 0.81,
                         "truth pion pt > 1 [GeV]");
    cut_label.DrawLatex(0.65, 0.77,
                         "|truth pion pdg| == 211");

    canvas.Update();
    canvas.SaveAs(output_path);

    TCanvas canvas2("matched_after_cuts", "dR cut thresholds", 1800, 1000);
    canvas2.SetLeftMargin(0.12);
    canvas2.SetBottomMargin(0.14);
    canvas2.SetTopMargin(0.09);
    canvas2.SetRightMargin(0.12);
    canvas2.SetTicks(1, 1);
    canvas2.SetGridy();
    double min_count = treshholds->GetBinContent(treshholds->GetMinimumBin());
    double max_count = treshholds->GetMaximum();
    double padding = 0.1 * (max_count - min_count);
    if (padding <= 0) padding = 0.05 * std::max(1.0, max_count);
    treshholds->SetMinimum(std::max(0.0, min_count - padding));
    treshholds->SetMaximum(max_count + padding);
    treshholds->Draw("HIST");
/*
    TLine l02(0.2, 20, 0.2, 3e5);
    TLine l015(0.15, 20, 0.15, 3e5);
    TLine l01(0.1, 20, 0.1, 3e5);
    TLine l005(0.05, 20, 0.05, 3e5);

    l02.Draw();
    l015.Draw();
    l01.Draw();
    l005.Draw();
*/
    TLatex cut_label2;
    cut_label2.SetNDC();
    cut_label2.SetTextFont(42);
    cut_label2.SetTextSize(0.029);
    cut_label2.DrawLatex(0.65, 0.81,
                         "truth pion pt > 1 [GeV]");
    cut_label2.DrawLatex(0.65, 0.77,
                         "|truth pion pdg| == 211");

    canvas2.Update();
    canvas2.SaveAs(output_dr_path);
}
