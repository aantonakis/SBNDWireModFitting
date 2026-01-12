


void plot_2d_splines(const char* infile, const char* prefix) {
    gROOT->SetBatch(kTRUE);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis); // or kBird, kRainBow, etc.
    gStyle->SetNumberContours(255);

    TFile* f = TFile::Open(infile, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "❌ Could not open file: " << infile << std::endl;
        return;
    }

    const int Nplanes = 3;
    const int NTPCs   = 2;

    for (int i = 0; i < Nplanes * NTPCs; ++i) {
        TString name = Form("splines_%d", i);
        auto g = (TGraph2DErrors*)f->Get(name);
        if (!g) {
            std::cerr << "⚠️  Graph " << name << " not found in file.\n";
            continue;
        }

        // Create 2D histogram from the TGraph2D for visualization
        TH2D* h2 = g->GetHistogram();
        if (!h2) {
            std::cerr << "⚠️  Could not create histogram for " << name << "\n";
            continue;
        }

        TCanvas* c = new TCanvas(Form("c_%d", i), name, 900, 800);
        c->SetRightMargin(0.15);
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);

/*
        h2->SetTitle(Form("2D Fit Ratios - %s;X coordinate;Y coordinate;Data/MC Ratio", name.Data()));
        h2->GetXaxis()->SetTitle(g->GetXaxis()->GetTitle());
        h2->GetYaxis()->SetTitle(g->GetYaxis()->GetTitle());
        h2->Draw("COLZ");
*/
        // Draw the graph points if desired
        //g->SetMarkerStyle(20);
        //g->SetMarkerSize(0.6);
        //g->SetMarkerColor(kBlack);
        //g->Draw("P SAME");
        //g->Draw("PCOL");
        //g->Draw("Colz");


        int n = g->GetN(); // number of points
        double* x = g->GetX();
        double* y = g->GetY();
    
        // Create a 2D scatter plot (TGraph)
        TGraph* g2 = new TGraph(n, x, y);
    
        g2->SetMarkerStyle(6); // dot style
        g2->SetMarkerColor(2); // dot style
        //g2->GetXaxis()->SetTitle("X");
        //g2->GetYaxis()->SetTitle("Y");
        g2->Draw("AP"); // A = axis, P = points

        TLine* ltop = new TLine(-200, 90, 0, 90);
        TLine* lbottom = new TLine(-200, -90, 0, -90);
        TLine* lleft = new TLine(-200, -90, -200, 90);
        TLine* lright = new TLine(200, -90, 200, 90);

        ltop->Draw("Same");
        lbottom->Draw("Same");
        lleft->Draw("Same");
        lright->Draw("Same");

        // Save as PNG and PDF
        TString pngname = Form("%s_%s.png", prefix, name.Data());
        //TString pdfname = Form("%s.pdf", name.Data());
        c->SaveAs(pngname);
        //c->SaveAs(pdfname);
    }

    f->Close();
}

