#include <TH2D.h>
#include <vector>
#include <algorithm>
#include <iostream>

struct Region {
    int x1, x2;
    int y1, y2;
};


double regionIntegral(TH2D* h, const Region& r) {
    return h->Integral(r.x1, r.x2, r.y1, r.y2);
}


void adaptiveSplit(TH2D* h, int target, std::vector<Region>& regions, Region region, int depth = 0) {
    double total = regionIntegral(h, region);
    int nx = region.x2 - region.x1 + 1;
    int ny = region.y2 - region.y1 + 1;

    // can't split further if only one bin left in each dimension
    if (nx <= 1 && ny <= 1) {
        regions.push_back(region);
        return;
    }

    // midpoint bin indices
    int midx = (region.x1 + region.x2) / 2;
    int midy = (region.y1 + region.y2) / 2;

    Region q1 = {region.x1, midx, midy+1, region.y2};
    Region q2 = {midx+1, region.x2, midy+1, region.y2};
    Region q3 = {region.x1, midx, region.y1, midy};
    Region q4 = {midx+1, region.x2, region.y1, midy};

    bool Splitq1 = regionIntegral(h, q1) > target;
    bool Splitq2 = regionIntegral(h, q2) > target;
    bool Splitq3 = regionIntegral(h, q3) > target;
    bool Splitq4 = regionIntegral(h, q4) > target;

    if (Splitq1 && Splitq2 && Splitq3 && Splitq4) {

      adaptiveSplit(h, target, regions, q1, depth+1);
      adaptiveSplit(h, target, regions, q2, depth+1);
      adaptiveSplit(h, target, regions, q3, depth+1);
      adaptiveSplit(h, target, regions, q4, depth+1);

    }
    else {
     regions.push_back(region);
     return;
    }

}

void makeAdaptiveBinningQuad(const char* input_file, const char* output_file, int target=100000, bool isX=true) {
    gROOT->SetBatch(kTRUE);

    TFile* f = TFile::Open(input_file, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);

    int tpc = 0;
    int plane = 2;

    // Get the histogram
    int idx = 2;
    std::cout << "Fitting Index: " << idx << std::endl;
    THnSparseD* hTrk = (THnSparseD*)f->Get(Form("hTrack%d", idx));
    if (!hTrk) {
      std::cerr << "Histogram not found: " << idx << std::endl;
      f->Close();
      return;
    }
    TH2D* h = hTrk->Projection(1, 0);

    int startx_bin = 1;
    int endx_bin = h->GetNbinsX()+1;

    if (isX) {
      if (tpc == 0) {
        startx_bin = 1;
        endx_bin = h->GetNbinsX()/2;
      }
      else {
        startx_bin = h->GetNbinsX()/2 + 1;
        endx_bin = h->GetNbinsX() + 1;
      }
    }

    std::vector<Region> regions;
    Region full = {startx_bin, endx_bin, 1, h->GetNbinsY()+1};
    adaptiveSplit(h, target, regions, full);

    // collect unique edges
    std::vector<double> xedges, yedges;
    //xedges.push_back(h->GetXaxis()->GetBinLowEdge(startx_bin));
    //yedges.push_back(h->GetYaxis()->GetBinLowEdge(1));

    for (const auto& r : regions) {
        xedges.push_back(h->GetXaxis()->GetBinLowEdge(r.x1));
        xedges.push_back(h->GetXaxis()->GetBinUpEdge(r.x2));
        yedges.push_back(h->GetYaxis()->GetBinLowEdge(r.y1));
        yedges.push_back(h->GetYaxis()->GetBinUpEdge(r.y2));
    }

    std::sort(xedges.begin(), xedges.end());
    std::sort(yedges.begin(), yedges.end());
    xedges.erase(std::unique(xedges.begin(), xedges.end()), xedges.end());
    yedges.erase(std::unique(yedges.begin(), yedges.end()), yedges.end());

    // build new histogram
    TH2D* hadapt = new TH2D(Form("%s_adapt", h->GetName()), 
                             Form("Adaptive binning of %s", h->GetTitle()),
                             xedges.size()-1, xedges.data(),
                             yedges.size()-1, yedges.data());

    // fill new histogram
    for (const auto& r : regions) {
        double content = regionIntegral(h, r);
        double xcen = 0.5*(h->GetXaxis()->GetBinCenter(r.x1) + h->GetXaxis()->GetBinCenter(r.x2));
        double ycen = 0.5*(h->GetYaxis()->GetBinCenter(r.y1) + h->GetYaxis()->GetBinCenter(r.y2));
        hadapt->Fill(xcen, ycen, content);
    }

    hadapt->SetStats(0);

    TH2D* h_test = (TH2D*)h->Clone("h_test");
    h_test->Reset();
    for (const auto& r : regions) {
        double content = regionIntegral(h, r);
        for (int i = r.x1; i <= r.x2; ++i) {
          for (int j = r.y1; j <= r.y2; ++j) {
            h_test->SetBinContent(i, j, content);
          }
        } 
    }


    outfile->cd();

    h->Write();
    hadapt->Write();
    h_test->Write();

    outfile->Close();

    delete h;
    delete hadapt;
    delete f;
}

