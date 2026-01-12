#include <TH2D.h>
#include <TNtuple.h>
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

void adaptiveSplit(TH2D* h, int target, std::vector<Region>& regions, Region region, int axis, int depth = 0) {
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

      adaptiveSplit(h, target, regions, q1, axis, depth+1);
      adaptiveSplit(h, target, regions, q2, axis, depth+1);
      adaptiveSplit(h, target, regions, q3, axis, depth+1);
      adaptiveSplit(h, target, regions, q4, axis, depth+1);

    }
    else {

      int x11; int x21; int y11; int y21;
      int x12; int x22; int y12; int y22;

      if (axis == 0) {
        x11 = region.x1; x21 = midx;
        x12 = midx + 1; x22 = region.x2;

        y11 = region.y1; y21 = region.y2;
        y12 = region.y1; y22 = region.y2;
      }
      else {
        x11 = region.x1; x21 = region.x2;
        x12 = region.x1; x22 = region.x2;

        y11 = region.y1; y21 = midy;
        y12 = midy + 1; y22 = region.y2;
      }

      Region s1 = {x11, x21, y11, y21};
      Region s2 = {x12, x22, y12, y22};

      bool Splits1 = regionIntegral(h, s1) > target;
      bool Splits2 = regionIntegral(h, s2) > target;
      if (Splits1 && Splits2) {
        regions.push_back(s1);
        regions.push_back(s2);
        return;
        //adaptiveSplit(h, target, regions, s1, axis, depth+1);
        //adaptiveSplit(h, target, regions, s2, axis, depth+1);

      }
      else {
       regions.push_back(region);
       return;
      }
    }
}

// In this version, try to make the binning symmetric for x vs theta_xw based on the positive angle binning

void makeAdaptiveBinningFineV3(const char* input_file, const char* output_file, int target=100000, bool isX=false, int axis=0) {
    gROOT->SetBatch(kTRUE);


    TNtuple* tree = new TNtuple("bin_tree", "bin_tree", "idx:x1:x2:y1:y2:counts");

    TFile* f = TFile::Open(input_file, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);

    for (int tpc = 0; tpc < 2; tpc++) {
      for (int plane = 0; plane < 3; plane++) {

        int idx = 3 * tpc + plane;
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
        int starty_bin = 1;
        int endy_bin = h->GetNbinsY()+1;

        if (isX) {

          starty_bin = h->GetYaxis()->FindBin(0.);
          //starty_bin = h->GetYaxis()->FindBin(-80.);
          endy_bin = h->GetYaxis()->FindBin(80.) - 1;
          if (tpc == 0) {
            startx_bin = 1;
            //endx_bin = h->GetNbinsX()/2;
            endx_bin = h->GetXaxis()->FindBin(0.)-1;
          }
          else {
            startx_bin = h->GetXaxis()->FindBin(0.)+1;
            //startx_bin = h->GetNbinsX()/2 + 1;
            endx_bin = h->GetNbinsX() + 1;
          }
        }

        std::vector<Region> regions;
        Region full = {startx_bin, endx_bin, starty_bin, endy_bin};
        adaptiveSplit(h, target, regions, full, axis);

        std::cout << "Made it through splitting algorithm!" << std::endl;


        TH2D* h_test = (TH2D*)h->Clone(Form("h%d", idx));
        TH2D* h_test_neg = (TH2D*)h->Clone(Form("h%d_neg", idx));
        h_test->Reset();
        h_test_neg->Reset();
        for (const auto& r : regions) {
          int y1 = r.y1;
          int y2 = r.y2;
          if (r.y2 < r.y1) {
            std::cout << "Swap Y" << std::endl;
            y1 = r.y2;
            y2 = r.y1;
          }
          double content = regionIntegral(h, r);
          //if (r.x1 == 1) std::cout << "Found x1 == 1" << std::endl;
          tree->Fill(idx, r.x1, r.x2, y1, y2, content);
          if (y2 < y1) {
            std::cout << "PROBLEM!!!" << std::endl;
            break;
          }
          std::cout << "x1 " << r.x1 << " x2 " << r.x2 << " y1 " << r.y1 << " y2 " << r.y2 << std::endl;
          int y1_neg;
          int y2_neg;
          // Now, we need to make negative theta regions for X
          if (isX) {
            double y1v = -1*h->GetYaxis()->GetBinCenter(y1);
            double y2v = -1*h->GetYaxis()->GetBinCenter(y2);
         
            if (y1v > y2v) {
              y2_neg = h->GetYaxis()->FindBin(y1v);    
              y1_neg = h->GetYaxis()->FindBin(y2v);
            }
            else {    
              y2_neg = h->GetYaxis()->FindBin(y2v);    
              y1_neg = h->GetYaxis()->FindBin(y1v);
            }
            tree->Fill(idx, r.x1, r.x2, y1_neg, y2_neg, content);
            std::cout << "x1 " << r.x1 << " x2 " << r.x2 << " y1_neg " << y1_neg << " y2_neg " << y2_neg << std::endl;
            std::cout << std::endl;
          }
          for (int i = r.x1; i <= r.x2; ++i) {
            for (int j = y1; j <= y2; ++j) {
              h_test->SetBinContent(i, j, content);
              h_test_neg->SetBinContent(i, j, content);
            }
          } 
          if (isX) {
            for (int i = r.x1; i <= r.x2; ++i) {
              for (int j = y1_neg; j <= y2_neg; ++j) {
                h_test_neg->SetBinContent(i, j, content);
              }
            }
          } 
        }
        outfile->cd();
        h_test->SetStats(0);
        h_test->Write();
        h_test->Delete();
        h_test_neg->SetStats(0);
        h_test_neg->Write();
        h_test_neg->Delete();
        h->Delete();
      }
    }

    tree->SetDirectory(0);
    outfile->cd();
    tree->Write();
    outfile->Close();

    delete f;
}

