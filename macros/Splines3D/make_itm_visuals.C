#include <iostream>
#include <vector>
#include <tuple>
#include <mutex>
#include <numeric>
#include <ROOT/TThreadExecutor.hxx>
#include "TFile.h"
#include "TNtuple.h"
#include "THnSparse.h"
#include "TH2.h"
#include "TH1.h"
#include "TGraph2DErrors.h"
#include "../../include_wire/Fitting.h"

const UInt_t kNplanes = 3;
const UInt_t kNTPCs   = 2;

// Struct for preloaded bin regions
struct BinRegion {
    int idx, X, x1, x2, y1, y2;
};

void make_itm_visuals(const char* mc_file,
                      const char* data_file,
                      const char* output_file,
                      const char* bin_file, int IDX, int num_plots=100)
{

    // Open bin file and preload TNtuple into memory
    TFile fbin(bin_file, "READ");
    TNtuple* t = (TNtuple*)fbin.Get("bin_tree");
    if (!t) {
        std::cerr << "No bin_tree found!" << std::endl;
        return;
    }

    float idx, X, x1, x2, y1, y2, counts;
    t->SetBranchAddress("idx", &idx);
    t->SetBranchAddress("X", &X);
    t->SetBranchAddress("x1", &x1);
    t->SetBranchAddress("x2", &x2);
    t->SetBranchAddress("y1", &y1);
    t->SetBranchAddress("y2", &y2);
    t->SetBranchAddress("counts", &counts);

    Long64_t nentries = t->GetEntries();
    std::vector<BinRegion> regions;
    //regions.reserve(nentries);
 

    int count = 0;   
    for (Long64_t i = 0; i < nentries; ++i) {

     
        t->GetEntry(i);
        if ((int)idx != IDX) continue;
        //if ((int)X > 10) continue; //TODO --> Fix this!
        BinRegion R = {(int)idx, (int)X, (int)x1, (int)x2, (int)y1, (int)y2};

        // Each thread opens its own files (safe for read)
        std::unique_ptr<TFile> fmc(TFile::Open(mc_file, "READ"));
        std::unique_ptr<TFile> fdata(TFile::Open(data_file, "READ"));

        auto hmc   = (THnSparseD*)fmc->Get(Form("hHit%d", R.idx));
        auto hdata = (THnSparseD*)fdata->Get(Form("hHit%d", R.idx));
        if (!hmc || !hdata) return;

        // Set the X-Range
        hmc->GetAxis(0)->SetRange(R.X, R.X);
        hdata->GetAxis(0)->SetRange(R.X, R.X);
        
        std::cout << "Made the X projections " <<  std::endl;

        auto h2D = hmc->Projection(2, 1);
        bool isX0 = (std::string(h2D->GetXaxis()->GetTitle()) == "x");
        bool isX1 = (std::string(h2D->GetYaxis()->GetTitle()) == "x");
        
        std::cout << "Made the YZ projection for axes info" << std::endl;

        // Apply ranges
        if (R.x2 > R.x1) {
          hmc->GetAxis(1)->SetRange(R.x1, R.x2);
          hdata->GetAxis(1)->SetRange(R.x1, R.x2);
        }
        else {
          hmc->GetAxis(1)->SetRange(R.x2, R.x1);
          hdata->GetAxis(1)->SetRange(R.x2, R.x1);
        }
        if (R.y2 > R.y1) {
          hmc->GetAxis(2)->SetRange(R.y1, R.y2);
          hdata->GetAxis(2)->SetRange(R.y1, R.y2);
        }
        else {
          hmc->GetAxis(2)->SetRange(R.y2, R.y1);
          hdata->GetAxis(2)->SetRange(R.y2, R.y1);
        }

        std::unique_ptr<TH1D> hmc_1d(hmc->Projection(3));
        std::unique_ptr<TH1D> hdata_1d(hdata->Projection(3));

        std::cout << "Made the hit dist projections" << std::endl;

        if (hmc_1d->Integral() < 1) {
          std::cout << "X " << R.X << " idx " << R.idx << std::endl;
          return;
        }

        double itm_result_mc[2];
        double itm_result_data[2];
        iterative_truncated_mean_std_err(hmc_1d.get(), -2, 1.75, 1e-4, itm_result_mc);
        iterative_truncated_mean_std_err(hdata_1d.get(), -2, 1.75, 1e-4, itm_result_data);

        std::cout << "Got ITM results" << std::endl;

        double xcen = 0.5*(h2D->GetXaxis()->GetBinCenter(R.x1) + h2D->GetXaxis()->GetBinCenter(R.x2));
        double ycen = 0.5*(h2D->GetYaxis()->GetBinCenter(R.y1) + h2D->GetYaxis()->GetBinCenter(R.y2));
        //double xerr = std::abs(h2D->GetXaxis()->GetBinCenter(R.x2) - h2D->GetXaxis()->GetBinCenter(R.x1)) / 2.0 + h2D->GetXaxis()->GetBinWidth(1)/2.0;
        //double yerr = std::abs(h2D->GetYaxis()->GetBinCenter(R.y2) - h2D->GetYaxis()->GetBinCenter(R.y1)) / 2.0 + h2D->GetYaxis()->GetBinWidth(1)/2.0;
        double ratio = itm_result_data[0] / itm_result_mc[0];

        std::cout << "computed ratio and position" << std::endl;


        // make plot
        TCanvas* c = new TCanvas(Form("c%d", count), "", 700, 500);
        hmc_1d->Scale(1./hmc_1d->Integral());  
        hdata_1d->Scale(1./hdata_1d->Integral());  
        hmc_1d->SetLineColor(4);
        hdata_1d->SetLineColor(1);
        hdata_1d->Draw("HISTE");
        hmc_1d->Draw("HISTE Same");
        

        c->Write();

        delete h2D;
        delete hmc;
        delete hdata;
        delete c;
        count += 1;

    }


}

