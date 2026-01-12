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
#include "../../../include_wire/Fitting.h"

const UInt_t kNplanes = 3;
const UInt_t kNTPCs   = 2;

// Struct for preloaded bin regions
struct BinRegion {
    int idx, x1, x2, y1, y2;
};

void make_2d_splines_prod_multi_v2(const char* mc_file,
                                const char* data_file,
                                const char* output_file,
                                const char* bin_file)
{
    gROOT->SetBatch(kTRUE);
    TH1::AddDirectory(kFALSE);
    ROOT::EnableImplicitMT();  // enable multithreading

    // Open bin file and preload TNtuple into memory
    TFile fbin(bin_file, "READ");
    TNtuple* t = (TNtuple*)fbin.Get("bin_tree");
    if (!t) {
        std::cerr << "No bin_tree found!" << std::endl;
        return;
    }

    float idx, x1, x2, y1, y2, counts;
    t->SetBranchAddress("idx", &idx);
    t->SetBranchAddress("x1", &x1);
    t->SetBranchAddress("x2", &x2);
    t->SetBranchAddress("y1", &y1);
    t->SetBranchAddress("y2", &y2);
    t->SetBranchAddress("counts", &counts);

    Long64_t nentries = t->GetEntries();
    std::vector<BinRegion> regions;
    regions.reserve(nentries);

    for (Long64_t i = 0; i < nentries; ++i) {
        t->GetEntry(i);
        regions.push_back({(int)idx, (int)x1, (int)x2, (int)y1, (int)y2});
    }
    std::cout << "Preloaded " << regions.size() << " regions." << std::endl;

    // Parallel section
    std::mutex result_mutex;
    std::vector<std::vector<std::tuple<double,double,double,double,double,double,double>>> results(6);

    auto process_region = [&](const BinRegion& R) {
        // Each thread opens its own files (safe for read)
        std::unique_ptr<TFile> fmc(TFile::Open(mc_file, "READ"));
        std::unique_ptr<TFile> fdata(TFile::Open(data_file, "READ"));

        auto hmc   = (THnSparseD*)fmc->Get(Form("hHit%d", R.idx));
        auto hdata = (THnSparseD*)fdata->Get(Form("hHit%d", R.idx));
        if (!hmc || !hdata) return;

        auto h2D = hmc->Projection(1, 0);
        bool isX0 = (std::string(h2D->GetXaxis()->GetTitle()) == "x");
        bool isX1 = (std::string(h2D->GetYaxis()->GetTitle()) == "x");
        
        // Apply ranges
        if (R.x2 > R.x1) {
          hmc->GetAxis(0)->SetRange(R.x1, R.x2);
          hdata->GetAxis(0)->SetRange(R.x1, R.x2);
        }
        else {
          hmc->GetAxis(0)->SetRange(R.x2, R.x1);
          hdata->GetAxis(0)->SetRange(R.x2, R.x1);
        }
        if (R.y2 > R.y1) {
          hmc->GetAxis(1)->SetRange(R.y1, R.y2);
          hdata->GetAxis(1)->SetRange(R.y1, R.y2);
        }
        else {
          hmc->GetAxis(1)->SetRange(R.y2, R.y1);
          hdata->GetAxis(1)->SetRange(R.y2, R.y1);
        }

        std::unique_ptr<TH1D> hmc_1d(hmc->Projection(2));
        std::unique_ptr<TH1D> hdata_1d(hdata->Projection(2));

        double itm_result_mc[2];
        double itm_result_data[2];
        iterative_truncated_mean_std_err(hmc_1d.get(), -2, 1.75, 1e-4, itm_result_mc);
        iterative_truncated_mean_std_err(hdata_1d.get(), -2, 1.75, 1e-4, itm_result_data);

        double xcen = 0.5*(h2D->GetXaxis()->GetBinCenter(R.x1) + h2D->GetXaxis()->GetBinCenter(R.x2));
        double ycen = 0.5*(h2D->GetYaxis()->GetBinCenter(R.y1) + h2D->GetYaxis()->GetBinCenter(R.y2));
        double xerr = std::abs(h2D->GetXaxis()->GetBinCenter(R.x2) - h2D->GetXaxis()->GetBinCenter(R.x1)) / 2.0 + h2D->GetXaxis()->GetBinWidth(1)/2.0;
        double yerr = std::abs(h2D->GetYaxis()->GetBinCenter(R.y2) - h2D->GetYaxis()->GetBinCenter(R.y1)) / 2.0 + h2D->GetYaxis()->GetBinWidth(1)/2.0;
        double ratio = itm_result_data[0] / itm_result_mc[0];

        // Determine if the region is on the border or is a corner
        bool onLeft   = ((R.x1 == 1) || (R.x2 == 1));
        //if (onLeft) std::cout << "Orig. On Left" << std::endl;
        bool onRight  = ((R.x2 == h2D->GetNbinsX()+1) || (R.x1 == h2D->GetNbinsX()+1));
        bool onBottom = ((R.y1 == 1) || (R.y2 == 1));
        bool onTop    = ((R.y2 == h2D->GetNbinsY()+1) || (R.y1 == h2D->GetNbinsY()+1));

        if (isX0) {
          int zero_bin = h2D->GetXaxis()->FindBin(0.);
          if (R.idx < 3) {
            onRight = ((R.x2 == zero_bin-1) || (R.x1 == zero_bin-1));
          }
          else {
            onLeft = ((R.x2 == zero_bin+1) || (R.x1 == zero_bin+1));
          }
          // Changed Theta XZ to cut out high-angle bins
          int top_bin = h2D->GetYaxis()->FindBin(79.5);
          int bot_bin = h2D->GetYaxis()->FindBin(-79.5);
          onBottom = ((R.y1 == bot_bin) || (R.y2 == bot_bin));
          onTop = ((R.y1 == top_bin) || (R.y2 == top_bin));
        }
        if (isX1) {
          std::cout << "Wrong! Put X label on the x-axis" <<std::endl;
          return;
        }

        bool corner_upleft = ( onLeft && onTop);
        bool corner_lowleft = ( onLeft && onBottom);
  
        bool corner_upright = ( onRight && onTop);
        bool corner_lowright = ( onRight && onBottom);

        {
            std::lock_guard<std::mutex> lock(result_mutex);
            results[R.idx].push_back({xcen, ycen,
                                      itm_result_mc[0], itm_result_data[0], ratio,
                                      itm_result_mc[1], itm_result_data[1]});

            // add points for edge cases
            if (onLeft) {
              //std::cout << "On Left" << std::endl;
              double xcen_t = xcen - xerr;
              //std::lock_guard<std::mutex> lock(result_mutex);
              results[R.idx].push_back({xcen_t, ycen,
                                        itm_result_mc[0], itm_result_data[0], ratio,
                                        itm_result_mc[1], itm_result_data[1]});
              
            }
            if (onRight) {
              //std::cout << "On Right" << std::endl;
              double xcen_t = xcen + xerr;
              //std::lock_guard<std::mutex> lock(result_mutex);
              results[R.idx].push_back({xcen_t, ycen,
                                        itm_result_mc[0], itm_result_data[0], ratio,
                                        itm_result_mc[1], itm_result_data[1]});

            }
            if (onTop) {
              //std::cout << "On Top" << std::endl;
              double ycen_t = ycen + yerr;
              //std::lock_guard<std::mutex> lock(result_mutex);
              results[R.idx].push_back({xcen, ycen_t,
                                        itm_result_mc[0], itm_result_data[0], ratio,
                                        itm_result_mc[1], itm_result_data[1]});
             
            }
            if (onBottom) {
              //std::cout << "On Bottom" << std::endl;
              double ycen_t = ycen - yerr;
              //std::lock_guard<std::mutex> lock(result_mutex);
              results[R.idx].push_back({xcen, ycen_t,
                                        itm_result_mc[0], itm_result_data[0], ratio,
                                        itm_result_mc[1], itm_result_data[1]});
            }
            if (corner_lowright) {
              double xcen_t = xcen + xerr;
              double ycen_t = ycen - yerr;
              //std::lock_guard<std::mutex> lock(result_mutex);
              results[R.idx].push_back({xcen_t, ycen_t,
                                        itm_result_mc[0], itm_result_data[0], ratio,
                                        itm_result_mc[1], itm_result_data[1]});

            }
            if (corner_upright) {
              double xcen_t = xcen + xerr;
              double ycen_t = ycen + yerr;
              //std::lock_guard<std::mutex> lock(result_mutex);
              results[R.idx].push_back({xcen_t, ycen_t,
                                        itm_result_mc[0], itm_result_data[0], ratio,
                                        itm_result_mc[1], itm_result_data[1]});

            }
            if (corner_lowleft) {
              double xcen_t = xcen - xerr;
              double ycen_t = ycen - yerr;
              //std::lock_guard<std::mutex> lock(result_mutex);
              results[R.idx].push_back({xcen_t, ycen_t,
                                        itm_result_mc[0], itm_result_data[0], ratio,
                                        itm_result_mc[1], itm_result_data[1]});

            }
            if (corner_upleft) {
              double xcen_t = xcen - xerr;
              double ycen_t = ycen + yerr;
              //std::lock_guard<std::mutex> lock(result_mutex);
              results[R.idx].push_back({xcen_t, ycen_t,
                                        itm_result_mc[0], itm_result_data[0], ratio,
                                        itm_result_mc[1], itm_result_data[1]});
     
            }

        }


        delete h2D;
        delete hmc;
        delete hdata;
    };

    ROOT::TThreadExecutor pool;
    pool.Foreach(process_region, regions);

    // Merge results into graphs
    TFile outfile(output_file, "RECREATE");
    TGraph2DErrors* mc_2d_fits[6];
    TGraph2DErrors* data_2d_fits[6];
    TGraph2DErrors* splines_2d[6];
    
    std::unique_ptr<TFile> f(TFile::Open(mc_file, "READ"));
    auto h = (THnSparseD*)f->Get(Form("hHit%d", 0));
    if (!h) return;
    auto haxes = h->Projection(1, 0);
    std::cout << "X-axis title " << haxes->GetXaxis()->GetTitle() << std::endl;
    std::cout << "Y-axis title " << haxes->GetYaxis()->GetTitle() << std::endl;
    for (int i = 0; i < 6; ++i) {
        mc_2d_fits[i]   = new TGraph2DErrors();
        data_2d_fits[i] = new TGraph2DErrors();
        splines_2d[i]   = new TGraph2DErrors();

        mc_2d_fits[i]->SetName(Form("mc_mpv_%d", i));
        data_2d_fits[i]->SetName(Form("data_mpv_%d", i));
        splines_2d[i]->SetName(Form("splines_%d", i));

        mc_2d_fits[i]->GetXaxis()->SetTitle(haxes->GetXaxis()->GetTitle());
        mc_2d_fits[i]->GetYaxis()->SetTitle(haxes->GetYaxis()->GetTitle());
        
        data_2d_fits[i]->GetXaxis()->SetTitle(haxes->GetXaxis()->GetTitle());
        data_2d_fits[i]->GetYaxis()->SetTitle(haxes->GetYaxis()->GetTitle());
        
        splines_2d[i]->GetXaxis()->SetTitle(haxes->GetXaxis()->GetTitle());
        splines_2d[i]->GetYaxis()->SetTitle(haxes->GetYaxis()->GetTitle());

        for (auto& r : results[i]) {
            double x, y, mpv_mc, mpv_data, ratio, err_mc, err_data;
            std::tie(x, y, mpv_mc, mpv_data, ratio, err_mc, err_data) = r;

            int n1 = mc_2d_fits[i]->GetN();
            mc_2d_fits[i]->SetPoint(n1, x, y, mpv_mc);
            mc_2d_fits[i]->SetPointError(n1, 0, 0, err_mc);

            int n2 = data_2d_fits[i]->GetN();
            data_2d_fits[i]->SetPoint(n2, x, y, mpv_data);
            data_2d_fits[i]->SetPointError(n2, 0, 0, err_data);

            int n3 = splines_2d[i]->GetN();
            splines_2d[i]->SetPoint(n3, x, y, ratio);
            splines_2d[i]->SetPointError(n3, 0, 0, 0);
        }

        outfile.cd();
        mc_2d_fits[i]->Write();
        data_2d_fits[i]->Write();
        splines_2d[i]->Write();
    }
    delete h;
    delete haxes;
    std::cout << "✅ Finished parallel fitting and writing results.\n";
}

