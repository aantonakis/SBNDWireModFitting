

#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>

#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFile.h"
#include "THnSparse.h"
#include "TNtuple.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVector3.h"
#include "TPRegexp.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "Math/Vector3D.h"

#include "../../../include_wire/Fitting.h"


const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;


void make_spline_1d(const char* input_mc, const char* input_data, const char* output_file) {
    gROOT->SetBatch(kTRUE);

    TFile* f_mc = TFile::Open(input_mc, "READ");
    if (!f_mc || f_mc->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }
    TFile* f_data = TFile::Open(input_data, "READ");
    if (!f_data || f_data->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);
    
    for (int tpc = 0; tpc < kNTPCs; tpc++) {
        for (int plane = 0; plane < kNplanes; plane++) {

            // Get the histogram
            int idx = 3 * tpc + plane;
            //std::string num_str = "hwidth"+std::to_string(idx); // convert int to string
            std::string num_str = "h1D"+std::to_string(idx); // convert int to string
            const char* cstr = num_str.c_str(); 
            
	    std::string mpv_str = "MPV"+std::to_string(idx); // convert int to string

            const char* cstr_mpv = mpv_str.c_str(); 
            TH1D* h_mc = (TH1D*)f_mc->Get(cstr_mpv); 
            TH1D* h_data = (TH1D*)f_data->Get(cstr_mpv); 

            if (!h_mc || !h_data) {
                std::cerr << "Histogram not found: " << idx << std::endl;
                f_mc->Close();
                f_data->Close();
                return;
            }

            TH1D* h_ratio = (TH1D*)h_data->Clone(Form("h_%s_ratio_%d", h_mc->GetYaxis()->GetTitle(), idx));
	
	    h_ratio->Divide(h_mc);
	   
	    h_ratio->GetXaxis()->SetTitle( h_mc->GetXaxis()->GetTitle());
	    
	    h_ratio->GetYaxis()->SetTitle("Data/MC");

	    // Create arrays for the spline input
    	    const int nbins = h_ratio->GetNbinsX();
    	    std::vector<double> x, y;

            for (int i = 1; i <= nbins; ++i) {
              double binCenter = h_ratio->GetBinCenter(i);
              double binContent = h_ratio->GetBinContent(i);
              if (binContent > 0) { // ignore empty bins if desired
                x.push_back(binCenter);
                y.push_back(binContent);
              }
            }

            // Build a TSpline3 (cubic spline)
            TSpline3 *spline = new TSpline3(Form("temp_%d", idx), &x[0], &y[0], x.size());
	    spline->SetName(Form("spline_%s_%d", h_mc->GetYaxis()->GetTitle(), idx));
            spline->SetLineColor(kRed);

            spline->SetLineWidth(2);

            outfile->cd();
           
            h_ratio->Write();
	    spline->Write(spline->GetName());	   


            
	    TCanvas* cmpv = new TCanvas(Form("cmpv_%d", idx), "", 700, 500);
            h_data->SetStats(0);
            h_data->GetYaxis()->SetRangeUser(0, 1.15*h_data->GetMaximum());
            h_data->SetLineColor(4);
            h_mc->SetLineColor(2);
            h_data->Draw("HISTE");
            h_mc->Draw("HISTE Same");
            TLegend* leg = new TLegend(0.3, 0.6, 0.7, 0.85);
            leg->AddEntry(h_data, "Data");
            leg->AddEntry(h_mc, "MC");
            leg->Draw("Same");
	    cmpv->Write();
	    
	    TCanvas* c = new TCanvas(Form("c_%d", idx), "", 700, 500);
	    h_ratio->GetYaxis()->SetRangeUser(0.85, 1.15);
	    h_ratio->Draw("HISTE");
	    spline->Draw("Same");
	    //cq->Update();
	    c->Write();	    

            delete h_mc;
            delete h_data;
            delete h_ratio;
	    delete spline;
            delete leg;
           
        }
    }
 
    delete f_mc;
    delete f_data;
    outfile->Close();

    printf("Finished processing.\n");

}


