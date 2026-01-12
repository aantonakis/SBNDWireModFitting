

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
#include "../../../include_wire/Corr.h"


const UInt_t kNplanes = 3;
const UInt_t kNTPCs = 2;


// This code takes a 3D hist where the 3rd dimension is a Hit width/charge and projects to one 
// of the first two dimensions for a 1D spline

void profile_1d_itm_std_proj(const char* input_file, const char* output_file, int dim) {
    
    gROOT->SetBatch(kTRUE);

    TFile* f = TFile::Open(input_file, "READ");
    if (!f || f->IsZombie()) {
      std::cerr << "File is a zombie mate" << std::endl;
      return;
    }

    TFile* outfile = new TFile(output_file, "RECREATE");

    TH1::AddDirectory(0);

    
    
    for (int tpc = 0; tpc < kNTPCs; tpc++) {
        for (int plane = 0; plane < kNplanes; plane++) {
		
            // Get the histogram
            int idx = 3 * tpc + plane;
	    std::cout << "Processing idx " << idx << std::endl;
            //std::string num_str = "hwidth"+std::to_string(idx); // convert int to string
            std::string num_str = "hHit"+std::to_string(idx); // convert int to string
            const char* cstr = num_str.c_str(); 
            THnSparseD* h = (THnSparseD*)f->Get(cstr);
            if (!h) {
                std::cerr << "Histogram not found: " << idx << std::endl;
                f->Close();
                return;
            }

            TH2D* h2D = h->Projection(2, dim);

            TH1D* h_result_temp = h->Projection(dim);
            h_result_temp->Reset();
            h_result_temp->SetName(Form("MPV%d", idx));
	    h_result_temp->GetYaxis()->SetTitle(h->GetAxis(dim)->GetTitle());
                      
	    TDirectory* P = outfile->mkdir(Form("Fits%d", idx));

    	    // Get ITM for each slice in the specified dimension
            for (int i = 1; i < h->GetAxis(dim)->GetNbins() + 1; i++) {
              h->GetAxis(dim)->SetRange(i, i);
              TH1D* h_1d = h->Projection(2); // charge
		
              Double_t itm_result[2];

              iterative_truncated_mean_std_err(h_1d, -2, 1.75, 1.0e-4, itm_result);
              
	      h_result_temp->SetBinContent(i, itm_result[0]);
              h_result_temp->SetBinError(i, itm_result[1]);
	
	      P->cd();
	      float edge_low = h->GetAxis(dim)->GetBinCenter(i) - (h->GetAxis(dim)->GetBinWidth(i))/2.;
	      float edge_high = h->GetAxis(dim)->GetBinCenter(i) + (h->GetAxis(dim)->GetBinWidth(i))/2;

	      h_1d->GetXaxis()->SetTitle(h->GetAxis(dim)->GetTitle());

	      h_1d->SetTitle(Form("Bin Range: [%f, %f]", edge_low, edge_high));

	      TCanvas* c = new TCanvas(Form("c_%d_%d", idx, i), "", 700, 500);
	      h_1d->Draw("HISTE");
	      TLine* lq = new TLine(itm_result[0], 0, itm_result[0], h_1d->GetMaximum());
	      lq->SetLineColor(2);
	      lq->SetLineWidth(2);
	      lq->Draw("Same");
	      c->Write();
              
            }
	
            outfile->cd();
	    TCanvas* c = new TCanvas(Form("c%d", idx), "", 700, 500);
            c->SetLogz();
            h2D->SetTitle(Form("TPC %d, Plane %d", tpc, plane));
	    h2D->SetStats(0);
	    h2D->Draw("Colz");
            h_result_temp->SetLineWidth(3);
            h_result_temp->SetLineColor(2);
            h_result_temp->Draw("HISTE Same");
            c->Write();

            h_result_temp->Write();

            delete h2D;
            delete h_result_temp;
            delete h;
           
        }
    }
 
    delete f;
    outfile->Close();

    printf("Finished processing.\n");

}


