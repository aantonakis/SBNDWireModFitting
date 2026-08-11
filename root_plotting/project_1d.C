#include "TFile.h"
#include "TH1D.h"
#include "THnSparse.h"
#include "TString.h"

#include <iostream>
#include <vector>

void project_1d(const char* input_file,
                       int axis_to_project,
                       const char* output_file = "projected_output.root")
{
  TFile* fin = TFile::Open(input_file, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Error: could not open input file: " << input_file << std::endl;
    return;
  }

  TFile* fout = TFile::Open(output_file, "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error: could not create output file: " << output_file << std::endl;
    fin->Close();
    return;
  }

  std::vector<TString> prefixes = {"hTrack", "hHit"};

  for (const auto& prefix : prefixes) {
    for (int i = 0; i <= 5; ++i) {
      std::cout << "Projecting " << prefix << " " << i << std::endl;
      TString hname = Form("%s%d", prefix.Data(), i);

      THnSparseD* hs = dynamic_cast<THnSparseD*>(fin->Get(hname));
      if (!hs) {
        std::cerr << "Warning: could not find THnSparseD named "
                  << hname << std::endl;
        continue;
      }

      int ndim = hs->GetNdimensions();
      if (axis_to_project < 0 || axis_to_project >= ndim) {
        std::cerr << "Error: requested axis " << axis_to_project
                  << " is invalid for " << hname
                  << " with ndim = " << ndim << std::endl;
        continue;
      }

      TH1D* hproj = hs->Projection(axis_to_project);

      if (!hproj) {
        std::cerr << "Warning: projection failed for " << hname << std::endl;
        continue;
      }

      TString outname = Form("%s_proj_axis%d", hname.Data(), axis_to_project);
      hproj->SetName(outname);
      hproj->SetTitle(Form("%s projected onto axis %d",
                           hname.Data(), axis_to_project));

      fout->cd();
      hproj->Write();

      std::cout << "Wrote " << outname << std::endl;
    }
  }

  fout->Close();
  fin->Close();

  std::cout << "Done. Output written to " << output_file << std::endl;
}
