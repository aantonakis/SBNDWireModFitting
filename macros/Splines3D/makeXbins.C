#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TString.h"

#include <vector>
#include <utility>
#include <iostream>

// Recursively split [b1,b2] (inclusive) if both halves pass threshold.
// Leaves are appended to `leaves` as (b1,b2) bin-index ranges.
static void SplitHalvesByCounts(const TH1* hx,
                                int b1, int b2,
                                double threshold,
                                std::vector<std::pair<int,int>>& leaves)
{

  //if (b1 >= b2) {
  //  leaves.emplace_back(b1, b2);
  //  return;
  //}

  if (b1 == b2) {
    leaves.emplace_back(b1, b2);
    return;
  }
  int mid = (b1 + b2) / 2; // floor
  //if (mid < b1 || mid >= b2) {
  //  leaves.emplace_back(b1, b2);
  //  return;
 // }

  const double cL = hx->Integral(b1, mid);
  const double cR = hx->Integral(mid + 1, b2);

  std::cout << "cL " << cL << std::endl; 
  std::cout << "cR " << cR << std::endl; 

  if (cL >= threshold && cR >= threshold) {
    std::cout << "Split again!" << std::endl;
    SplitHalvesByCounts(hx, b1, mid, threshold, leaves);
    SplitHalvesByCounts(hx, mid + 1, b2, threshold, leaves);
  } else {
    leaves.emplace_back(b1, b2);
  }
}

// Main entry point.
// Example call (ROOT prompt):
//   .x MakeAdaptiveXBinning.C("in.root","hx","out.root","hx_adapt",2000)
void makeXbins(const char* inFile,
                          const char* outFile, int tpc, int plane,
                          double threshold)
{
  TFile fin(inFile, "READ");
  if (fin.IsZombie()) {
    std::cerr << "ERROR: could not open input file " << inFile << "\n";
    return;
  }



  int idx = 3 * tpc + plane;
  std::cout << "Binning Index: " << idx << std::endl;
  THnSparseD* hTrk = (THnSparseD*)fin.Get(Form("hTrack%d", idx));
  if (!hTrk) {
    std::cerr << "Histogram not found: " << idx << std::endl;
    fin.Close();
    return;
  }


  TH1* hx = hTrk->Projection(0); 


  // Work with x bins only (exclude under/overflow): 1..NbinsX
  const int nb = hx->GetNbinsX();
  const TAxis* ax = hx->GetXaxis();
  std::cout << "Number of bins X " << nb << std::endl;
  // Build leaf bin-index ranges
  std::vector<std::pair<int,int>> leaves;
  //leaves.reserve(nb);

  int bin_zero = hx->FindBin(0.);
  if (tpc == 0) {

    SplitHalvesByCounts(hx, 1, bin_zero, threshold, leaves);
  }
  else {
    SplitHalvesByCounts(hx, bin_zero, nb+1, threshold, leaves);
  }

  std::cout << "Finished Splitting" << std::endl;

  // Convert leaf ranges to variable bin edges using the ORIGINAL axis edges.
  // Each leaf [b1,b2] becomes a new bin with edges:
  //   low = ax->GetBinLowEdge(b1)
  //   up  = ax->GetBinUpEdge(b2)
  std::vector<double> edges;
  //edges.reserve(leaves.size() + 1);

  // Sort is not necessary with this recursion (it preserves order),
  // but we'll assume it's already left-to-right.
  for (size_t i = 0; i < leaves.size(); ++i) {
    const int b1 = leaves[i].first;
    const int b2 = leaves[i].second;

    const double low = ax->GetBinLowEdge(b1);
    const double up  = ax->GetBinUpEdge(b2);

    if (i == 0) edges.push_back(low);
    edges.push_back(up);
  }

  // Create the new variable-bin histogram
  TH1D* hnew = new TH1D(Form("hx_%d", idx),
                        "",
                        (int)edges.size() - 1,
                        edges.data());

  hnew->Sumw2();

  // Fill the new histogram by integrating original hx over each leaf range.
  for (size_t i = 0; i < leaves.size(); ++i) {
    const int b1 = leaves[i].first;
    const int b2 = leaves[i].second;

    const double c = hx->Integral(b1, b2);

    // (Optional) propagate errors in quadrature by summing bin errors^2
    double err2 = 0.0;
    for (int b = b1; b <= b2; ++b) {
      const double e = hx->GetBinError(b);
      err2 += e * e;
    }

    const int newBin = (int)i + 1;
    hnew->SetBinContent(newBin, c);
    hnew->SetBinError(newBin, std::sqrt(err2));
  }

  // Print summary
  std::cout << "Original bins: " << nb << "\n";
  std::cout << "Adaptive bins: " << hnew->GetNbinsX() << "\n";
  std::cout << "Threshold:     " << threshold << "\n";
  for (size_t i = 0; i < leaves.size(); ++i) {
    int b1 = leaves[i].first, b2 = leaves[i].second;
    std::cout << "  newBin " << (i+1)
              << ": old bins [" << b1 << "," << b2 << "]"
              << "  x=[" << ax->GetBinLowEdge(b1) << "," << ax->GetBinUpEdge(b2) << "]"
              << "  counts=" << hx->Integral(b1,b2) << "\n";
  }

  // Save
  TFile fout(outFile, "RECREATE");
  hnew->Write();
  fout.Close();
  fin.Close();

}

