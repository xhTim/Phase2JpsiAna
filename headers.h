#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
#include <string>
#include <utility>
#include "math.h"
#include "sys/types.h"

#include "TROOT.h"
#include "TChain.h"
#include "TClass.h"
#include "TSystem.h"
#include "TUnixSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDLazy.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TAttFill.h"
#include "TString.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"
#include "TDatabasePDG.h"
#include "TParticlePDG.h"
#include "TPaletteAxis.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TSpectrum.h"
#include "TGaxis.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "THnSparse.h"
#include "TFitter.h"
#include "TFitResult.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TArrow.h"
#include "TArrayL.h"
#include "TArrayF.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TBrowser.h"
#include "TPad.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TLine.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TPDF.h"
#include "TPostScript.h"

//RooFit
//#include "RooGlobalFunc.h"
//#include "RooRealVar.h"
//#include "RooConstVar.h"
//#include "RooDataSet.h"
//#include "RooDataHist.h"
//#include "RooConstVar.h"
//#include "RooChebychev.h"
//#include "RooAddPdf.h"
//#include "RooSimultaneous.h"
//#include "RooCategory.h"
//#include "RooPlot.h"
//#include "RooGenericPdf.h"
//#include "RooCBShape.h"
//#include "RooGaussian.h"
//#include "RooExponential.h"
//#include "RooPolynomial.h"

#include "TVirtualFitter.h"
#include "TMinuit.h"

using namespace std;
//using namespace RooFit ;
