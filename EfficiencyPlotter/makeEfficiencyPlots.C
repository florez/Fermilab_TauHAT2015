// Load C++ libraries
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdio.h>
#include <string.h>

// Load ROOT library
 #include <TError.h>

// Load libraries (automatically included by ROOT)
#include <TFile.h>
#include <TH1.h>
#include <TFrame.h>
#include <TPad.h>
#include <TCollection.h>
#include <TDirectory.h>
#include <TKey.h>
#include <TClass.h>
#include <fstream>

void makeEfficiencyPlots(){
  
 ifstream inFile;
 inFile.open ("config.in", ios::in); 
 
 if (!inFile)
    {
      cerr << "ERROR: Can't open input file: " << endl;
      exit (1);
    }
 
  string inputType = ""; 

 // This set of lines are used to open and read the "config.in" file. 
 ///////////////////////////////////////////////////////////////////////    
 TEnv *params = new TEnv ("config_file");
 params->ReadFile ("config.in", kEnvChange);

 string inputROOTfile = params->GetValue ("inputFile", "inputROOTfile");
 string process       = params->GetValue ("process", "theprocess");      
 string dir_num       = params->GetValue ("dir_numerator", "fileNum");
 string dir_den       = params->GetValue ("dir_denominator", "fileDen");
 string histoName     = params->GetValue ("histoName", "massHisto"); 
 int marker           = params->GetValue ("markerStyle", 20);
 int setcolor         = params->GetValue ("colorStyle", 20);
 int rebin            = params->GetValue ("rebin", 1);
 int VariableBinSize  = params->GetValue ("variableBinSize", 0);
 int Nbins            = params->GetValue ("Nbins", 10);
 float binSize        = params->GetValue ("binSize", 100);
 //////////////////////////////////////////////////////////////////////

 TFile *file =  (TFile*) TFile::Open(inputROOTfile.c_str());
 file->cd(dir_num.c_str());

 // Create iterator for ROOT objects
 TIter nextKey_num (gDirectory->GetListOfKeys());
 TKey* key_num;
 TH1F* h_num;
 while ( (key_num = (TKey*) nextKey_num()) )  // Need double parentheses 
  {
     TObject* obj = key_num->ReadObj();
     TH1* histoObj = (TH1*) obj;
     // Get only 1D histograms from ROOT file and store in a TList
     if ( (histoObj->GetYaxis()->GetNbins() == 1) )
        {
          h_num = (TH1F*)histoObj;
          if (h_num->GetName() == histoName){h_num->Rebin(rebin); break;}
        }  
   }
 
  
 file->cd(dir_den.c_str());
 TIter nextKey_den (gDirectory->GetListOfKeys()); 
 TKey* key_den;
 TH1F* h_den;
 
 while ( (key_den = (TKey*) nextKey_den()) )  // Need double parentheses 
   {
    TObject* obj = key_num->ReadObj();
     TH1* histoObj = (TH1*) obj;
     if ( (histoObj->GetYaxis()->GetNbins() == 1) )
       {
          h_den = (TH1F*)histoObj;
          if (h_den->GetName() == histoName){h_den->Rebin(rebin); break;}
        }
     
  }

 // When you run out of statistics in the high mass or high pT tails
 // you might want to dump the events in falling in those regions, into 
 // one large bin. The lines below, show you how to do it. If you set the 
 // VariableBinSize to 1 in the config file, and make the plot for the 
 // di-lepton mass, you'll get an idea what this part of the code does.
  
 TH1F* H_bins_num;
 TH1F* H_bins_den;
 
 if (VariableBinSize==1){
   int h_Nbins = h_num->GetXaxis()->GetNbins();
   int arrayD = Nbins;
   float *BINS = new float [arrayD];

   for (int b = 0; b <= arrayD; b++ )
     {
       if (b < (arrayD-1))
	 {
	   BINS[b] = (float)b*(h_num->GetXaxis()->GetBinWidth(b));
	 }
       else {
	 BINS[b] = binSize;
       }
     }

   int bin_max = Nbins - 1;

   H_bins_num = new TH1F("h_NUM", "h_NUM", bin_max, BINS);
   H_bins_den = new TH1F("h_DEN", "h_DEN", bin_max, BINS);
   for (int b = 0; b <= (h_Nbins+1); b++ )
     {
       H_bins_num->SetBinContent(b, h_num->GetBinContent(b));
       H_bins_den->SetBinContent(b, h_den->GetBinContent(b));
     }
 } 

 TCanvas *c1 = new TCanvas("c1", "c1");
 c1->cd();

 if (VariableBinSize==1){ 
   TGraphAsymmErrors* gr1 = new TGraphAsymmErrors( H_bins_num, H_bins_den, "b(1,1) mode" ); 
   gr1->SetMarkerStyle(marker);
   gr1->SetMarkerColor(setcolor);
   gr1->Draw("AP");
 } else {
   TGraphAsymmErrors* gr1 = new TGraphAsymmErrors( h_num, h_den, "b(1,1) mode" );
   gr1->SetMarkerStyle(marker);
   gr1->SetMarkerColor(setcolor);
   gr1->Draw("AP");
 }
}

