#include <dlUtility.h>

int temp_draw_njl()
{
  TCanvas* can = new TCanvas("","",2000,1000);
  TFile* myfile = TFile::Open("../output/hists/all22njl.root");

  TH2D* thehist = (TH2D*)myfile->Get("njl_rn_jetet");
  TH1D* myhist = (TH1D*)myfile->Get("njl_rn_15");
  TH1D* projx = myhist;//thehist->ProjectionX("projx",15,100);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  projx->SetMarkerStyle(20);
  projx->SetMarkerColor(kBlack);
  projx->SetLineColor(kBlack);
  projx->GetYaxis()->SetRangeUser(5e4,1e5);
  for(int i=0; i<projx->GetNbinsX()+2; ++i)
    {
      if(std::isnan(projx->GetBinContent(i)))
	{
	  cout << "found nan " << i << " " << projx->GetBinContent(i) << endl;
	  projx->SetBinContent(i,0);
	  projx->SetBinError(i,0);
	}
    }
  for(int i=-1; i<12; ++i)
    {
      int lower = 47000+500*i;
      projx->GetXaxis()->SetRangeUser(lower,lower+500);
      if(i==-1) projx->GetXaxis()->SetRangeUser(47200,54000);
      projx->Draw("PE");
      sphenixtext(0.65,0.96);
      sqrt_s_text(0.65,0.92);
      float minet = 20;
      antikt_text(0.4,0.3,0.92);	  
      et_cut_text(minet,0.015,0.96);                                                                             
      dijet_cut_text(0.3,0.96);
      gPad->SaveAs(("../output/plots/temp"+to_string(i)+".png").c_str());
    }

  return 0;
}
