int temp_draw_njl()
{
  TFile* myfile = TFile::Open("../output/hists/all22njl.root");

  TH2D* thehist = (TH2D*)myfile->Get("njl_rn_jetet");

  TH1D* projx = thehist->ProjectionX("projx",20,100);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetOptStat(0);
  projx->SetMarkerStyle(20);
  projx->SetMarkerColor(kBlack);
  projx->SetLineColor(kBlack);
  projx->GetYaxis()->SetRangeUser(1e4,2e4);
  for(int i=0; i<12; ++i)
    {
      int lower = 47000+500*i;
      projx->GetXaxis()->SetRangeUser(lower,lower+500);
      projx->Draw("PE");
      gPad->SaveAs(("../output/plots/temp"+to_string(i)+".png").c_str());
    }

  return 0;
}
