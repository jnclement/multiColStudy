#include <dlUtility.h>
#include <TGraphErrors.h>

int set_axis_range_based_on_axis(TH3D* h, int axis, int from, int to)
{
  if(axis==0) h->GetXaxis()->SetRange(from,to);
  else if(axis==1) h->GetYaxis()->SetRange(from,to);
  else if(axis==2) h->GetZaxis()->SetRange(from,to);
  else return 1;
  return 0;
}

int get_bin_number_for_axis(TH3D* h, int axis, double value, int getlo)
{
  int bin = -1;
  if(axis==0) bin = h->GetXaxis()->FindBin(value);
  else if(axis==1) bin = h->GetYaxis()->FindBin(value);
  else if(axis==2) bin = h->GetZaxis()->FindBin(value);
  else return bin;
  bin += (getlo?0:1);
  return bin;
}
vector<TObject*> make_projections(TH3D* h, int axis, double from, double to)
{

  string projdir;
  vector<TObject*> outhists = {};

  
  if(axis==0) projdir = "yz";
  else if(axis==1) projdir = "xz";
  else if(axis==2) projdir = "xy";
  else return {};

  TH2D* proj1;
  TH2D* proj2;
  
  if(to > from)
    {
      set_axis_range_based_on_axis(h, axis, get_bin_number_for_axis(h, axis, from, 1), get_bin_number_for_axis(h, axis, to, 0));
      proj1 = (TH2D*)(h->Project3D(("proj1_"+projdir).c_str()));
    }
  else
    {
      set_axis_range_based_on_axis(h, axis, 0, get_bin_number_for_axis(h, axis, to, 0));
      proj1 = (TH2D*)(h->Project3D(("proj1_"+projdir).c_str()));
      //proj1->Draw("COLZ");
      set_axis_range_based_on_axis(h, axis, get_bin_number_for_axis(h, axis, from, 1), 9999999);
      proj2 = (TH2D*)(h->Project3D(("proj1_"+projdir).c_str()));
      proj1->Add(proj2);
    }

  proj1->GetZaxis()->SetTitle("Counts");

  
  outhists.push_back(proj1->Clone());
  ((TH2D*)outhists.at(0))->GetXaxis()->SetRange(1,((TH2D*)outhists.at(0))->GetNbinsX());
  ((TH2D*)outhists.at(0))->GetYaxis()->SetRange(1,((TH2D*)outhists.at(0))->GetNbinsY());
  outhists.push_back(((TH2D*)outhists.at(0))->ProjectionX());  
  outhists.push_back(((TH2D*)outhists.at(0))->ProjectionY());
  
  return outhists;
}

vector<TGraphErrors*> get_th2d_mean_tgraph(TH2D* h)
{
  TGraphErrors* means = new TGraphErrors();
  TGraphErrors* jes = new TGraphErrors();
  TGraphErrors* jer = new TGraphErrors();
  TGraphErrors* getrms = new TGraphErrors();
  getrms->SetMarkerStyle(20);
  getrms->SetMarkerColor(kRed);
  getrms->SetLineColor(kRed);
  getrms->SetMarkerSize(1);
  
  means->SetMarkerStyle(20);
  means->SetMarkerColor(kBlack);
  means->SetLineColor(kBlack);
  means->SetMarkerSize(1.5);
  means->SetLineWidth(2);
  
  int nx = h->GetNbinsX();
  int npoints = 0;
  //TCanvas* newcan = new TCanvas("newcan","",1000,1000);
  for(int i=1; i<=nx; ++i)
    {
      TH1D* projy = h->ProjectionY("_py",i,i);
      
      TF1* gaus = new TF1("gaus","gaus",0.3,projy->GetXaxis()->GetXmax());//projy->GetMean()-projy->GetStdDev(),projy->GetMean()+2*projy->GetStdDev());//
      projy->Fit(gaus,"QRIN");
      //projy->GetYaxis()->SetRangeUser(1e-12,1);
      //projy->Draw("PE");
      //gPad->SetLogy();
      //gPad->SaveAs(("../output/plots/"+std::string(h->GetName())+"_"+to_string(i)+".png").c_str());
      //gPad->SetLogy(0);
      double x = h->GetXaxis()->GetBinCenter(i);
      double mean = gaus->GetParameter(1);
      double err = gaus->GetParameter(2);
      if(mean == 0 || err == 0) continue;
      cout << x << " " <<  mean << " " << err << " " << gaus->GetParError(1) << " " << gaus->GetParError(2) << endl;
      jes->AddPoint(x,mean);
      jes->SetPointError(npoints,0,gaus->GetParError(1));
      jer->AddPoint(x,err/mean);
      jer->SetPointError(npoints,0,gaus->GetParError(2)/mean);
      means->AddPoint(x,mean);
      means->SetPointError(npoints,0,err);
      getrms->AddPoint(x,projy->GetMean());
      getrms->SetPointError(npoints,0,projy->GetStdDev());
      ++npoints;
      //means->AddPointError(x,mean,0,err);
      
      //delete gaus;
      //delete projy;
    }
  means->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  jes->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  jer->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  getrms->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  return {means,jes,jer,getrms};
}

int format_th1d(TH1D* hist, int color)
{
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.03);

  hist->SetMarkerSize(1.5);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  hist->SetLineWidth(2);
  hist->SetMarkerStyle(20);

  return 0;
}

int format_th2d(TH2D* hist)
{
  hist->GetXaxis()->SetLabelSize(0.03);
  hist->GetXaxis()->SetTitleSize(0.03);
  hist->GetYaxis()->SetLabelSize(0.03);
  hist->GetYaxis()->SetTitleSize(0.03);
  hist->GetZaxis()->SetLabelSize(0.03);
  hist->GetZaxis()->SetTitleSize(0.03);
  return 0;
}

int alltext(int hasz, string etcuttext = "", string zcuttext = "")
{
  sphenixtext(0.65,0.97);
  sqrt_s_text(0.65,0.93);
  antikt_text(0.4,0.3,0.93);
  //dijet_cut_text(0.6,0.96);
  drawText(etcuttext.c_str(),0.05,0.97);
  drawText(zcuttext.c_str(),0.05,0.93);
  drawText("MB,Jet10,Jet20,Jet30 PYTHIA",0.3,0.87);
  string hasztext = "";
  if(hasz == 0) hasztext = "z_{vtx} assumed 0";
  else if(hasz == 1) hasztext = "Reco z_{vtx} used";
  drawText(hasztext.c_str(),0.3,0.97);
  return 0;
}

TH1D* GraphToHist(const TGraphErrors* graph, const TString& name = "hFromGraph") {
    const int n = graph->GetN();
    if (n == 0) return nullptr;

    // Allocate bin edges based on the x-points
    std::vector<double> edges;
    edges.reserve(n + 1);

    // Compute bin edges between points (midpoints), extrapolate edges for first/last
    for (int i = 0; i < n; ++i) {
        double x, y;
        graph->GetPoint(i, x, y);
        if (i == 0) {
            double xNext;
            graph->GetPoint(i + 1, xNext, y);
            double width = xNext - x;
            edges.push_back(x - width / 2.0);
        } else {
            double xPrev;
            graph->GetPoint(i - 1, xPrev, y);
            double width = x - xPrev;
            edges.push_back(xPrev + width / 2.0);
        }
        if (i == n - 1) {
            double xPrev;
            graph->GetPoint(i - 1, xPrev, y);
            double width = x - xPrev;
            edges.push_back(x + width / 2.0);
        }
    }

    // Convert vector to ROOT-style array
    TH1D* h = new TH1D(name, graph->GetTitle(), n, &edges[0]);

    // Fill bin contents and errors
    for (int i = 0; i < n; ++i) {
        double x, y;
        graph->GetPoint(i, x, y);
        double ey = graph->GetErrorY(i);
        h->SetBinContent(i + 1, y);
        h->SetBinError(i + 1, ey);
    }

    return h;
}

TGraphErrors* getRatioGraphBinomial(TGraphErrors* gNum, TGraphErrors* gDen)
{
  int n = gNum->GetN();
  TGraphErrors* gRatio = new TGraphErrors();
  
  for (int i = 0; i < n; ++i) {
    double x1, y1, x2, y2;
    gNum->GetPoint(i, x1, y1);
    gDen->GetPoint(i, x2, y2);
    
    if (x1 != x2) {
      std::cerr << "Warning: X points do not match at i = " << i << std::endl;
      continue;
    }
    
    if (y2 == 0) continue;  // avoid divide by zero
    
    double ratio = y1 / y2;
    double binomError = sqrt(ratio * (1 - ratio) / y2);
    
    gRatio->AddPoint(x1, ratio);
    gRatio->SetPointError(i, 0, binomError);  // x error = 0
  }
  
  gRatio->SetMarkerStyle(20);
  gRatio->SetMarkerColor(kBlack);
  gRatio->SetLineColor(kBlack);
  gRatio->SetTitle("Ratio with Binomial Errors");
  
  return gRatio;
}

TGraphErrors* DivideGraphsMatchingX(TGraphErrors* gNum, TGraphErrors* gDen, bool useBinomialError = false) {
    TGraphErrors* gRatio = new TGraphErrors();
    gRatio->GetXaxis()->SetTitle(gNum->GetXaxis()->GetTitle());
    int iRatio = 0;

    for (int i = 0; i < gNum->GetN(); ++i) {
        double x1, y1;
        gNum->GetPoint(i, x1, y1);
        double ey1 = gNum->GetErrorY(i);

        // Try to find matching x in denominator graph
        bool foundMatch = false;
        for (int j = 0; j < gDen->GetN(); ++j) {
            double x2, y2;
            gDen->GetPoint(j, x2, y2);

            if (std::abs(x1 - x2) < 1e-8) {  // Matching x with tolerance
                if (y2 == 0) break;  // Skip division by zero

                double ey2 = gDen->GetErrorY(j);
                double r = y1 / y2;
                double er;

                if (useBinomialError) {
                    er = sqrt(r * (1 - r) / y2);  // binomial
                } else {
                    // propagated error
                    er = r * sqrt((ey1 / y1) * (ey1 / y1) + (ey2 / y2) * (ey2 / y2));
                }

                gRatio->SetPoint(iRatio, x1, r);
                gRatio->SetPointError(iRatio, 0.0, er);
                ++iRatio;
                foundMatch = true;
                break;
            }
        }

        if (!foundMatch) {
            std::cerr << "No matching x for numerator point x = " << x1 << std::endl;
        }
    }

    return gRatio;
}

int overlay_w_ratio_tgraph(TGraphErrors* wz, TGraphErrors* nz, TCanvas* can, TLegend* leg, string etcuttext, string zcuttext, string histtype, string xy)
{
  can->cd(1);
  TGraphErrors* ratio = DivideGraphsMatchingX(nz,wz);
  ratio->GetYaxis()->SetTitle("Ratio no-z / with-z");
  ratio->GetXaxis()->SetTitle(wz->GetXaxis()->GetTitle());
  ratio->GetXaxis()->SetTitleSize(0.075);
  ratio->GetYaxis()->SetTitleSize(0.075);
  ratio->GetXaxis()->SetLabelSize(0.075);
  ratio->GetYaxis()->SetLabelSize(0.075);
  //ratio->Divide(wz,nz,1,1,"B");

  wz->SetMarkerStyle(20);
  wz->SetMarkerColor(kBlack);
  wz->SetLineColor(kBlack);
  wz->SetMarkerSize(1.5);
  nz->SetMarkerStyle(20);
  nz->SetMarkerColor(kRed);
  nz->SetLineColor(kRed);
  nz->SetMarkerSize(1.5);

  ratio->SetMarkerStyle(20);
  ratio->SetMarkerColor(kRed);
  ratio->SetLineColor(kRed);
  ratio->SetMarkerSize(1.5);
  
  double ymax = wz->GetHistogram()->GetMaximum();
  double nzmax = nz->GetHistogram()->GetMaximum();
  if(nzmax > ymax) ymax = nzmax;
  ymax *= 1.5;

  //cout << wz->GetMaximum() << endl;
  //cout << ratio->GetXaxis()->GetXmin()<< " " << ratio->GetXaxis()->GetXmax() << endl;

  wz->Draw("APE");
  nz->Draw("SAME PE");
  wz->GetHistogram()->GetYaxis()->SetRangeUser(0,1.5);
  nz->GetHistogram()->GetYaxis()->SetRangeUser(0,1.5);
  wz->GetHistogram()->GetXaxis()->SetRangeUser(wz->GetHistogram()->GetXaxis()->GetXmin(),wz->GetHistogram()->GetXaxis()->GetXmax());
  nz->GetHistogram()->GetXaxis()->SetRangeUser(wz->GetHistogram()->GetXaxis()->GetXmin(),wz->GetHistogram()->GetXaxis()->GetXmax());  

  gPad->Update();
  
  can->cd(2);
  ratio->Draw("APE");
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->GetHistogram()->GetXaxis()->SetRangeUser(wz->GetHistogram()->GetXaxis()->GetXmin(),wz->GetHistogram()->GetXaxis()->GetXmax());

  gPad->Update();
  can->cd(0);
  alltext(-1, etcuttext, zcuttext);
  can->cd(1);
  leg->Draw();
  TLine* oneline = new TLine(wz->GetHistogram()->GetXaxis()->GetXmin(),1,wz->GetHistogram()->GetXaxis()->GetXmax(),1);
  can->cd(2);
  oneline->Draw();

  can->SaveAs(("../output/plots/"+histtype + "_"+xy+".png").c_str());

  //delete oneline;
  //delete ratio;
  
  return 0;
}

int overlay_w_ratio_th1d(TH1D* wz, TH1D* nz, TCanvas* can, TLegend* leg, string etcuttext, string zcuttext, string histtype, string xy)
{
  can->cd(1);
  TH1D* ratio = (TH1D*)wz->Clone();
  cout << ratio->GetNbinsX() << " " << ratio->GetBinLowEdge(201) << endl;
  
  ratio->Divide(wz,nz,1,1,"B");
  ratio->GetYaxis()->SetTitle("Ratio Reco z / Zero z");
  
  double ymax = wz->GetMaximum();
  double nzmax = nz->GetMaximum();
  if(nzmax > ymax) ymax = nzmax;
  ymax *= 1.5; 
  wz->GetYaxis()->SetRangeUser(0,ymax);
  nz->GetYaxis()->SetRangeUser(0,ymax);

  wz->Draw("PE");
  nz->Draw("SAME PE");
  can->cd(2);
  ratio->GetYaxis()->SetRangeUser(0,2);
  ratio->Draw("PE");
  can->cd(0);
  alltext(-1, etcuttext, zcuttext);
  can->cd(1);
  leg->Draw();
  TLine* oneline = new TLine(wz->GetXaxis()->GetXmin(),1,nz->GetXaxis()->GetXmax(),1);
  can->cd(2);
  oneline->Draw();

  can->SaveAs(("../output/plots/"+histtype + "_"+xy+".png").c_str());

  //delete oneline;

  return 0;
}

int FindLastBinAboveX(TH2D* h, double threshold) {
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    for (int i = nx; i >= 1; --i) {  // loop backward
        double sum = 0.0;
        for (int j = 1; j <= ny; ++j) {
            sum += h->GetBinContent(i, j);
        }
        if (sum > threshold) return i;
    }
    return -1;  // none found
}

int FindFirstBinAboveX(TH2D* h, double threshold) {
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    for (int i = 1; i <= nx; ++i) {  // loop backward
        double sum = 0.0;
        for (int j = 1; j <= ny; ++j) {
            sum += h->GetBinContent(i, j);
        }
        if (sum > threshold) return i;
    }
    return -1;  // none found
}

int FindLastBinAboveY(TH2D* h, double threshold) {
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    for (int j = ny; j >= 1; --j) {
        double sum = 0.0;
        for (int i = 1; i <= nx; ++i) {
            sum += h->GetBinContent(i, j);
        }
        if (sum > threshold) return j;
    }
    return -1;
}

int FindFirstBinAboveY(TH2D* h, double threshold) {
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    for (int j = 1; j <= ny; ++j) {
        double sum = 0.0;
        for (int i = 1; i <= nx; ++i) {
            sum += h->GetBinContent(i, j);
        }
        if (sum > threshold) return j;
    }
    return -1;
}

int draw_all(string histtype, vector<TObject*> wz, vector<TObject*> nz, string etcuttext = "", string zcuttext = "")
{

  TH2D* wz2 = (TH2D*)wz.at(0);
  TH1D* wzx = (TH1D*)wz.at(1);
  TH1D* wzy = (TH1D*)wz.at(2);

  TH2D* nz2 = (TH2D*)nz.at(0);
  TH1D* nzx = (TH1D*)nz.at(1);
  TH1D* nzy = (TH1D*)nz.at(2);

  //cout << wzy->GetNbinsX() << " " << wzy->GetBinLowEdge(200) << " " << nzy->GetNbinsX() << " " << nzy->GetBinLowEdge(200) << endl;
  
  format_th2d(wz2);
  format_th2d(nz2);
  format_th1d(wzx,kBlack);
  format_th1d(wzy,kBlack);
  format_th1d(nzx,kRed);
  format_th1d(nzy,kRed);

  vector<TGraphErrors*> wzgs = get_th2d_mean_tgraph(wz2);
  vector<TGraphErrors*> nzgs = get_th2d_mean_tgraph(nz2);
  TGraphErrors* wzg = wzgs.at(0);
  TGraphErrors* nzg = nzgs.at(0);
  //wzg->Draw("PE");
  TGraphErrors* jsw = wzgs.at(1);
  TGraphErrors* jsn = nzgs.at(1);
  TGraphErrors* jrw = wzgs.at(2);
  TGraphErrors* jrn = nzgs.at(2);

  for(int i=0; i<jsw->GetN(); ++i)
    {
      double x, y;
      jsn->GetPoint(i, x, y);
      cout << x << " " << y << endl;
    }

  TCanvas* can = new TCanvas("","",1000,1000);
  can->cd();
  gPad->SetTopMargin(0.15);
  gPad->SetBottomMargin(0.15);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.2);

  wz2->GetXaxis()->SetRangeUser(wz2->GetXaxis()->GetBinLowEdge(FindFirstBinAboveX(wz2,0)),wz2->GetXaxis()->GetBinLowEdge(FindLastBinAboveX(wz2,0)+1));
  wz2->GetYaxis()->SetRangeUser(wz2->GetYaxis()->GetBinLowEdge(FindFirstBinAboveY(wz2,0)),wz2->GetYaxis()->GetBinLowEdge(FindLastBinAboveY(wz2,0)+1));

  nz2->GetXaxis()->SetRangeUser(nz2->GetXaxis()->GetBinLowEdge(FindFirstBinAboveX(nz2,0)),nz2->GetXaxis()->GetBinLowEdge(FindLastBinAboveX(nz2,0)+1));
  nz2->GetYaxis()->SetRangeUser(nz2->GetYaxis()->GetBinLowEdge(FindFirstBinAboveY(nz2,0)),nz2->GetYaxis()->GetBinLowEdge(FindLastBinAboveY(nz2,0)+1));

  gPad->SetLogz();
  wz2->GetZaxis()->SetRangeUser(1e-12,wz2->GetMaximum()*2);
  wz2->Draw("COLZ");
  wzg->Draw("SAME PE");
  //wzgs.at(3)->Draw("SAME PE");
  alltext(1,etcuttext,zcuttext);
  can->SaveAs(("../output/plots/"+histtype + "_wz.png").c_str());
  nz2->GetZaxis()->SetRangeUser(1e-12,nz2->GetMaximum()*2);
  nz2->Draw("COLZ");
  nzg->Draw("SAME PE");
  //nzgs.at(3)->Draw("SAME PE");
  alltext(0,etcuttext,zcuttext);
  can->SaveAs(("../output/plots/"+histtype + "_nz.png").c_str());
  gPad->SetLogz(0);
  //delete can;
  can = new TCanvas("","",1000,1500);
  ratioPanelCanvas(can,0.3);
  can->cd(1);
  gPad->SetTopMargin(0.15);
  TLegend* leg = new TLegend(0.3,0.6,0.8,0.75);
  leg->AddEntry(wzx,"Reco z_{vtx} used for jets","p");
  leg->AddEntry(nzx,"z_{vtx} Assumed 0","p");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  overlay_w_ratio_th1d(wzx, nzx, can, leg, etcuttext, zcuttext, histtype, "x");
  overlay_w_ratio_th1d(wzy, nzy, can, leg, etcuttext, zcuttext, histtype, "y");
  overlay_w_ratio_tgraph(jsw, jsn, can, leg, etcuttext, zcuttext, histtype, "jes");
  overlay_w_ratio_tgraph(jrw, jrn, can, leg, etcuttext, zcuttext, histtype, "jer");
  //overlay_w_ratio_tgraph(wzg, nzg, can, leg, etcuttext, zcuttext, histtype, "mean");

  //delete leg;
  //delete can;
  //delete wzg;
  //delete nzg;
  
  return 0;
}

int get_and_draw(TH3D* hw, TH3D* hn, int axis, double from, double to, string etcuttext = "", string zcuttext = "", string histbase = "")
{

  string histtype = histbase+"_proj";
  if(axis == 0) histtype += "x";
  else if(axis==1) histtype += "y";
  else if(axis==2) histtype += "z";
  else return 1;

  histtype += "_";
  histtype += to_string(from);
  histtype += "_";
  histtype += to_string(to);
  
  vector<TObject*> wz = make_projections(hw, axis, from, to);
  vector<TObject*> nz = make_projections(hn, axis, from, to);
  
  draw_all(histtype,wz,nz,etcuttext,zcuttext);

  return 0;
}

int plot_th3d(string filename)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);

  TFile* file = TFile::Open(filename.c_str());
  if(!file) return 1;
  TH3D* hw = (TH3D*)file->Get("h3_resp_pT_zvtx");
  TH3D* hn = (TH3D*)file->Get("h3_resp_pT_zvtx_noz");

  TH3D* hwe = (TH3D*)file->Get("h3_resp_E_zvtx");
  TH3D* hne = (TH3D*)file->Get("h3_resp_E_zvtx_noz");

  hw->Rebin(4,1,1);
  hn->Rebin(4,1,1);
  hwe->Rebin(4,1,1);
  hwn->Rebin(4,1,1);
  
  TH2D* corre = (TH2D*)file->Get("noz_recoz_corret");

  corre->Draw("COLZ");
  gPad->SetLogz();
  gPad->Update();
  gPad->SaveAs("../output/plots/corre.png");
  gPad->SetLogz(0);
  if(!hw || !hn) return 2;
  hw->Scale(1./hw->Integral());
  hn->Scale(1./hn->Integral());
  get_and_draw(hw, hn, 2, 60, -60, "","|z_{vtx}^{truth}|>60 cm","gt60");
  get_and_draw(hw, hn, 2, -30, 30, "","|z_{vtx}^{truth}|<30 cm","lt30");
  get_and_draw(hwe, hne, 2, 60, -60, "","|z_{vtx}^{truth}|>60 cm","gr60e");
  get_and_draw(hwe, hne, 2, -30, 30, "","|z_{vtx}^{truth}|<30 cm","lt30e");  
  
  return 0;
}
