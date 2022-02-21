

{

  gStyle->SetNumberContours(99);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleStyle(2000);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetTitleYOffset(1.2);
  gStyle->SetTitleW(1);
  gStyle->SetTitleY(0.88);
  gStyle->SetTitleX(0.58);
  gStyle->SetTitleH(0.08);

  gStyle->SetLabelSize(0.06);
  gStyle->SetTitleXSize(0.06);
  gStyle->SetTitleYSize(0.06);
  gStyle->SetNdivisions(505);
  // gStyle->SetLegendTextSize(0.04);

  TFile *f1=new TFile("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/Cascade_Star_M_Sigma_Star_M_Sigma_Width_Output_041120_01.root");
  auto *hBE_dsss_Cascade_Star_Sigma_Star_Sigma  = (TH1D*) f1->Get("hBinding_dsss_M_Phs_WW_FF_noq");

  TFile *f2=new TFile("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/Cascade_Star_M_Sigma_Star_M_Lambda_Width_Output_041120_01.root");
  auto *hBE_dsss_Cascade_Star_Sigma_Star_Lambda  = (TH1D*) f2->Get("hBinding_dsss_M_Phs_WW_FF_noq");

  TFile *f3=new TFile("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/Omega_M_Delta_M_Width_Output_041120_01.root");
  auto *hBE_dsss_Omega_Delta  = (TH1D*) f3->Get("hBinding_dsss_M_Phs_WW_FF_noq");

  TFile *f4=new TFile("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/Cascade_M_Sigma_M_Width_Output_041120_01.root");
  auto *hBE_dsss_Cascade_Sigma  = (TH1D*) f4->Get("hBinding_dsss_M_Phs_WW_FF_noq");

  // Adjusting for factors
  hBE_dsss_Omega_Delta->Scale(0.5);
  hBE_dsss_Cascade_Star_Sigma_Star_Sigma->Scale(0.5);
  hBE_dsss_Cascade_Star_Sigma_Star_Lambda->Scale(0.5);

  // Printing out the widths
  Double_t dsss_Cascade_Star_Sigma_Star_lambda_width = hBE_dsss_Cascade_Star_Sigma_Star_Lambda->GetBinContent(hBE_dsss_Cascade_Star_Sigma_Star_Lambda->FindBin(-0.084));

  Double_t dsss_Cascade_Star_Sigma_Star_Sigma_width = hBE_dsss_Cascade_Star_Sigma_Star_Sigma->GetBinContent(hBE_dsss_Cascade_Star_Sigma_Star_Sigma->FindBin(-0.084));

  Double_t d_sss_Omega_Delta_width = hBE_dsss_Omega_Delta->GetBinContent(hBE_dsss_Omega_Delta->FindBin(-0.084));

  Double_t d_sss_Cascade_Sigma_width = hBE_dsss_Cascade_Sigma->GetBinContent(hBE_dsss_Cascade_Sigma->FindBin(-0.084);



  cout<<"d_{sss}->#Xi*#Sigma* (#Sigma*>#Lambda#pi):"<<dsss_Cascade_Star_Sigma_Star_lambda_width<<endl;
  cout<<"d_{sss}->#Xi*#Sigma* (#Sigma*>#Sigma#pi):"<<dsss_Cascade_Star_Sigma_Star_Sigma_width<<endl;
  cout<<"d_{sss}->#Omega#Delta: "<<d_sss_Omega_Delta_width<<endl;
  cout<<"d_{sss}->#Xi#Sigma: "<<d_sss_Cascade_Sigma_width<<endl;

  auto* hBE_dsss_total = new TH1D("hBE_dsss_total","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBE_dsss_total->GetXaxis()->SetNdivisions(511);
  hBE_dsss_total->GetXaxis()->SetLabelSize(0.05);
  hBE_dsss_total->GetXaxis()->SetTitleSize(0.05);
  hBE_dsss_total->GetYaxis()->SetNdivisions(510);
  hBE_dsss_total->GetYaxis()->SetLabelSize(0.05);
  hBE_dsss_total->GetYaxis()->SetTitleSize(0.05);

  auto* hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma = new TH1D("hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetXaxis()->SetNdivisions(511);
  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetXaxis()->SetLabelSize(0.05);
  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetXaxis()->SetTitleSize(0.05);
  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetYaxis()->SetNdivisions(510);
  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetYaxis()->SetLabelSize(0.05);
  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetYaxis()->SetTitleSize(0.05);

  auto* hBE_dsss_oct_plus_Omega_Delta = new TH1D("hBE_dsss_oct_plus_Omega_Delta","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBE_dsss_oct_plus_Omega_Delta->GetXaxis()->SetNdivisions(511);
  hBE_dsss_oct_plus_Omega_Delta->GetXaxis()->SetLabelSize(0.05);
  hBE_dsss_oct_plus_Omega_Delta->GetXaxis()->SetTitleSize(0.05);
  hBE_dsss_oct_plus_Omega_Delta->GetYaxis()->SetNdivisions(510);
  hBE_dsss_oct_plus_Omega_Delta->GetYaxis()->SetLabelSize(0.05);
  hBE_dsss_oct_plus_Omega_Delta->GetYaxis()->SetTitleSize(0.05);

  auto* hBR_dsss_total = new TH1D("hBR_dsss_total","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBR_dsss_total->GetXaxis()->SetNdivisions(511);
  hBR_dsss_total->GetXaxis()->SetLabelSize(0.05);
  hBR_dsss_total->GetXaxis()->SetTitleSize(0.05);
  hBR_dsss_total->GetYaxis()->SetNdivisions(510);
  hBR_dsss_total->GetYaxis()->SetLabelSize(0.05);
  hBR_dsss_total->GetYaxis()->SetTitleSize(0.05);

  auto* hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma = new TH1D("hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetXaxis()->SetNdivisions(511);
  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetXaxis()->SetLabelSize(0.05);
  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetXaxis()->SetTitleSize(0.05);
  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetYaxis()->SetNdivisions(510);
  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetYaxis()->SetLabelSize(0.05);
  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->GetYaxis()->SetTitleSize(0.05);

  auto* hBR_dsss_oct_plus_Omega_Delta = new TH1D("hBR_dsss_oct_plus_Omega_Delta","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBR_dsss_oct_plus_Omega_Delta->GetXaxis()->SetNdivisions(511);
  hBR_dsss_oct_plus_Omega_Delta->GetXaxis()->SetLabelSize(0.05);
  hBR_dsss_oct_plus_Omega_Delta->GetXaxis()->SetTitleSize(0.05);
  hBR_dsss_oct_plus_Omega_Delta->GetYaxis()->SetNdivisions(510);
  hBR_dsss_oct_plus_Omega_Delta->GetYaxis()->SetLabelSize(0.05);
  hBR_dsss_oct_plus_Omega_Delta->GetYaxis()->SetTitleSize(0.05);

  auto* hBR_dsss_Cascade_Star_Sigma_Star_Lambda = new TH1D("hBR_dsss_Cascade_Star_Sigma_Star_Lambda","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBR_dsss_Cascade_Star_Sigma_Star_Lambda->GetXaxis()->SetNdivisions(511);
  hBR_dsss_Cascade_Star_Sigma_Star_Lambda->GetXaxis()->SetLabelSize(0.05);
  hBR_dsss_Cascade_Star_Sigma_Star_Lambda->GetXaxis()->SetTitleSize(0.05);
  hBR_dsss_Cascade_Star_Sigma_Star_Lambda->GetYaxis()->SetNdivisions(510);
  hBR_dsss_Cascade_Star_Sigma_Star_Lambda->GetYaxis()->SetLabelSize(0.05);
  hBR_dsss_Cascade_Star_Sigma_Star_Lambda->GetYaxis()->SetTitleSize(0.05);

  auto* hBR_dsss_Cascade_Star_Sigma_Star_Sigma = new TH1D("hBR_dsss_Cascade_Star_Sigma_Star_Sigma","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBR_dsss_Cascade_Star_Sigma_Star_Sigma->GetXaxis()->SetNdivisions(511);
  hBR_dsss_Cascade_Star_Sigma_Star_Sigma->GetXaxis()->SetLabelSize(0.05);
  hBR_dsss_Cascade_Star_Sigma_Star_Sigma->GetXaxis()->SetTitleSize(0.05);
  hBR_dsss_Cascade_Star_Sigma_Star_Sigma->GetYaxis()->SetNdivisions(510);
  hBR_dsss_Cascade_Star_Sigma_Star_Sigma->GetYaxis()->SetLabelSize(0.05);
  hBR_dsss_Cascade_Star_Sigma_Star_Sigma->GetYaxis()->SetTitleSize(0.05);

  auto* hBR_dsss_Omega_Delta = new TH1D("hBR_dsss_Omega_Delta","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBR_dsss_Omega_Delta->GetXaxis()->SetNdivisions(511);
  hBR_dsss_Omega_Delta->GetXaxis()->SetLabelSize(0.05);
  hBR_dsss_Omega_Delta->GetXaxis()->SetTitleSize(0.05);
  hBR_dsss_Omega_Delta->GetYaxis()->SetNdivisions(510);
  hBR_dsss_Omega_Delta->GetYaxis()->SetLabelSize(0.05);
  hBR_dsss_Omega_Delta->GetYaxis()->SetTitleSize(0.05);

  auto* hBR_dsss_Cascade_Sigma = new TH1D("hBR_dsss_Cascade_Sigma","#Lambda = 0.16 GeV;B [GeV];#Gamma [GeV]",5000,-1.00062,4);
  hBR_dsss_Cascade_Sigma->GetXaxis()->SetNdivisions(511);
  hBR_dsss_Cascade_Sigma->GetXaxis()->SetLabelSize(0.05);
  hBR_dsss_Cascade_Sigma->GetXaxis()->SetTitleSize(0.05);
  hBR_dsss_Cascade_Sigma->GetYaxis()->SetNdivisions(510);
  hBR_dsss_Cascade_Sigma->GetYaxis()->SetLabelSize(0.05);
  hBR_dsss_Cascade_Sigma->GetYaxis()->SetTitleSize(0.05);

  //widths
  *hBE_dsss_total = *hBE_dsss_Omega_Delta + *hBE_dsss_Cascade_Sigma + *hBE_dsss_Cascade_Star_Sigma_Star_Sigma + *hBE_dsss_Cascade_Star_Sigma_Star_Lambda;

  //Branching ratios
  *hBR_dsss_Cascade_Sigma = *hBE_dsss_Cascade_Sigma;
  *hBR_dsss_Cascade_Sigma->Divide(hBE_dsss_Cascade_Sigma, hBE_dsss_total);

  *hBR_dsss_Omega_Delta = *hBE_dsss_Omega_Delta;
  hBR_dsss_Omega_Delta->Divide(hBE_dsss_Omega_Delta, hBE_dsss_total);

  *hBR_dsss_Cascade_Star_Sigma_Star_Sigma = *hBE_dsss_Cascade_Star_Sigma_Star_Sigma;
  hBR_dsss_Cascade_Star_Sigma_Star_Sigma->Divide(hBE_dsss_Cascade_Star_Sigma_Star_Sigma, hBE_dsss_total);

  *hBR_dsss_Cascade_Star_Sigma_Star_Lambda = *hBE_dsss_Cascade_Star_Sigma_Star_Lambda;
  hBR_dsss_Cascade_Star_Sigma_Star_Lambda->Divide(hBE_dsss_Cascade_Star_Sigma_Star_Lambda, hBE_dsss_total);

  *hBR_dsss_oct_plus_Omega_Delta = *hBR_dsss_Cascade_Sigma;
  *hBR_dsss_oct_plus_Omega_Delta = *hBR_dsss_Omega_Delta + *hBR_dsss_Cascade_Sigma;

  *hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma = *hBR_dsss_Cascade_Star_Sigma_Star_Sigma;
  *hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma = *hBR_dsss_Cascade_Star_Sigma_Star_Sigma + *hBR_dsss_oct_plus_Omega_Delta;

  *hBR_dsss_total = *hBR_dsss_Cascade_Star_Sigma_Star_Lambda;
  *hBR_dsss_total = *hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma + *hBR_dsss_Cascade_Star_Sigma_Star_Lambda;


  //set width colors
  hBE_dsss_Cascade_Sigma->SetLineColor(kViolet);
  hBE_dsss_Cascade_Sigma->SetFillColor(kViolet);

  hBE_dsss_oct_plus_Omega_Delta->SetLineColor(kBlue);
  hBE_dsss_oct_plus_Omega_Delta->SetFillColor(kBlue);

  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->SetLineColor(kGreen);
  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->SetFillColor(kGreen);

  hBE_dsss_total->SetLineColor(kRed);
  hBE_dsss_total->SetFillColor(kRed);

  //set branching ratio colors
  hBR_dsss_Cascade_Sigma->SetLineColor(kViolet);
  hBR_dsss_Cascade_Sigma->SetFillColor(kViolet);

  hBR_dsss_oct_plus_Omega_Delta->SetLineColor(kBlue);
  hBR_dsss_oct_plus_Omega_Delta->SetFillColor(kBlue);

  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->SetLineColor(kGreen);
  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->SetFillColor(kGreen);

  hBR_dsss_total->SetLineColor(kRed);
  hBR_dsss_total->SetFillColor(kRed);

  // Printing out the branching ratios
  Double_t dsss_Cascade_Star_Sigma_Star_lambda_BR = hBR_dsss_Cascade_Star_Sigma_Star_Lambda->GetBinContent(917);

  Double_t dsss_Cascade_Star_Sigma_Star_Sigma_BR = hBR_dsss_Cascade_Star_Sigma_Star_Sigma->GetBinContent(917);

  Double_t d_sss_Omega_Delta_BR = hBR_dsss_Omega_Delta->GetBinContent(917);

  Double_t d_sss_Cascade_Sigma_BR = hBR_dsss_Cascade_Sigma->GetBinContent(917);

  cout<<dsss_Cascade_Star_Sigma_Star_lambda_BR<<endl;
  cout<<dsss_Cascade_Star_Sigma_Star_Sigma_BR<<endl;
  cout<<d_sss_Omega_Delta_BR<<endl;
  cout<<d_sss_Cascade_Sigma_BR<<endl;



  TLine *line = new TLine(-0.084,0.001,-0.084,0.15);
  TLine *line2 = new TLine(-0.084,0.001,-0.084,1);

  line->SetLineColor(1);
  line->SetLineWidth(2);
  line->SetLineStyle(9);

  line2->SetLineColor(1);
  line2->SetLineWidth(2);
  line2->SetLineStyle(9);

  TBox *box = new TBox(-0.285,0.66,-0.205,0.72);
  box->SetLineColor(kBlack);
  box->SetLineWidth(2);
  box->SetFillColor(0);
  box->SetFillStyle(0);

  hBE_dsss_total->GetXaxis()->SetRangeUser(-0.3, 0.3);
  //SetNdivisions should be 505!
  hBE_dsss_total->GetXaxis()->SetNdivisions(505);
  hBE_dsss_total->GetXaxis()->SetTitle("B [GeV]");
  hBE_dsss_total->GetYaxis()->SetTitle("#Gamma [GeV]");

  hBE_dsss_total->SetTitle("d_{sss}#rightarrow#Xi*#Sigma* (#Sigma*#rightarrow#Lambda#pi)");
  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->SetTitle("d_{sss}#rightarrow#Xi*#Sigma* (#Sigma*#rightarrow#Sigma#pi)");
  hBE_dsss_oct_plus_Omega_Delta->SetTitle("d_{sss}#rightarrow#Omega#Delta");
  hBE_dsss_Cascade_Sigma->SetTitle("d_{sss}#rightarrow#Xi#Sigma");

  TCanvas *c1 = new TCanvas("c1","stacked hists",10,10,1000,1000);
  c1->cd();

  hBE_dsss_total->GetYaxis()->SetTitleOffset(1.1);

  hBE_dsss_total->DrawCopy("HIST logy");
  hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->DrawCopy("HIST same");
  hBE_dsss_oct_plus_Omega_Delta->DrawCopy("HIST same");
  hBE_dsss_Cascade_Sigma->DrawCopy("HIST same");

  line->Draw("same");

  // c1->BuildLegend(0.15,0.7,0.7,0.9,"","f");

  auto legend1 =  new TLegend(0.15,0.625,0.7,0.9);
  legend1->AddEntry(hBE_dsss_total,"d_{sss}#rightarrow#Xi*#Sigma* (#Sigma*#rightarrow#Lambda#pi)","f");
  legend1->AddEntry(hBE_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma,"d_{sss}#rightarrow#Xi*#Sigma* (#Sigma*#rightarrow#Sigma#pi)","f");
  legend1->AddEntry(hBE_dsss_oct_plus_Omega_Delta,"d_{sss}#rightarrow#Omega#Delta","f");
  legend1->AddEntry(hBE_dsss_Cascade_Sigma,"d_{sss}#rightarrow#Xi#Sigma","f");
  legend1->SetFillColor(0);
  legend1->SetFillStyle(0);
  legend1->SetLineWidth(0);
  legend1->SetTextSize(0.07);
  legend1->Draw();
  gPad->RedrawAxis();

  hBR_dsss_total->GetXaxis()->SetRangeUser(-0.3, 0.3);
  hBR_dsss_total->GetXaxis()->SetNdivisions(505);
  hBR_dsss_total->GetYaxis()->SetRangeUser(0, 1);
  hBR_dsss_total->GetYaxis()->SetNdivisions(505);
  hBR_dsss_total->GetXaxis()->SetTitle("B [GeV]");
  hBR_dsss_total->GetYaxis()->SetTitle("Branching Ratio");

  hBR_dsss_total->SetTitle("d_{sss}#rightarrow#Xi*#Sigma* (#Sigma*#rightarrow#Lambda#pi)");
  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->SetTitle("d_{sss}#rightarrow#Xi*#Sigma* (#Sigma*#rightarrow#Sigma#pi)");
  hBR_dsss_oct_plus_Omega_Delta->SetTitle("d_{sss}#rightarrow#Omega#Delta");
  hBR_dsss_Cascade_Sigma->SetTitle("d_{sss}#rightarrow#Xi#Sigma");

  TCanvas *c2 = new TCanvas("c2","Canvas",10,10,1000,1000);
  c2->Divide(1,1);
  c2->cd();
  hBR_dsss_total->GetYaxis()->SetTitleOffset(1.5);


  hBR_dsss_total->SetMinimum(0);
  hBR_dsss_total->SetMaximum(1);
  hBR_dsss_total->DrawCopy("HIST");
  hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma->DrawCopy("HIST same");
  hBR_dsss_oct_plus_Omega_Delta->DrawCopy("HIST same");
  hBR_dsss_Cascade_Sigma->DrawCopy("HIST same");


  auto legend2 =  new TLegend(0.15,0.625,0.7,0.9);
  legend2->AddEntry(hBR_dsss_total,"d_{sss}#rightarrow#Xi*#Sigma* (#Sigma*#rightarrow#Lambda#pi)","f");
  legend2->AddEntry(hBR_dsss_oct_plus_Omega_Delta_Sigma_Star_Sigma,"d_{sss}#rightarrow#Xi*#Sigma* (#Sigma*#rightarrow#Sigma#pi)","f");
  legend2->AddEntry(hBR_dsss_oct_plus_Omega_Delta,"d_{sss}#rightarrow#Omega#Delta","f");
  legend2->AddEntry(hBR_dsss_Cascade_Sigma,"d_{sss}#rightarrow#Xi#Sigma","f");
  legend2->SetFillColor(0);
  legend2->SetFillStyle(0);
  legend2->SetLineWidth(0);
  legend2->SetTextSize(0.07);
  legend2->SetTextColor(0);
  legend2->Draw();

  line2->Draw("same");
  gPad->RedrawAxis();


  box->Draw("same");


}
