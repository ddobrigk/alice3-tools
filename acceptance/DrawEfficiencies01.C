void DrawEfficiencies01(){ 
  Int_t fields[] = {5,10,15,20}; //in kGauss
  TString configs[] = {"def", "10layer", "12layer", "13layer"};
  const int nfields = sizeof(fields)/sizeof(Int_t); // number of field configurations
  const int nconfigs = 3; // number of layer configurations
  
  TFile *file[nfields][nconfigs];
  
  TH1D *hNum[nfields][nconfigs], *hGen[nfields][nconfigs];
  
  for(Int_t ic=0; ic<nconfigs; ic++){
    cout<<"Processing config: "<<configs[ic].Data()<<endl;
    for(Int_t fi=0; fi<nfields; fi++){
      cout<<"Processing field: "<<fields[fi]<<" kGauss"<<endl;
      
      file[fi][ic] = new TFile(Form("output_%iU_%s.root",fields[fi],configs[ic].Data()), "READ");
      hNum[fi][ic] = (TH1D*) file[fi][ic] -> Get("hReconstructedXiCC");
      hGen[fi][ic] = (TH1D*) file[fi][ic] -> Get("hGeneratedXiCC");
      hNum[fi][ic] -> SetName(Form("hNum_%i_%i", ic, fi));
      hGen[fi][ic] -> SetName(Form("hGen_%i_%i", ic, fi));
      hNum[fi][ic]->Sumw2();
      hGen[fi][ic]->Sumw2();
      hNum[fi][ic]->Divide(hGen[fi][ic]);
    }
  }
  
  const int drawField = 3;
  const int refConfig = 0;
  
  TCanvas *c1 = new TCanvas("c1", "", 800, 600);
  c1->Divide(1,2);
  c1->cd(1)->SetPad(0,0.5,1,1.0);
  c1->cd(2)->SetPad(0,0.0,1,0.5);
  
  c1->cd(1)->SetTicks(1,1);
  c1->cd(1)->SetLeftMargin(0.115);
  c1->cd(1)->SetBottomMargin(0.002);
  c1->cd(1)->SetRightMargin(0.03);
  c1->cd(1)->SetTopMargin(0.03);
  c1->cd(2)->SetTicks(1,1);
  c1->cd(2)->SetLeftMargin(0.115);
  c1->cd(2)->SetBottomMargin(0.172);
  c1->cd(2)->SetRightMargin(0.03);
  c1->cd(2)->SetTopMargin(0.002);
  
  c1->cd(1);
  
  Int_t markerStyles[] = {20, 21, 24, 25, 28, 33};
  Int_t markerColors[] = {kBlack, kGreen+1, kBlue+1, kRed+1, kViolet};
  
  for(Int_t ic=0; ic<nconfigs; ic++){
    hNum[drawField][ic]->SetMarkerStyle(markerStyles[ic]);
    hNum[drawField][ic]->SetMarkerColor(markerColors[ic]);
    hNum[drawField][ic]->SetLineColor(markerColors[ic]);
    hNum[drawField][ic]->SetMarkerSize(0.65);
  }
  
  gStyle->SetOptStat(0);
  hNum[drawField][0]->SetTitle("");
  hNum[drawField][0]->GetXaxis()->SetRangeUser(0,15);
  hNum[drawField][0]->GetYaxis()->SetRangeUser(0,.45);
  hNum[drawField][0]->GetYaxis()->SetTitle("Max acceptance");
  hNum[drawField][0]->GetYaxis()->SetTitleOffset(0.8);
  hNum[drawField][0]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  hNum[drawField][0]->GetYaxis()->SetTitleSize(0.07);
  hNum[drawField][0]->GetXaxis()->SetTitleSize(0.07);
  hNum[drawField][0]->GetYaxis()->SetLabelSize(0.052);
  hNum[drawField][0]->GetXaxis()->SetLabelSize(0.052);
  hNum[drawField][0]->Draw("");
  
  TH1D *hRatio[nconfigs];
  
  for(Int_t ic=0; ic<nconfigs; ic++){
    hNum[drawField][ic]->Draw("same");
    hRatio[ic] = (TH1D*) hNum[drawField][ic] -> Clone ( Form("hRatio_%i", ic) );
  }
  
  for(Int_t ic=0; ic<nconfigs; ic++){
    hRatio[ic] -> Divide( hNum[drawField][refConfig] );
  }

  
  TH1D *hDummy = (TH1D*) hNum[drawField][0]->Clone("hDummy");
  
  c1->cd(2);
  
  hDummy->GetYaxis()->SetTitle("Ratio to 11-layer");
  hDummy->GetYaxis()->SetTitleOffset(0.8);
  hDummy->GetYaxis()->SetRangeUser(0,1.6);
  hDummy->Reset();
  hDummy->Draw();
  for(Int_t ic=0; ic<nconfigs; ic++){
    if( ic != refConfig)
      hRatio[ic] -> Draw("same");
  }
  
  TH1D *hLayerArrangement[nconfigs];
  
//  for(Int_t ic=0; ic<nconfigs; ic++){
//    TString layerString = "";
//    hLayerArrangement[ic] = (TH1D*) file[ic]->Get("hLayerArrangement");
//    hLayerArrangement[ic] -> SetName( Form( "hLayerArrangement_%i", ic ) );
//
//    for(Int_t il=0; il<hLayerArrangement->GetNbinsX()+1; il++){
//      cout<<"Read in layer "<<il<<" at "<<hLayerArrangement->GetBinLowEdge(il+1)<<" cm"<<endl;
//      layerString += Form("%.1f", hLayerArrangement->GetBinLowEdge(il+1));
//      if( il != hLayerArrangement->GetNbinsX()) layerString += "-";
//    }
//  }
//
  TLegend *leg = new TLegend(0.712, 0.395, 0.929, 0.871);
  
  //TString configs[] = {"def", "10layer", "12layer", "13layer"};
  TString layerDescription[] = {"11-layer", "10-layer", "11-lay, red. OT rad.", "13-layer"};
  
  c1->cd(1);
  leg->SetBorderSize(0);
  leg->SetHeader(Form("Field: %i kGauss", fields[drawField]));
  for(Int_t fi=0; fi<nconfigs; fi++){
    cout<<"Processing field: "<<fields[fi]<<" kGauss"<<endl;
    leg->AddEntry(hNum[drawField][fi], layerDescription[fi], "LP");
  }
  leg->Draw();
  
  
//  for(Int_t fi=1; fi<nfields; fi++){
//    hNum[fi]->Draw("same");
//  }
//
//  TLatex *lat = new TLatex();
//  lat->SetNDC();
//  lat->SetTextAlign(32);
//  TString layerString = "";
//
//  TH1D *hLayerArrangement = (TH1D*) file[0]->Get("hLayerArrangement");
//
//  for(Int_t il=0; il<hLayerArrangement->GetNbinsX()+1; il++){
//    cout<<"Read in layer "<<il<<" at "<<hLayerArrangement->GetBinLowEdge(il+1)<<" cm"<<endl;
//    layerString += Form("%.1f", hLayerArrangement->GetBinLowEdge(il+1));
//    if( il != hLayerArrangement->GetNbinsX()) layerString += "-";
//
//  }
//
//  lat->SetTextSize(0.027);
//  lat->DrawLatex(0.93,0.915,"#Xi_{cc}^{++} efficiency with layer arrangement (cm):");
//
//  int minHits = hLayerArrangement->GetBinContent(1);
//  lat->DrawLatex(0.93,0.825,Form("N_{layers}: %i, N_{min}^{hits} = %i", hLayerArrangement->GetNbinsX()+1, minHits));
//  lat->SetTextFont(42);
//  lat->DrawLatex(0.93,0.87, layerString.Data());
//
//  TLegend *leg = new TLegend(0.75, 0.56, 0.93, 0.80);
//  leg->SetBorderSize(0);
//  for(Int_t fi=0; fi<nfields; fi++){
//    cout<<"Processing field: "<<fields[fi]<<" kGauss"<<endl;
//    leg->AddEntry(hNum[fi], Form("+%i kGauss",fields[fi]), "LP");
//  }
//  leg->Draw();
//
  c1->SaveAs("efficiency.pdf");
  
} 
