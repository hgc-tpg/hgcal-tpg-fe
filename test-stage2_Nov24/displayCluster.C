/**********************************************************************
 Created on : 15/07/2025
 Purpose    : Display the clusters of different passes
 Author     : Indranil Das, Research Associate, Imperial
 Email      : indranil.das@cern.ch | indra.ehep@gmail.com
**********************************************************************/

int displayCluster(Long64_t refEvent = 1)
{
  
  int GetHexEdges(bool hextype, std::pair<float,float> xyCentre, float R, std::pair<float,float> []);//std::vector<std::pair<float,float>>& edglst);
  //const char *infile16 = "/home/indra/temp/multiple_passes/TPGS2Emu_tree_ntuples_16_0_vbf.root";
  //const char *infile16 = "/home/indra/temp/multiple_passes/TPGReco_tree_ntuples_45_0_minbias.root";
  //const char *infile16 = "/home/indra/temp/multiple_passes/TPGReco_tree_ntuples_16_0.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter42/minbias_PU200_2kFiles/TPGReco_tree_ntuples_45_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter42/minbias_PU200_2kFiles/TPGReco_tree_ntuples_30_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter42/minbias_PU200_2kFiles/TPGReco_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter44/minbias_PU200_2kFiles/TPGReco_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter45/minbias_PU200_2kFiles/TPGReco_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter47/singlePion_PU0_Ideal/TPGReco_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter48/singlePion_PU0_Ideal/TPGReco_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter48/singlePion_PU0_Ideal/TPGReco_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter49/singlePion_PU0/TPGReco_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter50/singlePion_PU0_Realistic/TPGReco_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/ResultBkg_iter51/vbfHInv_0PU/TPGReco_tree_ntuples_16_merged.root";
  const char *infile16 = "/home/hep/idas/stage2_emulation_results/Result_iter53/vbfHInv_0PU/TPGS2Emu_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/stage2_emulation_results/Result_iter54/singlePion_PU0_Ideal/TPGS2Emu_tree_ntuples_16_merged.root";
  //const char *infile16 = "/home/hep/idas/temp/TPGReco_tree_ntuples_16_0_vbf.root";
  //const char *infile16 = "/home/hep/idas/temp/TPGReco_tree_ntuples_16_0_vbf_test-distance1.root";
  //const char *infile16 = "/home/hep/idas/temp/TPGReco_tree_ntuples_45_0_vbf.root";
  
  std::unique_ptr<TFile> fin16(TFile::Open(Form("%s",infile16)));
  std::unique_ptr<TTree> tr16((TTree*)fin16->Get("TPGS2Emu"));
  //std::unique_ptr<TTree> tr16((TTree*)fin16->Get("TPG_Reco"));
  
  tr16->SetBranchStatus("*",0);
  
  std::vector<float>  *genjet_pt = 0 ;
  tr16->SetBranchStatus("genjet_pt",1);
  tr16->SetBranchAddress("genjet_pt" , &genjet_pt);
  
  std::vector<float>  *genjet_eta = 0 ;
  tr16->SetBranchStatus("genjet_eta",1);
  tr16->SetBranchAddress("genjet_eta" , &genjet_eta);
  
  // std::vector<float>  *genjet_phi = 0 ;
  // tr16->SetBranchStatus("genjet_phi",1);
  // tr16->SetBranchAddress("genjet_phi" , &genjet_phi);
  
  std::vector<float>  *clus_pt = 0 ;
  // tr16->SetBranchStatus("clus_pt",1);
  // tr16->SetBranchAddress("clus_pt" , &clus_pt);
  tr16->SetBranchStatus("clus_pt_corr2D",1);
  tr16->SetBranchAddress("clus_pt_corr2D" , &clus_pt);
  
  std::vector<float>  *clus_eta = 0 ;
  // tr16->SetBranchStatus("clus_local_eta",1);
  // tr16->SetBranchAddress("clus_local_eta" , &clus_eta);
  tr16->SetBranchStatus("clus_eta",1);
  tr16->SetBranchAddress("clus_eta" , &clus_eta);
  
  std::vector<float>  *clus_phi = 0 ;
  // tr16->SetBranchStatus("clus_local_phi",1);
  // tr16->SetBranchAddress("clus_local_phi" , &clus_phi);
  tr16->SetBranchStatus("clus_phi",1);
  tr16->SetBranchAddress("clus_phi" , &clus_phi);
  
  // std::vector<float>  *clus_zcm = 0 ;
  // tr16->SetBranchStatus("clus_zcm",1);
  // tr16->SetBranchAddress("clus_zcm" , &clus_zcm);
  
  std::vector<unsigned int>  *clus_pass = 0 ;
  tr16->SetBranchStatus("clus_pass",1);
  tr16->SetBranchAddress("clus_pass" , &clus_pass);  
  
  TGraph *grGen = new TGraph(0);
  TGraph *grClus0 = new TGraph(0);
  TGraph *grClus1 = new TGraph(0);
  TGraph *grClus2 = new TGraph(0);
  TGraph *grClus3 = new TGraph(0);
  TGraph *grClus4 = new TGraph(0);
  TGraph *grClus5 = new TGraph(0);
  //Long64_t refEvent = 1;
  
  std::cout << std::fixed << std::setprecision(4) << std::setw(10) << std::endl;
  std::cout << "Nofevents: " << tr16->GetEntries() << std::endl;
  for (Long64_t ievent = 0 ; ievent < tr16->GetEntries() ; ievent++ ) {
  //for (Long64_t ievent = 0 ; ievent < 10000 ; ievent++ ) {
  //for (Long64_t ievent = 0 ; ievent < 100000 ; ievent++ ) {
  //for (Long64_t ievent = 0 ; ievent < 10 ; ievent++ ) {
    tr16->GetEntry(ievent) ;
    if(ievent%10000==0) std::cout << "Processing event: " << ievent << ", nof clusters: " << clus_pt->size() << ", nofgenjets: "<< genjet_pt->size() << std::endl;
    double zpos_ref = 400.0;
    
    //std::cout << "ievent: " << ievent << ", nof genjets : " << genjet_pt->size() << std::endl;
    double genpt_negeta = -1.0;
    for(int igenj = 0; igenj < genjet_pt->size() ; igenj++){
      if(genjet_eta->at(igenj)<0.) genpt_negeta = genjet_pt->at(igenj);
    //   double zpos = (genjet_eta->at(igenj)>0)?zpos_ref:-zpos_ref;
    //   double genjet_x = zpos * cos(genjet_phi->at(igenj))/sinh(genjet_eta->at(igenj));
    //   double genjet_y = zpos * sin(genjet_phi->at(igenj))/sinh(genjet_eta->at(igenj));
    //   //if(genjet_eta->at(igenj)<0.){
    //   if(genjet_eta->at(igenj)<0. and refEvent==ievent){
    //   	std::cout << "ievent: " << ievent
    // 		  << std::fixed << std::setprecision(5) << std::setw(10)
    // 		  << ", genjet (x,y): (" << genjet_x << ", " << genjet_y << ") "
    // 		  <<", pt: " << genjet_pt->at(igenj)
    // 		  << std::defaultfloat
    // 		  << std::endl;
    // 	//grGen->SetPoint(ipoint,genjet_x,genjet_y);
    // 	grGen->AddPoint(genjet_x,genjet_y);
    //   }
    }
    
    for(int iclus = 0; iclus < clus_pt->size() ; iclus++){
      double zpos = (clus_eta->at(iclus)>0)?zpos_ref:-zpos_ref;
      double clus_et = clus_pt->at(iclus);
      // double clus_x = zpos * cos(clus_phi->at(iclus))/sinh(clus_eta->at(iclus));
      // double clus_y = zpos * sin(clus_phi->at(iclus))/sinh(clus_eta->at(iclus));      
      double clus_x = cos(clus_phi->at(iclus))/sinh(clus_eta->at(iclus));
      double clus_y = sin(clus_phi->at(iclus))/sinh(clus_eta->at(iclus));      
      //if(clus_eta->at(iclus)<0. and refEvent==ievent){
      if(clus_eta->at(iclus)<0.){
      //if(clus_eta->at(iclus)<0. and clus_et>200.){
      //if(clus_eta->at(iclus)<0. and clus_et<100. and clus_et>10.){
      //if(clus_eta->at(iclus)<0. and genpt_negeta>200. and genpt_negeta>0.){
      	// std::cout << "ievent: " << ievent
	// 	  << std::fixed << std::setprecision(4) << std::setw(8)
	// 	  << ", clus (x,y,z): (" << clus_x << ", " << clus_y << ", " << zpos << ") "
	// 	  <<", pt: " << clus_pt->at(iclus) << ", pass: " << clus_pass->at(iclus)
	// 	  << std::defaultfloat
	// 	  << std::endl;
	if(clus_pass->at(iclus)==0){
	  grClus0->AddPoint(clus_x,clus_y);
	}else if(clus_pass->at(iclus)==1){
	  grClus1->AddPoint(clus_x,clus_y);
	}else if(clus_pass->at(iclus)==2){
	  grClus2->AddPoint(clus_x,clus_y);
	}else if(clus_pass->at(iclus)==3){
	  grClus3->AddPoint(clus_x,clus_y);
	}else if(clus_pass->at(iclus)==4){
	  grClus4->AddPoint(clus_x,clus_y);
	}else if(clus_pass->at(iclus)==5){
	  grClus5->AddPoint(clus_x,clus_y);
	}
      }
    }
    //std::cout << std::defaultfloat << std::endl;
    
    // genjet_pt->clear();
    // genjet_eta->clear();
    // genjet_phi->clear();
    clus_pt->clear();
    clus_eta->clear();
    clus_phi->clear();
    //clus_zcm->clear();
    
  }//event loop
  
  double smin_x = 0.1;
  double smax_x = 0.15;
  double smin_y = 0.0;
  double smax_y = 0.05;
  // double smax_x = 0.6;
  // double smin_x = -0.6;
  // double smax_y = 0.6;
  // double smin_y = -0.6;
  TGraph *grDraw = new TGraph(0);
  grDraw->AddPoint(smin_x,smin_y);
  grDraw->AddPoint(smax_x,smax_y);
  grDraw->AddPoint(smin_x,smax_y);
  grDraw->AddPoint(smax_x,smin_y);  
  grDraw->SetTitle("");
  grDraw->GetXaxis()->SetTitle("x/z (cm)");
  grDraw->GetYaxis()->SetTitle("y/z (cm)");
  grDraw->GetXaxis()->SetTitleOffset(1.2);
  grDraw->GetYaxis()->SetTitleOffset(1.4);
  
  grGen->SetTitle("");
  grGen->GetXaxis()->SetTitle("x (cm)");
  grGen->GetYaxis()->SetTitle("y (cm)");
  grGen->GetXaxis()->SetTitleOffset(1.2);
  grGen->GetYaxis()->SetTitleOffset(1.4);
  grGen->SetMinimum(-200);
  grGen->SetMaximum(200);
  grGen->GetXaxis()->SetRange(-200,200);
  grGen->SetLineColor(kBlack);
  grGen->SetMarkerColor(kBlack);
  grGen->SetMarkerSize(0.8);
  grGen->SetMarkerStyle(kFullCircle);
  
  grClus0->SetLineColor(kRed);
  grClus0->SetMarkerColor(kRed);
  grClus0->SetMarkerStyle(7);

  TH2D *h2Pass0 = new TH2D("h2Pass0","h2Pass0", 25,smin_x,smax_x, 25,smin_y,smax_y);
  for(int ip=0;ip<grClus0->GetN();ip++){
    if(grClus0->GetPointX(ip)>smin_x and grClus0->GetPointX(ip)<smax_x and grClus0->GetPointY(ip)>smin_y and grClus0->GetPointY(ip)<smax_y)
      h2Pass0->Fill(grClus0->GetPointX(ip),grClus0->GetPointY(ip));
  }
  h2Pass0->GetXaxis()->SetTitle("x/z (cm)");
  h2Pass0->GetYaxis()->SetTitle("y/z (cm)");
  h2Pass0->GetXaxis()->SetTitleOffset(1.2);
  h2Pass0->GetYaxis()->SetTitleOffset(1.4);
  
  grClus1->SetLineColor(kBlue);
  grClus1->SetMarkerColor(kBlue);
  grClus1->SetMarkerStyle(7);
  
  grClus2->SetLineColor(kMagenta);
  grClus2->SetMarkerColor(kMagenta);
  grClus2->SetMarkerStyle(7);
  
  grClus3->SetLineColor(kGreen+2);
  grClus3->SetMarkerColor(kGreen+2);
  grClus3->SetMarkerStyle(7);
  
  grClus4->SetLineColor(kOrange+2);
  grClus4->SetMarkerColor(kOrange+2);
  grClus4->SetMarkerStyle(7);
  
  grClus5->SetLineColor(kTeal+2);
  grClus5->SetMarkerColor(kTeal+2);
  grClus5->SetMarkerStyle(7);
  
  std::pair<float,float> xycentre = std::make_pair(0,0);
  float R = 180;
  std::pair<float,float> edges[6];
  TLine *l1[6];
  GetHexEdges(0, xycentre, R, edges);
  for(int ip=0;ip<6;ip++){
    std::cout << "ip : " << ip << ", (x,y): ("<< edges[ip].first << ", " << edges[ip].second << ") " <<std::endl;
    int ipnext = (ip==5)?0:ip+1;
    l1[ip] = new TLine(edges[ip].first,edges[ip].second,edges[ipnext].first,edges[ipnext].second);
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetTickx(1);
  c1->SetTicky(1);
  grDraw->Draw("AP");
  //grGen->Draw("P");
  grClus0->Draw("P");
  //grClus1->Draw("P");
  //grClus2->Draw("P");
  //grClus3->Draw("P");
  //grClus4->Draw("P");
  //grClus5->Draw("P");
  //for(int ip=0;ip<6;ip++) l1[ip]->Draw();
  c1->Update();
  //c1->SaveAs("/home/indra/temp/tt.pdf");
  c1->SaveAs("tt.pdf");
  
  TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
  c2->SetGridx();
  c2->SetGridy();
  c2->SetTickx(1);
  c2->SetTicky(1);
  h2Pass0->Draw("colz box");
  c2->Update();
  c2->SaveAs("tt1.pdf");

  std::unique_ptr<TFile> fout(TFile::Open("output.root","recreate"));
  h2Pass0->Write();
  grClus0->Write();
  grClus1->Write();
  
  // grDraw->Draw("AP");
  // grClus0->Draw("P");
  // grClus1->Draw("P");
  // grClus2->Draw("P");
  // grClus3->Draw("P");
  // grClus4->Draw("P");
  // grClus5->Draw("P");


  return 0;
}

int GetHexEdges(bool hextype, std::pair<float,float> xyCentre, float R, std::pair<float,float> edges[6])//std::vector<std::pair<float,float>>& edglst)
{
  
  //hextype 1 for regular, 0 truncated type;
  float r = TMath::Cos(TMath::Pi()/6.)*R;
  float aby2 = TMath::Sin(TMath::Pi()/6.)*R;
  if(hextype){
    edges[0] = std::make_pair(xyCentre.first, xyCentre.second+R);
    edges[1] = std::make_pair(xyCentre.first+r, xyCentre.second+aby2);
    edges[2] = std::make_pair(xyCentre.first+r, xyCentre.second-aby2);
    edges[3] = std::make_pair(xyCentre.first, xyCentre.second-R);
    edges[4] = std::make_pair(xyCentre.first-r, xyCentre.second-aby2);
    edges[5] = std::make_pair(xyCentre.first-r, xyCentre.second+aby2);    
  }else{
    edges[0] = std::make_pair(xyCentre.first+aby2, xyCentre.second+r);
    edges[1] = std::make_pair(xyCentre.first+R, xyCentre.second);
    edges[2] = std::make_pair(xyCentre.first+aby2, xyCentre.second-r);
    edges[3] = std::make_pair(xyCentre.first-aby2, xyCentre.second-r);
    edges[4] = std::make_pair(xyCentre.first-R, xyCentre.second);
    edges[5] = std::make_pair(xyCentre.first-aby2, xyCentre.second+r);    
  }
  
  // edglst.clear();
  // edglst.push_back(edges[0]);
  
  return true;
}
