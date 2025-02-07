double massL0 = 1.11568;
double xmin =1.09;
double xmax =1.15;
double background_gaus(double *x, double *par) {
    return par[0]*(1. + par[1]*x[0] + par[2]*TMath::Power(x[0],2) + par[3]*TMath::Power(x[0],3) + par[4]*TMath::Power(x[0],4)+ par[5]*TMath::Power(x[0],5))+par[6]*exp(-0.5*(TMath::Power((x[0]-par[7])/par[8],2)));
}
double background(double *x, double *par) {
    return par[0]*(1. + par[1]*x[0] + par[2]*TMath::Power(x[0],2) + par[3]*TMath::Power(x[0],3) + par[4]*TMath::Power(x[0],4)+ par[5]*TMath::Power(x[0],5));
}

double gaussian(double *x, double *par) {
    return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2)));
}
double double_gaussian(double *x, double *par) {
    return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2)))+par[3]*exp(-0.5*(TMath::Power((x[0]-par[4])/par[5],2)));
}

double fitting_function(double *x, double *par) {
    return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2))) + par[3]*(1. + par[4]*x[0] + par[5]*TMath::Power(x[0],2) + par[6]*TMath::Power(x[0],3) + par[7]*TMath::Power(x[0],4)+par[8]*TMath::Power(x[0],5));
}
double fitting_function_iter2(double *x, double *par) {
    return par[0]*exp(-0.5*(TMath::Power((x[0]-par[1])/par[2],2))) + par[3]*(1. + par[4]*x[0] + par[5]*TMath::Power(x[0],2) + par[6]*TMath::Power(x[0],3) + par[7]*TMath::Power(x[0],4)+ par[8]*TMath::Power(x[0],5))+par[9]*exp(-0.5*(TMath::Power((x[0]-par[10])/par[11],2)));
}
double func_non(double *x, double *par) {
    return par[0]*x[0]+par[1];
}
void v1_run8_sel_pt(){

	std::vector<char*> int_fitting_vect_pt={new char[100],new char[100],new char[100],new char[100],new char[100],new char[100],new char[100]};
  std::vector<TF1*> fit_vector_pt={new TF1(),new TF1(),new TF1(),new TF1(),new TF1(),new TF1(),new TF1()};
std::vector<char*> int_fitting_vect_eta={new char[100],new char[100],new char[100],new char[100],new char[100]};
std::vector<TF1*> fit_vector_eta={new TF1(),new TF1(),new TF1(),new TF1(),new TF1()};
//TFile *inFile_data = new TFile("/home/valeriy/bmn_lamb/JAM_lambda_qa_eff_050225_bckg_peak_rejected_3_less_minv_bins.root");
TFile *inFile_data = new TFile("/home/valeriy/bmn_lamb/JAM_lambda_qa_eff_no_weight_040225_protpion_bckg_study_8.root");
	TH1 *h_inv_mass_Psi;
	TH1 *h_inv_mass_sig;
	TH1 *h_inv_mass_bckg;
TH3D* invmass_pt_eta = (TH3D*)inFile_data->Get("h1_mass_pt_eta_basic_var");
TH3D* invmass_pt_eta_sig = (TH3D*)inFile_data->Get("h1_mass_pt_eta_selected_signal");
TH3D* invmass_pt_eta_bckg = (TH3D*)inFile_data->Get("h1_mass_pt_eta_selected_background");
TCanvas* c1 =new TCanvas("c1","A Simple Graph with Error bars",-1,50,2500,3500);
        c1->GetFrame()->SetBorderSize(12);
        c1->Divide(5,7,0,0,0);
	for( int i=1;i<8;i++){
        for(int j=1;j<6;j++){
		double gappts[8]={0,0.2,0.4,0.6,0.8,1.0,1.2,2.0};
double gapys[6]={-0.2,0,0.3,0.5,0.8,1.0};
               auto h_inv_mass_Psi2=invmass_pt_eta->ProjectionX(Form("invmass_pt%d_eta%d",i,j),i,i,j,j);
		h_inv_mass_Psi2->SetTitle(Form("%1.2f<pt<%1.2f %1.1f<y<%1.1f",gappts[i-1],gappts[i],gapys[j-1],gapys[j]));
		std::cout<<h_inv_mass_Psi2->GetTitle()<<std::endl;
		c1->cd(j+(i-1)*5);
                h_inv_mass_Psi2->Draw();
        }}
	c1->SaveAs("/home/valeriy/bmn_lamb/hardcut/lambda_invmass_map_NP_new_bins.png");
for( int i=1;i<8;i++){
	TF1 *fitting_fnc_full, *fitting_fnc_iter2, *fitting_bckg_full, *fitting_sig_full;
    char *int_fitting_fnc_full = new char[100];
    char *int_backFcn_full = new char[100];
    char *int_signalFcn_full = new char[100];
    char *int_fitting_fnc_iter2 = new char[100];
			h_inv_mass_Psi=invmass_pt_eta->ProjectionX(Form("invmass_pt%d_eta",i),i,i,3,5);
			h_inv_mass_sig=invmass_pt_eta_sig->ProjectionX(Form("invmass_pt%d_eta_sig",i),i,i,3,5);
			h_inv_mass_bckg=invmass_pt_eta_bckg->ProjectionX(Form("invmass_pt%d_eta_bckg",i),i,i,3,5);

			fitting_fnc_full = new TF1(int_fitting_fnc_full, fitting_function, xmin, xmax, 9);
        fitting_fnc_full->SetNpx(2000);
        fitting_fnc_full->SetLineWidth(4);
        fitting_fnc_full->SetLineColor(4);
        fitting_fnc_full->SetParNames("Strength", "Mean", "Sigma", "pol1", "pol2", "pol3", "pol4","pol5");
        int npar = fitting_fnc_full->GetNumberFreeParameters();
        for (int ip = 0; ip < npar; ++ip) fitting_fnc_full->SetParameter(ip, 0);
        fitting_fnc_full->SetParameter(0,h_inv_mass_Psi->GetMaximum());
fitting_fnc_full->SetParameter(1,massL0);
          fitting_fnc_full->SetParameter(2,0.002);

TCanvas *c2 = new TCanvas("c2","c2",-1,10,2500,500);
     h_inv_mass_Psi->Draw();
        h_inv_mass_Psi->Fit(int_fitting_fnc_full, "", "", xmin, xmax);
        double mass = fitting_fnc_full->GetParameter(npar - 8);
        double sigma = TMath::Abs(fitting_fnc_full->GetParameter(npar - 7));

        fitting_fnc_iter2 = new TF1("int_fitting_fnc_iter2",fitting_function_iter2,xmin,xmax,12);
        fitting_fnc_iter2->SetNpx(2000);
        fitting_fnc_iter2->SetLineWidth(4);
        fitting_fnc_iter2->SetLineColor(kGreen);
        fitting_fnc_iter2->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3");
        int npar2 = fitting_fnc_iter2->GetNumberFreeParameters();
        for (int ip = 0; ip < npar2; ++ip) fitting_fnc_iter2->SetParameter(ip,0);

        fitting_fnc_iter2->SetParameter(0,fitting_fnc_full->GetParameter(0));
        fitting_fnc_iter2->FixParameter(1,mass);
        fitting_fnc_iter2->SetParameter(2,sigma*1.5);
        fitting_fnc_iter2->SetParameter(9,fitting_fnc_full->GetParameter(0));
        fitting_fnc_iter2->FixParameter(10,mass);
        fitting_fnc_iter2->SetParameter(11,sigma);
        fitting_fnc_iter2->SetParLimits(2,sigma*1.5,sigma*2.5);
        fitting_fnc_iter2->SetParLimits(11,sigma/1.5,sigma*1.5);
        fitting_fnc_iter2->SetParLimits(9,0,fitting_fnc_full->GetParameter(0));
        h_inv_mass_Psi->Fit("int_fitting_fnc_iter2", "", "same", xmin, xmax);

	fitting_bckg_full = new TF1("int_backFcn_full",background,xmin,xmax,6);
fitting_sig_full = new TF1("int_signalFcn_full", double_gaussian,xmin,xmax,6);
        fitting_bckg_full->FixParameter(0,fitting_fnc_iter2->GetParameter(3));
        fitting_bckg_full->FixParameter(1,fitting_fnc_iter2->GetParameter(4));
        fitting_bckg_full->FixParameter(2,fitting_fnc_iter2->GetParameter(5));
        fitting_bckg_full->FixParameter(3,fitting_fnc_iter2->GetParameter(6));
        fitting_bckg_full->FixParameter(4,fitting_fnc_iter2->GetParameter(7));
        fitting_bckg_full->FixParameter(5,fitting_fnc_iter2->GetParameter(8));
//	fitting_bckg_full->FixParameter(6,fitting_fnc_iter2->GetParameter(0));
 //       fitting_bckg_full->FixParameter(7,fitting_fnc_iter2->GetParameter(1));
   //     fitting_bckg_full->FixParameter(8,fitting_fnc_iter2->GetParameter(2));


        h_inv_mass_Psi->Fit("int_backFcn_full","","same",xmin,xmax);
        fitting_sig_full->FixParameter(0,fitting_fnc_iter2->GetParameter(0));
        fitting_sig_full->FixParameter(1,fitting_fnc_iter2->GetParameter(1));
        fitting_sig_full->FixParameter(2,fitting_fnc_iter2->GetParameter(2));
        fitting_sig_full->FixParameter(3,fitting_fnc_iter2->GetParameter(9));
        fitting_sig_full->FixParameter(4,fitting_fnc_iter2->GetParameter(10));
        fitting_sig_full->FixParameter(5,fitting_fnc_iter2->GetParameter(11));
        h_inv_mass_Psi->Fit("int_signalFcn_full","","same",xmin,xmax);

	auto fit_signal_formula= [fitting_sig_full,fitting_bckg_full, fitting_fnc_iter2](double *x, double *par){
                return par[0]*(fitting_sig_full->Eval(x[0])/fitting_fnc_iter2->Eval(x[0]))+(par[1]+par[2]*x[0])*(fitting_bckg_full->Eval(x[0])/fitting_fnc_iter2->Eval(x[0]));

                                };
	fit_vector_pt.at(i-1)=new TF1(int_fitting_vect_pt.at(i-1), fit_signal_formula, 1.09, 1.15, 3);

}
for( int i=1;i<6;i++){
        TF1 *fitting_fnc_full, *fitting_fnc_iter2, *fitting_bckg_full, *fitting_sig_full;
    char *int_fitting_fnc_full = new char[100];
    char *int_backFcn_full = new char[100];
    char *int_signalFcn_full = new char[100];
    char *int_fitting_fnc_iter2 = new char[100];
                        h_inv_mass_Psi=invmass_pt_eta->ProjectionX(Form("invmass_pt_eta%d",i),3,6,i,i);
			h_inv_mass_sig=invmass_pt_eta_sig->ProjectionX(Form("invmass_pt_eta%d_sig",i),3,6,i,i);
			h_inv_mass_bckg=invmass_pt_eta_bckg->ProjectionX(Form("invmass_pt_eta%d_bckg",i),3,6,i,i);

                        fitting_fnc_full = new TF1(int_fitting_fnc_full, fitting_function, xmin, xmax, 9);
        fitting_fnc_full->SetNpx(2000);
        fitting_fnc_full->SetLineWidth(4);
        fitting_fnc_full->SetLineColor(kGreen);
        fitting_fnc_full->SetParNames("Strength", "Mean", "Sigma", "pol1", "pol2", "pol3", "pol4","pol5");
        int npar = fitting_fnc_full->GetNumberFreeParameters();
        for (int ip = 0; ip < npar; ++ip) fitting_fnc_full->SetParameter(ip, 0);
        fitting_fnc_full->SetParameter(0,h_inv_mass_Psi->GetMaximum());
fitting_fnc_full->SetParameter(1,massL0);
          fitting_fnc_full->SetParameter(2,0.002);

TCanvas *c2 = new TCanvas("c2","c2",-1,10,800,600);
     h_inv_mass_Psi->Draw();
        h_inv_mass_Psi->Fit(int_fitting_fnc_full, "", "", xmin, xmax);
        double mass = fitting_fnc_full->GetParameter(npar - 8);
        double sigma = TMath::Abs(fitting_fnc_full->GetParameter(npar - 7));

        fitting_fnc_iter2 = new TF1("int_fitting_fnc_iter2",fitting_function_iter2,xmin,xmax,12);
        fitting_fnc_iter2->SetNpx(2000);
        fitting_fnc_iter2->SetLineWidth(4);
        fitting_fnc_iter2->SetLineColor(kMagenta);
        fitting_fnc_iter2->SetParNames("Strength","Mean","Sigma","pol1","pol2","pol3");
        int npar2 = fitting_fnc_iter2->GetNumberFreeParameters();
        for (int ip = 0; ip < npar2; ++ip) fitting_fnc_iter2->SetParameter(ip,0);

        fitting_fnc_iter2->SetParameter(0,fitting_fnc_full->GetParameter(0));
        fitting_fnc_iter2->FixParameter(1,mass);
        fitting_fnc_iter2->SetParameter(2,sigma*1.5);
        fitting_fnc_iter2->SetParameter(9,fitting_fnc_full->GetParameter(0));
        fitting_fnc_iter2->FixParameter(10,mass);
        fitting_fnc_iter2->SetParameter(11,sigma);
        fitting_fnc_iter2->SetParLimits(2,sigma*1.5,sigma*2.5);
        fitting_fnc_iter2->SetParLimits(11,sigma/1.5,sigma*1.5);
        fitting_fnc_iter2->SetParLimits(9,0,fitting_fnc_full->GetParameter(0));
        h_inv_mass_Psi->Fit("int_fitting_fnc_iter2", "", "same", xmin, xmax);

        fitting_bckg_full = new TF1("int_backFcn_full",background,xmin,xmax,6);
	fitting_bckg_full->SetLineWidth(4);
        fitting_bckg_full->SetLineColor(kBlue);
fitting_sig_full = new TF1("int_signalFcn_full", double_gaussian,xmin,xmax,6);
        fitting_sig_full->SetLineWidth(4);
	fitting_sig_full->SetLineColor(kRed);
        fitting_bckg_full->FixParameter(0,fitting_fnc_iter2->GetParameter(3));
        fitting_bckg_full->FixParameter(1,fitting_fnc_iter2->GetParameter(4));
	fitting_bckg_full->FixParameter(2,fitting_fnc_iter2->GetParameter(5));
        fitting_bckg_full->FixParameter(3,fitting_fnc_iter2->GetParameter(6));
        fitting_bckg_full->FixParameter(4,fitting_fnc_iter2->GetParameter(7));
        fitting_bckg_full->FixParameter(5,fitting_fnc_iter2->GetParameter(8));
//	fitting_bckg_full->FixParameter(6,fitting_fnc_iter2->GetParameter(0));
  //      fitting_bckg_full->FixParameter(7,fitting_fnc_iter2->GetParameter(1));
   //     fitting_bckg_full->FixParameter(8,fitting_fnc_iter2->GetParameter(2));


        h_inv_mass_Psi->Fit("int_backFcn_full","","same",xmin,xmax);
        fitting_sig_full->FixParameter(0,fitting_fnc_iter2->GetParameter(0));
        fitting_sig_full->FixParameter(1,fitting_fnc_iter2->GetParameter(1));
        fitting_sig_full->FixParameter(2,fitting_fnc_iter2->GetParameter(2));
        fitting_sig_full->FixParameter(3,fitting_fnc_iter2->GetParameter(9));
        fitting_sig_full->FixParameter(4,fitting_fnc_iter2->GetParameter(10));
        fitting_sig_full->FixParameter(5,fitting_fnc_iter2->GetParameter(11));
        h_inv_mass_Psi->Fit("int_signalFcn_full","","same",xmin,xmax);


	fitting_fnc_full->Draw("same");
	fitting_fnc_iter2->Draw("same");
	fitting_bckg_full->Draw("same");

	

        auto fit_signal_formula= [fitting_sig_full,fitting_bckg_full, fitting_fnc_iter2](double *x, double *par){
                return par[0]*(fitting_sig_full->Eval(x[0])/fitting_fnc_iter2->Eval(x[0]))+(par[1]+par[2]*x[0])*(fitting_bckg_full->Eval(x[0])/fitting_fnc_iter2->Eval(x[0]));

                                };

        fit_vector_eta.at(i-1)=new TF1(int_fitting_vect_eta.at(i-1), fit_signal_formula, 1.09, 1.15, 3);
	 auto yield_signal_formula= [fitting_sig_full,fitting_bckg_full, fitting_fnc_iter2](double *x, double *par){
                return (fitting_sig_full->Eval(x[0])/fitting_fnc_iter2->Eval(x[0]));
                                };
	 auto yield_bckg_formula= [fitting_sig_full,fitting_bckg_full, fitting_fnc_iter2](double *x, double *par){
                return (fitting_bckg_full->Eval(x[0])/fitting_fnc_iter2->Eval(x[0]));
                                };

	 char *int_yield_sig = new char[100];
    char *int_yield_bckg = new char[100];
TF1* yield_fit_sig = new TF1(int_yield_sig,yield_signal_formula,1.09,1.15,0);
TF1* yield_fit_bckg = new TF1(int_yield_bckg,yield_bckg_formula,1.09,1.15,0);

TCanvas *c3 = new TCanvas("c3","c3",-1,10,800,600);
yield_fit_sig->SetLineColor(kRed);
yield_fit_bckg->SetLineColor(kBlue);
yield_fit_sig->Draw();
yield_fit_bckg->Draw("same");
c3->SaveAs(Form("/home/valeriy/bmn_lamb/hardcut/yields_func_lamb_eta_bin_%d.png",i));
TCanvas *c4 = new TCanvas("c4","c4",-1,10,800,600);
h_inv_mass_sig->SetLineColor(kRed-2);
h_inv_mass_bckg->SetLineColor(kBlue-2);
h_inv_mass_sig->Draw();
h_inv_mass_bckg->Draw("same");
fitting_bckg_full->Draw("same");
fitting_sig_full->Draw("same");
c4->SaveAs(Form("/home/valeriy/bmn_lamb/hardcut/sig_bckg_lamb_eta_bin_%d.png",i));



	c2->SaveAs(Form("/home/valeriy/bmn_lamb/hardcut/yields_lamb_eta_bin_%d.png",i));

}
  auto file = TFile::Open( "/home/valeriy/bmn_lamb/JAM_lambda_corr_bckg_peak_rejected_060225_1.root");
  std::vector<std::string> ep_vectors{ "F1_RESCALED", "F2_RESCALED", "F3_RESCALED" };
  std::vector<std::string> res_vectors{ "F1_RESCALED", "F2_RESCALED", "F3_RESCALED"};
  std::vector<std::string> sts_vectors{ "lambda_good_RESCALED"};
  std::vector<std::string> psi_vectors{ "psi_rp_PLAIN"};
  std::vector<std::string> sts_vectors2{ "lambda_good_RESCALED","tru_lambda_RESCALED"};
  std::array<std::string, 2> components{"x1x1centrality", "y1y1centrality"};

  auto file_out = TFile::Open( "~/bmn_lamb/JAM_lambda_correlation_test_290125.root", "recreate" );
  file_out->cd();
  file_out->mkdir("resolutions");
  file_out->mkdir("lambda");

  for( auto qa : ep_vectors ){
    Correlation<2> lambda_qa( file, "/", std::array{"lambda_good_RESCALED"s, qa}, components);
    auto res_v = Functions::VectorResolutions3S( file, "/", qa, res_vectors, components );
    for( auto R1 : res_v ) {

      auto v1_lambda = lambda_qa * 2 / R1;

      file_out->cd("resolutions");
      R1.Save("R1."+R1.Title());


      file_out->cd("lambda");
      v1_lambda.Save("v1."+R1.Title());
      

    }
  }
   for( auto psiqa : psi_vectors ){
  Correlation<2> lambda_qa( file, "/", std::array{"lambda_good_RESCALED"s,psiqa }, components);
  lambda_qa=lambda_qa*2;
   lambda_qa.Save("v1.psi_rp_PLAIN");
   Correlation<2> tru_lambda_qa( file, "/", std::array{"tru_lambda_RESCALED"s,psiqa }, components);
   tru_lambda_qa.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "simPt", 1, 0.4, 1.2 } } );
  tru_lambda_qa.Project({{"simProtonY"}} );
  tru_lambda_qa=tru_lambda_qa*2;
   tru_lambda_qa.Save("tru_v1.psi_rp_PLAIN");
   }
  for( auto psiqa : psi_vectors ){
  Correlation<2> lambda_qa_sig( file, "/", std::array{"signal_RESCALED"s,psiqa }, components);
  lambda_qa_sig=lambda_qa_sig*2;
  lambda_qa_sig.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 },/* { "candidate_rapidity", 1, 0.3, 1 }*//*, { "candidate_mass", 1, 1.105, 1.126 }*/ } );
//  lambda_qa.Project({{"candidate_mass"}} );
  lambda_qa_sig.Project({{"candidate_rapidity"}} );
  lambda_qa_sig.Save("v1_sig.psi_rp_PLAIN");

  Correlation<2> lambda_qa_bckg( file, "/", std::array{"background_RESCALED"s,psiqa }, components);
  lambda_qa_bckg=lambda_qa_bckg*2;
  lambda_qa_bckg.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 }/*, { "candidate_rapidity", 1, 0.3, 1 }*//*,{ "candidate_mass", 1, 1.105, 1.126 }*/ } );
  lambda_qa_bckg.Project({{"candidate_rapidity"}} );
 // lambda_qa.Project({{"candidate_mass"}} );
   lambda_qa_bckg.Save("v1_bckg.psi_rp_PLAIN");
}
  std::vector<std::string> vect{"v1.F1_RESCALED(F2_RESCALED,F3_RESCALED)","v1.F2_RESCALED(F1_RESCALED,F3_RESCALED)","v1.F3_RESCALED(F1_RESCALED,F2_RESCALED)","v1.psi_rp_PLAIN"};
  std::array<std::string, 1> comp{"y1y1centrality"};

for(int k=0;k<3;k++){

  for( int i=1;i<8;i++){

        for(int j=1;j<6;j++){
		Correlation<1> v1_lamb( file_out,
               "lambda",
               std::vector{vect.at(k)},
               comp);
                double gappt[8]={0,0.2,0.4,0.6,0.8,1.0,1.2,2.0};
double gapy[6]={-0.2,0,0.3,0.5,0.8,1.0};
                v1_lamb.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, gappt[i-1], gappt[i] },{"candidate_rapidity", 1,gapy[j-1],gapy[j] }} );
  v1_lamb.Project({{"candidate_mass"}} );
                v1_lamb.Save(((vect.at(k)+"y1y1centrality").append("1234567",i-1,1)).append("23456",j-1,1));
        }}
}

char* int_af_fit= new char[100];
std::vector<TF1*> fit_vector_eta_af={new TF1(int_af_fit,"[0]+x*[1]+x*x*x*[2]",0,1)};
for(int i=3;i<4;i++){


  Correlation<1> v1_lamb( file_out,
               "lambda",
               std::vector{vect.at(i)},
	       comp);
  v1_lamb.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_rapidity", 1, 0.3, 1 } } );
  v1_lamb.Project({{"candidate_mass"},{"candidate_pT"}} );
  Fitter v1_fit(v1_lamb,"candidate_mass","candidate_pT");
      v1_fit.Fit(fit_vector_pt,3);
      std::vector<Qn::DataContainerStatCalculate> param= v1_fit.GetFitResults();
      param.at(0).Write(("pT_"+vect.at(i)+".y1y1centrality").c_str());
v1_fit.DumpGraphs("~/bmn_lamb/calc2_test_fit_graphs_pt.root");
      Correlation<1> v1_lamb_eta( file_out,
               "lambda",
               std::vector{vect.at(i)},
               comp);
      v1_lamb_eta.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 }/*,{ "candidate_rapidity", 1, 0.3, 1 }*/ } );
  v1_lamb_eta.Project({{"candidate_mass"},{"candidate_rapidity"}} );
  v1_lamb_eta.Save(("mass_"+vect.at(i)+".y1y1centrality").c_str());
  std::cout<<"debug"<<std::endl;
  Fitter v1_fit_eta(v1_lamb_eta,"candidate_mass","candidate_rapidity");
      v1_fit_eta.Fit(fit_vector_eta,3);
      std::vector<Qn::DataContainerStatCalculate> param_eta= v1_fit_eta.GetFitResults();
      param_eta.at(0).Write(("eta_"+vect.at(i)+".y1y1centrality").c_str());

      v1_fit_eta.DumpGraphs("~/bmn_lamb/calc2_test_fit_graphs_eta_3.root");
      }
file_out->Close();

/*    Correlation v1_lamb_af("v1_lamb_af",std::array{vect.at(i)},std::array{param_eta.at(0)},vect);
    Fitter v1_fit_eta_2(v1_lamb_af,"candidate_rapidity","candidate_pT");
    std::cout<<"debug"<<std::endl;
      v1_fit_eta_2.Fit(fit_vector_eta_af,3);*/
//}
//      v1_fit.DumpGraphs("~/bmn_lamb/calc2_test_fit_graphs_pt.root");

//      v1_fit_eta.DumpGraphs("~/bmn_lamb/calc2_test_fit_graphs_eta_3.root");
auto file2 = TFile::Open( "/home/valeriy/bmn_lamb/JAM_lambda_corr_bckg_peak_rejected_060225_1.root");
auto file_out2 = TFile::Open( "~/bmn_lamb/JAM_lambda_correlation_test_290125.root", "update" );
file_out2->cd();
file_out2->cd("lambda");
for( auto psiqa : psi_vectors ){
  Correlation<2> lambda_qa_sig( file2, "/", std::array{"signal_RESCALED"s,psiqa }, components);
  lambda_qa_sig=lambda_qa_sig*2;
  lambda_qa_sig.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 }, { "candidate_rapidity", 1, 0.3, 0.5 }/*, { "candidate_mass", 1, 1.105, 1.126 }*/ } );
  lambda_qa_sig.Project({{"candidate_mass"}} );
  lambda_qa_sig.Save("v1_sig_3.psi_rp_PLAIN");

  Correlation<2> lambda_qa_bckg( file2, "/", std::array{"background_RESCALED"s,psiqa }, components);
  lambda_qa_bckg=lambda_qa_bckg*2;
  lambda_qa_bckg.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 }, { "candidate_rapidity", 1, 0.3, 0.5 }/*,{ "candidate_mass", 1, 1.105, 1.126 }*/ } );
  lambda_qa_bckg.Project({{"candidate_mass"}} );
   lambda_qa_bckg.Save("v1_bckg_3.psi_rp_PLAIN");

   Correlation<2> lambda_qa( file2, "/", std::array{"lambda_good_RESCALED"s,psiqa }, components);
  lambda_qa=lambda_qa*2;
  lambda_qa.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 },{ "candidate_rapidity", 1, 0.3, 0.5 } } );
  lambda_qa.Project({{"candidate_mass"}} );
   lambda_qa.Save("v1_3.psi_rp_PLAIN");

   char *int_sig = new char[100];
    char *int_bckg = new char[100];
fit_vector_eta.at(2)->Write("v1_total_fit_3");
TF1* sig = new TF1(int_sig,"[0]",1.09,1.15);
sig->FixParameter(0,fit_vector_eta.at(2)->GetParameter(0));
TF1* bckg = new TF1(int_bckg,"[1]*x+[0]",1.09,1.15);
bckg->FixParameter(0,fit_vector_eta.at(2)->GetParameter(1));
bckg->FixParameter(1,fit_vector_eta.at(2)->GetParameter(2));
sig->SetLineColor(kRed);
bckg->SetLineColor(kBlue);
sig->Write("v1_sig_fit_3");
bckg->Write("v1_bckg_fit_3");


}

for( auto psiqa : psi_vectors ){
  Correlation<2> lambda_qa_sig( file2, "/", std::array{"signal_RESCALED"s,psiqa }, components);
  lambda_qa_sig=lambda_qa_sig*2;
  lambda_qa_sig.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 }, { "candidate_rapidity", 1, 0.5, 0.8 }/*, { "candidate_mass", 1, 1.105, 1.126 }*/ } );
  lambda_qa_sig.Project({{"candidate_mass"}} );
  lambda_qa_sig.Save("v1_sig_4.psi_rp_PLAIN");

  Correlation<2> lambda_qa_bckg( file2, "/", std::array{"background_RESCALED"s,psiqa }, components);
  lambda_qa_bckg=lambda_qa_bckg*2;
  lambda_qa_bckg.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 }, { "candidate_rapidity", 1, 0.5, 0.8 }/*,{ "candidate_mass", 1, 1.105, 1.126 }*/ } );
  lambda_qa_bckg.Project({{"candidate_mass"}} );
   lambda_qa_bckg.Save("v1_bckg_4.psi_rp_PLAIN");

   Correlation<2> lambda_qa( file2, "/", std::array{"lambda_good_RESCALED"s,psiqa }, components);
  lambda_qa=lambda_qa*2;
  lambda_qa.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 },{ "candidate_rapidity", 1, 0.5, 0.8 } } );
  lambda_qa.Project({{"candidate_mass"}} );
   lambda_qa.Save("v1_4.psi_rp_PLAIN");

char *int_sig = new char[100];
    char *int_bckg = new char[100];
fit_vector_eta.at(3)->Write("v1_total_fit_4");
TF1* sig = new TF1(int_sig,"[0]",1.09,1.15);
sig->FixParameter(0,fit_vector_eta.at(3)->GetParameter(0));
TF1* bckg = new TF1(int_bckg,"[1]*x+[0]",1.09,1.15);
bckg->FixParameter(0,fit_vector_eta.at(3)->GetParameter(1));
bckg->FixParameter(1,fit_vector_eta.at(3)->GetParameter(2));
sig->SetLineColor(kRed);
bckg->SetLineColor(kBlue);
sig->Write("v1_sig_fit_4");
bckg->Write("v1_bckg_fit_4");
}

for( auto psiqa : psi_vectors ){
  Correlation<2> lambda_qa_sig( file2, "/", std::array{"signal_RESCALED"s,psiqa }, components);
  lambda_qa_sig=lambda_qa_sig*2;
  lambda_qa_sig.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 }, { "candidate_rapidity", 1, 0.8, 1 }/*, { "candidate_mass", 1, 1.105, 1.126 }*/ } );
  lambda_qa_sig.Project({{"candidate_mass"}} );
  lambda_qa_sig.Save("v1_sig_5.psi_rp_PLAIN");

  Correlation<2> lambda_qa_bckg( file2, "/", std::array{"background_RESCALED"s,psiqa }, components);
  lambda_qa_bckg=lambda_qa_bckg*2;
  lambda_qa_bckg.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 }, { "candidate_rapidity", 1, 0.8, 1 }/*,{ "candidate_mass", 1, 1.105, 1.126 }*/ } );
  lambda_qa_bckg.Project({{"candidate_mass"}} );
   lambda_qa_bckg.Save("v1_bckg_5.psi_rp_PLAIN");

   Correlation<2> lambda_qa( file2, "/", std::array{"lambda_good_RESCALED"s,psiqa }, components);
  lambda_qa=lambda_qa*2;
  lambda_qa.Rebin( std::vector <Qn::AxisD> { {"centrality", 1, 10, 40}, { "candidate_pT", 1, 0.4, 1.2 },{ "candidate_rapidity", 1, 0.8, 1 } } );
  lambda_qa.Project({{"candidate_mass"}} );
   lambda_qa.Save("v1_5.psi_rp_PLAIN");

   char *int_sig = new char[100];
    char *int_bckg = new char[100];
fit_vector_eta.at(4)->Write("v1_total_fit_5");
TF1* sig = new TF1(int_sig,"[0]",1.09,1.15);
sig->FixParameter(0,fit_vector_eta.at(4)->GetParameter(0));
TF1* bckg = new TF1(int_bckg,"[1]*x+[0]",1.09,1.15);
bckg->FixParameter(0,fit_vector_eta.at(4)->GetParameter(1));
bckg->FixParameter(1,fit_vector_eta.at(4)->GetParameter(2));
sig->SetLineColor(kRed);
bckg->SetLineColor(kBlue);
sig->Write("v1_sig_fit_5");
bckg->Write("v1_bckg_fit_5");

}



  file_out->Close();
  file->Close();
}
