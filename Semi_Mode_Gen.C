
// This macro calculates the corrections needed in the Semi-analytic
// scintillation light model. These corrections are parametrised as
// Gaisser-Hillas functions. These functions depend on: 1) distance,
// 2) offset-angle, and 3) distance to the center (border) of the
// LArTPC active volume. It allows to correct for the LAr light signals
// for optical detector with flat sensitive windows (circular and
// rectangular shepes) or semisphere (by a subsequent disks approach).

// Author: Diego Garcia-Gamez (dgamez@fnal.gov)
// Created: 20.02.2020


Double_t pol1(double *x,double *par)
{ 
  return (*x*par[1] + par[0]);
}
Double_t pol0(double *x,double *par)
{ 
  return par[0];
}

TGraphErrors *Profile_to_Graph(TProfile *p) {
  
  TAxis *xaxis = p->GetXaxis(); 
  Int_t nbins = xaxis->GetNbins(); 
  vector<double> vx, vy, vey;
  for (Int_t bin=0; bin <= nbins; bin++) {
    double y_value = p->GetBinContent(bin);
    if(y_value == 0) continue;
    vx.push_back(xaxis->GetBinCenter(bin));
    vy.push_back(y_value);
    vey.push_back(p->GetBinError(bin));
  }
  TGraphErrors *gr = new TGraphErrors(vx.size()-1, &vx[0], &vy[0], 0, &vey[0]);
  
  return gr;
  
}
///////////////////////////////////////////////////
Double_t GaisserHillas(double *x,double *par)
{
  //This is the Gaisser-Hillas function
  Double_t X_mu_0=par[3];
  Double_t Normalization=par[0];
  Double_t Diff=par[1]-X_mu_0;
  Double_t Term=pow((*x-X_mu_0)/Diff,Diff/par[2]);
  Double_t Exponential=TMath::Exp((par[1]-*x)/par[2]);
  
  return ( Normalization*Term*Exponential);
}

//Distance to the center in the Y-Z Plane
double GetDistanceCenter(const double center[2], double z, double y){
  z -= center[1];
  y -= center[0];
  
  return  sqrt(y*y + z*z);
}  

// solid angle of rectanglular aperture
// structure definition for solid angle of rectangle function
struct acc{
  // ax,ay,az = centre of rectangle; w = width; h = height
  double ax, ay, az, w, h;
};

double Rectangle_SolidAngle(double a, double b, double d){
  
  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));
  
}

double Rectangle_SolidAngle(acc& out, TVector3 v){
  
  //v is the position of the track segment with respect to
  //the center position of the arapuca window
  
  // arapuca plane fixed in x direction
  
  if( v.Y()==0.0 && v.Z()==0.0){
    return Rectangle_SolidAngle(out.w,out.h,v.X());
  }
  
  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.X());
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(B+b),d)-Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*(A+a),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }
  
  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.X());
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(b-B),d)+Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }
  
  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.X());
    double to_return = (Rectangle_SolidAngle(2*(A+a),2*(b-B),d)-Rectangle_SolidAngle(2*A,2*(b-B),d)+Rectangle_SolidAngle(2*(A+a),2*B,d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }
  
  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = std::abs(v.X());
    double to_return = (Rectangle_SolidAngle(2*(a-A),2*(B+b),d)-Rectangle_SolidAngle(2*(a-A),2*B,d)+Rectangle_SolidAngle(2*A,2*(B+b),d)-Rectangle_SolidAngle(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}

bool _mathmore_loaded_ = false;

double Omega(double* x, double *p) {
  const double d = x[0];
  const double h = x[1];
  const double b = p[0];
  if(b <= 0. || d < 0. || h <= 0.) return 0.; 
  const double aa = TMath::Sqrt(h*h/(h*h+(b+d)*(b+d)));
  if(d == 0) {
    return 2.*TMath::Pi()*(1.-aa);
  }
  const double bb = TMath::Sqrt(4*b*d/(h*h+(b+d)*(b+d)));
  const double cc = 4*b*d/((b+d)*(b+d));
  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  if(TMath::Abs(ROOT::Math::comp_ellint_1(bb) - bb) < 1e-10 && TMath::Abs(ROOT::Math::comp_ellint_3(cc,bb) - cc) <1e-10) {
    throw(std::runtime_error("please do gSystem->Load(\"libMathMore.so\") before running Omega for the first time!"));
  }
  if(d < b) {
    return 2.*TMath::Pi() - 2.*aa*(ROOT::Math::comp_ellint_1(bb) + TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb));
  }
  if(d == b) {
    return TMath::Pi() - 2.*aa*ROOT::Math::comp_ellint_1(bb);
  }
  if(d > b) {
    return 2.*aa*(TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb) - ROOT::Math::comp_ellint_1(bb));
  }
  return 0.;
}

double Omega(double d, double h, double b) {
  double x[2] = { d, h };
  double p[1] = { b };
  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  return Omega(x,p);
}

Double_t Omega_Dome(const int n, double d, double h0, double b0)
{
  // This function approaches a spheric section as a sequence of disks
  // where n is the number of disks used in the approach
  // d -> offset in z-y plane
  // h -> drift distance (in x)
  // b -> radius
  //distance between the different layers
  double step = b0/n;
  //The subsequent disk radios will be
  double b[n], h[n];
  for(int i=0; i<n; i++) {
    double gap = i*step;
    b[i] = sqrt(b0*b0 - gap*gap);
    h[i] = h0 - gap;   
  }

  double geo_factor = 0;
  if(n > 1)  {
    for(int i=0; i<n-1; i++) 
      geo_factor += (Omega(d, h[i], b[i]) - Omega(d, h[i], b[i+1]));
    geo_factor += Omega(d, h[n-1], b[n-1]);
  }
  else
    geo_factor = Omega(d, h0, b0);
  
  return geo_factor;
}

Double_t Omega_Dome_Model(double distance,double theta)
{
  // this function calculates the solid angle of a semi-sphere of radius b,
  // as a correction to the analytic formula of the on-axix solid angle,
  // as we move off-axis an angle theta. We have used 9-angular bins
  // with delta_theta width.
  
  // par0 = Radius correction close
  // par1 = Radius correction far
  // par2 = breaking distance betwween "close" and "far"

  double par0[9] = {0., 0., 0., 0., 0., 0.597542, 1.00872, 1.46993, 2.04221}; 
  double par1[9] = {0, 0, 0.19569, 0.300449, 0.555598, 0.854939, 1.39166, 2.19141, 2.57732};
  const double delta_theta = 10.;
  int j = int(theta/delta_theta);
  // 8" PMT radius
  const double b = 8*2.54/2.;
  // distance form which the model parameters break (empirical value)
  const double d_break = 5*b;//par2
  
  if(distance >= d_break) {
    double R_apparent_far = b - par1[j];
    return  (2*3.1416 * (1 - sqrt(1 - pow(R_apparent_far/distance,2))));
    
  }
  else {
    double R_apparent_close = b - par0[j];
    return (2*3.1416 * (1 - sqrt(1 - pow(R_apparent_close/distance,2))));
  }
  
}

void calcula(string positions, string path_files, string lista_files,
	     vector<double> &v1, vector<double> &v2, vector<double> &v3,
	     vector<double> &v4, vector<double> &v5,
	     bool IsRectangular, bool IsSphere, const int nLayers) {
  cout<<"calcula function ..."<<endl;

  double min_number_entries = 0;
  // width and height in cm of arapuca active window	
  double arapuca_w = 9.3;
  double arapuca_h = 46.8;
  // 8" PMT radius
  double b = 8*2.54/2.;
  // Y-Z coordinates of the active volume center
  //SBND geometry v01_04 values:
  const double centerYZ[2] = {0., 250.};
  //DUNE geometry dune1x2x6 values:
  //const double centerYZ[2] = {0., 700.};	
  // LAr absorption length in cm
  // This needs to match with the value used in the full-geant4 simulation
  const double L_abs = 2000.;
  
  gRandom->SetSeed(0);
  //getting the pmt positions (y and z)
  ifstream Traks_file1(positions.c_str());
  if(!Traks_file1) cerr << "WARNING:  Failed to open file with optical detector positions"<< endl;
  Traks_file1.seekg(0);
  vector<double> devx;
  vector<double> devy;
  vector<double> devz;
  
  double id, x, y, z;
  while(!(Traks_file1.eof())) { 
    Traks_file1 >> id >> x >> y >> z;
    devx.push_back(x);
    devy.push_back(y);
    devz.push_back(z);
  }

  //getting the file names
  ifstream Traks_file2(lista_files.c_str());
  if(!Traks_file2) cerr << "WARNING:  Failed to open file with Input file names"<< endl;
  Traks_file2.seekg(0);
  vector<string> names;
  string nombre;
  while(!(Traks_file2.eof())) { 
    Traks_file2 >> nombre;
    names.push_back(nombre);
  }


  const int numberDevices = devy.size() - 1;
  const int n_files = names.size() -1;
  cout<<"----> numberDevices: "<<numberDevices<<endl;
  cout<<"----> number of files: "<<names.size() -1<<endl;
  // LY used to estimate the "equivalent" dedx 
  double scint_light_yiels = 24000; //photons/MeV
  
  vector<double> v_hits, v_distance, v_rec_hits, v_offset_angle, v_d_center;
  //loop over files
  for(int n=0; n<n_files; n++) {
    
    int nDL[numberDevices];
    for(int i=0; i<numberDevices; i++)
      nDL[i] = 0;
      
    char inFile[120];
    sprintf(inFile,"%s%s",path_files.c_str(),names.at(n).c_str());    
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"Name of the input file: "<<inFile<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    TFile* f = new TFile(inFile);
    TTree *tree = (TTree *)f->Get("AllPhotons");   
    Int_t EventID, OpChannel;
    Float_t Wavelength, Time;
    tree->SetBranchAddress("EventID", &EventID);
    tree->SetBranchAddress("OpChannel", &OpChannel);
    tree->SetBranchAddress("Wavelength", &Wavelength);
    tree->SetBranchAddress("Time", &Time);

    TTree  *tree2 = (TTree *)f->Get("PhotonsGenerated");
    double X, Y, Z;
    int gen_photons;
    tree2->SetBranchAddress("X", &X);
    tree2->SetBranchAddress("Y", &Y);
    tree2->SetBranchAddress("Z", &Z);
    tree2->SetBranchAddress("entries", &gen_photons);
    tree2->GetEntry(0);
    double posSource[3]={X, Y, Z};
    double num_phot_generated = tree2->GetEntries()*gen_photons;
    double dedx = num_phot_generated/scint_light_yiels;
 
    cout<<"Energy Deposited [MeV] = "<<dedx<<" ---> N= "<<num_phot_generated<<endl;
    cout<<"scintillation point --> (x, y, z) = "<<posSource[0]<<"  "<<posSource[1]<<"  "<<posSource[2]<<endl;

    double distance_to_center = GetDistanceCenter(centerYZ, posSource[2], posSource[1]);
    
    for(int i=0; i!=tree->GetEntries(); ++i)
      {
	tree->GetEntry(i);
	if (Wavelength < 300) 
	  nDL[OpChannel]++;
      }
    
    //loop over the channels
    for(int i=0; i<numberDevices; i++) {
      // Next lines ONLY apply for SBND case to select 
      // only PMTs in TPC with x>0 for this test:	
      //if(devx.at(i) != 211.465) continue;
      
      int entries = nDL[i];
      if(entries < min_number_entries) continue;
      
      double distance_to_pmt = sqrt(pow(posSource[0] - devx.at(i),2) +
				    pow(posSource[1] - devy.at(i),2) +
				    pow(posSource[2] - devz.at(i),2));
      double coseno = sqrt(pow(posSource[0] - devx.at(i),2))/distance_to_pmt;
      double theta = acos(coseno)*180/3.1416;
      
      double geo_factor = -1;
      if(IsRectangular) {
	//------Rectangle---------
	acc detPoint;
	// centre coordinates of optical detector
	detPoint.ax = devx.at(i); detPoint.ay = devy.at(i); detPoint.az = devz.at(i);
	// width and height in cm of arapuca active window	
	detPoint.w = arapuca_w; detPoint.h = arapuca_h; 
	// get scintillation point coordinates relative to arapuca window centre
	TVector3 ScintPoint(posSource[0], posSource[1], posSource[2]);
	TVector3 OpDetPoint(devx.at(i), devy.at(i), devz.at(i));
	TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;	
	// calculate solid angle
	geo_factor  = Rectangle_SolidAngle(detPoint, ScintPoint_rel);
      }	
      else {
	if(IsSphere) {
	  geo_factor = Omega_Dome_Model(distance_to_pmt, theta);
	}
	else {
	//------Disk approach---------	
	double d,h;
	//offset in z-y plane
	d = sqrt((posSource[2] - devz.at(i))*(posSource[2] - devz.at(i)) + (posSource[1] - devy.at(i))*(posSource[1] - devy.at(i)));
	//drift distance (in x)
	h =  sqrt((devx.at(i) - posSource[0])*(devx.at(i) - posSource[0]));	
	// calculate solid angle
	if(theta<45) // semi-sphere approached as subsequent disks (for on-axis cases)
	  geo_factor = Omega_Dome(nLayers, d, h, b);
	else // semi-sphere approached as single disk (very off-axis cases)
	  geo_factor = Omega_Dome(1, d, h, b);
	//----------------------
      }
    }
      
      //pure geometric estimation of the number of arriving VUV photons
      double rec_N =  exp(-1.*distance_to_pmt/L_abs)*gRandom->Poisson(num_phot_generated*geo_factor/(4*3.1416));
      
      v_hits.push_back(entries);
      v_distance.push_back(distance_to_pmt);
      v_rec_hits.push_back(rec_N);
      v_offset_angle.push_back(theta);
      v_d_center.push_back(distance_to_center);
      
    }// channels
    delete f;
  }//files

  v1 = v_distance;
  v2 = v_hits;
  v3 = v_rec_hits;
  v4 = v_offset_angle;
  v5 = v_d_center;
}

void Semi_Mode_Gen() {

  //path to the directory where the Output files are 
  string path_files = "/PATH_TO_INPUT_FILES/";
  // path + file containing the LDS id and positions (x, y, z)
  string positions = path_files + "opch_info.txt";
  // path + file containing the list of Output file names
  string lista_files = path_files + "list_file_names.txt";

  // Two layer/disks are enough for an accurate correction
  // Keep in mind that you will use this approach if the estimation
  // of the corrections and also in the reconstruction of the signals
  // by the estimated correctio
  const int nLayers = 1;

  //If PMTs false, if (X-)ARAPUCAS true
  bool IsRectangular = false;
  //Option true for the PMTs semi-spheric
  bool IsSphere = true;	
	
  //offset- angle theta binning
  const int N = 9;
  double delta_angulo = 90./N;
  vector<double> angulo;

  double theta[N];
  for(int i=0; i<N; i++)
    theta[i] =  delta_angulo/2 + i*delta_angulo;
  
  //Distance range and step for the profile to fit with GH
  //SBND suggested values:
  double d_min = 0;
  double d_max = 500.;
  double step_d =25;
  //DUNE suggested values:
  /*double d_min = 0;
  double d_max = 1000.;
  double step_d =40;*/
	
  //Center distance bins
  //SBND suggested values:
  const int M = 8;
  double range_d = 320;
  //DUNE suggested values:
  /*const int M = 4;
  double range_d =1000;*/
  double delta_d = range_d/M;
  TH1D* h=new TH1D("","",range_d, 0, range_d);

  // Initializing parameters for GH fit
  string options = "W0Q";
  //SBND suggested values:
  double pars_ini[4] = {1., 100., 50, -250};//90cm <RS>
  //DUNE suggested values:
  //double pars_ini[4]= {1., 100., 50, -1000};//90cm <RS>
  TF1 *GH[N][M];
  TH1D* hd_centers[M];
  for(int k=0; k < M; k++) {
    hd_centers[k] = new TH1D("","",M,0,range_d);
    for(int j=0; j < N; j++) {
      GH[j][k] =  new TF1("GH",GaisserHillas,0.,d_max,4);
      GH[j][k]->SetParameters(pars_ini);
      GH[j][k]->SetParLimits(1,10,500);
      GH[j][k]->SetParLimits(2,10,2000);
      GH[j][k]->SetParLimits(3,-1000,0);
      //The user might need to modify the actual values of these
      //parameters depending on his/her particular case
      GH[j][k]->FixParameter(2,pars_ini[2]);
      GH[j][k]->FixParameter(3,pars_ini[3]);
      
    }
  }
  // Defining Profiles we want to fit
  TProfile* pdiff_d[N][M];
  for(int j=0; j < N; j++) {
    double theta_min = j*delta_angulo;
    double theta_max = theta_min + delta_angulo;
    angulo.push_back(theta_min + delta_angulo/2);
    for(int k=0; k < M; k++) {
      pdiff_d[j][k]  = new TProfile("","", int((d_max-d_min)/step_d), d_min, d_max, "s");
    }
  }
  // Calculating true and rec variables for the analysis
  vector<double> v_distance, v_hits, v_rec_hits, v_offset_angle, v_d_center;
  calcula(positions, path_files, lista_files, v_distance, v_hits, v_rec_hits, v_offset_angle, v_d_center, IsRectangular, IsSphere, nLayers);
  for(int i=0; i<v_distance.size(); i++) {
    double costheta = cos(3.1416*v_offset_angle.at(i)/180.);
    //which angulat bin
    int j = int(v_offset_angle.at(i)/delta_angulo);
    //which "center/crown" bin
    int k = int(v_d_center.at(i)/delta_d);
    pdiff_d[j][k]->Fill(v_distance.at(i), v_hits.at(i)/v_rec_hits.at(i)*costheta);
    hd_centers[k]->Fill(v_d_center.at(i));
    h->Fill(v_d_center.at(i));
  }

  // This plot is to see the distribution of "distances to the center" of our Y-Z plane
  // This will help us to optimise/verify our "crowns" definition is right (unbiased)
  TCanvas *canvas0 = new TCanvas("canvas0", "graph draw options",200,200,500,400);
  h->SetTitle("Help to choose the \"border-study\" distance bins");
  h->GetXaxis()->SetTitle("distance to centre in Y-Z plane [cm]");
  h->Draw("hist");
  canvas0->Update();
  canvas0->Modified();
  
  // Estimate range to fit in distance
  // Only wherever we have data points
  double max_x[N][M], min_x[N][M], n_entries[N][M];
  for(int j=0; j < N; j++) {
    for(int k=0; k < M; k++) {
      n_entries[j][k] = 0;
      TAxis *xaxis  =  pdiff_d[j][k]->GetXaxis(); 
      Int_t nbins  = xaxis->GetNbins(); 
      double min = d_max;
      double max = 0;
      for (Int_t bin=0;bin<=nbins;bin++) {
	n_entries[j][k] += pdiff_d[j][k]->GetBinContent(bin);
	
	if(pdiff_d[j][k]->GetBinContent(bin) == 0) continue;
	
	if(min > xaxis->GetBinCenter(bin))
	  min = xaxis->GetBinCenter(bin);
	if(max < xaxis->GetBinCenter(bin))
	  max = xaxis->GetBinCenter(bin);	
      }
      max_x[j][k] = max;
      min_x[j][k] = min;
    }
  }

  // This will be usefull for plotting
  vector<int> N_canvas;
  vector<double> d_center;
  string title[M];
  for(int k=0; k < M; k++) {
    for(int j=0; j < N; j++) 
      if(n_entries[j][k]>0) {
	N_canvas.push_back(k);
	break;
      }
    // save the center of our distance to center bins/crowns
    if(hd_centers[k]->GetEntries()>0) {
      d_center.push_back(hd_centers[k]->GetMean());
      title[k] =  "Distance to center " +  to_string(int(hd_centers[k]->GetMean())) + "cm";
    }
  }

  // Making plots and fits ...
  const int dim = N_canvas.size();
  TGraphErrors *gr[N][M];
  for(int k=0; k < M; k++) {
    for(int j=0; j < N; j++) {    
      if(n_entries[j][k] > 0) {
	gr[j][k] = Profile_to_Graph(pdiff_d[j][k]);
	gr[j][k]->SetLineColor(1 + j);
	gr[j][k]->SetMarkerColor(1 + j);
	gr[j][k]->SetMarkerStyle(20 + j);
	gr[j][k]->SetMarkerSize(0.6);
      }
    }
  }

  TCanvas *canvas1 = new TCanvas("canvas1", "graph draw options",200,200,1000,1000);
  if(N_canvas.size() > 1) {
    double nn = double(N_canvas.size())/2.;
    canvas1->Divide(2,ceil(nn));
  }
  double x_0[2] = {0, d_max};
  double y_0[2] = {0, 2};
  TGraph* gg0[dim];
  vector<double> p1[dim], p2[dim], p3[dim], p4[dim], ep1[dim], ep2[dim];
  TLegend *leg1=new TLegend(0.55, 0.4, 0.89, 0.89,NULL,"brNDC");
  leg1->SetHeader("");
  leg1->SetBorderSize(0);
  char label[N][20];
  for(int l=0; l < dim; l++) {
    int k = N_canvas.at(l);
    canvas1->cd(1 + l);
    gg0[l] = new TGraph(2,x_0,y_0);
    gg0[l]->SetTitle(title[k].c_str());
    gg0[l]->GetYaxis()->SetTitleSize(0.05);
    gg0[l]->GetYaxis()->SetTitleOffset(1.);
    gg0[l]->GetXaxis()->SetTitleSize(0.05);
    gg0[l]->GetXaxis()->SetTitleOffset(1.);
    gg0[l]->GetXaxis()->SetRangeUser(0,d_max);
    gg0[l]->GetYaxis()->SetRangeUser(0, 2);
    gg0[l]->GetXaxis()->SetTitle("distance [cm]");
    gg0[l]->GetYaxis()->SetTitle("N_{hit} / N_{#Omega}");
    gg0[l]->Draw("ap");
    
    for(int j=0; j < N; j++) {
      double pars_GH[4] = {-999, -999, -999, -999};
      double epars_GH[4]= {-999, -999, -999, -999};
      if(n_entries[j][k]>0) {
	gr[j][k]->Draw("p");
	// Fitting the simulation data
	gr[j][k]->Fit(GH[j][k],options.c_str(),"",min_x[j][k],max_x[j][k]);
	//Loading parameters	
	GH[j][k]->GetParameters(pars_GH);
	GH[j][k]->SetParameters(pars_GH);
	GH[j][k]->SetLineColor(1 + j);
	GH[j][k]->SetLineStyle(kDashed);
	GH[j][k]->SetRange(0*min_x[j][k], max_x[j][k]); 
	GH[j][k]->Draw("same");
	const double* epars_GH_ = GH[j][k]->GetParErrors();
	*epars_GH = *epars_GH_;
      }
      p1[l].push_back(pars_GH[0]);
      p2[l].push_back(pars_GH[1]);
      p3[l].push_back(pars_GH[2]); 
      p4[l].push_back(pars_GH[3]);
      ep1[l].push_back(epars_GH[0]);
      ep2[l].push_back(epars_GH[1]);
      if(l==0) {
	int a_min = j*delta_angulo;
	int a_max = (j+1)*delta_angulo;
	sprintf(label[j],"#theta #in [%i, %i] deg",a_min, a_max);
	leg1->AddEntry(gr[j][k],label[j],"p");
      }      
    }
    if(l==0)  leg1->Draw();
  }
  canvas1->Update();
  canvas1->Modified();
  canvas1->WaitPrimitive();
  // printing the GH parameters we get at our different conditions
  cout<<"  GH p1: "<<endl;
  for(int k=0; k < dim; k++){ 
    for(int j=0; j<N; j++)
      cout<<p1[k].at(j)<<", ";
    cout<<""<<endl;
  }
  cout<<"  GH p2: "<<endl;
  for(int k=0; k < dim; k++){
    for(int j=0; j<N; j++)
      cout<<p2[k].at(j)<<", ";
    cout<<""<<endl;
  }
  cout<<"  GH p3: "<<endl;
  for(int k=0; k < dim; k++){
    for(int j=0; j<N; j++)
      cout<<p3[k].at(j)<<", ";
    cout<<""<<endl;
  }
  cout<<"  GH p4: "<<endl;
  for(int k=0; k < dim; k++){ 
    for(int j=0; j<N; j++)
      cout<<p4[k].at(j)<<", ";
    cout<<""<<endl;
  }

  // For border effects
  // Studying the dependency of the corrections with how
  // far from the center (or close to the borders) we are
  //--------------------------------------------------------------------------
  TCanvas *c1 = new TCanvas("c1","A Simple Graph Example",200,10,1200,500);
  c1->Divide(2,1);
  vector<double> vNmax[N], vdmax[N], veNmax[N], vedmax[N], vd_center;
  TGraphErrors* gNmax[N];
  TGraphErrors* gdmax[N];
  vector<double> slopes1, slopes2, eslopes1, eslopes2, angles;
  TF1* f1[N];
  TF1* f2[N];
  for(int j=0; j < N; j++){  
    for(int k=0; k < dim; k++){
      vNmax[j].push_back(p1[k].at(j));
      vdmax[j].push_back(p2[k].at(j));
      veNmax[j].push_back(ep1[k].at(j));
      vedmax[j].push_back(ep2[k].at(j));
      if(j==0)
	vd_center.push_back(d_center.at(j));     
    }
    gNmax[j] = new TGraphErrors(vd_center.size(), &d_center[0], &vNmax[j][0], 0, &veNmax[j][0]);
    gdmax[j] = new TGraphErrors(vd_center.size(), &d_center[0], &vdmax[j][0], 0, &vedmax[j][0]);
    f1[j] =  new TF1("f1",pol1,0.,range_d,2);
    f1[j]->SetLineColor(1 + j);
    f1[j]->SetLineStyle(kDashed);
    f2[j] =  new TF1("f2",pol1,0.,range_d,2);
    f2[j]->SetLineColor(1 + j);
    f2[j]->SetLineStyle(kDashed);    
    gNmax[j]->Fit(f1[j],"Q","",0,range_d);
    gdmax[j]->Fit(f2[j],"Q","",0,range_d);
  
    slopes1.push_back(f1[j]->GetParameter(1));
    slopes2.push_back(f2[j]->GetParameter(1));   
    const double* tmp1 = f1[j]->GetParErrors();
    const double* tmp2 = f2[j]->GetParErrors();
    eslopes1.push_back(tmp1[1]);  
    eslopes2.push_back(tmp2[1]);
    angles.push_back(theta[j]);
    
    gNmax[j]->SetLineColor(1 + j);
    gNmax[j]->SetMarkerColor(1 + j);
    gNmax[j]->SetMarkerStyle(20 + j);
    gdmax[j]->SetLineColor(1 + j);
    gdmax[j]->SetMarkerColor(1 + j);
    gdmax[j]->SetMarkerStyle(20 + j);
    if(j==0) {
      gNmax[j]->SetTitle(0);
      gNmax[j]->GetXaxis()->SetTitle("distance to center [cm]");
      gNmax[j]->GetYaxis()->SetTitle("N_{max} Gaisser-Hillas");
      gNmax[j]->GetYaxis()->SetRangeUser(0,2.0);
      gdmax[j]->SetTitle(0);
      gdmax[j]->GetXaxis()->SetTitle("distance to center [cm]");
      gdmax[j]->GetYaxis()->SetTitle("d_{max} Gaisser-Hillas");
      gdmax[j]->GetYaxis()->SetRangeUser(-50,300);
      c1->cd(1);
      gNmax[j]->Draw("ap");
      c1->cd(2);
      gdmax[j]->Draw("ap");
    }  
    else {
      c1->cd(1);
      gNmax[j]->Draw("p same");
      c1->cd(2);
      gdmax[j]->Draw("p same");
    }
    c1->cd(1);
    f1[j]->Draw("l same");
    c1->cd(2);
    f2[j]->Draw("l same");
  }


  TGraphErrors* g1= new TGraphErrors(angles.size(),&angles[0],&slopes1[0],0, &eslopes1[0]);
  TF1* p0_m1 =  new TF1("p0_m1",pol0,0.,90,1);
  p0_m1->SetLineColor(kMagenta);
  p0_m1->SetLineStyle(kDashed);
  TGraphErrors* g2= new TGraphErrors(angles.size(),&angles[0],&slopes2[0],0, &eslopes2[0]);
  TF1* p0_m2 =  new TF1("p0_m2",pol0,0.,90,1);
  p0_m2->SetLineColor(kMagenta);
  p0_m2->SetLineStyle(kDashed);
  TCanvas *c3 = new TCanvas("c3","A Simple Graph Example",200,200,1000,500);
  c3->Divide(2,1);
  c3->cd(1);
  g1->SetTitle(0);
  g1->GetXaxis()->SetTitle("#theta [deg]");
  g1->GetYaxis()->SetTitle("slope N_{max} Gaisser-Hillas");
  g1->SetMarkerStyle(20);
  g1->Fit(p0_m1,"Q","",0,90);
  double m1 = p0_m1->GetParameter(0);
  g1->Draw("ap");
  c3->cd(2);
  g2->SetTitle(0);
  g2->GetXaxis()->SetTitle("#theta [deg]");
  g2->GetYaxis()->SetTitle("slope d_{max} Gaisser-Hillas");
  g2->SetMarkerStyle(20);
  g2->Fit(p0_m2,"Q","",0,90);
  double m2 = p0_m2->GetParameter(0);
  g2->Draw("ap");

  // Bringing the corrections to the Y-Z center
  double pGH[9][4];
  for(int i = 0; i < 9; i++) {
    pGH[i][0]= p1[0].at(i) + m1 * (0 - d_center.at(0));
    pGH[i][1]= p2[0].at(i) + m2 * (0 - d_center.at(0));
    pGH[i][2]= p3[0].at(i);
    pGH[i][3]= p4[0].at(i);
  }
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"Corrections to plug into LArSoft PhotonVisibilityServices"<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"  GH p1: ";
  for(int i = 0; i < 9; i++)
    cout<<pGH[i][0]<<", ";
  cout<<""<<endl;
  cout<<"  GH p2: ";
  for(int i = 0; i < 9; i++)
    cout<<pGH[i][1]<<", ";
  cout<<""<<endl;
  cout<<"  GH p3: ";
  for(int i = 0; i < 9; i++)
    cout<<pGH[i][2]<<", ";
  cout<<""<<endl;
  cout<<"  GH p4: ";
  for(int i = 0; i < 9; i++)
    cout<<pGH[i][3]<<", ";
  cout<<""<<endl;
  
  cout<<"m1 = "<<m1<<endl;
  cout<<"m2 = "<<m2<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  cout<<"-------------------------------------------------------"<<endl;
  
}
