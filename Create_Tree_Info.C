void Create_Tree_Info() {
  
  // *** User *** must write below the PATH/file_name with the list of files generated (probably in the Grid) containing the useful information
  string lista_files = "./List_of_Grid_files.txt";
  // *** User *** must write below the number of optical channels used in the list of files ("lista_files")
  const int numberDevices = 320;
  cout<<"----> numberDevices: "<<numberDevices<<endl;
  
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
  const int n_files = names.size() - 1;
  cout<<"----> number of files: "<<n_files<<endl;

  const Int_t kMaxDevices = 600;
  if(kMaxDevices < numberDevices) {
    cout<<"Warning: numberDevices larger than kMaxDevices, so change the last to continue!!!"<<endl;
    exit(0);
  }
 
 //Create a new tree file with all the information needed
  TFile *hfile = new TFile("MyNewFile.root","RECREATE");
  TTree *myTree = new TTree("myTree","A ROOT tree");
  int VUV_hits[kMaxDevices];
  int Vis_hits[kMaxDevices];
  double posX, posY, posZ;
  myTree->Branch("numberDevices",&numberDevices,"numberDevices/I");
  myTree->Branch("X", &posX, "X/D");
  myTree->Branch("Y", &posY, "Y/D");
  myTree->Branch("Z", &posZ, "Z/D");
  myTree->Branch("VUV_hits",VUV_hits,"VUV_hits[numberDevices]/I");
  myTree->Branch("Vis_hits",Vis_hits,"Vis_hits[numberDevices]/I");
  myTree->Branch("genPhotons", &genPhotons, "genPhotons/I");
  
  //loop over files
  for(int n=0; n<n_files; n++) {
   
    for(int i=0; i<numberDevices; i++) {
      VUV_hits[i] = 0;
      Vis_hits[i] = 0;
    }
    
    string input_file = names.at(n);
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"Name of the input file: "<<input_file<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    cout<<"-------------------------------------------------------"<<endl;
    TFile* f = new TFile(input_file.c_str());
    TTree *tree = (TTree *)f->Get("pmtresponse/AllPhotons");
    //tree->Print();   
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
    posX = X;
    posY = Y;
    posZ = Z;
    genPhotons = gen_photons;
  
    cout<<"x_value= "<<X<<"  y_value= "<<Y<<"  z_value= "<<Z<<endl;

    for(int i=0; i!=tree->GetEntries(); ++i)
      {
	tree->GetEntry(i);
	if (Wavelength < 300) 
	  VUV_hits[OpChannel]++;
	else
	  Vis_hits[OpChannel]++;
      }

    myTree->Fill();

    delete f;
  }//files

   hfile->cd();
   hfile->Write();
   hfile->Close();
  
}
