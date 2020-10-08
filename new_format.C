{

  //lista.txt is the list of the path/output_files we have generated for our study
  //An example line would be:
  ///pnfs/sbnd/scratch/users/gamez/v08_34_00/gen/semi_mode/28205516_0/372_scint.root

  ifstream Traks_file1("./lista.txt");
  if(!Traks_file1) cerr << "WARNING:  Failed to open file 1"<< endl;
  Traks_file1.seekg(0);

  int n = 0;
  // Ignore first n lines
  string line;
  for (int i=0;i<n;i++){
    getline(Traks_file1, line);
  }  
   
  string files;
  while(Traks_file1 >> files) {
    // The next two lines need to be checked by the users to verify they make sense in their cases 
    std::size_t first = files.find("_0/");
    std::size_t last = files.find("_sc");

    //cout<<first<<"  "<<last<<endl;
    std::string id_number = files.substr(first + 3,last - first - 3);

    int n_new = stoi(id_number);

    string names = "v04_tree";
 
    char inFile[120];
    sprintf(inFile,"%s",files.c_str());    
    //    cout<<"I: "<<inFile<<endl;

    char outFile[120];
    sprintf(outFile,"%s_%i.root",names.c_str(),n_new);    
    //n++;

    cout<<"Ouput name: "<<outFile<<endl;
    //Get old file, old tree and set top branch address
        TFile *oldfile = new TFile(inFile);
    
    TTree *oldtree = (TTree*)oldfile->Get("pmtresponse/AllPhotons");
    Int_t EventID, OpChannel;
    Float_t Wavelength, Time;
    oldtree->SetBranchAddress("EventID", &EventID);
    oldtree->SetBranchAddress("OpChannel", &OpChannel);
    oldtree->SetBranchAddress("Wavelength", &Wavelength);
    oldtree->SetBranchAddress("Time", &Time);
    
    TTree  *oldtree2 = (TTree *)oldfile->Get("generator/PhotonsGenerated");
    double X, Y, Z, T, PX, PY, PZ, PT;
    int Event_ID;
    oldtree2->SetBranchAddress("X", &X);
    oldtree2->SetBranchAddress("Y", &Y);
    oldtree2->SetBranchAddress("Z", &Z);
    oldtree2->SetBranchAddress("T", &T);
    oldtree2->SetBranchAddress("PX", &PX);
    oldtree2->SetBranchAddress("PY", &PY);
    oldtree2->SetBranchAddress("PZ", &PZ);
    oldtree2->SetBranchAddress("PT", &PT);
    oldtree2->SetBranchAddress("EventID", &Event_ID);
    oldtree2->SetBranchStatus("*",0);
    oldtree2->SetBranchStatus("X",1);
    oldtree2->SetBranchStatus("Y",1);
    oldtree2->SetBranchStatus("Z",1);
      
    //Create a new file + a clone of old trees in new file
    TFile *newfile = new TFile(outFile,"recreate");
    
    TTree *newtree = oldtree->CloneTree();
    
    TTree *newtree2 = oldtree2->CloneTree(0);
    int entries;
    TBranch *br = newtree2->Branch("entries", &entries, "entries/I");
    Long64_t nentries = oldtree2->GetEntries();
    for (Long64_t i=0;i<nentries; i++) {
      oldtree2->GetEntry(i);
      //newbranch
      entries = nentries;
      br->Fill();
      if (i < 1) {newtree2->Fill(); break;} 
    }
    
    
    newfile->cd();
    newtree->AutoSave();
    newtree2->AutoSave();

    delete oldfile;
    delete newfile;
    
  }
   
}
