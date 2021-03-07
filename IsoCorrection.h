using namespace std;

void rhoCorrection(Int_t type, Float_t* typeIso, TreeReader &data, vector<Int_t>& realpho, vector<Float_t>& rcIso){
  
  Float_t* phoEta = data.GetPtrFloat("phoEta");
  //Int_t    nPho = data.GetPtrInt("nPho");
  Float_t  rho                    = data.GetFloat("rho");

  
    
  Int_t npho = realpho.size();
  Int_t iEA = 0;
  
  //[type ch=0, pho=1, nh=2][eta section : iEA]
  Float_t EA[3][7] = {{0.0112, 0.0108, 0.0106, 0.01002, 0.0098, 0.0089, 0.000087},
		      {0.0668, 0.1054, 0.0786, 0.0223, 0.0078, 0.0028, 0.0137},
		      {0.1113, 0.0953, 0.0619, 0.0837, 0.1070, 0.1212, 0.1446}
  };

  for(Int_t ii=0; ii<npho; ii++){
    Int_t ipho = realpho[ii];
    if(fabs(phoEta[ipho] < 1.0)) iEA = 0;
    else if(fabs(phoEta[ipho]) > 1.0 && fabs(phoEta[ipho]) < 1.479) iEA = 1;
    else if(fabs(phoEta[ipho]) > 1.479 && fabs(phoEta[ipho]) < 2.0) iEA = 2;
    else if(fabs(phoEta[ipho]) > 2.0 && fabs(phoEta[ipho]) < 2.2) iEA = 3;
    else if(fabs(phoEta[ipho]) > 2.2 && fabs(phoEta[ipho]) < 2.3) iEA = 4;
    else if(fabs(phoEta[ipho]) > 2.3 && fabs(phoEta[ipho]) < 2.4) iEA = 5;
    else if(fabs(phoEta[ipho]) > 2.4) iEA = 6;

    rcIso.push_back(typeIso[ipho] - rho*EA[type][iEA]);
  }
}
