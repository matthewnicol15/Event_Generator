#include <TGenPhaseSpace.h>
#include <TComplex.h>
#include <TVector3.h>


void d2380_DeltaDelta_Width_Input()
{

  // Output file path and name
  TFile fileOutput1("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/d2380_Width_DeltaDelta_Mol_Output_1B_21022022_01.root","recreate");
  // TFile fileOutput1("/media/mn688/Elements1/PhD/Event_Generator/Output/d2380_Width_DeltaDelta_Output_10M_21022022_01.root","recreate");

  /////////////////////////////////////////////////////////////////////////////////
  /////  Create the histograms     //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  // Width against invariant mass of decuplet particles. Check gamma values correct
  TH2F* hInv_Delta_PP_Width=new TH2F("hInv_Delta_PP_Width","",200,1.0,1.5,200,0,0.2);
  hInv_Delta_PP_Width->GetXaxis()->SetNdivisions(505);
  hInv_Delta_PP_Width->GetXaxis()->SetLabelSize(0.05);
  hInv_Delta_PP_Width->GetXaxis()->SetTitleSize(0.05);
  hInv_Delta_PP_Width->GetYaxis()->SetNdivisions(505);
  hInv_Delta_PP_Width->GetYaxis()->SetLabelSize(0.05);
  hInv_Delta_PP_Width->GetYaxis()->SetTitleSize(0.05);

  TH2F* hInv_Delta_M_Width=new TH2F("hInv_Delta_M_Width","",200,1.0,1.5,200,0,0.2);
  hInv_Delta_M_Width->GetXaxis()->SetNdivisions(505);
  hInv_Delta_M_Width->GetXaxis()->SetLabelSize(0.05);
  hInv_Delta_M_Width->GetXaxis()->SetTitleSize(0.05);
  hInv_Delta_M_Width->GetYaxis()->SetNdivisions(505);
  hInv_Delta_M_Width->GetYaxis()->SetLabelSize(0.05);
  hInv_Delta_M_Width->GetYaxis()->SetTitleSize(0.05);

  // Creating invariant mass histograms
  TH1F* hInv_d2380_M_WW_FF=new TH1F("hInv_d2380_M_WW_FF","Invariant mass of p, #pi^{+}, n and #pi^{-};M(p #pi^{+} n #pi^{-}) [GeV]; Width",3000,2.00062,5);
  hInv_d2380_M_WW_FF->GetXaxis()->SetNdivisions(505);
  hInv_d2380_M_WW_FF->GetXaxis()->SetLabelSize(0.05);
  hInv_d2380_M_WW_FF->GetXaxis()->SetTitleSize(0.05);
  hInv_d2380_M_WW_FF->GetYaxis()->SetNdivisions(505);
  hInv_d2380_M_WW_FF->GetYaxis()->SetLabelSize(0.05);
  hInv_d2380_M_WW_FF->GetYaxis()->SetTitleSize(0.05);

  TH1F* hInv_d2380_M_Phs_WW_FF_noq=new TH1F("hInv_d2380_M_Phs_WW_FF_noq","Invariant mass of p, #pi^{+}, n and #pi^{-};M(p #pi^{+} n #pi^{-}) [GeV];Counts;Width",3000,2.00062,5);
  hInv_d2380_M_Phs_WW_FF_noq->GetXaxis()->SetNdivisions(505);
  hInv_d2380_M_Phs_WW_FF_noq->GetXaxis()->SetLabelSize(0.05);
  hInv_d2380_M_Phs_WW_FF_noq->GetXaxis()->SetTitleSize(0.05);
  hInv_d2380_M_Phs_WW_FF_noq->GetYaxis()->SetNdivisions(505);
  hInv_d2380_M_Phs_WW_FF_noq->GetYaxis()->SetLabelSize(0.05);
  hInv_d2380_M_Phs_WW_FF_noq->GetYaxis()->SetTitleSize(0.05);

  // Histogram of the width vs the binding energy of d2380
  TH1F* hBinding_d2380_M_WW_FF=new TH1F("hBinding_d2380_M_WW_FF","Binding energy of d2380;Binding Energy [GeV]; Width",5000,-1.00062,4);
  hBinding_d2380_M_WW_FF->GetXaxis()->SetNdivisions(505);
  hBinding_d2380_M_WW_FF->GetXaxis()->SetLabelSize(0.05);
  hBinding_d2380_M_WW_FF->GetXaxis()->SetTitleSize(0.05);
  hBinding_d2380_M_WW_FF->GetYaxis()->SetNdivisions(505);
  hBinding_d2380_M_WW_FF->GetYaxis()->SetLabelSize(0.05);
  hBinding_d2380_M_WW_FF->GetYaxis()->SetTitleSize(0.05);
  hInv_Delta_PP_Width->GetXaxis()->SetNdivisions(505);
  hInv_Delta_PP_Width->GetXaxis()->SetLabelSize(0.05);
  hInv_Delta_PP_Width->GetXaxis()->SetTitleSize(0.05);
  hInv_Delta_PP_Width->GetYaxis()->SetNdivisions(505);
  hInv_Delta_PP_Width->GetYaxis()->SetLabelSize(0.05);
  hInv_Delta_PP_Width->GetYaxis()->SetTitleSize(0.05);
  // Histogram of the width vs the binding energy of d2380
  TH1F* hBinding_d2380_M_Phs_WW_FF_noq=new TH1F("hBinding_d2380_M_Phs_WW_FF_noq","Binding energy of d2380;Binding Energy [GeV]; Width",5000,-1.00062,4);
  hBinding_d2380_M_Phs_WW_FF_noq->GetXaxis()->SetNdivisions(505);
  hBinding_d2380_M_Phs_WW_FF_noq->GetXaxis()->SetNdivisions(505);
  hBinding_d2380_M_Phs_WW_FF_noq->GetXaxis()->SetLabelSize(0.05);
  hBinding_d2380_M_Phs_WW_FF_noq->GetXaxis()->SetTitleSize(0.05);
  hBinding_d2380_M_Phs_WW_FF_noq->GetYaxis()->SetNdivisions(505);
  hBinding_d2380_M_Phs_WW_FF_noq->GetYaxis()->SetLabelSize(0.05);
  hBinding_d2380_M_Phs_WW_FF_noq->GetYaxis()->SetTitleSize(0.05);

  /////////////////////////////////////////////////////////////////////////////////
  /////  Create the particles     //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  //TLorentzVectors for the particles involved in the decay
  // Initial particle
  TLorentzVector d2380;
  // Particles from 1st decay vertex
  TLorentzVector *Delta_PP_Proton, *Delta_PP_Pip, *Delta_M_Neutron, *Delta_M_Pim;
  // TLorentzVectors for the Delta invariant masses
  TLorentzVector Delta_PP_Inv, Delta_M_Inv, Inv_d2380;
  //Binding Energies
  Double_t Binding_d2380;
  // TLorentzVectors for boosted pions in rest frame of Deltas
  TLorentzVector boosted_pim, boosted_pip;

  // Masses of particles
  Double_t Masses_1[4] = {0.93827, 0.13957, 0.93957, 0.13957}; // Masses of particles from 1st decay vertex
  Double_t m_d2380; // Mass of the d*2380 to be randomly generated later
  TRandom *random = new TRandom3();

  Double_t Threshold; // Sum of decay particles
  // Adds up the decay particle rest masses to get threshold
  for(Int_t j=0; j<4; j++){
    Threshold += Masses_1[j];
  }
  cout<<Threshold<<endl;

  /////////////////////////////////////////////////////////////////////////////////
  /////  Produce components for weights here    //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  // Phase weight
  Double_t Phasespace_Weight_1;
  // Momentum of pions in rest frame of parent Delta
  Double_t q_Pim, q_Pip, q_2del;
  // Breit-Wigner
  Double_t Delta_PP_Width, Delta_M_Width;
  // Nominal mass of Deltas
  Double_t m_Delta = 1.232000;
  Double_t R = 6.3; // s-channel resonance in the pn and SigmaCascade systems
  Double_t gamma_R = 100.012; // This contains the coupling constants and other constants,
  // fitted to yield a total width equal to binding energy
  Double_t Lambda = 0.16; // Cutoff parameter, this is adjusted to best reproduce the ABC effect (decup-decup only)
  Double_t FF; // Monopole formfactor, this accounts for potential barriers
  Double_t Width_Weight_noq, Width_Weight; // Overall weight to apply for the width
  // with (decup-decup only) and without q factor

  // TVector3s of Deltas for boosting
  TVector3 Delta_PP;
  TVector3 Delta_M;

  // Complex part for the propogators
  TComplex comp_P, comp_M;
  // Propogators for the two Deltas
  TComplex D_Delta_PP, D_Delta_M;
  /////////////////////////////////////////////////////////////////////////////////
  /////  Simulate events    //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////


  //How many events to simulate and percentage completed
  Int_t nevents=1000000000;
  Int_t Percentage=nevents/100;

  //Creating the different decay vertices
  TGenPhaseSpace Vertex_1;

  //Loop over simulated events
  for (Long64_t i=0;i<nevents;){

    //Counter, shows percentage of simulated events completed
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage);
      fflush (stderr);
    }
    // Random number for the d2380 mass
    m_d2380 = random->Uniform(Threshold,Threshold + 2);
    // Set the four-vector for the d2380
    d2380.SetPxPyPzE(0,0,0,m_d2380);

    //Set decay for 1st vertex and phasespace weights
    Vertex_1.SetDecay(d2380,4,Masses_1,"Fermi");
    Phasespace_Weight_1=Vertex_1.Generate();
    Delta_PP_Proton=Vertex_1.GetDecay(0);
    Delta_PP_Pip=Vertex_1.GetDecay(1);
    Delta_M_Neutron=Vertex_1.GetDecay(2);
    Delta_M_Pim=Vertex_1.GetDecay(3);

    Delta_PP_Inv = (TLorentzVector)*Delta_PP_Proton + (TLorentzVector)*Delta_PP_Pip;
    Delta_M_Inv = (TLorentzVector)*Delta_M_Neutron + (TLorentzVector)*Delta_M_Pim;

    Delta_PP.SetXYZ(Delta_PP_Inv.Px()/Delta_PP_Inv.E(),Delta_PP_Inv.Py()/Delta_PP_Inv.E(),Delta_PP_Inv.Pz()/Delta_PP_Inv.E());
    Delta_M.SetXYZ(Delta_M_Inv.Px()/Delta_M_Inv.E(),Delta_M_Inv.Py()/Delta_M_Inv.E(),Delta_M_Inv.Pz()/Delta_M_Inv.E());

    Inv_d2380 = (TLorentzVector)*Delta_PP_Proton + (TLorentzVector)*Delta_PP_Pip + (TLorentzVector)*Delta_M_Neutron + (TLorentzVector)*Delta_M_Pim;
    // Produce binding Energy
    Binding_d2380 = Inv_d2380.M() - (2*m_Delta);


    boosted_pip.SetPxPyPzE(Delta_PP_Pip->Px(), Delta_PP_Pip->Py(),Delta_PP_Pip->Pz(),Delta_PP_Pip->E());
    boosted_pim.SetPxPyPzE(Delta_M_Pim->Px(), Delta_M_Pim->Py(),Delta_M_Pim->Pz(),Delta_M_Pim->E());
    // Boosting pions in rest frame of parent particles
    boosted_pip.Boost(-Delta_PP);
    boosted_pim.Boost(-Delta_M);

    // Getting momentum of pions in rest frame of parent particles
    q_Pip = boosted_pip.Rho();
    q_Pim = boosted_pim.Rho();
    q_2del = 2*Delta_PP_Inv.Rho();

    Delta_PP_Width = (0.74*(pow(q_Pip,3))*pow(R,2)/(1 + pow(R,2)*pow(q_Pip,2)));
    Delta_M_Width = (0.74*(pow(q_Pim,3))*pow(R,2)/(1 + pow(R,2)*pow(q_Pim,2)));

    // Complex part of the propogators for Deltas
    comp_P=TComplex(0,m_Delta*Delta_PP_Width);
    comp_M=TComplex(0,m_Delta*Delta_M_Width);

    // Propogators for the two Deltas
    D_Delta_PP = (sqrt((m_Delta*Delta_PP_Width)/q_Pip))*(1.0/(pow(Delta_PP_Inv.M(),2) - pow(m_Delta,2) + comp_P));
    D_Delta_M = (sqrt((m_Delta*Delta_M_Width)/q_Pim)*(1.0/(pow(Delta_M_Inv.M(),2) - pow(m_Delta,2) + comp_M)));

    // Defining monopole formfactor
    // (octet-octet) use R in equation 3 in paper
    // (decup-decup) use Lambda using equation 5 in paper
    FF = (pow(Lambda,2)) / (pow(Lambda,2)+((pow(q_2del,2))/4));

    // Defining total weight equation for width
    Width_Weight_noq = pow(FF,2)*gamma_R*(D_Delta_PP*D_Delta_M).Rho2();

    hInv_Delta_PP_Width->Fill(Delta_PP_Inv.M(),Delta_PP_Width);
    hInv_Delta_M_Width->Fill(Delta_M_Inv.M(),Delta_M_Width);

    // Binding energies and invariant mass of d2380 with various weights applied
    hInv_d2380_M_Phs_WW_FF_noq->Fill(Inv_d2380.M(), Phasespace_Weight_1*Width_Weight_noq);
    hBinding_d2380_M_Phs_WW_FF_noq->Fill(Binding_d2380, Phasespace_Weight_1*Width_Weight_noq);


    i++;
  }

  // Normalise for number of events
  Double_t scale = 10.0/nevents;
  hInv_d2380_M_WW_FF->Scale(scale);
  hInv_d2380_M_Phs_WW_FF_noq->Scale(scale);
  hBinding_d2380_M_WW_FF->Scale(scale);
  hBinding_d2380_M_Phs_WW_FF_noq->Scale(scale);


  // Scale to the correct width of d2380->pn (8 MeV) for 84 MeV binding energy
  Double_t Before_Scaling = hBinding_d2380_M_Phs_WW_FF_noq->GetBinContent(hBinding_d2380_M_Phs_WW_FF_noq->FindBin(-0.084));
  Double_t Scale_Factor = 0.063 / Before_Scaling;
  cout<<Scale_Factor<<endl;
  // hInv_d2380_M_Phs_WW_FF_noq->Scale(Scale_Factor);
  // hBinding_d2380_M_Phs_WW_FF_noq->Scale(Scale_Factor);

  // Print out the width in MeV
  Double_t Width = 1000*hBinding_d2380_M_Phs_WW_FF_noq->GetBinContent(hBinding_d2380_M_Phs_WW_FF_noq->FindBin(-0.084));
  cout<<"Width = "<<Width<<" MeV"<<endl;

  // Scaling Factor determined to be 100.012 for decup-decup

  fileOutput1.Write();
  fileOutput1.Close();
}
