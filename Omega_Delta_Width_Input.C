#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TVirtualPad.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TRandom3.h>
#include <TComplex.h>

void Omega_Delta_Width_Input(){
  // Output file path and name
  TFile fileOutput1("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/Omega_M_Delta_M_Mol_Width_1B_Output_21022022_01.root","recreate");
  // TFile fileOutput1("../Output/Omega_M_Delta_M_Width_Output_21022022_01.root","recreate");

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Creating histograms here
  // Histogram of the width vs the mass of dsss for Lambda
  TH1F* hInv_dsss_M_Phs_WW_FF_noq=new TH1F("hInv_dsss_M_Phs_WW_FF_noq","Invariant mass of #Delta^{-}, n and #pi^{-}; M(#Delta^{-} n #pi^{-} [GeV]; Width",3000,2.00062,5);
  hInv_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetNdivisions(505);
  hInv_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetLabelSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetTitleSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetNdivisions(505);
  hInv_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetLabelSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetTitleSize(0.05);

  // Histogram of the width vs the binding energy of dsss for Lambda
  TH1F* hBinding_dsss_M_Phs_WW_FF_noq=new TH1F("hBinding_dsss_M_Phs_WW_FF_noq","Binding energy of dsss;Binding Energy [GeV]; Width",5000,-1.00062,4);
  hBinding_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetLabelSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetTitleSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetLabelSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetTitleSize(0.05);


  // Histogram of the invariant mass of Delta^{-}
  TH2F* hInv_Delta_M=new TH2F("hInv_Delta_M","Invariant mass of n and #pi^{-};M(n #pi^{-}) [GeV];Width",200,1,2,1000,0,1);
  hInv_Delta_M->GetXaxis()->SetNdivisions(505);
  hInv_Delta_M->GetXaxis()->SetLabelSize(0.05);
  hInv_Delta_M->GetXaxis()->SetTitleSize(0.05);
  hInv_Delta_M->GetYaxis()->SetNdivisions(505);
  hInv_Delta_M->GetYaxis()->SetLabelSize(0.05);
  hInv_Delta_M->GetYaxis()->SetTitleSize(0.05);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Produce particles here

  // Produce Vectors for particles here
  // Initial particle
  TLorentzVector dsss;
  // Particles from 1st decay vertex
  TLorentzVector *Omega_M, *Delta_M_Neutron, *Delta_M_Pim;
  // Invariant masses
  TLorentzVector Inv_Delta_M, Omega_Minus, Inv_dsss;
  //Binding Energies
  Double_t Binding_dsss;
  // Boosted pions in rest frame of Sigma^{*-} and Cascade^{*-}
  TLorentzVector Boosted_Delta_M_Pim;
  // TVector3s of Cascade^{*-} and Sigma^{*-}
  TVector3 Delta_M;

  // Masses of particles
  Double_t Masses_1[3] = {1.6725, 0.93957, 0.13957}; // Masses of particles from 1st decay vertex
  Double_t m_dsss; // Mass of dsss, to be randomly generated in event loop
  TRandom *random = new TRandom3(); // Random number generator to produce mass of dsss

  Double_t Threshold; // Sum of decay particles

  for(Int_t j=0; j<3; j++){
    Threshold += Masses_1[j];
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Produce components for weights here
  Double_t Phasespace_1; // Phase weight for decay vertices
  Double_t E_Delta_M_Pim; // Energy of pions in rest frame of parent particles
  Double_t q_Delta_M_Pim; // Momentum of pions in rest frame of parent particles
  Double_t q_Omega_Delta; // Combined momentum of Omega^{-} and Delta^{-}
  Double_t Delta_M_Width; // Breit Wigner widths
  Double_t m_Omega_M = 1.6725, m_Delta_M = 1.2320; // Nominal masses of parent particles
  Double_t gamma_BW_Delta = 0.74; // Gamma factor for Breit wigner weights
  Double_t R = 6.3; // s-channel resonance in the pn and SigmaCascade systems
  TComplex comp_Delta_M; // Complex part for the propogators
  TComplex D_Delta_M; // Propogators for parent particles with complex components
  // Double_t gamma_0 =  // Sum of decay widths of daughter particles
  Double_t gamma_R = 100.012; // This contains the coupling constants and other constants, fitted to yield a total width equal to binding energy
  Double_t Lambda = 0.16; // Cutoff parameter, this is adjusted to best reproduce the ABC effect
  Double_t FF; // Monopole formfactor, this accounts for potential barriers
  Double_t Width_Weight_noq; // Overall weight to apply for the width without q factor


  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

    // Randomly generate the mass of dsss from nominal masses of products up to an arbitrary number
    m_dsss = random->Uniform(Threshold,Threshold + 2);
    // Set the TLorentzVector for dsss now the mass has been generated
    dsss.SetPxPyPzE(0,0,0,m_dsss);

    //Set decay for 1st vertex
    Vertex_1.SetDecay(dsss,3,Masses_1,"Fermi");
    // Get the phasespace weight for decay vertex 1
    Phasespace_1=Vertex_1.Generate();
    // Assign the decay particles to the appropriate TLorentzVectors
    // These need to be in the same order as the array of masses
    Omega_M=Vertex_1.GetDecay(0);
    Delta_M_Neutron=Vertex_1.GetDecay(1);
    Delta_M_Pim=Vertex_1.GetDecay(2);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Produce invariant masses here
    Inv_Delta_M = (TLorentzVector)*Delta_M_Neutron + (TLorentzVector)*Delta_M_Pim; // Invariant mass for Delta^{-}
    Omega_Minus = (TLorentzVector)*Omega_M; // Non pointer TLorentzVector for Omega^{-}
    Inv_dsss = (TLorentzVector)*Delta_M_Neutron + (TLorentzVector)*Delta_M_Pim + (TLorentzVector)*Omega_M; // Invariant mass for dsss^{--}

    // Produce binding Energy
    Binding_dsss = Inv_dsss.M() - m_Omega_M - m_Delta_M;


    // Getting the momentum and energy of the pions for width equations
    // Delta Minus Pion
    E_Delta_M_Pim = (pow(Inv_Delta_M.M(),2) + pow(Masses_1[2],2) - pow(Masses_1[1],2))/(2*Inv_Delta_M.M());
    q_Delta_M_Pim = sqrt(abs(pow(E_Delta_M_Pim,2) - pow(Masses_1[2],2)));


    // Getting the momentum of the two parent particles
    q_Omega_Delta = Inv_Delta_M.Rho() + Omega_Minus.Rho();

    // Defining the Breit Wigner weights
    Delta_M_Width = (gamma_BW_Delta * (pow(q_Delta_M_Pim,3)) * pow(R,2) / (1 + pow(R,2) * pow(q_Delta_M_Pim,2)));

    // Defining the complex component of the propogators
    comp_Delta_M = TComplex(0,m_Delta_M * Delta_M_Width);

    // Defining propogator for Delta^{-}
    D_Delta_M = (sqrt((m_Delta_M*Delta_M_Width)/q_Delta_M_Pim))*(1.0/(pow(Inv_Delta_M.M(),2) - pow(m_Delta_M,2) + comp_Delta_M));

    // Defining monopole formfactor
    FF = (pow(Lambda,2)) / (pow(Lambda,2) + ((pow(q_Omega_Delta,2)) / 4));


    // Defining total weight equation for width
    Width_Weight_noq = /*Gamma_0 + */ pow(FF,2) * gamma_R * (D_Delta_M).Rho2();


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Filling histograms
    // Invariant mass of dsss with various weights applied
    hInv_dsss_M_Phs_WW_FF_noq->Fill(Inv_dsss.M(), Phasespace_1 * Width_Weight_noq); // All weights applied except phasespace and q factor

    // Binding energies of dss with various weights applied
    hBinding_dsss_M_Phs_WW_FF_noq->Fill(Binding_dsss, Phasespace_1 * Width_Weight_noq);

    i++;
  }

  // Determine scaling factor depending on number of events, this way histograms match no matter the statistics
  Double_t scale = 10.0 / nevents;

  // Scale histograms so they match if run with different number of events
  hInv_dsss_M_Phs_WW_FF_noq->Scale(scale);
  hBinding_dsss_M_Phs_WW_FF_noq->Scale(scale);

  // Print out the width in MeV
  Double_t Width = 1000*hBinding_dsss_M_Phs_WW_FF_noq->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq->FindBin(-0.084));
  cout<<"Width = "<<Width<<" MeV"<<endl;

  // Write the histograms to the output file
  fileOutput1.Write();
  // Close the output file
  fileOutput1.Close();
}
