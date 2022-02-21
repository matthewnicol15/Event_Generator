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

void Cascade_Star_Sigma_Star_Lambda_Width_Input(){
  // Output file path and name
  TFile fileOutput1("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/Cascade_Star_M_Sigma_Star_M_Lambda_Width_Mol_Output_1B_21022022_01.root","recreate");
  // TFile fileOutput1("../Output/Cascade_Star_M_Sigma_Star_M_Lambda_Width_Output_21022022_01.root","recreate");

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Creating histograms here
  // Histogram of the width vs the mass of dsss for Lambda
  TH1F* hInv_dsss_M_Phs_WW_FF_noq=new TH1F("hInv_dsss_M_Phs_WW_FF_noq","Invariant mass of 2p, 4#pi^{-} and #pi^{0};M(2p 4#pi^{-} #pi^{0}) [GeV]; Width",3000,2.00062,5);
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


  // Histogram of the invariant mass of Sigma^{*-}
  TH2F* hInv_Sigma_Star_M=new TH2F("hInv_Sigma_Star_M","Invariant mass of p and 2#pi^{-};M(p 2#pi^{-}) [GeV];Width",200,1,2,200,0.0,0.1);
  hInv_Sigma_Star_M->GetXaxis()->SetNdivisions(505);
  hInv_Sigma_Star_M->GetXaxis()->SetLabelSize(0.05);
  hInv_Sigma_Star_M->GetXaxis()->SetTitleSize(0.05);
  hInv_Sigma_Star_M->GetYaxis()->SetNdivisions(505);
  hInv_Sigma_Star_M->GetYaxis()->SetLabelSize(0.05);
  hInv_Sigma_Star_M->GetYaxis()->SetTitleSize(0.05);

  // Histogram of the invariant mass of Cascade^{*-}
  TH2F* hInv_Cascade_Star_M=new TH2F("hInv_Cascade_Star_M","Invariant mass of p, 2#pi^{-} and #pi^{0};M(p 2#pi^{-} #pi^{0}) [GeV];Width",200,1,2,200,0.0,0.1);
  hInv_Cascade_Star_M->GetXaxis()->SetNdivisions(505);
  hInv_Cascade_Star_M->GetXaxis()->SetLabelSize(0.05);
  hInv_Cascade_Star_M->GetXaxis()->SetTitleSize(0.05);
  hInv_Cascade_Star_M->GetYaxis()->SetNdivisions(505);
  hInv_Cascade_Star_M->GetYaxis()->SetLabelSize(0.05);
  hInv_Cascade_Star_M->GetYaxis()->SetTitleSize(0.05);

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Produce particles here

  // Produce Vectors for particles here
  // Initial particle
  TLorentzVector dsss;
  // Particles from 1st decay vertex
  TLorentzVector *Sigma_Star_M_Lambda, *Sigma_Star_M_Pim, *Cascade_Star_M_Cascade, *Cascade_Star_M_Pim;
  // Invariant masses
  TLorentzVector Inv_Sigma_Star_M, Inv_Cascade_Star_M, Inv_dsss;
  //Binding Energies
  Double_t Binding_dsss;
  // Boosted pions in rest frame of Sigma^{*-} and Cascade^{*-}
  TLorentzVector Boosted_Sigma_Star_M_Pim, Boosted_Cascade_Star_M_Pim;
  // TVector3s of Cascade^{*-} and Sigma^{*-}
  TVector3 Cascade_Star_M, Sigma_Star_M;

  // Masses of particles
  Double_t Masses_1[4] = {1.11568, 0.13957, 1.31486, 0.13957}; // Masses of particles from 1st decay vertex
  Double_t m_dsss; // Mass of dsss, to be randomly generated in event loop
  TRandom *random = new TRandom3(); // Random number generator to produce mass of dsss

  Double_t Threshold; // Sum of decay particles

  // This loop adds up the nominal masses of the array Masses_1 to give the
  // threshold mass
  for(Int_t j=0; j<4; j++){
    Threshold += Masses_1[j];
  }

  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Produce components for weights here
  Double_t Phasespace_1; // Phase weight for decay vertices
  Double_t E_Sigma_Star_M_Lambda_Pim, E_Sigma_Star_M_Sigma_Pim, E_Cascade_Star_M_Pim; // Energy of pions in rest frame of parent particles
  Double_t q_Sigma_Star_M_Lambda_Pim, q_Sigma_Star_M_Sigma_Pim, q_Cascade_Star_M_Pim; // Momentum of pions in rest frame of parent particles
  Double_t q_Sigma_Star_Cascade_Star; // Combined momentum of Sigma^{*-} and Cascade^{*-}
  Double_t Sigma_Star_M_Lambda_Width, Sigma_Star_M_Sigma_Width,Sigma_Star_M_Total_Width, Cascade_Star_M_Width; // Breit Wigner widths
  Double_t m_Sigma_Star_M = 1.3872, m_Cascade_Star_M = 1.5350; // Nominal masses of parent particles
  Double_t gamma_BW_Cascade_Star = 0.126, gamma_BW_Sigma_Star = 0.3; // Gamma factor for Breit wigner weights
  Double_t R = 6.3; // s-channel resonance in the pn and SigmaCascade systems
  TComplex comp_Sigma_Star_M, comp_Cascade_Star_M; // Complex part for the propogators
  TComplex D_Sigma_Star_M, D_Cascade_Star_M; // Propogators for parent particles with complex components
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
    Vertex_1.SetDecay(dsss,4,Masses_1,"Fermi");
    // Get the phasespace weight for decay vertex 1
    Phasespace_1=Vertex_1.Generate();
    // Assign the decay particles to the appropriate TLorentzVectors
    // These need to be in the same order as the array of masses
    Sigma_Star_M_Lambda=Vertex_1.GetDecay(0);
    Sigma_Star_M_Pim=Vertex_1.GetDecay(1);
    Cascade_Star_M_Cascade=Vertex_1.GetDecay(2);
    Cascade_Star_M_Pim=Vertex_1.GetDecay(3);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Produce invariant masses here
    Inv_Cascade_Star_M = (TLorentzVector)*Cascade_Star_M_Pim + (TLorentzVector)*Cascade_Star_M_Cascade; // Invariant mass for Cascade^{*-}
    Inv_Sigma_Star_M = (TLorentzVector)*Sigma_Star_M_Lambda + (TLorentzVector)*Sigma_Star_M_Pim; // Invariant mass for Sigma^{*-}
    Inv_dsss = (TLorentzVector)*Cascade_Star_M_Pim + (TLorentzVector)*Cascade_Star_M_Cascade + (TLorentzVector)*Sigma_Star_M_Lambda + (TLorentzVector)*Sigma_Star_M_Pim; // Invariant mass for dsss^{--}

    // Calculate binding Energy
    Binding_dsss = Inv_dsss.M() - m_Sigma_Star_M - m_Cascade_Star_M;


    // Getting the momentum and energy of the pions for width equations
    // Cascade Star Minus Pion
    E_Cascade_Star_M_Pim = (pow(Inv_Cascade_Star_M.M(),2) + pow(Masses_1[3],2) - pow(Masses_1[2],2))/(2*Inv_Cascade_Star_M.M());
    q_Cascade_Star_M_Pim = sqrt(abs(pow(E_Cascade_Star_M_Pim,2) - pow(Masses_1[3],2)));

    // Sigma Star Minus to Lambda Pion
    E_Sigma_Star_M_Lambda_Pim = (pow(Inv_Sigma_Star_M.M(),2) + pow(Masses_1[1],2) - pow(Masses_1[0],2))/(2*Inv_Sigma_Star_M.M());
    q_Sigma_Star_M_Lambda_Pim = sqrt(abs(pow(E_Sigma_Star_M_Lambda_Pim,2) - pow(Masses_1[3],2)));

    // Sigma Star Minus to Sigma Pion
    E_Sigma_Star_M_Sigma_Pim = (pow(Inv_Sigma_Star_M.M(),2) + pow(Masses_1[1],2) - pow(1.19745,2))/(2*Inv_Sigma_Star_M.M()); // Manually put in mass of sigma ground state because simulating lambda decay channel instead
    q_Sigma_Star_M_Sigma_Pim = sqrt(abs(pow(E_Sigma_Star_M_Sigma_Pim,2) - pow(Masses_1[1],2)));

    // Getting the momentum of the two parent particles
    q_Sigma_Star_Cascade_Star = Inv_Cascade_Star_M.Rho() + Inv_Sigma_Star_M.Rho();

    // Defining the Breit Wigner weights
    Cascade_Star_M_Width = (gamma_BW_Cascade_Star * (pow(q_Cascade_Star_M_Pim,3)) * R / (1 + R * pow(q_Cascade_Star_M_Pim,2)));
    Sigma_Star_M_Lambda_Width = (0.27 * (pow(q_Sigma_Star_M_Lambda_Pim,3)) * R / (1 + R * pow(q_Sigma_Star_M_Lambda_Pim,2)));
    Sigma_Star_M_Sigma_Width = (0.096 * (pow(q_Sigma_Star_M_Sigma_Pim,3)) * R / (1 + R * pow(q_Sigma_Star_M_Sigma_Pim,2)));
    Sigma_Star_M_Total_Width = Sigma_Star_M_Lambda_Width + Sigma_Star_M_Sigma_Width;

    // Defining the complex component of the propogators
    comp_Cascade_Star_M = TComplex(0,m_Cascade_Star_M * Cascade_Star_M_Width);
    comp_Sigma_Star_M = TComplex(0,m_Sigma_Star_M * Sigma_Star_M_Total_Width);

    // Defining propogators for Cascade^{*-} and Sigma^{*-}
    D_Cascade_Star_M = (sqrt((m_Cascade_Star_M*Cascade_Star_M_Width)/q_Cascade_Star_M_Pim))*(1.0/(pow(Inv_Cascade_Star_M.M(),2) - pow(m_Cascade_Star_M,2) + comp_Cascade_Star_M));
    D_Sigma_Star_M = (sqrt((m_Sigma_Star_M*Sigma_Star_M_Lambda_Width)/q_Sigma_Star_M_Lambda_Pim))*(1.0/(pow(Inv_Sigma_Star_M.M(),2) - pow(m_Sigma_Star_M,2) + comp_Sigma_Star_M));

    // Defining monopole formfactor
    FF = (pow(Lambda,2)) / (pow(Lambda,2) + ((pow(q_Sigma_Star_Cascade_Star,2)) / 4));


    // Defining total weight equation for width
    Width_Weight_noq = /*Gamma_0 + */ pow(FF,2) * gamma_R * (D_Sigma_Star_M * D_Cascade_Star_M).Rho2();


    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Filling histograms

    // Checking widths of 1st vertex decay particles
    hInv_Sigma_Star_M->Fill(Inv_Sigma_Star_M.M(),Sigma_Star_M_Total_Width);
    hInv_Cascade_Star_M->Fill(Inv_Cascade_Star_M.M(),Cascade_Star_M_Width);

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


  // Write the histograms to the output file
  fileOutput1.Write();
  // Close the output file
  fileOutput1.Close();
}
