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


void Cascade_Star_Sigma_Star_Sigma_Width_Input(){

   TFile *fin=new TFile("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/dsss_Width_Cascade_Star_Sigma_Star_Lambda_Output_10M_10032022_01.root");

   auto *hInv_Sigma_Star  = (TH2F*) fin->Get("hInv_Sigma_Star");

   Int_t nevents=1000000000; // how many events to simulate
   ostringstream Quantity;
   Quantity << "1B"; // Define quantity for output file name

   ostringstream File_Path;
   ostringstream Channel;
   ostringstream Date;
   ostringstream Output_File_Name;
   File_Path << "/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/";
   Channel << "dsss_Width_Cascade_Star_Sigma_Star_Sigma_Output";
   Date << "10032022_01";
   Output_File_Name << File_Path.str().c_str()<<Channel.str().c_str()<<"_"<<Quantity.str().c_str()<<"_"<<Date.str().c_str()<<".root";
   cout<<Output_File_Name.str().c_str()<<endl;


   // Output file path and name
   TFile fileOutput1(Output_File_Name.str().c_str(),"recreate");
  // TFile fileOutput1("../Output/dsss_Width_Cascade_Star_Sigma_Star_To_Sigma_Output_1B_21022022_01.root","recreate");

  /////////////////////////////////////////////////////////////////////////////////
  /////  Create the histograms     //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  // Molecular case
  // Histogram of the width vs the mass of dsss for Lambda
  TH1F* hInv_dsss_M_Phs_WW_FF_noq_Mol=new TH1F("hInv_dsss_M_Phs_WW_FF_noq_Mol","Invariant mass of 2p, 4#pi^{-} and #pi^{0};M(2p 4#pi^{-} #pi^{0}) [GeV]; Width",3000,2.00062,5);
  hInv_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetNdivisions(505);
  hInv_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetLabelSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetTitleSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetNdivisions(505);
  hInv_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetLabelSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetTitleSize(0.05);

  // Histogram of the width vs the binding energy of dsss for Lambda
  TH1F* hBinding_dsss_M_Phs_WW_FF_noq_Mol=new TH1F("hBinding_dsss_M_Phs_WW_FF_noq_Mol","Binding energy of dsss;Binding Energy [GeV]; Width",5000,-1.00062,4);
  hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetLabelSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetTitleSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetLabelSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetTitleSize(0.05);

  // Hexaquark case
  // Histogram of the width vs the mass of dsss for Lambda
  TH1F* hInv_dsss_M_Phs_WW_FF_noq_Hex=new TH1F("hInv_dsss_M_Phs_WW_FF_noq_Hex","Invariant mass of 2p, 4#pi^{-} and #pi^{0};M(2p 4#pi^{-} #pi^{0}) [GeV]; Width",3000,2.00062,5);
  hInv_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetNdivisions(505);
  hInv_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetLabelSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetTitleSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetNdivisions(505);
  hInv_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetLabelSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetTitleSize(0.05);

  // Histogram of the width vs the binding energy of dsss for Lambda
  TH1F* hBinding_dsss_M_Phs_WW_FF_noq_Hex=new TH1F("hBinding_dsss_M_Phs_WW_FF_noq_Hex","Binding energy of dsss;Binding Energy [GeV]; Width",5000,-1.00062,4);
  hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetLabelSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetTitleSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetLabelSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetTitleSize(0.05);

  // Histogram of the invariant mass of Cascade^{*-}
  TH2F* hInv_Cascade_Star=new TH2F("hInv_Cascade_Star","Invariant mass of p, 2#pi^{-} and #pi^{0};M(p 2#pi^{-} #pi^{0}) [GeV];Width",1100,1.43,1.54,200,0.0,0.01);
  hInv_Cascade_Star->GetXaxis()->SetNdivisions(505);
  hInv_Cascade_Star->GetXaxis()->SetLabelSize(0.05);
  hInv_Cascade_Star->GetXaxis()->SetTitleSize(0.05);
  hInv_Cascade_Star->GetYaxis()->SetNdivisions(505);
  hInv_Cascade_Star->GetYaxis()->SetLabelSize(0.05);
  hInv_Cascade_Star->GetYaxis()->SetTitleSize(0.05);

  /////////////////////////////////////////////////////////////////////////////////
  /////  Create the particles     //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  // Produce Vectors for particles here
  // Initial particle
  TLorentzVector dsss;
  // Particles from 1st decay vertex
  TLorentzVector *Sigma_Star_Sigma, *Sigma_Star_Pim, *Cascade_Star_Cascade, *Cascade_Star_Pim;
  // Invariant masses of decuplet baryons and dibaryons
  TLorentzVector Inv_Sigma_Star, Inv_Cascade_Star, Inv_dsss;
  //Binding Energies
  Double_t Binding_dsss_Mol, Binding_dsss_Hex;

  // Masses of particles
  Double_t Masses_1[4] = {1.1894, 0.13498, 1.31486, 0.13498}; // Masses of particles from 1st decay vertex
  Double_t m_dsss; // Mass of dsss, to be randomly generated in event loop
  TRandom *random = new TRandom3(); // Random number generator to produce mass of dsss

  Double_t Threshold; // Sum of decay particles
  Threshold = 1.1894 + 0.13498 + 1.31486 + 0.13498;

  cout<<Threshold<<endl;

  /////////////////////////////////////////////////////////////////////////////////
  /////  Define constants, variables and variables for widths  ////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  Double_t Phasespace_Weight_1; // Phase weight for decay vertices
  Double_t E_Sigma_Star_Lambda_Pim, E_Sigma_Star_Sigma_Pim, E_Cascade_Star_Pim; // Energy of pions in rest frame of parent particles
  Double_t q_Sigma_Star_Lambda_Pim, q_Sigma_Star_Sigma_Pim, q_Cascade_Star_Pim; // Momentum of pions in rest frame of parent particles
  Double_t q_Sigma_Star_Cascade_Star; // Combined momentum of Sigma^{*-} and Cascade^{*-}
  Double_t Sigma_Star_Lambda_Width, Sigma_Star_Sigma_Width,Sigma_Star_Total_Width, Cascade_Star_Width; // Breit Wigner widths
  Double_t m_Omega = 1.6725, m_Delta = 1.2320, m_Cascade_Star = 1.5318,m_Sigma_Star = 1.3828; // Nominal masses of parent particles
  Double_t gamma_BW_Cascade_Star = 0.1162, gamma_BW_Sigma_Star_Sigma = 0.1090; // Gamma factor for Breit wigner weights
  Double_t R = 6.3; // s-channel resonance in the pn and SigmaCascade systems
  TComplex comp_Sigma_Star, comp_Cascade_Star; // Complex part for the propogators
  TComplex D_Sigma_Star_Sigma, D_Sigma_Star_Lambda, D_Cascade_Star; // Propogators for parent particles with complex components
  // Double_t gamma_0 =  // Sum of decay widths of daughter particles
  Double_t gamma_R_Mol = 100.012; // This contains the coupling constants and other constants, fitted to yield a total width equal to binding energy
  Double_t gamma_R_Hex = 100.012; // This contains the coupling constants and other constants, fitted to yield a total width equal to binding energy
  Double_t Lambda = 0.16; // Cutoff parameter, this is adjusted to best reproduce the ABC effect
  Double_t FF; // Monopole formfactor, this accounts for potential barriers
  Double_t Width_Weight_noq_Mol, Width_Weight_noq_Hex; // Overall weight to apply for the width without q factor


  /////////////////////////////////////////////////////////////////////////////////
  /////  Simulate events  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  //How many events to simulate and percentage completed
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
    Phasespace_Weight_1=Vertex_1.Generate();
    // Assign the decay particles to the appropriate TLorentzVectors
    // These need to be in the same order as the array of masses
    Sigma_Star_Sigma=Vertex_1.GetDecay(0);
    Sigma_Star_Pim=Vertex_1.GetDecay(1);
    Cascade_Star_Cascade=Vertex_1.GetDecay(2);
    Cascade_Star_Pim=Vertex_1.GetDecay(3);

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Produce invariant masses here
    Inv_Cascade_Star = (TLorentzVector)*Cascade_Star_Pim + (TLorentzVector)*Cascade_Star_Cascade; // Invariant mass for Cascade^{*-}
    Inv_Sigma_Star = (TLorentzVector)*Sigma_Star_Sigma + (TLorentzVector)*Sigma_Star_Pim; // Invariant mass for Sigma^{*-}
    Inv_dsss = (TLorentzVector)*Cascade_Star_Pim + (TLorentzVector)*Cascade_Star_Cascade + (TLorentzVector)*Sigma_Star_Sigma + (TLorentzVector)*Sigma_Star_Pim; // Invariant mass for dsss^{--}

    // Produce binding Energy
    // Produce binding Energy
    // (Hexaquark)): generated mass of dibaryon - rest mass from all decup-decup pairs (use table 7 in paper for ratios)
    // (Molecule): generated mass of dibaryon - rest mass of lightest decup pair
      Binding_dsss_Hex = m_dsss - ((m_Sigma_Star + m_Cascade_Star + m_Omega + m_Delta)/2);
      Binding_dsss_Mol = m_dsss - (m_Omega + m_Delta);


    // Getting the momentum and energy of the pions for width equations
    // Cascade Star Minus Pion
    E_Cascade_Star_Pim = (pow(Inv_Cascade_Star.M(),2) + pow(Masses_1[3],2) - pow(Masses_1[2],2))/(2*Inv_Cascade_Star.M());
    q_Cascade_Star_Pim = sqrt(fabs(pow(E_Cascade_Star_Pim,2) - pow(Masses_1[3],2)));

    // Sigma Star Minus to Sigma Pion
    E_Sigma_Star_Sigma_Pim = (pow(Inv_Sigma_Star.M(),2) + pow(Masses_1[1],2) - pow(Masses_1[0],2))/(2*Inv_Sigma_Star.M());
    q_Sigma_Star_Sigma_Pim = sqrt(fabs(pow(E_Sigma_Star_Sigma_Pim,2) - pow(Masses_1[1],2)));

    // Getting the momentum of the two parent particles
    q_Sigma_Star_Cascade_Star = Inv_Cascade_Star.Rho() + Inv_Sigma_Star.Rho();

    // Defining the Breit Wigner weights
    Cascade_Star_Width = (gamma_BW_Cascade_Star * (pow(q_Cascade_Star_Pim,3)) * pow(R,2)) / (1 + pow(R,2) * pow(q_Cascade_Star_Pim,2));
    Sigma_Star_Sigma_Width = (gamma_BW_Sigma_Star_Sigma * (pow(q_Sigma_Star_Sigma_Pim,3)) * pow(R,2)) / (1 + pow(R,2) * pow(q_Sigma_Star_Sigma_Pim,2));

    // Find the total sigma* width from the histogram
    TH1F *project = (TH1F*)hInv_Sigma_Star->ProjectionY("project",hInv_Sigma_Star->GetXaxis()->FindBin(Inv_Sigma_Star.M()));
    Sigma_Star_Total_Width = project->GetBinCenter(project->GetMaximumBin());
    // Delete the projection, otherwise they add together
    delete project;

    // Defining the complex component of the propogators
    comp_Cascade_Star = TComplex(0,m_Cascade_Star * Cascade_Star_Width);
    comp_Sigma_Star = TComplex(0,m_Sigma_Star * Sigma_Star_Total_Width);

    // Defining propogators for Cascade^{*-} and Sigma^{*-}
    D_Cascade_Star = (sqrt((m_Cascade_Star*Cascade_Star_Width)/q_Cascade_Star_Pim))*(1.0/(pow(Inv_Cascade_Star.M(),2) - pow(m_Cascade_Star,2) + comp_Cascade_Star));
    D_Sigma_Star_Sigma = (sqrt((m_Sigma_Star*Sigma_Star_Sigma_Width)/q_Sigma_Star_Sigma_Pim))*(1.0/(pow(Inv_Sigma_Star.M(),2) - pow(m_Sigma_Star,2) + comp_Sigma_Star));

    // Defining monopole formfactor
    FF = (pow(Lambda,2)) / (pow(Lambda,2) + ((pow(q_Sigma_Star_Cascade_Star,2)) / 4));


    // Defining total weight equation for width (use equation 4 in paper)
    Width_Weight_noq_Mol = /*Gamma_0 + */ pow(FF,2) * gamma_R_Mol * (D_Cascade_Star*D_Sigma_Star_Sigma).Rho2();
    Width_Weight_noq_Hex = /*Gamma_0 + */ pow(FF,2) * gamma_R_Hex * (D_Cascade_Star*D_Sigma_Star_Sigma).Rho2();

    /////////////////////////////////////////////////////////////////////////////////
    /////  Filling histograms  //////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////

    // Invariant mass of dsss with weights applied
    // Invariant mass of dsss with weights applied
    hInv_dsss_M_Phs_WW_FF_noq_Mol->Fill(m_dsss, Phasespace_Weight_1 * Width_Weight_noq_Mol);
    hInv_dsss_M_Phs_WW_FF_noq_Hex->Fill(m_dsss, Phasespace_Weight_1 * Width_Weight_noq_Hex);

    // Binding energies of dss with weights applied
    hBinding_dsss_M_Phs_WW_FF_noq_Mol->Fill(Binding_dsss_Mol, Phasespace_Weight_1 * Width_Weight_noq_Mol);
    hBinding_dsss_M_Phs_WW_FF_noq_Hex->Fill(Binding_dsss_Hex, Phasespace_Weight_1 * Width_Weight_noq_Hex);

    // Invariant mass of decuplet particles to check widths and branching ratios are correct
    hInv_Cascade_Star->Fill(Inv_Cascade_Star.M(),Cascade_Star_Width);
    i++;
  }

  /////////////////////////////////////////////////////////////////////////////////
  /////  Determining scaling factors and widths  //////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////


  // Determine scaling factor depending on number of events, this way histograms match no matter the statistics
  Double_t scale = 10.0 / nevents;

  // Scale histograms so they match if run with different number of events
  hInv_dsss_M_Phs_WW_FF_noq_Mol->Scale(scale);
  hBinding_dsss_M_Phs_WW_FF_noq_Mol->Scale(scale);
  hInv_dsss_M_Phs_WW_FF_noq_Hex->Scale(scale);
  hBinding_dsss_M_Phs_WW_FF_noq_Hex->Scale(scale);

  // Print out the widths in MeV (expected binding energies from table 5 and 6)
  Double_t Width_Mol_84 = 1000*hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq_Mol->FindBin(-0.084));
  Double_t Width_Mol_Expected = 1000*hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq_Mol->FindBin(-0.002));
  Double_t Width_Hex_84 = 1000*hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq_Hex->FindBin(-0.084));
  Double_t Width_Hex_Expected = 1000*hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq_Hex->FindBin(-0.234));
  cout<<"Molecule Width at 84 MeV = "<<Width_Mol_84<<" MeV, Width at expected BE ="<<Width_Mol_Expected<<" MeV"<<endl;
  cout<<"Hexaquark Width at 84 MeV = "<<Width_Hex_84<<" MeV, Width at expected BE ="<<Width_Hex_Expected<<" MeV"<<endl;

  // Write the histograms to the output file
  fileOutput1.Write();
  // Close the output file
  fileOutput1.Close();
}
