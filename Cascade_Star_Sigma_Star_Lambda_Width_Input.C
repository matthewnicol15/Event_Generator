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

   Int_t nevents=10000000; // how many events to simulate
   ostringstream Quantity;
   Quantity << "10M"; // Define quantity for output file name

   ostringstream File_Path;
   ostringstream Channel;
   ostringstream Date;
   ostringstream Output_File_Name;
   File_Path << "/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/";
   Channel << "dsss_Width_Cascade_Star_Sigma_Star_Lambda_Output";
   Date << "10032022_01";
   Output_File_Name << File_Path.str().c_str()<<Channel.str().c_str()<<"_"<<Quantity.str().c_str()<<"_"<<Date.str().c_str()<<".root";
   cout<<Output_File_Name.str().c_str()<<endl;


   // Output file path and name
   TFile fileOutput1(Output_File_Name.str().c_str(),"recreate");
   // TFile fileOutput1("../Output/Cascade_Star_Sigma_Star_Lambda_Width_Output_21022022_01.root","recreate");

   /////////////////////////////////////////////////////////////////////////////////
   /////  Create the histograms     //////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   TH2F* hbranching_ratio=new TH2F("hbranching_ratio","",200,1.3,1.5,300,0.7,1);



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


   // Histogram of the invariant mass of Sigma^{*-}
   TH2F* hInv_Sigma_Star=new TH2F("hInv_Sigma_Star","Invariant mass of p and 2#pi^{-};M(p 2#pi^{-}) [GeV];Width",2300,1.2,3.5,600,0.0,0.6);
   hInv_Sigma_Star->GetXaxis()->SetNdivisions(505);
   hInv_Sigma_Star->GetXaxis()->SetLabelSize(0.05);
   hInv_Sigma_Star->GetXaxis()->SetTitleSize(0.05);
   hInv_Sigma_Star->GetYaxis()->SetNdivisions(505);
   hInv_Sigma_Star->GetYaxis()->SetLabelSize(0.05);
   hInv_Sigma_Star->GetYaxis()->SetTitleSize(0.05);

   TH2F* hInv_Sigma_Star_Sigma=new TH2F("hInv_Sigma_Star_Sigma","Invariant mass of p and 2#pi^{-};M(p 2#pi^{-}) [GeV];Width",1000,1.3,1.4,800,0.0,0.008);
   hInv_Sigma_Star_Sigma->GetXaxis()->SetNdivisions(505);
   hInv_Sigma_Star_Sigma->GetXaxis()->SetLabelSize(0.05);
   hInv_Sigma_Star_Sigma->GetXaxis()->SetTitleSize(0.05);
   hInv_Sigma_Star_Sigma->GetYaxis()->SetNdivisions(505);
   hInv_Sigma_Star_Sigma->GetYaxis()->SetLabelSize(0.05);
   hInv_Sigma_Star_Sigma->GetYaxis()->SetTitleSize(0.05);

   TH2F* hInv_Sigma_Star_Lambda=new TH2F("hInv_Sigma_Star_Lambda","Invariant mass of p and 2#pi^{-};M(p 2#pi^{-}) [GeV];Width",1500,1.25,1.4,400,0.0,0.04);
   hInv_Sigma_Star_Lambda->GetXaxis()->SetNdivisions(505);
   hInv_Sigma_Star_Lambda->GetXaxis()->SetLabelSize(0.05);
   hInv_Sigma_Star_Lambda->GetXaxis()->SetTitleSize(0.05);
   hInv_Sigma_Star_Lambda->GetYaxis()->SetNdivisions(505);
   hInv_Sigma_Star_Lambda->GetYaxis()->SetLabelSize(0.05);
   hInv_Sigma_Star_Lambda->GetYaxis()->SetTitleSize(0.05);

   // Histogram of the invariant mass of Cascade^{*-}
   TH2F* hInv_Cascade_Star=new TH2F("hInv_Cascade_Star","Invariant mass of p, 2#pi^{-} and #pi^{0};M(p 2#pi^{-} #pi^{0}) [GeV];Width",1100,1.43,1.54,200,0.0,0.01);
   hInv_Cascade_Star->GetXaxis()->SetNdivisions(505);
   hInv_Cascade_Star->GetXaxis()->SetLabelSize(0.05);
   hInv_Cascade_Star->GetXaxis()->SetTitleSize(0.05);
   hInv_Cascade_Star->GetYaxis()->SetNdivisions(505);
   hInv_Cascade_Star->GetYaxis()->SetLabelSize(0.05);
   hInv_Cascade_Star->GetYaxis()->SetTitleSize(0.05);

   /////////////////////////////////////////////////////////////////////////////////
   /////  Produce particles here     //////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   // Produce Vectors for particles here
   // Initial particle
   TLorentzVector dsss;
   // Particles from 1st decay vertex
   TLorentzVector *Sigma_Star_Lambda, *Sigma_Star_Pip, *Cascade_Star_Cascade, *Cascade_Star_Pi0;
   // Invariant masses
   TLorentzVector Inv_Sigma_Star, Inv_Cascade_Star, Inv_dsss;
   //Binding Energies
   Double_t Binding_dsss_Mol, Binding_dsss_Hex;

   // Masses of particles
   Double_t Masses_1[4] = {1.11568, 0.13957, 1.31486, 0.13498}; // Masses of particles from 1st decay vertex
   Double_t m_dsss; // Mass of dsss, to be randomly generated in event loop
   TRandom *random = new TRandom3(); // Random number generator to produce mass of dsss

   Double_t Threshold; // Sum of decay particles
   Threshold = 1.11568 + 0.13957 + 1.31486 + 0.13498;

   /////////////////////////////////////////////////////////////////////////////////
   /////  Define constants, variables and variables for widths  ////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   Double_t Phasespace_Weight_1; // Phase weight for decay vertices
   Double_t E_Sigma_Star_Lambda_Pip, E_Sigma_Star_Sigma_Pi0, E_Cascade_Star_Pi0; // Energy of pions in rest frame of parent particles
   Double_t q_Sigma_Star_Lambda_Pip,  q_Sigma_Star_Sigma_Pi0, q_Cascade_Star_Pi0; // Momentum of pions in rest frame of parent particles
   Double_t q_Sigma_Star_Cascade_Star; // Combined momentum of Sigma^{*-} and Cascade^{*-}
   Double_t Sigma_Star_Lambda_Width, Sigma_Star_Sigma_Width,Sigma_Star_Total_Width, Cascade_Star_Width; // Breit Wigner widths
   Double_t m_Omega = 1.6725, m_Delta = 1.2320, m_Cascade_Star = 1.5318,m_Sigma_Star = 1.3828; // Nominal masses of parent particles
   Double_t gamma_BW_Cascade_Star = 0.1162, gamma_BW_Sigma_Star_Sigma = 0.1090, gamma_BW_Sigma_Star_Lambda = 0.2678; // Gamma factor for Breit wigner weights
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
      Sigma_Star_Lambda=Vertex_1.GetDecay(0);
      Sigma_Star_Pip=Vertex_1.GetDecay(1);
      Cascade_Star_Cascade=Vertex_1.GetDecay(2);
      Cascade_Star_Pi0=Vertex_1.GetDecay(3);

      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Produce invariant masses here
      Inv_Cascade_Star = (TLorentzVector)*Cascade_Star_Pi0 + (TLorentzVector)*Cascade_Star_Cascade; // Invariant mass for Cascade^{*0}
      Inv_Sigma_Star = (TLorentzVector)*Sigma_Star_Lambda + (TLorentzVector)*Sigma_Star_Pip; // Invariant mass for Sigma^{*+}
      Inv_dsss = (TLorentzVector)*Cascade_Star_Pi0 + (TLorentzVector)*Cascade_Star_Cascade + (TLorentzVector)*Sigma_Star_Lambda + (TLorentzVector)*Sigma_Star_Pip; // Invariant mass for dsss^{+}

      // Produce binding Energy
      // Produce binding Energy
      // (Hexaquark)): generated mass of dibaryon - rest mass from all decup-decup pairs (use table 7 in paper for ratios)
      // (Molecule): generated mass of dibaryon - rest mass of lightest decup pair
      Binding_dsss_Hex = m_dsss - ((m_Sigma_Star + m_Cascade_Star + m_Omega + m_Delta)/2);
      Binding_dsss_Mol = m_dsss - (m_Omega + m_Delta);


      // Getting the momentum and energy of the pions for width equations
      // Cascade Star Minus Pion
      E_Cascade_Star_Pi0 = (pow(Inv_Cascade_Star.M(),2) + pow(Masses_1[3],2) - pow(Masses_1[2],2))/(2*Inv_Cascade_Star.M());
      q_Cascade_Star_Pi0 = sqrt(abs(pow(E_Cascade_Star_Pi0,2) - pow(Masses_1[3],2)));

      // Sigma Star Minus to Lambda Pion
      E_Sigma_Star_Lambda_Pip = (pow(Inv_Sigma_Star.M(),2) + pow(Masses_1[1],2) - pow(Masses_1[0],2))/(2*Inv_Sigma_Star.M());
      q_Sigma_Star_Lambda_Pip = sqrt(abs(pow(E_Sigma_Star_Lambda_Pip,2) - pow(Masses_1[1],2)));

      // Sigma Star Minus to Sigma Pion
      E_Sigma_Star_Sigma_Pi0 = (pow(Inv_Sigma_Star.M(),2) + pow(0.13498,2) - pow(1.19745,2))/(2*Inv_Sigma_Star.M()); // Manually put in mass of sigma ground state because simulating lambda decay channel instead
      q_Sigma_Star_Sigma_Pi0 = sqrt(abs(pow(E_Sigma_Star_Sigma_Pi0,2) - pow(0.13498,2)));


      // Getting the momentum of the two parent particles
      q_Sigma_Star_Cascade_Star = Inv_Cascade_Star.Rho() + Inv_Sigma_Star.Rho();

      // Defining the Breit Wigner weights
      Cascade_Star_Width = (gamma_BW_Cascade_Star * (pow(q_Cascade_Star_Pi0,3)) * pow(R,2)) / (1 + pow(R,2) * pow(q_Cascade_Star_Pi0,2));
      Sigma_Star_Lambda_Width = ( gamma_BW_Sigma_Star_Lambda * (pow(q_Sigma_Star_Lambda_Pip,3)) * pow(R,2)) / (1 + pow(R,2) * pow(q_Sigma_Star_Lambda_Pip,2));

      // If generated mass is below threshold for sigma*->sigmapi then
      // calculated momentum will be below zero. Therefore any that are below
      // zero should have a width of 0
      if(pow(E_Sigma_Star_Sigma_Pi0,2) - pow(0.13498,2) < 0){
         Sigma_Star_Sigma_Width = 0;
      }

      else {
         Sigma_Star_Sigma_Width = ( gamma_BW_Sigma_Star_Sigma * (pow(q_Sigma_Star_Sigma_Pi0,3)) * pow(R,2)) / (1 + pow(R,2) * pow(q_Sigma_Star_Sigma_Pi0,2));
      }

      // Calculating the total width from all possible decay branches
      Sigma_Star_Total_Width = Sigma_Star_Lambda_Width + Sigma_Star_Sigma_Width;

      // Defining the complex component of the propogators
      comp_Cascade_Star = TComplex(0,m_Cascade_Star * Cascade_Star_Width);
      comp_Sigma_Star = TComplex(0,m_Sigma_Star * Sigma_Star_Total_Width);

      // Defining propogators for Cascade^{*-} and Sigma^{*-}
      D_Cascade_Star = (sqrt((m_Cascade_Star*Cascade_Star_Width)/q_Cascade_Star_Pi0))*(1.0/(pow(Inv_Cascade_Star.M(),2) - pow(m_Cascade_Star,2) + comp_Cascade_Star));
      D_Sigma_Star_Lambda = (sqrt((m_Sigma_Star*Sigma_Star_Lambda_Width)/q_Sigma_Star_Lambda_Pip))*(1.0/(pow(Inv_Sigma_Star.M(),2) - pow(m_Sigma_Star,2) + comp_Sigma_Star));

      // Defining monopole formfactor
      FF = (pow(Lambda,2)) / (pow(Lambda,2) + ((pow(q_Sigma_Star_Cascade_Star,2)) / 4));


      // Defining total weight equation for width
      Width_Weight_noq_Mol = /*Gamma_0 + */ pow(FF,2) * gamma_R_Mol * (D_Cascade_Star*D_Sigma_Star_Lambda).Rho2();
      Width_Weight_noq_Hex = /*Gamma_0 + */ pow(FF,2) * gamma_R_Hex * (D_Cascade_Star*D_Sigma_Star_Lambda).Rho2();

      /////////////////////////////////////////////////////////////////////////////////
      /////  Filling histograms  //////////////////////////////////////////////////////
      /////////////////////////////////////////////////////////////////////////////////

      // Checking widths of 1st vertex decay particles
      hInv_Sigma_Star_Sigma->Fill(Inv_Sigma_Star.M(),Sigma_Star_Sigma_Width);
      hInv_Sigma_Star_Lambda->Fill(Inv_Sigma_Star.M(),Sigma_Star_Lambda_Width);
      hInv_Sigma_Star->Fill(Inv_Sigma_Star.M(),Sigma_Star_Total_Width);
      hInv_Cascade_Star->Fill(Inv_Cascade_Star.M(),Cascade_Star_Width);
      hbranching_ratio->Fill(Inv_Sigma_Star.M(),Sigma_Star_Lambda_Width / Sigma_Star_Total_Width);

      // Invariant mass of dsss with various weights applied
      hInv_dsss_M_Phs_WW_FF_noq_Mol->Fill(m_dsss, Phasespace_Weight_1 * Width_Weight_noq_Mol); // All weights applied except phasespace and q factor
      hInv_dsss_M_Phs_WW_FF_noq_Hex->Fill(m_dsss, Phasespace_Weight_1 * Width_Weight_noq_Hex); // All weights applied except phasespace and q factor

      // Binding energies of dss with various weights applied
      hBinding_dsss_M_Phs_WW_FF_noq_Mol->Fill(Binding_dsss_Mol, Phasespace_Weight_1 * Width_Weight_noq_Mol);
      hBinding_dsss_M_Phs_WW_FF_noq_Hex->Fill(Binding_dsss_Hex, Phasespace_Weight_1 * Width_Weight_noq_Hex);

      i++;
   }

   /////////////////////////////////////////////////////////////////////////////////
   /////  Determining scaling factors and widths  //////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   // Taking projections of the invariant mass plots at the mass of the baryon
   TH1F *project = (TH1F*)hInv_Sigma_Star_Sigma->ProjectionY("project",hInv_Sigma_Star_Sigma->GetXaxis()->FindBin(1.3828),hInv_Sigma_Star_Sigma->GetXaxis()->FindBin(1.3828));
   TH1F *project2 = (TH1F*)hInv_Sigma_Star_Lambda->ProjectionY("project2",hInv_Sigma_Star_Lambda->GetXaxis()->FindBin(1.3828),hInv_Sigma_Star_Lambda->GetXaxis()->FindBin(1.3828));
   TH1F *project3 = (TH1F*)hInv_Cascade_Star->ProjectionY("project3",hInv_Cascade_Star->GetXaxis()->FindBin(1.5318),hInv_Cascade_Star->GetXaxis()->FindBin(1.5318));
   TH1F *project4 = (TH1F*)hbranching_ratio->ProjectionY("project4",hbranching_ratio->GetXaxis()->FindBin(1.3828),hbranching_ratio->GetXaxis()->FindBin(1.3828));

   // Finding the width at the mass of the baryon
   Double_t Simga_Star_Sigma_Width_Initial = project->GetBinCenter(project->GetMaximumBin());
   Double_t Simga_Star_Lambda_Width_Initial = project2->GetBinCenter(project2->GetMaximumBin());
   Double_t Cascade_Star_Width_Initial = project3->GetBinCenter(project3->GetMaximumBin());
   Double_t Branching_Ratio_Sigma_Star_Lambda = project4->GetBinCenter(project4->GetMaximumBin());

   // These are checking if the current gamma is accurate for the specific
   // baryon decay. If not then multiply your gamma values by the scaling double
   // that is printed out.
   Double_t Scaling_Sigma_Star_Sigma = 0.00461 / Simga_Star_Sigma_Width_Initial;
   Double_t Scaling_Sigma_Star_Lambda = 0.03428 / Simga_Star_Lambda_Width_Initial;
   Double_t Scaling_Cascade_Star = 0.0091 / Cascade_Star_Width_Initial;
   Double_t Branching_Ratio_Error = (fabs(0.8815 - Branching_Ratio_Sigma_Star_Lambda)) / 0.8815;

   // Printing out the scaling factors and error in branching ratio
   cout<<"Scaling_Cascade_Star "<<Scaling_Cascade_Star<<endl;
   cout<<"Scaling_Sigma_Star_Sigma "<<Scaling_Sigma_Star_Sigma<<" Scaling_Sigma_Star_Lambda "<<Scaling_Sigma_Star_Lambda<<endl;
   cout<<"Branching ratio error "<<Branching_Ratio_Error<<endl;

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
   Double_t Width_Hex_Expected = 1000*hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq_Hex->FindBin(-0.236));
   cout<<"Molecule Width at 84 MeV = "<<Width_Mol_84<<" MeV, Width at expected BE ="<<Width_Mol_Expected<<" MeV"<<endl;
   cout<<"Hexaquark Width at 84 MeV = "<<Width_Hex_84<<" MeV, Width at expected BE ="<<Width_Hex_Expected<<" MeV"<<endl;


   // Write the histograms to the output file
   fileOutput1.Write();
   // Close the output file
   fileOutput1.Close();
}
