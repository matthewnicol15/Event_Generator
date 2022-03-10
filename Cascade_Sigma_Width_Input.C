#include <TGenPhaseSpace.h>
#include <TComplex.h>
#include <TVector3.h>

void Cascade_Sigma_Width_Input(){

   Int_t nevents=1000000000; // how many events to simulate
   ostringstream Quantity;
   Quantity << "1B"; // Define quantity for output file name

   ostringstream File_Path;
   ostringstream Channel;
   ostringstream Date;
   ostringstream Output_File_Name;
   File_Path << "/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/";
   Channel << "dsss_Width_Cascade_Sigma_Output";
   Date << "24022022_01";
   Output_File_Name << File_Path.str().c_str()<<Channel.str().c_str()<<"_"<<Quantity.str().c_str()<<"_"<<Date.str().c_str()<<".root";
   cout<<Output_File_Name.str().c_str()<<endl;


   // Wavefunction overlap file
   TFile *f1=new TFile("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Wavefunction_Overlap_500_bins_14022022_1.root");
   // TFile *f1=new TFile("/media/mn688/Elements1/PhD/Macros/Dibaryon_Paper/Wavefunction_Overlap_500_bins_14022022_1.root");

   // Output file path and name
   TFile fileOutput1(Output_File_Name.str().c_str(),"recreate");
   // TFile fileOutput1("/media/mn688/Elements1/PhD/Event_Generator/Output/d2380_Width_pn_Output_1B_21022022_01.root","recreate");

   /////////////////////////////////////////////////////////////////////////////////
   ///// Wavefunction Overlap      //////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   // State number of bins used in Density distributions macro
   Int_t Bins = 500;
   // binding energy determined from distance vs binding energy plot
   Double_t Calculated_BE[500];
   // Calculating the probability for this distance/binding energy
   Double_t Probability[500];

   // Setting the reduced mass for dsss to cascade sigma
   Double_t Mass_Red_dsss_DelOmg = (1.231 * 1.672) / (1.231 + 1.672);
   Double_t Mass_Red_dsss_XiSig = (1.383 * 1.532) / (1.383 + 1.532);
   Double_t Mass_Red_dsss_Tot = (Mass_Red_dsss_DelOmg + Mass_Red_dsss_XiSig) / 2;

   // Creating histograms for distance against binding energy
   TH1D h_d_BE("h_d_BE", "Distance vs binding energy;d [fm];BE", 500, 0, 5);
   // Creating histograms for probability against binding energy
   TH1D h_probability_BE("h_probability_BE", "Probability as a function of binding energy;BE ;Probability",1000,-2.0,0.0);
   // Getting the probability against distance plot from the wavefunction file
   auto *h_probability_d  = (TH1D*) f1->Get("h_probability_d");


   Double_t x[500], y[500];
   Int_t temp = 0;
   // Loop over the bins of probability against distance plot to make TGraph
   // of this for later use
   for(Int_t i=1; i<501; i++){

      if(h_probability_d->GetBinContent(i) > 0){

         x[temp] = h_probability_d->GetBinCenter(i);

         y[temp] = h_probability_d->GetBinContent(i);

         temp++;
      }
   }

   // Make TGraph of probability against distance
   auto prob_d_graph = new TGraph(temp, x, y);

   // Loop over bins
   for(Int_t i=1; i<501; i++){
      // Get the binding energy for each bin in the distance plot
      Calculated_BE [i-1] = -pow(0.197,2) / (2.0 * Mass_Red_dsss_Tot * pow(h_d_BE.GetBinCenter(i),2));
      // Get the probability for this distance and therefore binding energy
      Probability[i-1] = prob_d_graph->Eval(h_d_BE.GetBinCenter(i));
   }
   // Plot a TGraph of the probability against binding energy
   auto prob_BE = new TGraph(500, Calculated_BE, Probability);


   /////////////////////////////////////////////////////////////////////////////////
   /////  Create the histograms     //////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   // Molecular Model
   // Histogram of the width vs the mass of dibaryon
   TH1F* hInv_dsss_M_Phs_WW_FF_noq_Mol=new TH1F("hInv_dsss_M_Phs_WW_FF_noq_Mol","Invariant mass of #Cascade^{-} and #Sigma^{-};M(#Cascade^{-} #Sigma^{-}) [GeV]; Width",3000,2.00062,5);
   hInv_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetNdivisions(505);
   hInv_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetLabelSize(0.05);
   hInv_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetTitleSize(0.05);
   hInv_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetNdivisions(505);
   hInv_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetLabelSize(0.05);
   hInv_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetTitleSize(0.05);

   // Histogram of the width vs the binding energy of dibaryon
   TH1F* hBinding_dsss_M_Phs_WW_FF_noq_Mol=new TH1F("hBinding_dsss_M_Phs_WW_FF_noq_Mol","Binding energy of dsss;Binding Energy [GeV]; Width",5000,-1.00062,4);
   hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetNdivisions(505);
   hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetNdivisions(505);
   hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetLabelSize(0.05);
   hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetXaxis()->SetTitleSize(0.05);
   hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetNdivisions(505);
   hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetLabelSize(0.05);
   hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetYaxis()->SetTitleSize(0.05);

   // Hexaquark Model
   // Histogram of the width vs the mass of dibaryon
   TH1F* hInv_dsss_M_Phs_WW_FF_noq_Hex=new TH1F("hInv_dsss_M_Phs_WW_FF_noq_Hex","Invariant mass of #Cascade^{-} and #Sigma^{-};M(#Cascade^{-} #Sigma^{-}) [GeV]; Width",3000,2.00062,5);
   hInv_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetNdivisions(505);
   hInv_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetLabelSize(0.05);
   hInv_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetTitleSize(0.05);
   hInv_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetNdivisions(505);
   hInv_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetLabelSize(0.05);
   hInv_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetTitleSize(0.05);

   // Histogram of the width vs the binding energy of dibaryon
   TH1F* hBinding_dsss_M_Phs_WW_FF_noq_Hex=new TH1F("hBinding_dsss_M_Phs_WW_FF_noq_Hex","Binding energy of dsss;Binding Energy [GeV]; Width",5000,-1.00062,4);
   hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetNdivisions(505);
   hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetNdivisions(505);
   hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetLabelSize(0.05);
   hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetXaxis()->SetTitleSize(0.05);
   hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetNdivisions(505);
   hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetLabelSize(0.05);
   hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetYaxis()->SetTitleSize(0.05);

   /////////////////////////////////////////////////////////////////////////////////
   /////  Create the particles     //////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   // Produce Vectors for particles here
   // Initial particle
   TLorentzVector dsss;
   // Particles from 1st decay vertex
   TLorentzVector *Sigma, *Cascade;
   //Binding Energies
   Double_t Binding_dsss_Mol, Binding_dsss_Hex;


   // Masses of particles
   Double_t Masses_1[2] = {1.3149, 1.1894}; // Masses of particles from 1st decay vertex
   Double_t m_dsss; // Mass of dsss, to be randomly generated in event loop
   TRandom *random = new TRandom3(); // Random number generator to produce mass of dsss

   Double_t Threshold; // Sum of decay particles
   // Adds up the decay particle rest masses to get threshold
   Threshold = 1.3149 + 1.1894;
   cout<<Threshold<<endl;

   /////////////////////////////////////////////////////////////////////////////////
   /////  Produce components for weights here    //////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   Double_t Phasespace_1; // Phase weight for decay vertices
   Double_t E_Cascade, E_Sigma; // Energy in rest frame of parent particle
   Double_t q_Cascade, q_Sigma; // Momentum of particles
   Double_t m_Cascade = 1.3149, m_Sigma = 1.1894, m_Cascade_Star_M = 1.5350,m_Sigma_Star_M = 1.3872, m_Omega_M = 1.6725, m_Delta_M = 1.2320; // Nominal masses of parent particles
   Double_t R = 6.3; // s-channel resonance in the pn and SigmaCascade systems
   // Double_t gamma_0 =  // Sum of decay widths of daughter particles
   Double_t gamma_R_Mol = 2.68015; // This contains the coupling constants and other constants, fitted to yield a total width equal to binding energy
   Double_t gamma_R_Hex = 2.1964; // This contains the coupling constants and other constants, fitted to yield a total width equal to binding energy
   Double_t Lambda = 0.16; // Cutoff parameter, this is adjusted to best reproduce the ABC effect (decup-decup only)
   Double_t FF; // Monopole formfactor, this accounts for potential barriers
   Double_t Width_Weight_noq_Mol, Width_Weight_noq_Hex; // Overall weight to apply for the width
   // with (decup-decup only) and without q factor


   /////////////////////////////////////////////////////////////////////////////////
   /////  Simulate events    //////////////////////////////////////////////
   /////////////////////////////////////////////////////////////////////////////////

   //How many events to simulate and percentage completed
   Int_t Percentage=nevents/100;


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


      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Produce binding Energy
      // (Hexaquark)): generated mass of dibaryon - rest mass from all decup-decup pairs (use table 7 in paper for ratios)
      // (Molecule): generated mass of dibaryon - rest mass of lightest decup pair
      Binding_dsss_Hex = m_dsss - ((m_Sigma_Star_M + m_Cascade_Star_M + m_Omega_M + m_Delta_M)/2);
      Binding_dsss_Mol = m_dsss - (m_Omega_M + m_Delta_M);

      // Getting the momentum and energy of the particles for width equations
      // Cascade^{-}
      E_Cascade = (pow(m_dsss,2) + pow(Masses_1[0],2) - pow(Masses_1[1],2))/(2*m_dsss);
      q_Cascade = sqrt(abs(pow(E_Cascade,2) - pow(Masses_1[0],2)));

      // Sigma^{-}
      E_Sigma = (pow(m_dsss,2) + pow(Masses_1[1],2) - pow(Masses_1[0],2))/(2*m_dsss);
      q_Sigma = sqrt(abs(pow(E_Sigma,2) - pow(Masses_1[1],2)));

      // Defining monopole formfactor
      // (octet-octet) use R in equation 3 in paper
      // (decup-decup) use Lambda using equation 5 in paper
      FF = (pow(R,4))/(1+((pow(R,4))*(pow(q_Cascade,4))));


      // Defining total weight equation for width
      Width_Weight_noq_Mol = gamma_R_Mol * FF * pow(q_Cascade,5);
      Width_Weight_noq_Hex = gamma_R_Hex * FF * pow(q_Cascade,5);

      // Determine the Wavefunction probability based on the dibaryon binding energy
      Double_t Wavefunction_prob = prob_BE->Eval(-fabs(Binding_dsss_Mol));

      // Filling histograms
      // Binding energies and mass of dibaryon with weights applied
      hInv_dsss_M_Phs_WW_FF_noq_Mol->Fill(m_dsss, Width_Weight_noq_Mol*Wavefunction_prob);
      hBinding_dsss_M_Phs_WW_FF_noq_Mol->Fill(Binding_dsss_Mol, Width_Weight_noq_Mol*Wavefunction_prob);
      hInv_dsss_M_Phs_WW_FF_noq_Hex->Fill(m_dsss, Width_Weight_noq_Hex);
      hBinding_dsss_M_Phs_WW_FF_noq_Hex->Fill(Binding_dsss_Hex, Width_Weight_noq_Hex);
      i++;
   }

   // Normalise the histograms for the number of events simulated
   Double_t scale = 10.0/nevents;
   hInv_dsss_M_Phs_WW_FF_noq_Mol->Scale(scale);
   hBinding_dsss_M_Phs_WW_FF_noq_Mol->Scale(scale);
   hInv_dsss_M_Phs_WW_FF_noq_Hex->Scale(scale);
   hBinding_dsss_M_Phs_WW_FF_noq_Hex->Scale(scale);


   Double_t Width_Mol_84 = 1000*hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq_Mol->FindBin(-0.084));
   Double_t Width_Hex_84 = 1000*hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq_Hex->FindBin(-0.084));
   Double_t Width_Mol_Expected = 1000*hBinding_dsss_M_Phs_WW_FF_noq_Mol->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq_Mol->FindBin(-0.002));
   Double_t Width_Hex_Expected = 1000*hBinding_dsss_M_Phs_WW_FF_noq_Hex->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq_Hex->FindBin(-0.234));

   cout<<"Molecule Width at 84 MeV = "<<Width_Mol_84<<" MeV, Width at expected BE ="<<Width_Mol_Expected<<" MeV"<<endl;
   cout<<"Hexaquark Width at 84 MeV = "<<Width_Hex_84<<" MeV, Width at expected BE ="<<Width_Hex_Expected<<" MeV"<<endl;


   // Save and close the output file
   fileOutput1.Write();
   fileOutput1.Close();
}
