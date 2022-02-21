#include <TGenPhaseSpace.h>
#include <TComplex.h>
#include <TVector3.h>

void Cascade_Sigma_Width_Input(){

  // Wavefunction overlap file
  TFile *f1=new TFile("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Wavefunction_Overlap_500_bins_14022022_1.root");
  // TFile *f1=new TFile("/media/mn688/Elements1/PhD/Macros/Dibaryon_Paper/Wavefunction_Overlap_500_bins_14022022_1.root");

  // Output file path and name
  TFile fileOutput1("/shared/storage/physhad/JLab/mn688/Dibaryon_Paper/Output/dsss_Width_Cascade_M_Sigma_M_Hex_Output_1B_21022022_01.root","recreate");
  // TFile fileOutput1("/media/mn688/Elements1/PhD/Event_Generator/Output/dsss_Width_Cascade_M_Sigma_M_Output_1B_21022022_01.root","recreate");

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

  // Histogram of the width vs the mass of dibaryon
  TH1F* hInv_dsss_M_Phs_WW_FF_noq=new TH1F("hInv_dsss_M_Phs_WW_FF_noq","Invariant mass of #Cascade^{-} and #Sigma^{-};M(#Cascade^{-} #Sigma^{-}) [GeV]; Width",3000,2.00062,5);
  hInv_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetNdivisions(505);
  hInv_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetLabelSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetTitleSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetNdivisions(505);
  hInv_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetLabelSize(0.05);
  hInv_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetTitleSize(0.05);

  // Histogram of the width vs the binding energy of dibaryon
  TH1F* hBinding_dsss_M_Phs_WW_FF_noq=new TH1F("hBinding_dsss_M_Phs_WW_FF_noq","Binding energy of dsss;Binding Energy [GeV]; Width",5000,-1.00062,4);
  hBinding_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetLabelSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq->GetXaxis()->SetTitleSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetNdivisions(505);
  hBinding_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetLabelSize(0.05);
  hBinding_dsss_M_Phs_WW_FF_noq->GetYaxis()->SetTitleSize(0.05);


  /////////////////////////////////////////////////////////////////////////////////
  /////  Create the particles     //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  // Produce Vectors for particles here
  // Initial particle
  TLorentzVector dsss;
  // Particles from 1st decay vertex
  TLorentzVector *Sigma_M, *Cascade_M;
  //Binding Energies
  Double_t Binding_dsss;


  // Masses of particles
  Double_t Masses_1[2] = {1.3217, 1.1974}; // Masses of particles from 1st decay vertex
  Double_t m_dsss; // Mass of dsss, to be randomly generated in event loop
  TRandom *random = new TRandom3(); // Random number generator to produce mass of dsss

  Double_t Threshold; // Sum of decay particles
  // Adds up the decay particle rest masses to get threshold
  for(Int_t j=0; j<2; j++){
    Threshold += Masses_1[j];
  }
  cout<<Threshold<<endl;
  // Was giving a really high number for some reason so added in line below
  Threshold = 1.3217 + 1.1974;
  cout<<Threshold<<endl;

  /////////////////////////////////////////////////////////////////////////////////
  /////  Produce components for weights here    //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  Double_t Phasespace_1; // Phase weight for decay vertices
  Double_t E_Cascade_M, E_Sigma_M; // Energy in rest frame of parent particle
  Double_t q_Cascade_M, q_Sigma_M; // Momentum of particles
  Double_t m_Cascade_M = 1.3217, m_Sigma_M = 1.1974, m_Cascade_Star_M = 1.5350,m_Sigma_Star_M = 1.3872, m_Omega_M = 1.6725, m_Delta_M = 1.2320; // Nominal masses of parent particles
  Double_t R = 6.3; // s-channel resonance in the pn and SigmaCascade systems
  // Double_t gamma_0 =  // Sum of decay widths of daughter particles
  // Double_t gamma_R = 2.68015; // (molecular) This contains the coupling constants and other constants,
  Double_t gamma_R = 2.1964; // (compact hexaquark) This contains the coupling constants and other constants,
  Double_t Lambda = 0.16; // Cutoff parameter, this is adjusted to best reproduce the ABC effect (decup-decup only)
  Double_t FF; // Monopole formfactor, this accounts for potential barriers
  Double_t Width_Weight_noq; // Overall weight to apply for the width
  // with (decup-decup only) and without q factor


  /////////////////////////////////////////////////////////////////////////////////
  /////  Simulate events    //////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////

  //How many events to simulate and percentage completed
  Int_t nevents=1000000000;
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
    // (octet-octet): mass of dibaryon - contributions from decup-decup (use table 7 in paper for ratios)
    // (decup-decup): invariant mass of final state particles - rest mass of decup particles
    Binding_dsss = m_dsss - ((m_Sigma_Star_M + m_Cascade_Star_M + m_Omega_M + m_Delta_M)/2);


    // Getting the momentum and energy of the particles for width equations
    // Cascade^{-}
    E_Cascade_M = (pow(m_dsss,2) + pow(Masses_1[0],2) - pow(Masses_1[1],2))/(2*m_dsss);
    q_Cascade_M = sqrt(abs(pow(E_Cascade_M,2) - pow(Masses_1[0],2)));

    // Sigma^{-}
    E_Sigma_M = (pow(m_dsss,2) + pow(Masses_1[1],2) - pow(Masses_1[0],2))/(2*m_dsss);
    q_Sigma_M = sqrt(abs(pow(E_Sigma_M,2) - pow(Masses_1[1],2)));

    // Defining monopole formfactor
    // (octet-octet) use R in equation 3 in paper
    // (decup-decup) use Lambda using equation 5 in paper
    FF = (pow(R,4))/(1+((pow(R,4))*(pow(q_Cascade_M,4))));


    // Defining total weight equation for width
    Width_Weight_noq = FF * pow(q_Cascade_M,5);

    // Determine the Wavefunction probability based on the dibaryon binding energy
    Double_t Wavefunction_prob = prob_BE->Eval(-fabs(Binding_dsss));

    // Filling histograms
    // Binding energies and mass of dibaryon with weights applied
    hInv_dsss_M_Phs_WW_FF_noq->Fill(m_dsss, Width_Weight_noq/**Wavefunction_prob*/);
    hBinding_dsss_M_Phs_WW_FF_noq->Fill(Binding_dsss, Width_Weight_noq/**Wavefunction_prob*/);
    i++;
  }

  // Normalise the histograms for the number of events simulated
  Double_t scale = 10.0/nevents;
  hInv_dsss_M_Phs_WW_FF_noq->Scale(scale);
  hBinding_dsss_M_Phs_WW_FF_noq->Scale(scale);


  Double_t Width = 1000*hBinding_dsss_M_Phs_WW_FF_noq->GetBinContent(hBinding_dsss_M_Phs_WW_FF_noq->FindBin(-0.084));
  cout<<"Width = "<<Width<<" MeV"<<endl;

  // Save and close the output file
  fileOutput1.Write();
  fileOutput1.Close();
}
