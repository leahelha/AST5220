#include "Utils.h"
#include "BackgroundCosmology.h"
#include "RecombinationHistory.h"
#include "Perturbations.h"
#include "PowerSpectrum.h"
#include "SupernovaFitting.h"

int main(int argc, char **argv){
  Utils::StartTiming("Everything");

  //=========================================================================
  // Parameters
  //=========================================================================

  // // Background parameters
  double h           = 0.67;
  double OmegaB      = 0.05;
  double OmegaCDM    = 0.267;
  double OmegaK      = 0.0;
  double Neff        = 0.0; //3.046; //***
  double TCMB        = 2.7255;

  // Best fit background params
  //  h        OmegaM      OmegaK    
  // 0.701711   0.255027   0.0789514
  // double h           = 0.701711;
  // double OmegaB      = 0.05;
  // double OmegaCDM    = 0.255027-0.05;
  // double OmegaK      = 0.0789514;
  // double Neff        = 0.0; 
  // double TCMB        = 2.7255;


  // Recombination parameters
  double Yp          = 0; //= 0.245;

  // Power-spectrum parameters
  double A_s         = 2.1e-9;
  double n_s         = 0.965;
  double kpivot_mpc  = 0.05;

  //=========================================================================
  // Module I
  //=========================================================================

  // Set up and solve the background
  BackgroundCosmology cosmo(h, OmegaB, OmegaCDM, OmegaK, Neff, TCMB);
  cosmo.solve();
  cosmo.info();
  
  
  // Output background evolution quantities
  // cosmo.output("best_params_cosmology.txt");

  // Output background evolution quantities
  cosmo.output("cosmology.txt");

  // Do the supernova fits. Uncomment when you are ready to run this
  // Make sure you read the comments on the top of src/SupernovaFitting.h
  // mcmc_fit_to_supernova_data("data/supernovadata.txt", "results_supernovafitting.txt");

  
// return 0.0;
  //=========================================================================
  // Module II
  //=========================================================================

  // Solve the recombination history
  RecombinationHistory rec(&cosmo, Yp);
  rec.solve();
  rec.info();

  // Output recombination quantities
  rec.output("recombination.txt");
  

  //=========================================================================
  // Module III
  //=========================================================================
 
  // Solve the perturbations
  Perturbations pert(&cosmo, &rec);
  pert.solve();
  pert.info();
  
  // Output perturbation quantities
  double kvalue = 0.01 / Constants.Mpc;
  pert.output(kvalue, "perturbations_k0.01.txt");

  double kvalue2 = 0.1 / Constants.Mpc;
  pert.output(kvalue2, "perturbations_k0.1.txt");

  double kvalue3 = 0.001 / Constants.Mpc;
  pert.output(kvalue3, "perturbations_k0.001.txt");
  
  
  
  // Power spectrum parameters
  // double A_s = 2.1e-9;
  // double n_s = 0.965;
  // double kpivot_mpc = 0.05;
  //=========================================================================
  // Module IV
  //=========================================================================

  PowerSpectrum power(&cosmo, &rec, &pert, A_s, n_s, kpivot_mpc);
  power.solve();
  power.output("cells.txt");
  
  // Remove when module is completed
  return 0.0;
  Utils::EndTiming("Everything");
}
