#include"PowerSpectrum.h"

//====================================================
// Constructors
//====================================================

PowerSpectrum::PowerSpectrum(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec, 
    Perturbations *pert,
    double A_s,
    double n_s,
    double kpivot_mpc) : 
  cosmo(cosmo), 
  rec(rec), 
  pert(pert),
  A_s(A_s),
  n_s(n_s),
  kpivot_mpc(kpivot_mpc)
{}

//====================================================
// Do all the solving
//====================================================
void PowerSpectrum::solve(){

  //=========================================================================
  // TODO: Choose the range of k's and the resolution to compute Theta_ell(k)
  //=========================================================================


  Vector k_array(n_k);

  double k_min_log = log(Constants.k_min);
  double k_max_log = log(Constants.k_max);

  Vector log_k_array = Utils::linspace(k_min_log, k_max_log, n_k);
  for(int k = 0; k < n_k; k++){ 
    k_array[k] = exp(log_k_array[k]);
  }
  // std::cout << "k_min = " << k_array[0] << " k_max = " << k_array.back() << " k size is " << k_array.size() << "\n";
  std::cout << "SOLVE BEFORE BESSEL " <<"\n";
  //=========================================================================
  // TODO: Make splines for j_ell. 
  // Implement generate_bessel_function_splines
  //=========================================================================
  generate_bessel_function_splines();
  std::cout << "SOLVE AFTER BESSEL " <<"\n";
  //=========================================================================
  // TODO: Line of sight integration to get Theta_ell(k)
  // Implement line_of_sight_integration
  //=========================================================================
  line_of_sight_integration(k_array);
  std::cout << "SOLVE AFTER LOS " <<"\n";

  //=========================================================================
  // TODO: Integration to get Cell by solving dCell^f/dlogk = Delta(k) * f_ell(k)^2
  // Implement solve_for_cell
  //=========================================================================
  auto cell_TT = solve_for_cell(log_k_array, thetaT_ell_of_k_spline, thetaT_ell_of_k_spline);
  cell_TT_spline.create(ells, cell_TT, "Cell_TT_of_ell");
  

}

//====================================================
// Generate splines of j_ell(z) needed for LOS integration
//====================================================

void PowerSpectrum::generate_bessel_function_splines(){
  Utils::StartTiming("besselspline");
  
  // Make storage for the splines
  j_ell_splines = std::vector<Spline>(ells.size());
  double eta0 = cosmo->eta_of_x(0.0);
  
  //=============================================================================
  // TODO: Compute splines for bessel functions j_ell(z)
  // Choose a suitable range for each ell
  // NB: you don't want to go larger than z ~ 40000, then the bessel routines
  // might break down. Use j_ell(z) = Utils::j_ell(ell, z)
  //=============================================================================
  double z_start = 0;
  double z_end = k_max*eta0;
  double dz = 0.05;
  int n_z = abs(z_end-z_start)/dz; 
  
  
  Vector z_array = Utils::linspace(z_start, z_end, n_z);
   
  for(size_t i = 0; i < ells.size(); i++){
    const int ell = ells[i];
    
    Vector j_ell_array(z_array.size());

    for(int j = 0; j<z_array.size(); j ++){
      // std::cout << "SOLVE <3> BESSEL " <<  z_array[j] <<"\n";
      j_ell_array[j] = Utils::j_ell(ell, z_array[j]);
    }

    // Make the j_ell_splines[i] spline
  j_ell_splines[i].create(z_array, j_ell_array);
  }
  

  Utils::EndTiming("besselspline");
}

//====================================================
// Do the line of sight integration for a single
// source function
//====================================================

Vector2D PowerSpectrum::line_of_sight_integration_single(
    Vector & k_array, 
    std::function<double(double,double)> &source_function){
  Utils::StartTiming("lineofsight");
  
  double eta0 = cosmo->eta_of_x(0.0);

  // Make storage for the results
  Vector2D result = Vector2D(ells.size(), Vector(k_array.size()));

  double dx = 0.05; // Hans ***
  std::cout << " x_start_LOS = " << x_start_LOS << " x_end_LOS = " << x_end_LOS << "\n";
  // int n_x_LOS = abs(x_start_LOS-x_end_LOS)/dx; // Trying this out
  int n_x_LOS = abs(x_end_LOS-x_start_LOS)/dx;
  std::cout << " n_x_LOS = " << n_x_LOS << "\n";
  
  Vector x_array = Utils::linspace(x_start_LOS, x_end_LOS, n_x_LOS);
  
  double k_start = k_array[0];
  double k_end = k_array.back();


 

  for(size_t ik = 0; ik < k_array.size(); ik++){  // size_t ik = 0;

    double k = k_array[ik];
    //=============================================================================
    // TODO: Implement to solve for the general line of sight integral 
    // F_ell(k) = Int dx jell(k(eta-eta0)) * S(x,k) for all the ell values for the 
    // given value of k
    //=============================================================================
    for(int il = 0; il<ells.size(); il++){
      // std::cout << "il = " << il << " ik = " << ik << " eta0-cosmo->eta_of_x(x_start_LOS) = " <<  eta0-cosmo->eta_of_x(x_start_LOS) << "\n";

      int j_idx = il;

      double result_x = ( (source_function(x_array[0],k)*j_ell_splines[j_idx]( k*(eta0-cosmo->eta_of_x(x_start_LOS) ) )) //*dx)
                        + (source_function(x_array[n_x_LOS-1],k)*j_ell_splines[j_idx](k*(eta0-cosmo->eta_of_x(x_end_LOS)))) )*dx/2.0;
      
      
      for(int ix = 1; ix<n_x_LOS-1; ix++){

        double x = x_array[ix];
        double eta = cosmo->eta_of_x(x);
        double eta_ = cosmo->eta_of_x(x_array[ix-1]);

        double f_k = source_function(x_array[ix],k)*j_ell_splines[j_idx](k*(eta0-eta));
        double f_k_ = source_function(x_array[ix-1],k)*j_ell_splines[j_idx](k*(eta0-eta_));
        
        result_x += (f_k_+f_k)*dx/2.0;
        
      }
      
      result[il][ik] = result_x;
      
    }
    // // Store the result for Source_ell(k) in results[ell][ik]
    
    
  }

  Utils::EndTiming("lineofsight");
  return result;
}

//====================================================
// Do the line of sight integration
//====================================================
void PowerSpectrum::line_of_sight_integration(Vector & k_array){
  const int n_k        = k_array.size();
  const int n          = 100;
  const int nells      = ells.size();
  
  // Make storage for the splines we are to create
  thetaT_ell_of_k_spline = std::vector<Spline>(nells);
  std::cout << "SOLVE <3> LOS "  <<"\n";
  //============================================================================
  // TODO: Solve for Theta_ell(k) and spline the result
  //============================================================================

  // Make a function returning the source function
  std::function<double(double,double)> source_function_T = [&](double x, double k){
    return pert->get_Source_T(x,k);
  };
  std::cout << "SOLVE <3> LOS 2"  <<"\n";
  // Do the line of sight integration
  Vector2D thetaT_ell_of_k = line_of_sight_integration_single(k_array, source_function_T);

  std::cout << "SOLVE <3> LOS 3" << " k_min = " << k_array[0] << " k_max = " << k_array.back() << " n_k = " << n_k <<"\n";
  // Spline the result and store it in thetaT_ell_of_k_spline
  for(size_t il=0; il<ells.size(); il++){
    thetaT_ell_of_k_spline[il].create(k_array, thetaT_ell_of_k[il]);
    // std::cout << thetaT_ell_of_k_spline[il] << "\n";
  }
  
  // //============================================================================
  // // TODO: Solve for ThetaE_ell(k) and spline
  // //============================================================================
  // if(Constants.polarization){

  //   // ...
  //   // ...
  //   // ...
  //   // ...

  // }
}

//====================================================
// Compute Cell (could be TT or TE or EE) 
// Cell = Int_0^inf 4 * pi * P(k) f_ell g_ell dk/k
//====================================================
Vector PowerSpectrum::solve_for_cell(
    Vector & log_k_array,
    std::vector<Spline> & f_ell_spline,
    std::vector<Spline> & g_ell_spline){
  const int nells      = ells.size();
 
  double pi = M_PI;
  //============================================================================
  // TODO: Integrate Cell = Int 4 * pi * P(k) f_ell g_ell dk/k
  // or equivalently solve the ODE system dCell/dlogk = 4 * pi * P(k) * f_ell * g_ell
  //============================================================================
  // C_l = (4.0/pi)* int ( A_s *  ) 
  Vector result(ells.size());
  Vector k_array = exp(log_k_array);

  double dlogk = log_k_array[1] - log_k_array[0];
  std::cout << " SOLVE FOR CELL " << "\n";
  
  
  for(int il = 0; il<ells.size(); il++){
    double result_k = ( (4.0*pi*primordial_power_spectrum(k_array[0])*f_ell_spline[il](k_array[0])*g_ell_spline[il](k_array[0])*dlogk)
                    +  (4.0*pi*primordial_power_spectrum(k_array.back())*f_ell_spline[il](k_array.back())*g_ell_spline[il](k_array.back())*dlogk)) /2.0;
   
    for(int ik = 1; ik<k_array.size(); ik++){
      double k = k_array[ik]; //+ ik*dk;
      double k_ = k_array[ik-1];

      double Cell = 4.0*pi*primordial_power_spectrum(k)*f_ell_spline[il](k)*g_ell_spline[il](k)*dlogk;
      double Cell_ = 4.0*pi*primordial_power_spectrum(k_)*f_ell_spline[il](k_)*g_ell_spline[il](k_)*dlogk;

      result_k += (Cell_ + Cell)/2.0;
    }
    result[il] = result_k;
  }
  std::cout << " SOLVE FOR CELL  DONE" << "\n";
  return result;
}

//====================================================
// The primordial power-spectrum
//====================================================

double PowerSpectrum::primordial_power_spectrum(const double k) const{
  return A_s * pow( Constants.Mpc * k / kpivot_mpc , n_s - 1.0);
}

//====================================================
// P(k) in units of (Mpc)^3
//====================================================

double PowerSpectrum::get_matter_power_spectrum(const double x, const double k_mpc) const{
  double pofk = 0.0;

  double H0 = cosmo-> get_H0();
  double OmegaCDM0 = cosmo->get_OmegaCDM(0.0);
  double OmegaB0 = cosmo->get_OmegaB(0.0);

  double OmegaM0 = OmegaCDM0 + OmegaB0;

  double c = Constants.c;
  double k_SI = k_mpc ;
  double Phi = pert->get_Phi(x, k_SI);
  double a = exp(x); // Scale factor today

  //=============================================================================
  // TODO: Compute the matter power spectrum
  //=============================================================================
  
  double Delta_M = (c*c*k_SI*k_SI * Phi * a) / ((3.0/2.0) * OmegaM0 * H0*H0);
  double k3 = k_mpc * k_mpc * k_mpc;
  double P_primordial = 2.0*M_PI*M_PI * primordial_power_spectrum(k_SI) / k3; 

  pofk = abs(Delta_M*Delta_M) * P_primordial;

  return pofk;
}

//====================================================
// Get methods
//====================================================
double PowerSpectrum::get_cell_TT(const double ell) const{
  return cell_TT_spline(ell);
}
double PowerSpectrum::get_cell_TE(const double ell) const{
  return cell_TE_spline(ell);
}
double PowerSpectrum::get_cell_EE(const double ell) const{
  return cell_EE_spline(ell);
}



//====================================================
// Output the cells to file
//====================================================

void PowerSpectrum::output_cells(std::string filename) const{
  
  // MAKE EQUIVALENT FOR k's, you might have to pick one specific ell
  // make k_array vector 
  Vector k_array(n_k);

  double k_min_log = log(Constants.k_min);
  double k_max_log = log(Constants.k_max);

  Vector log_k_array = Utils::linspace(k_min_log, k_max_log, n_k);
  for(int k = 0; k < n_k; k++){ 
    k_array[k] = exp(log_k_array[k]);
  }


  // Output in standard units of muK^2

  std::ofstream fp(filename.c_str());

  const int ellmax = int(ells[ells.size()-1]);
  auto ellvalues = Utils::linspace(2, ellmax, ellmax-1);

  auto print_data = [&] (const double ell) {

    double normfactor  = (ell * (ell+1)) / (2.0 * M_PI) * pow(1e6 * cosmo->get_TCMB(), 2);
    double normfactorN = (ell * (ell+1)) / (2.0 * M_PI) 
      * pow(1e6 * cosmo->get_TCMB() *  pow(4.0/11.0, 1.0/3.0), 2);
    double normfactorL = (ell * (ell+1)) * (ell * (ell+1)) / (2.0 * M_PI);
    fp << ell                                 << " ";
    fp << cell_TT_spline( ell ) * normfactor  << " ";

    if(Constants.polarization){
      fp << cell_EE_spline( ell ) * normfactor  << " ";
      fp << cell_TE_spline( ell ) * normfactor  << " ";
    }
    fp << "\n";
  };
  std::for_each(ellvalues.begin(), ellvalues.end(), print_data);
  
}



void PowerSpectrum::output_transfer(std::string filename) const{
  
  // MAKE EQUIVALENT FOR k's, you might have to pick one specific ell
  // make k_array vector 
  Vector k_array(n_k);

  double k_min_log = log(Constants.k_min);
  double k_max_log = log(Constants.k_max);

  Vector log_k_array = Utils::linspace(k_min_log, k_max_log, n_k);
  for(int k = 0; k < n_k; k++){ 
    k_array[k] = exp(log_k_array[k]);
  }

  

  // Output in standard units of muK^2

  std::ofstream fp(filename.c_str());

  auto print_data = [&] (const double k) {

    fp << k  << " ";
    fp << get_matter_power_spectrum(0, k) << " ";
    fp << thetaT_ell_of_k_spline[0](k)    << " ";
    fp << thetaT_ell_of_k_spline[15](k)    << " ";
    fp << thetaT_ell_of_k_spline[19](k)    << " ";
    fp << thetaT_ell_of_k_spline[32](k)    << " ";
    fp << thetaT_ell_of_k_spline[42](k)    << " ";
    fp << thetaT_ell_of_k_spline[62](k)    << " ";
    fp << thetaT_ell_of_k_spline[5](k)    << " ";

    fp << "\n";
  };
  std::for_each(k_array.begin(), k_array.end(), print_data);
  
}




