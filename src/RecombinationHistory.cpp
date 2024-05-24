#include"RecombinationHistory.h"

//====================================================
// Constructors
//====================================================
   
RecombinationHistory::RecombinationHistory(
    BackgroundCosmology *cosmo, 
    double Yp) :
  cosmo(cosmo),
  Yp(Yp)
{}

//====================================================
// Do all the solving we need to do
//====================================================

void RecombinationHistory::solve(){
    
  // Compute and spline Xe, ne
  solve_number_density_electrons();
   
  // Compute and spline tau, dtaudx, ddtauddx, g, dgdx, ddgddx, ...
  solve_for_optical_depth_tau();
}

//====================================================
// Solve for X_e and n_e using Saha and Peebles and spline the result
//====================================================

void RecombinationHistory::solve_number_density_electrons(){
  
  
  const double OmegaB      = cosmo->get_OmegaB();
  const double h           = cosmo->get_h();
  const double H0 = Constants.H0_over_h*h;
  const double rho_c0 = (3.0*pow(H0, 2.0))/(8.0*Constants.pi*Constants.G);  

  std::ofstream outFile("Xe_test.txt");

  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
  Vector Xe_arr_Saha(npts_rec_arrays); // For finding out what the Saha approximation predicts
  Vector ne_arr(npts_rec_arrays);

  // Calculate recombination history
  bool saha_regime = true;
  for(int i = 0; i < npts_rec_arrays; i++){

    //==============================================================
    // TODO: Get X_e from solving the Saha equation so
    // implement the function electron_fraction_from_saha_equation
    //==============================================================
    auto Xe_ne_data = electron_fraction_from_saha_equation(x_array[i]);
    

    // Electron fraction and number density
    const double Xe_current = Xe_ne_data.first;
    const double ne_current = Xe_ne_data.second;

    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;
      Xe_arr_Saha[i] = Xe_current; 
      ne_arr[i] = ne_current; 
      

      outFile << Xe_arr[i] << " " <<  ne_arr[i] << " "  <<  x_array[i] << " "  << "\n";

      if (i % 100 == 0){
        // std::cout << " Saha = " << Xe_arr[i] << " for x = " << x_array[i] << " for i = " << i  <<"\n" ;
        
      }
      
    } else {
       // FINDING OUT WHAT THE SAHA EQUATION PREDICTS 
      // std::cout << "THIS ONE x = " << x_array[i] << "\n";
      auto Saha_temp = electron_fraction_from_saha_equation(x_array[i]);
      if (Saha_temp.first > 0.0001){
        Xe_arr_Saha[i] = Saha_temp.first;
      }
      

      // if (i % 100 == 0){
      // std::cout << "THIS ONE Saha = " << Xe_arr_Saha[i] << "\n";
      // }
     
      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================

      

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        
        return rhs_peebles_ode(x, Xe, dXedx);

      };
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      
      Vector Xe_init = {Xe_arr[i-1]};  // Trying to solve Peebles ODE from last Saha /last Peebles
     

      // int npts_peebles = npts_rec_arrays-(i-1);  // Relic from when I did Peebles all in one go
      double peebles_xstart = x_array[i-1];
      double peebles_xend = x_array[i]; //x_end;
      int npts_peebles = 2;


      Vector x_array_Peebs = Utils::linspace(peebles_xstart, peebles_xend, npts_peebles); 

      peebles_Xe_ode.solve(dXedx, x_array_Peebs, Xe_init); 
      auto solution_Xe = peebles_Xe_ode.get_data();
      

      const double a           = exp(x_array[i]);
      double nH = OmegaB*rho_c0 / (Constants.m_H*pow(a,3.0));
      Xe_arr[i] = solution_Xe[1][0];
      ne_arr[i] = solution_Xe[1][0]*nH;

      outFile  << Xe_arr[i] << " "  << ne_arr[i] << " "  << x_array[i] << " "  << "\n" ;
      
      if (i % 100 == 0){
        //std::cout << "Peebles solution = " << Xe_arr[i] << " ne = " <<  ne_arr[i] << " for x = " << x_array[i] << " for i = " << i << "\n" ;
      }

      /*
      for (size_t j = 0; j < solution_Xe.size(); ++j){
        int k = i+j;

        const double a           = exp(x_array[k]);
        double nH = OmegaB*rho_c0 / (Constants.m_H*pow(a,3)); 
        
        
        Xe_arr[i+j] = solution_Xe[j][0];
        ne_arr[i+j] = solution_Xe[j][0]*nH;
        
        
        // std::cout << "CHECK PEEBLES " << k << "  " ;

        if (i%1 == 0){
          std::cout << "Peebles solution = " << Xe_arr[i+j] << " ne = " <<  ne_arr[i+j] << " for x = " << x_array[k] << " for i = " << k << "\n" ;
          
          
          }

        } 
        */

      // break;  
     
      
     

    }

  }

  //=============================================================================
  // Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  
  Xe_of_x_spline.create(x_array, Xe_arr, "Xe_of_x");
  ne_of_x_spline.create(x_array, ne_arr, "ne_of_x");
  Xe_of_Saha.create(x_array, Xe_arr_Saha, "Xe_of_x_Saha");

  Utils::EndTiming("Xe");
}

//====================================================
// Solve the Saha equation to get ne and Xe
//====================================================
std::pair<double,double> RecombinationHistory::electron_fraction_from_saha_equation(double x) const{
  const double a           = exp(x);
 
  // Physical constants
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Fetch cosmological parameters
  const double OmegaB      = cosmo->get_OmegaB();
  const double h           = cosmo->get_h();
  const double OmegaCDM    = cosmo->get_OmegaCDM(); 
  const double OmegaK      = cosmo->get_OmegaK();
  const double Neff        = cosmo->get_Neff(); 
  const double TCMB        = cosmo->get_TCMB();

  const double H0 = H0_over_h*h;
  const double rho_c0 = (3.0*pow(H0, 2.0))/(8.0*Constants.pi*G);

  
  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;


  double Tb = TCMB/a;

  // Ignoring Yp, setting Yp = 0
  //double nb = (3.0*pow(H0, 2)*OmegaB)/(8.0*Constants.pi*G*m_H*pow(a, 3)); // ### Clean up: nb is unnecessarily defined here 
  double nH = OmegaB*rho_c0 / (m_H*pow(a,3.0));  
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  
  
  // we DONT set hbar and kb equal to 1
  double kb = k_b;
  double h_bar = hbar;


  double Xe_saha_Cfrac = (kb * m_e * Tb)/(2.0 * Constants.pi * pow(h_bar, 2));

  double Xe_saha_C = (1.0/nH)* pow(Xe_saha_Cfrac, (3.0/2.0)) * exp(-epsilon_0/(kb * Tb));  // Assuming nb = nH
  
  // Solving as quadratic function gives Xe = (-Xe_saha_C +- sqrt(pow(Xe_saha_C, 2)) + 4.0*Xe_saha_C))/2.0, only positive solution is:

  // std::cout << "C =  " << 4.0/Xe_saha_C << "\n";

  if (4.0/Xe_saha_C < 0.001){ 
    Xe = 1.0; //(Xe_saha_C/2.0) * (2.0/Xe_saha_C );  // This is just 1
    // std::cout << "C =  " << 4.0/Xe_saha_C << "\n";
  } 
  else{
    Xe = (Xe_saha_C/2.0) * ( -1.0+ sqrt( 1.0+ 4.0/Xe_saha_C ) );
    // std::cout << "C =  " << Xe_saha_C << "\n";
  }

  ne = Xe*nH;


  return std::pair<double,double>(Xe, ne);
}

//====================================================
// The right hand side of the dXedx Peebles ODE
//====================================================
int RecombinationHistory::rhs_peebles_ode(double x, const double *Xe, double *dXedx){

  // Current value of a and X_e
  const double X_e         = Xe[0];
  const double a           = exp(x);

  // Physical constants in SI units
  const double k_b         = Constants.k_b;
  const double G           = Constants.G;
  const double c           = Constants.c;
  const double m_e         = Constants.m_e;
  const double hbar        = Constants.hbar;
  const double m_H         = Constants.m_H;
  const double sigma_T     = Constants.sigma_T;
  const double lambda_2s1s = Constants.lambda_2s1s;
  const double epsilon_0   = Constants.epsilon_0;
  const double H0_over_h   = Constants.H0_over_h;

  // Cosmological parameters
  // const double OmegaB      = cosmo->get_OmegaB();
  const double OmegaB      = cosmo->get_OmegaB();
  const double h           = cosmo->get_h();
  const double OmegaCDM    = cosmo->get_OmegaCDM(); 
  const double OmegaK      = cosmo->get_OmegaK();
  const double Neff        = cosmo->get_Neff(); 
  const double TCMB        = cosmo->get_TCMB();
        double H           = cosmo -> H_of_x(x);
  
  const double H0          = H0_over_h*h;

  const double rho_c0 = (3.0*pow(H0, 2.0))/(8.0*Constants.pi*G);

  

  const double alpha = 1.0/(137.0359992);

  // we DONT set hbar, c and kb equal to 1
  double kb = k_b;
  double h_bar = hbar;
  double c_ = Constants.c;

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  double nH = (3.0*pow(H0, 2.0)*OmegaB)/(8.0*Constants.pi*G*m_H*pow(a, 3.0)); // 1/m^3
  double n1s = (1.0-X_e)*nH;  // 1/m^3

  double Tb = TCMB/a;
  double phi_2 = 0.448 * log(epsilon_0/(kb*Tb)); // dimensionless

  double Lambda_alpha = H * pow((3.0*epsilon_0), 3.0)/(pow((8.0*Constants.pi),2.0)*pow(c_, 3.0)*pow(h_bar, 3.0)*n1s); // s^-1
  double Lambda_2s_1s = lambda_2s1s; // s^-1

  double alpha_2 = (8.0*c_)/(sqrt(3.0*Constants.pi)) * Constants.sigma_T* sqrt(epsilon_0/(Tb*kb))*phi_2; // m^3/s

  double beta = alpha_2*pow(((m_e*Tb*kb)/(2.0*Constants.pi*pow(h_bar, 2.0))), (3.0/2.0))*exp(-epsilon_0/(Tb*kb)); // 1/s
  double beta_2 = alpha_2*pow(((m_e*Tb*kb)/(2.0*Constants.pi*pow(h_bar, 2.0))), (3.0/2.0))*exp(-epsilon_0/(4.0*Tb*kb)); // 1/s
  // double beta_2 = beta*exp( (3.0*epsilon_0)/(4.0*kb*Tb) ); // 1/s

  double Cr = (Lambda_2s_1s + Lambda_alpha)/(Lambda_2s_1s + Lambda_alpha + beta_2); // dimensionless



  double RHS = (Cr/H)*(beta * (1.0-X_e) - nH*alpha_2*X_e*X_e); //pow(X_e, 2));
  // std::cout << "RHS is " << RHS << " - nH*alpha_2*pow(X_e, 2)) = " <<  - nH*alpha_2*pow(X_e, 2) << "\n";
  // std::cout << "Xe is " << X_e << " Cr is " << Cr <<  "\n";
  dXedx[0] = RHS;

  return GSL_SUCCESS;
}

//====================================================
// Solve for the optical depth tau, compute the 
// visibility function and spline the result
//====================================================

void RecombinationHistory::solve_for_optical_depth_tau(){
  Utils::StartTiming("opticaldepth");

  // Set up x-arrays to integrate over. We split into three regions as we need extra points in reionisation
  const int npts =  npts_rec_arrays;
  Vector x_array_tau_reversed = Utils::linspace(x_end, x_start, npts);  // Reversing xarray so that we can place a boundary condition on tau(x=0)
  Vector x_array = Utils::linspace(x_start, x_end, npts);  

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODESolver tau_ode;
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    // Set the derivative for photon optical depth
    dtaudx[0] = -Constants.c * ne_of_x(x) * Constants.sigma_T  / (cosmo -> H_of_x(x)) ; //@@@
    return GSL_SUCCESS;
  };

  //=============================================================================
  // Set up and solve the ODE and make tau splines
  //=============================================================================
  Vector tau_init = {0.0};  // Tau(x=0) = 0
  tau_ode.solve(dtaudx,  x_array_tau_reversed, tau_init);  // @@@
  // tau_ode.solve(dtaudx,  x_array, tau_init); // @@@

  auto solution_tau = tau_ode.get_data();

   //=============================================================================
  // Compute visibility functions and spline everything
  //=============================================================================
  Vector g_tilde(npts);
  Vector tau_(npts);
  
  
  // Saving the tau solution in the correct x direction
  for (size_t j = 0; j < solution_tau.size(); ++j){

  int k = npts - j-1;

  g_tilde[j] =  (Constants.c*ne_of_x(x_array_tau_reversed[k])*Constants.sigma_T / (cosmo->H_of_x(x_array_tau_reversed[k]))) * exp(-solution_tau[k][0]); //@@@
  tau_[j] = solution_tau[k][0];

  // PRINTING OUT AS WE GO
  // if (j%4000 == 0){ 
  //   std::cout << "Check tau loop " << tau_[j] << "  g loop " << g_tilde[j]  << " x = " << x_array_tau_reversed[k] << " \n" ;
  //   std::cout << "exp(tau) " << exp(tau_[j]) << "  deriv " << - (Constants.c*ne_of_x(x_array_tau_reversed[k])*Constants.sigma_T / (cosmo->H_of_x(x_array_tau_reversed[k])))  << "\n" <<" \n" ;
  //   }

  } 
  
  g_tilde_of_x_spline.create(x_array, g_tilde, "g");
  tau_of_x_spline.create(x_array, tau_, "tau");

  // Saving spline of dtaudx
  Vector dtaudx_(npts);
  
  for(int i=0;i<npts;i++){
    dtaudx_[i] = -Constants.c*ne_of_x(x_array[i])*Constants.sigma_T/ (cosmo->H_of_x(x_array[i])); //@@@
  }

  dtaudx_of_x_spline.create(x_array, dtaudx_, "dtaudx"); 
  
  Utils::EndTiming("opticaldepth");


  // Generate dgtildedx 
  Vector dg_tildedx(npts);
  for(int i=0;i<npts;i++){
    dg_tildedx[i] = exp(-tau_of_x(x_array[i])) * (dtaudx_of_x(x_array[i])*dtaudx_of_x(x_array[i]) - ddtauddx_of_x(x_array[i])); // ***
  }
  dg_tildedx_of_x_spline.create(x_array, dg_tildedx, "dg_tildedx"); 


  // COMPUTING THE SOUND HORIZON
  Vector s(npts);
  
  
  for (size_t i = 0; i < npts +1; i++){

    double OmegaB = cosmo -> get_OmegaB(x_array[i]);
    double OmegaR = cosmo -> get_OmegaR(x_array[i]);
    double R = 4.0*OmegaR/(3.0*OmegaB*exp(x_array[i]));

    double cs = Constants.c * sqrt( R / ( 3.0 + 3.0*R )  );

    s[i] = cs/(cosmo -> Hp_of_x(x_array[i]));

  }
  
  sound_horizon_of_x_spline.create(x_array, s, "s");

}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{
  return dtaudx_of_x_spline(x);
}

double RecombinationHistory::ddtauddx_of_x(double x) const{ 
  return dtaudx_of_x_spline.deriv_x(x);   
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{
  return dg_tildedx_of_x_spline(x);
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{
  return dg_tildedx_of_x_spline.deriv_x(x);
}

double RecombinationHistory::Xe_of_x(double x) const{
  return Xe_of_x_spline(x);
}

double RecombinationHistory::ne_of_x(double x) const{
  return ne_of_x_spline(x);
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}


// I REALIZED WHILE WRITING THAT I HAVE MISINTERPRETEDED THIS FUNCTION AND NEED TO REDO THIS ***
double RecombinationHistory::sound_horizon_of_x(double x) const{
  return sound_horizon_of_x_spline(x);
}

double RecombinationHistory::Xe_of_x_Saha(double x) const{
  return Xe_of_Saha(x);
}
//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout <<  "Tau(0) = " <<  tau_of_x(0.0)  << "\n";
  std::cout <<  "Tau(-18) = " <<  tau_of_x(-18.0)  << "\n";
  std::cout <<  "gtilde(-7) = " <<  g_tilde_of_x(-7.0)  << "\n";
  std::cout <<  "gtilde(-6) = " <<  g_tilde_of_x(-6.0)  << "\n";
  std::cout <<  "s(-6.98462) = " <<  sound_horizon_of_x(-6.98462)  << "\n";
  //-6.98462


  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = npts_rec_arrays;
  const double x_min   = x_start;
  const double x_max   = x_end;

  Vector x_array = Utils::linspace(x_min, x_max, npts);
  auto print_data = [&] (const double x) {
    fp << x                    << " ";
    fp << Xe_of_x(x)           << " ";
    fp << ne_of_x(x)           << " ";

    fp << tau_of_x(x)          << " ";
    fp << dtaudx_of_x(x)       << " ";
    fp << ddtauddx_of_x(x)     << " ";

    fp << g_tilde_of_x(x)      << " ";
    fp << dgdx_tilde_of_x(x)   << " ";
    fp << ddgddx_tilde_of_x(x) << " ";

    fp << sound_horizon_of_x(x) << " ";

    fp << Xe_of_x_Saha(x) << " ";

    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

