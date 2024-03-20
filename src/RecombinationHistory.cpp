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
  Utils::StartTiming("Xe");
  
  //=============================================================================
  // TODO: Set up x-array and make arrays to store X_e(x) and n_e(x) on
  //=============================================================================
  Vector x_array = Utils::linspace(x_start, x_end, npts_rec_arrays);
  Vector Xe_arr(npts_rec_arrays);
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

    if (i % 100 == 0){
      std::cout << " Xe[i] current = " << Xe_current<< " ne[i] current = " << ne_current<<" x = " << x_array[i]<<" i = " << i <<"\n" ;
    }
    // Are we still in the Saha regime?
    if(Xe_current < Xe_saha_limit)
      saha_regime = false;

    if(saha_regime){
      
      //=============================================================================
      // TODO: Store the result we got from the Saha equation
      //=============================================================================
      Xe_arr[i] = Xe_current;  // *** 
      ne_arr[i] = ne_current; // 
      // std::cout << "Saha: Xe = " << Xe_arr[i] << " for x = " << x_array[i] << "\n";
      
    } else {

      //==============================================================
      // TODO: Compute X_e from current time til today by solving 
      // the Peebles equation (NB: if you solve all in one go remember to
      // exit the for-loop!)
      // Implement rhs_peebles_ode
      //==============================================================
      //...
      //...
      // std::cout << "CHECK PEEBLES" << "\n" ;

      // The Peebles ODE equation
      ODESolver peebles_Xe_ode;
      ODEFunction dXedx = [&](double x, const double *Xe, double *dXedx){
        return rhs_peebles_ode(x, Xe, dXedx);  
      };

      

      Vector Xe_init = {Xe_arr[i]};
      std::cout << "latest Xe after saha = " << Xe_init[0] << "\n" ;
     

      peebles_Xe_ode.solve(dXedx, x_array, Xe_init);
      auto solution_Xe = peebles_Xe_ode.get_data();

      for (size j = 0; j < solution_Xe.size(); ++j ){
        Xe_[j] = solution_Xe[j][0]
      }

      std::cout << "latest Solution_Xe" << solution_Xe.back()[0] <<  " size is" << solution_Xe.size() <<"\n" ;
      
      
      //=============================================================================
      // TODO: Set up IC, solve the ODE and fetch the result 
      //=============================================================================
      //...
      //...
    
    }
  }

  //=============================================================================
  // TODO: Spline the result. Implement and make sure the Xe_of_x, ne_of_x 
  // functions are working
  //=============================================================================
  //...
  //...

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
  const double rho_c0 = (3.0*pow(H0, 2)*OmegaB)/(8.0*Constants.pi*G);

  // Electron fraction and number density
  double Xe = 0.0;
  double ne = 0.0;


  double Tb = TCMB/a;

  // Ignoring Yp, setting Yp = 0
  //double nb = (3.0*pow(H0, 2)*OmegaB)/(8.0*Constants.pi*G*m_H*pow(a, 3)); // ### Clean up: nb is unnecessarily defined here 
  double nH = OmegaB*rho_c0 / (m_H*pow(a,2));  
  
  //=============================================================================
  // TODO: Compute Xe and ne from the Saha equation
  //=============================================================================
  
  
  // we DONT set hbar and kb equal to 1
  double kb = k_b;
  double h_bar = hbar;


  double Xe_saha_Cfrac = (kb * m_e * Tb)/(2.0*Constants.pi*pow(h_bar, 2));

  double Xe_saha_C = (1.0/nH)* pow(Xe_saha_Cfrac, (3.0/2.0)) *exp(-epsilon_0/(kb * Tb));  // Assuming nb = nH
  
  // Solving as quadratic function gives Xe = (-Xe_saha_C +- sqrt(pow(Xe_saha_C, 2)) + 4.0*Xe_saha_C))/2.0, only positive solution is:

  //std::cout << "C =  " << Xe_saha_C << "\n";

  if (4.0/Xe_saha_C < 0.01){ 
    Xe = -(Xe_saha_C / 2.0) + Xe_saha_C/2.0 * (1.0/2.0)*(1 + (Xe_saha_C/4.0)/2.0);
    Xe = 1.0; // Had to set an upper limit  *** is this allowed?
  } else {
    Xe = -(Xe_saha_C / 2.0) + Xe_saha_C/2.0 * sqrt((1+ 4.0/Xe_saha_C))/2.0;
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

  const double rho_c0 = (3.0*pow(H0, 2)*OmegaB)/(8.0*Constants.pi*G);

  

  const double alpha = 1.0/(137.0359992);

  // we DONT set hbar, c and kb equal to 1
  double kb = k_b;
  double h_bar = hbar;
  double c_ = Constants.c;

  //=============================================================================
  // TODO: Write the expression for dXedx
  //=============================================================================
  double nH = (3.0*pow(H0, 2)*OmegaB)/(8.0*Constants.pi*G*m_H*pow(a, 3)); // 1/m^3
  double n1s = (1.0-X_e)*nH;  // 1/m^3

  double Tb = TCMB/a;
  double phi_2 = 0.448 * log(epsilon_0/(kb*Tb)); // dimensionless

  double Lambda_alpha = H * pow((3.0*epsilon_0), 3)/(pow((8.0*Constants.pi),2)*pow(c_, 3)*pow(h_bar, 3)*n1s); // s^-1
  double Lambda_2s_1s = 8.227; // s^-1

  double alpha_2 = (64.0*Constants.pi)/(sqrt(27*Constants.pi)) * (pow(alpha, 2))/(pow(m_e, 2))* sqrt(epsilon_0/(Tb*kb))*phi_2; // m^3/s

  double beta = alpha_2*pow(((m_e*Tb*kb)/(2.0*Constants.pi*pow(h_bar, 2))), (3.0/2.0))*exp(-epsilon_0/(Tb*kb)); // 1/s
  double beta_2 = beta*exp( (3.0*epsilon_0)/(4.0*kb*Tb) ); // 1/s

  double Cr = (Lambda_2s_1s + Lambda_alpha)/(Lambda_2s_1s + Lambda_alpha + beta_2); // dimensionless



  double RHS = (Cr/H)*(beta * (1-X_e) - nH*alpha_2*pow(X_e, 2));
  
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
  const int npts = 1000;
  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE system dtau/dx, dtau_noreion/dx and dtau_baryon/dx
  ODEFunction dtaudx = [&](double x, const double *tau, double *dtaudx){

    //=============================================================================
    // TODO: Write the expression for dtaudx
    //=============================================================================
    //...
    //...

    // Set the derivative for photon optical depth
    dtaudx[0] = 0.0;

    return GSL_SUCCESS;
  };

  //=============================================================================
  // TODO: Set up and solve the ODE and make tau splines
  //=============================================================================
  //...
  //...

  //=============================================================================
  // TODO: Compute visibility functions and spline everything
  //=============================================================================
  //...
  //...

  Utils::EndTiming("opticaldepth");
}

//====================================================
// Get methods
//====================================================

double RecombinationHistory::tau_of_x(double x) const{
  return tau_of_x_spline(x);
}

double RecombinationHistory::dtaudx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement. Either from the tau-spline tau_of_x_spline.deriv_(x) or 
  // from a separate spline if you choose to do this
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddtauddx_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::g_tilde_of_x(double x) const{
  return g_tilde_of_x_spline(x);
}

double RecombinationHistory::dgdx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::ddgddx_tilde_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  return 0.0;
}

double RecombinationHistory::Xe_of_x(double x) const{

  //=============================================================================
  // TODO: Implement
  //=============================================================================
  //...
  //...

  double Xe_of_x = exp(log_Xe_of_x_spline(x)); // *** try this out

  return Xe_of_x;
}

double RecombinationHistory::ne_of_x(double x) const{

 
  //double ne_of_x = get_Xe_of_x(x)*(OmegaB*3*pow(H0, 2))/(Constants.m_H * pow(exp(x), 3) * 8 * Constants.pi * Constants.G)

  return 0.0;
}

double RecombinationHistory::get_Yp() const{
  return Yp;
}

//====================================================
// Print some useful info about the class
//====================================================
void RecombinationHistory::info() const{
  std::cout << "\n";
  std::cout << "Info about recombination/reionization history class:\n";
  std::cout << "Yp:          " << Yp          << "\n";
  std::cout << std::endl;
} 

//====================================================
// Output the data computed to file
//====================================================
void RecombinationHistory::output(const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts       = 5000;
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
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

