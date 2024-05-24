#include"Perturbations.h"

//===============================================================================
// Constructors
//===============================================================================

Perturbations::Perturbations(
    BackgroundCosmology *cosmo, 
    RecombinationHistory *rec) : 
  cosmo(cosmo), 
  rec(rec)
{}

//===============================================================================
// Do all the solving
//===============================================================================

void Perturbations::solve(){
  
  // Integrate all the perturbation equation and spline the result
  integrate_perturbations();


  
  // Compute source functions and spline the result
  std::cout << "\n" << "Computing source functions: TRUE" << "\n";
  compute_source_functions();

  
}

//===============================================================================
// The main work: integrate all the perturbations
// and spline the results
//===============================================================================

void Perturbations::integrate_perturbations(){
  Utils::StartTiming("integrateperturbation");

  //==============================================================================================
  // TODO: Set up the k-array for the k's we are going to integrate over
  // Start at k_min end at k_max with n_k points with either a
  // quadratic or a logarithmic spacing
  //==============================================================================================
  Vector k_array(n_k);
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  
  double k_min_log = log(k_min);
  double k_max_log = log(k_max);

  Vector log_k_array = Utils::linspace(k_min_log, k_max_log, n_k);
  for(int k = 0; k < n_k+1; k++){ // *** Added +1 
    k_array[k] = exp(log_k_array[k]);
  }


  Vector f_delta_cdm(n_k*n_x);
  Vector f_delta_b(n_k*n_x);
  Vector f_v_cdm(n_k*n_x);
  Vector f_v_b(n_k*n_x);
  Vector f_Phi(n_k*n_x);
  Vector f_Psi(n_k*n_x);
  std::vector<Vector> f_Theta(3, Vector(n_k*n_x));

  

  // Loop over all wavenumbers
  // #pragma omp parallel for schedule(dynamic, 1)
  for(int ik = 0; ik < n_k; ik++){
    // Progress bar...
    if( (10*ik) / n_k != (10*ik+10) / n_k ) {
      std::cout << (100*ik+100)/n_k << "% " << std::flush;
      if(ik == n_k-1) std::cout << std::endl;
    }
    
    // Current value of k
    double k = k_array[ik];

    // Find value to integrate to
    double x_end_tight = get_tight_coupling_time(k); 

    // Find index of x_end_tc in x_array
    int idx_end_tc;
    for (int i = 0; i < n_x; i++){ // Hans: it was +1 ? Why?
      if (x_array[i] == x_end_tight){
        idx_end_tc = i;
      }
    }
    

    //==============================================================================================
    // TODO: Tight coupling integration
    // Remember to implement the routines:
    // set_ic : The IC at the start
    // rhs_tight_coupling_ode : The dydx for our coupled ODE system
    //==============================================================================================

    // Set up initial conditions in the tight coupling regime
    auto y_tight_coupling_ini = set_ic(x_start, k);
    // double *Theta = &y_tight_coupling_ini[Constants.ind_start_theta_tc];

    // for(auto & v : y_tight_coupling_ini){
    //   std::cout << v << "\n";
    // }

    // The tight coupling ODE system
    ODEFunction dydx_tight_coupling = [&](double x, const double *y, double *dydx){
      return rhs_tight_coupling_ode(x, k, y, dydx);
    };

    // std::cout << "count tc regime" << count_tc << " idx_end_tc " << idx_end_tc << "\n";
   
    int n_x_tc = idx_end_tc;
    // Integrate from x_start -> x_end_tight
    Vector x_tc = Utils::linspace(x_start, x_array[idx_end_tc], n_x_tc); 


    Vector vector_y_tc_ini {y_tight_coupling_ini};
    ODESolver ode_y_tight_coupling;

    ode_y_tight_coupling.solve(dydx_tight_coupling, x_tc, vector_y_tc_ini);  
    auto solution_y_tight_coupling = ode_y_tight_coupling.get_data();

    // std::cout << "TIGHT COUPLING PHI INITIAL" << solution_y_tight_coupling[Constants.ind_Phi_tc][0]<< "\n";
    
    //===============================================================================================
    // TODO: Full equation integration
    // Remember to implement the routines:
    // set_ic_after_tight_coupling : The IC after tight coupling ends
    // rhs_full_ode : The dydx for our coupled ODE system
    //===============================================================================================
    int n_x_full = n_x - idx_end_tc + 1; //*** +1 for 1 index overlap 
    // Integrate from x_end_tight -> x_end
    Vector x_full = Utils::linspace(x_array[idx_end_tc], x_end, n_x_full);
    


    // Set up initial conditions (y_tight_coupling is the solution at the end of tight coupling)
    Vector y_tight_coupling = solution_y_tight_coupling.back(); 
    auto y_full_ini = set_ic_after_tight_coupling(y_tight_coupling, x_end_tight, k); 


    // The full ODE system
    ODEFunction dydx_full = [&](double x, const double *y, double *dydx){
      return rhs_full_ode(x, k, y, dydx);
      // std::cout << "THETA" << "\n";
      // std::cout << Theta[0] << " <- THETA full" << "\n";
    };
    
    ODESolver ode_y_full;
    ode_y_full.solve(dydx_full, x_full, y_full_ini);
    auto solution_y_full = ode_y_full.get_data();
    
    // std::cout << Theta[0] << " <- THETA full" << "\n";
    
    //==============================================================================================
    // TODO: remember to store the data found from integrating so we can
    // spline it below
    //==============================================================================================

    // Tight coupling regime
    for (int ix = 0; ix <idx_end_tc; ix++){  
      int idx = ix + n_x*ik;
      auto y = solution_y_tight_coupling[ix];

      // References to the quantities we are going to set
      double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
      double &delta_b         =  y[Constants.ind_deltab_tc];
      double &v_cdm           =  y[Constants.ind_vcdm_tc];
      double &v_b             =  y[Constants.ind_vb_tc];
      double &Phi             =  y[Constants.ind_Phi_tc];
      double *Theta           = &y[Constants.ind_start_theta_tc];
      double *Theta_p         = &y[Constants.ind_start_thetap_tc];
      double *Nu              = &y[Constants.ind_start_nu_tc];

      double x = x_array[ix];
      // std::cout << "x is " << x <<"\n";
      double k = k_array[ik];

      // Cosmological parameters and variables
      double a = exp(x);
      double c = Constants.c;
      double Hp = cosmo->Hp_of_x(x);
      double ck_Hp = c*k/Hp;

      double Omega_R = cosmo->get_OmegaR(0.0); // For Psi
      double h          = cosmo -> get_h();
      double H0          = Constants.H0_over_h*h;

      // Recombination variables
      double dtau = rec->dtaudx_of_x(x);

      // Other
      double Theta2 = - 20.0*ck_Hp/(45.0*dtau) *Theta[1];
      //***
      // SET: Photon temperature perturbations (Theta_ell)
      f_Theta[0][idx] = Theta[0];
      f_Theta[1][idx] = Theta[1];
      f_Theta[2][idx] = Theta2;//(-20.0*c*k)/(45.0*Hp*dtau)*Theta[1];
      // for( int l=3; l<Constants.n_ell_theta; l++){
      //   f_Theta[l][idx] = -l/(2.0*l+1.0) * ck_Hp/dtau *Theta[l-1];
      // }

      // SET: Scalar quantities (Gravitational potental, baryons and CDM)
    
      f_Phi[idx] = Phi;
      f_Psi[idx] = -Phi -((12.0*H0*H0)/(c*c*k*k*a*a)) * (Omega_R*Theta2); // *** THETA[2] became Theta2
      f_delta_cdm[idx] = delta_cdm;
      f_delta_b[idx] = delta_b;
      f_v_cdm[idx] = v_cdm;
      f_v_b[idx] = v_b;

      

      // if(k == 0.01){
      //   std::cout << "\n" << "TC in arrays Phi value"<< f_Phi[idx] << " IDK = " << idx << "\n";

      // }
      

    }



    // Full system 
    for(int ix = idx_end_tc; ix < n_x; ix ++){      
      int idx = ix + n_x*ik;
      
      //std::cout << "fetching yfull id = " << ix << "\n"; 
      auto y = solution_y_full[ix-idx_end_tc];

      // int real_idx = ix; //- idx_end_tc + 1; // +1 overlap

      // References to the quantities we are going to set
      double &delta_cdm       =  y[Constants.ind_deltacdm];
      double &delta_b         =  y[Constants.ind_deltab];
      double &v_cdm           =  y[Constants.ind_vcdm];
      double &v_b             =  y[Constants.ind_vb];
      double &Phi             =  y[Constants.ind_Phi];
      double *Theta           = &y[Constants.ind_start_theta];
      double *Theta_p         = &y[Constants.ind_start_thetap];
      double *Nu              = &y[Constants.ind_start_nu];

      double x = x_array[ix];
      double k = k_array[ik];

      // Cosmological parameters and variables
      //std::cout << "IDX " << idx  << " idx_end_tc " << idx_end_tc << " x IS " <<  x << "\n";
      // std::cout << "The length of the vector is: " << f_Theta[0].size() << "\n";
      
      double a = exp(x);
      double c = Constants.c;
      double Hp = cosmo->Hp_of_x(x);
      double ck_Hp = c*k/Hp;
      double Omega_R = cosmo->get_OmegaR(0.0); // For Psi
      double h          = cosmo -> get_h();
      double H0          = Constants.H0_over_h*h;
      
      
      // Recombination variables
      double dtau = rec->dtaudx_of_x(x); 
      
      
      // SET: Photon temperature perturbations (Theta_ell)
      
      f_Theta[0][idx] = Theta[0];
      f_Theta[1][idx] = Theta[1];
      f_Theta[2][idx] = Theta[2];//- 20.0*ck_Hp/(45.0*dtau) *Theta[1];// (-20.0*c*k)/(45.0*Hp*dtau)*Theta[1]; 

      // std::cout << "THETA[2] = " << Theta[2] << " THeta 2 = " << (- 20.0*ck_Hp/(45.0*dtau) *Theta[1]) << "\n";
      // for(int l=3; l<Constants.n_ell_theta; l++){
      //   f_Theta[l][idx] = -l/(2.0*l+1.0) * ck_Hp/dtau *Theta[l-1];
      // }

      // // SET: Scalar quantities (Gravitational potental, baryons and CDM)
      f_Phi[idx] = Phi;
      f_Psi[idx] = -Phi -((12.0*H0*H0)/(c*c*k*k*a*a)) * (Omega_R*Theta[2]); 
      


      f_delta_cdm[idx] = delta_cdm;
      f_delta_b[idx] = delta_b;
      f_v_cdm[idx] = v_cdm;
      f_v_b[idx] = v_b;
    }
  

    //===================================================================
    // To compute a 2D spline of a function f(x,k) the data must be given 
    // to the spline routine as a 1D array f_array with the points f(ix, ik) 
    // stored as f_array[ix + n_x * ik]
    // Example:
    // Vector x_array(n_x);
    // Vector k_array(n_k);
    // Vector f(n_x * n_k);
    // Spline2D y_spline;
    // f_spline.create(x_array, k_array, f_array);
    // We can now use the spline as f_spline(x, k)
    //
    // NB: If you use Theta_spline then you have to allocate it first,
    // before using it e.g.
    // Theta_spline = std::vector<Spline2D>(n_ell_theta);
    //
    //===================================================================
    //...
    //...

  }
  Utils::EndTiming("integrateperturbation");

  //=============================================================================
  // TODO: Make all splines needed: Theta0,Theta1,Theta2,Phi,Psi,...
  //=============================================================================
    // Spline2D delta_cdm_spline{"delta_cdm_spline"};
    // Spline2D delta_b_spline{"delta_b_spline"};
    // Spline2D v_cdm_spline{"v_cdm_spline"};
    // Spline2D v_b_spline{"v_b_spline"};
    // Spline2D Phi_spline{"Phi_spline"};
    // Spline2D Pi_spline{"Pi_spline"};
    // Spline2D Psi_spline{"Psi_spline"};

    // // e.g. Theta_spline = std::vector<Spline2D>(n_ell_theta); before using it
    // std::vector<Spline2D> Theta_spline;
    // std::vector<Spline2D> Theta_p_spline;
    // std::vector<Spline2D> Nu_spline;
   //=============================================================================
    // delta_cdm_spline.create

    delta_cdm_spline.create(x_array, k_array, f_delta_cdm, "delta_cdm_spline");
    delta_b_spline.create(x_array, k_array, f_delta_b, "delta_b_spline");

    v_cdm_spline.create(x_array, k_array, f_v_cdm, "v_cdm_spline");
    v_b_spline.create(x_array, k_array, f_v_b, "v_b_spline");

    Phi_spline.create(x_array, k_array, f_Phi, "Phi_spline");
    Psi_spline.create(x_array, k_array, f_Psi, "Psi_spline");

    

    // Spline2D Theta_spline;
    // // Theta_spline = std::vector<Spline2D>(Constants.n_ell_theta);
    // Theta_spline.create(k_array, x_array, f_Theta, "Theta");

    Theta_spline = std::vector<Spline2D>(3);
     for(int ell = 0; ell < 3; ell++){
        Theta_spline[ell].create(x_array, k_array, f_Theta[ell], "Theta_"+std::to_string(ell) );
    }
    
}

//====================================================
// Set IC at the start of the run (this is in the
// tight coupling regime)
//====================================================
Vector Perturbations::set_ic(const double x, const double k) const{

  // The vector we are going to fill
  Vector y_tc(Constants.n_ell_tot_tc);

  //=============================================================================
  // Compute where in the y_tc array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const int n_ell_tot_tc        = Constants.n_ell_tot_tc;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // References to the tight coupling quantities
  double &delta_cdm    =  y_tc[Constants.ind_deltacdm_tc];
  double &delta_b      =  y_tc[Constants.ind_deltab_tc];
  double &v_cdm        =  y_tc[Constants.ind_vcdm_tc];
  double &v_b          =  y_tc[Constants.ind_vb_tc];
  double &Phi          =  y_tc[Constants.ind_Phi_tc];
  double *Theta        = &y_tc[Constants.ind_start_theta_tc];
  double *Nu           = &y_tc[Constants.ind_start_nu_tc];

  double Hp            = cosmo -> Hp_of_x(x);
  double dtau          = rec -> dtaudx_of_x(x);
  //=============================================================================
  // TODO: Set the initial conditions in the tight coupling regime
  //=============================================================================
  double Psi = - 2.0/3.0;//- 1.0/ ( 3.0/2.0 + (2.0*f_nu)/5.0  );

  // SET: Scalar quantities (Gravitational potential, baryons and CDM)
  Phi = 2.0/3.0;//-Psi; //-(1.0 + (2.0*f_nu)/5.0)*Psi;

  delta_cdm = -3.0/2.0 * Psi;
  delta_b = -3.0/2.0 * Psi;

  v_cdm = -(Constants.c * k)/(2.0*Hp) * Psi;
  v_b = -(Constants.c * k)/(2.0*Hp) * Psi;


  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = -1.0/2.0 * Psi;
  Theta[1] = (Constants.c * k)/(6.0*Hp)*Psi;

  // We dont consider Theta polarization
  // Theta[2] =  -(20.0 * Constants.c)/(45.0 * Hp * dtau ) * Theta[1];

  return y_tc;
}

//====================================================
// Set IC for the full ODE system after tight coupling 
// regime ends
//====================================================

Vector Perturbations::set_ic_after_tight_coupling(
    const Vector &y_tc, 
    const double x, 
    const double k) const{

  // Make the vector we are going to fill
  Vector y(Constants.n_ell_tot_full);
  
  //=============================================================================
  // Compute where in the y array each component belongs and where corresponding
  // components are located in the y_tc array
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Number of multipoles we have in the full regime
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // Number of multipoles we have in the tight coupling regime
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;

  // References to the tight coupling quantities
  const double &delta_cdm_tc    =  y_tc[Constants.ind_deltacdm_tc];
  const double &delta_b_tc      =  y_tc[Constants.ind_deltab_tc];
  const double &v_cdm_tc        =  y_tc[Constants.ind_vcdm_tc];
  const double &v_b_tc          =  y_tc[Constants.ind_vb_tc];
  const double &Phi_tc          =  y_tc[Constants.ind_Phi_tc];
  const double *Theta_tc        = &y_tc[Constants.ind_start_theta_tc];
  const double *Nu_tc           = &y_tc[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set
  double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  double &delta_b         =  y[Constants.ind_deltab_tc];
  double &v_cdm           =  y[Constants.ind_vcdm_tc];
  double &v_b             =  y[Constants.ind_vb_tc];
  double &Phi             =  y[Constants.ind_Phi_tc];
  double *Theta           = &y[Constants.ind_start_theta_tc];
  double *Theta_p         = &y[Constants.ind_start_thetap_tc];
  double *Nu              = &y[Constants.ind_start_nu_tc];


  double c = Constants.c;
  double Hp = cosmo -> Hp_of_x(x);
  double ck_Hp = c*k/Hp;
  double dtau = rec -> dtaudx_of_x(x);

  //=============================================================================
  // TODO: fill in the initial conditions for the full equation system below
  // NB: remember that we have different number of multipoles in the two
  // regimes so be careful when assigning from the tc array
  //=============================================================================

  // SET: Scalar quantities (Gravitational potental, baryons and CDM)
  Phi = Phi_tc;
  delta_cdm = delta_cdm_tc;
  delta_b = delta_b_tc;
  v_cdm = v_cdm_tc;
  v_b = v_b_tc;

  // SET: Photon temperature perturbations (Theta_ell)
  Theta[0] = Theta_tc[0];
  Theta[1] = Theta_tc[1];
  Theta[2] = (-20.0*c*k)/(45.0*Hp*dtau)*Theta_tc[1];
  // Theta[2] = Theta_tc[2];
  // std::cout << "INITIAL CONDITION THETA[2] = " << Theta_tc[2] << "vs " << ((-20.0*c*k)/(45.0*Hp*dtau)*Theta[1]) << "\n";

  for( int l=3; l<n_ell_theta; l++){
    Theta[l] = -l/(2.0*l+1.0) * ck_Hp/dtau *Theta[l-1];

  }
  


  return y;
}

//====================================================
// The time when tight coupling end
//====================================================

double Perturbations::get_tight_coupling_time(const double k) const{
  double x_tight_coupling_end = 0.0;
  double c = Constants.c;
 
  

  

  //=============================================================================
  // TODO: compute and return x for when tight coupling ends
  // Remember all the three conditions in Callin

  // Condition 1: |k/(Hτ′)| < 1/10
  // Condition 2: |τ′| > 10
  // Condition 3: time x must be before Recombination
  //=============================================================================

  // x needs to be BEFORE the start of recombination. From analysis: x_recomb = −6.9882
  double x_rec = -8.3;
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  
  
  for (double x : x_array){
    if (x<x_rec){ 
      double Hp = cosmo -> Hp_of_x(x);
      double dtau = rec -> dtaudx_of_x(x);
      double ck_Hp = c*k/Hp;
      if(fabs(dtau)> 10.0){
        // if (fabs(k/(Hp * dtau)) < (1.0/10.0)){
        //   x_tight_coupling_end = x;        
        //   // std::cout << "x at the end of TC" << x_tight_coupling_end << "\n";
        // } 
        if (fabs(dtau)>(10.0*ck_Hp)){
          x_tight_coupling_end = x;        

        }
      }
    }
    else{
      break;
    }
  }

  // std::cout << "x_tight_coupling_end = " << x_tight_coupling_end << "\n";
  return x_tight_coupling_end;
}

//====================================================
// After integrsating the perturbation compute the
// source function(s)
//====================================================
void Perturbations::compute_source_functions(){
  Utils::StartTiming("source");

  //=============================================================================
  // Making the x and k arrays to evaluate over and use to make the splines
  //=============================================================================
  // Vector k_array(n_k);
  // Vector x_array(n_x);

  Vector k_array(n_k);
  Vector x_array = Utils::linspace(x_start, x_end, n_x);
  
  double k_min_log = log(k_min);
  double k_max_log = log(k_max);

  Vector log_k_array = Utils::linspace(k_min_log, k_max_log, n_k);
  for(int k = 0; k < n_k+1; k++){  // *** ADDED +1
    k_array[k] = exp(log_k_array[k]);
  }



  // Make storage for the source functions (in 1D array to be able to pass it to the spline)
  Vector ST_array(k_array.size() * x_array.size());
  // Vector SE_array(k_array.size() * x_array.size());

  Vector SW_array(k_array.size() * x_array.size());
  Vector ISW_array(k_array.size() * x_array.size());
  Vector Doppler_array(k_array.size() * x_array.size());
  Vector Polarization_array(k_array.size() * x_array.size());

  // Compute source functions
  for(auto ix = 0; ix < x_array.size(); ix++){
    const double x = x_array[ix];
    for(auto ik = 0; ik < k_array.size(); ik++){
      const double k = k_array[ik];

      // NB: This is the format the data needs to be stored 
      // in a 1D array for the 2D spline routine source(ix,ik) -> S_array[ix + nx * ik]
      const int index = ix + n_x * ik;

      //=============================================================================
      // TODO: Compute the source functions
      //=============================================================================
      // Fetch all the things we need...
      const double Hp       = cosmo->Hp_of_x(x);
      const double dHp      = cosmo->dHpdx_of_x(x);

      const double tau      = rec->tau_of_x(x);
      const double dtau     = rec->dtaudx_of_x(x);

      const double g_tilde  = rec->g_tilde_of_x(x);
      const double dg_tilde = rec-> dgdx_tilde_of_x(x); 
      const double ddg_tilde = rec->ddgddx_tilde_of_x(x);
      
      const double Psi = get_Psi(x,k);
      const double dPsi = Psi_spline.deriv_x(x,k);
      const double Phi = get_Phi(x,k);
      const double dPhi = Phi_spline.deriv_x(x, k);

      const double v_b = get_v_b(x,k);
      const double dv_b = v_b_spline.deriv_x(x, k);

      const double Theta_0 = get_Theta(x, k, 0);
      const double Theta_1 = get_Theta(x, k, 1);
      const double Theta_2 = get_Theta(x, k, 2);

      double Pi = Theta_2;
      double dPi = Theta_spline[2].deriv_x(x,k);
      double ddPi = Theta_spline[2].deriv_xx(x,k);

      const double c = Constants.c;
       
      SW_array[index] = g_tilde*(Theta_0 + Psi + 0.25*Pi) ;

      ISW_array[index] = exp(-tau)*(dPsi - dPhi);

      double d_Hgv_b_dx = dHp*g_tilde*v_b + Hp*dg_tilde*v_b + Hp*g_tilde*dv_b;
      Doppler_array[index] = - 1.0/(c*k) * (d_Hgv_b_dx);

      double dHpHpder_g_pi = dHp*dHp*g_tilde*Pi + Hp*dHp*g_tilde*Pi + Hp*dHp*dg_tilde*Pi + Hp*dHp*g_tilde*dPi;
      double deriv_rest = 3.0*Hp*dHp*(dg_tilde*Pi + g_tilde*dPi) + Hp*Hp*(ddg_tilde*Pi+2.0*dg_tilde*dPi+g_tilde*ddPi);
      Polarization_array[index] = 3.0/(4.0*c*c*k*k);

    

      // Temperatur source
      ST_array[index] = SW_array[index] + ISW_array[index] + Doppler_array[index] + Polarization_array[index];

      // Polarization source
      // if(Constants.polarization){
      //   SE_array[index] = 0.0;
      // }
    
    }
  }
  
  // Spline the source functions
  ST_spline.create(x_array, k_array, ST_array, "Source_Temp_x_k");
  // if(Constants.polarization){
  //   SE_spline.create (x_array, k_array, SE_array, "Source_Pol_x_k");
  // }

  Utils::EndTiming("source");
}

//====================================================
// The right hand side of the perturbations ODE
// in the tight coupling regime
//====================================================

// Derivatives in the tight coupling regime
int Perturbations::rhs_tight_coupling_ode(double x, double k, const double *y, double *dydx){

  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================
  
  // For integration of perturbations in tight coupling regime (Only 2 photon multipoles + neutrinos needed)
  const int n_ell_theta_tc      = Constants.n_ell_theta_tc;
  const int n_ell_neutrinos_tc  = Constants.n_ell_neutrinos_tc;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm_tc];
  const double &delta_b         =  y[Constants.ind_deltab_tc];
  const double &v_cdm           =  y[Constants.ind_vcdm_tc];
  const double &v_b             =  y[Constants.ind_vb_tc];
  const double &Phi             =  y[Constants.ind_Phi_tc];
  const double *Theta           = &y[Constants.ind_start_theta_tc];
  const double *Nu              = &y[Constants.ind_start_nu_tc];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm_tc];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab_tc];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm_tc];
  double &dv_bdx          =  dydx[Constants.ind_vb_tc];
  double &dPhidx          =  dydx[Constants.ind_Phi_tc];
  double *dThetadx        = &dydx[Constants.ind_start_theta_tc];
  double *dNudx           = &dydx[Constants.ind_start_nu_tc];
  // Cosmological parameters and variables
  double a = exp(x);
  double Hp = cosmo->Hp_of_x(x);
  double dHp = cosmo->dHpdx_of_x(x);
  double H = cosmo->H_of_x(x);
  double c = Constants.c;
  double h = cosmo -> get_h();

  const double H0  = Constants.H0_over_h*h;

  
  

  // double Omega_Nu = cosmo->get_OmegaNu; // Ignoring Neutrinos
  double Omega_B = cosmo->get_OmegaB(0.0);
  double Omega_CDM = cosmo->get_OmegaCDM(0.0);
  double Omega_R = cosmo->get_OmegaR(0.0); 

  // Recombination variables
  double tau = rec->tau_of_x(x);
  double dtau = rec->dtaudx_of_x(x);
  double ddtau = rec->ddtauddx_of_x(x);


  double ck_Hp = c*k/Hp;
  double R = (4.0*Omega_R)/(3.0*Omega_B*a);  // R is 1/R in Dodelson

  double Theta2 = - 20.0*ck_Hp/(45.0*dtau) *Theta[1];

  double Psi = -Phi -( (12.0*H0*H0)/(c*c*k*k*a*a) ) * (Omega_R*Theta2); 
              
  //=============================================================================
  // TODO: fill in the expressions for all the derivatives
  //=============================================================================

  // SET: Scalar quantities (Phi, delta, v, ...)

  double K = ((H0*H0)/(Hp*Hp))*( Omega_CDM*(1.0/a)*delta_cdm + Omega_B*(1.0/a)*delta_b + 4.0*Omega_R*(1.0/(a*a))*Theta[0] );  

  dPhidx = Psi - (1.0/3.0)*(ck_Hp*ck_Hp)*Phi + K/2.0; 

  dThetadx[0] = -ck_Hp*Theta[1] - dPhidx; 

  ddelta_cdmdx = ck_Hp * v_cdm - 3.0*dPhidx;
  ddelta_bdx = ck_Hp * v_b - 3.0*dPhidx;

  dv_cdmdx = -v_cdm - ck_Hp*Psi;
  


  
  // Tight coupling specific 
  double q1 = -( (1.0 - R)*dtau + (1.0 + R)*ddtau ) * (3.0*Theta[1] + v_b) -ck_Hp*Psi + (1.0 - (dHp/Hp) )*ck_Hp*(-Theta[0]+2.0*Theta2) - ck_Hp*dThetadx[0];
  double q2 = (1.0 + R)*dtau + (dHp/Hp) -1.0;
  double q = q1/q2;

  dv_bdx = (1.0/(1.0 + R)) * ( -v_b - ck_Hp*Psi +R*( q + ck_Hp*(-Theta[0] + 2.0*Theta2) -ck_Hp*Psi) ) ;

  dThetadx[1] = (q - dv_bdx)/3.0 ;

  // SET: Photon multipoles (Theta_ell)
  // """""""""""""""""""""""""""" GREEN CHANGE *** """"""""""""""""""""""""""""""""""""""""""""""
  // Theta[2] = - 20.0*ck_Hp/(45.0*dtau) *Theta[1];
  // if (n_ell_theta_tc>2){ 
  //   // double l = n_ell_neutrinos_tc + 2; // *** Might need to define l as a double
  //   for(int l=3; l<n_ell_theta_tc; l++){
  //     Theta[l] = -l/(2.0*l+1.0) *ck_Hp/dtau * Theta[l-1];
  //    }
  // }

  return GSL_SUCCESS;
}

//====================================================
// The right hand side of the full ODE
//====================================================

int Perturbations::rhs_full_ode(double x, double k, const double *y, double *dydx){
  
  //=============================================================================
  // Compute where in the y / dydx array each component belongs
  // This is just an example of how to do it to make it easier
  // Feel free to organize the component any way you like
  //=============================================================================

  // Index and number of the different quantities
  const int n_ell_theta         = Constants.n_ell_theta;
  const int n_ell_thetap        = Constants.n_ell_thetap;
  const int n_ell_neutrinos     = Constants.n_ell_neutrinos;
  const bool polarization       = Constants.polarization;
  const bool neutrinos          = Constants.neutrinos;

  // The different quantities in the y array
  const double &delta_cdm       =  y[Constants.ind_deltacdm];
  const double &delta_b         =  y[Constants.ind_deltab];
  const double &v_cdm           =  y[Constants.ind_vcdm];
  const double &v_b             =  y[Constants.ind_vb];
  const double &Phi             =  y[Constants.ind_Phi];
  const double *Theta           = &y[Constants.ind_start_theta];
  const double *Theta_p         = &y[Constants.ind_start_thetap];
  const double *Nu              = &y[Constants.ind_start_nu];

  // References to the quantities we are going to set in the dydx array
  double &ddelta_cdmdx    =  dydx[Constants.ind_deltacdm];
  double &ddelta_bdx      =  dydx[Constants.ind_deltab];
  double &dv_cdmdx        =  dydx[Constants.ind_vcdm];
  double &dv_bdx          =  dydx[Constants.ind_vb];
  double &dPhidx          =  dydx[Constants.ind_Phi];
  double *dThetadx        = &dydx[Constants.ind_start_theta];
  double *dTheta_pdx      = &dydx[Constants.ind_start_thetap];
  double *dNudx           = &dydx[Constants.ind_start_nu];

  // Cosmological parameters and variables
  double a = exp(x);
  double c = Constants.c;

  double Hp = cosmo->Hp_of_x(x);
  double H = cosmo->H_of_x(x);
  
  double eta = cosmo->eta_of_x(x);

  // double Omega_Nu = cosmo->get_OmegaNu; // Ignoring Neutrinos
  double Omega_B = cosmo->get_OmegaB(0.0);
  double Omega_CDM = cosmo->get_OmegaCDM(0.0);
  double Omega_R = cosmo->get_OmegaR(0.0); 

  double h          = cosmo -> get_h();
  const double H0          = Constants.H0_over_h*h;


  // Recombination variables
  double tau = rec->tau_of_x(x);
  double dtau = rec->dtaudx_of_x(x);
  double ddtau = rec->ddtauddx_of_x(x);


  double ck_Hp = c*k/Hp;
  double R = (4.0*Omega_R)/(3.0*Omega_B*a);  // R is 1/R in Dodelson
  
  // SET: Scalar quantities (Phi, delta, v, ...)
  double Psi = -Phi -( (12.0*H0*H0)/(c*c*k*k*a*a) ) * (Omega_R*Theta[2]); 

  double K = ((H0*H0)/(Hp*Hp))*( Omega_CDM*(1.0/a)*delta_cdm + Omega_B*(1.0/a)*delta_b + 4.0*Omega_R*(1.0/(a*a))*Theta[0] );  
  dPhidx = Psi - (1.0/3.0)*(ck_Hp*ck_Hp)*Phi + K/2.0;

  ddelta_cdmdx = ck_Hp * v_cdm - 3.0*dPhidx;
  ddelta_bdx = ck_Hp * v_b - 3.0*dPhidx;

  dv_cdmdx = -v_cdm - ck_Hp*Psi;
  dv_bdx  = -v_b - ck_Hp*Psi + dtau*R * (3.0*Theta[1] + v_b) ;

   

  // SET: Photon multipoles (Theta_ell)

  double Pi = Theta[2]; 
  

  dThetadx[0] = -ck_Hp*Theta[1] - dPhidx;
  dThetadx[1] = ck_Hp/3.0 *Theta[0] - (2.0*ck_Hp/3.0)*Theta[2] + ck_Hp/3.0 *Psi + dtau*(Theta[1] + (1.0/3.0)*v_b);
  dThetadx[2] = (2.0/(2.0*2.0+1.0))*ck_Hp*Theta[2-1] - (2.0 +1)/(2.0*2.0 + 1.0)*ck_Hp*Theta[2+1] + dtau*(Theta[2]-(1.0/10.0)*Pi*1.0); // Kroenecker delta = 1


  for(int l=3; l<n_ell_theta-1; l++){
    // dThetadx[l] = (l/(2.0*l+1.0))*ck_Hp*Theta[l-1] - (l +1)/(2.0*l + 1.0)*ck_Hp*Theta[l+1] + dtau*(Theta[l]-(1.0/10.0)*Pi*0.0); // Kroenecker delta = 0
    dThetadx[l] = (l/(2.0*l+1.0))*ck_Hp*Theta[l-1] - (l +1)/(2.0*l + 1.0)*ck_Hp*Theta[l+1] + dtau*(Theta[l]); // Kroenecker delta = 0
  }

  dThetadx[n_ell_theta-1] =  ck_Hp*Theta[n_ell_theta-1-1] - c*( (n_ell_theta)/(Hp*eta) ) *Theta[n_ell_theta -1 ] + dtau * Theta[n_ell_theta-1]; 
  // dThetadx[n_ell_theta] =  ck_Hp*Theta[n_ell_theta-1] - c*( (n_ell_theta+1.0)/(Hp*eta)) *Theta[n_ell_theta ] + dtau * Theta[n_ell_theta];
  

  return GSL_SUCCESS;
}

//====================================================
// Get methods
//====================================================

double Perturbations::get_delta_cdm(const double x, const double k) const{
  return delta_cdm_spline(x,k);
}
double Perturbations::get_delta_b(const double x, const double k) const{
  return delta_b_spline(x,k);
}
double Perturbations::get_v_cdm(const double x, const double k) const{
  return v_cdm_spline(x,k);
}
double Perturbations::get_v_b(const double x, const double k) const{
  return v_b_spline(x,k);
}
double Perturbations::get_Phi(const double x, const double k) const{
  return Phi_spline(x,k);
}
double Perturbations::get_Psi(const double x, const double k) const{
  return Psi_spline(x,k);
}


double Perturbations::get_Pi(const double x, const double k) const{
  return 0.0; //Pi_spline(x,k);
}
double Perturbations::get_Source_T(const double x, const double k) const{
  return ST_spline(x,k);
}
double Perturbations::get_Source_E(const double x, const double k) const{
  return SE_spline(x,k);
}
double Perturbations::get_Theta(const double x, const double k, const int ell) const{
  return Theta_spline[ell](x,k);
}
double Perturbations::get_Theta_p(const double x, const double k, const int ell) const{
  return 0.0;//Theta_p_spline[ell](x,k);
}
double Perturbations::get_Nu(const double x, const double k, const int ell) const{
  return 0.0;//Nu_spline[ell](x,k);
}

//====================================================
// Print some useful info about the class
//====================================================

void Perturbations::info() const{
  std::cout << "\n";
  std::cout << "Info about perturbations class:\n";
  std::cout << "x_start:       " << x_start                << "\n";
  std::cout << "x_end:         " << x_end                  << "\n";
  std::cout << "n_x:           " << n_x                    << "\n";
  std::cout << "k_min (1/Mpc): " << k_min * Constants.Mpc  << "\n";
  std::cout << "k_max (1/Mpc): " << k_max * Constants.Mpc  << "\n";
  std::cout << "n_k:           " << n_k                    << "\n";
  if(Constants.polarization)
    std::cout << "We include polarization\n";
  else
    std::cout << "We do not include polarization\n";
  if(Constants.neutrinos)
    std::cout << "We include neutrinos\n";
  else
    std::cout << "We do not include neutrinos\n";

  std::cout << "Information about the perturbation system:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm         << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab           << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm             << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb               << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi              << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta      << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta          << "\n";
  if(Constants.polarization){
    std::cout << "ind_start_thetap:   " << Constants.ind_start_thetap   << "\n";
    std::cout << "n_ell_thetap:       " << Constants.n_ell_thetap       << "\n";
  }
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu       << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos    << "\n";
  }
  std::cout << "n_ell_tot_full:     " << Constants.n_ell_tot_full       << "\n";

  std::cout << "Information about the perturbation system in tight coupling:\n";
  std::cout << "ind_deltacdm:       " << Constants.ind_deltacdm_tc      << "\n";
  std::cout << "ind_deltab:         " << Constants.ind_deltab_tc        << "\n";
  std::cout << "ind_v_cdm:          " << Constants.ind_vcdm_tc          << "\n";
  std::cout << "ind_v_b:            " << Constants.ind_vb_tc            << "\n";
  std::cout << "ind_Phi:            " << Constants.ind_Phi_tc           << "\n";
  std::cout << "ind_start_theta:    " << Constants.ind_start_theta_tc   << "\n";
  std::cout << "n_ell_theta:        " << Constants.n_ell_theta_tc       << "\n";
  if(Constants.neutrinos){
    std::cout << "ind_start_nu:       " << Constants.ind_start_nu_tc    << "\n";
    std::cout << "n_ell_neutrinos     " << Constants.n_ell_neutrinos_tc << "\n";
  }
  std::cout << "n_ell_tot_tc:       " << Constants.n_ell_tot_tc         << "\n";
  std::cout << std::endl;
}

//====================================================
// Output some results to file for a given value of k
//====================================================

void Perturbations::output(const double k, const std::string filename) const{
  std::ofstream fp(filename.c_str());
  const int npts = 5000;
  auto x_array = Utils::linspace(x_start, x_end, npts);

  auto print_data = [&] (const double x) {
    double arg = k * (cosmo->eta_of_x(0.0) - cosmo->eta_of_x(x));
    fp << x                  << " ";
    fp << get_Theta(x,k,0)   << " ";
    fp << get_Theta(x,k,1)   << " ";
    fp << get_Theta(x,k,2)   << " ";
    fp << get_Phi(x,k)       << " ";
    fp << get_Psi(x,k)       << " ";
    

    fp << get_delta_cdm(x,k)   << " ";
    fp << get_v_cdm(x,k)       << " ";
    fp << get_delta_b(x,k)       << " ";
    fp << get_v_b(x,k)        << " ";


    // fp << get_Source_T(x,k)  << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(5,   arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(50,  arg)           << " ";
    // fp << get_Source_T(x,k) * Utils::j_ell(500, arg)           << " ";

    // fp << get_Pi(x,k)        << " ";
    
    fp << "\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);

}

