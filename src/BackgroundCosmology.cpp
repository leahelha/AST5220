#include "BackgroundCosmology.h"


//====================================================
// Constructors
//====================================================
    
BackgroundCosmology::BackgroundCosmology(
    double h, 
    double OmegaB, 
    double OmegaCDM, 
    double OmegaK,
    double Neff, 
    double TCMB) :
  h(h),
  OmegaB(OmegaB),
  OmegaCDM(OmegaCDM),
  OmegaK(OmegaK),
  Neff(Neff), 
  TCMB(TCMB)
{

  double pi = M_PI;
  


  H0 = Constants.H0_over_h * h;
  OmegaR = 2*pow(pi,2)/30 * pow((Constants.k_b* TCMB), 4)/(pow(Constants.hbar,3)*pow(Constants.c, 5))*(8*pi*Constants.G)/(3*pow(H0,2));
                                  

  OmegaNu = Neff * 7.0/8.0 * pow((4.0/11.0), 4/3) * OmegaR; 
  OmegaLambda = 1 - (OmegaK + OmegaB + OmegaCDM + OmegaR + OmegaNu);

}


// Solve the background
void BackgroundCosmology::solve(){
  Utils::StartTiming("Eta");
    

  double x_start = Constants.x_start;
  double x_end = Constants.x_end;
  double npts = Constants.npts ;//1e4;


  Vector x_array = Utils::linspace(x_start, x_end, npts);

  // The ODE for deta/dx
  ODEFunction detadx = [&](double x, const double *eta, double *detadx){

    detadx[0] = Constants.c * (1.0/Hp_of_x(x)); 

    return GSL_SUCCESS;
  };




  Vector eta_init = {0.00};
  ODESolver ode;
  ode.solve(detadx, x_array, eta_init);
  auto solution = ode.get_data();

  // eta(today) i enheter av meter eta/c i Gyr
  // std::cout << solution.back()[0] / Constants.c / (60.*60.*24*365*1e9) << " Gyr";

  size_t num_rows = solution.size();
  size_t num_columns = (num_rows > 0) ? solution[0].size() : 0;

  // Declare and initialize eta vector
  Vector eta(num_rows, 0.0); // Initialize with num_rows elements, all set to 0.0

  // Populate the eta vector with the first column of the solution
  for (size_t i = 0; i < num_rows; ++i) {
    eta[i] = solution[i][0];
  }
  // Print the shape
  //printf("Shape of the solution: %zu rows x %zu columns\n", num_rows, num_columns);
  

  
  eta_of_x_spline.create(x_array, eta, "eta_of_x");  

  Utils::EndTiming("Eta");

  // The ODE for deta/dx
  ODEFunction dtdx = [&](double x, const double *t, double *dtdx){

    dtdx[0] = 1.0/H_of_x(x);

    return GSL_SUCCESS;
  };

  Vector t_init = {0.0}; // *** Just trying something @@@

  ODESolver odet;
  odet.solve(dtdx, x_array, t_init);
  auto solution_time = odet.get_data();
  // Declare and initialize t vector
  Vector t(num_rows, 0.0); // Initialize with num_rows elements, all set to 0.0

  // Populate the eta vector with the first column of the solution
  for (size_t i = 0; i < solution_time.size(); ++i) {
    t[i] = solution_time[i][0];
  }

  t_of_x_spline.create(x_array, t, "t_of_x");
  std::cout << "cosmic time = " << solution_time.back()[0] / (60.*60.*24*365*1e9) << " Gyr";

  


}



double BackgroundCosmology::H_of_x(double x) const{
  
  double H_of_x = get_H0()*sqrt((OmegaB+OmegaCDM)*exp(-3*x)+ (OmegaR + OmegaNu)*exp(-4*x)+OmegaK*exp(-2*x)+OmegaLambda);

  return H_of_x;
}

double BackgroundCosmology::Hp_of_x(double x) const{

  // x = ln a -> a = e^x
  double Hp_of_x = exp(x) * H_of_x(x); 

  return Hp_of_x;
}

double BackgroundCosmology::dHpdx_of_x(double x) const{

  double dHdx = -(H0*exp(-4.0*x)*(2.0*OmegaK*exp(2.0*x)+(3.0*OmegaCDM+3.0*OmegaB)*exp(x)+4.0*OmegaR+4.0*OmegaNu))/(2.0*sqrt(OmegaK*exp(-2.0*x)+(OmegaCDM+OmegaB)*exp(-3.0*x)+(OmegaR+OmegaNu)*exp(-4.0*x)+OmegaLambda));
  double dHpdx_of_x = exp(x) * H_of_x(x) + exp(x) * dHdx; 

  return dHpdx_of_x;
}


double BackgroundCosmology::ddHpddx_of_x(double x) const{

  double dHdx = -(H0*exp(-4.0*x)*(2.0*OmegaK*exp(2.0*x)+(3.0*OmegaCDM+3.0*OmegaB)*exp(x)+4.0*OmegaR+4.0*OmegaNu))/(2.0*sqrt(OmegaK*exp(-2.0*x)+(OmegaCDM+OmegaB)*exp(-3.0*x)+(OmegaR+OmegaNu)*exp(-4.0*x)+OmegaLambda));
  double ddHddx = (H0*exp(-8.0*x)*(8.0*OmegaK*OmegaLambda*exp(6.0*x)+(18.0*OmegaCDM+18.0*OmegaB)*OmegaLambda*exp(5.0*x)+(32.0*OmegaLambda*OmegaR+32.0*OmegaLambda*OmegaNu+4.0*pow(OmegaK,2.0))*exp(4.0*x)+(14.0*OmegaCDM+14.0*OmegaB)*OmegaK*exp(3.0*x)+(24.0*OmegaK*OmegaR+24.0*OmegaK*OmegaNu+9.0*pow(OmegaCDM,2)+18.0*OmegaB*OmegaCDM+9.0*pow(OmegaB,2))*exp(2.0*x)+((26.0*OmegaCDM+26.0*OmegaB)*OmegaR+(26.0*OmegaCDM+26.0*OmegaB)*OmegaNu)*exp(x)+16.0*pow(OmegaR,2)+32.0*OmegaNu*OmegaR+16.0*pow(OmegaNu,2)))/(4.0*pow((OmegaK*exp(-2.0*x)+(OmegaCDM+OmegaB)*exp(-3.0*x)+(OmegaR+OmegaNu)*exp(-4.0*x)+OmegaLambda),(3.0/2.0)));

  double ddHpddx_of_x = exp(x) * (H_of_x(x) + 2.0 * dHdx + ddHddx) ;  //*** 

  return ddHpddx_of_x;
}


double BackgroundCosmology::get_OmegaB(double x) const{ 
  if(x == 0.0) return OmegaB;

  double Hx_over_H0 = H_of_x(x)/get_H0();
  double f_Hx_H0 = 1.0 / (exp(3.0*x)*(pow(Hx_over_H0, 2)));  // function of H(x) and H0 in the expression for OmegaB
  
  return OmegaB*(f_Hx_H0);
}


double BackgroundCosmology::get_OmegaR(double x) const{ 
  if(x == 0.0) return OmegaR;

  double Hx_over_H0 = H_of_x(x)/get_H0();
  double f_Hx_H0 = 1.0 / (exp(4.0*x)*(pow(Hx_over_H0, 2)));  // function of H(x) and H0 in the expression for OmegaR
  
  return OmegaR*(f_Hx_H0);
}

double BackgroundCosmology::get_OmegaNu(double x) const{ 
  if(x == 0.0) return OmegaNu;

  double Hx_over_H0 = H_of_x(x)/get_H0();
  double f_Hx_H0 = 1.0 / (exp(4.0*x)*(pow(Hx_over_H0, 2)));  // function of H(x) and H0 in the expression for OmegaNu
  
  return OmegaNu*(f_Hx_H0);
}

double BackgroundCosmology::get_OmegaCDM(double x) const{ 
  if(x == 0.0) return OmegaCDM;


  double Hx_over_H0 = H_of_x(x)/get_H0();
  double f_Hx_H0 = 1.0 / (exp(3.0*x)*(pow(Hx_over_H0, 2)));  // function of H(x) and H0 in the expression for OmegaCDM
  
  return OmegaCDM*(f_Hx_H0);
}

double BackgroundCosmology::get_OmegaLambda(double x) const{ 
  if(x == 0.0) return OmegaLambda;

  double Hx_over_H0 = H_of_x(x)/get_H0();
  double f_Hx_H0 = 1.0 / ((pow(Hx_over_H0, 2)));  // function of H(x) and H0 in the expression for OmegaLambda
  
  return OmegaLambda*(f_Hx_H0);
}

double BackgroundCosmology::get_OmegaK(double x) const{ 
  if(x == 0.0) return OmegaK;

  double Hx_over_H0 = H_of_x(x)/get_H0();
  double f_Hx_H0 = 1.0 / (exp(2.0*x)*(pow(Hx_over_H0, 2)));  // function of H(x) and H0 in the expression for OmegaK
  
  return OmegaK*(f_Hx_H0);


}
// double BackgroundCosmology::get_angular_distance_of_x(double x)const{

//   double A_sin = sin(sqrt(fabs(OmegaK)) * get_H0() * get_comoving_distance_of_x(x) / Constants.c); 
//   double B = (sqrt(fabs(OmegaK)) * get_H0() * get_comoving_distance_of_x(x) / Constants.c);
//   double A_sinh = sinh(sqrt(fabs(OmegaK)) * get_H0() * get_comoving_distance_of_x(x) / Constants.c); 
   

//   double r = 0.0; // Declare r outside of the if blocks

//   if (get_OmegaK(x) == 0)
//       r = get_comoving_distance_of_x(x);

//   else if (get_OmegaK(x) < 0)
//       r = get_comoving_distance_of_x(x) * A_sin / B;

//   else if (get_OmegaK(x) > 0)
//       r = get_comoving_distance_of_x(x) * A_sinh / B;

//   //double dL = r / exp(x);

//   double dA = exp(x)*r;
//   //std::cout << "r is " << r ;
//   return dA;
// }
    
double BackgroundCosmology::get_luminosity_distance_of_x(double x) const{
  //=============================================================================
  // My MCMC method didnt want to run after I made the dA function 
  //=============================================================================
  //double dA = get_angular_distance_of_x(x);
  
  double A_sin = sin(sqrt(fabs(OmegaK)) * get_H0() * get_comoving_distance_of_x(x) / Constants.c); 
  double B = (sqrt(fabs(OmegaK)) * get_H0() * get_comoving_distance_of_x(x) / Constants.c);
  double A_sinh = sinh(sqrt(fabs(OmegaK)) * get_H0() * get_comoving_distance_of_x(x) / Constants.c); 
   

  double r = 0.0; // Declare r outside of the if blocks

  if (get_OmegaK(x) == 0)
      r = get_comoving_distance_of_x(x);

  else if (get_OmegaK(x) < 0)
      r = get_comoving_distance_of_x(x) * A_sin / B;

  else if (get_OmegaK(x) > 0)
      r = get_comoving_distance_of_x(x) * A_sinh / B;

  double dL = r/(exp(x));

  return dL;
}

double BackgroundCosmology::get_comoving_distance_of_x(double x) const{

  double eta_0 = eta_of_x(0.0);
  double eta_ = eta_of_x(x);

  double Chi = eta_0 - eta_ ;

  return Chi;
}

double BackgroundCosmology::eta_of_x(double x) const{
  return eta_of_x_spline(x);
}

//*** I added this
double BackgroundCosmology::t_of_x(double x) const{
  return t_of_x_spline(x);
}

double BackgroundCosmology::get_H0() const{ 
  return H0; 
}

double BackgroundCosmology::get_h() const{ 
  return h; 
}

double BackgroundCosmology::get_Neff() const{ 
  return Neff; 
}

double BackgroundCosmology::get_TCMB(double x) const{ 
  if(x == 0.0) return TCMB;
  return TCMB * exp(-x); 
}

//====================================================
// Print out info about the class
//====================================================
void BackgroundCosmology::info() const{ 
  std::cout << "\n";
  std::cout << "Info about cosmology class:\n";
  std::cout << "OmegaB:      " << OmegaB      << "\n";
  std::cout << "OmegaCDM:    " << OmegaCDM    << "\n";
  std::cout << "OmegaLambda: " << OmegaLambda << "\n";
  std::cout << "OmegaK:      " << OmegaK      << "\n";
  std::cout << "OmegaNu:     " << OmegaNu     << "\n";
  std::cout << "OmegaR:      " << OmegaR      << "\n";
  std::cout << "OmegaSum:    " << OmegaB+OmegaCDM+OmegaLambda+OmegaK+OmegaR<< "\n"; //***
  std::cout << "Neff:        " << Neff        << "\n";
  std::cout << "h:           " << h           << "\n";
  std::cout << "TCMB:        " << TCMB        << "\n";

  std::cout << std::endl;
} 

//====================================================
// Output some data to file
//====================================================
void BackgroundCosmology::output(const std::string filename) const{
  const double x_min = Constants.x_start; //original -10.0 //*** @@@  changed
  const double x_max =  Constants.x_end;  //original 0.0 //*** @@@
  const int    n_pts =  100;  //original 100  //*** @@@
  
  Vector x_array = Utils::linspace(x_min, x_max, n_pts);

  std::ofstream fp(filename.c_str());
  auto print_data = [&] (const double x) {
    fp << x                  << " ";
    fp << eta_of_x(x)        << " ";
    fp << Hp_of_x(x)         << " ";
    fp << dHpdx_of_x(x)      << " ";
    fp << get_OmegaB(x)      << " ";
    fp << get_OmegaCDM(x)    << " ";
    fp << get_OmegaLambda(x) << " ";
    fp << get_OmegaR(x)      << " ";
    fp << get_OmegaNu(x)     << " ";
    fp << get_OmegaK(x)      << " ";
    fp << ddHpddx_of_x(x)        << " "; //***
    fp << t_of_x(x)        << " "; //***
    fp << get_luminosity_distance_of_x(x) << " ";
    fp <<"\n";
  };
  std::for_each(x_array.begin(), x_array.end(), print_data);
}

