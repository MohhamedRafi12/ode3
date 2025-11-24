///
/// @file 
/// @brief Test of Runge-Kutta solver for series of ODEs
/// @author Bob Hirosky
/// @date 31 Dec 2019 
/// 
/// Use the Rk4 solver for coupled ODEs to solve for projectile 
/// motion with air resistance
///
/// Definition of our variables
/// x    = time <br>
/// y[0] = position along i axis  ; f_ri = dri/dt => velocity along i axis  <br>
/// y[1] = velocity along i axis  ; f_vi = dvi/dt => acceleration along i axis <br>
/// y[2] = position along j axis  ; f_rj = drj/dt => velocity along j axis <br>
/// y[3] = velocity along j axis  ; f_vj = dvj/dt => acceleration along j axis <br>


#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include <unistd.h> 
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;


struct Params {
  double g;   ///< acceleration [m/s^2]
  double m;   ///< mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double air_k;  ///< constant for air resistance. mass DOES matter with air resistance
} ;

// functions to describe simple projectile motion
// here use use ri,rj,rk to define directions to prevent confusion with
// standard ODE notation, where x=independent variable, \vec y=dependent variable(s)


/// \brief Change in position along \f$\hat i\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_ri(double x, const vector<double> &y, void *params=0){ 
  (void) x;   // prevent unused variable warning
  return y[1];
}

/// \brief Change in velocity along  \f$\hat i\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_vi(double x, const vector<double> &y, void *params=0){ 
  (void) x;
  Params *p = (Params*)params;
  return -p->air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[1] / p->m;
  // return 0;  // if no air, no forces/acceleration along i direction in this problem
}

/// \brief Change in position along \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Air resistance model: F= \f$k v^2\f$
///
double f_rj(double x, const vector<double> &y, void *params=0){  
  (void) x;   // prevent unused variable warning
  return y[3];
}

/// Change in velocity along  \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_vj(double x, const vector<double> &y, void *params=0){  
  (void) x;
  Params *p = (Params*)params;
  return -p->air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[3] / p->m - p->g;
  // return -g;    // if no air constant acceleration along -j direction: F/m = -g
}

/// \brief Stopping condition
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Returns 0(1) to flag continuation(termination) of calculation 
double f_stop(double x, const vector<double> &y, void *params=0){
  (void) x;
  if (y[2]<0) return 1;  // stop calulation if the current step takes height to negative value
  return 0;  // continue calculation
}
/// \brief Use RK4 method to describe simple projectile motion.
int main(int argc, char **argv){

  // setup default parameters
  Params pars;
  pars.g=9.81;
  pars.m=10.0;
  pars.air_k=0.1;
  void *p_par = (void*) &pars;

  double theta=45;   // initial angle degrees
  double v0=100;     // m/s
  
  int c;
  while ((c = getopt (argc, argv, "v:t:m:k:")) != -1)
    switch (c) {
    case 'v':
      v0 = atof(optarg);
      break;
    case 't':
      theta = atof(optarg);
      break;
    case 'm':
      pars.m = atof(optarg);
      break;
    case 'k':
      pars.air_k = atof(optarg);
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  
  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 1.1*dh;
  // ******************************************************************************

  // *** test 2: Use RK4SolveN to calculate simple projectile motion
  vector<pfunc_t> v_fun(4);   // 4 element vector of function pointers
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;

  vector<double> y(4);
  // initial conditions are starting position, velocity and angle, equivalently ri,rj,vi,vj
  y[0]=0;   // init position on i-axis
  y[1]=v0*cos(theta*3.14159/180);  // init velocity along i axis
  y[2]=0;   // repeat for j-axis
  y[3]=v0*sin(theta*3.14159/180);
  cout << "Vinit: " << v0 << " m/s" << endl;
  cout << "Angle: " << theta << " deg" << endl;
  cout << "(vx,vy) " << y[1] << " , "  <<  y[3] << " m/s" << endl;

  
  double x=0;           // t0
  double xmax=20;  // tmax
  int nsteps=200;
  // fixed step size algorithm
  auto tgN = RK4SolveN(v_fun, y, nsteps, x, xmax, p_par, f_stop);
  // example of variable step algorithm, here the estimate accuracy is limited to 1e-4
  // in the plot you will see the change in step size thoughout the time interval
  //auto tgN = RK4SolveNA(v_fun, y, nsteps, x, xmax, p_par, f_stop,1e-4);
  
  TCanvas *c2 = new TCanvas("c2","ODE solutions 2",dw,dh);
  tgN[2].Draw("al*");
  c2->Draw();

  cout << "Final velocity = " << sqrt(y[1]*y[1]+y[3]*y[3]) << endl;
  
  // save our graphs
  TFile *tf=new TFile("RKnDemo.root","recreate");
  for (unsigned i=0; i<v_fun.size(); i++){
    tgN[i].Write();
  }
  tf->Close();

  
  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
}

