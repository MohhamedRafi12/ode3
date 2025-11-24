///
/// @file 
/// @brief Basic example of using RK solver in the GSL
/// Use GSL to solve the projectile motion problem
/// Below we will use the 8th order Runge-Kutta algorithm implemented in the GSL
///

#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TGClient.h"
#include <unistd.h> 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include <iostream>
using std::cout;
using std::endl;


// functions to describe simple projectile motion
// here use use ri,rj,rk to define directions to prevent confusion with
// standard ODE notation, where x=independent variable, \vec y=dependent variable(s)

// first we define our vaiables
// x    = time
// y[0] = position along i axis   ; f_ri = dri/dt => velocity along i axis  
// y[1] = velocity along i axis   ; f_vi = dvi/dt => acceleration along i axis
// y[2] = position along j axis   ; f_rj = drj/dt => velocity along j axis  
// y[3] = velocity along j axis   ; f_vj = dvj/dt => acceleration along j axis


// the function for the gsl solver holds the equations of motion of each coordinate
int func (double t, const double y[], double f[], void *params) {
  const double g = 9.81;
  double m = ((double *)params)[0];
  double air_k = ((double *)params)[1];
  double v = sqrt(y[1]*y[1] + y[3]*y[3]);
  f[0] = y[1];                         // f_ri
  f[1] = -air_k * v * y[1] / m;        // f_vi
  f[2] = y[3];                         // f_rj
  f[3] = -air_k * v * y[3] / m - g;    // f_vj
  return GSL_SUCCESS;
}



int main(int argc, char **argv){

  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  // ******************************************************************************
  // ** this block is useful for supporting both high and std resolution screens **
  UInt_t dh = gClient->GetDisplayHeight()/2;   // fix plot to 1/2 screen height  
  //UInt_t dw = gClient->GetDisplayWidth();
  UInt_t dw = 2.2*dh;
  // ******************************************************************************

  
  // m, air_k specified as parameters
  double params[]={1.0, 0.12};

  // Init the ODE solver, by defining the system of equations.
  // The second parameter below is the Jacobian (a function if provided to calcualte
  // partial derivatives of the functions).  This is not used for RK algorithms,
  // so we set the function pointer to 0.
  // The params array allows us to set parameters in the functions.
  // {func, jacoboan fcn, dimension, params (can be 0 if not used)}
  gsl_odeiv2_system sys = {func, 0, 4, params};
  

  // Choose the algorithm
  // "The driver object is a high level wrapper that combines the
  // evolution, control and stepper objects for easy use."
  // see https://www.gnu.org/software/gsl/doc/html/ode-initval.html
  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, 1e-9, 1e-9, 0.0);

  // starting coordinates
  double y[4]={0,10,0,10};   // x=0, vx=10 m/s, y=0, vy=10 m/s

  // make two graphs from the trajectory, using small time steps
  TCanvas *c2 = new TCanvas("c2","Projectile GSL",dw,dh);
  c2->Divide(2,1);
  TGraph *tgYvX = new TGraph();
  tgYvX->SetTitle("Height vs distance;Horizontal distance [m];Vertical distance [m]");
  TGraph *tgKEvt = new TGraph();
  tgKEvt->SetTitle("Kinetic E vs time;time [s];KE [J]");

  // calculate a trajectory
  double t=0;
  double dt=0.001;  // 1 millisecond time steps

  // here we integrate over successive time steps until our defined stopping condition is satisfied
  while (y[2]>=0) {
    tgYvX->SetPoint(tgYvX->GetN(),y[0],y[2]);
    double KE=0.5*params[0]*(y[1]*y[1]+y[3]*y[3]);
    tgKEvt->SetPoint(tgKEvt->GetN(),t,KE);
    // take a step from t to t+dt.  Note that t and y are updated by the step
    int status = gsl_odeiv2_driver_apply (d, &t, t+dt, y);
    if (status != GSL_SUCCESS) {
      printf ("error, return value=%d\n", status);
      return -1;
    }
  }
    
  // free the driver
  gsl_odeiv2_driver_free (d);

  c2->cd(1);
  tgYvX->Draw("AC");
  c2->cd(2);
  tgKEvt->Draw("AC");
  
  
  // save our graphs
  TFile *tf=new TFile("projGSL.root","recreate");
  tgYvX->Write();
  tgKEvt->Write();
  tf->Close();

  cout << "Press ^c to exit" << endl;
  theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
  theApp.Run();
  
  return 0;
}

