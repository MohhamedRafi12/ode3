///
/// @file 
/// @brief Test of Runge-Kutta solver for series of ODEs
/// @author Bob Hirosky
/// @date 31 Dec 2019 
/// 
/// Use the Rk4 solver for coupled ODEs to solve for projectile 
/// motion with air resistance
///
/// This is similar to the test program RKnDemo, but a low level interface is used
/// control the steps in the solution.  The coordiante values are returned at each step
/// and not TGraphs are returned.
///
/// Definition of our variables
/// x    = time <br>
/// y[0] = position along i axis  ; f_ri = dri/dt => velocity along i axis  <br>
/// y[1] = velocity along i axis  ; f_vi = dvi/dt => acceleration along i axis <br>
/// y[2] = position along j axis  ; f_rj = drj/dt => velocity along j axis <br>
/// y[3] = velocity along j axis  ; f_vj = dvj/dt => acceleration along j axis <br>


#include "RKn.hpp"
#include <unistd.h> 
#include <cstdio>
#include <iostream>

using namespace std;
const double g=9.81;    ///< acceleration [m/s^2]
const double m=1.0;     ///< mass of object [kg], nb for simple projectile motion does not depend on the mass
const double air_k=0.1; ///< constant for air resistance. mass DOES matter with air resistance

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
double f_vi(double x, const vector<double> &y, void *params){ 
  (void) x;
  return -air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[1] / m;
  // return 0;  // if no air, no forces/acceleration along i direction in this problem
}

/// \brief Change in position along \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
double f_rj(double x, const vector<double> &y, void *params=0){  
  (void) x;   // prevent unused variable warning
  return y[3];
}

/// Change in velocity along  \f$\hat j\f$ axis
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Air resistance model: F= \f$k v^2\f$
///
double f_vj(double x, const vector<double> &y, void *params){  
  (void) x;
  return -air_k * sqrt(y[1]*y[1] + y[3]*y[3]) * y[3] / m - g;
  // return g;    // if no air constant acceleration along -j direction: F/m = -g
}

/// \brief Stopping condition
/// \param[in] x independent variable
/// \param[in] y dependent variables
///
/// Returns 0(1) to flag continuation(termination) of calculation 
double f_stop(double x, const vector<double> &y){
  (void) x;
  if (y[2]<0) return 1;  // stop calulation if the current step takes height to negative value
  return 0;  // continue calculation
}
/// \brief Use RK4 method to describe simple projectile motion.
int main(int argc, char **argv){

  // *** test 2: Use RK4SolveN to calculate simple projectile motion
  vector<pfunc_t> v_fun(4);   // 4 element vector of function pointers
  v_fun[0]=f_ri;
  v_fun[1]=f_vi;
  v_fun[2]=f_rj;
  v_fun[3]=f_vj;
  vector<double> y(4);
  // initial conditions are starting position, velocity and angle, equivalently ri,rj,vi,vj
  y[0]=0;   // init position on i-axis
  y[1]=70;  // init velocity along i axis
  y[2]=0;   // repeat for j-axis
  y[3]=70;

  double t=0;
  double h=0.01;  // step size
  do {
    printf("%10.3lf\t%10.3lf%10.3lf\n",t,y[0],y[2]);
    y=RK4StepN(v_fun,y,t,h);
    t+=h;
  } while (f_stop(t,y)==0);
  
}

