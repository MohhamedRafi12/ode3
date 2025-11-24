#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"

#include <unistd.h>
#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

struct Params {
    double g;
    double m;
    double d;
    double b;
    double c;
};

/// Drag force magnitude = b v + c v^2
/// Return components of acceleration
static inline void drag_accel(double vx, double vy, double vz,
                              const Params* p,
                              double& ax, double& ay, double& az)
{
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    if (v == 0){
        ax = ay = az = 0;
        return;
    }

    double Fdrag = p->b * v + p->c * v*v;   // magnitude of drag
    double a = Fdrag / p->m;

    ax = -a * (vx / v);
    ay = -a * (vy / v);
    az = -a * (vz / v);
}

// ODE components for 3D projectile
double f_rx(double t, const vector<double>& y, void* params){ (void)t; return y[1]; }
double f_vx(double t, const vector<double>& y, void* params){
    (void)t;
    Params* p = (Params*)params;
    double ax, ay, az;
    drag_accel(y[1], y[3], y[5], p, ax, ay, az);
    return ax;
}
double f_ry(double t, const vector<double>& y, void* params){ (void)t; return y[3]; }
double f_vy(double t, const vector<double>& y, void* params){
    (void)t;
    Params* p = (Params*)params;
    double ax, ay, az;
    drag_accel(y[1], y[3], y[5], p, ax, ay, az);
    return ay;
}
double f_rz(double t, const vector<double>& y, void* params){ (void)t; return y[5]; }
double f_vz(double t, const vector<double>& y, void* params){
    (void)t;
    Params* p = (Params*)params;
    double ax, ay, az;
    drag_accel(y[1], y[3], y[5], p, ax, ay, az);
    return az - p->g;
}

/// Stop integration when ball goes below strike zone floor
double f_stop(double t, const vector<double>& y, void* params){
    (void)t; (void)params;
    return (y[4] < 0.0);  // stop when z becomes negative
}

// Simulate pitch for a given v0; return final z at the moment x passes xend
double simulate_pitch(double v0, double xend, double z0, double theta0_deg, const Params& p)
{
    double theta = theta0_deg * M_PI/180.0;

    // state vector: [x, vx, y, vy, z, vz]
    vector<double> y(6);
    y[0] = 0;
    y[1] = v0 * cos(theta);
    y[2] = 0;
    y[3] = 0;
    y[4] = z0;
    y[5] = v0 * sin(theta);

    // bundle function pointers
    vector<pfunc_t> f(6);
    f[0] = f_rx;  f[1] = f_vx;
    f[2] = f_ry;  f[3] = f_vy;
    f[4] = f_rz;  f[5] = f_vz;

    double t = 0;
    double h = 0.001;          // fine resolution
    double tmax = 2.0;         // should reach plate within ~0.4 s for ~50 m/s

    double x_prev = y[0];
    double z_prev = y[4];

    while (t < tmax)
    {
        // advance one RK4 step
        vector<double> y_next = RK4StepN(f, y, t, h, (void*)&p);

        double x_curr = y_next[0];
        double z_curr = y_next[4];

        // if we hit the ground before getting to the plate -> definitely too low
        if (z_curr < 0.0 && x_curr < xend) {
            return -1e6;  // sentinel "way too low"
        }

        // check if we passed x = xend
        if (x_prev <= xend && x_curr >= xend)
        {
            // linear interpolation in x to find z at exactly x = xend
            double alpha = (xend - x_prev) / (x_curr - x_prev);
            double z_hit = z_prev + alpha * (z_curr - z_prev);
            return z_hit;
        }

        // update state
        y = y_next;
        x_prev = x_curr;
        z_prev = z_curr;
        t += h;
    }

    // if we never reached xend within tmax, treat as too low / short
    return -1e6;
}


int main(int argc, char** argv)
{
    Params pars;
    pars.g = 9.81;
    pars.m = 0.145;
    pars.d = 0.075;
    pars.b = 1.6e-4 * pars.d;          // b = (1.6e-4) * d
    pars.c = 0.25    * pars.d * pars.d; // c = 0.25 * d^2

    double xend = 18.5;
    double z0   = 1.4;
    double theta0 = 1.0;

    bool showPlot = false;

    int c;
    while ((c = getopt(argc, argv, "x:z:t:p")) != -1) {
        switch (c) {
        case 'x': xend   = atof(optarg); break;
        case 'z': z0     = atof(optarg); break;
        case 't': theta0 = atof(optarg); break;
        case 'p': showPlot = true; break;
        }
    }

    TApplication theApp("App", &argc, argv);

    double targetZ = 0.9;
    double vlow = 20;
    double vhigh = 60;

    // binary search for v0
    for (int i = 0; i < 30; i++){
        double mid = 0.5*(vlow + vhigh);
        double zhit = simulate_pitch(mid, xend, z0, theta0, pars);

        // if pitch is too low at plate (or never makes plate), increase speed
        if (zhit < targetZ)
            vlow = mid;
        else
            vhigh = mid;
    }

    double vPitch = 0.5*(vlow + vhigh);

    // required output block (do not change)
    printf("********************************\n");
    printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n", xend, z0, theta0);
    printf("v_pitch = %lf m/s\n", vPitch);
    printf("********************************\n");

    if (showPlot){
        cout << "Press ^c to exit\n";
        theApp.SetIdleTimer(30, ".q");
        theApp.Run();
    }

    return 0;
}
