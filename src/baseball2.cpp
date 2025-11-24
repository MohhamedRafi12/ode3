///
/// baseball2.cpp – Simulations of Baseball Pitches (Fitzpatrick style)
///
/// Pitches:
///   ip = 0 : slider
///   ip = 1 : curveball
///   ip = 2 : screwball
///   ip = 3 : fastball
///
/// All pitches start at x=y=z=0 (m), with x toward the plate,
/// y = horizontal (left/right), z = vertical (height).
///
/// We integrate in SI units and convert to feet only for printing/plotting.
///

#include "RKn.hpp"

#include "TLine.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1F.h"

#include <unistd.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace std;

struct Params {
    double g;        // gravity [m/s^2]
    double m;        // mass [kg]
    double Cdrag;    // quadratic drag coefficient [kg/m]
    double S0;       // Magnus coefficient [kg/m]
    double omega;    // spin magnitude [rad/s]
    double phi;      // spin-axis angle (deg) from +z toward +y
    double x_end_m;  // plate distance [m]
};

static const double FT_TO_M  = 0.3048;
static const double M_TO_FT  = 1.0 / FT_TO_M;
static const double MPH_TO_MPS = 0.44704;

// -------------------------------------------------------------
// Helper: compute accelerations with drag + Magnus + gravity
// State ordering: y[0]=x, y[1]=vx, y[2]=y, y[3]=vy, y[4]=z, y[5]=vz
// -------------------------------------------------------------
static inline void compute_accel(const vector<double>& y,
                                 const Params* p,
                                 double& ax, double& ay, double& az)
{
    double vx = y[1];
    double vy = y[3];
    double vz = y[5];

    // velocity magnitude
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    if (v == 0.0) {
        ax = ay = az = 0.0;
        return;
    }

    // --- Quadratic drag: F_D = -Cdrag * v * vhat * v
    //     (i.e. magnitude Cdrag * v^2, direction -vhat)
    double fdx = -p->Cdrag * v * vx;
    double fdy = -p->Cdrag * v * vy;
    double fdz = -p->Cdrag * v * vz;

    // --- Magnus force: F_M = S0 * (omega x v)
    double phi_rad = p->phi * M_PI / 180.0;
    // spin axis has no x-component, lies in yz-plane
    double wy = p->omega * sin(phi_rad);
    double wz = p->omega * cos(phi_rad);
    double wx = 0.0;

    // cross product w x v
    double w_cross_v_x = wy * vz - wz * vy;
    double w_cross_v_y = wz * vx - wx * vz;
    double w_cross_v_z = wx * vy - wy * vx;

    double fmx = p->S0 * w_cross_v_x;
    double fmy = p->S0 * w_cross_v_y;
    double fmz = p->S0 * w_cross_v_z;

    // gravity
    double fgx = 0.0;
    double fgy = 0.0;
    double fgz = -p->m * p->g;

    // total force and resulting acceleration
    double Fx = fdx + fmx + fgx;
    double Fy = fdy + fmy + fgy;
    double Fz = fdz + fmz + fgz;

    ax = Fx / p->m;
    ay = Fy / p->m;
    az = Fz / p->m;
}

// -------------------------------------------------------------
// ODE components
// -------------------------------------------------------------
double f_x(double t, const vector<double>& y, void* params) { return y[1]; }
double f_vx(double t, const vector<double>& y, void* params)
{
    Params* p = (Params*)params;
    double ax, ay, az;
    compute_accel(y, p, ax, ay, az);
    return ax;
}
double f_y(double t, const vector<double>& y, void* params) { return y[3]; }
double f_vy(double t, const vector<double>& y, void* params)
{
    Params* p = (Params*)params;
    double ax, ay, az;
    compute_accel(y, p, ax, ay, az);
    return ay;
}
double f_z(double t, const vector<double>& y, void* params) { return y[5]; }
double f_vz(double t, const vector<double>& y, void* params)
{
    Params* p = (Params*)params;
    double ax, ay, az;
    compute_accel(y, p, ax, ay, az);
    return az;
}

// -------------------------------------------------------------
// Set up initial conditions and pitch parameters for ip = 0..3
// -------------------------------------------------------------
void SetupPitch(int ip, vector<double>& y0, Params& pars, TString& title)
{
    // common physical constants
    pars.g = 9.81;
    pars.m = 0.145;      // kg
    double R = 0.0366;   // radius [m] ~ 1.44 in
    double rho = 1.2;    // air density [kg/m^3]
    double Cd  = 0.35;   // drag coefficient
    double A   = M_PI * R * R;

    pars.Cdrag = 0.5 * rho * Cd * A;   // ~7.5e-4 kg/m
    pars.S0    = 4.1e-4;               // tuned Magnus coefficient

    pars.x_end_m = 60.0 * FT_TO_M;     // 60 ft to plate

    // state vector: x, vx, y, vy, z, vz
    y0.assign(6, 0.0);

    double v0_mph;
    double theta_deg = 1.0;
    double phi_deg;
    double spin_rpm = 1800.0;

    if (ip == 0) {
        title = "Slider";
        v0_mph  = 85.0;
        phi_deg =   0.0;
    }
    else if (ip == 1) {
        title = "Curveball";
        v0_mph  = 85.0;
        phi_deg =  45.0;
    }
    else if (ip == 2) {
        title = "Screwball";
        v0_mph  = 85.0;
        phi_deg = 135.0;
    }
    else { // ip == 3
        title = "Fastball";
        v0_mph  = 95.0;
        phi_deg = 225.0;
    }

    double v0   = v0_mph * MPH_TO_MPS;
    double theta = theta_deg * M_PI / 180.0;

    pars.phi   = phi_deg;
    pars.omega = spin_rpm * 2.0 * M_PI / 60.0; // rad/s

    // initial conditions
    y0[0] = 0.0;              // x
    y0[2] = 0.0;              // y (lateral)
    y0[4] = 0.0;              // z (vertical)

    y0[1] = v0 * cos(theta);  // vx (toward plate)
    y0[3] = 0.0;              // vy
    y0[5] = v0 * sin(theta);  // vz (upwards)
}

// -------------------------------------------------------------
// main
// -------------------------------------------------------------
int main(int argc, char **argv)
{
    vector<double> y0(6);

    bool showPlot = true;
    int ip = 1;  // default: curveball

    int opt;
    while ((opt = getopt(argc, argv, "p:n")) != -1) {
        switch(opt) {
        case 'p':
            ip = atoi(optarg);
            break;
        case 'n':
            showPlot = false;
            break;
        }
    }

    Params pars;
    TString title;
    SetupPitch(ip, y0, pars, title);

    TApplication theApp("App", &argc, argv); // ROOT app

    // ODE system
    vector<pfunc_t> f(6);
    f[0] = f_x;  f[1] = f_vx;
    f[2] = f_y;  f[3] = f_vy;
    f[4] = f_z;  f[5] = f_vz;

    // integration variables
    double t    = 0.0;
    double h    = 1.0e-4;     // Fitzpatrick's time step
    double tmax = 1.0;        // should reach plate before this

    vector<double> y = y0;
    vector<double> y_prev = y0;
    double t_prev = t;

    // graphs for vertical (z) and horizontal (y) displacement vs x
    TGraph* g_vert = new TGraph();
    TGraph* g_lat  = new TGraph();

    int npts = 0;

    // integrate until x reaches plate or time runs out
    while (t < tmax && y[0] < pars.x_end_m) {
        // store current point (in feet)
        double x_ft = y[0] * M_TO_FT;
        double y_ft = y[2] * M_TO_FT;  // lateral
        double z_ft = y[4] * M_TO_FT;  // vertical

        g_vert->SetPoint(npts, x_ft, z_ft);
        g_lat ->SetPoint(npts, x_ft, y_ft);
        ++npts;

        // step
        y_prev = y;
        t_prev = t;
        y = RK4StepN(f, y, t, h, (void*)&pars);
        t += h;
    }

    // Interpolate to exactly x_end if we over-shot
    double xend_m  = pars.x_end_m;
    double xend_ft = xend_m * M_TO_FT;

    double x_prev = y_prev[0];
    double x_curr = y[0];
    double alpha  = 0.0;

    if (x_curr != x_prev) {
        alpha = (xend_m - x_prev) / (x_curr - x_prev);
        if (alpha < 0.0) alpha = 0.0;
        if (alpha > 1.0) alpha = 1.0;
    }

    vector<double> y_end(6);
    for (int i = 0; i < 6; ++i) {
        y_end[i] = y_prev[i] + alpha * (y[i] - y_prev[i]);
    }

    double yend_ft  = y_end[2] * M_TO_FT;
    double zend_ft  = y_end[4] * M_TO_FT;
    double vxend_ft = y_end[1] * M_TO_FT;
    double vyend_ft = y_end[3] * M_TO_FT;
    double vzend_ft = y_end[5] * M_TO_FT;

    // also add the final point to the graphs
    g_vert->SetPoint(npts, xend_ft, zend_ft);
    g_lat ->SetPoint(npts, xend_ft, yend_ft);

    // print results (in feet / ft/s)
    printf("********************************\n");
    printf("Pitch type ip = %d (%s)\n", ip, title.Data());
    printf("Coordinates when x = 60 ft\n");
    printf("(x,y,z)   = (%lf,%lf,%lf)   [ft]\n",
           xend_ft, yend_ft, zend_ft);
    printf("(vx,vy,vz) = (%lf,%lf,%lf) [ft/s]\n",
           vxend_ft, vyend_ft, vzend_ft);
    printf("********************************\n");

    // ---------------------------------------------------------
    // Plot: both vertical z(x) and lateral y(x) on same axes
    // ---------------------------------------------------------
    if (showPlot) {
        gStyle->SetOptStat(0);

        TCanvas* c1 = new TCanvas("c1", title, 800, 600);

        // set a frame roughly matching Fitzpatrick's plots
        TH1F* frame = c1->DrawFrame(0.0, -4.0, 60.0, 2.0);
        frame->SetTitle(Form("%s trajectory; x (ft); y,z (ft)", title.Data()));

        // style the curves: solid = vertical (z), dashed = lateral (y)
        g_vert->SetLineColor(kBlack);
        g_vert->SetLineWidth(3);

        g_lat->SetLineColor(kBlack);
        g_lat->SetLineStyle(2);
        g_lat->SetLineWidth(3);

        g_vert->Draw("L SAME");
        g_lat->Draw("L SAME");

        // vertical dotted lines at release and plate (x=0,60)
        TLine* l0  = new TLine(0.0, -4.0, 0.0, 2.0);
        TLine* l60 = new TLine(60.0,-4.0,60.0,2.0);
        l0->SetLineStyle(3);
        l60->SetLineStyle(3);
        l0->Draw();
        l60->Draw();

        // legend
        TLegend* leg = new TLegend(0.60, 0.80, 0.89, 0.90);
        leg->AddEntry(g_vert, "vertical z(x)", "l");
        leg->AddEntry(g_lat,  "horizontal y(x)", "l");
        leg->Draw();

        c1->SaveAs(Form("pitch_%d.pdf", ip));

        cout << "Press ^C to exit" << endl;
        theApp.SetIdleTimer(30, ".q");
        theApp.Run();
    }

    return 0;
}
