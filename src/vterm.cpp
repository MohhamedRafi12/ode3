#include "RKn.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

struct Params {
    double g;
    double m;
    double air_k;
};

/// ri' = vx
double f_ri(double t, const vector<double>& y, void* params){
    return y[1];
}

/// vx' = drag_x / m
double f_vi(double t, const vector<double>& y, void* params){
    Params* p = (Params*)params;
    double vx = y[1];
    double vy = y[3];
    double v = sqrt(vx*vx + vy*vy);
    return -p->air_k * v * vx / p->m;
}

/// rj' = vy
double f_rj(double t, const vector<double>& y, void* params){
    return y[3];
}

/// vy' = drag_y/m − g
double f_vj(double t, const vector<double>& y, void* params){
    Params* p = (Params*)params;
    double vx = y[1];
    double vy = y[3];
    double v = sqrt(vx*vx + vy*vy);
    return -p->air_k * v * vy / p->m - p->g;
}

/// stop when projectile hits ground
double f_stop(double t, const vector<double>& y, void* params){
    return (y[2] < 0); 
}

/// Compute terminal velocity by letting projectile fall vertically
double compute_terminal_velocity(double m, double air_k){
    Params p; 
    p.g = 9.81;
    p.m = m;
    p.air_k = air_k;

    vector<pfunc_t> fn(4);
    fn[0] = f_ri;
    fn[1] = f_vi;
    fn[2] = f_rj;
    fn[3] = f_vj;

    // start high so we can actually reach vt
    vector<double> y(4);
    y[0] = 0;      // x
    y[1] = 0;      // vx
    y[2] = 3000;   // drop from 3 km
    y[3] = 0;      // vy

    double t = 0;
    double h = 0.01;

    double theoretical_vt = sqrt(p.g * p.m / p.air_k);
    double last_v = 0;

    for (int i = 0; i < 500000; ++i) {
        y = RK4StepN(fn, y, t, h, &p);
        t += h;
        double v = sqrt(y[1]*y[1] + y[3]*y[3]);

        if (fabs(v - last_v) < 1e-5 &&
            fabs(v - theoretical_vt) < 1e-3 * theoretical_vt)
            return v;

        last_v = v;
    }
    return last_v;
}

double energy_error_no_air(double h, const char* outfile = nullptr)
{
    Params p;
    p.g = 9.81;
    p.m = 1.0;
    p.air_k = 0.0;  // no air

    vector<pfunc_t> fn(4);
    fn[0]=f_ri;
    fn[1]=f_vi;
    fn[2]=f_rj;
    fn[3]=f_vj;

    vector<double> y(4);
    y[0]=0;                   // x
    y[1]=50*cos(M_PI/4);      // vx
    y[2]=0;                   // y
    y[3]=50*sin(M_PI/4);      // vy

    double t = 0.0;

    ofstream fout;
    if (outfile) {
        fout.open(outfile);
        fout << "# h = " << h << "\n";
        fout << "# t   KE   PE   E_total\n";
    }

    double m = p.m;
    double KE0 = 0.5*m*(y[1]*y[1] + y[3]*y[3]);
    double PE0 = m*p.g*y[2];
    double E0 = KE0 + PE0;

    double max_rel_err = 0.0;

    while (true) {
        double KE = 0.5*m*(y[1]*y[1] + y[3]*y[3]);
        double PE = m*p.g*y[2];
        double E  = KE + PE;

        double rel_err = fabs(E - E0) / fabs(E0);
        if (rel_err > max_rel_err) max_rel_err = rel_err;

        if (fout.is_open()) {
            fout << t << " " << KE << " " << PE << " " << E << "\n";
        }

        y = RK4StepN(fn, y, t, h, &p);
        t += h;
        if (y[2] < 0) break;
    }

    if (fout.is_open()) fout.close();
    return max_rel_err;
}

int main() {

    {
        double h_ref = 0.01;
        // write full time-series for plotting
        energy_error_no_air(h_ref, "energy_no_air.txt");
    }

    {
        vector<double> hlist = {0.2, 0.1, 0.05, 0.02, 0.01, 0.005};

        ofstream ferr("energy_error_vs_h.txt");
        ferr << "# h    max_rel_energy_error\n";

        for (double h : hlist) {
            double err = energy_error_no_air(h);
            ferr << h << "  " << err << "\n";
            cout << "h = " << h << "  max rel energy error = " << err << "\n";
        }

        ferr.close();
        cout << "\nEnergy error vs h written to energy_error_vs_h.txt\n";
    }

    double default_m = 1.0;   // kg
    double default_k = 0.1;   // drag
    double vt = compute_terminal_velocity(default_m, default_k);

    cout << "Default terminal velocity (m=1kg, k=0.1) = " << vt << " m/s\n";

    ofstream fvt("vt_vs_mass.txt");
    fvt << "# mass(kg)   vt(m/s)\n";

    for(double m = 0.001; m <= 10.0; m*=1.2){
        double vt_m = compute_terminal_velocity(m, default_k);
        fvt << m << "  " << vt_m << "\n";
        cout << "m = " << m << " kg → vt = " << vt_m << "\n";
    }

    cout << "\nData written to:\n"
         << "  energy_no_air.txt\n"
         << "  energy_error_vs_h.txt\n"
         << "  vt_vs_mass.txt\n";

    return 0;
}
