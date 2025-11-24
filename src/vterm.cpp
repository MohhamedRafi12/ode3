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

    vector<double> y(4);
    y[0] = 0;     // x pos
    y[1] = 0;     // vx
    y[2] = 100;   // drop from 100 m
    y[3] = 0;     // vy = 0

    double t = 0;
    double h = 0.01;

    double last_v = 0;
    for(int i=0; i<20000; i++){
        y = RK4StepN(fn, y, t, h, &p);
        t += h;
        double v = sqrt(y[1]*y[1] + y[3]*y[3]);
        if (fabs(v - last_v) < 1e-5) return v;
        last_v = v;
    }
    return last_v;
}

int main(){

    /*******************************************************
     *  PART A1 — Projectile motion WITHOUT air resistance
     *******************************************************/
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
        y[0]=0;             // x
        y[1]=50*cos(M_PI/4); // vx
        y[2]=0;             // y
        y[3]=50*sin(M_PI/4); // vy

        double t=0;
        double h=0.01;

        ofstream fout("energy_no_air.txt");
        fout << "# t  KE  PE  E_total\n";

        double m = p.m;
        while(true){
            double KE = 0.5*m*(y[1]*y[1] + y[3]*y[3]);
            double PE = m*p.g*y[2];
            fout << t << " " << KE << " " << PE << " " << KE+PE << "\n";

            y = RK4StepN(fn, y, t, h, &p);
            t += h;
            if(y[2] < 0) break;
        }
    }

    /*******************************************************
     *  PART A2 — Terminal velocity for default parameters
     *******************************************************/
    double default_m = 1.0;   // kg
    double default_k = 0.1;   // drag
    double vt = compute_terminal_velocity(default_m, default_k);

    cout << "Default terminal velocity (m=1kg, k=0.1) = " << vt << " m/s\n";

    /*******************************************************
     *  PART A3 — Sweep: vt vs mass  (1 g to 10 kg)
     *******************************************************/
    ofstream fvt("vt_vs_mass.txt");
    fvt << "# mass(kg)   vt(m/s)\n";

    for(double m = 0.001; m <= 10.0; m*=1.2){
        double vt_m = compute_terminal_velocity(m, default_k);
        fvt << m << "  " << vt_m << "\n";
        cout << "m = " << m << " kg → vt = " << vt_m << "\n";
    }

    cout << "\nData written to:\n"
         << "  energy_no_air.txt\n"
         << "  vt_vs_mass.txt\n";

    return 0;
}
