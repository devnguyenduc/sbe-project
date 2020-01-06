#ifndef BLOCH_H
#define BLOCH_H

#include <iostream>
#include <iomanip>
#include <complex>
#include <math.h>
using namespace std;
//=========================================================================================
//parameters
const double PI = 3.14159265358979;
const complex<double> i_ = 1i;
const double N_phonon = 1.0;
//=========================================================================================
//=========================================================================================
class SBE
{
private:
    const double hbar = 658.5;
    const double Er = 4.2;
    double Chi = 0;
    double phase = 0;
    double dephasing_time = 0;

public:
    double var_energy = 0;
    int N = 0;
    SBE()
    {
    }
    //=======================================================================================
    SBE(int N, double Energy, double Chi, double dephasing_time)
    {
        this->N = N;
        this->var_energy = 2.0 * (double)(Energy / N);
        this->Chi = Chi;
        this->phase = sqrt(this->Er) * this->var_energy / PI;
        this->dephasing_time = dephasing_time;
    }
    //---------------------------------------------------------------------------------------
    //define function g
    double fun_g(double a, double b)
    {
        if (a == b)
            return 0;
        return log((sqrt(a) + sqrt(b))) - log(abs(sqrt(a) - sqrt(b)));
    }
    //---------------------------------------------------------------------------------------
    //define function E
    complex<double> fun_E(int i, int j, int k, complex<double> **old)
    {
        complex<double> sum = 0.0;
        for (int a = 0; a < (this->N); a++)
        {
            sum += this->fun_g(k, a) * 2.0 * old[0][a];
        }
        return this->phase * sum;
    }
    //---------------------------------------------------------------------------------------
    virtual complex<double> solve(int i, int j, int k, complex<double> **old) {
        return 0;
    }
    //---------------------------------------------------------------------------------------
    //define function Omega
    complex<double> fun_Omega(int i, int j, int k, complex<double> **old)
    {
        complex<double> sum = 0;
        double time_x = -(double)(i * i) / (25. * 25.);
        //cout << time_x << endl;
        //cout << "Chi " << this->Chi << endl;
        double time_solution = sqrt(PI) * this->Chi * exp(time_x) * (1 / 50.);
        //cout << "time so" << time_solution << endl;
        // cout << time_solution << endl;
        for (int a = 0; a < (this->N); a++)
        {
            sum += this->fun_g(k, a) * old[1][a];
        }
        return time_solution + this->phase * sum / this->hbar;
    }
    //---------------------------------------------------------------------------------------
    //define function f
    double fun_f(int i, int j, int k, complex<double> **old)
    {
        // complex<double> sum = 0.0;
        // for (int a = 0; a < this->N; a++)
        // {
        //     sum += fun_g(k, a) * (H1[k][a] * old[k + 2][a] + H2[k][a] * old[k + N + 2][a]);
        // }
        // return (-2.0 * imag(this->fun_Omega(i, j, k, old) * conj(old[1][k])) + imag(1.0 * sum));
        return (-2.0 * imag(this->fun_Omega(i, j, k, old) * conj(old[1][k])));
    }
    //---------------------------------------------------------------------------------------
    //define function p
    complex<double> fun_p(int i, int j, int k, complex<double> **old)
    {
        return (
            (-i_ / this->hbar) * ((k)*var_energy - 30.0 - this->fun_E(i, j, k, old)) * old[1][k] + i_ * (1.0 - 2.0 * old[0][k]) * this->fun_Omega(i, j, k, old) - old[1][k] / dephasing_time);
    }
    //---------------------------------------------------------------------------------------
    //define function G1
    complex<double> fun_G1(int i, int j, int k, complex<double> **old)
    {
        return (i_ / hbar) * H1[j][k] * (old[0][j] * old[0][k] + N_phonon * old[0][k] - (N_phonon + 1) * old[0][j]);
    }
    //---------------------------------------------------------------------------------------
    //define function G2
    complex<double> fun_G2(int i, int j, int k, complex<double> **old)
    {
        return (i_ / hbar) * H2[j][k] * (-old[0][j] * old[0][k] - N_phonon * old[0][j] + (N_phonon + 1) * old[0][k]);
    }
    ~SBE()
    {
       // cout << "SBE is deleted!" << endl;
    }
};
//==========================================================================================
class density : public SBE
{
public:
    using SBE::SBE;
    virtual complex<double> solve(int i, int j, int k, complex<double> **old)
    {
        return this->fun_f(i, j, k, old);
    };
    ~density() {
       // cout << "delete density" << endl;    
    }
};
//==========================================================================================
class polarize : public SBE
{
public:
    using SBE::SBE;
    virtual complex<double> solve(int i, int j, int k, complex<double> **old)
    {
        return this->fun_p(i, j, k, old);
    }
    ~polarize() {
       // cout << "delete polarize" << endl;    
    }
};
//==========================================================================================
class G1 : public SBE
{
public:
    using SBE::SBE;
    virtual complex<double> solve(int i, int j, int k, complex<double> **old)
    {
        return fun_G1(i, j, k, old);
    };
    ~G1() {
        //cout << "delete G1" << endl;    
    }
};
//==========================================================================================
class G2 : public SBE
{
public:
    using SBE::SBE;
    virtual complex<double> solve(int i, int j, int k, complex<double> **old)
    {
        return fun_G2(i, j, k, old);
    };
    ~G2() {
        //cout << "delete G2" << endl;    
    }
};
#endif
