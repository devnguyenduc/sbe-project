#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>
const int Nenergy = 500, Ntime = 3, Nomega = 1000;
const double hbar = 659.1; //meV.fs
const double Er = 4.2;     //meV
const double Deltat = 50;  //fs
const double Chi = 0.1;
const double Detuning = 50; //meV
const double emax = 300;    //meV
double de;
const double tmax = 0;
double tmin;
const double T2 = 100; //fs
const double alpha = 0.069;
const double phonon = 36; //meV
const double mRatio = 6.86;
double re;
double rh;
const double gam = 1;
double C1, C2, C3, C4, C5;
double Nphonon;
double G[Nenergy + 1][Nenergy + 1];
/////////////////////////////////////////
typedef struct
{
    double Fe[Nenergy + 1];
    double Fh[Nenergy + 1];
    double complex P[Nenergy + 1];
    double complex Ge1[Nenergy + 1][Nenergy + 1];
    double complex Ge2[Nenergy + 1][Nenergy + 1];
    double complex Gh1[Nenergy + 1][Nenergy + 1];
    double complex Gh2[Nenergy + 1][Nenergy + 1];
    double complex Gp1[Nenergy + 1][Nenergy + 1];
    double complex Gp2[Nenergy + 1][Nenergy + 1];
} Data1;
typedef struct
{
    double TIME[Ntime + 1];
    double ENERGY[Nenergy + 1];
    double Fe[Nenergy + 1][Ntime + 1];
    double Fh[Nenergy + 1][Ntime + 1];
    double complex P[Nenergy + 1][Ntime + 1];
    double complex TotalP[Ntime + 1];
    double TotalFe[Ntime + 1];
    double TotalFh[Ntime + 1];
    double OMEGA[Nomega + 1], Alpha[Nomega + 1];
    double complex Spectrum[Nomega + 1];
} Data2;
/////////////////////////////////////////
void InitConstants()
{
    tmin = -4 * Deltat;
    re = 1 / (1 + 1 / mRatio);
    rh = 1 / (1 + mRatio);
    de = emax / Nenergy;
    C1 = alpha / M_PI / hbar * pow(phonon, 1.5) * sqrt(de);
    C2 = sqrt(Er * de) / M_PI / hbar;
    C3 = 0.5 * sqrt(M_PI) / Deltat * Chi;
    C4 = I * C1 / 2;
    C5 = C2 * hbar;
    Nphonon = 1 / (exp(36 / 25) - 1);
}
void CreateGmatrix(double G[Nenergy + 1][Nenergy + 1])
{
    double a, b;
    int i, j;
    for (i = 1; i <= Nenergy + 1; i++)
    {
        for (j = i + 1; j <= Nenergy + 1; j++)
        {
            a = sqrt(i);
            b = sqrt(j);
            G[i][j] = 1 / sqrt(i) * log(fabs((a + b) / (a - b)));
            G[j][i] = G[i][j];
        }
    }
}
/////////////////////////////////////////
void Zero(Data1 *R)
{
    int i, j;
    for (i = 0; i <= Nenergy; i++)
    {
        R->Fe[i] = 0;
        R->Fh[i] = 0;
        R->P[i] = 0;
        for (j = 0; j <= Nenergy; j++)
        {
            R->Ge1[i][j] = 0;
            R->Ge2[i][j] = 0;
            R->Gh1[i][j] = 0;
            R->Gh2[i][j] = 0;
            R->Gp1[i][j] = 0;
            R->Gp2[i][j] = 0;
        }
    }
}
void Copy(Data1 *R, Data1 *X)
{
    int i, j;
    for (i = 0; i <= Nenergy; i++)
    {
        R->Fe[i] = X->Fe[i];
        R->Fh[i] = X->Fh[i];
        R->P[i] = X->P[i];
        for (j = 0; j <= Nenergy; j++)
        {
            R->Ge1[i][j] = X->Ge1[i][j];
            R->Ge2[i][j] = X->Ge2[i][j];
            R->Gh1[i][j] = X->Gh1[i][j];
            R->Gh2[i][j] = X->Gh2[i][j];
            R->Gp1[i][j] = X->Gp1[i][j];
            R->Gp2[i][j] = X->Gp2[i][j];
        }
    }
}
void Combine(Data1 *X1, double x1, Data1 *X2, double x2, Data1 *R)
{
    int i, j;
    for (i = 0; i <= Nenergy; i++)
    {
        R->Fe[i] = x1 * X1->Fe[i] + x2 * X2->Fe[i];
        R->Fh[i] = x1 * X1->Fh[i] + x2 * X2->Fh[i];
        R->P[i] = x1 * X1->P[i] + x2 * X2->P[i];
        for (j = 0; j <= Nenergy; j++)
        {
            R->Ge1[i][j] = x1 * X1->Ge1[i][j] + x2 * X2->Ge1[i][j];
            R->Ge2[i][j] = x1 * X1->Ge2[i][j] + x2 * X2->Ge2[i][j];
            R->Gh1[i][j] = x1 * X1->Gh1[i][j] + x2 * X2->Gh1[i][j];
            R->Gh2[i][j] = x1 * X1->Gh2[i][j] + x2 * X2->Gh2[i][j];
            R->Gp1[i][j] = x1 * X1->Gp1[i][j] + x2 * X2->Gp1[i][j];
            R->Gp2[i][j] = x1 * X1->Gp2[i][j] + x2 * X2->Gp2[i][j];
        }
    }
}
void CreateTime(Data2 *R)
{
    int i;
    double h = (tmax - tmin) / Ntime;
    for (i = 0; i <= Ntime; i++)
    {
        R->TIME[i] = tmin + i * h;
    }
}
void CreateEnergy(Data2 *R)
{
    int i;
    for (i = 0; i <= Nenergy; i++)
    {
        R->ENERGY[i] = i * de;
    }
}
void CreateOmega(Data2 *R)
{
    int i;
    double step = 60. / Nomega;
    for (i = 0; i <= Nomega; i++)
    {
        R->OMEGA[i] = -60 + i * step; //meV
    }
}
/////////////////////////////////////////
double complex Omega(int n, double t, Data1 *X)
{
    double complex R = 0;
    int l;
    for (l = 1; l <= Nenergy; l++)
    {
        if (l != n)
        {
            R += G[n + 1][l + 1] * X->P[l];
        }
    }
    R *= C2;
    R += C3 * exp(-pow(t / Deltat, 2));
    return R;
}
double En(int n, Data1 *X)
{
    double R = 0;
    int l;
    for (l = 0; l < Nenergy; l++)
    {
        //considering removing if, since G[i][i] = 0
        if (l != n)
        {
            R += G[n + 1][l + 1] * (X->Fe[l] + X->Fh[l]);
        }
    }
    R *= C5;
    return R;
}
double complex Ke1(int n, int n1, double t)
{
    double complex R;
    R = cexp((I * (n - n1) * de * re - gam + I * phonon) * t / hbar);
    return R;
}
double complex Ke2(int n, int n1, double t)
{
    double complex R;
    R = cexp((I * (n - n1) * de * re - gam - I * phonon) * t / hbar);
    return R;
}
double complex Kh1(int n, int n1, double t)
{
    double complex R;
    R = cexp((I * (n - n1) * de * rh - gam + I * phonon) * t / hbar);
    return R;
}
double complex Kh2(int n, int n1, double t)
{
    double complex R;
    R = cexp((I * (n - n1) * de * rh - gam - I * phonon) * t / hbar);
    return R;
}
double FeSum(int n, double t, Data1 *X)
{
    double R = 0;
    int n1;
    for (n1 = 0; n1 < Nenergy; n1++)
    {
        R += cimag(Ke1(n, n1, t) * X->Ge1[n][n1] + Ke2(n, n1, t) * X->Ge2[n][n1]) * G[n + 1][n1 + 1];
    }
    return R;
}
double FhSum(int n, double t, Data1 *X)
{
    double R = 0;
    int n1;
    for (n1 = 0; n1 < Nenergy; n1++)
    {
        R += cimag(Kh1(n, n1, t) * X->Gh1[n][n1] + Kh2(n, n1, t) * X->Gh2[n][n1]) * G[n + 1][n1 + 1];
    }
    return R;
}
double complex Kp1(int n, int n1, double t)
{
    double complex R;
    R = cexp((I * (n * rh - n1 * re) * de - gam - I * Detuning - I * phonon) * t / hbar);
    return R;
}
double complex Kp2(int n, int n1, double t)
{
    double complex R;
    R = cexp((I * (n * rh - n1 * re) * de - gam - I * Detuning + I * phonon) * t / hbar);
    return R;
}
double FpSum(int n, double t, Data1 *X)
{
    double R = 0;
    int n1;
    for (n1 = 0; n1 < Nenergy; n1++)
    {
        R += (Kp1(n, n1, t) * X->Gp1[n][n1] + Kp2(n, n1, t) * X->Gp2[n][n1]) * G[n + 1][n1 + 1];
    }
    return R;
}
void Derivative(Data1 *X0, Data1 *R, double t)
{
    int n, n1;
    for (n = 0; n <= Nenergy; n++)
    {
        double complex On = Omega(n, t, X0);
        double E = En(n, X0);
        R->Fe[n] = -2 * cimag(On * conj(X0->P[n])) + C1 * FeSum(n, t, X0);
        R->Fh[n] = -2 * cimag(On * conj(X0->P[n])) + C1 * FhSum(n, t, X0);
        R->P[n] = -I / hbar * (n * de - Detuning - E) * X0->P[n] + I * (1 - X0->Fe[n] - X0->Fh[n]) * On - C4 * FpSum(n, t, X0);
        for (n1 = 0; n1 <= Nenergy; n1++)
        {
            R->Ge1[n][n1] = -1 / Ke1(n, n1, t) * I / hbar * (X0->Fe[n] * X0->Fe[n1] + Nphonon * X0->Fe[n] - (Nphonon + 1) * X0->Fe[n1] + conj(X0->P[n]) * X0->P[n1]);
            R->Ge2[n][n1] = 1 / Ke2(n, n1, t) * I / hbar * (X0->Fe[n] * X0->Fe[n1] + Nphonon * X0->Fe[n] - (Nphonon + 1) * X0->Fe[n1] + conj(X0->P[n]) * X0->P[n1]);
            R->Gh1[n][n1] = -1 / Kh1(n, n1, t) * I / hbar * (X0->Fh[n] * X0->Fh[n1] + Nphonon * X0->Fh[n] - (Nphonon + 1) * X0->Fh[n1] + conj(X0->P[n]) * X0->P[n1]);
            R->Gh2[n][n1] = 1 / Kh2(n, n1, t) * I / hbar * (X0->Fh[n] * X0->Fh[n1] + Nphonon * X0->Fh[n] - (Nphonon + 1) * X0->Fh[n1] + conj(X0->P[n]) * X0->P[n1]);
            R->Gp1[n][n1] = 1 / Kp1(n, n1, t) * I / hbar * (X0->P[n] * (X0->Fh[n1] - X0->Fe[n1] + X0->P[n1] * (X0->Fh[n] - X0->Fe[n])));
            R->Gp2[n][n1] = -1 / Kp2(n, n1, t) * I / hbar * (X0->P[n] * (X0->Fh[n1] - X0->Fe[n1] + X0->P[n1] * (X0->Fh[n] - X0->Fe[n])));
        }
    }
}
/////////////////////////////////////////
void Rk4step(Data1 *X0, Data1 *R, double t0, double t1)
{
    double h = t1 - t0;
    static Data1 K1, K2, K3, K4, TEMP;

    Derivative(X0, &K1, t0);

    Combine(X0, 1, &K1, h / 2., &TEMP);
    Derivative(&TEMP, &K2, t0 + h / 2.);

    Combine(X0, 1, &K2, h / 2., &TEMP);
    Derivative(&TEMP, &K3, t0 + h / 2.);

    Combine(X0, 1, &K3, h, &TEMP);
    Derivative(&TEMP, &K3, t0 + h);

    Combine(&K1, 1, &K2, 2, &TEMP);
    Combine(&TEMP, 1, &K3, 2, &TEMP);
    Combine(&TEMP, 1, &K4, 1, &TEMP);
    Combine(X0, 1, &TEMP, h / 6., R);
}
void WriteData2(Data1 *X, int ColIndex, Data2 *R)
{
    int i = Nenergy;
    R->TotalP[ColIndex] = -sqrt(i + 1) * X->P[i];
    R->TotalFe[ColIndex] = -sqrt(i + 1) * X->Fe[i];
    R->TotalFh[ColIndex] = -sqrt(i + 1) * X->Fh[i];
    for (i = 0; i <= Nenergy; i++)
    {
        R->Fe[i][ColIndex] = X->Fe[i];
        R->Fh[i][ColIndex] = X->Fh[i];
        R->P[i][ColIndex] = X->P[i];
        R->TotalP[ColIndex] += sqrt(i + 1) * X->P[i];
        R->TotalFe[ColIndex] += sqrt(i + 1) * X->Fe[i];
        R->TotalFh[ColIndex] += sqrt(i + 1) * X->Fh[i];
    }
}
void PrintData2(Data2 *R)
{
    int i, j;
    FILE *file;

    // Density function Fe
    file = fopen("Fe.txt", "w");
    //First line: Energies axis
    fprintf(file, "%+e \t", (double)Nenergy);
    for (i = 0; i < Nenergy; i++)
    {
        fprintf(file, "%+e \t", R->ENERGY[i]);
    }
    fprintf(file, "%+e\n", R->ENERGY[Nenergy]);
    //Other lines and time axis
    for (i = 0; i <= Ntime; i++)
    {
        fprintf(file, "%+e \t", R->TIME[i]);
        for (j = 0; j < Nenergy; j++)
        {
            fprintf(file, "%+e \t", R->Fe[j][i]);
        }
        fprintf(file, "%+e\n", R->Fe[Nenergy][i]);
        fclose(file);

        // Density function Fh
        file = fopen("Fh.txt", "w");
        //First line: Energies axis
        fprintf(file, "%+e \t", (double)Nenergy);
        for (i = 0; i < Nenergy; i++)
        {
            fprintf(file, "%+e \t", R->ENERGY[i]);
        }
        fprintf(file, "%+e\n", R->ENERGY[Nenergy]);
        //Other lines and time axis
        for (i = 0; i <= Ntime; i++)
        {
            fprintf(file, "%+e \t", R->TIME[i]);
            for (j = 0; j < Nenergy; j++)
            {
                fprintf(file, "%+e \t", R->Fh[j][i]);
            }
            fprintf(file, "%+e\n", R->Fh[Nenergy][i]);
        }
        fclose(file);

        // Polarization P
        file = fopen("P.txt", "w");
        //First line: Energies axis
        fprintf(file, "%+e \t", (double)Nenergy);
        for (i = 0; i < Nenergy; i++)
        {
            fprintf(file, "%+e \t", R->ENERGY[i]);
        }
        fprintf(file, "%+e\n", R->ENERGY[Nenergy]);

        //Other lines and time axis
        for (i = 0; i <= Ntime; i++)
        {
            fprintf(file, "%+e \t", R->TIME[i]);
            for (j = 0; j < Nenergy; j++)
            {
                fprintf(file, "%+e \t", (cabs(R->P[j][i])));
            }
            fprintf(file, "%+e\n", cabs(R->P[Nenergy][i]));
        }
        fclose(file);

        // Total Fe
        file = fopen("TotalFe.txt", "w");
        for (i = 0; i <= Ntime; i++)
        {
            fprintf(file, "%+e\t%+e\n", R->TIME[i], R->TotalFe[i]);
        }
        fclose(file);

        // Total Fh
        file = fopen("TotalFh.txt", "w");
        for (i = 0; i <= Ntime; i++)
        {
            fprintf(file, "%+e\t%+e\n", R->TIME[i], R->TotalFh[i]);
        }
        fclose(file);

        // Total P
        file = fopen("TotalP.txt", "w");
        for (i = 0; i <= Ntime; i++)
        {
            fprintf(file, "%+e\t%+e\n", R->TIME[i], cabs(R->TotalP[i]));
        }
        fclose(file);

        // Alpha and Spectrum
        file = fopen("Alpha.txt", "w");
        for (i = 0; i <= Nomega; i++)
        {
            fprintf(file, "%+e\t%+e\t%e\n", R->OMEGA[i], creal(R->Alpha[i]), cabs(R->Spectrum[i]));
        }
        fclose(file);
    }
}
void Alpha(Data2 *R)
{
    int i, j;
    for (i = 0; i <= Nomega; i++)
    {
        double complex E = 0, P = 0;
        for (j = 1; j <= Ntime; j++)
        {
            E += cexp(-pow(R->TIME[j] / Deltat, 2) + I * R->OMEGA[i] * R->TIME[j] / hbar);
            P += R->TotalP[j] * cexp(I * R->OMEGA[i] * R->TIME[j] / hbar);
        }
        R->Alpha[i] = cimag(P / E);
        R->Spectrum[i] = E;
    }
}
void Rk4run(Data2 *R)
{
    int i;
    static Data1 X0, X1;
    Zero(&X0);
    WriteData2(&X0, 0, R);
    for (i = 1; i <= Ntime; i++)
    {
        Rk4step(&X0, &X1, R->TIME[i - 1], R->TIME[i]);
        WriteData2(&X1, i, R);
        Copy(&X0, &X1);
    }
    Alpha(R);
}
int main()
{
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    InitConstants();
    CreateGmatrix(G);

    static Data2 R;
    CreateTime(&R);
    CreateEnergy(&R);
    CreateOmega(&R);
    Rk4run(&R);
    PrintData2(&R);

    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time = %f\n", cpu_time_used);
    return 0;
}
