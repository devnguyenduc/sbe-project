#ifndef CLASS_H
#define CLASS_H

#include <iostream>
#include <math.h>
#include <string>
#include <complex>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "bloch.h"
#include "tools.h"
using namespace std;
//------------------------------------------------------------------------------------------
void addVector(int dim, complex<double> *vec1, complex<double> *vec2, complex<double> *&vec)
{
	for (int i = 0; i < dim; i++)
	{
		vec[i] = vec1[i] + vec2[i];
	}
}
//------------------------------------------------------------------------------------------
void matMul(int dim, double factor, complex<double> *vec_init, complex<double> *&vec)
{
	//complex<double> *temp = new complex<double>[dim];
	for (int i = 0; i < dim; i++)
	{
	    vec[i] = factor * vec_init[i];
	}
	//return temp;
}
//------------------------------------------------------------------------------------------
void zeros(int num, complex<double> *&vec)
{
	for (int i = 0; i < num; i++)
	{
		vec[i] = {0, 0};
	}
}
//------------------------------------------------------------------------------------------

class RungeKutta4
{
public:
	SBE *dev;
	int h = 2;
	complex<double> **solution;
	complex<double> *result;
	RungeKutta4(SBE *dev)
	{
		this->dev = dev;
		this->result = new complex<double>[this->dev->N];
	}
	~RungeKutta4()
	{
                //cout << "delete RungeKutta 4" << endl;
		delete this->dev;
		for (int j = 0; j < (2 * this->dev->N + 2); j++)
		{
			delete []this->solution[j];
		}
		delete[] this->solution;
		delete[] this->result;
		//	cout << "RungeKutta4 is freed" << endl;
	}
	//------------------------------------------------------------------------------------------
	void add_solution(complex<double> **solution)
	{
		// Copy a array 2 dimension with rows = 2 * N + 2, cols = N;
		int N = this->dev->N;
		this->solution = new complex<double>*[2 * N + 2];
		for(int j = 0; j < 2 * N + 2; j++){
			this->solution[j] = new complex<double>[N];
			for(int k = 0; k < N; k++){
				this->solution[j][k] = 0.0;
			}
		}
		copy_array(this->solution, solution, 2 * N + 2, N);
	}
	//------------------------------------------------------------------------------------------

	void oneStep(int i, int j)
	{
		// Get N Collums on functions
		int N = this->dev->N;
		complex<double> *k1 = new complex<double>[N];
		complex<double> *k2 = new complex<double>[N];
		complex<double> *k3 = new complex<double>[N];
		complex<double> *k4 = new complex<double>[N];

		// copy a array 2 dimension save originize data on this->solutions 2 * N + 2 rows and N cols
		complex<double> **temp = new complex<double> *[2 * N + 2];
		for (int j = 0; j < 2 * N + 2; j++)
		{
			temp[j] = new complex<double>[N];
			for(int k = 0; k < N; k++){
				temp[j][k] = 0.0;
			}
		}
		copy_array(temp, this->solution, 2 * (this->dev->N) + 2, this->dev->N);
		for (int a = 0; a < N; a++)
		{
			k1[a] = (double)this->h * this->dev->solve(i, j, a, this->solution);
		}
		for (int a = 0; a < N; a++)
		{
			//cout << "k2" << endl;
			complex<double> *nsolutions = new complex<double>[N];
                        zeros(N, nsolutions);
                        complex<double> *mat_mul = new complex<double>[N];
                        zeros(N, nsolutions);
                        matMul(N, 0.5, k1, mat_mul);
                        addVector(N, this->solution[j], mat_mul, nsolutions);
			// temp[j] will change: temp[j] copy a array 1 dimension with N cols
			// WARNINGS: this step don't change variable in this->solutions !
			copy_array(temp[j], nsolutions, N);
			delete []nsolutions;
                        delete []mat_mul;
			k2[a] = (double)this->h * this->dev->solve(i + this->h * (1 / 2), j, a, temp);
		};

		for (int a = 0; a < N; a++)
		{
			//cout << "k3" << endl;
			complex<double> *nsolutions = new complex<double>[N];
                        zeros(N, nsolutions);
                        complex<double> *mat_mul = new complex<double>[N];
                        zeros(N, mat_mul);
                        matMul(N, 0.5, k2, mat_mul);
                        addVector(N, this->solution[j], mat_mul, nsolutions);
			// temp[j] will change: temp[j] copy a array 1 dimension with N cols
			// WARNINGS: this step don't change variable in this->solutions !
			copy_array(temp[j], nsolutions, N);
			delete []nsolutions;
                        delete []mat_mul;
			k3[a] = (double)this->h * this->dev->solve(i + this->h * (1 / 2), j, a, temp);
		};

		for (int a = 0; a < N; a++)
		{
			//cout << "k4" << endl;
			complex<double> *nsolutions = new complex<double>[N];
                        zeros(N, nsolutions);
                        addVector(N, this->solution[j], k3, nsolutions);
			copy_array(temp[j], nsolutions, N);
			delete []nsolutions;
			// temp[j] will change: temp[j] copy a array 1 dimension with N cols
			// WARNINGS: this step don't change variable in this->solutions !
			k4[a] = (double)this->h * this->dev->solve(i + this->h * (1), j, a, temp);
		};

		for (int a = 0; a < N; a++)
		{
			this->result[a] = this->solution[j][a] + (k1[a] + 2.0 * k2[a] + 2.0 * k3[a] + k4[a]) / 6.0;
			//cout << this->result[i][j] << endl;
		}

		// free k1 -> k4;
		delete[] k1;
		delete[] k2;
		delete[] k3;
		delete[] k4;
		// free temp
		for (int j = 0; j < (2 * N + 2); j++)
		{
			delete []temp[j];
		}
		delete[] temp;
	}
	//------------------------------------------------------------------------------------------
	void print()
	{
		for (int i = 0; i < this->dev->N; i++)
		{
			std::cout << result[i] << '\n';
		}
	}
};
#endif
