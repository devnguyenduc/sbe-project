#include <complex>
using namespace std;
complex<double> **H1;
complex<double> **H2;
#include "class.h"
#define input_len 60
//------------------------------------------------------------------------------------------
complex<double> *zeros(int num)
{
	complex<double> *x = new complex<double>[num];
	for (int i = 0; i < num; i++)
	{
		x[i] = {0, 0};
	}
	return x;
}
//------------------------------------------------------------------------------------------
double write(complex<double> *array, int m)
{
	complex<double> tot = 0.0;
	for (int j = 0; j < m; j++)
	{
		tot += sqrt(j) * array[j];
	}
	return abs(tot);
}
//------------------------------------------------------------------------------------------
void shrink(RungeKutta4 *&Rk4, int i, int j, SBE *dev, complex<double> **solution)
{
	Rk4 = new RungeKutta4(dev);
	//call method add_solution
	Rk4->add_solution(solution);
	Rk4->oneStep(i, j);
}
//------------------------------------------------------------------------------------------
int main()
{
	//number of iteration
	//so khoang chia nang luong
	cout << "Printing date and time ..." << endl;
	time_t my_time = time(NULL);
	printf("%s", ctime(&my_time));
	int N{0};
	int loop{0};
	double Chi{0.0};
	double T2_0{0.0};
	double gamma{0.0};
	double n{0.0};
	double E_max{0.0};
	double Detune{0.0};
	double omega_0{0.0};
	double T{0.0};
	double alpha{0.0};
	double h_e{0.0};

	string input[input_len];
	int a = 0;
	ifstream infile;
	infile.open("input");
	cout << "Reading input file ..." << endl;
	while (true)
	{
		string data;
		infile >> data;
		input[a] = data;
		a++;
		if (infile.eof())
			break;
	}
	infile.close();
	for (int b = 0; b < input_len; b++)
	{
		if (input[b] == "N")
			N = stoi(input[b + 1]);
		if (input[b] == "loop")
			loop = stoi(input[b + 1]);
		if (input[b] == "Chi")
			Chi = stod(input[b + 1]);
		if (input[b] == "T2_0")
			T2_0 = stod(input[b + 1]);
		if (input[b] == "gamma")
			gamma = stod(input[b + 1]);
		if (input[b] == "n")
			n = stod(input[b + 1]);
		if (input[b] == "E_max")
			E_max = stod(input[b + 1]);
		if (input[b] == "Detune")
			Detune = stod(input[b + 1]);
		if (input[b] == "omega_0")
			omega_0 = stod(input[b + 1]);
		if (input[b] == "T")
			T = stod(input[b + 1]);
		if (input[b] == "alpha")
			alpha = stod(input[b + 1]);
		if (input[b] == "h_e")
			h_e = stod(input[b + 1]);
	}
	//cout << omega_0 << endl;
	double var_energy = E_max / N;
	const double hbar = 658.5;
	const double kB = 0.086;
	double T2 = 1 / ((1.0 / T2_0) + gamma * pow(n, 1.0 / 3.0));
	double N_phonon = 1.0 / (exp(omega_0 / (kB * T)) - 1.0);
	//generate object rk4 of RungKutta class
	//create mang 2 chieu old Nx(2N+2)
	complex<double> **old = new complex<double> *[2 * N + 2];
	complex<double> **Fn = new complex<double> * [loop];
	complex<double> **Pn = new complex<double> * [loop];
	for(int i = 0; i < loop; i++){
		Fn[i] = new complex<double> [N];
		Pn[i] = new complex<double> [N];
	}
	complex<double> *Nt = new complex<double>[loop];
	complex<double> *Pt = new complex<double>[loop];
	for (int j = 0; j < (2 * N + 2); j++)
	{
		old[j] = new complex<double>[N];
		for (int k = 0; k < (N); k++)
		{
			old[j][k] = 0.0;
		}
	}
	//create H1,H2
	H1 = new complex<double> *[N];
	for (int j = 0; j < (N); j++)
	{
		H1[j] = new complex<double>[N];
		for (int k = 0; k < N; k++)
		{
			H1[j][k] = 0.0;
		}
	}
	H2 = new complex<double> *[N];
	for (int j = 0; j < (N); j++)
	{
		H2[j] = new complex<double>[N];
		for (int k = 0; k < N; k++)
		{
			H2[j][k] = 0.0;
		}
	}
	//------------------------------------------------------------------------------------------
	int save_data = 0;
	RungeKutta4 *Rk4 = NULL;
	complex<double> **neo = NULL;
	//------------------------------------------------------------------------------------------
	cout << "RUNNING CALCULATION..." << endl;
	for (int i = -75; i < 501; i += 2)
	{
		neo = new complex<double> *[2 * N + 2];
		for (int j = 0; j < (2 * N + 2); j++)
		{
			neo[j] = new complex<double>[N];
			for (int k = 0; k < (N); k++)
			{
				neo[j][k] = 0.0;
			}
		}
		//tinh H1
		//==========================================================================================
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				complex<double> power = i_ * (var_energy * (j - k) - hbar * omega_0 + i_) * (double)(-i) / hbar;
				H1[j][k] = exp(power);
			}
		}
		//tinh H2
		//------------------------------------------------------------------------------------------
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				complex<double> power = i_ * (var_energy * (j - k) + hbar * omega_0 + i_) * (double)(-i) / hbar;
				H2[j][k] = exp(power);
			}
		}
		//tinh toan chinh
		//------------------------------------------------------------------------------------------
		//cho j chay theo 2N+2 hang cua mang old
		for (int j = 0; j < 2; j++)
		{
			//==========================================================================================
			if (j == 0)
			{
				density *fn = new density(N, E_max, Chi, T2);
				shrink(Rk4, i, j, fn, old);
				copy_array(neo[j], Rk4->result, N);
				// Save data Fn
				copy_array(Fn[save_data], neo[j], N);
				// Save data Nt
				Nt[save_data] = write(neo[j], N);
				delete Rk4;
				Rk4 = NULL;
            
			}
			//==========================================================================================
			else if (j == 1)
			{
				polarize *pn = new polarize(N, E_max, Chi, T2);
				shrink(Rk4, i, j, pn, old);
				copy_array(neo[j], Rk4->result, N);
				// save data Pn
				copy_array(Pn[save_data], neo[j], N);
				// save data Pt
				Pt[save_data] = write(neo[j], N);
				delete Rk4;
				Rk4 = NULL;
				
			}
			//==========================================================================================
			//for j chay tu 2 den N+1
			else if (1 < j && j < (N + 2))
			{
				G1 *G_one = new G1(N, E_max, Chi, T2);
				shrink(Rk4, i, j - 2, G_one, old);
				copy_array(neo[j], Rk4->result, N);
				delete Rk4;
				Rk4 = NULL;
			}
			//==========================================================================================
			//for j chay tu N+2 toi 2N+1
			//final[j]=G2(), goi rk4 cho moi dong cua G2
			else if ((N + 1) < j && j < (2 * N + 2))
			{
				G2 *G_two = new G2(N, E_max, Chi, T2);
				shrink(Rk4, i, j - (N + 2), G_two, old);
				copy_array(neo[j], Rk4->result, N);
				delete Rk4;
				Rk4 = NULL;
			}
			//==========================================================================================
			//ghi gia tri sau khi tinh vong lap i vao mang final Nx(2N+2)
		}
		for (int j = 0; j < 2 * N + 2; j++)
		{
			for (int k = 0; k < N; k++)
			{
				old[j][k] = neo[j][k];
			}
		}
		//==========================================================================================
		//xoa "neo" sau vong lap i
		for (int j = 0; j < (2 * N + 2); j++)
		{
			delete[] neo[j];
		}
		delete[] neo;
		neo = NULL;
		save_data++;
	}
	cout << "END CALCULATION!" <<endl;
	//==========================================================================================
	cout << "FREE MEMORY..." << endl;
 	//free H1,H2
	for (int j = 0; j < N; j++)
	{
		delete[] H1[j];
	}
	delete[] H1;
	H1 = NULL;
	//free H2
	for (int j = 0; j < N; j++)
	{
		delete[] H2[j];
	}
	delete[] H2;
	H2 = NULL;
	//free mang old
	for (int j = 0; j < (2 * N + 2); j++)
	{
		delete[] old[j];
	}
	delete[] old;
	old = NULL;
	cout << "FREE MEMORY SUCCESSFULL!" << endl;
	cout << "WRITE DATA ..." << endl;
	string path = "output/detuning_" + to_string((int)Detune);
	save_csv(path + "_Nt.csv", Nt, loop, "");
	delete []Nt;
	save_csv(path + "_Pt.csv", Pt, loop, "");
	delete []Pt;
	save_csv(path + "_Pn.csv", Pn, loop, N , "");
	for(int i = 0; i < loop; i++){
		delete[] Pn[i];
	}
	delete[] Pn;
	save_csv(path + "_Fn.csv", Fn, loop, N, "");
	for(int i = 0; i < loop; i++){
		delete[] Fn[i];
	}
	delete[] Fn;
	cout << "WRITE DATA SUCCESSFULL!" << endl;
	cout << "END PROCESSING!" << endl;
	return 0;
}
