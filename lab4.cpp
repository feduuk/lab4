#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

const int N = 10;

#include <iostream>

void inversion(double **A, int N)
{
	double temp;

	double **E = new double *[N];

	for (int i = 0; i < N; i++)
		E[i] = new double[N];

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			E[i][j] = 0.0;

			if (i == j)
				E[i][j] = 1.0;
		}

	for (int k = 0; k < N; k++)
	{
		temp = A[k][k];

		for (int j = 0; j < N; j++)
		{
			A[k][j] /= temp;
			E[k][j] /= temp;
		}

		for (int i = k + 1; i < N; i++)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int k = N - 1; k > 0; k--)
	{
		for (int i = k - 1; i >= 0; i--)
		{
			temp = A[i][k];

			for (int j = 0; j < N; j++)
			{
				A[i][j] -= A[k][j] * temp;
				E[i][j] -= E[k][j] * temp;
			}
		}
	}

	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
			A[i][j] = E[i][j];

	for (int i = 0; i < N; i++)
		delete[] E[i];

	delete[] E;
}


void gaus(double **array, double *farr)
{
	double **array1 = new double*[N];
	for (int i = 0; i < N; i++)
	{
		array1[i] = new double[N];
	}
	double *farr1 = new double[N];

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			array1[i][j] = array[i][j];
		}
	}
	for (int i = 0; i < N; i++)
	{
		farr1[i] = farr[i];
	}
	double buf;
	for (int k = 0; k < N - 1; k++)
	{
		for (int i = k+1; i < N; i++)
		{ 
			buf = array1[i][k] / array1[k][k];
			for (int j = k; j < N; j++)
			{
				array1[i][j] -= buf * array1[k][j];
			}
			farr1[i] -= buf * farr1[k];
			//array1[i][k] = 0;
			//array1[i][k] = array1[i][k] - array1[i][k] * array1[k][k] / array1[k][k];
		}

	}
	

	cout << endl << "triangular system: " << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << array1[i][j] << "    ";
		}
		cout << endl;
	}
	cout << endl;
	for (int i = 0; i < N; i++)
	{
		cout << farr1[i] << endl;
	}

	double *solution = new double[N];


	double buff = 0;
	cout << endl  << "Solution: " << endl;
	for (int i = N - 1; i >= 0; i--)
	{
		buff = 0;
		for (int j = i+1; j < N; j++)
		{
			buff += array1[i][j] * solution[j];
		}

		solution[i] = (farr1[i] - buff) / array1[i][i];
	}


	for (int i = 0; i < N; i++)
	{
		cout << solution[i] << endl;
	}

	cout << endl << "Nevzka: " << endl;
	double xbuff = 0;
	double nevzka[N];
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			xbuff += array1[i][j] * solution[j];
		}
		nevzka[i] = -farr1[i] + xbuff;
		xbuff = 0;
	}
	for (int i = 0; i < N; i++)
	{
		cout << nevzka[i] << endl;
	}
	
	double norm = 0;
	for (int i = 0; i < N; i++)
	{
		if (norm < abs(nevzka[i]))
		{
			norm = abs(nevzka[i]);
		}
	}
	cout << endl << "Norma nevzki: " << norm << endl;
}
void gauss(double **a, double *y, int n)
{
	double *x, max;
	int k, index;
	const double eps = 0.00001;  // точность
	x = new double[n];
	k = 0;
	while (k < n)
	{
		// Поиск строки с максимальным a[i][k]
		max = abs(a[k][k]);
		index = k;
		for (int i = k + 1; i < n; i++)
		{
			if (abs(a[i][k]) > max)
			{
				max = abs(a[i][k]);
				index = i;
			}
		}
		// Перестановка строк
		if (max < eps)
		{
			// нет ненулевых диагональных элементов
			cout << "Решение получить невозможно из-за нулевого столбца ";
			cout << index << " матрицы A" << endl;
			return;
		}
		for (int j = 0; j < n; j++)
		{
			double temp = a[k][j];
			a[k][j] = a[index][j];
			a[index][j] = temp;
		}
		double temp = y[k];
		y[k] = y[index];
		y[index] = temp;
		// Нормализация уравнений
		for (int i = k; i < n; i++)
		{
			double temp = a[i][k];
			if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] / temp;
			y[i] = y[i] / temp;
			if (i == k)  continue; // уравнение не вычитать само из себя
			for (int j = 0; j < n; j++)
				a[i][j] = a[i][j] - a[k][j];
			y[i] = y[i] - y[k];
		}
		k++;
	}
	// обратная подстановка
	for (k = n - 1; k >= 0; k--)
	{
		x[k] = y[k];
		for (int i = 0; i < k; i++)
			y[i] = y[i] - a[i][k] * x[k];
	}
	for (int i = 0; i < N; i++)
	{
		cout << x[i] << endl;
	}
}
void pvr(double **array, double *farr)
{
	int n, i, j, k = 0;
	double eps;
	double norma;
	double w;
	double A[10][10], B[10], x[10], xn[10];
	n = N;
	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			A[i][j] = array[i][j];
		}
	}
	for (i = 0; i<n; i++)
	{
		B[i] = farr[i];
	}
	cout << ("Enter accuracy:");
	cin >> eps;
	cout << "Enter w: ";
	cin >> w;

	for (i = 0; i<10; i++)
	{
		xn[i] = 0;
		x[i] = xn[i];
	}
	do
	{
		k++;
		norma = 0;

		for (i = 0; i<n; i++)
		{
			x[i] = B[i];
			for (j = 0; j<n; j++)
			{
				if (i != j)
					x[i] = x[i] - A[i][j] * x[j];
			}
			x[i] /= A[i][i];

			x[i] = w * x[i] + (1 - w)*xn[i];

			if (fabs(x[i] - xn[i]) > norma)
				norma = fabs(x[i] - xn[i]);
			xn[i] = x[i];
		}
	} while (norma > eps);
	cout << "Number of iterations: ";
	cout << k << endl;
	cout << "Solution: " << endl;
	for (i = 0; i<n; i++)
		cout << x[i] << endl;

	
	cout << endl << "Nevzka: " << endl;
	double xbuff = 0;
	double nevzka[N];
	for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					xbuff += A[i][j] * x[j];
				}
				nevzka[i] = -B[i] + xbuff;
				xbuff = 0;
			}
	for (int i = 0; i < N; i++)
	{
		cout << nevzka[i] << endl;
	}
	double norm = 0;
	for (int i = 0; i < N; i++)
	{
		if (norm < abs(nevzka[i]))
		{
			norm = abs(nevzka[i]);
			
		}
	}
	cout << endl << "Norma nevzki: " << norm << endl;
}

int main(int argc, char *argv[])
{
	setlocale(0, "");
	double **array = new double*[N];
	for (int i = 0; i < N; i++)
	{
		array[i] = new double[N];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			int i2 = i + 1;
			int j2 = j + 1;
			if (i == j)
				array[i][j] = 1;
			else
				array[i][j] = 1.0 / (i2 + j2);
		}
	}

	cout << "matr A: " << endl;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			cout << array[i][j] << "    ";
		}
		cout << endl;
	}


	double *farr = new double[N];
	for (int i = 0; i < N; i++)
	{
		int i2 = i + 1;
		farr[i] = 1.0 / i2;
	}

	cout << endl;
	cout << "matr b: " << endl;
	for (int i = 0; i < N; i++)
	{
		cout << farr[i] << endl;
	}

	//число обусловленности
	double **matr = new double *[N];
	for (int i = 0; i < N; i++)
	{
		matr[i] = new double[N];
	}
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			matr[i][j] = array[i][j];
		}
	}
	
	inversion(matr, N);

	double norma1 = 0;
	double buff = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{ 
			buff += abs(array[i][j]);
		}
		if (norma1 < buff)
		{
			norma1 = buff;
		}
		buff = 0;
	}
	double norma2 = 0;
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			buff += abs(matr[i][j]);
		}
		if (norma2 < buff)
		{
			norma2 = buff;
		}
		buff = 0;
	}
	double obus = norma1 * norma2;
	cout << endl << "Number obuslovlenosti: " << obus << endl;






	//метод гаусса
	cout << endl << "method Gaussa: " << endl;
	gaus(array, farr);
	//gauss(array, farr, N);

	//метод пвр
	cout << endl << "method pvr: " << endl << endl;
	pvr(array, farr);





	



	for (int i = 0; i < N; i++)
	{
		delete[] array[i];
	}
	delete[] array;
	delete[] farr;
	system("pause");
	return 0;
}
