#include "Header.h"


Matrix::Matrix(int x, ifstream& vect)
{
	this->n = x;
	col = 0;
	vector< vector< real > > B(n);
	for (int i = 0; i < n; i++)
	{
		vector<real> buf(n);
		B[i] = buf;
	}
	matrix = B;
	F.resize(n);
	for (int i = 0; i < n; i++)
	{
		vect >> F[i];
	}
}

Matrix::Matrix(int x)
{
	this->n = x;
	col = 0;
	vector< vector< real > > B(n);
	for (int i = 0; i < n; i++)
	{
		vector<real> buf(n);
		B[i] = buf;
	}
	matrix = B;
	F.resize(n);
}

void Matrix::setMatrix(vector< vector<real> > A, int x, vector<real> B)
{
	n = x;
	matrix = A;
	F = B;
}

Matrix::Matrix(void)
{
}


Matrix::~Matrix(void)
{
}

vector<real> Matrix::Gauss() // метод гауса обычный
{
	for (int i = 1; i < n; i++)
		for (int j = i; j < n; j++)
		{
			real m = matrix[j][i - 1] / matrix[i - 1][i - 1];
			for (int k = 0; k < n; k++)
				matrix[j][k] = matrix[j][k] - m * matrix[i - 1][k];
			F[j] = F[j] - m * F[i - 1];
		}
	for (int k = n - 1; k >= 0; k--)
	{
		dubl buf = 0;
		for (int j = k + 1; j < n; j++)
		{
			buf += matrix[k][j] * F[j];
		}
		F[k] = F[k] - buf;
		F[k] = F[k] / matrix[k][k];
	}
	return F;
}

vector<real> Matrix::GoodGauss() // метод гаусса с ведущим элементом не работает
{
	real b = 0;
	for (int i = 0; i < n; i++)
	{
		dubl max = matrix[i][i];
		int k = i;
		for (int c = i; c < n; c++)
		{
			if (matrix[c][i] > max)
			{
				max = matrix[c][i];
				k = c;
			}
		}
		matrix[i].swap(matrix[k]);
		b = F[i];
		F[i] = F[k];
		F[k] = b;
		for (int j = i + 1; j < n; j++)
		{
			dubl m = matrix[j][i] / matrix[i][i];
			for (int k = 0; k < n; k++)
				matrix[j][k] = matrix[j][k] - m * matrix[i][k];
			F[j] = F[j] - m * F[i];
		}
	}
	for (int k = n - 1; k >= 0; k--)
	{
		dubl buf = 0;
		for (int j = k + 1; j < n; j++)
		{
			buf += matrix[k][j] * F[j];
		}
		F[k] = F[k] - buf;
		F[k] = F[k] / matrix[k][k];
	}
	return F;
}

void Matrix::Gilbert() // построение матрицы гильберта
{
	for (int i = 0; i < n; i++)
	{
		F[i] = i + 1;
		for (int j = 0; j < n; j++)
			matrix[i][j] = (double)1.0 / (i + j + 1);
	}
	GilbertVect();
}

void Matrix::ToProf(MatrixProf* A) // перевод матрицы в профильный формат
{
	getCol();
	vector <real> bdi(n);
	vector <int> bia(n + 1);
	vector <real> bau(col);
	vector <real> bal(col);
	int s = 0;
	int flag;
	for (int i = 0; i < n; i++)
	{
		bdi[i] = matrix[i][i];
		bia[i] = s;
		flag = 0;
		for (int j = 0; j < i; j++)
		{
			if (matrix[i][j] != 0 || matrix[j][i] != 0)
			{
				flag = 1;
			}
			if (flag == 1)
			{
				bau[s] = matrix[i][j];
				bal[s] = matrix[j][i];
				s++;
			}
		}
	}
	bia[n] = s;
	A->setProf(n, col, bdi, bia, bal, bau, F);
}

void Matrix::getCol() // подсчет значащих элементов в нижнем и верхнем треугольнике матрицы
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < i; j++)
		{
			if (matrix[i][j] != 0 || matrix[j][i] != 0)
				col++;
		}
}

void Matrix::GilbertVect()
{
	for (int i = 0; i < n; i++)
	{
		dubl sum = 0;
		for (int k = 0; k < n; k++)
		{
			sum += matrix[i][k] * (k + 1);
		}
		F[i] = sum;
	}
}

void Matrix::Plus(real x)
{
	matrix[0][0] += x;
	F[0] += x;
}