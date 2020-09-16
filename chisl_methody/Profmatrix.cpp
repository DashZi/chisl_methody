#include "Header.h"


MatrixProf::MatrixProf(void)
{
	n = 0;
	col = 0;
}

MatrixProf::MatrixProf(int x, int c)
{
	n = x;
	col = c;
	di.resize(n);
	ia.resize(n + 1);
	al.resize(col);
	au.resize(col);
	F.resize(n);
}


MatrixProf::~MatrixProf(void)
{
}

void MatrixProf::load(ifstream& matrix, ifstream& vect)
{
	for (int i = 0; i < n + 1; i++)
	{
		matrix >> ia[i];
	}
	for (int i = 0; i < n; i++)
	{
		matrix >> di[i];
	}
	for (int i = 0; i < col; i++)
	{
		matrix >> al[i];
	}
	for (int i = 0; i < col; i++)
	{
		matrix >> au[i];
	}
	for (int i = 0; i < n; i++)
	{
		vect >> F[i];
	}
}

void MatrixProf::LUDec() //  LU разложение
{
	for (int i = 0; i < n; i++) {
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (ia[i + 1] - ia[i]);
		dubl bdi = 0;
		for (int k = ia[i]; k < ia[i + 1]; k++, j++)
		{
			int ki = ia[i];
			int kj = ia[j];
			int dif = k - ia[i] - ia[j + 1] + ia[j];
			if (dif < 0)
				kj += abs(dif);
			else
				ki += dif;
			dubl bal = 0;
			dubl bau = 0;
			for (ki; ki < k; ki++, kj++)
			{
				bal += al[ki] * au[kj];
				bau += au[ki] * al[kj];
			}
			al[k] = al[k] - bal;
			au[k] = au[k] - bau;
			au[k] = au[k] / di[j];
			bdi += al[k] * au[k];
		}
		di[i] -= bdi;
	}
}

void MatrixProf::ToTight(Matrix* A) // перевод матрицы в плотный вид
{
	vector< vector< real > > B(n);
	for (int i = 0; i < n; i++)
	{
		vector<real> buf(n);
		B[i] = buf;
	}
	for (int i = 0; i < n; i++)
	{
		int j = i - (ia[i + 1] - ia[i]);
		for (int k = ia[i]; k < ia[i + 1]; k++, j++)
		{
			B[i][j] = al[k];
			B[j][i] = au[k];
		}
		B[i][i] = di[i];
	}
	A->setMatrix(B, n, F);
}

void MatrixProf::Direct() // прямой ход Ly=F
{
	for (int i = 0; i < n; i++)
	{
		int j = i - (ia[i + 1] - ia[i]);
		dubl sum = 0;
		for (int k = ia[i]; k < ia[i + 1]; k++, j++)
		{
			sum += F[j] * al[k];
		}
		F[i] = (F[i] - sum) / di[i];
	}
}

void MatrixProf::Reverse() // обратный ход Ux=y
{
	for (int i = n - 1; i >= 0; i--)
	{
		int j = i - (ia[i + 1] - ia[i]);
		for (int k = ia[i]; k < ia[i + 1]; k++, j++)
		{
			F[j] -= F[i] * au[k];
		}
	}
}

vector <real> MatrixProf::SLAU() // функция, которая решает
{
	LUDec();
	Direct();
	Reverse();
	return F;
}

void MatrixProf::Plus(real x)
{
	di[0] += x;
	F[0] += x;
}

void MatrixProf::setProf(int bn, int bcol, vector<real> bdi, vector<int> bia, vector<real> bal, vector<real> bau, vector<real> bF) // сет функция. мне стыдно
{
	di.resize(bn);
	ia.resize(bn + 1);
	al.resize(bcol);
	au.resize(bcol);
	F.resize(bn);
	n = bn;
	col = bcol;
	di = bdi;
	ia = bia;
	al = bal;
	au = bau;
	F = bF;
}