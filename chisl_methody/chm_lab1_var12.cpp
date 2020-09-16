#include "Header.h"

void main()
{
	setlocale(LC_CTYPE, "Russian");
	ifstream size("size.txt");
	ifstream matr("matrix.txt");
	ifstream vect("vector.txt");
	ofstream proffile("prof.txt");
	ofstream tightfile("tight.txt");
	int n, col;
	size >> n >> col;
	MatrixProf myMat(n, col);
	vector<real> X(n);
	vector<real> Y(n);
	Matrix tightMat(n);
	int flag = 0;
	try
	{
		cout << "1 - LU из файла. 2 - гильберт в плотной и в LU. 3" << endl;
		cin >> flag;
		switch (flag)
		{
		case 1:
		{
			real x = 0;
			for (int i = 0; i < 20; i++)
			{

				ifstream matr("matrix.txt");
				ifstream vect("vector.txt");
				myMat.load(matr, vect);
				x = pow(10.0, -i);
				//myMat.Plus(x);
				myMat.ToTight(&tightMat);
				X = myMat.SLAU();
				//Y = tightMat.Gauss();
				tightfile.precision(7);
				proffile.precision(7);
				for (int i = 0; i < n; i++)
				{
					//tightfile << Y[i] << endl;
					proffile << X[i] << endl;
				}
				tightfile << endl;
				proffile << endl;
			}
			break;
		}
		case 2:
		{
			for (int i = 10; i < 11; i++)
			{
				Matrix tightMat(i);
				tightMat.Gilbert();
				MatrixProf myMat(1, 1);
				tightMat.ToProf(&myMat);
				X = myMat.SLAU();
				Y = tightMat.Gauss();
				tightfile.precision(15);
				proffile.precision(15);
				for (int j = 0; j < i; j++)
				{
					tightfile << Y[j] << endl;
					proffile << X[j] << endl;
				}
				tightfile << endl;
				proffile << endl;
			}
			break;
		}
		default:
		{
			cout << "другого не дано" << endl;
			break;
		}
		}
	}
	catch (int error)
	{
		switch (error)
		{
		case 1:
		{
			cout << "что-то пошло не так";
			system("pause");
			break;
		}
		}
	}
}
