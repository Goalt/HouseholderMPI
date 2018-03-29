#include <mpi.h>
#include <iostream> 
#include "functions.h"

using namespace std;

int main(int argc, char *argv[])
{
	int rankInGroup, sizeGroup;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rankInGroup);
	MPI_Comm_size(MPI_COMM_WORLD,&sizeGroup);

	if ((sizeGroup == 1) && (argc == 3)) {
		int N = atoi(argv[1]);
		vector<Column> columns;
		vector<double> b(N);
		vector< pair<int, double> > res;

		allocMemmory(columns, N, rankInGroup, sizeGroup);
		readMatrix(columns, b, N, string(argv[2]));
		double t1, t2, t3, t4, tmp;

		t1 = MPI_Wtime();
		makeTriangle(columns, b, N, rankInGroup, sizeGroup);
		if(checkMatrix(columns, N, rankInGroup, sizeGroup)) {
			cout << "det = 0" << endl;
			freeMemmory(columns);
			MPI_Finalize();
			return 0;
		}
		t1 = MPI_Wtime() - t1;
		MPI_Reduce(&t1, &tmp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); t1 = tmp;

		t2 = MPI_Wtime();
		backGaus(columns, b, res, N, rankInGroup, sizeGroup);
		t2 = MPI_Wtime() - t2;
		MPI_Reduce(&t2, &tmp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); t2 = tmp;

		gatherResults(res, N, rankInGroup, sizeGroup);
		readMatrix(columns, b, N, string(argv[2]));

		t3 = MPI_Wtime();
		double t =  residual(columns, b, res, N, rankInGroup, sizeGroup);
		t3 = MPI_Wtime() - t3;
		MPI_Reduce(&t3, &tmp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); t3 = tmp;
		
		cout << "Кол-во узлов: " << sizeGroup << endl;
		cout << "Размер матрицы: " << N << endl;
		cout << "Невязка: " << t << endl;
		cout << "Время прямого хода: " << t1 << endl;
		cout << "Время обратного хода: " << t2 << endl;
		cout << "Время вычисления невязки: " << t3 << endl << endl;
		// printResults(res, rankInGroup);

		freeMemmory(columns);
		MPI_Finalize();
		return 0;
	}

	if (argc != 2) {
		cout << "Usage: " << argv[0] << " [N]" << endl;
		MPI_Finalize();
		return 0;
	}

	int N = atoi(argv[1]);

	if (sizeGroup > N) {
		cout << "Too many procs" << endl;
		MPI_Finalize();
		return 0;
	}	

	vector<Column> columns;
	vector<double> b(N);
	vector< pair<int, double> > res;

	allocMemmory(columns, N, rankInGroup, sizeGroup);
	fillMas(columns, N);
	fillB(b, columns, N);

	double t1, t2, t3, t4, tmp;

	t1 = MPI_Wtime();
	makeTriangle(columns, b, N, rankInGroup, sizeGroup);
	if(checkMatrix(columns, N, rankInGroup, sizeGroup)) {
			if (rankInGroup == 0)
				cout << "det = 0" << endl;

			freeMemmory(columns);
			MPI_Finalize();
			return 0;
	}
	t1 = MPI_Wtime() - t1;
	MPI_Reduce(&t1, &tmp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); t1 = tmp;

	t2 = MPI_Wtime();
	backGaus(columns, b, res, N, rankInGroup, sizeGroup);
	t2 = MPI_Wtime() - t2;
	MPI_Reduce(&t2, &tmp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); t2 = tmp;

	gatherResults(res, N, rankInGroup, sizeGroup);
	fillMas(columns, N);
	fillB(b, columns, N);
	// printMatrix(columns, N, rankInGroup, sizeGroup);

	t3 = MPI_Wtime();
	double t =  residual(columns, b, res, N, rankInGroup, sizeGroup);
	t3 = MPI_Wtime() - t3;
	MPI_Reduce(&t3, &tmp, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); t3 = tmp;

	if (rankInGroup == 0) {
		cout << "Кол-во узлов: " << sizeGroup << endl;
		cout << "Размер матрицы: " << N << endl;
		cout << "Невязка: " << t << endl;
		cout << "Время прямого хода: " << t1 << endl;
		cout << "Время обратного хода: " << t2 << endl;
		cout << "Время вычисления невязки: " << t3 << endl << endl;
	}
	// printResults(res, rankInGroup);

	freeMemmory(columns);
	MPI_Finalize();
	return 0;
}