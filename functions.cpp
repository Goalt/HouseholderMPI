#ifndef FUNCTIONS_CPP
#define FUNCTIONS_CPP
#include "functions.h"

using namespace std;

void fillMas(std::vector<Column> &columns, int N) {
	for(int i = 0; i < columns.size(); i++) {
		for(int k = 0; k < N; k++) {
			if(columns[i].j == k)
				columns[i].mas[k] = 1.0/double(1 + columns[i].j + k) + 0.1;
			else
				columns[i].mas[k] = 1.0/double(1 + columns[i].j + k);
			// if(columns[i].j == k)
			// 	columns[i].mas[k] = 1;
			// else
			// 	columns[i].mas[k] = 0;
		}
	}
}

void fillB(std::vector<double> &b, std::vector<Column> &columns, int N) {
	for(int j = 0; j < N; j++) {
		double sum = 0;
		for(int i = 0; i < columns.size(); i++)
			if(columns[i].j % 2)
				sum += columns[i].mas[j];
		double res;
		MPI_Allreduce(&sum, &res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		b[j] = res;

	}
	// for(int k = 0; k < b.size(); k++) {
	// 	b[k] = k;
	// }
}

bool cmpColumns (Column a, Column b) { 
	return a.j < b.j;
}

void printMatrix(std::vector<Column> columns, int N, int rankInGroup, int sizeGroup) {
	if(rankInGroup == 0) {
		for(int j = 1; j < sizeGroup; j++) {
			int columnsRec;
			MPI_Recv(&columnsRec, 1, MPI_INT, j, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			for(int k = 0; k < columnsRec; k++) {
				Column tmp;
				tmp.mas = (double*)malloc(N * sizeof(double));
				MPI_Recv(tmp.mas, N, MPI_DOUBLE, j, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&(tmp.j), 1, MPI_INT, j, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				columns.push_back(tmp);
			}
		}

		std::sort (columns.begin(), columns.end(), cmpColumns);
		for(int i = 0; i < N; i++) {
			for(int j = 0; j < columns.size(); j++)
				if(fabs(columns[j].mas[i]) < EPS)
					cout << 0 << " ";
				else
					cout << columns[j].mas[i] << " ";
			cout << endl;
		}
	} else {
		int size = columns.size();
		MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		for(int i = 0; i < size; i++) {
			MPI_Send(columns[i].mas, N, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&(columns[i].j), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		}
	}
}

void makeTriangle(std::vector<Column> &columns, std::vector<double> &b, int N, int rankInGroup, int sizeGroup) {
	int cntProc = 0;
	int cntColumn = 0;

	for (int i = 0; i < N-1; ++i) {
		if(cntProc == rankInGroup) {
			// Составить x
			double S = 0.0;
			for (int j = i+1; j < N; ++j)
				S += columns[cntColumn].mas[j]*columns[cntColumn].mas[j];

			double a = sqrt(S + columns[cntColumn].mas[i]*columns[cntColumn].mas[i]);
			
			double* x;
			x = (double*) malloc(N*sizeof(double));
			memset(x, 0, N*sizeof(double));

			x[i] = columns[cntColumn].mas[i] - a;
			for (int j = i+1; j < N; ++j)
				x[j] = columns[cntColumn].mas[j];

			double normX = sqrt(S + x[i]*x[i]);

			if(fabs(normX) > EPS)
				for (int j = 0; j < N; ++j)
					x[j] /= normX;

			MPI_Bcast(x, N, MPI_DOUBLE, cntProc, MPI_COMM_WORLD);
			
			// Произвести перемножение
			for(int j = 0; j < columns.size(); j++) {
				double sum = 0;
				for(int k = 0; k < N; ++k)
					sum += x[k] * columns[j].mas[k];

				for(int k = 0; k < N; k++)
					columns[j].mas[k] -= 2 * x[k] * sum;
			}

			double sum = 0;
			for(int j = 0; j < N; j++)
				sum += x[j] * b[j];

			for(int j = 0; j < N; j++)
				b[j] -= 2 * sum * x[j];

			free(x);
			cntColumn++;
		} else {
			double* x;
			x = (double*) malloc(N*sizeof(double));
			memset(x, 0, N*sizeof(double));
			
			MPI_Bcast(x, N, MPI_DOUBLE, cntProc, MPI_COMM_WORLD);
			
			for(int j = 0; j < columns.size(); j++) {
				double sum = 0;
				for(int k = 0; k < N; ++k)
					sum += x[k] * columns[j].mas[k];

				for(int k = 0; k < N; k++)
					columns[j].mas[k] -= 2 * x[k] * sum;
			}

			double sum = 0;
			for(int j = 0; j < N; j++)
				sum += x[j] * b[j];

			for(int j = 0; j < N; j++)
				b[j] -= 2 * sum * x[j];

			free(x);
		}

		cntProc = (cntProc + 1) % sizeGroup;
	}
}

// void backGaus(std::vector<Column> &columns, std::vector<double> &b, std::vector<std::pair<int, double> > &res, int N, int rankInGroup, int sizeGroup) {
// 	int cntColumn = columns.size()-1;
// 	int reciever = rankInGroup - 1 < 0 ? sizeGroup - 1 : rankInGroup - 1;
	
// 	std::vector< std::pair<int, double> > results;

// 	if(N-1 == columns[cntColumn].j) {
// 		double x = b[columns[cntColumn].j] / columns[cntColumn].mas[columns[cntColumn].j];
// 		results.push_back(std::pair<int, double>(columns[cntColumn].j, x));
// 		for(int j = 0; j < N; j++)
// 			b[j] -= x * columns[cntColumn].mas[j];
// 		cntColumn--;
// 		MPI_Send(b.data(), b.size(), MPI_DOUBLE, reciever, 0, MPI_COMM_WORLD);
// 	}

// 	while(cntColumn >= 0) {
// 		MPI_Recv(b.data(), b.size(), MPI_DOUBLE, (rankInGroup + 1) % sizeGroup, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
// 		double x = b[columns[cntColumn].j] / columns[cntColumn].mas[columns[cntColumn].j];
// 		results.push_back(std::pair<int, double>(columns[cntColumn].j, x));
// 		for(int j = 0; j < N; j++)
// 			b[j] -= x * columns[cntColumn].mas[j];
// 		cntColumn--;

// 		if (columns[cntColumn+1].j == 0)
// 			break;
// 		MPI_Send(b.data(), b.size(), MPI_DOUBLE, reciever, 0, MPI_COMM_WORLD);
// 	}

// 	if(rankInGroup == 0) {
// 		for (int i = 1; i < sizeGroup; ++i){
// 			int cntX;
// 			MPI_Recv(&cntX, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
// 			for (int j = 0; j < cntX; ++j) {
// 				int numx;
// 				double x;
// 				MPI_Recv(&numx, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
// 				MPI_Recv(&x, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
// 				results.push_back(std::pair<int, double>(numx, x));	
// 			}
// 		}
// 	} else {
// 		int size = results.size();
// 		MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
// 		for(int i = 0; i < results.size(); i++) {
// 			MPI_Send(&(results[i].first), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
// 			MPI_Send(&(results[i].second), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
// 		}
// 	}

// 	res = results;
// }


void gatherResults(std::vector<std::pair<int, double> > &res, int N, int rankInGroup, int sizeGroup) {
	if(rankInGroup == 0) {
		for (int i = 1; i < sizeGroup; ++i){
			int cntX;
			MPI_Recv(&cntX, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (int j = 0; j < cntX; ++j) {
				int numx;
				double x;
				MPI_Recv(&numx, 1, MPI_INT, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&x, 1, MPI_DOUBLE, i, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				res.push_back(std::pair<int, double>(numx, x));	
			}
		}
	} else {
		int size = res.size();
		MPI_Send(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		for(int i = 0; i < res.size(); i++) {
			MPI_Send(&(res[i].first), 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&(res[i].second), 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
}

void backGaus(std::vector<Column> &columns, std::vector<double> &b, std::vector<std::pair<int, double> > &res, int N, int rankInGroup, int sizeGroup) {
	int cntStep = N-1;
	int cntColumn = columns.size() - 1;
	int reciever = cntStep%sizeGroup;

	while(cntStep >= 0) {

		if((cntColumn >= 0) && (columns[cntColumn].j == cntStep)) {
			double x = b[cntStep];
			for(int k = 0; k < res.size(); k++)
				x -= res[k].second * columns[res[k].first/sizeGroup].mas[cntStep];
			
			double tmpX;
			MPI_Reduce(&x, &tmpX, 1, MPI_DOUBLE, MPI_SUM, reciever, MPI_COMM_WORLD);
			x = tmpX;

			x /= columns[cntColumn].mas[cntStep];

			res.push_back(std::pair<int, double>(cntStep, x));
			cntColumn--;
		} else if(res.size()){
			double sendX = 0;
			for(int k = 0; k < res.size(); k++)
				sendX += res[k].second * columns[res[k].first/sizeGroup].mas[cntStep];
			sendX *= -1;

			double tmpSendX;
			MPI_Reduce(&sendX, &tmpSendX, 1, MPI_DOUBLE, MPI_SUM, reciever, MPI_COMM_WORLD);
		} else {
			double null = 0.0;
			double tmpNull = 0.0;
			MPI_Reduce(&null, &tmpNull, 1, MPI_DOUBLE, MPI_SUM, reciever, MPI_COMM_WORLD);
		}

		cntStep--;
		reciever = cntStep%sizeGroup;
	}

	// MPI_Comm_free(&tmp);
	// MPI_Group_free(&reduceGroup);
	// MPI_Group_free(&world_group);
}

void readMatrix(std::vector<Column> &columns, std::vector<double> &b, int N, std::string fileName) {
	ifstream file(fileName.data());

	for(int i = 0; i < N; ++i)
		file >> b[i];

	for (int i = 0; i < N; ++i)
		for (int j = 0; j < N; ++j)
			file >> columns[j].mas[i];

	file.close();
}

void printResults(const std::vector<std::pair<int, double> > &res, int rankInGroup) {
	if(rankInGroup == 0) {
		for (int i = 0; i < res.size(); ++i)
			cout << "x" << res[i].first << " = " << res[i].second << endl;
	}
}

double residual(std::vector<Column> &columns, std::vector<double> &b, std::vector<std::pair<int, double> > &res, int N, int rankInGroup, int sizeGroup) {
	map<int, double> resMap;
	for (int i = 0; i < res.size(); ++i)
		resMap[res[i].first] = res[i].second;

	if (rankInGroup == 0) {
		double* tmp = (double *)malloc(N * sizeof(double));
		memset(tmp, 0, N * sizeof(double));
		for(int i = 0; i < N; i++) {
			for (int j = 0; j < columns.size(); j++)
				tmp[i] += columns[j].mas[i] * resMap[columns[j].j];
		}

		double* tmp2 = (double *)malloc(N * sizeof(double));
		MPI_Reduce(tmp, tmp2, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		double res = 0.0;
		for(int i = 0; i < N; i++) 
			res += (tmp2[i] - b[i]) * (tmp2[i] - b[i]);
		res = sqrt(res); 
		free(tmp);
		free(tmp2);
		return res;
	} else {
		double* tmp = (double *)malloc(N * sizeof(double));
		memset(tmp, 0, N * sizeof(double));
		for(int i = 0; i < N; i++) {
			for (int j = 0; j < columns.size(); j++)
				tmp[i] += columns[j].mas[i] * resMap[columns[j].j];
		}

		MPI_Reduce(tmp, NULL, N, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		free(tmp);
	}

	return 0;
}

void allocMemmory(std::vector<Column> &columns, int N, int rankInGroup, int sizeGroup) {
	int size = N / sizeGroup + (rankInGroup < (N % sizeGroup) ? 1 : 0);
	columns.resize(size);
	int k = rankInGroup;
	for(int i = 0; i < columns.size(); i++) {
		columns[i].mas = (double*)malloc(N*sizeof(double));
		memset(columns[i].mas, 0, N*sizeof(double));
		columns[i].j = k;
		k += sizeGroup; 
	}
}

void freeMemmory(std::vector<Column> &columns) {
	for(int i = 0; i < columns.size(); i++)
		free(columns[i].mas);
}

int checkMatrix(std::vector<Column> &columns, int N, int rankInGroup, int sizeGroup) {
	int error = 0;
	for (int i = 0; i < columns.size(); ++i)
		if (fabs(columns[i].mas[columns[i].j]) < EPS) {
			error = 1;
			break;
		}

	int tmp;
	MPI_Allreduce(&error, &tmp, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	return tmp;
}
#endif