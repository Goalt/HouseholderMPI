#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <vector>
#include <mpi.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <math.h>
#include <string>
#include <fstream>
#include <map>

#define EPS 0.0000000000000000001

struct Column {
	double* mas;
	int j;
};

void gatherResults(std::vector<std::pair<int, double> > &res, int N, int rankInGroup, int sizeGroup);

void allocMemmory(std::vector<Column> &columns, int N, int rankInGroup, int sizeGroup);

void freeMemmory(std::vector<Column> &columns);

void fillMas(std::vector<Column> &columns, int N);

void fillB(std::vector<double> &b, std::vector<Column> &columns, int N);

void printMatrix(std::vector<Column> columns, int N, int rankInGroup, int sizeGroup);

void makeTriangle(std::vector<Column> &columns, std::vector<double> &b, int N, int rankInGroup, int sizeGroup);

void backGaus(std::vector<Column> &columns, std::vector<double> &b, std::vector<std::pair<int, double> > &res, int N, int rankInGroup, int sizeGroup);

void readMatrix(std::vector<Column> &columns, std::vector<double> &b, int N, std::string fileName);

void printResults(const std::vector<std::pair<int, double> > &res, int rankInGroup);

double residual(std::vector<Column> &columns, std::vector<double> &b, std::vector<std::pair<int, double> > &res, int N, int rankInGroup, int sizeGroup);

int checkMatrix(std::vector<Column> &columns, int N, int rankInGroup, int sizeGroup);
#endif