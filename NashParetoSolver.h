#pragma once
#include "BimatrixGame.h"

using Matrix = std::vector<std::vector<double>>;

class NashParetoSolver {
private:
	BimatrixGame	m_gameMatrix;
	Matrix		m_firstPlayerMatrix;
	Matrix		m_secondPlayerMatrix;
	Matrix		m_x;
	Matrix		m_y;
	
	std::vector<std::pair<double, double>> m_nashEquilibrium;
	std::vector<std::pair<double, double>> m_paretoOptimality;

	bool FindNashEquilibrium();
	bool FindParetoOptimality();
	bool FindEquilibriumSituation();

	void PrintNashEquilibrium();
	void PrintParetoOptimality();
	void PrintNashParetoIntersection();

	double GetGameValue(Matrix matrix);
	Matrix GetInverse(const Matrix vect);
	Matrix GetCofactor(const Matrix vect);
	Matrix GetTranspose(const Matrix matrix1);
	Matrix MultMatrix(const Matrix lhs, const Matrix rhs);
	double GetDeterminant(const Matrix vect);

public:
	NashParetoSolver();
	NashParetoSolver(BimatrixGame gm);
	void iSolve();
	void iSolveTheorem();
	void iPrint();
	void iPrintXY();
};
