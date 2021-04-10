#pragma once
#include "BimatrixGame.h"

using MatrixGame = std::vector<std::vector<double>>;

class NashParetoSolver {
private:
	BimatrixGame	m_gameMatrix;
	MatrixGame		m_firstPlayerMatrix;
	MatrixGame		m_secondPlayerMatrix;
	MatrixGame		m_x;
	MatrixGame		m_y;
	
	std::vector<std::pair<double, double>> m_nashEquilibrium;
	std::vector<std::pair<double, double>> m_paretoOptimality;

	bool FindNashEquilibrium();
	bool FindParetoOptimality();
	bool FindEquilibriumSituation();

	void PrintNashEquilibrium();
	void PrintParetoOptimality();
	void PrintNashParetoIntersection();

	double GetGameValue(MatrixGame matrix);
	MatrixGame GetInverse(const MatrixGame vect);
	MatrixGame GetCofactor(const MatrixGame vect);
	MatrixGame GetTranspose(const MatrixGame matrix1);
	MatrixGame MultMatrix(const MatrixGame lhs, const MatrixGame rhs);
	double GetDeterminant(const MatrixGame vect);

public:
	NashParetoSolver();
	NashParetoSolver(BimatrixGame gm);
	void iSolve();
	void iSolveTheorem();
	void iPrint();
	void iPrintXY();
};
