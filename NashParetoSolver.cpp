#include "NashParetoSolver.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

bool NashParetoSolver::FindNashEquilibrium()
{
	for (size_t i = 0; i < m_gameMatrix.size(); ++i) {
		for (size_t j = 0; j < m_gameMatrix[i].size(); ++j) {		
			bool isEquilibrium = true;
			for (size_t k = 0; k < m_gameMatrix.size(); ++k) {
				if (i == k) continue;
				if (m_gameMatrix[i][j].first < m_gameMatrix[k][j].first) {
					isEquilibrium = false;
					break;
				}
			}
			if (isEquilibrium) {
				for (size_t k = 0; k < m_gameMatrix[i].size(); ++k) {
					if (j == k) continue;
					if (m_gameMatrix[i][j].second < m_gameMatrix[i][k].second) {
						isEquilibrium = false;
						break;
					}
				}
			}
			if (isEquilibrium) {
				m_nashEquilibrium.push_back(m_gameMatrix[i][j]);
			}
		}
	}
	return true;
}

bool NashParetoSolver::FindParetoOptimality()
{
	for (size_t i = 0; i < m_gameMatrix.size(); ++i) {
		for (size_t j = 0; j < m_gameMatrix[i].size(); ++j) {

			bool isOptimal = true;
			for (size_t k = 0; k < m_gameMatrix.size(); ++k) {
				for (size_t l = 0; l < m_gameMatrix[i].size(); ++l) {
					if (i == k && j == l) continue;
					if (m_gameMatrix[k][l].first > m_gameMatrix[i][j].first) {
						if (m_gameMatrix[k][l].second > m_gameMatrix[i][j].second) {
							isOptimal = false;
							break;
						}
					}
					if (isOptimal) {
						if (m_gameMatrix[k][l].second > m_gameMatrix[i][j].second) {
							if (m_gameMatrix[k][l].first > m_gameMatrix[i][j].first) {
								isOptimal = false;
								break;
							}
						}
					}
				}
			}
			if (isOptimal)
				m_paretoOptimality.push_back(m_gameMatrix[i][j]);
		}
	}
	return true;
}

bool NashParetoSolver::FindEquilibriumSituation()
{
	MatrixGame v2 = { {GetGameValue(m_secondPlayerMatrix)} };
	MatrixGame v1 = { {GetGameValue(m_firstPlayerMatrix)} };
	MatrixGame ident = { {1, 1} };

	m_x = MultMatrix(MultMatrix(v2, ident), GetInverse(m_secondPlayerMatrix));
	m_y = GetTranspose(MultMatrix(MultMatrix(GetInverse(m_firstPlayerMatrix), GetTranspose(ident)), v1));

	return true;
}

void NashParetoSolver::PrintNashEquilibrium()
{
	std::cout << "Nash: {";
	for (size_t i = 0; i < m_nashEquilibrium.size(); ++i) {
		std::cout << "(" << m_nashEquilibrium.at(i).first << "/"
			<< m_nashEquilibrium.at(i).second << ")";
		if (i != m_nashEquilibrium.size() - 1) std::cout << ", ";
	}
	std::cout << "}";
}

void NashParetoSolver::PrintParetoOptimality()
{
	std::cout << "Pareto: {";
	for (size_t i = 0; i < m_paretoOptimality.size(); ++i) {
		std::cout << "(" << m_paretoOptimality.at(i).first << "/"
			<< m_paretoOptimality.at(i).second << ")";
		if (i != m_paretoOptimality.size() - 1) std::cout << ", ";
	}
	std::cout << "}";
}

void NashParetoSolver::PrintNashParetoIntersection()
{
	std::vector<std::pair<double, double>> intersection;
	std::sort(m_paretoOptimality.begin(), m_paretoOptimality.end());
	std::sort(m_nashEquilibrium.begin(), m_nashEquilibrium.end());
	std::set_intersection(
		m_paretoOptimality.begin(), m_paretoOptimality.end(),
		m_nashEquilibrium.begin(), m_nashEquilibrium.end(),
		std::back_inserter(intersection));
	
	std::cout << "Pareto-Nash Intersection: {";
	for (size_t i = 0; i < intersection.size(); ++i) {
		std::cout << "(" << intersection.at(i).first << "/"
			<< intersection.at(i).second << ")";
		if (i != intersection.size() - 1) std::cout << ", ";
	}
	std::cout << "}";
}

double NashParetoSolver::GetGameValue(MatrixGame matrix)
{
	MatrixGame ident = { {1, 1} };
	MatrixGame identTrans = GetTranspose(ident);

	MatrixGame ansMatrix = MultMatrix(MultMatrix(ident, GetInverse(matrix)), identTrans);

	return (1 / ansMatrix[0][0]);
}

MatrixGame NashParetoSolver::GetInverse(const MatrixGame vect)
{
	if (GetDeterminant(vect) == 0) {
		throw std::runtime_error("Determinant is 0");
	}

	double d = 1.0 / GetDeterminant(vect);
	MatrixGame solution(vect.size(), std::vector<double>(vect.size()));

	for (size_t i = 0; i < vect.size(); i++) {
		for (size_t j = 0; j < vect.size(); j++) {
			solution[i][j] = vect[i][j];
		}
	}

	solution = GetTranspose(GetCofactor(solution));

	for (size_t i = 0; i < vect.size(); i++) {
		for (size_t j = 0; j < vect.size(); j++) {
			solution[i][j] *= d;
		}
	}

	return solution;
}

MatrixGame NashParetoSolver::GetCofactor(const MatrixGame vect)
{
	if (vect.size() != vect[0].size()) {
		throw std::runtime_error("Matrix is not quadratic");
	}

	MatrixGame solution(vect.size(), std::vector<double>(vect.size()));
	MatrixGame subVect(vect.size() - 1, std::vector<double>(vect.size() - 1));

	for (std::size_t i = 0; i < vect.size(); i++) {
		for (std::size_t j = 0; j < vect[0].size(); j++) {

			int p = 0;
			for (size_t x = 0; x < vect.size(); x++) {
				if (x == i) {
					continue;
				}
				int q = 0;

				for (size_t y = 0; y < vect.size(); y++) {
					if (y == j) {
						continue;
					}

					subVect[p][q] = vect[x][y];
					q++;
				}
				p++;
			}
			solution[i][j] = pow(-1, i + j) * GetDeterminant(subVect);
		}
	}
	return solution;
}

MatrixGame NashParetoSolver::GetTranspose(const MatrixGame matrix1)
{
	MatrixGame solution(matrix1[0].size(), std::vector<double>(matrix1.size()));

	for (size_t i = 0; i < matrix1.size(); i++) {
		for (size_t j = 0; j < matrix1[0].size(); j++) {
			solution[j][i] = matrix1[i][j];
		}
	}
	return solution;
}

MatrixGame NashParetoSolver::MultMatrix(const MatrixGame lhs, const MatrixGame rhs)
{
	const int lhsRowsCount = lhs.size();
	const int lhsColsCount = lhs[0].size();
	const int rhsColsCount = rhs[0].size();

	MatrixGame solution(lhsRowsCount, std::vector<double>(rhsColsCount, 0));
	for (auto j = 0; j < rhsColsCount; ++j)
	{
		for (auto k = 0; k < lhsColsCount; ++k)
		{
			for (auto i = 0; i < lhsRowsCount; ++i)
			{
				solution[i][j] += lhs[i][k] * rhs[k][j];
			}
		}
	}
	return solution;
}

double NashParetoSolver::GetDeterminant(const MatrixGame vect)
{
	if (vect.size() != vect[0].size()) {
		throw std::runtime_error("Matrix is not quadratic");
	}
	int dimension = vect.size();

	if (dimension == 0) {
		return 1;
	}

	if (dimension == 1) {
		return vect[0][0];
	}

	if (dimension == 2) {
		return vect[0][0] * vect[1][1] - vect[0][1] * vect[1][0];
	}

	double result = 0;
	int sign = 1;
	for (int i = 0; i < dimension; i++) {

		MatrixGame subVect(dimension - 1, std::vector<double>(dimension - 1));
		for (int m = 1; m < dimension; m++) {
			int z = 0;
			for (int n = 0; n < dimension; n++) {
				if (n != i) {
					subVect[m - 1][z] = vect[m][n];
					z++;
				}
			}
		}

		result = result + sign * vect[0][i] * GetDeterminant(subVect);
		sign = -sign;
	}

	return result;
}

NashParetoSolver::NashParetoSolver() {}

NashParetoSolver::NashParetoSolver(BimatrixGame gm)
{
	m_gameMatrix = gm;
	for (const auto &row : m_gameMatrix) {
		std::vector<double> firstPlayerRow;
		std::vector<double> secondPlayerRow;
		for (const auto &element : row) {
			firstPlayerRow.push_back(element.first);
			secondPlayerRow.push_back(element.second);
		}
		m_firstPlayerMatrix.push_back(firstPlayerRow);
		m_secondPlayerMatrix.push_back(secondPlayerRow);
	}
}

void NashParetoSolver::iSolve()
{
	FindNashEquilibrium();
	FindParetoOptimality();
}

void NashParetoSolver::iSolveTheorem() 
{
	FindEquilibriumSituation();
	std::cout << std::endl;
	iPrintXY();
}

void NashParetoSolver::iPrint()
{
	PrintNashEquilibrium();
	std::cout << std::endl;
	PrintParetoOptimality();
	std::cout << std::endl;
	PrintNashParetoIntersection();
}

void NashParetoSolver::iPrintXY()
{
	std::cout << "x = [";
	for (auto i = 0; i < m_x.at(0).size(); ++i) {
		std::cout << std::setprecision(3) << m_x.at(0).at(i);
		if (i != m_x.at(0).size() - 1)
			std::cout << ", ";
	}
	std::cout << "]" << std::endl;

	std::cout << "y = [";
	for (auto i = 0; i < m_y.at(0).size(); ++i) {
		std::cout << std::setprecision(3) << m_y.at(0).at(i);
		if (i != m_y.at(0).size() - 1)
			std::cout << ", ";
	}
	std::cout << "]" << std::endl;
}
