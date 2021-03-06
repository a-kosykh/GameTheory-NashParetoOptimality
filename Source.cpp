#include <iostream>
#include "NashParetoSolver.h"



BimatrixGame familyDispute = {
	{{4, 1}, {0, 0}},
	{{0, 0}, {1, 4}}
};

BimatrixGame prisonersDilemma = {
	{{-5.0, -5.0}, {0.0, -10.0}},
	{{-10.0, 0.0}, {-1.0, -1.0}}
};

BimatrixGame test = {
	{{0.0, 1.0}, {11.0, 4.0}},
	{{7.0, 8.0}, {6.0, 3.0}}
};
BimatrixGame test2 = {
	{{5.0, 8.0}, {7.0, 4.0}},
	{{11.0, 7.0}, {6.0, 9.0}}
};

int main() {
	std::cout << "Kosykh A.M" << std::endl;
	std::cout << "Variant 5" << std::endl;
	std::cout << "Lab 3" << std::endl;
	std::cout << "Random 10x10 game:" << std::endl;
	BimatrixGame randomGame = createRandomGame(10);
	printGame(randomGame);
	NashParetoSolver randomGameSolver(createRandomGame(10));
	randomGameSolver.iSolve();
	randomGameSolver.iPrint();
	std::cout << std::endl << std::endl;

	std::cout << "Family Dispute Game:" << std::endl;
	printGame(familyDispute);
	NashParetoSolver familyDisputeGS(familyDispute);
	familyDisputeGS.iSolve();
	familyDisputeGS.iPrint();
	std::cout << std::endl << std::endl;

	std::cout << "CrossRoad Game:" << std::endl;
	double e1 = 0.5;
	double e2 = 1.2;
	BimatrixGame crossRoad = {
		{{1, 1}, {1 - e1, 2}},
		{{2, 1 - e2}, {0, 0}}
	};
	printGame(crossRoad);
	NashParetoSolver crossRoadGS(crossRoad);
	crossRoadGS.iSolve();
	crossRoadGS.iPrint();
	std::cout << std::endl << std::endl;

	std::cout << "Prisoner's Dilemma Game:" << std::endl;
	printGame(prisonersDilemma);
	NashParetoSolver prisonersDilemmaGS(prisonersDilemma);
	prisonersDilemmaGS.iSolve();
	prisonersDilemmaGS.iPrint();
	std::cout << std::endl << std::endl;

	std::cout << "theorem test Game:" << std::endl;
	printGame(test2);
	NashParetoSolver testGS(test2);
	testGS.iSolve();
	testGS.iPrint();
	testGS.iSolveTheorem();
	std::cout << std::endl << std::endl;

	return 0;
}