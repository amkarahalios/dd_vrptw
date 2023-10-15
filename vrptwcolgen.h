#ifndef VRPTWCG
#define VRPTWCG

#include <vector>
#include <chrono>

#include "vrptwdecisiondiagram.h"
#include "vrptw.h"
#include "cvrpsep/cnstrmgr.h"

struct ColGenStats
{
  public:
    ColGenStats() : numIterations(0), millisecondsCompilingDD(0), millisecondsSolvingPricing(0), millisecondsSolvingMaster(0), lowerBound (0.0), solutionValue(0.0) {};

    void print() const
    {
      std::cout << "STATS - iterations[" << numIterations << "] ddTime[";
      std::cout << int(millisecondsCompilingDD / 1000) << "] pricingTime[";
      std::cout << int(millisecondsSolvingPricing / 1000) << "] masterTime[" << int(millisecondsSolvingMaster/1000);
      std::cout << "] lb[" << lowerBound;
      std::cout << "] sol[" << solutionValue << "]";
      std::cout << " time: [" << getNumSeconds() << "]" << std::endl;
    }

    int getNumSeconds() const
    {
      return int((millisecondsCompilingDD+millisecondsSolvingPricing+millisecondsSolvingMaster) / 1000);
    }

    int numIterations;
    float millisecondsCompilingDD;
    float millisecondsSolvingPricing;
    float millisecondsSolvingMaster;
    double lowerBound;
    double solutionValue;
};

class VRPTWColGen
{
  public:
    VRPTWColGen(VRPTW vrptw, VRPTWDDParameters params, PricingProblemType pricingProblemType, InitialStateSpace initialStateSpace, int s);

    bool solve();
    ColGenStats getStats() { return stats; }

  private:
    void addColumn(std::vector<int> route);
    void routeToColumn(const std::vector<int>& route, std::vector<int>& column);
    void initializeColumns();
    bool setupAndSolveRMP();

    bool solvePricingProblem();
    bool solvePricingProblemDD();
    bool solvePricingProblemFukasawa();

    bool separateFractionalSolution();

    VRPTW vrptw;
    PricingProblemType pricingProblemType;
    int s;
    VRPTWDecisionDiagram routeDD;

    std::vector<std::vector<int>> columns;
    std::vector<double> columnCosts;

    double bestLpDistance;
    std::vector<double> lpSolution;
    std::vector<double> setCoverDualVariables;
    double singlePathDual;
    double numTrucksDualVariable;

    //CnstrMgrPointer allCutsCMP;
    //std::vector<std::vector<int>> lhsCuts;
    //std::vector<int> rhsCuts;

    ColGenStats stats;
};

#endif
