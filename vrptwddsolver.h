#ifndef VRPTWDDSOLVER
#define VRPTWDDSOLVER

#include <vector>
#include <chrono>

#include "vrptwdecisiondiagram.h"
#include "vrptw.h"

#include "cnstrmgr.h"
#include "capsep.h"
#include "combsep.h"

struct DDStats
{
  public:
    DDStats() : lpIterations(0), numLagIterations(0), numSSPIterations(0), numSeparations(0), millisecondsCompiling(0), millisecondsSolvingSSP(0), millisecondsSolvingLP(0), lowerBound(0), upperBound(INF), numCuts(0)
    {
      startTime = std::chrono::high_resolution_clock::now();
    };
    
    void print(int ddSize) const
    {
      std::cout << "STATS - lpIterations[" << lpIterations << "] lagIterations[";
      std::cout << numLagIterations << "] sspIterations[" << numSSPIterations;
      std::cout << "] numSeparations[" << numSeparations << "] compileTime[" << int(millisecondsCompiling / 1000);
      std::cout << "] sspSolveTime[" << int(millisecondsSolvingSSP / 1000);
      std::cout << "] lpSolveTime[" << int(millisecondsSolvingLP / 1000);
      std::cout << "] lb[" << lowerBound << "] ub[" << upperBound << "]";
      std::cout << " size: [" << ddSize << "]";
      std::cout << " time: [" << getNumSeconds() << "]" << std::endl;
    }

    int getNumSeconds() const
    {
      auto currTime = std::chrono::high_resolution_clock::now();
      return int(std::chrono::duration_cast<std::chrono::milliseconds>(currTime - startTime).count() / 1000);
    }

    int lpIterations;
    int numLagIterations;
    int numSSPIterations;
    int numSeparations;
    float millisecondsCompiling;
    float millisecondsSolvingSSP;
    float millisecondsSolvingLP;
    std::chrono::high_resolution_clock::time_point startTime;
    double lowerBound;
    double upperBound;
    int numCuts;
};

class VRPTWDDSolver
{
  public:
    VRPTWDDSolver(VRPTW vrptw, LPSolveType lpSolveType, InitialStateSpace initialStateSpace, int s, int maxS, bool useCuts, int timeoutSeconds);

    bool solve();
    DDStats getStats() { return stats; }

    // public for testing purpose only
    void addRCCs(const std::vector<int>& edgeTail, const std::vector<int>& edgeHead, const std::vector<double>& edgeFlow, bool& cutAdded);
    void addCombs(std::vector<int>& edgeTail, std::vector<int>& edgeHead, std::vector<double>& edgeFlow, bool& cutAdded);
    void convertArcIndicesForVRPTWSep(const std::vector<double>& routeFlows,
                                     const std::vector<std::set<int>>& decomposedRoutes,
                                     std::vector<int>& edgeTail,
                                     std::vector<int>& edgeHead,
                                     std::vector<double>& edgeFlow);


  private:
    void addSRCCuts(std::vector<double>& srcDuals);

    bool solveLP(std::vector<double>& lambda, std::vector<double>& capDuals, std::vector<double>& combDuals, std::vector<double>& srcDuals);
    bool solveIP(std::vector<double>& lambda, std::vector<double>& capDuals, std::vector<double>& combDuals, std::vector<double>& srcDuals);
    void addColumn();
    void initializeColumns();
    bool solvePricingProblem(std::vector<double>& lambda);
    bool solveLPCG(std::vector<double>& lambda, std::vector<double>& capDuals, std::vector<double>& combDuals, std::vector<double>& srcDuals);
    bool solveLagrangeanRelaxation(std::vector<double>& lambda, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& muSRC);

    void updateMultipliers(std::vector<double>& lambda, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals, std::vector<double>& stepSizes, const std::set<int>& solutionArcs, double lagrangeanLowerBound, int iteration);

    VRPTW vrptw;
    VRPTWDecisionDiagram routeDD;
    int maxS;
    bool useCuts;
    int timeoutSeconds;
    LPSolveType lpSolveType;
    std::vector<std::vector<int>> bestSolution;

    std::vector<double> bestLambdaArcFixing;
    double bestLambdaPercentFixed = 0.0;
    double bestLambdaLowerBound = 0.0;

    int Dim = 100;
    double epsForIntegrality = 0.0001;
    CnstrMgrPointer MyCutsCMP;
    CnstrMgrPointer MyOldCutsCMP;

    DDStats stats;
};

#endif
