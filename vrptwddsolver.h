#ifndef VRPTWDDSOLVER
#define VRPTWDDSOLVER

#include <vector>
#include <chrono>

#include "vrptwdecisiondiagram.h"
#include "vrptw.h"

#include "cnstrmgr.h"
#include "capsep.h"
#include "combsep.h"

enum PhaseType
{
  INITIAL_DUAL = 0,
  SEPARATION = 1,
  ROBUST_PRIMAL = 2,
  ROBUST_CUT_DUALS = 3,
  NONROBUST_PRIMAL = 4,
  NONROBUST_CUT_DUALS = 5
};

struct DDStats
{
  public:
    DDStats() : lpIterations(0), numLagIterations(0), numLagIterationsWithResets(0), numSSPIterations(0), numSeparations(0), millisecondsCompiling(0), millisecondsSolvingSSP(0), millisecondsSolvingLP(0), lowerBound(0), upperBound(INF), numCuts(0)
    {
      startTime = std::chrono::high_resolution_clock::now();
    };
    
    void print(int ddNumArcs, int ddNumFixedArcs) const
    {
      std::cout << "STATS - lpIterations[" << lpIterations << "] lagIterations[";
      std::cout << numLagIterations << "] sspIterations[" << numSSPIterations;
      std::cout << "] numSeparations[" << numSeparations << "] numCuts[" << numCuts;
      std::cout << "] compileTime[" << int(millisecondsCompiling / 1000);
      std::cout << "] sspSolveTime[" << int(millisecondsSolvingSSP / 1000);
      std::cout << "] lpSolveTime[" << int(millisecondsSolvingLP / 1000);
      std::cout << "] lb[" << lowerBound << "] ub[" << upperBound << "]";
      std::cout << " numArcs: [" << ddNumArcs << "]";
      std::cout << " numFixed: [" << ddNumFixedArcs << "]";
      std::cout << " time: [" << getNumSeconds() << "]" << std::endl;
    }

    int getNumSeconds() const
    {
      auto currTime = std::chrono::high_resolution_clock::now();
      return int(std::chrono::duration_cast<std::chrono::milliseconds>(currTime - startTime).count() / 1000);
    }

    int lpIterations;
    int numLagIterations;
    int numLagIterationsWithResets;
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
    VRPTWDDSolver(VRPTW vrptw, VRPTWDDParameters);

    bool solve(bool solveIP);
    DDStats getStats() { return stats; }

    // public for testing purpose only
    void addRCCs(const std::vector<int>& edgeTail, const std::vector<int>& edgeHead, const std::vector<double>& edgeFlow, std::vector<int>& rccArcs, int maxNumCuts, bool& cutAdded, Dual& dual);
    void addCombs(std::vector<int>& edgeTail, std::vector<int>& edgeHead, std::vector<double>& edgeFlow, bool& cutAdded);
    void convertArcIndicesForVRPTWSep(const Primal& primal,
                                      std::vector<int>& edgeTail,
                                      std::vector<int>& edgeHead,
                                      std::vector<double>& edgeFlow,
                                      std::vector<int>& rccArcs,
                                      std::vector<double>& rccArcFlows);


  private:
    bool addCutsUsingCurrentPrimal(Dual& duals);
    void addSRCCuts(std::vector<double>& srcDuals, const std::vector<double>& violations);

    bool solveLP(Dual& duals);
    bool solveIP(Dual& duals);

    void initializeDual(Dual& dual);
    bool solveLagrangeanRelaxation(Dual& duals);
    bool solveLagrangeanRelaxationVolumeAlgorithm(Dual& duals);

    void addColumn();
    void initializeColumns();
    bool solvePricingProblem(std::vector<double>& lambda);
    bool solveLPCG(Dual& duals);

    void updateMultipliers(Dual& duals, double lagrangeanLowerBound, int iteration);
    void updateMultipliersVolumeAlgorithm(Dual& duals, Primal& primal, double lagrangeanLowerBound, int iteration);
    void printMultipliers(Dual& duals);
    void repairMultipliers(Dual& repairedDual, LPSolveType lpSolveType);
    void resizeMultipliers(const Dual& dual1, Dual& dual2);
 
    void constructNextPrimal(double alpha, const std::vector<std::vector<int>>& decomposedRoutes, const std::vector<std::vector<int>>& decomposedRouteArcs, Primal& primal);
    void getGradient(const Primal& primal, const Dual& dual, std::vector<double>& gradient);
    double calculateTwoNorm(const std::vector<double>& gradient);
    double calculateDotProduct(const std::vector<double>& vector1, const std::vector<double>& vector2);

    VRPTW vrptw;
    VRPTWDecisionDiagram routeDD;
    VRPTWDDParameters params;
    std::vector<std::vector<int>> bestSolution;
    PhaseType phaseType;
    double muPercentImproved;

    // for arc fixing
    Dual bestDualArcFixing;
    double bestDualArcFixingPercent;

    // for subgradient descent / volume algorithm
    Dual bestDual;
    double bestDualValue;
    std::vector<double> previousGradient;

    // for volume algorithm
    double stepSizeMultiplier;
    int stepSizeMultiplierIteration;
    int stepSizeMultiplierIterationCutoff;
    double alphaLowerBound;
    int alphaLowerBoundIteration;
    double alphaLowerBoundCheckValue;
    double targetLowerBound;

    // for primal
    Primal primal;

    // for cuts
    int Dim = 100;
    double epsForIntegrality = 0.0001;
    CnstrMgrPointer MyCutsCMP;
    CnstrMgrPointer MyOldCutsCMP;

    // used to deactivate cuts
    std::vector<int> capCutTooSmallCounters;
    std::vector<int> cliqueCutTooSmallCounters;

    DDStats stats;
};

#endif
