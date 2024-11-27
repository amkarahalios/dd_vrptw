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
    DDStats() : lpIterations(0), numLagIterations(0), numLagIterationsWithResets(0), numSSPIterations(0), numSeparations(0), millisecondsCompiling(0), millisecondsSolvingSSP(0), millisecondsSolvingLP(0), millisecondsRepairingLAG(0), millisecondsFindingCuts(0), millisecondsTryingAlpha(0), millisecondsDecompose(0), millisecondsYellow(0), millisecondsUpdateDual(0), millisecondsFix(0), lowerBound(0), upperBound(INF), numCuts(0)
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

      std::cout << "STATS1 - lpIterations[" << lpIterations << "] lagIterations[";
      std::cout << numLagIterations << "] sspIterations[" << numSSPIterations;
      std::cout << "] numSeparations[" << numSeparations << "] numCuts[" << numCuts;
      std::cout << "] compileTime[" << int(millisecondsCompiling / 1000);
      std::cout << "] sspSolveTime[" << int(millisecondsSolvingSSP / 1000);
      std::cout << "] lpSolveTime[" << int(millisecondsSolvingLP / 1000);
      std::cout << "] repairTime[" << int(millisecondsRepairingLAG / 1000);
      std::cout << "] cutsTime[" << int(millisecondsFindingCuts / 1000);
      std::cout << "] alphaTime[" << int(millisecondsTryingAlpha / 1000);
      std::cout << "] decompTime[" << int(millisecondsDecompose / 1000);
      std::cout << "] yellowTime[" << int(millisecondsYellow / 1000);
      std::cout << "] updateTime[" << int(millisecondsUpdateDual / 1000);
      std::cout << "] fixTime[" << int(millisecondsFix/ 1000);
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

    float millisecondsRepairingLAG;
    float millisecondsFindingCuts;
    float millisecondsTryingAlpha;
    float millisecondsDecompose;
    float millisecondsYellow;
    float millisecondsUpdateDual;
    float millisecondsFix;
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

    bool primalHeuristicMIP(std::vector<std::vector<int>>& routesByLocationPrimalHeuristic);

    void initializeDual(Dual& dual);
    bool solveLagrangeanRelaxation(Dual& duals, SGDAlgorithm& sgdAlgo);

    void addColumn();
    void initializeColumns();
    bool solvePricingProblem(std::vector<double>& lambda);
    bool solveLPCG(Dual& duals);

    void updateMultipliers(Dual& duals, double lagrangeanLowerBound, int iteration);
    void updateMultipliersVolumeAlgorithm(Dual& duals, Primal& primal, int iteration);
    void printMultipliers(Dual& duals);
    void repairMultipliers(Dual& repairedDual, LPSolveType lpSolveType);
    void resizeMultipliers(const Dual& dual1, Dual& dual2);
    void resizeMultipliersAndCopy(const Dual& dual1, Dual& dual2);
    void resizeMultipliers(const Dual& dual1, std::vector<Dual>& dual2);
 
    void constructNextPrimal(double alpha, const std::vector<std::vector<int>>& decomposedRoutes, const std::vector<std::vector<int>>& decomposedRouteArcs, Primal& primal);
    void getGradient(const Primal& primal, const Dual& dual, std::vector<double>& gradient);
    double calculateTwoNorm(const std::vector<double>& gradient);
    double calculateDotProduct(const std::vector<double>& vector1, const std::vector<double>& vector2);

    VRPTW vrptw;
    VRPTWDecisionDiagram routeDD;
    VRPTWDDParameters params;
    std::vector<std::vector<int>> bestSolution;
    double muPercentImproved;

    // for arc fixing
    std::vector<Dual> currDualsArcFixing;
    std::vector<Dual> bestDualsArcFixingLag;
    double bestDualsArcFixingLagPercent;
    Dual bestDualArcFixingLP;
    double bestDualArcFixingLPPercent;

    // for subgradient descent / volume algorithm
    Dual bestDual;
    double bestDualValue;
    std::vector<double> previousGradient;

    // for volume algorithm
    double stepSizeMultiplier;
    int stepSizeMultiplierIteration;
    int stepSizeMultiplierIterationCutoff;
    double minStepSizeMultiplier;
    double alphaLowerBound;
    int alphaLowerBoundIteration;
    double targetLowerBound;
    double percentGapToStartCuts;
    int numCapCutUpdateGroups;
    double targetBoundIncrease;
    int numIterationsRestart;
    int numIterationsDeactivations;

    // static params
    int infeasibleRoutesBatchSize;
    double deactivateCutValueThreshold;
    int deactivateCutIterThreshold;
    double lagOptimalityGapToStartRepairing;
    double percentFixedToChangeToCPLEX;
    int numArcsToChangeToCPLEX;
    int numLagCuts;
    int kappaIterations;
    bool shouldRestart;
    bool useScaling;

    // for primal
    Primal primal;
    std::vector<std::vector<int>> primalRoutes;
    std::vector<double> primalRouteCosts;
    std::map<int,std::vector<int>> primalRouteLocationIndices;

    // for cuts
    int Dim = 100;
    double epsForIntegrality = 0.0001;
    CnstrMgrPointer MyCutsCMP;
    CnstrMgrPointer MyOldCutsCMP;

    // used to deactivate cuts
    std::vector<int> capCutTooSmallCounters;
    std::vector<int> srcCutTooSmallCounters;

    DDStats stats;
};

#endif
