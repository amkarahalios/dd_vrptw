#ifndef VRPTWDDSOLVER
#define VRPTWDDSOLVER

#include <vector>
#include <deque>
#include <chrono>

#include "vrptwdecisiondiagram.h"
#include "vrptw.h"

#include "cnstrmgr.h"
#include "capsep.h"
#include "combsep.h"

struct DDStats
{
  public:
    DDStats() : lpIterations(0), numLagIterations(0), numLagIterationsWithResets(0), numSSPIterations(0), numSeparations(0), millisecondsCompiling(0), millisecondsSolvingSSP(0), millisecondsSolvingLP(0), millisecondsRepairingLAG(0), millisecondsFindingCuts(0), millisecondsTryingAlpha(0), millisecondsDecompose(0), millisecondsYellow(0), millisecondsUpdateDual(0), millisecondsFix(0), lowerBound(0), numHeuristicIPs(0), numHeuristicLNSs(0), numPrimalLNSRepairs(0), upperBound(INF), numCuts(0)
    {
      startTime = std::chrono::high_resolution_clock::now();
    };
    
    void print(int ddNumArcs, int ddNumFixedArcs) const;
    int getNumSeconds() const;

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
    int numHeuristicIPs;
    int numHeuristicLNSs;
    int numPrimalLNSRepairs;
};

class VRPTWDDSolver
{
  public:
    VRPTWDDSolver(VRPTW vrptw, VRPTWDDParameters);

    bool solve(bool solveIP);
    DDStats getStats() { return stats; }

  private:
    bool solveLP(Dual& duals);
    bool solveIP(Dual& duals);
 
    void initializeDual(Dual& dual);
    bool solveLagrangeanRelaxation(Dual& duals, SGDAlgorithm& sgdAlgo);
 
    void constructNextPrimal(double alpha, const std::vector<std::vector<int>>& decomposedRoutes, const std::vector<std::vector<int>>& decomposedRouteArcs, Primal& primal);
    void getGradient(const Primal& primal, const Dual& dual, std::vector<double>& gradient);
    double calculateTwoNorm(const std::vector<double>& gradient);
    double calculateDotProduct(const std::vector<double>& vector1, const std::vector<double>& vector2);
 
    void updateMultipliers(Dual& duals, double lagrangeanLowerBound, int iteration);
    void updateMultipliersVolumeAlgorithm(Dual& duals, Primal& primal, int iteration);
    void printMultipliers(Dual& duals);
    void repairMultipliers(Dual& repairedDual, LPSolveType lpSolveType);
    void resizeMultipliers(const Dual& dual1, Dual& dual2);
    void resizeMultipliers(const Dual& dual1, std::vector<Dual>& dual2);
    void resizeMultipliers(const Dual& dual1, std::deque<Dual>& dual2);

    bool largeNeighborhoodSearch(const std::vector<std::vector<int>>& feasibleSolution, const Dual& dual, int timeoutSeconds);
    bool primalHeuristicMIP(FlowType flowType, std::vector<std::vector<int>>& routesByLocationPrimalHeuristic, std::set<int>& returnRouteIndices);
    void addRouteToPrimalRoutes(std::vector<int> route);
    bool intraRouteSwaps(std::vector<int>& route, double& routeCost);

    void destroyBySingleRoute(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& destroyedElements, int numElementsToDestroy);
    void destroyByShaw(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& destroyedElements, int numElementsToDestroy);
    double computeShawRelatedness(const std::vector<int>& route1, int location1, const std::vector<int>& route2, int location2);
    void destroyByWorst(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& destroyedElements, int numElementsToDestroy);
    void destroyRandomly(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& destroyedElements, int numElementsToDestroy);
    void chooseRandomLocationsFromRoutes(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& randomElementsChosen);
    void removeElements(std::set<int>& removedElements, const std::vector<std::vector<int>>& feasibleSolution, int numElementsDestroy, int extraElementsDestroy, int sequenceIndex);

    bool addCutsUsingCurrentPrimal(Dual& dual, const std::vector<std::vector<int>>& decomposedRoutes);
    void separateSequencesAndTruncate(const std::vector<std::vector<int>>& inputSequences, std::vector<std::vector<int>>& outputSequences, std::vector<std::vector<int>>& outputSequencesArcs);

    void addSRCCuts(std::vector<double>& srcDuals, const std::vector<double>& violations);
    void addRCCs(const std::vector<int>& edgeTail, const std::vector<int>& edgeHead, const std::vector<double>& edgeFlow, std::vector<int>& rccArcs, int maxNumCuts, bool& cutAdded, Dual& dual);
    void addCombs(std::vector<int>& edgeTail, std::vector<int>& edgeHead, std::vector<double>& edgeFlow, bool& cutAdded);
    void convertArcIndicesForVRPTWSep(const Primal& primal,
                                      std::vector<int>& edgeTail,
                                      std::vector<int>& edgeHead,
                                      std::vector<double>& edgeFlow,
                                      std::vector<int>& rccArcs,
                                      std::vector<double>& rccArcFlows);

    VRPTW vrptw;
    VRPTWDecisionDiagram routeDD;
    VRPTWDDParameters params;
    std::vector<std::vector<int>> bestSolution;

    // for arc fixing
    std::deque<Dual> bestDualsArcFixingLag;
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
    double targetBoundIncrease;

    // static params
    int infeasibleRoutesBatchSize;
    double percentFixedToChangeToCPLEX;
    int numArcsToChangeToCPLEX;
    int numLagCuts;
    int kappaIterations;
    double muPercentImproved;
    double percentGapToStartCuts;

    // for primal
    Primal primal;
    std::vector<std::vector<int>> primalRoutes;
    std::vector<double> primalRouteCosts;
    std::map<int,std::vector<int>> primalRouteLocationIndices;
    std::set<int> mipPoolSolutionIndices;

    // for cuts
    int Dim = 100;
    double epsForIntegrality = 0.0001;
    CnstrMgrPointer MyCutsCMP;
    CnstrMgrPointer MyOldCutsCMP;

    DDStats stats;
};

#endif
