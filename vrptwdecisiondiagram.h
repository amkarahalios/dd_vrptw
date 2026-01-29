#ifndef VRPTWDD
#define VRPTWDD

#include <vector>
#include <deque>
#include <set>
#include <unordered_map>
#include <map>

#include "vrptw.h"

struct Dual
{
  public:
    Dual() {}
    Dual(const Dual& dual)
    {
      lambda = dual.lambda;
      capDuals = dual.capDuals;
      combDuals = dual.combDuals;
      srcDuals = dual.srcDuals;
      fixedPathDual = dual.fixedPathDual;
    }

    std::vector<double> lambda;
    std::vector<double> capDuals;
    std::vector<double> combDuals;
    std::vector<double> srcDuals;
    double fixedPathDual;
};

struct Primal
{
  public:
    Primal() {}
    Primal(Primal& primal)
    {
      xDecompositions = primal.xDecompositions;
      xDecompositionArcs = primal.xDecompositionArcs;
      xDecompositionFlows = primal.xDecompositionFlows;
    }

    std::vector<std::vector<int>> xDecompositions;
    std::vector<std::vector<int>> xDecompositionArcs;
    std::vector<double> xDecompositionFlows;
};

struct VRPTWNodeState
{
  VRPTWNodeState(int _counter, int _load, int _timeWithMultiplier, int _lastVisited, bool _isExact, std::set<int> _visited) : counter(_counter), load(_load), timeWithMultiplier(_timeWithMultiplier), lastVisited(_lastVisited), isExact(_isExact), visited(_visited) {};
  VRPTWNodeState(const VRPTWNodeState& state) : counter(state.counter), load(state.load), timeWithMultiplier(state.timeWithMultiplier), lastVisited(state.lastVisited), isExact(state.isExact), visited(state.visited) {};

  bool operator==(const VRPTWNodeState & rhs) const
  {
    if ((counter == rhs.counter) && (load == rhs.load) && (timeWithMultiplier == rhs.timeWithMultiplier) && (lastVisited == rhs.lastVisited) && (visited == rhs.visited) && (isExact == rhs.isExact))
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  // does not need full ordering, more of equivalence classes in nodeOrdering
  bool operator<(const VRPTWNodeState & rhs) const
  {
    if (counter < rhs.counter)
    {
      return true;
    }
    else if (counter > rhs.counter)
    {
      return false;
    }
    else
    {
      if (load < rhs.load)
      {
        return true;
      }
      else if (load > rhs.load)
      {
        return false;
      }
      else
      {
        if (timeWithMultiplier < rhs.timeWithMultiplier)
        {
          return true;
        }
        else
        {
          return false;
        }
      }
    }
  }

  int counter;
  int load;
  int timeWithMultiplier;
  int lastVisited;
  bool isExact;
  std::set<int> visited;
};

struct hash_state
{
  std::size_t operator()(const VRPTWNodeState& state) const
  {
    std::size_t h1 = std::hash<int>()(state.counter);
    std::size_t h2 = std::hash<int>()(state.load);
    std::size_t h3 = std::hash<int>()(state.timeWithMultiplier);
    std::size_t h4 = std::hash<int>()(state.lastVisited);
    std::size_t h5 = 0;
    if (!state.visited.empty())
    {
      h5 = std::hash<int>()(*state.visited.begin());
    }
    return ((((h1 ^ h2) ^ h3) ^ h4) ^ h5);
  }
};

class VRPTWNode
{
  public:
    VRPTWNode(VRPTWNodeState _state) : state(_state), potential(0), shortestPathDistance(INF) {};

    void print() const;
    void printForDD(int nodeIndex) const;

    VRPTWNodeState state;
    std::vector<int> outArcs;
    std::vector<int> inArcs;

    double potential;
    double shortestPathDistance;
};

struct VRPTWArc
{
  public:
    VRPTWArc(int _fromNodeIndex, int _toNodeIndex, int _location, double _distance) : fromNodeIndex(_fromNodeIndex), toNodeIndex(_toNodeIndex), location(_location), distance(_distance), decompositionFlow(0.0), heuristicFlow(0.0), coeff(0.0), cijPi(0.0), isReverseArc(false) {};
    VRPTWArc(int _fromNodeIndex, int _toNodeIndex, int _location, double _distance, double _coeff, double _cijPi) : fromNodeIndex(_fromNodeIndex), toNodeIndex(_toNodeIndex), location(_location), distance(_distance), decompositionFlow(0.0), heuristicFlow(0.0), coeff(_coeff), cijPi(_cijPi), isReverseArc(true) {};

    bool operator==(const VRPTWArc& rhs) const
    {
      if ((fromNodeIndex == rhs.fromNodeIndex) && (toNodeIndex == rhs.toNodeIndex))
      {
        return true;
      }
      else
      {
        return false;
      }
    }

    void print() const;

    int fromNodeIndex;
    int toNodeIndex;

    int location;
    double distance;

    double decompositionFlow;
    double heuristicFlow;
    double coeff;
    double cijPi;

    bool isReverseArc;
};

class VRPTWDecisionDiagram
{
  public:
    VRPTWDecisionDiagram(const VRPTW& _vrptw, VRPTWDDParameters _params);

    void setTimeStepSize();
    void updateVRPTW(const VRPTW& _vrptw)
    {
      vrptw = _vrptw;
    }
    void updateTimeStateMultiplierByTen();

    void print() const;
    int numNodes() const { return nodes.size(); }
    double evaluateRouteCost(const std::vector<int>& routeByArc);
    double getPercentFixedArcs() const
    {
      return fixedArcs.size() * 100.0 / (arcs.size() - removedArcs.size() - arcReverseArc.size());
    }
    int getNumArcsNotRemovedOrReverse() const
    {
      return arcs.size() - removedArcs.size() - arcReverseArc.size();
    }
    int getNumArcsNotRemovedOrReverseOrFixed() const
    {
      return arcs.size() - removedArcs.size() - arcReverseArc.size() - getNumFixedArcs();
    }
    int getNumFixedArcs() const
    {
      return fixedArcs.size();
    }
    bool isArcAlive(int arcIndex) const
    {
      return (fixedArcs.find(arcIndex) == fixedArcs.end()) && (removedArcs.find(arcIndex) == removedArcs.end());
    }

   std::pair<int,int> getFromAndToLocations(int arcIndex) { return std::make_pair(nodes[arcs[arcIndex].fromNodeIndex].state.lastVisited,arcs[arcIndex].location); }

    // Compilation
    void compileNgRoute(int s);
    void compileEmpty();
    void setupNgSets(int s);

    // Set arc coeffs to different values
    void setCoeffsAsDistances();
    void setCoeffsAsPreciseDistances();
    void setCoeffsAsDistancesMinusLagrangean(const std::vector<double>& lambda);
    void setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(const Dual& dual, LPSolveType lpSolveType);

    // Arc aggregation calculations
    void getNumberOfTimesLocationsCovered(std::unordered_map<int,double>& locationsCovered);
    void getCutSetValues(std::vector<double>& cutValues);
    void getCombValues(std::vector<double>& combValues);
    void getSrcCutValues(std::vector<double>& srcCutValues);

    void getNumberOfTimesLocationsCoveredRoutes(const Primal& primal, std::unordered_map<int,double>& locationsCovered);
    void getCutSetValuesRoutes(const Primal& primal, std::vector<double>& cutValues);
    void getCombValuesRoutes(const Primal& primal, std::vector<double>& combValues);
    void getSrcCutValuesRoutes(const Primal& primal, std::vector<double>& srcCutValues);

    int findSRCs(const Primal& primal, int limit, std::vector<double>& violations);
    void convertSolutionForVRPTWSep(std::vector<int>& edgeTail,
                                   std::vector<int>& edgeHead,
                                   std::vector<double>& edgeFlow);
    void getFeasiblePrimalIndices(const Primal& primal, std::vector<int>& feasibleIndices);

    // SSP - allows negative weights unlike dijkstra
    double computeShortestPathBFS(ShortestPathMode mode, std::vector<int>& routeByLocation);
    double computeShortestPathBFSWang(std::vector<int>& treeByParentArcs, std::vector<int>& routeByArcs, double& longestShortestPathLength);

    // Network flow for lp/ip
    double setupAndSolveFlowModel(FlowType flowType, IncludeCoverConstraints includeCoverConstraints, UseColumnGeneration useCg, const std::set<int>& initialPrimalArcIndices, Dual& dual, bool removeConstraintsForTesting, int timeoutSeconds);
    double getDualObjectiveValue(const Dual& dual, LPSolveType lpSolveType);
    double fixArcs(const Dual& dual, LPSolveType lpSolveType, double upperBound);
    double fixArcs(const std::vector<Dual>& dual, LPSolveType lpSolveType, double upperBound);
    double fixArcs(const std::deque<Dual>& dual, LPSolveType lpSolveType, double upperBound);

    // Primal heuristics
    double repairSolution(const std::vector<std::vector<int>>& feasibleSolution, const std::set<int>& destroyedElements, const Dual& dual, std::vector<std::vector<int>>& newBestRoutes, int timeoutSeconds);
 
    int selectArcWithLargestFlowFromNode(int nodeIndex);
    bool addLocationToRoute(int location, std::vector<int>& route);
    bool primalHeuristicGreedy(std::vector<std::vector<int>>& routesByLocation);
    void createTruncatedRoute(const std::vector<int>& route, std::vector<int>& truncatedRoute);
    bool prefixIntraRouteSwaps(std::vector<int>& route, double& routeCost);
    bool intraRouteSwaps(std::vector<int>& route, double& routeCost);

    int addPrimalHeuristicSuffixArc(int fromNodeIndex, const std::vector<int>& suffix, int suffixIndex);
    void generateHeuristicRoutesGreedy(std::vector<std::vector<int>>& routes);
    bool generateHeuristicRoutesLiterature(std::vector<std::vector<int>>& routes);
    double calculateInsertionCost(std::vector<int>& bestRoute, InsertionCriteria insertionCriteria, double gamma, int location);
    double insertLocationIntoRoute(std::vector<int>& bestRoute, int location);
    void populateEmptyRouteWithRandomLocation(std::vector<int>& route, std::vector<int>& candidateList);
    void sequentialInsertion(std::vector<std::vector<int>>& routes, std::vector<int> candidateList, InsertionCriteria insertionCriteria);
    void sequentialInsertionSimple(std::vector<std::vector<int>>& routes, std::vector<int> candidateList);
    void parallelInsertion(std::vector<std::vector<int>>& routes, std::vector<int> candidateList, InsertionCriteria insertionCriteria);
    void createMaximalSequence(const std::vector<int>& sequence, std::vector<int>& maximalSequence);

    // Separation methods
    bool doesRouteExistByArcs(const std::vector<int>& routeArcs, std::vector<int>& locationArcs) const;
    bool doesRouteExistByLocations(const std::vector<int>& routeByLocation, std::vector<int>& routeArcs) const;
    void decomposeRoutes(std::vector<int>& routeArcs, std::vector<double>& flows, std::vector<std::vector<int>>& routeDecomposition, std::vector<std::vector<int>>& decomposedArcs, DecompositionReason dr);
    void getSolutionArcs(std::set<int>& solutionArcs);
    void separateRoute(const std::vector<int>& routeArcs);
    bool isRouteFeasible(const std::vector<int>& route);
    void repairRoute(const std::vector<int>& route, std::vector<int>& feasibleRoute, std::vector<int>& feasibleRouteArcs);

    // Merge methods
    void findMergeNodesReducedCost(const Dual& dual, LPSolveType solveType, double limitToMerge);
    bool canMergeNodes(int nodeIndex1, int nodeIndex2);
    void mergeNodes(int nodeIndex1, int nodeIndex2);

    // Min cost flow for lagrangean
    void runDijkstra(ShortestPathMode mode, std::vector<int>& shortestPathByArc);
    int reverseArc(int arcIndex);
    void updateResidualGraph(const std::vector<int>& shortestPathByArc);
    double createSolutionFromReverseArcsAndReset();
    double solveMinCostFlowModel(const std::vector<double>& duals, std::vector<std::vector<int>>& shortestPathsByArc, bool& isDualFeasible, double& minReducedCost);

    // Min cost flow Wang for lagrangean
    void getTreeByChildArcsFromTreeByParentArcs(const std::vector<int>& treeByParentArcs, std::vector<std::vector<int>>& treeByChildArcs);
    void identifyNodesForUpdate(int branchNodeIndex, const std::vector<std::vector<int>>& treeByChildArcs, std::set<int>& nodesToUpdate);
    int findMultiPathNode(const std::set<int>& nodesToUpdate, int currentNumPaths);
    void findMultiPath(int newPathArcIndex, const std::vector<int>& treeByParentArcs, std::vector<int>& newShortestPathByArcs);
    void clipPermanentArc(int arcIndex);
    double createSolutionFromReverseArcsAndResetWang(const std::set<int>& clippedArcs, const std::vector<std::vector<int>>& shortestPathsByArc);
    void updateResidualGraphWang(const std::vector<int>& shortestPathByArc, std::vector<int>& treeByParentArc);
    void dijkstraWithBatchProc(std::vector<int>& treeByParentArcs, std::vector<std::vector<int>>& treeByChildArcs, std::set<int>& nodesToUpdate);
    double solveMinCostFlowModelWang(const Dual& dual, std::vector<std::vector<int>>& shortestPathsByArc, bool& isDualFeasible, double& shortestPathLength, double& longestShortestPathLength);
    bool checkFeasibleDual(const Dual& dual, LPSolveType lpSolveType);
    void getLocationsOnArcPaths(const std::vector<std::vector<int>>& shortestPathsByArc, std::set<int>& locations);

    // Debugging / testing
    const std::vector<VRPTWNode>& getNodes() const { return nodes; }
    const std::vector<VRPTWArc>& getArcs() const { return arcs; }
    bool checkAn32k5SolutionPossible() const;
    bool checkAn37k5SolutionPossible() const;
    bool checkAn33k6SolutionPossible() const;
    bool checkAn36k5SolutionPossible() const;
    bool checkAn48k7SolutionPossible() const;
    bool checkXn181k23SolutionPossible() const;
    bool checkC141SolutionPossible() const;
    bool checkLC121SolutionPossible() const;
    bool checkLRC121SolutionPossible() const;
    bool checkLC145SolutionCost();
    bool checkSinglePathPossible() const;
    void printCuts() const;
    void printFixedArcs() const;

    // Cuts
    bool addCapCutSet(const std::vector<int>& cutSet, const std::vector<int>& sequenceArcs, double rhs, RCCType rccType, bool useScaling);
    std::vector<int> getCapCutSet(int index) const { return capCutSets[index]; }
    double getCapCutSetRHS(int index) const { return capCutSetsRHS[index]; }
    double getCapCutSetScaling(int index) const { return capCutSetsScaling[index]; }
    int getNumCapCuts() const { return capCutSetsRHS.size(); }
    RCCType getRCCType(int index) const { return rccTypes[index]; }
    double getRCCCoeff(std::vector<int> route, int rccIndex);
    bool calculateRCCCoeff(int arcIndex, int rccIndex, double& coeff);

    void addCombCutTeeth(const std::vector<std::set<int>>& teeth);
    void addCombCutRHS(const double& rhs) { combRHS.push_back(rhs); }
    int getCombCutRHS(int index) const { return combRHS[index]; }
    int getNumCombCuts() const { return combRHS.size(); }

    int getNumSrcCuts() const { return srcCuts.size(); }
    int getSRCCoeff(int numTimesVisited, SRCType srcType);
    SRCType getSRCType(int index) { return srcCutTypes[index]; }
    int getSRCRHS(SRCType srcType);

  private:
    bool generateNewStateExact(VRPTWNodeState& state, int action);
    bool generateNewStateRelaxation(VRPTWNodeState& state, int action, int nodeIndex);

    int addNode(const VRPTWNodeState& state);
    int addArc(int fromNodeIndex, int toNodeIndex);
    bool moveArc(int arcIndex, int newToNodeIndex);
    void removeArc(int arcIndex);
    int addReverseArc(int forwardArcIndex);
 
    double getSRCFlowViolation(const Primal& primal, const std::set<int>& cutSet, SRCType srcType);
    void getSRCArcsAndCoeffs(const Primal& primal, const std::set<int>& cutSet, SRCType srcType, std::vector<int>& srcArcs, std::vector<int>& srcCoeffs);
    bool checkExistingSrcCuts(const std::set<int>& testSet, std::vector<int>& bestArcSet, std::vector<int>& bestCoeffs, SRCType srcType);
    int getNumTimesSetVisited(std::vector<int> route, const std::set<int>& cutSet);
    void trimRouteToLocations(std::vector<int>& routeArcs, const std::set<int>& cutSet);
    bool isArcIncrementalSRC(int numTimesVisited, SRCType srcType);
    void getRouteSRCArcs(const std::vector<int>& routeArcs, const std::set<int>& cutSet, SRCType srcType, std::set<int>& arcSet);
    int getSRC3Coeff(int numTimesVisited);
    int getSRC4Coeff(int numTimesVisited);
    int getSRC5V1Coeff(int numTimesVisited);
    int getSRC5V2Coeff(int numTimesVisited);

    VRPTW vrptw;
    VRPTWDDParameters params;
    std::vector<VRPTWNode> nodes;
    std::vector<VRPTWArc> arcs;
    std::set<int> fixedArcs;
    std::set<int> removedArcs;
    std::vector<std::set<int>> ngSets;

    std::vector<std::set<int>> compilationAllVisitedDown;

    int rootNodeIndex;
    int terminalNodeIndex;

    std::unordered_map<VRPTWNodeState,int,hash_state> stateToNodes;
    std::unordered_map<int,std::vector<int>> locationToArcs;
    std::unordered_map<int,int> arcReverseArc;
    std::map<VRPTWNodeState,std::vector<int>> nodeOrdering;

    // Discretization parameters
    int timeStepSize;
    int loadStepSize;
    int timeWindowBinary;
    int capacityBinary;
    int counterBinary;

    // Rounded Capacity Cuts
    std::vector<std::vector<int>> capCutSets;
    std::vector<double> capCutSetsRHS;
    std::vector<double> capCutSetsScaling;
    std::vector<std::vector<int>> capCutSetArcs;
    std::vector<std::vector<double>> capCutSetCoeffs;
    std::vector<RCCType> rccTypes;

    // SRC Cuts (SRC)
    std::vector<std::set<int>> srcCuts;
    std::vector<std::vector<int>> srcCutSeparatedArcs;
    std::vector<std::vector<int>> srcCutSeparatedCoeffs;
    std::vector<SRCType> srcCutTypes;
    int separatedFeasibleRouteCounter;

    // Strengthened Comb Cuts
    std::vector<std::vector<std::set<int>>> teeths;
    std::vector<int> combRHS;
    std::vector<std::vector<int>> combCutArcs;
};

#endif
