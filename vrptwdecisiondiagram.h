#ifndef VRPTWDD
#define VRPTWDD

#include <vector>
#include <deque>
#include <set>
#include <unordered_map>
#include <map>

#include "vrptw.h"

struct VRPTWNodeState
{
  VRPTWNodeState(int _counter, int _capacity, int _timeWithMultiplier, int _lastVisited, std::set<int> _visited) : counter(_counter), capacity(_capacity), timeWithMultiplier(_timeWithMultiplier), lastVisited(_lastVisited), visited(_visited) {};
  VRPTWNodeState(const VRPTWNodeState& state) : counter(state.counter), capacity(state.capacity), timeWithMultiplier(state.timeWithMultiplier), lastVisited(state.lastVisited), visited(state.visited) {};

  bool operator==(const VRPTWNodeState & rhs) const
  {
    if ((counter == rhs.counter) && (capacity == rhs.capacity) && (timeWithMultiplier == rhs.timeWithMultiplier) && (lastVisited == rhs.lastVisited) && (visited == rhs.visited))
    {
      return true;
    }
    else
    {
      return false;
    }
  }

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
      if (capacity < rhs.capacity)
      {
        return true;
      }
      else if (capacity > rhs.capacity)
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
  int capacity;
  int timeWithMultiplier;
  int lastVisited;
  std::set<int> visited;
};

struct hash_state
{
  std::size_t operator()(const VRPTWNodeState& state) const
  {
    std::size_t h1 = std::hash<int>()(state.counter);
    std::size_t h2 = std::hash<int>()(state.capacity);
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
      return fixedArcs.size() * 100.0 / arcs.size();
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

    // compilation
    void compileExactFukasawa(int s);
    void compileNgRoute(int s);

    // set arc coeffs to different values
    void setCoeffsAsDistances();
    void setCoeffsAsDistancesMinusLagrangean(const std::vector<double>& lambda);
    void setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(const std::vector<double>& lambda, const std::vector<double>& capDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, LPSolveType solveType);

    // arc aggregation calculations
    void getNumberOfTimesLocationsCovered(std::unordered_map<int,double>& locationsCovered);
    void getCutSetValues(std::vector<double>& cutValues);
    void getCombValues(std::vector<double>& combValues);
    void getCliqueCutValues(std::vector<double>& cliqueCutValues);
    void addConnectedNodesToBlacklist(int nodeIndex, int demandLimit, std::set<int>& blacklist, const std::set<int>& nodesUsed);
    bool areNodesConnected(int nodeIndex1, int nodeIndex2);
    //void findSRCThree(const std::vector<std::set<int>>& decomposition, const std::vector<double>& stepSizes, bool& cutAdded);
    void convertSolutionForVRPTWSep(std::vector<int>& edgeTail,
                                   std::vector<int>& edgeHead,
                                   std::vector<double>& edgeFlow,
                                   std::set<int>& rccArcs);

    // can be used for col gen or ssp, allows negative weights unlike dijkstra
    double computeShortestPathBFS(ShortestPathMode mode, std::vector<int>& routeByLocation);
    double computeShortestPathBFSWang(std::vector<int>& treeByParentArcs, std::vector<int>& routeByArcs);
    void addColumnForLPCG(const std::vector<int>& route);

    // network flow for lp/ip
    void initializeColumnsByLPDecomp();
    double setupAndSolveFlowModel(FlowType flowType, IncludeCoverConstraints includeCoverConstraints, UseColumnGeneration useCg, std::vector<double>& duals, double& singlePathDual, std::vector<double>& capDuals, std::vector<double>& combDuals, std::vector<double>& srcDuals);
    double fixArcs(const std::vector<double>& lambda, double singlePathDual, const std::vector<double>& capDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, double lowerBound, LPSolveType solveType);

    // primal heuristic and separation methods
    void primalHeuristic(std::vector<std::vector<int>>& routesByLocation);
    bool doesRouteExistByArcs(const std::vector<int>& routeArcs) const;
    bool doesRouteExistByLocations(const std::vector<int>& routeByLocation) const;
    void decomposeRoutes(std::vector<int>& routeArcs, std::vector<double>& flows, std::vector<std::vector<int>>& routeDecomposition, int maxS, DecompositionReason dr);
    void getSolutionArcs(std::set<int>& solutionArcs);
    void separateInfeasibleRoute(const std::vector<int>& routeArcs, int maxS);

    // min cost flow for lagrangean
    void runDijkstra(ShortestPathMode mode, std::vector<int>& shortestPathByArc);
    int reverseArc(int arcIndex);
    void updateResidualGraph(const std::vector<int>& shortestPathByArc);
    double createSolutionFromReverseArcsAndReset();
    double solveMinCostFlowModel(const std::vector<double>& duals, std::vector<std::vector<int>>& shortestPathsByArc, bool& isDualFeasible, double& minReducedCost);

    // min cost flow Wang for lagrangean
    void getTreeByChildArcsFromTreeByParentArcs(const std::vector<int>& treeByParentArcs, std::vector<std::vector<int>>& treeByChildArcs);
    void identifyNodesForUpdate(int branchNodeIndex, const std::vector<std::vector<int>>& treeByChildArcs, std::set<int>& nodesToUpdate);
    int findMultiPathNode(const std::set<int>& nodesToUpdate, int currentNumPaths);
    void findMultiPath(int newPathArcIndex, const std::vector<int>& treeByParentArcs, std::vector<int>& newShortestPathByArcs);
    void clipPermanentArc(int arcIndex);
    double createSolutionFromReverseArcsAndResetWang(const std::set<int>& clippedArcs, const std::vector<std::vector<int>>& shortestPathsByArc);
    void updateResidualGraphWang(const std::vector<int>& shortestPathByArc, std::vector<int>& treeByParentArc);
    void dijkstraWithBatchProc(std::vector<int>& treeByParentArcs, std::vector<std::vector<int>>& treeByChildArcs, std::set<int>& nodesToUpdate);
    double solveMinCostFlowModelWang(const std::vector<double>& duals, const std::vector<double>& capDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, std::vector<std::vector<int>>& shortestPathsByArc, bool& isDualFeasible, double& shortestPathLength);
    bool checkFeasibleDual(const std::vector<double>& lambda, const std::vector<double>& rccDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, LPSolveType solveType);
    void getLocationsOnArcPaths(const std::vector<std::vector<int>>& shortestPathsByArc, std::set<int>& locations);

    // for testing
    const std::vector<VRPTWNode>& getNodes() const { return nodes; }
    const std::vector<VRPTWArc>& getArcs() const { return arcs; }
    bool checkC141SolutionPossible() const;
    bool checkLC121SolutionPossible() const;
    bool checkLRC121SolutionPossible() const;
    bool checkSinglePathPossible() const;

    // for cuts
    void addCapCutSet(const std::vector<int>& cutSet, const std::set<int>& rccArcs);
    void addCapCutSetRHS(const double& rhs) { capCutSetsRHS.push_back(rhs); }
    int getCapCutSetRHS(int index) const { return capCutSetsRHS[index]; }
    int getNumCapCuts() const { return capCutSetsRHS.size(); }

    int getNumCliqueCuts() const { return cliqueCuts.size(); }
    bool isCliqueCutActive(int index);

    void addCombCutTeeth(const std::vector<std::set<int>>& teeth);
    void addCombCutRHS(const double& rhs) { combRHS.push_back(rhs); }
    int getCombCutRHS(int index) const { return combRHS[index]; }
    int getNumCombCuts() const { return combRHS.size(); }

  private:
    int addNode(const VRPTWNodeState& state);
    int addArc(int fromNodeIndex, int toNodeIndex);
    int addReverseArc(int forwardArcIndex);
    bool generateNewStateFromExact(VRPTWNodeState& state, int action);
    bool moveArc(int arcIndex, int newToNodeIndex);
    void removeArc(int arcIndex);
 
    VRPTW vrptw;
    VRPTWDDParameters params;
    std::vector<VRPTWNode> nodes;
    std::vector<VRPTWArc> arcs;
    std::set<int> fixedArcs;
    std::set<int> removedArcs;

    int rootNodeIndex;
    int terminalNodeIndex;

    std::unordered_map<VRPTWNodeState,int,hash_state> stateToNodes;
    std::unordered_map<int,std::vector<int>> locationToArcs;
    std::unordered_map<int,int> arcReverseArc;
    std::map<VRPTWNodeState,std::vector<int>> nodeOrdering;

    std::vector<bool> arcsUsed;

    // bucket parameters
    int timeStepSize;
    int capacityStepSize;

    // data structures for cuts
    // Rounded Capacity Cuts
    std::vector<std::vector<int>> capCutSets;
    std::vector<double> capCutSetsRHS;
    std::vector<std::vector<std::vector<int>>> capCutSetRoutes;
    std::vector<std::vector<int>> capCutSetArcs;

    // Clique Cuts
    std::set<std::set<int>> cliqueCuts;
    std::vector<std::vector<int>> cliqueCutArcs;
    std::vector<bool> cliqueCutActive;

    // Strengthened Comb Cuts
    std::vector<std::vector<std::set<int>>> teeths;
    std::vector<int> combRHS;
    std::vector<std::vector<int>> combCutArcs;
};

#endif
