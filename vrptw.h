#ifndef VRPTW_H
#define VRPTW_H

#define DBG(x)

#include <string>
#include <cmath>
#include <set>
#include <regex>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "vrptwstatic.h"

const long INF = 1e9;

enum ProblemType
{
  CVRP = 0,
  TW = 1,
  TSPTW = 2,
  SOP = 3,
  PDP = 4
};

enum FixedNumPaths
{
  FIXED_NUM_PATHS = 0,
  FLEXIBLE_NUM_PATHS = 1
};

enum CircuitOrPath
{
  CIRCUIT = 0,
  PATH = 1
};

enum StateSpace
{
  NG = 0,
  Q = 1
};

enum PricingProblemType
{
  DD = 0,
  DP = 1
};

enum LPSolveType
{
  LPSolver = 0,
  LAGSolver = 1
};

enum DecompositionReason
{
  SEPARATE = 0,
  DECOMPOSE = 1
};

enum IncludeCoverConstraints
{
  Y = 0,
  N = 1
};

enum UseColumnGeneration
{
  USE_CG = 0,
  NO_CG = 1
};

enum FlowType
{
  LP = 0,
  IP = 1
};

enum AllowMultipleVisits
{
  YES = 0,
  NO = 1
};

enum SRCType
{
  SRC3 = 0,
  SRC4 = 1,
  SRC5V1 = 2,
  SRC5V2 = 2
};

enum RCCType
{
  Type1 = 0,
  Type3 = 1
};


enum SGDAlgorithm
{
  SGD = 0,
  VA = 1
};

enum ShortestPathMode
{
  UPDATE_POTENTIALS = 0,
  SHORTEST_PATH = 1
};

enum VRPTWCapacityType
{
  RELAX_CAPACITY = 0,
  NO_RELAX_CAPACITY = 1,
};

enum VRPTWCounterType
{
  USE_COUNTER = 0,
  NO_USE_COUNTER = 1
};

enum VRPTWTimeWindowType
{
  TIME_WINDOWS = 0,
  NO_TIME_WINDOWS = 1
};

enum PrimalHeuristic
{
  STANDALONE = 0,
  CE_MIP = 1,
  CE_MIP_LNS = 2,
  CE_GREEDY = 3,
  BEST_KNOWN = 4
};

enum PrimalHeuristicRemovalStrategy
{
  RANDOM = 0,
  SHAW = 1,
  WORST = 2,
  SEQUENCE_RANDOM = 3,
  SEQUENCE_SHAW = 4,
  SEQUENCE_WORST = 5
};

enum PrimalHeuristicInsertionAblation
{
  NO_ABLATION = 0,
  INTRA_SWAP = 1,
  TRUNCATED = 2
};

enum InsertionCriteria
{
  MCFIC = 0,
  NFIC = 1
};

enum LNSNodeLimitMethod
{
  ALL = 0,
  NUM_DESTROYED_REDUCED_COST = 1,
  NUM_DESTROYED_DISTANCE = 2
};

struct VRPTWDDParameters
{
  int timeoutSeconds;

  LPSolveType lpSolveType;
  StateSpace stateSpace;
  int ngSetSize;
  bool changeToLP;
  bool useSeparations;

  bool useVariableFixing;
  bool useMuSSP;
  bool repairDuals;

  int numLagCuts;
  bool useRCCs;
  bool useSRC3s;
  bool useSRC4s;
  bool useSRC5V1s;
  bool useSRC5V2s;
  bool useVolumeAlgorithm;
  bool limitRCCs;
  bool useScaling;

  PrimalHeuristic primalHeuristic;
  int lnsTimeoutSeconds;
  PrimalHeuristicRemovalStrategy removalStrategy;
  PrimalHeuristicInsertionAblation insertionAblation;
};

struct VRPTWSolution
{
  VRPTWSolution(std::vector<std::vector<int>> _routes, double _totalDistance) : routes(_routes), totalDistance(_totalDistance) {}

  std::vector<std::vector<int>> routes;
  double totalDistance;
};

struct VRPTW
{
  public:
    VRPTW(std::string _fileName);

    double evaluateRouteDistance(const std::vector<int>& routeByLocation);
    double evaluateSolutionCost(const std::vector<std::vector<int>>& routesByLocation);
    double calculateLocationPickupTime(const std::vector<int>& route, int routeIndex);

    int getLongestPossibleRoute();
    void recomputeDistancesPDPTW();

    void calculateRouteLengthUpperBound();
    void calculateMaxValues();

    void writeVrpFormat();

    std::string fileName;
    bool isInstanceClosed;

    int capacity;
    std::vector<int> xCoords;
    std::vector<int> yCoords;
    std::vector<std::vector<double> > distances;
    std::vector<std::vector<double> > preciseDistances;
    std::vector<std::set<int> > precedences;
    std::vector<std::set<int> > reliances;
    std::vector<int> demands;
    std::vector<double> demandsForSeparation;
    std::vector<int> demandsForCombs;

    std::vector<int> startTimes;
    std::vector<int> endTimes;
    std::vector<int> serviceTimes;

    int depot;
    int numLocations;
    int numVehicles;

    VRPTWCapacityType vrptwCapacityType;
    VRPTWTimeWindowType vrptwTimeWindowType;
    VRPTWCounterType counterType;
    FixedNumPaths fixedNumPaths;
    CircuitOrPath circuitOrPath;
    ProblemType problemType;
    int timeStateMultiplier;
    int timeStateDiscretization;
    int finalTimeStateMultiplier;
    int loadDiscretization;
    double instanceUpperBound;
    int routeLengthUpperBound;

    int maxStartTime;
    int maxDemand;
    double maxDistance;
    int fixedRouteCost;
};

#endif
