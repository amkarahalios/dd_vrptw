#include <bits/stdc++.h>
#include <cmath>
#include <random>
#include <algorithm>
#include <chrono>

#include <ilcplex/ilocplex.h>
#include "vrptwdecisiondiagram.h"

void VRPTWNode::printForDD(int nodeIndex) const
{
  std::cout << nodeIndex << ": {";
  for (int loc : state.visited)
  {
    std::cout << loc << " ";
  }
  std::cout << "c: ";
  std::cout << "cnt: " << state.counter;
  std::cout << " cap: " << state.load;
  std::cout << " time: " << state.timeWithMultiplier;
  std::cout << " curr: " << state.lastVisited;
  std::cout << " pot: " << potential;
  std::cout << " sp: " << shortestPathDistance;
  std::cout << " }";
};

void VRPTWNode::print() const
{
  std::cout << "counter: " << state.counter << " ";
  std::cout << "load: " << state.load << " ";
  std::cout << "time: " << state.timeWithMultiplier << " ";
  std::cout << "location: " << state.lastVisited << " ";
  std::cout << "visited: ";
  for (int location : state.visited)
  {
    std::cout << location << " ";
  }
  std::cout << "out arc indices: ";
  for (int arcIndex : outArcs)
  {
    std::cout << arcIndex << ",";
  }
  std::cout << std::endl;
}

void VRPTWArc::print() const
{
  std::cout << "fromNodeIndex: " << fromNodeIndex << " ";
  std::cout << "toNodeIndex: " << toNodeIndex << " ";
  std::cout << "location: " << location << " ";
  std::cout << "distance: " << distance << " ";
  std::cout << "isReverseArc: " << isReverseArc << " ";
  std::cout << std::endl;
}

VRPTWDecisionDiagram::VRPTWDecisionDiagram(const VRPTW& _vrptw, const VRPTWDDParameters _params) : vrptw(_vrptw), params(_params)
{
  // default time and load step sizes
  loadStepSize = 1 * (*std::min_element(vrptw.demands.begin()+1, vrptw.demands.end()));
  if (vrptw.vrptwCapacityType == NO_RELAX_CAPACITY)
  {
    loadStepSize = 1;
  }

  if (vrptw.loadDiscretization != 0)
  {
    loadStepSize = vrptw.loadDiscretization;
  }

  setTimeStepSize();
  separatedFeasibleRouteCounter = 0;
 
  // create some binary values for initial relaxation to use different state information
  timeWindowBinary = 1;
  if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::NO_TIME_WINDOWS)
  {
    timeWindowBinary = 0;
  }

  capacityBinary = 1;
  if (vrptw.vrptwCapacityType == VRPTWCapacityType::RELAX_CAPACITY)
  {
    capacityBinary = 0;
  }

  counterBinary = 1;
  if (vrptw.counterType == VRPTWCounterType::NO_USE_COUNTER)
  {
    counterBinary = 0;
  }
};

void VRPTWDecisionDiagram::setTimeStepSize()
{
  // default time and load step sizes
    timeStepSize = vrptw.timeStateMultiplier * vrptw.timeStateDiscretization;
};

void VRPTWDecisionDiagram::updateTimeStateMultiplierByTen()
{
  // update time step
  setTimeStepSize();

  // update node states and arc distances
  for (int nodeIndex=0; nodeIndex<nodes.size(); ++nodeIndex)
  {
    nodes[nodeIndex].state.timeWithMultiplier = nodes[nodeIndex].state.timeWithMultiplier * 10;
  }

  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    int fromLocation = nodes[arcs[arcIndex].fromNodeIndex].state.lastVisited;
    int toLocation = nodes[arcs[arcIndex].toNodeIndex].state.lastVisited;
    arcs[arcIndex].distance = vrptw.distances[fromLocation][toLocation];
  }
};

void VRPTWDecisionDiagram::print() const
{
  std::cout << std::endl << "digraph D {" << std::endl;
  for (int nodeIndex=0; nodeIndex<nodes.size(); ++nodeIndex)
  {
    for (int arcIndex : nodes[nodeIndex].outArcs)
    {
      nodes[nodeIndex].printForDD(nodeIndex);
      std::cout << "\"";
      std::cout << " -> " << arcIndex << ": [loc: " << arcs[arcIndex].location << ", fl: " << arcs[arcIndex].heuristicFlow << ", d-l: ";
      std::cout << arcs[arcIndex].coeff << ", rc: ";
      std::cout << arcs[arcIndex].cijPi << ", d: ";
      std::cout << arcs[arcIndex].distance << ", rev: " << arcs[arcIndex].isReverseArc << "] -> ";
      int toNodeIndex = arcs[arcIndex].toNodeIndex;
      nodes[toNodeIndex].printForDD(toNodeIndex);
      std::cout << "\"";
      std::cout << std::endl;
    }
  }
  std::cout << std::endl << "}" << std::endl;
};

double VRPTWDecisionDiagram::evaluateRouteCost(const std::vector<int>& routeByArc)
{
  double cost = 0;
  for (auto arcIndex : routeByArc)
  {
    cost = cost + arcs[arcIndex].coeff;
  }

  return cost;
};

int VRPTWDecisionDiagram::addNode(const VRPTWNodeState& state)
{
  int countStatesInMap = stateToNodes.count(state);
  bool stateExists = (countStatesInMap == 1);
  if (!stateExists)
  {
    VRPTWNode newNode(state);
    int newNodeIndex = nodes.size();

    nodes.push_back(newNode);
    stateToNodes.insert(std::make_pair(state, newNodeIndex));
    nodeOrdering[state].push_back(newNodeIndex);
    return newNodeIndex;
  }

  return stateToNodes[state];
};

int VRPTWDecisionDiagram::addArc(int fromNodeIndex, int toNodeIndex)
{
  int fromLocation = nodes[fromNodeIndex].state.lastVisited;
  int toLocation = nodes[toNodeIndex].state.lastVisited;
  double distance = vrptw.distances[fromLocation][toLocation];

  VRPTWArc newArc(fromNodeIndex, toNodeIndex, toLocation, distance);
  int newArcIndex = arcs.size();
  arcs.push_back(newArc);

  nodes[fromNodeIndex].outArcs.push_back(newArcIndex);
  nodes[toNodeIndex].inArcs.push_back(newArcIndex);

  locationToArcs[toLocation].push_back(newArcIndex);

  if (!arcs[newArcIndex].isReverseArc)
  {
    // add to appropriate cut sets
    for (int capCutIndex=0; capCutIndex<capCutSets.size(); ++capCutIndex)
    {
      RCCType rccType = rccTypes[capCutIndex];

      double coeff = 0;
      bool nonZeroCoeff = calculateRCCCoeff(newArcIndex, capCutIndex, coeff);
      if (nonZeroCoeff)
      {
        double scaling = capCutSetsScaling[capCutIndex];
        capCutSetArcs[capCutIndex].push_back(newArcIndex);
        capCutSetCoeffs[capCutIndex].push_back(coeff * scaling);
      }
    }

    // use duals for combs
    for (int combIndex=0; combIndex<combRHS.size(); ++combIndex)
    {
      auto teeth = teeths[combIndex];
      int fromLoc = nodes[arcs[newArcIndex].fromNodeIndex].state.lastVisited;
      int toLoc = nodes[arcs[newArcIndex].toNodeIndex].state.lastVisited;
      for (int toothIndex=0; toothIndex<teeth.size(); ++toothIndex)
      {
        auto tooth = teeth[toothIndex];
        bool fromLocInSet = (std::find(tooth.begin(), tooth.end(), fromLoc) != tooth.end());
        bool toLocInSet = (std::find(tooth.begin(), tooth.end(), toLoc) != tooth.end());
        if ((fromLocInSet && !toLocInSet) || (!fromLocInSet && toLocInSet))
        {
          // keep this collection as a vector because this can be called multiple times
          combCutArcs[combIndex].push_back(newArcIndex);
        }
      }
    }

    // src cuts
    if (nodes[toNodeIndex].state.isExact)
    {
      auto visitedLocationsWithArc = nodes[toNodeIndex].state.visited;
      for (int srcIndex=0; srcIndex<srcCuts.size(); ++srcIndex)
      {
        auto cutSet = srcCuts[srcIndex];
        auto srcType = srcCutTypes[srcIndex];
        if (cutSet.find(toLocation) != cutSet.end())
        {
          std::set<int> locationsOverlap;
          for (int location : cutSet)
          {
            if (visitedLocationsWithArc.find(location) != visitedLocationsWithArc.end())
            {
              locationsOverlap.insert(location);
            }
          }

          bool isIncrementalArc = isArcIncrementalSRC(static_cast<int>(locationsOverlap.size()), srcType);
          if (isIncrementalArc)
          {
            srcCutSeparatedArcs[srcIndex].push_back(newArcIndex);
            srcCutSeparatedCoeffs[srcIndex].push_back(1);
          }
        }
      }
    }
  }

  return newArcIndex;
}
      
int VRPTWDecisionDiagram::addPrimalHeuristicSuffixArc(int fromNodeIndex, const std::vector<int>& suffix, int suffixIndex)
{
  int loc1 = nodes[fromNodeIndex].state.lastVisited;
  double distance = 0;
  for (int index=0; index<suffix.size(); ++index)
  {
    int loc2 = suffix[index];
    distance += vrptw.distances[loc1][loc2];
    loc1 = loc2;
  }

  // suffix indicator, change index by 1 so 0 isn't like going to the depot
  int suffixIndicator = (-1 * suffixIndex) - 1;
  VRPTWArc newArc(fromNodeIndex, terminalNodeIndex, suffixIndicator, distance);
  int newArcIndex = arcs.size();
  arcs.push_back(newArc);

  nodes[fromNodeIndex].outArcs.push_back(newArcIndex);
  nodes[terminalNodeIndex].inArcs.push_back(newArcIndex);

  for (int loc : suffix)
  {
    locationToArcs[loc].push_back(newArcIndex);
  }

  return newArcIndex;
}

// Tries to insert pickup/delivery at best indices in route
// Maybe this should build the transition graph, might save time?
// Basically build the route nodes and try the pickup/delivery transition at each iteration?
double VRPTWDecisionDiagram::calculateInsertionCost(std::vector<int>& route, InsertionCriteria insertionCriteria, double gamma, int location)
{
  int pickup = -1;
  int delivery = -1;
  if (!vrptw.reliances.empty() && !vrptw.precedences.empty())
  {
    if (!vrptw.precedences[location].empty())
    {
      return INF;
    }
    else
    {
      delivery = *(vrptw.reliances[location].begin());
      pickup = location;
    }
  }
  else
  {
    pickup = location;
  }

  // Kind of like repair solution with one route and two locations for extra transitions?
  double bestInsertionCost = INF;
  std::vector<int> bestRoute;
  for (int insertionIndex1=1; insertionIndex1<route.size(); ++insertionIndex1)
  {
    std::vector<int> newRoute;
    for (int routeIndex=0; routeIndex<route.size(); ++routeIndex)
    {
      if (insertionIndex1 == routeIndex)
      {
        newRoute.push_back(pickup);
      }
      newRoute.push_back(route[routeIndex]);
    }

    newRoute.pop_back();
    if (!isRouteFeasible(newRoute))
    {
      continue;
    }
    newRoute.push_back(0);
 
    std::vector<int> newRouteFinal;
    if (delivery != -1)
    {
      for (int insertionIndex2=1; insertionIndex2<newRoute.size(); ++insertionIndex2)
      {
        newRouteFinal = newRoute;
        newRouteFinal.insert(newRouteFinal.begin()+insertionIndex2, delivery);
      }
    }
    else
    {
      newRouteFinal = newRoute;
    }

    if (isRouteFeasible(newRouteFinal))
    {
      double routeCost = vrptw.evaluateRouteDistance(newRouteFinal);
      if (routeCost < bestInsertionCost)
      {
        bestInsertionCost = routeCost;
        bestRoute = newRouteFinal;
      }
    }
  }

  route.clear();
  for (int loc : bestRoute)
  {
    route.push_back(loc);
  }

  return bestInsertionCost;
}

void VRPTWDecisionDiagram::parallelInsertion(std::vector<std::vector<int>>& routes, std::vector<int> candidateList, InsertionCriteria insertionCriteria)
{
  double gamma = 0.0;
  if (insertionCriteria == InsertionCriteria::MCFIC)
  {
    int randomInteger = std::rand();
    int randomMultiplier = randomInteger % 18;
    gamma = 0.1 * randomMultiplier;
  }

  while (!candidateList.empty())
  {
    bool insertedLocation = false;
    double bestInsertionCost = INF;
    int bestInsertionLocation = -1;
    int bestInsertionIndex = -1;
    int bestInsertionRouteIndex = -1;
    for (int routeIndex=0; routeIndex<routes.size(); ++routeIndex)
    {
      for (int insertionLocation : candidateList)
      {
        std::vector<int> insertionRoute = routes[routeIndex];
        double insertionCost = calculateInsertionCost(insertionRoute, insertionCriteria, gamma, insertionLocation);
      }
    }

    // Add new vehicle if needed
    if (!insertedLocation)
    {
      std::vector<int> emptyRoute;
      emptyRoute.push_back(0);
      emptyRoute.push_back(0);
      routes.push_back(emptyRoute);
    }
  }
}

double VRPTWDecisionDiagram::insertLocationIntoRoute(std::vector<int>& route, int location)
{
  int pickup = -1;
  int delivery = -1;
  if (!vrptw.reliances.empty() && !vrptw.precedences.empty())
  {
    if (!vrptw.precedences[location].empty())
    {
      return INF;
    }
    else
    {
      delivery = *(vrptw.reliances[location].begin());
      pickup = location;
    }
  }
  else
  {
    pickup = location;
  }

  for (int insertionIndex1=1; insertionIndex1<route.size(); ++insertionIndex1)
  {
    std::vector<int> newRoute;
    for (int routeIndex=0; routeIndex<route.size(); ++routeIndex)
    {
      if (insertionIndex1 == routeIndex)
      {
        newRoute.push_back(pickup);
      }
      newRoute.push_back(route[routeIndex]);
    }

    newRoute.pop_back();
    if (!isRouteFeasible(newRoute))
    {
      continue;
    }
    newRoute.push_back(0);
 
    std::vector<int> newRouteFinal;
    if (delivery != -1)
    {
      for (int insertionIndex2=1; insertionIndex2<newRoute.size(); ++insertionIndex2)
      {
        newRouteFinal = newRoute;
        newRouteFinal.insert(newRouteFinal.begin()+insertionIndex2, delivery);
      }
    }
    else
    {
      newRouteFinal = newRoute;
    }
 
    if (isRouteFeasible(newRouteFinal))
    {
      double routeCost = vrptw.evaluateRouteDistance(newRouteFinal);
      route.clear();
      for (int loc : newRouteFinal)
      {
        route.push_back(loc);
      }
      return routeCost;
    }
  }

  route.clear();
  return INF;
}

void VRPTWDecisionDiagram::populateEmptyRouteWithRandomLocation(std::vector<int>& route, std::vector<int>& candidateList)
{
  // Randomly assign it locations
  std::vector<int> newRoute = {0, 0};
  while (newRoute.size() == 2)
  {
    int randomInteger0 = std::rand();
    int randomCandidateListIndex = randomInteger0 % candidateList.size();
    int randomLocation = candidateList[randomCandidateListIndex];
    newRoute.insert(newRoute.begin()+1, randomLocation);
    candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), randomLocation), candidateList.end());
    if (!vrptw.reliances.empty())
    {
      for (int relianceLocation : vrptw.reliances[randomLocation])
      {
        newRoute.insert(newRoute.end()-1, relianceLocation);
        candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), relianceLocation), candidateList.end());
      }
    }
    if (!isRouteFeasible(newRoute))
    {
      for (int newRouteLocation : newRoute)
      {
        if (newRouteLocation != 0)
        {
          candidateList.push_back(newRouteLocation);
        }
      }
      newRoute = {0,0};
    }
  }

  route.clear();
  for (int loc : newRoute)
  {
    route.push_back(loc);
  }
}

void VRPTWDecisionDiagram::sequentialInsertionSimple(std::vector<std::vector<int>>& routes, std::vector<int> candidateList)
{
  // make first route
  std::vector<int> newRoute;
  populateEmptyRouteWithRandomLocation(newRoute, candidateList);
  routes.push_back(newRoute);

  int currRouteIndex = 0;
  while (!candidateList.empty())
  {
    bool isLocationInserted = false;
    std::vector<int> routeWithInsertion;
    int locationToInsert;
    for (int insertionLocation : candidateList)
    {
      std::vector<int> insertionRoute = routes[currRouteIndex];
      double insertionCost = insertLocationIntoRoute(insertionRoute, insertionLocation);
      if (insertionCost < INF-1)
      {
        routeWithInsertion = insertionRoute;
        locationToInsert = insertionLocation;
        break;
      }
    }

    if (!routeWithInsertion.empty())
    {
      isLocationInserted = true;
      candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), locationToInsert), candidateList.end());
      if (!vrptw.precedences.empty() && !vrptw.precedences[locationToInsert].empty())
      {
        int secondLocation = *(vrptw.precedences[locationToInsert].begin());
        candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), secondLocation), candidateList.end());
      }
      if (!vrptw.reliances.empty() && !vrptw.reliances[locationToInsert].empty())
      {
        int secondLocation = *(vrptw.reliances[locationToInsert].begin());
        candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), secondLocation), candidateList.end());
      }
      routes[currRouteIndex] = routeWithInsertion;
    }

    // Add new vehicle if needed
    if (!isLocationInserted)
    {
      std::vector<int> newRoute;
      populateEmptyRouteWithRandomLocation(newRoute, candidateList);
      routes.push_back(newRoute);
      ++currRouteIndex;
    }
  }
}

void VRPTWDecisionDiagram::sequentialInsertion(std::vector<std::vector<int>>& routes, std::vector<int> candidateList, InsertionCriteria insertionCriteria)
{
  double gamma = 0.0;
  if (insertionCriteria == InsertionCriteria::MCFIC)
  {
    int randomInteger = std::rand();
    int randomMultiplier = randomInteger % 18;
    gamma = 0.1 * randomMultiplier;
  }

  while (!candidateList.empty())
  {
    // For each route, try to insert a location in the best index possible
    bool insertedLocation = false;
    for (int routeIndex=0; routeIndex<routes.size(); ++routeIndex)
    {
      double bestInsertionCost = INF;
      std::vector<int> bestNewRoute;
      int bestLocation = -1;
      for (int insertionLocation : candidateList)
      {
        if (!vrptw.precedences.empty() && !vrptw.precedences[insertionLocation].empty())
        {
          continue;
        }
        std::vector<int> insertionRoute = routes[routeIndex];
        double insertionCost = calculateInsertionCost(insertionRoute, insertionCriteria, gamma, insertionLocation);
        if (insertionCost < bestInsertionCost)
        {
          bestInsertionCost = insertionCost;
          bestNewRoute = insertionRoute;
          bestLocation = insertionLocation;
          break;
        }
      }

      if (!bestNewRoute.empty())
      {
        insertedLocation = true;
        candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), bestLocation), candidateList.end());
        if (!vrptw.precedences.empty() && !vrptw.precedences[bestLocation].empty())
        {
          int secondLocation = *(vrptw.precedences[bestLocation].begin());
          candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), secondLocation), candidateList.end());
        }
        if (!vrptw.reliances.empty() && !vrptw.reliances[bestLocation].empty())
        {
          int secondLocation = *(vrptw.reliances[bestLocation].begin());
          candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), secondLocation), candidateList.end());
        }
        routes[routeIndex] = bestNewRoute;
      }
    }

    // Add new vehicle if needed
    if (!insertedLocation)
    {
      std::vector<int> newRoute;
      newRoute.push_back(0);
      newRoute.push_back(0);

      // Randomly assign it locations
      while (static_cast<int>(newRoute.size()) == 2)
      {
        int randomInteger0 = std::rand();
        int randomCandidateListIndex = randomInteger0 % candidateList.size();
        int randomLocation = candidateList[randomCandidateListIndex];
        newRoute.insert(newRoute.begin()+1, randomLocation);
        candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), randomLocation), candidateList.end());
        if (!vrptw.reliances.empty())
        {
          for (int relianceLocation : vrptw.reliances[randomLocation])
          {
            newRoute.insert(newRoute.end()-1, relianceLocation);
            candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), relianceLocation), candidateList.end());
          }
        }
        if (!isRouteFeasible(newRoute))
        {
          for (int newRouteLocation : newRoute)
          {
            if (newRouteLocation != 0)
            {
              candidateList.push_back(newRouteLocation);
            }
          }
          newRoute = {0,0};
        }
      }
 
      routes.push_back(newRoute);
    }
  }
}

bool VRPTWDecisionDiagram::generateHeuristicRoutesLiterature(std::vector<std::vector<int>>& routes)
{
  while (true)
  {
    // initialize candidate list
    std::vector<int> candidateList;
    for (int location=1; location<vrptw.numLocations; ++location)
    {
      candidateList.push_back(location);
    }

    // initialize empty routes
    routes.clear();
    for (int vehicleIndex=0; vehicleIndex<vrptw.numVehicles; ++vehicleIndex)
    {
      std::vector<int> emptyRoute;
      emptyRoute.push_back(0);
      emptyRoute.push_back(0);
      routes.push_back(emptyRoute);
    }

    // put random location in each route
    for (int routeIndex=0; routeIndex<routes.size(); ++routeIndex)
    {
      std::vector<int> newRoute = routes[routeIndex];
      while (newRoute.size() == 2)
      {
        int randomInteger0 = std::rand();
        int randomCandidateListIndex = randomInteger0 % candidateList.size();
        int randomLocation = candidateList[randomCandidateListIndex];
        newRoute.insert(newRoute.begin()+1, randomLocation);
        candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), randomLocation), candidateList.end());
        if (!vrptw.reliances.empty())
        {
          for (int relianceLocation : vrptw.reliances[randomLocation])
          {
            newRoute.insert(newRoute.end()-1, relianceLocation);
            candidateList.erase(std::remove(candidateList.begin(), candidateList.end(), relianceLocation), candidateList.end());
          }
        }
        if (!isRouteFeasible(newRoute))
        {
          for (int newRouteLocation : newRoute)
          {
            if (newRouteLocation != 0)
            {
              candidateList.push_back(newRouteLocation);
            }
          }
          newRoute = {0,0};
        }
      }
      routes[routeIndex] = newRoute;
    }

    // choose random insertion criteria and strategy
    int randomInteger1 = std::rand();
    bool insertionStrategyBoolean = ((randomInteger1 % 2) == 0);

    InsertionCriteria insertionCriteria = InsertionCriteria::MCFIC; 

    // run insertion heuristic
    if (vrptw.problemType == ProblemType::PDP)
    {
      insertionStrategyBoolean = true;

      candidateList.clear();
      for (int location=1; location<vrptw.numLocations; ++location)
      {
        candidateList.push_back(location);
      }

      routes.clear();
    }

    std::random_shuffle(candidateList.begin(), candidateList.end());
    if (insertionStrategyBoolean)
    {
      std::cout << "using sequential insertion strategy" << std::endl;
      sequentialInsertion(routes, candidateList, insertionCriteria);
    }
    else
    {
      std::cout << "using sequential insertion strategy" << std::endl;
      sequentialInsertion(routes, candidateList, insertionCriteria);
    }

    if (vrptw.fixedNumPaths == FixedNumPaths::FLEXIBLE_NUM_PATHS)
    {
      break;
    }
    else
    {
      if (routes.size() == vrptw.numVehicles)
      {
        break;
      }
    }
  }

  for (auto route : routes)
  {
    if (!isRouteFeasible(route))
    {
      for (int loc : route)
      {
        std::cout << loc << std::endl;
      }
      return false;
    }
  }

  for (auto primalRoute : routes)
  {
    std::cout << "initial primalRoute: ";
    for (int loc : primalRoute)
    {
      std::cout << loc << " ";
    }
    std::cout << std::endl;
  }
  return true;
}

void VRPTWDecisionDiagram::generateHeuristicRoutesGreedy(std::vector<std::vector<int>>& routes)
{
  std::set<int> locationsUsed;
  locationsUsed.insert(0);
  while (locationsUsed.size() < vrptw.numLocations)
  {
    std::vector<int> newRoute;
    newRoute.push_back(0);
    newRoute.push_back(0);

    std::vector<int> locationsToTry;
    for (int location=1; location<vrptw.numLocations; ++location)
    {
      if (locationsUsed.find(location) == locationsUsed.end())
      {
        locationsToTry.push_back(location);
      }
    }

    for (auto location : locationsToTry)
    {
      addLocationToRoute(location, newRoute);
    }

    if (newRoute.size() > 2)
    {
      if (isRouteFeasible(newRoute))
      {
        routes.push_back(newRoute);
        for (int loc : newRoute)
        {
          locationsUsed.insert(loc);
        }
      }
    }
  }

  for (auto primalRoute : routes)
  {
    std::cout << "initial primalRoute: ";
    for (int loc : primalRoute)
    {
      std::cout << loc << " ";
    }
    std::cout << std::endl;
  }
}
 
void VRPTWDecisionDiagram::createMaximalSequence(const std::vector<int>& sequence, std::vector<int>& maximalSequence)
{
  std::vector<int> locationsToTry;
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    locationsToTry.push_back(location);
  }
  std::random_shuffle(locationsToTry.begin(), locationsToTry.end());

  std::vector<int> updatedSequence = sequence;
  while (true)
  {
    bool isUpdated = false;
    for (int location : locationsToTry)
    {
      for (int index=1; index<updatedSequence.size(); ++index)
      {
        std::vector<int> testSequence = updatedSequence;
        testSequence.insert(testSequence.begin() + index, location);
        bool isFeasible = isRouteFeasible(testSequence);
        if (isFeasible)
        {
          isUpdated = true;
          updatedSequence = testSequence;
          break;
        }
      }

      if (isUpdated)
      {
        break;
      }
    }

    if (!isUpdated)
    {
      break;
    }
  }

  maximalSequence = updatedSequence;
}

int VRPTWDecisionDiagram::addReverseArc(int forwardArcIndex)
{
  VRPTWArc forwardArc = arcs[forwardArcIndex];
  VRPTWArc reverseArc(forwardArc.toNodeIndex, forwardArc.fromNodeIndex, nodes[forwardArc.fromNodeIndex].state.lastVisited, -forwardArc.distance, -forwardArc.coeff, -forwardArc.cijPi);
  int reverseArcIndex = arcs.size();
  arcs.push_back(reverseArc);

  nodes[forwardArc.toNodeIndex].outArcs.push_back(reverseArcIndex);
  nodes[forwardArc.fromNodeIndex].inArcs.push_back(reverseArcIndex);

  // store for next time
  arcReverseArc[forwardArcIndex] = reverseArcIndex;
  arcReverseArc[reverseArcIndex] = forwardArcIndex;

  return reverseArcIndex;
}

// create next relaxed state from current one
bool VRPTWDecisionDiagram::generateNewStateRelaxation(VRPTWNodeState& newState, int location, int nodeIndex)
{
  newState.isExact = false;

  newState.counter = newState.counter + counterBinary;
  if (newState.counter > vrptw.routeLengthUpperBound)
  {
    return false;
  }

  newState.load = newState.load + vrptw.demands[location];
  if (newState.load > vrptw.capacity)
  {
    return false;
  }

  if (std::find(newState.visited.begin(), newState.visited.end(), location) != newState.visited.end())
  {
    return false;
  }

  // precedence only considered for root node, otherwise need to learn by separation
  if ((nodeIndex == 0) && !vrptw.precedences[location].empty())
  {
    if (vrptw.precedences[location].size() > 0)
    {
      return false;
    }
  }

  // if next location is a precedence... can't go there
  if (!vrptw.precedences.empty())
  {
    if (vrptw.precedences[newState.lastVisited].find(location) != vrptw.precedences[newState.lastVisited].end())
    {
      return false;
    }
  }

  // check if we can make it to next location in time
  double earliestStartTime = 0.0;
  if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
  {
    int lastVisitedLocation = newState.lastVisited;
    earliestStartTime = (newState.timeWithMultiplier / (vrptw.timeStateMultiplier * 1.0)) + vrptw.distances[lastVisitedLocation][location] + vrptw.serviceTimes[lastVisitedLocation];
    if (earliestStartTime > vrptw.endTimes[location])
    {
      return false;
    }
  }

  newState.lastVisited = location;

  // check if we can make it back in time
  if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
  {
    if (earliestStartTime + vrptw.distances[location][0] + vrptw.serviceTimes[location] > vrptw.endTimes[0])
    {
      return false;
    }
  }

  // use ng set to get next visited set
  std::set<int> visitedSet;
  for (int loc : newState.visited)
  {
    visitedSet.insert(loc);
  }
  std::set<int> ngSet = ngSets[location];
  std::set<int> newVisited;
  std::set_intersection(visitedSet.begin(), visitedSet.end(), ngSet.begin(), ngSet.end(), std::inserter(newVisited, newVisited.begin()));
  newVisited.insert(location);
  newState.visited = newVisited;

  int newTimeWithMultiplier = (int)(earliestStartTime * vrptw.timeStateMultiplier);

  // testing: do not discretize if the distance + service time is under the step size
  // use more dynamic bucketing if otherwise it'll make a bad loop
  // bucketing
  int timeQuotient = newTimeWithMultiplier / timeStepSize;
  int newTimeWithMultiplierDiscretized = timeStepSize * timeQuotient * timeWindowBinary;
  newTimeWithMultiplierDiscretized = std::max(newTimeWithMultiplierDiscretized, vrptw.startTimes[location] * vrptw.timeStateMultiplier) * timeWindowBinary;
  newState.timeWithMultiplier = newTimeWithMultiplierDiscretized;

  int loadQuotient = newState.load / loadStepSize;
  int newCapacityDiscretized = loadStepSize * loadQuotient * capacityBinary;
  newState.load = newCapacityDiscretized;

  return true;
}

// state should be exact for true/false return to be true
bool VRPTWDecisionDiagram::generateNewStateExact(VRPTWNodeState& newState, int location)
{
  newState.isExact = true;

  if (location == 0)
  {
    std::set<int> initialDeque = {};
    VRPTWNodeState rootNodeState(1,0,0,0,true,initialDeque);
    if (newState == rootNodeState)
    {
      return true;
    }

    if (vrptw.problemType == ProblemType::PDP)
    {
      const auto visitedSet = newState.visited;
      for (int loc : newState.visited)
      {
        if (!vrptw.reliances.empty())
        {
          for (int relianceLoc : vrptw.reliances[loc])
          {
            if (visitedSet.find(relianceLoc) == visitedSet.end())
            {
              return false;
            }
          }
        }
      }
    }
    else
    {
      return true;
    }
  }

  newState.counter = newState.counter + 1;
  if (vrptw.counterType == VRPTWCounterType::NO_USE_COUNTER)
  {
    newState.counter = newState.counter - 1;
  }

  if (newState.counter > vrptw.routeLengthUpperBound + 1)
  {
    return false;
  }

  auto inserted = newState.visited.insert(location);
  if (!inserted.second)
  {
    return false;
  }

  newState.timeWithMultiplier = (int)(((newState.timeWithMultiplier / (vrptw.timeStateMultiplier * 1.0)) + vrptw.serviceTimes[newState.lastVisited] + vrptw.distances[newState.lastVisited][location]) * vrptw.timeStateMultiplier);
  newState.timeWithMultiplier = std::max(newState.timeWithMultiplier, vrptw.startTimes[location] * vrptw.timeStateMultiplier);
  if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::NO_TIME_WINDOWS)
  {
    newState.timeWithMultiplier = 0;
  }
  if ((newState.timeWithMultiplier / (vrptw.timeStateMultiplier * 1.0)) > vrptw.endTimes[location])
  {
    return false;
  }

  newState.lastVisited = location;

  newState.load = newState.load + vrptw.demands[location];
  if ((vrptw.fixedNumPaths == FixedNumPaths::FIXED_NUM_PATHS) && (vrptw.numVehicles == 1))
  {
    newState.load = 0;
  }
  if (newState.load > vrptw.capacity)
  {
    return false;
  }
 
  if (!vrptw.precedences.empty())
  {
    for (int precLoc : vrptw.precedences[location])
    {
      if (newState.visited.find(precLoc) == newState.visited.end())
      {
        return false;
      }
    }
  }

  // check that everything can be dropped off in time
  if (!vrptw.reliances.empty())
  {
    for (int seenLocation : newState.visited)
    {
      auto reliances = vrptw.reliances[seenLocation];
      for (int relianceLocation : reliances)
      {
        if (newState.visited.find(relianceLocation) == newState.visited.end())
        {
          if ((newState.timeWithMultiplier / (vrptw.timeStateMultiplier * 1.0)) + vrptw.distances[location][relianceLocation] > vrptw.endTimes[relianceLocation])
          {
            return false;
          }
        }
      }
    }
  }

  return true;
}

bool VRPTWDecisionDiagram::moveArc(int arcIndex, int newToNodeIndex)
{
  // update arc info
  VRPTWArc& arc = arcs[arcIndex];
  if (arc.toNodeIndex == newToNodeIndex)
  {
    return false;
  }

  nodes[arc.toNodeIndex].inArcs.erase(std::remove(nodes[arc.toNodeIndex].inArcs.begin(), nodes[arc.toNodeIndex].inArcs.end(), arcIndex), nodes[arc.toNodeIndex].inArcs.end());
  arc.toNodeIndex = newToNodeIndex;
  nodes[newToNodeIndex].inArcs.push_back(arcIndex);

  // remove from reverse arc pairs
  int reverseArcIndex = -1;
  auto reverseArcIt = arcReverseArc.find(arcIndex);
  if (reverseArcIt != arcReverseArc.end())
  {
    reverseArcIndex = reverseArcIt->second;
  }

  // update reverse arc map
  if (arcReverseArc.find(arcIndex) != arcReverseArc.end())
  {
    arcReverseArc.erase(arcIndex);
  }
  if (reverseArcIndex != -1)
  {
    if (arcReverseArc.find(reverseArcIndex) != arcReverseArc.end())
    {
      arcReverseArc.erase(reverseArcIndex);
    }
  }
 
  // src cuts
  if (!arcs[arcIndex].isReverseArc)
  {
    if (nodes[arc.toNodeIndex].state.isExact)
    {
      auto visitedLocationsWithArc = nodes[arc.toNodeIndex].state.visited;
      for (int srcIndex=0; srcIndex<srcCuts.size(); ++srcIndex)
      {
        auto cutSet = srcCuts[srcIndex];
        auto srcType = srcCutTypes[srcIndex];
        if (cutSet.find(arc.location) != cutSet.end())
        {
          std::set<int> locationsOverlap;
          for (int location : cutSet)
          {
            if (visitedLocationsWithArc.find(location) != visitedLocationsWithArc.end())
            {
              locationsOverlap.insert(location);
            }
          }

          bool isIncrementalArc = isArcIncrementalSRC(static_cast<int>(locationsOverlap.size()), srcType);
          if (isIncrementalArc)
          {
            if (std::find(srcCutSeparatedArcs[srcIndex].begin(), srcCutSeparatedArcs[srcIndex].end(), arcIndex) == srcCutSeparatedArcs[srcIndex].end())
            {
              srcCutSeparatedArcs[srcIndex].push_back(arcIndex);
              srcCutSeparatedCoeffs[srcIndex].push_back(1);
            }
          }
        }
      }
    }
  }

  return true;
}

void VRPTWDecisionDiagram::removeArc(int arcIndex)
{
  const VRPTWArc& arc = arcs[arcIndex];
  auto fromNodeOutArcIt = std::find(nodes[arc.fromNodeIndex].outArcs.begin(), nodes[arc.fromNodeIndex].outArcs.end(), arcIndex);
  if (fromNodeOutArcIt != std::end(nodes[arc.fromNodeIndex].outArcs))
  {
    nodes[arc.fromNodeIndex].outArcs.erase(fromNodeOutArcIt);
  }
  auto toNodeInArcIt = std::find(nodes[arc.toNodeIndex].inArcs.begin(), nodes[arc.toNodeIndex].inArcs.end(), arcIndex);
  if (toNodeInArcIt != std::end(nodes[arc.toNodeIndex].inArcs))
  {
    nodes[arc.toNodeIndex].inArcs.erase(toNodeInArcIt);
  }

  removedArcs.insert(arcIndex);
}

bool sort_by_second(const std::pair<int,int>& a, const std::pair<int,int>& b)
{
  return (a.second < b.second);
}

void VRPTWDecisionDiagram::setupNgSets(int s)
{
  std::set<int> emptySet;
  ngSets.push_back(emptySet);
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    std::vector<std::pair<int,double>> indexDistances;
    for (int otherLoc=1; otherLoc<vrptw.numLocations; ++otherLoc)
    {
      indexDistances.push_back(std::make_pair(otherLoc,vrptw.distances[location][otherLoc]));
    }
    std::sort(indexDistances.begin(), indexDistances.end(), sort_by_second);

    std::set<int> ngSet;
    ngSet.insert(location);
    for (int index=0; index<s; ++index)
    {
      ngSet.insert(indexDistances[index].first);
    }

    ngSets.push_back(ngSet);
  }
}

void VRPTWDecisionDiagram::compileEmpty()
{
  // reserve but do not resize
  if ((vrptw.fixedNumPaths != FIXED_NUM_PATHS) || (vrptw.numVehicles > 1))
  {
    nodes.reserve(vrptw.capacity * (std::pow(vrptw.numLocations,2)));
  }
  else
  {
    nodes.reserve(std::pow(vrptw.numLocations,2));
  }
  arcs.reserve(std::pow(vrptw.numLocations,2)/100);

  // root node r
  std::set<int> initialDeque = {};
  VRPTWNodeState rootNodeState(1,0,0,0,true,initialDeque);
  addNode(rootNodeState);
  rootNodeIndex = 0;

  // terminal node t
  int terminalNodeLoc = 0;
  int largestLocationIndex = vrptw.numLocations;
  if (vrptw.circuitOrPath == CircuitOrPath::PATH)
  {
    terminalNodeLoc = vrptw.numLocations-1;
    largestLocationIndex = vrptw.numLocations-1;
  }
  double terminalEndTime = 0;
  if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
  {
    terminalEndTime = vrptw.endTimes[0]*(vrptw.timeStateMultiplier);
  }
  VRPTWNodeState terminalNodeState(vrptw.numLocations+2,vrptw.capacity+1,terminalEndTime+1,terminalNodeLoc,false,initialDeque);
  addNode(terminalNodeState);
  terminalNodeIndex = 1;
}

void VRPTWDecisionDiagram::compileNgRoute(int s)
{
  std::cout << "time step mult: " << vrptw.timeStateMultiplier << std::endl;
  std::cout << "time step: " << timeStepSize << std::endl;
  std::cout << "cap step: " << loadStepSize << std::endl;

  // reserve but do not resize
  if ((vrptw.fixedNumPaths != FIXED_NUM_PATHS) || (vrptw.numVehicles > 1))
  {
    nodes.reserve(std::pow(vrptw.numLocations,2));
  }
  else
  {
    nodes.reserve(std::pow(vrptw.numLocations,2));
  }
  arcs.reserve(std::pow(vrptw.numLocations,2)/100);

  // setup ng sets
  setupNgSets(s);

  // setup queue as dictionary
  std::map<int,std::vector<int>> priorityQueue;

  // root node r
  std::set<int> initialDeque = {};
  VRPTWNodeState rootNodeState(1,0,0,0,true,initialDeque);
  addNode(rootNodeState);
  rootNodeIndex = 0;
  priorityQueue.insert(std::make_pair(0,std::vector<int>(1,0)));
  compilationAllVisitedDown.push_back(initialDeque);

  // terminal node t
  int terminalNodeLoc = 0;
  int largestLocationIndex = vrptw.numLocations;
  if (vrptw.circuitOrPath == CircuitOrPath::PATH)
  {
    terminalNodeLoc = vrptw.numLocations-1;
    largestLocationIndex = vrptw.numLocations-1;
  }
  double terminalEndTime = 0;
  if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
  {
    terminalEndTime = vrptw.endTimes[0]*(vrptw.timeStateMultiplier);
  }
  VRPTWNodeState terminalNodeState(vrptw.numLocations+2,vrptw.capacity+1,terminalEndTime+1,terminalNodeLoc,false,initialDeque);
  addNode(terminalNodeState);
  terminalNodeIndex = 1;

  // Use layers for reasoning about compilationAllVisitedDown
  int maxPriorityValue = 0;
  if (vrptw.counterType == VRPTWCounterType::USE_COUNTER)
  {
    priorityQueue.insert(std::make_pair(vrptw.routeLengthUpperBound+1,std::vector<int>(1,1)));
    maxPriorityValue = vrptw.routeLengthUpperBound+1;
  }
  else if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
  {
    priorityQueue.insert(std::make_pair(terminalEndTime,std::vector<int>(1,1)));
    maxPriorityValue = terminalEndTime;
  }
  else
  {
    priorityQueue.insert(std::make_pair(vrptw.capacity+1,std::vector<int>(1,1)));
    maxPriorityValue = vrptw.capacity + 1;
  }
  compilationAllVisitedDown.push_back(initialDeque);

  // create nodes except r/t
  int priorityIteratorValue = -1;
  while (priorityIteratorValue <= maxPriorityValue) 
  {
    ++priorityIteratorValue;
    if (priorityQueue.find(priorityIteratorValue) == priorityQueue.end())
    {
      continue;
    }

    std::vector<int> priorityItems = priorityQueue.find(priorityIteratorValue)->second;
    for (int nodeIndex : priorityItems)
    {
      const VRPTWNodeState stateToCheck = nodes[nodeIndex].state;

      const auto allVisited = compilationAllVisitedDown[nodeIndex];
      for (int location=1; location<largestLocationIndex; ++location)
      {
        if (allVisited.find(location) != allVisited.end())
        {
          continue;
        }

        VRPTWNodeState newState(stateToCheck);
        bool newNodeFeasible = generateNewStateRelaxation(newState, location, nodeIndex);
        if (newNodeFeasible)
        {
          int numNodes = static_cast<int>(nodes.size());
          int newNodeIndex = addNode(newState);
          int newNumNodes = static_cast<int>(nodes.size());
          if (newNumNodes > numNodes)
          {
            int priorityQueueKeyValue = -1;
            if (vrptw.counterType == VRPTWCounterType::USE_COUNTER)
            {
              priorityQueueKeyValue = newState.counter;
            }
            else if (vrptw.vrptwCapacityType == VRPTWCapacityType::NO_RELAX_CAPACITY)
            {
              priorityQueueKeyValue = newState.load;
            }
            else if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
            {
              priorityQueueKeyValue = newState.timeWithMultiplier;
            }

            if (priorityQueue.find(priorityQueueKeyValue) != priorityQueue.end())
            {
              priorityQueue[priorityQueueKeyValue].push_back(newNodeIndex);
            }
            else
            {
              const auto newElement = std::make_pair(priorityQueueKeyValue, std::vector<int>(1,newNodeIndex));
              priorityQueue.insert(newElement);
            }

            compilationAllVisitedDown.push_back(compilationAllVisitedDown[nodeIndex]);
            compilationAllVisitedDown[newNodeIndex].insert(location);
          }
          else
          {
            std::set<int> intersection;
            std::set<int> allDownCurrNode = compilationAllVisitedDown[nodeIndex];
            allDownCurrNode.insert(location);
            std::set_intersection(compilationAllVisitedDown[newNodeIndex].begin(), compilationAllVisitedDown[newNodeIndex].end(), allDownCurrNode.begin(), allDownCurrNode.end(), std::inserter(intersection,intersection.begin()));
            compilationAllVisitedDown[newNodeIndex] = intersection;
          }

          addArc(nodeIndex, newNodeIndex);
        }
      }
    }
  }

  // Add arcs to terminal node
  int nodeIndex = 2;
  while (nodeIndex < nodes.size())
  {
    if (vrptw.circuitOrPath == CircuitOrPath::PATH)
    {
      if (nodes[nodeIndex].state.counter == vrptw.numLocations - 1)
      {
        addArc(nodeIndex, terminalNodeIndex);
      }
    }
    else
    {
      addArc(nodeIndex, terminalNodeIndex);
    }
    nodeIndex = nodeIndex + 1;
  }

  std::cout << "num nodes: " << nodes.size() << std::endl;
  std::cout << "num arcs: " << arcs.size() << std::endl;
};

void VRPTWDecisionDiagram::setCoeffsAsDistances()
{
  for (VRPTWArc& arc : arcs)
  {
    arc.coeff = arc.distance;
    arc.cijPi = arc.distance;
  }
};

void VRPTWDecisionDiagram::setCoeffsAsPreciseDistances()
{
  for (VRPTWArc& arc : arcs)
  {
    double preciseDistance = vrptw.preciseDistances[arc.location][nodes[arc.fromNodeIndex].state.lastVisited];
    arc.coeff = preciseDistance;
    arc.cijPi = preciseDistance;
  }
};

void VRPTWDecisionDiagram::setCoeffsAsDistancesMinusLagrangean(const std::vector<double>& lambda)
{
  for (VRPTWArc& arc : arcs)
  {
    arc.coeff = arc.distance - lambda[arc.location];
    arc.cijPi = arc.distance - lambda[arc.location];
  }
}

void VRPTWDecisionDiagram::setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(const Dual& dual, LPSolveType solveType)
{
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    VRPTWArc& arc = arcs[arcIndex];
    arc.coeff = arc.distance - dual.lambda[arc.location];
    arc.cijPi = arc.distance - dual.lambda[arc.location];

    if (arc.fromNodeIndex == rootNodeIndex)
    {
      arc.coeff = arc.coeff + vrptw.fixedRouteCost;
      arc.cijPi = arc.cijPi + vrptw.fixedRouteCost;
    }
  }

  // use duals for size as we may have added more cuts, but are using previous solution
  for (int capCutIndex=0; capCutIndex<dual.capDuals.size(); ++capCutIndex)
  {
    for (int index=0; index<capCutSetArcs[capCutIndex].size(); ++index)
    {
      int arcIndex = capCutSetArcs[capCutIndex][index];
      double coeff = capCutSetCoeffs[capCutIndex][index];
      VRPTWArc& arc = arcs[arcIndex];
      if (solveType == LPSolveType::LPSolver)
      {
        arc.coeff = arc.coeff - dual.capDuals[capCutIndex] * coeff;
        arc.cijPi = arc.cijPi - dual.capDuals[capCutIndex] * coeff;
      }
      else
      {
        arc.coeff = arc.coeff + dual.capDuals[capCutIndex] * coeff;
        arc.cijPi = arc.cijPi + dual.capDuals[capCutIndex] * coeff;
      }
    }
  }

  // use duals for combs
  for (int combIndex=0; combIndex<dual.combDuals.size(); ++combIndex)
  {
    for (int arcIndex : combCutArcs[combIndex])
    {
      VRPTWArc& arc = arcs[arcIndex];
      arc.coeff = arc.coeff - dual.combDuals[combIndex];
      arc.cijPi = arc.cijPi - dual.combDuals[combIndex];
    }
  }

  // use duals for src
  for (int srcCutIndex=0; srcCutIndex<dual.srcDuals.size(); ++srcCutIndex)
  {
    auto currSrcCutSeparatedArcs = srcCutSeparatedArcs[srcCutIndex];
    auto currSrcCutSeparatedCoeffs = srcCutSeparatedCoeffs[srcCutIndex];
    for (int cutIndex=0; cutIndex<currSrcCutSeparatedArcs.size(); ++cutIndex)
    {
      int arcIndex = currSrcCutSeparatedArcs[cutIndex];
      int coeff = currSrcCutSeparatedCoeffs[cutIndex];

      VRPTWArc& arc = arcs[arcIndex];
      if (solveType == LPSolveType::LPSolver)
      {
        arc.coeff = arc.coeff - dual.srcDuals[srcCutIndex]*coeff;
        arc.cijPi = arc.cijPi - dual.srcDuals[srcCutIndex]*coeff;
      }
      else
      {
        arc.coeff = arc.coeff + dual.srcDuals[srcCutIndex]*coeff;
        arc.cijPi = arc.cijPi + dual.srcDuals[srcCutIndex]*coeff;
      }
    }
  }
}

void VRPTWDecisionDiagram::getNumberOfTimesLocationsCovered(std::unordered_map<int,double>& locationsCovered)
{
  for (const VRPTWArc& arc : arcs)
  {
    locationsCovered[arc.location] += arc.decompositionFlow;
  }
};

void VRPTWDecisionDiagram::getNumberOfTimesLocationsCoveredRoutes(const Primal& primal, std::unordered_map<int,double>& locationsCovered)
{
  for (int index=0; index<primal.xDecompositionFlows.size(); ++index)
  {
    double flow = primal.xDecompositionFlows[index];
    auto route = primal.xDecompositions[index];
    for (int loc : route)
    {
      locationsCovered[loc] += flow;
    }
  }
};

void VRPTWDecisionDiagram::getCutSetValues(std::vector<double>& cutValues)
{
  cutValues.clear();
  for (int capCutIndex=0; capCutIndex<capCutSets.size(); ++capCutIndex)
  {
    double cutSetValue = 0.0;
    for (int index=0; index<capCutSetArcs[capCutIndex].size(); ++index)
    {
      int arcIndex = capCutSetArcs[capCutIndex][index];
      double coeff = capCutSetCoeffs[capCutIndex][index];
      cutSetValue += arcs[arcIndex].heuristicFlow * coeff;
    }
    cutValues.push_back(cutSetValue);
  }
};

void VRPTWDecisionDiagram::getCutSetValuesRoutes(const Primal& primal, std::vector<double>& cutValues)
{
  cutValues.clear();
  for (int capCutIndex=0; capCutIndex<capCutSets.size(); ++capCutIndex)
  {
    double cutSetValue = 0.0;

    for (int index=0; index<primal.xDecompositionFlows.size(); ++index)
    {
      double flow = primal.xDecompositionFlows[index];
      auto route = primal.xDecompositions[index];
      double coeff = getRCCCoeff(route, capCutIndex);
      cutSetValue += (flow * coeff);
    }

    cutValues.push_back(cutSetValue);
  }
};

void VRPTWDecisionDiagram::getCombValues(std::vector<double>& combValues)
{
  combValues.clear();
  for (int combIndex=0; combIndex<teeths.size(); ++combIndex)
  {
    double combValue = 0.0;
    for (int arcIndex : combCutArcs[combIndex])
    {
      combValue += arcs[arcIndex].heuristicFlow;
    }
    combValues.push_back(combValue);
  }
};

void VRPTWDecisionDiagram::getCombValuesRoutes(const Primal& primal, std::vector<double>& combValues)
{
  combValues.clear();
  for (int combIndex=0; combIndex<teeths.size(); ++combIndex)
  {
    double combValue = 0.0;

    auto teeth = teeths[combIndex];
    for (auto tooth : teeth)
    {
      for (int index=0; index<primal.xDecompositionFlows.size(); ++index)
      {
        double flow = primal.xDecompositionFlows[index];
        auto route = primal.xDecompositions[index];

        int fromLoc = route[0];
        int coeff = 0;
        for (int routeIndex=1; routeIndex<route.size(); ++routeIndex)
        {
          int toLoc = route[routeIndex];
          bool fromLocInSet = (std::find(tooth.begin(), tooth.end(), fromLoc) != tooth.end());
          bool toLocInSet = (std::find(tooth.begin(), tooth.end(), toLoc) != tooth.end());
          if ((fromLocInSet && !toLocInSet) || (!fromLocInSet && toLocInSet))
          {
            coeff += 1;
          }

          fromLoc = toLoc;
        }

        combValue += (flow * coeff);
      }
    }

    combValues.push_back(combValue);
  }
};

void VRPTWDecisionDiagram::getSrcCutValues(std::vector<double>& srcCutValues)
{
  srcCutValues.clear();
  for (int srcCutIndex=0; srcCutIndex<srcCuts.size(); ++srcCutIndex)
  {
    // get all disjoint path nodes possible greedily
    double value = 0.0;
    auto currSrcCutSeparatedArcs = srcCutSeparatedArcs[srcCutIndex];
    auto currSrcCutSeparatedCoeffs = srcCutSeparatedCoeffs[srcCutIndex];
    for (int cutIndex=0; cutIndex<currSrcCutSeparatedArcs.size(); ++cutIndex)
    {
      int arcIndex = currSrcCutSeparatedArcs[cutIndex];
      int coeff = currSrcCutSeparatedCoeffs[cutIndex];
      value += arcs[arcIndex].heuristicFlow * coeff;
    }

    srcCutValues.push_back(value);
  }
}

void VRPTWDecisionDiagram::getSrcCutValuesRoutes(const Primal& primal, std::vector<double>& srcCutValues)
{
  srcCutValues.clear();
  for (int srcCutIndex=0; srcCutIndex<srcCuts.size(); ++srcCutIndex)
  {
    auto cutSet = srcCuts[srcCutIndex];
    auto srcType = srcCutTypes[srcCutIndex];
    auto currSrcCutSeparatedArcs = srcCutSeparatedArcs[srcCutIndex];
    auto currSrcCutSeparatedCoeffs = srcCutSeparatedCoeffs[srcCutIndex];

    double feasibleRoutesFlow = 0.0;
    for (int index=0; index<primal.xDecompositionFlows.size(); ++index)
    {
      auto routeArcs = primal.xDecompositionArcs[index];
      double flow = primal.xDecompositionFlows[index];
      for (int cutIndex=0; cutIndex<currSrcCutSeparatedArcs.size(); ++cutIndex)
      {
        int currCutArcIndex = currSrcCutSeparatedArcs[cutIndex];
        if (std::find(routeArcs.begin(), routeArcs.end(), currCutArcIndex) != routeArcs.end())
        {
          int coeff = currSrcCutSeparatedCoeffs[cutIndex];
          feasibleRoutesFlow += flow * coeff;
        }
      }
    }

    srcCutValues.push_back(feasibleRoutesFlow);
  }
}

int VRPTWDecisionDiagram::getSRCRHS(SRCType srcType)
{
  if ((srcType == SRCType::SRC3) || (srcType == SRCType::SRC5V1))
  {
    return 1;
  }
  else if ((srcType == SRCType::SRC4) || (srcType == SRCType::SRC5V2))
  {
    return 2;
  }

  return 0;
};

double VRPTWDecisionDiagram::getSRCFlowViolation(const Primal& primal, const std::set<int>& cutSet, SRCType srcType)
{
  double feasibleRoutesFlow = 0.0;
  for (int index=0; index<primal.xDecompositions.size(); ++index)
  {
    auto route = primal.xDecompositions[index];
    int numTimesVisited = getNumTimesSetVisited(route, cutSet);
    int srcCoeff = getSRCCoeff(numTimesVisited, srcType);
    if (srcCoeff > 0)
    {
      double coeffTimesFlow = srcCoeff * primal.xDecompositionFlows[index];
      feasibleRoutesFlow = feasibleRoutesFlow + coeffTimesFlow;
    }
  }

  int rhs = getSRCRHS(srcType);
  double violation = feasibleRoutesFlow - rhs;

  return violation;
};

void VRPTWDecisionDiagram::trimRouteToLocations(std::vector<int>& routeArcs, const std::set<int>& cutSet)
{
  int count = 0;
  for (int index=0; index<routeArcs.size(); ++index)
  {
    if (cutSet.find(arcs[routeArcs[index]].location) != cutSet.end())
    {
      ++count;
    }

    if (count == cutSet.size())
    {
      routeArcs.resize(index+1);
      return;
    }
  }
}

void VRPTWDecisionDiagram::getRouteSRCArcs(const std::vector<int>& routeArcs, const std::set<int>& cutSet, SRCType srcType, std::set<int>& arcSet)
{
  std::set<int> locations1;
  for (int index=0; index<routeArcs.size(); ++index)
  {
    int arcIndex = routeArcs[index];
    int location = arcs[arcIndex].location;
    if (cutSet.find(location) != cutSet.end())
    {
      locations1.insert(location);
      bool isIncrementalArc = isArcIncrementalSRC(locations1.size(), srcType);
      if (isIncrementalArc)
      {
        arcSet.insert(arcIndex);
      }
    }
  }
}

void VRPTWDecisionDiagram::getSRCArcsAndCoeffs(const Primal& primal, const std::set<int>& cutSet, SRCType srcType, std::vector<int>& srcArcs, std::vector<int>& srcCoeffs)
{
  for (int index=0; index<primal.xDecompositionArcs.size(); ++index)
  {
    auto routeArcs = primal.xDecompositionArcs[index];
    auto route = primal.xDecompositions[index];

    int numTimesVisited = getNumTimesSetVisited(route, cutSet);
    int srcCoeff = getSRCCoeff(numTimesVisited, srcType);
    if (srcCoeff > 0)
    {
      std::set<int> srcArcsRoute;
      getRouteSRCArcs(routeArcs, cutSet, srcType, srcArcsRoute);
      for (int arcIndex : srcArcsRoute)
      {
        if (std::find(srcArcs.begin(), srcArcs.end(), arcIndex) == srcArcs.end())
        {
          srcArcs.push_back(arcIndex);
          srcCoeffs.push_back(1);
        }
      }
    }
  }
};

// add to src cuts, might already have one with same test set
bool VRPTWDecisionDiagram::checkExistingSrcCuts(const std::set<int>& testSet, std::vector<int>& srcArcs, std::vector<int>& srcCoeffs, SRCType srcType)
{
  for (int index=0; index<srcCuts.size(); ++index)
  {
    auto existingTestSet = srcCuts[index];
    auto existingSRCType = srcCutTypes[index];
    if ((testSet == existingTestSet) && (existingSRCType == srcType))
    {
      auto existingArcs = srcCutSeparatedArcs[index];
      auto existingCoeffs = srcCutSeparatedCoeffs[index];
      for (int newArcIndex=0; newArcIndex<srcArcs.size(); ++newArcIndex)
      {
        int newArc = srcArcs[newArcIndex];
        int newArcCoeff = srcCoeffs[newArcIndex];
        if (std::find(existingArcs.begin(), existingArcs.end(), newArc) == existingArcs.end())
        {
          srcCutSeparatedArcs[index].push_back(newArc);
          srcCutSeparatedCoeffs[index].push_back(newArcCoeff);
        }
        std::cout << "added arc " << newArc << " to existing src cut" << std::endl;
      }
      return false;
    }
  }

  return true;
};

bool sort_by_second3(const std::pair<std::pair<SRCType,std::set<int>>,double>& a, const std::pair<std::pair<SRCType,std::set<int>>,double>& b)
{
  return (a.second > b.second);
}

int VRPTWDecisionDiagram::findSRCs(const Primal& primal, int limit, std::vector<double>& violations)
{
  double violation = -1;
  std::vector<std::pair<std::pair<SRCType,std::set<int>>,double>> testSetViolations;
  for (int i=1; i<vrptw.numLocations-2; ++i)
  {
    for (int j=i+1; j<vrptw.numLocations-1; ++j)
    {
      for (int k=j+1; k<vrptw.numLocations; ++k)
      {
        // SRC3
        if (params.useSRC3s)
        {
          std::set<int> testSetSRC3 = {i,j,k};
          violation = getSRCFlowViolation(primal, testSetSRC3, SRCType::SRC3);
          if (violation > 0.1)
          {
            testSetViolations.push_back(std::make_pair(std::make_pair(SRCType::SRC3,testSetSRC3),violation));
            continue;
          }
        }

        // SRC4
        if (params.useSRC4s)
        {
          for (int l=k+1; l<vrptw.numLocations; ++l)
          {
            std::set<int> testSetSRC4 = {i,j,k,l};
            violation = getSRCFlowViolation(primal, testSetSRC4, SRCType::SRC4);
            if (violation > 0.1)
            {
              testSetViolations.push_back(std::make_pair(std::make_pair(SRCType::SRC4,testSetSRC4),violation));
              continue;
            }

            // SRC5V1
            if (params.useSRC5V1s)
            {
              for (int m=l+1; m<vrptw.numLocations; ++m)
              {
                std::set<int> testSetSRC5V1 = {i,j,k,l,m};
                violation = getSRCFlowViolation(primal, testSetSRC5V1, SRCType::SRC5V1);
                if (violation > 0.1)
                {
                  testSetViolations.push_back(std::make_pair(std::make_pair(SRCType::SRC5V1,testSetSRC5V1),violation));
                  continue;
                }
              }
            }
          }
        }
      }
    }
  }
 
  // order distances
  std::sort(testSetViolations.begin(), testSetViolations.end(), sort_by_second3);
  testSetViolations.resize(std::min(limit,static_cast<int>(testSetViolations.size())));
  for (auto typeTestSetViolation : testSetViolations)
  {
    auto srcType = typeTestSetViolation.first.first;
    auto testSet = typeTestSetViolation.first.second;
    auto violation = typeTestSetViolation.second;

    std::vector<int> srcArcs;
    std::vector<int> srcCoeffs;
    getSRCArcsAndCoeffs(primal, testSet, srcType, srcArcs, srcCoeffs);
    bool shouldAddNewCut = checkExistingSrcCuts(testSet, srcArcs, srcCoeffs, srcType);
    if (shouldAddNewCut)
    {
      violations.push_back(violation);
      std::vector<int> emptyVector;
      std::cout << "add src cut type: " << srcType << " with arcs: ";
      srcCuts.push_back(testSet);
      srcCutTypes.push_back(srcType);

      // add for all arcs instead of just the ones in the violation
      // would need to ensure that they are exact
      /*
      srcArcs.clear();
      srcCoeffs.clear();
      for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
      {
        if (testSet.find(arcs[arcIndex].location) != testSet.end())
        {
          auto visitedLocationsWithArc = nodes[arcs[arcIndex].toNodeIndex].state.visited;
          std::set<int> locationsOverlap;
          for (int location : testSet)
          {
            if (visitedLocationsWithArc.find(location) != visitedLocationsWithArc.end())
            {
              locationsOverlap.insert(location);
            }
          }

          bool isIncrementalArc = isArcIncrementalSRC(locationsOverlap.size(), srcType);
          if (isIncrementalArc)
          {
            srcArcs.push_back(arcIndex);
            srcCoeffs.push_back(1);
          }
        }
      }
      */

      srcCutSeparatedArcs.push_back(srcArcs);
      srcCutSeparatedCoeffs.push_back(srcCoeffs);
      for (int a : srcArcs)
      {
        std::cout << a << " ";
      }

      std::cout << "for test set: ";
      for (int loc : testSet)
      {
        std::cout << loc << " ";
      }
      std::cout << std::endl;

      std::cout << "index: " << srcCuts.size() << " violation: " << violation << std::endl;
      std::cout << "arcs: ";
      for (int arcIndex : srcArcs)
      {
        std::cout << arcIndex << ",";
      }
      std::cout << std::endl;
    }
    else
    {
      std::cout << "cut existed for test set: " << std::endl;
      for (int loc : testSet)
      {
        std::cout << loc << " ";
      }
      std::cout << std::endl;
    }
  }

  return static_cast<int>(testSetViolations.size());
}

double VRPTWDecisionDiagram::computeShortestPathBFS(ShortestPathMode mode, std::vector<int>& routeByLocation)
{
  // setup data structures to store shortestpaths and the path itself
  for (VRPTWNode& node : nodes)
  {
    node.shortestPathDistance = INF;
  }
  nodes[rootNodeIndex].shortestPathDistance = 0;

  std::vector<int> priorNodeForShortestPath;
  std::vector<int> priorLocationForShortestPath;
  std::vector<int> priorArcIndexShortestPath;
  priorNodeForShortestPath.resize(nodes.size());
  priorLocationForShortestPath.resize(nodes.size());
  priorArcIndexShortestPath.resize(nodes.size());

  // dp to find shortest path
  for (auto capNodeIndices : nodeOrdering)
  {
    for (int nodeIndex : capNodeIndices.second)
    {
      for (int arcIndex : nodes[nodeIndex].outArcs)
      {
        int toNodeIndex = arcs[arcIndex].toNodeIndex;

        double distanceFromNode = nodes[nodeIndex].shortestPathDistance + arcs[arcIndex].coeff;
        if (distanceFromNode < nodes[toNodeIndex].shortestPathDistance)
        {
          priorLocationForShortestPath[toNodeIndex] = arcs[arcIndex].location;
          priorNodeForShortestPath[toNodeIndex] = nodeIndex;
          priorArcIndexShortestPath[toNodeIndex] = arcIndex;
          nodes[toNodeIndex].shortestPathDistance = distanceFromNode;
        }
      }
    }
  }

  // construct route
  if (mode == ShortestPathMode::SHORTEST_PATH)
  {
    routeByLocation.push_back(priorLocationForShortestPath[terminalNodeIndex]);
    int currentNodeIndex = priorNodeForShortestPath[terminalNodeIndex];
    int currentArcIndex = priorArcIndexShortestPath[terminalNodeIndex];
    while (true)
    {
      currentArcIndex = priorArcIndexShortestPath[currentNodeIndex];
      routeByLocation.push_back(priorLocationForShortestPath[currentNodeIndex]);
      currentNodeIndex = priorNodeForShortestPath[currentNodeIndex];
      if (currentNodeIndex == rootNodeIndex)
      {
        routeByLocation.push_back(0);
        std::reverse(routeByLocation.begin(), routeByLocation.end());
        break;
      }
    }
  }
 
  // update potentials and reduced costs
  if (mode == ShortestPathMode::UPDATE_POTENTIALS)
  {
    for (VRPTWNode& node : nodes)
    {
      // pi = pi - d;
      node.potential = node.potential - node.shortestPathDistance;
    }

    for (VRPTWNode& node : nodes)
    {
      for (int arcIndex : node.outArcs)
      {
        // c^pi_ij = c_ij + pi(i) - pi(j)
        const VRPTWNode& toNode = nodes[arcs[arcIndex].toNodeIndex];
        arcs[arcIndex].cijPi = arcs[arcIndex].coeff - node.potential + toNode.potential;
      }
    }
  }

  return nodes[terminalNodeIndex].shortestPathDistance;
};

double VRPTWDecisionDiagram::computeShortestPathBFSWang(std::vector<int>& treeByParentArcs, std::vector<int>& routeByArc, double& maxShortestPath)
{
  // setup data structures to store shortestpaths and the path itself
  for (VRPTWNode& node : nodes)
  {
    node.shortestPathDistance = INF;
  }
  nodes[rootNodeIndex].shortestPathDistance = 0.0;

  // keep track to guide arc fixing / repair process
  maxShortestPath = 0;

  std::vector<int> priorNodeForShortestPath;
  std::vector<int> priorLocationForShortestPath;
  priorNodeForShortestPath.resize(nodes.size());
  priorLocationForShortestPath.resize(nodes.size());

  // dp to find shortest path
  for (auto capNodeIndices : nodeOrdering)
  {
    for (int nodeIndex : capNodeIndices.second)
    {
      for (int arcIndex : nodes[nodeIndex].outArcs)
      {
        int toNodeIndex = arcs[arcIndex].toNodeIndex;

        double distanceFromNode = nodes[nodeIndex].shortestPathDistance + arcs[arcIndex].coeff;
        if (distanceFromNode < nodes[toNodeIndex].shortestPathDistance)
        {
          priorLocationForShortestPath[toNodeIndex] = arcs[arcIndex].location;
          priorNodeForShortestPath[toNodeIndex] = nodeIndex;
          treeByParentArcs[toNodeIndex] = arcIndex;
          nodes[toNodeIndex].shortestPathDistance = distanceFromNode;
          if (distanceFromNode < vrptw.instanceUpperBound)
          {
            maxShortestPath = std::max(maxShortestPath, distanceFromNode);
          }
        }
      }
    }
  }

  // construct route
  routeByArc.push_back(treeByParentArcs[terminalNodeIndex]);
  int currentNodeIndex = priorNodeForShortestPath[terminalNodeIndex];
  int currentArcIndex = treeByParentArcs[terminalNodeIndex];
  while (true)
  {
    currentArcIndex = treeByParentArcs[currentNodeIndex];
    routeByArc.push_back(currentArcIndex);
    currentNodeIndex = priorNodeForShortestPath[currentNodeIndex];
    if (currentNodeIndex == rootNodeIndex)
    {
      std::reverse(routeByArc.begin(), routeByArc.end());
      break;
    }
  }

  return nodes[terminalNodeIndex].shortestPathDistance;
};

double VRPTWDecisionDiagram::setupAndSolveFlowModel(FlowType flowType, IncludeCoverConstraints includeCoverConstraints, UseColumnGeneration useCg, const std::set<int>& initialPrimalArcIndices, Dual& duals, bool removeForTesting, int timeoutSeconds)
{
  // Setup model
  IloEnv env;
  IloModel flowModel(env);
  IloRangeArray coverConstraints(env);
  IloRangeArray flowConservationConstraints(env);
  IloRangeArray fixedPathConstraint(env);
  IloRangeArray capConstraints(env);
  IloRangeArray combConstraints(env);
  IloRangeArray srcConstraints(env);
  IloNumVarArray x(env, arcs.size());
  IloNumArray initialPrimal(env);
  IloExpr objective(env);

  // Setup variables and objective
  int numNonZeroVars = 0;
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    if (!arcs[arcIndex].isReverseArc && isArcAlive(arcIndex))
    {
      if (flowType == FlowType::IP)
      {
        // NOTE(akarahal) try seeing if upper bounds get dual values
        //x[arcIndex] = IloNumVar(env, 0, IloInfinity, ILOINT);
        x[arcIndex] = IloNumVar(env, 0, 1, ILOINT);

        // PDP minimize number of vehicles used as first objective
        // Additional term in objective for smaller routes
        if (vrptw.problemType == ProblemType::PDP)
        {
          objective += x[arcIndex] * arcs[arcIndex].coeff;
          if (arcs[arcIndex].fromNodeIndex == rootNodeIndex)
          {
            objective += x[arcIndex] * vrptw.fixedRouteCost;
          }
          /*
          else if (arcs[arcIndex].toNodeIndex == terminalNodeIndex)
          {
            objective -= x[arcIndex] * std::pow(nodes[arcs[arcIndex].fromNodeIndex].state.counter,2);
          }
          */
        }
        else
        {
          objective += x[arcIndex] * arcs[arcIndex].coeff;
        }

        // Initial solution values
        if (initialPrimalArcIndices.find(arcIndex) != initialPrimalArcIndices.end())
        {
          initialPrimal.add(1);
        }
        else
        {
          initialPrimal.add(0);
        }

        ++numNonZeroVars;
      }
      else
      {
        if (useCg == UseColumnGeneration::USE_CG)
        {
          x[arcIndex] = IloNumVar(env, 0, 0);
        }
        else
        {
          ++numNonZeroVars;
          x[arcIndex] = IloNumVar(env, 0, 1);
          objective += x[arcIndex] * arcs[arcIndex].coeff;
        }
      }
    }
    else
    {
      x[arcIndex] = IloNumVar(env, 0, 0);
    }
  }

  std::cout << "num vars: " << numNonZeroVars << std::endl;
  flowModel.add(IloMinimize(env, objective));
  objective.end();

  // Cover Constraints
  if (includeCoverConstraints == IncludeCoverConstraints::Y)
  {
    for (int location=0; location<vrptw.numLocations; ++location)
    {
      IloExpr sumLocationArcs(env);
      for (int arcIndex : locationToArcs[location])
      {
        if (!arcs[arcIndex].isReverseArc)
        {
          sumLocationArcs += x[arcIndex];
        }
      }

      if ((location == 0) && (vrptw.circuitOrPath == CircuitOrPath::PATH))
      {
        coverConstraints.add(x[0] >= 0);
      }
      else
      {
        coverConstraints.add(sumLocationArcs >= 1);
      }
    }
    flowModel.add(coverConstraints);
  }

  // Flow Conservation Constraints
  for (int nodeIndex=2; nodeIndex<nodes.size(); ++nodeIndex)
  {
    IloExpr sumInMinusOutNode(env);
    for (int arcIndex : nodes[nodeIndex].outArcs)
    {
      sumInMinusOutNode += x[arcIndex];
    }

    for (int arcIndex : nodes[nodeIndex].inArcs)
    {
      sumInMinusOutNode -= x[arcIndex];
    }

    flowConservationConstraints.add(sumInMinusOutNode == 0);
  }
  flowModel.add(flowConservationConstraints);

  // Fixed Number of Vehicles
  if (vrptw.fixedNumPaths == FIXED_NUM_PATHS)
  {
    IloExpr fixedNumPaths(env);
    for (auto arcIndex : nodes[rootNodeIndex].outArcs)
    {
      fixedNumPaths += x[arcIndex];
    }

    fixedPathConstraint.add(fixedNumPaths == vrptw.numVehicles);
    flowModel.add(fixedPathConstraint);
  }

  // Cuts for LP only
  if (flowType == FlowType::LP)
  {
    if (includeCoverConstraints == IncludeCoverConstraints::Y)
    {
      // RCC - rounded capacity cuts, might even need to use for IPs
      for (int capCutIndex=0; capCutIndex<capCutSets.size(); ++capCutIndex)
      {
        if (!removeForTesting)
        {
          IloExpr cutSetSum(env);
          for (int index=0; index<capCutSetArcs[capCutIndex].size(); ++index)
          {
            int arcIndex = capCutSetArcs[capCutIndex][index];
            double coeff = capCutSetCoeffs[capCutIndex][index];
            cutSetSum += x[arcIndex] * coeff;
          }

          capConstraints.add(cutSetSum <= capCutSetsRHS[capCutIndex]);
        }
        else
        {
          capConstraints.add(x[0] <= 1);
        }
      }
      flowModel.add(capConstraints);

      // Comb - strengthened comb cuts
      for (int combIndex=0; combIndex<teeths.size(); ++combIndex)
      {
        IloExpr combSum(env);
        for (int arcIndex : combCutArcs[combIndex])
        {
          combSum += x[arcIndex];
        }
     
        combConstraints.add(combSum >= combRHS[combIndex]);
      }
      flowModel.add(combConstraints);

      // SRC - subset row cuts
      for (int srcCutIndex=0; srcCutIndex<srcCuts.size(); ++srcCutIndex)
      {
        IloExpr srcCutSum(env);
        auto currSrcCutSeparatedArcs = srcCutSeparatedArcs[srcCutIndex];
        auto currSrcCutSeparatedCoeffs = srcCutSeparatedCoeffs[srcCutIndex];
        for (int cutIndex=0; cutIndex<currSrcCutSeparatedArcs.size(); ++cutIndex)
        {
          int arcIndex = currSrcCutSeparatedArcs[cutIndex];
          int coeff = currSrcCutSeparatedCoeffs[cutIndex];
          srcCutSum += x[arcIndex] * coeff;
        }

        SRCType srcType = srcCutTypes[srcCutIndex];
        int rhs = getSRCRHS(srcType);
        srcConstraints.add(srcCutSum <= rhs);
      }
      flowModel.add(srcConstraints);
    }
  }

  // Solver
  IloCplex solver(flowModel);
  solver.setOut(env.getNullStream());
  solver.setWarning(env.getNullStream());
  solver.setError(env.getNullStream());
  solver.setParam(IloCplex::Param::TimeLimit, std::max(5, timeoutSeconds));
  //solver.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Barrier);
  //solver.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
  solver.setParam(IloCplex::Param::Threads, 1);
  solver.exportModel("LPflowmodel.lp");

  // Warm starts
  if ((flowType == FlowType::LP) && (includeCoverConstraints == IncludeCoverConstraints::Y))
  {
    IloNumArray startDuals(env);
    for (auto dual : duals.lambda)
    {
      startDuals.add(dual);
    }
    solver.setStart(NULL, NULL, x, NULL, startDuals, coverConstraints);
  }
  else
  {
    if (!initialPrimalArcIndices.empty())
    {
      solver.addMIPStart(x, initialPrimal);
    }
  }
  solver.solve();

  // Get results
  IloAlgorithm::Status solverStatus = solver.getStatus();
  if ((solverStatus == IloAlgorithm::Optimal) || ((solverStatus == IloAlgorithm::Feasible) && (!initialPrimalArcIndices.empty())))
  {
    for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
    {
      if (!arcs[arcIndex].isReverseArc && isArcAlive(arcIndex))
      {
        arcs[arcIndex].decompositionFlow = solver.getValue(x[arcIndex]);
        arcs[arcIndex].heuristicFlow= solver.getValue(x[arcIndex]);
      }
      /*
      if (arcs[arcIndex].decompositionFlow > 0)
      {
        std::cout << "primal [" << arcIndex << "]: " << arcs[arcIndex].decompositionFlow << std::endl;
      }
      */
    }

    // Get duals for LP
    double objValue = solver.getObjValue();
    if (flowType == FlowType::LP)
    {
      double dualValue = 0.0;
      if (includeCoverConstraints == IncludeCoverConstraints::Y)
      {
        IloNumArray lpDuals(env);
        solver.getDuals(lpDuals, coverConstraints);
        for (int dualIndex=0; dualIndex<vrptw.numLocations; ++dualIndex)
        {
          duals.lambda[dualIndex] = lpDuals[dualIndex];
          dualValue += lpDuals[dualIndex];
        }
        std::cout << "DUALS: ";
        for (int dualIndex=0; dualIndex<vrptw.numLocations; ++dualIndex)
        {
          std::cout << lpDuals[dualIndex] << ",";
        }
        std::cout << std::endl;
      }

      if (vrptw.fixedNumPaths == FIXED_NUM_PATHS)
      {
        IloNumArray fixedPathDualFromLP(env);
        solver.getDuals(fixedPathDualFromLP, fixedPathConstraint);
        duals.fixedPathDual = fixedPathDualFromLP[0];
        std::cout << "fixed path dual: " << duals.fixedPathDual << std::endl;
      }

      if (includeCoverConstraints == IncludeCoverConstraints::Y)
      {
        IloNumArray lpCapDuals(env);
        solver.getDuals(lpCapDuals, capConstraints);
        duals.capDuals.resize(capCutSets.size());
        for (int dualIndex=0; dualIndex<capCutSets.size(); ++dualIndex)
        {
          duals.capDuals[dualIndex] = lpCapDuals[dualIndex];
          dualValue += (lpCapDuals[dualIndex] * capCutSetsRHS[dualIndex]);
          if (duals.capDuals[dualIndex] < 0)
          {
            std::cout << "cap dual [" << dualIndex << "]: " << duals.capDuals[dualIndex] << std::endl;
          }
        }
   
        IloNumArray lpCombDuals(env);
        solver.getDuals(lpCombDuals, combConstraints);
        duals.combDuals.resize(teeths.size());
        for (int dualIndex=0; dualIndex<teeths.size(); ++dualIndex)
        {
          duals.combDuals[dualIndex] = lpCombDuals[dualIndex];
          dualValue += lpCombDuals[dualIndex] * combRHS[dualIndex];
        }
        for (int dualIndex=0; dualIndex<teeths.size(); ++dualIndex)
        {
          std::cout << "comb dual [" << dualIndex << "]: " << duals.combDuals[dualIndex] << std::endl;
        }

        IloNumArray lpSrcDuals(env);
        solver.getDuals(lpSrcDuals, srcConstraints);
        duals.srcDuals.resize(srcCuts.size());
        for (int dualIndex=0; dualIndex<srcCuts.size(); ++dualIndex)
        {
          duals.srcDuals[dualIndex] = lpSrcDuals[dualIndex];
          if (duals.srcDuals[dualIndex] < 0)
          {
            std::cout << "src dual [" << dualIndex << "]: " << duals.srcDuals[dualIndex] << std::endl;
          }
          dualValue += lpSrcDuals[dualIndex];
        }
      }

      dualValue = dualValue + duals.fixedPathDual*vrptw.numVehicles;
      std::cout << "lp dual value: " << dualValue << std::endl;
      std::cout << "obj val: " << objValue << std::endl;
      if ((dualValue <= objValue - 0.00001) || (objValue <= dualValue - 0.00001))
      {
        std::cout << "WARNING dual values for unit upper bounds on arcs" << std::endl;
      }
    }

    env.end();
    return objValue;
  }
  else
  {
    std::cout << "results not optimal" << std::endl;
    std::cout << solverStatus << std::endl;
    env.end();
    return -1;
  }

  env.end();
  return -1;
};

double VRPTWDecisionDiagram::getDualObjectiveValue(const Dual& dual, LPSolveType solveType)
{
  double lowerBound = 0.0;
  for (int dualIndex=0; dualIndex<vrptw.numLocations; ++dualIndex)
  {
    lowerBound += dual.lambda[dualIndex];
  }

  for (int dualIndex=0; dualIndex<capCutSets.size(); ++dualIndex)
  {
    if (solveType == LPSolveType::LPSolver)
    {
      lowerBound += (dual.capDuals[dualIndex] * capCutSetsRHS[dualIndex]);
    }
    else
    {
      lowerBound += (-1 * dual.capDuals[dualIndex] * capCutSetsRHS[dualIndex]);
    }
  }

  for (int dualIndex=0; dualIndex<teeths.size(); ++dualIndex)
  {
    lowerBound += dual.combDuals[dualIndex] * combRHS[dualIndex];
  }

  for (int dualIndex=0; dualIndex<srcCuts.size(); ++dualIndex)
  {
    SRCType srcType = srcCutTypes[dualIndex];
    int rhs = getSRCRHS(srcType);
    if (solveType == LPSolveType::LPSolver)
    {
      lowerBound += dual.srcDuals[dualIndex] * rhs;
    }
    else
    {
      lowerBound -= dual.srcDuals[dualIndex] * rhs;
    }
  }

  lowerBound = lowerBound + dual.fixedPathDual*vrptw.numVehicles;

  return lowerBound;
}

void VRPTWDecisionDiagram::findMergeNodesReducedCost(const Dual& dual, LPSolveType solveType, double limitToMerge)
{
  std::cout << "merge nodes" << std::endl;

  setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(dual, solveType);

  // check over nodes for merge opportunities
  bool keepChecking = true;
  while (keepChecking)
  {
    std::vector<double> shortestPathDown(nodes.size(), INF);
    shortestPathDown[rootNodeIndex] = 0;
    std::vector<double> shortestPathUp(nodes.size(), INF);
    shortestPathUp[terminalNodeIndex] = 0;
    std::vector<std::set<int>> allVisitedDown(nodes.size());
    std::vector<bool> nodeSeen(nodes.size(), false);

    // dp to find shortest path down
    for (auto capNodeIndices : nodeOrdering)
    {
      for (int nodeIndex : capNodeIndices.second)
      {
        for (int arcIndex : nodes[nodeIndex].outArcs)
        {
          int toNodeIndex = arcs[arcIndex].toNodeIndex;

          double distanceFromNode = shortestPathDown[nodeIndex] + arcs[arcIndex].coeff;
          if (distanceFromNode < shortestPathDown[toNodeIndex])
          {
            shortestPathDown[toNodeIndex] = distanceFromNode;
          }

          int location = arcs[arcIndex].location;
          if (nodeSeen[toNodeIndex])
          {
            std::set<int> intersection;
            std::set<int> allDownCurrNode = allVisitedDown[nodeIndex];
            allDownCurrNode.insert(location);
            std::set_intersection(allVisitedDown[toNodeIndex].begin(), allVisitedDown[toNodeIndex].end(), allDownCurrNode.begin(), allDownCurrNode.end(), std::inserter(intersection,intersection.begin()));
            allVisitedDown[toNodeIndex] = intersection;
          }
          else
          {
            allVisitedDown[toNodeIndex] = allVisitedDown[nodeIndex];
            allVisitedDown[toNodeIndex].insert(location);
            nodeSeen[toNodeIndex] = true;
          }
        }
      }
    }

    // dp to find shortest path up
    for (auto it=nodeOrdering.rbegin(); it!=nodeOrdering.rend(); ++it)
    {
      for (int nodeIndex : it->second)
      {
        for (int arcIndex : nodes[nodeIndex].inArcs)
        {
          int fromNodeIndex = arcs[arcIndex].fromNodeIndex;

          double distanceFromNode = shortestPathUp[nodeIndex] + arcs[arcIndex].coeff;
          if (distanceFromNode < shortestPathUp[fromNodeIndex])
          {
            shortestPathUp[fromNodeIndex] = distanceFromNode;
          }
        }
      }
    }

    // check all pairs for big enough (up + down) and compatible merging
    std::vector<int> nodesToCheck;
    for (int nodeIndex1=0; nodeIndex1<nodes.size(); ++nodeIndex1)
    {
      double upPlusDown1 = shortestPathUp[nodeIndex1] + shortestPathDown[nodeIndex1];
      if ((upPlusDown1 > limitToMerge) && (shortestPathUp[nodeIndex1] < INF) && (shortestPathDown[nodeIndex1] < INF))
      {
        nodesToCheck.push_back(nodeIndex1);
      }
    }
    
    keepChecking = false;
    for (int nodeIndexToCheck1=0; nodeIndexToCheck1<nodesToCheck.size(); ++nodeIndexToCheck1)
    {
      int nodeIndex1 = nodesToCheck[nodeIndexToCheck1];
      for (int nodeIndexToCheck2=nodeIndexToCheck1+1; nodeIndexToCheck2<nodesToCheck.size(); ++nodeIndexToCheck2)
      {
        int nodeIndex2 = nodesToCheck[nodeIndexToCheck2];
        double down1up2 = shortestPathDown[nodeIndex1] + shortestPathUp[nodeIndex2];
        double down2up1 = shortestPathDown[nodeIndex2] + shortestPathUp[nodeIndex1];
        if ((down1up2 > limitToMerge) && (down2up1 > limitToMerge))
        {
          if (canMergeNodes(nodeIndex1, nodeIndex2))
          {
            mergeNodes(nodeIndex1, nodeIndex2);
            std::cout << "merged: " << nodeIndex1 << " and " << nodeIndex2 << std::endl;
            keepChecking = true;
          }
        }

        if (keepChecking)
        {
          break;
        }
      }
 
      if (keepChecking)
      {
        break;
      }
    }
  }
}

bool VRPTWDecisionDiagram::canMergeNodes(int nodeIndex1, int nodeIndex2)
{
  const VRPTWNodeState& state1 = nodes[nodeIndex1].state;
  const VRPTWNodeState& state2 = nodes[nodeIndex2].state;

  if (state1.lastVisited != state2.lastVisited)
  {
    return false;
  }

  // check for DAG to be maintained
  if (vrptw.counterType == VRPTWCounterType::USE_COUNTER)
  {
    if (state1.counter != state2.counter)
    {
      return false;
    }
  }
  else
  {
    for (int arcIndex : nodes[nodeIndex2].inArcs)
    {
      if (nodes[arcs[arcIndex].fromNodeIndex].state.load >= state1.load)
      {
        return false;
      }
    }
  }

  if (state1.load > state2.load)
  {
    return false;
  }
 
  if (state1.timeWithMultiplier > state2.timeWithMultiplier)
  {
    return false;
  }

  if (!std::includes(state2.visited.begin(), state2.visited.end(),
                     state1.visited.begin(), state1.visited.end()))
  {
    return false;
  }

  return true;
}

void VRPTWDecisionDiagram::mergeNodes(int nodeIndex1, int nodeIndex2)
{
  // move all in arcs from nodeIndex2 to nodeIndex1
  for (int arcIndex : nodes[nodeIndex2].inArcs)
  {
    moveArc(arcIndex, nodeIndex1);
  }
}

double VRPTWDecisionDiagram::fixArcs(const std::vector<Dual>& duals, LPSolveType solveType, double upperBound)
{
  for (const Dual& dual : duals)
  {
    fixArcs(dual, solveType, upperBound);
  }

  return getPercentFixedArcs();
};

double VRPTWDecisionDiagram::fixArcs(const std::deque<Dual>& duals, LPSolveType solveType, double upperBound)
{
  for (const Dual& dual : duals)
  {
    fixArcs(dual, solveType, upperBound);
  }

  return getPercentFixedArcs();
};

double VRPTWDecisionDiagram::fixArcs(const Dual& dual, LPSolveType solveType, double upperBound)
{
  // get value of dual
  double lowerBound = getDualObjectiveValue(dual, solveType);
  std::cout << "fixing. lower bound used: " << lowerBound << std::endl;

  setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(dual, solveType);

  std::vector<double> shortestPathDown(nodes.size(), INF);
  shortestPathDown[rootNodeIndex] = 0;
  std::vector<double> shortestPathUp(nodes.size(), INF);
  shortestPathUp[terminalNodeIndex] = 0;

  // dp to find shortest path down
  for (auto capNodeIndices : nodeOrdering)
  {
    for (int nodeIndex : capNodeIndices.second)
    {
      for (int arcIndex : nodes[nodeIndex].outArcs)
      {
        int toNodeIndex = arcs[arcIndex].toNodeIndex;

        double distanceFromNode = shortestPathDown[nodeIndex] + arcs[arcIndex].coeff;
        if (distanceFromNode < shortestPathDown[toNodeIndex])
        {
          shortestPathDown[toNodeIndex] = distanceFromNode;
        }
      }
    }
  }

  // dp to find shortest path up
  for (auto it=nodeOrdering.rbegin(); it!=nodeOrdering.rend(); ++it)
  {
    for (int nodeIndex : it->second)
    {
      for (int arcIndex : nodes[nodeIndex].inArcs)
      {
        int fromNodeIndex = arcs[arcIndex].fromNodeIndex;

        double distanceFromNode = shortestPathUp[nodeIndex] + arcs[arcIndex].coeff;
        if (distanceFromNode < shortestPathUp[fromNodeIndex])
        {
          shortestPathUp[fromNodeIndex] = distanceFromNode;
        }
      }
    }
  }

  // only check arcs that are not removed ... dont loop over arcs.size()
  std::set<double> setToCheck;
  for (int nodeIndex=0; nodeIndex<nodes.size(); ++nodeIndex)
  {
    const std::vector<int> arcIndices = nodes[nodeIndex].outArcs;
    for (int arcIndex : arcIndices)
    {
      // fix arcs based on lb + rc > ub
      const VRPTWArc& arcToCheck = arcs[arcIndex];
      double bestPossibleReducedCost = shortestPathDown[arcToCheck.fromNodeIndex] + shortestPathUp[arcToCheck.toNodeIndex] + arcToCheck.coeff - dual.fixedPathDual;

      bool removeReducedCost = (lowerBound + bestPossibleReducedCost) > (upperBound + 0.00001);
      if (removeReducedCost)
      {
        // remove from graph if there
        bool removed = false;
        VRPTWArc& arc = arcs[arcIndex];
        auto fromNodeOutArcIt = std::find(nodes[arc.fromNodeIndex].outArcs.begin(), nodes[arc.fromNodeIndex].outArcs.end(), arcIndex);
        if (fromNodeOutArcIt != std::end(nodes[arc.fromNodeIndex].outArcs))
        {
          removed = true;
          nodes[arc.fromNodeIndex].outArcs.erase(fromNodeOutArcIt);
        }
        auto toNodeInArcIt = std::find(nodes[arc.toNodeIndex].inArcs.begin(), nodes[arc.toNodeIndex].inArcs.end(), arcIndex);
        if (toNodeInArcIt != std::end(nodes[arc.toNodeIndex].inArcs))
        {
          removed = true;
          nodes[arc.toNodeIndex].inArcs.erase(toNodeInArcIt);
        }

        if (removed)
        {
          auto elementInserted = fixedArcs.insert(arcIndex);
          arc.heuristicFlow = 0.0;
          arc.decompositionFlow = 0.0;
        }
      }
    }
  }

  DBG(std::cout << "number of fixed arcs: " << fixedArcs.size() << std::endl;
  std::cout << "percent fixed arcs: " << getPercentFixedArcs() << std::endl;)
  //std::cout << "shortest path down = " << shortestPathDown[terminalNodeIndex] << std::endl;
  return getPercentFixedArcs();
}

int VRPTWDecisionDiagram::selectArcWithLargestFlowFromNode(int nodeIndex)
{
  int arcIndexLargestFlow = -1;
  double largestFlow = 0;
  for (int outArcIndex : nodes[nodeIndex].outArcs)
  {
    double flow = arcs[outArcIndex].heuristicFlow;
    if (flow > 0.00001)
    {
      if (flow > largestFlow)
      {
        largestFlow = flow;
        arcIndexLargestFlow = outArcIndex;
      }
    }
  }

  return arcIndexLargestFlow;
}
 
bool VRPTWDecisionDiagram::addLocationToRoute(int location, std::vector<int>& route)
{
  int bestInsertionIndex = -1;
  double bestDistance = INF;
  for (int insertionIndex=1; insertionIndex<route.size(); ++insertionIndex)
  {
    std::vector<int> newRoute;
    for (int routeIndex=0; routeIndex<route.size(); ++routeIndex)
    {
      if (insertionIndex == routeIndex)
      {
        newRoute.push_back(location);
        if (!vrptw.reliances.empty())
        {
          for (int relianceLoc : vrptw.reliances[location])
          {
            newRoute.push_back(relianceLoc);
          }
        }
      }
      newRoute.push_back(route[routeIndex]);
    }
    bool isFeasible = isRouteFeasible(newRoute);
    if (isFeasible)
    {
      double distance = vrptw.evaluateRouteDistance(newRoute);

      if (distance < bestDistance)
      {
        bestDistance = distance;
        bestInsertionIndex = insertionIndex;
      }
    }
  }

  if (bestInsertionIndex == -1)
  {
    return false;
  }
  else
  {
    route.insert(route.begin() + bestInsertionIndex, location);
    int index = 1;
    if (!vrptw.reliances.empty())
    {
      for (int relianceLoc : vrptw.reliances[location])
      {
        route.insert(route.begin() + bestInsertionIndex + index, relianceLoc); 
        ++index;
      }
    }
    return true;
  }
}

// greedily take paths with flow
// choose highest flow value
// follow path down until cannot go further
// after this, try adding elements to current sequences in loop
// can try adding element anywhere in each sequence, and use the best one
bool VRPTWDecisionDiagram::primalHeuristicGreedy(std::vector<std::vector<int>>& routesByLocationPrimalHeuristic)
{
  bool rootFlowExists = true;
  std::set<int> locationsCovered;
  locationsCovered.insert(0);
  while (rootFlowExists && (locationsCovered.size() < vrptw.numLocations))
  {
    int currentNodeIndex = rootNodeIndex;
    bool continueRoute = true;
    std::vector<int> route;
    route.push_back(0);
    while (continueRoute)
    {
      int arcIndex = selectArcWithLargestFlowFromNode(currentNodeIndex);

      if (arcIndex == -1)
      {
        continueRoute = false;
        if (currentNodeIndex == rootNodeIndex)
        {
          rootFlowExists = false;
        }
        break;
      }
      else
      {
        int nextArcLocation = arcs[arcIndex].location;
        int nextNodeIndex = arcs[arcIndex].toNodeIndex;
        if (locationsCovered.find(nextArcLocation) == locationsCovered.end())
        {
          route.push_back(nextArcLocation);
          locationsCovered.insert(nextArcLocation);
        }

        if (nextNodeIndex == terminalNodeIndex)
        {
          continueRoute = false;
        }
        else
        {
          currentNodeIndex = nextNodeIndex;
        }
 
        arcs[arcIndex].heuristicFlow = 0;
      }
    }

    // add route if at least one location was added
    if (route.size() > 1)
    {
      if (route.back() != 0)
      {
        route.push_back(0);
      }
      routesByLocationPrimalHeuristic.push_back(route);
    }
  }
 
  // when done with flow decomposition, greedily add other locations
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    if (locationsCovered.find(location) == locationsCovered.end())
    {
      for (int routeIndex=0; routeIndex<routesByLocationPrimalHeuristic.size(); ++routeIndex)
      {
        std::vector<int> route = routesByLocationPrimalHeuristic[routeIndex];
        bool addedLocation = addLocationToRoute(location, route);
        if (addedLocation)
        {
          locationsCovered.insert(location);
          routesByLocationPrimalHeuristic[routeIndex] = route;
          break;
        }
      }
    }
  }

  // if not all locations covered, try to make new sequences to cover them all
  while (true)
  {
    bool addedNewLocation = false;

    std::vector<int> newRoute;
    newRoute.push_back(0);
    newRoute.push_back(0);
    for (int location=1; location<vrptw.numLocations; ++location)
    {
      if (locationsCovered.find(location) == locationsCovered.end())
      {
        bool addedLocation = addLocationToRoute(location, newRoute);
        if (addedLocation)
        {
          locationsCovered.insert(location);
          addedNewLocation = true;
        }
      }
    }

    if (newRoute.size() > 2)
    {
      routesByLocationPrimalHeuristic.push_back(newRoute);
    }

    if (!addedNewLocation)
    {
      break;
    }
  }

  // if all locations covered, return solution and true
  if (locationsCovered.size() == vrptw.numLocations)
  {
    for (auto primalRoute : routesByLocationPrimalHeuristic)
    {
      std::cout << "primalRoute: ";
      for (int loc : primalRoute)
      {
        std::cout << loc << " ";
      }
      std::cout << std::endl;
    }
    return true;
  }
  else
  {
    std::cout << "only " << locationsCovered.size() << " locations used" << std::endl;
    routesByLocationPrimalHeuristic.clear();
    return false;
  }
};

void VRPTWDecisionDiagram::createTruncatedRoute(const std::vector<int>& route, std::vector<int>& truncatedRoute)
{
  VRPTWNodeState currState = nodes[rootNodeIndex].state;
  truncatedRoute.push_back(0);
  for (int index=0; index<(int)route.size(); ++index)
  {
    int nextLocation = route[index];
    if (nextLocation == 0)
    {
      continue;
    }

    VRPTWNodeState newState(currState);
    bool newNodeFeasible = generateNewStateExact(newState, nextLocation);
    if (!newNodeFeasible)
    {
      truncatedRoute.push_back(0);
      return;
    }
    else
    {
      truncatedRoute.push_back(nextLocation);
    }

    currState = newState;
  }
 
  truncatedRoute.push_back(0);
}

bool VRPTWDecisionDiagram::prefixIntraRouteSwaps(std::vector<int>& route, double& routeCost)
{
  std::cout << "step 1" << std::endl;
  for (int removalIndex=1; removalIndex<(int)route.size(); ++removalIndex)
  {
    int removedElement = route[removalIndex];
    for (int insertionIndex=1; insertionIndex<route.size()-1; ++insertionIndex)
    {
      std::vector<int> newRoute = route;
      newRoute.erase(newRoute.begin() + removalIndex);
      std::cout << "erase" << std::endl;
      newRoute.insert(newRoute.begin() + insertionIndex, removedElement);
      std::cout << "insert" << std::endl;
      if (isRouteFeasible(newRoute))
      {
        std::cout << "feasible route" << std::endl;
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        std::cout << "cost calculated" << std::endl;
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "prefix improved with intra-route reinsertion" << std::endl;
          return true;
        }
      }
    }
  }

  std::cout << "step 2" << std::endl;
  for (int removalIndex=1; removalIndex<(int)route.size()-1; ++removalIndex)
  {
    int removedElement1 = route[removalIndex];
    int removedElement2 = route[removalIndex+1];
    for (int insertionIndex=1; insertionIndex<route.size()-1; ++insertionIndex)
    {
      std::vector<int> newRoute = route;
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.insert(newRoute.begin() + insertionIndex, removedElement1);
      newRoute.insert(newRoute.begin() + insertionIndex+1, removedElement2);
      if (isRouteFeasible(newRoute))
      {
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "prefix improved with intra-route Or-opt2" << std::endl;
          return true;
        }
      }
    }
  }

  std::cout << "step 3" << std::endl;
  for (int removalIndex=1; removalIndex<static_cast<int>(route.size())-2; ++removalIndex)
  {
    int removedElement1 = route[removalIndex];
    int removedElement2 = route[removalIndex+1];
    int removedElement3 = route[removalIndex+2];
    for (int insertionIndex=1; insertionIndex<route.size()-2; ++insertionIndex)
    {
      std::vector<int> newRoute = route;
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.insert(newRoute.begin() + insertionIndex, removedElement1);
      newRoute.insert(newRoute.begin() + insertionIndex+1, removedElement2);
      newRoute.insert(newRoute.begin() + insertionIndex+2, removedElement3);
      if (isRouteFeasible(newRoute))
      {
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "prefix improved with intra-route Or-opt3" << std::endl;
          return true;
        }
      }
    }
  }

  std::cout << "step 4" << std::endl;
  for (int index1=1; index1<(int)route.size(); ++index1)
  {
    int element1 = route[index1];
    for (int index2=1; index2<(int)route.size(); ++index2)
    {
      int element2 = route[index2];

      std::vector<int> newRoute = route;
      newRoute[index2] = element1;
      newRoute[index1] = element2;
      if (isRouteFeasible(newRoute))
      {
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "prefix improved with intra-route exchange" << std::endl;
          return true;
        }
      }
    }
  }

  std::cout << "step 5" << std::endl;
  for (int index1=1; index1<(int)route.size()-2; ++index1)
  {
    for (int index2=index1+1; index2<(int)route.size(); ++index2)
    {
      std::vector<int> newRoute;
      for (int i=0; i<=index1; ++i)
      {
        newRoute.push_back(route[i]);
      }
      for (int i=index2; i>index1; --i)
      {
        newRoute.push_back(route[i]);
      }
      for (int i=index2+1; i<(int)route.size(); ++i)
      {
        newRoute.push_back(route[i]);
      }

      if (isRouteFeasible(newRoute))
      {
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "prefix improved with intra-route 2-opt" << std::endl;
          return true;
        }
      }
    }
  }

  return false;
}

bool VRPTWDecisionDiagram::intraRouteSwaps(std::vector<int>& route, double& routeCost)
{
  for (int removalIndex=1; removalIndex<route.size()-1; ++removalIndex)
  {
    int removedElement = route[removalIndex];
    for (int insertionIndex=1; insertionIndex<route.size()-2; ++insertionIndex)
    {
      std::vector<int> newRoute = route;
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.insert(newRoute.begin() + insertionIndex, removedElement);
      if (isRouteFeasible(newRoute))
      {
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "improved with intra-route reinsertion" << std::endl;
          return true;
        }
      }
    }
  }

  for (int removalIndex=1; removalIndex<route.size()-2; ++removalIndex)
  {
    int removedElement1 = route[removalIndex];
    int removedElement2 = route[removalIndex+1];
    for (int insertionIndex=1; insertionIndex<route.size()-2; ++insertionIndex)
    {
      std::vector<int> newRoute = route;
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.insert(newRoute.begin() + insertionIndex, removedElement1);
      newRoute.insert(newRoute.begin() + insertionIndex+1, removedElement2);
      if (isRouteFeasible(newRoute))
      {
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "improved with intra-route Or-opt2" << std::endl;
          return true;
        }
      }
    }
  }

  for (int removalIndex=1; removalIndex<static_cast<int>(route.size())-3; ++removalIndex)
  {
    int removedElement1 = route[removalIndex];
    int removedElement2 = route[removalIndex+1];
    int removedElement3 = route[removalIndex+2];
    for (int insertionIndex=1; insertionIndex<route.size()-3; ++insertionIndex)
    {
      std::vector<int> newRoute = route;
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.erase(newRoute.begin() + removalIndex);
      newRoute.insert(newRoute.begin() + insertionIndex, removedElement1);
      newRoute.insert(newRoute.begin() + insertionIndex+1, removedElement2);
      newRoute.insert(newRoute.begin() + insertionIndex+2, removedElement3);
      if (isRouteFeasible(newRoute))
      {
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "improved with intra-route Or-opt3" << std::endl;
          return true;
        }
      }
    }
  }

  for (int index1=1; index1<route.size()-1; ++index1)
  {
    int element1 = route[index1];
    for (int index2=1; index2<route.size()-1; ++index2)
    {
      int element2 = route[index2];

      std::vector<int> newRoute = route;
      newRoute[index2] = element1;
      newRoute[index1] = element2;
      if (isRouteFeasible(newRoute))
      {
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "improved with intra-route exchange" << std::endl;
          return true;
        }
      }
    }
  }

  for (int index1=1; index1<static_cast<int>(route.size())-3; ++index1)
  {
    for (int index2=index1+1; index2<route.size()-1; ++index2)
    {
      std::vector<int> newRoute;
      for (int i=0; i<=index1; ++i)
      {
        newRoute.push_back(route[i]);
      }
      for (int i=index2; i>index1; --i)
      {
        newRoute.push_back(route[i]);
      }
      for (int i=index2+1; i<route.size(); ++i)
      {
        newRoute.push_back(route[i]);
      }

      if (isRouteFeasible(newRoute))
      {
        double newRouteCost = vrptw.evaluateRouteDistance(newRoute);
        if (newRouteCost < routeCost)
        {
          routeCost = newRouteCost;
          route = newRoute;
          //std::cout << "improved with intra-route 2-opt" << std::endl;
          return true;
        }
      }
    }
  }
 
  return false;
}

bool sort_by_second2(const std::pair<std::pair<int,int>,double>& a, const std::pair<std::pair<int,int>,double>& b)
{
  return (a.second < b.second);
}

double VRPTWDecisionDiagram::repairSolution(const std::vector<std::vector<int>>& feasibleSolution, const std::set<int>& destroyedElements, const Dual& dual, std::vector<std::vector<int>>& newBestRoutes, int timeoutSeconds)
{
  double feasibleSolutionValue = vrptw.evaluateSolutionCost(feasibleSolution);
  int feasibleSolutionNumVehicles = static_cast<int>(feasibleSolution.size());
  feasibleSolutionValue -= vrptw.fixedRouteCost * feasibleSolutionNumVehicles;
  double routeCostUpperBound = (1.2 * feasibleSolutionValue) / (feasibleSolutionNumVehicles-1);
  int numSkipped = 0;

  // 1. Build a partial route one-at-a-time, and consider adding destroyed elements
  nodes.reserve(std::pow(vrptw.numLocations,2));
  arcs.reserve(std::pow(vrptw.numLocations,2)/100);

  std::set<int> initialDeque = {};
  VRPTWNodeState rootNodeState(1,0,0,0,true,initialDeque);
  addNode(rootNodeState);
  rootNodeIndex = 0;

  int terminalNodeLoc = 0;
  if (vrptw.circuitOrPath == CircuitOrPath::PATH)
  {
    terminalNodeLoc = vrptw.numLocations-1;
  }
  double terminalEndTime = 0;
  if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
  {
    terminalEndTime = vrptw.endTimes[0]*(vrptw.timeStateMultiplier);
  }
  VRPTWNodeState terminalNodeState(vrptw.numLocations+2,vrptw.capacity+1,terminalEndTime+1,terminalNodeLoc,false,initialDeque);
  addNode(terminalNodeState);
  terminalNodeIndex = 1;
  int maxPriorityValue = 1;
  if (vrptw.counterType == VRPTWCounterType::USE_COUNTER)
  {
    maxPriorityValue = vrptw.numLocations + 2;
  }
  else if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
  {
    maxPriorityValue = terminalEndTime + 2;
  }
  else
  {
    maxPriorityValue = vrptw.capacity + 2;
  }

  // Keep the initial solution to warm start the MIP
  std::set<int> initialPrimalArcIndices;
  std::set<std::pair<int,int>> generatedArcs;

  // Loop over routes
  for (int routeIndex=0; routeIndex<feasibleSolution.size(); ++routeIndex)
  {
    // Create destroyed route
    std::vector<int> destroyedRoute;
    for (int loc : feasibleSolution[routeIndex])
    {
      if (destroyedElements.find(loc) == destroyedElements.end())
      {
        destroyedRoute.push_back(loc);
      }
    }

    // Get costs
    double initialRouteCost = vrptw.evaluateRouteDistance(feasibleSolution[routeIndex]) - vrptw.fixedRouteCost;
    double destroyedRouteCost= vrptw.evaluateRouteDistance(destroyedRoute) - vrptw.fixedRouteCost;

    // Start creating DP
    std::map<int,std::vector<int>> priorityQueue;
    std::map<int,int> nodeRouteIndex;
    std::map<int,int> nodeNumNewPickups;

    // Start with root node
    priorityQueue.insert(std::make_pair(0,std::vector<int>(1,0)));
    nodeRouteIndex[rootNodeIndex] = 0;
    nodeNumNewPickups[rootNodeIndex] = 0;
 
    // Create the initial full route
    VRPTWNodeState currState = nodes[rootNodeIndex].state;
    int currNodeIndex = rootNodeIndex;
    for (int location : feasibleSolution[routeIndex])
    {
      if (location == 0)
      {
        continue;
      }

      // Create next state/node
      VRPTWNodeState newState(currState);
      bool newNodeFeasible = generateNewStateExact(newState, location);
      if (newNodeFeasible)
      {
        int numNodes = nodes.size();
        int nextNodeIndex = addNode(newState);
        int newNumNodes = nodes.size();

        // Create next arc
        int arcIndex = addArc(currNodeIndex, nextNodeIndex);
        auto nodePair = std::make_pair(currNodeIndex, nextNodeIndex);
        generatedArcs.insert(nodePair);
        initialPrimalArcIndices.insert(arcIndex);

        // Update current values
        currNodeIndex = nextNodeIndex;
        currState = newState;
      }
    }

    // Terminal arc
    int fullRouteArcIndex = addArc(currNodeIndex, terminalNodeIndex);
    auto nodePair1 = std::make_pair(currNodeIndex, terminalNodeIndex);
    generatedArcs.insert(nodePair1);
    initialPrimalArcIndices.insert(fullRouteArcIndex);

    // Allow totally destroyed routes to have removals
    if (destroyedRoute.size() == 2)
    {
      std::map<int,std::vector<int>> priorityQueue1;
      std::map<int,int> nodeRouteIndex1;

      // Start with root node
      priorityQueue1.insert(std::make_pair(0,std::vector<int>(1,0)));
      nodeRouteIndex1[rootNodeIndex] = 0;

      int priorityIteratorValue = -1;
      while (priorityIteratorValue <= maxPriorityValue)
      {
        ++priorityIteratorValue;
        if (priorityQueue1.find(priorityIteratorValue) == priorityQueue1.end())
        {
          continue;
        }

        std::vector<int> priorityItems = priorityQueue1.find(priorityIteratorValue)->second;
        for (int nodeIndex : priorityItems)
        {
          const VRPTWNodeState stateToCheck = nodes[nodeIndex].state;

          // Transition with destroyed elements
          auto route = feasibleSolution[routeIndex];
          int currIndex = nodeRouteIndex1[nodeIndex];
          for (int index=currIndex+1; index<route.size(); ++index)
          {
            int location = route[index];
            if (location == 0)
            {
              continue;
            }

            VRPTWNodeState newState(stateToCheck);
            bool newNodeFeasible = generateNewStateExact(newState, location);
            if (newNodeFeasible)
            {
              int numNodes = nodes.size();
              int nextNodeIndex = addNode(newState);
              int newNumNodes = nodes.size();
              int priorityQueueKeyValue = -1;
              if (vrptw.counterType == VRPTWCounterType::USE_COUNTER)
              {
                priorityQueueKeyValue = newState.counter;
              }
              else if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
              {
                priorityQueueKeyValue = newState.timeWithMultiplier;
              }
              else
              {
                priorityQueueKeyValue = newState.load;
              }

              if (priorityQueue1.find(priorityQueueKeyValue) != priorityQueue1.end())
              {
                priorityQueue1[priorityQueueKeyValue].push_back(nextNodeIndex);
              }
              else
              {
                const auto newElement = std::make_pair(priorityQueueKeyValue, std::vector<int>(1,nextNodeIndex));
                priorityQueue1.insert(newElement);
              }

              nodeRouteIndex1[nextNodeIndex] = index;

              // Create arc from state to state
              std::pair<int,int> nodePair = std::make_pair(nodeIndex, nextNodeIndex);
              if (generatedArcs.find(nodePair) == generatedArcs.end())
              {
                generatedArcs.insert(nodePair);
                addArc(nodePair.first, nodePair.second);
              }

              // Create arcs to terminal node
              if (newState.counter % 2 == 1)
              {
                auto nodePair = std::make_pair(nextNodeIndex, terminalNodeIndex);
                if (generatedArcs.find(nodePair) == generatedArcs.end())
                {
                  generatedArcs.insert(nodePair);
                  addArc(nextNodeIndex, terminalNodeIndex);
                }
              }
            }
          }
        }
      }
      continue;
    }
 
    // Create the initial destroyed route
    currState = nodes[rootNodeIndex].state;
    currNodeIndex = rootNodeIndex;
    int currCost = 0;
    for (int location : destroyedRoute)
    {
      if (location == 0)
      {
        continue;
      }

      // Create next state/node
      VRPTWNodeState newState(currState);
      bool newNodeFeasible = generateNewStateExact(newState, location);
      if (newNodeFeasible)
      {
        currCost += vrptw.distances[currState.lastVisited][location];
        int numNodes = nodes.size();
        int nextNodeIndex = addNode(newState);
        int newNumNodes = nodes.size();
        nodeRouteIndex[nextNodeIndex] = nodeRouteIndex[currNodeIndex] + 1;
        nodeNumNewPickups[nextNodeIndex] = 0;
        int priorityQueueKeyValue = currCost;
    
        if (priorityQueue.find(priorityQueueKeyValue) != priorityQueue.end())
        {
          priorityQueue[priorityQueueKeyValue].push_back(nextNodeIndex);
        }
        else
        {
          const auto newElement = std::make_pair(priorityQueueKeyValue, std::vector<int>(1,nextNodeIndex));
          priorityQueue.insert(newElement);
        }

        // Create next arc
        auto nodePair = std::make_pair(currNodeIndex, nextNodeIndex);
        if (generatedArcs.find(nodePair) == generatedArcs.end())
        {
          generatedArcs.insert(nodePair);
          addArc(currNodeIndex, nextNodeIndex);
        }

        // Update current values
        currNodeIndex = nextNodeIndex;
        currState = newState;
      }
    }

    // Terminal arc
    auto nodePair2 = std::make_pair(currNodeIndex, terminalNodeIndex);
    if (generatedArcs.find(nodePair2) == generatedArcs.end())
    {
      generatedArcs.insert(nodePair2);
      addArc(nodePair2.first, nodePair2.second);
    }

    // Limit time and number of nodes
    auto startTime = std::chrono::high_resolution_clock::now();

    // Create restricted dp with possible insertions of destroyed elements in the route
    int priorityIteratorValue = -1;
    while (priorityIteratorValue <= maxPriorityValue)
    {
      auto currTime = std::chrono::high_resolution_clock::now();
      auto numSeconds = std::chrono::duration_cast<std::chrono::seconds>(currTime - startTime).count();
      if (numSeconds > timeoutSeconds)
      {
        break;
      }

      ++priorityIteratorValue;
      if (priorityQueue.find(priorityIteratorValue) == priorityQueue.end())
      {
        continue;
      }

      std::vector<int> priorityItems = priorityQueue.find(priorityIteratorValue)->second;
      for (int nodeIndex : priorityItems)
      {
        int index = nodeRouteIndex[nodeIndex];
        int numPickups = nodeNumNewPickups[nodeIndex];
        const VRPTWNodeState stateToCheck = nodes[nodeIndex].state;
 
        // Transition with destroyed elements
        for (int location : destroyedElements)
        {
          VRPTWNodeState newState(stateToCheck);
          bool newNodeFeasible = generateNewStateExact(newState, location);
          if (newNodeFeasible)
          {
            // Need to check that the rest of the route works
            int delivery = -1;
            int pickup = -1;
            if (!vrptw.reliances.empty() && !vrptw.reliances[location].empty())
            {
              pickup = location;
              delivery = *(vrptw.reliances[location].begin());
            }
            else
            {
              pickup = location;
              if (vrptw.problemType == ProblemType::PDP)
              {
                continue;
              }
            }

            // For pickups, check all possible delivery indices, add the new full route if possible
            int newStateCost = priorityIteratorValue + vrptw.distances[newState.lastVisited][location];
            for (int insertDeliveryIndex=index+1; insertDeliveryIndex<destroyedRoute.size(); ++insertDeliveryIndex)
            {
              if ((vrptw.problemType != ProblemType::PDP) && (insertDeliveryIndex > index+1))
              {
                continue;
              }

              std::vector<VRPTWNodeState> newStates = {newState};
              std::vector<int> newStateRouteIndices = {index};
              std::vector<int> newStateCosts = {newStateCost};
              VRPTWNodeState suffixState(newState);
              bool isSuffixFeasible = true;
              for (int suffixIndex=index+1; suffixIndex<destroyedRoute.size(); ++suffixIndex)
              {
                // Insert delivery index
                int suffixLocation = destroyedRoute[suffixIndex];
                if (vrptw.problemType == ProblemType::PDP)
                {
                  if (insertDeliveryIndex == suffixIndex)
                  {
                    bool newNodeFeasible1 = generateNewStateExact(suffixState, delivery);
                    if (newNodeFeasible1)
                    {
                      newStates.push_back(suffixState);
                      newStateRouteIndices.push_back(suffixIndex-1);
                      newStateCost += vrptw.distances[suffixLocation][delivery];
                      newStateCosts.push_back(newStateCost);
                    }
                    else
                    {
                      isSuffixFeasible = false;
                      break;
                    }
                  }
                }

                // Insert suffix location
                if (suffixLocation != 0)
                {
                  bool newNodeFeasible2 = generateNewStateExact(suffixState, suffixLocation);
                  if (newNodeFeasible2)
                  {
                    newStates.push_back(suffixState);
                    newStateRouteIndices.push_back(suffixIndex);
                    newStateCost += vrptw.distances[suffixLocation][delivery];
                    newStateCosts.push_back(newStateCost);
                  }
                  else
                  {
                    isSuffixFeasible = false;
                    break;
                  }
                }
              }

              // For feasible suffix, create nodes and arcs
              if (!isSuffixFeasible)
              {
                break;
              }
              else
              {
                // Limit with cost analysis during route creation (check UB vs. # routes)
                /*
                if (newStateCost > routeCostUpperBound)
                {
                  ++numSkipped;
                  continue;
                }
                */

                // Limit options with randomness
                // Use idea from Pisiginer, and check for how much the insertion changes obj?
                // Granular neighborhood idea?
                /*
                int limitingRandomInteger = std::rand();
                int limitingRandomIntegerModTen = limitingRandomInteger % 10;
                if (limitingRandomIntegerModTen <= 0)
                {
                  continue;
                }
                */

                int previousNodeIndex = -1;
                for (int newStateIndex=0; newStateIndex<newStates.size(); ++newStateIndex)
                {
                  auto state = newStates[newStateIndex];
                  int numNodes = nodes.size();
                  int newNodeIndex = addNode(state);
                  int newNumNodes = nodes.size();
                  int stateRouteIndex = newStateRouteIndices[newStateIndex];
                  nodeRouteIndex[newNodeIndex] = stateRouteIndex;
                  nodeNumNewPickups[newNodeIndex] = numPickups + 1;
                  int priorityQueueKeyValue = newStateCosts[newStateIndex];

                  if (priorityQueue.find(priorityQueueKeyValue) != priorityQueue.end())
                  {
                    priorityQueue[priorityQueueKeyValue].push_back(newNodeIndex);
                  }
                  else
                  {
                    const auto newElement = std::make_pair(priorityQueueKeyValue, std::vector<int>(1,newNodeIndex));
                    priorityQueue.insert(newElement);
                  }

                  // Create arc from state to state
                  std::pair<int,int> nodePair;
                  if (previousNodeIndex > -1)
                  {
                    nodePair = std::make_pair(previousNodeIndex, newNodeIndex);
                  }
                  else
                  {
                    nodePair = std::make_pair(nodeIndex, newNodeIndex);
                  }
                  if (generatedArcs.find(nodePair) == generatedArcs.end())
                  {
                    generatedArcs.insert(nodePair);
                    addArc(nodePair.first, nodePair.second);
                  }

                  // Create arcs to terminal node
                  if (stateRouteIndex == static_cast<int>(destroyedRoute.size()-2))
                  {
                    auto nodePair = std::make_pair(newNodeIndex, terminalNodeIndex);
                    if (generatedArcs.find(nodePair) == generatedArcs.end())
                    {
                      generatedArcs.insert(nodePair);
                      addArc(newNodeIndex, terminalNodeIndex);
                    }
                  }
 
                  previousNodeIndex = newNodeIndex;
                }
              }
            }
          }
        }
      }
    }
  }

/*
  // Add some short routes with only destroyed elements
  std::map<std::pair<int,std::set<int>>, int> visitedSetBestCosts;
  std::map<std::pair<int,std::set<int>>, int> visitedSetBestState;
  std::map<int, int> nodeIndexCosts;
  std::set<int> dominatedNodeIndices;

  std::map<int,std::vector<int>> priorityQueue2;
  priorityQueue2.insert(std::make_pair(0,std::vector<int>(1,0)));
  nodeIndexCosts[rootNodeIndex] = 0;
  int priorityIteratorValue = -1;
  auto startTime = std::chrono::high_resolution_clock::now();
  while (priorityIteratorValue <= maxPriorityValue)
  {
    auto currTime = std::chrono::high_resolution_clock::now();
    auto numSeconds = std::chrono::duration_cast<std::chrono::seconds>(currTime - startTime).count();
    if (numSeconds > timeoutSeconds)
    {
      break;
    }

    ++priorityIteratorValue;
    if (priorityQueue2.find(priorityIteratorValue) == priorityQueue2.end())
    {
      continue;
    }

    std::vector<int> priorityItems = priorityQueue2.find(priorityIteratorValue)->second;
    for (int nodeIndex : priorityItems)
    {
      if (dominatedNodeIndices.find(nodeIndex) != dominatedNodeIndices.end())
      {
        continue;
      }
      const VRPTWNodeState stateToCheck = nodes[nodeIndex].state;
      int stateToCheckCost = nodeIndexCosts[nodeIndex];

      // Transition with destroyed elements
      for (int location : destroyedElements)
      {
        VRPTWNodeState newState(stateToCheck);
        bool newNodeFeasible = generateNewStateExact(newState, location);
        if (newNodeFeasible)
        {
          int numNodes = nodes.size();
          int nextNodeIndex = addNode(newState);
          int newNumNodes = nodes.size();

          // Check cost and dominance
          int newCost = stateToCheckCost + vrptw.distances[stateToCheck.lastVisited][location];
          nodeIndexCosts[nextNodeIndex] = newCost;
          auto visitedLastVisited = std::make_pair(newState.lastVisited, newState.visited);
          auto visitedSetBestCostIter = visitedSetBestCosts.find(visitedLastVisited);
          if (visitedSetBestCostIter != visitedSetBestCosts.end())
          {
            int currBest = visitedSetBestCostIter->second;
            if (newCost < currBest)
            {
              visitedSetBestCosts[visitedLastVisited] = newCost;
              dominatedNodeIndices.insert(visitedSetBestState[visitedLastVisited]);
              visitedSetBestState[visitedLastVisited] = nextNodeIndex;
            }
            else
            {
              continue;
            }
          }

          int priorityQueueKeyValue = -1;
          if (vrptw.counterType == VRPTWCounterType::USE_COUNTER)
          {
            priorityQueueKeyValue = newState.counter;
          }
          else if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
          {
            priorityQueueKeyValue = newState.timeWithMultiplier;
          }
          else
          {
            priorityQueueKeyValue = newState.load;
          }

          if (priorityQueue2.find(priorityQueueKeyValue) != priorityQueue2.end())
          {
            priorityQueue2[priorityQueueKeyValue].push_back(nextNodeIndex);
          }
          else
          {
            const auto newElement = std::make_pair(priorityQueueKeyValue, std::vector<int>(1,nextNodeIndex));
            priorityQueue2.insert(newElement);
          }

          // Create arc from state to state
          std::pair<int,int> nodePair = std::make_pair(nodeIndex, nextNodeIndex);
          if (generatedArcs.find(nodePair) == generatedArcs.end())
          {
            generatedArcs.insert(nodePair);
            addArc(nodePair.first, nodePair.second);
          }

          // Create arcs to terminal node
          if (newState.counter % 2 == 1)
          {
            auto nodePair = std::make_pair(nextNodeIndex, terminalNodeIndex);
            if (generatedArcs.find(nodePair) == generatedArcs.end())
            {
              generatedArcs.insert(nodePair);
              addArc(nextNodeIndex, terminalNodeIndex);
            }
          }
        }
      }
    }
  }
*/
  //print();
  //std::cout << "num skipped: " << numSkipped << std::endl;

  // 2. Solve the Arc Flow Formulation as a MIP
  if (vrptw.problemType == ProblemType::PDP)
  {
    setCoeffsAsPreciseDistances();
  }
  else
  {
    setCoeffsAsDistances();
  }
  Dual duals;
  double mipObjective = setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, initialPrimalArcIndices, duals, false, std::max(1,std::min(5,timeoutSeconds)));
  if (mipObjective < 0)
  {
    std::cout << "repair solution issue" << std::endl;
    print();
    return mipObjective;
  }

  // 3. Return the solution
  std::vector<int> infeasibleRoute;
  std::vector<std::vector<int>> decomposedArcs;
  std::vector<double> routeFlows;
  std::vector<std::vector<int>> decomposedRoutes;
  decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, DecompositionReason::DECOMPOSE);

  newBestRoutes.clear();
  for (auto route : decomposedRoutes)
  {
    std::vector<int> newRoute;
    for (int loc : route)
    {
      newRoute.push_back(loc);
    }

    newBestRoutes.push_back(newRoute);
  }

  return mipObjective;
}

bool VRPTWDecisionDiagram::doesRouteExistByArcs(const std::vector<int>& routeArcs, std::vector<int>& routeLocations) const
{
  routeLocations.push_back(0);
  for (int index=0; index<(routeArcs.size()-1); ++index)
  {
    int currNodeIndex = arcs[routeArcs[index]].toNodeIndex;
    int nextNodeIndex = arcs[routeArcs[index+1]].fromNodeIndex;
    if (currNodeIndex != nextNodeIndex)
    {
      return false;
    }
    else
    {
      routeLocations.push_back(arcs[routeArcs[index]].location);
    }
  }

  return true;
};

bool VRPTWDecisionDiagram::doesRouteExistByLocations(const std::vector<int>& routeLocations, std::vector<int>& routeArcs) const
{
  // check regular r-t paths
  bool regularPathExists = true;

  int routeIndex = 0;
  if (routeLocations[0] == 0)
  {
    routeIndex = 1;
  }

  int currNodeIndex = rootNodeIndex;
  while (currNodeIndex != terminalNodeIndex)
  {
    bool nextLocationPossible = false;
    for (int arcIndex : nodes[currNodeIndex].outArcs)
    {
      if (arcs[arcIndex].location == routeLocations[routeIndex])
      {
        routeArcs.push_back(arcIndex);
        currNodeIndex = arcs[arcIndex].toNodeIndex;
        nextLocationPossible = true;
        break;
      }
    }

    if (!nextLocationPossible)
    {
      regularPathExists = false;
      break;
    }

    routeIndex = routeIndex + 1;
  }

  // check all special 0-starting SRC paths
  bool specialSRCPathExists = false;
  if (!regularPathExists)
  {
    for (int arcIndex : nodes[rootNodeIndex].outArcs)
    {
      routeArcs.clear();
      routeArcs.push_back(arcIndex);
      if (arcs[arcIndex].location == 0)
      {
        // check special SRC arc
        int routeIndex = 0;
        if (routeLocations[0] == 0)
        {
          routeIndex = 1;
        }

        int currNodeIndex = arcs[arcIndex].toNodeIndex;
        while (currNodeIndex != terminalNodeIndex)
        {
          bool nextLocationPossible = false;
          for (int currNodeArcIndex: nodes[currNodeIndex].outArcs)
          {
            if (arcs[currNodeArcIndex].location == routeLocations[routeIndex])
            {
              routeArcs.push_back(currNodeArcIndex);
              currNodeIndex = arcs[currNodeArcIndex].toNodeIndex;
              nextLocationPossible = true;
              break;
            }
          }

          if (!nextLocationPossible)
          {
            break;
          }

          routeIndex = routeIndex + 1;
          if (routeIndex == routeLocations.size())
          {
            specialSRCPathExists = true;
            break;
          }
        }
      }

      if (specialSRCPathExists)
      {
        break;
      }
    }
  }

  if (!regularPathExists && !specialSRCPathExists)
  {
    routeArcs.clear();
    return false;
  }

  return true;
};

void VRPTWDecisionDiagram::decomposeRoutes(std::vector<int>& routeArcs, std::vector<double>& flows, std::vector<std::vector<int>>& routeDecomposition, std::vector<std::vector<int>>& decomposedArcs, DecompositionReason decompositionReason)
{
  // store decomposition flows to restore after
  std::vector<double> decompositionFlows;
  if (decompositionReason == DecompositionReason::DECOMPOSE)
  {
    for (int index=0; index<arcs.size(); ++index)
    {
      decompositionFlows.push_back(arcs[index].decompositionFlow);
    }
  }

  std::vector<int> route;
  std::vector<int> routeByLoc;
  std::deque<int> lastLocationsVisited;
  std::deque<int> lastArcsVisited;
  std::deque<int> lastNodesVisited;
  bool continueDecomposing = true;
  while (continueDecomposing)
  {
    route.clear();
    routeByLoc.clear();
    routeByLoc.push_back(0);
    lastLocationsVisited.clear();
    lastArcsVisited.clear();
    lastNodesVisited.clear();

    bool continueRoute = true;
    int currentNodeIndex = rootNodeIndex;
    VRPTWNodeState rootNodeState(nodes[rootNodeIndex].state);
    VRPTWNodeState currState(rootNodeState);
    double routeFlow = INF;
    while (continueRoute)
    {
      bool arcWithFlowFound = false;
      for (int arcIndex : nodes[currentNodeIndex].outArcs)
      {
        if (arcs[arcIndex].decompositionFlow > 0.00001)
        {
          arcWithFlowFound = true;
          currentNodeIndex = arcs[arcIndex].toNodeIndex;
          routeFlow = std::min(routeFlow, arcs[arcIndex].decompositionFlow);
          route.push_back(arcIndex);
          routeByLoc.push_back(arcs[arcIndex].location);

          if (currentNodeIndex == terminalNodeIndex)
          {
            continueRoute = false;
          }
          else
          {
            int currLocation = routeByLoc[routeByLoc.size()-1];
            int lastLocation = routeByLoc[routeByLoc.size()-2];

            VRPTWNodeState nextState(currState);
            bool nextNodeFeasible = generateNewStateExact(nextState, currLocation);
            currState = nextState;
            if (decompositionReason == DecompositionReason::SEPARATE)
            {
              if (!nextNodeFeasible)
              {
                for (int index=0; index<route.size(); ++index)
                {
                  routeArcs.push_back(route[index]);
                }
                for (int routeArcIndex : route)
                {
                  arcs[routeArcIndex].decompositionFlow = arcs[routeArcIndex].decompositionFlow - routeFlow;
                }
                return;
              }
            }
          }

          break;
        }
      }

      if (!arcWithFlowFound && (currentNodeIndex == rootNodeIndex))
      {
        if (decompositionReason == DecompositionReason::DECOMPOSE)
        {
          for (int index=0; index<arcs.size(); ++index)
          {
            arcs[index].decompositionFlow = decompositionFlows[index];
          }
        }
        return;
      }
    }

    flows.push_back(routeFlow);
    routeDecomposition.push_back(routeByLoc);
    decomposedArcs.push_back(route);
    
    for (int arcIndex : route)
    {
      arcs[arcIndex].decompositionFlow = arcs[arcIndex].decompositionFlow - routeFlow;
    }

    /*
    std::cout << "[loc] route: ";
    for (int arcIndex : route)
    {
      std::cout << arcs[arcIndex].location << " ";
    }
    std::cout << std::endl;
 
    std::cout << "[arc index] route: ";
    for (int arcIndex : route)
    {
      std::cout << arcIndex << " ";
      arcs[arcIndex].decompositionFlow = arcs[arcIndex].decompositionFlow - routeFlow;
    }
    std::cout << std::endl;
    */
  }
};

void VRPTWDecisionDiagram::getSolutionArcs(std::set<int>& solutionArcs)
{
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    if (arcs[arcIndex].heuristicFlow > 0.000001)
    {
      solutionArcs.insert(arcIndex);
    }
  }
};

void VRPTWDecisionDiagram::separateRoute(const std::vector<int>& routeArcs)
{
  // clear fixed arcs when dd changes
  for (int arcIndex : fixedArcs)
  {
    const VRPTWArc& arc = arcs[arcIndex];
    nodes[arc.fromNodeIndex].outArcs.push_back(arcIndex);
    nodes[arc.toNodeIndex].inArcs.push_back(arcIndex);
  }
  fixedArcs.clear();

  std::cout << "separating route arcs: ";
  for (int arcIndex : routeArcs)
  {
    std::cout << arcIndex << " ";
  }
  std::cout << std::endl;

  std::cout << "separating route locations: ";
  for (int arcIndex : routeArcs)
  {
    std::cout << arcs[arcIndex].location << " ";
  }
  std::cout << std::endl;

  // index to routeArcs.size()-1 because the last arc will be removed because it causes a conflict
  bool haveMovedFirstArc = false;
  int currNodeIndex = arcs[routeArcs[0]].fromNodeIndex;
  for (int index=0; index<(routeArcs.size()-1); ++index)
  {
    int nextLocation = arcs[routeArcs[index]].location;

    // new node on next layer
    int nextNodeIndex = -1;
    bool newNodeCreated = false;
    int prevNumNodes = nodes.size();
    VRPTWNodeState newState(nodes[currNodeIndex].state);
    bool newNodeFeasible = generateNewStateExact(newState, nextLocation);

    // if feasible, create new state if necessary, etc.
    if (newNodeFeasible)
    {
      nextNodeIndex = addNode(newState);
      if (prevNumNodes < nodes.size())
      {
        newNodeCreated = true;
      }

      // if new node created, copy outgoing arcs
      if (newNodeCreated)
      {
        DBG(std::cout << "new node created: " << nodes.size() - 1 << std::endl;)
        for (int oldArcIndex : nodes[arcs[routeArcs[index]].toNodeIndex].outArcs)
        {
          int oldArcLocation = arcs[oldArcIndex].location;
          int oldArcToNodeIndex = arcs[oldArcIndex].toNodeIndex;
          VRPTWNodeState copyArcState(newState);
          bool copyArcFeasible = generateNewStateExact(copyArcState, oldArcLocation);
          if (copyArcFeasible)
          {
            int newArcIndex = addArc(nextNodeIndex, oldArcToNodeIndex);
            DBG(std::cout << "copy arc by adding arc index: " << newArcIndex << std::endl;)
          }
        }
      }
      else
      {
        // if node existed, remove arcs if possible, because we now know it's exact
        for (int existingNodeArcIndex : nodes[nextNodeIndex].outArcs)
        {
          int existingNodeArcLocation = arcs[existingNodeArcIndex].location;
          VRPTWNodeState copyArcState(newState);
          bool existingArcFeasible = generateNewStateExact(copyArcState, existingNodeArcLocation);
          if (!existingArcFeasible)
          {
            //std::cout << "remove arc " << existingNodeArcIndex << std::endl;
            removeArc(existingNodeArcIndex);
          }
        }
      }

      // move arc from curr node to new node, note the location identifies the action
      for (int currNodeArcIndex : nodes[currNodeIndex].outArcs)
      {
        if (arcs[currNodeArcIndex].location == nextLocation)
        {
          if (arcs[currNodeArcIndex].toNodeIndex == nextNodeIndex)
          {
            break;
          }
          else
          {
            DBG(std::cout << "move arc " << currNodeArcIndex << " to node " << nextNodeIndex << std::endl;)
            moveArc(currNodeArcIndex, nextNodeIndex);
          }
        }
      }

      currNodeIndex = nextNodeIndex;
    }
    else
    {
      // if new node not feasible, remove the arc and return
      for (int currNodeArcIndex : nodes[currNodeIndex].outArcs)
      {
        if (arcs[currNodeArcIndex].location == nextLocation)
        {
          //std::cout << "remove arc " << currNodeArcIndex << std::endl;
          removeArc(currNodeArcIndex);
          return;
        }
      }
      return;
    }
  }
 
  DBG(std::cout << "done separating routes" << std::endl;)
};

bool VRPTWDecisionDiagram::isRouteFeasible(const std::vector<int>& route)
{
  VRPTWNodeState currState = nodes[rootNodeIndex].state;
  for (int index=0; index<route.size(); ++index)
  {
    int nextLocation = route[index];
    VRPTWNodeState newState(currState);
    bool newNodeFeasible = generateNewStateExact(newState, nextLocation);
    if (!newNodeFeasible)
    {
      return false;
    }

    currState = newState;
  }

  return true;
};

void VRPTWDecisionDiagram::repairRoute(const std::vector<int>& route, std::vector<int>& feasibleRoute, std::vector<int>& feasibleRouteArcs)
{
  VRPTWNode currNode = nodes[rootNodeIndex];
  VRPTWNodeState currState = nodes[rootNodeIndex].state;
  for (int index=0; index<route.size(); ++index)
  {
    int nextLocation = route[index];

    VRPTWNodeState newState(currState);
    bool newNodeFeasible = generateNewStateExact(newState, nextLocation);
    if (newNodeFeasible)
    {
      feasibleRoute.push_back(nextLocation);
      for (int arcIndex : currNode.outArcs)
      {
        if (arcs[arcIndex].location == nextLocation)
        {
          feasibleRouteArcs.push_back(arcIndex);
          break;
        }
      }

      currState = newState;
    }
  }
};

// note need to be careful with potentials, coeffs, cijPi
// shortest distance -> ( (nodeIndex, arcIn arcIndex)
typedef std::pair<double,std::pair<int,int>> pqKeyType;
void VRPTWDecisionDiagram::runDijkstra(ShortestPathMode mode, std::vector<int>& shortestPathByArc)
{
  DBG(std::cout << "run dijkstra" << std::endl;)
  for (VRPTWNode& node : nodes)
  {
    node.shortestPathDistance = INF;
  }

  // implement as min cost priority queue pq
  std::priority_queue<pqKeyType, std::vector<pqKeyType>, std::greater<pqKeyType>> pq;
  std::set<int> permanent;
  std::vector<int> shortestPathsArcIn;
  shortestPathsArcIn.resize(nodes.size());

  // add source vertex into pq
  pq.push(std::make_pair(0, std::make_pair(rootNodeIndex, rootNodeIndex)));
  nodes[rootNodeIndex].shortestPathDistance = 0;
  nodes[rootNodeIndex].potential = 0;

  // while not all distances set or pq not empty
  while (!pq.empty())
  {
    // extract min distance vertex u from pq
    pqKeyType u = pq.top();
    pq.pop();
    int nodeIndex = u.second.first;
    int fromArcIndex = u.second.second;
    DBG(std::cout << "popped node: " << nodeIndex << " from arc: " << fromArcIndex << std::endl;)
    if (permanent.find(nodeIndex) == permanent.end())
    {
      permanent.insert(nodeIndex);
      shortestPathsArcIn[nodeIndex] = fromArcIndex;
    }
    else
    {
      continue;
    }

    // relax arcs out of popped node
    for (int arcIndex : nodes[nodeIndex].outArcs)
    {
      // if shorter path to v through u, update distance v and add to pq
      int toNodeIndex = arcs[arcIndex].toNodeIndex;
      double possibleDistance = (nodes[nodeIndex].shortestPathDistance + arcs[arcIndex].cijPi);
      if (possibleDistance <= nodes[toNodeIndex].shortestPathDistance)
      {
        nodes[toNodeIndex].shortestPathDistance = possibleDistance;
        if (permanent.find(toNodeIndex) == permanent.end())
        {
          pq.push(std::make_pair(possibleDistance, std::make_pair(toNodeIndex, arcIndex)));
        }
      }
    }
  }

  // create shortest path
  bool endRoute = false;
  int currNodeIndex = terminalNodeIndex;
  if (mode == ShortestPathMode::SHORTEST_PATH)
  {
    while (!endRoute)
    {
      int arcIndex = shortestPathsArcIn[currNodeIndex];
      currNodeIndex = arcs[arcIndex].fromNodeIndex;
      shortestPathByArc.push_back(arcIndex);
      if (currNodeIndex == rootNodeIndex)
      {
        endRoute = true;
      }
    }

    std::reverse(shortestPathByArc.begin(), shortestPathByArc.end());
    DBG(
      std::cout << "shortest path by arc: ";
      for (int a : shortestPathByArc)
      {
        std::cout << a << ", ";
      }
      std::cout << std::endl;
    )
  }
 
  // update potentials and reduced costs
  if (mode == ShortestPathMode::UPDATE_POTENTIALS)
  {
    for (VRPTWNode& node : nodes)
    {
      // pi = pi - d;
      node.potential = node.potential - node.shortestPathDistance;
    }

    for (VRPTWNode& node : nodes)
    {
      for (int arcIndex : node.outArcs)
      {
        // c^pi_ij = c_ij + pi(i) - pi(j)
        const VRPTWNode& toNode = nodes[arcs[arcIndex].toNodeIndex];
        arcs[arcIndex].cijPi = arcs[arcIndex].coeff - node.potential + toNode.potential;
      }
    }
  }
}

int VRPTWDecisionDiagram::reverseArc(int forwardArcIndex)
{
  VRPTWArc forwardArc = arcs[forwardArcIndex];

  // add reverse arc
  int reverseArcIndex = -1;
  auto reverseArcIt = arcReverseArc.find(forwardArcIndex);
  if (reverseArcIt == arcReverseArc.end())
  {
    reverseArcIndex = addReverseArc(forwardArcIndex);
  }
  else
  {
    reverseArcIndex = reverseArcIt->second;

    arcs[reverseArcIndex].coeff = -forwardArc.coeff;
    arcs[reverseArcIndex].cijPi = -forwardArc.cijPi;
    nodes[forwardArc.fromNodeIndex].inArcs.push_back(reverseArcIndex);
    nodes[forwardArc.toNodeIndex].outArcs.push_back(reverseArcIndex);
  }

  // remove forward arc
  auto fromNodeOutArcIt = std::find(nodes[forwardArc.fromNodeIndex].outArcs.begin(), nodes[forwardArc.fromNodeIndex].outArcs.end(), forwardArcIndex);
  nodes[forwardArc.fromNodeIndex].outArcs.erase(fromNodeOutArcIt);
  auto toNodeInArcIt = std::find(nodes[forwardArc.toNodeIndex].inArcs.begin(), nodes[forwardArc.toNodeIndex].inArcs.end(), forwardArcIndex);
  nodes[forwardArc.toNodeIndex].inArcs.erase(toNodeInArcIt);

  return reverseArcIndex;
};

void VRPTWDecisionDiagram::updateResidualGraph(const std::vector<int>& shortestPathByArc)
{
  for (int arcIndex : shortestPathByArc)
  {
    reverseArc(arcIndex);
  }
}

void VRPTWDecisionDiagram::updateResidualGraphWang(const std::vector<int>& shortestPathByArc, std::vector<int>& treeByParentArc)
{
  for (int arcNum=0; arcNum<shortestPathByArc.size(); ++arcNum)
  {
    int arcIndex = shortestPathByArc[arcNum];
    int reverseArcIndex = reverseArc(arcIndex);
    const VRPTWArc& arc = arcs[arcIndex];

    // update tree too
    if (arc.fromNodeIndex != rootNodeIndex)
    {
      treeByParentArc[arc.fromNodeIndex] = reverseArcIndex;
    }
  }
}

double VRPTWDecisionDiagram::createSolutionFromReverseArcsAndReset()
{
  DBG(std::cout << "creating solution" << std::endl;)
  DBG(print();)
  // only go over arcs currently in the graph
  double solutionValue = 0.0;
  for (VRPTWNode node : nodes)
  {
    for (int arcIndex : node.outArcs)
    {
      if (arcs[arcIndex].isReverseArc)
      {
        int solutionArcIndex = reverseArc(arcIndex);
        arcs[solutionArcIndex].decompositionFlow = 1;
        arcs[solutionArcIndex].heuristicFlow = 1;
        solutionValue = solutionValue + arcs[solutionArcIndex].coeff;
      }
    }
  }

  DBG(print();)
  DBG(std::cout << "solution created" << std::endl;)
  return solutionValue;
};

// 1 be careful - paths by arc may contain reverse arcs
// 2 be careful - not all arcs were flipped and/or clipped if terminated early
double VRPTWDecisionDiagram::createSolutionFromReverseArcsAndResetWang(const std::set<int>& clippedArcs, const std::vector<std::vector<int>>& shortestPathsByArc)
{
  std::cout << "creating solution" << std::endl;
  //print();
  // Reverse arcs
  for (const VRPTWNode& node : nodes)
  {
    std::vector<int> outArcs = node.outArcs;
    for (int arcIndex : outArcs)
    {
      if (arcs[arcIndex].isReverseArc)
      {
        reverseArc(arcIndex);
      }
    }
  }

  // Clipped arcs
  for (int arcIndex : clippedArcs)
  {
    VRPTWArc& arc = arcs[arcIndex];
    nodes[arc.fromNodeIndex].outArcs.push_back(arcIndex);
    nodes[arc.toNodeIndex].inArcs.push_back(arcIndex);
  }

  // Create solution
  double solutionValue = 0.0;
  for (auto path : shortestPathsByArc)
  {
    for (int solutionArcIndex : path)
    {
      DBG(std::cout << arcs[solutionArcIndex].fromNodeIndex << "->" << arcs[solutionArcIndex].toNodeIndex << "->";)
      solutionValue = solutionValue + arcs[solutionArcIndex].coeff;
      if (!arcs[solutionArcIndex].isReverseArc)
      {
        arcs[solutionArcIndex].decompositionFlow = 1;
        arcs[solutionArcIndex].heuristicFlow = 1;
      }
      else
      {
        arcs[arcReverseArc[solutionArcIndex]].decompositionFlow = 0;
        arcs[arcReverseArc[solutionArcIndex]].heuristicFlow = 0;
      }
    }
    DBG(std::cout << std::endl;)
  }

  DBG(print();)
  DBG(std::cout << "solution created" << std::endl;)
  return solutionValue;
};

double VRPTWDecisionDiagram::solveMinCostFlowModel(const std::vector<double>& duals, std::vector<std::vector<int>>& shortestPathsByArc, bool& isDualFeasible, double& minReducedCost)
{
  // start potentials and flows at 0
  for (VRPTWNode& node : nodes)
  {
    node.potential = 0;
  }
  for (VRPTWArc& arc : arcs)
  {
    arc.decompositionFlow = 0.0;
    arc.heuristicFlow = 0.0;
  }

  // while negative s-t path exists, send flow and make updates
  // take at most k (number of trucks) negative s-t paths
  int numPaths = 0;
  bool continueFindingPaths = true;
  while (continueFindingPaths)
  {
    // cannot run dijkstras for shortest path with negative lengths
    // distances should be reduced costs c_ij^pi = c_ij - pi_i + pi_j
    if (numPaths == 0)
    {
      setCoeffsAsDistancesMinusLagrangean(duals);
      std::vector<int> empty;
      minReducedCost = computeShortestPathBFS(ShortestPathMode::UPDATE_POTENTIALS, empty);
      DBG(std::cout << "computed shortest path bfs" << std::endl;)
      DBG(print();)

      // shortest path tells us dual feasibility
      if (minReducedCost < 0.00001)
      {
        isDualFeasible = false;
        DBG(std::cout << "dual not feasible: " << shortestPathLength << std::endl;)
      }
      else
      {
        isDualFeasible = true;
      }
    }
    else
    {
      // need to update potentials to all nodes not just find shortest path...
      std::vector<int> empty;
      runDijkstra(ShortestPathMode::UPDATE_POTENTIALS, empty);
    }
    DBG(std::cout << "potentials updated" << std::endl;)
    DBG(print();)

    // calculate shortest path
    std::vector<int> shortestPathByArc;
    runDijkstra(ShortestPathMode::SHORTEST_PATH, shortestPathByArc);

    // convert back to real distance to check if it is negative cost
    double shortestPathDistance = 0;
    for (int arcIndex : shortestPathByArc)
    {
      shortestPathDistance = shortestPathDistance + arcs[arcIndex].coeff;
    }

    // terminate if shortest path distance is positive
    if ((shortestPathDistance >= -0.000000001) && (numPaths != 0))
    {
      break;
    }

    // get residual graph and continue or decide to stop
    numPaths = numPaths + 1;

    DBG(std::cout << "updating residual graph" << std::endl;)
    DBG(print();)
    updateResidualGraph(shortestPathByArc);
    DBG(std::cout << "updated residual graph" << std::endl;)
    DBG(print();)
    shortestPathsByArc.push_back(shortestPathByArc);

    DBG(
    std::cout << "path: ";
    for (int arcIndex : shortestPathByArc)
    {
      std::cout << arcs[arcIndex].fromNodeIndex << ", ";
    }
    std::cout << std::endl;)

    if (nodes[rootNodeIndex].outArcs.empty())
    {
      continueFindingPaths = false;
    }
  }

  DBG(
  for (auto path : shortestPathsByArc)
  {
    for (int solutionArcIndex : path)
    {
      std::cout << arcs[solutionArcIndex].fromNodeIndex << "->" << arcs[solutionArcIndex].toNodeIndex << "->";
    }
    std::cout << std::endl;
  }
  )
  // create solution from reverse arcs and clean up residual graph backward arcs
  DBG(std::cout << "finished running, now create solution" << std::endl;)
  return createSolutionFromReverseArcsAndReset();
};

void VRPTWDecisionDiagram::getTreeByChildArcsFromTreeByParentArcs(const std::vector<int>& treeByParentArcs, std::vector<std::vector<int>>& treeByChildArcs)
{
  treeByChildArcs.clear();
  treeByChildArcs.resize(treeByParentArcs.size());
  for (int nodeIndex=0; nodeIndex<treeByParentArcs.size(); ++nodeIndex)
  {
    if (nodeIndex != rootNodeIndex)
    {
      int arcIndex = treeByParentArcs[nodeIndex];
      const VRPTWArc& arc = arcs[arcIndex];
      treeByChildArcs[arc.fromNodeIndex].push_back(arcIndex);
    }
  }
}
 
void VRPTWDecisionDiagram::identifyNodesForUpdate(int branchNode, const std::vector<std::vector<int>>& treeByChildArcs, std::set<int>& nodesToUpdate)
{
  std::vector<int> newNodesToUpdate;
  newNodesToUpdate.push_back(branchNode);
  for (int newNodeIndex=0; newNodeIndex<newNodesToUpdate.size(); ++newNodeIndex)
  {
    int nodeIndexToUpdate = newNodesToUpdate[newNodeIndex];
    for (int arcIndex : treeByChildArcs[nodeIndexToUpdate])
    {
      const VRPTWArc& arc = arcs[arcIndex];
      newNodesToUpdate.push_back(arc.toNodeIndex);
    }
  }

  for (int nodeIndex : newNodesToUpdate)
  {
    nodesToUpdate.insert(nodeIndex);
  }
}

bool sort_by_second1(const std::pair<int,double>& a, const std::pair<int,double>& b)
{
  return (a.second < b.second);
}

int VRPTWDecisionDiagram::findMultiPathNode(const std::set<int>& nodesToUpdate, int currentNumPaths)
{
  // get all distances
  std::vector<std::pair<int,double>> nodeIndexDistanceToTerminal;
  for (int arcIndex : nodes[terminalNodeIndex].inArcs)
  {
    const VRPTWArc& arc = arcs[arcIndex];
    int nodeIndex = arc.fromNodeIndex;
    double shortestPathLength = nodes[nodeIndex].shortestPathDistance + arc.cijPi;
    nodeIndexDistanceToTerminal.push_back(std::make_pair(arcIndex, shortestPathLength));
  }

  // order distances
  std::sort(nodeIndexDistanceToTerminal.begin(), nodeIndexDistanceToTerminal.end(), sort_by_second1);

  // check if currentNumPath-th shortest path is in the nodesToUpdate set
  int newArcIndex = nodeIndexDistanceToTerminal[currentNumPaths].first;
  int newNodeIndex = arcs[newArcIndex].fromNodeIndex;
  if (nodesToUpdate.find(newNodeIndex) != nodesToUpdate.end())
  {
    return -1;
  }
  else
  {
    return newArcIndex;
  }
}

void VRPTWDecisionDiagram::findMultiPath(int newPathArcIndex, const std::vector<int>& treeByParentArcs, std::vector<int>& newShortestPathByArcs)
{
  newShortestPathByArcs.clear();
  newShortestPathByArcs.push_back(newPathArcIndex);
  int nextNodeIndex = arcs[newPathArcIndex].fromNodeIndex;
  while (nextNodeIndex != rootNodeIndex)
  {
    int arcIndex = treeByParentArcs[nextNodeIndex];
    newShortestPathByArcs.push_back(arcIndex);
    nextNodeIndex = arcs[arcIndex].fromNodeIndex;
  }

  std::reverse(newShortestPathByArcs.begin(), newShortestPathByArcs.end());
}

void VRPTWDecisionDiagram::clipPermanentArc(int arcIndex)
{
  const VRPTWArc& arc = arcs[arcIndex];
  auto fromNodeOutArcIt = std::find(nodes[arc.fromNodeIndex].outArcs.begin(), nodes[arc.fromNodeIndex].outArcs.end(), arcIndex);
  nodes[arc.fromNodeIndex].outArcs.erase(fromNodeOutArcIt);
  auto toNodeInArcIt = std::find(nodes[arc.toNodeIndex].inArcs.begin(), nodes[arc.toNodeIndex].inArcs.end(), arcIndex);
  nodes[arc.toNodeIndex].inArcs.erase(toNodeInArcIt);
}

void VRPTWDecisionDiagram::dijkstraWithBatchProc(std::vector<int>& treeByParentArcs, std::vector<std::vector<int>>& treeByChildArcs, std::set<int>& nodesToUpdate)
{
  // get temporary distances
  std::vector<std::pair<int,double>> tempDistances(nodes.size(), std::make_pair(-1,INF));
  DBG(std::cout << "nodes to update: ";)
  for (int nodeIndex : nodesToUpdate)
  {
    DBG(std::cout << nodeIndex << ", ";)
    const VRPTWNode& node = nodes[nodeIndex];
    for (int arcIndex : node.inArcs)
    {
      // only relax outside arcs at first
      const VRPTWArc& arc = arcs[arcIndex];
      if (nodesToUpdate.find(arc.fromNodeIndex) == nodesToUpdate.end())
      {
        double possibleTempDistance = nodes[arc.fromNodeIndex].shortestPathDistance + arc.cijPi;
        if (tempDistances[nodeIndex].second > possibleTempDistance)
        {
          tempDistances[nodeIndex] = std::make_pair(arcIndex, possibleTempDistance);
        }
      }
    }
  }
  DBG(std::cout << std::endl;)

  // setup queue for dijkstra
  std::priority_queue<pqKeyType, std::vector<pqKeyType>, std::greater<pqKeyType>> pq;
  std::set<int> permanent;

  // only add appropriate nodes to queue
  for (int nodeIndex : nodesToUpdate)
  {
    if (nodeIndex == terminalNodeIndex)
    {
      continue;
    }

    // NOTE: (else if) added edge case - if node parent is terminal and can still be reached
    // NOTE CONT: then it also needs to be added :)
    int parentNodeIndex = arcs[treeByParentArcs[nodeIndex]].fromNodeIndex;
    if (tempDistances[nodeIndex].second < tempDistances[parentNodeIndex].second)
    {
      pq.push(std::make_pair(tempDistances[nodeIndex].second, std::make_pair(nodeIndex, tempDistances[nodeIndex].first)));
    }
    else if (tempDistances[nodeIndex].first != -1)
    {
      if (parentNodeIndex == terminalNodeIndex)
      {
        pq.push(std::make_pair(tempDistances[nodeIndex].second, std::make_pair(nodeIndex, tempDistances[nodeIndex].first)));
      }
    }
  }

  // start running dijkstras
  while (!pq.empty())
  {
    // extract min distance vertex u from pq
    pqKeyType u = pq.top();
    pq.pop();
    int nodeIndex = u.second.first;
    int fromArcIndex = u.second.second;
    if (permanent.find(nodeIndex) == permanent.end())
    {
      permanent.insert(nodeIndex);
      treeByParentArcs[nodeIndex] = fromArcIndex;
      nodes[nodeIndex].shortestPathDistance = u.first;
    }
    else
    {
      continue;
    }

    // relax arcs out of popped node and downhill nodes, update tree too
    std::vector<int> downhillNodes;
    downhillNodes.push_back(nodeIndex);
    for (int index=0; index<downhillNodes.size(); ++index)
    {
      int downhillNodeIndex = downhillNodes[index];
      // relax downhill children
      for (int arcIndex : treeByChildArcs[downhillNodeIndex])
      {
        int toNodeIndex = arcs[arcIndex].toNodeIndex;

        // label downhill nodes
        if (nodesToUpdate.find(toNodeIndex) != nodesToUpdate.end())
        {
          downhillNodes.push_back(toNodeIndex);
          if (permanent.find(toNodeIndex) == permanent.end())
          {
            // NOTE: make sure to label downhill across branches too!
            permanent.insert(toNodeIndex);
            treeByParentArcs[toNodeIndex] = arcIndex;
            nodes[toNodeIndex].shortestPathDistance = u.first;
            for (int arcIndex : nodes[toNodeIndex].outArcs)
            {
              int relaxNodeIndex = arcs[arcIndex].toNodeIndex;
              // NOTE: funky thing with child/parent from residual update...
              // NOTE CONT: we don't want terminal node in queue b/c has "child" nodes
              if (relaxNodeIndex == terminalNodeIndex)
              {
                continue;
              }
              if ((nodesToUpdate.find(relaxNodeIndex) != nodesToUpdate.end()) && (permanent.find(relaxNodeIndex) == permanent.end()))
              {
                double newPossibleDistance = u.first + arcs[arcIndex].cijPi;
                pq.push(std::make_pair(newPossibleDistance, std::make_pair(relaxNodeIndex, arcIndex)));
              }
            }
          }
        }
      }
    }

    // relax across 0-tree arcs
    for (int arcIndex : nodes[nodeIndex].outArcs)
    {
      int toNodeIndex = arcs[arcIndex].toNodeIndex;
      // NOTE: funky thing with child/parent from residual update...
      // NOTE CONT: we don't want terminal node in queue b/c has "child" nodes
      if (toNodeIndex == terminalNodeIndex)
      {
        continue;
      }
      if ((nodesToUpdate.find(toNodeIndex) != nodesToUpdate.end()) && (permanent.find(toNodeIndex) == permanent.end()))
      {
        double newPossibleDistance = u.first + arcs[arcIndex].cijPi;
        pq.push(std::make_pair(newPossibleDistance, std::make_pair(toNodeIndex, arcIndex)));
      }
    }
  }

  // some nodes become unreachable
  for (int nodeIndex : nodesToUpdate)
  {
    if (permanent.find(nodeIndex) == permanent.end())
    {
      nodes[nodeIndex].shortestPathDistance = INF;
    }
  }

  // terminal node
  VRPTWNode& node = nodes[terminalNodeIndex];
  node.shortestPathDistance = INF;
  for (int arcIndex : node.inArcs)
  {
    const VRPTWArc& arc = arcs[arcIndex];
    double possibleTempDistance = nodes[arc.fromNodeIndex].shortestPathDistance + arc.cijPi;
    if (node.shortestPathDistance > possibleTempDistance)
    {
      node.shortestPathDistance = possibleTempDistance;
      treeByParentArcs[terminalNodeIndex] = arcIndex;
    }
  }
  DBG(print();)
}

double VRPTWDecisionDiagram::solveMinCostFlowModelWang(const Dual& dual, std::vector<std::vector<int>>& shortestPathsByArc, bool& isDualFeasible, double& shortestPathLength, double& longestShortestPathLength)
{
  // start potentials and flows at 0
  for (VRPTWNode& node : nodes)
  {
    node.potential = 0.0;
    node.shortestPathDistance = INF;
  }
  nodes[rootNodeIndex].shortestPathDistance = 0.0;
  for (VRPTWArc& arc : arcs)
  {
    arc.decompositionFlow = 0.0;
    arc.heuristicFlow = 0.0;
  }

  // setup some data structures (follow Wang,Wang,Wang NeurIps paper)
  std::vector<int> treeByParentArcs;
  treeByParentArcs.resize(nodes.size());
  std::vector<int> shortestPathByArc;
  setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(dual, LPSolveType::LAGSolver);
  shortestPathLength = computeShortestPathBFSWang(treeByParentArcs, shortestPathByArc, longestShortestPathLength);
  // NOTE: when fixing arcs, if there is no path to terminal, no paths exist
  if (shortestPathLength == INF)
  {
    return 0.0;
  }

  // shortest path tells us dual feasibility
  if (shortestPathLength < 0)
  {
    isDualFeasible = false;
    DBG(std::cout << "dual not feasible: " << shortestPathLength << std::endl;)
  }
  else
  {
    isDualFeasible = true;
    if (vrptw.fixedNumPaths != FixedNumPaths::FIXED_NUM_PATHS)
    {
      return 0.0;
    }
  }

  std::vector<std::vector<int>> treeByChildArcs;
  treeByChildArcs.resize(nodes.size());
  getTreeByChildArcsFromTreeByParentArcs(treeByParentArcs, treeByChildArcs);

  // while negative s-t path exists, send flow and make updates
  // take at most k (number of trucks) negative s-t paths
  std::set<int> nodesToUpdate;
  std::vector<std::vector<int>> allCurrentPaths;
  allCurrentPaths.push_back(shortestPathByArc);
  shortestPathsByArc.push_back(shortestPathByArc);
  std::set<int> clippedArcs;
  int currentNumPaths = 1;
  while (true)
  {
    // identify nodes 4 update
    int branchNodeIndex = arcs[shortestPathByArc[0]].toNodeIndex;
    identifyNodesForUpdate(branchNodeIndex, treeByChildArcs, nodesToUpdate);

    // find multi path
    int newPathArcIndex = findMultiPathNode(nodesToUpdate, currentNumPaths);
    if ((newPathArcIndex != -1) && (nodesToUpdate.find(arcs[newPathArcIndex].fromNodeIndex) == nodesToUpdate.end()) && !((vrptw.fixedNumPaths == FixedNumPaths::FIXED_NUM_PATHS) && (shortestPathsByArc.size() == vrptw.numVehicles)))
    {
      DBG(std::cout << "find multi path" << std::endl;)
      findMultiPath(newPathArcIndex, treeByParentArcs, shortestPathByArc);
      double currRouteCost = evaluateRouteCost(shortestPathByArc);
      DBG(std::cout << "r: " << currRouteCost << std::endl;)
      if ((currRouteCost >= -0.000000001) && (vrptw.fixedNumPaths == FixedNumPaths::FLEXIBLE_NUM_PATHS))
      {
        break;
      }
      else
      {
        allCurrentPaths.push_back(shortestPathByArc);
        shortestPathsByArc.push_back(shortestPathByArc);
        ++currentNumPaths;
      }
    }
    else
    {
      // convert edge cost and shortest distances to residuals
      for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
      {
        VRPTWArc& arc = arcs[arcIndex];
        arc.cijPi = arc.cijPi + nodes[arc.fromNodeIndex].shortestPathDistance - nodes[arc.toNodeIndex].shortestPathDistance;
      }
      for (int nodeIndex=0; nodeIndex<nodes.size(); ++nodeIndex)
      {
        nodes[nodeIndex].shortestPathDistance = 0.0;
      }

      // residual graph and update tree
      for (auto pathByArc : allCurrentPaths)
      {
        updateResidualGraphWang(pathByArc, treeByParentArcs);
      }
      getTreeByChildArcsFromTreeByParentArcs(treeByParentArcs, treeByChildArcs);
      // remove terminal as child
      int terminalParentIndex = arcs[treeByParentArcs[terminalNodeIndex]].fromNodeIndex;
      auto childArcsIt = std::find(treeByChildArcs[terminalParentIndex].begin(), treeByChildArcs[terminalParentIndex].end(), treeByParentArcs[terminalNodeIndex]);
      treeByChildArcs[terminalParentIndex].erase(childArcsIt);

      // clip permanent edges
      for (auto pathByArc : allCurrentPaths)
      {
        // clip the reverse because we updated residual
        clipPermanentArc(arcReverseArc[pathByArc.front()]);
        clipPermanentArc(arcReverseArc[pathByArc.back()]);
        clippedArcs.insert(pathByArc.front());
        clippedArcs.insert(pathByArc.back());
      }

      // fixed number of paths reached - stop trying for more
      if ((vrptw.fixedNumPaths == FixedNumPaths::FIXED_NUM_PATHS) && (vrptw.numVehicles == shortestPathsByArc.size()))
      {
        break;
      }

      // dijkstra with batch proc
      DBG(std::cout << "batch proc" << std::endl;)
      dijkstraWithBatchProc(treeByParentArcs, treeByChildArcs, nodesToUpdate);
      DBG(std::cout << "finish batch proc" << std::endl;)
      getTreeByChildArcsFromTreeByParentArcs(treeByParentArcs, treeByChildArcs);

      // NOTE: when fixing arcs, may not be able to find terminal node inArcs
      if (nodes[terminalNodeIndex].shortestPathDistance == INF)
      {
        break;
      }

      // update new shortest path
      shortestPathByArc.clear();
      int currentArcIndex = treeByParentArcs[terminalNodeIndex];
      int currentNodeIndex = arcs[currentArcIndex].fromNodeIndex;
      shortestPathByArc.push_back(currentArcIndex);
      while (true)
      {
        currentArcIndex = treeByParentArcs[currentNodeIndex];
        shortestPathByArc.push_back(currentArcIndex);
        currentNodeIndex = arcs[currentArcIndex].fromNodeIndex;
        if (currentNodeIndex == rootNodeIndex)
        {
          std::reverse(shortestPathByArc.begin(), shortestPathByArc.end());
          break;
        }
      }
      DBG(std::cout << "finish shortest path update" << std::endl;)

      double currRouteCost = evaluateRouteCost(shortestPathByArc);
      DBG(std::cout << "r: " << currRouteCost << std::endl;)
      if ((currRouteCost >= -0.000000001) && (vrptw.fixedNumPaths == FixedNumPaths::FLEXIBLE_NUM_PATHS))
      {
        break;
      }

      nodesToUpdate.clear();
      allCurrentPaths.clear();
      allCurrentPaths.push_back(shortestPathByArc);
      shortestPathsByArc.push_back(shortestPathByArc);
      currentNumPaths = 1;
    }
  }

  // create solution from reverse arcs and clean up residual graph backward arcs
  /*
  for (auto path : shortestPathsByArc)
  {
    std::cout << "path by loc: ";
    for (int arcIndex : path)
    {
      std::cout << arcs[arcIndex].location << ", ";
    }
    std::cout << std::endl;
    std::cout << "by arc: ";
    for (int arcIndex : path)
    {
      std::cout << arcIndex << ", ";
    }
    std::cout << std::endl;
    std::cout << "by is reverse: ";
    for (int arcIndex : path)
    {
      std::cout << arcs[arcIndex].isReverseArc << ", ";
    }
    std::cout << std::endl;
    std::cout << "by reverse arcs: ";
    for (int arcIndex : path)
    {
      if (arcReverseArc.find(arcIndex) != arcReverseArc.end())
      {
        std::cout << arcReverseArc[arcIndex] << ", ";
      }
      else
      {
        std::cout << "-" << ", ";
      }
    }
    std::cout << std::endl;
  }
  */
  DBG(std::cout << "finished running, now create solution" << std::endl;)
  return createSolutionFromReverseArcsAndResetWang(clippedArcs, shortestPathsByArc);
};

bool VRPTWDecisionDiagram::checkFeasibleDual(const Dual& dual, LPSolveType solveType)
{
  // get shortest path
  setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(dual, solveType);
  std::vector<int> emptyRoute;
  double shortestPathDistance = computeShortestPathBFS(ShortestPathMode::SHORTEST_PATH, emptyRoute);

  // dual is infeasible if there exists a negative reduced cost route
  if (shortestPathDistance < 0)
  {
    return false;
  }

  return true;
};

void VRPTWDecisionDiagram::getLocationsOnArcPaths(const std::vector<std::vector<int>>& shortestPathsByArc, std::set<int>& locations)
{
  for (auto path : shortestPathsByArc)
  {
    for (int arcIndex : path)
    {
      int loc = arcs[arcIndex].location;
      if (loc != 0)
      {
        locations.insert(loc);
      }
    }
  }
};

bool VRPTWDecisionDiagram::checkC141SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({5,160,337,330,143,215,46,147,153,89,166,258,0});
  routesByLocation.push_back({17,326,214,355,117,247,162,320,168,283,289,158,104,227,0});
  routesByLocation.push_back({23,371,38,333,111,194,363,251,316,276,142,0});
  routesByLocation.push_back({33,255,29,157,212,169,237,28,0});
  routesByLocation.push_back({49,298,87,70,302,201,25,96,353,0});
  routesByLocation.push_back({51,181,389,1,116,127,245,242,259,0});
  routesByLocation.push_back({56,141,252,370,274,138,170,4,190,275,0});
  routesByLocation.push_back({63,213,121,233,173,161,50,383,261,131,0});
  routesByLocation.push_back({73,265,380,85,328,177,95,30,0});
  routesByLocation.push_back({80,315,221,284,351,58,356,134,107,0});
  routesByLocation.push_back({92,266,165,308,146,268,79,0});
  routesByLocation.push_back({94,54,180,179,386,379,312,246,81,197,84,0});
  routesByLocation.push_back({97,375,149,101,7,257,300,13,15,0});
  routesByLocation.push_back({102,69,222,159,304,71,113,211,59,253,0});
  routesByLocation.push_back({105,155,395,361,272,306,377,124,2,347,376,0});
  routesByLocation.push_back({114,44,348,235,388,281,129,263,176,137,384,0});
  routesByLocation.push_back({126,115,331,398,66,277,325,368,37,20,0});
  routesByLocation.push_back({186,154,381,202,342,309,118,21,346,45,0});
  routesByLocation.push_back({189,174,390,271,244,57,48,345,98,250,234,0});
  routesByLocation.push_back({267,364,340,373,167,148,378,164,367,391,0});
  routesByLocation.push_back({286,88,217,209,185,10,152,336,0});
  routesByLocation.push_back({287,31,362,108,248,171,243,239,192,128,27,0});
  routesByLocation.push_back({290,182,236,61,218,11,369,110,35,223,240,0});
  routesByLocation.push_back({291,52,16,314,136,317,22,399,132,0});
  routesByLocation.push_back({292,226,109,172,26,350,282,352,40,183,123,75,228,0});
  routesByLocation.push_back({301,139,280,392,36,191,349,344,91,130,372,0});
  routesByLocation.push_back({303,144,24,93,229,204,295,103,220,122,0});
  routesByLocation.push_back({322,151,270,193,19,55,120,294,0});
  routesByLocation.push_back({324,156,133,210,382,338,256,205,0});
  routesByLocation.push_back({332,310,339,262,400,8,319,260,43,230,41,0});
  routesByLocation.push_back({334,216,358,195,285,18,238,297,329,0});
  routesByLocation.push_back({335,385,199,82,296,359,254,293,60,53,313,86,0});
  routesByLocation.push_back({343,269,64,208,100,145,318,106,163,0});
  routesByLocation.push_back({357,77,74,394,68,14,288,6,178,327,0});
  routesByLocation.push_back({360,198,12,9,47,299,305,219,249,0});
  routesByLocation.push_back({366,150,83,354,206,99,184,175,3,0});
  routesByLocation.push_back({374,32,62,39,67,65,396,307,264,0});
  routesByLocation.push_back({387,278,78,231,188,125,225,323,34,224,0});
  routesByLocation.push_back({393,112,135,140,279,90,196,321,76,311,241,0});
  routesByLocation.push_back({397,72,203,365,273,341,207,187,200,119,232,42,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR C141: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  return true;
}

bool VRPTWDecisionDiagram::checkLC121SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({32,171,65,86,115,94,51,174,136,189,0});
  routesByLocation.push_back({177,3,88,8,186,127,98,157,137,183,0});
  routesByLocation.push_back({21,23,182,75,163,194,145,195,52,92,0});
  routesByLocation.push_back({161,104,18,54,185,132,7,181,117,49,0});
  routesByLocation.push_back({60,211,82,180,84,191,125,4,72,17,0});
  routesByLocation.push_back({148,103,197,203,124,141,69,200,0});
  routesByLocation.push_back({170,134,50,156,112,168,79,205,29,87,42,123,0});
  routesByLocation.push_back({114,159,38,150,22,151,16,140,204,187,142,111,63,56,0});
  routesByLocation.push_back({190,5,10,193,46,128,106,167,207,34,95,158,0});
  routesByLocation.push_back({57,118,83,143,176,36,206,33,121,165,188,108,0});
  routesByLocation.push_back({93,55,135,58,202,184,199,37,81,138,0});
  routesByLocation.push_back({133,48,26,152,40,153,169,89,105,15,59,198,0});
  routesByLocation.push_back({164,210,66,147,160,47,91,70,0});
  routesByLocation.push_back({101,144,119,166,35,126,71,9,1,99,53,201,0});
  routesByLocation.push_back({30,120,19,192,196,97,14,96,130,28,74,149,0});
  routesByLocation.push_back({20,41,85,80,31,25,172,77,110,162,0});
  routesByLocation.push_back({73,116,12,129,11,6,122,139,0});
  routesByLocation.push_back({62,131,44,102,146,208,68,76,0});
  routesByLocation.push_back({45,178,27,173,154,209,24,61,100,64,179,109,0});
  routesByLocation.push_back({113,155,78,175,13,43,2,90,67,39,107,212,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR LC121: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  std::cout << "clear" << std::endl;
  return true;
}

bool VRPTWDecisionDiagram::checkLRC121SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({128,168,53,151,160,109,189,212,125,17,0});
  routesByLocation.push_back({163,28,205,6,49,14,132,78,171,130,187,150,0});
  routesByLocation.push_back({92,195,31,52,74,101,192,50,33,157,0});
  routesByLocation.push_back({172,197,5,170,73,135,2,76,97,64,0});
  routesByLocation.push_back({153,47,200,103,90,79,204,96,34,175,0});
  routesByLocation.push_back({162,63,149,26,62,159,20,46,138,186,120,39,0});
  routesByLocation.push_back({24,124,43,180,136,3,110,178,93,193,19,146,0});
  routesByLocation.push_back({142,152,167,61,198,91,202,145,116,148,0});
  routesByLocation.push_back({88,155,144,71,137,123,209,13,8,1,60,75,0});
  routesByLocation.push_back({42,181,35,164,104,25,10,211,77,122,141,185,0});
  routesByLocation.push_back({140,54,89,87,121,208,115,11,30,67,0});
  routesByLocation.push_back({165,85,45,56,108,38,32,203,23,117,55,154,0});
  routesByLocation.push_back({4,182,57,133,147,40,16,22,106,70,0});
  routesByLocation.push_back({84,199,127,118,105,69,176,113,191,58,134,174,0});
  routesByLocation.push_back({100,9,12,59,143,173,188,29,139,21,166,82,0});
  routesByLocation.push_back({169,98,119,161,72,129,41,194,210,158,0});
  routesByLocation.push_back({131,156,102,196,183,190,111,44,95,83,201,36,0});
  routesByLocation.push_back({48,15,179,7,68,18,184,27,51,80,126,206,0});
  routesByLocation.push_back({37,177,65,99,94,112,114,86,81,107,66,207,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR LRC121: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  std::cout << "clear" << std::endl;
  return true;
}

bool VRPTWDecisionDiagram::checkAn32k5SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({21,31,19,17,13,7,26,0});
  routesByLocation.push_back({12,1,16,30,0});
  routesByLocation.push_back({27,24,0});
  routesByLocation.push_back({29,18,8,9,22,15,10,25,5,20,0});
  routesByLocation.push_back({14,28,11,4,23,3,2,6,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR An32k5: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  std::cout << "clear" << std::endl;
  return true;
}

bool VRPTWDecisionDiagram::checkAn37k5SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({22,13,10,6,5,33,4,7,0});
  routesByLocation.push_back({1,12,2,19,20,23,14,17,0});
  routesByLocation.push_back({36,29,32,28,31,30,15,0});
  routesByLocation.push_back({3,24,9,11,27,8,25,35,18,26,34,0});
  routesByLocation.push_back({21,16,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR An37k5: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  std::cout << "clear" << std::endl;
  return true;
}

bool VRPTWDecisionDiagram::checkAn33k6SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  //routesByLocation.push_back({21,31,19,17,13,7,26,0});
  routesByLocation.push_back({5,2,20,15,9,3,8,4,0});
  routesByLocation.push_back({31,24,23,26,22,0});
  routesByLocation.push_back({17,11,29,19,7,0});
  routesByLocation.push_back({10,12,21,0});
  routesByLocation.push_back({28,27,30,16,25,32,0});
  routesByLocation.push_back({13,6,18,1,14,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR An33k6: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  std::cout << "clear" << std::endl;
  return true;
}

bool VRPTWDecisionDiagram::checkAn36k5SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({9,6,3,4,19,31,12,0});
  routesByLocation.push_back({28,14,34,23,2,35,8,15,0});
  routesByLocation.push_back({16,11,24,27,25,5,20,0});
  routesByLocation.push_back({10,7,26,0});
  routesByLocation.push_back({1,22,32,13,17,30,29,33,18,21,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR An36k5: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  std::cout << "clear" << std::endl;
  return true;
}

bool VRPTWDecisionDiagram::checkAn48k7SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({34,9,24,4,11,42,15,27,45,0});
  routesByLocation.push_back({32,36,38,19,25,22,6,0});
  routesByLocation.push_back({40,7,8,39,26,20,3,37,0});
  routesByLocation.push_back({41,2,10,47,17,14,0});
  routesByLocation.push_back({44,35,18,0});
  routesByLocation.push_back({28,29,21,30,13,46,33,16,0});
  routesByLocation.push_back({23,43,31,1,5,12,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR An48k7: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  std::cout << "clear" << std::endl;
  return true;
}

bool VRPTWDecisionDiagram::checkXn181k23SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({114,67,130,73,124,144,165,104,0});
  routesByLocation.push_back({8,76,18,159,171,96,80,58,0});
  routesByLocation.push_back({88,78,93,3,102,21,98,176,0});
  routesByLocation.push_back({49,137,107,20,161,109,68,19,0});
  routesByLocation.push_back({14,36,158,79,30,9,128,60,0});
  routesByLocation.push_back({52,51,29,27,122,157,47,143,0});
  routesByLocation.push_back({142,11,42,84,26,90,163,95,0});
  routesByLocation.push_back({132,140,63,85,89,2,149,152,0});
  routesByLocation.push_back({116,92,131,169,61,65,32,44,0});
  routesByLocation.push_back({45,148,180,145,0});
  routesByLocation.push_back({59,106,154,38,129,136,77,13,0});
  routesByLocation.push_back({133,173,119,28,69,100,56,99,0});
  routesByLocation.push_back({139,160,172,141,81,66,25,72,0});
  routesByLocation.push_back({1,16,50,94,146,166,74,46,0});
  routesByLocation.push_back({123,174,48,6,126,120,101,112,0});
  routesByLocation.push_back({153,82,103,125,83,134,7,91,0});
  routesByLocation.push_back({31,57,110,175,162,168,41,22,0});
  routesByLocation.push_back({75,113,4,53,33,151,35,117,0});
  routesByLocation.push_back({39,167,87,147,37,164,15,127,0});
  routesByLocation.push_back({178,156,177,40,150,97,34,118,0});
  routesByLocation.push_back({17,71,24,179,115,138,155,170,0});
  routesByLocation.push_back({55,135,5,105,64,86,62,23,0});
  routesByLocation.push_back({10,43,111,54,70,121,12,108,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR Xn181k23: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  std::cout << "clear" << std::endl;
  return true;
}

bool VRPTWDecisionDiagram::checkLC145SolutionCost()
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({290, 182, 236, 61, 218, 11, 369, 406, 110, 35, 223, 240,0});
  routesByLocation.push_back({49, 298, 87, 70, 302, 201, 25, 96, 413, 353,0});
  routesByLocation.push_back({5, 160, 337, 330, 143, 215, 46, 147, 153, 89, 166, 258,0});
  routesByLocation.push_back({17, 326, 214, 355, 117, 247, 162, 320, 168, 283, 289, 158, 104, 227,0});
  routesByLocation.push_back({23, 371, 38, 333, 111, 420, 194, 363, 251, 316, 276, 142,0});
  routesByLocation.push_back({33, 255, 29, 157, 212, 169, 237, 28,0});
  routesByLocation.push_back({51, 403, 181, 389, 1, 116, 127, 245, 242, 259,0});
  routesByLocation.push_back({56, 141, 252, 370, 274, 138, 170, 4, 190, 275,0});
  routesByLocation.push_back({63, 213, 121, 233, 173, 161, 50, 383, 261, 131,0});
  routesByLocation.push_back({73, 265, 380, 85, 328, 177, 95, 30,0});
  routesByLocation.push_back({80, 315, 221, 284, 351, 58, 356, 134, 107, 401,0});
  routesByLocation.push_back({92, 266, 165, 308, 146, 268, 418, 79,0});
  routesByLocation.push_back({94, 54, 180, 179, 386, 379, 312, 414, 246, 81, 197, 84,0});
  routesByLocation.push_back({97, 422, 375, 149, 101, 7, 257, 300, 13, 15,0});
  routesByLocation.push_back({102, 69, 222, 159, 304, 71, 113, 211, 59, 253,0});
  routesByLocation.push_back({105, 155, 395, 361, 272, 306, 412, 377, 124, 2, 347, 376,0});
  routesByLocation.push_back({114, 44, 348, 235, 388, 281, 129, 263, 176, 137, 384, 409,0});
  routesByLocation.push_back({126, 115, 331, 398, 66, 277, 325, 368, 37, 20,0});
  routesByLocation.push_back({186, 154, 381, 202, 342, 309, 118, 21, 346, 45,0});
  routesByLocation.push_back({189, 174, 390, 271, 244, 57, 402, 48, 345, 98, 250, 234,0});
  routesByLocation.push_back({267, 364, 340, 373, 167, 148, 378, 164, 367, 391,0});
  routesByLocation.push_back({286, 88, 217, 209, 185, 10, 152, 336,0});
  routesByLocation.push_back({287, 31, 362, 108, 248, 171, 243, 239, 192, 128, 27, 408,0});
  routesByLocation.push_back({291, 52, 407, 16, 314, 136, 317, 22, 399, 132, 311, 415,0});
  routesByLocation.push_back({292, 226, 109, 172, 26, 350, 282, 352, 40, 405, 183, 123, 75, 228,0});
  routesByLocation.push_back({301, 139, 280, 392, 36, 419, 191, 349, 344, 91, 130, 372,0});
  routesByLocation.push_back({303, 144, 24, 93, 229, 204, 295, 103, 220, 122,0});
  routesByLocation.push_back({322, 151, 270, 193, 19, 55, 120, 294,0});
  routesByLocation.push_back({324, 156, 133, 210, 382, 338, 256, 205,0});
  routesByLocation.push_back({332, 310, 339, 262, 400, 8, 417, 319, 260, 43, 230, 41,0});
  routesByLocation.push_back({334, 216, 358, 195, 285, 18, 238, 297, 329, 411, 163, 410,0});
  routesByLocation.push_back({335, 385, 199, 82, 296, 359, 254, 293, 60, 53, 313, 86,0});
  routesByLocation.push_back({343, 269, 64, 208, 100, 145, 318, 106,0});
  routesByLocation.push_back({357, 77, 74, 394, 68, 14, 288, 6, 178, 327,0});
  routesByLocation.push_back({360, 198, 12, 9, 47, 299, 305, 219, 416, 249,0});
  routesByLocation.push_back({366, 150, 83, 354, 206, 99, 184, 175, 3, 421,0});
  routesByLocation.push_back({374, 32, 62, 39, 67, 404, 65, 396, 307, 264,0});
  routesByLocation.push_back({387, 278, 78, 231, 188, 125, 225, 323, 34, 224,0});
  routesByLocation.push_back({393, 112, 135, 140, 279, 90, 196, 321, 76, 241,0});
  routesByLocation.push_back({397, 72, 203, 365, 273, 341, 207, 187, 200, 119, 232, 42,0});
  double totalCost = 0.0;
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    double cost = vrptw.evaluateRouteDistance(routesByLocation[routeIndex]);
  std::cout << "cost: " << cost << std::endl;
    totalCost += cost;
  }
  std::cout << "total cost: " << totalCost << std::endl;

  return true;
}

bool VRPTWDecisionDiagram::checkSinglePathPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({14,8,11,21,5,20,19,22,24,25,10,9,12,13,28,2,1,16,15,23,17,4,3,26,6,7,27,18,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    std::vector<int> updatedRouteArcs;
    if (!doesRouteExistByLocations(routesByLocation[routeIndex], updatedRouteArcs))
    {
      std::cout << "ERROR RC_202.3.txt: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  return true;
}

void VRPTWDecisionDiagram::printCuts() const
{
  for (int index=0; index<srcCuts.size(); ++index)
  {
    std::cout << "index: " << index << std::endl;
    auto cutSet = srcCuts[index];
    std::cout << "locations: " << std::endl;
    for (int v : cutSet)
    {
      std::cout << v << ",";
    }
    std::cout << std::endl;
    
    auto arcIndices = srcCutSeparatedArcs[index];
    std::cout << "arc indices: " << std::endl;
    for (int v : arcIndices)
    {
      std::cout << v << ",";
    }
    std::cout << std::endl;

    auto coeffs = srcCutSeparatedCoeffs[index];
    std::cout << "coeffs: " << std::endl;
    for (int v : coeffs)
    {
      std::cout << v << ",";
    }
    std::cout << std::endl;
  }
}

void VRPTWDecisionDiagram::printFixedArcs() const
{
  std::cout << "fixed arcs: " << std::endl;
  for (int arcIndex : fixedArcs)
  {
    std::cout << arcIndex << ",";
  }
  std::cout << std::endl;
}

void VRPTWDecisionDiagram::convertSolutionForVRPTWSep(std::vector<int>& edgeTail,
                                             std::vector<int>& edgeHead,
                                             std::vector<double>& edgeFlow)
{
  std::map<std::pair<int,int>,double> edgeFlows;
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    const VRPTWArc& arc = arcs[arcIndex];
    if (arc.heuristicFlow > 0.0001)
    {
      int currLoc = (nodes[arc.fromNodeIndex].state.lastVisited == 0) ? vrptw.numLocations : nodes[arc.fromNodeIndex].state.lastVisited;
      int nextLoc = (arc.location == 0) ? vrptw.numLocations : arc.location;

      // special src paths have 0,0 start
      if ((currLoc == vrptw.numLocations) && (currLoc == nextLoc))
      {
        continue;
      }

      if (currLoc < nextLoc)
      {
        edgeFlows[std::make_pair(currLoc,nextLoc)] = edgeFlows[std::make_pair(currLoc,nextLoc)] + arc.heuristicFlow;
      }
      else
      {
        edgeFlows[std::make_pair(nextLoc,currLoc)] = edgeFlows[std::make_pair(nextLoc,currLoc)] + arc.heuristicFlow;
      }
    }
  }

  edgeTail.push_back(0);
  edgeHead.push_back(0);
  edgeFlow.push_back(0);
  for (auto edge : edgeFlows)
  {
    edgeTail.push_back(edge.first.first);
    edgeHead.push_back(edge.first.second);
    edgeFlow.push_back(edge.second);
  }
};

bool VRPTWDecisionDiagram::calculateRCCCoeff(int arcIndex, int rccIndex, double& coeff)
{
  auto cutSet = capCutSets[rccIndex];
  auto rccType = rccTypes[rccIndex];
  if (rccType == RCCType::Type1)
  {
    if (!arcs[arcIndex].isReverseArc)
    {
      int fromLoc = nodes[arcs[arcIndex].fromNodeIndex].state.lastVisited;
      int toLoc = nodes[arcs[arcIndex].toNodeIndex].state.lastVisited;
      bool fromLocInSet = (std::find(cutSet.begin(), cutSet.end(), fromLoc) != cutSet.end());
      bool toLocInSet = (std::find(cutSet.begin(), cutSet.end(), toLoc) != cutSet.end());
      if (fromLocInSet && toLocInSet)
      {
        coeff = 1;
        return true;
      }
    }
  }
  else if (rccType == RCCType::Type3)
  {
    if (!arcs[arcIndex].isReverseArc)
    {
      int fromLoc = nodes[arcs[arcIndex].fromNodeIndex].state.lastVisited;
      int toLoc = nodes[arcs[arcIndex].toNodeIndex].state.lastVisited;
      bool fromLocInSet = (std::find(cutSet.begin(), cutSet.end(), fromLoc) != cutSet.end());
      bool toLocInSet = (std::find(cutSet.begin(), cutSet.end(), toLoc) != cutSet.end());
      if ((fromLoc != 0) && (toLoc != 0) && !fromLocInSet && !toLocInSet)
      {
        coeff = 1;
        return true;
      }

      if (((fromLoc == 0) && toLocInSet) || ((toLoc == 0) && fromLocInSet))
      {
        coeff = -0.5;
        return true;
      }
 
      if (((fromLoc == 0) && !toLocInSet) || ((toLoc == 0) && !fromLocInSet))
      {
        coeff = 0.5;
        return true;
      }
    }
  }

  coeff = 0;
  return false;
}

// RCCType1 is x(S:S) <= |S| - k(S)
// RCCType3 is x(Sbar:Sbar) + 0.5x({0}:Sbar) - 0.5x({0}:S) <= |Sbar| - k(S)
bool VRPTWDecisionDiagram::addCapCutSet(const std::vector<int>& cutSet, const std::vector<int>& sequenceArcs, double rhs, RCCType rccType, bool useScaling)
{
  // check existing cap cut sets
  int existingCapCutSetIndex = -1;
  for (int capCutSetIndex=0; capCutSetIndex<capCutSets.size(); ++capCutSetIndex)
  {
    auto existingCutSet = capCutSets[capCutSetIndex];
    if (existingCutSet == cutSet)
    {
      std::cout << "cut set exists: " << capCutSetIndex << ", add arcs" << std::endl;
      existingCapCutSetIndex = capCutSetIndex;
      break;
    }
  }

  // if none exists already, make a new one
  double scaling = 1.0;
  if (existingCapCutSetIndex == -1)
  {
    capCutSets.push_back(cutSet);
    capCutSetArcs.resize(capCutSets.size());
    capCutSetCoeffs.resize(capCutSets.size());
    if (useScaling)
    {
      if ((rhs > 0.0001) || (rhs < -0.0001))
      {
        scaling = 1.0 / std::abs(rhs);
      }
    }
    capCutSetsRHS.push_back(rhs * scaling);
    capCutSetsScaling.push_back(scaling);
    rccTypes.push_back(rccType);
  }

  // add arcs to cut sets
  int currCutSetIndex = static_cast<int>(capCutSets.size())-1;
  if (existingCapCutSetIndex != -1)
  {
    currCutSetIndex = existingCapCutSetIndex;
    scaling = capCutSetsScaling[currCutSetIndex];
  }
  if (sequenceArcs.empty())
  {
    for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
    {
      double coeff = 0;
      bool nonZeroCoeff = calculateRCCCoeff(arcIndex, currCutSetIndex, coeff);
      if (nonZeroCoeff)
      {
        if (std::find(capCutSetArcs[currCutSetIndex].begin(), capCutSetArcs[currCutSetIndex].end(), arcIndex) == capCutSetArcs[currCutSetIndex].end())
        {
          capCutSetArcs[currCutSetIndex].push_back(arcIndex);
          capCutSetCoeffs[currCutSetIndex].push_back(coeff * scaling);
        }
      }
    }
  }
  else
  {
    for (int arcIndex : sequenceArcs)
    {
      double coeff = 0;
      bool nonZeroCoeff = calculateRCCCoeff(arcIndex, currCutSetIndex, coeff);
      if (nonZeroCoeff)
      {
        if (std::find(capCutSetArcs[currCutSetIndex].begin(), capCutSetArcs[currCutSetIndex].end(), arcIndex) == capCutSetArcs[currCutSetIndex].end())
        {
          capCutSetArcs[currCutSetIndex].push_back(arcIndex);
          capCutSetCoeffs[currCutSetIndex].push_back(coeff * scaling);
        }
      }
    }
  }

  if (existingCapCutSetIndex != -1)
  {
    return true;
  }
  else
  {
    return false;
  }
};

void VRPTWDecisionDiagram::addCombCutTeeth(const std::vector<std::set<int>>& teeth)
{
  teeths.push_back(teeth);
  combCutArcs.resize(teeths.size());
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    if (!arcs[arcIndex].isReverseArc)
    {
      int fromLoc = nodes[arcs[arcIndex].fromNodeIndex].state.lastVisited;
      int toLoc = nodes[arcs[arcIndex].toNodeIndex].state.lastVisited;
      for (int toothIndex=0; toothIndex<teeth.size(); ++toothIndex)
      {
        auto tooth = teeth[toothIndex];
        bool fromLocInSet = (std::find(tooth.begin(), tooth.end(), fromLoc) != tooth.end());
        bool toLocInSet = (std::find(tooth.begin(), tooth.end(), toLoc) != tooth.end());
        if (((fromLocInSet && !toLocInSet) || (!fromLocInSet && toLocInSet)) && ((fromLoc != 0) || (toLoc != 0)))
        {
          combCutArcs[teeths.size()-1].push_back(arcIndex);
        }
      }
    }
  }
};

double VRPTWDecisionDiagram::getRCCCoeff(std::vector<int> route, int rccIndex)
{
  RCCType rccType = rccTypes[rccIndex];
  double scaling = capCutSetsScaling[rccIndex];
  if (rccType == RCCType::Type1)
  {
    double coeff = 0;
    int fromLoc = 0;
    auto cutSet = capCutSets[rccIndex];
    for (int toLoc : route)
    {
      bool fromLocInSet = (std::find(cutSet.begin(), cutSet.end(), fromLoc) != cutSet.end());
      bool toLocInSet = (std::find(cutSet.begin(), cutSet.end(), toLoc) != cutSet.end());
      if (fromLocInSet && toLocInSet)
      {
        coeff = coeff + 1;
      }

      fromLoc = toLoc;
    }

    return coeff * scaling;
  }
  else if (rccType == RCCType::Type3)
  {
    double coeff = 0;
    int fromLoc = 0;
    auto cutSet = capCutSets[rccIndex];
    for (int toLoc : route)
    {
      bool fromLocInSet = (std::find(cutSet.begin(), cutSet.end(), fromLoc) != cutSet.end());
      bool toLocInSet = (std::find(cutSet.begin(), cutSet.end(), toLoc) != cutSet.end());
      if ((fromLoc != 0) && (toLoc != 0) && !fromLocInSet && !toLocInSet)
      {
        coeff = coeff + 1;
      }

      if (((fromLoc == 0) && toLocInSet) || ((toLoc == 0) && fromLocInSet))
      {
        coeff = coeff - 0.5;
      }
 
      if (((fromLoc == 0) && !toLocInSet && (toLoc != 0)) || ((toLoc == 0) && !fromLocInSet && (fromLoc != 0)))
      {
        coeff = coeff + 0.5;
      }

      fromLoc = toLoc;
    }

    return coeff * scaling;
  }

  return 0;
};

int VRPTWDecisionDiagram::getNumTimesSetVisited(std::vector<int> route, const std::set<int>& set)
{
  int numTimesVisited = 0;
  for (int toLoc : route)
  {
    if (toLoc == 0)
    {
      continue;
    }

    if (set.find(toLoc) != set.end())
    {
      numTimesVisited = numTimesVisited + 1;
    }
  }

  return numTimesVisited;
}

// |S| = 3, p=1/2
// |S| = 4, p=2/3
// |S| = 5, p=1/3
// |S| = 5, p=1/2
int VRPTWDecisionDiagram::getSRCCoeff(int numTimesVisited, SRCType srcType)
{
  if (srcType == SRCType::SRC3)
  {
    return getSRC3Coeff(numTimesVisited);
  }
  else if (srcType == SRCType::SRC4)
  {
    return getSRC4Coeff(numTimesVisited);
  }
  else if (srcType == SRCType::SRC5V1)
  {
    return getSRC5V1Coeff(numTimesVisited);
  }
  else if (srcType == SRCType::SRC5V2)
  {
    return getSRC5V2Coeff(numTimesVisited);
  }

  return 0;
};
    
bool VRPTWDecisionDiagram::isArcIncrementalSRC(int numTimesVisited, SRCType srcType)
{
  if (numTimesVisited <= 1)
  {
    return false;
  }

  if (srcType == SRCType::SRC3)
  {
    return std::floor(numTimesVisited / 2.0) > std::floor((numTimesVisited - 1) / 2.0);
  }
  else if (srcType == SRCType::SRC4)
  {
    return std::floor(numTimesVisited * 2.0 / 3.0) > std::floor((numTimesVisited - 1) * 2.0 / 3.0);
  }
  else if (srcType == SRCType::SRC5V1)
  {
    return std::floor(numTimesVisited / 3.0) > std::floor((numTimesVisited - 1) / 3.0);
  }
  else if (srcType == SRCType::SRC5V2)
  {
    return std::floor(numTimesVisited / 2.0) > std::floor((numTimesVisited - 1) / 2.0);
  }

  return false;
}

int VRPTWDecisionDiagram::getSRC3Coeff(int numTimesVisited)
{
  return static_cast<int>(std::floor(numTimesVisited / 2.0));
};

int VRPTWDecisionDiagram::getSRC4Coeff(int numTimesVisited)
{
  return static_cast<int>(std::floor(numTimesVisited * 2.0 / 3.0));
};

int VRPTWDecisionDiagram::getSRC5V1Coeff(int numTimesVisited)
{
  return static_cast<int>(std::floor(numTimesVisited / 3.0));
};

int VRPTWDecisionDiagram::getSRC5V2Coeff(int numTimesVisited)
{
  return static_cast<int>(std::floor(numTimesVisited / 2.0));
};
