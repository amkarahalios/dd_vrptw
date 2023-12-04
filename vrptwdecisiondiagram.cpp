#include <bits/stdc++.h>
#include <cmath>

extern "C" {
#include "graph.h"
#include "cliquer.h"
}
#include "vrptwdecisiondiagram.h"
#include <ilcplex/ilocplex.h>

boolean record_clique_func(set_t s,graph_t *g,clique_options *opts);
boolean record_clique_func(set_t s,graph_t *g,clique_options *opts) {
/*
        if (clique_count>=clique_list_size) {
                clique_list=(setelement**)realloc(clique_list,(clique_list_size+512) *
                                    sizeof(set_t));
                clique_list_size+=512;
        }
        clique_list[clique_count]=set_duplicate(s);
        clique_count++;
*/
        return TRUE;
}

boolean clique_print_time(int level, int i, int n, int max,
                          double cputime, double realtime,
                          clique_options *opts) {
        static float prev_time=100;
        static int prev_i=100;
        static int prev_max=100;
        static int prev_level=0;
        static double total_time=-100;
        FILE *fp=opts->output;
        int j;
        double clique_time = ABS(prev_time-realtime);
        total_time = total_time + clique_time;
        if (total_time > 1)
        {
          total_time=0;
          return FALSE;
        }

        if (fp==NULL)
                fp=stdout;

        if (ABS(prev_time-realtime)>0.1 || i==n || i<prev_i || max!=prev_max ||
            level!=prev_level) {
                for (j=1; j<level; j++)
                        fprintf(fp,"  ");
                if (realtime-prev_time < 0.01 || i<=prev_i)
                        fprintf(fp,"%3d/%d (max %2d)  %2.2f s  "
                                "(0.00 s/round)\n",i,n,max,
                                realtime);
                else
                        fprintf(fp,"%3d/%d (max %2d)  %2.2f s  "
                                "(%2.2f s/round)\n",
                                i,n,max,realtime,
                                (realtime-prev_time)/(i-prev_i));
                prev_time=realtime;
                prev_i=i;
                prev_max=max;
                prev_level=level;
        }
        return TRUE;
}

void VRPTWNode::printForDD(int nodeIndex) const
{
  std::cout << nodeIndex << ": {";
  for (int loc : state.visited)
  {
    std::cout << loc << " ";
  }
  std::cout << "c: ";
  std::cout << "cnt: " << state.counter;
  std::cout << " cap: " << state.capacity;
  std::cout << " time: " << state.timeWithMultiplier;
  std::cout << " curr: " << state.lastVisited;
  std::cout << " pot: " << potential;
  std::cout << " sp: " << shortestPathDistance;
  std::cout << " }";
};

void VRPTWNode::print() const
{
  std::cout << "counter: " << state.counter << " ";
  std::cout << "capacity: " << state.capacity << " ";
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
  // default time and capacity step sizes
  if (params.bucketsPerVertex <= 0)
  {
    capacityStepSize = 1 * (*std::min_element(vrptw.demands.begin()+1, vrptw.demands.end()));
    if (vrptw.vrptwCapacityType == NO_RELAX_CAPACITY)
    {
      capacityStepSize = 1;
    }

    if (vrptw.capacityDiscretization != 0)
    {
      capacityStepSize = vrptw.capacityDiscretization;
    }
  }
  else
  {
    if (vrptw.vrptwCapacityType == NO_RELAX_CAPACITY)
    {
      capacityStepSize = 1;
    }
    else
    {
      double bucketStepSize1 = vrptw.capacity * 1.0 / std::sqrt(params.bucketsPerVertex);
      capacityStepSize = bucketStepSize1;
    }
  }

  setTimeStepSize();

  //TODO(akarahal) compute better bound on number of locations in longest route
  //if (vrptw.counterType == VRPTWCounterType::COUNTER)
  //{
    // compile DD to set max counter value as longest path
    // use default step sizes
    //  compileNg(8);
  //}
};

void VRPTWDecisionDiagram::setTimeStepSize()
{
  // default time and capacity step sizes
  if (params.bucketsPerVertex <= 0)
  {
    timeStepSize = vrptw.timeStateMultiplier * vrptw.timeStateDiscretization;
  }
  else
  {
    if (vrptw.vrptwCapacityType == NO_RELAX_CAPACITY)
    {
      double bucketStepSize = (vrptw.endTimes[0] - vrptw.startTimes[0]) * 1.0 / params.bucketsPerVertex;
      timeStepSize = vrptw.timeStateMultiplier * bucketStepSize;
    }
    else
    {
      double bucketStepSize = (vrptw.endTimes[0] - vrptw.startTimes[0]) * 1.0 / std::sqrt(params.bucketsPerVertex);
      timeStepSize = vrptw.timeStateMultiplier * bucketStepSize;

    }
  }
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
    //std::cout << nodes.size() << std::endl;
    //newNode.print();

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
  arcsUsed.push_back(false);

  nodes[fromNodeIndex].outArcs.push_back(newArcIndex);
  nodes[toNodeIndex].inArcs.push_back(newArcIndex);

  locationToArcs[toLocation].push_back(newArcIndex);

  if (!arcs[newArcIndex].isReverseArc)
  {
    // add to appropriate cut sets
    capCutSetArcs.resize(capCutSets.size());
    for (int capCutIndex=0; capCutIndex<capCutSets.size(); ++capCutIndex)
    {
      auto capCutSet = capCutSets[capCutIndex];
      int fromLoc = nodes[arcs[newArcIndex].fromNodeIndex].state.lastVisited;
      int toLoc = nodes[arcs[newArcIndex].toNodeIndex].state.lastVisited;
      bool fromLocInSet = (std::find(capCutSet.begin(), capCutSet.end(), fromLoc) != capCutSet.end());
      bool toLocInSet = (std::find(capCutSet.begin(), capCutSet.end(), toLoc) != capCutSet.end());
      if (fromLocInSet && toLocInSet)
      {
        capCutSetArcs[capCutIndex].push_back(newArcIndex);
      }
    }

    // use duals for combs
    combCutArcs.resize(combRHS.size());
    for (int combIndex=0; combIndex<combRHS.size(); ++combIndex)
    {
      std::set<int> combCutArcSet;
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

    // use duals for src
    cliqueCutArcs.resize(cliqueCuts.size());
    int cliqueCutIndex = 0;
    for (auto cliqueCut : cliqueCuts)
    {
      for (auto nodeIndex : cliqueCut)
      {
        const auto node = nodes[nodeIndex];
        if (std::find(node.outArcs.begin(), node.outArcs.end(), newArcIndex) != node.outArcs.end())
        {
          cliqueCutArcs[cliqueCutIndex].push_back(newArcIndex);
        }
      }

      ++cliqueCutIndex;
    }
  }

  return newArcIndex;
}

int VRPTWDecisionDiagram::addReverseArc(int forwardArcIndex)
{
  VRPTWArc forwardArc = arcs[forwardArcIndex];
  VRPTWArc reverseArc(forwardArc.toNodeIndex, forwardArc.fromNodeIndex, nodes[forwardArc.fromNodeIndex].state.lastVisited, -forwardArc.distance, -forwardArc.coeff, -forwardArc.cijPi);
  int reverseArcIndex = arcs.size();
  arcs.push_back(reverseArc);
  arcsUsed.push_back(false);

  nodes[forwardArc.toNodeIndex].outArcs.push_back(reverseArcIndex);
  nodes[forwardArc.fromNodeIndex].inArcs.push_back(reverseArcIndex);

  // store for next time
  arcReverseArc[forwardArcIndex] = reverseArcIndex;
  arcReverseArc[reverseArcIndex] = forwardArcIndex;

  return reverseArcIndex;
}

// TODO(akarahal) figure out PDP, need to update capacities, can decrease when drop off.
// state should be exact for true/false return to be true
bool VRPTWDecisionDiagram::generateNewStateFromExact(VRPTWNodeState& newState, int location)
{
  newState.counter = newState.counter + 1;
  if (vrptw.counterType == VRPTWCounterType::NO_USE_COUNTER)
  {
    newState.counter = newState.counter - 1;
  }
  if (newState.counter > vrptw.numLocations + 1)
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
  newState.capacity = newState.capacity + vrptw.demands[location];
  if (vrptw.oneOrMorePaths == OneOrMorePaths::ONE_PATH)
  {
    newState.capacity = 0;
  }
  if (newState.capacity > vrptw.capacity)
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

// only allow up to s = 2
void VRPTWDecisionDiagram::compileExactFukasawa(int s)
{
  // reserve but do not resize
  nodes.reserve(vrptw.capacity * (std::pow(vrptw.numLocations,s)));
  arcs.reserve(std::pow(vrptw.numLocations,2));
  arcsUsed.reserve(std::pow(vrptw.numLocations,2));

  // root node r and terminal node t
  std::set<int> initialDeque = {};
  VRPTWNodeState rootNodeState(0,0,0,0,initialDeque);
  addNode(rootNodeState);
  rootNodeIndex = 0;
  VRPTWNodeState terminalNodeState(vrptw.numLocations*1000,vrptw.capacity,vrptw.endTimes[0],0,initialDeque);
  addNode(terminalNodeState);
  terminalNodeIndex = 1;

  // create nodes except r/t
  int nodeIndex = 0;
  while (nodeIndex < nodes.size())
  {
    for (int location=1; location<vrptw.numLocations; ++location)
    {
      // check capacity
      int newCapacity = nodes[nodeIndex].state.capacity + vrptw.demands[location];
      if (newCapacity > vrptw.capacity)
      {
        continue;
      }

      // check visited set
      if (std::find(nodes[nodeIndex].state.visited.begin(), nodes[nodeIndex].state.visited.end(), location) != nodes[nodeIndex].state.visited.end())
      {
        continue;
      }

      // check time to begin isnt past endTime
      int lastVisitedLocation = nodes[nodeIndex].state.lastVisited;
      double earliestStartTime = (nodes[nodeIndex].state.timeWithMultiplier / (vrptw.timeStateMultiplier * 1.0)) + vrptw.distances[lastVisitedLocation][location] + vrptw.serviceTimes[lastVisitedLocation];
      if (earliestStartTime > vrptw.endTimes[location])
      {
        continue;
      }

      std::set<int> newVisited = nodes[nodeIndex].state.visited;
      newVisited.insert(location);
      if (s == 1)
      {
        newVisited.erase(nodes[nodeIndex].state.lastVisited);
      }
      else if (s == 2)
      {
        for (int loc : nodes[nodeIndex].state.visited)
        {
          if (loc != nodes[nodeIndex].state.lastVisited)
          {
            newVisited.erase(loc);
          }
        }
      }
      else
      {
        std::cout << "ERROR - dont use s > 2" << std::endl;
        return;
      }

      int newTimeWithMultiplier = (int)(earliestStartTime * vrptw.timeStateMultiplier);
      newTimeWithMultiplier = std::max(newTimeWithMultiplier, vrptw.startTimes[location] * vrptw.timeStateMultiplier);

      // bucketing
      // use more dynamic bucketing if otherwise it'll make a bad loop
      int timeQuotient = newTimeWithMultiplier / timeStepSize;
      int newTimeWithMultiplierDiscretized = timeStepSize * timeQuotient;
      if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::NO_TIME_WINDOWS)
      {
        newTimeWithMultiplierDiscretized = 0;
      }
      else if (newTimeWithMultiplierDiscretized <= newTimeWithMultiplier)
      {
        newTimeWithMultiplierDiscretized = newTimeWithMultiplier;
      }
 
      int loadQuotient = newCapacity / capacityStepSize;
      int newCapacityDiscretized = capacityStepSize * loadQuotient;
      if (vrptw.vrptwCapacityType == VRPTWCapacityType::RELAX_CAPACITY)
      {
        newCapacityDiscretized = 0;
      }

      int newCounter = nodes[nodeIndex].state.counter + 1;
      if (vrptw.counterType == VRPTWCounterType::NO_USE_COUNTER)
      {
        newCounter = nodes[nodeIndex].state.counter;
      }
      VRPTWNodeState newState(newCounter, newCapacityDiscretized, newTimeWithMultiplierDiscretized, location, newVisited);
      int newNodeIndex = addNode(newState);
      addArc(nodeIndex, newNodeIndex);
    }

    nodeIndex = nodeIndex + 1;
  }

  // add arcs to terminal node
  nodeIndex = 2;
  while (nodeIndex < nodes.size())
  {
    // need to be able to make it back in time
    int lastVisitedLocation = nodes[nodeIndex].state.lastVisited;
    double earliestStartTime = (nodes[nodeIndex].state.timeWithMultiplier / (vrptw.timeStateMultiplier * 1.0)) + vrptw.distances[lastVisitedLocation][0] + vrptw.serviceTimes[lastVisitedLocation];
    if (earliestStartTime > vrptw.endTimes[0])
    {
      continue;
    }
    addArc(nodeIndex, terminalNodeIndex);
    nodeIndex = nodeIndex + 1;
  }

  DBG(print();)
};

bool sort_by_second(const std::pair<int,int>& a, const std::pair<int,int>& b)
{
  return (a.second < b.second);
}

void VRPTWDecisionDiagram::compileNgRoute(int s)
{
  std::cout << "time step mult: " << vrptw.timeStateMultiplier << std::endl;
  std::cout << "time step: " << timeStepSize << std::endl;
  std::cout << "cap step: " << capacityStepSize << std::endl;
  // reserve but do not resize
  if (vrptw.oneOrMorePaths != ONE_PATH)
  {
    nodes.reserve(vrptw.capacity * (std::pow(vrptw.numLocations,2)));
  }
  else
  {
    nodes.reserve(std::pow(vrptw.numLocations,2));
  }
  arcs.reserve(std::pow(vrptw.numLocations,2)/100);

  // get ng sets as closest neighbors
  std::vector<std::set<int>> ngSets;
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

  int timeWindowBinary = 1;
  if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::NO_TIME_WINDOWS)
  {
    timeWindowBinary = 0;
  }

  int capacityBinary = 1;
  if (vrptw.vrptwCapacityType == VRPTWCapacityType::RELAX_CAPACITY)
  {
    capacityBinary = 0;
  }

  int counterBinary = 1;
  if (vrptw.counterType == VRPTWCounterType::NO_USE_COUNTER)
  {
    counterBinary = 0;
  }

  // root node r and terminal node t
  std::set<int> initialDeque = {};
  VRPTWNodeState rootNodeState(1,0,0,0,initialDeque);
  addNode(rootNodeState);
  rootNodeIndex = 0;
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
  VRPTWNodeState terminalNodeState(vrptw.numLocations,vrptw.capacity,terminalEndTime,terminalNodeLoc,initialDeque);
  addNode(terminalNodeState);
  terminalNodeIndex = 1;

  // create nodes except r/t
  int nodeIndex = 0;
  while (nodeIndex < nodes.size())
  {
    const VRPTWNode& nodeToCheck = nodes[nodeIndex];
    if (nodeToCheck.state.counter >= vrptw.numLocations)
    {
      nodeIndex = nodeIndex + 1;
      continue;
    }
    for (int location=1; location<largestLocationIndex; ++location)
    {
      // adding node in this loop could change its location, and we want to go by reference
      const VRPTWNode& node = nodes[nodeIndex];
      int newCapacity = node.state.capacity + vrptw.demands[location];
      if (newCapacity > vrptw.capacity)
      {
        continue;
      }

      if (std::find(node.state.visited.begin(), node.state.visited.end(), location) != node.state.visited.end())
      {
        continue;
      }

      // precedence only considered for root node, otherwise need to learn by separation
      bool precedenceIssue = false;
      if ((nodeIndex == 0) && !vrptw.precedences[location].empty())
      {
        if (vrptw.precedences[location].size() > 0)
        {
          precedenceIssue = true;
        }
      }

      // if next location is a precedence... can't go there
      if (vrptw.precedences[node.state.lastVisited].find(location) != vrptw.precedences[node.state.lastVisited].end())
      {
        precedenceIssue = true;
      }

      if (precedenceIssue)
      {
        continue;
      }

      int lastVisitedLocation = node.state.lastVisited;
      double earliestStartTime = 0.0;
      if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
      {
        lastVisitedLocation = node.state.lastVisited;
        earliestStartTime = (node.state.timeWithMultiplier / (vrptw.timeStateMultiplier * 1.0)) + vrptw.distances[lastVisitedLocation][location] + vrptw.serviceTimes[lastVisitedLocation];
        if (earliestStartTime > vrptw.endTimes[location])
        {
          continue;
        }
      }

      std::set<int> visitedSet;
      for (int loc : node.state.visited)
      {
        visitedSet.insert(loc);
      }
      std::set<int> ngSet = ngSets[location];
      std::set<int> newVisited;
      std::set_intersection(visitedSet.begin(), visitedSet.end(), ngSet.begin(), ngSet.end(), std::inserter(newVisited, newVisited.begin()));
      newVisited.insert(location);

      int newCounter = node.state.counter + counterBinary;

      int newTimeWithMultiplier = (int)(earliestStartTime * vrptw.timeStateMultiplier);
      newTimeWithMultiplier = std::max(newTimeWithMultiplier, vrptw.startTimes[location] * vrptw.timeStateMultiplier);

      // bucketing
      // use more dynamic bucketing if otherwise it'll make a bad loop
      int timeQuotient = newTimeWithMultiplier / timeStepSize;
      int newTimeWithMultiplierDiscretized = timeStepSize * timeQuotient * timeWindowBinary;
      /*
      if ((vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS) && (newTimeWithMultiplierDiscretized <= node.state.timeWithMultiplier))
      {
        newTimeWithMultiplierDiscretized = newTimeWithMultiplier;
      }
      */
 
      int loadQuotient = newCapacity / capacityStepSize;
      int newCapacityDiscretized = capacityStepSize * loadQuotient * capacityBinary;

      // time window removals - need to be able to make it back in time
      if (vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS)
      {
        if (earliestStartTime + vrptw.distances[location][terminalNodeLoc] + vrptw.serviceTimes[location] > vrptw.endTimes[terminalNodeLoc])
        {
          continue;
        }
      }

      VRPTWNodeState newState(newCounter, newCapacityDiscretized, newTimeWithMultiplierDiscretized, location, newVisited);
      int newNodeIndex = addNode(newState);
      addArc(nodeIndex, newNodeIndex);
    }

    nodeIndex = nodeIndex + 1;
  }

  DBG(print();)
  // add arcs to terminal node
  nodeIndex = 2;
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

  //print();
  //checkLC121SolutionPossible();
  //checkLRC121SolutionPossible();
};

void VRPTWDecisionDiagram::setCoeffsAsDistances()
{
  for (VRPTWArc& arc : arcs)
  {
    arc.coeff = arc.distance;
    arc.cijPi = arc.distance;
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

void VRPTWDecisionDiagram::setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(const std::vector<double>& lambda, const std::vector<double>& capDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, LPSolveType solveType)
{
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    VRPTWArc& arc = arcs[arcIndex];
    arc.coeff = arc.distance - lambda[arc.location];
    arc.cijPi = arc.distance - lambda[arc.location];
  }

  // use duals for size as we may have added more cuts, but are using previous solution
  for (int capCutIndex=0; capCutIndex<capDuals.size(); ++capCutIndex)
  {
    for (int arcIndex : capCutSetArcs[capCutIndex])
    {
      VRPTWArc& arc = arcs[arcIndex];
      if (solveType == LPSolveType::LPSolver)
      {
        arc.coeff = arc.coeff - capDuals[capCutIndex];
        arc.cijPi = arc.cijPi - capDuals[capCutIndex];
      }
      else
      {
        arc.coeff = arc.coeff + capDuals[capCutIndex];
        arc.cijPi = arc.cijPi + capDuals[capCutIndex];
      }
    }
  }

  // use duals for combs
  for (int combIndex=0; combIndex<combDuals.size(); ++combIndex)
  {
    for (int arcIndex : combCutArcs[combIndex])
    {
      VRPTWArc& arc = arcs[arcIndex];
      arc.coeff = arc.coeff - combDuals[combIndex];
      arc.cijPi = arc.cijPi - combDuals[combIndex];
    }
  }

  // use duals for src
  for (int cliqueCutIndex=0; cliqueCutIndex<srcDuals.size(); ++cliqueCutIndex)
  {
    if (isCliqueCutActive(cliqueCutIndex))
    {
      for (int arcIndex : cliqueCutArcs[cliqueCutIndex])
      {
        VRPTWArc& arc = arcs[arcIndex];
        if (solveType == LPSolveType::LPSolver)
        {
          arc.coeff = arc.coeff - srcDuals[cliqueCutIndex];
          arc.cijPi = arc.cijPi - srcDuals[cliqueCutIndex];
        }
        else
        {
          arc.coeff = arc.coeff + srcDuals[cliqueCutIndex];
          arc.cijPi = arc.cijPi + srcDuals[cliqueCutIndex];
        }
      }
    }
    else
    {
      srcDuals[cliqueCutIndex] = 0.0;
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

void VRPTWDecisionDiagram::getCutSetValues(std::vector<double>& cutValues)
{
  cutValues.clear();
  for (int capCutIndex=0; capCutIndex<capCutSets.size(); ++capCutIndex)
  {
    double cutSetValue = 0.0;
    for (int arcIndex : capCutSetArcs[capCutIndex])
    {
      cutSetValue += arcs[arcIndex].heuristicFlow;
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

void VRPTWDecisionDiagram::getCliqueCutValues(std::vector<double>& cliqueCutValues)
{
  cliqueCutValues.clear();
  int cliqueCutIndex = 0;
  for (auto cliqueCut : cliqueCuts)
  {
    if (isCliqueCutActive(cliqueCutIndex))
    {
      // get all disjoint path nodes possible greedily
      double value = 0.0;
      for (auto arcIndex : cliqueCutArcs[cliqueCutIndex])
      {
        value = value + arcs[arcIndex].heuristicFlow;
      }

      cliqueCutValues.push_back(value);
      ++cliqueCutIndex;
    }
  }
}

void VRPTWDecisionDiagram::addConnectedNodesToBlacklist(int nodeIndex, int demandLimit, std::set<int>& blacklist, const std::set<int>& nodesUsed)
{
  int load = nodes[nodeIndex].state.capacity;
  std::vector<bool> seenNode(nodes.size(), false);
  for (auto capNodeIndices : nodeOrdering)
  {
    // mark seen nodes and try to find nodeIndex2
    for (int nodeIndex : capNodeIndices.second)
    {
      // skip to load range
      if (nodes[nodeIndex].state.capacity < load)
      {
        break;
      }
      else if (nodes[nodeIndex].state.capacity > demandLimit)
      {
        return;
      }

      if (!seenNode[nodeIndex])
      {
        continue;
      }

      for (int arcIndex : nodes[nodeIndex].outArcs)
      {
        int toNodeIndex = arcs[arcIndex].toNodeIndex;
        seenNode[toNodeIndex] = true;
        if (nodesUsed.find(toNodeIndex) != nodesUsed.end())
        {
          blacklist.insert(toNodeIndex);
        }
      }
    }
  }
}

bool VRPTWDecisionDiagram::areNodesConnected(int nodeIndex1, int nodeIndex2)
{
  int load1 = nodes[nodeIndex1].state.capacity;
  int load2 = nodes[nodeIndex2].state.capacity;
  std::vector<bool> seenNode(nodes.size(), false);
  if (load1 < load2)
  {
    seenNode[nodeIndex1] = true;
    for (auto capNodeIndices : nodeOrdering)
    {
      // mark seen nodes and try to find nodeIndex2
      for (int nodeIndex : capNodeIndices.second)
      {
        // skip to load range
        if (nodes[nodeIndex].state.capacity < load1)
        {
          break;
        }
        else if (nodes[nodeIndex].state.capacity > load2)
        {
          return false;
        }

        if (!seenNode[nodeIndex])
        {
          continue;
        }
        else if (nodeIndex == nodeIndex2)
        {
          return true;
        }

        for (int arcIndex : nodes[nodeIndex].outArcs)
        {
          int toNodeIndex = arcs[arcIndex].toNodeIndex;
          seenNode[toNodeIndex] = true;
        }
      }
    }

    return false;
  }
  else if (load1 > load2)
  {
    seenNode[nodeIndex2] = true;
    for (auto capNodeIndices : nodeOrdering)
    {
      // mark seen nodes and try to find nodeIndex2
      for (int nodeIndex : capNodeIndices.second)
      {
        // skip to load range
        if (nodes[nodeIndex].state.capacity < load2)
        {
          break;
        }
        else if (nodes[nodeIndex].state.capacity > load1)
        {
          return false;
        }

        if (!seenNode[nodeIndex])
        {
          continue;
        }
        else if (nodeIndex == nodeIndex1)
        {
          return true;
        }

        for (int arcIndex : nodes[nodeIndex].outArcs)
        {
          int toNodeIndex = arcs[arcIndex].toNodeIndex;
          seenNode[toNodeIndex] = true;
        }
      }
    }

    return false;
  }
  else
  {
    return false;
  }
}

/*
void VRPTWDecisionDiagram::findSRCThree(const std::vector<std::set<int>>& decomposition, const std::vector<double>& stepSizes, bool& cutAdded)
{
  // find nodes in solution and sort by load
  std::unordered_map<int,int> nodeFlowWeight;
  std::set<int> nodesUsed;
  if (decomposition.empty())
  {
    for (auto arc : arcs)
    {
      if (arc.heuristicFlow > 0.01)
      {
        int fromNodeIndex = arc.fromNodeIndex;
        int toNodeIndex = arc.toNodeIndex;
        if (fromNodeIndex != rootNodeIndex)
        {
          nodesUsed.insert(fromNodeIndex);
        }
        nodesUsed.insert(toNodeIndex);
        nodeFlowWeight[fromNodeIndex] = nodeFlowWeight[fromNodeIndex] + std::ceil(arc.heuristicFlow * 100.0);
      }
    }
  }
  else
  {
    double totalFlow = 0.0;
    for (double flow : stepSizes)
    {
      totalFlow += flow;
    }

    std::map<std::pair<int,int>,double> edgeFlows;
    for (int routeIndex=0; routeIndex<decomposition.size(); ++routeIndex)
    {
      auto route = decomposition[routeIndex];
      double flow = (stepSizes[routeIndex] / totalFlow);
      if (flow > 0.001)
      {
        for (int arcIndex : route)
        {
          int fromNodeIndex = arcs[arcIndex].fromNodeIndex;
          int toNodeIndex = arcs[arcIndex].toNodeIndex;
          if (fromNodeIndex != rootNodeIndex)
          {
            nodesUsed.insert(fromNodeIndex);
          }
          nodesUsed.insert(toNodeIndex);
          nodeFlowWeight[fromNodeIndex] = nodeFlowWeight[fromNodeIndex] + std::ceil(flow * 100.0);
        }
      }
    }
  }

  std::cout << "finding max clique" << std::endl;
  std::vector<int> nodesUsedVector(nodesUsed.begin(), nodesUsed.end());
  graph_t* g;
  g=(graph_t*)calloc(1,sizeof(graph_t));
  g->n = (int)nodesUsedVector.size();
  g->edges = (set_t*)calloc(g->n,sizeof(set_t));
  g->weights = (int*)calloc(g->n,sizeof(int));
  for (int i=0; i<g->n; ++i)
  {
    g->edges[i] = set_new(g->n);
    if (nodeFlowWeight.find(nodesUsedVector[i]) != nodeFlowWeight.end())
    {
      g->weights[i] = nodeFlowWeight[nodesUsedVector[i]];
    }
    else
    {
      g->weights[i] = 0.0;
    }
  }
  for (int index1=0; index1<nodesUsedVector.size()-1; ++index1)
  {
    if (nodeFlowWeight[nodesUsedVector[index1]] >= 1)
    {
      for (int index2=index1+1; index2<nodesUsedVector.size(); ++index2)
      {
        if (nodeFlowWeight[nodesUsedVector[index2]] >= 1)
        {
          int nodeIndex1 = nodesUsedVector[index1];
          int nodeIndex2 = nodesUsedVector[index2];
          auto nodeState1 = std::set<int>(nodes[nodeIndex1].state.visited.begin(), nodes[nodeIndex1].state.visited.end());
          auto nodeState2 = std::set<int>(nodes[nodeIndex2].state.visited.begin(), nodes[nodeIndex2].state.visited.end());
          std::set<int> intersection;
          std::set_intersection(nodeState1.begin(), nodeState1.end(),
                                nodeState2.begin(), nodeState2.end(),
                                std::inserter(intersection, intersection.begin()));
          if (!intersection.empty())
          {
            if (!areNodesConnected(nodeIndex1, nodeIndex2))
            {
              GRAPH_ADD_EDGE(g,index1,index2);
            }
          }
        }
      }
    }
  }
  std::cout << "graph setup complete" << std::endl;
  clique_options* opts;
  opts = (clique_options*)malloc(sizeof(clique_options));
  opts->time_function=clique_print_time;
  opts->output=stderr;
  opts->reorder_function = NULL;
  opts->reorder_map = NULL;
  opts->user_function=record_clique_func;
  opts->user_data=NULL;
  opts->clique_list=NULL;
  opts->clique_list_length=0;
  //set_t clique = clique_find_single(g,103,200,FALSE,opts);
  set_t clique = clique_find_single(g,0,0,FALSE,opts);
  if (!clique)
  {
    std::cout << "no clique returned" << std::endl;
    return;
  }
  int cliqueWeight = graph_subgraph_weight(g,clique);
  std::cout << "got max clique with weight: " << cliqueWeight << std::endl;
  std::set<int> cliqueSet;
  for (int i=0; i<SET_MAX_SIZE(clique); ++i)
  {
    if (SET_CONTAINS(clique,i))
    {
      cliqueSet.insert(nodesUsedVector[i]);
      std::cout << nodesUsedVector[i] << " has state: ";
      for (auto n1 : nodes[nodesUsedVector[i]].state.visited)
      {
        std::cout << n1 << ",";
      }
      std::cout << std::endl;
    }
  }

  if ((cliqueWeight > (100 + cliqueSet.size())) && (cliqueSet.size() > 1))
  {
    // NOTE(akarahal) should try to maximize nodes in the clique set, even if not in flow
    int beforeSize = cliqueCuts.size();
    cliqueCuts.insert(cliqueSet);
    int afterSize = cliqueCuts.size();
    if (beforeSize < afterSize)
    {
      cutAdded = true;
      cliqueCutActive.push_back(true);

      // add relevant arcs
      cliqueCutArcs.resize(cliqueCuts.size());
      for (auto nodeIndex : cliqueSet)
      {
        const auto node = nodes[nodeIndex];
        for (int arcIndex : node.outArcs)
        {
          cliqueCutArcs.back().push_back(arcIndex);
        }
      }
    }
  }
}
*/

void VRPTWDecisionDiagram::addColumnForLPCG(const std::vector<int>& route)
{
  int currentNodeIndex = rootNodeIndex;
  for (int routeIndex=0; routeIndex<route.size(); ++routeIndex)
  {
    for (int arcIndex : nodes[currentNodeIndex].outArcs)
    {
      if (arcs[arcIndex].location == route[routeIndex])
      {
        currentNodeIndex = arcs[arcIndex].toNodeIndex;
        arcsUsed[arcIndex] = true;
        break;
      }
    }

    if (currentNodeIndex == terminalNodeIndex)
    {
      return;
    }
  }
};

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
    arcsUsed[priorArcIndexShortestPath[terminalNodeIndex]] = true;
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

double VRPTWDecisionDiagram::computeShortestPathBFSWang(std::vector<int>& treeByParentArcs, std::vector<int>& routeByArc)
{
  // setup data structures to store shortestpaths and the path itself
  for (VRPTWNode& node : nodes)
  {
    node.shortestPathDistance = INF;
  }
  nodes[rootNodeIndex].shortestPathDistance = 0.0;

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
        }
      }
    }
  }

  // construct route
  routeByArc.push_back(treeByParentArcs[terminalNodeIndex]);
  arcsUsed[treeByParentArcs[terminalNodeIndex]] = true;
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

void VRPTWDecisionDiagram::initializeColumnsByLPDecomp()
{
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    if (arcs[arcIndex].decompositionFlow > 0.0000001)
    {
      arcsUsed[arcIndex] = true;
    }
  }
}

// can we merge nodes to reduce problem size?
// can we warm start with a solution? dual solution?
// can we reduce memory?
double VRPTWDecisionDiagram::setupAndSolveFlowModel(FlowType flowType, IncludeCoverConstraints includeCoverConstraints, UseColumnGeneration useCg, std::vector<double>& duals, double& singlePathDual, std::vector<double>& capDuals, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  // setup model
  IloEnv env;
  IloModel flowModel(env);
  IloRangeArray coverConstraints(env);
  IloRangeArray flowConservationConstraints(env);
  IloRangeArray singlePathConstraint(env);
  IloRangeArray capConstraints(env);
  IloRangeArray combConstraints(env);
  IloRangeArray srcConstraints(env);
  IloNumVarArray x(env, arcs.size());
  IloExpr objective(env);

  // setup variables and objective
  int numNonZeroVars = 0;
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    if (!arcs[arcIndex].isReverseArc && isArcAlive(arcIndex))
    {
      if (flowType == FlowType::IP)
      {
        x[arcIndex] = IloNumVar(env, 0, 1, ILOINT);
        objective += x[arcIndex] * arcs[arcIndex].coeff;
        //if ((arcs[arcIndex].fromNodeIndex == rootNodeIndex) && (vrptw.problemType == ProblemType::PDP))
        //{
        //  objective += x[arcIndex] * 1;
        //}
        ++numNonZeroVars;
      }
      else
      {
        if ((useCg == UseColumnGeneration::USE_CG) && !arcsUsed[arcIndex])
        {
          x[arcIndex] = IloNumVar(env, 0, 0);
        }
        else
        {
          ++numNonZeroVars;
          x[arcIndex] = IloNumVar(env, 0, 1);
          objective += x[arcIndex] * arcs[arcIndex].coeff;
          //if ((arcs[arcIndex].fromNodeIndex == rootNodeIndex) && (vrptw.problemType == ProblemType::PDP))
          //{
          //  objective += x[arcIndex] * 1;
          //}
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

  // setup coverConstraints
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
  }
  flowModel.add(coverConstraints);

  // setup flow conservation constraints (skip root and terminal)
  // deal with fixed nodes
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

  // one path for TSPs, should this even be solved with an LP? lol
  if (vrptw.oneOrMorePaths == ONE_PATH)
  {
    IloExpr onePath(env);
    for (auto arcIndex : nodes[rootNodeIndex].outArcs)
    {
      onePath += x[arcIndex];
    }
    singlePathConstraint.add(onePath == 1);
    flowModel.add(singlePathConstraint);
  }

  if (flowType == FlowType::LP)
  {
    // RCC - rounded capacity cuts, might even need to use for IPs
    for (int capCutIndex=0; capCutIndex<capCutSets.size(); ++capCutIndex)
    {
      IloExpr cutSetSum(env);
      for (int arcIndex : capCutSetArcs[capCutIndex])
      {
        cutSetSum += x[arcIndex];
      }

      capConstraints.add(cutSetSum <= capCutSetsRHS[capCutIndex]);
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
    int cliqueCutIndex = 0;
    for (auto cliqueCut : cliqueCuts)
    {
      if (isCliqueCutActive(cliqueCutIndex))
      {
        IloExpr cliqueCutSum(env);
        for (auto arcIndex : cliqueCutArcs[cliqueCutIndex])
        {
          cliqueCutSum += x[arcIndex];
        }

        srcConstraints.add(cliqueCutSum <= 1);
        ++cliqueCutIndex;
      }
      else
      {
        srcConstraints.add(x[0] <= 1);
      }
    }
    flowModel.add(srcConstraints);
  }

  // solve model
  IloCplex solver(flowModel);
  solver.setOut(env.getNullStream());
  solver.setWarning(env.getNullStream());
  solver.setError(env.getNullStream());
  //solver.setParam(IloCplex::Param::TimeLimit, timelimit);
  //solver.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Barrier);
  solver.setParam(IloCplex::Param::Threads, 1);
  //solver.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
  solver.exportModel("LPflowmodel.lp");

  // x is var[i], dj[i] reduced cost, pi[i] specifes starting dual for rng[i]
  //solver.setVectors(x, dj, var, slack, pi, rng)
  IloNumArray startDuals(env);
  for (auto dual : duals)
  {
    startDuals.add(dual);
  }
  solver.setStart(NULL, NULL, x, NULL, startDuals, coverConstraints);
  solver.solve();

  // get results
  IloAlgorithm::Status solverStatus = solver.getStatus();
  if (solverStatus == IloAlgorithm::Optimal)
  {
    // store results
    for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
    {
      if (!arcs[arcIndex].isReverseArc)
      {
        arcs[arcIndex].decompositionFlow = solver.getValue(x[arcIndex]);
        arcs[arcIndex].heuristicFlow= solver.getValue(x[arcIndex]);
      }
      DBG(
        if (arcs[arcIndex].decompositionFlow > 0)
        {
          std::cout << "primal [" << arcIndex << "]: " << arcs[arcIndex].decompositionFlow << std::endl;
        }
      )
    }

    double objValue = solver.getObjValue();
    /*
    if (vrptw.problemType == ProblemType::PDP)
    {
      for (auto arcIndex : nodes[rootNodeIndex].outArcs)
      {
        objValue = objValue - (arcs[arcIndex].decompositionFlow * 1);
      }
    }
    */

    if (flowType == FlowType::LP)
    {
      if (includeCoverConstraints == IncludeCoverConstraints::Y)
      {
        IloNumArray lpDuals(env);
        solver.getDuals(lpDuals, coverConstraints);
        for (int dualIndex=0; dualIndex<vrptw.numLocations; ++dualIndex)
        {
          duals[dualIndex] = lpDuals[dualIndex];
        }
        DBG(std::cout << "DUALS: ";
        for (int dualIndex=0; dualIndex<vrptw.numLocations; ++dualIndex)
        {
          //std::cout << "cover dual [" << dualIndex << "]: " << lpDuals[dualIndex] << std::endl;
          std::cout << lpDuals[dualIndex] << ",";
        }
        std::cout << std::endl;)
      }

      if (vrptw.oneOrMorePaths == ONE_PATH)
      {
        IloNumArray singlePathDualFromLP(env);
        solver.getDuals(singlePathDualFromLP, singlePathConstraint);
        singlePathDual = singlePathDualFromLP[0];
        DBG(std::cout << "single path dual: " << singlePathDual << std::endl;)
      }

      IloNumArray lpCapDuals(env);
      solver.getDuals(lpCapDuals, capConstraints);
      capDuals.resize(capCutSets.size());
      for (int dualIndex=0; dualIndex<capCutSets.size(); ++dualIndex)
      {
        capDuals[dualIndex] = lpCapDuals[dualIndex];
      }
      DBG(
      for (int dualIndex=0; dualIndex<capCutSets.size(); ++dualIndex)
      {
        std::cout << "cap dual [" << dualIndex << "]: " << capDuals[dualIndex] << std::endl;
      })
 
      IloNumArray lpCombDuals(env);
      solver.getDuals(lpCombDuals, combConstraints);
      combDuals.resize(teeths.size());
      for (int dualIndex=0; dualIndex<teeths.size(); ++dualIndex)
      {
        combDuals[dualIndex] = lpCombDuals[dualIndex];
      }
      DBG(
      for (int dualIndex=0; dualIndex<teeths.size(); ++dualIndex)
      {
        std::cout << "comb dual [" << dualIndex << "]: " << combDuals[dualIndex] << std::endl;
      })

      IloNumArray lpSrcDuals(env);
      solver.getDuals(lpSrcDuals, srcConstraints);
      srcDuals.resize(cliqueCuts.size());
      for (int dualIndex=0; dualIndex<cliqueCuts.size(); ++dualIndex)
      {
        srcDuals[dualIndex] = lpSrcDuals[dualIndex];
        DBG(std::cout << "src dual [" << dualIndex << "]: " << srcDuals[dualIndex] << std::endl;)
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
    return 0;
  }

  env.end();
  return 0;
};

double VRPTWDecisionDiagram::fixArcs(const std::vector<double>& lambda, double singlePathDual, const std::vector<double>& capDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, double lowerBound, LPSolveType solveType)
{
  DBG(std::cout << "fixing" << std::endl;)
  setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(lambda, capDuals, combDuals, srcDuals, solveType);

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

  // fix arcs based on lb + rc > ub
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    const VRPTWArc& arcToCheck = arcs[arcIndex];
    double bestPossibleReducedCost = shortestPathDown[arcToCheck.fromNodeIndex] + shortestPathUp[arcToCheck.toNodeIndex] + arcToCheck.coeff - singlePathDual;
    if ((lowerBound + bestPossibleReducedCost) > (vrptw.instanceUpperBound + 0.00001))
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
  DBG(std::cout << "number of fixed arcs: " << fixedArcs.size() << std::endl;
  std::cout << "percent fixed arcs: " << getPercentFixedArcs() << std::endl;)
  return getPercentFixedArcs();
}

// greedily take paths with flow
void VRPTWDecisionDiagram::primalHeuristic(std::vector<std::vector<int>>& routesByLocationPrimalHeuristic)
{
  DBG(std::cout << "running primal heuristic" << std::endl;)
  std::set<int> locationsCovered;
  locationsCovered.insert(0);
  bool rootFlowExists = true;
  while (locationsCovered.size() != vrptw.numLocations)
  {
    std::vector<int> route;
    route.push_back(0);
    int currentNodeIndex = rootNodeIndex;
    bool continueRoute = true;
    while (continueRoute && rootFlowExists)
    {
      int arcIndex = -1;
      for (int outArcIndex : nodes[currentNodeIndex].outArcs)
      {
        if (arcs[outArcIndex].heuristicFlow > 0.00001)
        {
          arcIndex = outArcIndex;
          break;
        }
      }

      if (arcIndex == -1)
      {
        continueRoute = false;
        if (currentNodeIndex == rootNodeIndex)
        {
          rootFlowExists = false;
        }
        break;
      }

      arcs[arcIndex].heuristicFlow = 0;
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
    }

    // add route if at least one location was added
    route.push_back(0);
    if (route.size() > 2)
    {
      routesByLocationPrimalHeuristic.push_back(route);
    }
 
    // add the rest of the locations
    if (!rootFlowExists)
    {
      std::vector<int> capacities;
      for (auto currRoute : routesByLocationPrimalHeuristic)
      {
        int cap = 0;
        for (int loc : currRoute)
        {
          cap = cap + vrptw.demands[loc];
        }
        capacities.push_back(cap);
      }

      for (int location=1; location<vrptw.numLocations; ++location)
      {
        if (locationsCovered.find(location) == locationsCovered.end())
        {
          for (int routeIndex=0; routeIndex<capacities.size(); ++routeIndex)
          {
            if (capacities[routeIndex] + vrptw.demands[location] <= vrptw.capacity)
            {
              capacities[routeIndex] = capacities[routeIndex] + vrptw.demands[location];
              auto beforeLast = routesByLocationPrimalHeuristic[routeIndex].end()-1;
              routesByLocationPrimalHeuristic[routeIndex].insert(beforeLast, location);
              locationsCovered.insert(location);
              break;
            }
          }
        }
      }

      std::vector<int> extraRoute;
      extraRoute.push_back(0);
      int extraRouteCap = 0;
      for (int location=1; location<vrptw.numLocations; ++location)
      {
        if (locationsCovered.find(location) == locationsCovered.end())
        {
          if ((extraRouteCap + vrptw.demands[location]) <= vrptw.capacity)
          {
            extraRouteCap = extraRouteCap + vrptw.demands[location];
            extraRoute.push_back(location);
            locationsCovered.insert(location);
          }
        }
      }
      extraRoute.push_back(0);
      if (extraRoute.size() > 2)
      {
        routesByLocationPrimalHeuristic.push_back(extraRoute);
      }
    }
  }

  DBG(
    for (auto primalRoute: routesByLocationPrimalHeuristic)
    {
      std::cout << "primalRoute: ";
      for (int loc : primalRoute)
      {
        std::cout << loc << " ";
      }
      std::cout << std::endl;
    }
  )
};

bool VRPTWDecisionDiagram::doesRouteExistByArcs(const std::vector<int>& routeArcs) const
{
  for (int index=0; index<(routeArcs.size()-1); ++index)
  {
    int currNodeIndex = arcs[routeArcs[index]].toNodeIndex;
    int nextNodeIndex = arcs[routeArcs[index+1]].fromNodeIndex;
    if (currNodeIndex != nextNodeIndex)
    {
      return false;
    }
  }

  return true;
};

bool VRPTWDecisionDiagram::doesRouteExistByLocations(const std::vector<int>& routeLocations) const
{
  int routeIndex = 0;
  int currNodeIndex = rootNodeIndex;
  while (currNodeIndex != terminalNodeIndex)
  {
    bool nextLocationPossible = false;
    for (int arcIndex : nodes[currNodeIndex].outArcs)
    {
      if (arcs[arcIndex].location == routeLocations[routeIndex])
      {
        currNodeIndex = arcs[arcIndex].toNodeIndex;
        nextLocationPossible = true;
        break;
      }
    }

    if (!nextLocationPossible)
    {
      std::cout << "not possible loc: " << routeLocations[routeIndex] << std::endl;
      print();
      return false;
    }
    routeIndex = routeIndex + 1;
  }

  return true;
};

void VRPTWDecisionDiagram::decomposeRoutes(std::vector<int>& routeArcs, std::vector<double>& flows, std::vector<std::vector<int>>& routeDecomposition, int maxS, DecompositionReason decompositionReason)
{
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
    double routeFlow = INF;
    double routeTime = 0;
    int routeLoad = 0;
    while (continueRoute)
    {
      bool arcWithFlowFound = false;
      for (int arcIndex : nodes[currentNodeIndex].outArcs)
      {
        if (arcs[arcIndex].decompositionFlow > 0.0000001)
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
          // get time of current route decomposing
          int currLocation = routeByLoc[routeByLoc.size()-1];
          int lastLocation = routeByLoc[routeByLoc.size()-2];
          routeTime = std::max(routeTime + vrptw.serviceTimes[lastLocation] + vrptw.distances[lastLocation][currLocation], vrptw.startTimes[currLocation]*1.0);
          routeLoad = routeLoad + vrptw.demands[currLocation];
          if (decompositionReason == DecompositionReason::SEPARATE)
          {
            for (int lastLocationIndex=0; lastLocationIndex<lastLocationsVisited.size(); ++lastLocationIndex)
            {
              // if we have already visited, return the route of arcs causing the infeasibility
              if (lastLocationsVisited[lastLocationIndex] == currLocation)
              {
                // can be smarter and only return subsection, but this complicates the separation routine
                //for (int index=lastLocationIndex+1; index<lastLocationsVisited.size(); ++index)
                for (int index=0; index<route.size(); ++index)
                {
                  //routeArcs.push_back(lastArcsVisited[index]);
                  routeArcs.push_back(route[index]);
                }
                //routeArcs.push_back(arcIndex);
                for (int routeArcIndex : route)
                {
                  arcs[routeArcIndex].decompositionFlow = arcs[routeArcIndex].decompositionFlow - routeFlow;
                }
                DBG(std::cout << "infeas already visited, separate " << currLocation << std::endl;)
                return;
              }
            }

            // if we are past the time window, return the route of arcs causing this
            if ((vrptw.vrptwTimeWindowType == VRPTWTimeWindowType::TIME_WINDOWS) && (vrptw.endTimes[currLocation] < routeTime))
            {
              DBG(std::cout << "infeas time window, separate " << vrptw.endTimes[currLocation] << " " << routeTime << std::endl;)
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

            // if we are over capacity, return the route of arcs causing this
            if ((vrptw.vrptwCapacityType == VRPTWCapacityType::RELAX_CAPACITY) && (routeLoad > vrptw.capacity))
            {
              DBG(std::cout << "infeas capacity, separate " << vrptw.capacity << " " << routeLoad << std::endl;)
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

            // if precedence not followed, return route of arcs causing this
            if (!vrptw.precedences.empty())
            {
              bool precedenceBroke = false;
              for (int index=0; index<route.size()-1; ++index)
              {
                if (vrptw.precedences[routeByLoc[index]].find(currLocation) != vrptw.precedences[routeByLoc[index]].end())
                {
                  precedenceBroke = true;
                  DBG(std::cout << "1 infeas precedence" << std::endl;)
                }
              }

              for (int precLoc : vrptw.precedences[currLocation])
              {
                if (std::find(routeByLoc.begin(), routeByLoc.end()-1, precLoc) == routeByLoc.end()-1)
                {
                  precedenceBroke = true;
                  DBG(std::cout << "2 infeas precedence " << precLoc << " " << currLocation << std::endl;)
                }
              }

              if (precedenceBroke)
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

            if (lastLocationsVisited.size() == maxS)
            {
              lastLocationsVisited.push_back(arcs[arcIndex].location);
              lastLocationsVisited.pop_front();
              lastNodesVisited.push_back(currentNodeIndex);
              lastNodesVisited.pop_front();
              lastArcsVisited.push_back(arcIndex);
              lastArcsVisited.pop_front();
            }
            else
            {
              lastLocationsVisited.push_back(arcs[arcIndex].location);
              lastNodesVisited.push_back(currentNodeIndex);
              lastArcsVisited.push_back(arcIndex);
            }
          }

          break;
        }
      }

      if (!arcWithFlowFound && (currentNodeIndex == rootNodeIndex))
      {
        return;
      }
    }

    flows.push_back(routeFlow);
    routeDecomposition.push_back(routeByLoc);

    DBG(std::cout << "decomposed route: ";)
    for (int arcIndex : route)
    {
      DBG(std::cout << arcs[arcIndex].location << " ";
      std::cout << arcIndex << " ";)
      arcs[arcIndex].decompositionFlow = arcs[arcIndex].decompositionFlow - routeFlow;
    }
    DBG(std::cout << std::endl;)
  }

  //print();
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

void VRPTWDecisionDiagram::separateInfeasibleRoute(const std::vector<int>& routeArcs, int maxS)
{
  // clear fixed arcs when dd changes
  for (int arcIndex : fixedArcs)
  {
    const VRPTWArc& arc = arcs[arcIndex];
    nodes[arc.fromNodeIndex].outArcs.push_back(arcIndex);
    nodes[arc.toNodeIndex].inArcs.push_back(arcIndex);
  }
  fixedArcs.clear();

  DBG(std::cout << "separating route arcs: ";
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
  std::cout << std::endl;)

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
    bool newNodeFeasible = generateNewStateFromExact(newState, nextLocation);

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
          bool copyArcFeasible = generateNewStateFromExact(copyArcState, oldArcLocation);
          if (copyArcFeasible)
          {
            int newArcIndex = addArc(nextNodeIndex, oldArcToNodeIndex);
            DBG(std::cout << "copy arc by adding arc index: " << newArcIndex << std::endl;)
            arcsUsed[newArcIndex] = true;
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
          bool existingArcFeasible = generateNewStateFromExact(copyArcState, existingNodeArcLocation);
          if (!existingArcFeasible)
          {
            DBG(std::cout << "remove arc " << existingNodeArcIndex << std::endl;)
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
          if (arcs[currNodeArcIndex].toNodeIndex == nextNodeIndex)
          {
            break;
          }
          else
          {
            DBG(std::cout << "remove arc " << currNodeArcIndex << std::endl;)
            removeArc(currNodeArcIndex);
            return;
          }
        }
      }
      return;
    }
  }
 
  DBG(std::cout << "done separating routes" << std::endl;)
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
  DBG(std::cout << "creating solution" << std::endl;)
  DBG(print();)
  // only go over arcs currently in the graph
  for (VRPTWNode node : nodes)
  {
    for (int arcIndex : node.outArcs)
    {
      if (arcs[arcIndex].isReverseArc)
      {
        int solutionArcIndex = reverseArc(arcIndex);
      }
    }
  }

  for (int arcIndex : clippedArcs)
  {
    VRPTWArc& arc = arcs[arcIndex];
    nodes[arc.fromNodeIndex].outArcs.push_back(arcIndex);
    nodes[arc.toNodeIndex].inArcs.push_back(arcIndex);
  }

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
      if (minReducedCost < 0)
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

double VRPTWDecisionDiagram::solveMinCostFlowModelWang(const std::vector<double>& duals, const std::vector<double>& capDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, std::vector<std::vector<int>>& shortestPathsByArc, bool& isDualFeasible, double& shortestPathLength)
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
  setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(duals, capDuals, combDuals, srcDuals, LPSolveType::LAGSolver);
  shortestPathLength = computeShortestPathBFSWang(treeByParentArcs, shortestPathByArc);
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
  int totalNumPaths = 0;
  int currentNumPaths = 1;
  while (true)
  {
    // identify nodes 4 update
    int branchNodeIndex = arcs[shortestPathByArc[0]].toNodeIndex;
    identifyNodesForUpdate(branchNodeIndex, treeByChildArcs, nodesToUpdate);

    // find multi path
    int newPathArcIndex = findMultiPathNode(nodesToUpdate, currentNumPaths);
    if ((newPathArcIndex != -1) && (nodesToUpdate.find(arcs[newPathArcIndex].fromNodeIndex) == nodesToUpdate.end()) && (vrptw.oneOrMorePaths != OneOrMorePaths::ONE_PATH))
    {
      DBG(std::cout << "find multi path" << std::endl;)
      findMultiPath(newPathArcIndex, treeByParentArcs, shortestPathByArc);
      if (evaluateRouteCost(shortestPathByArc) >= -0.000000001)
      {
        break;
      }
      else
      {
        allCurrentPaths.push_back(shortestPathByArc);
        shortestPathsByArc.push_back(shortestPathByArc);
      }
      ++currentNumPaths;
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

      // one path - stop trying for more
      if (vrptw.oneOrMorePaths == OneOrMorePaths::ONE_PATH)
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

      nodesToUpdate.clear();
      allCurrentPaths.clear();
      if (evaluateRouteCost(shortestPathByArc) >= -0.000000001)
      {
        break;
      }
      allCurrentPaths.push_back(shortestPathByArc);
      shortestPathsByArc.push_back(shortestPathByArc);
      totalNumPaths = totalNumPaths + currentNumPaths;
      currentNumPaths = 1;
    }
  }

  // create solution from reverse arcs and clean up residual graph backward arcs
  DBG(
  for (auto path : shortestPathsByArc)
  {
    std::cout << "path: ";
    for (int arcIndex : path)
    {
      std::cout << arcs[arcIndex].fromNodeIndex << ", ";
    }
    std::cout << std::endl;
  })
  DBG(std::cout << "finished running, now create solution" << std::endl;)
  return createSolutionFromReverseArcsAndResetWang(clippedArcs, shortestPathsByArc);
};

bool VRPTWDecisionDiagram::checkFeasibleDual(const std::vector<double>& lambda, const std::vector<double>& rccDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, LPSolveType solveType)
{
  // get shortest path
  setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(lambda, rccDuals, combDuals, srcDuals, solveType);
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
    if (!doesRouteExistByLocations(routesByLocation[routeIndex]))
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
    if (!doesRouteExistByLocations(routesByLocation[routeIndex]))
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
    if (!doesRouteExistByLocations(routesByLocation[routeIndex]))
    {
      std::cout << "ERROR LRC121: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  std::cout << "clear" << std::endl;
  return true;
}

bool VRPTWDecisionDiagram::checkSinglePathPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({14,8,11,21,5,20,19,22,24,25,10,9,12,13,28,2,1,16,15,23,17,4,3,26,6,7,27,18,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    if (!doesRouteExistByLocations(routesByLocation[routeIndex]))
    {
      std::cout << "ERROR RC_202.3.txt: failed at index: " << routeIndex << std::endl;
      return false;
    }
  }

  return true;
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

void VRPTWDecisionDiagram::addCapCutSet(const std::vector<int>& cutSet)
{
  capCutSets.push_back(cutSet);
  capCutSetArcs.resize(capCutSets.size());
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    if (!arcs[arcIndex].isReverseArc)
    {
      int fromLoc = nodes[arcs[arcIndex].fromNodeIndex].state.lastVisited;
      int toLoc = nodes[arcs[arcIndex].toNodeIndex].state.lastVisited;
      bool fromLocInSet = (std::find(cutSet.begin(), cutSet.end(), fromLoc) != cutSet.end());
      bool toLocInSet = (std::find(cutSet.begin(), cutSet.end(), toLoc) != cutSet.end());
      if (fromLocInSet && toLocInSet)
      {
        capCutSetArcs.back().push_back(arcIndex);
      }
    }
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
        if ((fromLocInSet && !toLocInSet) || (!fromLocInSet && toLocInSet))
        {
          combCutArcs.back().push_back(arcIndex);
        }
      }
    }
  }
};

bool VRPTWDecisionDiagram::isCliqueCutActive(int index)
{
  if (cliqueCutActive[index])
  {
    auto it = cliqueCuts.begin();
    std::advance(it, index);
    auto cliqueSet = *it;
    std::vector<int> cliqueSetVec(cliqueSet.begin(), cliqueSet.end());
    for (int index1=0; index1<cliqueSet.size()-1; ++index1)
    {
      for (int index2=index1+1; index2<cliqueSet.size(); ++index2)
      {
        int nodeIndex1 = cliqueSetVec[index1];
        int nodeIndex2 = cliqueSetVec[index2];
        if (!areNodesConnected(nodeIndex1, nodeIndex2))
        {
          cliqueCutActive[index] = false;
          return false;
        }
      }
    }
  }

  return true;
};

