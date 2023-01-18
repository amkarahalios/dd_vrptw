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
  std::cout << "cap: " << state.capacity;
  std::cout << " time: " << state.timeTimesTen;
  std::cout << " curr: " << state.lastVisited;
  std::cout << " pot: " << potential;
  std::cout << " sp: " << shortestPathDistance;
  std::cout << " }";
};

void VRPTWNode::print() const
{
  std::cout << "capacity: " << state.capacity << " ";
  std::cout << "time: " << state.timeTimesTen << " ";
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
    capacityNodeIndices[state.capacity].push_back(newNodeIndex);
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

void VRPTWDecisionDiagram::moveArc(int arcIndex, int newToNodeIndex)
{
  // update arc info
  VRPTWArc& arc = arcs[arcIndex];
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
  VRPTWNodeState rootNodeState(0,0,0,initialDeque);
  addNode(rootNodeState);
  rootNodeIndex = 0;
  VRPTWNodeState terminalNodeState(vrptw.capacity,vrptw.endTimes[0],0,initialDeque);
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
      double earliestStartTime = (nodes[nodeIndex].state.timeTimesTen / 10.0) + vrptw.distances[location][lastVisitedLocation] + vrptw.serviceTimes[lastVisitedLocation];
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

      int newTimeTimesTen = (int)(earliestStartTime * 10);
      newTimeTimesTen = std::max(newTimeTimesTen, vrptw.startTimes[location] * 10);

      VRPTWNodeState newState(newCapacity, newTimeTimesTen, location, newVisited);
      int newNodeIndex = addNode(newState);
      addArc(nodeIndex, newNodeIndex);
    }

    nodeIndex = nodeIndex + 1;
  }

  // add arcs to terminal node
  nodeIndex = 2;
  while (nodeIndex < nodes.size())
  {
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
  // reserve but do not resize
  nodes.reserve(vrptw.capacity * (std::pow(vrptw.numLocations,2)));
  arcs.reserve(std::pow(vrptw.numLocations,2));

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
    for (int index=0; index<s; ++index)
    {
      ngSet.insert(indexDistances[index].first);
    }

    ngSets.push_back(ngSet);
  }

  // root node r and terminal node t
  std::set<int> initialDeque = {};
  VRPTWNodeState rootNodeState(0,0,0,initialDeque);
  addNode(rootNodeState);
  rootNodeIndex = 0;
  VRPTWNodeState terminalNodeState(vrptw.capacity,vrptw.endTimes[0],0,initialDeque);
  addNode(terminalNodeState);
  terminalNodeIndex = 1;

  // create nodes except r/t
  int nodeIndex = 0;
  while (nodeIndex < nodes.size())
  {
    const VRPTWNode& node = nodes[nodeIndex];
    for (int location=1; location<vrptw.numLocations; ++location)
    {
      int newCapacity = node.state.capacity + vrptw.demands[location];
      if (newCapacity > vrptw.capacity)
      {
        continue;
      }
      if (std::find(node.state.visited.begin(), node.state.visited.end(), location) != node.state.visited.end())
      {
        continue;
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

      int newTimeTimesTen = 0; // fix this
      VRPTWNodeState newState(newCapacity, newTimeTimesTen, location, newVisited);
      int newNodeIndex = addNode(newState);
      addArc(nodeIndex, newNodeIndex);
    }

    nodeIndex = nodeIndex + 1;
  }

  // add arcs to terminal node
  nodeIndex = 2;
  while (nodeIndex < nodes.size())
  {
    addArc(nodeIndex, terminalNodeIndex);
    nodeIndex = nodeIndex + 1;
  }

  DBG(print();)
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
  for (auto capNodeIndices : capacityNodeIndices)
  {
    // skip to load range
    if (capNodeIndices.first < load)
    {
      continue;
    }
    else if (capNodeIndices.first > demandLimit)
    {
      return;
    }

    // mark seen nodes and try to find nodeIndex2
    for (int nodeIndex : capNodeIndices.second)
    {
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
    for (auto capNodeIndices : capacityNodeIndices)
    {
      // skip to load range
      if (capNodeIndices.first < load1)
      {
        continue;
      }
      else if (capNodeIndices.first > load2)
      {
        return false;
      }

      // mark seen nodes and try to find nodeIndex2
      for (int nodeIndex : capNodeIndices.second)
      {
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
    for (auto capNodeIndices : capacityNodeIndices)
    {
      // skip to load range
      if (capNodeIndices.first < load2)
      {
        continue;
      }
      else if (capNodeIndices.first > load1)
      {
        return false;
      }

      // mark seen nodes and try to find nodeIndex2
      for (int nodeIndex : capNodeIndices.second)
      {
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
  for (auto capNodeIndices : capacityNodeIndices)
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
  for (auto capNodeIndices : capacityNodeIndices)
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

double VRPTWDecisionDiagram::setupAndSolveFlowModel(FlowType flowType, IncludeCoverConstraints includeCoverConstraints, UseColumnGeneration useCg, std::vector<double>& duals, double& truckDual, std::vector<double>& capDuals, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  // setup model
  IloEnv env;
  IloModel flowModel(env);
  IloRangeArray coverConstraints(env);
  IloRangeArray flowConservationConstraints(env);
  IloRangeArray capConstraints(env);
  IloRangeArray combConstraints(env);
  IloRangeArray srcConstraints(env);
  IloNumVarArray x(env, arcs.size());
  IloExpr objective(env);

  // setup variables and objective
  int numNonZeroVars = 0;
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    if (!arcs[arcIndex].isReverseArc && (fixedArcs.find(arcIndex) == fixedArcs.end()))
    {
      if (flowType == FlowType::IP)
      {
        x[arcIndex] = IloNumVar(env, 0, 1, ILOINT);
        objective += x[arcIndex] * arcs[arcIndex].coeff;
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

      coverConstraints.add(sumLocationArcs >= 1);
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

  // fix arcs
  for (auto arcIndex : fixedArcs)
  {
    flowModel.add(x[arcIndex] == 0);
  }

  if (flowType == FlowType::LP)
  {
    // RCC - rounded capacity cuts
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
        DBG(
          for (int dualIndex=0; dualIndex<vrptw.numLocations; ++dualIndex)
          {
            std::cout << "cover dual [" << dualIndex << "]: " << lpDuals[dualIndex] << std::endl;
          }
        )
      }

      IloNumArray lpCapDuals(env);
      solver.getDuals(lpCapDuals, capConstraints);
      capDuals.resize(capCutSets.size());
      for (int dualIndex=0; dualIndex<capCutSets.size(); ++dualIndex)
      {
        capDuals[dualIndex] = lpCapDuals[dualIndex];
      }
      for (int dualIndex=0; dualIndex<capCutSets.size(); ++dualIndex)
      {
        std::cout << "cap dual [" << dualIndex << "]: " << capDuals[dualIndex] << std::endl;
      }
   
      IloNumArray lpCombDuals(env);
      solver.getDuals(lpCombDuals, combConstraints);
      combDuals.resize(teeths.size());
      for (int dualIndex=0; dualIndex<teeths.size(); ++dualIndex)
      {
        combDuals[dualIndex] = lpCombDuals[dualIndex];
      }
      for (int dualIndex=0; dualIndex<teeths.size(); ++dualIndex)
      {
        std::cout << "comb dual [" << dualIndex << "]: " << combDuals[dualIndex] << std::endl;
      }

      IloNumArray lpSrcDuals(env);
      solver.getDuals(lpSrcDuals, srcConstraints);
      srcDuals.resize(cliqueCuts.size());
      for (int dualIndex=0; dualIndex<cliqueCuts.size(); ++dualIndex)
      {
        srcDuals[dualIndex] = lpSrcDuals[dualIndex];
        std::cout << "src dual [" << dualIndex << "]: " << srcDuals[dualIndex] << std::endl;
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

double VRPTWDecisionDiagram::fixArcs(const std::vector<double>& lambda, const double& truckDual, const std::vector<double>& capDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, double lowerBound, LPSolveType solveType)
{
  setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(lambda, capDuals, combDuals, srcDuals, solveType);

  std::vector<double> shortestPathDown(nodes.size(), INF);
  shortestPathDown[rootNodeIndex] = 0;
  std::vector<double> shortestPathUp(nodes.size(), INF);
  shortestPathUp[terminalNodeIndex] = 0;

  // dp to find shortest path down
  for (auto capNodeIndices : capacityNodeIndices)
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
  for (auto it=capacityNodeIndices.rbegin(); it!=capacityNodeIndices.rend(); ++it)
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

  // TODO: try saving some of these, or using multiple duals
  // fix arcs based on lb + rc > ub
  for (int arcIndex=0; arcIndex<arcs.size(); ++arcIndex)
  {
    const VRPTWArc& arc = arcs[arcIndex];
    double bestPossibleReducedCost = shortestPathDown[arc.fromNodeIndex] + shortestPathUp[arc.toNodeIndex] - truckDual + arc.coeff;
    if ((lowerBound + bestPossibleReducedCost) > (vrptw.hgsUpperBound + 0.00001))
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
  std::cout << "number of fixed arcs: " << fixedArcs.size() << std::endl;
  std::cout << "percent fixed arcs: " << getPercentFixedArcs() << std::endl;
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
  int firstLocation = nodes[arcs[routeArcs[0]].fromNodeIndex].state.lastVisited;
  int lastLocation = arcs[routeArcs.back()].location;
  if (firstLocation != lastLocation)
  {
    return false;
  }

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
    routeIndex = routeIndex + 1;
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
      return false;
    }
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

          if (decompositionReason == DecompositionReason::SEPARATE)
          {
            for (int lastLocationIndex=0; lastLocationIndex<lastLocationsVisited.size(); ++lastLocationIndex)
            {
              // if we have already visited, return the route of arcs causing the infeasibility
              if (lastLocationsVisited[lastLocationIndex] == arcs[arcIndex].location)
              {
                for (int index=lastLocationIndex+1; index<lastLocationsVisited.size(); ++index)
                {
                  routeArcs.push_back(lastArcsVisited[index]);
                }
                routeArcs.push_back(arcIndex);
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
      DBG(std::cout << arcIndex << " ";)
      arcs[arcIndex].decompositionFlow = arcs[arcIndex].decompositionFlow - routeFlow;
    }
    DBG(std::cout << std::endl;)
  }

  print();
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
  // clear when dd changes
  for (int arcIndex : fixedArcs)
  {
    VRPTWArc& arc = arcs[arcIndex];
    nodes[arc.fromNodeIndex].outArcs.push_back(arcIndex);
    nodes[arc.toNodeIndex].inArcs.push_back(arcIndex);
  }
  fixedArcs.clear();

  DBG(std::cout << "separating route: ";
  for (int arcIndex : routeArcs)
  {
    std::cout << arcIndex << " ";
  }
  std::cout << std::endl;)

  int currNodeIndex = arcs[routeArcs[0]].fromNodeIndex;
  for (int index=0; index<(routeArcs.size()-1); ++index)
  {
    int nextLocation = arcs[routeArcs[index]].location;

    // new node on next layer
    int nextNodeIndex = -1;
    bool newNodeCreated = false;
    int prevNumNodes = nodes.size();
    VRPTWNodeState newState(nodes[currNodeIndex].state);
    newState.visited.insert(nextLocation);
    /*
    if (newState.visited.size() > maxS)
    {
      newState.visited.pop_front();
    }
    */
    // NOTE(akarahal) special for experiments with Q
    if (maxS == 2)
    {
      if (newState.visited.size() > maxS)
      {
        for (int loc : nodes[currNodeIndex].state.visited)
        {
          if (loc != nodes[currNodeIndex].state.lastVisited)
          {
            newState.visited.erase(loc);
          }
        }
      }
    }
    newState.lastVisited = nextLocation;
    newState.capacity = newState.capacity + vrptw.demands[nextLocation];
    nextNodeIndex = addNode(newState);
    if (prevNumNodes < nodes.size())
    {
      newNodeCreated = true;
    }
 
    // cases (new node created, or already existed)
    if (newNodeCreated)
    {
      // copy possible outgoing arcs except next one on route that we know we'll switch
      for (int oldArcIndex : nodes[arcs[routeArcs[index]].toNodeIndex].outArcs)
      {
        VRPTWArc& oldArc = arcs[oldArcIndex];
        if ((std::find(newState.visited.begin(), newState.visited.end(), oldArc.location) == newState.visited.end()) && (oldArc.location != arcs[routeArcs[index+1]].location))
        {
          int newArcIndex = addArc(nextNodeIndex, oldArc.toNodeIndex);
          arcsUsed[newArcIndex] = true;
        }
      }
    }

    // move first arc
    if (index == 0)
    {
      int arcToMoveIndex = routeArcs[index];
      moveArc(arcToMoveIndex, nextNodeIndex);
    }
    else
    {
      if (newNodeCreated)
      {
        // add arc from curr state to this new node
        int newArcIndex = addArc(currNodeIndex, nextNodeIndex);
        arcsUsed[newArcIndex] = true;
      }
      else
      {
        // arc exists or may need to move arc
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
              moveArc(currNodeArcIndex, nextNodeIndex);
            }
          }
        }
      }
    }

    currNodeIndex = nextNodeIndex;
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

double VRPTWDecisionDiagram::solveMinCostFlowModel(const std::vector<double>& duals, std::vector<std::vector<int>>& shortestPathsByArc)
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
      computeShortestPathBFS(ShortestPathMode::UPDATE_POTENTIALS, empty);
      DBG(std::cout << "computed shortest path bfs" << std::endl;)
      DBG(print();)
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

double VRPTWDecisionDiagram::solveMinCostFlowModelWang(const std::vector<double>& duals, const std::vector<double>& capDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, std::vector<std::vector<int>>& shortestPathsByArc, bool& isDualFeasible, double& truckDual)
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

  // setup some data structures (fllow Wang,Wang,Wang NeurIps paper)
  std::vector<int> treeByParentArcs;
  treeByParentArcs.resize(nodes.size());
  std::vector<int> shortestPathByArc;
  setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(duals, capDuals, combDuals, srcDuals, LPSolveType::LAGSolver);
  double shortestPathLength = computeShortestPathBFSWang(treeByParentArcs, shortestPathByArc);
  // NOTE: when fixing arcs, if there is no path to terminal, no paths exist
  if (shortestPathLength == INF)
  {
    return 0.0;
  }

  // shortest path tells us dual feasibility
  if (shortestPathLength < 0)
  {
    isDualFeasible = false;
  }
  else
  {
    isDualFeasible = true;
    truckDual = shortestPathLength;
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
    if ((newPathArcIndex != -1) && (nodesToUpdate.find(arcs[newPathArcIndex].fromNodeIndex) == nodesToUpdate.end()))
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

bool VRPTWDecisionDiagram::checkFeasibleDual(const std::vector<double>& lambda, const std::vector<double>& rccDuals, const std::vector<double>& combDuals, std::vector<double>& srcDuals, LPSolveType solveType, double& truckDual)
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

  // otherwise find associated truck dual
  truckDual = shortestPathDistance;
  return true;
};

bool VRPTWDecisionDiagram::checkAn33k5SolutionPossible() const
{
  std::vector<std::vector<int>> routesByLocation;
  routesByLocation.push_back({0,2,0});
  routesByLocation.push_back({0,12,5,27,30,5,27,12,0});
  routesByLocation.push_back({0,20,32,13,26,7,8,13,32,20,0});
  routesByLocation.push_back({0,11,19,1,14,21,1,14,21,1,14,19,11,0});
  routesByLocation.push_back({0,16,3,9,17,3,16,29,0});
  routesByLocation.push_back({0,24,6,19,21,1,31,18,0});
  routesByLocation.push_back({0,5,25,30,27,25,5,4,0});
  routesByLocation.push_back({0,23,18,31,29,15,22,0});
  routesByLocation.push_back({0,11,19,14,21,1,14,21,1,14,21,1,29,0}); // index 8
  routesByLocation.push_back({0,12,10,17,9,16,15,0});
  routesByLocation.push_back({0,4,5,27,25,30,10,0});
  routesByLocation.push_back({0,4,5,27,25,30,12,4,0});
  routesByLocation.push_back({0,11,31,1,21,14,19,6,24,0});
  routesByLocation.push_back({0,24,6,19,14,1,21,14,19,11,0});
  routesByLocation.push_back({0,20,4,12,10,17,9,29,0});
  routesByLocation.push_back({0,30,25,27,30,25,12,0});
  routesByLocation.push_back({0,22,15,29,31,18,28,0});
  routesByLocation.push_back({0,22,23,18,28,23,22,0});
  routesByLocation.push_back({0,2,32,13,8,7,26,5,12,0});
  routesByLocation.push_back({0,32,13,8,7,26,8,7,26,20,0}); // index 19
  routesByLocation.push_back({0,11,19,21,14,1,21,19,6,24,0});
  routesByLocation.push_back({0,29,3,9,17,3,16,29,0});
  routesByLocation.push_back({0,22,23,11,6,24,2,0});
  routesByLocation.push_back({0,23,28,18,31,28,23,0});
  routesByLocation.push_back({0,2,24,6,11,18,28,0});
  routesByLocation.push_back({0,16,3,9,17,3,15,0});
  routesByLocation.push_back({0,12,27,25,30,27,12,0});
  routesByLocation.push_back({0,4,5,30,25,27,5,4,0});
  for (int routeIndex=0; routeIndex<routesByLocation.size(); ++routeIndex)
  {
    if (!doesRouteExistByLocations(routesByLocation[routeIndex]))
    {
      std::cout << "failed at index: " << routeIndex << std::endl;
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

