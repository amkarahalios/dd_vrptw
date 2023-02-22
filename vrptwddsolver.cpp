#include "vrptwddsolver.h"
#include <math.h>

VRPTWDDSolver::VRPTWDDSolver(VRPTW _vrptw, LPSolveType _lpSolveType, InitialStateSpace initialStateSpace, int _s, int _maxS, bool _useCuts, int _timeoutSeconds) : vrptw(_vrptw), routeDD(_vrptw, 100, 1 * (*std::min_element(vrptw.demands.begin()+1, vrptw.demands.end()))), maxS(_maxS), useCuts(_useCuts), timeoutSeconds(_timeoutSeconds), lpSolveType(_lpSolveType)
{
  for (int location=0; location<vrptw.numLocations; ++location)
  {
    //stats.upperBound += (vrptw.distances[location][0] * 2);
    stats.upperBound = vrptw.hgsUpperBound;
  }

  auto startCompileTime = std::chrono::high_resolution_clock::now();
  if (initialStateSpace == InitialStateSpace::Q)
  {
    routeDD.compileExactFukasawa(_s);
  } 
  else if (initialStateSpace == InitialStateSpace::NG)
  {
    routeDD.compileNgRoute(_s);
  }

  auto endCompileTime = std::chrono::high_resolution_clock::now();
  auto compileSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endCompileTime - startCompileTime).count();
  stats.millisecondsCompiling = stats.millisecondsCompiling + compileSolveTime;
  std::cout << "done compiling dd" << std::endl;

  CMGR_CreateCMgr(&MyCutsCMP,Dim);
  CMGR_CreateCMgr(&MyOldCutsCMP,Dim);
};

void VRPTWDDSolver::convertArcIndicesForVRPTWSep(const std::vector<double>& routeFlows,
                                               const std::vector<std::set<int>>& decomposedRoutes,
                                               std::vector<int>& edgeTail,
                                               std::vector<int>& edgeHead,
                                               std::vector<double>& edgeFlow)
{
  double totalFlow = 0.0;
  for (double flow : routeFlows)
  {
    totalFlow += flow;
  }

  std::map<std::pair<int,int>,double> edgeFlows;
  for (int routeIndex=0; routeIndex<decomposedRoutes.size(); ++routeIndex)
  {
    auto route = decomposedRoutes[routeIndex];
    for (int arcIndex : route)
    {
      std::pair<int,int> fromAndToIndices = routeDD.getFromAndToLocations(arcIndex);
      int currLoc = (fromAndToIndices.first == 0) ? vrptw.numLocations : fromAndToIndices.first;
      int nextLoc = (fromAndToIndices.second == 0) ? vrptw.numLocations : fromAndToIndices.second;
      if (currLoc < nextLoc)
      {
        edgeFlows[std::make_pair(currLoc,nextLoc)] = edgeFlows[std::make_pair(currLoc,nextLoc)] + (routeFlows[routeIndex]) / totalFlow;
      }
      else
      {
        edgeFlows[std::make_pair(nextLoc,currLoc)] = edgeFlows[std::make_pair(nextLoc,currLoc)] + (routeFlows[routeIndex]) / totalFlow;
      }
    }
  }

  edgeTail.push_back(0);
  edgeHead.push_back(0);
  edgeFlow.push_back(0);
  for (auto edge : edgeFlows)
  {
    DBG(std::cout << "t: " << edge.first.first << " h: " << edge.first.second << " flow: " << edge.second << std::endl;)
    edgeTail.push_back(edge.first.first);
    edgeHead.push_back(edge.first.second);
    edgeFlow.push_back(edge.second);
  }
};

void VRPTWDDSolver::addRCCs(const std::vector<std::vector<int>>& routes, bool& cutAdded)
{
  // Rounded Capacity Cuts when capacity is relaxed
  for (auto route : routes)
  {
    std::set<int> locations;
    int load = 0;
    for (int loc : route)
    {
      auto inserted = locations.insert(loc);
      if (inserted.second)
      {
        load = load + vrptw.demands[loc];
      }
      if (load > vrptw.capacity)
      {
        std::cout << "new capacity cutset: ";
        std::vector<int> cutSet;
        for (int loc : locations)
        {
          if (loc != 0)
          {
            cutSet.push_back(loc);
            std::cout << loc << " ";
          }
        }
        int RHS = cutSet.size() - ceil(load * 1.0 / vrptw.capacity);
        std::cout << "<= " << RHS << std::endl;
        routeDD.addCapCutSet(cutSet);
        routeDD.addCapCutSetRHS(RHS);
        stats.numCuts = stats.numCuts + 1;
        cutAdded = true;
        break;
      }
    }
  }

  if (!cutAdded)
  {
    std::cout << "no rcc cuts found" << std::endl;
  }
}

void VRPTWDDSolver::addRCCs(const std::vector<int>& edgeTail, const std::vector<int>& edgeHead, const std::vector<double>& edgeFlow, bool& cutAdded)
{
  // Rounded Capacity Cuts
  char integerAndFeas = '0';
  double epsForViolation = 0.1;
  double maxViolation = 0.0;
  int maxNoCuts = 10;
  CAPSEP_SeparateCapCuts(vrptw.numLocations-1,
                         &(vrptw.demandsForSeparation[0]),
                         vrptw.capacity,
                         static_cast<int>(edgeTail.size()-1), // num edges in solution
                         &(edgeTail[0]), // tails of edges
                         &(edgeHead[0]), // heads of edges
                         &(edgeFlow[0]), // flow value of edges
                         MyOldCutsCMP,
                         maxNoCuts, // max cuts to be returned
                         epsForIntegrality,
                         epsForViolation,
                         &integerAndFeas, // returned by method, 1 means int sol
                         &maxViolation, // violation of the cut with largest violation
                         MyCutsCMP); // contains cut
  if (integerAndFeas)
  {
    std::cout << "integral solution" << std::endl;

    // may still need to add cut if we relaxed capacity constraint
    // want to check each route that it's under capacity
    return;
  }

  int numCuts = MyCutsCMP->Size;
  if (numCuts > 0)
  {
    std::cout << "max violation: " << maxViolation << std::endl;
    cutAdded = true;
    for (int cutIndex=0; cutIndex<numCuts; ++cutIndex)
    {
      std::cout << "new capacity cutset: ";
      std::vector<int> cutSet;
      int cutsetSize = MyCutsCMP->CPL[cutIndex]->IntListSize;
      for (int locationIndex=1; locationIndex<=cutsetSize; ++locationIndex)
      {
        int location = MyCutsCMP->CPL[cutIndex]->IntList[locationIndex];
        std::cout << location << " ";
        cutSet.push_back(location);
      }
      double RHS = MyCutsCMP->CPL[cutIndex]->RHS;
      std::cout << "<= " << RHS << std::endl;

      routeDD.addCapCutSet(cutSet);
      routeDD.addCapCutSetRHS(RHS);
    }
    stats.numCuts = stats.numCuts + numCuts;

    for (int cutIndex=0; cutIndex<numCuts; ++cutIndex)
    {
      CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP,cutIndex,0);
    }
    MyCutsCMP->Size = 0;
  }
  else
  {
    std::cout << "no rcc cuts found" << std::endl;
  }
}

void VRPTWDDSolver::addCombs(std::vector<int>& edgeTail, std::vector<int>& edgeHead, std::vector<double>& edgeFlow, bool& cutAdded)
{
  // Add Comb Cuts
  int qMin = 0;
  int maxNoCuts = 5;
  double maxViolation = 0.0;
  COMBSEP_SeparateCombs(vrptw.numLocations-1,
                         &(vrptw.demandsForCombs[0]),
			 vrptw.capacity,
                         qMin,
                         static_cast<int>(edgeTail.size()-1), // num edges in solution
                         &(edgeTail[0]), // tails of edges
                         &(edgeHead[0]), // heads of edges
                         &(edgeFlow[0]), // flow value of edges
                         maxNoCuts, // max cuts to be returned
                         &maxViolation, // violation of the cut with largest violation
                         MyCutsCMP); // contains cut

  int numCuts = MyCutsCMP->Size;
  if (numCuts > 0)
  {
    std::cout << "max violation: " << maxViolation << std::endl;
    cutAdded = true;
    for (int cutIndex=0; cutIndex<numCuts; ++cutIndex)
    {
      int numTeeth = MyCutsCMP->CPL[cutIndex]->Key;
      std::vector<std::set<int>> teeth;
      std::cout << "new strengthened comb with " << numTeeth << " teeth" << std::endl;

      std::set<int> handle;
      std::cout << "handle: ";
      for (int rowIndex=1; rowIndex<=MyCutsCMP->CPL[cutIndex]->IntListSize; ++rowIndex)
      {
        int row = MyCutsCMP->CPL[cutIndex]->IntList[rowIndex];
        if (row == vrptw.numLocations)
        {
          row = 0;
        }
	handle.insert(row);
        std::cout << row << ",";
      }
      std::cout << std::endl;
      teeth.push_back(handle);

      for (int toothIndex=1; toothIndex<=numTeeth; ++toothIndex)
      {
        std::set<int> tooth;
        std::cout << "tooth: ";
        int minIndex = MyCutsCMP->CPL[cutIndex]->ExtList[toothIndex];
	int maxIndex = -1;
	if (toothIndex == numTeeth)
	{
          maxIndex = MyCutsCMP->CPL[cutIndex]->ExtListSize;
	}
	else
	{
          maxIndex = MyCutsCMP->CPL[cutIndex]->ExtList[toothIndex+1] - 1;
	}

	for (int nodeIndex=minIndex; nodeIndex<=maxIndex; ++nodeIndex)
	{
          int row = MyCutsCMP->CPL[cutIndex]->ExtList[nodeIndex];
          if (row == vrptw.numLocations)
          {
            row = 0;
          }
	  tooth.insert(row);
          std::cout << row << ",";
	}
        std::cout << std::endl;

	teeth.push_back(tooth);
      }
      routeDD.addCombCutTeeth(teeth);

      int RHS = MyCutsCMP->CPL[cutIndex]->RHS;
      std::cout << ">= " << RHS << std::endl;
      routeDD.addCombCutRHS(RHS);
    }
    stats.numCuts = stats.numCuts + numCuts;

    for (int cutIndex=0; cutIndex<numCuts; ++cutIndex)
    {
      CMGR_MoveCnstr(MyCutsCMP, MyOldCutsCMP,cutIndex,0);
    }
    MyCutsCMP->Size = 0;
  }
  else
  {
    std::cout << "no comb cuts found" << std::endl;
  }
}

void VRPTWDDSolver::addSRCCuts(std::vector<double>& srcDuals)
{
  for (int index=srcDuals.size(); index<routeDD.getNumCliqueCuts(); ++index)
  {
    srcDuals.push_back(0);
  }
}

bool VRPTWDDSolver::solve()
{
  // for Lagrangean method, the choice of initial labmda is important
  // choice 1: simple all 0
  // choice 2: optimal value will be min distance from depot
  // choice 3: optimal value will be distances from depot
  //std::vector<double> lambda(vrptw.numLocations, 0);
  //std::vector<double> lambda(vrptw.numLocations, *std::min_element(vrptw.distances[0].begin()+1, vrptw.distances[0].end()));
  bool lpFlowType = FlowType::LP;
  std::vector<double> lambda;
  lambda.push_back(0);
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    lambda.push_back(2 * vrptw.distances[0][location] * vrptw.demands[location] / vrptw.capacity);
  }
  // so arc fixing does not try to use IP achieved lb
  double bestLPValue = 0.0;
  bestLambdaArcFixing.resize(lambda.size());
  std::vector<double> mu;
  std::vector<double> combDuals;
  std::vector<double> srcDuals;

  bool changedLagToLP = false;
  bool finishedSolving = false;
  while (!finishedSolving)
  {
    bool solved = false;
    if (lpSolveType == LPSolveType::LPSolver)
    {
      if (stats.lpIterations > 1)
      {
        double percentArcsFixed = 0.0;
        if (!changedLagToLP)
        {
          if (vrptw.oneOrMorePaths != OneOrMorePaths::ONE_PATH)
          {
            percentArcsFixed = routeDD.fixArcs(lambda, mu, combDuals, srcDuals, bestLPValue, lpSolveType);
          }
        }
        if ((percentArcsFixed < bestLambdaPercentFixed) || changedLagToLP)
        {
          if (changedLagToLP)
          {
            changedLagToLP = false;
          }
          std::vector<double> emptyMu;
          std::vector<double> emptyComb;
          std::vector<double> emptySrc;
          routeDD.fixArcs(bestLambdaArcFixing, emptyMu, emptyComb, emptySrc, bestLambdaLowerBound, lpSolveType);
        }
      }

      if (lpFlowType == FlowType::LP)
      {
        solved = solveLP(lambda, mu, combDuals, srcDuals);
        bestLPValue = stats.lowerBound;
      }
      else
      {
        solved = solveIP(lambda, mu, combDuals, srcDuals);
      }
    }
    else if (lpSolveType == LPSolveType::LAGSolver)
    {
      solved = solveLagrangeanRelaxation(lambda, mu, combDuals, srcDuals);
      // if switching to LP, go right to next iter and use cuts
      if (lpSolveType == LPSolveType::LPSolver)
      {
        changedLagToLP = true;
        useCuts = true;
        continue;
      }
    }

    if (!solved)
    {
      std::cout << "finished solving - error" << std::endl;
      return false;
    }

    // get primal heuristic flow decomp or just a flow decomp for the solution
    std::vector<std::vector<int>> routesByLocationPrimalHeuristic;
    //routeDD.primalHeuristic(routesByLocationPrimalHeuristic);
    //double currentUpperBound = vrptw.evaluateSolutionCost(routesByLocationPrimalHeuristic);
    //if (currentUpperBound < stats.upperBound)
    //{
    //  stats.upperBound = currentUpperBound;
    //}

    // infeasibilities for LP solve
    if (lpSolveType == LPSolveType::LPSolver)
    {
      bool infeasibilityFound = false;
      bool cutAdded = false;

      // do cuts first in case separations mess up dd structure
      if (useCuts && (lpFlowType == FlowType::LP))
      {
        // RCC - rounded capacity cuts
        std::vector<int> edgeTail;
        std::vector<int> edgeHead;
        std::vector<double> edgeFlow;
        routeDD.convertSolutionForVRPTWSep(edgeTail, edgeHead, edgeFlow);
        addRCCs(edgeTail, edgeHead, edgeFlow, cutAdded);
      }

      bool stopFindingInfeasibilities = false;
      std::vector<std::vector<int>> infeasibilities;
      std::vector<double> routeFlows;
      std::vector<std::vector<int>> decomposedRoutes;
      while (!stopFindingInfeasibilities)
      {
        std::vector<int> infeasibleRoute;
        routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, maxS, DecompositionReason::SEPARATE);
        if (!infeasibleRoute.empty())
        {
          infeasibilities.push_back(infeasibleRoute);
          infeasibilityFound = true;
        }
        else
        {
          stopFindingInfeasibilities = true;
        }
      }

      if (infeasibilityFound)
      {
        // also check capacity feasibility
        bool allRoutesCapacityFeasible = true;
        for (auto route : decomposedRoutes)
        {
          int load = 0;
          for (int index=0; index<route.size()-1; ++index)
          {
            int loc = route[index];
            load = load + vrptw.demands[loc];
          }
          if (load > vrptw.capacity)
          {
            allRoutesCapacityFeasible = false;
          }
        }
        // RCC - rounded capacity cuts when dd relaxes capacity
        if (useCuts && (vrptw.vrptwType != VRPTWType::NO_RELAX_CAPACITY) && !allRoutesCapacityFeasible)
        {
          addRCCs(decomposedRoutes, cutAdded);
        }

        for (auto infeasibleRoute : infeasibilities)
        {
          if (routeDD.doesRouteExistByArcs(infeasibleRoute))
          {
            stats.numSeparations = stats.numSeparations + 1;
            routeDD.separateInfeasibleRoute(infeasibleRoute, maxS);
          }
        }
      }
      else
      {
        // print decomposed routes and check capacity feasibility
        std::cout << "routes: " << std::endl;
        bool allRoutesCapacityFeasible = true;
        int index = 0;
        for (auto route : decomposedRoutes)
        {
          int load = 0;
          std::cout << "flow {" << routeFlows[index] << "}: ";
          for (int index=0; index<route.size()-1; ++index)
          {
            int loc = route[index];
            int nextLoc = route[index+1];
            load = load + vrptw.demands[loc];
            std::cout << loc << " ";
          }
          std::cout << std::endl;
          index = index + 1;
          if (load > vrptw.capacity)
          {
            allRoutesCapacityFeasible = false;
          }
        }
        std::cout << "no more separations possible" << std::endl;

        // RCC - rounded capacity cuts when dd relaxes capacity
        if (useCuts && (vrptw.vrptwType != VRPTWType::NO_RELAX_CAPACITY) && !allRoutesCapacityFeasible)
        {
          addRCCs(decomposedRoutes, cutAdded);
        }
      }

      if (!cutAdded && !infeasibilityFound)
      {
        std::cout << "no more separations or cuts possible" << std::endl;
        if (stats.lowerBound != stats.upperBound)
        {
          if (lpFlowType == FlowType::LP)
          {
            std::cout << "LP solved - switching to IP" << std::endl;
            lpFlowType = FlowType::IP;
          }
          else
          {
            stats.upperBound = stats.lowerBound;
          }
        }
        else
        {
          finishedSolving = true;
        }
      }
    }

    if (stats.getNumSeconds() >= timeoutSeconds)
    {
      finishedSolving = true;
    }

    if (stats.lowerBound == stats.upperBound)
    {
      finishedSolving = true;
    }
  }

  stats.print(routeDD.size());
  std::cout << "finished solving - complete" << std::endl;
  CMGR_FreeMemCMgr(&MyCutsCMP);
  CMGR_FreeMemCMgr(&MyOldCutsCMP);
  return true;
}

bool VRPTWDDSolver::solveIP(std::vector<double>& lambda, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  routeDD.setCoeffsAsDistances();
  auto startLPTime = std::chrono::high_resolution_clock::now();
  std::cout << "solving IP" << std::endl;
  stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, lambda, mu, combDuals, srcDuals);
  auto endLPTime = std::chrono::high_resolution_clock::now();
  auto lpSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endLPTime - startLPTime).count();
  stats.millisecondsSolvingLP = stats.millisecondsSolvingLP + lpSolveTime;

  stats.lpIterations = stats.lpIterations + 1;
  stats.print(routeDD.size());
  return true;
};

bool VRPTWDDSolver::solveLP(std::vector<double>& lambda, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  routeDD.setCoeffsAsDistances();
  auto startLPTime = std::chrono::high_resolution_clock::now();
  bool MIPOn = false;
  if ((routeDD.getPercentFixedArcs() >= 95) && MIPOn)
  {
    std::cout << "solving IP" << std::endl;
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, lambda, mu, combDuals, srcDuals);
  }
  else
  {
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, lambda, mu, combDuals, srcDuals);
  }
  auto endLPTime = std::chrono::high_resolution_clock::now();
  auto lpSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endLPTime - startLPTime).count();
  stats.millisecondsSolvingLP = stats.millisecondsSolvingLP + lpSolveTime;

  stats.lpIterations = stats.lpIterations + 1;
  stats.print(routeDD.size());
  return true;
};

bool VRPTWDDSolver::solvePricingProblem(std::vector<double>& lambda)
{
  routeDD.setCoeffsAsDistancesMinusLagrangean(lambda);

  std::vector<int> newRoute;
  double shortestPathLength = routeDD.computeShortestPathBFS(ShortestPathMode::SHORTEST_PATH, newRoute);
  shortestPathLength = shortestPathLength;
  //std::cout << "spl: " << shortestPathLength << std::endl;
  if (shortestPathLength < -0.00000001)
  {
    std::cout << "adding route: ";
    for (int loc : newRoute)
    {
      std::cout << loc << ",";
    }
    std::cout << std::endl;
    routeDD.addColumnForLPCG(newRoute);
    return true;
  }

  return false;
}

void VRPTWDDSolver::initializeColumns()
{
  // all single stop routes
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    std::vector<int> route;
    route.push_back(vrptw.depot);
    route.push_back(location);
    route.push_back(vrptw.depot);
    routeDD.addColumnForLPCG(route);
  }

  // some full routes
  std::set<int> locationsAdded;
  while (locationsAdded.size() < (vrptw.numLocations-1))
  {
    std::vector<int> route;
    route.push_back(0);
    int currentDemand = 0;
    for (int location=1; location<vrptw.numLocations; ++location)
    {
      if (locationsAdded.find(location) != locationsAdded.end())
      {
        continue;
      }

      if (currentDemand + vrptw.demands[location] <= vrptw.capacity)
      {
        route.push_back(location);
        locationsAdded.insert(location);
        currentDemand = currentDemand + vrptw.demands[location];
      }
    }
    route.push_back(0);
    DBG(std::cout << "adding initial route: ";)
    for (int loc : route)
    {
      std::cout << loc << ",";
    }
    std::cout << std::endl;
    routeDD.addColumnForLPCG(route);
  }
};

bool VRPTWDDSolver::solveLPCG(std::vector<double>& lambda, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  // column generation - RMP <-> pricing problem
  //initializeColumns();

  auto startLPTime = std::chrono::high_resolution_clock::now();
  bool solved = false;
  double flowObj = stats.lowerBound;
  while (!solved)
  {
    routeDD.setCoeffsAsDistances();
    flowObj = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::USE_CG, lambda, mu, combDuals, srcDuals);
    std::cout << "flowobj: " << flowObj << std::endl;

    bool addedColumn = solvePricingProblem(lambda);
    if (!addedColumn)
    {
      solved = true;
    }
  }

  auto endLPTime = std::chrono::high_resolution_clock::now();
  auto lpSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endLPTime - startLPTime).count();
  stats.millisecondsSolvingLP = stats.millisecondsSolvingLP + lpSolveTime;
  stats.lowerBound = flowObj;
  stats.lpIterations = stats.lpIterations + 1;
  stats.print(routeDD.size());
  return true;
};

void VRPTWDDSolver::updateMultipliers(std::vector<double>& lambda, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals, std::vector<double>& stepSizes, const std::set<int>& solutionArcs, double lagrangeanLowerBound, int k)
{
  // get values of LHS for all dualized inequalities
  std::unordered_map<int,double> locationsCovered;
  routeDD.getNumberOfTimesLocationsCovered(locationsCovered);

  std::vector<double> cutValues;
  routeDD.getCutSetValues(cutValues);
 
  std::vector<double> combValues;
  routeDD.getCombValues(combValues);

  std::vector<double> cliqueCutValues;
  routeDD.getCliqueCutValues(cliqueCutValues);

  // gamma_(k) = b - Ax_(k)
  // Beasley - when multiplier is already 0 and step direction is negative, don't include in ||gamma||^2
  double normGammaSquared = 0;
  std::vector<double> gamma(vrptw.numLocations + mu.size() + srcDuals.size() + combDuals.size(), 0);
  for (int i=1; i<vrptw.numLocations; ++i)
  {
    gamma[i] = 1 - locationsCovered[i];
    if ((lambda[i] <= 0.000001) && (gamma[i] <= 0.000001))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[i], 2);
    }
  }

  // cap cuts are in the form Cx <= r so change to -Cx >= -r
  for (int i=0; i<mu.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations;
    gamma[gammaIndex] = (-1 * routeDD.getCapCutSetRHS(i)) + cutValues[i];
    if ((mu[i] <= 0.01) && (gamma[gammaIndex] <= 0.01))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
  }

  // src cuts are in the form Cx <= r so chang eto -Cx >= -r
  for (int i=0; i<srcDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations + mu.size();
    gamma[gammaIndex] = -1 + cliqueCutValues[i];
    if ((srcDuals[i] <= 0.01) && (gamma[gammaIndex] <= 0.01))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
  }

  for (int i=0; i<combDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations + mu.size() + srcDuals.size();
    gamma[gammaIndex] = routeDD.getCombCutRHS(i) - combValues[i];
    if ((combDuals[i] <= 0.01) && (gamma[gammaIndex] <= 0.01))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
  }

  // eta_(k) = 0.05 * 100 / (100 + k)
  //double eta = 0.05 * 100 / (100 + stats.numLagIterations);
  double eta = 0.05 * 100 / (100 + k);

  // psi_(star) = psi_(best) * (1 + eta_(k))
  double psiStar = stats.lowerBound * (1 + eta);

  // alpha_(k) = (psi_(star) - psi(lambda(k))) / ||gamma_(k)||_(2)^2
  //double alpha = std::max(0.0001, (psiStar - lagrangeanLowerBound) / normGammaSquared);
  //std::cout << "lag iter: " << k << std::endl;
  //std::cout << "lag iter: " << stats.numLagIterations << std::endl;
  //std::cout << "lag bound: " << lagrangeanLowerBound << std::endl;
  double alpha = (psiStar - lagrangeanLowerBound) / normGammaSquared;
  stepSizes.push_back(alpha);

  DBG(
  std::cout << "alpha: " << alpha << std::endl;
  for (double g : gamma)
  {
    std::cout << g << ",";
  }
  std::cout << std::endl;

  for (double l : lambda)
  {
    std::cout << l << ",";
  }
  std::cout << std::endl;

  for (double m : mu)
  {
    std::cout << m << ",";
  }
  std::cout << std::endl;

  for (double s : srcDuals)
  {
    std::cout << s << ",";
  }
  std::cout << std::endl;
  for (double c : combDuals)
  {
    std::cout << c << ",";
  }
  std::cout << std::endl;)

  // lambda_(k+1) = lambda_(k) + alpha_(k) * gamma_(k)
  for (int i=1; i<vrptw.numLocations; ++i)
  {
    lambda[i] = std::max(0.0, lambda[i] + alpha * gamma[i]);
  }

  // same for mu
  for (int i=0; i<mu.size(); ++i)
  {
    mu[i] = std::max(0.0, mu[i] + alpha * gamma[i+vrptw.numLocations]);
  }
 
  // same for srcDuals
  for (int i=0; i<srcDuals.size(); ++i)
  {
    srcDuals[i] = std::max(0.0, srcDuals[i] + alpha * gamma[i+vrptw.numLocations+mu.size()]);
  }

  // same for combs
  for (int i=0; i<combDuals.size(); ++i)
  {
    combDuals[i] = std::max(0.0, combDuals[i] + alpha * gamma[i+vrptw.numLocations+mu.size()+srcDuals.size()]);
  }
};

bool VRPTWDDSolver::solveLagrangeanRelaxation(std::vector<double>& lambda, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  stats.lpIterations = stats.lpIterations + 1;
  bool shouldTerminate = false;
  int kappaIterations = 100;
  double muPercentImproved = 0.001;
  int lastMuImprovedIteration = 0;
  double muLowerBound = 0.0;
  double currIterLowerBound = 0.0;
  int numLagIterations = 0;
  std::vector<std::vector<int>> infeasibleRoutes;
  int infeasibleRoutesSizeLimit = 100;
  std::vector<double> stepSizes;
  std::vector<std::set<int>> xDecompositions;
  stats.print(routeDD.size());
  while (!shouldTerminate)
  {
    while (!shouldTerminate && (infeasibleRoutes.size() < infeasibleRoutesSizeLimit))
    {
      ++stats.numLagIterations;
      ++numLagIterations;

      std::vector<std::vector<int>> shortestPaths;
      auto startSSPTime = std::chrono::high_resolution_clock::now();
      DBG(std::cout << "start" << std::endl;)
      bool isDualFeasible = false;
      double minReducedCost = 0.0;
      //double lagrangeanLowerBound = routeDD.solveMinCostFlowModel(lambda, shortestPaths);
      double lagrangeanLowerBound = routeDD.solveMinCostFlowModelWang(lambda, mu, combDuals, srcDuals, shortestPaths, isDualFeasible, minReducedCost);
      DBG(std::cout << "finish" << std::endl;)
      stats.numSSPIterations = stats.numSSPIterations + shortestPaths.size();
      auto endSSPTime = std::chrono::high_resolution_clock::now();
      auto sspSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endSSPTime - startSSPTime).count();
      stats.millisecondsSolvingSSP = stats.millisecondsSolvingSSP + sspSolveTime;
      DBG(
      for (int index=0; index<lambda.size(); ++index)
      {
        std::cout << lambda[index] << ",";
      }
      std::cout << std::endl;)

      // dTx - lambda_T(Ax - b), so add sum of lambdas
      for (int index=0; index<lambda.size(); ++index)
      {
        lagrangeanLowerBound += lambda[index];
      }

      // - mu_T(-Cx + r) because capCuts are Cx <= r (not >=)
      for (int index=0; index<mu.size(); ++index)
      {
        lagrangeanLowerBound += mu[index] * routeDD.getCapCutSetRHS(index) * -1;
      }
 
      // - combDual_T(Ax - RHS)
      for (int index=0; index<combDuals.size(); ++index)
      {
        lagrangeanLowerBound += (combDuals[index] * routeDD.getCombCutRHS(index));
      }

      // - src_Duals_T(-Cx + r) because srcDuals are Cx <= r (not >=)
      for (int index=0; index<srcDuals.size(); ++index)
      {
        lagrangeanLowerBound += srcDuals[index] * -1;
      }

      DBG(
      for (int index=0; index<lambda.size(); ++index)
      {
        std::cout << index << ": " << lambda[index] << std::endl;
      }
      std::cout << "curr lb: " << lagrangeanLowerBound << std::endl;
      )

      // keep track of this iteration and overall
      if (currIterLowerBound < lagrangeanLowerBound)
      {
        if (((1 + muPercentImproved) * muLowerBound) < lagrangeanLowerBound)
        {
          lastMuImprovedIteration = numLagIterations;
          muLowerBound = lagrangeanLowerBound;
        }
        currIterLowerBound = lagrangeanLowerBound;
      }
      if (stats.lowerBound < lagrangeanLowerBound)
      {
        stats.lowerBound = lagrangeanLowerBound;
        stats.print(routeDD.size());
      }

      // get primal solution
      std::set<int> solutionArcs;
      routeDD.getSolutionArcs(solutionArcs);
      if ((stats.lpIterations >= 2) && useCuts)
      {
        xDecompositions.push_back(solutionArcs);
      }

      // need decomposition flow to work
      updateMultipliers(lambda, mu, combDuals, srcDuals, stepSizes, solutionArcs, lagrangeanLowerBound, numLagIterations);

      // check for cycles up to certain size and add to be separated
      //if ((maxS > 1) && (stats.lpIterations >= 2) && isSeparationRound)
      if (maxS > 1)
      {
        bool infeasibleRouteFound = true;
        while (infeasibleRouteFound)
        {
          std::vector<int> infeasibleRoute;
          std::vector<double> flows;
          std::vector<std::vector<int>> routes;
          routeDD.decomposeRoutes(infeasibleRoute, flows, routes, maxS, DecompositionReason::SEPARATE);
          if (!infeasibleRoute.empty())
          {
            infeasibleRoutes.push_back(infeasibleRoute);
          }
          else
          {
            infeasibleRouteFound = false;
          }
        }
      }

      // termination criteria for this iteration
      if ((numLagIterations - lastMuImprovedIteration) > kappaIterations)
      {
        shouldTerminate = true;
      }

      // overall termination criteria
      if (stats.getNumSeconds() >= timeoutSeconds)
      {
        shouldTerminate = true;
      }
 
      // should fix after decomposing in case we fix an arc that is in the current solution?
      // see if we can fix arcs based on a feasible dual
      // some rounds let's force dual feasibility to find the bound and fix some arcs?
      double percentFixed = 0.0;
      if (isDualFeasible)
      {
        if (vrptw.oneOrMorePaths != OneOrMorePaths::ONE_PATH)
        {
          percentFixed = routeDD.fixArcs(lambda, mu, combDuals, srcDuals, lagrangeanLowerBound, lpSolveType);
          if (percentFixed > bestLambdaPercentFixed)
          {
            bestLambdaPercentFixed = percentFixed;
            bestLambdaLowerBound = lagrangeanLowerBound;
            for (int index=0; index<lambda.size(); ++index)
            {
              bestLambdaArcFixing[index] = lambda[index];
            }
          }

          double checkLambdaLB = 0.0;
          for (double dual : lambda)
          {
            checkLambdaLB += dual;
          }
          std::cout << "check lb sum lambda: " << checkLambdaLB << std::endl;
        }
      }
      else
      {
        // can repair dual to make it feasible, and then try to fix arcs. Min reduced cost would be negative
        std::vector<double> repairedLambda(lambda);
        while (true)
        {
          routeDD.setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(repairedLambda, mu, combDuals, srcDuals, LPSolveType::LAGSolver);

          std::vector<int> treeByParentArcs;
          treeByParentArcs.resize(routeDD.getNodes().size());
          std::vector<int> shortestPathByArc;
          double shortestPathLength = routeDD.computeShortestPathBFSWang(treeByParentArcs, shortestPathByArc);
          //std::cout << "spl: " << shortestPathLength << std::endl;
          if (shortestPathLength > 0)
          {
            break;
          }
          else
          {
            std::vector<std::vector<int>> shortestPathsByArc;
            shortestPathsByArc.push_back(shortestPathByArc);
            std::set<int> locations;
            routeDD.getLocationsOnArcPaths(shortestPathsByArc, locations);
            for (int loc : locations)
            {
              repairedLambda[loc] = std::max(0.0, repairedLambda[loc] + (shortestPathLength / (shortestPathByArc.size() - 2)) - 0.000001);
            }
          }
        }

        double repairedBound = 0.0;
        for (int index=1; index<repairedLambda.size(); ++index)
        {
          repairedBound += repairedLambda[index];
        }
        std::cout << "repaired lb: " << repairedBound << std::endl;
        percentFixed = routeDD.fixArcs(repairedLambda, mu, combDuals, srcDuals, repairedBound, lpSolveType);
        if (percentFixed > bestLambdaPercentFixed)
        {
          bestLambdaPercentFixed = percentFixed;
          bestLambdaLowerBound = repairedBound;
          for (int index=0; index<lambda.size(); ++index)
          {
            bestLambdaArcFixing[index] = repairedLambda[index];
          }
        }

        if (percentFixed > 95)
        {
          for (int index=0; index<lambda.size(); ++index)
          {
            lambda[index] = repairedLambda[index];
          }
        }
      }

      if (percentFixed > 95)
      {
        lpSolveType = LPSolveType::LPSolver;
        stats.lpIterations = 1;
        shouldTerminate = true;
        routeDD.print();
      }
    }

    // add separations each round
    //if ((maxS > 1) && !shouldTerminate && isSeparationRound)
    if ((maxS > 1) && !shouldTerminate)
    {
      for (auto infeasibleRouteToSeparate : infeasibleRoutes)
      {
        if (routeDD.doesRouteExistByArcs(infeasibleRouteToSeparate))
        {
          stats.numSeparations = stats.numSeparations + 1;
          routeDD.separateInfeasibleRoute(infeasibleRouteToSeparate, maxS);
        }
      }
    }

    // add cutting planes at the end of 2
    if (shouldTerminate && useCuts && (stats.lpIterations >= 2))
    {
      // for Rounded Capacity Cuts
      bool cutAdded = false;
      std::vector<int> edgeTail;
      std::vector<int> edgeHead;
      std::vector<double> edgeFlow;
      convertArcIndicesForVRPTWSep(stepSizes, xDecompositions, edgeTail, edgeHead, edgeFlow);
      addRCCs(edgeTail, edgeHead, edgeFlow, cutAdded);
      mu.resize(routeDD.getNumCapCuts());

      // Strengthened Combs
      addCombs(edgeTail, edgeHead, edgeFlow, cutAdded);
      combDuals.resize(routeDD.getNumCombCuts());

      // for Subset Row Cuts add right away during cut rounds
      // dont add for X was taking too long
      // TOOD: maybe only add for small instances?
      routeDD.findSRCThree(xDecompositions, stepSizes, cutAdded);
      addSRCCuts(srcDuals);
      xDecompositions.clear();
      stepSizes.clear();

      // should we reset iteration count when we add more dual variables?
      numLagIterations = 0;
    }
    infeasibleRoutes.clear();
  }

  return true;
}
