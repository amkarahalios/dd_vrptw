#include "vrptwddsolver.h"
#include <math.h>

VRPTWDDSolver::VRPTWDDSolver(VRPTW _vrptw, VRPTWDDParameters _params) : vrptw(_vrptw), routeDD(_vrptw, _params), params(_params)
{
  for (int location=0; location<vrptw.numLocations; ++location)
  {
    //stats.upperBound += (vrptw.distances[location][0] * 2);
    stats.upperBound = vrptw.instanceUpperBound;
  }

  std::cout << "compiling DD" << std::endl;
  auto startCompileTime = std::chrono::high_resolution_clock::now();
  if (params.initialStateSpace == InitialStateSpace::Q)
  {
    routeDD.compileExactFukasawa(params.s);
  } 
  else if (params.initialStateSpace == InitialStateSpace::NG)
  {
    routeDD.compileNgRoute(params.s);
    //routeDD.checkLRC121SolutionPossible();
  }

  auto endCompileTime = std::chrono::high_resolution_clock::now();
  auto compileSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endCompileTime - startCompileTime).count();
  stats.millisecondsCompiling = stats.millisecondsCompiling + compileSolveTime;
  std::cout << "done compiling DD" << std::endl;

  // set and log parameters
  std::cout << "batch size for lag: " << params.infeasibleRoutesBatchSize << std::endl;
  std::cout << "iteration delay to start separating: " << params.lagIterationDelayToStartSeparating << std::endl;
  std::cout << "optimality gap to start repairing: " << params.lagOptimalityGapToStartRepairing << std::endl;
  std::cout << "percent fixed to change to CPLEX: " << params.percentFixedToChangeToCPLEX << std::endl;
  std::cout << "num arcs to change to CPLEX: " << params.numArcsToChangeToCPLEX << std::endl;
  std::cout << "num arcs to change to LAG: " << params.numArcsToChangeToLAG << std::endl;
  std::cout << "num lag iters to cut: " << params.numLagItersForCuts << std::endl;
  std::cout << "num cuts lag: " << params.numLagCuts << std::endl;
  std::cout << "cut phase: " << params.cutPhase << std::endl;
  std::cout << "deactivate cut value threshold: " << params.deactivateCutValueThreshold << std::endl;
  std::cout << "deactivate cut iter threshold: " << params.deactivateCutIterThreshold << std::endl;
  std::cout << "momentum beta: " << params.momentumBeta << std::endl;
  std::cout << "stall alpha factor: " << params.stallAlphaFactor << std::endl;

  if ((params.lpSolveType == LPSolveType::LPSolver) && (routeDD.getArcs().size() >= params.numArcsToChangeToLAG))
  {
    std::cout << "changing solve type to LAG" << std::endl;
    params.lpSolveType = LPSolveType::LAGSolver;
  }
  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());

  // set best lambda lower bound to 0
  bestLambdaLowerBound = 0.0;
  bestLambdaPercentFixed = 0.0;
  bestSinglePathDualFixing = 0.0;
  alphaFactor = 1.0;
  lambdaStore.resize(lambdaStoreSize);
  singlePathStore.resize(lambdaStoreSize);
  lambdaLowerBoundStore.resize(lambdaStoreSize);

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
        std::cout << "new capacity cutset for relaxation: ";
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

void VRPTWDDSolver::addRCCs(const std::vector<int>& edgeTail, const std::vector<int>& edgeHead, const std::vector<double>& edgeFlow, int maxNumCuts, bool& cutAdded)
{
  // Rounded Capacity Cuts
  char integerAndFeas = '0';
  double epsForViolation = 0.1;
  double maxViolation = 0.0;
  CAPSEP_SeparateCapCuts(vrptw.numLocations-1,
                         &(vrptw.demandsForSeparation[0]),
                         vrptw.capacity,
                         static_cast<int>(edgeTail.size()-1), // num edges in solution
                         &(edgeTail[0]), // tails of edges
                         &(edgeHead[0]), // heads of edges
                         &(edgeFlow[0]), // flow value of edges
                         MyOldCutsCMP,
                         maxNumCuts, // max cuts to be returned
                         epsForIntegrality,
                         epsForViolation,
                         &integerAndFeas, // returned by method, 1 means int sol
                         &maxViolation, // violation of the cut with largest violation
                         MyCutsCMP); // contains cut
  if (integerAndFeas)
  {
    std::cout << "integral solution - no cut added" << std::endl;

    // may still need to add cut if we relaxed capacity constraint
    // want to check each route that it's under capacity
    return;
  }

  int numCuts = MyCutsCMP->Size;
  if (numCuts > 0)
  {
    std::cout << "checking cuts, max violation: " << maxViolation << ", num cuts: " << numCuts << std::endl;
    cutAdded = true;
    std::vector<std::set<int>> cutSets;
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

      // ensure we don't add families all together, as this negatively affects subgradient descent
      bool inFamily = false;
      for (int index=0; index<cutSets.size(); ++index)
      {
        auto referenceCutSet = cutSets[index];
        if (std::includes(cutSet.begin(), cutSet.end(), referenceCutSet.begin(), referenceCutSet.end()) || std::includes(referenceCutSet.begin(), referenceCutSet.end(), cutSet.begin(), cutSet.end()))
        {
          inFamily = true;
        }
      }

      if (inFamily)
      {
        std::cout << "nested set, do not add" << std::endl;
      }
      else
      {
        std::set<int> cutSetAsSet(cutSet.begin(), cutSet.end());
        cutSets.push_back(cutSetAsSet);
        routeDD.addCapCutSet(cutSet);
        routeDD.addCapCutSetRHS(RHS);
        stats.numCuts = stats.numCuts + 1;
      }
    }

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

bool VRPTWDDSolver::solve(bool shouldSolveIP)
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
    if (vrptw.oneOrMorePaths == OneOrMorePaths::ONE_PATH)
    {
      if (vrptw.instanceUpperBound < INF)
      {
        double splitUpperBound = vrptw.instanceUpperBound * 1.0 / vrptw.numLocations;
        lambda.push_back(splitUpperBound);
      }
      else
      {
        lambda.push_back(0);
      }
    }
    else
    {
      lambda.push_back(2 * vrptw.distances[0][location] * std::abs(vrptw.demands[location]) / vrptw.capacity);
    }
  }

  // so arc fixing does not try to use IP achieved lb
  double bestLPValue = 0.0;
  bestLambdaArcFixing.resize(lambda.size());
  bestLambda.resize(lambda.size());
  previousLambdaMomentum = lambda;

  double singlePathDual = 0.0;
  std::vector<double> mu;
  std::vector<double> combDuals;
  std::vector<double> srcDuals;

  bool changedLagToLP = false;
  bool finishedSolving = false;
  while (!finishedSolving)
  {
    bool solved = false;
    if (params.lpSolveType == LPSolveType::LPSolver)
    {
      if (stats.lpIterations > 1)
      {
        double percentArcsFixed = 0.0;
        if (!changedLagToLP)
        {
          if (params.useVariableFixing)
          {
            percentArcsFixed = routeDD.fixArcs(lambda, singlePathDual, mu, combDuals, srcDuals, bestLPValue, params.lpSolveType);
          }

          // we don't currently, but we should still save as best
        }
        else
        {
          // TODO(akarahal) if using a counter to avoid loops for LAG solve...
          // ...can now get rid of them by merging nodes!
          // TODO(akarahal) may be helpful to save even more duals for fixing!
          changedLagToLP = false;
        }

        // right now best is only kept for LAG, without cuts
        std::vector<double> emptyMu;
        std::vector<double> emptyComb;
        std::vector<double> emptySrc;
        if (params.useVariableFixing)
        {
          routeDD.fixArcs(bestLambdaArcFixing, bestSinglePathDualFixing, bestMuArcFixing, emptyComb, emptySrc, bestLambdaLowerBound, params.lpSolveType);
        }
      }

      if (lpFlowType == FlowType::LP)
      {
        //routeDD.checkLC121SolutionPossible();
        //routeDD.checkLRC121SolutionPossible();
        solved = solveLP(lambda, singlePathDual, mu, combDuals, srcDuals);
        DBG(double best = 0.0;
        best = best + singlePathDual;
        for (auto l : lambda)
        {
          best = best + l;
        }
        std::cout << "best? " << best << std::endl;
        if ((best < (stats.lowerBound - 0.001)) || (best > (stats.lowerBound + 0.001)))
        {
          std::cout << "OOOOOOOOOOO" << std::endl;
        })
        bestLPValue = stats.lowerBound;
      }
      else
      {
        solved = solveIP(lambda, singlePathDual, mu, combDuals, srcDuals);
      }
    }
    else if (params.lpSolveType == LPSolveType::LAGSolver)
    {
      // maybe restart Lag with the best lambda so far?
      solved = solveLagrangeanRelaxation(lambda, singlePathDual, mu, combDuals, srcDuals);
      // if switching to LP, go right to next iter and use cuts
      if (params.lpSolveType == LPSolveType::LPSolver)
      {
        changedLagToLP = true;
        if ((vrptw.problemType == ProblemType::CVRP) || (vrptw.problemType == ProblemType::TW))
        {
          params.useCuts = true;
        }
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
    if (params.lpSolveType == LPSolveType::LPSolver)
    {
      bool infeasibilityFound = false;
      bool cutAdded = false;

      // do cuts first in case separations mess up dd structure
      if (params.useCuts && (lpFlowType == FlowType::LP))
      {
        // RCC - rounded capacity cuts
        std::vector<int> edgeTail;
        std::vector<int> edgeHead;
        std::vector<double> edgeFlow;
        routeDD.convertSolutionForVRPTWSep(edgeTail, edgeHead, edgeFlow);
        addRCCs(edgeTail, edgeHead, edgeFlow, 100, cutAdded);
      }

      bool stopFindingInfeasibilities = false;
      std::vector<std::vector<int>> infeasibilities;
      std::vector<double> routeFlows;
      std::vector<std::vector<int>> decomposedRoutes;
      while (!stopFindingInfeasibilities)
      {
        std::vector<int> infeasibleRoute;
        routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, params.maxS, DecompositionReason::SEPARATE);
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
        if (params.useCuts && (vrptw.vrptwCapacityType == VRPTWCapacityType::RELAX_CAPACITY) && !allRoutesCapacityFeasible)
        {
          addRCCs(decomposedRoutes, cutAdded);
        }

        for (auto infeasibleRoute : infeasibilities)
        {
          if (routeDD.doesRouteExistByArcs(infeasibleRoute))
          {
            stats.numSeparations = stats.numSeparations + 1;
            routeDD.separateInfeasibleRoute(infeasibleRoute, params.maxS);
            //routeDD.checkLRC121SolutionPossible();
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
        if (params.useCuts && (vrptw.vrptwCapacityType == VRPTWCapacityType::RELAX_CAPACITY) && !allRoutesCapacityFeasible)
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
            if (shouldSolveIP)
            {
              std::cout << "LP solved - switching to IP" << std::endl;
              lpFlowType = FlowType::IP;
            }
            else
            {
              finishedSolving = true;
            }
          }
          else
          {
            // right now, update precision for pdptw only
            if (vrptw.timeStateMultiplier == vrptw.finalTimeStateMultiplier)
	    {
              stats.upperBound = stats.lowerBound;
	    }
	    else if (vrptw.problemType == ProblemType::PDP)
	    {
              std::cout << "updating precision using multiplier" << std::endl;
              vrptw.timeStateMultiplier = vrptw.timeStateMultiplier * 10;
              vrptw.recomputeDistancesPDPTW();
              std::cout << "new timeStateMultiplier: " << vrptw.timeStateMultiplier << std::endl;
              routeDD.updateVRPTW(vrptw);
              routeDD.updateTimeStateMultiplierByTen();
	    }
          }
        }
        else
        {
          finishedSolving = true;
        }
      }
    }

    if (stats.getNumSeconds() >= params.timeoutSeconds)
    {
      finishedSolving = true;
    }

    // allow lag to finish if close enough
    if ((params.lpSolveType == LPSolveType::LAGSolver) && (stats.lowerBound + 0.01 > stats.upperBound))
    {
      finishedSolving = true;
    }

    if (stats.lowerBound == stats.upperBound)
    {
      finishedSolving = true;
    }
  }

  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
  std::cout << "finished solving - complete" << std::endl;
  CMGR_FreeMemCMgr(&MyCutsCMP);
  CMGR_FreeMemCMgr(&MyOldCutsCMP);
  return true;
}

bool VRPTWDDSolver::solveIP(std::vector<double>& lambda, double& singlePathDual, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  routeDD.setCoeffsAsDistances();
  auto startLPTime = std::chrono::high_resolution_clock::now();
  std::cout << "solving IP" << std::endl;
  stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, lambda, singlePathDual, mu, combDuals, srcDuals);
  auto endLPTime = std::chrono::high_resolution_clock::now();
  auto lpSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endLPTime - startLPTime).count();
  stats.millisecondsSolvingLP = stats.millisecondsSolvingLP + lpSolveTime;

  stats.lpIterations = stats.lpIterations + 1;
  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
  return true;
};

bool VRPTWDDSolver::solveLP(std::vector<double>& lambda, double& singlePathDual, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  routeDD.setCoeffsAsDistances();
  auto startLPTime = std::chrono::high_resolution_clock::now();
  bool MIPOn = false;
  if ((routeDD.getPercentFixedArcs() >= 97.5) && MIPOn)
  {
    std::cout << "solving IP" << std::endl;
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, lambda, singlePathDual, mu, combDuals, srcDuals);
  }
  else
  {
    double oldLb = stats.lowerBound;
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, lambda, singlePathDual, mu, combDuals, srcDuals);
    DBG(std::cout << "DUALS " << stats.getNumSeconds() << "," << stats.lowerBound << "," << routeDD.getNumArcsNotRemovedOrReverse() << "," << routeDD.getNumFixedArcs() << ",";
    for (int dualIndex=0; dualIndex<vrptw.numLocations; ++dualIndex)
    {
      std::cout << lambda[dualIndex] << ",";
    }
    std::cout << std::endl;)
    if (oldLb > stats.lowerBound + 0.01)
    {
      std::cout << "ERROR possible - lower bound not monotonically increasing, but this can happen!" << std::endl;
    }
  }
  auto endLPTime = std::chrono::high_resolution_clock::now();
  auto lpSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endLPTime - startLPTime).count();
  stats.millisecondsSolvingLP = stats.millisecondsSolvingLP + lpSolveTime;

  stats.lpIterations = stats.lpIterations + 1;
  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
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

bool VRPTWDDSolver::solveLPCG(std::vector<double>& lambda, double& singlePathDual, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  // column generation - RMP <-> pricing problem
  //initializeColumns();

  auto startLPTime = std::chrono::high_resolution_clock::now();
  bool solved = false;
  double flowObj = stats.lowerBound;
  while (!solved)
  {
    routeDD.setCoeffsAsDistances();
    flowObj = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::USE_CG, lambda, singlePathDual, mu, combDuals, srcDuals);
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
  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
  return true;
};

void VRPTWDDSolver::printMultipliers(std::vector<double>& lambda, std::vector<double>& mu)
{
  std::cout << "lambda: ";
  for (double l : lambda)
  {
    std::cout << l << ",";
  }
  std::cout << std::endl;

  std::cout << "mu: ";
  for (double m : mu)
  {
    std::cout << m << ",";
  }
  std::cout << std::endl;
}

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
    if ((lambda[i] <= 0.000001) && (gamma[i] <= 0))
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
    if (deactivatedCuts.find(i) == deactivatedCuts.end())
    {
      // problem: even with one cut introduced, the RHS - cutValues can be huge
      // problem cont: so the mu goes up a bunch, then back to 0, then up a bunch, etc.
      // idea: cap the gradient size
      double gammaCap = 2.0;
      gamma[gammaIndex] = (-1 * routeDD.getCapCutSetRHS(i)) + cutValues[i];
      std::cout << "rhs: " << routeDD.getCapCutSetRHS(i) << " cuts: " << cutValues[i] << std::endl;
      if (gamma[gammaIndex] > 0)
      {
        gamma[gammaIndex] = std::min(gamma[gammaIndex], gammaCap);
        std::cout << "gamma for index " << i << " capped to " << gammaCap << std::endl;
      }
      else
      {
        gamma[gammaIndex] = std::max(gamma[gammaIndex], -1 * gammaCap);
        std::cout << "gamma for index " << i << " capped to -" << gammaCap << std::endl;
      }

      if ((mu[i] <= 0.1) && (gamma[gammaIndex] <= 0))
      {
        continue;
      }
      else
      {
        normGammaSquared += std::pow(gamma[gammaIndex], 2);
      }
    }
    else
    {
      mu[i] = 0;
      gamma[gammaIndex] = 0;
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
  double alpha = (psiStar - lagrangeanLowerBound) / normGammaSquared;

  // half step size when no progress
  alpha = alpha * alphaFactor;
  stepSizes.push_back(alpha);
  std::cout << "alpha: " << alpha << std::endl;

  DBG(printMultipliers(lambda, mu);)

  // lambda_(k+1) = lambda_(k) + alpha_(k) * gamma_(k)
  for (int i=1; i<vrptw.numLocations; ++i)
  {
    previousLambdaMomentum[i] = lambda[i];
    lambda[i] = std::max(0.0, lambda[i] + alpha * gamma[i]);

    // beta for momentum / heavy ball method
    if (params.cutPhase && (stats.lpIterations > 1))
    {
      lambda[i] = std::max(0.0, lambda[i] + params.momentumBeta * (lambda[i] - previousLambdaMomentum[i]));
    }
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

  // remove if too small for too long
  for (int muIndex=0; muIndex<mu.size(); ++muIndex)
  {
    if (deactivatedCuts.find(muIndex) == deactivatedCuts.end())
    {
      if (mu[muIndex] < params.deactivateCutValueThreshold)
      {
        cutTooSmallCounters[muIndex] = cutTooSmallCounters[muIndex] + 1;
        if (cutTooSmallCounters[muIndex] > params.deactivateCutIterThreshold)
        {
          std::cout << "deactivated cut at index: " << muIndex << std::endl;
          deactivatedCuts.insert(muIndex);
        }
      }
      else
      {
        cutTooSmallCounters[muIndex] = 0;
      }
    }
  }
};

bool VRPTWDDSolver::solveLagrangeanRelaxation(std::vector<double>& lambda, double& singlePathDual, std::vector<double>& mu, std::vector<double>& combDuals, std::vector<double>& srcDuals)
{
  stats.lpIterations = stats.lpIterations + 1;
  bool shouldTerminate = false;
  int kappaIterations = 100;
  double muPercentImproved = 0.001;
  int lastMuImprovedIteration = 0;
  double muLowerBound = 0.0;
  double currIterLowerBound = 0.0;
  int numLagIterations = 0;

  double startingLowerBound = stats.lowerBound;

  // Need to decide when step size has gotten too small and need to 'restart'
  // When there is little progress, this method terminates and gets called again.
  // So, then we can increase the 'psiStar' estimation with a smaller iteration value to increase the 'alpha' step size
  stats.numLagIterationsWithResets = std::ceil(stats.numLagIterationsWithResets / stats.lpIterations);

  std::vector<std::vector<int>> infeasibleRoutes;
  std::vector<double> stepSizes;
  std::vector<std::set<int>> xDecompositions;
  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());

  while (!shouldTerminate)
  {
    while (!shouldTerminate && (infeasibleRoutes.size() < params.infeasibleRoutesBatchSize))
    {
      ++stats.numLagIterations;
      ++numLagIterations;
      ++stats.numLagIterationsWithResets;

      std::vector<std::vector<int>> shortestPaths;
      auto startSSPTime = std::chrono::high_resolution_clock::now();
      DBG(std::cout << "start" << std::endl;)
      bool isDualFeasible = false;
      double minReducedCost = 0.0;
      std::vector<double> repairedLambda(lambda);
      double percentFixed = 0.0;
      if (params.useVariableFixing && (bestLambdaLowerBound > 0.001))
      {
        std::vector<double> emptyComb;
        std::vector<double> emptySrc;
        routeDD.fixArcs(bestLambdaArcFixing, bestSinglePathDualFixing, bestMuArcFixing, emptyComb, emptySrc, bestLambdaLowerBound, params.lpSolveType);
      }

      int notMuSSPSeconds = 0;
      double notMuSSPLowerBound = 0.0;
      int notMuSSPNumPaths = 0;
      if (!params.useMuSSP)
      {
        notMuSSPLowerBound = routeDD.solveMinCostFlowModel(lambda, shortestPaths, isDualFeasible, minReducedCost);
        auto endNotMuSSPTime = std::chrono::high_resolution_clock::now();
        auto notMuSSPTime = std::chrono::duration_cast<std::chrono::milliseconds>(endNotMuSSPTime - startSSPTime).count();
        notMuSSPSeconds = notMuSSPTime / 1000.0;
        minReducedCost = 0.0;
        notMuSSPNumPaths = shortestPaths.size();
        shortestPaths.clear();
        isDualFeasible = false;
      }

      double lagrangeanLowerBound = routeDD.solveMinCostFlowModelWang(lambda, mu, combDuals, srcDuals, shortestPaths, isDualFeasible, minReducedCost);
      //routeDD.checkLC121SolutionPossible();
      //routeDD.checkLRC121SolutionPossible();
      //routeDD.checkC141SolutionPossible();
      DBG(std::cout << "finish" << std::endl;)
      stats.numSSPIterations = stats.numSSPIterations + shortestPaths.size();
      auto endMuSSPTime = std::chrono::high_resolution_clock::now();
      auto sspSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endMuSSPTime - startSSPTime).count();
      stats.millisecondsSolvingSSP = stats.millisecondsSolvingSSP + sspSolveTime;

      if (!params.useMuSSP)
      {
        // ensure they achieve the same value
        if ((notMuSSPLowerBound - lagrangeanLowerBound > 0.00001) || (lagrangeanLowerBound - notMuSSPLowerBound > 0.00001))
        {
          std::cout << "ERROR, muSSP: " << lagrangeanLowerBound << " vs. SSP: " << notMuSSPLowerBound << std::endl;
          std::cout << "ERROR, muSSP num paths: " << shortestPaths.size() << " vs. SSP: " << notMuSSPNumPaths << std::endl;
        }
        else
        {
          std::cout << "COMPARISON size[" << routeDD.getNumArcsNotRemovedOrReverseOrFixed() << "] muSSP[" << (sspSolveTime / 1000.0) - notMuSSPSeconds << "] SSP[" << notMuSSPSeconds << "]" << std::endl;
        }
      }
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

      DBG(std::cout << "DUALS " << stats.getNumSeconds() << "," << stats.lowerBound << "," << routeDD.getNumArcsNotRemovedOrReverse() << "," << routeDD.getNumFixedArcs() << ",";
      for (int dualIndex=0; dualIndex<vrptw.numLocations; ++dualIndex)
      {
        std::cout << lambda[dualIndex] << ",";
      }
      std::cout << std::endl;)

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
        stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
        for (int index=0; index<lambda.size(); ++index)
        {
          bestLambda[index] = lambda[index];
        }
        printMultipliers(lambda, mu);
      }

      // get primal solution
      std::set<int> solutionArcs;
      routeDD.getSolutionArcs(solutionArcs);
      if (params.useCuts)
      {
        xDecompositions.push_back(solutionArcs);
      }

      // one path can only separate one at a time and likely needs many more iterations
      // so reset when progress slows down too much, because likely step size is too small
      // uses solution arcs
      //if (vrptw.oneOrMorePaths == OneOrMorePaths::ONE_PATH)
      //{
      //  updateMultipliers(lambda, mu, combDuals, srcDuals, stepSizes, solutionArcs, lagrangeanLowerBound, numLagIterations);
      //}
      //else
      //{
      updateMultipliers(lambda, mu, combDuals, srcDuals, stepSizes, solutionArcs, lagrangeanLowerBound, stats.numLagIterations);
      //}

      // check for cycles up to certain size and add to be separated
      //if ((params.maxS > 1) && (stats.lpIterations >= 2) && isSeparationRound)
      if (params.maxS > 1)
      {
        bool infeasibleRouteFound = true;
        while (infeasibleRouteFound)
        {
          std::vector<int> infeasibleRoute;
          std::vector<double> flows;
          std::vector<std::vector<int>> routes;
          routeDD.decomposeRoutes(infeasibleRoute, flows, routes, params.maxS, DecompositionReason::SEPARATE);
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
      if ((stats.getNumSeconds() >= params.timeoutSeconds) || (stats.lowerBound + 0.01 > stats.upperBound))
      {
        shouldTerminate = true;
      }

      // use repaired lambda as a copy of lambda for now
      // should fix after decomposing in case we fix an arc that is in the current solution
      // see if we can fix arcs based on a feasible dual
      // some rounds let's force dual feasibility to find the bound and fix some arcs?
      if (isDualFeasible & (vrptw.oneOrMorePaths != OneOrMorePaths::ONE_PATH))
      {
        if (params.useVariableFixing)
        {
          percentFixed = routeDD.fixArcs(repairedLambda, singlePathDual, mu, combDuals, srcDuals, lagrangeanLowerBound, params.lpSolveType);
          //int indexToStore = stats.lagIterations % lambdaStoreSize;
          //lambdaStore[] = repairedLambda;
          //singlePathStore[] = singlePathDual;
          //lambdaLowerBoundStore[] = lagrangeanLowerBound;
        }
        if (percentFixed > bestLambdaPercentFixed)
        {
          std::cout << "updating best lambda" << std::endl;
          bestLambdaPercentFixed = percentFixed;
          bestLambdaLowerBound = lagrangeanLowerBound;
          for (int index=0; index<lambda.size(); ++index)
          {
            bestLambdaArcFixing[index] = repairedLambda[index];
          }
          for (int index=0; index<mu.size(); ++index)
          {
            bestMuArcFixing[index] = mu[index];
          }
        }

        double checkLambdaLB = 0.0;
        for (double dual : repairedLambda)
        {
          checkLambdaLB += dual;
        }
        DBG(std::cout << "check lb sum lambda: " << checkLambdaLB << std::endl;)
      }
      else
      {
        // can repair dual to make it feasible, and then try to fix arcs. Min reduced cost would be negative
        // with one path can repair with the onePathDual
        // repair can be costly so only try when within 2% optimality gap
        if (params.useVariableFixing && ((stats.upperBound - lagrangeanLowerBound) * 100.0 / stats.upperBound < params.lagOptimalityGapToStartRepairing))
        {
          if (vrptw.oneOrMorePaths == OneOrMorePaths::ONE_PATH)
          {
            // "reducedCost" is without the single dual variable, dual feasible means c - uTA >= v
            // minReducedCost gives min c - uTA
            singlePathDual = minReducedCost;
            double repairedBound = 0.0;
            for (auto l : repairedLambda)
            {
              repairedBound += l;
            }
            repairedBound = repairedBound + singlePathDual;

            DBG(std::cout << "single path dual: " << singlePathDual << std::endl;
            std::cout << "repaired lb: " << repairedBound << std::endl;)
            percentFixed = routeDD.fixArcs(repairedLambda, singlePathDual, mu, combDuals, srcDuals, repairedBound, params.lpSolveType);
            if (percentFixed > bestLambdaPercentFixed)
            {
              std::cout << "updating best lambda" << std::endl;
              bestLambdaPercentFixed = percentFixed;
              bestLambdaLowerBound = repairedBound;
              bestSinglePathDualFixing = singlePathDual;
              for (int index=0; index<lambda.size(); ++index)
              {
                bestLambdaArcFixing[index] = repairedLambda[index];
              }
              for (int index=0; index<mu.size(); ++index)
              {
                bestMuArcFixing[index] = mu[index];
              }
            }
          }
          else
          {
            while (true)
            {
              routeDD.setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(repairedLambda, mu, combDuals, srcDuals, LPSolveType::LAGSolver);

              std::vector<int> treeByParentArcs;
              treeByParentArcs.resize(routeDD.getNodes().size());
              std::vector<int> shortestPathByArc;
              double shortestPathLength = routeDD.computeShortestPathBFSWang(treeByParentArcs, shortestPathByArc);
              //double shortestPathLength = routeDD.computeShortestPathBFS(treeByParentArcs, shortestPathByArc);
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
            // - mu_T(-Cx + r) because capCuts are Cx <= r (not >=)
            for (int index=0; index<mu.size(); ++index)
            {
              repairedBound += mu[index] * routeDD.getCapCutSetRHS(index) * -1;
            }

            // - combDual_T(Ax - RHS)
            for (int index=0; index<combDuals.size(); ++index)
            {
              repairedBound += (combDuals[index] * routeDD.getCombCutRHS(index));
            }
            DBG(std::cout << "repaired lb: " << repairedBound << std::endl;)
            // The repaired bound might be the best lb yet!
            if (stats.lowerBound < repairedBound)
            {
              stats.lowerBound = repairedBound;
              stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
              printMultipliers(lambda, mu);
            }

            percentFixed = routeDD.fixArcs(repairedLambda, singlePathDual, mu, combDuals, srcDuals, repairedBound, params.lpSolveType);
            if (percentFixed > bestLambdaPercentFixed)
            {
              std::cout << "updating best lambda" << std::endl;
              bestLambdaPercentFixed = percentFixed;
              bestLambdaLowerBound = repairedBound;
              for (int index=0; index<lambda.size(); ++index)
              {
                bestLambdaArcFixing[index] = repairedLambda[index];
              }
              for (int index=0; index<mu.size(); ++index)
              {
                bestMuArcFixing[index] = mu[index];
              }
            }
          }
        }
      }

      if ((percentFixed > params.percentFixedToChangeToCPLEX) || (routeDD.getNumArcsNotRemovedOrReverseOrFixed() < params.numArcsToChangeToCPLEX))
      {
        params.lpSolveType = LPSolveType::LPSolver;
        stats.lpIterations = 1;
        shouldTerminate = true;
        std::cout << "switching from LAG to LP solver" << std::endl;
      }
    }

    // add separations each round
    //if ((params.maxS > 1) && !shouldTerminate && isSeparationRound)
    if ((params.maxS > 1) && !shouldTerminate)
    {
      // delay lagrangean separations to try to keep DD small
      if (stats.numLagIterations < params.lagIterationDelayToStartSeparating)
      {
        infeasibleRoutes.clear();
      }

      // stop separating after first lag termination
      if (params.cutPhase && (stats.lpIterations > 1))
      {
        infeasibleRoutes.clear();
      }

      for (auto infeasibleRouteToSeparate : infeasibleRoutes)
      {
        if (routeDD.doesRouteExistByArcs(infeasibleRouteToSeparate))
        {
          stats.numSeparations = stats.numSeparations + 1;
          routeDD.separateInfeasibleRoute(infeasibleRouteToSeparate, params.maxS);
          //routeDD.checkLRC121SolutionPossible();
        }
      }
    }

    // add cuts if specified
    if (!shouldTerminate && (params.numLagItersForCuts != 0) && (stats.numLagIterations % params.numLagItersForCuts == 0))
    {
      // for Rounded Capacity Cuts
      bool cutAdded = false;
      std::vector<int> edgeTail;
      std::vector<int> edgeHead;
      std::vector<double> edgeFlow;
      convertArcIndicesForVRPTWSep(stepSizes, xDecompositions, edgeTail, edgeHead, edgeFlow);
      addRCCs(edgeTail, edgeHead, edgeFlow, params.numLagCuts, cutAdded);

      mu.resize(routeDD.getNumCapCuts());
      bestMuArcFixing.resize(mu.size());
      cutTooSmallCounters.resize(mu.size());

      xDecompositions.clear();
      stepSizes.clear();
    }

    if (shouldTerminate && params.useCuts && (!params.cutPhase || (stats.lpIterations > 1)))
    {
      // for Rounded Capacity Cuts
      bool cutAdded = false;
      std::vector<int> edgeTail;
      std::vector<int> edgeHead;
      std::vector<double> edgeFlow;
      convertArcIndicesForVRPTWSep(stepSizes, xDecompositions, edgeTail, edgeHead, edgeFlow);
      addRCCs(edgeTail, edgeHead, edgeFlow, params.numLagCuts, cutAdded);
      mu.resize(routeDD.getNumCapCuts());
      bestMuArcFixing.resize(mu.size());
      cutTooSmallCounters.resize(mu.size());

      // Strengthened Combs
      //addCombs(edgeTail, edgeHead, edgeFlow, cutAdded);
      //combDuals.resize(routeDD.getNumCombCuts());

      // for Subset Row Cuts add right away during cut rounds
      // dont add for X was taking too long
      // TOOD: maybe only add for small instances?
      /*
      routeDD.findSRCThree(xDecompositions, stepSizes, cutAdded);
      addSRCCuts(srcDuals);
      */
      xDecompositions.clear();
      stepSizes.clear();
    }
    infeasibleRoutes.clear();
  }

  if (stats.lowerBound < startingLowerBound + 0.01)
  {
    alphaFactor = alphaFactor / 2;
    std::cout << "halving alphaFactor: " << alphaFactor << std::endl;
  }

  return true;
}
