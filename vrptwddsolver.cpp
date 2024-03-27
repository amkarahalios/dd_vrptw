#include "vrptwddsolver.h"
#include <math.h>

VRPTWDDSolver::VRPTWDDSolver(VRPTW _vrptw, VRPTWDDParameters _params) : vrptw(_vrptw), routeDD(_vrptw, _params), params(_params), phaseType(PhaseType::INITIAL_DUAL)
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

  muPercentImproved = 0.001;

  bestDualArcFixingPercent = 0.0;
  bestDualValue = 0.0;
  stepSizeMultiplier = 1.0;
  stepSizeMultiplierIteration = 0;
  stepSizeMultiplierIterationCutoff = 20;
  alphaLowerBound = 0.1;
  alphaLowerBoundIteration = 0;
  alphaLowerBoundCheckValue = 0;
  targetLowerBound = 1.0;

  CMGR_CreateCMgr(&MyCutsCMP,Dim);
  CMGR_CreateCMgr(&MyOldCutsCMP,Dim);
};

void VRPTWDDSolver::convertArcIndicesForVRPTWSep(const Primal& currPrimal,
                                                std::vector<int>& edgeTail,
                                                std::vector<int>& edgeHead,
                                                std::vector<double>& edgeFlow,
                                                std::vector<int>& rccArcs,
                                                std::vector<double>& rccArcFlows)
{
  std::map<std::pair<int,int>,double> edgeFlows;
  for (int routeIndex=0; routeIndex<currPrimal.xDecompositionArcs.size(); ++routeIndex)
  {
    auto route = currPrimal.xDecompositions[routeIndex];
    auto routeArcs = currPrimal.xDecompositionArcs[routeIndex];
    // Note: Could check and only use feasible if we want
    for (int arcIndex : routeArcs)
    {
      std::pair<int,int> fromAndToIndices = routeDD.getFromAndToLocations(arcIndex);
      int currLoc = (fromAndToIndices.first == 0) ? vrptw.numLocations : fromAndToIndices.first;
      int nextLoc = (fromAndToIndices.second == 0) ? vrptw.numLocations : fromAndToIndices.second;
      if (currLoc < nextLoc)
      {
        edgeFlows[std::make_pair(currLoc,nextLoc)] = edgeFlows[std::make_pair(currLoc,nextLoc)] + currPrimal.xDecompositionFlows[routeIndex];
      }
      else
      {
        edgeFlows[std::make_pair(nextLoc,currLoc)] = edgeFlows[std::make_pair(nextLoc,currLoc)] + currPrimal.xDecompositionFlows[routeIndex];
      }

      rccArcs.push_back(arcIndex);
      rccArcFlows.push_back(currPrimal.xDecompositionFlows[routeIndex]);
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

void VRPTWDDSolver::addRCCs(const std::vector<int>& edgeTail, const std::vector<int>& edgeHead, const std::vector<double>& edgeFlow, std::vector<int>& rccArcs, int maxNumCuts, bool& cutAdded, Dual& dual)
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
      /*
      bool inFamily = false;
      for (int index=0; index<cutSets.size(); ++index)
      {
        if (deactivatedCuts.find(index) == deactivatedCuts.end())
        {
          auto referenceCutSet = cutSets[index];
          if (std::includes(cutSet.begin(), cutSet.end(), referenceCutSet.begin(), referenceCutSet.end()) || std::includes(referenceCutSet.begin(), referenceCutSet.end(), cutSet.begin(), cutSet.end()))
          {
            inFamily = true;
          }
        }
      }

      if (inFamily)
      {
        std::cout << "nested set, do not add" << std::endl;
      }
      */
      //else
      //{
      std::set<int> cutSetAsSet(cutSet.begin(), cutSet.end());
      cutSets.push_back(cutSetAsSet);
      routeDD.addCapCutSet(cutSet, rccArcs, RHS, params.lpSolveType);
      stats.numCuts = stats.numCuts + 1;
      //}

      dual.capDuals.push_back(maxViolation);
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

void VRPTWDDSolver::addSRCCuts(std::vector<double>& srcDuals, const std::vector<double>& violations)
{
  int newViolationIndex = 0;
  for (int index=srcDuals.size(); index<routeDD.getNumCliqueCuts(); ++index)
  {
    srcDuals.push_back(violations[newViolationIndex]);
    ++newViolationIndex;
  }
}

void VRPTWDDSolver::resizeMultipliers(const Dual& dual1, Dual& dual2)
{
  dual2.lambda.resize(dual1.lambda.size());
  dual2.capDuals.resize(dual1.capDuals.size());
  dual2.combDuals.resize(dual1.combDuals.size());
  dual2.srcDuals.resize(dual1.srcDuals.size());
}

void VRPTWDDSolver::initializeDual(Dual& dual)
{
  dual.lambda.resize(vrptw.numLocations);

  dual.lambda[0] = 0;
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    if ((vrptw.fixedNumPaths == FixedNumPaths::FIXED_NUM_PATHS) && (vrptw.numVehicles == 1))
    {
      if (vrptw.instanceUpperBound < INF)
      {
        double splitUpperBound = vrptw.instanceUpperBound * 1.0 / vrptw.numLocations;
        dual.lambda[location] = splitUpperBound;
      }
      else
      {
        dual.lambda[location] = 0;
      }
    }
    else
    {
      dual.lambda[location] = 2 * vrptw.distances[0][location] * std::abs(vrptw.demands[location]) / vrptw.capacity;
    }
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

  Dual dual;
  initializeDual(dual);

  // so arc fixing does not try to use IP achieved lb
  bestDual.lambda.resize(dual.lambda.size());
  bestDualArcFixing.lambda.resize(dual.lambda.size());

  int averageRouteLength = 0;
  int numRoutesInAverage = 0;
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
            repairMultipliers(dual, LPSolveType::LPSolver);
            percentArcsFixed = routeDD.fixArcs(dual, LPSolveType::LPSolver);
          }
        }
        else
        {
          // TODO(akarahal) if using a counter to avoid loops for LAG solve...
          // ...can now get rid of them by merging nodes!
          // TODO(akarahal) may be helpful to save even more duals for fixing!
          changedLagToLP = false;
        }
      }

      // TODO(akarahal) store best duals for arc fixing better
      if (params.useVariableFixing)
      {
        repairMultipliers(bestDualArcFixing, LPSolveType::LAGSolver);
        routeDD.fixArcs(bestDualArcFixing, LPSolveType::LAGSolver);
      }

      if (lpFlowType == FlowType::LP)
      {
        routeDD.strengthenSRCs(averageRouteLength);
        solved = solveLP(dual);
      }
      else
      {
        solved = solveIP(dual);
      }
    }
    else if (params.lpSolveType == LPSolveType::LAGSolver)
    {
      // maybe restart Lag with the best lambda so far?
      //solved = solveLagrangeanRelaxation(dual);
      solved = solveLagrangeanRelaxationVolumeAlgorithm(dual);

      // if switching to LP, go right to next iter and use cuts
      if (params.lpSolveType == LPSolveType::LPSolver)
      {
        changedLagToLP = true;
        params.useCuts = true;
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
      bool cutAdded = false;
      int numSrcAdded = 0;
      std::vector<std::vector<int>> decomposedRoutes;
      std::vector<std::vector<int>> decomposedArcs;
      std::vector<double> routeFlows;

      // do cuts first in case separations mess up dd structure
      if (params.useCuts && (lpFlowType == FlowType::LP))
      {
        // RCC - rounded capacity cuts
        std::vector<int> edgeTail;
        std::vector<int> edgeHead;
        std::vector<double> edgeFlow;
        std::vector<int> rccArcs;
        std::vector<double> rccArcFlows;
        routeDD.convertSolutionForVRPTWSep(edgeTail, edgeHead, edgeFlow, rccArcs, rccArcFlows);
        addRCCs(edgeTail, edgeHead, edgeFlow, rccArcs, 100, cutAdded, dual);
        dual.capDuals.resize(routeDD.getNumCapCuts());

        // Subset Row Cuts
        std::vector<int> infeasibleRoute;
        routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, params.maxS, DecompositionReason::DECOMPOSE);
        averageRouteLength = averageRouteLength * numRoutesInAverage;
        for (auto route : decomposedRoutes)
        {
          averageRouteLength = averageRouteLength + route.size();
        }
        numRoutesInAverage += decomposedRoutes.size();
        averageRouteLength = averageRouteLength / numRoutesInAverage;

        // currently adding by full separation, can strengthen with all up / all down
        //routeDD.print();
        std::vector<double> violations;
        numSrcAdded = routeDD.findSRC3s(primal, 10, violations);
        numSrcAdded += routeDD.findSRC4s(primal, 10, violations);
        numSrcAdded += routeDD.findSRC5V1s(primal, 10, violations);
        numSrcAdded += routeDD.findSRC5V2s(primal, 10, violations);
        if (numSrcAdded > 0)
        {
          cutAdded = true;
          stats.numCuts = stats.numCuts + numSrcAdded;
        }

        addSRCCuts(dual.srcDuals, violations);

        resizeMultipliers(dual, bestDual);
        resizeMultipliers(dual, bestDualArcFixing);
      }

      bool stopFindingInfeasibilities = false;
      std::vector<std::vector<int>> infeasibilities;
      if (numSrcAdded == 0)
      {
        while (!stopFindingInfeasibilities)
        {
          std::vector<int> infeasibleRoute;
          routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, params.maxS, DecompositionReason::SEPARATE);
          if (!infeasibleRoute.empty())
          {
            infeasibilities.push_back(infeasibleRoute);
          }
          else
          {
            stopFindingInfeasibilities = true;
          }
        }
      }
      else
      {
        // adding cuts could separate graph / ruin flow structure
        // in this case, check the decomposed paths directly
        for (int index=0; index<decomposedRoutes.size(); ++index)
        {
          auto route = decomposedRoutes[index];
          if (!routeDD.isRouteFeasible(route))
          {
            auto routeArcs = decomposedArcs[index];
            infeasibilities.push_back(routeArcs);
          }
        }
      }

      if (infeasibilities.size() > 0)
      {
        for (auto infeasibleRoute : infeasibilities)
        {
          if (routeDD.doesRouteExistByArcs(infeasibleRoute))
          {
            stats.numSeparations = stats.numSeparations + 1;
            routeDD.separateInfeasibleRoute(infeasibleRoute, params.maxS);
          }
        }
      }
      else
      {
        // print decomposed routes and check capacity feasibility
        std::cout << "routes: " << std::endl;
        int index = 0;
        for (auto route : decomposedRoutes)
        {
          int load = 0;
          std::cout << "flow {" << routeFlows[index] << "}: ";
          for (int index=0; index<route.size()-1; ++index)
          {
            std::cout << route[index] << " ";
          }
          std::cout << std::endl;
          index = index + 1;
        }
        std::cout << "no more separations possible" << std::endl;
      }

      if (!cutAdded && (infeasibilities.size() == 0))
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

bool VRPTWDDSolver::solveIP(Dual& duals)
{
  routeDD.setCoeffsAsDistances();
  auto startLPTime = std::chrono::high_resolution_clock::now();
  std::cout << "solving IP" << std::endl;
  stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, duals);
  auto endLPTime = std::chrono::high_resolution_clock::now();
  auto lpSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endLPTime - startLPTime).count();
  stats.millisecondsSolvingLP = stats.millisecondsSolvingLP + lpSolveTime;

  stats.lpIterations = stats.lpIterations + 1;
  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
  return true;
};

bool VRPTWDDSolver::solveLP(Dual& duals)
{
  routeDD.setCoeffsAsDistances();
  auto startLPTime = std::chrono::high_resolution_clock::now();
  bool MIPOn = false;
  if ((routeDD.getPercentFixedArcs() >= 97.5) && MIPOn)
  {
    std::cout << "solving IP" << std::endl;
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, duals);
  }
  else
  {
    double oldLb = stats.lowerBound;
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, duals);
    if (oldLb > stats.lowerBound + 0.01)
    {
      std::cout << "WARNING possible - lower bound not monotonically increasing, but this can happen!" << std::endl;
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

bool VRPTWDDSolver::solveLPCG(Dual& duals)
{
  // column generation - RMP <-> pricing problem
  //initializeColumns();

  auto startLPTime = std::chrono::high_resolution_clock::now();
  bool solved = false;
  double flowObj = stats.lowerBound;
  while (!solved)
  {
    routeDD.setCoeffsAsDistances();
    flowObj = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::USE_CG, duals);
    std::cout << "flowobj: " << flowObj << std::endl;

    bool addedColumn = solvePricingProblem(duals.lambda);
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

void VRPTWDDSolver::repairMultipliers(Dual& repairedDual, LPSolveType solveType)
{
  routeDD.clearRelaxedSrcs();

  while (true)
  {
    routeDD.setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(repairedDual, solveType);

    std::vector<int> treeByParentArcs;
    treeByParentArcs.resize(routeDD.getNodes().size());
    std::vector<int> shortestPathByArc;
    double shortestPathLength = routeDD.computeShortestPathBFSWang(treeByParentArcs, shortestPathByArc);
    shortestPathLength = shortestPathLength - repairedDual.fixedPathDual;
    if (shortestPathLength >= 0)
    {
      break;
    }
    else
    {
      bool updated = false;
      std::vector<std::vector<int>> shortestPathsByArc;
      shortestPathsByArc.push_back(shortestPathByArc);
      std::set<int> locations;
      routeDD.getLocationsOnArcPaths(shortestPathsByArc, locations);
      double updateAmount = std::min(shortestPathLength, -0.001);
      if (shortestPathByArc.size() - 1 > 0)
      {
        for (int loc : locations)
        {
          if (repairedDual.lambda[loc] > 0.00001)
          {
            updated = true;
            repairedDual.lambda[loc] = std::max(0.0, repairedDual.lambda[loc] + (updateAmount / (shortestPathByArc.size() - 1)) - 0.000001);
          }
        }
      }

      if (!updated)
      {
        if (vrptw.fixedNumPaths == FixedNumPaths::FIXED_NUM_PATHS)
        {
          repairedDual.fixedPathDual = repairedDual.fixedPathDual + shortestPathLength - 0.001;
        }

        for (int index=1; index<vrptw.numLocations; ++index)
        {
          repairedDual.lambda[index] = std::max(0.0, repairedDual.lambda[index] - 0.01);
        }
      }
    }
  }
};

void VRPTWDDSolver::printMultipliers(Dual& dual)
{
  std::cout << "lambda: ";
  for (double l : dual.lambda)
  {
    std::cout << l << ",";
  }
  std::cout << std::endl;

  std::cout << "mu: ";
  for (int index=0; index<dual.capDuals.size(); ++index)
  {
    if (routeDD.isCapCutActive(index))
    {
      std::cout << dual.capDuals[index] << ",";
    }
  }
  std::cout << std::endl;

  std::cout << "srcDuals: ";
  for (int index=0; index<dual.srcDuals.size(); ++index)
  {
    if (routeDD.isCliqueCutActive(index))
    {
      std::cout << dual.srcDuals[index] << ",";
    }
  }
  std::cout << std::endl;
}

void VRPTWDDSolver::updateMultipliers(Dual& dual, double lagrangeanLowerBound, int iteration)
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
  std::vector<double> gamma(vrptw.numLocations + dual.capDuals.size() + dual.srcDuals.size() + dual.combDuals.size(), 0);
  for (int i=1; i<vrptw.numLocations; ++i)
  {
    gamma[i] = 1 - locationsCovered[i];
    if ((dual.lambda[i] <= 0.000001) && (gamma[i] <= 0))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[i], 2);
    }
  }

  // cap cuts are in the form Cx <= r so change to -Cx >= -r
  for (int i=0; i<dual.capDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations;
    if (routeDD.isCapCutActive(i))
    {
      gamma[gammaIndex] = (-1 * routeDD.getCapCutSetRHS(i)) + cutValues[i];

      // problem: even with one cut introduced, the RHS - cutValues can be huge
      // problem cont: so the mu goes up a bunch, then back to 0, then up a bunch, etc.
      // idea: cap the gradient size
      /*
      double gammaCap = 1.0;
      std::cout << "rhs: " << routeDD.getCapCutSetRHS(i) << " cuts: " << cutValues[i] << std::endl;
      if ((gamma[gammaIndex] > 0) && (gamma[gammaIndex] > gammaCap))
      {
        gamma[gammaIndex] = std::min(gamma[gammaIndex], gammaCap);
        std::cout << "gamma for mu index " << i << " capped to " << gamma[gammaIndex] << std::endl;
      }
      else if ((gamma[gammaIndex] < 0) && (gamma[gammaIndex] < -1 * gammaCap))
      {
        gamma[gammaIndex] = std::max(gamma[gammaIndex], -1 * gammaCap);
        std::cout << "gamma for mu index " << i << " capped to " << gamma[gammaIndex] << std::endl;
      }
      */

      if ((dual.capDuals[i] <= 0.1) && (gamma[gammaIndex] <= 0))
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
      dual.capDuals[i] = 0;
      gamma[gammaIndex] = 0;
    }
  }

  // src cuts are in the form Cx <= r so chang eto -Cx >= -r
  for (int i=0; i<dual.srcDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations + dual.capDuals.size();
    SRCType srcType = routeDD.getSRCType(i);
    int rhs = routeDD.getSRCRHS(srcType);
    gamma[gammaIndex] = (-1*rhs) + cliqueCutValues[i];

    // problem: even with one cut introduced, the RHS - cutValues can be huge
    // problem cont: so the mu goes up a bunch, then back to 0, then up a bunch, etc.
    // idea: cap the gradient size
    /*
    double gammaCap = 1.0;
    if ((gamma[gammaIndex] > 0) && (gamma[gammaIndex] > gammaCap))
    {
      gamma[gammaIndex] = std::min(gamma[gammaIndex], gammaCap);
      std::cout << "gamma for src index " << i << " capped to " << gamma[gammaIndex] << std::endl;
    }
    else if ((gamma[gammaIndex] < 0) && (gamma[gammaIndex] < gammaCap))
    {
      gamma[gammaIndex] = std::max(gamma[gammaIndex], -1 * gammaCap);
      std::cout << "gamma for src index " << i << " capped to " << gamma[gammaIndex] << std::endl;
    }
    */

    if ((dual.srcDuals[i] <= 0.01) && (gamma[gammaIndex] <= 0.01))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
  }

  for (int i=0; i<dual.combDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations + dual.capDuals.size() + dual.srcDuals.size();
    gamma[gammaIndex] = routeDD.getCombCutRHS(i) - combValues[i];
    if ((dual.combDuals[i] <= 0.01) && (gamma[gammaIndex] <= 0.01))
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
  double eta = 0.05 * 100 / (100 + iteration);

  // psi_(star) = psi_(best) * (1 + eta_(k))
  double psiStar = stats.lowerBound * (1 + eta);

  // alpha_(k) = (psi_(star) - psi(lambda(k))) / ||gamma_(k)||_(2)^2
  double alpha = (psiStar - lagrangeanLowerBound) / normGammaSquared;

  std::cout << "alpha: " << alpha << std::endl;
  std::cout << "psiStar: " << psiStar << std::endl;
  std::cout << "lagLB: " << lagrangeanLowerBound << std::endl;
  std::cout << "||gamma||^2: " << normGammaSquared << std::endl;

  // lambda_(k+1) = lambda_(k) + alpha_(k) * gamma_(k)
  for (int i=1; i<vrptw.numLocations; ++i)
  {
    dual.lambda[i] = std::max(0.0, dual.lambda[i] + alpha * gamma[i]);
  }

  // same for mu
  for (int i=0; i<dual.capDuals.size(); ++i)
  {
    dual.capDuals[i] = std::max(0.0, dual.capDuals[i] + alpha * gamma[i+vrptw.numLocations]);
  }
 
  // same for srcDuals
  for (int i=0; i<dual.srcDuals.size(); ++i)
  {
    dual.srcDuals[i] = std::max(0.0, dual.srcDuals[i] + alpha * gamma[i+vrptw.numLocations+dual.capDuals.size()]);
  }

  // same for combs
  for (int i=0; i<dual.combDuals.size(); ++i)
  {
    dual.combDuals[i] = std::max(0.0, dual.combDuals[i] + alpha * gamma[i+vrptw.numLocations+dual.capDuals.size()+dual.srcDuals.size()]);
  }

  // remove if too small for too long
  for (int muIndex=0; muIndex<dual.capDuals.size(); ++muIndex)
  {
    if (routeDD.isCapCutActive(muIndex))
    {
      if (dual.capDuals[muIndex] < params.deactivateCutValueThreshold)
      {
        capCutTooSmallCounters[muIndex] = capCutTooSmallCounters[muIndex] + 1;
        if (capCutTooSmallCounters[muIndex] > params.deactivateCutIterThreshold)
        {
          std::cout << "deactivated cap cut at index: " << muIndex << std::endl;
          routeDD.deactivateCapCut(muIndex);
          stats.numCuts = stats.numCuts - 1;
        }
      }
      else
      {
        capCutTooSmallCounters[muIndex] = 0;
      }
    }
  }
 
  // remove if too small for too long
  for (int srcIndex=0; srcIndex<dual.srcDuals.size(); ++srcIndex)
  {
    if (routeDD.isCliqueCutActive(srcIndex))
    {
      if (dual.srcDuals[srcIndex] < params.deactivateCutValueThreshold)
      {
        cliqueCutTooSmallCounters[srcIndex] = cliqueCutTooSmallCounters[srcIndex] + 1;
        if (cliqueCutTooSmallCounters[srcIndex] > params.deactivateCutIterThreshold)
        {
          std::cout << "deactivated src cut at index: " << srcIndex << std::endl;
          routeDD.deactivateCliqueCut(srcIndex);
          stats.numCuts = stats.numCuts - 1;
        }
      }
      else
      {
        cliqueCutTooSmallCounters[srcIndex] = 0;
      }
    }
  }
};

// Subgradient Descent
bool VRPTWDDSolver::solveLagrangeanRelaxation(Dual& dual)
{
  bool shouldTerminate = false;
  int kappaIterations = 100;
  int lastMuImprovedIteration = 0;
  double muLowerBound = 0.0;
  double currIterLowerBound = 0.0;
  int numLagIterations = 0;

  primal = Primal();
  std::vector<std::vector<int>> infeasibleRoutes;

  stats.lpIterations = stats.lpIterations + 1;
  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
  while (!shouldTerminate)
  {
    while (!shouldTerminate && (infeasibleRoutes.size() < params.infeasibleRoutesBatchSize))
    {
      ++stats.numLagIterations;
      ++numLagIterations;
      ++stats.numLagIterationsWithResets;

      // Arc fixing
      bool isDualFeasible = false;
      double minReducedCost = 0.0;
      Dual repairedDual(repairedDual);
      double percentFixed = 0.0;
      if (params.useVariableFixing && (bestDualArcFixingPercent > 0.001))
      {
        repairMultipliers(bestDualArcFixing, LPSolveType::LAGSolver);
        routeDD.fixArcs(bestDualArcFixing, LPSolveType::LAGSolver);
      }

      // Run muSSP (or SSP to test differences in time)
      std::vector<std::vector<int>> shortestPaths;
      auto startSSPTime = std::chrono::high_resolution_clock::now();
      int notMuSSPSeconds = 0;
      double notMuSSPLowerBound = 0.0;
      int notMuSSPNumPaths = 0;
      if (!params.useMuSSP)
      {
        notMuSSPLowerBound = routeDD.solveMinCostFlowModel(dual.lambda, shortestPaths, isDualFeasible, minReducedCost);
        auto endNotMuSSPTime = std::chrono::high_resolution_clock::now();
        auto notMuSSPTime = std::chrono::duration_cast<std::chrono::milliseconds>(endNotMuSSPTime - startSSPTime).count();
        notMuSSPSeconds = notMuSSPTime / 1000.0;
        minReducedCost = 0.0;
        notMuSSPNumPaths = shortestPaths.size();
        shortestPaths.clear();
        isDualFeasible = false;
      }

      double lagrangeanLowerBound = routeDD.solveMinCostFlowModelWang(dual, shortestPaths, isDualFeasible, minReducedCost);
      stats.numSSPIterations = stats.numSSPIterations + shortestPaths.size();
      auto endMuSSPTime = std::chrono::high_resolution_clock::now();
      auto sspSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endMuSSPTime - startSSPTime).count();
      stats.millisecondsSolvingSSP = stats.millisecondsSolvingSSP + sspSolveTime;

      // Test that SSP and muSSP get the same values
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

      // Impute the dual bound value and value of fixedPathDual
      // dTx - lambda_T(Ax - b), so add sum of lambdas
      double dualBoundWithoutFixedPathDual = 0.0;
      for (int index=0; index<dual.lambda.size(); ++index)
      {
        lagrangeanLowerBound += dual.lambda[index];
        dualBoundWithoutFixedPathDual += dual.lambda[index];
      }

      // - mu_T(-Cx + r) because capCuts are Cx <= r (not >=)
      for (int index=0; index<dual.capDuals.size(); ++index)
      {
        lagrangeanLowerBound += dual.capDuals[index] * routeDD.getCapCutSetRHS(index) * -1;
        dualBoundWithoutFixedPathDual += dual.capDuals[index] * routeDD.getCapCutSetRHS(index) * -1;
      }
 
      // - combDual_T(Ax - RHS)
      for (int index=0; index<dual.combDuals.size(); ++index)
      {
        lagrangeanLowerBound += (dual.combDuals[index] * routeDD.getCombCutRHS(index));
        dualBoundWithoutFixedPathDual += (dual.combDuals[index] * routeDD.getCombCutRHS(index));
      }

      // - src_Duals_T(-Cx + r) because srcDuals are Cx <= r (not >=)
      for (int index=0; index<dual.srcDuals.size(); ++index)
      {
        lagrangeanLowerBound += dual.srcDuals[index] * -1;
        dualBoundWithoutFixedPathDual += dual.srcDuals[index] * -1;
      }

      if (vrptw.fixedNumPaths == FixedNumPaths::FIXED_NUM_PATHS)
      {
        dual.fixedPathDual = (lagrangeanLowerBound - dualBoundWithoutFixedPathDual) / vrptw.numVehicles;
        std::cout << "fixed path dual: " << dual.fixedPathDual << std::endl;
        isDualFeasible = false;
      }
      else
      {
        dual.fixedPathDual = 0.0;
      }

      // Keep track of best-known so far
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
        printMultipliers(dual);
      }

      // Compute primal solution
      std::vector<int> infeasibleRoute;
      std::vector<double> routeFlows;
      std::vector<std::vector<int>> decomposedRoutes;
      std::vector<std::vector<int>> decomposedArcs;
      routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, params.maxS, DecompositionReason::DECOMPOSE);
      if (params.useCuts)
      {
        // use weighting 95% previous, 5% current
        double alpha = 0.05;
        int firstNonZeroIndex = INF;
        for (int flowIndex=0; flowIndex<primal.xDecompositionFlows.size(); ++flowIndex)
        {
          double newFlow = (1 - alpha) * primal.xDecompositionFlows[flowIndex];
          if (newFlow > 0.0001)
          {
            primal.xDecompositionFlows[flowIndex] = (1 - alpha) * primal.xDecompositionFlows[flowIndex];
            firstNonZeroIndex = std::min(firstNonZeroIndex, flowIndex);
          }
          else
          {
            primal.xDecompositionFlows[flowIndex] = 0;
          }
        }

        // erase any 0 flow if too small
        if ((firstNonZeroIndex > 0) && (firstNonZeroIndex < INF))
        {
          primal.xDecompositions.erase(primal.xDecompositions.begin(), primal.xDecompositions.begin() + firstNonZeroIndex);
          primal.xDecompositionArcs.erase(primal.xDecompositionArcs.begin(), primal.xDecompositionArcs.begin() + firstNonZeroIndex);
          primal.xDecompositionFlows.erase(primal.xDecompositionFlows.begin(), primal.xDecompositionFlows.begin() + firstNonZeroIndex);
        }

        for (int index=0; index<decomposedRoutes.size(); ++index)
        {
          auto route = decomposedRoutes[index];
          auto routeArcs = decomposedArcs[index];
          primal.xDecompositions.push_back(route);
          primal.xDecompositionArcs.push_back(routeArcs);
          primal.xDecompositionFlows.push_back(alpha);
        }
      }

      updateMultipliers(dual, lagrangeanLowerBound, stats.numLagIterations);

      // Check for infeasibilities
      if (params.maxS > 1)
      {
        bool infeasibleRouteFound = true;
        while (infeasibleRouteFound)
        {
          std::vector<int> infeasibleRoute;
          std::vector<double> flows;
          std::vector<std::vector<int>> routes;
          std::vector<std::vector<int>> decomposedArcs;
          routeDD.decomposeRoutes(infeasibleRoute, flows, routes, decomposedArcs, params.maxS, DecompositionReason::SEPARATE);
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

      // Check termination criteria
      if ((numLagIterations - lastMuImprovedIteration) > kappaIterations)
      {
        shouldTerminate = true;
      }

      if ((stats.getNumSeconds() >= params.timeoutSeconds) || (stats.lowerBound + 0.01 > stats.upperBound))
      {
        shouldTerminate = true;
      }

      // Repair dual to be feasible if necessary and maintain best dual for arc fixing
      // use repaired lambda as a copy of lambda for now
      // should fix after decomposing in case we fix an arc that is in the current solution
      if (isDualFeasible & (vrptw.fixedNumPaths != FixedNumPaths::FIXED_NUM_PATHS))
      {
        if (params.useVariableFixing)
        {
          repairMultipliers(repairedDual, LPSolveType::LAGSolver);
          percentFixed = routeDD.fixArcs(repairedDual, params.lpSolveType);
          std::cout << "percent fixed: " << percentFixed << std::endl;
        }
        if (percentFixed > bestDualArcFixingPercent)
        {
          std::cout << "updating best lambda arc fixing" << std::endl;
          bestDualArcFixing = repairedDual;
          bestDualArcFixingPercent = percentFixed;
        }
      }
      else
      {
        // Repair can be costly so only try when within a percent of the optimality gap
        if (params.useVariableFixing && ((stats.upperBound - lagrangeanLowerBound) * 100.0 / stats.upperBound < params.lagOptimalityGapToStartRepairing))
        {
          repairMultipliers(repairedDual, LPSolveType::LAGSolver);

          double repairedBound = routeDD.getDualObjectiveValue(repairedDual, LPSolveType::LAGSolver);
          if (stats.lowerBound < repairedBound)
          {
            stats.lowerBound = repairedBound;
            stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
            printMultipliers(dual);
          }

          percentFixed = routeDD.fixArcs(repairedDual, params.lpSolveType);
          std::cout << "percent fixed: " << percentFixed << std::endl;
          if (percentFixed > bestDualArcFixingPercent)
          {
            std::cout << "updating best lambda arc fixing" << std::endl;
            bestDualArcFixing = repairedDual;
            bestDualArcFixingPercent = percentFixed;
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

    // Separations
    if ((params.maxS > 1) && !shouldTerminate)
    {
      // delay lagrangean separations to try to keep DD small
      if (stats.numLagIterations < params.lagIterationDelayToStartSeparating)
      {
        infeasibleRoutes.clear();
      }

      for (auto infeasibleRouteToSeparate : infeasibleRoutes)
      {
        // when separating changes, can use locations to re-find the route
        if (routeDD.doesRouteExistByArcs(infeasibleRouteToSeparate))
        {
          stats.numSeparations = stats.numSeparations + 1;
          routeDD.separateInfeasibleRoute(infeasibleRouteToSeparate, params.maxS);
        }
      }
    }

    // Cuts
    if (!shouldTerminate && (params.numLagItersForCuts != 0) && (stats.numLagIterations % params.numLagItersForCuts == 0))
    {
      addCutsUsingCurrentPrimal(dual);
    }

    if (shouldTerminate && params.useCuts)
    {
      addCutsUsingCurrentPrimal(dual);
    }
    infeasibleRoutes.clear();
  }

  return true;
}

// Barahona Volume Algorithm
// Always move from pi_{bar}, the best dual solution found so far
// Use subgradient 1 - Ax_{bar}
// Use red, yellow, green iterations that update step size dynamically
bool VRPTWDDSolver::solveLagrangeanRelaxationVolumeAlgorithm(Dual& dual)
{
  bool shouldTerminate = false;
  int kappaIterations = 100;
  int lastMuImprovedIteration = 0;
  double muLowerBound = 0.0;
  double currIterLowerBound = 0.0;
  int numLagIterations = 0;

  std::vector<std::vector<int>> infeasibleRoutesByArc;
  stats.lpIterations = stats.lpIterations + 1;
  std::vector<std::vector<int>> infeasibleRoutes;

  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
  while (!shouldTerminate)
  {
    while (!shouldTerminate && (infeasibleRoutesByArc.size() < params.infeasibleRoutesBatchSize))
    {
      ++stats.numLagIterations;
      ++numLagIterations;
      ++stats.numLagIterationsWithResets;

      // Arc fixing
      bool isDualFeasible = false;
      double minReducedCost = 0.0;
      Dual repairedDual(dual);
      double percentFixed = 0.0;
      if (params.useVariableFixing && (bestDualArcFixingPercent > 0.001))
      {
        repairMultipliers(bestDualArcFixing, LPSolveType::LAGSolver);
        routeDD.fixArcs(bestDualArcFixing, LPSolveType::LAGSolver);
      }

      // Run muSSP
      std::vector<std::vector<int>> shortestPaths;
      auto startSSPTime = std::chrono::high_resolution_clock::now();
      double lagrangeanLowerBound = routeDD.solveMinCostFlowModelWang(dual, shortestPaths, isDualFeasible, minReducedCost);
      stats.numSSPIterations = stats.numSSPIterations + shortestPaths.size();
      auto endMuSSPTime = std::chrono::high_resolution_clock::now();
      auto sspSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endMuSSPTime - startSSPTime).count();
      stats.millisecondsSolvingSSP = stats.millisecondsSolvingSSP + sspSolveTime;

      // Impute the dual bound value and value of fixedPathDual
      double dualBoundWithoutFixedPathDual = 0.0;
      for (int index=0; index<dual.lambda.size(); ++index)
      {
        lagrangeanLowerBound += dual.lambda[index];
        dualBoundWithoutFixedPathDual += dual.lambda[index];
      }

      for (int index=0; index<dual.capDuals.size(); ++index)
      {
        lagrangeanLowerBound += dual.capDuals[index] * routeDD.getCapCutSetRHS(index) * -1;
        dualBoundWithoutFixedPathDual += dual.capDuals[index] * routeDD.getCapCutSetRHS(index) * -1;
      }
 
      for (int index=0; index<dual.combDuals.size(); ++index)
      {
        lagrangeanLowerBound += (dual.combDuals[index] * routeDD.getCombCutRHS(index));
        dualBoundWithoutFixedPathDual += (dual.combDuals[index] * routeDD.getCombCutRHS(index));
      }

      for (int index=0; index<dual.srcDuals.size(); ++index)
      {
        lagrangeanLowerBound += dual.srcDuals[index] * -1;
        dualBoundWithoutFixedPathDual += dual.srcDuals[index] * -1;
      }

      if (vrptw.fixedNumPaths == FixedNumPaths::FIXED_NUM_PATHS)
      {
        dual.fixedPathDual = (lagrangeanLowerBound - dualBoundWithoutFixedPathDual) / vrptw.numVehicles;
        std::cout << "fixed path dual: " << dual.fixedPathDual << std::endl;
        isDualFeasible = false;
      }
      else
      {
        dual.fixedPathDual = 0.0;
      }

      // Keep track of best-known solution
      if (currIterLowerBound < lagrangeanLowerBound)
      {
        if (((1 + muPercentImproved) * muLowerBound) < lagrangeanLowerBound)
        {
          lastMuImprovedIteration = numLagIterations;
          muLowerBound = lagrangeanLowerBound;
        }
        currIterLowerBound = lagrangeanLowerBound;
      }

      // Update step sizes
      bool isImproved = false;
      if (bestDualValue < lagrangeanLowerBound)
      {
        isImproved = true;
        stats.lowerBound = std::max(stats.lowerBound, lagrangeanLowerBound);
        stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
        printMultipliers(dual);

        bestDual = dual;
        bestDualValue = lagrangeanLowerBound;

        stepSizeMultiplierIteration = 0;
        if (lagrangeanLowerBound > 0.95 * targetLowerBound)
        {
          targetLowerBound = std::min(lagrangeanLowerBound * 1.05, stats.upperBound);
        }

        alphaLowerBoundIteration = 0;
        alphaLowerBoundCheckValue = bestDualValue;
      }
      else
      {
        if (stepSizeMultiplierIteration == stepSizeMultiplierIterationCutoff)
        {
          stepSizeMultiplier = std::max(0.0005, stepSizeMultiplier * 0.66);
          stepSizeMultiplierIteration = 0;
        }
        else
        {
          ++stepSizeMultiplierIteration;
        }
 
        ++alphaLowerBoundIteration;
        if (alphaLowerBoundIteration == 100)
        {
          alphaLowerBoundIteration = 0;
          alphaLowerBound = alphaLowerBound / 2.0;
        }
      }

      // Compute primal solution
      std::vector<int> infeasibleRoute;
      std::vector<double> routeFlows;
      std::vector<std::vector<int>> decomposedRoutes;
      std::vector<std::vector<int>> decomposedArcs;
      routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, params.maxS, DecompositionReason::DECOMPOSE);

      // for yellow, check v^{t} dot (1 - Ax^{t})
      if (isImproved)
      {
        Primal latestPrimal(primal);
        constructNextPrimal(1.0, decomposedRoutes, decomposedArcs, latestPrimal);

        std::vector<double> latestGradient;
        getGradient(latestPrimal, dual, latestGradient);

        previousGradient.resize(latestGradient.size());
        double dotProduct = calculateDotProduct(previousGradient, latestGradient);
        if (dotProduct > 0)
        {
          stepSizeMultiplier = std::min(stepSizeMultiplier * 1.1, 2.0);
        }
      }

      // tune the weighting of the primals
      // min ||b-A(alpha x^{t} + (1-alpha)x^{bar}|| s.t. u/10 <= alpha <= u
      // currently using weighting 95% previous, 5% current
      double alphaTry = alphaLowerBound / 10;
      double bestAlpha = alphaTry;
      double bestAlphaValue = INF;
      if (!primal.xDecompositionFlows.empty())
      {
        while (alphaTry < alphaLowerBound)
        {
          Primal nextPrimal(primal);
          constructNextPrimal(alphaTry, decomposedRoutes, decomposedArcs, nextPrimal);

          std::vector<double> gradient;
          getGradient(nextPrimal, dual, gradient);

          double alphaValue = calculateTwoNorm(gradient);
          if (alphaValue < bestAlphaValue)
          {
            bestAlphaValue = alphaValue;
            bestAlpha = alphaTry;
            previousGradient = gradient;
          }
          alphaTry = std::min(alphaTry + (alphaLowerBound / 10), alphaLowerBound);
        }
      }
      else
      {
        bestAlpha = 1.0;
      }
      std::cout << "best alpha: " << bestAlpha << std::endl;

      Primal nextPrimal(primal);
      constructNextPrimal(bestAlpha, decomposedRoutes, decomposedArcs, nextPrimal);
      primal = nextPrimal;

      // Update dual
      updateMultipliersVolumeAlgorithm(dual, primal, lagrangeanLowerBound, stats.numLagIterations);

      // Check for infeasibilities
      if (phaseType == PhaseType::SEPARATION)
      {
        bool infeasibleRouteFound = true;
        while (infeasibleRouteFound)
        {
          std::vector<int> infeasibleRouteByArcs;
          std::vector<double> flows;
          std::vector<std::vector<int>> routes;
          std::vector<std::vector<int>> separationDecomposedArcs;
          routeDD.decomposeRoutes(infeasibleRouteByArcs, flows, routes, separationDecomposedArcs, params.maxS, DecompositionReason::SEPARATE);
          if (!infeasibleRouteByArcs.empty())
          {
            infeasibleRoutesByArc.push_back(infeasibleRouteByArcs);
          }
          else
          {
            infeasibleRouteFound = false;
          }
        }
      }

      // Check termination criteria
      if ((numLagIterations - lastMuImprovedIteration) > kappaIterations)
      {
        shouldTerminate = true;
      }

      if ((stats.getNumSeconds() >= params.timeoutSeconds) || (stats.lowerBound + 0.01 > stats.upperBound))
      {
        shouldTerminate = true;
      }

      // Repair Dual
      if (isDualFeasible & (vrptw.fixedNumPaths != FixedNumPaths::FIXED_NUM_PATHS))
      {
        if (params.useVariableFixing)
        {
          repairMultipliers(repairedDual, LPSolveType::LAGSolver);
          percentFixed = routeDD.fixArcs(repairedDual, LPSolveType::LAGSolver);
          std::cout << "percent fixed: " << percentFixed << std::endl;
        }
        if (percentFixed > bestDualArcFixingPercent)
        {
          std::cout << "updating best lambda arc fixing" << std::endl;
          bestDualArcFixing = repairedDual;
          bestDualArcFixingPercent = percentFixed;
        }
      }
      else
      {
        // Repair can be costly so only try when within a percent of the optimality gap
        if (params.useVariableFixing && ((stats.upperBound - lagrangeanLowerBound) * 100.0 / stats.upperBound < params.lagOptimalityGapToStartRepairing))
        {
          repairMultipliers(repairedDual, LPSolveType::LAGSolver);
          double repairedBound = routeDD.getDualObjectiveValue(repairedDual, LPSolveType::LAGSolver);

          DBG(std::cout << "repaired lb: " << repairedBound << std::endl;)
          if (bestDualValue < repairedBound)
          {
            stats.lowerBound = std::max(stats.lowerBound, repairedBound);
            stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
            printMultipliers(repairedDual);

            bestDual = repairedDual;
            bestDualValue = repairedBound;
          }

          percentFixed = routeDD.fixArcs(repairedDual, LPSolveType::LAGSolver);
          std::cout << "percent fixed: " << percentFixed << std::endl;
          if (percentFixed > bestDualArcFixingPercent)
          {
            std::cout << "updating best lambda arc fixing" << std::endl;
            bestDualArcFixingPercent = percentFixed;
            bestDualArcFixing = repairedDual;
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

    // Separations
    if ((phaseType == PhaseType::SEPARATION) && !shouldTerminate)
    {
      // delay lagrangean separations to try to keep DD small
      if (stats.numLagIterations < params.lagIterationDelayToStartSeparating)
      {
        infeasibleRoutesByArc.clear();
      }

      for (auto infeasibleRouteToSeparateByArc : infeasibleRoutesByArc)
      {
        // when separating changes, can use locations to re-find the route
        if (routeDD.doesRouteExistByArcs(infeasibleRouteToSeparateByArc))
        {
          stats.numSeparations = stats.numSeparations + 1;
          routeDD.separateInfeasibleRoute(infeasibleRouteToSeparateByArc, params.maxS);
        }
      }
    }

    // Cuts
    if (!shouldTerminate && (params.numLagItersForCuts != 0) && (stats.numLagIterations % params.numLagItersForCuts == 0))
    {
      addCutsUsingCurrentPrimal(dual);
    }

    if (shouldTerminate && params.useCuts)
    {
      addCutsUsingCurrentPrimal(dual);
    }

    // Phase changes
    if (shouldTerminate)
    {
      if (phaseType == PhaseType::INITIAL_DUAL)
      {
        std::cout << "start separation phase" << std::endl;
        phaseType = PhaseType::SEPARATION;
      }
      else if (phaseType == PhaseType::SEPARATION)
      {
        std::cout << "start robust primal phase" << std::endl;
        phaseType = PhaseType::ROBUST_PRIMAL;

        // fresh primal and use primal to add cuts
        primal = Primal();
        params.useCuts = true;
      }
      else if (phaseType == PhaseType::ROBUST_PRIMAL)
      {
        std::cout << "start robust cut dual phase" << std::endl;
        phaseType = PhaseType::ROBUST_CUT_DUALS;

        // fresh step sizes and duals
        /*
        stepSizeMultiplier = 1.0;
        stepSizeMultiplierIteration = 0;
        stepSizeMultiplierIterationCutoff = 100;
        alphaLowerBound = 0.1;
        alphaLowerBoundIteration = 0;
        alphaLowerBoundCheckValue = 0;

        initializeDual(dual);
        repairMultipliers(dual, LPSolveType::LAGSolver);
        bestDual = dual;
        bestDualValue = routeDD.getDualObjectiveValue(bestDual, LPSolveType::LAGSolver);
        targetLowerBound = bestDualValue * 1.05;
        */

        params.useCuts = false;
      }
      else if (phaseType == PhaseType::ROBUST_CUT_DUALS)
      {
        std::cout << "start nonrobust primal phase" << std::endl;
        phaseType = PhaseType::NONROBUST_PRIMAL;
        params.useCuts = true;
      }
      else if (phaseType == PhaseType::NONROBUST_PRIMAL)
      {
        std::cout << "start nonrobust cut dual phase" << std::endl;
        phaseType = PhaseType::NONROBUST_CUT_DUALS;

        // fresh step sizes and duals
        /*
        stepSizeMultiplier = 1.0;
        stepSizeMultiplierIteration = 0;
        stepSizeMultiplierIterationCutoff = 20;
        alphaLowerBound = 0.1;
        alphaLowerBoundIteration = 0;
        alphaLowerBoundCheckValue = 0;
        bestDualValue = 0.0;
        muPercentImproved = muPercentImproved / 2.0;

        initializeDual(dual);
        */
        params.useCuts = false;
      }
      else if (phaseType == PhaseType::NONROBUST_CUT_DUALS)
      {
        std::cout << "start nonrobust primal phase again" << std::endl;
        phaseType = PhaseType::NONROBUST_PRIMAL;
        params.useCuts = true;
      }
    }

    infeasibleRoutesByArc.clear();
  }

  return true;
}

void VRPTWDDSolver::constructNextPrimal(double alphaTry, const std::vector<std::vector<int>>& decomposedRoutes, const std::vector<std::vector<int>>& decomposedRouteArcs, Primal& nextPrimal)
{
  for (int flowIndex=0; flowIndex<nextPrimal.xDecompositionFlows.size(); ++flowIndex)
  {
    nextPrimal.xDecompositionFlows[flowIndex] = (1 - alphaTry) * nextPrimal.xDecompositionFlows[flowIndex];
  }

  // add feasible versions if infeasible
  for (int index=0; index<decomposedRoutes.size(); ++index)
  {
    auto route = decomposedRoutes[index];

    // check to add flow to current set of routes
    bool alreadyExists = false;
    for (int currIndex=0; currIndex<nextPrimal.xDecompositionFlows.size(); ++currIndex)
    {
      if (route == nextPrimal.xDecompositions[currIndex])
      {
        nextPrimal.xDecompositionFlows[currIndex] += alphaTry;
        alreadyExists = true;
        break;
      }
    }

    // for infeasible routes, add relevant / closests feasible route
    if (!alreadyExists)
    {
      //if (routeDD.isRouteFeasible(route))
      //{
      auto routeArcs = decomposedRouteArcs[index];
      nextPrimal.xDecompositions.push_back(route);
      nextPrimal.xDecompositionArcs.push_back(routeArcs);
      nextPrimal.xDecompositionFlows.push_back(alphaTry);
      //}
      /*
      else
      {
        // repair route and then add
        std::vector<int> feasibleRoute;
        std::vector<int> feasibleRouteArcs;
        routeDD.repairRoute(route, feasibleRoute, feasibleRouteArcs);

        // check to add flow to current set of routes
        bool repairAlreadyExists = false;
        for (int currIndex=0; currIndex<nextPrimal.xDecompositionFlows.size(); ++currIndex)
        {
          if (feasibleRoute == nextPrimal.xDecompositions[currIndex])
          {
            nextPrimal.xDecompositionFlows[currIndex] += alphaTry;
            repairAlreadyExists = true;
            break;
          }
        }

        if (!repairAlreadyExists)
        {
          nextPrimal.xDecompositions.push_back(feasibleRoute);
          nextPrimal.xDecompositionArcs.push_back(feasibleRouteArcs);
          nextPrimal.xDecompositionFlows.push_back(alphaTry);
        }
      }
      */
    }
  }
};

double VRPTWDDSolver::calculateTwoNorm(const std::vector<double>& gamma)
{
  double total = 0.0;
  for (double g : gamma)
  {
    total = total + g*g;
  }

  return std::sqrt(total);
};

double VRPTWDDSolver::calculateDotProduct(const std::vector<double>& vector1, const std::vector<double>& vector2)
{
  double total = 0.0;
  for (int index=0; index<vector1.size(); ++index)
  {
    total = total + vector1[index] * vector2[index];
  }

  return total;
};

void VRPTWDDSolver::getGradient(const Primal& currPrimal, const Dual& dual, std::vector<double>& gradient)
{
  // get values of LHS for all dualized inequalities
  std::unordered_map<int,double> locationsCovered;
  routeDD.getNumberOfTimesLocationsCoveredRoutes(currPrimal, locationsCovered);

  std::vector<double> cutValues;
  routeDD.getCutSetValuesRoutes(currPrimal, cutValues);
 
  std::vector<double> combValues;
  routeDD.getCombValuesRoutes(currPrimal, combValues);

  std::vector<double> cliqueCutValues;
  routeDD.getCliqueCutValuesRoutes(currPrimal, cliqueCutValues);

  // gradient_(k) = b - Ax_(k)
  gradient.resize(vrptw.numLocations + dual.capDuals.size() + dual.srcDuals.size() + dual.combDuals.size(), 0);
  for (int i=1; i<vrptw.numLocations; ++i)
  {
    gradient[i] = 1 - locationsCovered[i];
  }

  // cap cuts are in the form Cx <= r so change to -Cx >= -r
  for (int i=0; i<dual.capDuals.size(); ++i)
  {
    int gradientIndex = i + vrptw.numLocations;
    if (routeDD.isCapCutActive(i))
    {
      gradient[gradientIndex] = (-1 * routeDD.getCapCutSetRHS(i)) + cutValues[i];
    }
    else
    {
      gradient[gradientIndex] = 0;
    }
  }

  // src cuts are in the form Cx <= r so chang eto -Cx >= -r
  for (int i=0; i<dual.srcDuals.size(); ++i)
  {
    int gradientIndex = i + vrptw.numLocations + dual.capDuals.size();
    SRCType srcType = routeDD.getSRCType(i);
    int rhs = routeDD.getSRCRHS(srcType);
    gradient[gradientIndex] = (-1*rhs) + cliqueCutValues[i];
  }

  for (int i=0; i<dual.combDuals.size(); ++i)
  {
    int gradientIndex = i + vrptw.numLocations + dual.capDuals.size() + dual.srcDuals.size();
    gradient[gradientIndex] = routeDD.getCombCutRHS(i) - combValues[i];
  }
};

void VRPTWDDSolver::updateMultipliersVolumeAlgorithm(Dual& dual, Primal& currPrimal, double lagrangeanLowerBound, int iteration)
{
  // get values of LHS for all dualized inequalities
  std::unordered_map<int,double> locationsCovered;
  routeDD.getNumberOfTimesLocationsCoveredRoutes(currPrimal, locationsCovered);

  std::vector<double> cutValues;
  routeDD.getCutSetValuesRoutes(currPrimal, cutValues);
 
  std::vector<double> combValues;
  routeDD.getCombValuesRoutes(currPrimal, combValues);

  std::vector<double> cliqueCutValues;
  routeDD.getCliqueCutValuesRoutes(currPrimal, cliqueCutValues);

  // gamma_(k) = b - Ax_(k)
  // Beasley - when multiplier is already 0 and step direction is negative, don't include in ||gamma||^2
  double normGammaSquared = 0;
  std::vector<double> gamma(vrptw.numLocations + dual.capDuals.size() + dual.srcDuals.size() + dual.combDuals.size(), 0);
  for (int i=1; i<vrptw.numLocations; ++i)
  {
    gamma[i] = 1 - locationsCovered[i];
    if ((bestDual.lambda[i] <= 0.000001) && (gamma[i] <= 0))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[i], 2);
    }
  }

  // cap cuts are in the form Cx <= r so change to -Cx >= -r
  // scale down to RHS 1 like the other constraints
  for (int i=0; i<dual.capDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations;
    if (routeDD.isCapCutActive(i))
    {
      gamma[gammaIndex] = (-1 * routeDD.getCapCutSetRHS(i)) + cutValues[i];

      if ((bestDual.capDuals[i] <= 0.1) && (gamma[gammaIndex] <= 0))
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
      bestDual.capDuals[i] = 0;
      gamma[gammaIndex] = 0;
    }
  }

  // src cuts are in the form Cx <= r so chang eto -Cx >= -r
  for (int i=0; i<dual.srcDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations + dual.capDuals.size();
    SRCType srcType = routeDD.getSRCType(i);
    int rhs = routeDD.getSRCRHS(srcType);
    gamma[gammaIndex] = (-1*rhs) + cliqueCutValues[i];

    if ((bestDual.srcDuals[i] <= 0.01) && (gamma[gammaIndex] <= 0.01))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
  }

  for (int i=0; i<dual.combDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations + dual.capDuals.size() + dual.srcDuals.size();
    gamma[gammaIndex] = routeDD.getCombCutRHS(i) - combValues[i];
    if ((bestDual.combDuals[i] <= 0.01) && (gamma[gammaIndex] <= 0.01))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
  }

  for (int index=0; index<gamma.size(); ++index)
  {
    std::cout << index << ":" << gamma[index] << std::endl;
  }

  double alpha = stepSizeMultiplier * (targetLowerBound - bestDualValue) / normGammaSquared;
  std::cout << "alpha: " << alpha << std::endl;
  std::cout << "lagLB: " << lagrangeanLowerBound << std::endl;
  std::cout << "stepSizeMultiplier: " << stepSizeMultiplier << std::endl;
  std::cout << "targetLowerBound: " << targetLowerBound << std::endl;
  std::cout << "||gamma||^2: " << normGammaSquared << std::endl;

  printMultipliers(dual);
  printMultipliers(bestDual);

  // lambda_(k+1) = lambda_(k) + alpha_(k) * gamma_(k)
  for (int i=1; i<vrptw.numLocations; ++i)
  {
    dual.lambda[i] = std::max(0.0, bestDual.lambda[i] + alpha * gamma[i]);
  }

  // same for mu
  for (int i=0; i<dual.capDuals.size(); ++i)
  {
    dual.capDuals[i] = std::max(0.0, bestDual.capDuals[i] + alpha * gamma[i+vrptw.numLocations]);
  }
 
  // same for srcDuals
  for (int i=0; i<dual.srcDuals.size(); ++i)
  {
    dual.srcDuals[i] = std::max(0.0, bestDual.srcDuals[i] + alpha * gamma[i+vrptw.numLocations+dual.capDuals.size()]);
  }

  // same for combs
  for (int i=0; i<dual.combDuals.size(); ++i)
  {
    dual.combDuals[i] = std::max(0.0, bestDual.combDuals[i] + alpha * gamma[i+vrptw.numLocations+dual.capDuals.size()+dual.srcDuals.size()]);
  }

  // remove if too small for too long
  for (int muIndex=0; muIndex<dual.capDuals.size(); ++muIndex)
  {
    if (routeDD.isCapCutActive(muIndex))
    {
      if (dual.capDuals[muIndex] < params.deactivateCutValueThreshold)
      {
        capCutTooSmallCounters[muIndex] = capCutTooSmallCounters[muIndex] + 1;
        if (capCutTooSmallCounters[muIndex] > params.deactivateCutIterThreshold)
        {
          std::cout << "deactivated cut at index: " << muIndex << std::endl;
          routeDD.deactivateCapCut(muIndex);
          stats.numCuts = stats.numCuts - 1;
        }
      }
      else
      {
        capCutTooSmallCounters[muIndex] = 0;
      }
    }
  }
 
  // remove if too small for too long
  for (int srcIndex=0; srcIndex<dual.srcDuals.size(); ++srcIndex)
  {
    if (routeDD.isCliqueCutActive(srcIndex))
    {
      if (dual.srcDuals[srcIndex] < params.deactivateCutValueThreshold)
      {
        cliqueCutTooSmallCounters[srcIndex] = cliqueCutTooSmallCounters[srcIndex] + 1;
        if (cliqueCutTooSmallCounters[srcIndex] > params.deactivateCutIterThreshold)
        {
          std::cout << "deactivated src cut at index: " << srcIndex << std::endl;
          routeDD.deactivateCliqueCut(srcIndex);
          stats.numCuts = stats.numCuts - 1;
        }
      }
      else
      {
        cliqueCutTooSmallCounters[srcIndex] = 0;
      }
    }
  }
};

bool VRPTWDDSolver::addCutsUsingCurrentPrimal(Dual& dual)
{
  // Rounded Capacity Cuts
  bool cutAdded = false;
  std::vector<int> edgeTail;
  std::vector<int> edgeHead;
  std::vector<double> edgeFlow;
  std::vector<int> rccArcs;
  std::vector<double> rccArcFlows;
  if ((phaseType == PhaseType::ROBUST_PRIMAL) || (phaseType == PhaseType::NONROBUST_PRIMAL))
  {
    convertArcIndicesForVRPTWSep(primal, edgeTail, edgeHead, edgeFlow, rccArcs, rccArcFlows);
    addRCCs(edgeTail, edgeHead, edgeFlow, rccArcs, params.numLagCuts, cutAdded, dual);

    // Strengthened Combs
    //addCombs(edgeTail, edgeHead, edgeFlow, cutAdded);
    //combDuals.resize(routeDD.getNumCombCuts());
  }

  // Subset Row Cuts
  double totalFlow = 0;
  std::vector<double> flowByVertex(vrptw.numLocations, 0);
  for (int index=0; index<primal.xDecompositionFlows.size(); ++index)
  {
    double flow = primal.xDecompositionFlows[index];
    totalFlow += flow;
    auto route = primal.xDecompositions[index];
    std::cout << "route with flow " << flow << " : ";
    for (int loc : route)
    {
      flowByVertex[loc] += flow;
      std::cout << loc << ",";
    }
    std::cout << " ";
    auto routeArcs = primal.xDecompositionArcs[index];
    for (int arcIndex : routeArcs)
    {
      std::cout << arcIndex << ",";
    }
    std::cout << std::endl;
  }
  std::cout << "total flow: " << totalFlow << std::endl;
  for (int index=0; index<vrptw.numLocations; ++index)
  {
    std::cout << "loc: " << index << " flow: " << flowByVertex[index] << std::endl;
  }

  if (phaseType == PhaseType::NONROBUST_PRIMAL)
  {
    std::vector<double> violations;
    std::cout << "finding src3s" << std::endl;
    int numSrc3Added = routeDD.findSRC3s(primal, 2, violations);
    std::cout << "finding src4s" << std::endl;
    int numSrc4Added = routeDD.findSRC4s(primal, 1, violations);

    int numSrcAdded = numSrc3Added + numSrc4Added;
    if (numSrcAdded > 0)
    {
      cutAdded = true;
      stats.numCuts += numSrcAdded;
      std::cout << "adding src cuts" << std::endl;
      addSRCCuts(dual.srcDuals, violations);

      // strengthen src using average route length
      int averageRouteLength = 0;
      for (int index=0; index<primal.xDecompositions.size(); ++index)
      {
        averageRouteLength += primal.xDecompositions[index].size();
      }
      averageRouteLength = averageRouteLength / primal.xDecompositions.size();
      routeDD.strengthenSRCs(averageRouteLength);
    }
  }

  resizeMultipliers(dual, bestDual);
  resizeMultipliers(dual, bestDualArcFixing);
  capCutTooSmallCounters.resize(dual.capDuals.size());
  cliqueCutTooSmallCounters.resize(dual.srcDuals.size());

  return cutAdded;
};
