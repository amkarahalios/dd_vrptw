#include "vrptwddsolver.h"
#include <algorithm>
#include <chrono>
#include <random>
#include <math.h>
#include <ilcplex/ilocplex.h>

void DDStats::print(int ddNumArcs, int ddNumFixedArcs) const
{
  std::cout << "STATS - lpIterations[" << lpIterations << "] lagIterations[";
  std::cout << numLagIterations << "] sspIterations[" << numSSPIterations;
  std::cout << "] numSeparations[" << numSeparations << "] numCuts[" << numCuts;
  std::cout << "] compileTime[" << int(millisecondsCompiling / 1000);
  std::cout << "] sspSolveTime[" << int(millisecondsSolvingSSP / 1000);
  std::cout << "] lpSolveTime[" << int(millisecondsSolvingLP / 1000);
  std::cout << "] lb[" << lowerBound << "] ub[" << upperBound << "]";
  std::cout << " numArcs: [" << ddNumArcs << "]";
  std::cout << " numFixed: [" << ddNumFixedArcs << "]";
  std::cout << " numHeuristicIPs: [" << numHeuristicIPs << "]";
  std::cout << " numHeuristicLNS: [" << numHeuristicLNSs << "]";
  std::cout << " numPrimalLNSRepairs: [" << numPrimalLNSRepairs << "]";
  std::cout << " time: [" << getNumSeconds() << "]" << std::endl;

  std::cout << "STATS1 - lpIterations[" << lpIterations << "] lagIterations[";
  std::cout << numLagIterations << "] sspIterations[" << numSSPIterations;
  std::cout << "] numSeparations[" << numSeparations << "] numCuts[" << numCuts;
  std::cout << "] compileTime[" << int(millisecondsCompiling / 1000);
  std::cout << "] sspSolveTime[" << int(millisecondsSolvingSSP / 1000);
  std::cout << "] lpSolveTime[" << int(millisecondsSolvingLP / 1000);
  std::cout << "] repairTime[" << int(millisecondsRepairingLAG / 1000);
  std::cout << "] cutsTime[" << int(millisecondsFindingCuts / 1000);
  std::cout << "] alphaTime[" << int(millisecondsTryingAlpha / 1000);
  std::cout << "] decompTime[" << int(millisecondsDecompose / 1000);
  std::cout << "] yellowTime[" << int(millisecondsYellow / 1000);
  std::cout << "] updateTime[" << int(millisecondsUpdateDual / 1000);
  std::cout << "] fixTime[" << int(millisecondsFix/ 1000);
  std::cout << "] lb[" << lowerBound << "] ub[" << upperBound << "]";
  std::cout << " numArcs: [" << ddNumArcs << "]";
  std::cout << " numFixed: [" << ddNumFixedArcs << "]";
  std::cout << " time: [" << getNumSeconds() << "]" << std::endl;
}
 
int DDStats::getNumSeconds() const
{
  auto currTime = std::chrono::high_resolution_clock::now();
  return int(std::chrono::duration_cast<std::chrono::milliseconds>(currTime - startTime).count() / 1000);
}

VRPTWDDSolver::VRPTWDDSolver(VRPTW _vrptw, VRPTWDDParameters _params) : vrptw(_vrptw), routeDD(_vrptw, _params), params(_params)
{
  if (params.primalHeuristic == PrimalHeuristic::BEST_KNOWN)
  {
    stats.upperBound = vrptw.instanceUpperBound + 1;
  }
  else
  {
    stats.upperBound = INF;
  }

  // Seed for random components
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::srand(seed);
  std::cout << "Random seed:" << seed << std::endl;

  // Compile DD
  std::cout << "compiling DD" << std::endl;
  if (params.primalHeuristic != PrimalHeuristic::STANDALONE)
  {
    auto startCompileTime = std::chrono::high_resolution_clock::now();
    if (params.stateSpace == StateSpace::NG)
    {
      routeDD.compileNgRoute(params.ngSetSize);
    }
    auto endCompileTime = std::chrono::high_resolution_clock::now();
    auto compileSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endCompileTime - startCompileTime).count();
    stats.millisecondsCompiling = stats.millisecondsCompiling + compileSolveTime;
  }
  else
  {
    routeDD.compileEmpty();
  }
  std::cout << "done compiling DD" << std::endl;

  // Primal heuristic
  if ((params.primalHeuristic != PrimalHeuristic::BEST_KNOWN) && (params.primalHeuristic != PrimalHeuristic::STANDALONE))
  {
    std::cout << "running primal heuristic" << std::endl;

    // Standalone solver runs Heuristic + MIP Pool + LNS in a loop
    int timeToRunHeuristics = params.lnsTimeoutSeconds;
    if (params.primalHeuristic == PrimalHeuristic::STANDALONE)
    {
      timeToRunHeuristics = params.timeoutSeconds;
    }

    // Enter LNS process
    int heuristicTimeSeconds = 0;
    auto heuristicStartTime = std::chrono::high_resolution_clock::now();
    while (heuristicTimeSeconds < timeToRunHeuristics)
    {
      // Heuristic for Initial Routes
      std::vector<std::vector<int>> heuristicRoutes;
      bool isHeuristicFeasible = routeDD.generateHeuristicRoutesLiterature(heuristicRoutes);
      if (!isHeuristicFeasible)
      {
        continue;
      }
      for (auto route : heuristicRoutes)
      {
        addRouteToPrimalRoutes(route);
      }
      double heuristicRoutesUpperBound = vrptw.evaluateSolutionCost(heuristicRoutes);

      // MIP Pool
      std::vector<std::vector<int>> mipPoolSolution;
      bool primalHeuristicFeasible = primalHeuristicMIP(FlowType::IP, mipPoolSolution, mipPoolSolutionIndices);
      if (primalHeuristicFeasible)
      {
        double mipPoolUpperBound = vrptw.evaluateSolutionCost(mipPoolSolution);
        if (mipPoolUpperBound < stats.upperBound)
        {
          std::cout << "improved ub from MIP pool: " << mipPoolUpperBound << std::endl;
          stats.upperBound = mipPoolUpperBound;
          stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
          heuristicRoutes = mipPoolSolution;
          heuristicStartTime = std::chrono::high_resolution_clock::now();
        }
      }
   
      // LNS
      Dual dual;
      dual.lambda.resize(vrptw.numLocations);
      if (params.primalHeuristic == PrimalHeuristic::CE_MIP_LNS)
      {
        bool isImproved = largeNeighborhoodSearch(heuristicRoutes, dual, params.lnsTimeoutSeconds);
        if (isImproved && (params.primalHeuristic == PrimalHeuristic::STANDALONE))
        {
          heuristicStartTime = std::chrono::high_resolution_clock::now();
        }
      }
   
      auto heuristicCurrTime = std::chrono::high_resolution_clock::now();
      heuristicTimeSeconds = std::chrono::duration_cast<std::chrono::seconds>(heuristicCurrTime- heuristicStartTime).count();
    }
  }

  // Set and Log Parameters
  infeasibleRoutesBatchSize = 1;
  percentFixedToChangeToCPLEX = 95;
  numArcsToChangeToCPLEX = 100000;

  std::cout << "batch size for lag: " << infeasibleRoutesBatchSize << std::endl;
  std::cout << "percent arcs fixed to change to CPLEX: " << percentFixedToChangeToCPLEX << std::endl;
  std::cout << "number arcs fixed to change to CPLEX: " << numArcsToChangeToCPLEX << std::endl;

  // best duals for fixing
  bestDualsArcFixingLagPercent = 0.0;
  bestDualArcFixingLPPercent = -1;
  bestDualValue = 0.0;

  // va
  stepSizeMultiplier = 1.0;
  stepSizeMultiplierIteration = 0;
  stepSizeMultiplierIterationCutoff = 5;
  alphaLowerBound = 0.1;
  alphaLowerBoundIteration = 0;
  targetLowerBound = 1.0;
  targetBoundIncrease = 1.05;
  minStepSizeMultiplier = 0.0001;

  // va parameters
  std::cout << "step size iteration cutoff: " << stepSizeMultiplierIterationCutoff << std::endl;
  std::cout << "target bound increase: " << targetBoundIncrease << std::endl;
  std::cout << "min step size multiplier: " << minStepSizeMultiplier << std::endl;

  // dynamic parameters
  std::cout << "ngSetSize: " << params.ngSetSize << std::endl;
  std::cout << "change to LP: " << params.changeToLP << std::endl;
  std::cout << "var fixing: " << params.useVariableFixing << std::endl;
  std::cout << "useSeparations: " << params.useSeparations << std::endl;
  std::cout << "muSSP: " << params.useMuSSP << std::endl;

  // cuts parameters
  std::cout << "numLagCuts: " << params.numLagCuts << std::endl;
  std::cout << "RCCs: " << params.useRCCs << std::endl;
  std::cout << "SRC3s: " << params.useSRC3s << std::endl;
  std::cout << "SRC4s: " << params.useSRC4s << std::endl;
  std::cout << "SRC5V1s: " << params.useSRC5V1s << std::endl;
  std::cout << "SRC5V2s: " << params.useSRC5V2s << std::endl;
  std::cout << "useVolumeAlgorithm: " << params.useVolumeAlgorithm << std::endl;
  std::cout << "limit rccs: " << params.limitRCCs << std::endl;
  std::cout << "useScaling: " << params.useScaling << std::endl;
  std::cout << "lowImprovementIterationsToNextRound: " << params.lowImprovementIterationsToNextRound << std::endl;
  std::cout << "percentImprovementThreshold: " << params.percentImprovementThreshold << std::endl;

  // primal heuristic params
  std::cout << "primal heuristic: " << params.primalHeuristic << std::endl;
  std::cout << "lns timeout: " << params.lnsTimeoutSeconds << std::endl;
  std::cout << "removal strategy: " << params.removalStrategy << std::endl;
  std::cout << "insertion ablation: " << params.insertionAblation << std::endl;

  CMGR_CreateCMgr(&MyCutsCMP,Dim);
  CMGR_CreateCMgr(&MyOldCutsCMP,Dim);

  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
};

bool VRPTWDDSolver::solve(bool shouldSolveIP)
{
  bool lpFlowType = FlowType::LP;

  Dual dual;
  initializeDual(dual);
  bestDual = dual;
  bestDualArcFixingLP = dual;

  SGDAlgorithm sgdAlgo = SGDAlgorithm::SGD;
  if (params.useVolumeAlgorithm)
  {
    sgdAlgo = SGDAlgorithm::VA;
  }

  int averageRouteLength = 0;
  int numRoutesInAverage = 0;
  bool changedLagToLP = false;
  bool finishedSolving = false;
  while (!finishedSolving)
  {
    bool solved = false;
    if (params.lpSolveType == LPSolveType::LPSolver)
    {
      double percentArcsFixed = 0.0;
      if (!changedLagToLP)
      {
        if (params.useVariableFixing)
        {
          if (bestDualArcFixingLPPercent > 0)
          {
            repairMultipliers(bestDualArcFixingLP, LPSolveType::LPSolver);
            percentArcsFixed = routeDD.fixArcs(bestDualArcFixingLP, LPSolveType::LPSolver, stats.upperBound);
          }

          repairMultipliers(dual, LPSolveType::LPSolver);
          percentArcsFixed = routeDD.fixArcs(dual, LPSolveType::LPSolver, stats.upperBound);
          if (percentArcsFixed > bestDualArcFixingLPPercent)
          {
            bestDualArcFixingLPPercent = percentArcsFixed;
            bestDualArcFixingLP = dual;
            std::cout << "update best LP dual for fixing. Percent: " << bestDualArcFixingLPPercent << std::endl;
          }
        }
      }
      else
      {
        // TODO(akarahal) if using a counter to avoid loops for LAG solve...
        // ...can now get rid of the counters by merging nodes!
        changedLagToLP = false;
      }

      // TODO(akarahal) store best duals for arc fixing better
      if (params.useVariableFixing)
      {
        std::cout << "fix with best duals from LAG. Using: " << bestDualsArcFixingLag.size() << " dual(s)" << std::endl;
        for (int dualMultiplierIndex=0; dualMultiplierIndex<bestDualsArcFixingLag.size(); ++dualMultiplierIndex)
        {
          Dual dualToRepair = bestDualsArcFixingLag[dualMultiplierIndex];
          repairMultipliers(dualToRepair, LPSolveType::LAGSolver);
          bestDualsArcFixingLag[dualMultiplierIndex] = dualToRepair;
        }
        routeDD.fixArcs(bestDualsArcFixingLag, LPSolveType::LAGSolver, stats.upperBound);
      }

      if (lpFlowType == FlowType::LP)
      {
        solved = solveLP(dual);
      }
      else
      {
        solved = solveIP(dual);
      }
    }
    else if (params.lpSolveType == LPSolveType::LAGSolver)
    {
      solved = solveLagrangeanRelaxation(dual, sgdAlgo);

      if (params.lpSolveType == LPSolveType::LPSolver)
      {
        changedLagToLP = true;
        continue;
      }
    }

    if (!solved)
    {
      std::cout << "finished solving - error" << std::endl;
      return false;
    }

    // Primal Heuristic
    std::vector<std::vector<int>> routesByLocationPrimalHeuristic;
    bool primalHeuristicFeasible = false;
    if (params.primalHeuristic == PrimalHeuristic::CE_GREEDY)
    {
      //primalHeuristicFeasible = routeDD.primalHeuristicGreedy(routesByLocationPrimalHeuristic);
    }
    //else if ((params.primalHeuristic == PrimalHeuristic::CE_MIP) || (params.primalHeuristic == PrimalHeuristic::CE_MIP_LNS) || (params.primalHeuristic == PrimalHeuristic::BEST_KNOWN))
    else if ((params.primalHeuristic == PrimalHeuristic::CE_MIP) || (params.primalHeuristic == PrimalHeuristic::CE_MIP_LNS))
    {
      primalHeuristicFeasible = primalHeuristicMIP(FlowType::IP, routesByLocationPrimalHeuristic, mipPoolSolutionIndices);
    }

    // Evaluate Solution
    if (primalHeuristicFeasible)
    {
      double heuristicUpperBound = vrptw.evaluateSolutionCost(routesByLocationPrimalHeuristic);
      if (heuristicUpperBound < stats.upperBound)
      {
        stats.upperBound = heuristicUpperBound;
        std::cout << "improved ub from MIP primal heuristic: " << heuristicUpperBound << std::endl;
        stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
      }
 
      // Run LNS
      if (params.primalHeuristic == PrimalHeuristic::CE_MIP_LNS)
      {
        largeNeighborhoodSearch(routesByLocationPrimalHeuristic, dual, params.lnsTimeoutSeconds);
      }
    }

    // Infeasibilities for LP solve
    if (params.lpSolveType == LPSolveType::LPSolver)
    {
      bool cutAdded = false;
      std::vector<std::vector<int>> decomposedRoutes;
      std::vector<std::vector<int>> decomposedArcs;
      std::vector<double> routeFlows;

      if (lpFlowType == FlowType::LP)
      {
        if (params.useRCCs)
        {
          auto startTimeCut = std::chrono::high_resolution_clock::now();
          // RCC - Rounded Capacity Cuts
          std::vector<int> edgeTail;
          std::vector<int> edgeHead;
          std::vector<double> edgeFlow;
          routeDD.convertSolutionForVRPTWSep(edgeTail, edgeHead, edgeFlow);

          std::vector<int> sequenceArcs;
          if (vrptw.problemType != ProblemType::PDP)
          {
            addRCCs(edgeTail, edgeHead, edgeFlow, sequenceArcs, 100, cutAdded, dual);
          }

          resizeMultipliers(dual, bestDual);
          resizeMultipliers(dual, bestDualsArcFixingLag);
          resizeMultipliers(dual, bestDualArcFixingLP);

          // Strengthened Combs
          //addCombs(edgeTail, edgeHead, edgeFlow, cutAdded);
          //dual.combDuals.resize(routeDD.getNumCombCuts());
          auto endTimeCut = std::chrono::high_resolution_clock::now();
          auto totalTimeCut = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeCut - startTimeCut).count();
          stats.millisecondsFindingCuts += totalTimeCut;
        }

        // SRC - Subset Row Cuts
        // Part 1 - calculate primals, don't separate until after decomposition separations happen
        if (params.useSRC3s || params.useSRC4s || params.useSRC5V1s || params.useSRC5V2s)
        {
          auto startTimeCut = std::chrono::high_resolution_clock::now();
          std::vector<int> infeasibleRoute;
          routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, DecompositionReason::DECOMPOSE);
          primal = Primal();
          for (int index=0; index<decomposedRoutes.size(); ++index)
          {
            auto route = decomposedRoutes[index];
            auto routeArcs = decomposedArcs[index];
            auto routeFlow = routeFlows[index];
            primal.xDecompositions.push_back(route);
            primal.xDecompositionArcs.push_back(routeArcs);
            primal.xDecompositionFlows.push_back(routeFlow);
          }

          auto endTimeCut = std::chrono::high_resolution_clock::now();
          auto totalTimeCut = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeCut - startTimeCut).count();
          stats.millisecondsFindingCuts += totalTimeCut;
        }
      }

      bool stopFindingInfeasibilities = false;
      std::vector<std::vector<int>> infeasibilities;
      int numPrimalRoutes = primalRoutes.size();
      if (params.useSeparations)
      {
        while (!stopFindingInfeasibilities)
        {
          std::vector<int> infeasibleRoute;
          routeFlows.clear();
          decomposedRoutes.clear();
          decomposedArcs.clear();
          routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, DecompositionReason::DECOMPOSE);
          if ((params.primalHeuristic != PrimalHeuristic::BEST_KNOWN) && (params.primalHeuristic != PrimalHeuristic::STANDALONE))
          {
            for (auto route : decomposedRoutes)
            {
              if (routeDD.isRouteFeasible(route))
              {
                addRouteToPrimalRoutes(route);
              }
              else
              {
                if (params.insertionAblation != PrimalHeuristicInsertionAblation::TRUNCATED)
                {
                  std::vector<int> truncatedRoute;
                  routeDD.createTruncatedRoute(route, truncatedRoute);
                  addRouteToPrimalRoutes(truncatedRoute);
                }
              }
            }
          }

          infeasibleRoute.clear();
          routeFlows.clear();
          decomposedRoutes.clear();
          decomposedArcs.clear();
          routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, DecompositionReason::SEPARATE);
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

      if (infeasibilities.size() > 0)
      {
        for (auto infeasibleRoute : infeasibilities)
        {
          std::vector<int> emptyRoute;
          if (routeDD.doesRouteExistByArcs(infeasibleRoute, emptyRoute))
          {
            stats.numSeparations = stats.numSeparations + 1;
            routeDD.separateRoute(infeasibleRoute);
          }
        }
      }
      else
      {
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
 
      // SRC - Subset Row Cuts
      // Part 2 - these may require separations, which may alter the structure
      if (params.useSRC3s || params.useSRC4s || params.useSRC5V1s || params.useSRC5V2s)
      {
        // Should check decomposition first for violated cuts... if exist, then separate!

        // Convert primal routes in case some infeasibilities were separated
        // Now we have all feasible or truncated routes
        std::vector<std::vector<int>> truncatedExactSequences;
        std::vector<std::vector<int>> truncatedExactSequencesArcs;
        separateSequencesAndTruncate(primal.xDecompositions, truncatedExactSequences, truncatedExactSequencesArcs);
        Primal cutPrimal;
        cutPrimal.xDecompositions = truncatedExactSequences;
        cutPrimal.xDecompositionArcs = truncatedExactSequencesArcs;
        cutPrimal.xDecompositionFlows = primal.xDecompositionFlows;

        auto startTimeCut = std::chrono::high_resolution_clock::now();
        std::vector<double> violations;
        int numSrcAdded = 0;
        numSrcAdded = routeDD.findSRCs(cutPrimal, 100, violations);
        if (numSrcAdded > 0)
        {
          cutAdded = true;
          stats.numCuts = stats.numCuts + numSrcAdded;
        }

        std::vector<double> zeroViolations;
        zeroViolations.resize(violations.size());
        addSRCCuts(dual.srcDuals, zeroViolations);

        resizeMultipliers(dual, bestDual);
        resizeMultipliers(dual, bestDualsArcFixingLag);
        resizeMultipliers(dual, bestDualArcFixingLP);

        auto endTimeCut = std::chrono::high_resolution_clock::now();
        auto totalTimeCut = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeCut - startTimeCut).count();
        stats.millisecondsFindingCuts += totalTimeCut;
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

    if (stats.lowerBound + 0.01 >= stats.upperBound)
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
  const std::set<int> initialPrimalArcIndices;
  stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, initialPrimalArcIndices, duals, false, params.timeoutSeconds);
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
    const std::set<int> initialPrimalArcIndices;
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, initialPrimalArcIndices, duals, false, params.timeoutSeconds);
  }
  else
  {
    double oldLb = stats.lowerBound;
    const std::set<int> initialPrimalArcIndices;
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, initialPrimalArcIndices, duals, false, params.timeoutSeconds);
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

  //routeDD.print();
  return true;
};

void VRPTWDDSolver::addRouteToPrimalRoutes(std::vector<int> route)
{
  // ensure it's a real route to add it
  if (route.size() <= 2)
  {
    return;
  }
  else
  {
    if ((route[0] != 0) || (route[-1] != 0))
    {
      return;
    }
  }

  if (routeDD.isRouteFeasible(route))
  {
    // Try to improve with intra-route swaps
    double routeCost = vrptw.evaluateRouteDistance(route);

    if (params.insertionAblation != PrimalHeuristicInsertionAblation::INTRA_SWAP)
    {
      bool improved = true;
      while (improved)
      {
        improved = routeDD.intraRouteSwaps(route, routeCost);
      }
    }

    // Add route
    if (std::find(primalRoutes.begin(), primalRoutes.end(), route) == primalRoutes.end())
    {
      primalRoutes.push_back(route);
      primalRouteCosts.push_back(routeCost);
      int primalRouteIndex = primalRoutes.size()-1;
      for (int loc : route)
      {
        primalRouteLocationIndices[loc].push_back(primalRouteIndex);
      }

      std::cout << "added route to primal routes: ";
      for (int loc : route)
      {
        std::cout << loc << ",";
      }
      std::cout << std::endl;
    }
  }
}

// primal heuristic with MIP
bool VRPTWDDSolver::primalHeuristicMIP(FlowType flowType, std::vector<std::vector<int>>& returnRoutes, std::set<int>& returnRouteIndices)
{
  ++stats.numHeuristicIPs;

  // setup model
  IloEnv env;
  IloModel columnModel(env);
  IloRangeArray coverConstraints(env);
  IloRangeArray fixedPathConstraint(env);
  IloNumVarArray x(env, primalRoutes.size());
  IloNumArray initialPrimal(env);
  IloExpr objective(env);

  // setup variables and objective
  int varNum = 0;
  for (auto route : primalRoutes)
  {
    if (flowType == FlowType::IP)
    {
      x[varNum] = IloNumVar(env, 0, 1, ILOINT);
    }
    else
    {
      x[varNum] = IloNumVar(env, 0, 1);
    }
    objective += x[varNum] * primalRouteCosts[varNum];

    if (returnRouteIndices.find(varNum) != returnRouteIndices.end())
    {
      initialPrimal.add(1);
    }
    else
    {
      initialPrimal.add(0);
    }

    ++varNum;
  }

  // setup cover constraints
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    IloExpr sumLocationRoutes(env);

    for (int routeIndex : primalRouteLocationIndices[location])
    {
      sumLocationRoutes += x[routeIndex];
    }

    coverConstraints.add(sumLocationRoutes >= 1);
  }
  columnModel.add(coverConstraints);

  // set fixed number of vehicles
  if (vrptw.fixedNumPaths == FIXED_NUM_PATHS)
  {
    IloExpr fixedNumPaths(env);
    int routeIndex = 0;
    for (auto route : primalRoutes)
    {
      fixedNumPaths += x[routeIndex];
      ++routeIndex;
    }

    fixedPathConstraint.add(fixedNumPaths == vrptw.numVehicles);
    columnModel.add(fixedPathConstraint);
  }

  std::cout << "primal heuristic num vars: " << varNum << std::endl;
  columnModel.add(IloMinimize(env, objective));
  objective.end();

  // setup solver
  IloCplex solver(columnModel);
  solver.setOut(env.getNullStream());
  solver.setWarning(env.getNullStream());
  solver.setError(env.getNullStream());
  solver.setParam(IloCplex::Param::Threads, 1);
  solver.setParam(IloCplex::Param::TimeLimit, 5);
  //solver.exportModel("MIPHeuristicModel.lp");

  // Warm start
  if (!returnRouteIndices.empty())
  {
    solver.addMIPStart(x, initialPrimal);
  }

  // solve model
  solver.solve();

  // get results
  IloAlgorithm::Status solverStatus = solver.getStatus();
  if ((solverStatus == IloAlgorithm::Optimal) || ((solverStatus == IloAlgorithm::Feasible) && (!returnRouteIndices.empty())))
  {
    // store results
    returnRoutes.clear();
    returnRouteIndices.clear();
    for (int routeIndex=0; routeIndex<primalRoutes.size(); ++routeIndex)
    {
      double value = solver.getValue(x[routeIndex]);
      if (value >= 0.999)
      {
        returnRoutes.push_back(primalRoutes[routeIndex]);
        returnRouteIndices.insert(routeIndex);
      }
    }

    double objectiveValue = solver.getObjValue();
    std::cout << "primal heuristic MIP obj: " << objectiveValue << std::endl;

/*
    std::cout << "feasible primal routes: " << std::endl;
    for (auto primalRoute : returnRoutes)
    {
      std::cout << "route: " << std::endl;
      for (int loc : primalRoute)
      {
        std::cout << loc << " ";
      }
      std::cout << std::endl;
    }
*/

    return true;
  }

  return false;
}

void VRPTWDDSolver::destroyBySingleRoute(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& destroyedElements, int numElementsDestroy)
{
  bool routeSizeExists = false;
  while (!routeSizeExists)
  {
    for (auto route : feasibleSolution)
    {
      if (route.size() == numElementsDestroy)
      {
        routeSizeExists = true;
      }
    }
    if (!routeSizeExists)
    {
      numElementsDestroy + 2;
    }
  }

  while (destroyedElements.size() <= numElementsDestroy)
  {
    int currNumDestroyed = destroyedElements.size();
    int randomIntegerDestroy = std::rand();
    int randomIndexDestroy = randomIntegerDestroy % feasibleSolution.size();
    auto route = feasibleSolution[randomIndexDestroy];
    if (route.size() == numElementsDestroy)
    {
      for (int location : route)
      {
        if (location != 0)
        {
          destroyedElements.insert(location);
        }
      }
      break;
    }
  }
}

void VRPTWDDSolver::chooseRandomLocationsFromRoutes(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& randomElementsChosen)
{
  std::vector<int> locationsRemaining;
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    if (randomElementsChosen.find(location) == randomElementsChosen.end())
    {
      locationsRemaining.push_back(location);
    }
  }

  int randomIndex = std::rand() % static_cast<int>(locationsRemaining.size());
  int locationToDestroy = locationsRemaining[randomIndex];
  randomElementsChosen.insert(locationToDestroy);
  if (!vrptw.reliances.empty() && !vrptw.precedences.empty())
  {
    for (int relianceLocation : vrptw.reliances[locationToDestroy])
    {
      randomElementsChosen.insert(relianceLocation);
    }
    for (int precedenceLocation : vrptw.precedences[locationToDestroy])
    {
      randomElementsChosen.insert(precedenceLocation);
    }
  }
}

void VRPTWDDSolver::destroyRandomly(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& destroyedElements, int numElementsDestroy)
{
  int limit = std::min(numElementsDestroy,vrptw.numLocations-2);
  while (static_cast<int>(destroyedElements.size()) <= limit)
  {
    chooseRandomLocationsFromRoutes(feasibleSolution, destroyedElements);
  }
}

double VRPTWDDSolver::computeShawRelatedness(const std::vector<int>& route1, int index1, const std::vector<int>& route2, int index2)
{
  int location1 = route1[index1];
  int location2 = route2[index2];

  // Get pickup/delivery for PDPTW
  int pickup1 = -1;
  int pickup2 = -1;
  int delivery1 = -1;
  int delivery2 = -1;
  if (!vrptw.reliances.empty() and !vrptw.precedences.empty())
  {
    if (vrptw.reliances[location1].empty())
    {
      pickup1 = location1;
      delivery1 = *vrptw.reliances[location1].begin();
    }
    else
    {
      delivery1 = location1;
      pickup1 = *vrptw.precedences[location1].begin();
    }

    if (vrptw.reliances[location1].empty())
    {
      pickup2 = location2;
      delivery2 = *vrptw.reliances[location2].begin();
    }
    else
    {
      delivery2 = location2;
      pickup2 = *vrptw.precedences[location2].begin();
    }
  }
  else
  {
    pickup1 = location1;
    pickup2 = location2;
  }

  // Calculate metrics
  double distance = vrptw.distances[pickup1][pickup2];
  double loadDifference = std::abs(vrptw.demands[pickup1]  - vrptw.demands[pickup2]);

  double startTime1 = vrptw.calculateLocationPickupTime(route1, index1);
  double startTime2 = vrptw.calculateLocationPickupTime(route2, index2);
  double timeDifference = std::abs(startTime2 - startTime1);

  if (!vrptw.reliances.empty() and !vrptw.precedences.empty())
  {
    distance += vrptw.distances[delivery1][delivery2];

    // Get indices of deliveries
    int deliveryIndex1 = std::distance(route1.begin(), std::find(route1.begin(), route1.end(), delivery1));
    int deliveryIndex2 = std::distance(route2.begin(), std::find(route2.begin(), route2.end(), delivery2));
    double deliveryStartTime1 = vrptw.calculateLocationPickupTime(route1, deliveryIndex1);
    double deliveryStartTime2 = vrptw.calculateLocationPickupTime(route2, deliveryIndex2);
    timeDifference += std::abs(deliveryStartTime2 - deliveryStartTime1);
  }

  double normalizedDistanceValue = distance / (2*vrptw.maxDistance);
  double normalizedTimeValue = timeDifference / vrptw.maxStartTime;
  double normalizedLoadValue = loadDifference / (2*vrptw.maxDemand);

  double shawValue = normalizedDistanceValue + normalizedTimeValue + normalizedLoadValue;
  return shawValue;
}

void VRPTWDDSolver::destroyByShaw(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& destroyedElements, int numElementsDestroy)
{
  while (destroyedElements.empty())
  {
    chooseRandomLocationsFromRoutes(feasibleSolution, destroyedElements);
  }

  while (destroyedElements.size() <= numElementsDestroy)
  {
    int randomDestroyedLocationIndex = std::rand() % static_cast<int>(destroyedElements.size());
    int randomDestroyedLocation = *std::next(destroyedElements.begin(), randomDestroyedLocationIndex);
    int randomDestroyedLocationRouteIndex = -1;
    int randomDestroyedLocationRouteIndexIndex = -1;
    for (int routeIndex=0; routeIndex<feasibleSolution.size(); ++routeIndex)
    {
      auto route = feasibleSolution[routeIndex];
      auto it = std::find(route.begin(), route.end(), randomDestroyedLocation);
      if (it != route.end())
      {
        randomDestroyedLocationRouteIndex = routeIndex;
        randomDestroyedLocationRouteIndexIndex = std::distance(route.begin(), it);
        break;
      }
    }
    int mostRelatedLocation = -1;
    double mostRelatedScore = INF;
    for (int candidateLocation=1; candidateLocation<vrptw.numLocations; ++candidateLocation)
    for (int routeIndex=0; routeIndex<feasibleSolution.size(); ++routeIndex)
    {
      auto candidateLocationRoute = feasibleSolution[routeIndex];
      for (int candidateLocationIndex=1; candidateLocationIndex<candidateLocationRoute.size()-1; ++candidateLocationIndex)
      {
        int candidateLocation = candidateLocationRoute[candidateLocationIndex];
        if (destroyedElements.find(candidateLocation) == destroyedElements.end())
        {
          double shawValue = computeShawRelatedness(candidateLocationRoute, candidateLocationIndex, feasibleSolution[randomDestroyedLocationRouteIndex], randomDestroyedLocationRouteIndexIndex);
          if (shawValue < mostRelatedScore)
          {
            mostRelatedScore = shawValue;
            mostRelatedLocation = candidateLocation;
          }
        }
      }
    }

    destroyedElements.insert(mostRelatedLocation);
    if (!vrptw.reliances.empty() && !vrptw.precedences.empty())
    {
      for (int relianceLocation : vrptw.reliances[mostRelatedLocation])
      {
        destroyedElements.insert(relianceLocation);
      }
      for (int precedenceLocation : vrptw.precedences[mostRelatedLocation])
      {
        destroyedElements.insert(precedenceLocation);
      }
    }
  }
}

void VRPTWDDSolver::destroyByWorst(const std::vector<std::vector<int>>& feasibleSolution, std::set<int>& destroyedElements, int numElementsDestroy)
{
  while (destroyedElements.size() <= numElementsDestroy)
  {
    int worstLocation = -1;
    double worstDistance = -1;
    for (int routeIndex=0; routeIndex<feasibleSolution.size(); ++routeIndex)
    {
      auto route = feasibleSolution[routeIndex];
      int currLocation = 0;
      for (int routeIndexIndex=1; routeIndexIndex<route.size()-1; ++routeIndexIndex)
      {
        int routeLocation = route[routeIndexIndex];
        int nextRouteLocation = route[routeIndexIndex+1];
        if (destroyedElements.find(routeLocation) != destroyedElements.end())
        {
          continue;
        }

        double distance = vrptw.distances[currLocation][routeLocation] + vrptw.distances[routeLocation][nextRouteLocation];
        if (distance > worstDistance)
        {
          worstDistance = distance;
          worstLocation = routeLocation;
        }
        currLocation = routeLocation;
      }
    }

    destroyedElements.insert(worstLocation);
    if (!vrptw.reliances.empty() && !vrptw.precedences.empty())
    {
      for (int relianceLocation : vrptw.reliances[worstLocation])
      {
        destroyedElements.insert(relianceLocation);
      }
      for (int precedenceLocation : vrptw.precedences[worstLocation])
      {
        destroyedElements.insert(precedenceLocation);
      }
    }
  }
}

void VRPTWDDSolver::removeElements(std::set<int>& destroyedElements, const std::vector<std::vector<int>>& currSolution, int numElementsDestroy, int extraElementsDestroy, int sequenceIndex)
{
  // Use params PrimalHeuristicRemovalStrategy
  if (params.removalStrategy == PrimalHeuristicRemovalStrategy::RANDOM)
  {
    destroyRandomly(currSolution, destroyedElements, numElementsDestroy);
  }
  else if (params.removalStrategy == PrimalHeuristicRemovalStrategy::SHAW)
  {
    destroyByShaw(currSolution, destroyedElements, numElementsDestroy);
  }
  else if (params.removalStrategy == PrimalHeuristicRemovalStrategy::WORST)
  {
    destroyByWorst(currSolution, destroyedElements, numElementsDestroy);
  }
  else
  {
    auto route = currSolution[sequenceIndex];
    if ((static_cast<int>(route.size()) - 2) == numElementsDestroy)
    {
      //std::cout << "destroy by single route index: " << sequenceIndex << std::endl;
      for (int loc : route)
      {
        if (loc != 0)
        {
          destroyedElements.insert(loc);
        }
      }

      if (params.removalStrategy == PrimalHeuristicRemovalStrategy::SEQUENCE_RANDOM)
      {
        destroyRandomly(currSolution, destroyedElements, extraElementsDestroy);
      }
      else if (params.removalStrategy == PrimalHeuristicRemovalStrategy::SEQUENCE_SHAW)
      {
        destroyByShaw(currSolution, destroyedElements, extraElementsDestroy);
      }
      else if (params.removalStrategy == PrimalHeuristicRemovalStrategy::SEQUENCE_WORST)
      {
        destroyByWorst(currSolution, destroyedElements, extraElementsDestroy);
      }
    }
  }
}

// PDPTW needs to store second obj that isn't the real obj but the MIP with fake weights.
bool VRPTWDDSolver::largeNeighborhoodSearch(const std::vector<std::vector<int>>& feasibleSolution, const Dual& dual, int timeoutSeconds)
{
  // Stats setup
  ++stats.numHeuristicLNSs;
  bool isImproved = false;
  std::vector<std::vector<int>> currSolution = feasibleSolution;
  double currUpperBound = vrptw.evaluateSolutionCost(feasibleSolution);

  // Parameter setup
  int numElementsDestroy = 2;
  int extraElementsDestroy = 1;
  int sequenceIndex = 0;

  int maxSequenceSize = 0;
  for (auto sequence : feasibleSolution)
  {
    maxSequenceSize = std::max(static_cast<int>(sequence.size()), maxSequenceSize);
  }

  // Run LNS in loop for timeout
  bool shouldContinue = true;
  auto startTime = std::chrono::high_resolution_clock::now();
  std::set<std::set<int>> destroyedSetsTried;
  while (shouldContinue)
  {
    // Setup timeout
    auto currTime = std::chrono::high_resolution_clock::now();
    auto numSeconds = std::chrono::duration_cast<std::chrono::seconds>(currTime - startTime).count();
    int remainingSeconds = timeoutSeconds - numSeconds;
    if (remainingSeconds < 0)
    {
      shouldContinue = false;
      break;
    }

    // Removal
    std::set<int> destroyedElements;
    removeElements(destroyedElements, currSolution, numElementsDestroy, extraElementsDestroy, sequenceIndex);

    // Repair with LNS
    std::vector<std::vector<int>> newBestRoutes;
    double lnsHeuristicUpperBound = INF; 
    if (destroyedElements.size() > 0)
    {
      VRPTWDecisionDiagram heuristicDD(vrptw, params);
      double mipObjective = heuristicDD.repairSolution(currSolution, destroyedElements, dual, newBestRoutes, remainingSeconds);
      ++stats.numPrimalLNSRepairs;
      if (mipObjective > 0)
      {
        lnsHeuristicUpperBound = vrptw.evaluateSolutionCost(newBestRoutes);
      }
    }

    // Handle solution
    if (lnsHeuristicUpperBound < currUpperBound)
    {
      isImproved = true;
      std::cout << "improved iteration ub from lns primal heuristic: " << lnsHeuristicUpperBound << std::endl;
      std::cout << "improved primal routes: " << std::endl;
      for (auto primalRoute : newBestRoutes)
      {
        std::cout << "route: " << std::endl;
        for (int loc : primalRoute)
        {
          std::cout << loc << " ";
        }
        std::cout << std::endl;
      }

      // Add to primal routes
      int previousNumRoutes = static_cast<int>(primalRoutes.size());
      for (auto route : newBestRoutes)
      {
        addRouteToPrimalRoutes(route);
      }
      int afterNumRoutes = static_cast<int>(primalRoutes.size());

      // Update current stats
      currSolution = newBestRoutes;
      currUpperBound = lnsHeuristicUpperBound;

      if (lnsHeuristicUpperBound < stats.upperBound)
      {
        stats.upperBound = std::min(stats.upperBound, lnsHeuristicUpperBound);
        stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
      }

      // Run MIP Pool if any added - could help
      if (previousNumRoutes < afterNumRoutes)
      {
        std::cout << "routes added, try running MIP Pool" << std::endl;
        std::vector<std::vector<int>> mipPoolSolution;
        bool primalHeuristicFeasible = primalHeuristicMIP(FlowType::IP, mipPoolSolution, mipPoolSolutionIndices);
        if (primalHeuristicFeasible)
        {
          double mipPoolUpperBound = vrptw.evaluateSolutionCost(mipPoolSolution);
          if (mipPoolUpperBound < stats.upperBound)
          {
            std::cout << "improved ub from MIP pool: " << mipPoolUpperBound << std::endl;
            currSolution = mipPoolSolution;
            stats.upperBound = mipPoolUpperBound;
            stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
            startTime = std::chrono::high_resolution_clock::now();
          }
        }
      }

      // Update parameters
      startTime = std::chrono::high_resolution_clock::now();
      destroyedSetsTried.clear();
      numElementsDestroy = 2;
      extraElementsDestroy = 1;
      sequenceIndex = 0;
    }
    else
    {
      // Not improved, make updates
      if ((params.removalStrategy == PrimalHeuristicRemovalStrategy::SEQUENCE_RANDOM) || (params.removalStrategy == PrimalHeuristicRemovalStrategy::SEQUENCE_SHAW) || (params.removalStrategy == PrimalHeuristicRemovalStrategy::SEQUENCE_WORST))
      {
        ++sequenceIndex;
        if (sequenceIndex == currSolution.size())
        {
          sequenceIndex = 0;
          ++numElementsDestroy;
        }

        if (numElementsDestroy > maxSequenceSize)
        {
          ++extraElementsDestroy;
        }
      }
      else
      {
        numElementsDestroy = numElementsDestroy + 1;
      }
 
      if (numElementsDestroy == vrptw.numLocations - 2)
      {
        shouldContinue = false;
        break;
      }
    }
  }

  return isImproved;
}

void VRPTWDDSolver::repairMultipliers(Dual& repairedDual, LPSolveType solveType)
{
  auto startTimeRepair = std::chrono::high_resolution_clock::now();
  while (true)
  {
    routeDD.setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(repairedDual, solveType);

    std::vector<int> treeByParentArcs;
    treeByParentArcs.resize(routeDD.getNodes().size());
    std::vector<int> shortestPathByArc;
    double longestShortestPathLength = 0;
    double shortestPathLength = routeDD.computeShortestPathBFSWang(treeByParentArcs, shortestPathByArc, longestShortestPathLength);
    shortestPathLength = shortestPathLength - repairedDual.fixedPathDual;
    //std::cout << "spl: " << shortestPathLength << std::endl;
    if (shortestPathLength >= 0.0000000001)
    {
      break;
    }
    else
    {
      std::cout << "Dual needs repair, spl: " << shortestPathLength << std::endl;
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
/*
      std::cout << "Dual needs repair, spl: " << shortestPathLength << std::endl;
      bool updated = false;
      std::vector<std::vector<int>> shortestPathsByArc;
      shortestPathsByArc.push_back(shortestPathByArc);
      std::set<int> locations;
      routeDD.getLocationsOnArcPaths(shortestPathsByArc, locations);
      double updateAmount = std::min(shortestPathLength, -1.0);

      if (vrptw.fixedNumPaths == FixedNumPaths::FIXED_NUM_PATHS)
      {
        repairedDual.fixedPathDual = repairedDual.fixedPathDual + shortestPathLength - 0.001;
      }

      for (int index=1; index<vrptw.numLocations; ++index)
      {
        repairedDual.lambda[index] = std::max(0.0, repairedDual.lambda[index] + (updateAmount / std::max(1,((int)shortestPathByArc.size() - 1))) - 1);
      }
*/
    }
  }

  //printMultipliers(repairedDual);

  auto endTimeRepair = std::chrono::high_resolution_clock::now();
  auto totalTimeRepair = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeRepair - startTimeRepair).count();
  stats.millisecondsRepairingLAG += totalTimeRepair;
};

void VRPTWDDSolver::resizeMultipliers(const Dual& dual1, std::vector<Dual>& duals2)
{
  for (int index=0; index<duals2.size(); ++index)
  {
    Dual& dual2 = duals2[index];
    dual2.lambda.resize(dual1.lambda.size());
    dual2.capDuals.resize(dual1.capDuals.size());
    dual2.combDuals.resize(dual1.combDuals.size());
    dual2.srcDuals.resize(dual1.srcDuals.size());
  }
}

void VRPTWDDSolver::resizeMultipliers(const Dual& dual1, std::deque<Dual>& duals2)
{
  for (int index=0; index<duals2.size(); ++index)
  {
    Dual& dual2 = duals2[index];
    dual2.lambda.resize(dual1.lambda.size());
    dual2.capDuals.resize(dual1.capDuals.size());
    dual2.combDuals.resize(dual1.combDuals.size());
    dual2.srcDuals.resize(dual1.srcDuals.size());
  }
}

void VRPTWDDSolver::resizeMultipliers(const Dual& dual1, Dual& dual2)
{
  dual2.lambda.resize(dual1.lambda.size());
  dual2.capDuals.resize(dual1.capDuals.size());
  dual2.combDuals.resize(dual1.combDuals.size());
  dual2.srcDuals.resize(dual1.srcDuals.size());
}

// for Lagrangean method, the choice of initial labmda is important
// choice 1: simple all 0
// choice 2: optimal value will be min distance from depot
// choice 3: optimal value will be distances from depot
//std::vector<double> lambda(vrptw.numLocations, 0);
//std::vector<double> lambda(vrptw.numLocations, *std::min_element(vrptw.distances[0].begin()+1, vrptw.distances[0].end()));
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
      else if (stats.upperBound < INF)
      {
        double splitUpperBound = stats.upperBound * 1.0 / vrptw.numLocations;
        dual.lambda[location] = splitUpperBound;
      }
      else
      {
        dual.lambda[location] = 0;
      }
    }
    else
    {
      if (vrptw.problemType == ProblemType::PDP)
      {
        dual.lambda[location] = vrptw.distances[0][location] * std::abs(vrptw.demands[location]) / vrptw.capacity;
      }
      else
      {
        dual.lambda[location] = 2 * vrptw.distances[0][location] * std::abs(vrptw.demands[location]) / vrptw.capacity;
      }
    }
  }
}


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
    std::cout << dual.capDuals[index] << ",";
  }
  std::cout << std::endl;
 
  std::cout << "combDuals: ";
  for (int index=0; index<dual.combDuals.size(); ++index)
  {
    std::cout << dual.combDuals[index] << ",";
  }
  std::cout << std::endl;

  std::cout << "srcDuals: ";
  for (int index=0; index<dual.srcDuals.size(); ++index)
  {
    std::cout << dual.srcDuals[index] << ",";
  }
  std::cout << std::endl;

  std::cout << "fixed dual: " << dual.fixedPathDual << std::endl;
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

  std::vector<double> srcCutValues;
  routeDD.getSrcCutValues(srcCutValues);

  // gamma_(k) = b - Ax_(k)
  // Beasley - when multiplier is already 0 and step direction is negative, don't include in ||gamma||^2
  double normGammaSquared = 0;
  std::vector<double> gamma(vrptw.numLocations + dual.capDuals.size() + dual.srcDuals.size() + dual.combDuals.size(), 0);
  for (int i=1; i<vrptw.numLocations; ++i)
  {
    gamma[i] = 1 - locationsCovered[i];
    //std::cout << "gamma: " << i << " val: " << gamma[i] << std::endl;
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
    gamma[gammaIndex] = (-1 * routeDD.getCapCutSetRHS(i)) + cutValues[i];
    if ((dual.capDuals[i] <= 0.001) && (gamma[gammaIndex] <= 0.001))
    {
      gamma[gammaIndex] = 0;
    }
    else
    {
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
    //std::cout << "rcc gamma: " << gammaIndex << " val: " << gamma[gammaIndex] << std::endl;
  }

  // src cuts are in the form Cx <= r so chang eto -Cx >= -r
  for (int i=0; i<dual.srcDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations + dual.capDuals.size();
    SRCType srcType = routeDD.getSRCType(i);
    int rhs = routeDD.getSRCRHS(srcType);
    gamma[gammaIndex] = (-1*rhs) + srcCutValues[i];
    //std::cout << "src gamma: " << gammaIndex << " val: " << gamma[gammaIndex] << std::endl;

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
  double lbGoal = std::max(lagrangeanLowerBound * (1 + eta), lagrangeanLowerBound * (1 - eta));
  double psiStar = std::min(stats.upperBound, lbGoal);

  // alpha_(k) = (psi_(star) - psi(lambda(k))) / ||gamma_(k)||_(2)^2
  double alpha = (psiStar - lagrangeanLowerBound) / normGammaSquared;
 
  std::cout << "alpha: " << alpha << std::endl;
  std::cout << "psiStar: " << psiStar << std::endl;
  std::cout << "lagLB: " << lagrangeanLowerBound << std::endl;
  std::cout << "||gamma||^2: " << normGammaSquared << std::endl;

  //printMultipliers(dual);
  //printMultipliers(bestDual);

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
};

// SGD Algorithm or
// Volume Algorithm (Barahona et al.)
// Always move from pi_{bar}, the best dual solution found so far
// Use subgradient 1 - Ax_{bar}
// Use red, yellow, green iterations that update step size dynamically
bool VRPTWDDSolver::solveLagrangeanRelaxation(Dual& dual, SGDAlgorithm& sgdAlgo)
{
  // initialize values
  bool shouldTerminate = false;
  int lastMuImprovedIteration = 0;
  double muLowerBound = 0.0;
  double currIterLowerBound = 0.0;
  int numLagIterations = 0;
  primal = Primal();
  alphaLowerBound = 0.1;
  alphaLowerBoundIteration = 0;

  // warmstart new iterations
  dual = bestDual;

  // for cuts
  std::vector<std::vector<int>> routesForCuts;
 
  // run algorithm
  stats.lpIterations = stats.lpIterations + 1;
  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
  while (!shouldTerminate)
  {
    std::vector<std::vector<int>> infeasibleRoutes;
    std::vector<std::vector<int>> infeasibleRoutesByArc;
    while (!shouldTerminate && (infeasibleRoutesByArc.size() < infeasibleRoutesBatchSize))
    {
      stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
      ++stats.numLagIterations;
      ++numLagIterations;
      ++stats.numLagIterationsWithResets;

      // Arc fixing
      bool isDualFeasible = false;
      bool tryArcFixingOrRepair = false;
      double longestShortestPathLength = 0;
      double minReducedCost = 0.0;
      Dual repairedDual(dual);
      double percentFixed = 0.0;
      if (params.useVariableFixing && (bestDualsArcFixingLagPercent > 0.001))
      {
        auto startTimeFix = std::chrono::high_resolution_clock::now();
        std::cout << "fix with best duals from LAG. Using: " << bestDualsArcFixingLag.size() << " duals" << std::endl;
        std::vector<Dual> dualsToUseForFixing;
        for (int dualMultiplierIndex=0; dualMultiplierIndex<bestDualsArcFixingLag.size(); ++dualMultiplierIndex)
        {
          Dual dualToRepair = bestDualsArcFixingLag[dualMultiplierIndex];
          repairMultipliers(dualToRepair, LPSolveType::LAGSolver);
          //bestDualsArcFixingLag[dualMultiplierIndex] = dualToRepair;
          dualsToUseForFixing.push_back(dualToRepair);
        }
        routeDD.fixArcs(dualsToUseForFixing, LPSolveType::LAGSolver, stats.upperBound);
        auto endTimeFix = std::chrono::high_resolution_clock::now();
        auto totalTimeFix = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeFix - startTimeFix).count();
        stats.millisecondsFix += totalTimeFix;
      }

      // Run muSSP
      std::vector<std::vector<int>> shortestPaths;
      auto startSSPTime = std::chrono::high_resolution_clock::now();
      double lagrangeanLowerBound = routeDD.solveMinCostFlowModelWang(dual, shortestPaths, isDualFeasible, minReducedCost, longestShortestPathLength);

/*
      const std::set<int> initialPrimalArcIndices;
      Dual newDual = dual;
      routeDD.setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(newDual, LPSolveType::LAGSolver);
      double bestBound = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::N, UseColumnGeneration::NO_CG, initialPrimalArcIndices, newDual, false, params.timeoutSeconds);
      std::cout << "best bound LP: " << bestBound << std::endl;
      if ((bestBound >= lagrangeanLowerBound + 0.001) || (bestBound <= lagrangeanLowerBound - 0.001))
      {
        std::cout << "ERROR: " << bestBound << " vs. " << lagrangeanLowerBound << std::endl;
        routeDD.print();
        return false;
      }
*/

      stats.numSSPIterations = stats.numSSPIterations + shortestPaths.size();
      auto endMuSSPTime = std::chrono::high_resolution_clock::now();
      auto sspSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endMuSSPTime - startSSPTime).count();
      stats.millisecondsSolvingSSP = stats.millisecondsSolvingSSP + sspSolveTime;
 
      // Heuristic for trying arc fixing or repair
      if (stats.lowerBound + longestShortestPathLength > stats.upperBound)
      {
        if (stats.numLagIterations % 10 == 0)
        {
          tryArcFixingOrRepair = true;
        }
      }

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
        repairedDual.fixedPathDual = dual.fixedPathDual;

        minReducedCost = minReducedCost - dual.fixedPathDual;
        if (minReducedCost >= 0)
        {
          isDualFeasible = true;
        }
        else
        {
          isDualFeasible = false;
        }
      }
      else
      {
        dual.fixedPathDual = 0.0;
        repairedDual.fixedPathDual = dual.fixedPathDual;
      }

      std::cout << "Lagrangean bound solution value: " << lagrangeanLowerBound << std::endl;

      // Keep track of best-known solution
      if (currIterLowerBound < lagrangeanLowerBound)
      {
        if (((1 + params.percentImprovementThreshold) * muLowerBound) < lagrangeanLowerBound)
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
        //printMultipliers(dual);

        bestDual = dual;
        bestDualValue = lagrangeanLowerBound;
        std::cout << "new best dual, value: " << bestDualValue << std::endl;
        /*
        double percentGap = (stats.upperBound - bestDualValue) * 100.0 / stats.upperBound;
        if (percentGap < 10)
        {
          double limitToMerge = (stats.upperBound - bestDualValue) / 2.0;
          routeDD.findMergeNodesReducedCost(bestDual, LPSolveType::LAGSolver, limitToMerge);
          stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
        }
        */

        stepSizeMultiplierIteration = 0;
        if (lagrangeanLowerBound > 0.95 * targetLowerBound)
        {
          targetLowerBound = std::min(lagrangeanLowerBound * targetBoundIncrease, stats.upperBound);
        }

        alphaLowerBoundIteration = 0;
      }
      else
      {
        if (stepSizeMultiplierIteration == stepSizeMultiplierIterationCutoff)
        {
          stepSizeMultiplier = std::max(minStepSizeMultiplier, stepSizeMultiplier * 0.66);
          stepSizeMultiplierIteration = 0;
        }
        else
        {
          ++stepSizeMultiplierIteration;
        }
 
        ++alphaLowerBoundIteration;
        if (alphaLowerBoundIteration == 20)
        {
          alphaLowerBoundIteration = 0;
          alphaLowerBound = std::max(alphaLowerBound / 2.0, 0.00001);
        }
      }

      // Compute primal solution
      auto startTimeDecompose = std::chrono::high_resolution_clock::now();
      std::vector<int> infeasibleRoute;
      std::vector<double> routeFlows;
      std::vector<std::vector<int>> decomposedRoutes;
      std::vector<std::vector<int>> decomposedArcs;
      routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, DecompositionReason::DECOMPOSE);
      routesForCuts = decomposedRoutes;
      auto endTimeDecompose = std::chrono::high_resolution_clock::now();
      auto totalTimeDecompose = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeDecompose - startTimeDecompose).count();
      stats.millisecondsDecompose += totalTimeDecompose;

      // Add to set for primal heuristic
      int numPrimalRoutes = primalRoutes.size();
      if ((params.primalHeuristic != PrimalHeuristic::BEST_KNOWN) && (params.primalHeuristic != PrimalHeuristic::STANDALONE))
      {
        for (auto route : decomposedRoutes)
        {
          // include feasible routes, and truncated infeasible routes
          if (routeDD.isRouteFeasible(route))
          {
            addRouteToPrimalRoutes(route);
          }
          else
          {
            if (params.insertionAblation != PrimalHeuristicInsertionAblation::TRUNCATED)
            {
              std::vector<int> truncatedRoute;
              routeDD.createTruncatedRoute(route, truncatedRoute);
              addRouteToPrimalRoutes(truncatedRoute);
            }
          }
        }
      }

      // Primal Heuristic
      std::vector<std::vector<int>> routesByLocationPrimalHeuristic;
      bool primalHeuristicFeasible = false;
      if (params.primalHeuristic == PrimalHeuristic::CE_GREEDY)
      {
        //primalHeuristicFeasible = routeDD.primalHeuristicGreedy(routesByLocationPrimalHeuristic);
      }
      else if ((params.primalHeuristic == PrimalHeuristic::CE_MIP) || (params.primalHeuristic == PrimalHeuristic::CE_MIP_LNS))
      {
        primalHeuristicFeasible = primalHeuristicMIP(FlowType::IP, routesByLocationPrimalHeuristic, mipPoolSolutionIndices);
      }

      // Evaluate solution
      if (primalHeuristicFeasible)
      {
        double heuristicUpperBound = vrptw.evaluateSolutionCost(routesByLocationPrimalHeuristic);
        if (heuristicUpperBound < stats.upperBound)
        {
          stats.upperBound = heuristicUpperBound;
          std::cout << "improved ub from MIP primal heuristic: " << heuristicUpperBound << std::endl;
          stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());

          // Run LNS
          if (params.primalHeuristic == PrimalHeuristic::CE_MIP_LNS)
          {
            largeNeighborhoodSearch(routesByLocationPrimalHeuristic, dual, params.lnsTimeoutSeconds);
          }
        }
      }

      // for yellow, check v^{t} dot (1 - Ax^{t})
      if (sgdAlgo == SGDAlgorithm::VA)
      {
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
        auto startTimeTryingAlpha = std::chrono::high_resolution_clock::now();
        double alphaTry = alphaLowerBound / 10;
        double bestAlpha = alphaTry;
        double bestAlphaValue = INF;
        if (!primal.xDecompositionFlows.empty())
        {
          bestAlpha = (alphaLowerBound / 10 + alphaLowerBound) / 2;
        }
        else
        {
          bestAlpha = 1.0;
        }
        std::cout << "best alpha: " << bestAlpha << std::endl;

        Primal nextPrimal(primal);
        constructNextPrimal(bestAlpha, decomposedRoutes, decomposedArcs, nextPrimal);
        primal = nextPrimal;
        auto endTimeTryingAlpha = std::chrono::high_resolution_clock::now();
        auto totalTimeTryingAlpha = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeTryingAlpha - startTimeTryingAlpha).count();
        stats.millisecondsTryingAlpha += totalTimeTryingAlpha;
      }

      // Update dual
      auto startTimeUpdateDual = std::chrono::high_resolution_clock::now();
      if (sgdAlgo == SGDAlgorithm::VA)
      {
        updateMultipliersVolumeAlgorithm(dual, primal, stats.numLagIterations);
      }
      else
      {
        //updateMultipliers(dual, lagrangeanLowerBound, stats.numLagIterations);
        updateMultipliers(dual, lagrangeanLowerBound, numLagIterations);
      }
      auto endTimeUpdateDual = std::chrono::high_resolution_clock::now();
      auto totalTimeUpdateDual = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeUpdateDual - startTimeUpdateDual).count();
      stats.millisecondsUpdateDual += totalTimeUpdateDual;

      // Check for infeasibilities
      if (params.useSeparations)
      {
        bool infeasibleRouteFound = true;
        while (infeasibleRouteFound)
        {
          std::vector<int> infeasibleRouteByArcs;
          std::vector<double> flows;
          std::vector<std::vector<int>> routes;
          std::vector<std::vector<int>> separationDecomposedArcs;
          routeDD.decomposeRoutes(infeasibleRouteByArcs, flows, routes, separationDecomposedArcs, DecompositionReason::SEPARATE);
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
      if ((numLagIterations - lastMuImprovedIteration) > params.lowImprovementIterationsToNextRound)
      {
        shouldTerminate = true;
      }

      if ((stats.getNumSeconds() >= params.timeoutSeconds) || (stats.lowerBound + 0.01 > stats.upperBound))
      {
        shouldTerminate = true;
      }

      // Repair Dual - separations can change the paths available
      if (tryArcFixingOrRepair)
      {
        repairMultipliers(repairedDual, LPSolveType::LAGSolver);
        double repairedBound = routeDD.getDualObjectiveValue(repairedDual, LPSolveType::LAGSolver);
        if (bestDualValue < repairedBound)
        {
          std::cout << "repaired lb: " << repairedBound << std::endl;
          stats.lowerBound = std::max(stats.lowerBound, repairedBound);
          stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
          //printMultipliers(repairedDual);

          if (repairedBound > 0.95 * targetLowerBound)
          {
            targetLowerBound = std::min(repairedBound * targetBoundIncrease, stats.upperBound);
          }
        }
      }

      // Arc Fixing
      if (params.useVariableFixing && tryArcFixingOrRepair)
      {
        percentFixed = routeDD.fixArcs(repairedDual, LPSolveType::LAGSolver, stats.upperBound);
        std::cout << "percent fixed: " << percentFixed << std::endl;
      }

      if (percentFixed > bestDualsArcFixingLagPercent)
      {
        std::cout << "updating best lambda arc fixing" << std::endl;
        bestDualsArcFixingLag.push_back(repairedDual);
        if (bestDualsArcFixingLag.size() > 5)
        {
          bestDualsArcFixingLag.pop_front();
        }
        bestDualsArcFixingLagPercent = percentFixed;
      }

      if (params.changeToLP)
      {
        if ((percentFixed > percentFixedToChangeToCPLEX) || (routeDD.getNumArcsNotRemovedOrReverseOrFixed() < numArcsToChangeToCPLEX))
        {
          params.lpSolveType = LPSolveType::LPSolver;
          stats.lpIterations = 1;
          shouldTerminate = true;
          std::cout << "switching from LAG to LP solver" << std::endl;
        }
      }
    }

    // Separations
    if (!shouldTerminate)
    {
      for (auto infeasibleRouteToSeparateByArc : infeasibleRoutesByArc)
      {
        // when separating changes, can use locations to re-find the route
        std::vector<int> emptyRoute;
        if (routeDD.doesRouteExistByArcs(infeasibleRouteToSeparateByArc, emptyRoute))
        {
          stats.numSeparations = stats.numSeparations + 1;
          routeDD.separateRoute(infeasibleRouteToSeparateByArc);
        }
      }
    }

    // Consider restarting
    bool shouldCut = true;
    if ((shouldTerminate && (params.useRCCs || params.useSRC3s || params.useSRC4s || params.useSRC5V1s || params.useSRC5V2s)) && shouldCut)
    {
      addCutsUsingCurrentPrimal(dual, routesForCuts);
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
      auto routeArcs = decomposedRouteArcs[index];
      nextPrimal.xDecompositions.push_back(route);
      nextPrimal.xDecompositionArcs.push_back(routeArcs);
      nextPrimal.xDecompositionFlows.push_back(alphaTry);
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

  std::vector<double> srcCutValues;
  routeDD.getSrcCutValuesRoutes(currPrimal, srcCutValues);

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
    gradient[gradientIndex] = (-1 * routeDD.getCapCutSetRHS(i)) + cutValues[i];
  }

  // src cuts are in the form Cx <= r so change to -Cx >= -r
  for (int i=0; i<dual.srcDuals.size(); ++i)
  {
    int gradientIndex = i + vrptw.numLocations + dual.capDuals.size();
    SRCType srcType = routeDD.getSRCType(i);
    int rhs = routeDD.getSRCRHS(srcType);
    gradient[gradientIndex] = (-1*rhs) + srcCutValues[i];
  }

  for (int i=0; i<dual.combDuals.size(); ++i)
  {
    int gradientIndex = i + vrptw.numLocations + dual.capDuals.size() + dual.srcDuals.size();
    gradient[gradientIndex] = routeDD.getCombCutRHS(i) - combValues[i];
  }
};

void VRPTWDDSolver::updateMultipliersVolumeAlgorithm(Dual& dual, Primal& currPrimal, int iteration)
{
  // get values of LHS for all dualized inequalities
  std::unordered_map<int,double> locationsCovered;
  routeDD.getNumberOfTimesLocationsCoveredRoutes(currPrimal, locationsCovered);

  std::vector<double> cutValues;
  routeDD.getCutSetValuesRoutes(currPrimal, cutValues);
 
  std::vector<double> combValues;
  routeDD.getCombValuesRoutes(currPrimal, combValues);

  std::vector<double> srcCutValues;
  routeDD.getSrcCutValuesRoutes(currPrimal, srcCutValues);

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
 
    std::cout << "gamma at lambda " << i << ":" << gamma[i] << std::endl;
  }

  // cap cuts are in the form Cx <= r so change to -Cx >= -r
  for (int i=0; i<dual.capDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations;

    gamma[gammaIndex] = (-1 * routeDD.getCapCutSetRHS(i)) + cutValues[i];
    if ((dual.capDuals[i] <= 0.001) && (gamma[gammaIndex] <= 0.001))
    {
      gamma[gammaIndex] = 0;
    }
    else
    {
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
  }

  // src cuts are in the form Cx <= r so change to -Cx >= -r
  for (int i=0; i<dual.srcDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations + dual.capDuals.size();
    SRCType srcType = routeDD.getSRCType(i);
    int rhs = routeDD.getSRCRHS(srcType);
    gamma[gammaIndex] = (-1*rhs) + srcCutValues[i];

    if ((bestDual.srcDuals[i] <= 0.01) && (gamma[gammaIndex] <= 0.01))
    {
      continue;
    }
    else
    {
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
 
    std::cout << "gamma at src " << i << ":" << gamma[gammaIndex] << std::endl;
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
 
    std::cout << "gamma at comb " << i << ":" << gamma[gammaIndex] << std::endl;
  }

  double eta = 0.05 * 100 / (100 + iteration);
  double psiStar = stats.lowerBound + (1 + eta);
  double alpha = stepSizeMultiplier * (psiStar - bestDualValue) / normGammaSquared;
  std::cout << "alpha: " << alpha << std::endl;
  std::cout << "bestDualValue: " << bestDualValue << std::endl;
  std::cout << "stepSizeMultiplier: " << stepSizeMultiplier << std::endl;
  std::cout << "targetLowerBound: " << psiStar << std::endl;
  std::cout << "||gamma||^2: " << normGammaSquared << std::endl;

  //printMultipliers(dual);
  //printMultipliers(bestDual);

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
};

void VRPTWDDSolver::separateSequencesAndTruncate(const std::vector<std::vector<int>>& inputSequences, std::vector<std::vector<int>>& outputSequences, std::vector<std::vector<int>>& outputSequencesArcs)
{
  // Separate and add to output sequences
  for (auto sequence : inputSequences)
  {
    if (routeDD.isRouteFeasible(sequence))
    {
      outputSequences.push_back(sequence);
    }
    else
    {
      std::vector<int> truncatedSequence;
      routeDD.createTruncatedRoute(sequence, truncatedSequence);
      outputSequences.push_back(truncatedSequence);
    }
  }
 
  std::cout << "sequences truncated and separated:" << std::endl;
  for (auto route : outputSequences)
  {
    for (int loc : route)
    {
      std::cout << loc << ",";
    }
    std::cout << " ";
  }
  std::cout << std::endl;

  for (auto sequence : outputSequences)
  {
    std::vector<int> sequenceArcs;
    routeDD.doesRouteExistByLocations(sequence, sequenceArcs);
    routeDD.separateRoute(sequenceArcs);
    outputSequencesArcs.push_back(sequenceArcs);
  }
}

bool VRPTWDDSolver::addCutsUsingCurrentPrimal(Dual& dual, const std::vector<std::vector<int>>& decomposedRoutes)
{
  std::cout << "Adding cuts for routes:" << std::endl;
  for (auto route : decomposedRoutes)
  {
    for (int loc : route)
    {
      std::cout << loc << ",";
    }
    std::cout << " ";
  }
  std::cout << std::endl;

  auto startTimeCut = std::chrono::high_resolution_clock::now();

  // Convert primal routes in case some infeasibilities were separated
  // Now we have all feasible or truncated routes
  std::vector<std::vector<int>> truncatedExactSequences;
  std::vector<std::vector<int>> truncatedExactSequencesArcs;
  separateSequencesAndTruncate(decomposedRoutes, truncatedExactSequences, truncatedExactSequencesArcs);

  std::cout << "separated sequences" << std::endl;
 
  // Use 0.99 for separation package to not see as integral solution
  Primal cutPrimal;
  cutPrimal.xDecompositions = truncatedExactSequences;
  cutPrimal.xDecompositionArcs = truncatedExactSequencesArcs;
  cutPrimal.xDecompositionFlows = std::vector<double>(truncatedExactSequences.size(),0.99);

  // Print the solution
  for (int index=0; index<cutPrimal.xDecompositionFlows.size(); ++index)
  {
    double flow = cutPrimal.xDecompositionFlows[index];
    auto route = cutPrimal.xDecompositions[index];
    std::cout << "route with flow " << flow << " : ";
    for (int loc : route)
    {
      std::cout << loc << ",";
    }
    std::cout << " ";
    auto routeArcs = cutPrimal.xDecompositionArcs[index];
    for (int arcIndex : routeArcs)
    {
      std::cout << arcIndex << ",";
    }
    std::cout << std::endl;
  }

  // Robust Cuts
  bool cutAdded = false;
  std::vector<int> edgeTail;
  std::vector<int> edgeHead;
  std::vector<double> edgeFlow;
  std::vector<int> rccArcs;
  std::vector<double> rccArcFlows;
  if (params.useRCCs)
  {
    // Convert to Sep format
    convertArcIndicesForVRPTWSep(cutPrimal, edgeTail, edgeHead, edgeFlow, rccArcs, rccArcFlows);

    // Find cuts based on these feasible and truncated routes
    // RCC need to use the package, so update that method, currently not sparse at all...
    // Much more sparse, helps with convergence of subgradient descent
    std::vector<int> allCutArcs;
    for (auto sequence : truncatedExactSequencesArcs)
    {
      for (int arcIndex : sequence)
      {
        allCutArcs.push_back(arcIndex);
      }
    }

    if (vrptw.problemType != ProblemType::PDP)
    {
      addRCCs(edgeTail, edgeHead, edgeFlow, allCutArcs, params.numLagCuts, cutAdded, dual);
    }

    // Strengthened Combs
    //addCombs(edgeTail, edgeHead, edgeFlow, cutAdded);
    //dual.combDuals.resize(routeDD.getNumCombCuts());
  }

  // Non-Robust Cuts
  if (params.useSRC3s || params.useSRC4s || params.useSRC5V1s || params.useSRC5V2s)
  {
    double totalFlow = 0;
    std::vector<double> flowByVertex(vrptw.numLocations, 0);
    for (int index=0; index<cutPrimal.xDecompositionFlows.size(); ++index)
    {
      double flow = cutPrimal.xDecompositionFlows[index];
      totalFlow += flow;
      auto route = cutPrimal.xDecompositions[index];
      for (int loc : route)
      {
        flowByVertex[loc] += flow;
      }
    }
    std::cout << "total flow: " << totalFlow << std::endl;
    for (int index=0; index<vrptw.numLocations; ++index)
    {
      std::cout << "loc: " << index << " flow: " << flowByVertex[index] << std::endl;
    }

    std::vector<double> violations;
    int numSrcAdded = 0;
    if (stats.lpIterations > 10)
    {
      numSrcAdded += routeDD.findSRCs(cutPrimal, params.numLagCuts, violations);
    }
    if (numSrcAdded > 0)
    {
      cutAdded = true;
      stats.numCuts += numSrcAdded;
      std::cout << "adding src cuts" << std::endl;
      std::vector<double> zeroViolations;
      zeroViolations.resize(violations.size());
      addSRCCuts(dual.srcDuals, zeroViolations);
    }
  }

  resizeMultipliers(dual, bestDual);
  resizeMultipliers(dual, bestDualsArcFixingLag);
  resizeMultipliers(dual, bestDualArcFixingLP);

  auto endTimeCut = std::chrono::high_resolution_clock::now();
  auto totalTimeCut = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeCut - startTimeCut).count();
  stats.millisecondsFindingCuts += totalTimeCut;

  return cutAdded;
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

      // special src paths have 0,0 start
      if ((currLoc == vrptw.numLocations) && (currLoc == nextLoc))
      {
        continue;
      }

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

void VRPTWDDSolver::addRCCs(const std::vector<int>& edgeTail, const std::vector<int>& edgeHead, const std::vector<double>& edgeFlow, std::vector<int>& sequenceArcs, int maxNumCuts, bool& cutAdded, Dual& dual)
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
    return;
  }

  int numCuts = MyCutsCMP->Size;
  if (numCuts > 0)
  {
    std::cout << "checking cuts, max violation: " << maxViolation << ", num cuts: " << numCuts << std::endl;
    for (int cutIndex=0; cutIndex<std::min(numCuts,maxNumCuts); ++cutIndex)
    {
      std::vector<int> cutSet;
      int cutsetSize = MyCutsCMP->CPL[cutIndex]->IntListSize;
      for (int locationIndex=1; locationIndex<=cutsetSize; ++locationIndex)
      {
        int location = MyCutsCMP->CPL[cutIndex]->IntList[locationIndex];
        cutSet.push_back(location);
      }
      double RHS = MyCutsCMP->CPL[cutIndex]->RHS;

      // try to get type 3 rccs working
      bool newCut = false;
      bool cutExisted = false;
      cutExisted = routeDD.addCapCutSet(cutSet, sequenceArcs, RHS, RCCType::Type1, params.useScaling);
      newCut = true;
      /*
      if (static_cast<int>(cutSet.size()) <= vrptw.numLocations / 2)
      {
      }
      else
      {
        cutExisted = routeDD.addCapCutSet(cutSet, sequenceArcs, RHS, RCCType::Type3, params.useScaling);
      }
      */
      if (!cutExisted && newCut)
      {
        stats.numCuts = stats.numCuts + 1;
        cutAdded = true;
        dual.capDuals.push_back(0);
        std::cout << "add capacity cutset type: ";
        for (int location : cutSet)
        {
          std::cout << location << " ";
        }
        std::cout << "<= " << RHS << std::endl;
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
  double maxViolation = 0.0;
  COMBSEP_SeparateCombs(vrptw.numLocations-1,
                         &(vrptw.demandsForCombs[0]),
			 vrptw.capacity,
                         qMin,
                         static_cast<int>(edgeTail.size()-1), // num edges in solution
                         &(edgeTail[0]), // tails of edges
                         &(edgeHead[0]), // heads of edges
                         &(edgeFlow[0]), // flow value of edges
                         params.numLagCuts, // max cuts to be returned
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
      std::cout << "new comb with " << numTeeth << " teeth" << std::endl;

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
      teeth.push_back(handle);

      for (int toothIndex=1; toothIndex<=numTeeth; ++toothIndex)
      {
        std::set<int> tooth;
        std::cout << " tooth: ";
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

	teeth.push_back(tooth);
      }
      routeDD.addCombCutTeeth(teeth);

      int RHS = MyCutsCMP->CPL[cutIndex]->RHS;
      std::cout << " >= " << RHS << std::endl;
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
  for (int index=srcDuals.size(); index<routeDD.getNumSrcCuts(); ++index)
  {
    srcDuals.push_back(violations[newViolationIndex]);
    ++newViolationIndex;
  }
}

