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
  if (params.stateSpace == StateSpace::Q)
  {
    routeDD.compileExactFukasawa(params.ngSetSize);
  } 
  else if (params.stateSpace == StateSpace::NG)
  {
    routeDD.compileNgRoute(params.ngSetSize);
    //routeDD.checkLRC121SolutionPossible();
  }

  auto endCompileTime = std::chrono::high_resolution_clock::now();
  auto compileSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endCompileTime - startCompileTime).count();
  stats.millisecondsCompiling = stats.millisecondsCompiling + compileSolveTime;
  std::cout << "done compiling DD" << std::endl;

  // set and log parameters
  // static parameters
  infeasibleRoutesBatchSize = 1;
  deactivateCutValueThreshold = 0.01;
  deactivateCutIterThreshold = 1000000000;
  percentFixedToChangeToCPLEX = 97.5;
  numArcsToChangeToCPLEX = 100000;
  numLagCuts = 20;
  kappaIterations = 40;
  muPercentImproved = 0.001;
  percentGapToStartCuts = 0.95;

  std::cout << "batch size for lag: " << infeasibleRoutesBatchSize << std::endl;
  std::cout << "deactivate cut value threshold: " << deactivateCutValueThreshold << std::endl;
  std::cout << "deactivate cut iter threshold: " << deactivateCutIterThreshold << std::endl;
  std::cout << "percent arcs fixed to change to CPLEX: " << percentFixedToChangeToCPLEX << std::endl;
  std::cout << "number arcs fixed to change to CPLEX: " << numArcsToChangeToCPLEX << std::endl;
  std::cout << "kappa iterations: " << kappaIterations << std::endl;
  std::cout << "num lag cuts: " << numLagCuts << std::endl;
  std::cout << "mu percent improved: " << muPercentImproved << std::endl;
  std::cout << "percent to start cuts: " << percentGapToStartCuts << std::endl;

  bestDualsArcFixingPercent = 0.0;
  bestDualValue = 0.0;
  stepSizeMultiplier = 1.0;
  stepSizeMultiplierIteration = 0;
  stepSizeMultiplierIterationCutoff = 5;
  alphaLowerBound = 0.1;
  alphaLowerBoundIteration = 0;
  targetLowerBound = 1.0;
  targetBoundIncrease = 1.05;
  minStepSizeMultiplier = 0.0001;

  shouldRestart = false;
  numCapCutUpdateGroups = 5;
  numIterationsRestart = 5;
  numIterationsDeactivations = 2;

  // va parameters
  std::cout << "step size iteration cutoff: " << stepSizeMultiplierIterationCutoff << std::endl;
  std::cout << "target bound increase: " << targetBoundIncrease << std::endl;
  std::cout << "min step size multiplier: " << minStepSizeMultiplier << std::endl;
  std::cout << "num cap cut update groups: " << numCapCutUpdateGroups << std::endl;
  std::cout << "num iterations restart: " << numIterationsRestart << std::endl;
  std::cout << "num iterations deactivation: " << numIterationsDeactivations << std::endl;

  // dynamic parameters
  std::cout << "ngSetSize: " << params.ngSetSize << std::endl;
  std::cout << "change to LP: " << params.changeToLP << std::endl;
  std::cout << "var fixing: " << params.useVariableFixing << std::endl;
  std::cout << "useSeparations: " << params.useSeparations << std::endl;
  std::cout << "muSSP: " << params.useMuSSP << std::endl;
  std::cout << "repairDuals: " << params.repairDuals << std::endl;
  std::cout << "robustCuts: " << params.useRobustCuts << std::endl;
  std::cout << "nonRobustCuts: " << params.useNonRobustCuts << std::endl;
  std::cout << "volume algo: " << params.useVolumeAlgorithm << std::endl;
  std::cout << "use scaling: " << params.useScaling << std::endl;
  std::cout << "use sparse RCCs: " << params.useSparseRCCs << std::endl;
  std::cout << "use restarts: " << params.useRestarts << std::endl;

  CMGR_CreateCMgr(&MyCutsCMP,Dim);
  CMGR_CreateCMgr(&MyOldCutsCMP,Dim);

  stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
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

      // do not add same cut if encountered twice
      bool alreadyExists = false;
      bool tooSimilarToExisting = false;
      auto cutSetAsSet = std::set<int>(cutSet.begin(), cutSet.end());
      for (int index=0; index<routeDD.getNumCapCuts(); ++index)
      {
        if (routeDD.isCapCutActive(index))
        {
          auto existingCutSet = routeDD.getCapCutSet(index);
          auto existingCutSetAsSet = std::set<int>(existingCutSet.begin(), existingCutSet.end());
          if (cutSetAsSet == existingCutSetAsSet)
          {
            alreadyExists = true;
            break;
          }
/*
          if (std::includes(cutSetAsSet.begin(), cutSetAsSet.end(), existingCutSetAsSet.begin(), existingCutSetAsSet.end()))
          {
            tooSimilarToExisting = true;
            break;
          }
*/
        }
      }

      if (alreadyExists)
      {
        std::cout << "rcc already exists, do not add" << std::endl;
      }
      else if (tooSimilarToExisting)
      {
        std::cout << "rcc too similar to existing, do not add" << std::endl;
      }
      else if (cutSet.size() > vrptw.numLocations - 5)
      {
        std::cout << "not adding, size of cut set too large" << std::endl;
      }
      else
      {
        std::set<int> cutSetAsSet(cutSet.begin(), cutSet.end());
        cutSets.push_back(cutSetAsSet);
        if (params.lpSolveType == LPSolveType::LPSolver)
        {
          dual.capDuals.push_back(0);
          dual.capCutUsage.push_back(false);
        }
        else
        {
          dual.capDuals.push_back(1);
          dual.capCutUsage.push_back(true);
        }
        stats.numCuts = stats.numCuts + 1;

        if ((cutSet.size() <= vrptw.numLocations / 2) || !params.useSparseRCCs)
        {
          routeDD.addCapCutSet(cutSet, rccArcs, RHS, RCCType::Type1, params.lpSolveType, params.useScaling);
          std::cout << "using type 1 with rhs: " << RHS << std::endl;
        }
        else
        {
          double newRHS = RHS - cutSet.size() + ((vrptw.numLocations - 1) - cutSet.size());
          routeDD.addCapCutSet(cutSet, rccArcs, newRHS, RCCType::Type3, params.lpSolveType, params.useScaling);
          std::cout << "using type 3 with rhs: " << newRHS << std::endl;
        }
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
                         numLagCuts, // max cuts to be returned
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
  for (int index=srcDuals.size(); index<routeDD.getNumCliqueCuts(); ++index)
  {
    srcDuals.push_back(violations[newViolationIndex]);
    ++newViolationIndex;
  }
}

void VRPTWDDSolver::resizeMultipliers(const Dual& dual1, std::vector<Dual>& duals2)
{
  for (int index=0; index<duals2.size(); ++index)
  {
    Dual& dual2 = duals2[index];
    dual2.lambda.resize(dual1.lambda.size());
    dual2.capDuals.resize(dual1.capDuals.size());
    dual2.capCutUsage.resize(dual1.capCutUsage.size());
    dual2.combDuals.resize(dual1.combDuals.size());
    dual2.srcDuals.resize(dual1.srcDuals.size());
  }
}

void VRPTWDDSolver::resizeMultipliers(const Dual& dual1, Dual& dual2)
{
  dual2.lambda.resize(dual1.lambda.size());
  dual2.capDuals.resize(dual1.capDuals.size());
  dual2.capCutUsage.resize(dual1.capCutUsage.size());
  dual2.combDuals.resize(dual1.combDuals.size());
  dual2.srcDuals.resize(dual1.srcDuals.size());
}

void VRPTWDDSolver::resizeMultipliersAndCopy(const Dual& dual1, Dual& dual2)
{
  dual2.lambda.resize(dual1.lambda.size());

  if (dual1.capDuals.size() > dual2.capDuals.size())
  {
    int startIndex = dual2.capDuals.size();
    dual2.capDuals.resize(dual1.capDuals.size());
    dual2.capCutUsage.resize(dual1.capCutUsage.size());
 
    for (int index = startIndex; index<dual1.capDuals.size(); ++index)
    {
      dual2.capDuals[index] = dual1.capDuals[index];
    }
  }

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

/*
  for (int index=0; index<dual.capDuals.size(); ++index)
  {
    dual.capDuals[index] = 0;
  }
  for (int index=0; index<dual.combDuals.size(); ++index)
  {
    dual.combDuals[index] = 0;
  }
  for (int index=0; index<dual.srcDuals.size(); ++index)
  {
    dual.srcDuals[index] = 0;
  }

  dual.fixedPathDual = 0;
*/
}

bool VRPTWDDSolver::solve(bool shouldSolveIP)
{
  bool lpFlowType = FlowType::LP;

  Dual dual;
  initializeDual(dual);
  bestDual = dual;

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
      if (stats.lpIterations > 1)
      {
        double percentArcsFixed = 0.0;
        if (!changedLagToLP)
        {
          if (params.useVariableFixing)
          {
            //std::cout << "initial value: " << routeDD.getDualObjectiveValue(dual, LPSolveType::LPSolver) << std::endl;
            //repairMultipliers(dual, LPSolveType::LPSolver);
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
        //repairMultipliers(bestDualArcFixing, LPSolveType::LAGSolver);
        routeDD.fixArcs(bestDualsArcFixing, LPSolveType::LAGSolver);
      }

      if (lpFlowType == FlowType::LP)
      {
        //routeDD.strengthenSRCs(averageRouteLength);
        //routeDD.checkAn32k5SolutionPossible();
        //routeDD.checkAn36k5SolutionPossible();
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

      // if switching to LP, go right to next iter and use cuts
      if (params.lpSolveType == LPSolveType::LPSolver)
      {
        changedLagToLP = true;
        params.useRobustCuts = true;
        params.useNonRobustCuts = true;
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
      if (lpFlowType == FlowType::LP)
      {
        auto startTimeCut = std::chrono::high_resolution_clock::now();

        if (params.useRobustCuts)
        {
          // RCC - rounded capacity cuts
          std::vector<int> edgeTail;
          std::vector<int> edgeHead;
          std::vector<double> edgeFlow;
          std::vector<int> rccArcs;
          std::vector<double> rccArcFlows;
          routeDD.convertSolutionForVRPTWSep(edgeTail, edgeHead, edgeFlow, rccArcs, rccArcFlows);
          addRCCs(edgeTail, edgeHead, edgeFlow, rccArcs, 100, cutAdded, dual);

          // Strengthened Combs
          addCombs(edgeTail, edgeHead, edgeFlow, cutAdded);
          dual.combDuals.resize(routeDD.getNumCombCuts());
        }

        if (params.useNonRobustCuts)
        {
          // Subset Row Cuts
          std::vector<int> infeasibleRoute;
          routeDD.decomposeRoutes(infeasibleRoute, routeFlows, decomposedRoutes, decomposedArcs, DecompositionReason::DECOMPOSE);
          averageRouteLength = averageRouteLength * numRoutesInAverage;
          for (auto route : decomposedRoutes)
          {
            averageRouteLength = averageRouteLength + route.size();
          }
          numRoutesInAverage += decomposedRoutes.size();
          averageRouteLength = averageRouteLength / numRoutesInAverage;

          // currently adding by full separation, can strengthen with all up / all down
          std::vector<double> violations;
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
          numSrcAdded = routeDD.findSRC3s(primal, 100, violations);
          numSrcAdded += routeDD.findSRC4s(primal, 100, violations);
          //numSrcAdded += routeDD.findSRC5V1s(primal, 100, violations);
          //numSrcAdded += routeDD.findSRC5V2s(primal, 100, violations);

          if (numSrcAdded > 0)
          {
            cutAdded = true;
            stats.numCuts = stats.numCuts + numSrcAdded;
          }

          std::vector<double> zeroViolations;
          zeroViolations.resize(violations.size());
          addSRCCuts(dual.srcDuals, zeroViolations);

          currDualsArcFixing.clear();
        }

        resizeMultipliers(dual, bestDual);
        resizeMultipliers(dual, bestDualsArcFixing);
        resizeMultipliers(dual, currDualsArcFixing);
 
        auto endTimeCut = std::chrono::high_resolution_clock::now();
        auto totalTimeCut = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeCut - startTimeCut).count();
        stats.millisecondsFindingCuts += totalTimeCut;
      }

      bool stopFindingInfeasibilities = false;
      std::vector<std::vector<int>> infeasibilities;
      if ((numSrcAdded == 0) && params.useSeparations)
      {
        while (!stopFindingInfeasibilities)
        {
          std::vector<int> infeasibleRoute;
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
      /*
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
            //infeasibilities.push_back(routeArcs);
          }
        }
      }
      */

      if (infeasibilities.size() > 0)
      {
        for (auto infeasibleRoute : infeasibilities)
        {
          if (routeDD.doesRouteExistByArcs(infeasibleRoute))
          {
            stats.numSeparations = stats.numSeparations + 1;
            routeDD.separateInfeasibleRoute(infeasibleRoute);
            currDualsArcFixing.clear();
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
  stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, duals, false);
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
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::IP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, duals, false);
  }
  else
  {
    double oldLb = stats.lowerBound;
    stats.lowerBound = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, duals, false);
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
    flowObj = routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::USE_CG, duals, false);
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
  auto startTimeRepair = std::chrono::high_resolution_clock::now();
  routeDD.clearRelaxedSrcs();

  while (true)
  {
    routeDD.setCoeffsAsDistancesMinusLagrangeanPlusCapDualsPlusSrcDualsPlusCombDuals(repairedDual, solveType);

    std::vector<int> treeByParentArcs;
    treeByParentArcs.resize(routeDD.getNodes().size());
    std::vector<int> shortestPathByArc;
    double longestShortestPathLength = 0;
    double shortestPathLength = routeDD.computeShortestPathBFSWang(treeByParentArcs, shortestPathByArc, longestShortestPathLength);
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
 
  auto endTimeRepair = std::chrono::high_resolution_clock::now();
  auto totalTimeRepair = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeRepair - startTimeRepair).count();
  stats.millisecondsRepairingLAG += totalTimeRepair;
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
 
  std::cout << "combDuals: ";
  for (int index=0; index<dual.combDuals.size(); ++index)
  {
    std::cout << dual.combDuals[index] << ",";
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
    std::cout << "gamma: " << i << " val: " << gamma[i] << std::endl;
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
      gamma[gammaIndex] = (-1 * routeDD.getCapCutSetRHS(i) * routeDD.getCapCutSetScaling(i)) + cutValues[i];
      if ((i % numCapCutUpdateGroups != stats.numLagIterations % numCapCutUpdateGroups) || ((dual.capDuals[i] <= 0.001) && (gamma[gammaIndex] <= 0.001)))
      {
        gamma[gammaIndex] = 0;
      }
      else
      {
        std::cout << "only consider one gamma: " << i << " val: " << gamma[gammaIndex] << std::endl;
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
  double psiStar = std::min(stats.upperBound, lagrangeanLowerBound * (1 + eta));

  // alpha_(k) = (psi_(star) - psi(lambda(k))) / ||gamma_(k)||_(2)^2
  double alpha = (psiStar - lagrangeanLowerBound) / normGammaSquared;
 
  std::cout << "alpha: " << alpha << std::endl;
  std::cout << "psiStar: " << psiStar << std::endl;
  std::cout << "lagLB: " << lagrangeanLowerBound << std::endl;
  std::cout << "||gamma||^2: " << normGammaSquared << std::endl;

  printMultipliers(dual);
  printMultipliers(bestDual);

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
      if (dual.capDuals[muIndex] < deactivateCutValueThreshold)
      {
        capCutTooSmallCounters[muIndex] = capCutTooSmallCounters[muIndex] + 1;
        if (capCutTooSmallCounters[muIndex] > deactivateCutIterThreshold)
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
      if (dual.srcDuals[srcIndex] < deactivateCutValueThreshold)
      {
        cliqueCutTooSmallCounters[srcIndex] = cliqueCutTooSmallCounters[srcIndex] + 1;
        if (cliqueCutTooSmallCounters[srcIndex] > deactivateCutIterThreshold)
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

  // restart
  if (shouldRestart)
  {
    std::cout << "restarting" << std::endl;
    initializeDual(dual);
    initializeDual(bestDual);
    bestDualValue = 1.0;
    targetLowerBound = 1.0;
    targetBoundIncrease = 1.01;
    stepSizeMultiplier = 1.0;
    stepSizeMultiplierIteration = 0;
    shouldRestart = false;

    //kappaIterations = 100;
    //numLagCuts = 20;
    muPercentImproved = 0.0001;
    minStepSizeMultiplier = 0.0001;
    numCapCutUpdateGroups = 1;

    // switch to use SGD
    if ((sgdAlgo == SGDAlgorithm::VA) && (stats.lpIterations > numIterationsRestart))
    {
      std::cout << "switch to SGD" << std::endl;
      sgdAlgo = SGDAlgorithm::SGD;
      params.useRobustCuts = false;
      params.useNonRobustCuts = false;
    }
    else if (sgdAlgo == SGDAlgorithm::SGD)
    {
      std::cout << "switch to VA" << std::endl;
      sgdAlgo = SGDAlgorithm::VA;
      params.useRobustCuts = true;
      params.useNonRobustCuts = true;
    }
  }

  // check best LP
  /*
  Dual testDual;
  testDual = dual;
  routeDD.setCoeffsAsDistances();
  routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, testDual, true);
  routeDD.setupAndSolveFlowModel(FlowType::LP, IncludeCoverConstraints::Y, UseColumnGeneration::NO_CG, testDual, false);
  */

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
      if (params.useVariableFixing && (bestDualsArcFixingPercent > 0.001))
      {
        //repairMultipliers(bestDualsArcFixing, LPSolveType::LAGSolver);
        auto startTimeFix = std::chrono::high_resolution_clock::now();
        routeDD.fixArcs(bestDualsArcFixing, LPSolveType::LAGSolver);
        auto endTimeFix = std::chrono::high_resolution_clock::now();
        auto totalTimeFix = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeFix - startTimeFix).count();
        stats.millisecondsFix += totalTimeFix;
      }

      // Run muSSP
      std::vector<std::vector<int>> shortestPaths;
      auto startSSPTime = std::chrono::high_resolution_clock::now();
      double lagrangeanLowerBound = routeDD.solveMinCostFlowModelWang(dual, shortestPaths, isDualFeasible, minReducedCost, longestShortestPathLength);
      stats.numSSPIterations = stats.numSSPIterations + shortestPaths.size();
      auto endMuSSPTime = std::chrono::high_resolution_clock::now();
      auto sspSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endMuSSPTime - startSSPTime).count();
      stats.millisecondsSolvingSSP = stats.millisecondsSolvingSSP + sspSolveTime;
 
      // Heuristic for trying arc fixing or repair
      if (stats.lowerBound + longestShortestPathLength > stats.upperBound)
      {
        tryArcFixingOrRepair = true;
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
        lagrangeanLowerBound += dual.capDuals[index] * routeDD.getCapCutSetRHS(index) * -1 * routeDD.getCapCutSetScaling(index);
        dualBoundWithoutFixedPathDual += dual.capDuals[index] * routeDD.getCapCutSetRHS(index) * -1 * routeDD.getCapCutSetScaling(index);
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
      }

      std::cout << "lag: " << lagrangeanLowerBound << std::endl;

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
        std::cout << "new best dual, value: " << bestDualValue << std::endl;
        for (int index=0; index<bestDual.capCutUsage.size(); ++index)
        {
          if (bestDual.capDuals[index] > 0.000001)
          {
            bestDual.capCutUsage[index] = true;
          }
        }
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
      auto endTimeDecompose = std::chrono::high_resolution_clock::now();
      auto totalTimeDecompose = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeDecompose - startTimeDecompose).count();
      stats.millisecondsDecompose += totalTimeDecompose;

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
      /*
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
          alphaTry = std::min(alphaTry + (alphaLowerBound / 2), alphaLowerBound);
          //alphaTry = std::min(alphaTry  alphaLowerBound);
        }
      }
      */
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

      // Update dual
      auto startTimeUpdateDual = std::chrono::high_resolution_clock::now();
      if (sgdAlgo == SGDAlgorithm::VA)
      {
        updateMultipliersVolumeAlgorithm(dual, primal, stats.numLagIterations);
      }
      else
      {
        updateMultipliers(dual, lagrangeanLowerBound, stats.numLagIterations);
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
      if ((numLagIterations - lastMuImprovedIteration) > kappaIterations)
      {
        shouldTerminate = true;
      }

      if ((stats.getNumSeconds() >= params.timeoutSeconds) || (stats.lowerBound + 0.01 > stats.upperBound))
      {
        shouldTerminate = true;
      }

      // Repair Dual
      //if (isDualFeasible & (vrptw.fixedNumPaths != FixedNumPaths::FIXED_NUM_PATHS))
      if (isDualFeasible)
      {
        if (params.useVariableFixing && tryArcFixingOrRepair)
        {
          //repairMultipliers(repairedDual, LPSolveType::LAGSolver);
          double previousPercentFixed = percentFixed;
          percentFixed = routeDD.fixArcs(repairedDual, LPSolveType::LAGSolver);
          std::cout << "percent fixed: " << percentFixed << std::endl;
          if (((percentFixed - previousPercentFixed) > 1) && (currDualsArcFixing.size() < 5))
          {
            currDualsArcFixing.push_back(repairedDual);
          }
        }

        if (percentFixed > bestDualsArcFixingPercent)
        {
          std::cout << "updating best lambda arc fixing" << std::endl;
          bestDualsArcFixing = currDualsArcFixing;
          bestDualsArcFixingPercent = percentFixed;
        }
      }
      else
      {
        if (params.repairDuals)
        {
          if (params.useVariableFixing && tryArcFixingOrRepair)
          {
            repairMultipliers(repairedDual, LPSolveType::LAGSolver);
            double repairedBound = routeDD.getDualObjectiveValue(repairedDual, LPSolveType::LAGSolver);

            DBG(std::cout << "repaired lb: " << repairedBound << std::endl;)
            if (bestDualValue < repairedBound)
            {
              stats.lowerBound = std::max(stats.lowerBound, repairedBound);
              stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
              printMultipliers(repairedDual);

              // NOTE(akarahal) don't allow best dual from repair, because might not have same with wang calculation
              bestDual = repairedDual;
              bestDualValue = repairedBound;
              std::cout << "new best dual from repair, value: " << bestDualValue << std::endl;
              for (int index=0; index<bestDual.capCutUsage.size(); ++index)
              {
                if (bestDual.capDuals[index] > 0.0001)
                {
                  bestDual.capCutUsage[index] = true;
                }
              }
              /*
              double percentGap = (stats.upperBound - bestDualValue) * 100.0 / stats.upperBound;
              if (percentGap < 10)
              {
                double limitToMerge = (stats.upperBound - bestDualValue) / 2.0;
                routeDD.findMergeNodesReducedCost(bestDual, LPSolveType::LAGSolver, limitToMerge);
                stats.print(routeDD.getNumArcsNotRemovedOrReverse(), routeDD.getNumFixedArcs());
              }
              */

              if (repairedBound > 0.95 * targetLowerBound)
              {
                targetLowerBound = std::min(repairedBound * 1.05, stats.upperBound);
              }
            }

            double previousPercentFixed = percentFixed;
            percentFixed = routeDD.fixArcs(repairedDual, LPSolveType::LAGSolver);
            if (((percentFixed - previousPercentFixed) > 1) && (currDualsArcFixing.size() < 5))
            {
              currDualsArcFixing.push_back(repairedDual);
            }
            std::cout << "percent fixed: " << percentFixed << std::endl;
            if (percentFixed > bestDualsArcFixingPercent)
            {
              std::cout << "updating from repair - best lambda arc fixing" << std::endl;
              bestDualsArcFixingPercent = percentFixed;
              bestDualsArcFixing = currDualsArcFixing;
            }
          }
        }
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
    if ((phaseType == PhaseType::SEPARATION) && !shouldTerminate)
    {
      for (auto infeasibleRouteToSeparateByArc : infeasibleRoutesByArc)
      {
        // when separating changes, can use locations to re-find the route
        if (routeDD.doesRouteExistByArcs(infeasibleRouteToSeparateByArc))
        {
          stats.numSeparations = stats.numSeparations + 1;
          routeDD.separateInfeasibleRoute(infeasibleRouteToSeparateByArc);
          currDualsArcFixing.clear();
        }
      }
    }

    // Cuts
    bool shouldCut = false;
    if (stats.lowerBound > stats.upperBound * percentGapToStartCuts)
    {
      shouldCut = true;
    }

    // Consider restarting
    if ((shouldTerminate && (params.useRobustCuts || params.useNonRobustCuts)) && shouldCut)
    {
      addCutsUsingCurrentPrimal(dual);
      if (params.useRestarts)
      {
        shouldRestart = true;
      }
    }

    // Switch back to VA to add more cuts?
    if (shouldTerminate && (sgdAlgo == SGDAlgorithm::SGD) && (stats.lpIterations % numIterationsRestart == 0))
    {
      if (params.useRestarts)
      {
        shouldRestart = true;
      }
    }

    infeasibleRoutesByArc.clear();
  }

  return true;
}
 
      // deactivate unused cut duals
      /*
      if (stats.lpIterations % numIterationsDeactivations == 0)
      {
        for (int capIndex=0; capIndex<bestDual.capCutUsage.size(); ++capIndex)
        {
          if (routeDD.isCapCutActive(capIndex))
          {
            if (!bestDual.capCutUsage[capIndex])
            {
              std::cout << "deactivated unused cap cut at index: " << capIndex << std::endl;
              routeDD.deactivateCapCut(capIndex);
            }
          }
        }
      }
      */

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
      gradient[gradientIndex] = (-1 * routeDD.getCapCutSetRHS(i) * routeDD.getCapCutSetScaling(i)) + cutValues[i];
    }
    else
    {
      gradient[gradientIndex] = 0;
    }
  }

  // src cuts are in the form Cx <= r so change to -Cx >= -r
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

void VRPTWDDSolver::updateMultipliersVolumeAlgorithm(Dual& dual, Primal& currPrimal, int iteration)
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
 
    std::cout << "gamma at lambda " << i << ":" << gamma[i] << std::endl;
  }

  // cap cuts are in the form Cx <= r so change to -Cx >= -r
  //int keepIndex = 0;
  //if (stats.numCuts > 0)
  //{
  //  keepIndex = stats.numLagIterations % dual.capDuals.size();
  //}
  for (int i=0; i<dual.capDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations;
    if (!routeDD.isCapCutActive(i))
    {
      bestDual.capDuals[i] = 0;
      gamma[gammaIndex] = 0;
      continue;
    }

    //if (i != keepIndex)
    gamma[gammaIndex] = (-1 * routeDD.getCapCutSetRHS(i) * routeDD.getCapCutSetScaling(i)) + cutValues[i];
    if ((i % numCapCutUpdateGroups != stats.numLagIterations % numCapCutUpdateGroups) || ((dual.capDuals[i] <= 0.001) && (gamma[gammaIndex] <= 0.001)))
    {
      gamma[gammaIndex] = 0;
    }
    else
    {
      std::cout << "only consider one gamma: " << i << " val: " << gamma[gammaIndex] << std::endl;
      normGammaSquared += std::pow(gamma[gammaIndex], 2);
    }
  }

  // src cuts are in the form Cx <= r so change to -Cx >= -r
  for (int i=0; i<dual.srcDuals.size(); ++i)
  {
    int gammaIndex = i + vrptw.numLocations + dual.capDuals.size();
    if (routeDD.isCliqueCutActive(i))
    {
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
    else
    {
      bestDual.srcDuals[i] = 0;
      gamma[gammaIndex] = 0;
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

  double alpha = stepSizeMultiplier * (targetLowerBound - bestDualValue) / normGammaSquared;
  std::cout << "alpha: " << alpha << std::endl;
  std::cout << "bestDualValue: " << bestDualValue << std::endl;
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
      if (dual.capDuals[muIndex] < deactivateCutValueThreshold)
      {
        capCutTooSmallCounters[muIndex] = capCutTooSmallCounters[muIndex] + 1;
        if (capCutTooSmallCounters[muIndex] > deactivateCutIterThreshold)
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
      if (dual.srcDuals[srcIndex] < deactivateCutValueThreshold)
      {
        cliqueCutTooSmallCounters[srcIndex] = cliqueCutTooSmallCounters[srcIndex] + 1;
        if (cliqueCutTooSmallCounters[srcIndex] > deactivateCutIterThreshold)
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
  auto startTimeCut = std::chrono::high_resolution_clock::now();

  // Rounded Capacity Cuts
  bool cutAdded = false;
  std::vector<int> edgeTail;
  std::vector<int> edgeHead;
  std::vector<double> edgeFlow;
  std::vector<int> rccArcs;
  std::vector<double> rccArcFlows;
  if (params.useRobustCuts)
  {
    convertArcIndicesForVRPTWSep(primal, edgeTail, edgeHead, edgeFlow, rccArcs, rccArcFlows);
    addRCCs(edgeTail, edgeHead, edgeFlow, rccArcs, numLagCuts, cutAdded, dual);

    // Strengthened Combs
    //addCombs(edgeTail, edgeHead, edgeFlow, cutAdded);
    //dual.combDuals.resize(routeDD.getNumCombCuts());
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

  if (params.useNonRobustCuts)
  {
    std::vector<double> violations;
    //std::cout << "finding src3s" << std::endl;
    //int numSrc3Added = routeDD.findSRC3s(primal, numLagCuts, violations);
    //std::cout << "finding src4s" << std::endl;
    //int numSrc4Added = routeDD.findSRC4s(primal, numLagCuts, violations);
    int numSrc3Added = 0;
    int numSrc4Added = 0;

    int numSrcAdded = numSrc3Added + numSrc4Added;
    if (numSrcAdded > 0)
    {
      cutAdded = true;
      stats.numCuts += numSrcAdded;
      std::cout << "adding src cuts" << std::endl;
      std::vector<double> zeroViolations;
      zeroViolations.resize(violations.size());
      addSRCCuts(dual.srcDuals, violations);

      // strengthen src using average route length
      int averageRouteLength = 0;
      for (int index=0; index<primal.xDecompositions.size(); ++index)
      {
        averageRouteLength += primal.xDecompositions[index].size();
      }
      averageRouteLength = averageRouteLength / primal.xDecompositions.size();
      //routeDD.strengthenSRCs(averageRouteLength);
 
      currDualsArcFixing.clear();
    }
  }

  //resizeMultipliersAndCopy(dual, bestDual);
  resizeMultipliers(dual, bestDual);
  resizeMultipliers(dual, bestDualsArcFixing);
  resizeMultipliers(dual, currDualsArcFixing);
  capCutTooSmallCounters.resize(dual.capDuals.size());
  cliqueCutTooSmallCounters.resize(dual.srcDuals.size());

  auto endTimeCut = std::chrono::high_resolution_clock::now();
  auto totalTimeCut = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeCut - startTimeCut).count();
  stats.millisecondsFindingCuts += totalTimeCut;

  return cutAdded;
};
