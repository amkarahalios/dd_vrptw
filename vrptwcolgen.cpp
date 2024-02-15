#include <ilcplex/ilocplex.h>

#include "vrptwcolgen.h"
#include "cvrpsep/capsep.h"

VRPTWColGen::VRPTWColGen(VRPTW _vrptw, VRPTWDDParameters _params, PricingProblemType _pricingProblemType, InitialStateSpace initialStateSpace, int _s) : vrptw(_vrptw), routeDD(_vrptw, _params), bestLpDistance(0.0), singlePathDual(0.0), pricingProblemType(_pricingProblemType), s(_s)
{
  if (pricingProblemType == PricingProblemType::DD)
  {
    auto startCompileTime = std::chrono::high_resolution_clock::now();
    if (initialStateSpace == InitialStateSpace::Q)
    {
      routeDD.compileExactFukasawa(_s);
    }
    else if (initialStateSpace == InitialStateSpace::NG)
    {
      std::cout << "begin compiling" << std::endl;
      routeDD.compileNgRoute(_s);
      //routeDD.print();
    }
    std::cout << "DD size: " << routeDD.getNumArcsNotRemovedOrReverse() << std::endl;
    auto endCompileTime = std::chrono::high_resolution_clock::now();
    auto compileSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endCompileTime - startCompileTime).count();
    stats.millisecondsCompilingDD = stats.millisecondsCompilingDD + compileSolveTime;
  }

  initializeColumns();

  int dim = 100;
  //CMGR_CreateCMgr(&allCutsCMP,dim);
};

void VRPTWColGen::routeToColumn(const std::vector<int>& route, std::vector<int>& column)
{
  column.resize(vrptw.numLocations);
  for (int location : route)
  {
    ++column[location];
  }
}

void VRPTWColGen::addColumn(std::vector<int> route)
{
  // need to put into column format because a q-route can visit same location twice
  std::vector<int> column;
  routeToColumn(route, column);
  columns.push_back(column);

  double cost = 0.0;
  for (int index=0; index<(route.size()-1); ++index)
  {
    cost = cost + vrptw.distances[route[index]][route[index+1]];
  }
  columnCosts.push_back(cost);
}

void VRPTWColGen::initializeColumns()
{
  // all single stop routes
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    std::vector<int> route;
    route.push_back(0);
    route.push_back(location);
    route.push_back(0);
    addColumn(route);
  }

  // some full routes
  if (vrptw.problemType = ProblemType::CVRP)
  {
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
      std::cout << "adding initial route: ";
      for (int loc : route)
      {
        std::cout << loc << ",";
      }
      std::cout << std::endl;
      addColumn(route);
    }

    locationsAdded.clear();
    while (locationsAdded.size() < (vrptw.numLocations-1))
    {
      std::vector<int> route;
      route.push_back(0);
      int currentDemand = 0;
      for (int location=vrptw.numLocations-1; location>0; --location)
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
      addColumn(route);
    }
  }
};

bool VRPTWColGen::solve()
{
  // column generation - RMP <-> pricing problem
  bool solved = false;
  while (!solved)
  {
    stats.numIterations = stats.numIterations + 1;
    stats.print();

    auto startRMPTime = std::chrono::high_resolution_clock::now();
    bool rmpSolved = setupAndSolveRMP();
    auto endRMPTime = std::chrono::high_resolution_clock::now();
    auto rmpSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endRMPTime - startRMPTime).count();
    stats.millisecondsSolvingMaster = stats.millisecondsSolvingMaster + rmpSolveTime;

    if (!rmpSolved)
    {
      return false;
    }

    auto startPricingTime = std::chrono::high_resolution_clock::now();
    bool addedColumn = solvePricingProblem();
    auto endPricingTime = std::chrono::high_resolution_clock::now();
    auto pricingSolveTime = std::chrono::duration_cast<std::chrono::milliseconds>(endPricingTime - startPricingTime).count();
    stats.millisecondsSolvingPricing = stats.millisecondsSolvingPricing + pricingSolveTime;
    if (!addedColumn)
    {
      solved = true;
    }
  }

  // cut generation
  bool isCutGenerated = separateFractionalSolution();
  stats.print();

  return true;
};

bool VRPTWColGen::separateFractionalSolution()
{/*
  // need proper solution format to give to cvrpsep code
  char integerAndFeasible;
  double maxViolation;
  int numberOfEdges = 0;
  std::vector<double> demands;
  for (int demandIndex=1; demandIndex<vrptw.demands.size(); ++demandIndex)
  {
    demands.push_back(vrptw.demands[demandIndex]);
  }

  std::vector<int> edgeTails;
  std::vector<int> edgeHeads;
  std::vector<double> edgeValues;
  CnstrMgrPointer thisRoundCutsCMP;
  CAPSEP_SeparateCapCuts(vrptw.numLocations-1,
                         demands.data(),
                         vrptw.capacity,
                         numberOfEdges,
                         edgeTails.data(),
                         edgeHeads.data(),
                         edgeValues.data(),
                         allCutsCMP,
                         100,
                         0.0001,
                         0.0001,
                         &integerAndFeasible,
                         &maxViolation,
                         thisRoundCutsCMP);

  // get cuts back
  for (int cutIndex=rhsCuts.size(); cutIndex<allCutsCMP->Size; ++cutIndex)
  {
    std::vector<int> lhsList;
    for (int lhsIndex=0; lhsIndex<=allCutsCMP->CPL[cutIndex]->IntListSize; ++lhsIndex)
    {
      lhsList.push_back(allCutsCMP->CPL[cutIndex]->IntList[lhsIndex]);
    }

    int rhs = allCutsCMP->CPL[cutIndex]->RHS;
    lhsCuts.push_back(lhsList);
    rhsCuts.push_back(rhs);
  }
*/
  return false;
};

bool VRPTWColGen::setupAndSolveRMP()
{
  // setup model
  IloEnv env;
  IloModel setCoverModel(env);
  IloRangeArray coverConstraints(env);
  IloRangeArray singlePathConstraint(env);

  // setup variables
  IloNumVarArray x(env, columns.size());
  for (int columnIndex=0; columnIndex<columns.size(); ++columnIndex)
  {
    x[columnIndex] = IloNumVar(env, 0, 1);
  }

  // setup constraints
  for (int location=0; location<vrptw.numLocations; ++location)
  {
    IloExpr lhs(env);

    for (int columnIndex=0; columnIndex<columns.size(); ++columnIndex)
    {
      lhs += columns[columnIndex][location] * x[columnIndex];
    }

    if (location == 0)
    {
      coverConstraints.add(lhs >= 0);
    }
    else
    {
      coverConstraints.add(lhs >= 1);
    }
  }
  setCoverModel.add(coverConstraints);

  // one path for TSPs
  if (vrptw.fixedNumPaths == FIXED_NUM_PATHS)
  {
    IloExpr onePath(env);
    for (int columnIndex=0; columnIndex<columns.size(); ++columnIndex)
    {
      onePath += x[columnIndex];
    }
    singlePathConstraint.add(onePath == 1);
    setCoverModel.add(singlePathConstraint);
  }

  // setup objective
  IloExpr objective(env);
  for (int columnIndex=0; columnIndex<columns.size(); ++columnIndex)
  {
    objective += columnCosts[columnIndex] * x[columnIndex];
  }
  setCoverModel.add(IloMinimize(env, objective));
  objective.end();

  // solve model
  IloCplex solver(setCoverModel);
  solver.setOut(env.getNullStream());
  solver.setWarning(env.getNullStream());
  solver.setError(env.getNullStream());
  //solver.setParam(IloCplex::Param::TimeLimit, timelimit);
  solver.setParam(IloCplex::Param::Threads, 1);
  solver.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
  //solver.exportModel("LPflowmodel.lp");
  solver.solve();

  // get results
  stats.solutionValue = solver.getObjValue();
  IloAlgorithm::Status solverStatus = solver.getStatus();
  if (solverStatus == IloAlgorithm::Optimal)
  {
    // store results
    lpSolution.clear();
    lpSolution.resize(columns.size());
    for (int columnIndex=0; columnIndex<columns.size(); ++columnIndex)
    {
      lpSolution[columnIndex] = solver.getValue(x[columnIndex]);
      DBG(
        if (lpSolution[columnIndex] > 0)
        {
          std::cout << "primal [" << columnIndex << "]: " << lpSolution[columnIndex] << std::endl;
          std::cout << "col: ";
          for (int loc : columns[columnIndex])
          {
            std::cout << loc << ",";
          }
          std::cout << std::endl;
          std::cout << "cost: " << columnCosts[columnIndex] << std::endl;
        }
      )
    }

    // store duals
    IloNumArray setCoverDuals(env);
    solver.getDuals(setCoverDuals, coverConstraints);
    setCoverDualVariables.clear();
    setCoverDualVariables.resize(vrptw.numLocations);
    for (int dualIndex=0; dualIndex<vrptw.numLocations; ++dualIndex)
    {
      setCoverDualVariables[dualIndex] = static_cast<double>(setCoverDuals[dualIndex]);
      std::cout << "dual [" << dualIndex << "]: " << setCoverDualVariables[dualIndex] << std::endl;
    }

    if (vrptw.fixedNumPaths == FIXED_NUM_PATHS)
    {
      IloNumArray singlePathDualFromLP(env);
      solver.getDuals(singlePathDualFromLP, singlePathConstraint);
      singlePathDual = singlePathDualFromLP[0];
      std::cout << "single path dual: " << singlePathDual << std::endl;
    }

    env.end();
    return true;
  }
  else
  {
    std::cout << "results not optimal" << std::endl;
    std::cout << solverStatus << std::endl;
    env.end();
    return false;
  }

  env.end();
}

bool VRPTWColGen::solvePricingProblem()
{
  bool addedColumnNegativeReducedCost = false;
  if (pricingProblemType == PricingProblemType::DD)
  {
    addedColumnNegativeReducedCost = solvePricingProblemDD();
  }
  else if (pricingProblemType == PricingProblemType::DP)
  {
    addedColumnNegativeReducedCost = solvePricingProblemFukasawa();
  }

  if (!addedColumnNegativeReducedCost)
  {
    std::cout << "no negative reduced cost paths" << std::endl;
    stats.lowerBound = stats.solutionValue;
    return false;
  }

  return true;
}

struct DPEntry
{
  DPEntry(double _cost, int _priorEntryRow, int _priorEntryCol, int _priorEntryIndex) : cost(_cost), priorEntryRow(_priorEntryRow), priorEntryCol(_priorEntryCol), priorEntryIndex(_priorEntryIndex), isSet(true) {}
  DPEntry() : cost(100000), priorEntryRow(0), priorEntryCol(0), priorEntryIndex(0), isSet(false) {}

  double cost;
  int priorEntryRow;
  int priorEntryCol;
  int priorEntryIndex;
  bool isSet;
};

bool VRPTWColGen::solvePricingProblemFukasawa()
{
  // solve using dynamic programming per Fukasawa et al. 2006
  // M is a (capacity x n) matrix
  // remove 2-cycles by keeping two shortest paths at each entry
  // solve from top left to bottom right
  std::vector<std::vector<DPEntry>> M;
  M.resize(vrptw.capacity+1);
  for (int demand=0; demand<=vrptw.capacity; ++demand)
  {
    M[demand].resize(vrptw.numLocations);
  }

  // initialize with one location routes
  for (int location=1; location<vrptw.numLocations; ++location)
  {
    double cost = vrptw.distances[0][location] - setCoverDualVariables[location];
    DPEntry dpEntry(cost, 0, 0, 0);
    M[vrptw.demands[location]][location] = dpEntry;
  }

  // loop over the matrix from top left to bottom right
  for (int demand=0; demand<=vrptw.capacity; ++demand)
  {
    for (int location=1; location<vrptw.numLocations; ++location)
    {
      if (M[demand][location].isSet)
      {
        for (int toLocation=1; toLocation<vrptw.numLocations; ++toLocation)
        {
          if (location == toLocation)
          {
            continue;
          }

          int newDemand = demand + vrptw.demands[toLocation];
          if (newDemand <= vrptw.capacity)
          {
            double newCost = M[demand][location].cost + vrptw.distances[location][toLocation] - setCoverDualVariables[toLocation];
            if (!M[newDemand][toLocation].isSet || (M[newDemand][toLocation].cost > newCost))
            {
              DPEntry dpEntry(newCost, demand, location, 0);
              M[newDemand][toLocation] = dpEntry;
            }
          }
        }
      }
    }
  }

  // extend to depot and choose negative cost if there is one
  double minCost = 1;
  int minCostRow = 0;
  int minCostCol = 0;
  for (int demand=0; demand<=vrptw.capacity; ++demand)
  {
    for (int location=1; location<vrptw.numLocations; ++location)
    {
      double cost = M[demand][location].cost + vrptw.distances[location][0];
      if (cost < minCost)
      {
        minCost = cost;
        minCostRow = demand;
        minCostCol = location;
      }
    }
  }

  // build route backwards then reverse
  bool continueRoute = true;
  std::vector<int> newRoute;
  newRoute.push_back(0);
  newRoute.push_back(minCostCol);
  while (continueRoute)
  {
    int previousRow = M[minCostRow][minCostCol].priorEntryRow;
    int previousCol = M[minCostRow][minCostCol].priorEntryCol;
    if (M[previousRow][previousCol].isSet)
    {
      newRoute.push_back(previousCol);
      minCostRow = previousRow;
      minCostCol = previousCol;
    }
    else
    {
      continueRoute = false;
    }
  }
  newRoute.push_back(0);
  std::reverse(newRoute.begin(), newRoute.end());

  if (minCost < -0.000001)
  {
    addColumn(newRoute);
    return true;
  }

  return false;
}

bool VRPTWColGen::solvePricingProblemDD()
{
  routeDD.setCoeffsAsDistancesMinusLagrangean(setCoverDualVariables);

  std::vector<int> newRoute;
  double shortestPathLength = routeDD.computeShortestPathBFS(ShortestPathMode::SHORTEST_PATH, newRoute);
  shortestPathLength = shortestPathLength - singlePathDual;
  DBG(
    std::cout << "spl: " << shortestPathLength << std::endl;
    std::cout << "adding route: ";
    for (int loc : newRoute)
    {
      std::cout << loc << ",";
    }
    std::cout << std::endl;
  )
  if (shortestPathLength < -0.00000001)
  {
    addColumn(newRoute);
    // lower bound needs revamp without num trucks bound
    //stats.lowerBound = std::max(stats.lowerBound, stats.solutionValue + (vrptw.
    stats.print();
    return true;
  }

  return false;
}
