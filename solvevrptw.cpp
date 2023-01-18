#include <iostream>
#include <chrono>

#include "vrptw.h"
#include "vrptwdecisiondiagram.h"
#include "vrptwcolgen.h"
#include "vrptwddsolver.h"
#include "cvrpsep/cnstrmgr.h"

void usage()
{
  std::cout << std::endl;
  std::cout << "USAGE: dd_vrp -fileName <fileName> -solver <solver> [options]" << std::endl;
  std::cout << "USAGE: dd_vrp -fileName <fileName> -solver COL_GEN DD/DP Q/NG s/k maxS cuts timeout" << std::endl;
  std::cout << "USAGE: dd_vrp -fileName <fileName> -solver COL_ELIM LP/LAG Q/NG s/k maxS cuts timeout" << std::endl;
}

int main(int argc, char** argv)
{
  if (argc < 3)
  {
    usage();
    return -1;
  }

  std::string fileName = argv[1];
  VRPTW vrptw(fileName);

  auto start = std::chrono::high_resolution_clock::now();
  std::string solverName = argv[2];
  if (solverName == "COL_GEN")
  {
    std::string pricingProblemTypeString = argv[3];
    PricingProblemType pricingProblemType = PricingProblemType::DD;
    if (pricingProblemTypeString == "DP")
    {
      pricingProblemType = PricingProblemType::DP;
    }
    else if (pricingProblemTypeString == "DD")
    {
      pricingProblemType = PricingProblemType::DD;
    }
    else
    {
      return -1;
    }

    std::string initialStateSpaceString = argv[4];
    InitialStateSpace initialStateSpace = InitialStateSpace::NG;
    if (initialStateSpaceString == "NG")
    {
      initialStateSpace = InitialStateSpace::NG;
    }
    else if (initialStateSpaceString == "Q")
    {
      initialStateSpace = InitialStateSpace::Q;
    }
    else
    {
      return -1;
    }

    std::string sValue = argv[5];
    int s = std::stoi(sValue);

    VRPTWColGen colGenModel(vrptw, pricingProblemType, initialStateSpace, s);
    bool solved = colGenModel.solve();
  }
  else if (solverName == "COL_ELIM")
  {
    std::string solveTypeString = argv[3];
    LPSolveType solveType = LPSolveType::LPSolver;
    if (solveTypeString == "LAG")
    {
      solveType = LPSolveType::LAGSolver;
    }
    else if (solveTypeString == "LP")
    {
      solveType = LPSolveType::LPSolver;
    }
    else
    {
      return -1;
    }

    std::string initialStateSpaceString = argv[4];
    InitialStateSpace initialStateSpace = InitialStateSpace::NG;
    if (initialStateSpaceString == "NG")
    {
      initialStateSpace = InitialStateSpace::NG;
    }
    else if (initialStateSpaceString == "Q")
    {
      initialStateSpace = InitialStateSpace::Q;
    }
    else
    {
      return -1;
    }

    std::string sValue = argv[5];
    std::string maxSValue = argv[6];
    std::string useCutsString = argv[7];
    std::string timeoutString = argv[8];
    int s = std::stoi(sValue);
    int maxS = std::stoi(maxSValue);

    bool useCuts = false;
    if (useCutsString == "Y")
    {
      useCuts = true;
    }

    int timeout = std::stoi(timeoutString);

    VRPTWDDSolver ddSolver(vrptw, solveType, initialStateSpace, s, maxS, useCuts, timeout);
    bool solved = ddSolver.solve();
  }
  else
  {
    usage();
    return -1;
  }

  auto stop = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
  std::cout << "time elapsed: " << duration.count() << std::endl;

  return 0;
}
