#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
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
  std::cout << "USAGE: ./solver filePath COL_ELIM paramFilePath" << std::endl;
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

    std::string stateSpaceString = argv[4];
    StateSpace stateSpace = StateSpace::NG;
    if (stateSpaceString == "NG")
    {
      stateSpace = StateSpace::NG;
    }
    else if (stateSpaceString == "Q")
    {
      stateSpace = StateSpace::Q;
    }
    else
    {
      return -1;
    }

    std::string sValue = argv[5];
    int s = std::stoi(sValue);

    VRPTWDDParameters params;
    VRPTWColGen colGenModel(vrptw, params, pricingProblemType, stateSpace, s);
    bool solved = colGenModel.solve();
  }
  else if (solverName == "COL_ELIM")
  {
    VRPTWDDParameters params;
    std::string paramFilePath = argv[3];
    std::ifstream paramIfs(paramFilePath);
    std::string line;
    int lineIndex = 0;
    while (std::getline(paramIfs, line))
    {
      std::istringstream iss(line);
      std::string paramValueString;
      int wordIndex = 0;
      while (std::getline(iss, paramValueString, ':'))
      {
        if (wordIndex == 1)
        {
          if (lineIndex == 0)
          {
            params.timeoutSeconds = std::stoi(paramValueString);
          }
          else if (lineIndex == 2)
          {
            if (paramValueString == "LAG")
            {
              params.lpSolveType = LPSolveType::LAGSolver;
            }
            else if (paramValueString == "LP")
            {
              params.lpSolveType = LPSolveType::LPSolver;
            }
            else if (paramValueString == "LAG-LP")
            {
              params.lpSolveType = LPSolveType::LPSolver;
            }
            else
            {
              return -1;
            }
          }
          else if (lineIndex == 3)
          {
            if (paramValueString == "NG")
            {
              params.stateSpace = StateSpace::NG;
            }
            else if (paramValueString == "Q")
            {
              params.stateSpace = StateSpace::Q;
            }
            else
            {
              return -1;
            }
          }
          else if (lineIndex == 4)
          {
            params.ngSetSize = std::stoi(paramValueString);
          }
          else if (lineIndex == 5)
          {
            params.changeToLP = false;
            if (paramValueString == "Y")
            {
              params.changeToLP = true;
            }
          }
          else if (lineIndex == 6)
          {
            params.useSeparations = false;
            if (paramValueString == "Y")
            {
              params.useSeparations = true;
            }
          }
          else if (lineIndex == 8)
          {
            params.useVariableFixing = false;
            if (paramValueString == "Y")
            {
              params.useVariableFixing = true;
            }
          }
          else if (lineIndex == 9)
          {
            params.useMuSSP = false;
            if (paramValueString == "Y")
            {
              params.useMuSSP = true;
            }
          }
          else if (lineIndex == 10)
          {
            params.repairDuals = false;
            if (paramValueString == "Y")
            {
              params.repairDuals = true;
            }
          }
          else if (lineIndex == 12)
          {
            params.useRobustCuts = false;
            if (paramValueString == "Y")
            {
              params.useRobustCuts = true;
            }
          }
          else if (lineIndex == 13)
          {
            params.useNonRobustCuts = false;
            if (paramValueString == "Y")
            {
              params.useNonRobustCuts = true;
            }
          }
          else if (lineIndex == 14)
          {
            params.useVolumeAlgorithm = false;
            if (paramValueString == "Y")
            {
              params.useVolumeAlgorithm = true;
            }
          }
          else if (lineIndex == 15)
          {
            params.useScaling = false;
            if (paramValueString == "Y")
            {
              params.useScaling = true;
            }
          }
          else if (lineIndex == 16)
          {
            params.useSparseRCCs = false;
            if (paramValueString == "Y")
            {
              params.useSparseRCCs = true;
            }
          }
          else if (lineIndex == 17)
          {
            params.useRestarts = false;
            if (paramValueString == "Y")
            {
              params.useRestarts = true;
            }
          }
          else if (lineIndex == 18)
          {
            params.switchSepToCuts = false;
            if (paramValueString == "Y")
            {
              params.switchSepToCuts = true;
            }
          }
          else if (lineIndex == 20)
          {
            params.primalHeuristic = PrimalHeuristic::BEST_UB;
            if (paramValueString == "GREEDY")
            {
              params.primalHeuristic = PrimalHeuristic::GREEDY;
            }
            else if (paramValueString == "MIP")
            {
              params.primalHeuristic = PrimalHeuristic::MIP;
            }
          }
          else if (lineIndex == 21)
          {
            params.primalHeuristicLNS = PrimalHeuristicLNS::NONE;
            if (paramValueString == "LOCAL_SEARCH")
            {
              params.primalHeuristicLNS = PrimalHeuristicLNS::LOCAL_SEARCH;
            }
            else if (paramValueString == "DESTROY_REPAIR")
            {
              params.primalHeuristicLNS = PrimalHeuristicLNS::DESTROY_REPAIR;
            }
            else if (paramValueString == "BEAM")
            {
              params.primalHeuristicLNS = PrimalHeuristicLNS::BEAM;
            }
          }
          else if (lineIndex == 22)
          {
            params.lnsTimeout = std::stoi(paramValueString);
          }
          else if (lineIndex == 23)
          {
            params.lnsRandomPercent = std::stoi(paramValueString);
          }
          else if (lineIndex == 24)
          {
            params.lnsIterationsToUpdateParams = std::stoi(paramValueString);
          }
        }
        wordIndex = wordIndex + 1;
      }
      lineIndex = lineIndex + 1;
    }

    VRPTWDDSolver ddSolver(vrptw, params);
    bool solved = ddSolver.solve(true);
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
