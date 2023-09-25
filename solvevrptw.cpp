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
            if (paramValueString == "LAG")
            {
              params.lpSolveType = LPSolveType::LAGSolver;
            }
            else if (paramValueString == "LP")
            {
              params.lpSolveType = LPSolveType::LPSolver;
            }
            else
            {
              return -1;
            }
          }
          else if (lineIndex == 1)
          {
            if (paramValueString == "NG")
            {
              params.initialStateSpace = InitialStateSpace::NG;
            }
            else if (paramValueString == "Q")
            {
              params.initialStateSpace = InitialStateSpace::Q;
            }
            else
            {
              return -1;
            }
          }
          else if (lineIndex == 2)
          {
            params.s = std::stoi(paramValueString);
          }
          else if (lineIndex == 3)
          {
            params.maxS = std::stoi(paramValueString);
          }
          else if (lineIndex == 4)
          {
            params.useCuts = false;
            if (paramValueString == "Y")
            {
              params.useCuts = true;
            }
          }
          else if (lineIndex == 5)
          {
            params.useVariableFixing = false;
            if (paramValueString == "Y")
            {
              params.useVariableFixing = true;
            }
          }
          else if (lineIndex == 6)
          {
            params.useMuSSP = false;
            if (paramValueString == "Y")
            {
              params.useMuSSP = true;
            }
          }
          else if (lineIndex == 7)
          {
            params.timeoutSeconds = std::stoi(paramValueString);
          }
          else if (lineIndex == 8)
          {
            params.infeasibleRoutesBatchSize = std::stoi(paramValueString);
          }
          else if (lineIndex == 9)
          {
            params.lagIterationDelayToStartSeparating = std::stod(paramValueString);
          }
          else if (lineIndex == 10)
          {
            params.lagOptimalityGapToStartRepairing = std::stod(paramValueString);
          }
          else if (lineIndex == 11)
          {
            params.percentFixedToChangeToCPLEX = std::stod(paramValueString);
          }
          else if (lineIndex == 12)
          {
            params.numArcsToChangeToCPLEX = std::stoi(paramValueString);
          }
          else if (lineIndex == 13)
          {
            params.numArcsToChangeToLAG = std::stoi(paramValueString);
          }
          else if (lineIndex == 14)
          {
            params.numLagItersForCuts = std::stoi(paramValueString);
          }
          else if (lineIndex == 15)
          {
            params.numLagCuts = std::stoi(paramValueString);
          }
          else if (lineIndex == 16)
          {
            params.cutPhase = false;
            if (paramValueString == "Y")
            {
              params.cutPhase = true;
            }
          }
          else if (lineIndex == 17)
          {
            params.deactivateCutValueThreshold = std::stod(paramValueString);
          }
          else if (lineIndex == 18)
          {
            params.deactivateCutIterThreshold = std::stoi(paramValueString);
          }
          else if (lineIndex == 19)
          {
            params.momentumBeta = std::stod(paramValueString);
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
