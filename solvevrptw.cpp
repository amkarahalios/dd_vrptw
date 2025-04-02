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

  //vrptw.writeVrpFormat();

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
            params.numLagCuts = std::stoi(paramValueString);
          }
          else if (lineIndex == 13)
          {
            params.useRCCs = false;
            if (paramValueString == "Y")
            {
              params.useRCCs = true;
            }
          }
          else if (lineIndex == 14)
          {
            params.useSRC3s = false;
            if (paramValueString == "Y")
            {
              params.useSRC3s = true;
            }
          }
          else if (lineIndex == 15)
          {
            params.useSRC4s = false;
            if (paramValueString == "Y")
            {
              params.useSRC4s = true;
            }
          }
          else if (lineIndex == 16)
          {
            params.useSRC5V1s = false;
            if (paramValueString == "Y")
            {
              params.useSRC5V1s = true;
            }
          }
          else if (lineIndex == 17)
          {
            params.useSRC5V2s = false;
            if (paramValueString == "Y")
            {
              params.useSRC5V2s = true;
            }
          }
          else if (lineIndex == 18)
          {
            params.useVolumeAlgorithm = false;
            if (paramValueString == "Y")
            {
              params.useVolumeAlgorithm = true;
            }
          }
          else if (lineIndex == 19)
          {
            params.limitRCCs = false;
            if (paramValueString == "Y")
            {
              params.limitRCCs = true;
            }
          }
          else if (lineIndex == 20)
          {
            params.useScaling = false;
            if (paramValueString == "Y")
            {
              params.useScaling = true;
            }
          }
          else if (lineIndex == 22)
          {
            params.primalHeuristic = PrimalHeuristic::STANDALONE;
            if (paramValueString == "CE_MIP")
            {
              params.primalHeuristic = PrimalHeuristic::CE_MIP;
            }
            else if (paramValueString == "CE_MIP_LNS")
            {
              params.primalHeuristic = PrimalHeuristic::CE_MIP_LNS;
            }
            else if (paramValueString == "CE_GREEDY")
            {
              params.primalHeuristic = PrimalHeuristic::CE_GREEDY;
            }
            else if (paramValueString == "BEST_KNOWN")
            {
              params.primalHeuristic = PrimalHeuristic::BEST_KNOWN;
            }
          }
          else if (lineIndex == 23)
          {
            params.lnsTimeoutSeconds = std::stoi(paramValueString);
          }
          else if (lineIndex == 24)
          {
            if (paramValueString == "RANDOM")
            {
              params.removalStrategy = PrimalHeuristicRemovalStrategy::RANDOM;
            }
            else if (paramValueString == "SHAW")
            {
              params.removalStrategy = PrimalHeuristicRemovalStrategy::SHAW;
            }
            else if (paramValueString == "WORST")
            {
              params.removalStrategy = PrimalHeuristicRemovalStrategy::WORST;
            }
            else if (paramValueString == "SEQUENCE_RANDOM")
            {
              params.removalStrategy = PrimalHeuristicRemovalStrategy::SEQUENCE_RANDOM;
            }
            else if (paramValueString == "SEQUENCE_SHAW")
            {
              params.removalStrategy = PrimalHeuristicRemovalStrategy::SEQUENCE_SHAW;
            }
            else if (paramValueString == "SEQUENCE_WORST")
            {
              params.removalStrategy = PrimalHeuristicRemovalStrategy::SEQUENCE_WORST;
            }
          }
          else if (lineIndex == 25)
          {
            if (paramValueString == "NONE")
            {
              params.insertionAblation = PrimalHeuristicInsertionAblation::NO_ABLATION;
            }
            else if (paramValueString == "INTRA_SWAP")
            {
              params.insertionAblation = PrimalHeuristicInsertionAblation::INTRA_SWAP;
            }
            else if (paramValueString == "TRUNCATED")
            {
              params.insertionAblation = PrimalHeuristicInsertionAblation::TRUNCATED;
            }
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
