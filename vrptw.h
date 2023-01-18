#ifndef VRPTW_H
#define VRPTW_H

#define DBG(x)

#include <string>
#include <cmath>
#include <regex>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

const long INF = 1e10;

enum InitialStateSpace
{
  NG = 0,
  Q = 1
};

enum PricingProblemType
{
  DD = 0,
  DP = 1
};

enum LPSolveType
{
  LPSolver = 0,
  LAGSolver = 1
};

enum DecompositionReason
{
  SEPARATE = 0,
  DECOMPOSE = 1
};

enum IncludeCoverConstraints
{
  Y = 0,
  N = 1
};

enum UseColumnGeneration
{
  USE_CG = 0,
  NO_CG = 1
};

enum FlowType
{
  LP = 0,
  IP = 1
};

enum AllowMultipleVisits
{
  YES = 0,
  NO = 1
};

enum ShortestPathMode
{
  UPDATE_POTENTIALS = 0,
  SHORTEST_PATH = 1
};

struct VRPTWSolution
{
  VRPTWSolution(std::vector<std::vector<int>> _routes, double _totalDistance) : routes(_routes), totalDistance(_totalDistance) {}

  std::vector<std::vector<int>> routes;
  double totalDistance;
};

const static std::map<std::string,double> instanceOptimalSolutions =
{
{"X-n979-k58.vrp",118976}
};

struct VRPTW
{
  public:
    VRPTW(std::string fileName)
    {
      std::vector<std::pair<double,double> > coordinates;

      bool vehicleCapacitySection = false;
      bool customerSection = false;
      std::ifstream infile(fileName);
      std::string line;
      while (std::getline(infile, line))
      {
        std::cout << line << std::endl;
        if (line.empty())
        {
          continue;
        }

        if (line.find("NUMBER") != std::string::npos)
        {
          vehicleCapacitySection = true;
          customerSection = false;
          continue;
        }
        else if (line.find("CUST NO.") != std::string::npos)
        {
          vehicleCapacitySection = false;
          customerSection = true;
          continue;
        }
 
        if (vehicleCapacitySection)
        {
          std::istringstream iss(line);
          int number, truckCapacity;
          if (!(iss >> number >> truckCapacity))
          {
            continue;
          }
          else
          {
            capacity = truckCapacity;
          }
        }


        if (customerSection)
        {
          std::istringstream iss(line);
          int custNo, x, y, demand, startTime, endTime, serviceTime;
          if (!(iss >> custNo >> x >> y >> demand >> startTime >> endTime >> serviceTime))
          {
            continue;
          }
          else
          {
            demands.push_back(demand);
            demandsForSeparation.push_back(demand);
            demandsForCombs.push_back(demand);

            startTimes.push_back(startTime);
            endTimes.push_back(endTime);
            serviceTimes.push_back(serviceTime);
            coordinates.push_back(std::make_pair(x, y));
          }
        }
      }

      distances.resize(demands.size());
      for (int i=0; i < demands.size(); ++i)
      {
        distances[i].resize(demands.size());
      }

      for (int i=0; i < demands.size(); ++i)
      {
        for (int j=0; j < demands.size(); ++j)
        {
          double distance = (int)(std::sqrt(std::pow((coordinates[i].first - coordinates[j].first),2) + std::pow((coordinates[i].second - coordinates[j].second),2))+0.5);
          distances[i][j] = distance;
          distances[j][i] = distance;
        }
      }

      numLocations = demands.size();

      // formally add if this works
      std::string instanceName = fileName.substr(fileName.find_last_of("/") + 1);
      if (instanceOptimalSolutions.find(instanceName) != instanceOptimalSolutions.end())
      {
        hgsUpperBound = instanceOptimalSolutions.find(instanceName)->second;
      }
      else
      {
        hgsUpperBound = INF;
      }
      std::cout << "hgs upper bound: " << hgsUpperBound << std::endl;
    };

    double evaluateSolutionCost(const std::vector<std::vector<int>>& routesByLocation)
    {
      double cost = 0;
      for (auto route : routesByLocation)
      {
        int previousLoc = 0;
        for (int loc : route)
        {
          cost = cost + distances[loc][previousLoc];
          previousLoc = loc;
        }
      }

      return cost;
    };

    int getLongestPossibleRoute()
    {
      std::vector<int> demandsCopy;
      for (int demand : demands)
      {
        demandsCopy.push_back(demand);
      }

      int totalDemand = 0;
      std::sort(demandsCopy.begin(), demandsCopy.end());
      for (int index=0; index<demandsCopy.size(); ++index)
      {
        totalDemand = totalDemand + demandsCopy[index];
        if (totalDemand > capacity)
        {
          return index;
        }
      }

      return demandsCopy.size();
    };

    int isRouteInfeasible(const std::vector<int>& routeByLocation)
    {
      std::vector<int> locationsVisitedFreq;
      locationsVisitedFreq.resize(numLocations);
      for (int location : routeByLocation)
      {
        ++locationsVisitedFreq[location];
        if (locationsVisitedFreq[location] == 2)
        {
          if (location != 0)
          {
            return location;
          }
        }
      }

      return -1;
    }

    int capacity;
    std::vector<std::vector<double> > distances;
    std::vector<int> demands;
    std::vector<double> demandsForSeparation;
    std::vector<int> demandsForCombs;

    std::vector<int> startTimes;
    std::vector<int> endTimes;
    std::vector<int> serviceTimes;

    int depot;
    int numLocations;

    double hgsUpperBound;
};

#endif
