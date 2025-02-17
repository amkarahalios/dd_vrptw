#include "vrptw.h"

VRPTW::VRPTW(std::string _fileName)
{
  fileName = _fileName;

  // VRPTW instance format
  std::regex vrptwInstanceRegex("[^_]vrptw");
  std::smatch vrptwInstanceMatch;
  if (std::regex_search(fileName, vrptwInstanceMatch, vrptwInstanceRegex))
  {
    fixedNumPaths = FixedNumPaths::FLEXIBLE_NUM_PATHS;
    circuitOrPath = CircuitOrPath::CIRCUIT;
    problemType = ProblemType::TW;
    vrptwTimeWindowType = VRPTWTimeWindowType::TIME_WINDOWS;
    timeStateMultiplier = 10;
    finalTimeStateMultiplier = timeStateMultiplier;
    timeStateDiscretization = 1;
    loadDiscretization = 0;
    std::vector<std::pair<double,double> > coordinates;

    std::regex noCapacityRelaxRegex(".*HG.*C1|R1|RC1.*txt");
    std::smatch noCapacityRelaxMatch;
    vrptwCapacityType = VRPTWCapacityType::RELAX_CAPACITY;
    counterType = VRPTWCounterType::USE_COUNTER;
    if (std::regex_search(fileName, noCapacityRelaxMatch, noCapacityRelaxRegex))
    {
      std::cout << "do not relax capacity constraint" << std::endl;
      vrptwCapacityType = VRPTWCapacityType::NO_RELAX_CAPACITY;
      counterType = VRPTWCounterType::NO_USE_COUNTER;
    }
    else
    {
      std::cout << "relax capacity constraint" << std::endl;
    }

    bool vehicleCapacitySection = false;
    bool customerSection = false;
    std::ifstream infile(fileName);
    std::string line;
    while (std::getline(infile, line))
    {
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
        int numberOfVehicles, truckCapacity;
        if (!(iss >> numberOfVehicles >> truckCapacity))
        {
          continue;
        }
        else
        {
          numVehicles = numberOfVehicles;
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

    precedences.resize(demands.size());
    distances.resize(demands.size());
    for (int i=0; i < demands.size(); ++i)
    {
      distances[i].resize(demands.size());
    }

    for (int i=0; i < demands.size(); ++i)
    {
      for (int j=0; j < demands.size(); ++j)
      {
        double distance = std::sqrt(std::pow((coordinates[i].first - coordinates[j].first),2) + std::pow((coordinates[i].second - coordinates[j].second),2));
        distance = (int)( timeStateMultiplier * distance ) / (timeStateMultiplier * 1.0);
        distances[i][j] = distance;
        distances[j][i] = distance;
      }
    }

    int minServiceTime = *std::min_element(serviceTimes.begin()+1, serviceTimes.end());
    int minDistanceDepot = (*std::min_element(distances[0].begin()+1, distances[0].end())) + 1;
    timeStateDiscretization = std::max(timeStateDiscretization, minServiceTime);

    numLocations = demands.size();
    std::cout << "num locations: " << numLocations << std::endl;
    routeLengthUpperBound = numLocations;

    std::string instanceName = fileName.substr(fileName.find_last_of("/") + 1);
    if (instanceOptimalSolutions.find(instanceName) != instanceOptimalSolutions.end())
    {
      instanceUpperBound = instanceOptimalSolutions.find(instanceName)->second;
      if (closedInstances.find(instanceName) != closedInstances.end())
      {
        isInstanceClosed = true;
      }
      else
      {
        isInstanceClosed = false;
      }
    }
    else
    {
      instanceUpperBound = INF;
    }
    std::cout << "instance upper bound: " << instanceUpperBound << std::endl;
  }

  // CVRP instance format, relax num trucks constraint
  std::regex cvrpInstanceRegex("([A-Z]).*-k([0-9]*).*");
  std::smatch cvrpInstanceMatch;
  if (std::regex_search(fileName, cvrpInstanceMatch, cvrpInstanceRegex))
  {
    timeStateMultiplier = 1;
    finalTimeStateMultiplier = timeStateMultiplier;
    timeStateDiscretization = 10;
    loadDiscretization = 0;

    if (cvrpInstanceMatch[1] == "X")
    {
      fixedNumPaths = FixedNumPaths::FLEXIBLE_NUM_PATHS;
      std::cout << "flexible number of vehicles" << std::endl;
    }
    else
    {
      fixedNumPaths = FixedNumPaths::FIXED_NUM_PATHS;
      std::cout << "fixed number of vehicles" << std::endl;
    }
    numVehicles = std::stoi(cvrpInstanceMatch[2]);
    std::cout << "num vehicles: " << numVehicles << std::endl;

    circuitOrPath = CircuitOrPath::CIRCUIT;
    problemType = ProblemType::CVRP;
    vrptwCapacityType = VRPTWCapacityType::NO_RELAX_CAPACITY;
    counterType = VRPTWCounterType::NO_USE_COUNTER;
    vrptwTimeWindowType = VRPTWTimeWindowType::NO_TIME_WINDOWS;
    std::vector<std::pair<double,double> > coordinates;

    bool nodeCoordSection = false;
    bool demandSection = false;
    bool depotSection = false;
    std::ifstream infile(fileName);
    std::string line;
    while (std::getline(infile, line))
    {
      if (line.find("NODE_COORD_SECTION") != std::string::npos)
      {
        nodeCoordSection = true;
        demandSection = false;
        depotSection = false;
        continue;
      }
      else if (line.find("DEMAND_SECTION") != std::string::npos)
      {
        nodeCoordSection = false;
        demandSection = true;
        depotSection = false;
        continue;
      }
      else if (line.find("DEPOT_SECTION") != std::string::npos)
      {
        nodeCoordSection = false;
        demandSection = false;
        depotSection = true;
        continue;
      }
      else if (line.find("CAPACITY : ") != std::string::npos)
      {
        std::regex capacityRegex("CAPACITY :[ \t]+([0-9]*)");
        std::smatch match;
        if (std::regex_search(line, match, capacityRegex))
        {
          capacity = std::stoi(match[1]);
          std::cout << "capacity: " << capacity << std::endl;
        }
      }

      if (nodeCoordSection)
      {
        std::istringstream iss(line);
        double a, b, c;
        if (!(iss >> a >> b >> c))
        {
          break;
        }
        else
        {
          coordinates.push_back(std::make_pair(b, c));
        }
      }
 
      if (demandSection)
      {
        std::istringstream iss(line);
        int a, b;
        if (!(iss >> a >> b))
        {
          break;
        }
        else
        {
          demands.push_back(b);
          demandsForSeparation.push_back(b);
          demandsForCombs.push_back(b);

          // cvrp no constraint
          startTimes.push_back(0);
          endTimes.push_back(100000);
          serviceTimes.push_back(0);
        }
      }

      if (depotSection)
      {
        std::istringstream iss(line);
        int a;
        if (!(iss >> a))
        {
          break;
        }
        else
        {
          depot = a - 1;
          depotSection = false;
        }
      }
    }

    precedences.resize(demands.size());
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
    std::cout << "num locations: " << numLocations << std::endl;
    routeLengthUpperBound = numLocations;

    DBG(
      std::cout << "distances: " << std::endl;
      for (int loc1=0; loc1<numLocations; ++loc1)
      {
        for (int loc2=0; loc2<numLocations; ++loc2)
        {
          std::cout << "[" << loc1 << "][" << loc2 << "]: " << distances[loc1][loc2] << std::endl;
        }
      }

      std::cout << "demands: " << std::endl;
      for (int loc=0; loc<numLocations; ++loc)
      {
        std::cout << "[" << loc << "]: " << demands[loc] << std::endl;
      }
    )

    std::string instanceName = fileName.substr(fileName.find_last_of("/") + 1);
    if (instanceOptimalSolutions.find(instanceName) != instanceOptimalSolutions.end())
    {
      instanceUpperBound = instanceOptimalSolutions.find(instanceName)->second;
      if (closedInstances.find(instanceName) != closedInstances.end())
      {
        isInstanceClosed = true;
      }
      else
      {
        isInstanceClosed = false;
      }
    }
    else
    {
      instanceUpperBound = INF;
    }
    std::cout << "instance upper bound: " << instanceUpperBound << std::endl;
  }

  // TSPTW instance format
  std::regex tsptwInstanceRegex(".*tsptw");
  std::smatch tsptwInstanceMatch;
  if (std::regex_search(fileName, tsptwInstanceMatch, tsptwInstanceRegex))
  {
    // afg, dumas, gendreau-dumas, ohlmann thomas, integer
    timeStateMultiplier = 1;
    finalTimeStateMultiplier = timeStateMultiplier;
    timeStateDiscretization = 1;
    loadDiscretization = 0;

    // solomon potvin bengio, solomon pesant, 0.0000
    std::regex solomonRegex(".*Solomon.*");
    std::smatch solomonMatch;
    if (std::regex_search(fileName, solomonMatch, solomonRegex))
    {
      timeStateMultiplier = 10000;
      finalTimeStateMultiplier = timeStateMultiplier;
      timeStateDiscretization = 100;
    }

    // langevin 0.1
    std::regex langevinRegex(".*Langevin.*");
    std::smatch langevinMatch;
    if (std::regex_search(fileName, langevinMatch, langevinRegex))
    {
      timeStateMultiplier = 10;
      finalTimeStateMultiplier = timeStateMultiplier;
      timeStateDiscretization = 10;
    }

    fixedNumPaths = FixedNumPaths::FIXED_NUM_PATHS;
    numVehicles = 1;
    circuitOrPath = CircuitOrPath::CIRCUIT;
    problemType = ProblemType::TSPTW;
    std::vector<std::pair<double,double> > coordinates;
    vrptwCapacityType = VRPTWCapacityType::RELAX_CAPACITY;
    counterType = VRPTWCounterType::USE_COUNTER;
    vrptwTimeWindowType = VRPTWTimeWindowType::TIME_WINDOWS;

    bool numLocationsSection = true;
    bool distanceMatrixSection = false;
    int distanceMatrixRow = 0;
    bool timeWindowSection = false;
    int timeWindowRow = 0;
    depot = 0;
    capacity = 10000;
    std::ifstream infile(fileName);
    std::string line;
    while (std::getline(infile, line))
    {
      if (line.empty())
      {
        continue;
      }

      if (timeWindowSection)
      {
        std::istringstream iss(line);
        int startTime, endTime;
        if (!(iss >> startTime >> endTime))
        {
          continue;
        }
        else
        {
          startTimes.push_back(startTime);
          endTimes.push_back(endTime);
        }
        
        timeWindowRow = timeWindowRow + 1;
        if (timeWindowRow == numLocations)
        {
          timeWindowSection = false;
        }
      }

      if (distanceMatrixSection)
      {
        std::istringstream iss(line);
        int column = 0;
        double distance;
        while (iss >> distance)
        {
          distances[distanceMatrixRow][column] = distance;
          /*
          if ((distanceMatrixRow != column) && (distance == 0.0) && (distanceMatrixRow != 0))
          {
            std::cout << "0 distance will be bad for input, row: " << distanceMatrixRow << " col: " << column << std::endl;
            return;
          }
          */
          column = column + 1;

          demands.push_back(1);
          demandsForSeparation.push_back(1);
          demandsForCombs.push_back(1);
          serviceTimes.push_back(0);
        }

        distanceMatrixRow = distanceMatrixRow + 1;
        if (distanceMatrixRow == numLocations)
        {
          distanceMatrixSection = false;
          timeWindowSection = true;
        }
      }

      if (numLocationsSection)
      {
        std::istringstream iss(line);
        int numLocationsInput;
        if (!(iss >> numLocationsInput))
        {
          continue;
        }
        else
        {
          numLocations = numLocationsInput;
          routeLengthUpperBound = numLocations;
          std::cout << "num loc: " << numLocations << std::endl;
        }

        numLocationsSection = false;
        distanceMatrixSection = true;
        distances.resize(numLocations);
        precedences.resize(numLocations);
        for (int i=0; i < distances.size(); ++i)
        {
          distances[i].resize(distances.size());
        }
      }
    }

    std::string instanceName = fileName.substr(fileName.find_last_of("/") + 1);
    if (instanceOptimalSolutions.find(instanceName) != instanceOptimalSolutions.end())
    {
      instanceUpperBound = instanceOptimalSolutions.find(instanceName)->second;
      if (closedInstances.find(instanceName) != closedInstances.end())
      {
        isInstanceClosed = true;
      }
      else
      {
        isInstanceClosed = false;
      }

      // rounding in TSPTW
      instanceUpperBound += 0.01;
    }
    else
    {
      instanceUpperBound = INF;
    }
    std::cout << "instance upper bound: " << instanceUpperBound << std::endl;
  }

  // SOP instance format
  std::regex sopInstanceRegex(".*sop");
  std::smatch sopInstanceMatch;
  if (std::regex_search(fileName, sopInstanceMatch, sopInstanceRegex))
  {
    timeStateMultiplier = 1;
    finalTimeStateMultiplier = timeStateMultiplier;
    timeStateDiscretization = 1;
    loadDiscretization = 0;

    fixedNumPaths = FixedNumPaths::FIXED_NUM_PATHS;
    numVehicles = 1;
    circuitOrPath = CircuitOrPath::PATH;
    problemType = ProblemType::SOP;
    std::vector<std::pair<double,double> > coordinates;
    vrptwCapacityType = VRPTWCapacityType::RELAX_CAPACITY;
    counterType = VRPTWCounterType::USE_COUNTER;
    vrptwTimeWindowType = VRPTWTimeWindowType::NO_TIME_WINDOWS;

    bool numLocationsSection = true;
    bool distanceMatrixSection = false;
    int distanceMatrixRow = 0;
    depot = 0;
    capacity = 10000;
    std::ifstream infile(fileName);
    std::string line;
    while (std::getline(infile, line))
    {
      if (line.find("EDGE_WEIGHT_SECTION") != std::string::npos)
      {
        numLocationsSection = true;
        continue;
      }

      if (distanceMatrixSection)
      {
        std::istringstream iss(line);
        int column = 0;
        double distance;
        while (iss >> distance)
        {
          // -1 indicates that vertex j must precede vertex i
          // so based on vertices seen we know which can come next?
          // each vertex as a set that all needs to be in the state to be an arc
          if ((column != 0) && (distance == -1))
          {
            precedences[distanceMatrixRow].insert(column);
            reliances[column].insert(distanceMatrixRow);
          }

          distances[distanceMatrixRow][column] = distance;
          column = column + 1;

          demands.push_back(1);
          demandsForSeparation.push_back(1);
          demandsForCombs.push_back(1);
          serviceTimes.push_back(0);
          startTimes.push_back(0);
          endTimes.push_back(1000000);
        }

        distanceMatrixRow = distanceMatrixRow + 1;
        if (distanceMatrixRow == numLocations)
        {
          distanceMatrixSection = false;
        }
      }

      if (numLocationsSection)
      {
        std::istringstream iss(line);
        int numLocationsInput;
        if (!(iss >> numLocationsInput))
        {
          continue;
        }
        else
        {
          numLocations = numLocationsInput;
          routeLengthUpperBound = numLocations;
          std::cout << "num loc: " << numLocations << std::endl;
        }

        numLocationsSection = false;
        distanceMatrixSection = true;
        distances.resize(numLocations);
        for (int i=0; i < distances.size(); ++i)
        {
          distances[i].resize(distances.size());
        }
        precedences.resize(numLocations);
        reliances.resize(numLocations);
      }
    }

    DBG(for (int i=0; i<precedences.size(); ++i)
    {
      for (auto j : precedences[i])
      {
        std::cout << i << "," << j << std::endl;
      }
    })

    std::string instanceName = fileName.substr(fileName.find_last_of("/") + 1);
    if (instanceOptimalSolutions.find(instanceName) != instanceOptimalSolutions.end())
    {
      instanceUpperBound = instanceOptimalSolutions.find(instanceName)->second;
      if (closedInstances.find(instanceName) != closedInstances.end())
      {
        isInstanceClosed = true;
      }
      else
      {
        isInstanceClosed = false;
      }
    }
    else
    {
      instanceUpperBound = INF;
    }
    std::cout << "instance upper bound: " << instanceUpperBound << std::endl;
  }

  // PDP instance format
  std::regex pdpInstanceRegex(".*pdp");
  std::smatch pdpInstanceMatch;
  if (std::regex_search(fileName, pdpInstanceMatch, pdpInstanceRegex))
  {
    //timeStateMultiplier = 1;
    //timeStateDiscretization = 10;
    //timeStateMultiplier = 10;
    //timeStateDiscretization = 1;
    timeStateMultiplier = 10;
    timeStateDiscretization = 10;
    //timeStateMultiplier = 1000;
    //timeStateDiscretization = 100;
    //timeStateMultiplier = 10000;
    //timeStateDiscretization = 1000;
    finalTimeStateMultiplier = 1000000;
    loadDiscretization = 1;

    fixedNumPaths = FixedNumPaths::FLEXIBLE_NUM_PATHS;
    circuitOrPath = CircuitOrPath::CIRCUIT;
    problemType = ProblemType::PDP;
    std::vector<std::pair<double,double> > coordinates;
    vrptwCapacityType = VRPTWCapacityType::RELAX_CAPACITY;
    counterType = VRPTWCounterType::USE_COUNTER;
    vrptwTimeWindowType = VRPTWTimeWindowType::TIME_WINDOWS;

    bool dataSection = true;
    bool locationDataSection = false;
    depot = 0;
    precedences.resize(2000);
    reliances.resize(2000);
    std::ifstream infile(fileName);
    std::string line;
    while (std::getline(infile, line))
    {
      if (line.empty())
      {
        continue;
      }

      if (dataSection)
      {
        std::istringstream iss(line);
        int capacityInput, speed;
        if (!(iss >> numVehicles >> capacityInput >> speed))
        {
          continue;
        }
        else
        {
          capacity = capacityInput;
        }
        dataSection = false;
        locationDataSection = true;
      }

      if (locationDataSection)
      {
        std::istringstream iss(line);
        int location, x, y, demand, tw1, tw2, serviceTime, pickup, delivery;
        while (iss >> location >> x >> y >> demand >> tw1 >> tw2 >> serviceTime >> pickup >> delivery)
        {
          coordinates.push_back(std::make_pair(x, y));
          demands.push_back(demand);
          demandsForSeparation.push_back(demand);
          demandsForCombs.push_back(demand);
          startTimes.push_back(tw1);
          endTimes.push_back(tw2);
          serviceTimes.push_back(serviceTime);
          if (location != 0)
          {
            if (pickup == 0)
            {
              precedences[delivery].insert(location);
              reliances[location].insert(delivery);
            }
            if (delivery == 0)
            {
              precedences[location].insert(pickup);
              reliances[pickup].insert(location);
            }
          }
        }
      }
    }

    std::cout << "num vehicles: " << numVehicles << std::endl;

    numLocations = demands.size();
    routeLengthUpperBound = numLocations;
    distances.resize(numLocations);
    preciseDistances.resize(numLocations);
    for (int i=0; i < distances.size(); ++i)
    {
      distances[i].resize(distances.size());
      preciseDistances[i].resize(distances.size());
    }
    for (int i=0; i < demands.size(); ++i)
    {
      for (int j=0; j < demands.size(); ++j)
      {
        double distance = std::sqrt(std::pow((coordinates[i].first - coordinates[j].first),2) + std::pow((coordinates[i].second - coordinates[j].second),2));
        distance = (int)( timeStateMultiplier * distance ) / (timeStateMultiplier * 1.00);
        distances[i][j] = distance;
        distances[j][i] = distance;
        double preciseDistance = std::sqrt(std::pow((coordinates[i].first - coordinates[j].first),2) + std::pow((coordinates[i].second - coordinates[j].second),2));
        preciseDistance = (int)( finalTimeStateMultiplier * preciseDistance ) / (finalTimeStateMultiplier * 1.00);
        preciseDistances[i][j] = preciseDistance;
        preciseDistances[j][i] = preciseDistance;
      }
    }

    int minServiceTime = 1000000;
    for (int loc1=0; loc1<numLocations; ++loc1)
    {
      if (serviceTimes[loc1] != 0)
      {
        minServiceTime = std::min(serviceTimes[loc1], minServiceTime);
      }
    }
    std::cout << "min service time: " << minServiceTime << std::endl;
    timeStateDiscretization = std::max(timeStateDiscretization, minServiceTime);

    std::string instanceName = fileName.substr(fileName.find_last_of("/") + 1);
    if (instanceOptimalSolutions.find(instanceName) != instanceOptimalSolutions.end())
    {
      instanceUpperBound = instanceOptimalSolutions.find(instanceName)->second;
      if (closedInstances.find(instanceName) != closedInstances.end())
      {
        isInstanceClosed = true;
      }
      else
      {
        isInstanceClosed = false;
      }

      // rounding in TSPTW
      //instanceUpperBound += 0.01;
    }
    else
    {
      instanceUpperBound = INF;
    }
    std::cout << "instance upper bound: " << instanceUpperBound << std::endl;
    infile.close();
  }

  calculateMaxValues();
  calculateRouteLengthUpperBound();
  std::cout << "route length upper bound: " << routeLengthUpperBound << std::endl;
};

double VRPTW::evaluateRouteDistance(const std::vector<int>& routeByLocation)
{
  double cost = 0;
  int previousLoc = 0;
  for (int loc : routeByLocation)
  {
    if (problemType == ProblemType::PDP)
    {
      cost = cost + preciseDistances[loc][previousLoc];
    }
    else
    {
      cost = cost + distances[loc][previousLoc];
    }
    previousLoc = loc;
  }

  return cost;
}

double VRPTW::evaluateSolutionCost(const std::vector<std::vector<int>>& routesByLocation)
{
  double cost = 0;
  for (auto route : routesByLocation)
  {
    cost += evaluateRouteDistance(route);
  }

  return cost;
};

int VRPTW::getLongestPossibleRoute()
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

void VRPTW::recomputeDistancesPDPTW()
{
  bool dataSection = true;
  bool locationDataSection = false;
  std::vector<std::pair<double,double> > coordinates;
  std::ifstream infile(fileName);
  std::string line;
  while (std::getline(infile, line))
  {
    if (line.empty())
    {
      continue;
    }

    if (dataSection)
    {
      std::istringstream iss(line);
      int capacityInput, speed;
      if (!(iss >> numVehicles >> capacityInput >> speed))
      {
        continue;
      }
      else
      {
        capacity = capacityInput;
      }
      dataSection = false;
      locationDataSection = true;
    }

    if (locationDataSection)
    {
      std::istringstream iss(line);
      int location, x, y, demand, tw1, tw2, serviceTime, pickup, delivery;
      while (iss >> location >> x >> y >> demand >> tw1 >> tw2 >> serviceTime >> pickup >> delivery)
      {
        coordinates.push_back(std::make_pair(x, y));
      }
    }
  }

  distances.clear();
  distances.resize(numLocations);
  for (int i=0; i < distances.size(); ++i)
  {
    distances[i].resize(distances.size());
  }
  for (int i=0; i < demands.size(); ++i)
  {
    for (int j=0; j < demands.size(); ++j)
    {
      double distance = std::sqrt(std::pow((coordinates[i].first - coordinates[j].first),2) + std::pow((coordinates[i].second - coordinates[j].second),2));
      distance = (int)( timeStateMultiplier * distance ) / (timeStateMultiplier * 1.00);
      distances[i][j] = distance;
      distances[j][i] = distance;
    }
  }
};

void VRPTW::calculateRouteLengthUpperBound()
{
  // capacity
  if ((problemType == ProblemType::TW) || (problemType == ProblemType::CVRP))
  {
    std::vector<int> sortedDemands = demands;
    std::sort(sortedDemands.begin(), sortedDemands.end());
    int greedyLocations = 0;
    int greedyLoad = 0;
    while (greedyLoad < capacity)
    {
      for (int index=0; index<sortedDemands.size(); ++index)
      {
        greedyLocations = greedyLocations + 1;
        greedyLoad = greedyLoad + sortedDemands[index];
        if (greedyLoad >= capacity)
        {
          break;
        }
      }
    }

    routeLengthUpperBound = std::min(routeLengthUpperBound, greedyLocations);
  }

  // distances
  if ((problemType == ProblemType::PDP) || (problemType == ProblemType::TW))
  {
    std::set<int> greedyLocations;
    int greedyTime = 0;
    while (greedyTime < endTimes[0])
    {
      double minDistance = 1e9;
      int minLoc = -1;
      for (int index1=0; index1<numLocations; ++index1)
      {
        if (greedyLocations.find(index1) != greedyLocations.end())
        {
          continue;
        }
        for (int index2=0; index2<numLocations; ++index2)
        {
          if (index1 == index2)
          {
            continue;
          }

          double distance = distances[index1][index2] + serviceTimes[index1];
          if (distance < minDistance)
          {
            minDistance = distance;
            minLoc = index1;
          }
        }
      }

      greedyLocations.insert(minLoc);
      greedyTime = greedyTime + minDistance;
    }

    routeLengthUpperBound = std::min(routeLengthUpperBound, (int)greedyLocations.size());
  }
}

void VRPTW::calculateMaxValues()
{
  maxDistance = -1;
  maxDemand = -1;
  maxStartTime = -1;
  for (int i=0; i<numLocations; ++i)
  {
    maxDemand = std::max(demands[i], maxDemand);
    maxStartTime = std::max(startTimes[i], maxStartTime);
    for (int j=0; j<numLocations; ++j)
    {
      maxDistance = std::max(distances[i][j], maxDistance);
    }
  }
}
