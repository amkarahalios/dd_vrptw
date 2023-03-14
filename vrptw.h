#ifndef VRPTW_H
#define VRPTW_H

#define DBG(x)

#include <string>
#include <cmath>
#include <set>
#include <regex>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

const long INF = 1e10;

enum OneOrMorePaths
{
  ONE_PATH = 0,
  MORE_PATHS = 1
};

enum CircuitOrPath
{
  CIRCUIT = 0,
  PATH = 1
};

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

enum VRPTWCapacityType
{
  RELAX_CAPACITY = 0,
  NO_RELAX_CAPACITY = 1,
};

enum VRPTWTimeWindowType
{
  TIME_WINDOWS = 0,
  NO_TIME_WINDOWS = 1,
};

struct VRPTWSolution
{
  VRPTWSolution(std::vector<std::vector<int>> _routes, double _totalDistance) : routes(_routes), totalDistance(_totalDistance) {}

  std::vector<std::vector<int>> routes;
  double totalDistance;
};

const static std::map<std::string,double> instanceOptimalSolutions =
{
{"C101.vrptw",827.3},
{"C102.vrptw",827.3},
{"C103.vrptw",826.3},
{"C104.vrptw",822.9},
{"C105.vrptw",827.3},
{"C106.vrptw",827.3},
{"C107.vrptw",827.3},
{"C108.vrptw",827.3},
{"C109.vrptw",827.3},
{"C201.vrptw",589.1},
{"C202.vrptw",589.1},
{"C203.vrptw",588.7},
{"C204.vrptw",588.1},
{"C205.vrptw",586.4},
{"C206.vrptw",586},
{"C207.vrptw",585.8},
{"C208.vrptw",585.8},
{"R101.vrptw",1637.7},
{"R102.vrptw",1466.6},
{"R103.vrptw",1208.7},
{"R104.vrptw",971.5},
{"R105.vrptw",1355.3},
{"R106.vrptw",1234.6},
{"R107.vrptw",1064.6},
{"R108.vrptw",932.1},
{"R109.vrptw",1146.9},
{"R110.vrptw",1068},
{"R111.vrptw",1048.7},
{"R112.vrptw",948.6},
{"R201.vrptw",1143.2},
{"R202.vrptw",1029.6},
{"R203.vrptw",870.8},
{"R204.vrptw",731.3},
{"R205.vrptw",949.8},
{"R206.vrptw",875.9},
{"R207.vrptw",794},
{"R208.vrptw",701},
{"R209.vrptw",854.8},
{"R210.vrptw",900.5},
{"R211.vrptw",746.7},
{"RC101.vrptw",1619.8},
{"RC102.vrptw",1457.4},
{"RC103.vrptw",1258},
{"RC104.vrptw",1132.3},
{"RC105.vrptw",1513.7},
{"RC106.vrptw",1372.7},
{"RC107.vrptw",1207.8},
{"RC108.vrptw",1114.2},
{"RC201.vrptw",1261.8},
{"RC202.vrptw",1092.3},
{"RC203.vrptw",923.7},
{"RC204.vrptw",783.5},
{"RC205.vrptw",1154},
{"RC206.vrptw",1051.1},
{"RC207.vrptw",962.9},
{"RC208.vrptw",776.1},
{"C1_2_1.vrptw",2698.6},
{"C1_2_2.vrptw",2694.3},
{"C1_2_3.vrptw",2675.8},
{"C1_2_4.vrptw",2625.6},
{"C1_2_5.vrptw",2694.9},
{"C1_2_6.vrptw",2694.9},
{"C1_2_7.vrptw",2694.9},
{"C1_2_8.vrptw",2684},
{"C1_2_9.vrptw",2639.6},
{"C1_2_10.vrptw",2624.7},
{"C2_2_1.vrptw",1922.1},
{"C2_2_2.vrptw",1851.4},
{"C2_2_3.vrptw",1763.4},
{"C2_2_4.vrptw",1695},
{"C2_2_5.vrptw",1869.6},
{"C2_2_6.vrptw",1844.8},
{"C2_2_7.vrptw",1842.2},
{"C2_2_8.vrptw",1813.7},
{"C2_2_9.vrptw",1815},
{"C2_2_10.vrptw",1791.2},
{"R1_2_1.vrptw",4667.2},
{"R1_2_2.vrptw",3919.9},
{"R1_2_3.vrptw",3373.9},
{"R1_2_4.vrptw",3047.6},
{"R1_2_5.vrptw",4053.2},
{"R1_2_6.vrptw",3559.1},
{"R1_2_7.vrptw",3141.9},
{"R1_2_8.vrptw",2938.4},
{"R1_2_9.vrptw",3734.7},
{"R1_2_10.vrptw",3293.1},
{"R2_2_1.vrptw",3468},
{"R2_2_2.vrptw",3008.2},
{"R2_2_3.vrptw",2537.5},
{"R2_2_4.vrptw",1928.5},
{"R2_2_5.vrptw",3061.1},
{"R2_2_6.vrptw",2675.4},
{"R2_2_7.vrptw",2304.7},
{"R2_2_8.vrptw",1842.4},
{"R2_2_9.vrptw",2843.3},
{"R2_2_10.vrptw",2549.4},
{"RC1_2_1.vrptw",3516.9},
{"RC1_2_2.vrptw",3221.6},
{"RC1_2_3.vrptw",3001.4},
{"RC1_2_4.vrptw",2845.2},
{"RC1_2_5.vrptw",3325.6},
{"RC1_2_6.vrptw",3300.7},
{"RC1_2_7.vrptw",3177.8},
{"RC1_2_8.vrptw",3060},
{"RC1_2_9.vrptw",3073.3},
{"RC1_2_10.vrptw",2990.5},
{"RC2_2_1.vrptw",2797.4},
{"RC2_2_2.vrptw",2481.6},
{"RC2_2_3.vrptw",2227.7},
{"RC2_2_4.vrptw",1854.8},
{"RC2_2_5.vrptw",2491.4},
{"RC2_2_6.vrptw",2495.1},
{"RC2_2_7.vrptw",2287.7},
{"RC2_2_8.vrptw",2151.2},
{"RC2_2_9.vrptw",2086.6},
{"RC2_2_10.vrptw",1989.2},
{"C1_4_1.vrptw",7138.8},
{"C1_4_2.vrptw",7113.3},
{"C1_4_3.vrptw",6929.9},
{"C1_4_4.vrptw",6777.7},
{"C1_4_5.vrptw",7138.8},
{"C1_4_6.vrptw",7140.1},
{"C1_4_7.vrptw",7136.2},
{"C1_4_8.vrptw",7083},
{"C1_4_9.vrptw",6927.8},
{"C1_4_10.vrptw",6825.4},
{"C2_4_1.vrptw",4100.3},
{"C2_4_2.vrptw",3914.1},
{"C2_4_3.vrptw",3755.2},
{"C2_4_4.vrptw",3523.7},
{"C2_4_5.vrptw",3923.2},
{"C2_4_6.vrptw",3860.1},
{"C2_4_7.vrptw",3870.9},
{"C2_4_8.vrptw",3773.7},
{"C2_4_9.vrptw",3842.1},
{"C2_4_10.vrptw",3665.1},
{"R1_4_1.vrptw",10305.8},
{"R1_4_2.vrptw",8873.3},
{"R1_4_3.vrptw",7784.3},
{"R1_4_4.vrptw",7266.2},
{"R1_4_5.vrptw",9184.6},
{"R1_4_6.vrptw",8340.4},
{"R1_4_7.vrptw",7599.8},
{"R1_4_8.vrptw",7240.5},
{"R1_4_9.vrptw",8677.5},
{"R1_4_10.vrptw",8077.8},
{"R2_4_1.vrptw",7520.7},
{"R2_4_2.vrptw",6482.8},
{"R2_4_3.vrptw",5372.9},
{"R2_4_4.vrptw",4211.2},
{"R2_4_5.vrptw",6567.9},
{"R2_4_6.vrptw",5813.5},
{"R2_4_7.vrptw",4893.5},
{"R2_4_8.vrptw",4000.1},
{"R2_4_9.vrptw",6067.8},
{"R2_4_10.vrptw",5645.9},
{"RC1_4_1.vrptw",8522.9},
{"RC1_4_2.vrptw",7878.2},
{"RC1_4_3.vrptw",7516.9},
{"RC1_4_4.vrptw",7292.9},
{"RC1_4_5.vrptw",8152.3},
{"RC1_4_6.vrptw",8148},
{"RC1_4_7.vrptw",7932.5},
{"RC1_4_8.vrptw",7757.2},
{"RC1_4_9.vrptw",7717.7},
{"RC1_4_10.vrptw",7581.2},
{"RC2_4_1.vrptw",6147.3},
{"RC2_4_2.vrptw",5407.5},
{"RC2_4_3.vrptw",4573},
{"RC2_4_4.vrptw",3597.9},
{"RC2_4_5.vrptw",5392.3},
{"RC2_4_6.vrptw",5324.6},
{"RC2_4_7.vrptw",4987.8},
{"RC2_4_8.vrptw",4693.3},
{"RC2_4_9.vrptw",4510.4},
{"RC2_4_10.vrptw",4252.3},
{"C1_6_1.vrptw",14076.6},
{"C1_6_2.vrptw",13948.3},
{"C1_6_3.vrptw",13757},
{"C1_6_4.vrptw",13538.6},
{"C1_6_5.vrptw",14066.8},
{"C1_6_6.vrptw",14070.9},
{"C1_6_7.vrptw",14066.8},
{"C1_6_8.vrptw",13991.2},
{"C1_6_9.vrptw",13664.5},
{"C1_6_10.vrptw",13617.5},
{"C2_6_1.vrptw",7752.2},
{"C2_6_2.vrptw",7471.5},
{"C2_6_3.vrptw",7215},
{"C2_6_4.vrptw",6877},
{"C2_6_5.vrptw",7553.8},
{"C2_6_6.vrptw",7449.8},
{"C2_6_7.vrptw",7491.3},
{"C2_6_8.vrptw",7303.7},
{"C2_6_9.vrptw",7303.2},
{"C2_6_10.vrptw",7123.9},
{"R1_6_1.vrptw",21274.2},
{"R1_6_2.vrptw",18558.7},
{"R1_6_3.vrptw",16874.9},
{"R1_6_4.vrptw",15721.4},
{"R1_6_5.vrptw",19294.9},
{"R1_6_6.vrptw",17763.7},
{"R1_6_7.vrptw",16496.2},
{"R1_6_8.vrptw",15584.3},
{"R1_6_9.vrptw",18474.1},
{"R1_6_10.vrptw",17583.7},
{"R2_6_1.vrptw",15145.3},
{"R2_6_2.vrptw",12976.3},
{"R2_6_3.vrptw",10455.3},
{"R2_6_4.vrptw",7915.1},
{"R2_6_5.vrptw",13790.2},
{"R2_6_6.vrptw",11847.8},
{"R2_6_7.vrptw",9777.9},
{"R2_6_8.vrptw",7512.3},
{"R2_6_9.vrptw",12736.8},
{"R2_6_10.vrptw",11837},
{"RC1_6_1.vrptw",16960.1},
{"RC1_6_2.vrptw",15890.6},
{"RC1_6_3.vrptw",15181.3},
{"RC1_6_4.vrptw",14753.2},
{"RC1_6_5.vrptw",16536.3},
{"RC1_6_6.vrptw",16473.3},
{"RC1_6_7.vrptw",16055.3},
{"RC1_6_8.vrptw",15891.8},
{"RC1_6_9.vrptw",15803.5},
{"RC1_6_10.vrptw",15651.3},
{"RC2_6_1.vrptw",11966.1},
{"RC2_6_2.vrptw",10336.9},
{"RC2_6_3.vrptw",8894.9},
{"RC2_6_4.vrptw",6967.5},
{"RC2_6_5.vrptw",11080.7},
{"RC2_6_6.vrptw",10830.5},
{"RC2_6_7.vrptw",10289.4},
{"RC2_6_8.vrptw",9779},
{"RC2_6_9.vrptw",9436},
{"RC2_6_10.vrptw",8974.7},
{"C1_8_1.vrptw",25156.9},
{"C1_8_2.vrptw",24974.1},
{"C1_8_3.vrptw",24156.1},
{"C1_8_4.vrptw",23797.3},
{"C1_8_5.vrptw",25138.6},
{"C1_8_6.vrptw",25133.3},
{"C1_8_7.vrptw",25127.3},
{"C1_8_8.vrptw",24809.7},
{"C1_8_9.vrptw",24200.4},
{"C1_8_10.vrptw",24026.7},
{"C2_8_1.vrptw",11631.9},
{"C2_8_2.vrptw",11394.5},
{"C2_8_3.vrptw",11138.1},
{"C2_8_4.vrptw",10650},
{"C2_8_5.vrptw",11395.6},
{"C2_8_6.vrptw",11316.3},
{"C2_8_7.vrptw",11332.9},
{"C2_8_8.vrptw",11133.9},
{"C2_8_9.vrptw",11140.4},
{"C2_8_10.vrptw",10946},
{"R1_8_1.vrptw",36345},
{"R1_8_2.vrptw",32277.6},
{"R1_8_3.vrptw",29304.5},
{"R1_8_4.vrptw",27734.7},
{"R1_8_5.vrptw",33494.2},
{"R1_8_6.vrptw",30872.4},
{"R1_8_7.vrptw",28789},
{"R1_8_8.vrptw",27609.4},
{"R1_8_9.vrptw",32257.3},
{"R1_8_10.vrptw",30918.4},
{"R2_8_1.vrptw",24969.8},
{"R2_8_2.vrptw",21312.2},
{"R2_8_3.vrptw",17234.8},
{"R2_8_4.vrptw",13160.8},
{"R2_8_5.vrptw",22801.6},
{"R2_8_6.vrptw",19740.5},
{"R2_8_7.vrptw",16357.5},
{"R2_8_8.vrptw",12611.7},
{"R2_8_9.vrptw",21282.7},
{"R2_8_10.vrptw",19984.8},
{"RC1_8_1.vrptw",29978.9},
{"RC1_8_2.vrptw",28290.1},
{"RC1_8_3.vrptw",27447.7},
{"RC1_8_4.vrptw",26557.2},
{"RC1_8_5.vrptw",29219.9},
{"RC1_8_6.vrptw",29194.2},
{"RC1_8_7.vrptw",28788.6},
{"RC1_8_8.vrptw",28418.1},
{"RC1_8_9.vrptw",28347.1},
{"RC1_8_10.vrptw",28168.5},
{"RC2_8_1.vrptw",19201.3},
{"RC2_8_2.vrptw",16709.5},
{"RC2_8_3.vrptw",14013.6},
{"RC2_8_4.vrptw",10969.4},
{"RC2_8_5.vrptw",17466.1},
{"RC2_8_6.vrptw",17195.1},
{"RC2_8_7.vrptw",16362.2},
{"RC2_8_8.vrptw",15528.8},
{"RC2_8_9.vrptw",15183},
{"RC2_8_10.vrptw",14370.9},
{"C1_10_1.vrptw",42444.8},
{"C1_10_2.vrptw",41337.8},
{"C1_10_3.vrptw",40064.4},
{"C1_10_4.vrptw",39434.1},
{"C1_10_5.vrptw",42434.8},
{"C1_10_6.vrptw",42437},
{"C1_10_7.vrptw",42420.4},
{"C1_10_8.vrptw",41652.1},
{"C1_10_9.vrptw",40288.4},
{"C1_10_10.vrptw",39816.8},
{"C2_10_1.vrptw",16841.1},
{"C2_10_2.vrptw",16462.6},
{"C2_10_3.vrptw",16036.5},
{"C2_10_4.vrptw",15459.5},
{"C2_10_5.vrptw",16521.3},
{"C2_10_6.vrptw",16290.7},
{"C2_10_7.vrptw",16378.4},
{"C2_10_8.vrptw",16029.1},
{"C2_10_9.vrptw",16075.4},
{"C2_10_10.vrptw",15728.6},
{"R1_10_1.vrptw",53046.5},
{"R1_10_2.vrptw",48263.1},
{"R1_10_3.vrptw",44677.1},
{"R1_10_4.vrptw",42440.7},
{"R1_10_5.vrptw",50406.7},
{"R1_10_6.vrptw",46930.3},
{"R1_10_7.vrptw",43997.4},
{"R1_10_8.vrptw",42279.3},
{"R1_10_9.vrptw",49162.8},
{"R1_10_10.vrptw",47364.6},
{"R2_10_1.vrptw",36881},
{"R2_10_2.vrptw",31241.9},
{"R2_10_3.vrptw",24399},
{"R2_10_4.vrptw",17811.5},
{"R2_10_5.vrptw",34132.8},
{"R2_10_6.vrptw",29124.7},
{"R2_10_7.vrptw",23102.2},
{"R2_10_8.vrptw",17403.8},
{"R2_10_9.vrptw",31990.6},
{"R2_10_10.vrptw",29840.5},
{"RC1_10_1.vrptw",45790.8},
{"RC1_10_2.vrptw",43678.3},
{"RC1_10_3.vrptw",42122},
{"RC1_10_4.vrptw",41357.4},
{"RC1_10_5.vrptw",45028.1},
{"RC1_10_6.vrptw",44903.6},
{"RC1_10_7.vrptw",44417.1},
{"RC1_10_8.vrptw",43916.5},
{"RC1_10_9.vrptw",43858.1},
{"RC1_10_10.vrptw",43533.7},
{"RC2_10_1.vrptw",28122.6},
{"RC2_10_2.vrptw",24248.6},
{"RC2_10_3.vrptw",19618.1},
{"RC2_10_4.vrptw",15657},
{"RC2_10_5.vrptw",25797.5},
{"RC2_10_6.vrptw",25782.5},
{"RC2_10_7.vrptw",24395.8},
{"RC2_10_8.vrptw",23280.2},
{"RC2_10_9.vrptw",22731.6},
{"RC2_10_10.vrptw",21736.1},
{"A-n32-k5.vrp",784},
{"A-n33-k5.vrp",661},
{"A-n33-k6.vrp",742},
{"A-n34-k5.vrp",778},
{"A-n36-k5.vrp",799},
{"A-n37-k5.vrp",669},
{"A-n37-k6.vrp",949},
{"A-n38-k5.vrp",730},
{"A-n39-k5.vrp",822},
{"A-n39-k6.vrp",831},
{"A-n44-k6.vrp",937},
{"A-n45-k6.vrp",944},
{"A-n45-k7.vrp",1146},
{"A-n46-k7.vrp",914},
{"A-n48-k7.vrp",1073},
{"A-n53-k7.vrp",1010},
{"A-n54-k7.vrp",1167},
{"A-n55-k9.vrp",1073},
{"A-n60-k9.vrp",1354},
{"A-n61-k9.vrp",1034},
{"A-n62-k8.vrp",1288},
{"A-n63-k10.vrp",1314},
{"A-n63-k9.vrp",1616},
{"A-n64-k9.vrp",1401},
{"A-n65-k9.vrp",1174},
{"A-n69-k9.vrp",1159},
{"A-n80-k10.vrp",1763},
{"B-n31-k5.vrp",672},
{"B-n34-k5.vrp",788},
{"B-n35-k5.vrp",955},
{"B-n38-k6.vrp",805},
{"B-n39-k5.vrp",549},
{"B-n41-k6.vrp",829},
{"B-n43-k6.vrp",742},
{"B-n44-k7.vrp",909},
{"B-n45-k5.vrp",751},
{"B-n45-k6.vrp",678},
{"B-n50-k7.vrp",741},
{"B-n50-k8.vrp",1312},
{"B-n51-k7.vrp",1032},
{"B-n52-k7.vrp",747},
{"B-n56-k7.vrp",707},
{"B-n57-k7.vrp",1153},
{"B-n57-k9.vrp",1598},
{"B-n63-k10.vrp",1496},
{"B-n64-k9.vrp",861},
{"B-n66-k9.vrp",1316},
{"B-n67-k10.vrp",1032},
{"B-n68-k9.vrp",1272},
{"B-n78-k10.vrp",1221},
{"E-n101-k14.vrp",1071},
{"E-n101-k8.vrp",817},
{"E-n13-k4.vrp",247},
{"E-n22-k4.vrp",375},
{"E-n23-k3.vrp",569},
{"E-n30-k3.vrp",534},
{"E-n31-k7.vrp",379},
{"E-n33-k4.vrp",835},
{"E-n51-k5.vrp",521},
{"E-n76-k10.vrp",830},
{"E-n76-k14.vrp",1021},
{"E-n76-k7.vrp",682},
{"E-n76-k8.vrp",735},
{"F-n135-k7.vrp",1162},
{"F-n45-k4.vrp",724},
{"F-n72-k4.vrp",237},
{"P-n101-k4.vrp",681},
{"P-n16-k8.vrp",450},
{"P-n19-k2.vrp",212},
{"P-n20-k2.vrp",216},
{"P-n21-k2.vrp",211},
{"P-n22-k2.vrp",216},
{"P-n22-k8.vrp",603},
{"P-n23-k8.vrp",529},
{"P-n40-k5.vrp",458},
{"P-n45-k5.vrp",510},
{"P-n50-k10.vrp",696},
{"P-n50-k7.vrp",554},
{"P-n50-k8.vrp",631},
{"P-n51-k10.vrp",741},
{"P-n55-k10.vrp",694},
{"P-n55-k15.vrp",989},
{"P-n55-k7.vrp",568},
{"P-n55-k8.vrp",588},
{"P-n60-k10.vrp",744},
{"P-n60-k15.vrp",968},
{"P-n65-k10.vrp",792},
{"P-n70-k10.vrp",827},
{"P-n76-k4.vrp",593},
{"P-n76-k5.vrp",627},
{"M-n101-k10.vrp",820},
{"M-n121-k7.vrp",1034},
{"M-n151-k12.vrp",1053},
{"M-n200-k16.vrp",1373},
{"M-n200-k17.vrp",1373},
{"X-n1001-k43.vrp",72355},
{"X-n101-k25.vrp",27591},
{"X-n106-k14.vrp",26362},
{"X-n110-k13.vrp",14971},
{"X-n115-k10.vrp",12747},
{"X-n120-k6.vrp",13332},
{"X-n125-k30.vrp",55539},
{"X-n129-k18.vrp",28940},
{"X-n134-k13.vrp",10916},
{"X-n139-k10.vrp",13590},
{"X-n143-k7.vrp",15700},
{"X-n148-k46.vrp",43448},
{"X-n153-k22.vrp",21220},
{"X-n157-k13.vrp",16876},
{"X-n162-k11.vrp",14138},
{"X-n167-k10.vrp",20557},
{"X-n172-k51.vrp",45607},
{"X-n176-k26.vrp",47812},
{"X-n181-k23.vrp",25569},
{"X-n186-k15.vrp",24145},
{"X-n190-k8.vrp",16980},
{"X-n195-k51.vrp",44225},
{"X-n200-k36.vrp",58578},
{"X-n204-k19.vrp",19565},
{"X-n209-k16.vrp",30656},
{"X-n214-k11.vrp",10856},
{"X-n219-k73.vrp",117595},
{"X-n223-k34.vrp",40437},
{"X-n228-k23.vrp",25742},
{"X-n233-k16.vrp",19230},
{"X-n237-k14.vrp",27042},
{"X-n242-k48.vrp",82751},
{"X-n247-k50.vrp",37274},
{"X-n251-k28.vrp",38684},
{"X-n256-k16.vrp",18839},
{"X-n261-k13.vrp",26558},
{"X-n266-k58.vrp",75478},
{"X-n270-k35.vrp",35291},
{"X-n275-k28.vrp",21245},
{"X-n280-k17.vrp",33503},
{"X-n284-k15.vrp",20215},
{"X-n289-k60.vrp",95151},
{"X-n294-k50.vrp",47161},
{"X-n298-k31.vrp",34231},
{"X-n303-k21.vrp",21736},
{"X-n308-k13.vrp",25859},
{"X-n313-k71.vrp",94043},
{"X-n317-k53.vrp",78355},
{"X-n322-k28.vrp",29834},
{"X-n327-k20.vrp",27532},
{"X-n331-k15.vrp",31102},
{"X-n336-k84.vrp",139111},
{"X-n344-k43.vrp",42050},
{"X-n351-k40.vrp",25896},
{"X-n359-k29.vrp",51505},
{"X-n367-k17.vrp",22814},
{"X-n376-k94.vrp",147713},
{"X-n384-k52.vrp",65928},
{"X-n393-k38.vrp",38260},
{"X-n401-k29.vrp",66154},
{"X-n411-k19.vrp",19712},
{"X-n420-k130.vrp",107798},
{"X-n429-k61.vrp",65449},
{"X-n439-k37.vrp",36391},
{"X-n449-k29.vrp",55233},
{"X-n459-k26.vrp",24139},
{"X-n469-k138.vrp",221824},
{"X-n480-k70.vrp",89449},
{"X-n491-k59.vrp",66483},
{"X-n502-k39.vrp",69226},
{"X-n513-k21.vrp",24201},
{"X-n524-k153.vrp",154593},
{"X-n536-k96.vrp",94846},
{"X-n548-k50.vrp",86700},
{"X-n561-k42.vrp",42717},
{"X-n573-k30.vrp",50673},
{"X-n586-k159.vrp",190316},
{"X-n599-k92.vrp",108451},
{"X-n613-k62.vrp",59353},
{"X-n627-k43.vrp",62164},
{"X-n641-k35.vrp",63682},
{"X-n655-k131.vrp",106780},
{"X-n670-k130.vrp",146332},
{"X-n685-k75.vrp",68205},
{"X-n701-k44.vrp",81923},
{"X-n716-k35.vrp",43373},
{"X-n733-k159.vrp",136187},
{"X-n749-k98.vrp",77269},
{"X-n766-k71.vrp",114417},
{"X-n783-k48.vrp",72386},
{"X-n801-k40.vrp",73305},
{"X-n819-k171.vrp",158121},
{"X-n837-k142.vrp",193737},
{"X-n856-k95.vrp",88965},
{"X-n876-k59.vrp",99299},
{"X-n895-k37.vrp",53860},
{"X-n916-k207.vrp",329179},
{"X-n936-k151.vrp",132715},
{"X-n957-k87.vrp",85465},
{"X-n979-k58.vrp",118976},
{"lou-n88-k2.vrp",41},
{"rbg010a.tsptw",671},
{"rbg016a.tsptw",938},
{"rbg016b.tsptw",1304},
{"rbg017.2.tsptw",852},
{"rbg017.tsptw",893},
{"rbg017a.tsptw",4296},
{"rbg019a.tsptw",1262},
{"rbg019b.tsptw",1866},
{"rbg019c.tsptw",4536},
{"rbg019d.tsptw",1356},
{"rbg020a.tsptw",4689},
{"rbg021.2.tsptw",4528},
{"rbg021.3.tsptw",4528},
{"rbg021.4.tsptw",4525},
{"rbg021.5.tsptw",4515},
{"rbg021.6.tsptw",4480},
{"rbg021.7.tsptw",4479},
{"rbg021.8.tsptw",4478},
{"rbg021.9.tsptw",4478},
{"rbg021.tsptw",4536},
{"rbg027a.tsptw",5091},
{"rbg031a.tsptw",1863},
{"rbg033a.tsptw",2069},
{"rbg034a.tsptw",2222},
{"rbg035a.2.tsptw",2056},
{"rbg035a.tsptw",2144},
{"rbg038a.tsptw",2480},
{"rbg040a.tsptw",2378},
{"rbg041a.tsptw",2598},
{"rbg042a.tsptw",2772},
{"rbg048a.tsptw",9383},
{"rbg049a.tsptw",10018},
{"rbg050a.tsptw",2953},
{"rbg050b.tsptw",9863},
{"rbg050c.tsptw",10024},
{"rbg055a.tsptw",3761},
{"rbg067a.tsptw",4625},
{"rbg086a.tsptw",8400},
{"rbg092a.tsptw",7158},
{"rbg125a.tsptw",7936},
{"rbg132.2.tsptw",8192},
{"rbg132.tsptw",8468},
{"rbg152.3.tsptw",9788},
{"rbg152.tsptw",10032},
{"rbg172a.tsptw",10950},
{"rbg193.2.tsptw",12142},
{"rbg193.tsptw",12535},
{"rbg201a.tsptw",12948},
{"rbg233.2.tsptw",14495},
{"rbg233.tsptw",14992},
{"n100w20.001.tsptw",738},
{"n100w20.002.tsptw",715},
{"n100w20.003.tsptw",762},
{"n100w20.004.tsptw",799},
{"n100w20.005.tsptw",774},
{"n100w40.001.tsptw",770},
{"n100w40.002.tsptw",653},
{"n100w40.003.tsptw",736},
{"n100w40.004.tsptw",651},
{"n100w40.005.tsptw",699},
{"n100w60.001.tsptw",655},
{"n100w60.002.tsptw",659},
{"n100w60.003.tsptw",744},
{"n100w60.004.tsptw",764},
{"n100w60.005.tsptw",661},
{"n150w20.001.tsptw",925},
{"n150w20.002.tsptw",864},
{"n150w20.003.tsptw",834},
{"n150w20.004.tsptw",873},
{"n150w20.005.tsptw",846},
{"n150w40.001.tsptw",918},
{"n150w40.002.tsptw",941},
{"n150w40.003.tsptw",727},
{"n150w40.004.tsptw",764},
{"n150w40.005.tsptw",824},
{"n150w60.001.tsptw",859},
{"n150w60.002.tsptw",782},
{"n150w60.003.tsptw",793},
{"n150w60.004.tsptw",819},
{"n150w60.005.tsptw",840},
{"n200w20.001.tsptw",1019},
{"n200w20.002.tsptw",972},
{"n200w20.003.tsptw",1050},
{"n200w20.004.tsptw",984},
{"n200w20.005.tsptw",1020},
{"n200w40.001.tsptw",1023},
{"n200w40.002.tsptw",948},
{"n200w40.003.tsptw",933},
{"n200w40.004.tsptw",980},
{"n200w40.005.tsptw",1037},
{"n20w100.001.tsptw",237},
{"n20w100.002.tsptw",222},
{"n20w100.003.tsptw",310},
{"n20w100.004.tsptw",349},
{"n20w100.005.tsptw",258},
{"n20w20.001.tsptw",378},
{"n20w20.002.tsptw",286},
{"n20w20.003.tsptw",394},
{"n20w20.004.tsptw",396},
{"n20w20.005.tsptw",352},
{"n20w40.001.tsptw",254},
{"n20w40.002.tsptw",333},
{"n20w40.003.tsptw",317},
{"n20w40.004.tsptw",388},
{"n20w40.005.tsptw",288},
{"n20w60.001.tsptw",335},
{"n20w60.002.tsptw",244},
{"n20w60.003.tsptw",352},
{"n20w60.004.tsptw",280},
{"n20w60.005.tsptw",338},
{"n20w80.001.tsptw",329},
{"n20w80.002.tsptw",338},
{"n20w80.003.tsptw",320},
{"n20w80.004.tsptw",304},
{"n20w80.005.tsptw",264},
{"n40w100.001.tsptw",429},
{"n40w100.002.tsptw",358},
{"n40w100.003.tsptw",364},
{"n40w100.004.tsptw",357},
{"n40w100.005.tsptw",377},
{"n40w20.001.tsptw",500},
{"n40w20.002.tsptw",552},
{"n40w20.003.tsptw",478},
{"n40w20.004.tsptw",404},
{"n40w20.005.tsptw",499},
{"n40w40.001.tsptw",465},
{"n40w40.002.tsptw",461},
{"n40w40.003.tsptw",474},
{"n40w40.004.tsptw",452},
{"n40w40.005.tsptw",453},
{"n40w60.001.tsptw",494},
{"n40w60.002.tsptw",470},
{"n40w60.003.tsptw",408},
{"n40w60.004.tsptw",382},
{"n40w60.005.tsptw",328},
{"n40w80.001.tsptw",395},
{"n40w80.002.tsptw",431},
{"n40w80.003.tsptw",412},
{"n40w80.004.tsptw",417},
{"n40w80.005.tsptw",344},
{"n60w100.001.tsptw",515},
{"n60w100.002.tsptw",538},
{"n60w100.003.tsptw",560},
{"n60w100.004.tsptw",510},
{"n60w100.005.tsptw",451},
{"n60w20.001.tsptw",551},
{"n60w20.002.tsptw",605},
{"n60w20.003.tsptw",533},
{"n60w20.004.tsptw",616},
{"n60w20.005.tsptw",603},
{"n60w40.001.tsptw",591},
{"n60w40.002.tsptw",621},
{"n60w40.003.tsptw",603},
{"n60w40.004.tsptw",597},
{"n60w40.005.tsptw",539},
{"n60w60.001.tsptw",609},
{"n60w60.002.tsptw",566},
{"n60w60.003.tsptw",485},
{"n60w60.004.tsptw",571},
{"n60w60.005.tsptw",569},
{"n60w80.001.tsptw",458},
{"n60w80.002.tsptw",498},
{"n60w80.003.tsptw",550},
{"n60w80.004.tsptw",566},
{"n60w80.005.tsptw",468},
{"n80w20.001.tsptw",616},
{"n80w20.002.tsptw",737},
{"n80w20.003.tsptw",667},
{"n80w20.004.tsptw",615},
{"n80w20.005.tsptw",748},
{"n80w40.001.tsptw",606},
{"n80w40.002.tsptw",618},
{"n80w40.003.tsptw",674},
{"n80w40.004.tsptw",557},
{"n80w40.005.tsptw",695},
{"n80w60.001.tsptw",554},
{"n80w60.002.tsptw",633},
{"n80w60.003.tsptw",651},
{"n80w60.004.tsptw",619},
{"n80w60.005.tsptw",575},
{"n80w80.001.tsptw",624},
{"n80w80.002.tsptw",592},
{"n80w80.003.tsptw",589},
{"n80w80.004.tsptw",594},
{"n80w80.005.tsptw",570},
{"n100w100.001.tsptw",643},
{"n100w100.002.tsptw",619},
{"n100w100.003.tsptw",685},
{"n100w100.004.tsptw",684},
{"n100w100.005.tsptw",572},
{"n100w120.001.tsptw",629},
{"n100w120.002.tsptw",540},
{"n100w120.003.tsptw",617},
{"n100w120.004.tsptw",663},
{"n100w120.005.tsptw",537},
{"n100w140.001.tsptw",604},
{"n100w140.002.tsptw",615},
{"n100w140.003.tsptw",481},
{"n100w140.004.tsptw",533},
{"n100w140.005.tsptw",509},
{"n100w160.001.tsptw",582},
{"n100w160.002.tsptw",532},
{"n100w160.003.tsptw",495},
{"n100w160.004.tsptw",580},
{"n100w160.005.tsptw",586},
{"n100w80.001.tsptw",670},
{"n100w80.002.tsptw",668},
{"n100w80.003.tsptw",691},
{"n100w80.004.tsptw",700},
{"n100w80.005.tsptw",603},
{"n20w120.001.tsptw",267},
{"n20w120.002.tsptw",218},
{"n20w120.003.tsptw",303},
{"n20w120.004.tsptw",300},
{"n20w120.005.tsptw",240},
{"n20w140.001.tsptw",176},
{"n20w140.002.tsptw",272},
{"n20w140.003.tsptw",236},
{"n20w140.004.tsptw",255},
{"n20w140.005.tsptw",225},
{"n20w160.001.tsptw",241},
{"n20w160.002.tsptw",201},
{"n20w160.003.tsptw",201},
{"n20w160.004.tsptw",203},
{"n20w160.005.tsptw",245},
{"n20w180.001.tsptw",253},
{"n20w180.002.tsptw",265},
{"n20w180.003.tsptw",271},
{"n20w180.004.tsptw",201},
{"n20w180.005.tsptw",193},
{"n20w200.001.tsptw",233},
{"n20w200.002.tsptw",203},
{"n20w200.003.tsptw",249},
{"n20w200.004.tsptw",293},
{"n20w200.005.tsptw",227},
{"n40w120.001.tsptw",434},
{"n40w120.002.tsptw",445},
{"n40w120.003.tsptw",357},
{"n40w120.004.tsptw",303},
{"n40w120.005.tsptw",350},
{"n40w140.001.tsptw",328},
{"n40w140.002.tsptw",383},
{"n40w140.003.tsptw",398},
{"n40w140.004.tsptw",342},
{"n40w140.005.tsptw",371},
{"n40w160.001.tsptw",348},
{"n40w160.002.tsptw",337},
{"n40w160.003.tsptw",346},
{"n40w160.004.tsptw",288},
{"n40w160.005.tsptw",315},
{"n40w180.001.tsptw",337},
{"n40w180.002.tsptw",347},
{"n40w180.003.tsptw",279},
{"n40w180.004.tsptw",354},
{"n40w180.005.tsptw",335},
{"n40w200.001.tsptw",330},
{"n40w200.002.tsptw",303},
{"n40w200.003.tsptw",339},
{"n40w200.004.tsptw",301},
{"n40w200.005.tsptw",296},
{"n60w120.001.tsptw",384},
{"n60w120.002.tsptw",427},
{"n60w120.003.tsptw",407},
{"n60w120.004.tsptw",490},
{"n60w120.005.tsptw",547},
{"n60w140.001.tsptw",423},
{"n60w140.002.tsptw",462},
{"n60w140.003.tsptw",427},
{"n60w140.004.tsptw",488},
{"n60w140.005.tsptw",460},
{"n60w160.001.tsptw",560},
{"n60w160.002.tsptw",423},
{"n60w160.003.tsptw",434},
{"n60w160.004.tsptw",401},
{"n60w160.005.tsptw",502},
{"n60w180.001.tsptw",411},
{"n60w180.002.tsptw",399},
{"n60w180.003.tsptw",445},
{"n60w180.004.tsptw",456},
{"n60w180.005.tsptw",395},
{"n60w200.001.tsptw",410},
{"n60w200.002.tsptw",414},
{"n60w200.003.tsptw",455},
{"n60w200.004.tsptw",431},
{"n60w200.005.tsptw",427},
{"n80w100.001.tsptw",565},
{"n80w100.002.tsptw",567},
{"n80w100.003.tsptw",580},
{"n80w100.004.tsptw",649},
{"n80w100.005.tsptw",532},
{"n80w120.001.tsptw",498},
{"n80w120.002.tsptw",577},
{"n80w120.003.tsptw",540},
{"n80w120.004.tsptw",501},
{"n80w120.005.tsptw",591},
{"n80w140.001.tsptw",512},
{"n80w140.002.tsptw",470},
{"n80w140.003.tsptw",580},
{"n80w140.004.tsptw",423},
{"n80w140.005.tsptw",545},
{"n80w160.001.tsptw",506},
{"n80w160.002.tsptw",549},
{"n80w160.003.tsptw",521},
{"n80w160.004.tsptw",509},
{"n80w160.005.tsptw",439},
{"n80w180.001.tsptw",551},
{"n80w180.002.tsptw",479},
{"n80w180.003.tsptw",524},
{"n80w180.004.tsptw",479},
{"n80w180.005.tsptw",470},
{"n80w200.001.tsptw",490},
{"n80w200.002.tsptw",488},
{"n80w200.003.tsptw",464},
{"n80w200.004.tsptw",526},
{"n80w200.005.tsptw",439},
{"N20ft301.tsptw",661.6},
{"N20ft302.tsptw",684.2},
{"N20ft303.tsptw",746.4},
{"N20ft304.tsptw",817},
{"N20ft305.tsptw",716.5},
{"N20ft306.tsptw",727.8},
{"N20ft307.tsptw",691.8},
{"N20ft308.tsptw",788.2},
{"N20ft309.tsptw",730.7},
{"N20ft310.tsptw",683},
{"N20ft401.tsptw",660.8},
{"N20ft402.tsptw",684.2},
{"N20ft403.tsptw",746.4},
{"N20ft404.tsptw",817},
{"N20ft405.tsptw",716.5},
{"N20ft406.tsptw",727.8},
{"N20ft407.tsptw",691.8},
{"N20ft408.tsptw",757.3},
{"N20ft409.tsptw",730.7},
{"N20ft410.tsptw",683},
{"N40ft201.tsptw",1100.6},
{"N40ft202.tsptw",1010.4},
{"N40ft203.tsptw",876.8},
{"N40ft204.tsptw",885.8},
{"N40ft205.tsptw",940.9},
{"N40ft206.tsptw",1054.2},
{"N40ft207.tsptw",867.5},
{"N40ft208.tsptw",1050.7},
{"N40ft209.tsptw",1013.9},
{"N40ft210.tsptw",1026.3},
{"N40ft401.tsptw",1085},
{"N40ft402.tsptw",995.6},
{"N40ft403.tsptw",845.8},
{"N40ft404.tsptw",868},
{"N40ft405.tsptw",936.5},
{"N40ft406.tsptw",969.1},
{"N40ft407.tsptw",831.2},
{"N40ft408.tsptw",1002.7},
{"N40ft409.tsptw",1000.5},
{"N40ft410.tsptw",983.8},
{"N60ft201.tsptw",1353.5},
{"N60ft202.tsptw",1161.6},
{"N60ft203.tsptw",1182.9},
{"N60ft204.tsptw",1257.5},
{"N60ft205.tsptw",1184.1},
{"N60ft206.tsptw",1199.6},
{"N60ft207.tsptw",1299},
{"N60ft208.tsptw",1113},
{"N60ft209.tsptw",1171.3},
{"N60ft210.tsptw",1234.3},
{"N60ft301.tsptw",1337},
{"N60ft302.tsptw",1089.5},
{"N60ft303.tsptw",1179},
{"N60ft304.tsptw",1230},
{"N60ft305.tsptw",1151.6},
{"N60ft306.tsptw",1167.9},
{"N60ft307.tsptw",1220.1},
{"N60ft308.tsptw",1097.6},
{"N60ft309.tsptw",1140.6},
{"N60ft310.tsptw",1219.2},
{"N60ft401.tsptw",1335},
{"N60ft402.tsptw",1088.1},
{"N60ft403.tsptw",1173.7},
{"N60ft404.tsptw",1184.7},
{"N60ft405.tsptw",1146.2},
{"N60ft406.tsptw",1140.2},
{"N60ft407.tsptw",1198.9},
{"N60ft408.tsptw",1029.4},
{"N60ft409.tsptw",1121.4},
{"N60ft410.tsptw",1189.6},
{"n150w120.001.tsptw",734},
{"n150w120.002.tsptw",677},
{"n150w120.003.tsptw",747},
{"n150w120.004.tsptw",763},
{"n150w120.005.tsptw",689},
{"n150w140.001.tsptw",762},
{"n150w140.002.tsptw",755},
{"n150w140.003.tsptw",613},
{"n150w140.004.tsptw",676},
{"n150w140.005.tsptw",663},
{"n150w160.001.tsptw",706},
{"n150w160.002.tsptw",711},
{"n150w160.003.tsptw",608},
{"n150w160.004.tsptw",672},
{"n150w160.005.tsptw",658},
{"n200w120.001.tsptw",799},
{"n200w120.002.tsptw",721},
{"n200w120.003.tsptw",880},
{"n200w120.004.tsptw",777},
{"n200w120.005.tsptw",841},
{"n200w140.001.tsptw",834},
{"n200w140.002.tsptw",760},
{"n200w140.003.tsptw",758},
{"n200w140.004.tsptw",816},
{"n200w140.005.tsptw",822},
{"rc201.0.tsptw",628.62},
{"rc201.1.tsptw",654.7},
{"rc201.2.tsptw",707.65},
{"rc201.3.tsptw",422.54},
{"rc202.0.tsptw",496.22},
{"rc202.1.tsptw",426.53},
{"rc202.2.tsptw",611.77},
{"rc202.3.tsptw",627.85},
{"rc203.0.tsptw",727.45},
{"rc203.1.tsptw",726.99},
{"rc203.2.tsptw",617.46},
{"rc204.0.tsptw",541.45},
{"rc204.1.tsptw",485.37},
{"rc204.2.tsptw",778.4},
{"rc205.0.tsptw",511.65},
{"rc205.1.tsptw",491.22},
{"rc205.2.tsptw",714.69},
{"rc205.3.tsptw",601.24},
{"rc206.0.tsptw",835.23},
{"rc206.1.tsptw",664.73},
{"rc206.2.tsptw",655.37},
{"rc207.0.tsptw",806.69},
{"rc207.1.tsptw",726.36},
{"rc207.2.tsptw",546.41},
{"rc208.0.tsptw",820.56},
{"rc208.1.tsptw",509.04},
{"rc208.2.tsptw",503.92},
{"rc_201.1.tsptw",444.54},
{"rc_201.2.tsptw",711.54},
{"rc_201.3.tsptw",790.61},
{"rc_201.4.tsptw",793.64},
{"rc_202.1.tsptw",771.78},
{"rc_202.2.tsptw",304.14},
{"rc_202.3.tsptw",837.72},
{"rc_202.4.tsptw",793.03},
{"rc_203.1.tsptw",453.48},
{"rc_203.2.tsptw",784.16},
{"rc_203.3.tsptw",817.53},
{"rc_203.4.tsptw",314.29},
{"rc_204.1.tsptw",878.64},
{"rc_204.2.tsptw",662.16},
{"rc_204.3.tsptw",455.03},
{"rc_205.1.tsptw",343.21},
{"rc_205.2.tsptw",755.93},
{"rc_205.3.tsptw",825.06},
{"rc_205.4.tsptw",760.47},
{"rc_206.1.tsptw",117.85},
{"rc_206.2.tsptw",828.06},
{"rc_206.3.tsptw",574.42},
{"rc_206.4.tsptw",831.67},
{"rc_207.1.tsptw",732.68},
{"rc_207.2.tsptw",701.25},
{"rc_207.3.tsptw",682.4},
{"rc_207.4.tsptw",119.64},
{"rc_208.1.tsptw",789.25},
{"rc_208.2.tsptw",533.78},
{"rc_208.3.tsptw",634.44},
{"ESC07.sop",2125},
{"ESC12.sop",1675},
{"ESC25.sop",1681},
{"ESC47.sop",1288},
{"ESC63.sop",62},
{"ESC78.sop",18230},
{"br17.10.sop",55},
{"br17.12.sop",55},
{"ft53.1.sop",7531},
{"ft53.2.sop",8026},
{"ft53.3.sop",10262},
{"ft53.4.sop",14425},
{"ft70.1.sop",39313},
{"ft70.2.sop",40419},
{"ft70.3.sop",42535},
{"ft70.4.sop",53530},
{"kro124p.1.sop",39420},
{"kro124p.2.sop",41336},
{"kro124p.3.sop",49499},
{"kro124p.4.sop",76103},
{"p43.1.sop",28140},
{"p43.2.sop",28480},
{"p43.3.sop",28835},
{"p43.4.sop",83005},
{"prob.42.sop",243},
{"prob.100.sop",1163},
{"rbg048a.sop",351},
{"rbg050c.sop",467},
{"rbg109a.sop",1038},
{"rbg150a.sop",1750},
{"rbg174a.sop",2033},
{"rbg253a.sop",2950},
{"rbg323a.sop",3140},
{"rbg341a.sop",2568},
{"rbg358a.sop",2545},
{"rbg378a.sop",2816},
{"ry48p.1.sop",15805},
{"ry48p.2.sop",16666},
{"ry48p.3.sop",19894},
{"ry48p.4.sop",31446}
};

struct VRPTW
{
  public:
    VRPTW(std::string fileName)
    {
      // VRPTW instance format
      std::regex vrptwInstanceRegex("[^_]vrptw");
      std::smatch vrptwInstanceMatch;
      if (std::regex_search(fileName, vrptwInstanceMatch, vrptwInstanceRegex))
      {
        oneOrMorePaths = OneOrMorePaths::MORE_PATHS;
        circuitOrPath = CircuitOrPath::CIRCUIT;
        timeStateMultiplier = 10;
        timeStateDiscretization = 10;
        std::vector<std::pair<double,double> > coordinates;

        std::regex noCapacityRelaxRegex(".*HG.*C1|R1|RC1.*txt");
        std::smatch noCapacityRelaxMatch;
        vrptwCapacityType = VRPTWCapacityType::RELAX_CAPACITY;
        if (std::regex_search(fileName, noCapacityRelaxMatch, noCapacityRelaxRegex))
        {
          std::cout << "do not relax capacity constraint" << std::endl;
          vrptwCapacityType = VRPTWCapacityType::NO_RELAX_CAPACITY;
        }
        else
        {
          std::cout << "relax capacity constraint" << std::endl;
        }
        vrptwTimeWindowType = VRPTWTimeWindowType::TIME_WINDOWS;

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

        int minServiceTime = *std::min_element(serviceTimes.begin()+1, serviceTimes.end());
        timeStateDiscretization = std::max(timeStateDiscretization, minServiceTime);

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
            distance = (int)( 10 * distance ) / 10.0;
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
      }

      // CVRP instance format, relax num trucks constraint
      std::regex cvrpInstanceRegex(".*-k.*");
      std::smatch cvrpInstanceMatch;
      if (std::regex_search(fileName, cvrpInstanceMatch, cvrpInstanceRegex))
      {
        timeStateMultiplier = 1;
        timeStateDiscretization = 10;
        oneOrMorePaths = OneOrMorePaths::MORE_PATHS;
        circuitOrPath = CircuitOrPath::CIRCUIT;
        vrptwCapacityType = VRPTWCapacityType::NO_RELAX_CAPACITY;
        vrptwTimeWindowType = VRPTWTimeWindowType::NO_TIME_WINDOWS;
        std::cout << "relax time window constraint" << std::endl;
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
              std::cout << "obtain capacity: " << std::endl;
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
      }
 
      // TSPTW instance format
      std::regex tsptwInstanceRegex(".*tsptw");
      std::smatch tsptwInstanceMatch;
      if (std::regex_search(fileName, tsptwInstanceMatch, tsptwInstanceRegex))
      {
        // afg, dumas, gendreau-dumas, ohlmann thomas, integer
        timeStateMultiplier = 1;
        timeStateDiscretization = 1;

        // solomon potvin bengio, solomon pesant, 0.0000
        std::regex solomonRegex(".*Solomon.*");
        std::smatch solomonMatch;
        if (std::regex_search(fileName, solomonMatch, solomonRegex))
        {
          timeStateMultiplier = 10000;
          timeStateDiscretization = 100;
        }

        // langevin 0.1
        std::regex langevinRegex(".*Langevin.*");
        std::smatch langevinMatch;
        if (std::regex_search(fileName, langevinMatch, langevinRegex))
        {
          timeStateMultiplier = 10;
          timeStateDiscretization = 10;
        }

        oneOrMorePaths = OneOrMorePaths::ONE_PATH;
        circuitOrPath = CircuitOrPath::CIRCUIT;
        std::vector<std::pair<double,double> > coordinates;
        vrptwCapacityType = VRPTWCapacityType::RELAX_CAPACITY;
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
              std::cout << "num loc: " << numLocations << std::endl;
            }

            numLocationsSection = false;
            distanceMatrixSection = true;
            distances.resize(numLocations);
            for (int i=0; i < distances.size(); ++i)
            {
              distances[i].resize(distances.size());
            }
          }
        }

        // formally add if this works
        std::string instanceName = fileName.substr(fileName.find_last_of("/") + 1);
        if (instanceOptimalSolutions.find(instanceName) != instanceOptimalSolutions.end())
        {
          hgsUpperBound = instanceOptimalSolutions.find(instanceName)->second;

          // rounding in TSPTW
          hgsUpperBound += 0.01;
        }
        else
        {
          hgsUpperBound = INF;
        }
        std::cout << "hgs upper bound: " << hgsUpperBound << std::endl;
      }

      // SOP instance format
      std::regex sopInstanceRegex(".*sop");
      std::smatch sopInstanceMatch;
      if (std::regex_search(fileName, sopInstanceMatch, sopInstanceRegex))
      {
        timeStateMultiplier = 1;
        timeStateDiscretization = 1;

        oneOrMorePaths = OneOrMorePaths::ONE_PATH;
        circuitOrPath = CircuitOrPath::PATH;
        std::vector<std::pair<double,double> > coordinates;
        vrptwCapacityType = VRPTWCapacityType::RELAX_CAPACITY;
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
              }

              distances[distanceMatrixRow][column] = distance;
              column = column + 1;

              demands.push_back(1);
              demandsForSeparation.push_back(1);
              demandsForCombs.push_back(1);
              serviceTimes.push_back(0);
              startTimes.push_back(0);
              endTimes.push_back(0);
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
          }
        }

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
      }

      // PDP instance format
      std::regex pdpInstanceRegex(".*pdp");
      std::smatch pdpInstanceMatch;
      if (std::regex_search(fileName, pdpInstanceMatch, pdpInstanceRegex))
      {
        timeStateMultiplier = 100;
        timeStateDiscretization = 10;

        oneOrMorePaths = OneOrMorePaths::MORE_PATHS;
        circuitOrPath = CircuitOrPath::CIRCUIT;
        std::vector<std::pair<double,double> > coordinates;
        vrptwCapacityType = VRPTWCapacityType::NO_RELAX_CAPACITY;
        vrptwTimeWindowType = VRPTWTimeWindowType::TIME_WINDOWS;

        bool dataSection = true;
        bool locationDataSection = false;
        depot = 0;
        precedences.resize(2000);
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
            int numVehicles, capacityInput, speed;
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
              if (pickup == 0)
              {
                precedences[delivery].insert(location);
              }
              if (delivery == 0)
              {
                precedences[location].insert(pickup);
              }
            }
          }
        }

        numLocations = demands.size();
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
            distance = (int)( 100 * distance ) / 100.00;
            distances[i][j] = distance;
            distances[j][i] = distance;
          }
        }

        // formally add if this works
        std::string instanceName = fileName.substr(fileName.find_last_of("/") + 1);
        if (instanceOptimalSolutions.find(instanceName) != instanceOptimalSolutions.end())
        {
          hgsUpperBound = instanceOptimalSolutions.find(instanceName)->second;

          // rounding in TSPTW
          hgsUpperBound += 0.01;
        }
        else
        {
          hgsUpperBound = INF;
        }
        std::cout << "hgs upper bound: " << hgsUpperBound << std::endl;
      }

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
    std::vector<std::set<int> > precedences;
    std::vector<int> demands;
    std::vector<double> demandsForSeparation;
    std::vector<int> demandsForCombs;

    std::vector<int> startTimes;
    std::vector<int> endTimes;
    std::vector<int> serviceTimes;

    int depot;
    int numLocations;

    VRPTWCapacityType vrptwCapacityType;
    VRPTWTimeWindowType vrptwTimeWindowType;
    OneOrMorePaths oneOrMorePaths;
    CircuitOrPath circuitOrPath;
    int timeStateMultiplier;
    int timeStateDiscretization;
    double hgsUpperBound;
};

#endif
