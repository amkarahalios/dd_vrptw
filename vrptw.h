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

enum OneOrMorePaths
{
  ONE_PATH = 0,
  MORE_PATHS = 1
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

enum VRPTWType
{
  RELAX_CAPACITY = 0,
  NO_RELAX_CAPACITY = 1,
  RELAX_TIME_WINDOWS = 2
};

struct VRPTWSolution
{
  VRPTWSolution(std::vector<std::vector<int>> _routes, double _totalDistance) : routes(_routes), totalDistance(_totalDistance) {}

  std::vector<std::vector<int>> routes;
  double totalDistance;
};

const static std::map<std::string,double> instanceOptimalSolutions =
{
{"C101.txt",827.3},
{"C102.txt",827.3},
{"C103.txt",826.3},
{"C104.txt",822.9},
{"C105.txt",827.3},
{"C106.txt",827.3},
{"C107.txt",827.3},
{"C108.txt",827.3},
{"C109.txt",827.3},
{"C201.txt",589.1},
{"C202.txt",589.1},
{"C203.txt",588.7},
{"C204.txt",588.1},
{"C205.txt",586.4},
{"C206.txt",586},
{"C207.txt",585.8},
{"C208.txt",585.8},
{"R101.txt",1637.7},
{"R102.txt",1466.6},
{"R103.txt",1208.7},
{"R104.txt",971.5},
{"R105.txt",1355.3},
{"R106.txt",1234.6},
{"R107.txt",1064.6},
{"R108.txt",932.1},
{"R109.txt",1146.9},
{"R110.txt",1068},
{"R111.txt",1048.7},
{"R112.txt",948.6},
{"R201.txt",1143.2},
{"R202.txt",1029.6},
{"R203.txt",870.8},
{"R204.txt",731.3},
{"R205.txt",949.8},
{"R206.txt",875.9},
{"R207.txt",794},
{"R208.txt",701},
{"R209.txt",854.8},
{"R210.txt",900.5},
{"R211.txt",746.7},
{"RC101.txt",1619.8},
{"RC102.txt",1457.4},
{"RC103.txt",1258},
{"RC104.txt",1132.3},
{"RC105.txt",1513.7},
{"RC106.txt",1372.7},
{"RC107.txt",1207.8},
{"RC108.txt",1114.2},
{"RC201.txt",1261.8},
{"RC202.txt",1092.3},
{"RC203.txt",923.7},
{"RC204.txt",783.5},
{"RC205.txt",1154},
{"RC206.txt",1051.1},
{"RC207.txt",962.9},
{"RC208.txt",776.1},
{"C1_2_1.txt",2698.6},
{"C1_2_2.txt",2694.3},
{"C1_2_3.txt",2675.8},
{"C1_2_4.txt",2625.6},
{"C1_2_5.txt",2694.9},
{"C1_2_6.txt",2694.9},
{"C1_2_7.txt",2694.9},
{"C1_2_8.txt",2684},
{"C1_2_9.txt",2639.6},
{"C1_2_10.txt",2624.7},
{"C2_2_1.txt",1922.1},
{"C2_2_2.txt",1851.4},
{"C2_2_3.txt",1763.4},
{"C2_2_4.txt",1695},
{"C2_2_5.txt",1869.6},
{"C2_2_6.txt",1844.8},
{"C2_2_7.txt",1842.2},
{"C2_2_8.txt",1813.7},
{"C2_2_9.txt",1815},
{"C2_2_10.txt",1791.2},
{"R1_2_1.txt",4667.2},
{"R1_2_2.txt",3919.9},
{"R1_2_3.txt",3373.9},
{"R1_2_4.txt",3047.6},
{"R1_2_5.txt",4053.2},
{"R1_2_6.txt",3559.1},
{"R1_2_7.txt",3141.9},
{"R1_2_8.txt",2938.4},
{"R1_2_9.txt",3734.7},
{"R1_2_10.txt",3293.1},
{"R2_2_1.txt",3468},
{"R2_2_2.txt",3008.2},
{"R2_2_3.txt",2537.5},
{"R2_2_4.txt",1928.5},
{"R2_2_5.txt",3061.1},
{"R2_2_6.txt",2675.4},
{"R2_2_7.txt",2304.7},
{"R2_2_8.txt",1842.4},
{"R2_2_9.txt",2843.3},
{"R2_2_10.txt",2549.4},
{"RC1_2_1.txt",3516.9},
{"RC1_2_2.txt",3221.6},
{"RC1_2_3.txt",3001.4},
{"RC1_2_4.txt",2845.2},
{"RC1_2_5.txt",3325.6},
{"RC1_2_6.txt",3300.7},
{"RC1_2_7.txt",3177.8},
{"RC1_2_8.txt",3060},
{"RC1_2_9.txt",3073.3},
{"RC1_2_10.txt",2990.5},
{"RC2_2_1.txt",2797.4},
{"RC2_2_2.txt",2481.6},
{"RC2_2_3.txt",2227.7},
{"RC2_2_4.txt",1854.8},
{"RC2_2_5.txt",2491.4},
{"RC2_2_6.txt",2495.1},
{"RC2_2_7.txt",2287.7},
{"RC2_2_8.txt",2151.2},
{"RC2_2_9.txt",2086.6},
{"RC2_2_10.txt",1989.2},
{"C1_4_1.txt",7138.8},
{"C1_4_2.txt",7113.3},
{"C1_4_3.txt",6929.9},
{"C1_4_4.txt",6777.7},
{"C1_4_5.txt",7138.8},
{"C1_4_6.txt",7140.1},
{"C1_4_7.txt",7136.2},
{"C1_4_8.txt",7083},
{"C1_4_9.txt",6927.8},
{"C1_4_10.txt",6825.4},
{"C2_4_1.txt",4100.3},
{"C2_4_2.txt",3914.1},
{"C2_4_3.txt",3755.2},
{"C2_4_4.txt",3523.7},
{"C2_4_5.txt",3923.2},
{"C2_4_6.txt",3860.1},
{"C2_4_7.txt",3870.9},
{"C2_4_8.txt",3773.7},
{"C2_4_9.txt",3842.1},
{"C2_4_10.txt",3665.1},
{"R1_4_1.txt",10305.8},
{"R1_4_2.txt",8873.3},
{"R1_4_3.txt",7784.3},
{"R1_4_4.txt",7266.2},
{"R1_4_5.txt",9184.6},
{"R1_4_6.txt",8340.4},
{"R1_4_7.txt",7599.8},
{"R1_4_8.txt",7240.5},
{"R1_4_9.txt",8677.5},
{"R1_4_10.txt",8077.8},
{"R2_4_1.txt",7520.7},
{"R2_4_2.txt",6482.8},
{"R2_4_3.txt",5372.9},
{"R2_4_4.txt",4211.2},
{"R2_4_5.txt",6567.9},
{"R2_4_6.txt",5813.5},
{"R2_4_7.txt",4893.5},
{"R2_4_8.txt",4000.1},
{"R2_4_9.txt",6067.8},
{"R2_4_10.txt",5645.9},
{"RC1_4_1.txt",8522.9},
{"RC1_4_2.txt",7878.2},
{"RC1_4_3.txt",7516.9},
{"RC1_4_4.txt",7292.9},
{"RC1_4_5.txt",8152.3},
{"RC1_4_6.txt",8148},
{"RC1_4_7.txt",7932.5},
{"RC1_4_8.txt",7757.2},
{"RC1_4_9.txt",7717.7},
{"RC1_4_10.txt",7581.2},
{"RC2_4_1.txt",6147.3},
{"RC2_4_2.txt",5407.5},
{"RC2_4_3.txt",4573},
{"RC2_4_4.txt",3597.9},
{"RC2_4_5.txt",5392.3},
{"RC2_4_6.txt",5324.6},
{"RC2_4_7.txt",4987.8},
{"RC2_4_8.txt",4693.3},
{"RC2_4_9.txt",4510.4},
{"RC2_4_10.txt",4252.3},
{"C1_6_1.txt",14076.6},
{"C1_6_2.txt",13948.3},
{"C1_6_3.txt",13757},
{"C1_6_4.txt",13538.6},
{"C1_6_5.txt",14066.8},
{"C1_6_6.txt",14070.9},
{"C1_6_7.txt",14066.8},
{"C1_6_8.txt",13991.2},
{"C1_6_9.txt",13664.5},
{"C1_6_10.txt",13617.5},
{"C2_6_1.txt",7752.2},
{"C2_6_2.txt",7471.5},
{"C2_6_3.txt",7215},
{"C2_6_4.txt",6877},
{"C2_6_5.txt",7553.8},
{"C2_6_6.txt",7449.8},
{"C2_6_7.txt",7491.3},
{"C2_6_8.txt",7303.7},
{"C2_6_9.txt",7303.2},
{"C2_6_10.txt",7123.9},
{"R1_6_1.txt",21274.2},
{"R1_6_2.txt",18558.7},
{"R1_6_3.txt",16874.9},
{"R1_6_4.txt",15721.4},
{"R1_6_5.txt",19294.9},
{"R1_6_6.txt",17763.7},
{"R1_6_7.txt",16496.2},
{"R1_6_8.txt",15584.3},
{"R1_6_9.txt",18474.1},
{"R1_6_10.txt",17583.7},
{"R2_6_1.txt",15145.3},
{"R2_6_2.txt",12976.3},
{"R2_6_3.txt",10455.3},
{"R2_6_4.txt",7915.1},
{"R2_6_5.txt",13790.2},
{"R2_6_6.txt",11847.8},
{"R2_6_7.txt",9777.9},
{"R2_6_8.txt",7512.3},
{"R2_6_9.txt",12736.8},
{"R2_6_10.txt",11837},
{"RC1_6_1.txt",16960.1},
{"RC1_6_2.txt",15890.6},
{"RC1_6_3.txt",15181.3},
{"RC1_6_4.txt",14753.2},
{"RC1_6_5.txt",16536.3},
{"RC1_6_6.txt",16473.3},
{"RC1_6_7.txt",16055.3},
{"RC1_6_8.txt",15891.8},
{"RC1_6_9.txt",15803.5},
{"RC1_6_10.txt",15651.3},
{"RC2_6_1.txt",11966.1},
{"RC2_6_2.txt",10336.9},
{"RC2_6_3.txt",8894.9},
{"RC2_6_4.txt",6967.5},
{"RC2_6_5.txt",11080.7},
{"RC2_6_6.txt",10830.5},
{"RC2_6_7.txt",10289.4},
{"RC2_6_8.txt",9779},
{"RC2_6_9.txt",9436},
{"RC2_6_10.txt",8974.7},
{"C1_8_1.txt",25156.9},
{"C1_8_2.txt",24974.1},
{"C1_8_3.txt",24156.1},
{"C1_8_4.txt",23797.3},
{"C1_8_5.txt",25138.6},
{"C1_8_6.txt",25133.3},
{"C1_8_7.txt",25127.3},
{"C1_8_8.txt",24809.7},
{"C1_8_9.txt",24200.4},
{"C1_8_10.txt",24026.7},
{"C2_8_1.txt",11631.9},
{"C2_8_2.txt",11394.5},
{"C2_8_3.txt",11138.1},
{"C2_8_4.txt",10650},
{"C2_8_5.txt",11395.6},
{"C2_8_6.txt",11316.3},
{"C2_8_7.txt",11332.9},
{"C2_8_8.txt",11133.9},
{"C2_8_9.txt",11140.4},
{"C2_8_10.txt",10946},
{"R1_8_1.txt",36345},
{"R1_8_2.txt",32277.6},
{"R1_8_3.txt",29304.5},
{"R1_8_4.txt",27734.7},
{"R1_8_5.txt",33494.2},
{"R1_8_6.txt",30872.4},
{"R1_8_7.txt",28789},
{"R1_8_8.txt",27609.4},
{"R1_8_9.txt",32257.3},
{"R1_8_10.txt",30918.4},
{"R2_8_1.txt",24969.8},
{"R2_8_2.txt",21312.2},
{"R2_8_3.txt",17234.8},
{"R2_8_4.txt",13160.8},
{"R2_8_5.txt",22801.6},
{"R2_8_6.txt",19740.5},
{"R2_8_7.txt",16357.5},
{"R2_8_8.txt",12611.7},
{"R2_8_9.txt",21282.7},
{"R2_8_10.txt",19984.8},
{"RC1_8_1.txt",29978.9},
{"RC1_8_2.txt",28290.1},
{"RC1_8_3.txt",27447.7},
{"RC1_8_4.txt",26557.2},
{"RC1_8_5.txt",29219.9},
{"RC1_8_6.txt",29194.2},
{"RC1_8_7.txt",28788.6},
{"RC1_8_8.txt",28418.1},
{"RC1_8_9.txt",28347.1},
{"RC1_8_10.txt",28168.5},
{"RC2_8_1.txt",19201.3},
{"RC2_8_2.txt",16709.5},
{"RC2_8_3.txt",14013.6},
{"RC2_8_4.txt",10969.4},
{"RC2_8_5.txt",17466.1},
{"RC2_8_6.txt",17195.1},
{"RC2_8_7.txt",16362.2},
{"RC2_8_8.txt",15528.8},
{"RC2_8_9.txt",15183},
{"RC2_8_10.txt",14370.9},
{"C1_10_1.txt",42444.8},
{"C1_10_2.txt",41337.8},
{"C1_10_3.txt",40064.4},
{"C1_10_4.txt",39434.1},
{"C1_10_5.txt",42434.8},
{"C1_10_6.txt",42437},
{"C1_10_7.txt",42420.4},
{"C1_10_8.txt",41652.1},
{"C1_10_9.txt",40288.4},
{"C1_10_10.txt",39816.8},
{"C2_10_1.txt",16841.1},
{"C2_10_2.txt",16462.6},
{"C2_10_3.txt",16036.5},
{"C2_10_4.txt",15459.5},
{"C2_10_5.txt",16521.3},
{"C2_10_6.txt",16290.7},
{"C2_10_7.txt",16378.4},
{"C2_10_8.txt",16029.1},
{"C2_10_9.txt",16075.4},
{"C2_10_10.txt",15728.6},
{"R1_10_1.txt",53046.5},
{"R1_10_2.txt",48263.1},
{"R1_10_3.txt",44677.1},
{"R1_10_4.txt",42440.7},
{"R1_10_5.txt",50406.7},
{"R1_10_6.txt",46930.3},
{"R1_10_7.txt",43997.4},
{"R1_10_8.txt",42279.3},
{"R1_10_9.txt",49162.8},
{"R1_10_10.txt",47364.6},
{"R2_10_1.txt",36881},
{"R2_10_2.txt",31241.9},
{"R2_10_3.txt",24399},
{"R2_10_4.txt",17811.5},
{"R2_10_5.txt",34132.8},
{"R2_10_6.txt",29124.7},
{"R2_10_7.txt",23102.2},
{"R2_10_8.txt",17403.8},
{"R2_10_9.txt",31990.6},
{"R2_10_10.txt",29840.5},
{"RC1_10_1.txt",45790.8},
{"RC1_10_2.txt",43678.3},
{"RC1_10_3.txt",42122},
{"RC1_10_4.txt",41357.4},
{"RC1_10_5.txt",45028.1},
{"RC1_10_6.txt",44903.6},
{"RC1_10_7.txt",44417.1},
{"RC1_10_8.txt",43916.5},
{"RC1_10_9.txt",43858.1},
{"RC1_10_10.txt",43533.7},
{"RC2_10_1.txt",28122.6},
{"RC2_10_2.txt",24248.6},
{"RC2_10_3.txt",19618.1},
{"RC2_10_4.txt",15657},
{"RC2_10_5.txt",25797.5},
{"RC2_10_6.txt",25782.5},
{"RC2_10_7.txt",24395.8},
{"RC2_10_8.txt",23280.2},
{"RC2_10_9.txt",22731.6},
{"RC2_10_10.txt",21736.1},
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
{"rc201.0.tw",628.62},
{"rc201.1.tw",654.7},
{"rc201.2.tw",707.65},
{"rc201.3.tw",422.54},
{"rc202.0.tw",496.22},
{"rc202.1.tw",426.53},
{"rc202.2.tw",611.77},
{"rc202.3.tw",627.85},
{"rc203.0.tw",727.45},
{"rc203.1.tw",726.99},
{"rc203.2.tw",617.46},
{"rc204.0.tw",541.45},
{"rc204.1.tw",485.37},
{"rc204.2.tw",778.4},
{"rc205.0.tw",511.65},
{"rc205.1.tw",491.22},
{"rc205.2.tw",714.69},
{"rc205.3.tw",601.24},
{"rc206.0.tw",835.23},
{"rc206.1.tw",664.73},
{"rc206.2.tw",655.37},
{"rc207.0.tw",806.69},
{"rc207.1.tw",726.36},
{"rc207.2.tw",546.41},
{"rc208.0.tw",820.56},
{"rc208.1.tw",509.04},
{"rc208.2.tw",503.92},
{"rc_201.1.tw",444.54},
{"rc_201.2.tw",711.54},
{"rc_201.3.tw",790.61},
{"rc_201.4.tw",793.64},
{"rc_202.1.tw",771.78},
{"rc_202.2.tw",304.14},
{"rc_202.3.tw",837.72},
{"rc_202.4.tw",793.03},
{"rc_203.1.tw",453.48},
{"rc_203.2.tw",784.16},
{"rc_203.3.tw",817.53},
{"rc_203.4.tw",314.29},
{"rc_204.1.tw",878.64},
{"rc_204.2.tw",662.16},
{"rc_204.3.tw",455.03},
{"rc_205.1.tw",343.21},
{"rc_205.2.tw",755.93},
{"rc_205.3.tw",825.06},
{"rc_205.4.tw",760.47},
{"rc_206.1.tw",117.85},
{"rc_206.2.tw",828.06},
{"rc_206.3.tw",574.42},
{"rc_206.4.tw",831.67},
{"rc_207.1.tw",732.68},
{"rc_207.2.tw",701.25},
{"rc_207.3.tw",682.4},
{"rc_207.4.tw",119.64},
{"rc_208.1.tw",789.25},
{"rc_208.2.tw",533.78},
{"rc_208.3.tw",634.44},
{"rbg010a.tw",671},
{"rbg016a.tw",938},
{"rbg016b.tw",1304},
{"rbg017.2.tw",852},
{"rbg017.tw",893},
{"rbg017a.tw",4296},
{"rbg019a.tw",1262},
{"rbg019b.tw",1866},
{"rbg019c.tw",4536},
{"rbg019d.tw",1356},
{"rbg020a.tw",4689},
{"rbg021.2.tw",4528},
{"rbg021.3.tw",4528},
{"rbg021.4.tw",4525},
{"rbg021.5.tw",4515},
{"rbg021.6.tw",4480},
{"rbg021.7.tw",4479},
{"rbg021.8.tw",4478},
{"rbg021.9.tw",4478},
{"rbg021.tw",4536},
{"rbg027a.tw",5091},
{"rbg031a.tw",1863},
{"rbg033a.tw",2069},
{"rbg034a.tw",2222},
{"rbg035a.2.tw",2056},
{"rbg035a.tw",2144},
{"rbg038a.tw",2480},
{"rbg040a.tw",2378},
{"rbg041a.tw",2598},
{"rbg042a.tw",2772},
{"rbg048a.tw",9383},
{"rbg049a.tw",10018},
{"rbg050a.tw",2953},
{"rbg050b.tw",9863},
{"rbg050c.tw",10024},
{"rbg055a.tw",3761},
{"rbg067a.tw",4625},
{"rbg086a.tw",8400},
{"rbg092a.tw",7158},
{"rbg125a.tw",7936},
{"rbg132.2.tw",8192},
{"rbg132.tw",8468},
{"rbg152.3.tw",9788},
{"rbg152.tw",10032},
{"rbg172a.tw",10950},
{"rbg193.2.tw",12142},
{"rbg193.tw",12535},
{"rbg201a.tw",12948},
{"rbg233.2.tw",14495},
{"rbg233.tw",14992},
{"n100w20.001.tw",738},
{"n100w20.002.tw",715},
{"n100w20.003.tw",762},
{"n100w20.004.tw",799},
{"n100w20.005.tw",774},
{"n100w40.001.tw",770},
{"n100w40.002.tw",653},
{"n100w40.003.tw",736},
{"n100w40.004.tw",651},
{"n100w40.005.tw",699},
{"n100w60.001.tw",655},
{"n100w60.002.tw",659},
{"n100w60.003.tw",744},
{"n100w60.004.tw",764},
{"n100w60.005.tw",661},
{"n150w20.001.tw",925},
{"n150w20.002.tw",864},
{"n150w20.003.tw",834},
{"n150w20.004.tw",873},
{"n150w20.005.tw",846},
{"n150w40.001.tw",918},
{"n150w40.002.tw",941},
{"n150w40.003.tw",727},
{"n150w40.004.tw",764},
{"n150w40.005.tw",824},
{"n150w60.001.tw",859},
{"n150w60.002.tw",782},
{"n150w60.003.tw",793},
{"n150w60.004.tw",819},
{"n150w60.005.tw",840},
{"n200w20.001.tw",1019},
{"n200w20.002.tw",972},
{"n200w20.003.tw",1050},
{"n200w20.004.tw",984},
{"n200w20.005.tw",1020},
{"n200w40.001.tw",1023},
{"n200w40.002.tw",948},
{"n200w40.003.tw",933},
{"n200w40.004.tw",980},
{"n200w40.005.tw",1037},
{"n20w100.001.tw",237},
{"n20w100.002.tw",222},
{"n20w100.003.tw",310},
{"n20w100.004.tw",349},
{"n20w100.005.tw",258},
{"n20w20.001.tw",378},
{"n20w20.002.tw",286},
{"n20w20.003.tw",394},
{"n20w20.004.tw",396},
{"n20w20.005.tw",352},
{"n20w40.001.tw",254},
{"n20w40.002.tw",333},
{"n20w40.003.tw",317},
{"n20w40.004.tw",388},
{"n20w40.005.tw",288},
{"n20w60.001.tw",335},
{"n20w60.002.tw",244},
{"n20w60.003.tw",352},
{"n20w60.004.tw",280},
{"n20w60.005.tw",338},
{"n20w80.001.tw",329},
{"n20w80.002.tw",338},
{"n20w80.003.tw",320},
{"n20w80.004.tw",304},
{"n20w80.005.tw",264},
{"n40w100.001.tw",429},
{"n40w100.002.tw",358},
{"n40w100.003.tw",364},
{"n40w100.004.tw",357},
{"n40w100.005.tw",377},
{"n40w20.001.tw",500},
{"n40w20.002.tw",552},
{"n40w20.003.tw",478},
{"n40w20.004.tw",404},
{"n40w20.005.tw",499},
{"n40w40.001.tw",465},
{"n40w40.002.tw",461},
{"n40w40.003.tw",474},
{"n40w40.004.tw",452},
{"n40w40.005.tw",453},
{"n40w60.001.tw",494},
{"n40w60.002.tw",470},
{"n40w60.003.tw",408},
{"n40w60.004.tw",382},
{"n40w60.005.tw",328},
{"n40w80.001.tw",395},
{"n40w80.002.tw",431},
{"n40w80.003.tw",412},
{"n40w80.004.tw",417},
{"n40w80.005.tw",344},
{"n60w100.001.tw",515},
{"n60w100.002.tw",538},
{"n60w100.003.tw",560},
{"n60w100.004.tw",510},
{"n60w100.005.tw",451},
{"n60w20.001.tw",551},
{"n60w20.002.tw",605},
{"n60w20.003.tw",533},
{"n60w20.004.tw",616},
{"n60w20.005.tw",603},
{"n60w40.001.tw",591},
{"n60w40.002.tw",621},
{"n60w40.003.tw",603},
{"n60w40.004.tw",597},
{"n60w40.005.tw",539},
{"n60w60.001.tw",609},
{"n60w60.002.tw",566},
{"n60w60.003.tw",485},
{"n60w60.004.tw",571},
{"n60w60.005.tw",569},
{"n60w80.001.tw",458},
{"n60w80.002.tw",498},
{"n60w80.003.tw",550},
{"n60w80.004.tw",566},
{"n60w80.005.tw",468},
{"n80w20.001.tw",616},
{"n80w20.002.tw",737},
{"n80w20.003.tw",667},
{"n80w20.004.tw",615},
{"n80w20.005.tw",748},
{"n80w40.001.tw",606},
{"n80w40.002.tw",618},
{"n80w40.003.tw",674},
{"n80w40.004.tw",557},
{"n80w40.005.tw",695},
{"n80w60.001.tw",554},
{"n80w60.002.tw",633},
{"n80w60.003.tw",651},
{"n80w60.004.tw",619},
{"n80w60.005.tw",575},
{"n80w80.001.tw",624},
{"n80w80.002.tw",592},
{"n80w80.003.tw",589},
{"n80w80.004.tw",594},
{"n80w80.005.tw",570},
{"n100w100.001.tw",643},
{"n100w100.002.tw",619},
{"n100w100.003.tw",685},
{"n100w100.004.tw",684},
{"n100w100.005.tw",572},
{"n100w120.001.tw",629},
{"n100w120.002.tw",540},
{"n100w120.003.tw",617},
{"n100w120.004.tw",663},
{"n100w120.005.tw",537},
{"n100w140.001.tw",604},
{"n100w140.002.tw",615},
{"n100w140.003.tw",481},
{"n100w140.004.tw",533},
{"n100w140.005.tw",509},
{"n100w160.001.tw",582},
{"n100w160.002.tw",532},
{"n100w160.003.tw",495},
{"n100w160.004.tw",580},
{"n100w160.005.tw",586},
{"n100w80.001.tw",670},
{"n100w80.002.tw",668},
{"n100w80.003.tw",691},
{"n100w80.004.tw",700},
{"n100w80.005.tw",603},
{"n20w120.001.tw",267},
{"n20w120.002.tw",218},
{"n20w120.003.tw",303},
{"n20w120.004.tw",300},
{"n20w120.005.tw",240},
{"n20w140.001.tw",176},
{"n20w140.002.tw",272},
{"n20w140.003.tw",236},
{"n20w140.004.tw",255},
{"n20w140.005.tw",225},
{"n20w160.001.tw",241},
{"n20w160.002.tw",201},
{"n20w160.003.tw",201},
{"n20w160.004.tw",203},
{"n20w160.005.tw",245},
{"n20w180.001.tw",253},
{"n20w180.002.tw",265},
{"n20w180.003.tw",271},
{"n20w180.004.tw",201},
{"n20w180.005.tw",193},
{"n20w200.001.tw",233},
{"n20w200.002.tw",203},
{"n20w200.003.tw",249},
{"n20w200.004.tw",293},
{"n20w200.005.tw",227},
{"n40w120.001.tw",434},
{"n40w120.002.tw",445},
{"n40w120.003.tw",357},
{"n40w120.004.tw",303},
{"n40w120.005.tw",350},
{"n40w140.001.tw",328},
{"n40w140.002.tw",383},
{"n40w140.003.tw",398},
{"n40w140.004.tw",342},
{"n40w140.005.tw",371},
{"n40w160.001.tw",348},
{"n40w160.002.tw",337},
{"n40w160.003.tw",346},
{"n40w160.004.tw",288},
{"n40w160.005.tw",315},
{"n40w180.001.tw",337},
{"n40w180.002.tw",347},
{"n40w180.003.tw",279},
{"n40w180.004.tw",354},
{"n40w180.005.tw",335},
{"n40w200.001.tw",330},
{"n40w200.002.tw",303},
{"n40w200.003.tw",339},
{"n40w200.004.tw",301},
{"n40w200.005.tw",296},
{"n60w120.001.tw",384},
{"n60w120.002.tw",427},
{"n60w120.003.tw",407},
{"n60w120.004.tw",490},
{"n60w120.005.tw",547},
{"n60w140.001.tw",423},
{"n60w140.002.tw",462},
{"n60w140.003.tw",427},
{"n60w140.004.tw",488},
{"n60w140.005.tw",460},
{"n60w160.001.tw",560},
{"n60w160.002.tw",423},
{"n60w160.003.tw",434},
{"n60w160.004.tw",401},
{"n60w160.005.tw",502},
{"n60w180.001.tw",411},
{"n60w180.002.tw",399},
{"n60w180.003.tw",445},
{"n60w180.004.tw",456},
{"n60w180.005.tw",395},
{"n60w200.001.tw",410},
{"n60w200.002.tw",414},
{"n60w200.003.tw",455},
{"n60w200.004.tw",431},
{"n60w200.005.tw",427},
{"n80w100.001.tw",565},
{"n80w100.002.tw",567},
{"n80w100.003.tw",580},
{"n80w100.004.tw",649},
{"n80w100.005.tw",532},
{"n80w120.001.tw",498},
{"n80w120.002.tw",577},
{"n80w120.003.tw",540},
{"n80w120.004.tw",501},
{"n80w120.005.tw",591},
{"n80w140.001.tw",512},
{"n80w140.002.tw",470},
{"n80w140.003.tw",580},
{"n80w140.004.tw",423},
{"n80w140.005.tw",545},
{"n80w160.001.tw",506},
{"n80w160.002.tw",549},
{"n80w160.003.tw",521},
{"n80w160.004.tw",509},
{"n80w160.005.tw",439},
{"n80w180.001.tw",551},
{"n80w180.002.tw",479},
{"n80w180.003.tw",524},
{"n80w180.004.tw",479},
{"n80w180.005.tw",470},
{"n80w200.001.tw",490},
{"n80w200.002.tw",488},
{"n80w200.003.tw",464},
{"n80w200.004.tw",526},
{"n80w200.005.tw",439},
{"N20ft301.tw",661.6},
{"N20ft302.tw",684.2},
{"N20ft303.tw",746.4},
{"N20ft304.tw",817},
{"N20ft305.tw",716.5},
{"N20ft306.tw",727.8},
{"N20ft307.tw",691.8},
{"N20ft308.tw",788.2},
{"N20ft309.tw",730.7},
{"N20ft310.tw",683},
{"N20ft401.tw",660.8},
{"N20ft402.tw",684.2},
{"N20ft403.tw",746.4},
{"N20ft404.tw",817},
{"N20ft405.tw",716.5},
{"N20ft406.tw",727.8},
{"N20ft407.tw",691.8},
{"N20ft408.tw",757.3},
{"N20ft409.tw",730.7},
{"N20ft410.tw",683},
{"N40ft201.tw",1100.6},
{"N40ft202.tw",1010.4},
{"N40ft203.tw",876.8},
{"N40ft204.tw",885.8},
{"N40ft205.tw",940.9},
{"N40ft206.tw",1054.2},
{"N40ft207.tw",867.5},
{"N40ft208.tw",1050.7},
{"N40ft209.tw",1013.9},
{"N40ft210.tw",1026.3},
{"N40ft401.tw",1085},
{"N40ft402.tw",995.6},
{"N40ft403.tw",845.8},
{"N40ft404.tw",868},
{"N40ft405.tw",936.5},
{"N40ft406.tw",969.1},
{"N40ft407.tw",831.2},
{"N40ft408.tw",1002.7},
{"N40ft409.tw",1000.5},
{"N40ft410.tw",983.8},
{"N60ft201.tw",1353.5},
{"N60ft202.tw",1161.6},
{"N60ft203.tw",1182.9},
{"N60ft204.tw",1257.5},
{"N60ft205.tw",1184.1},
{"N60ft206.tw",1199.6},
{"N60ft207.tw",1299},
{"N60ft208.tw",1113},
{"N60ft209.tw",1171.3},
{"N60ft210.tw",1234.3},
{"N60ft301.tw",1337},
{"N60ft302.tw",1089.5},
{"N60ft303.tw",1179},
{"N60ft304.tw",1230},
{"N60ft305.tw",1151.6},
{"N60ft306.tw",1167.9},
{"N60ft307.tw",1220.1},
{"N60ft308.tw",1097.6},
{"N60ft309.tw",1140.6},
{"N60ft310.tw",1219.2},
{"N60ft401.tw",1335},
{"N60ft402.tw",1088.1},
{"N60ft403.tw",1173.7},
{"N60ft404.tw",1184.7},
{"N60ft405.tw",1146.2},
{"N60ft406.tw",1140.2},
{"N60ft407.tw",1198.9},
{"N60ft408.tw",1029.4},
{"N60ft409.tw",1121.4},
{"N60ft410.tw",1189.6},
{"n150w120.001.tw",734},
{"n150w120.002.tw",677},
{"n150w120.003.tw",747},
{"n150w120.004.tw",763},
{"n150w120.005.tw",689},
{"n150w140.001.tw",762},
{"n150w140.002.tw",755},
{"n150w140.003.tw",613},
{"n150w140.004.tw",676},
{"n150w140.005.tw",663},
{"n150w160.001.tw",706},
{"n150w160.002.tw",711},
{"n150w160.003.tw",608},
{"n150w160.004.tw",672},
{"n150w160.005.tw",658},
{"n200w120.001.tw",799},
{"n200w120.002.tw",721},
{"n200w120.003.tw",880},
{"n200w120.004.tw",777},
{"n200w120.005.tw",841},
{"n200w140.001.tw",834},
{"n200w140.002.tw",760},
{"n200w140.003.tw",758},
{"n200w140.004.tw",816},
{"n200w140.005.tw",822},
{"rc201.0.tw",628.62},
{"rc201.1.tw",654.7},
{"rc201.2.tw",707.65},
{"rc201.3.tw",422.54},
{"rc202.0.tw",496.22},
{"rc202.1.tw",426.53},
{"rc202.2.tw",611.77},
{"rc202.3.tw",627.85},
{"rc203.0.tw",727.45},
{"rc203.1.tw",726.99},
{"rc203.2.tw",617.46},
{"rc204.0.tw",541.45},
{"rc204.1.tw",485.37},
{"rc204.2.tw",778.4},
{"rc205.0.tw",511.65},
{"rc205.1.tw",491.22},
{"rc205.2.tw",714.69},
{"rc205.3.tw",601.24},
{"rc206.0.tw",835.23},
{"rc206.1.tw",664.73},
{"rc206.2.tw",655.37},
{"rc207.0.tw",806.69},
{"rc207.1.tw",726.36},
{"rc207.2.tw",546.41},
{"rc208.0.tw",820.56},
{"rc208.1.tw",509.04},
{"rc208.2.tw",503.92},
{"rc_201.1.tw",444.54},
{"rc_201.2.tw",711.54},
{"rc_201.3.tw",790.61},
{"rc_201.4.tw",793.64},
{"rc_202.1.tw",771.78},
{"rc_202.2.tw",304.14},
{"rc_202.3.tw",837.72},
{"rc_202.4.tw",793.03},
{"rc_203.1.tw",453.48},
{"rc_203.2.tw",784.16},
{"rc_203.3.tw",817.53},
{"rc_203.4.tw",314.29},
{"rc_204.1.tw",878.64},
{"rc_204.2.tw",662.16},
{"rc_204.3.tw",455.03},
{"rc_205.1.tw",343.21},
{"rc_205.2.tw",755.93},
{"rc_205.3.tw",825.06},
{"rc_205.4.tw",760.47},
{"rc_206.1.tw",117.85},
{"rc_206.2.tw",828.06},
{"rc_206.3.tw",574.42},
{"rc_206.4.tw",831.67},
{"rc_207.1.tw",732.68},
{"rc_207.2.tw",701.25},
{"rc_207.3.tw",682.4},
{"rc_207.4.tw",119.64},
{"rc_208.1.tw",789.25},
{"rc_208.2.tw",533.78},
{"rc_208.3.tw",634.44}
};

struct VRPTW
{
  public:
    VRPTW(std::string fileName)
    {
      // VRPTW instance format
      std::regex vrptwInstanceRegex(".*txt");
      std::smatch vrptwInstanceMatch;
      if (std::regex_search(fileName, vrptwInstanceMatch, vrptwInstanceRegex))
      {
        oneOrMorePaths = OneOrMorePaths::MORE_PATHS;
        std::vector<std::pair<double,double> > coordinates;

        std::regex noCapacityRelaxRegex(".*HG.*C1|R1|RC1.*txt");
        std::smatch noCapacityRelaxMatch;
        vrptwType = VRPTWType::RELAX_CAPACITY;
        if (std::regex_search(fileName, noCapacityRelaxMatch, noCapacityRelaxRegex))
        {
          std::cout << "do not relax capacity constraint" << std::endl;
          vrptwType = VRPTWType::NO_RELAX_CAPACITY;
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
        oneOrMorePaths = OneOrMorePaths::MORE_PATHS;
        vrptwType = VRPTWType::RELAX_TIME_WINDOWS;
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
      std::regex tsptwInstanceRegex(".*tw");
      std::smatch tsptwInstanceMatch;
      if (std::regex_search(fileName, tsptwInstanceMatch, tsptwInstanceRegex))
      {
        oneOrMorePaths = OneOrMorePaths::ONE_PATH;
        std::vector<std::pair<double,double> > coordinates;
        vrptwType = VRPTWType::RELAX_CAPACITY;

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
              if ((distanceMatrixRow != column) && (distance == 0.0) && (distanceMatrixRow != 0))
              {
                std::cout << "0 distance will be bad for input, row: " << distanceMatrixRow << " col: " << column << std::endl;
                return;
              }
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
    std::vector<int> demands;
    std::vector<double> demandsForSeparation;
    std::vector<int> demandsForCombs;

    std::vector<int> startTimes;
    std::vector<int> endTimes;
    std::vector<int> serviceTimes;

    int depot;
    int numLocations;

    VRPTWType vrptwType;
    OneOrMorePaths oneOrMorePaths;
    double hgsUpperBound;
};

#endif
