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

enum VRPTWType
{
  RELAX_CAPACITY = 0,
  NO_RELAX_CAPACITY = 1
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
{"lou-n88-k2.vrp",41}
};

struct VRPTW
{
  public:
    VRPTW(std::string fileName)
    {
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
    double hgsUpperBound;
};

#endif
