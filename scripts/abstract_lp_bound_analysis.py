# Create table to compare CE with SOA
import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections

instance_upper_bounds = {
"A-n32-k5.vrp":784,
"A-n33-k5.vrp":661,
"A-n33-k6.vrp":742,
"A-n34-k5.vrp":778,
"A-n36-k5.vrp":799,
"A-n37-k5.vrp":669,
"A-n37-k6.vrp":949,
"A-n38-k5.vrp":730,
"A-n39-k5.vrp":822,
"A-n39-k6.vrp":831,
"A-n44-k6.vrp":937,
"A-n45-k6.vrp":944,
"A-n45-k7.vrp":1146,
"A-n46-k7.vrp":914,
"A-n48-k7.vrp":1073,
"A-n53-k7.vrp":1010,
"A-n54-k7.vrp":1167,
"A-n55-k9.vrp":1073,
"A-n60-k9.vrp":1354,
"A-n61-k9.vrp":1034,
"A-n62-k8.vrp":1288,
"A-n63-k10.vrp":1314,
"A-n63-k9.vrp":1616,
"A-n64-k9.vrp":1401,
"A-n65-k9.vrp":1174,
"A-n69-k9.vrp":1159,
"A-n80-k10.vrp":1763,
"B-n31-k5.vrp":672,
"B-n34-k5.vrp":788,
"B-n35-k5.vrp":955,
"B-n38-k6.vrp":805,
"B-n39-k5.vrp":549,
"B-n41-k6.vrp":829,
"B-n43-k6.vrp":742,
"B-n44-k7.vrp":909,
"B-n45-k5.vrp":751,
"B-n45-k6.vrp":678,
"B-n50-k7.vrp":741,
"B-n50-k8.vrp":1312,
"B-n51-k7.vrp":1032,
"B-n52-k7.vrp":747,
"B-n56-k7.vrp":707,
"B-n57-k7.vrp":1153,
"B-n57-k9.vrp":1598,
"B-n63-k10.vrp":1496,
"B-n64-k9.vrp":861,
"B-n66-k9.vrp":1316,
"B-n67-k10.vrp":1032,
"B-n68-k9.vrp":1272,
"B-n78-k10.vrp":1221,
"E-n101-k14.vrp":1071,
"E-n101-k8.vrp":817,
"E-n13-k4.vrp":247,
"E-n22-k4.vrp":375,
"E-n23-k3.vrp":569,
"E-n30-k3.vrp":534,
"E-n31-k7.vrp":379,
"E-n33-k4.vrp":835,
"E-n51-k5.vrp":521,
"E-n76-k10.vrp":830,
"E-n76-k14.vrp":1021,
"E-n76-k7.vrp":682,
"E-n76-k8.vrp":735,
"F-n135-k7.vrp":1162,
"F-n45-k4.vrp":724,
"F-n72-k4.vrp":237,
"P-n101-k4.vrp":681,
"P-n16-k8.vrp":450,
"P-n19-k2.vrp":212,
"P-n20-k2.vrp":216,
"P-n21-k2.vrp":211,
"P-n22-k2.vrp":216,
"P-n22-k8.vrp":603,
"P-n23-k8.vrp":529,
"P-n40-k5.vrp":458,
"P-n45-k5.vrp":510,
"P-n50-k10.vrp":696,
"P-n50-k7.vrp":554,
"P-n50-k8.vrp":631,
"P-n51-k10.vrp":741,
"P-n55-k10.vrp":694,
"P-n55-k15.vrp":989,
"P-n55-k7.vrp":568,
"P-n55-k8.vrp":588,
"P-n60-k10.vrp":744,
"P-n60-k15.vrp":968,
"P-n65-k10.vrp":792,
"P-n70-k10.vrp":827,
"P-n76-k4.vrp":593,
"P-n76-k5.vrp":627,
"M-n101-k10.vrp":820,
"M-n121-k7.vrp":1034,
"M-n151-k12.vrp":1053,
"M-n200-k16.vrp":1373,
"M-n200-k17.vrp":1373
}

# Log lines for DDSolver
#STATS - lpIterations[1] lagIterations[32] sspIterations[186] numSeparations[0] compileTime[240] sspSolveTime[349] lpSolveTime[0] lb[794.025] ub[1e+10] numArcs: [87254] numFixed: [0] time: [590]

colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\] numCuts\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] numArcs: \[([0-9]+)\] numFixed: \[([0-9]+)\] time: \[([0-9]+)\]")
#colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] size: \[([0-9]+)\] time: \[([0-9]+)\]")

logs_dir = "/Users/akarahal/Desktop/dd_vrptw/logs/"
test_set = ["CVRP_LAG_NG2_src_phase", "CVRP_LAG_NG2_rcc_phase"]

instances = []
instance_dir = "/Users/akarahal/Desktop/dd_vrptw/instances/CVRP/"
instances_to_lp_opt = []

for instance in os.listdir(instance_dir):
  instances.append(instance)

time_results = []
results = []
for test in test_set:
  for instance in instances:
    log_file_name = logs_dir + test + '/' + instance + '.log'
    if not os.path.exists(log_file_name):
      continue

    col_elim = False
    log_file = open(log_file_name, "r")
    for line in log_file:
      # only go until end of LP solver
      if "switching to IP" in line:
        #print("stop")
        instances_to_lp_opt.append((instance,test))
        break

      colelim_match = colelim_pattern.match(line)
      if colelim_match:
        col_elim = True
        lpIterations = int(colelim_match.group(1))
        lagIterations = int(colelim_match.group(2))
        sspIterations = colelim_match.group(3)
        numSeparations = int(colelim_match.group(4))
        numCuts= int(colelim_match.group(5))
        compileTime = colelim_match.group(6)
        sspSolveTime = colelim_match.group(7)
        lpSolveTime = colelim_match.group(8)
        lb = float(colelim_match.group(9))
        ub = float(colelim_match.group(10))
        numArcs = int(colelim_match.group(11))
        numFixed = int(colelim_match.group(12))
        time = float(colelim_match.group(13))
        time_result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'lagIterations': lagIterations, 'numSep' : numSeparations, 'numCut': numCuts, 'lb' : lb, 'ub' : ub, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
        time_results.append(time_result)
 
    if col_elim:
      result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'lagIterations': lagIterations, 'numSep' : numSeparations, 'numCut' : numCuts, 'lb' : lb, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
      results.append(result)
    else:
      result = {'instance': instance, 'test': test, 'iterations' : numpy.nan, 'lagIterations': numpy.nan, 'numSep' : numpy.nan, 'numCut': numpy.nan, 'lb' : numpy.nan, 'numArcs': numpy.nan, 'numFixed': numpy.nan, 'time' : numpy.nan}
      results.append(result)

# give table of results
results_df = pandas.DataFrame(results)
results_df.sort_values(by=['instance','lb'],inplace=True)
results_df.reset_index(drop=True,inplace=True)
table_results_df = results_df[['instance','test','lb','time','iterations','numSep','numCut','lagIterations','numArcs']]
print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))

# analysis
for instance in instances:
  if (instance, "CVRP_LAG_NG2_rcc_phase") in instances_to_lp_opt:
    instance_df = table_results_df[table_results_df['instance'] == instance]
    src_df = instance_df[instance_df['test'] == "CVRP_LAG_NG2_src_phase"]
    rcc_df = instance_df[instance_df['test'] == "CVRP_LAG_NG2_rcc_phase"]
    src_lb = src_df['lb'].values[0]
    rcc_lb = rcc_df['lb'].values[0]
    if src_lb > rcc_lb:
      print(instance, src_lb, rcc_lb)
