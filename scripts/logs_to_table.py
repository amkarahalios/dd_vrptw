import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections

instance_upper_bounds = {
"C101.txt":827.3,
"C102.txt":827.3,
"C103.txt":826.3,
"C104.txt":822.9,
"C105.txt":827.3,
"C106.txt":827.3,
"C107.txt":827.3,
"C108.txt":827.3,
"C109.txt":827.3,
"C201.txt":589.1,
"C202.txt":589.1,
"C203.txt":588.7,
"C204.txt":588.1,
"C205.txt":586.4,
"C206.txt":586,
"C207.txt":585.8,
"C208.txt":585.8,
"R101.txt":1637.7,
"R102.txt":1466.6,
"R103.txt":1208.7,
"R104.txt":971.5,
"R105.txt":1355.3,
"R106.txt":1234.6,
"R107.txt":1064.6,
"R108.txt":932.1,
"R109.txt":1146.9,
"R110.txt":1068,
"R111.txt":1048.7,
"R112.txt":948.6,
"R201.txt":1143.2,
"R202.txt":1029.6,
"R203.txt":870.8,
"R204.txt":731.3,
"R205.txt":949.8,
"R206.txt":875.9,
"R207.txt":794,
"R208.txt":701,
"R209.txt":854.8,
"R210.txt":900.5,
"R211.txt":746.7,
"RC101.txt":1619.8,
"RC102.txt":1457.4,
"RC103.txt":1258,
"RC104.txt":1132.3,
"RC105.txt":1513.7,
"RC106.txt":1372.7,
"RC107.txt":1207.8,
"RC108.txt":1114.2,
"RC201.txt":1261.8,
"RC202.txt":1092.3,
"RC203.txt":923.7,
"RC204.txt":783.5,
"RC205.txt":1154,
"RC206.txt":1051.1,
"RC207.txt":962.9,
"RC208.txt":776.1,
"C1_2_1.txt":2698.6,
"C1_2_2.txt":2694.3,
"C1_2_3.txt":2675.8,
"C1_2_4.txt":2625.6,
"C1_2_5.txt":2694.9,
"C1_2_6.txt":2694.9,
"C1_2_7.txt":2694.9,
"C1_2_8.txt":2684,
"C1_2_9.txt":2639.6,
"C1_2_10.txt":2624.7,
"C2_2_1.txt":1922.1,
"C2_2_2.txt":1851.4,
"C2_2_3.txt":1763.4,
"C2_2_4.txt":1695,
"C2_2_5.txt":1869.6,
"C2_2_6.txt":1844.8,
"C2_2_7.txt":1842.2,
"C2_2_8.txt":1813.7,
"C2_2_9.txt":1815,
"C2_2_10.txt":1791.2,
"R1_2_1.txt":4667.2,
"R1_2_2.txt":3919.9,
"R1_2_3.txt":3373.9,
"R1_2_4.txt":3047.6,
"R1_2_5.txt":4053.2,
"R1_2_6.txt":3559.1,
"R1_2_7.txt":3141.9,
"R1_2_8.txt":2938.4,
"R1_2_9.txt":3734.7,
"R1_2_10.txt":3293.1,
"R2_2_1.txt":3468,
"R2_2_2.txt":3008.2,
"R2_2_3.txt":2537.5,
"R2_2_4.txt":1928.5,
"R2_2_5.txt":3061.1,
"R2_2_6.txt":2675.4,
"R2_2_7.txt":2304.7,
"R2_2_8.txt":1842.4,
"R2_2_9.txt":2843.3,
"R2_2_10.txt":2549.4,
"RC1_2_1.txt":3516.9,
"RC1_2_2.txt":3221.6,
"RC1_2_3.txt":3001.4,
"RC1_2_4.txt":2845.2,
"RC1_2_5.txt":3325.6,
"RC1_2_6.txt":3300.7,
"RC1_2_7.txt":3177.8,
"RC1_2_8.txt":3060,
"RC1_2_9.txt":3073.3,
"RC1_2_10.txt":2990.5,
"RC2_2_1.txt":2797.4,
"RC2_2_2.txt":2481.6,
"RC2_2_3.txt":2227.7,
"RC2_2_4.txt":1854.8,
"RC2_2_5.txt":2491.4,
"RC2_2_6.txt":2495.1,
"RC2_2_7.txt":2287.7,
"RC2_2_8.txt":2151.2,
"RC2_2_9.txt":2086.6,
"RC2_2_10.txt":1989.2,
"C1_4_1.txt":7138.8,
"C1_4_2.txt":7113.3,
"C1_4_3.txt":6929.9,
"C1_4_4.txt":6777.7,
"C1_4_5.txt":7138.8,
"C1_4_6.txt":7140.1,
"C1_4_7.txt":7136.2,
"C1_4_8.txt":7083,
"C1_4_9.txt":6927.8,
"C1_4_10.txt":6825.4,
"C2_4_1.txt":4100.3,
"C2_4_2.txt":3914.1,
"C2_4_3.txt":3755.2,
"C2_4_4.txt":3523.7,
"C2_4_5.txt":3923.2,
"C2_4_6.txt":3860.1,
"C2_4_7.txt":3870.9,
"C2_4_8.txt":3773.7,
"C2_4_9.txt":3842.1,
"C2_4_10.txt":3665.1,
"R1_4_1.txt":10305.8,
"R1_4_2.txt":8873.3,
"R1_4_3.txt":7784.3,
"R1_4_4.txt":7266.2,
"R1_4_5.txt":9184.6,
"R1_4_6.txt":8340.4,
"R1_4_7.txt":7599.8,
"R1_4_8.txt":7240.5,
"R1_4_9.txt":8677.5,
"R1_4_10.txt":8077.8,
"R2_4_1.txt":7520.7,
"R2_4_2.txt":6482.8,
"R2_4_3.txt":5372.9,
"R2_4_4.txt":4211.2,
"R2_4_5.txt":6567.9,
"R2_4_6.txt":5813.5,
"R2_4_7.txt":4893.5,
"R2_4_8.txt":4000.1,
"R2_4_9.txt":6067.8,
"R2_4_10.txt":5645.9,
"RC1_4_1.txt":8522.9,
"RC1_4_2.txt":7878.2,
"RC1_4_3.txt":7516.9,
"RC1_4_4.txt":7292.9,
"RC1_4_5.txt":8152.3,
"RC1_4_6.txt":8148,
"RC1_4_7.txt":7932.5,
"RC1_4_8.txt":7757.2,
"RC1_4_9.txt":7717.7,
"RC1_4_10.txt":7581.2,
"RC2_4_1.txt":6147.3,
"RC2_4_2.txt":5407.5,
"RC2_4_3.txt":4573,
"RC2_4_4.txt":3597.9,
"RC2_4_5.txt":5392.3,
"RC2_4_6.txt":5324.6,
"RC2_4_7.txt":4987.8,
"RC2_4_8.txt":4693.3,
"RC2_4_9.txt":4510.4,
"RC2_4_10.txt":4252.3,
"C1_6_1.txt":14076.6,
"C1_6_2.txt":13948.3,
"C1_6_3.txt":13757,
"C1_6_4.txt":13538.6,
"C1_6_5.txt":14066.8,
"C1_6_6.txt":14070.9,
"C1_6_7.txt":14066.8,
"C1_6_8.txt":13991.2,
"C1_6_9.txt":13664.5,
"C1_6_10.txt":13617.5,
"C2_6_1.txt":7752.2,
"C2_6_2.txt":7471.5,
"C2_6_3.txt":7215,
"C2_6_4.txt":6877,
"C2_6_5.txt":7553.8,
"C2_6_6.txt":7449.8,
"C2_6_7.txt":7491.3,
"C2_6_8.txt":7303.7,
"C2_6_9.txt":7303.2,
"C2_6_10.txt":7123.9,
"R1_6_1.txt":21274.2,
"R1_6_2.txt":18558.7,
"R1_6_3.txt":16874.9,
"R1_6_4.txt":15721.4,
"R1_6_5.txt":19294.9,
"R1_6_6.txt":17763.7,
"R1_6_7.txt":16496.2,
"R1_6_8.txt":15584.3,
"R1_6_9.txt":18474.1,
"R1_6_10.txt":17583.7,
"R2_6_1.txt":15145.3,
"R2_6_2.txt":12976.3,
"R2_6_3.txt":10455.3,
"R2_6_4.txt":7915.1,
"R2_6_5.txt":13790.2,
"R2_6_6.txt":11847.8,
"R2_6_7.txt":9777.9,
"R2_6_8.txt":7512.3,
"R2_6_9.txt":12736.8,
"R2_6_10.txt":11837,
"RC1_6_1.txt":16960.1,
"RC1_6_2.txt":15890.6,
"RC1_6_3.txt":15181.3,
"RC1_6_4.txt":14753.2,
"RC1_6_5.txt":16536.3,
"RC1_6_6.txt":16473.3,
"RC1_6_7.txt":16055.3,
"RC1_6_8.txt":15891.8,
"RC1_6_9.txt":15803.5,
"RC1_6_10.txt":15651.3,
"RC2_6_1.txt":11966.1,
"RC2_6_2.txt":10336.9,
"RC2_6_3.txt":8894.9,
"RC2_6_4.txt":6967.5,
"RC2_6_5.txt":11080.7,
"RC2_6_6.txt":10830.5,
"RC2_6_7.txt":10289.4,
"RC2_6_8.txt":9779,
"RC2_6_9.txt":9436,
"RC2_6_10.txt":8974.7,
"C1_8_1.txt":25156.9,
"C1_8_2.txt":24974.1,
"C1_8_3.txt":24156.1,
"C1_8_4.txt":23797.3,
"C1_8_5.txt":25138.6,
"C1_8_6.txt":25133.3,
"C1_8_7.txt":25127.3,
"C1_8_8.txt":24809.7,
"C1_8_9.txt":24200.4,
"C1_8_10.txt":24026.7,
"C2_8_1.txt":11631.9,
"C2_8_2.txt":11394.5,
"C2_8_3.txt":11138.1,
"C2_8_4.txt":10650,
"C2_8_5.txt":11395.6,
"C2_8_6.txt":11316.3,
"C2_8_7.txt":11332.9,
"C2_8_8.txt":11133.9,
"C2_8_9.txt":11140.4,
"C2_8_10.txt":10946,
"R1_8_1.txt":36345,
"R1_8_2.txt":32277.6,
"R1_8_3.txt":29304.5,
"R1_8_4.txt":27734.7,
"R1_8_5.txt":33494.2,
"R1_8_6.txt":30872.4,
"R1_8_7.txt":28789,
"R1_8_8.txt":27609.4,
"R1_8_9.txt":32257.3,
"R1_8_10.txt":30918.4,
"R2_8_1.txt":24969.8,
"R2_8_2.txt":21312.2,
"R2_8_3.txt":17234.8,
"R2_8_4.txt":13160.8,
"R2_8_5.txt":22801.6,
"R2_8_6.txt":19740.5,
"R2_8_7.txt":16357.5,
"R2_8_8.txt":12611.7,
"R2_8_9.txt":21282.7,
"R2_8_10.txt":19984.8,
"RC1_8_1.txt":29978.9,
"RC1_8_2.txt":28290.1,
"RC1_8_3.txt":27447.7,
"RC1_8_4.txt":26557.2,
"RC1_8_5.txt":29219.9,
"RC1_8_6.txt":29194.2,
"RC1_8_7.txt":28788.6,
"RC1_8_8.txt":28418.1,
"RC1_8_9.txt":28347.1,
"RC1_8_10.txt":28168.5,
"RC2_8_1.txt":19201.3,
"RC2_8_2.txt":16709.5,
"RC2_8_3.txt":14013.6,
"RC2_8_4.txt":10969.4,
"RC2_8_5.txt":17466.1,
"RC2_8_6.txt":17195.1,
"RC2_8_7.txt":16362.2,
"RC2_8_8.txt":15528.8,
"RC2_8_9.txt":15183,
"RC2_8_10.txt":14370.9,
"C1_10_1.txt":42444.8,
"C1_10_2.txt":41337.8,
"C1_10_3.txt":40064.4,
"C1_10_4.txt":39434.1,
"C1_10_5.txt":42434.8,
"C1_10_6.txt":42437,
"C1_10_7.txt":42420.4,
"C1_10_8.txt":41652.1,
"C1_10_9.txt":40288.4,
"C1_10_10.txt":39816.8,
"C2_10_1.txt":16841.1,
"C2_10_2.txt":16462.6,
"C2_10_3.txt":16036.5,
"C2_10_4.txt":15459.5,
"C2_10_5.txt":16521.3,
"C2_10_6.txt":16290.7,
"C2_10_7.txt":16378.4,
"C2_10_8.txt":16029.1,
"C2_10_9.txt":16075.4,
"C2_10_10.txt":15728.6,
"R1_10_1.txt":53046.5,
"R1_10_2.txt":48263.1,
"R1_10_3.txt":44677.1,
"R1_10_4.txt":42440.7,
"R1_10_5.txt":50406.7,
"R1_10_6.txt":46930.3,
"R1_10_7.txt":43997.4,
"R1_10_8.txt":42279.3,
"R1_10_9.txt":49162.8,
"R1_10_10.txt":47364.6,
"R2_10_1.txt":36881,
"R2_10_2.txt":31241.9,
"R2_10_3.txt":24399,
"R2_10_4.txt":17811.5,
"R2_10_5.txt":34132.8,
"R2_10_6.txt":29124.7,
"R2_10_7.txt":23102.2,
"R2_10_8.txt":17403.8,
"R2_10_9.txt":31990.6,
"R2_10_10.txt":29840.5,
"RC1_10_1.txt":45790.8,
"RC1_10_2.txt":43678.3,
"RC1_10_3.txt":42122,
"RC1_10_4.txt":41357.4,
"RC1_10_5.txt":45028.1,
"RC1_10_6.txt":44903.6,
"RC1_10_7.txt":44417.1,
"RC1_10_8.txt":43916.5,
"RC1_10_9.txt":43858.1,
"RC1_10_10.txt":43533.7,
"RC2_10_1.txt":28122.6,
"RC2_10_2.txt":24248.6,
"RC2_10_3.txt":19618.1,
"RC2_10_4.txt":15657,
"RC2_10_5.txt":25797.5,
"RC2_10_6.txt":25782.5,
"RC2_10_7.txt":24395.8,
"RC2_10_8.txt":23280.2,
"RC2_10_9.txt":22731.6,
"RC2_10_10.txt":21736.1
}

held_bounds = {
  "zeroin.i.3.col":30
}

# Log lines for DDSolver
#STATS - lpIterations[1] lagIterations[32] sspIterations[186] numSeparations[0] compileTime[240] sspSolveTime[349] lpSolveTime[0] lb[794.025] ub[1e+10] size: [87254] time: [590]

# Log lines for ColGenSolver
#STATS - iterations[402] ddTime[50] pricingTime[112] masterTime[97] lb[1581.24] time: [260]

colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\] compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] size: \[([0-9]+)\] time: \[([0-9]+)\]")
colgen_pattern = re.compile("STATS - iterations\[([0-9]+)\] ddTime\[([0-9]+)\] pricingTime\[([0-9]+)\] masterTime\[([0-9]+)\] lb\[([0-9]+.*)\] time: \[([0-9]+)\]")
colgen_pattern1 = re.compile("STATS - iterations\[([0-9]+)\] ddTime\[([0-9]+)\] pricingTime\[([0-9]+)\] masterTime\[([0-9]+)\] lb\[([0-9]+.*)\] sol\[([0-9]+.*)\] time: \[([0-9]+)\]")
colgen_complete_pattern = re.compile("no negative reduced cost paths")
colgen_pattern_size = re.compile("DD size: ([0-9]+)")
held_initial_pattern = re.compile("Finished initial bounds: LB ([0-9]+) and UB ([0-9]+) in.*([0-9]+)\.0000.*seconds.")
held_update_time_pattern = re.compile("Compute_coloring took ([0-9]+)..* seconds.*")
held_improved_lb_pattern = re.compile("Lower bound improved: LB ([0-9]+) and UB ([0-9]+).*")
held_improved_ub_pattern = re.compile("Upper bound improved: LB ([0-9]+) and UB ([0-9]+).*")
bdd_pattern = re.compile(".*Time elapsed: ([0-9]+).*BDD.*Lower bound: ([0-9]+).*Upper bound: ([0-9]+).*")

root_node_pattern = re.compile(".*Time elapsed: ([0-9]+).*BDD.*Lower bound: ([0-9]+).*Upper bound: ([0-9]+).*")
time_pattern = re.compile("Current Bounds - LB.*UB.* Time: ([0-9]+).*")
lb_pattern = re.compile("Current Bounds - LB: ([0-9]+) UB:.*")
ub_pattern = re.compile("Current Bounds - LB:.*UB: ([0-9]+).*")
finish_pattern = re.compile("Done solving.*time:\[([0-9]+).*\] LB:\[([0-9]+)\] UB:\[([0-9]+)\]")

logs_dir = "/Users/akarahal/Desktop/dd_vrptw/logs/"
test_set = ["col_elim_hg_lag_ng_4_20_N_3600"]

instances = []
instance_dir = "/Users/akarahal/Desktop/dd_vrptw/instances/"
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
    col_gen = False
    col_gen_complete = False
    held = False
    log_file = open(log_file_name, "r")
    for line in log_file:
      colelim_match = colelim_pattern.match(line)
      if colelim_match:
        col_elim = True
        lpIterations = colelim_match.group(1)
        lagIterations = colelim_match.group(2)
        sspIterations = colelim_match.group(3)
        numSeparations = colelim_match.group(4)
        compileTime = colelim_match.group(5)
        sspSolveTime = colelim_match.group(6)
        lpSolveTime = colelim_match.group(7)
        lb = math.ceil(float(colelim_match.group(8)))
        ub = int(colelim_match.group(9))
        size = colelim_match.group(10)
        time = colelim_match.group(11)
        time_result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'numSep' : numSeparations, 'lb' : lb, 'ub' : ub, 'size': size, 'time' : time}
        time_results.append(time_result)
 
      colgen_match = colgen_pattern.match(line)
      colgen_match1 = colgen_pattern1.match(line)
      if colgen_match and not colgen_match1:
        col_gen = True
        iterations = colgen_match.group(1)
        ddTime = colgen_match.group(2)
        pricingTime = colgen_match.group(3)
        masterTime = colgen_match.group(4)
        lb = math.ceil(float(colgen_match.group(5)))
        time = colgen_match.group(6)
        # remove DD time?
        time = float(time) - float(ddTime)
 
      if colgen_match1:
        col_gen = True
        iterations = colgen_match1.group(1)
        ddTime = colgen_match1.group(2)
        pricingTime = colgen_match1.group(3)
        masterTime = colgen_match1.group(4)
        lb = math.ceil(float(colgen_match1.group(5)))
        sol = float(colgen_match1.group(6))
        time = colgen_match1.group(7)
        # remove DD time?
        time = float(time) - float(ddTime)
        time_result = {'instance': instance, 'test': test, 'iterations' : iterations, 'lb' : lb, 'size' : size, 'time' : time}
        time_results.append(time_result)

      colgen_size_match = colgen_pattern_size.match(line)
      if colgen_size_match:
        size = colgen_size_match.group(1)

      colgen_complete_match = colgen_complete_pattern.match(line)
      if colgen_complete_match:
        col_gen_complete = True
 
      held_initial_match = held_initial_pattern.match(line)
      if held_initial_match:
        held = True
        lb = float(held_initial_match.group(1))
        ub = float(held_initial_match.group(2))
        time = float(held_initial_match.group(3))
        time_result = {'instance': instance, 'test': test, 'lb' : lb, 'ub' : ub, 'time' : time}
        time_results.append(time_result)
 
      held_update_time_match = held_update_time_pattern.match(line)
      if held_update_time_match:
        time = float(held_update_time_match.group(1))

      held_improved_lb_match = held_improved_lb_pattern.match(line)
      if held_improved_lb_match:
        lb = float(held_improved_lb_match.group(1))
        ub = float(held_improved_lb_match.group(2))
        time_result = {'instance': instance, 'test': test, 'lb' : lb, 'ub' : ub, 'time' : time}
        time_results.append(time_result)
 
      held_improved_ub_match = held_improved_ub_pattern.match(line)
      if held_improved_ub_match:
        lb = float(held_improved_ub_match.group(1))
        ub = float(held_improved_ub_match.group(2))
        time_result = {'instance': instance, 'test': test, 'lb' : lb, 'ub' : ub, 'time' : time}
        time_results.append(time_result)

      bdd_match = bdd_pattern.match(line)
      if bdd_match:
        time = float(bdd_match.group(1))
        lb = float(bdd_match.group(2))
        ub = float(bdd_match.group(3))
        time_result = {'instance': instance, 'test': test, 'lb' : lb, 'ub' : ub, 'time' : time}
        time_results.append(time_result)

    if col_gen and col_gen_complete:
      result = {'instance': instance, 'test': test, 'iterations' : iterations, 'lb' : lb, 'size': size, 'time' : time}
      results.append(result)
      time_results.append(result)
    elif col_elim:
      result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'numSep' : numSeparations, 'lb' : lb, 'size': size, 'time' : time}
      results.append(result)
    elif held:
      result = {'instance': instance, 'test': test, 'lb' : lb, 'time' : time}
      results.append(result)
    else:
      result = {'instance': instance, 'test': test, 'iterations' : numpy.nan, 'numSep' : numpy.nan, 'lb' : numpy.nan, 'size': numpy.nan, 'time' : numpy.nan}
      results.append(result)

# setup results df
results_df = pandas.DataFrame(results)
results_df.sort_values(by=['instance','lb'],inplace=True)
results_df.reset_index(drop=True,inplace=True)
time_results_df = pandas.DataFrame(time_results)
time_results_df.sort_values(by=['instance','lb'],inplace=True)
time_results_df.reset_index(drop=True,inplace=True)

# print table
print_table = True
if print_table:
  table_results_df = results_df[['instance','test','lb','time']]
  print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))

# experiment 1 - Exact DD Size vs. DD Elim Size
experiment1 = False
if experiment1:
  x = []
  y = []
  label = []
  for instance in set(results_df['instance']):
    instance_results = results_df[results_df['instance'] == instance]
    instance_tests = set(instance_results['test'])
    #if (("test_A_colelim_lag_q_2_2_1800" in instance_tests) and ("test_A_colelim_lp_q_1_2_1800" in instance_tests) and ("test_A_colgen_dd_q_2_1800" in instance_tests)):
    #  colelim_q_2_2 = instance_results[instance_results['test'] == "test_A_colelim_lag_q_2_2_1800"].iloc[0]
    #  colelim_q_1_2 = instance_results[instance_results['test'] == "test_A_colelim_lp_q_1_2_1800"].iloc[0]
    #  colgen_q_2_2 = instance_results[instance_results['test'] == "test_A_colgen_dd_q_2_1800"].iloc[0]
    #if (("A_colelim_lp_q_1_2_1800" in instance_tests) and ("A_colelim_lp_q_2_2_1800" in instance_tests) and ("A_colgen_dd_q_2_1800" in instance_tests)):
    if (("A_colelim_lp_q_1_2_1800" in instance_tests) and ("A_colgen_dd_q_2_1800" in instance_tests)):
      colelim_q_1_2 = instance_results[instance_results['test'] == "A_colelim_lp_q_1_2_1800"].iloc[0]
      colgen_q_2_2 = instance_results[instance_results['test'] == "A_colgen_dd_q_2_1800"].iloc[0]
      if (colgen_q_2_2['lb'] == colelim_q_1_2['lb']):
        exact_size = int(colgen_q_2_2['size'])
        elim_size = int(colelim_q_1_2['size'])
        x.append(exact_size)
        y.append(elim_size)
        label.append(instance)
        plt.plot([exact_size],[elim_size],marker='o',label=instance)
    if (("B_colelim_lp_q_1_2_1800" in instance_tests) and ("B_colgen_dd_q_2_1800" in instance_tests)):
      colelim_q_1_2 = instance_results[instance_results['test'] == "B_colelim_lp_q_1_2_1800"].iloc[0]
      colgen_q_2_2 = instance_results[instance_results['test'] == "B_colgen_dd_q_2_1800"].iloc[0]
      if (colgen_q_2_2['lb'] == colelim_q_1_2['lb']):
        exact_size = int(colgen_q_2_2['size'])
        elim_size = int(colelim_q_1_2['size'])
        x.append(exact_size)
        y.append(elim_size)
        label.append(instance)
        plt.plot([exact_size],[elim_size],marker='o',label=instance)
  plt.xlabel('exact size')
  plt.ylabel('elim size at opt')
  plt.xscale('log')
  plt.yscale('log')
  plt.xlim([10**2,10**5])
  plt.ylim([10**2,10**5])
  #plt.legend()
  plt.title('q=1->2 vs. q=2')
  plt.show()

# experiment 2 - plot the lb versus time for each method
experiment2 = True
if experiment2:
  # used these for paper exp2
  # using these to test X instances
  tests_to_compare = ["col_elim_hg_lag_ng_4_20_N_3600"]
  #instances_to_consider = ["X-n327-k20.vrp","X-n344-k43.vrp","X-n359-k29.vrp","X-n367-k17.vrp","X-n480-k70.vrp","X-n502-k39.vrp","X-n561-k42.vrp","X-n573-k30.vrp","X-n801-k40.vrp","X-n957-k87.vrp"]
  instances_to_consider = instances
  #instances_to_consider = instances['X']
  for instance in instances_to_consider:
    if 'X' in instance:
      time = list(range(0,50000,100))
    else:
      time = list(range(0,3600,100))
    lbs = []
    labels = []
    instance_results = time_results_df[time_results_df['instance'] == instance]
    for test in tests_to_compare:
      instance_test_results = instance_results[instance_results['test'] == test]
      if not instance_test_results.empty:
        instance_test_lb = []
        for t in time:
          instance_test_results_time = instance_test_results[instance_test_results['time'].astype(int) <= t]
          if (instance_test_results_time.empty):
            lb = 0
          else:
            lb = max(instance_test_results_time['lb'])
          instance_test_lb.append(lb)
        lbs.append(instance_test_lb)
        labels.append(test)

    min_val = 10000000
    max_val = 0
    for (label, lb) in zip(labels, lbs):
      if len([l for l in lb if l >0]) == 0:
        continue
      to_plot = zip(time,lb)
      to_plot = [data for data in to_plot if data[1] > 0]
      times_to_plot = [d[0] for d in to_plot]
      lbs_to_plot = [d[1] for d in to_plot]
      plt.plot(times_to_plot,lbs_to_plot,label=label)
      min_value = min([l for l in lb if l > 0])
      max_value = max(lb)
      min_val = min(min_val, min_value)
      max_val = max(max_val, max_value)
    plt.title(instance)
    plt.xlabel('time (s)')
    plt.ylabel('lb')
    if instance in instance_upper_bounds:
      plt.axhline(y=instance_upper_bounds[instance],linewidth=2,label="opt",color="magenta")
      plt.ylim(min_val-10,instance_upper_bounds[instance]+10)
    else:
      plt.ylim(min_val,max_val)
    if instance in pecin_root_lbs:
      if 'X' in instance:
        plt.scatter([pecin_root_lbs[instance][1]*60], [pecin_root_lbs[instance][0]],label="pecin_root_lb")
      else:
        plt.scatter([pecin_root_lbs[instance][1]], [pecin_root_lbs[instance][0]],label="pecin_root_lb")
    plt.legend()
    plt.show()

# experiment 3 - lp runtimes
experiment3 = False
if experiment3:
  tests_to_include = ["test_A_colelim_lp_ng_2_10_1800", "test_A_colelim_lp_ng_5_10_1800","test_M_colgen_dd_q_1_1800"]
  instances_to_consider = ["A-n32-k5.vrp", "A-n33-k5.vrp","A-n36-k5.vrp","A-n37-k6.vrp","A-n38-k5.vrp","A-n80-k10.vrp","M-n101-k10.vrp","M-n121-k7.vrp","M-n151-k12.vrp"]
  #instances_to_consider = ["A-n80-k10.vrp","M-n101-k10.vrp","M-n121-k7.vrp","M-n151-k12.vrp"]
  ddSizes = []
  lpSolveTimes = []
  labels = []
  for instance in instances_to_consider:
    instance_results = time_results_df[time_results_df['instance'] == instance]
    for test in tests_to_include:
      instance_test_results = instance_results[instance_results['test'] == test]
      if not instance_test_results.empty:
        instance_test_lb = []
        curr_time = 0
        instance_test_results.dropna(inplace=True)
        for index, row in instance_test_results.iterrows():
          solve_time = int(row['time']) - curr_time
          curr_time = int(row['time'])
          dd_size = int(row['size'])
          ddSizes.append(dd_size)
          lpSolveTimes.append(solve_time)
          if 'M' in instance:
            labels.append('blue')
          elif 'A' in instance:
            labels.append('orange')

  plt.scatter(ddSizes, lpSolveTimes, color=labels)
  plt.title('LP time vs. DD size')
  plt.xlabel('DD size')
  plt.ylabel('time (s)')
  plt.show()


'''                   
                    'CELAGmuSSP' : ['test_A_colelim_lag_q_1_2_N_3600',
                                    'test_B_colelim_lag_q_1_2_N_3600',
                                    'test_M_colelim_lag_q_1_2_N_3600'],
'''
# experiment Baseline
# We want this to be a performance plot
# First showing how many instances are solved within d% at times 0-3600
# Then showing how many instances are solved within x% of the optimal value
experiment4 = False
if experiment4:
  methodsToTests = {'CE_MIP' : ['all_col_elim_lp_ub_3600'],
                    'CE_BANDB' : ['all_col_elim_lp_bandb_3600'],
                    'Held' : ['held_3600']}
  # d% gap
  # List of times
  # List of percent gaps
  # x_axis names?
  d = 5
  times_list = list(range(0,3700,100))
  gap_list = list(range(d+1,21,1))
  line_styles = ['solid', 'dotted', 'dashed', 'dashdot']
  line_colors = ['b', 'r', 'g', 'k']

  # For each method build:
  # - num solved within d% at each time
  # - num solved within each percent gaps
  method_num = 0
  for method in methodsToTests:
    count = 0
    method_time_list = [0] * len(times_list)
    method_gap_list = [0] * len(gap_list)
    for test in methodsToTests[method]:
      for instance in instances:
        count = count + 1
        optimal = instance_upper_bounds[instance]
        test_instance_results = time_results_df[(time_results_df['instance'] == instance) & (time_results_df['test'] == test)]
        if not test_instance_results.empty:
          print(method)
          print(test)
          # time info
          for t_iter in range(len(times_list)):
            t = times_list[t_iter]
            time_result = test_instance_results[test_instance_results['time'].astype(int) <= t]
            if time_result.empty:
              lb = 0
            else:
              lb = max(time_result['lb'])
            if lb >= (optimal * (1-(d*1.0/100))):
              method_time_list[t_iter] = method_time_list[t_iter] + 1
          # gap info
          best_lb = max(test_instance_results['lb'])
          for gap_iter in range(len(gap_list)):
            gap = gap_list[gap_iter]
            if best_lb >= (optimal * (1- (0.01*gap))):
              method_gap_list[gap_iter] = method_gap_list[gap_iter] + 1

    # plot each of these time then reverse gaps?
    x_plot = times_list.copy()
    for gap_iter in range(len(gap_list)):
        x_plot.append(3600 + (3600 / len(gap_list) * (gap_iter+1)))
    y_plot = method_time_list + method_gap_list
    plt.plot(x_plot,y_plot,label=method,linestyle=line_styles[method_num],color=line_colors[method_num])
    method_num = method_num + 1
    print("method count")
    print(count)

  # plot
  plt.xlabel('time (s) | optimality gap (%)')
  plt.ylabel('# instances')

  # update x axis
  x_ticks = [100,1000,1800,3600]
  x_labels = [100,1000,1800,3600]
  #x_ticks_raw = times_list.copy()
  #x_labels_raw = times_list.copy()

  for gap_iter in range(len(gap_list)):
    if gap_iter % 2 == 0 and gap_iter > 0:
      x_ticks.append(3600 + (3600 / len(gap_list) * (gap_iter+1)))
      x_labels.append(str(gap_list[gap_iter]) + '%')

  plt.xticks(x_ticks, x_labels)
  plt.legend()
  plt.axvline(x=3600,color='k')
  plt.show()

# Get the average optimality gap for each class of instances
experiment5 = False
if experiment5:
  instance_type_tests = {'A' : ['test_A_colelim_lp_ng_2_20_Y_3600', 'test_A_colelim_lag_ng_2_20_Y_3600'],
                         'B' : ['test_B_colelim_lp_ng_2_20_Y_3600', 'test_B_colelim_lag_ng_2_20_Y_3600'],
                         'E' : ['test_E_colelim_lp_ng_2_20_Y_3600', 'test_E_colelim_lag_ng_2_20_Y_3600'],
                         'F' : ['test_F_colelim_lp_ng_2_20_Y_3600', 'test_F_colelim_lag_ng_2_20_Y_3600'],
                         'M' : ['test_M_colelim_lp_ng_2_20_Y_3600', 'test_M_colelim_lag_ng_2_20_Y_3600'],
                         'P' : ['test_P_colelim_lp_ng_2_20_Y_3600', 'test_P_colelim_lag_ng_2_20_Y_3600'],
                         'X' : ['test_X_colelim_lp_ng_2_20_Y_7200', 'test_X_colelim_lag_ng_2_20_Y_7200']}

  best_gaps = {}
  instance_set_list = ['A','B','E','F','M','P','X']
  for instance_type in instance_set_list:
    instance_type_gaps = []
    for instance in instances[instance_type]:
      # Remove instances not in Pecin
      q2_non_compile = ['A-n32','A-n33','A-n34','A-n36','B-n31','B-n34','B-n35']
      shouldSkip = False
      for q2_no in q2_non_compile:
        if q2_no in instance:
          shouldSkip = True
          break
      if shouldSkip:
        continue

      optimal = instance_optimal_values[instance]
      instance_best_gap = 10000000000
      for test in instance_type_tests[instance_type]:
        test_instance_results = time_results_df[(time_results_df['instance'] == instance) & (time_results_df['test'] == test)]
        if not test_instance_results.empty:
          lb = max(test_instance_results['lb'])
        else:
          lb = 0

        gap = (optimal - lb) * 100.0 / optimal
        instance_best_gap = min(instance_best_gap, gap)
      if instance_best_gap != 100.0:
        instance_type_gaps.append(instance_best_gap)
      else:
        print(instance)
    print(instance_type_gaps)
    instance_type_key = instance_type
    if instance_type == 'M':
      instance_type_key = 'E'
    best_gaps[instance_type_key] = 1.0 * sum(instance_type_gaps) / len(instance_type_gaps)
  print(best_gaps)

# experiment Baseline
# We want this to be a performance plot
# Don't use best ub like experiment4 though
# First showing how many instances are solved within d% at times 0-3600
# Then showing how many instances are solved within x% of the optimal value
experiment6 = False
if experiment6:
  methodsToTests = {'CE_MIP' : ['all_col_elim_lp_ub_3600'],
                    'CE_BANDB' : ['all_col_elim_lp_bandb_ub_3600'],
                    'CE_BDD' : ['all_col_elim_bdd_3600'],
                    'Held' : ['held_3600']}
  # List of times
  # List of percent gaps
  # x_axis names?
  d = 0
  times_list = list(range(0,3700,100))
  gap_list = list(range(d,21,1))
  line_styles = ['solid', 'dotted', 'dashed', 'dashdot']
  line_colors = ['b', 'r', 'g', 'k']

  # For each method build:
  # - num solved within d% at each time
  # - num solved within each percent gaps
  method_num = 0
  for method in methodsToTests:
    count = 0
    method_time_list = [0] * len(times_list)
    method_gap_list = [0] * len(gap_list)
    for test in methodsToTests[method]:
      for instance in instances:
        count = count + 1
        optimal = instance_upper_bounds[instance]
        test_instance_results = time_results_df[(time_results_df['instance'] == instance) & (time_results_df['test'] == test)]
        if not test_instance_results.empty:
          print(method)
          print(test)
          # time info
          for t_iter in range(len(times_list)):
            t = times_list[t_iter]
            time_result = test_instance_results[test_instance_results['time'].astype(int) <= t]
            if time_result.empty:
              lb = 0
              ub = 1000000
            else:
              lb = max(time_result['lb'])
              ub = min(time_result['ub'])
            if lb >= ub * (1-(d*1.0/100)):
              method_time_list[t_iter] = method_time_list[t_iter] + 1
          # gap info
          best_lb = max(test_instance_results['lb'])
          best_ub = min(test_instance_results['ub'])
          for gap_iter in range(len(gap_list)):
            gap = gap_list[gap_iter]
            if best_lb >= (best_ub * (1- (0.01*gap))):
              method_gap_list[gap_iter] = method_gap_list[gap_iter] + 1

    # plot each of these time then reverse gaps?
    x_plot = times_list.copy()
    for gap_iter in range(len(gap_list)):
        x_plot.append(3600 + (3600 / len(gap_list) * (gap_iter+1)))
    y_plot = method_time_list + method_gap_list
    plt.plot(x_plot,y_plot,label=method,linestyle=line_styles[method_num],color=line_colors[method_num])
    method_num = method_num + 1
    print("method count")
    print(count)

  # plot
  plt.xlabel('time (s) | optimality gap (%)')
  plt.ylabel('# instances')

  # update x axis
  x_ticks = [100,1000,1800,3600]
  x_labels = [100,1000,1800,3600]
  #x_ticks_raw = times_list.copy()
  #x_labels_raw = times_list.copy()

  for gap_iter in range(len(gap_list)):
    if gap_iter % 2 == 0 and gap_iter > 0:
      x_ticks.append(3600 + (3600 / len(gap_list) * (gap_iter+1)))
      x_labels.append(str(gap_list[gap_iter]) + '%')

  plt.xticks(x_ticks, x_labels)
  plt.legend()
  plt.axvline(x=3600,color='k')
  plt.show()
