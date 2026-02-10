# Create table to compare CE with SOA
import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections
import gzip

instance_upper_bounds = {
"LC1_2_1.pdp":2704.6,
"LC1_2_2.pdp":2764.6,
"LC1_2_3.pdp":2772.2,
"LC1_2_4.pdp":2661.4,
"LC1_2_5.pdp":2702.0,
"LC1_2_6.pdp":2701.0,
"LC1_2_7.pdp":2701.0,
"LC1_2_8.pdp":2689.8,
"LC1_2_9.pdp":2724.2,
"LC1_2_10.pdp":2741.6,
"LC2_2_1.pdp":1931.44,
"LC2_2_2.pdp":1881.40,
"LC2_2_3.pdp":1844.33,
"LC2_2_4.pdp":1767.12,
"LC2_2_5.pdp":1891.21,
"LC2_2_6.pdp":1857.78,
"LC2_2_7.pdp":1850.13,
"LC2_2_8.pdp":1824.34,
"LC2_2_9.pdp":1854.21,
"LC2_2_10.pdp":1817.45,
"LR1_2_1.pdp":4819.1,
"LR1_2_2.pdp":4093.1,
"LR1_2_3.pdp":3486.8,
"LR1_2_4.pdp":2830.7,
"LR1_2_5.pdp":4221.6,
"LR1_2_6.pdp":3763.0,
"LR1_2_7.pdp":3112.9,
"LR1_2_8.pdp":2645.5,
"LR1_2_9.pdp":3953.5,
"LR1_2_10.pdp":3386.3,
"LR2_2_1.pdp":4073.10,
"LR2_2_2.pdp":3796.00,
"LR2_2_3.pdp":3098.36,
"LR2_2_4.pdp":2486.00,
"LR2_2_5.pdp":3438.39,
"LR2_2_6.pdp":4457.95,
"LR2_2_7.pdp":3098.35,
"LR2_2_8.pdp":2449.36,
"LR2_2_9.pdp":3922.11,
"LR2_2_10.pdp":3254.83,
"LRC1_2_1.pdp":3606.1,
"LRC1_2_2.pdp":3292.4,
"LRC1_2_3.pdp":3079.5,
"LRC1_2_4.pdp":2525.8,
"LRC1_2_5.pdp":3715.8,
"LRC1_2_6.pdp":3360.9,
"LRC1_2_7.pdp":3317.7,
"LRC1_2_8.pdp":3086.5,
"LRC1_2_9.pdp":3053.8,
"LRC1_2_10.pdp":2837.5,
"LRC2_2_1.pdp":3595.18,
"LRC2_2_2.pdp":3158.25,
"LRC2_2_3.pdp":2881.99,
"LRC2_2_4.pdp":2835.40,
"LRC2_2_5.pdp":2776.93,
"LRC2_2_6.pdp":2707.96,
"LRC2_2_7.pdp":3010.68,
"LRC2_2_8.pdp":2399.89,
"LRC2_2_9.pdp":2208.49,
"LRC2_2_10.pdp":2437.88,
"LC1_4_1.pdp":7152.06,
"LC1_4_2.pdp":8007.79,
"LC1_4_3.pdp":8678.23,
"LC1_4_4.pdp":6451.68,
"LC1_4_5.pdp":7150.00,
"LC1_4_6.pdp":7154.02,
"LC1_4_7.pdp":7149.43,
"LC1_4_8.pdp":8305.42,
"LC1_4_9.pdp":7451.20,
"LC1_4_10.pdp":7850.22,
"LC2_4_1.pdp":4116.33,
"LC2_4_2.pdp":4144.29,
"LC2_4_3.pdp":4401.08,
"LC2_4_4.pdp":3743.95,
"LC2_4_5.pdp":4030.63,
"LC2_4_6.pdp":3900.29,
"LC2_4_7.pdp":3962.51,
"LC2_4_8.pdp":3844.45,
"LC2_4_9.pdp":4188.93,
"LC2_4_10.pdp":3828.44,
"LR1_4_1.pdp":10639.75,
"LR1_4_2.pdp":11009.51,
"LR1_4_3.pdp":9245.48,
"LR1_4_4.pdp":7007.90,
"LR1_4_5.pdp":11374.06,
"LR1_4_6.pdp":11334.11,
"LR1_4_7.pdp":8675.36,
"LR1_4_8.pdp":6164.82,
"LR1_4_9.pdp":9859.47,
"LR1_4_10.pdp":8192.65,
"LR2_4_1.pdp":9726.88,
"LR2_4_2.pdp":9405.40,
"LR2_4_3.pdp":10176.94,
"LR2_4_4.pdp":6201.84,
"LR2_4_5.pdp":9894.46,
"LR2_4_6.pdp":8946.91,
"LR2_4_7.pdp":7993.16,
"LR2_4_8.pdp":5260.42,
"LR2_4_9.pdp":7926.07,
"LR2_4_10.pdp":7596.62,
"LRC1_4_1.pdp":9124.52,
"LRC1_4_2.pdp":8346.06,
"LRC1_4_3.pdp":7805.16,
"LRC1_4_4.pdp":5803.31,
"LRC1_4_5.pdp":8847.40,
"LRC1_4_6.pdp":8394.47,
"LRC1_4_7.pdp":8037.87,
"LRC1_4_8.pdp":7930.15,
"LRC1_4_9.pdp":8004.24,
"LRC1_4_10.pdp":7064.36,
"LRC2_4_1.pdp":9738.95,
"LRC2_4_2.pdp":7159.95,
"LRC2_4_3.pdp":6426.47,
"LRC2_4_4.pdp":5238.77,
"LRC2_4_5.pdp":7309.54,
"LRC2_4_6.pdp":6337.08,
"LRC2_4_7.pdp":6292.23,
"LRC2_4_8.pdp":5767.72,
"LRC2_4_9.pdp":6272.92,
"LRC2_4_10.pdp":5395.71,
"LC1_6_1.pdp":14095.64,
"LC1_6_2.pdp":15048.16,
"LC1_6_3.pdp":14874.68,
"LC1_6_4.pdp":13300.55,
"LC1_6_5.pdp":14086.30,
"LC1_6_6.pdp":14090.79,
"LC1_6_7.pdp":14083.76,
"LC1_6_8.pdp":14880.70,
"LC1_6_9.pdp":15207.58,
"LC1_6_10.pdp":15432.55,
"LC2_6_1.pdp":7977.98,
"LC2_6_2.pdp":9900.48,
"LC2_6_3.pdp":8512.94,
"LC2_6_4.pdp":7860.38,
"LC2_6_5.pdp":9051.53,
"LC2_6_6.pdp":8775.55,
"LC2_6_7.pdp":9376.58,
"LC2_6_8.pdp":7579.63,
"LC2_6_9.pdp":8714.22,
"LC2_6_10.pdp":7946.60,
"LR1_6_1.pdp":22821.65,
"LR1_6_2.pdp":20137.22,
"LR1_6_3.pdp":17846.17,
"LR1_6_4.pdp":13127.56,
"LR1_6_5.pdp":23623.52,
"LR1_6_6.pdp":23084.50,
"LR1_6_7.pdp":16963.37,
"LR1_6_8.pdp":11917.70,
"LR1_6_9.pdp":21835.87,
"LR1_6_10.pdp":19298.25,
"LR2_6_1.pdp":21759.33,
"LR2_6_2.pdp":21289.70,
"LR2_6_3.pdp":17480.06,
"LR2_6_4.pdp":10639.08,
"LR2_6_5.pdp":22625.89,
"LR2_6_6.pdp":19099.49,
"LR2_6_7.pdp":14645.32,
"LR2_6_8.pdp":12341.90,
"LR2_6_9.pdp":18259.20,
"LR2_6_10.pdp":16390.64,
"LRC1_6_1.pdp":18288.90,
"LRC1_6_2.pdp":16515.41,
"LRC1_6_3.pdp":13975.98,
"LRC1_6_4.pdp":10800.44,
"LRC1_6_5.pdp":17463.94,
"LRC1_6_6.pdp":19025.36,
"LRC1_6_7.pdp":15914.85,
"LRC1_6_8.pdp":15317.63,
"LRC1_6_9.pdp":15344.64,
"LRC1_6_10.pdp":13963.66,
"LRC2_6_1.pdp":14578.92,
"LRC2_6_2.pdp":13850.22,
"LRC2_6_3.pdp":12432.91,
"LRC2_6_4.pdp":9833.92,
"LRC2_6_5.pdp":14380.37,
"LRC2_6_6.pdp":14962.08,
"LRC2_6_7.pdp":14094.26,
"LRC2_6_8.pdp":13179.58,
"LRC2_6_9.pdp":16309.78,
"LRC2_6_10.pdp":13054.50,
"LC1_8_1.pdp":25184.38,
"LC1_8_2.pdp":29735.85,
"LC1_8_3.pdp":27287.80,
"LC1_8_4.pdp":22686.08,
"LC1_8_5.pdp":25211.22,
"LC1_8_6.pdp":25164.25,
"LC1_8_7.pdp":25158.38,
"LC1_8_8.pdp":26688.17,
"LC1_8_9.pdp":26975.61,
"LC1_8_10.pdp":27099.59,
"LC2_8_1.pdp":11687.06,
"LC2_8_2.pdp":13713.13,
"LC2_8_3.pdp":12981.37,
"LC2_8_4.pdp":12917.34,
"LC2_8_5.pdp":12298.33,
"LC2_8_6.pdp":12645.71,
"LC2_8_7.pdp":14041.47,
"LC2_8_8.pdp":12924.10,
"LC2_8_9.pdp":11629.41,
"LC2_8_10.pdp":12226.42,
"LR1_8_1.pdp":39291.32,
"LR1_8_2.pdp":34074.49,
"LR1_8_3.pdp":29225.36,
"LR1_8_4.pdp":20675.70,
"LR1_8_5.pdp":39809.68,
"LR1_8_6.pdp":36244.90,
"LR1_8_7.pdp":27078.09,
"LR1_8_8.pdp":20705.45,
"LR1_8_9.pdp":38499.00,
"LR1_8_10.pdp":31091.07,
"LR2_8_1.pdp":41873.12,
"LR2_8_2.pdp":36549.68,
"LR2_8_3.pdp":26391.07,
"LR2_8_4.pdp":20618.03,
"LR2_8_5.pdp":34816.50,
"LR2_8_6.pdp":28814.65,
"LR2_8_7.pdp":25502.39,
"LR2_8_8.pdp":18327.23,
"LR2_8_9.pdp":30515.50,
"LR2_8_10.pdp":30136.86,
"LRC1_8_1.pdp":32252.28,
"LRC1_8_2.pdp":27878.89,
"LRC1_8_3.pdp":24371.95,
"LRC1_8_4.pdp":18208.51,
"LRC1_8_5.pdp":31169.16,
"LRC1_8_6.pdp":28961.66,
"LRC1_8_7.pdp":28768.40,
"LRC1_8_8.pdp":26902.93,
"LRC1_8_9.pdp":24854.96,
"LRC1_8_10.pdp":24622.59,
"LRC2_8_1.pdp":23074.22,
"LRC2_8_2.pdp":22220.95,
"LRC2_8_3.pdp":20379.26,
"LRC2_8_4.pdp":14711.80,
"LRC2_8_5.pdp":23602.44,
"LRC2_8_6.pdp":22591.26,
"LRC2_8_7.pdp":25436.79,
"LRC2_8_8.pdp":22604.36,
"LRC2_8_9.pdp":22594.69,
"LRC2_8_10.pdp":19604.12,
"LC1_10_1.pdp":42488.66,
"LC1_10_2.pdp":44548.51,
"LC1_10_3.pdp":45906.30,
"LC1_10_4.pdp":37782.17,
"LC1_10_5.pdp":42477.40,
"LC1_10_6.pdp":42838.39,
"LC1_10_7.pdp":42854.99,
"LC1_10_8.pdp":42949.56,
"LC1_10_9.pdp":43564.01,
"LC1_10_10.pdp":42929.15,
"LC2_10_1.pdp":16879.24,
"LC2_10_2.pdp":20764.09,
"LC2_10_3.pdp":19299.48,
"LC2_10_4.pdp":17886.97,
"LC2_10_5.pdp":17137.53,
"LC2_10_6.pdp":19387.75,
"LC2_10_7.pdp":18389.37,
"LC2_10_8.pdp":17015.03,
"LC2_10_9.pdp":18225.30,
"LC2_10_10.pdp":17043.64,
"LR1_10_1.pdp":56744.91,
"LR1_10_2.pdp":49349.81,
"LR1_10_3.pdp":41483.89,
"LR1_10_4.pdp":30792.87,
"LR1_10_5.pdp":59053.68,
"LR1_10_6.pdp":51219.61,
"LR1_10_7.pdp":38563.20,
"LR1_10_8.pdp":29613.27,
"LR1_10_9.pdp":54582.25,
"LR1_10_10.pdp":45510.71,
"LR2_10_1.pdp":62859.29,
"LR2_10_2.pdp":51357.30,
"LR2_10_3.pdp":42908.33,
"LR2_10_4.pdp":26595.39,
"LR2_10_5.pdp":55038.21,
"LR2_10_6.pdp":46253.37,
"LR2_10_7.pdp":39648.54,
"LR2_10_8.pdp":26744.49,
"LR2_10_9.pdp":50781.83,
"LR2_10_10.pdp":44888.99,
"LRC1_10_1.pdp":49111.78,
"LRC1_10_2.pdp":45547.38,
"LRC1_10_3.pdp":35616.85,
"LRC1_10_4.pdp":27211.87,
"LRC1_10_5.pdp":50323.04,
"LRC1_10_6.pdp":45115.22,
"LRC1_10_7.pdp":41560.52,
"LRC1_10_8.pdp":40770.10,
"LRC1_10_9.pdp":40934.27,
"LRC1_10_10.pdp":36539.03,
"LRC2_10_1.pdp":34463.46,
"LRC2_10_2.pdp":38619.13,
"LRC2_10_3.pdp":27218.08,
"LRC2_10_4.pdp":23212.46,
"LRC2_10_5.pdp":40638.13,
"LRC2_10_6.pdp":30910.65,
"LRC2_10_7.pdp":33275.24,
"LRC2_10_10.pdp":29085.28
}

baldacci_bounds = {
"LC1_2_1.pdp":(2704.6,2704.6,3.3),
"LC1_2_2.pdp":(2764.6,2764.6,21.5),
"LC1_2_3.pdp":(2772.2,2772.2,114.9),
"LC1_2_4.pdp":(2395.8,2661.4,454.2),
"LC1_2_5.pdp":(2702.0,2702.0,4.8),
"LC1_2_6.pdp":(2701.0,2701.0,7.4),
"LC1_2_7.pdp":(2701.0,2701.0,7.7),
"LC1_2_8.pdp":(2689.8,2689.8,16.0),
"LC1_2_9.pdp":(2724.2,2724.2,55.3),
"LC1_2_10.pdp":(2741.6,2741.6,137.1),
"LR1_2_1.pdp":(4819.1,4819.1,1.6),
"LR1_2_2.pdp":(4093.1,4093.1,20.6),
"LR1_2_3.pdp":(3486.8,3486.8,3690.8),
"LR1_2_4.pdp":(2341.8,2830.7,1809.6),
"LR1_2_5.pdp":(4221.6,4221.6,2.6),
"LR1_2_6.pdp":(3763.0,3763.0,180.9),
"LR1_2_7.pdp":(2761.8,3112.9,1320.4),
"LR1_2_8.pdp":(2150.8,2645.5,566.9),
"LR1_2_9.pdp":(3953.3,3953.5,15.4),
"LR1_2_10.pdp":(3386.3,3386.3,1376.7),
"LRC1_2_1.pdp":(3606.1,3606.1,3.1),
"LRC1_2_2.pdp":(3292.4,3292.4,322.3),
"LRC1_2_3.pdp":(2497.8,3079.5,304.3),
"LRC1_2_4.pdp":(1981.0,2525.8,188.2),
"LRC1_2_5.pdp":(3715.8,3715.8,42.1),
"LRC1_2_6.pdp":(3360.9,3360.9,7.0),
"LRC1_2_7.pdp":(3317.7,3317.7,408.2),
"LRC1_2_8.pdp":(3086.5,3086.5,1562.7),
"LRC1_2_9.pdp":(3053.8,3053.8,1757.2),
"LRC1_2_10.pdp":(2335.5,2837.5,217.4),
"LC1_10_1.pdp":(42488.7,42488.7,79.5),
"LC1_10_5.pdp":(42477.4,42477.4,118.7),
"LR1_10_1.pdp":(56744.9,56744.9,233.1),
"LR1_10_5.pdp":(52536.3,52901.3,4068.8),
"LRC1_10_1.pdp":(48398.8,48666.5,2533.3),
"LRC1_10_5.pdp":(38177.8,49287.1,1650.3)
}

# Log lines for DDSolver
#STATS - lpIterations[1] lagIterations[32] sspIterations[186] numSeparations[0] compileTime[240] sspSolveTime[349] lpSolveTime[0] lb[794.025] ub[1e+10] numArcs: [87254] numFixed: [0] time: [590]

colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] numArcs: \[([0-9]+)\] numFixed: \[([0-9]+)\].*time: \[([0-9]+)\]")
colelim_finish_pattern = re.compile("time elapsed: ([0-9]+)")

logs_dir = "/Users/akarahal/Desktop/dd_vrptw/new_new_logs/"
#test_set = ["pdp_200_LAG_NG2_upper_bound_counter", "pdp_400_LAG_NG2_upper_bound_counter", "pdp_600_LAG_NG2_upper_bound_counter", "pdptw800_LAG_NG2_upper_bound_counter", "pdptw1000_LAG_NG2_upper_bound_counter"]
test_set = ["pdp_200_LAG_NG2_3600_fix_0", "pdp_400_LAG_NG2_3600_fix_0", "pdp_600_LAG_NG2_3600_fix_0", "pdptw800_LAG_NG2_3600_fix_2", "pdptw1000_LAG_NG2_3600_fix_2"]

instances = []
base_instance_dir = "/Users/akarahal/Desktop/dd_vrptw/instances/"
dir_list = ["pdp_200","pdp_400","pdp_600",'pdptw800','pdptw1000']

for d in dir_list:
  instance_dir = base_instance_dir + d + '/'
  for instance in os.listdir(instance_dir):
    instances.append(instance)

ce_optimal_instances = []
time_results = []
results = []
for test in test_set:
  for instance in instances:
    log_file_name = logs_dir + test + '/' + instance + '.log.gz'
    if not os.path.exists(log_file_name):
      continue

    col_elim = False
    with gzip.open(log_file_name, "rt", encoding='utf-8') as log_file:
      for line in log_file:
        colelim_match = colelim_pattern.match(line)
        if colelim_match:
          col_elim = True
          lpIterations = int(colelim_match.group(1))
          lagIterations = int(colelim_match.group(2))
          sspIterations = colelim_match.group(3)
          numSeparations = int(colelim_match.group(4))
          compileTime = colelim_match.group(5)
          sspSolveTime = colelim_match.group(6)
          lpSolveTime = colelim_match.group(7)
          lb = float(colelim_match.group(8))
          ub = float(colelim_match.group(9))
          numArcs = int(colelim_match.group(10))
          numFixed = int(colelim_match.group(11))
          time = float(colelim_match.group(12))
          time_result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'lagIterations': lagIterations, 'numSep' : numSeparations, 'lb' : lb, 'ub' : ub, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
          time_results.append(time_result)
 
        finish_match = colelim_finish_pattern.match(line)
        if finish_match:
          finish_time = float(finish_match.group(1))
          if finish_time < 3599:
            ce_optimal_instances.append(instance)
   
      if col_elim:
        result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'lagIterations': lagIterations, 'numSep' : numSeparations, 'lb' : lb, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
        results.append(result)
      else:
        result = {'instance': instance, 'test': test, 'iterations' : numpy.nan, 'lagIterations': numpy.nan, 'numSep' : numpy.nan, 'lb' : numpy.nan, 'numArcs': numpy.nan, 'numFixed': numpy.nan, 'time' : numpy.nan}
        results.append(result)

# add VrpSolver to results and time_results
#statistics_cols: instance & :Optimal & cutoff & :bcRecRootDb & :bcTimeRootEval & :bcCountNodeProc & :bcRecBestDb & :bcRecBestInc & :bcTimeMain \\
#statistics: C1_4_10 & 0 & 6825.5 & 6820.19 & 3420.03 & 1 & 6820.19 & -- & 3600.35 \\
#statistics: C2_2_1 & 1 & 1922.2 & 1922.10 & 228.49 & 1 & 1922.10 & 1922.10 & 228.50 \\
vrpsolver_pattern = re.compile('statistics: ([A-Z]+.*[0-9]) & [0-9] & [0-9]+.* & [0-9]+.* & [0-9]+.* & ([0-9]+) & ([0-9]+.*) & ([0-9\--]+.*) & ([0-9]+.*) \\.*')

vrpsolver_results = {}
logs_dir = "/Users/akarahal/Desktop/dd_vrptw/logs/VrpSolver/"
log_names = ['lilim200-results', 'lilim400-results', 'lilim600-results', 'lilim800-results', 'lilim1000-results']
for test in test_set:
  for log_name in log_names:
    log_file_name = logs_dir + '/' + log_name + '.log'
    if not os.path.exists(log_file_name):
      continue

    log_file = open(log_file_name, "r")
    for line in log_file:
      match = vrpsolver_pattern.match(line)
      if match:
        # might need to update some names
        instance = match.group(1)
        instance = instance + '.vrptw'
        instance = instance.replace('210','2_10')
        instance = instance.replace('410','4_10')
        num_nodes = int(match.group(2))
        lb = float(match.group(3))
        ub = match.group(4)
        if ub != "--":
          ub = float(match.group(4))
        time = float(match.group(5))
        vrpsolver_results[instance] = (num_nodes, lb, ub, time)

# give table of results
results_df = pandas.DataFrame(results)
results_df.sort_values(by=['instance','lb'],inplace=True)
results_df.reset_index(drop=True,inplace=True)
table_results_df = results_df[['instance','test','lb','time','iterations','numSep','lagIterations','numArcs']]
print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))

# get instance attributes: number locations
instance_num_locations = {}
for instance in instances:
  for d in dir_list:
    instance_dir = base_instance_dir + d + '/'
    instance_file_name = instance_dir + '/' + instance
    if os.path.exists(instance_file_name):
      instance_file = open(instance_file_name, "r")
      for line in instance_file:
        split_line = line.split()
        if len(split_line) > 0:
          num_locations = split_line[0]
      instance_num_locations[instance] = num_locations

soa_no_bound_count = 0
soa_optimal_count = 0
soa_gaps = []

ce_no_bound_count = 0
ce_optimal_count = 0
ce_gaps = []

# create output table for SOA comparison
num_printed = 0
for i, row in table_results_df.iterrows():
  instance = row['instance']

  baldacci_compare = False
  if baldacci_compare:
    if instance not in baldacci_bounds.keys():
      continue

  instance_name = instance.split('.')[0]
  instance_name = instance_name.replace("_","\\_")
  lb_value = row['lb']
  if not math.isnan(lb_value):
    #lb_value = math.ceil(float(lb_value) * 10) / 10.0
    lb_value = float(lb_value)
  else:
    lb_value = 0

  numArcs = row['numArcs']
  if not math.isnan(numArcs):
    numArcs = int(numArcs)

  time = row['time']
  if not math.isnan(time):
    time = int(time)

  if (not instance in ce_optimal_instances) and (lb_value < instance_upper_bounds[instance]):
    time = 3600

  numLpIterations = row['iterations']
  if not math.isnan(numLpIterations):
    numLpIterations = int(numLpIterations)
  else:
    numLpIterations = '-'
 
  numLagIterations = row['lagIterations']
  if not math.isnan(numLagIterations):
    numLagIterations = int(numLagIterations)
  else:
    numLagIterations = '-'

  numSeparations = row['numSep']
  if not math.isnan(numSeparations):
    numSeparations = int(numSeparations)
  else:
    numSeparations = '-'

  #gap = round((instance_upper_bounds[instance] - lb_value) * 100.0 / instance_upper_bounds[instance], 1)
  if instance in vrpsolver_results.keys():
    vrpsolver_result = vrpsolver_results[instance]
    vrpsolver_num_nodes = vrpsolver_result[0]
    vrpsolver_lb = vrpsolver_result[1]
    #vrpsolver_ub = vrpsolver_result[2]
    vrpsolver_time = int(vrpsolver_result[3])
  else:
    vrpsolver_num_nodes = '-'
    vrpsolver_lb = '-'
    #vrpsolver_ub = vrpsolver_result[1]

  if (vrpsolver_lb == '-') or (vrpsolver_lb < instance_upper_bounds[instance]):
    vrpsolver_time = 3600

  if baldacci_compare:
    if instance in baldacci_bounds.keys():
      baldacci_result = baldacci_bounds[instance]
      baldacci_lb = baldacci_result[0]
      baldacci_ub = baldacci_result[1]
      baldacci_time = baldacci_result[2]
    else:
      continue
  elif instance in baldacci_bounds.keys():
    continue
   
  if (lb_value == 0) or (lb_value == '-'):
    lb_value = '-'
    numLpIterations = '-'
    numLagIterations = '-'
    numSeparations = '-'

  if baldacci_compare:
    if (lb_value == '-') and (baldacci_lb == '-') and (vrpsolver_lb == '-'):
      continue
    else:
      num_printed = num_printed + 1
  else:
    if (lb_value == '-'):
      continue
    else:
      num_printed = num_printed + 1

  if (num_printed == 25) or ((num_printed > 27) and ((num_printed - 25) % 27 == 0)):
    print("")

  if baldacci_compare:
    if (lb_value != '-') and (baldacci_lb != '-') and (lb_value > baldacci_lb):
      print(f"{instance_name} & {instance_upper_bounds[instance]} & & {baldacci_lb} & {baldacci_ub} & {baldacci_time} & & {vrpsolver_lb} & {vrpsolver_time} & & {lb_value} & {numLpIterations} & {numLagIterations} & {numSeparations} & {time} \\\\")
    else:
      print(f"{instance_name} & {instance_upper_bounds[instance]} & & {baldacci_lb} & {baldacci_ub} & {baldacci_time} & & {vrpsolver_lb} & {vrpsolver_time} & & {lb_value} & {numLpIterations} & {numLagIterations} & {numSeparations} & {time} \\\\")

    if baldacci_lb == '-':
      soa_no_bound_count = soa_no_bound_count + 1
    elif baldacci_lb == instance_upper_bounds[instance]:
      if baldacci_lb != '-':
        soa_optimal_count = soa_optimal_count + 1
    else:
      gap = 100.0 * (instance_upper_bounds[instance] - baldacci_lb) / instance_upper_bounds[instance]
      soa_gaps.append(gap)

    if lb_value == '-':
      ce_no_bound_count = ce_no_bound_count + 1
    #elif lb_value == instance_upper_bounds[instance]:
    elif (instance in ce_optimal_instances):
      ce_optimal_count = ce_optimal_count + 1
    else:
      gap = 100.0 * (instance_upper_bounds[instance] - lb_value) / instance_upper_bounds[instance]
      ce_gaps.append(gap)
  else:
    #if (lb_value != '-') and (lb_value == instance_upper_bounds[instance]):
    if (lb_value != '-') and (instance in ce_optimal_instances):
      print(f"{instance_name} & {instance_upper_bounds[instance]} & & \\textbf\u007b{lb_value}\u007d & {numLpIterations} & {numLagIterations} & {numSeparations} & {time} \\\\")
    else:
      print(f"{instance_name} & {instance_upper_bounds[instance]} & & {lb_value} & {numLpIterations} & {numLagIterations} & {numSeparations} & {time} \\\\")
 
    if lb_value == '-':
      ce_no_bound_count = ce_no_bound_count + 1
    #elif lb_value == instance_upper_bounds[instance]:
    elif (instance in ce_optimal_instances):
      ce_optimal_count = ce_optimal_count + 1
    else:
      gap = 100.0 * (instance_upper_bounds[instance] - lb_value) / instance_upper_bounds[instance]
      ce_gaps.append(gap)

#print("soa optimal:")
#print(soa_optimal_count)
#print("soa no bound:")
#print(soa_no_bound_count)
#print("soa gap:")
#print(sum(soa_gaps) / len(soa_gaps))
print("ce optimal:")
print(ce_optimal_count)
print("ce no bound:")
print(ce_no_bound_count)
print("ce gap:")
print(sum(ce_gaps) / len(ce_gaps))
