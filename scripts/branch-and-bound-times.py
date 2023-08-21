# Plot the difference in solve time using variable ordering or not
import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections

# Log lines for DDSolver
#STATS - lpIterations[1] lagIterations[32] sspIterations[186] numSeparations[0] compileTime[240] sspSolveTime[349] lpSolveTime[0] lb[794.025] ub[1e+10] numArcs: [87254] numFixed: [0] time: [590]

colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] numArcs: \[([0-9]+)\] numFixed: \[([0-9]+)\] time: \[([0-9]+)\]")
#colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] size: \[([0-9]+)\] time: \[([0-9]+)\]")

logs_dir = "/Users/akarahal/Desktop/dd_graph_color/logs/"
test_set = ["col_elim_gc_lp_0_0_MIP_3600", "col_elim_gc_lp_0_0_BANDB_3600"]

instances = []
instance_dir = "/Users/akarahal/Desktop/dd_graph_color/instances/GC/"

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
        lb = int(colelim_match.group(8))
        ub = int(colelim_match.group(9))
        numArcs = int(colelim_match.group(10))
        numFixed = int(colelim_match.group(11))
        time = float(colelim_match.group(12))
        time_result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'numSep' : numSeparations, 'lb' : lb, 'ub' : ub, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
        time_results.append(time_result)
 
    if col_elim:
      result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'numSep' : numSeparations, 'lb' : lb, 'ub': ub, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
      results.append(result)
    else:
      result = {'instance': instance, 'test': test, 'iterations' : numpy.nan, 'numSep' : numpy.nan, 'lb' : numpy.nan, 'ub': ub, 'numArcs': numpy.nan, 'numFixed': numpy.nan, 'time' : numpy.nan}
      results.append(result)

# give table of results
results_df = pandas.DataFrame(results)
results_df.sort_values(by=['instance','lb'],inplace=True)
results_df.reset_index(drop=True,inplace=True)
table_results_df = results_df[['instance','test','lb','ub','time']]
print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))

# setup time_results df
time_results_df = pandas.DataFrame(time_results)
time_results_df.sort_values(by=['instance','lb'],inplace=True)
time_results_df.reset_index(drop=True,inplace=True)

method_pairs = [["col_elim_gc_lp_0_0_BANDB_3600", "col_elim_gc_lp_0_0_MIP_3600"]]

# gather times to solve
for method_pair in method_pairs:
  bandb_times = []
  no_bandb_times = []
  for instance in instances:
    bandb_method = method_pair[0]
    bandb_test_instance_results = time_results_df[(time_results_df['instance'] == instance) & (time_results_df['test'] == bandb_method)]
    bandb_solve = True
    if bandb_test_instance_results.empty:
      bandb_solve = False
    bandb_best_lb = max(bandb_test_instance_results['lb'])
    bandb_best_ub = min(bandb_test_instance_results['ub'])
    if bandb_best_lb != bandb_best_ub:
      bandb_solve = False

    no_bandb_method = method_pair[1]
    no_bandb_test_instance_results = time_results_df[(time_results_df['instance'] == instance) & (time_results_df['test'] == no_bandb_method)]
    not_bandb_solve = True
    if no_bandb_test_instance_results.empty:
      not_bandb_solve = False
    no_bandb_best_lb = max(no_bandb_test_instance_results['lb'])
    no_bandb_best_ub = min(no_bandb_test_instance_results['ub'])
    if no_bandb_best_lb != no_bandb_best_ub:
      not_bandb_solve = False

    if bandb_solve and not_bandb_solve:
      bandb_time = max(bandb_test_instance_results['time'])
      bandb_times.append(bandb_time)

      no_bandb_time = max(no_bandb_test_instance_results['time'])
      no_bandb_times.append(no_bandb_time)
    elif bandb_solve:
      print(f'bandb solves: {instance}')
    elif not_bandb_solve:
      print(f'not-bandb solves: {instance}')

  diffs = [bandb - no_bandb for (bandb,no_bandb) in zip(bandb_times,no_bandb_times)]
  avg_diff = sum(diffs) / len(diffs)
  print(f'bandb avg diff: {avg_diff}')

  mults = [no_bandb / bandb for (bandb,no_bandb) in zip(bandb_times,no_bandb_times) if bandb != 0]
  avg_mult = sum(mults) / len(mults)
  print(f'bandb avg mult: {avg_mult}')
  plt.scatter(bandb_times, no_bandb_times, marker='x', color='b')

plt.xlabel('time to solve with branch-and-bound (s)')
plt.xscale('log')

plt.ylabel('time to solve without branch-and-bound (s)')
plt.yscale('log')

plt.plot([0,10000],[0,10000])
plt.show()
