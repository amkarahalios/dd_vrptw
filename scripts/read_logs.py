import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections
import static
 
def get_instances(instance_set_name, root_directory):
  instance_dir = root_directory + "/instances/"

  instances = []
  if instance_set_name == "Solomon":
    instance_dir = instance_dir + "/Vrp-Set-Solomon/"
  elif instance_set_name == "HG":
    instance_dir = instance_dir + "/Vrp-Set-HG/"
  elif instance_set_name == "CVRP":
    instance_dir = instance_dir + "/CVRP/"
  elif instance_set_name == "X":
    instance_dir = instance_dir + "/X/"
  elif instance_set_name == "PDPTW":
    instance_dir = instance_dir + "/pdp_200/"

  for instance in os.listdir(instance_dir):
    if instance in static.instance_upper_bounds:
      instances.append(instance)

  instances.sort()
  return instances

def get_column_elimination_results(instances, parameter_set_names, root_directory):
  regex_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] numArcs: \[([0-9]+)\] numFixed: \[([0-9]+)\] numHeuristicIPs: \[([0-9]+)\] numHeuristicLNS: \[([0-9]+)\] numPrimalLNSRepairs: \[([0-9]+)\] time: \[([0-9]+)\]")

  results = []
  for parameter_name in parameter_set_names:
    for instance in instances:
      log_file_name = root_directory + "/new_logs/" + parameter_name + '/' + instance + '.log'
      if not os.path.exists(log_file_name):
        continue

      log_file = open(log_file_name, "r")
      for line in log_file:
        match = regex_pattern.match(line)
        if match:
          lpIterations = match.group(1)
          lagIterations = match.group(2)
          sspIterations = match.group(3)
          numSeparations = match.group(4)
          compileTime = match.group(5)
          sspSolveTime = match.group(6)
          lpSolveTime = match.group(7)
          lb = float(match.group(8))
          ub = float(match.group(9))
          numArcs = int(match.group(10))
          numFixed = int(match.group(11))
          numIPs = int(match.group(12))
          numLNS = int(match.group(13))
          numRepairs = int(match.group(14))
          time = float(match.group(15))
          time_result = {'instance': instance, 'parameter': parameter_name, 'lpIt' : lpIterations, 'lagIt' : lagIterations, 'numSep' : numSeparations, 'lb' : lb, 'ub' : ub, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
          results.append(time_result)

  results_df = pandas.DataFrame(results)
  results_df.sort_values(by=['instance','lb','ub'],inplace=True)
  results_df.reset_index(drop=True,inplace=True)

  print_table = True
  if print_table:
    table_results_df = results_df[['instance','parameter','lb','ub','lpIt','lagIt','time']]
    print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))
 
  # In case it went over time
  results_df = results_df[results_df['time'] <= 3600]

  return results_df

def get_vrpsolver_results(instances, vrpsolver_file_names, instance_set_name, root_directory):
  regex_instance = re.compile('  instance  =>  "/CVRP/data/X/(.*)"')
  regex_stats = re.compile('statistics: ([A-Z]+.*[0-9]) & [0-9] & [0-9]+.* & [0-9]+.* & [0-9]+.* & ([0-9]+) & ([0-9]+.*) & ([0-9\--]+.*) & ([0-9]+.*) \\.*')
  regex_pb = re.compile('<DWph.*et=([0-9]+).*PB=([0-9]+).*')

  results = []
  for file_name in vrpsolver_file_names:
    log_file_name = root_directory + "/new_logs/VRPSolver/" + file_name
    if not os.path.exists(log_file_name):
      continue

    log_file = open(log_file_name, "r")
    for line in log_file:
      file_name_match = regex_instance.match(line)
      if file_name_match:
        instance = file_name_match.group(1)

      stats_match = regex_stats.match(line)
      if stats_match:
        # update some names
        instance = stats_match.group(1)
        if instance_set_name == 'Solomon':
          instance = instance + '.vrptw'
        else:
          instance = instance + '.vrp'
        instance = instance.replace('210','2_10')
        instance = instance.replace('410','4_10')
        num_nodes = int(stats_match.group(2))
        lb = float(stats_match.group(3))
        ub = stats_match.group(4)
        if ub != "--":
          ub = float(stats_match.group(4))
        time = float(stats_match.group(5))
        if instance in instances:
          result = {'parameter': 'VRPSolver', 'instance' : instance, 'lb' : lb, 'ub' : ub, 'time' : time} 
          results.append(result)

      pb_match = regex_pb.match(line)
      if pb_match:
        time = float(pb_match.group(1))
        pb = float(pb_match.group(2))
        if instance in instances:
          result = {'parameter': 'VRPSolver', 'instance' : instance, 'lb' : 0, 'ub' : pb, 'time' : time} 
          results.append(result)

  results_df = pandas.DataFrame(results)
  if results_df.empty:
    return results_df
  else:
    results_df.sort_values(by=['instance','lb','ub'],inplace=True)
    results_df.reset_index(drop=True,inplace=True)

  # In case it went over time
  results_df = results_df[results_df['time'] <= 3600]
 
  print_table = True
  if print_table:
    table_results_df = results_df[['instance','parameter','lb','ub','time']]
    print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))
 
  vrpsolver_instances = set(results_df['instance'])
  for instance in instances:
    if not instance in vrpsolver_instances:
      result = {'parameter': 'VRPSolver', 'instance' : instance, 'lb' : 0, 'ub' : 1e9, 'time' : time} 
      results.append(result)
      print(f"PyVRP missing instance {instance}")

  return results_df

def get_pyvrp_results(instances, pyvrp_file_names, instance_set_name, root_directory):
  regex_file_name = re.compile("Solving the file: (.*)")
  regex_stats = re.compile("[H \|]*[0-9]+[ \|]*([0-9].*)s[ \|]*[0-9]+[ \|]*[0-9]+[ \|]*([0-9]+)[ \|]*[0-9]+[ \|]*[0-9]+[ \|]*[0-9]+[ \|]")

  results = []
  for file_name in pyvrp_file_names:
    log_file_name = root_directory + "/new_logs/PyVRP/" + file_name
    if not os.path.exists(log_file_name):
      continue

    log_file = open(log_file_name, "r")
    for line in log_file:
      file_match = regex_file_name.match(line)
      if file_match:
        instance = file_match.group(1)

      stats_match = regex_stats.match(line)
      if stats_match:
        lb = 1
        time = float(stats_match.group(1))
        ub = float(stats_match.group(2))
        if instance in instances:
          result = {'parameter': 'PyVRP', 'instance' : instance, 'lb' : lb, 'ub' : ub, 'time' : time} 
          results.append(result)

  results_df = pandas.DataFrame(results)
  if results_df.empty:
    return results_df
  else:
    results_df.sort_values(by=['instance','lb','ub'],inplace=True)
    results_df.reset_index(drop=True,inplace=True)

  # In case it went over time
  results_df = results_df[results_df['time'] <= 3600]
 
  print_table = True
  if print_table:
    table_results_df = results_df[['instance','parameter','lb','ub','time']]
    print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))

  pyvrp_instances = set(results_df['instance'])
  for instance in instances:
    if not instance in pyvrp_instances:
      result = {'parameter': 'PyVRP', 'instance' : instance, 'lb' : 0, 'ub' : 1e9, 'time' : time} 
      results.append(result)
      print(f"PyVRP missing instance {instance}")

  return results_df
