import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections
import static

def calculate_gaps(ce_results_df, times_for_metrics):
  instances = set(ce_results_df['instance'])
  ce_parameter_set = set(ce_results_df['parameter'])

  gap_results = []
  for instance in instances:
    best_known = static.instance_upper_bounds[instance]
    if ".pdp" in instance:
      best_known = best_known + 10000 * static.instance_num_vehicles[instance]

    ce_instance_results = ce_results_df[ce_results_df['instance'] == instance]
    if not ce_instance_results.empty:
      for parameter in ce_parameter_set:
        parameter_primal_gaps = {}
        parameter_rel_gaps = {}
        instance_parameter_results = ce_instance_results[ce_instance_results['parameter'] == parameter]
        for end_time in times_for_metrics:
          if not instance_parameter_results.empty:
            instance_parameter_time_results = instance_parameter_results[instance_parameter_results['time'] <= end_time]
            if not instance_parameter_time_results.empty:
              instance_parameter_ubs = instance_parameter_time_results['ub']
              ub = min(instance_parameter_ubs)
              primal_gap = min(1,(ub - best_known) / max(best_known,ub))

              instance_parameter_lbs = instance_parameter_time_results['lb']
              lb = max(instance_parameter_lbs)
              rel_gap = (ub - lb) / ub

              gap_info = {'instance': instance, 'parameter' : parameter, 'rel_gap' : rel_gap, 'primal_gap' : primal_gap, 'time' : end_time} 
            else:
              gap_info = {'instance': instance, 'parameter' : parameter, 'rel_gap' : 1, 'primal_gap' : 1, 'time' : end_time} 
          else:
            gap_info = {'instance': instance, 'parameter' : parameter, 'rel_gap' : 1, 'primal_gap' : 1, 'time' : end_time} 
          gap_results.append(gap_info)

  gap_results_df = pandas.DataFrame(gap_results)
  return gap_results_df

def calculate_average_gaps(gap_results_df, times_for_metrics):
  average_results = []
  parameter_set = set(gap_results_df['parameter'])
  for parameter in parameter_set:
    parameter_df = gap_results_df[gap_results_df['parameter'] == parameter]
    for time in times_for_metrics:
      parameter_time_df = parameter_df[parameter_df['time'] == time]
      average_primal_gap = sum(parameter_time_df['primal_gap']) / len(parameter_time_df['primal_gap'])
      average_rel_gap = sum(parameter_time_df['rel_gap']) / len(parameter_time_df['rel_gap'])
      average_info = {'parameter' : parameter, 'time' : time, 'rel_gap' : average_rel_gap, 'primal_gap' : average_primal_gap}
      average_results.append(average_info)

  average_results_df = pandas.DataFrame(average_results)
  return average_results_df

def plot_average_gaps(instance_set_name, parameter_set_names_strings, average_gap_results_df):
  #parameter_set = set(average_gap_results_df['parameter'])
  parameter_set = [x[0] for x in parameter_set_names_strings.items()]

  markers = ['.','1','x','s','v','p','d','_']
  i = 0
  for parameter in parameter_set:
    parameter_df = average_gap_results_df[average_gap_results_df['parameter'] == parameter]
    plt.plot(list(parameter_df['time']), list(parameter_df['primal_gap']), label=parameter_set_names_strings[parameter], markersize=4, marker=markers[i])
    i = i + 1

  if instance_set_name == "Solomon":
    plt.title("Average p(t) for VRPTW Solomon")
  elif instance_set_name == "HG":
    plt.title("Average p(t) for VRPTW HG")
  elif instance_set_name == "CVRP":
    plt.title("Average p(t) for CVRP")
  elif instance_set_name == "X":
    plt.title("Average p(t) for CVRP X")
  elif instance_set_name == "PDPTW":
    plt.title("Average p(t) for PDPTW LiLim")

  plt.ylabel('p(t)')
  plt.xlabel('time (s)')
  plt.legend()
  plt.show()

def plot_gaps_by_instance(parameter_set_names_strings, gap_results):
  #parameter_set = set(gap_results['parameter'])
  parameter_set = [x[0] for x in parameter_set_names_strings.items()]
  instances = list(set(gap_results['instance']))
  instances.sort()

  markers = ['.','1','x','s','v','p','d','_']
  for instance in instances:
    i = 0
    instance_df = gap_results[gap_results['instance'] == instance]
    for parameter in parameter_set:
      parameter_df = instance_df[instance_df['parameter'] == parameter]
      plt.plot(list(parameter_df['time']), list(parameter_df['rel_gap']), label=parameter_set_names_strings[parameter], marker=markers[i])
      i = i + 1

    plt.ylim(0,0.1)

    plt.title(f"UB and LB for {instance}")
    plt.ylabel('objective value')
    plt.xlabel('time (s)')
    plt.legend()
    plt.show()

def print_instance_table(parameter_set_names_strings, gap_results_df):
  instances = list(set(gap_results_df['instance']))
  instances.sort()
  parameters = [parameter[0] for parameter in parameter_set_names_strings.items()]
  for instance in instances:
    best_known = static.instance_upper_bounds[instance]
    instance_results_df = gap_results_df[gap_results_df['instance'] == instance]
    result_string = f"{instance}"
    for parameter in parameters:
      instance_parameter_results_df = instance_results_df[instance_results_df['parameter'] == parameter]
      if instance_parameter_results_df.empty:
        primal_gap = 1
        rel_gap = 1
      else:
        primal_gap = round(min(instance_parameter_results_df['primal_gap']), 3)
        rel_gap = round(min(instance_parameter_results_df['rel_gap']), 3)
      #result_string = result_string + f" & & {primal_gap} & {rel_gap}"
      result_string = result_string + f" & {rel_gap}"

    result_string = result_string + "\\\\"
    print(result_string)

  header_string = "instance & & "
  for parameter in parameters:
    header_string = header_string + f"{parameter_set_names_strings[parameter]} & & "

  print(header_string)

def print_aggregated_table(instance_set_name, parameter_set_names_strings, gap_results_df):
  instances = list(set(gap_results_df['instance']))
  instances.sort()
  parameters = [parameter[0] for parameter in parameter_set_names_strings.items()]

  vrptw_aggregations = {
    "C1" : "^C1_",
    "C2" : "^C2_",
    "R" : "^R[12]_",
    "RC" : "^RC[12]_"
  }
  x_aggregations = {
    "X100" : "X-n1..-",
    "X200" : "X-n2..-",
    "X300-500" : "X-n[345]",
    "X600-1000" : "X-n[6789]"
  }
  pdp_aggregations = {
    "LC1 200" : "LC1_2",
    "LR1 200" : "LR1_2",
    "LC1 400" : "LC1_4",
    "LR1 400" : "LR1_4",
  }
  cvrp_aggregations = {
    "A" : "A-",
    "B" : "B-",
    "E" : "E-",
    "F" : "F-",
    "M" : "M-",
    "P" : "P-"
  }

  if instance_set_name == "Solomon":
    aggregations = {"Solomon" : ""}
  elif instance_set_name == "HG" or instance_set_name == "HGC1C2":
    aggregations = vrptw_aggregations
  elif instance_set_name == "CVRP":
    aggregations = cvrp_aggregations
  elif instance_set_name == "X":
    aggregations = x_aggregations
  elif instance_set_name == "PDPTW":
    aggregations = pdp_aggregations

  for aggregation_name, aggregation_regex in aggregations.items():
    aggregation_gap_results_df = gap_results_df[gap_results_df['instance'].str.contains(aggregation_regex, regex=True)]
    aggregation_instances = set(aggregation_gap_results_df['instance'])

    averages_string = f"{aggregation_name}"
    for parameter in parameters:
      parameter_results_df = aggregation_gap_results_df[aggregation_gap_results_df['parameter'] == parameter]
      if parameter_results_df.empty:
        continue
      parameter_primal_gaps = []
      parameter_rel_gaps = []
      for instance in aggregation_instances:
        parameter_instance_results_df = parameter_results_df[parameter_results_df['instance'] == instance]
        if parameter_instance_results_df.empty:
          primal_gap = 1
          rel_gap = 1
        else:
          primal_gap = round(min(parameter_instance_results_df['primal_gap']), 3)
          rel_gap = round(min(parameter_instance_results_df['rel_gap']), 3)
        parameter_primal_gaps.append(primal_gap)
        parameter_rel_gaps.append(rel_gap)
      average_primal_gap = round(sum(parameter_primal_gaps) / len(parameter_primal_gaps), 3)
      average_rel_gap = round(sum(parameter_rel_gaps) / len(parameter_rel_gaps), 3)
      averages_string = averages_string + f" & & {average_primal_gap} & {average_rel_gap}"

    header_string = "instance_class & "
    for parameter in parameters:
      header_string = header_string + f"{parameter_set_names_strings[parameter]} & & "
    print(averages_string)
  print(header_string)

def print_times_table(parameter_set_names_strings, gap_results_df):
  instances = list(set(gap_results_df['instance']))
  instances.sort()
  parameters = [parameter[0] for parameter in parameter_set_names_strings.items()]

  time_aggregations = [300, 3600]
  for time_aggregation in time_aggregations:
    aggregation_gap_results_df = gap_results_df[gap_results_df['time'] <= time_aggregation]

    averages_string = f"{time_aggregation}s & "
    for parameter in parameters:
      parameter_results_df = aggregation_gap_results_df[aggregation_gap_results_df['parameter'] == parameter]
      parameter_primal_gaps = []
      parameter_rel_gaps = []
      for instance in instances:
        parameter_instance_results_df = parameter_results_df[parameter_results_df['instance'] == instance]
        if parameter_instance_results_df.empty:
          primal_gap = 1
          rel_gap = 1
        else:
          primal_gap = round(min(parameter_instance_results_df['primal_gap']), 3)
          rel_gap = round(min(parameter_instance_results_df['rel_gap']), 3)
        parameter_primal_gaps.append(primal_gap)
        parameter_rel_gaps.append(rel_gap)
      average_primal_gap = round(sum(parameter_primal_gaps) / len(parameter_primal_gaps), 3)
      average_rel_gap = round(sum(parameter_rel_gaps) / len(parameter_rel_gaps), 3)
      averages_string = averages_string + f"{average_primal_gap} & {average_rel_gap} & & "

    header_string = "instance & "
    for parameter in parameters:
      header_string = header_string + f"{parameter_set_names_strings[parameter]} & & "
    print(averages_string)
  print(header_string)


def print_cuts_info(parameter_set_names_strings, ce_results_df):
  special_instances = ['A-n33-k6.vrp','A-n34-k5.vrp','A-n37-k5.vrp','A-n37-k6.vrp','A-n38-k5.vrp','A-n39-k5.vrp','A-n39-k6.vrp','A-n45-k6.vrp','A-n45-k7.vrp','A-n48-k7.vrp','A-n53-k7.vrp','A-n55-k9.vrp','B-n31-k5.vrp','B-n34-k5.vrp','E-n51-k5.vrp','P-n16-k8.vrp','P-n20-k2.vrp','P-n22-k2.vrp','P-n40-k5.vrp','P-n45-k5.vrp','P-n50-k10.vrp','P-n50-k7.vrp','P-n50-k8.vrp','P-n51-k10.vrp','P-n55-k10.vrp','P-n55-k15.vrp','P-n55-k7.vrp','P-n55-k8.vrp','P-n60-k10.vrp','P-n60-k15.vrp','P-n65-k10.vrp']

  size_diff_rcc = []
  size_diff_src3 = []
  size_diff_src34 = []
  size_diff_src345 = []

  rcc_rcc_cuts = []

  src3_rcc_cuts = []
  src3_src3_cuts = []
 
  src34_rcc_cuts = []
  src34_src3_cuts = []
  src34_src4_cuts = []
 
  src345_rcc_cuts = []
  src345_src3_cuts = []
  src345_src4_cuts = []
  src345_src5_cuts = []

  for instance in special_instances:
    print(instance)
    instance_df = ce_results_df[ce_results_df['instance'] == instance]
    print(instance_df)

    none_df = instance_df[instance_df['parameter'].str.contains("CE_None_0")]
    print(none_df)
    size_none = max(none_df['numArcs'])

    rcc_df = instance_df[instance_df['parameter'].str.contains("CE_RCC_0")]
    size_rcc = max(rcc_df['numArcs'])
    size_diff_rcc.append((size_rcc - size_none) / size_none)
    rcc_rcc_cuts.append(max(rcc_df['rcc']))

    src3_df = instance_df[instance_df['parameter'].str.contains("CE_RCCSRC3_0")]
    size_src3 = max(src3_df['numArcs'])
    size_diff_src3.append((size_src3 - size_none) / size_none)
    src3_rcc_cuts.append(max(src3_df['rcc']))
    src3_src3_cuts.append(max(src3_df['src3']))

    src34_df = instance_df[instance_df['parameter'].str.contains("CE_RCCSRC3SRC4_0")]
    size_src34 = max(src34_df['numArcs'])
    size_diff_src34.append((size_src34 - size_none) / size_none)
    src34_rcc_cuts.append(max(src34_df['rcc']))
    src34_src3_cuts.append(max(src34_df['src3']))
    src34_src4_cuts.append(max(src34_df['src4']))

    src345_df = instance_df[instance_df['parameter'].str.contains("CE_RCCSRC3SRC4SRC5_0")]
    size_src345 = max(src345_df['numArcs'])
    size_diff_src345.append((size_src345 - size_none) / size_none)
    src345_rcc_cuts.append(max(src345_df['rcc']))
    src345_src3_cuts.append(max(src345_df['src3']))
    src345_src4_cuts.append(max(src345_df['src4']))
    src345_src5_cuts.append(max(src345_df['src5']))

  print(f'rcc size diff: {sum(size_diff_rcc) / len(size_diff_rcc)}')
  print(f'src3 size diff: {sum(size_diff_src3) / len(size_diff_src3)}')
  print(f'src34 size diff: {sum(size_diff_src34) / len(size_diff_src34)}')
  print(f'src345 size diff: {sum(size_diff_src345) / len(size_diff_src345)}')
  
  print(f'rcc rcc: {sum(rcc_rcc_cuts) / len(rcc_rcc_cuts)}')
  
  print(f'src3 rcc: {sum(src3_rcc_cuts) / len(src3_rcc_cuts)}')
  print(f'src3 src3: {sum(src3_src3_cuts) / len(src3_src3_cuts)}')
 
  print(f'src34 rcc: {sum(src34_rcc_cuts) / len(src34_rcc_cuts)}')
  print(f'src34 src3: {sum(src34_src3_cuts) / len(src34_src3_cuts)}')
  print(f'src34 src4: {sum(src34_src4_cuts) / len(src34_src4_cuts)}')
 
  print(f'src345 rcc: {sum(src345_rcc_cuts) / len(src345_rcc_cuts)}')
  print(f'src345 src3: {sum(src345_src3_cuts) / len(src345_src3_cuts)}')
  print(f'src345 src4: {sum(src345_src4_cuts) / len(src345_src4_cuts)}')
  print(f'src345 src4: {sum(src345_src5_cuts) / len(src345_src5_cuts)}')
