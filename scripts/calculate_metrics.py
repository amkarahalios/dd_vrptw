import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections
import static

def calculate_gaps(ce_results_df, vrpsolver_results_df, times_for_metrics):
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
        if not instance_parameter_results.empty:
          for end_time in times_for_metrics:
              instance_parameter_time_results = instance_parameter_results[instance_parameter_results['time'] < end_time]
              if not instance_parameter_time_results.empty:
                instance_parameter_ubs = instance_parameter_time_results['ub']
                ub = min(instance_parameter_ubs)
                primal_gap = abs(ub - best_known) / max(abs(best_known), abs(ub))

                instance_parameter_lbs = instance_parameter_time_results['lb']
                lb = max(instance_parameter_lbs)
                rel_gap = (ub-lb) / max(ub,lb)

                gap_info = {'instance': instance, 'parameter' : parameter, 'rel_gap' : rel_gap, 'primal_gap' : primal_gap, 'time' : end_time} 
              else:
                gap_info = {'instance': instance, 'parameter' : parameter, 'rel_gap' : 1, 'primal_gap' : 1, 'time' : end_time} 
              gap_results.append(gap_info)

    if instance in static.pyvrp_solutions:
      ub = static.pyvrp_solutions[instance]
      primal_gap = abs(ub - best_known) / max(abs(best_known), abs(ub))
      rel_gap = 1
      for end_time in times_for_metrics:
        gap_info = {'instance': instance, 'parameter' : "PyVRP", 'rel_gap' : rel_gap, 'primal_gap' : primal_gap, 'time' : end_time} 
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
  parameter_set = set(average_gap_results_df['parameter'])

  for parameter in parameter_set:
    parameter_df = average_gap_results_df[average_gap_results_df['parameter'] == parameter]
    plt.plot(list(parameter_df['time']), list(parameter_df['primal_gap']), label=parameter_set_names_strings[parameter])

  if instance_set_name == "Solomon":
    plt.title("Average p(t) for VRPTW Solomon")
  elif instance_set_name == "HG":
    plt.title("Average p(t) for VRPTW HG")
  elif instance_set_name == "CVRP":
    plt.title("Average p(t) for CVRP")
  elif instance_set_name == "PDPTW":
    plt.title("Average p(t) for PDPTW")

  plt.ylabel('p(t)')
  plt.xlabel('time (s)')
  plt.legend()
  plt.show()

def print_instance_table(parameter_set_names_strings, gap_results_df):
  instances = list(set(gap_results_df['instance']))
  instances.sort()
  parameters = set(gap_results_df['parameter'])
  for instance in instances:
    best_known = static.instance_upper_bounds[instance]
    instance_results_df = gap_results_df[gap_results_df['instance'] == instance]
    result_string = f"{instance} & {best_known} & & "
    for parameter in parameters:
      instance_parameter_results_df = instance_results_df[instance_results_df['parameter'] == parameter]
      primal_gap = round(min(instance_parameter_results_df['primal_gap']), 3)
      rel_gap = round(min(instance_parameter_results_df['rel_gap']), 3)
      result_string = result_string + f"{primal_gap} & {rel_gap} & & "

    print(result_string)

  header_string = "instance & best_known & & "
  for parameter in parameters:
    header_string = header_string + f"{parameter_set_names_strings[parameter]} & & "

  print(header_string)

def print_aggregated_table(parameter_set_names_strings, gap_results_df):
  instances = list(set(gap_results_df['instance']))
  instances.sort()
  parameters = set(gap_results_df['parameter'])
  averages_string = "instance & & "
  for parameter in parameters:
    parameter_results_df = gap_results_df[gap_results_df['parameter'] == parameter]
    parameter_primal_gaps = []
    parameter_rel_gaps = []
    for instance in instances:
      parameter_instance_results_df = parameter_results_df[parameter_results_df['instance'] == instance]
      primal_gap = round(min(parameter_instance_results_df['primal_gap']), 3)
      rel_gap = round(min(parameter_instance_results_df['rel_gap']), 3)
      parameter_primal_gaps.append(primal_gap)
      parameter_rel_gaps.append(rel_gap)
    average_primal_gap = round(sum(parameter_primal_gaps) / len(parameter_primal_gaps), 3)
    average_rel_gap = round(sum(parameter_rel_gaps) / len(parameter_rel_gaps), 3)
    averages_string = averages_string + f"{average_primal_gap} & {average_rel_gap} & & "

  header_string = "instance & & "
  for parameter in parameters:
    header_string = header_string + f"{parameter_set_names_strings[parameter]} & & "
  print(header_string)
