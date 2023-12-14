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
"COG-10teams.mcol":0,
"COG-air04.mcol":377,
"COG-air05.mcol":0,
"COG-atlanta-ip.mcol":15,
"COG-cap6000.mcol":14,
"COG-ds.mcol":500,
"COG-gesa2-o.mcol":13,
"COG-misc07.mcol":39,
"COG-mkc.mcol":169,
"COG-mod011.mcol":13,
"COG-mzzv11.mcol":101,
"COG-mzzv42z.mcol":91,
"COG-net12.mcol":17,
"COG-nsrand-ipx.mcol":0,
"COG-opt1217.mcol":0,
"COG-rd-rplusc-21.mcol":109,
"COG-rout.mcol":32,
"COG-swath.mcol":0
}

instance_lower_bounds = {
"COG-10teams.mcol":0,
"COG-air04.mcol":377,
"COG-air05.mcol":0,
"COG-atlanta-ip.mcol":15,
"COG-cap6000.mcol":14,
"COG-ds.mcol":500,
"COG-gesa2-o.mcol":12,
"COG-misc07.mcol":36,
"COG-mkc.mcol":169,
"COG-mod011.mcol":12,
"COG-mzzv11.mcol":101,
"COG-mzzv42z.mcol":91,
"COG-net12.mcol":17,
"COG-nsrand-ipx.mcol":0,
"COG-opt1217.mcol":0,
"COG-rd-rplusc-21.mcol":109,
"COG-rout.mcol":30,
"COG-swath.mcol":0
}

instance_times = {
"COG-10teams.mcol":3600,
"COG-air04.mcol":1.8,
"COG-air05.mcol":3600,
"COG-atlanta-ip.mcol":1844,
"COG-cap6000.mcol":304,
"COG-ds.mcol":6.5,
"COG-gesa2-o.mcol":3600,
"COG-misc07.mcol":3600,
"COG-mkc.mcol":0.1,
"COG-mod011.mcol":3600,
"COG-mzzv11.mcol":0.1,
"COG-mzzv42z.mcol":0.1,
"COG-net12.mcol":1301,
"COG-nsrand-ipx.mcol":3600,
"COG-opt1217.mcol":3600,
"COG-rd-rplusc-21.mcol":0.0,
"COG-rout.mcol":3600,
"COG-swath.mcol":3600
}

max_weight_cliques = {
"COG-10teams.mcol":73,
"COG-air04.mcol":377,
"COG-air05.mcol":413,
"COG-atlanta-ip.mcol":15,
"COG-cap6000.mcol":14,
"COG-ds.mcol":1,
"COG-fiber.mcol":19,
"COG-gesa2.mcol":18,
"COG-gesa2-o.mcol":12,
"COG-glass4.mcol":22,
"COG-harp2.mcol":41,
"COG-misc07.mcol":36,
"COG-mkc.mcol":169,
"COG-mod011.mcol":12,
"COG-momentum1.mcol":50,
"COG-momentum2.mcol":58,
"COG-momentum3.mcol":66,
"COG-msc98-ip.mcol":20,
"COG-mzzv11.mcol":101,
"COG-mzzv42z.mcol":91,
"COG-net12.mcol":17,
"COG-nsrand-ipx.mcol":30,
"COG-opt1217.mcol":26,
"COG-p2756.mcol":24,
"COG-rd-rplusc-21.mcol":109,
"COG-roll3000.mcol":18,
"COG-rout.mcol":30,
"COG-seymour.mcol":18,
"COG-swath.mcol":317
}

# Log lines for DDSolver
#STATS - lpIterations[1] lagIterations[32] sspIterations[186] numSeparations[0] compileTime[240] sspSolveTime[349] lpSolveTime[0] lb[794.025] ub[1e+10] numArcs: [87254] numFixed: [0] time: [590]

colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] numArcs: \[([0-9]+)\] numFixed: \[([0-9]+)\] time: \[([0-9]+)\]")
#colelim_pattern = re.compile("STATS - lpIterations\[([0-9]+)\] lagIterations\[([0-9]+)\] sspIterations\[([0-9]+)\] numSeparations\[([0-9]+)\].*compileTime\[([0-9]+)\] sspSolveTime\[([0-9]+)\] lpSolveTime\[([0-9]+)\] lb\[([0-9]+.*)\] ub\[([0-9]+.*)\] size: \[([0-9]+)\] time: \[([0-9]+)\]")

logs_dir = "/Users/akarahal/Desktop/dd_graph_color/logs/"
test_set = ["col_elim_COGs_lp_x_0_MIP_3600"]

instances = []
instance_dir = "/Users/akarahal/Desktop/dd_graph_color/instances/COGs/"

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
        lpIterations = int(colelim_match.group(1))
        lagIterations = int(colelim_match.group(2))
        sspIterations = colelim_match.group(3)
        numSeparations = int(colelim_match.group(4))
        compileTime = colelim_match.group(5)
        sspSolveTime = colelim_match.group(6)
        lpSolveTime = colelim_match.group(7)
        lb = int(colelim_match.group(8))
        ub = int(colelim_match.group(9))
        numArcs = int(colelim_match.group(10))
        numFixed = int(colelim_match.group(11))
        time = float(colelim_match.group(12))
        time_result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'lagIterations': lagIterations, 'numSep' : numSeparations, 'lb' : lb, 'ub' : ub, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
        time_results.append(time_result)
 
    if col_elim:
      result = {'instance': instance, 'test': test, 'iterations' : lpIterations, 'lagIterations': lagIterations, 'numSep' : numSeparations, 'lb' : lb, 'ub' : ub, 'numArcs': numArcs, 'numFixed': numFixed, 'time' : time}
      results.append(result)
    else:
      result = {'instance': instance, 'test': test, 'iterations' : numpy.nan, 'lagIterations': numpy.nan, 'numSep' : numpy.nan, 'lb' : numpy.nan, 'ub' : numpy.nan, 'numArcs': numpy.nan, 'numFixed': numpy.nan, 'time' : numpy.nan}
      results.append(result)

# give table of results
results_df = pandas.DataFrame(results)
results_df.sort_values(by=['instance','lb'],inplace=True)
results_df.reset_index(drop=True,inplace=True)
table_results_df = results_df[['instance','test','lb','ub','time','iterations','numSep','lagIterations','numArcs']]
print(tabulate.tabulate(table_results_df, headers=table_results_df.columns))

# get sizes for sop instances
pattern = re.compile("p edges ([0-9]+) ([0-9]+)")
instance_num_vertices = {}
instance_num_edges = {}
for instance in instances:
  instance_file_name = instance_dir + '/' + instance
  instance_file = open(instance_file_name, "r")
  for line in instance_file:
    match = pattern.match(line)
    if match:
      num_vertices = int(match.group(1))
      instance_num_vertices[instance] = num_vertices
 
      num_edges = int(match.group(2))
      instance_num_edges[instance] = num_edges
      break

# create output table for SOA comparison
for i, row in table_results_df.iterrows():
  instance = row['instance']
  if instance not in instance_lower_bounds:
    continue
  instance_name = instance.rsplit('.',1)[0]
  instance_name = instance_name.replace("_","\\_")
  lb_value = row['lb']
  if not math.isnan(lb_value):
    lb_value = int(lb_value)
  else:
    lb_value = '-'

  ub_value = row['ub']
  if not math.isnan(ub_value):
    ub_value = int(ub_value)
  else:
    ub_value = '-'

  numArcs = row['numArcs']
  if numArcs > 100000:
    numArcs = '{:.2e}'.format(numArcs) 

  time = row['time']
  if not math.isnan(time):
    time = int(time)
  else:
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

  instance_lower_bound = instance_lower_bounds[instance]
  if instance_lower_bound == 0:
    instance_lower_bound = '-'
 
  instance_upper_bound = instance_upper_bounds[instance]
  if instance_upper_bound == 0:
    instance_upper_bound = '-'

  instance_upper_bound = instance_upper_bounds[instance]
  if instance_upper_bound == 0:
    instance_upper_bound = '-'
 
  if instance_time == 0 and instance_upper_bound == '-' and instance_lower_bound == '-':
    instance_time = '-'
 
  if instance_upper_bound != '-' and ub_value != '-' and instance_upper_bound > ub_value:
    print(f"{instance_name} & {instance_num_vertices[instance]} & {instance_num_edges[instance]} & {max_weight_cliques[instance]} & & {instance_lower_bound} & {instance_upper_bound} & {instance_time} & & {lb_value} & \\textbf\u007b{ub_value}\u007d & {numLpIterations} & {numSeparations} & {time} \\\\")
  elif ub_value == lb_value and instance_upper_bound == '-' and lb_value != '-':
    print(f"{instance_name} & {instance_num_vertices[instance]} & {instance_num_edges[instance]} & {max_weight_cliques[instance]} & & {instance_lower_bound} & {instance_upper_bound} & {instance_time} & & {lb_value} & \\textbf\u007b{ub_value}\u007d & {numLpIterations} & {numSeparations} & {time} \\\\")
  else:
    print(f"{instance_name} & {instance_num_vertices[instance]} & {instance_num_edges[instance]} & {max_weight_cliques[instance]} & & {instance_lower_bound} & {instance_upper_bound} & {instance_time} & & {lb_value} & {ub_value} & {numLpIterations} & {numSeparations} & {time} \\\\")
