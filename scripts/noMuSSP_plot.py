import os
import math
import tabulate
import pandas
import numpy
import matplotlib.pyplot as plt
import re
import collections
import statistics

logs_dir = "/Users/akarahal/Desktop/dd_vrptw/logs/"
instance_dir = "/Users/akarahal/Desktop/dd_vrptw/instances/X123456789/"

instances = []
for instance in os.listdir(instance_dir):
  instances.append(instance)

sizes = []
muSSP_runtimes = []
SSP_runtimes = []
colelim_pattern = re.compile("COMPARISON size\[([0-9]+)\] muSSP\[([0-9]+.*)\] SSP\[([0-9]+.*)\]")
for instance in instances:
  log_file_name = logs_dir + 'noMuSSP' + '/' + instance + '.log'
  if not os.path.exists(log_file_name):
    continue

  log_file = open(log_file_name, "r")
  for line in log_file:
    colelim_match = colelim_pattern.match(line)
    if colelim_match:
      sizes.append(int(colelim_match.group(1)))
      muSSP_runtimes.append(float(colelim_match.group(2)))
      SSP_runtimes.append(float(colelim_match.group(3)))

diffs = [ssp / muSSP for (ssp,muSSP) in zip(SSP_runtimes, muSSP_runtimes)]
print(statistics.mean(diffs))
plt.scatter(SSP_runtimes, muSSP_runtimes, marker = 'x')
 
plt.xlabel('SSP Runtime (s)')
plt.xscale('log')
plt.xlim([1,max(max(SSP_runtimes),max(muSSP_runtimes))])

plt.ylabel('muSSP Runtime (s)')
plt.yscale('log')
plt.ylim([1,max(max(SSP_runtimes),max(muSSP_runtimes))])

plt.plot([0,10000],[0,10000])

plt.show()
