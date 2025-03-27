import static
import read_logs
import calculate_metrics
import pandas

# Input Information
root_directory = "/Users/akarahal/Desktop/dd_vrptw/"

'''
parameter_set_names_strings = {
"Vrp-Set-HG-C12-124_primal_param_files/CE_5_1_0" : 'CE_5',
"Vrp-Set-HG-C12-124_primal_param_files/CE_30_1_0": 'CE_30',
"Vrp-Set-HG-C12-124_primal_param_files/CE_MIP_5_1_0" : "CE_MIP_5",
"Vrp-Set-HG-C12-124_primal_param_files/CE_MIP_30_1_0" : "CE_MIP_30",
"Vrp-Set-HG-C12-124_primal_param_files/CE_MIP_LNS_5_1_0" : "CE_MIP_LNS_5",
"Vrp-Set-HG-C12-124_primal_param_files/CE_MIP_LNS_30_1_0" : "CE_MIP_LNS_30"
#"PyVRP" : "PyVRP"
}
'''

parameter_set_names_strings = {
"X_primal_param_files/CE_MIP_LNS_5_1_0" : "CE_MIP_LNS_5",
"PyVRP" : "PyVRP"
}

parameter_set_names_strings = {
"Vrp-Set-HG_primal_param_files/CE_MIP_LNS_5_1_0" : "CE_MIP_LNS_5",
"PyVRP" : "PyVRP"
}

parameter_set_names_strings = {
"X_primal_param_files/CE_MIP_LNS_SRC_5_1_0" : "CE_MIP_LNS_SRC_5_1_0",
"X_primal_param_files/CE_MIP_LNS_VA_5_1_0" : "CE_MIP_LNS_VA_5_1_0",
"X_primal_param_files/CE_MIP_LNS_VA_SRC_5_1_0" : "CE_MIP_LNS_VA_SRC_5_1_0",
"X_primal_param_files/CE_MIP_LNS_VA_SRC_All_5_1_0" : "CE_MIP_LNS_VA_SRC_All_5_1_0",
"X_primal_param_files/CE_MIP_LNS_VA_SRC_Restarts_5_1_0" : "CE_MIP_LNS_VA_SRC_Restarts_5_1_0",
"PyVRP" : "PyVRP"
}

'''
parameter_set_names_strings = {
"Vrp-Set-HG-C12-124_primal_param_files/removals/CE_Random_0/" : "Random",
"Vrp-Set-HG-C12-124_primal_param_files/removals/CE_Worst_0/" : "Worst",
"Vrp-Set-HG-C12-124_primal_param_files/removals/CE_Shaw_0/" : "Shaw",
"Vrp-Set-HG-C12-124_primal_param_files/removals/CE_SequenceRandom_0/" : "Sequence+Random",
"Vrp-Set-HG-C12-124_primal_param_files/removals/CE_SequenceWorst_0/" : "Sequence+Worst",
"Vrp-Set-HG-C12-124_primal_param_files/removals/CE_SequenceShaw_0/" : "Sequence+Shaw",
}

parameter_set_names_strings = {
"Vrp-Set-HG-C12-124_primal_param_files/ablations/CE_None_0/" : "None",
"Vrp-Set-HG-C12-124_primal_param_files/ablations/CE_Intra_0/" : "IntraSequenceSwaps",
"Vrp-Set-HG-C12-124_primal_param_files/ablations/CE_Truncated_0/" : "TruncatedSequences",
}

parameter_set_names_strings = {
"Vrp-Set-HG-C12-124_primal_param_files/times/CE_5_0/" : "5s",
"Vrp-Set-HG-C12-124_primal_param_files/times/CE_30_0/" : "30s",
"Vrp-Set-HG-C12-124_primal_param_files/times/CE_60_0/" : "60s",
}
'''

'''
parameter_set_names_strings = {
"Vrp-Set-HG-C12-124_primal_param_files/CE_5_1_0" : 'CE',
"Vrp-Set-HG-C12-124_primal_param_files/CE_MIP_5_1_0" : "CE+MIP",
"Vrp-Set-HG-C12-124_primal_param_files/CE_MIP_LNS_5_1_0" : "CE+MIP+LNS",
}

'''
parameter_set_names_strings = {
  "X_primal_param_files/best/CE_Random_0/" : "CE",
  "PyVRP" : "PyVRP",
  "VRPSolver" : "VRPSolver"
}

parameter_set_names_strings = {
  "Vrp-Set-HG_primal_param_files/best/CE_Random_0/" : "CE",
  "PyVRP" : "PyVRP",
  "VRPSolver" : "VRPSolver"
}

'''
parameter_set_names_strings = {
  "pdp_600_primal_param_files/best/CE_Random_0/" : "CE_Random",
  "pdp_600_primal_param_files/best/CE_SequenceRandom_0/" : "CE_Sequence",
  "PyVRP" : "PyVRP",
  "VRPSolver" : "VRPSolver"
}
'''

parameter_set_names = [x[0] for x in parameter_set_names_strings.items()]

instance_set_names = ["Solomon", "HG", "CVRP", "X", "PDPTW"]
instance_set_name = instance_set_names[1]

vrpsolver_file_names = ['cvrp-X-noub1.log','cvrp-X-noub2.log','cvrp-X-noub3.log','cvrp-X-noub4.log','cvrp-X-noub5.log']
pyvrp_file_names = ['cvrp_instances1.log','cvrp_instances2.log','cvrp_instances3.log','cvrp_instances4.log','vrptw_instances1.log','vrptw_instances2.log','vrptw_instances3.log','vrptw_instances4.log','vrptw_instances5.log','vrptw_instances6.log']

times_for_metrics = range(0,3600,100)

# Run Evaluation for the given Experiments and Instance Set
instances = read_logs.get_instances(instance_set_name, root_directory)

ce_results_df = read_logs.get_column_elimination_results(instances, parameter_set_names, root_directory)
ce_results_df = ce_results_df[['instance','parameter','lb','ub','time']]

if "VRPSolver" in parameter_set_names_strings:
  vrpsolver_results_df = read_logs.get_vrpsolver_results(instances, vrpsolver_file_names, instance_set_name, root_directory)
  ce_results_df = pandas.concat([ce_results_df, vrpsolver_results_df])

if "PyVRP" in parameter_set_names_strings:
  pyvrp_results_df = read_logs.get_pyvrp_results(instances, pyvrp_file_names, instance_set_name, root_directory)
  ce_results_df = pandas.concat([ce_results_df, pyvrp_results_df])

gap_results_df = calculate_metrics.calculate_gaps(ce_results_df, times_for_metrics)

average_gap_results_df = calculate_metrics.calculate_average_gaps(gap_results_df, times_for_metrics)

#calculate_metrics.plot_average_gaps(instance_set_name, parameter_set_names_strings, average_gap_results_df)

#calculate_metrics.plot_gaps_by_instance(parameter_set_names_strings, gap_results_df)

calculate_metrics.print_instance_table(parameter_set_names_strings, gap_results_df)

calculate_metrics.print_aggregated_table(instance_set_name, parameter_set_names_strings, gap_results_df)

calculate_metrics.print_times_table(parameter_set_names_strings, gap_results_df)
