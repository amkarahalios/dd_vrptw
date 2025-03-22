import static
import read_logs
import calculate_metrics

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

'''
parameter_set_names_strings = {
"Vrp-Set-HG_primal_param_files/CE_MIP_LNS_5_1_0" : "CE_MIP_LNS_5",
"PyVRP" : "PyVRP"
}
'''

parameter_set_names = [x[0] for x in parameter_set_names_strings.items()]

instance_set_names = ["Solomon", "HG", "CVRP", "X", "PDPTW"]
instance_set_name = instance_set_names[3]

vrpsolver_file_names = []

times_for_metrics = range(0,3600,20)

# Run Evaluation for the given Experiments and Instance Set
instances = read_logs.get_instances(instance_set_name, root_directory)

ce_results_df = read_logs.get_column_elimination_results(instances, parameter_set_names, root_directory)

vrpsolver_results_df = read_logs.get_vrpsolver_results(instances, vrpsolver_file_names, instance_set_name, root_directory)

gap_results_df = calculate_metrics.calculate_gaps(ce_results_df, vrpsolver_results_df, times_for_metrics)

average_gap_results_df = calculate_metrics.calculate_average_gaps(gap_results_df, times_for_metrics)

calculate_metrics.plot_average_gaps(instance_set_name, parameter_set_names_strings, average_gap_results_df)

calculate_metrics.plot_gaps_by_instance(parameter_set_names_strings, gap_results_df)

calculate_metrics.print_instance_table(parameter_set_names_strings, gap_results_df)

calculate_metrics.print_aggregated_table(parameter_set_names_strings, gap_results_df)
