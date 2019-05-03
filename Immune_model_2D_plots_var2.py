# Derek Park
# April2019
import numpy
import scipy
from scipy import integrate
#import matplotlib
#from matplotlib import pyplot as plt
import pandas
import os
import datetime
import multiprocessing
import MLImmuneModel
import canonicalparams
import time
import sys

class MyException(Exception):

    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

# The method to evaluate the model with a specific parameter combination
def eval_model(vars):



    # Create the parameter set
    the_params = canonicalparams.CanonicalModelParameters()

    # Set the parameters according to vars
    for variable in vars.keys():
        setattr(the_params, variable, vars[variable])


    # Create the model
    model = MLImmuneModel.ImmuneModel(the_params, None, None)

    # Leadup growth
    model.initialize_model()
    # Treat with the level of chemo for 10 cycles
    this_run_chemo_strength = vars['chemo_strength']

    for i in range(10):

        model.treat(None, override=True, manual_chemo=this_run_chemo_strength)


    # Get the results
    simulation_results = model.get_results()

    simulation_results.to_csv(vars['save_name'])

    final_tumor_size = simulation_results['Ts'].tolist()[-1]

    percent_change_in_size = ((final_tumor_size)/10**8.0 - 1.0) * 100

    vars['Result'] = percent_change_in_size

    RECIST_category = -100

    if percent_change_in_size <= -99.0:
        RECIST_category = -2
    elif percent_change_in_size <= -30.0:
        RECIST_category = -1
    elif percent_change_in_size <= 20.0:
        RECIST_category = 0
    else:
        RECIST_category = 1

    vars['RECIST'] = RECIST_category

    vars['final size'] = final_tumor_size
    return vars

# def generate_2D_data(var1_name, var1_lower_bound, var1_upper_bound, var1_range_length, var1_scale,
#                     var2_name, var2_lower_bound, var2_upper_bound, var2_range_length, var2_scale,
#                     save_location, save_title):
def generate_spec_data():
    #The objective of this method is to generate 2D parameter scanning plots of the immune model.



    # print 'First variable name is: ' + str(var1_name) + '. It is to be varied in the range [' + str(var1_lower_bound) + ', ' + str(var1_upper_bound) + '] in a ' + var1_scale + ' scale with ' + str(var1_range_length) + ' data points.'
    # print 'First variable name is: ' + str(var2_name) + '. It is to be varied in the range [' + str(var2_lower_bound) + ', ' + str(var2_upper_bound) + '] in a ' + var2_scale + ' scale with ' + str(var2_range_length) + ' data points.'
    #
    # # Check to see that these are valid variables of the immune model
    # print 'Checking to see that these are valid parameters of the model...'
    #
    # parameters_found = 0
    #
    # canon_params = canonicalparams.CanonicalModelParameters()
    # model_params = dir(canon_params)
    #
    #
    # if var1_name in model_params:
    #     print 'Parameter 1: \'' + var1_name + '\' is found.'
    #     parameters_found += 1
    # if var2_name in model_params:
    #     print 'Parameter2: \'' + var2_name + '\' is found.'
    #     parameters_found += 1
    #
    # if parameters_found != 2:
    #
    #     raise MyException('Both parameters were not found')
    #
    # # Generate the variable 1 list of values:
    #
    # var1_range = []
    #
    # if var1_scale == 'log':
    #
    #     var1_range = numpy.logspace(numpy.log10(var1_lower_bound), numpy.log10(var1_upper_bound), num=var1_range_length, endpoint=False)
    # else:
    #     var1_range = numpy.linspace(var1_lower_bound, var1_upper_bound, num=var1_range_length, endpoint=False)
    #
    #
    # # Generate the variable 2 list of values:
    #
    # var2_range = []
    #
    # if var2_scale == 'log':
    #
    #     var2_range = numpy.logspace(numpy.log10(var2_lower_bound), numpy.log10(var2_upper_bound), num=var2_range_length, endpoint=False)
    # else:
    #     var2_range = numpy.linspace(var2_lower_bound, var2_upper_bound, num=var2_range_length, endpoint=False)
    #
    #
    # # Create a list of dictionaries that are to be passed to the function to evaluate the model
    #
    # print 'Generating the list of parameters to be passed to the model.'

    parameter_dictionary_list = []


    # for v1 in var1_range:
    #
    #     for v2 in var2_range:
    #
    #         curr_dict = {}
    #         curr_dict[var1_name] = v1
    #         curr_dict[var2_name] = v2
    #         ###################################################### PARAMETERS FOR SPECIFIC VACCINE TRIALS #######################################################
    #         #curr_dict['vaccine_strength'] = 10.0
    #         #curr_dict['M_0'] = 10**6.0
    #         #curr_dict['vaccine_half_life'] = 3.0
    #         ###################################################### PARAMETERS FOR SPECIFIC VACCINE TRIALS #######################################################
    #         parameter_dictionary_list.append(curr_dict)

    param_names = [ 'rho',#
                    'c',#
                    'gamma',#
                    'delta_R',#
                    'r_M',#
                    'omega']
    param_vals = [[0, 0.05, 1.0],
                  [0.001, 0.01, 0.1],
                  [10., 100., 1000.],
                  [0.01, 0.1, 1.0],
                  [0.001, 0.01, 0.1],
                  [0.001, 0.01, 0.1]]
    constant_chemo = 0.25

    run = True
    for name, val_set in zip(param_names, param_vals):

        if run:
            for val in val_set:
                curr_dict = {}
                curr_dict['chemo_strength'] = constant_chemo
                curr_dict[name] = val
                curr_dict['save_name'] = 'chemo_strength=%f_%s=%f.csv' %(constant_chemo, name, val)
                parameter_dictionary_list.append(curr_dict)
            

    task_size = len(parameter_dictionary_list)
    #print parameter_dictionary_list

    print str(task_size) + ' simulations to be completed.'


    start_time = datetime.datetime.now()

    print 'Simulations started at: ' + str(start_time)



    parameter_variation_results = []

    p = multiprocessing.Pool()
    rs = p.imap_unordered(eval_model, parameter_dictionary_list)
    num_CPUs = multiprocessing.cpu_count()
    p.close() # No more work
    while (True):
        completed = rs._index
        if (completed == task_size): break
        curr_time = datetime.datetime.now()
        if completed != 0:
            time_per_sim = (curr_time - start_time)/completed
            output_string = ""
            output_string = "Waiting for " + str(task_size - completed) + " tasks to complete on " + str(num_CPUs) + " cores. Projected finish is: " + str(curr_time + time_per_sim*(task_size-completed))
            #f = open(save_title+"_progress.txt",'w')
            f = open("specific_runs_progress.txt",'w')
            print(output_string)
            try:
                f.write(output_string)
            except:
                print("Error writing")
            f.close()

        time.sleep(5)

    for i in rs:
        if i:
            parameter_variation_results.append(i)

    end_time = datetime.datetime.now()

    print 'Simulations finished at ' + str(end_time)
    print 'Time per simulation for ' + str(task_size) + ' simulations was: ' + str((end_time - start_time)/task_size) + '. \nTotal time was: ' + str(end_time- start_time) + '.'

    # Now organize the output into a table for formatting
    results_list_form = []

    for result in parameter_variation_results:
        curr_result_formatted = []
        curr_result_formatted.append(result[var1_name])
        curr_result_formatted.append(result[var2_name])
        curr_result_formatted.append(result['Result'])
        curr_result_formatted.append(result['RECIST'])
        curr_result_formatted.append(result['final size'])
        #print result
        #print '\n'
        results_list_form.append(curr_result_formatted)

    results_frame = pandas.DataFrame(results_list_form, columns=[var1_name, var2_name, 'Percent Change in Tumor Size', 'RECIST', 'final size'])
    try:
        results_frame.to_csv(os.path.join(save_location, save_title))
    except IOError:
        results_frame.to_csv('%s_error.csv' %(save_title))



if __name__ == '__main__':

    # var1_name = 'chemo_strength'
    # var2_name = 'r_T'
    #
    # var1_lower_bound = 0
    # var1_upper_bound = 0.9
    # var1_scale = 'linear'
    # var1_range_length = 150
    #
    # var2_lower_bound = 500
    # var2_upper_bound = 1500
    # var2_scale = 'linear'
    # var2_range_length = 150
    #
    # save_location = '/home/parkds/PhD_projects/Chemo-Immune_model/Vaccine_Trials/constant_regime_vac_strength=10/chemo_vs_offset_vac=10.csv'
    # save_location = '.'
    # save_title = 'toff8_chemo_vs_rt.csv'

    # if len(sys.argv) > 1:
    #     print(sys.argv)
    #     var1_name = str(sys.argv[1])
    #     var1_lower_bound = float(sys.argv[2])
    #     var1_upper_bound = float(sys.argv[3])
    #     var1_range_length = int(sys.argv[4])
    #     var1_scale = str(sys.argv[5])
    #
    #     var2_name = str(sys.argv[6])
    #     var2_lower_bound = float(sys.argv[7])
    #     var2_upper_bound = float(sys.argv[8])
    #     var2_range_length = int(sys.argv[9])
    #     var2_scale = str(sys.argv[10])
    #
    #     save_location = str(sys.argv[11])
    #     save_title = str(sys.argv[12])
    #
    #     generate_2D_data(var1_name, var1_lower_bound, var1_upper_bound, var1_range_length, var1_scale,
    #                     var2_name, var2_lower_bound, var2_upper_bound, var2_range_length, var2_scale,
    #                     save_location, save_title)
    # else:
    #     print("Error for command line args:\n ")
    #     print(sys.argv)

    generate_spec_data()
    # print(eval_model({'chemo_strength':0.15,
    #                   'r_Ts':1250}))
