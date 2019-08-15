# Derek Park
# September, 2016
import numpy
import scipy
from scipy import integrate
import matplotlib
from matplotlib import pyplot as plt
import pandas
import os
import datetime
import multiprocessing
import ImmuneModel_vaccine as ImmuneModel
import time


# The purpose is to generate specific plots of the patient


class MyException(Exception):

    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def eval_model(vars):

    mod = ImmuneModel.ImmuneModel_vaccine('specific')
    mod.set_params(vars)
    var1 = vars.keys()[0]
    var1_val = vars[var1]
    var2 = vars.keys()[1]
    var2_val = vars[var2]

    simulation_results = mod.simulate(10) # Simulate and get all of the results

    #simulation_results.to_csv('RECIST_outcomes/Individual_runs/' + var1 + '_=_' + str(var1_val) + '_' + var2 + '_=_' + str(var2_val) + '.csv')
    final_tumor_size = simulation_results['T'].tolist()[-1]

    percent_change_in_size = ((final_tumor_size) / 10 ** 8.0 - 1.0) * 100

    vars['Result'] = percent_change_in_size

    maximum_effector_cell_count = simulation_results['E'].max()
    vars['max_E_cell_count'] = maximum_effector_cell_count
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
    vars['final_size'] = final_tumor_size

    return vars


def single_plot(input_vals, output=False):

    var1_name, var1_value, var2_name, var2_value = input_vals

    vars = {}
    vars[var1_name] = var1_value
    vars[var2_name] = var2_value
    vars['vaccine_strength'] = 0.5
    vars['vaccine_start_offset'] = 0
    vars['vaccine_dosing_interval'] = 14
    mod = ImmuneModel.ImmuneModel_vaccine('specific')
    print vars
    mod.set_params(vars)


    simulation_results = mod.simulate(10) # Simulate and get all of the results

    # Generate the series of alpha values during the course of therapy
    start_of_therapy_time = mod.get_start_of_therapy_time()
    time = simulation_results.index
    alpha_series = []
    for t in time:
        alpha_series.append(mod.current_alpha(mod.alpha, t))

    print len(alpha_series)
    print time

    #simulation_results.to_csv(var1_name + '_=_' + str(var1_value) + '_' + var2_name + '_=_' + str(var2_value) + '.csv')
    return simulation_results
    if output:
        matplotlib.rcParams.update({'font.size': 30})
        plt.figure(1)
        #plt.title(var1_name + '_=_' + str(var1_value) + '_' + var2_name + '_=_' + str(var2_value))
        #plt.title('x')
        ax = plt.subplot(111)
        plt.yscale('log')
        plt.xlabel('Days')
        plt.ylabel('Cellular Compartment Size')
        plt.ylim([1, 10**12.5])
        time = simulation_results.index
        ax.plot(time, simulation_results['T'], label='T', linewidth=5)
        ax.plot(time, simulation_results['E'], label='E', linewidth=5)
        ax.plot(time, simulation_results['R'], label='R', linewidth=5)
        ax.plot(time, simulation_results['M'], label='M', linewidth=5)
        ax.plot(time, simulation_results['N'], label='N', linewidth=5)

        # Shrink current axis's height by 10% on the bottom

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height*0.2,
                         box.width, box.height * 0.8])

        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10),
                  fancybox=True, shadow=True, ncol=5)

        plt.figure(2)

        plt.title("Alpha during the course of therapy")
        plt.ylabel('Alpha value')
        plt.xlabel('Days')
        ax = plt.subplot(111)
        ax.plot(time, alpha_series, label='alpha', linewidth=5)
        plt.show(block=True)

option = 0



# Plot a single treatment
if option == 0:
    res = single_plot(['M_0', 10**5.0, 'chemo_strength', 0.6], output=True)
    res.to_csv('M0_low.csv')

# Results for different vaccine strengths versus chemo
if option == 0.5:

    M_0 = 1000000.0

    chemo_range = numpy.linspace(0.2, .9, num=100)
    vaccine_strengths = [0, 1.0, 10.0, 100.0, 1000.0]
    results_series = []

    k = 1
    for vac_strength in vaccine_strengths:
        vac_series = []

        i = 1
        for chemo in chemo_range:
            print str(k) + "." + str(i)
            vars = {}
            vars['M_0'] = M_0
            vars['chemo_strength'] = chemo
            vars['vaccine_start_offset'] = 0
            vars['vaccine_dosing_interval'] = 14
            vars['vaccine_half_life'] = 3.0
            vars['vaccine_strength'] = vac_strength
            vars['chemo_effect_on_immune'] = 1.0
            final_change = eval_model(vars)['final_size']
            vac_series.append(final_change)
            i += 1

        results_series.append(vac_series)
        k += 1

    if True:
        # Subtract to create relative differences:
        baseline = results_series[0]
        i = 0
        for list in results_series:
            results_series[i] = numpy.subtract(baseline, list)
            i += 1
    matplotlib.rcParams.update({'font.size': 40})
    plt.figure(1)
    # plt.title(var1_name + '_=_' + str(var1_value) + '_' + var2_name + '_=_' + str(var2_value))
    plt.title('')
    ax = plt.subplot(111)
    plt.xlabel('Chemotherapy strength')
    plt.ylabel('Improvement in Tumor Response (cells)')
    #plt.yscale('log')
    #plt.ylim([1,10**9])

    i = 1
    while i < len(vaccine_strengths):
        series_label = str(vaccine_strengths[i])
        ax.plot(chemo_range, results_series[i], label=series_label, linewidth=8)
        i += 1
    # Shrink current axis's height by 10% on the bottom

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2,
                     box.width, box.height * 0.8])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10),
              fancybox=True, shadow=True, ncol=5)

    plt.show(block=True)

    transposed_results = []

    i = 0
    while i < len(results_series[0]):
        curr_line = []
        k = 0

        while k < len(results_series):
            curr_line.append(results_series[k][i])
            k += 1

        transposed_results.append(curr_line)
        i +=1

    indices = chemo_range
    results_frame = pandas.DataFrame(transposed_results, columns=vaccine_strengths, index=indices)
    results_frame.to_csv('vaccine_improvement_vs_chemo.csv', index_label='Chemo_strength')

# Plot the vaccine improvement results
if option == 0.6:
        data = pandas.read_csv('vaccine_improvement_vs_chemo.csv')
        matplotlib.rcParams.update({'font.size': 40})
        plt.figure(1)
        # plt.title(var1_name + '_=_' + str(var1_value) + '_' + var2_name + '_=_' + str(var2_value))
        plt.title('')
        ax = plt.subplot(111)
        plt.xlabel('Chemotherapy strength')
        #plt.gca().axes.get_yaxis().set_ticks([])
        plt.ylabel('Tumor Size Reduction (cells)')
        # plt.yscale('log')

        plt.plot(data.iloc[:, 0], data.iloc[:, 2], label='1',linewidth=8)
        plt.plot(data.iloc[:, 0], data.iloc[:, 3], label='10',linewidth=8)
        plt.plot(data.iloc[:, 0], data.iloc[:, 4], label='100',linewidth=8)
        plt.plot(data.iloc[:, 0], data.iloc[:, 5], label='1000',linewidth=8)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.2,
                         box.width, box.height * 0.8])
        # Put a legend below current axis
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20),
                  fancybox=True, shadow=True, ncol=5)
        plt.show(block=True)

# RECIST outcomes for different vaccine strengths vs chemo
if option == 0.75:

    M_0 = 1000000.0

    chemo_range = numpy.linspace(0.2, .9, num=100)
    vaccine_strengths = [0, 1.0, 10.0, 100.0, 1000.0]
    results_series = []

    k = 1
    for vac_strength in vaccine_strengths:
        vac_series = []

        i = 1
        for chemo in chemo_range:
            print str(k) + "." + str(i)
            vars = {}
            vars['M_0'] = M_0
            vars['chemo_strength'] = chemo
            vars['vaccine_start_offset'] = 0
            vars['vaccine_dosing_interval'] = 14
            vars['vaccine_half_life'] = 3.0
            vars['vaccine_strength'] = vac_strength
            vars['chemo_effect_on_immune'] = 1.0
            final_change = eval_model(vars)['RECIST'] - k*0.1
            vac_series.append(final_change)
            i += 1

        results_series.append(vac_series)
        k += 1


    matplotlib.rcParams.update({'font.size': 40})
    plt.figure(1)
    # plt.title(var1_name + '_=_' + str(var1_value) + '_' + var2_name + '_=_' + str(var2_value))
    plt.title('')
    ax = plt.subplot(111)
    plt.xlabel('Chemotherapy strength')
    plt.gca().axes.get_yaxis().set_ticks([])
    # plt.ylabel('RECIST response')
    # plt.yscale('log')
    # plt.ylim([1,10**9])

    i = 0
    while i < len(vaccine_strengths):
        series_label = str(vaccine_strengths[i])
        ax.plot(chemo_range, results_series[i], label=series_label, linewidth=8)
        i += 1
    # Shrink current axis's height by 10% on the bottom

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2,
                     box.width, box.height * 0.8])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.20),
              fancybox=True, shadow=True, ncol=5)

    plt.show(block=True)

    transposed_results = []

    i = 0
    while i < len(results_series[0]):
        curr_line = []
        k = 0

        while k < len(results_series):
            curr_line.append(results_series[k][i])
            k += 1

        transposed_results.append(curr_line)
        i += 1

    indices = chemo_range
    results_frame = pandas.DataFrame(transposed_results, columns=vaccine_strengths, index=indices)
    results_frame.to_csv('vaccine_improvement_vs_chemo_RECIST.csv', index_label='Chemo_strength')

# Plot the vaccine vs RECIST results
if option == 0.8:

    data = pandas.read_csv('vaccine_improvement_vs_chemo_RECIST.csv')
    matplotlib.rcParams.update({'font.size': 40})
    plt.figure(1)
    # plt.title(var1_name + '_=_' + str(var1_value) + '_' + var2_name + '_=_' + str(var2_value))
    plt.title('')
    ax = plt.subplot(111)
    plt.xlabel('Chemotherapy strength')
    plt.gca().axes.get_yaxis().set_ticks([])
    # plt.ylabel('RECIST response')
    # plt.yscale('log')
    plt.ylim([-3, 1.5])
    plt.axhline(y=2.0, color='black', linestyle='dashed', linewidth=3)
    plt.axhline(y=1.0, color='black', linestyle='dashed', linewidth=3)
    plt.axhline(y=0, color='black', linestyle='dashed', linewidth=3)
    plt.axhline(y= -1.0, color='black', linestyle='dashed', linewidth=3)
    plt.axhline(y=-2.0, color='black', linestyle='dashed', linewidth=3)
    plt.plot(data.iloc[:, 0], data.iloc[:, 1], color='black', marker='s', linewidth=0, markersize=15)
    plt.plot(data.iloc[:,0], data.iloc[:,2:],  marker='s', linewidth=0, markersize=15)
    plt.show(block=True)


if option == 1:

    M_0 = 1000000.0

    chemo_range = numpy.linspace(0.2,.9,num=100)
    #chemo_range = [0.25, 0.6, 0.9]

    #immune_diffs = [1.0, 0.0]
    immune_diffs = numpy.linspace(0.1, 1.5, num=100)
    results_series = []

    k = 1
    for diff in immune_diffs:
        vac_series = []

        i = 1
        for chemo in chemo_range:
            print str(k) + "." + str(i)
            vars = {}
            vars['M_0'] = M_0
            vars['chemo_strength'] = chemo
            vars['vaccine_start_offset'] = 0
            vars['vaccine_dosing_interval'] = 14
            vars['vaccine_half_life'] = 3.0
            vars['vaccine_strength'] = 0.0
            vars['chemo_effect_on_immune'] = diff
            final_change = eval_model(vars)['final_size']
            vac_series.append(final_change)
            i += 1
        results_series.append(vac_series)
        k += 1

    if False:
        # Subtract to create relative differences:
        baseline = results_series[0]
        i = 0
        for list in results_series:

            results_series[i] = numpy.subtract(baseline, list)
            i += 1
    matplotlib.rcParams.update({'font.size': 22})
    plt.figure(1)
    # plt.title(var1_name + '_=_' + str(var1_value) + '_' + var2_name + '_=_' + str(var2_value))
    plt.title('')
    ax = plt.subplot(111)
    plt.xlabel('Chemotherapy strength')
    plt.ylabel('Final tumor size')
    #plt.yscale('log')
    #plt.ylim([-10**8,10**8])

    i = 0
    while i < len(immune_diffs):
        series_label = "Rel. Immune Effect = " + str(immune_diffs[i])
        ax.plot(chemo_range, results_series[i], label=series_label, linewidth=3)
        i += 1
    # Shrink current axis's height by 10% on the bottom

    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2,
                     box.width, box.height * 0.8])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10),
              fancybox=True, shadow=True, ncol=5)


    plt.show(block=True)

    indices = map(str, immune_diffs)
    indices = immune_diffs
    results_frame = pandas.DataFrame(results_series, columns= map(str, chemo_range), index=indices)
    results_frame.to_csv('final_tumor_sizes_rel_immune_effect_changing_chemo.csv', index_label='Chemo_strength')

if option == 2:

    df = pandas.read_csv('final_tumor_sizes_rel_immune_effect_changing_chemo.csv', index_col=0, header=0)
    #df = df.T[0]
    #df.to_csv('final_tumor_sizes_rel_immune_effect_changing_chemo=0.25.csv')
    print df.index

    #print df
    matplotlib.rcParams.update({'font.size': 40})
    plt.figure(1)
    # plt.title(var1_name + '_=_' + str(var1_value) + '_' + var2_name + '_=_' + str(var2_value))
    plt.title('')
    ax = plt.subplot(111)
    plt.xlabel('Relative chemotherapeutic Effect')
    plt.ylabel('Final Tumor Size (cells)')
    plt.axhline(y=10**8, color='black', linestyle='dashed', linewidth='3')
    plt.yscale('log')

    start_col = 0
    #labels = df.columns
    labels = ['0.25', '0.6', '0.9']
    #labels = ['0.25']
    print
    for i in range(len(labels)):
        if i >= 0:
            plt.plot(df.index, df.iloc[:,i+start_col], label=labels[i], linewidth=8)




    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2,
                     box.width, box.height * 0.8])
    # Put a legend below current axis
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
              fancybox=True, shadow=True, ncol=5)

    plt.show(block=True)


#Run a simple MC simulation of treatment with different attenuation effects
if option == 3:
    #os.chdir('//home//parkds//PhD_projects//Chemo-Immune_model//Vaccine_Trials//')
    num_patients = 500
    chemo_strength = [0.9, 0.8, 0.7, 0.6]
    #chemo_strength = [0.6]
    attenuation_mean = 1.0
    attenuation_sigma_percent = 25


    for chemo in chemo_strength:
        percent_changes = []
        final_sizes = []

        for i in range(num_patients):

            curr_attenuation = numpy.random.normal(loc=attenuation_mean, scale=(attenuation_mean * attenuation_sigma_percent/100.0))

            vars = {}
            vars['M_0'] = 10**6
            vars['chemo_strength'] = chemo
            vars['vaccine_start_offset'] = 0
            vars['vaccine_dosing_interval'] = 14
            vars['vaccine_half_life'] = 3.0
            vars['vaccine_strength'] = 0.0
            vars['chemo_effect_on_immune'] = curr_attenuation
            results = eval_model(vars)

            final_size = results['final_size']
            final_sizes.append(final_size)
            final_percent_change = (final_size - 10**8)/10**8 *100.0
            percent_changes.append(final_percent_change)
            print str(chemo) + '.' + str(i)

        dframe = pandas.DataFrame([percent_changes, final_sizes], index = ['Percent changes', 'final size']).T
        file_name = 'MC_' + str(num_patients) + '_patients_chemo=' + str(chemo) + '_atten_mean=' + str(attenuation_mean)
        file_name += '_atten_sigma=' + str(attenuation_sigma_percent) + '.csv'
        dframe.to_csv(file_name)

#Make the waterfall plot
if option == 4:
    f_name = 'MC_500_patients_chemo=0.9_atten_mean=0.8_atten_sigma=25'

    dframe = pandas.read_csv(f_name + '.csv')

    percent_changes = dframe['Percent changes'].sort_values(ascending=False).tolist()
    print percent_changes

    matplotlib.rcParams.update({'font.size': 30})
    plt.figure(1)
    # plt.title(var1_name + '_=_' + str(var1_value) + '_' + var2_name + '_=_' + str(var2_value))
    #plt.title('x')
    ax = plt.subplot(111)
    #plt.yscale('log')
    #plt.xlabel('Days')
    plt.ylabel('Percent change in Tumor Size')
    plt.ylim([-100, 100])
    print len(percent_changes)
    plt.bar(numpy.arange(500), percent_changes)
    #plt.savefig(f_name + '.png')
    plt.show(block=True)
