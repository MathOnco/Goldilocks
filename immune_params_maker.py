import numpy as np

# variable_combos = [['chemo_strength', 'alpha'],
#                     ['chemo_strength', 'rho'],
#                     ['chemo_strength', 'c'],
#                     ['chemo_strength', 'gamma'],
#                     ['chemo_strength', 'r_M'],
#                     ['chemo_strength', 'omega'],
#                     ['chemo_strength', 'delta_R'],]
variable_combos = [['chemo_strength', 'r_Ts'],
                    ['chemo_strength', 'k_0'],
                    ['chemo_strength', 'rho']]


# variable_bounds = [[[0, 0.9],[-1, 3]], #Alpha
#                    [[0, 0.9],[0., 1.]], #rho
#                    [[0, 0.9],[-3, 1.]], #c
#                    [[0, 0.9],[1., 3.]], # gamma
#                    [[0, 0.9],[-3., -1]], # r_M
#                    [[0, 0.9],[-3., -1]], #omega
#                    [[0, 0.9],[-2, 1]]] # delta_R
variable_bounds = [[[0.05, 0.9],[500., 1500.]], #k_0
                    [[0.05, 0.9],[0.5,1.5]], #r_T
                  [[0.05, 0.9],[0, 1.]]] #rho

# variable_scales = [['linear', 'log'],
#                     ['linear', 'linear'],
#                     ['linear', 'log'],
#                     ['linear', 'log'],
#                     ['linear', 'log'],
#                     ['linear', 'log'],
#                     ['linear', 'log'],]
variable_scales = [['linear', 'linear'],
                    ['linear', 'linear'],
                    ['linear', 'linear']]


# variable_lengths = [[150, 150],
#                     [150, 150],
#                     [150, 150],
#                     [150, 150],
#                     [150, 150],
#                     [150, 150],
#                     [150, 150]]
variable_lengths = [[150, 150],
                    [150, 150],
                    [150, 150]]



param_file_string = ""
for combo, bounds, scales, lengths in zip(variable_combos, variable_bounds, variable_scales, variable_lengths):

    var1_name = combo[0]

    var1_lower_bound = bounds[0][0]
    var1_upper_bound = bounds[0][1]
    var1_range_length = lengths[0]
    var1_scale = scales[0]

    var2_name = combo[1]

    var2_lower_bound = bounds[1][0]
    var2_upper_bound = bounds[1][1]

    var2_scale = scales[1]

    var2_breakpoints = np.linspace(var2_lower_bound, var2_upper_bound, 6, endpoint = False)

    for k in range(len(var2_breakpoints)):
        this_var2_lower = var2_breakpoints[k]
        if k + 1 < len(var2_breakpoints):
            this_var2_upper = var2_breakpoints[k+1]
        else:
            this_var2_upper = var2_upper_bound

        save_location = '.'
        save_title = 'fig3_params-%s_vs_%s_%i.csv' % (var1_name, var2_name, k)
        param_file_string += "%s %f %f %i %s %s %f %f %i %s %s %s\n" %(var1_name, var1_lower_bound, var1_upper_bound, var1_range_length, var1_scale,
                                                                var2_name, this_var2_lower, this_var2_upper, 25, var2_scale,
                                                                save_location, save_title)

k = 0
for line in param_file_string.split('\n'):
    print(line)
    with open('fig3_param_files/fig3_params-%i.txt' %(k), 'w') as f:
        f.write(line)
        f.close()
    k += 1
