# Derek Park
# July 2019
import numpy
import scipy
from scipy import integrate
import matplotlib
from matplotlib import pyplot as plt
import pandas
import os
import datetime

class MyException(Exception):

    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class ImmuneModel_vaccine:

    #Starting cell populations
    T_0 = 10**7   #Tumor cells
    E_0 = 100   #Effector cells
    R_0 = 100   #Treg cells
    M_0 = 10**6.5  #Memory cells
    N_0 = 10**12 - 10**6.5   #Naive cells



    #Variables for the dT/dt equation:

    r_T = 1000      #the Tumor growth rate
    k_0 = 0.1      #The killing efficiency of E cells on suppressing tumor growth
    B = 0.75        #The scaling factor for the effectiveness of R cells on suppressing E-mediated tumor death as the R/E ratio increases
    T_trans = 10**6 # The size at which the tumor transitions from exponential to power law growth
    m = 1/2.    #the power of the power-law growth rate of the tumor
    P = 3.0     # A smoothing term b/w exponential and power law growth

    #Global variables
    T_off = 5   #The time at which CTL expansion turns of and that compartment begins to contract
    E_max = 10**12  #The maximum number of M + N cell



    #Variables for the dE/dt equaton:
    alpha = 1   #The conversion efficiency of memory cells into effector cells
    rho = 0    #The proportional death rate of E cell
    delta_E = 0.13  #The death rate of the E cells during the contraction phase
    c = 0.01        #The proportional increase in death rate of E cells during contraction due to R cells
    gamma = 100  #The M cell multiplier since multiple E cells come from 1 M cell



    #Variables for the dR/dt equaton:
    sigma = 0.01    #The proportional increase in R due to T (i.e. R recruitment by T)
    delta_R = 0.1  #The proportionl death rate of R



    #Variables for the dM/dt equaton
    r_M = 0.01      #The growth rate of the M cells
    omega = 0.01     #The fraction of lost E cells that are converted back into M cells

    #Variables for the dN/dt equaton
    r_N = 0.1      #The growth rate for the N cell


    chemo_strength = 0.75


    treatment_interval = 14  #The interval between treatments

    number_cycles = 5 # The number of chemotherapy cycles

    model_UID = ""

    # These are the parameters to be varied:
    param_1_label = "k0"
    param_1_value = 0.0
    param_2_label = "alpha"
    param_2_value = 0.0

    start_time = 100

    resistance=True


    vaccine_dosing_interval = 14.0
    vaccine_strength = 0
    vaccine_start_offset = 3.0
    vaccine_half_life = 3.0

    chemo_effect_on_immune = 1.0

    time_start_therapy = -1

    def __init__(self, UID, resistance=True):
        self.resistance = resistance
        self.model_UID = UID

        # Starting cell populations
        self.T_0 = 10 ** 7  # Tumor cells
        self.E_0 = 100  # Effector cells
        self.R_0 = 100  # Treg cells
        self.M_0 = 10 ** 6  # Memory cells
        self.N_0 = 10 ** 12 - 10 ** 6  # Naive cells




        #Variables for the dT/dt equation:

        self.r_T = 1000      #the Tumor growth rate
        self.k_0 = 1.0      #The killing efficiency of E cells on suppressing tumor growth
        self.B = 0.75        #The scaling factor for the effectiveness of R cells on suppressing E-mediated tumor death as the R/E ratio increases
        self.T_trans = 10**6 # The size at which the tumor transitions from exponential to power law growth
        self.m = 1/2.    #the power of the power-law growth rate of the tumor
        self.P = 3.0     # A smoothing term b/w exponential and power law growth

        #Global variables
        self.T_off = 5   #The time at which CTL expansion turns of and that compartment begins to contract
        self.E_max = 10**12  #The maximum number of M + N cell



        #Variables for the dE/dt equaton:
        self.alpha = 1   #The conversion efficiency of memory cells into effector cells
        self.rho = 0    #The proportional death rate of E cell
        self.delta_E = 0.13  #The death rate of the E cells during the contraction phase
        self.c = 0.01        #The proportional increase in death rate of E cells during contraction due to R cells
        self.gamma = 100  #The M cell multiplier since multiple E cells come from 1 M cell



        #Variables for the dR/dt equaton:
        self.sigma = 0.01    #The proportional increase in R due to T (i.e. R recruitment by T)
        self.delta_R = 0.1  #The proportionl death rate of R



        #Variables f    or the dM/dt equaton
        self.r_M = 0.01      #The growth rate of the M cells
        self.omega = 0.01     #The fraction of lost E cells that are converted back into M cells

        #Variables for the dN/dt equaton
        self.r_N = 0.1      #The growth rate for the N cell


        self.chemo_strength = 0.25


        self.treatment_interval = 14  #The interval between treatments

        self.number_cycles = 5 # The number of chemotherapy cycles

        self.vaccine_dosing_interval = 14.0
        self.vaccine_strength = 0.5
        self.vaccine_start_offset = 3.0
        self.vaccine_half_life = 3.0

        self.chemo_effect_on_immune = 1.0

    # Return the time that therapy starts
    def get_start_of_therapy_time(self):

        return self.time_start_therapy

    # Set the parameters of this model
    def set_params(self, dict):

        for key in dict:

            if hasattr(self, key):
                setattr(self, key, dict[key])
            else:
                print "Wrong variable: " + key


    # Return the current numerical parameters of this model
    def return_params(self, separate=False):

        x = {key:value for key, value in self.__dict__.items() if not key.startswith('__') and not callable(key) and key != 'model_UID'}

        banned_parameters = ['T_0', 'T_off', 'number_cycles']

        x = {}
        for key,value in self.__dict__.items() :
            #print key
            if key not in banned_parameters  and not key.startswith('__') and not callable(key) and key != 'model_UID':
                x[key] = value



        if not separate:
            return x
        else:
            return x.keys(), x.values()

    # Return a dynamically changing alpha based on a simulated immunotherapy dose
    def current_alpha(self, base_alpha, time):


        dose_times = [self.time_start_therapy + self.vaccine_start_offset]
        dose_times.append(dose_times[-1] + self.vaccine_dosing_interval)
        dose_times.append(dose_times[-1] + self.vaccine_dosing_interval)

        #The half life of the dose in days:
        half_life = self.vaccine_half_life

        curr_alpha = base_alpha

        # Calculate the immunostimulatory effect of each dose
        for indiv_dose in dose_times:
            time_since_dose = time - indiv_dose

            # i.e. the dose should be given already
            if time_since_dose >= 0:
                current_dose_strength = base_alpha * self.vaccine_strength * 0.5 **(time_since_dose/half_life)

                curr_alpha += current_dose_strength




        return curr_alpha

    # The Equations that govern expansion of the system.
    def expansion_equations(self, t0, y):

        T = y[0]
        E = y[1]
        R = y[2]
        M = y[3]
        N = y[4]

        Z = (M + N)/self.E_max

        T_growth_denom = ((1/(self.T_trans**(self.m - 1)*self.r_T))**self.P + (T**(1-self.m)/self.r_T)**self.P)**(1/3.0)

        dT = T/T_growth_denom - self.k_0 * (T * E)/(T + E) * (1 - self.B * R/(R + E))
        #dT = T/((1/(0.001*self.r_T))^3 + (T^(1-self.m)/self.r_T)^3)^(1/3) - self.k_0 * (T * E)/(T + E) * (1 - self.B * R/(R + E))
        dE = (1 - Z) * self.gamma * self.alpha * (T * M) / (T + M) - self.rho * E

        dR = self.sigma * T - self.delta_R * R
        dM = -1 * ((1 - Z) * self.alpha * T * M / (T + M) ) + self.r_M * M * (1 - Z)
        dN = self.r_N * N * (1 - Z)


        return [dT, dE, dR, dM, dN]

    # The equations that govern contraction of the immune system after attenuation
    def contraction_equations(self, t0, y):
        #print "y is: " + str(y)
        #print "t0 is: " + str(t0)
        T = y[0]
        E = y[1]
        R = y[2]
        M = y[3]
        N = y[4]

        Z = (M + N)/self.E_max

        T_growth_denom = ((1/(self.T_trans**(self.m - 1)*self.r_T))**self.P + (T**(1-self.m)/self.r_T)**self.P)**(1/3.0)

        dT = T/T_growth_denom - self.k_0 * (T * E)/(T + E) * (1 - self.B * R/(R + E))
        #dT = T/((1/(0.001*self.r_T))**3. + (T**(1-self.m)/self.r_T)**3.)**(1/3.) - self.k_0 * (T * E)/(T + E) * (1 - self.B * R/(R + E))
        dE = -1 *  self.delta_E * E * (1.0 + self.c * R/(R + E))
        # In case of the E < 1, then we just don't have any change since this is effectively 0. If we do not implement this, there are artefacts due to the numeric solver driving E into negative values
        if E < 1.0:
            dE = 0
        dR = self.sigma * T - self.delta_R * R
        dM = self.delta_E * E * self.omega + self.r_M * M * (1 - Z)
        dN = self.r_N * N * (1 - Z)

        return [dT, dE, dR, dM, dN]


    # Experimental try using the heavyside:

    def heavyside_growth(self, t0, y, current_Toff, simulation_start):

        curr_runtime = datetime.datetime.now() - simulation_start

        if curr_runtime > datetime.timedelta(minutes=3):
            print 'Integration timeout'
            raise Exception('integration timeout')

        T = y[0]
        E = y[1]
        R = y[2]
        M = y[3]
        N = y[4]


        curr_alpha = self.current_alpha(self.alpha, t0)

        heavyside_Toff_minus_t = 0
        if t0 < current_Toff:
            heavyside_Toff_minus_t = 1

        heavyside_t_minus_Toff = 0
        if t0 > current_Toff:
            heavyside_t_minus_Toff = 1

        Z = (M + N)/self.E_max

        T_growth_denom = ((1/(self.T_trans**(self.m - 1)*self.r_T))**self.P + (T**(1-self.m)/self.r_T)**self.P)**(1/3.0)

        dT = T/T_growth_denom - self.k_0 * (T * E)/(T + E) * (1 - self.B * R/(R + E))
        #dT = T/((1/(0.001*self.r_T))^3 + (T^(1-self.m)/self.r_T)^3)^(1/3) - self.k_0 * (T * E)/(T + E) * (1 - self.B * R/(R + E))

        dE = heavyside_Toff_minus_t * ((1 - Z) * self.gamma * curr_alpha * (T * M) / (T + M) - self.rho * E) - heavyside_t_minus_Toff * (self.delta_E * E * (1.0 + self.c * R/(R + E)))

        dR = self.sigma * T - self.delta_R * R
        dM = heavyside_Toff_minus_t * (-1 * ((1 - Z) * curr_alpha * T * M / (T + M) )) + heavyside_t_minus_Toff * (self.delta_E * E * self.omega)  + self.r_M * M * (1 - Z)
        dN = self.r_N * N * (1 - Z)



        #Exception conditions

        if T < 10**-1 and dT < 0:
            dT = 0 #If it is less than .1 tumor cells and it looks like it is giong to decrease, then make no change. However, it can increase!

        if E < 10**-1 and dE < 0:
            dE = 0

        if R < 10**-1 and dR < 0:
            dR = 0

        if M < 10**-1 and dM < 0:
            dM = 0
        if N < 10**-1 and dN < 0:
            dN = 0

       # print [T, E, R, M, N]
        return [dT, dE, dR, dM, dN]

    # See if the integration has been going too long
    def going_too_long(self):

        current_time = datetime.datetime.now()

        if current_time - self.start_time > datetime.timedelta(minutes=0.5):
            return -1
        else:
            return None

    # A function to describe the system going through multiple cycles of chemotherapy
    def cycle(self, time_start, initial_pops, simulation_start_time):

        T0 = initial_pops[0]
        E0 = initial_pops[1]
        R0 = initial_pops[2]
        M0 = initial_pops[3]
        N0 = initial_pops[4]

        def integration_timeout(t,y):
            curr_time = datetime.datetime.now()
            if (curr_time - simulation_start_time) > datetime.timedelta(minutes=0.5):
                print 'Integration Timeout'
                return -1
            return 0



        total_time = [time_start]
        total_growth = [initial_pops]

        heavyside = True
        # if not heavyside:
        #     #Simulate immune expansion
        #
        #     solver = scipy.integrate.ode(self.expansion_equations)
        #     # solver.set_integrator('dopri5', nsteps=500)
        #     solver.set_integrator('vode', nsteps=50000, method='bdf')
        #     solver.set_initial_value(initial_pops, time_start)
        #     solver.set_solout(integration_timeout)
        #
        #
        #
        #
        #     dt = 0.001
        #
        #     while solver.successful() and solver.t < (time_start + self.T_off):
        #
        #         solver.integrate(solver.t + dt)
        #         total_time.append(solver.t)
        #         total_growth.append(solver.y)
        #
        #
        #     # Simulate immune contraction
        #
        #     solver = scipy.integrate.ode(self.contraction_equations)
        #     # solver.set_integrator('dopri5', nsteps=500)
        #     solver.set_integrator('vode', nsteps=50000, method='bdf')
        #     solver.set_initial_value(total_growth[-1], total_time[-1])
        #     solver.set_solout(integration_timeout)
        #
        #     while solver.successful() and solver.t < (time_start + self.treatment_interval):
        #
        #         solver.integrate(solver.t + dt)
        #         #print "y values: " + str(solver.y)
        #         total_time.append(solver.t)
        #         total_growth.append(solver.y)
        # else:
        try:
            self.start_time = datetime.datetime.now()
            # Simulate the whole cycle using the heavyside equations
            solver = scipy.integrate.ode(self.heavyside_growth)
            # solver.set_integrator('dopri5', nsteps=500)
            solver.set_integrator('vode', nsteps=5000000, method='bdf')
            solver.set_f_params(time_start + self.T_off, simulation_start_time)
            solver.set_initial_value(initial_pops, time_start)


            dt = 0.001

            while solver.successful() and solver.t < (time_start + self.treatment_interval):

                solver.integrate(solver.t + dt)
                #print "y values: " + str(solver.y)
                total_time.append(solver.t)
                total_growth.append(solver.y)
        except:
            print 'Returning None'
            return None, None

        return total_time, total_growth



    #Simulate the leadup growth

    def leadup(self, time_start, initial_pops):

        T0 = initial_pops[0]
        E0 = initial_pops[1]
        R0 = initial_pops[2]
        M0 = initial_pops[3]
        N0 = initial_pops[4]

        #print time_series[0:end_index].tolist()

        solver = scipy.integrate.ode(self.contraction_equations)
        # solver.set_integrator('dopri5', nsteps=500)
        solver.set_integrator('zvode', nsteps=5000, method='bdf')
        solver.set_initial_value(initial_pops, time_start)


        total_time = [time_start]
        total_growth = [initial_pops]

        dt = 0.01

        do_integration = True
        self.start_time = datetime.datetime.now()
        while do_integration:

            if solver.y[0] > 10**8:
                do_integration = False

            solver.integrate(solver.t + dt)
            total_time.append(solver.t)
            total_growth.append(solver.y)

            if datetime.datetime.now()-self.start_time > datetime.timedelta(seconds=3.0):
                raise Exception('Tumor did not grow at the start')

        return total_time, total_growth

    # Simulate therapy
    def simulate(self, num_cycles):

        T0 = self.T_0
        E0 = self.E_0
        R0 = self.R_0
        M0 = self.M_0
        N0 = self.N_0

        #print "Feedback is: " + str((M0 + N0)/self.E_max)

        # Simulate the leadup in growth
        try:
            total_time, total_growth = self.leadup(0, [T0, E0, R0, M0, N0])
        except:
            return None
        curr_cycle = 0

        time_simulation_start = datetime.datetime.now()
        self.time_start_therapy = total_time[-1]
        while curr_cycle < num_cycles:
            #print curr_cycle

            curr_start_time = total_time[-1]
            #print curr_start_time
            # Treat it with chemo
            curr_start_growth = total_growth[-1] * (1 - self.chemo_strength)

            if self.resistance:
                new_pops = []
                new_T = total_growth[-1][0]* (1 - (self.chemo_strength - (self.chemo_strength*0.25 * curr_cycle)/num_cycles))
                new_E = total_growth[-1][1]* (1 - self.chemo_strength * self.chemo_effect_on_immune)
                new_R = total_growth[-1][2]* (1 - self.chemo_strength * self.chemo_effect_on_immune)
                new_M = total_growth[-1][3]* (1 - self.chemo_strength * self.chemo_effect_on_immune)
                new_N = total_growth[-1][4]* (1 - self.chemo_strength * self.chemo_effect_on_immune)

                new_pops = [new_T, new_E, new_R, new_M, new_N]
                curr_start_growth = new_pops

            cycle_time, cycle_growth =  self.cycle(curr_start_time, curr_start_growth, time_simulation_start)

            if not cycle_growth:

                return None

            total_time = total_time + cycle_time
            total_growth = numpy.vstack((total_growth,cycle_growth))
            #print "Tot growth: " + str(total_growth[-1])
            curr_cycle += 1


        model_file_name = "test/" + self.model_UID + "_" + str(self.param_1_label) + "_" + str(self.param_1_value) + "_vs_" + str(self.param_2_label) + "_" + str(self.param_2_value)


        

        return dataframe


mod = ImmuneModel_vaccine('test')
print mod.return_params(separate=True)[0]
