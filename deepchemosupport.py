import pandas as pd

class ModelParameters:

    #Starting cell populations
    Ts_0 = 10**7   #Tumor cells
    Tr_0 = 0 #10**5 # Resistant Tumor cells
    E_0 = 100   #Effector cells
    R_0 = 100   #Treg cells
    M_0 = 10**6.5  #Memory cells
    N_0 = 10**12   #Naive cells




    #Variables for the dT/dt equation:

    r_Ts = 1000      #the Tumor growth rate
    r_Tr = 100
    k_0s = 1.0     #The killing efficiency of E cells on suppressing tumor growth
    k_0r = 0.25
    k_0 = 0.15
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
        self.initial_pops = [self.Ts_0, self.Tr_0, self.E_0, self.R_0, self.M_0, self.N_0]
        self.resistance = resistance
        self.model_UID = UID

    



class Schedule:

    def __init__(self):

        self.intervals = [14, 14*2,14*3, 14*4]
