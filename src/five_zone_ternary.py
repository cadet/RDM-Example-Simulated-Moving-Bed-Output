# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# (Example_3_SMB)=
# # Ternary separation

# %% [markdown]
# Case study III examines the separation of the three components 2’-deoxycytidine `A`, 2’-deoxyguanosine `B` and 2’-deoxyadenosine `C`. A standard five-zone SMB system with two extract ports and without partial-feeding or partial-closing is assumed. <br> This experiment was originally simulated and published by Mun et al. in "Improving performance of a five-zone simulated moving bed chromatography for ternary separation by simultaneous use of partial-feeding and partial-closing of the product port in charge of collecting the intermediate-affinity solute molecules" (Sungyong Mun, Journal of Chromatography A 2011;  1218(44):8060-8074) <br> https://doi.org/10.1016/j.chroma.2011.09.015.
#

# %% [markdown]
# ```{figure} ./figures/case_study3.jpg
# :width: 600px
# <div style="text-align: center">
# (Fig. 4, He et al.) Integrated five-zone SMB scheme with two extract ports
# <div>

# %% [markdown]
# (Ternary separation with a five zone system is particularly effective when the target component is present in much higher abundance than the more strongly retained component.)

# %% [markdown]
# ## Differences in Parameters He's paper / He's Matlab code / Mun's paper
#    
#     1.)  feed concentrations match in He paper text and He table cF = 1.0 g/cm^-3 
#     -> Do not match original Mun paper: Mun Table 1 cF = 1.0 g/L -> would be 1.0e3 g/cm^3
#     feed concentrations are too small
#     Matlab code:     
#     concentrationFeed 	= [1.0, 1.0, 1.0];    % g/cm^3 [concentration_compA, concentration_compB]  => wrong unit 
#     opt.molMass         = [227.217, 267.24, 251.24192]; % The molar mass of each components
#     
#     He's paper: feed.c = [4.41e3, 3.75e3, 3.98e3]  # c_in [mol / m^3]  
#     Muns paper c + Matlab code M = [4.401079144606257, 3.74195479718605, 3.9802275034357324]  [mol / m^3]
#
#     2.) column.particle_porosity = 1.0e-5  # ε_p [-] 
#     Matlab code: opt.porosityParticle    = 0.00000001;   % e_p very small to ensure e_t = e_c  => would be 1.0e-8
#
#     3.) column.film_diffusion = component_system.n_comp * [1.6e4]  # k_f [m / s]  
#     Matlab code: opt.filmDiffusion             = [5.0e-5, 2.5e-5, 5.0e-5];  % K_f
#
#     4.) column.pore_diffusion = component_system.n_comp * [5e-5]  # D_p [m² / s]
#     Matlab code: opt.diffusionParticle         = [1.6e4, 1.6e4, 1.6e4];  % D_p
#
# # volumetric flow rates, time checked
#
#     % Purities archieved at operation point a for components A,B,C: 99.62  71.16	96.78
#     flowRate.recycle    = 2.9230e-7;      % m^3/s  == QI, not QIV or QV!
#
# # Parameters missing
#     opt.tolIter         = 1e-4;  % tolerance of the SMB stopping criterion
#     opt.nMaxIter        = 1000;  % the maximum iteration step in SMB
#     opt.nThreads        = 8;     % threads of CPU, up to your computer
#
#     opt.timePoints      = 1000;  % the observed time-points
#     -> For Feed concentration setup
#     Feed.time = linspace(0, opt.switch, opt.timePoints);  # 1000 evaluated points during 1 switching time
#     Feed.concentration = zeros(length(Feed.time), opt.nComponents);
#
#     
#     opt.Purity_extract1_limit   = 0.95;  % used for constructing constraints
#     opt.Purity_extract2_limit   = 0.65;  % used for constructing constraints
#     opt.Purity_raffinate_limit  = 0.99;  % used for constructing constraints
#     opt.Penalty_factor          = 10;    % penalty factor in penalty function
#
#     opt.yLim            = max(concentrationFeed ./ opt.molMass) * 1.1; % the magnitude for plotting
#
#
#  #   Capable of placing a CSTR or DPFR apparatues before and after the calculated column
#
#     Continuous Stirred Tank Reactor
#     opt.enable_CSTR = false;
#     opt.CSTR_length = 0.01;
#
#     Dispersive Plug Flow Reactor
#     opt.enable_DPFR = false;
#
#     opt.DPFR_length = 0.0066;
#     opt.DPFR_nCells = 50;
#
#     opt.DPFR_velocity   = 0.00315;
#     opt.DPFR_dispersion = 2.5e-20;
#

# %%
# Transformation of datatypes for adsorption/desorption rates with mass transfer
import numpy as np
mass_transfer = [1, 0.5, 0.1]
arr1 = np.divide([3.15, 7.40, 23.0], mass_transfer)# Henry_1 = 3.15; Henry_2 = 7.40, Henry_3 = 23.0, second ternary separation system (Mun et al.) ;  list(np.divide([1, 1, 1], mass_transfer)) 
arr2 = np.divide([1, 1, 1], mass_transfer) 
k_a = [float(x) for x in arr1]
k_d = [float(x) for x in arr2]

# %%
from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Linear
from CADETProcess.processModel import Inlet, Outlet, GeneralRateModel
import numpy as np

# Component System
component_system = ComponentSystem(['A', 'B', 'C'])
#Components at extraction from:     
#    opt.comp_raf_ID  = 1; % the target component withdrawn from the raffinate ports
#    opt.comp_ext1_ID = 3; % the target component withdrawn from the extract_1 ports
#    opt.comp_ext2_ID = 2; % the target component withdrawn from the extract_2 ports

# Binding Model
binding_model = Linear(component_system)
binding_model.is_kinetic = True
mass_transfer = [1, 0.5, 0.1]
#binding_model.adsorption_rate = [3.15, 7.40, 23.0]# Henry_1 = 3.15; Henry_2 = 7.40, Henry_3 = 23.0, second ternary separation system (Mun et al.) ; 
binding_model.adsorption_rate = k_a
#binding_model.desorption_rate = [1, 1, 1] # k_kin = Mass-transfer coefficient (ap km ), 1/s, kd = 1/k_kin
binding_model.desorption_rate = k_d

# Column
column = GeneralRateModel(component_system, name='column')
column.binding_model = binding_model
column.length = 0.150 # L [m]
column.diameter = 1.0e-2  # d [m]
column.bed_porosity = 0.80  # ε_c [-]
column.particle_porosity = 1.0e-8  # ε_p [-] 
#Matlab code: opt.porosityParticle    = 0.00000001;   % e_p very small to ensure e_t = e_c  => would be 1.0e-8
column.particle_radius = 1.50e-5  # r_p [m]

#column.film_diffusion = component_system.n_comp * [1.6e4]  # k_f [m / s]  
column.film_diffusion = [5.0e-5, 5.0e-5, 5.0e-5]
#Matlab code: opt.filmDiffusion             = [5.0e-5, 2.5e-5, 5.0e-5];  % K_f   Componente B fits better with 5.0e-5

#column.pore_diffusion = component_system.n_comp * [5e-5]  # D_p [m² / s]
column.pore_diffusion = [1.6e4, 1.6e4, 1.6e4]
#Matlab code: opt.diffusionParticle         = [1.6e4, 1.6e4, 1.6e4];  % D_p


column.axial_dispersion = 3.81e-10  # D_ax [m² / s  #Matlab code: 3.8148e-10;
column.discretization.npar = 1  # N_r
column.discretization.ncol = 40  # N_z

column.solution_recorder.write_solution_bulk = True

eluent = Inlet(component_system, name='eluent')
eluent.c = [0, 0, 0]  # c_in_D [mol / m^3]
eluent.flow_rate = 2.34e-7  # Q_D [m^3 / s]  Matlab code: flowRate.desorbent  = 2.3412e-7;      % m^3/s



feed = Inlet(component_system, name='feed')
#feed.c = [4.41e3, 3.75e3, 3.98e3]  # c_in [mol / m^3]  
feed.c = [4.40, 3.74, 3.98]  # c_in [mol / m^3] 
# Muns paper c + Matlab code M = [4.401079144606257, 3.74195479718605, 3.9802275034357324]
feed.flow_rate = 1.67e-8  # Q_F [m^3 / s]  Matlab code: flowRate.feed       = 1.6667e-8;      % m^3/s

# %% [markdown]
# All zones are connected to each other in series `SerialZone`
# ```
# zone_I -> extract_1 + zone_II
# Q_I = Q_E + Q_II 
#
# zone_II -> extract_2 + zone_III
# Q_II = Q_E2 + Q_III 
#
# zoneIV -> raffinate + zone_V
# Q_IV = Q_R + Q_V 
#
# w_e1 = Q_E / Q_I = 0.6438
# w_e2 = Q_E2 / Q_II = 0.4419
# w_r = Q_R / Q_IV = 0.224

# %%
#Q_R = 1.68e-8  # Operating point a, Table 3 Mun et al.
#Q_E = 1.88e-7  # Operating point a, Table 3
#Q_E2 = 4.64e-8  
#Q_I = 2.92e-7
#Q_II = 1.05e-7
#Q_IV =7.50e-8

#Q_E / Q_I 
#Q_E2 / Q_II 
#Q_R / Q_IV 

# %%
extract_1 = Outlet(component_system, name = 'extract_1')
extract_2 = Outlet(component_system, name = 'extract_2')
raffinate = Outlet(component_system, name = 'raffinate')
from CADETProcess.modelBuilder import SerialZone

zone_I = SerialZone(component_system, 'zone_I', n_columns=1, valve_dead_volume=1e-9)
zone_II = SerialZone(component_system, 'zone_II', n_columns=1, valve_dead_volume = 1e-9)
zone_III = SerialZone(component_system, 'zone_III', n_columns=1, valve_dead_volume = 1e-9)
zone_IV = SerialZone(component_system, 'zone_IV', n_columns=1, valve_dead_volume = 1e-9)
zone_V = SerialZone(component_system, 'zone_V', n_columns=1, valve_dead_volume = 1e-9)

from CADETProcess.modelBuilder import CarouselBuilder

builder = CarouselBuilder(component_system, 'smb')
builder.valve_dead_volume = 1e-9
builder.column = column
builder.add_unit(eluent)
builder.add_unit(feed)
builder.add_unit(extract_1)
builder.add_unit(extract_2)
builder.add_unit(raffinate)

builder.add_unit(zone_I)
builder.add_unit(zone_II)
builder.add_unit(zone_III)
builder.add_unit(zone_IV)
builder.add_unit(zone_V)

builder.add_connection(eluent, zone_I)

builder.add_connection(zone_I, extract_1)
builder.add_connection(zone_I, zone_II)
w_e1 = 0.6438
builder.set_output_state(zone_I, [w_e1, 1 - w_e1])

builder.add_connection(zone_II, extract_2)
builder.add_connection(zone_II, zone_III)
w_e2 = 0.4419
builder.set_output_state(zone_II, [w_e2, 1 - w_e2])

builder.add_connection(zone_III, zone_IV)
builder.add_connection(feed, zone_IV)

builder.add_connection(zone_IV, raffinate)
builder.add_connection(zone_IV, zone_V)
w_r = 0.224
builder.set_output_state(zone_IV, [w_r, 1 - w_r])

builder.add_connection(zone_V, zone_I)

builder.switch_time = 264  # Fig 8 - standard mode at 200.99 steps n = step number/switching number

process = builder.build_process()

# %%
from CADETProcess.simulator import Cadet
process_simulator = Cadet()
#process_simulator.evaluate_stationarity = True
process_simulator.n_cycles = 41  #200.99 steps = 200.99 switch times -> 40.198
process_simulator.use_dll = True
#process_simulator.timeout = 15*60
#simulate first 8 switch times (1 iteration), conc bis 1mol

process_simulator.time_integrator_parameters.abstol = 1e-10
process_simulator.time_integrator_parameters.reltol = 1e-6  # Not in Matlab code!, not in Klatt paper 
process_simulator.time_integrator_parameters.init_step_size = 1e-14
process_simulator.time_integrator_parameters.max_step_size = 5e6

#process_simulator.time_resolution = builder.switch_time / 1000  # default value is 1 second.  Matlab code: Feed.time = linspace(0, opt.switch, opt.timePoints);

simulation_results = process_simulator.simulate(process)
cycle = 41
_ = simulation_results.solution.raffinate.inlet.plot()
_ = simulation_results.solution.extract_1.inlet.plot()
_ = simulation_results.solution.extract_2.inlet.plot()
_ = simulation_results.solution.raffinate.inlet.plot(start = cycle * 0, end = 200.99 * builder.switch_time)

# %%
import numpy as np
raff_mM = simulation_results.solution.raffinate.inlet.solution
ext1_mM = simulation_results.solution.extract_1.inlet.solution
ext2_mM = simulation_results.solution.extract_2.inlet.solution
t = simulation_results.time_complete

# Transformation from mM to g/L
molar_mass = [227.217, 267.24, 251.24192]  # molar mass of each component
raff = np.multiply(raff_mM, molar_mass) * 1e-3
ext_1 = np.multiply(ext1_mM, molar_mass) * 1e-3
ext_2 = np.multiply(ext2_mM, molar_mass) * 1e-3

# %%
import matplotlib.pyplot as plt

fig, axs = plt.subplots(2, 2, figsize=(20, 17))
ax1 = axs[0, 0]  # Raff
ax2 = axs[0, 1]  # Ext2
ax3 = axs[1, 0]  # Ext1

ax1.plot(t / builder.switch_time, raff)
ax1.set_title("Raffinate")
ax1.set_xlabel("Switches")
ax1.set_ylabel("c [g / L]")
ax1.set_ylim(0, 1)
#ax1.set_xlim(0, 8)

ax2.plot(t / builder.switch_time, ext_2)
ax2.set_title("Extract 2")
ax2.set_xlabel("Switches")
ax2.set_ylabel("c [g / L]")
ax2.set_ylim(0, 0.8)
#ax2.set_xlim(0, 8)

ax3.plot(t / builder.switch_time, ext_1)
ax3.set_title("Extract 1")
ax3.set_xlabel("Switches")
ax3.set_ylabel("c [g / L]")
ax3.set_ylim(0, 0.8)
#ax3.set_xlim(0, 40)

# %%
#class CarouselSolutionBulk(SolutionBase): 
#from CADETProcess.modelBuilder.carouselBuilder import CarouselSolutionBulk
#axial_conc = CarouselSolutionBulk(builder, simulation_results)
#axial_conc.component_system
#axial_conc.solution
#axial_conc.axial_coordinates
#axial_conc.time
#simulation_results.solution
#ax1, _ = axial_conc.plot_at_time(t=200.01 * builder.switch_time) #, ax = _plot_solution_1D(ax = ) )
#ax2, _ = axial_conc.plot_at_time(t=200.99 * builder.switch_time)

# CADET output in mM = mol / m^3
# mM -> g / L 

# component C correct 


# %%
#dir(ax1[0])
#fig = ax1[0].figure  # Hole die Figure vom Axes-Objekt
    # (In vielen Umgebungen funktioniert das, aber nicht überall)
   # Funktioniert global für aktive Figuren


#fig = ax1[0].figure
#fig.show()
#ax = ax1[0]
#ax.set_title("zone_I")
#plt.show()


#print(ax.get_title())        # → "zone_I"
#print(ax.get_xlim())         # → z. B. (0.0, 1.0)
#print(dir(ax1[0]))               # → alle Attribute und Methoden
#ax1[0].plot()
#print(ax1[0].lines)  # zeigt alle Linien-Objekte


# %% [markdown]
# ```{figure} ./figures/ternary.png
# :width: 800px
# <div style="text-align: center">
# (Fig. 8, Mun et al.) internal concentration profiles of the five-zone SMBs, (a) Standard (at 200.01 steps), (b) Standard (at 200.99 steps), Blue line: component A, red line: component B, green line: component C.
# <div>

# %% [markdown]
# ```{figure} ./figures/ternary_separation_Mun.png
# :width: 800px
# <div style="text-align: center">
# (Fig. 8, Mun et al.) internal concentration profiles of the five-zone SMBs, Operation mode: Standard mode (at 200.99 steps)
# <div>
