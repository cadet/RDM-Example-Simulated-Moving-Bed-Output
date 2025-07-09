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
# (Example_SMB)=
# # Simulated Moving Bed

# %% [markdown]
# The following example is a reproduction of part of the research results published in "Efficient numerical simulation of simulated moving bed chromatography with a single-column solver" (Qiao-Le He, Samuel Leweke, Eric von Lieres, Computers & Chemical Engineering 2018; 111:183-198. doi:10.1016/j.compchemeng.2017.12.022.) <br>
# https://www.sciencedirect.com/science/article/pii/S0098135417304520
#

# %% [markdown]
# The first of the five case studies depicted in the paper evaluates the separation of the two components fructose `A` and glucose `B` in a four-zone simulated moving bed (SMB) system with eight columns. The binding behavior follows a linear isotherm. 
# This experiment was originally simulated and published by Klatt et al. in "Model-based control of a simulated moving bed chromatographic process for the separation of fructose and glucose" (Karsten-Ulrich Klatt, Felix Hanisch, Guido Dünnebier, Journal of Process Control 2002; 12(2):203-219. https://doi.org/10.1016/S0959-1524(01)00005-1.) <br>
# https://www.sciencedirect.com/science/article/abs/pii/S0959152401000051
#
#
# Continuous processes generally have many benefits over batch processes, like higher efficiency and throughput. However, a truly continuous chromatography process is not practically feasable. This so called true moving bed chromatography would entail the solid phase of the chromatography column (the bed) moving in the opposite direction of the mobile phase. This would induce the local separation of components in the feed solution based on their column binding behaviour. Those components could then by retrieved by different outlet streams located upstream and downstream of the feed inlet. 
#
# ```{figure} ./figures/true_MB.png
# :width: 400px
# <div style="text-align: center">
#
#  [Link to Youtube Video](https://www.youtube.com/watch?v=xhhJxb48tgc)
# <div>
#
#

# %% [markdown]
# Simulated Moving Bed Chromatography is a way to approach such a continuous process in practice. This is realized by having multiple chromatography columns connected to each other in a caroussel. By periodically switching the inlet and outlet valves connected to the columns ("column switching") the movement of the solid phase can be mimicked. 
#
# A general SMB system contains four different external units connected to the columns:
# 1. Feed (Inlet): Component mixture
# 2. Raffinate (Outlet): Faster eluting component (elutes before feed plug flow)
# 3. Extract (Outlet): Slower eluting component (elutes after feed plug flow), interacts more strongly with the column solid phase<br>
# There can be multiple extract or raffinate outlets in a SMB system depending on the number of components to be separated.
# 4. Desorbant / Eluent / Solvent (Inlet): Solution to elute extract from column before raffinate plug flow enters the column again to prevent mixing of the separated components
#
# The position of each column relative to these external units determine their specific “Zone” in the SMB system. Based on thier external units, every zone has a different function in the seperation process and a different flow rate within the columns. The most basic SMB system for binary separation is made up of four zones with one chromatography column in each zone. 
#

# %% [markdown]
# ```{figure} ./figures/four_zone.jpg
# :width: 600px
# <div style="text-align: center">
# (Fig. 1, He et al.) Schematic of four-zone SMB chromatography for binary separations. Column positions are periodically switched in opposite direction of liquid flow.

# %%
# Four-zone binary SMB

# %% [markdown]
# ```{figure} ./figures/case_study1_practical_setup.jpg
# :width: 600px
# <div style="text-align: center">
# (Fig. 5, He et al.) Four-zone SMB schemes with eight columns indicating positions of the associated hold-up volumes
# <div>

# %% [markdown]
# As seen in Fig. 5, hold up-volumes between all columns and external units exist and should ideally be considered in thier effect on the SMB elution. (triangle theory) They generally increase retention time and dispersion. (The hold-up volume on either side of a column, i.e., tubing and frits, can be described by a CSTR that is moved through the network together with that column. )(four categories: (1) tubing between multi-port valve and column inlet plus frit before packed bed, (2) frit after packed bed plus tubing between column outlet and multi-port valve, (3) tubing between injection point and multi-port valve, and (4) tubing between multi-port valve and detector. Each of these categories can be modeled as one or more PFR, CSTR and/or DPFR in series. )
# -> CADET-SMB allows to consider hold-up volumes in the column network. This is demonstrated by introducing CSTR models, Eq. (10), in case study I as illustrated by Fig. 5. The residence time, τ
# CSTR, is varied between 0s, 5s and 10s. Fig. 14 shows the impact of these hold-up volumes on the column states in CSS.

# %% [markdown]
# # Differences in He's Matlab code / He's paper / original case study in [Klatt's paper](https://www.sciencedirect.com/science/article/pii/S0959152401000051?ref=pdf_download&fr=RR-2&rr=94b07706292368ec#TBL1):
#
# ## Parameters from Matlab code (getParameters_binary_case2.m) not in CarouselBuilder Example:
#
#         % The parameter setting for simulator
#         opt.tolIter         = 1e-4;
#
#         % The parameter settting for the SMB
#         opt.Purity_limit    = [0.99, 0.99];
#         opt.Penalty_factor  = 10;
#
#         % opt.compTargID = [2, 1];
#         opt.structID    = [2 2 2 2];
#         opt.diffusionParticleSurface  = [0.0 0.0];  # probably automatically 0 if not defined
#         opt.enableDebug = true;
#
#         % The parameter setting for simulator
#         opt.nThreads        = 4;(/8)
#         opt.timePoints      = 1000         !!
#
#         opt.yLim            = max(concentrationFeed ./ opt.molMass);
#         
#         (Viscosity in Klatt paper, not in He paper)
#
# ## Interstitial velocities missing (Calculation already checked for all flow rates Q): 
#         Interstitial velocities calculated explicitly in Matlab code, done automatically by CarouselBuilder?
#         
#         % Interstitial velocity = flow_rate / (across_area * opt.porosityColumn)
#         interstVelocity.recycle   = flowRate.recycle / (crossArea*opt.porosityColumn);      % m/s
#         interstVelocity.feed      = flowRate.feed / (crossArea*opt.porosityColumn);         % m/s
#         interstVelocity.raffinate = flowRate.raffinate / (crossArea*opt.porosityColumn);    % m/s
#         interstVelocity.desorbent = flowRate.desorbent / (crossArea*opt.porosityColumn);    % m/s
#         interstVelocity.extract   = flowRate.extract / (crossArea*opt.porosityColumn);      % m/s
#
#
#         process_simulator.time_integrator_parameters.reltol = 1e-6  # Not in Matlab code, not in Klatt paper, but in He paper
#         

# %% [markdown]
# ##  Matlab code capable of placing a CSTR or DPFR apparatues before and after the calculated column:
#
# %   Continuous Stirred Tank Reactor
#     opt.enable_CSTR = false;
#     opt.CSTR_length = 0.01;
#
# %   Dispersive Plug Flow Reactor
#     opt.enable_DPFR = false;
#
#     opt.DPFR_length = 0.0066;
#     opt.DPFR_nCells = 50;
#
#     opt.DPFR_velocity   = 0.00315;
#     opt.DPFR_dispersion = 2.5e-20;
#
#
# ## Differences in parameters:
#     
#     1.)  He paper: "The inlet concentrations are converted from 0.55 g/cm^3 assuming that fructose and glucose have the same molar mass of 180 g/mol." 
#     => Would be: (0.55/180)*1e6 = 3055.555 = 3.06e3 mol/m^3
#     
#     Text does not match table in paper
#     
#     # He paper Table 1: 2.78e3 
#     # original Klatt paper:  "cF = 0.5 g/cm3" => would also be 2.78e3
#
#     Matlab code: 
#         concentrationFeed 	= [0.5, 0.5];   % g/m^3 [concentration_compA, concentration_compB]   => wrong unit, is actually g/cm^3
#         opt.molMass         = [180.16, 180.16];
#         opt.yLim            = max(concentrationFeed ./ opt.molMass);
#
#         Feed.time = linspace(0, opt.switch, opt.timePoints);
#         Feed.concentration = zeros(length(Feed.time), opt.nComponents);  
#     
#     2.) Film diffusion and particle diffusion exchanged in He's paper and Matlab code
#     Matlab code: opt.filmDiffusion             = [5e-5 5e-5];
#     Matlab code: opt.diffusionParticle         = [1.6e4 1.6e4];
#
#
#     3.) Recycle flow rate nomenclature:
#     Matlab code: Recycle flow rate QI: flowRate.recycle    = 0.1395e-6;      % m^3/s 
#     Klatt/He paper: Recycle flow rate QIV = 0.0981 cm 3 /s = 9.81e-8 m^3/s 
#     
#     4.) Order of Henry coefficients exchanged:
#     # Matlab code: opt.KA = [0.28 0.54]; % [comp_A, comp_B], A for raffinate, B for extract
#     opt.comp_raf_ID = 1;  % the target component withdrawn from the raffinate ports
#     opt.comp_ext_ID = 2;  % the target component withdrawn from the extract ports    

# %% [markdown]
# The first of the five case studies depicted in the paper evaluates the separation of the two components glucose `A` and fructose `B` in a four-zone simulated moving bed (SMB) system with eight columns. The binding behavior follows a linear isotherm. 
# This experiment was originally simulated and published by Klatt et al. in "Model-based control of a simulated moving bed chromatographic process for the separation of fructose and glucose" (Karsten-Ulrich Klatt, Felix Hanisch, Guido Dünnebier, Journal of Process Control 2002; 12(2):203-219. https://doi.org/10.1016/S0959-1524(01)00005-1.) <br>
# https://www.sciencedirect.com/science/article/abs/pii/S0959152401000051
#
# To simulate a SMB process, first the physical properties of the columns and the Inlet are defined. The mass transfer within the column is characterized by the equilibrium-dispersive model (EDM) which can be derived from the `GeneralRateModel` by defining the spatial discretization `column.discretization.npar` as 1. In the finite volume method, only one radial cell is assumed. The axial column dimension `column.discretization.ncol` is set to 40 axial cells. ((All numerical values are taken from Table 1.(4. Case Studies))). The Henry coefficient can be assumed to equal the equilibrium constant under ideal, linear conditions. 
# The axial concentration of every column can later be visualized at a specific time by plotting the `column.solution_recorder.write_solution_bulk` concentration. 
#

# %%
from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Linear
from CADETProcess.processModel import Inlet, Outlet, GeneralRateModel

# Component System
component_system = ComponentSystem(['A', 'B'])
# Matlab code: opt.comp_raf_ID = 1;  % the target component withdrawn from the raffinate ports
# Matlab code: opt.comp_ext_ID = 2;  % the target component withdrawn from the extract ports

# Binding Model
binding_model = Linear(component_system)
binding_model.is_kinetic = False
#binding_model.adsorption_rate = [0.54, 0.28]  # Henry_1 = 	0.54 = fructose; Henry_2 = 0.28 = glucose
binding_model.adsorption_rate = [0.28, 0.54]
# Matlab code: opt.KA = [0.28 0.54]; % [comp_A, comp_B], A for raffinate = glucose, B for extract = fructose
binding_model.desorption_rate = [1, 1]
 

# Column
column = GeneralRateModel(component_system, name='column')
column.binding_model = binding_model

column.length = 0.536  # L [m]
column.diameter = 2.6e-2  # d [m]
column.bed_porosity = 0.38  # ε_c [-] 
# Matlab: porosityColumn

column.particle_porosity = 1.0e-5  # ε_p [-] 
column.particle_radius = 1.63e-3  # r_p [m]
#Matlab code:  opt.particleRadius      = 0.325e-2 /2 = 0,001625; 
#column.film_diffusion = component_system.n_comp * [1.6e4]  # k_f [m / s]
column.film_diffusion = component_system.n_comp * [5e-5]
#Matlab code: opt.filmDiffusion             = [5e-5 5e-5];

#column.pore_diffusion = component_system.n_comp * [5e-5]  # D_p [m² / s]
column.pore_diffusion = component_system.n_comp * [1.6e4]
#Matlab code: opt.diffusionParticle         = [1.6e4 1.6e4];

column.axial_dispersion = 3.81e-6  # D_ax [m² / s]
#Matlab code: opt.dispersionColumn          = ones(1,opt.nZone) .* 3.8148e-6;
column.discretization.npar = 1  # N_r
column.discretization.ncol = 40  # N_z

column.solution_recorder.write_solution_bulk = True

eluent = Inlet(component_system, name='eluent') #Name in paper = "desorbent"
eluent.c = [0, 0]  # c_in_D [mol / m^3]
#Matlab code: Desorbent.concentration = zeros(length(Feed.time), opt.nComponents);
eluent.flow_rate = 4.14e-8  # Q_D [m^3 / s] 
#

feed = Inlet(component_system, name='feed')
feed.c = [2.78e3, 2.78e3]  # c_in [mol / m^3] => He Matlab, Klatt, NOT HE PAPER
#He paper: "The inlet concentrations are converted from 0.55 g/m^3 assuming that fructose and glucose have the same molar mass of 180 g/mol." => [3052.84, 3052.84])

#Matlab code: Feed.time = linspace(0, opt.switch, opt.timePoints)
#Matlab code: Feed.concentration = zeros(length(Feed.time), opt.nComponents)
feed.flow_rate = 2.0e-8  # Q_F [m^3 / s]

# %% [markdown]
# The unit system of the SMB is implemented using the `CarouselBuilder` from CADETProcess. Four zones with two columns in each zone and the connections to their respective external units are implemented as seen in Fig. 5. For more information please refer: [here](https://cadet-process.readthedocs.io/en/stable/user_guide/tools/carousel_builder.html#). (Not using SMB builder because `n_columns` = 2) <br>
# The percentile of the volume flow that leaves `zone_I` for `extract`(`w_e`) or `zone_II` (`1-w_e`) can be deducted by examining the differences in the volumetric flow rate as all columns have the same cross section `A`(Table 1.). The same can be done for the percentile of the volume flow that leaves `zone_III` for `raffinate`(`w_r`) and `zone_IV` (`1-w_r`) <br>
# The continuity equation for laminar flow is assumed:
# `A_1 * v_1 = A_2 * v_2`
# ```
# zone_I -> extract + zone_II
# w_e = Q_E / Q_I 
# w_e = (3.48e-8) / (1.4e-7) = 0.249  (= 0.24857142857142855)
#
# zoneIII -> raffinate + zone_IV
# w_r = Q_R / Q_III 
# w_r = (2.66e-8 ) / (1.25e-7) = 0.213 (= 0.2128000000000000)
# ``` 

# %%
extract = Outlet(component_system, name='extract')
raffinate = Outlet(component_system, name='raffinate')
from CADETProcess.modelBuilder import SerialZone

zone_I = SerialZone(component_system, 'zone_I', n_columns = 2, valve_dead_volume=1e-9)
zone_II = SerialZone(component_system, 'zone_II', n_columns = 2, valve_dead_volume=1e-9)
zone_III = SerialZone(component_system, 'zone_III', n_columns = 2, valve_dead_volume=1e-9)
zone_IV = SerialZone(component_system, 'zone_IV', n_columns = 2, valve_dead_volume=1e-9)

from CADETProcess.modelBuilder import CarouselBuilder

builder = CarouselBuilder(component_system, 'smb')
builder.valve_dead_volume = 1e-9
builder.column = column
builder.add_unit(eluent)
builder.add_unit(feed)

builder.add_unit(extract)
builder.add_unit(raffinate)

builder.add_unit(zone_I)
builder.add_unit(zone_II)
builder.add_unit(zone_III)
builder.add_unit(zone_IV)

builder.add_connection(eluent, zone_I)

builder.add_connection(zone_I, extract)
builder.add_connection(zone_I, zone_II)
w_e = 0.249  
builder.set_output_state(zone_I, [w_e, 1-w_e])

builder.add_connection(zone_II, zone_III)

builder.add_connection(feed, zone_III)

builder.add_connection(zone_III, raffinate)
builder.add_connection(zone_III, zone_IV)
w_r = 0.213
builder.set_output_state(zone_III, [w_r, 1-w_r])

builder.add_connection(zone_IV, zone_I)

builder.switch_time = 1552

process = builder.build_process()

# %% [markdown]
# ((CSS - cyclic steady state: dynamic trajectory is repeated after every switch => all columns have same state after specific time laps ))

# %%
from CADETProcess.simulator import Cadet
process_simulator = Cadet()
#process_simulator.evaluate_stationarity = True
process_simulator.n_cycles = 13
process_simulator.use_dll = True
#Upper panels show the first 40 switching times (five iterations), 

process_simulator.time_integrator_parameters.abstol = 1e-10
process_simulator.time_integrator_parameters.reltol = 1e-6  # Not in Matlab code!, not in Klatt paper 
process_simulator.time_integrator_parameters.init_step_size = 1e-14
process_simulator.time_integrator_parameters.max_step_size = 5e6

simulation_results = process_simulator.simulate(process) #509s

# %%
import matplotlib.pyplot as plt
import numpy as np

raff = simulation_results.solution.raffinate.inlet.solution
ext = simulation_results.solution.extract.inlet.solution
t = simulation_results.time_complete

# n = c * V = pi * L * (d/2)**2 
#volume = np.pi * column.length * (column.diameter / 2) ** 2

fig, axs = plt.subplots(2, 2, figsize=(20, 8))
ax1 = axs[1, 0]  # Extract 8 swt
ax2 = axs[1, 1]  # Raffinate 8 swt
ax3 = axs[0, 0]  # Extract 40 swt
ax4 = axs[0, 1]  # Raffinate 40 swt

# Extrakt 8 swt
ax1.plot(t / builder.switch_time, ext)
ax1.set_title("Extract")
ax1.set_xlabel("Switches")
ax1.set_ylabel("c [mM]")
ax1.set_ylim(0, 1500)
ax1.set_xlim(0, 8)


# Raffinat 8 swt
ax2.plot(t / builder.switch_time, raff)
ax2.set_title("Raffinate")
ax2.set_xlabel("Switches")
ax2.set_ylabel("c [mM]")
ax2.set_ylim(0, 1500)
ax2.set_xlim(0, 8)

# Plot links: Extrakt 40swt
ax3.plot(t / builder.switch_time, ext)
ax3.set_title("Extract")
ax3.set_xlabel("Switches")
ax3.set_ylabel("c [mM]")
ax3.set_ylim(0, 3000)
ax3.set_xlim(0, 40)

# Plot rechts: Raffinat 40 swt
ax4.plot(t / builder.switch_time, raff)
ax4.set_title("Raffinate")
ax4.set_xlabel("Switches")
ax4.set_ylabel("c [mM]")
ax4.set_ylim(0, 3000)
ax4.set_xlim(0, 40)

plt.suptitle("Concentration profiles at extract and raffinate ports for 40 and 8 switch times", fontsize = 18)
plt.tight_layout()
plt.show()

# %% [markdown]
# Component A glucose is recovered in the raffinate and fructose in the extract outlet. The establishment of the cyclic steady state (CSS) can be seen in the upper graphs. At around the 30th switching time, the concentration profiles start to not change noticibly anymore.  
#
# To visualize the axial concentration of every column at the CSS, the **bulk** concentrations are plotted at the end of the period before the 104th switch, the start of the 13th cycle. 

# %%
from CADETProcess.modelBuilder.carouselBuilder import CarouselSolutionBulk
axial_conc = CarouselSolutionBulk(builder, simulation_results)
axial_conc.component_system
axial_conc.solution
axial_conc.axial_coordinates
axial_conc.time
simulation_results.solution
axial_conc.plot_at_time(t = 104 * builder.switch_time - 1)
axial_conc.plot_at_time(t = 104 * builder.switch_time)
#He paper: CSS of case I computed by STD-FPI. The combined axial concentration profiles in all columns are displayed at multiples of the switching time. In this example, a total of 104 switches (= 13 cycles) was required to fall below an error of ɛ
# Klatt paper: Fig. 4. Axial profile (CSS, end of period) of the SMB fructose/glucose separation during nominal operation and for perturbed values of the isotherm parameters.

# t = 48*switchtime = 6 cycles 
#for t = switching time -> only 1 column switch, have to at least switch once for every column
#8x switching time = 1 cycle
#needs a few cycles to get to CSS => columns are filled completely 
