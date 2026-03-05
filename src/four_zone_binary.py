# ---
# jupyter:
#   jupytext:
#     formats: py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Four-Zone Binary Separation

# %% [markdown]
# This case study depicts the separation of the two components glucose `A` and fructose `B` in a four-zone simulated moving bed (SMB) system with **eight columns**. The binding behavior follows a **linear isotherm**. 
# This experiment was originally simulated and published by Klatt et al. in [Model-based control of a simulated moving bed chromatographic process for the separation of fructose and glucose](https://doi.org/10.1016/S0959-1524(01)00005-1) (Karsten-Ulrich Klatt, Felix Hanisch, Guido Dünnebier, Journal of Process Control 2002; 12(2):203-219.)
#
# Later publications have used this experiment as a case study and presented additional results. A part of the results replicated here were published by He et al. in [Efficient numerical simulation of simulated moving bed chromatography with a single-column solver](https://www.sciencedirect.com/science/article/pii/S0098135417304520) (Qiao-Le He, Samuel Leweke, Eric von Lieres, Computers & Chemical Engineering 2018; 111:183-198. doi:10.1016/j.compchemeng.2017.12.022.)
#

# %% [markdown]
# ```{figure} ./figures/case_study1_practical_setup.jpg
# :width: 600px
# :align: center
#
# Four-zone SMB schemes with eight columns indicating positions of the associated hold-up volumes
# [Fig. 5, He et al.](https://www.sciencedirect.com/science/article/pii/S0098135417304520#fig0004)

# %% [markdown]
# As seen in Fig. 5, **hold up-volumes** between all columns and external units within the SMB system exist and should ideally be considered in thier effect on the elution. They generally increase retention time and dispersion and could be described by **CSTRs** within the [**flow sheet**](https://cadet-process.readthedocs.io/en/v0.10.1/user_guide/process_model/flow_sheet.html). In the following examples hold-up volumes are not considered.

# %% [markdown]
# To simulate a SMB process, first the physical properties of the columns and the Inlet are defined. The mass transfer within the column is characterized by the **equilibrium-dispersive model (EDM)** which can be derived from the `GeneralRateModel` by defining the spatial discretization `column.discretization.npar` as 1. In the finite volume method, only one radial cell is assumed. The axial column dimension `column.discretization.ncol` is set to 40 axial cells. [(4.Case studies, He et al.)](https://www.sciencedirect.com/science/article/pii/S0098135417304520#sec0020)
# <br> 
# The process parameters are based on [Table 1, Klatt et al.](https://www.sciencedirect.com/science/article/pii/S0959152401000051?via%3Dihub#TBL1) and [Table 1 CaseI, He et al.](https://www.sciencedirect.com/science/article/pii/S0098135417304520#tbl0001), [4.1 Case study I, He et al.](https://www.sciencedirect.com/science/article/pii/S0098135417304520#sec0021). To avoid rate-limiting`film_diffusion` and `pore_diffusion`, their numerical values were swapped compared to those in the publication of He et al.. The feed concentration is converted from 0.5 g/m^3 in accordance with Klatt et al., assuming that fructose and glucose have the same molar mass of 180 g/mol. The Henry coefficient can be assumed to equal the equilibrium constant under ideal, linear conditions. 
# The axial concentration of every column can later be visualized at a specific time by plotting the `column.solution_recorder.write_solution_bulk` concentration. 

# %% [markdown]
# ## Setup

# %%
import numpy as np
from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Linear
from CADETProcess.processModel import Inlet, Outlet, GeneralRateModel

# %% [markdown]
# ### Component System

# %%
component_system = ComponentSystem(['A', 'B'])

# %% [markdown]
# ### Binding Model

# %%
binding_model = Linear(component_system)
binding_model.is_kinetic = False
binding_model.adsorption_rate = [0.28, 0.54]
binding_model.desorption_rate = [1, 1]

# %% [markdown]
# ### Unit Operations

# %%
# Transport Model
column = GeneralRateModel(component_system, name='column')
column.binding_model = binding_model
column.length = 0.536  # L [m]
column.diameter = 2.6e-2  # d [m]
column.bed_porosity = 0.38  # ε_c [-] 
column.particle_porosity = 1.0e-5  # ε_p [-] 
column.particle_radius = 1.63e-3  # r_p [m]
column.film_diffusion = component_system.n_comp * [5e-5]  # k_f [m / s]
column.pore_diffusion = component_system.n_comp * [1.6e4]  # D_p [m² / s]
column.axial_dispersion = 3.81e-6  # D_ax [m² / s]
column.discretization.npar = 1  # N_r
column.discretization.ncol = 40  # N_z
column.solution_recorder.write_solution_bulk = True

eluent = Inlet(component_system, name='eluent') #Name in paper = "desorbent"
eluent.c = [0, 0]  # c_in_D [mol / m^3]
eluent.flow_rate = 4.14e-8  # Q_D [m^3 / s] 

feed = Inlet(component_system, name='feed')
feed.c = [2.78e3, 2.78e3]  # c_in [mol / m^3] => He Matlab, Klatt, NOT HE PAPER
feed.flow_rate = 2.0e-8  # Q_F [m^3 / s]

# %% [markdown]
# ### SMB Flow Sheet

# %% [markdown]
# The unit system of the SMB is implemented using the `CarouselBuilder` from CADET-Process. Four zones with two columns in each zone and the connections to their respective external units are defined as seen in [Fig. 1, Klatt et al.](https://www.sciencedirect.com/science/article/pii/S0959152401000051?via%3Dihub#FIGGR1). For more information on the `CarouselBuilder`, please refer: **[here](https://cadet-process.readthedocs.io/en/stable/user_guide/tools/carousel_builder.html#)**. The SMB builder is not used in the examples depicted here because the systems contain multiple `n_columns` per zone. <br>
# The percentile of the volume flow that leaves `zone_I` for `extract`(`w_e`) or `zone_II` (`1-w_e`) can be calculated by dividing the volumetric flow rates of the outlet with the inlet. The same can be done for the percentile of the volume flow that leaves `zone_III` for `raffinate`(`w_r`) and `zone_IV` (`1-w_r`) <br>
#
# ```
# zone_I -> extract + zone_II
# w_e = Q_E / Q_I = 0.249 
#
# zoneIII -> raffinate + zone_IV
# w_r = Q_R / Q_III = 0.213
# ``` 

# %%
extract = Outlet(component_system, name='extract')
raffinate = Outlet(component_system, name='raffinate')
from CADETProcess.modelBuilder import SerialZone

zone_I = SerialZone(component_system, 'zone_I', n_columns = 2, valve_parameters={"valve_dead_volume":1e-9})
zone_II = SerialZone(component_system, 'zone_II', n_columns = 2, valve_parameters={"valve_dead_volume":1e-9})
zone_III = SerialZone(component_system, 'zone_III', n_columns = 2, valve_parameters={"valve_dead_volume":1e-9})
zone_IV = SerialZone(component_system, 'zone_IV', n_columns = 2, valve_parameters={"valve_dead_volume":1e-9})

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
# ### Process
# The SMB process is simulated for 13 `cycles`. During this period, the **column switching** will have been performed 104 times and every column will have been at every possible position within the four zones 12 times. The `time integrator parameters` are set according to [4. Case Studies, He et al.](https://www.sciencedirect.com/science/article/pii/S0098135417304520#sec0020), with the relative tolerance `reltol` set to a sensible value.

# %%
from CADETProcess.simulator import Cadet
process_simulator = Cadet()
#process_simulator.evaluate_stationarity = True
process_simulator.n_cycles = 13
process_simulator.use_dll = True

process_simulator.time_integrator_parameters.abstol = 1e-10
process_simulator.time_integrator_parameters.reltol = 1e-6 
process_simulator.time_integrator_parameters.init_step_size = 1e-14
process_simulator.time_integrator_parameters.max_step_size = 5e6

simulation_results = process_simulator.simulate(process)

# %% [markdown]
# ## Results
# The **extract** and **raffinate** outlet concentrations are plotted for the first 40 and eight switching times respectively, replicating [Fig. 8, He et al.](https://www.sciencedirect.com/science/article/pii/S0098135417304520#fig0009).

# %%
import matplotlib.pyplot as plt

raff = simulation_results.solution.raffinate.inlet.solution
ext = simulation_results.solution.extract.inlet.solution
t = simulation_results.time_complete

fig, axs = plt.subplots(2, 2, figsize=(20, 8))
ax1 = axs[1, 0]  # Extract 8 swt
ax2 = axs[1, 1]  # Raffinate 8 swt
ax3 = axs[0, 0]  # Extract 40 swt
ax4 = axs[0, 1]  # Raffinate 40 swt

# Extract at 8 switching times
ax1.plot(t / builder.switch_time, ext)
ax1.set_title("Extract")
ax1.set_xlabel("Switches")
ax1.set_ylabel("c [mM]")
ax1.set_ylim(0, 1500)
ax1.set_xlim(0, 8)

# Raffinate at 8 switching times
ax2.plot(t / builder.switch_time, raff)
ax2.set_title("Raffinate")
ax2.set_xlabel("Switches")
ax2.set_ylabel("c [mM]")
ax2.set_ylim(0, 1500)
ax2.set_xlim(0, 8)

# Extract at 40 switching times
ax3.plot(t / builder.switch_time, ext)
ax3.set_title("Extract")
ax3.set_xlabel("Switches")
ax3.set_ylabel("c [mM]")
ax3.set_ylim(0, 3000)
ax3.set_xlim(0, 40)

# Raffinate at 40 switching times
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
# Component A (glucose, blue) is recovered in the **raffinate** and fructose (red) in the **extract** outlet. The establishment of the **cyclic steady state (CSS)** can be seen in the upper graphs. At around the 30th switching time, the concentration profiles during a switch period start to not change noticibly anymore.  
#
#
# ### Axial concentrations
# To visualize the axial concentrations of every column at a timepoint during the CSS, the **bulk concentrations** are plotted at the end of the period before the 104th switch [Fig. 2, Klatt et al.](https://www.sciencedirect.com/science/article/pii/S0959152401000051?via%3Dihub#FIGGR2) and at the start of the 13th cycle [Fig. 8, He et al.](https://www.sciencedirect.com/science/article/pii/S0098135417304520#fig0008). The effect of a **column switch** can be observed as the local concentration profile is shifted to the column upstream of the current one by exactly half a zone.

# %%
from CADETProcess.modelBuilder.carouselBuilder import CarouselSolutionBulk
axial_conc = CarouselSolutionBulk(builder, simulation_results)
axial_conc.plot_at_time(t = 104 * builder.switch_time - 1)
axial_conc.plot_at_time(t = 104 * builder.switch_time)
