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
# (Example_3_SMB)=
# # Five-Zone Ternary Separation

# %% [markdown]
# The following case study examines the separation of the three components 2’-deoxycytidine `A`, 2’-deoxyguanosine `B` and 2’-deoxyadenosine `C`. A standard **five-zone SMB** system with **two extract ports** and without partial-feeding or partial-closing is assumed. A `LumpedRateModelWithoutPores` is used as the unit operation model.
#
# This experiment was originally simulated and published in [Improving performance of a five-zone simulated moving bed chromatography for ternary separation by simultaneous use of partial-feeding and partial-closing of the product port in charge of collecting the intermediate-affinity solute molecules](https://doi.org/10.1016/j.chroma.2011.09.015) (Sungyong Mun, Journal of Chromatography A 2011;  1218(44):8060-8074).
#

# %% [markdown]
# ```{figure} ./figures/case_study3.jpg
# :width: 600px
# :align: center
#
# Integrated five-zone SMB scheme with two extract ports
# [Fig. 4, He et al.](https://www.sciencedirect.com/science/article/pii/S0098135417304520#fig0004)

# %% [markdown]
# The following process parameters are taken from [Table 1](https://www.sciencedirect.com/science/article/pii/S002196731101363X?via%3Dihub#tbl0005). The feed concentrations `feed.c` (mol/m^3) are calculated based on the respective `molar_mass` for each component and a feed of 1 g/L. The `axial dispersion` is not explicitly given in the paper. A small numerical value is chosen as to not inhibit the elution, the effect of the axial dispersion is negligible in this case. <br>The mass transfer follows a **linear binding model**, in which the `adsorption_rate` is given as the product of the `mass-transfer coefficient` and the `Henry constant` for each component. The Henry coefficients of the three components are made sure to be sufficiently different to assure their separation [4. Results and discussion](https://www.sciencedirect.com/science/article/pii/S002196731101363X?via%3Dihub#sec0040). The `desorption_rate` is equal to the mass-transfer coefficient of each component [Equation 9b, 9j](https://www.sciencedirect.com/science/article/pii/S002196731101363X?via%3Dihub#sec0035). 
#

# %% [markdown]
# ## Setup

# %%
import numpy as np
from CADETProcess.processModel import ComponentSystem
from CADETProcess.processModel import Linear
from CADETProcess.processModel import Inlet, Outlet, LumpedRateModelWithoutPores

# %% [markdown]
# ### Component System

# %%
component_system = ComponentSystem(['A', 'B', 'C'])

# %% [markdown]
# ### Binding Model

# %%
binding_model = Linear(component_system)
binding_model.is_kinetic = True
binding_model.adsorption_rate = [3.15, 7.40*0.5, 23.0*0.1]  # k_a = apkm * H [1/s]
binding_model.desorption_rate = [1, 0.5, 0.1]  # k_d = apkm [1/s]

# %% [markdown]
# ### Unit Operations

# %%
# Transport Model
column = LumpedRateModelWithoutPores(component_system, name='column')
column.binding_model = binding_model
column.length = 0.150  # L_c [m]
column.diameter = 1.0e-2  # d_c [m]
column.total_porosity = 0.80  # ε [-]
column.axial_dispersion = 1e-7  # E_b [m² / s]
column.discretization.ncol = 40  # N_z
column.solution_recorder.write_solution_bulk = True

eluent = Inlet(component_system, name='eluent')
eluent.c = [0, 0, 0]  # c_in_D [mol / m^3]
eluent.flow_rate = 1.908e-7  # Q_des [m^3 / s] operating point c

feed = Inlet(component_system, name='feed')
feed.c = [4.40, 3.74, 3.98]  # c_in [mol / m^3]
feed.flow_rate = 1.67e-8  # Q_feed [m^3 / s]

# %% [markdown]
# ### SMB Flow Sheet
#
# All zones are connected to each other in series. The zones are set up using the `SerialZone` class. The flow rates of each zone are taken from [Table 3, Point c](https://www.sciencedirect.com/science/article/pii/S002196731101363X?via%3Dihub#tbl0015) and the fractions of the flow going into extract port I `w_e1`, extract port 2 `w_e2` and the raffinate port `w_r` are calculated. Using the `CarouselBuilder` all connections of the zones are set up to resemble the [five-zone SMB scheme](https://www.sciencedirect.com/science/article/pii/S002196731101363X?via%3Dihub#fig0005).
# ```
# zone_I -> extract_1 + zone_II
# Q_I = Q_E1 + Q_II 
#
# zone_II -> extract_2 + zone_III
# Q_II = Q_E2 + Q_III 
#
# zoneIV -> raffinate + zone_V
# Q_IV = Q_R + Q_V 
#
# w_e1 = Q_E1 / Q_I = 0.595
# w_e2 = Q_E2 / Q_II = 0.395
# w_r = Q_R / Q_IV = 0.368

# %%
extract_1 = Outlet(component_system, name = 'extract_1')
extract_2 = Outlet(component_system, name = 'extract_2')
raffinate = Outlet(component_system, name = 'raffinate')
from CADETProcess.modelBuilder import SerialZone

zone_I = SerialZone(component_system, 'zone_I', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})
zone_II = SerialZone(component_system, 'zone_II', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})
zone_III = SerialZone(component_system, 'zone_III', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})
zone_IV = SerialZone(component_system, 'zone_IV', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})
zone_V = SerialZone(component_system, 'zone_V', n_columns=1, valve_parameters={"valve_dead_volume":1e-9})

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
w_e1 = 0.595
builder.set_output_state(zone_I, [w_e1, 1 - w_e1])

builder.add_connection(zone_II, extract_2)
builder.add_connection(zone_II, zone_III)
w_e2 = 0.395
builder.set_output_state(zone_II, [w_e2, 1 - w_e2])

builder.add_connection(zone_III, zone_IV)

builder.add_connection(feed, zone_IV)
builder.add_connection(zone_IV, raffinate)
builder.add_connection(zone_IV, zone_V)
w_r = 0.368
builder.set_output_state(zone_IV, [w_r, 1 - w_r])

builder.add_connection(zone_V, zone_I)

builder.switch_time = 324  

process = builder.build_process()

# %% [markdown]
# ### Process

# %%
from CADETProcess.simulator import Cadet
process_simulator = Cadet()
process_simulator.n_cycles = 41 
process_simulator.use_dll = True
process_simulator.time_integrator_parameters.abstol = 1e-10
process_simulator.time_integrator_parameters.reltol = 1e-6  
process_simulator.time_integrator_parameters.init_step_size = 1e-14
process_simulator.time_integrator_parameters.max_step_size = 5e6

simulation_results = process_simulator.simulate(process)

# %% [markdown]
# ## Results

# %% [markdown]
# The process above simulates the SMB during 41 cycles. There are five switching times for every cycle, one for each column within the system. CADET-Process generates the concentrations in the `simulation_results` in SI-units (mM) by default. To compare the results to the publication, the concentrations are converted to g/L. 

# %%
raff_mM = simulation_results.solution.raffinate.inlet.solution
ext1_mM = simulation_results.solution.extract_1.inlet.solution
ext2_mM = simulation_results.solution.extract_2.inlet.solution
t = simulation_results.time_complete

# Conversion from mM to g/L
molar_mass = [227.22, 267.24, 251.24]  
raff = np.multiply(raff_mM, molar_mass) * 1e-3
ext_1 = np.multiply(ext1_mM, molar_mass) * 1e-3
ext_2 = np.multiply(ext2_mM, molar_mass) * 1e-3

# %% [markdown]
# ### Concentrations at the outlet ports
#
# To compare the simulation results to those of Mun ([Fig. 11 a, b, c](https://www.sciencedirect.com/science/article/pii/S002196731101363X?via%3Dihub#fig0055)), the concentration of every component is averaged over one switching period. This results in a new average every 324s. As there are 5 columns that switch a total of 41 times, this results in 205 total average concentrations for every component. Dividing the total simulation time by the switch time yields the same number of 205 steps for `n_averages`. The averaging is done for the **raffinate**, **extract 1** and **extract 2** ports.

# %%
n_averages = int(len(raff) // builder.switch_time)  

raff_average = (
    raff[: int(n_averages * builder.switch_time)]
    .reshape(n_averages, int(builder.switch_time), raff.shape[1])
    .mean(axis=1)
)
ext1_average = (
    ext_1[: int(n_averages * builder.switch_time)]
    .reshape(n_averages, int(builder.switch_time), raff.shape[1])
    .mean(axis=1)
)
ext2_average = (
    ext_2[: int(n_averages * builder.switch_time)]
    .reshape(n_averages, int(builder.switch_time), raff.shape[1])
    .mean(axis=1)
)

# %%
import matplotlib.pyplot as plt

n_steps = range(0,n_averages)

fig, axs = plt.subplots(2, 2, figsize=(20, 17))
ax1 = axs[0, 0]  # Raff
ax2 = axs[0, 1]  # Ext2
ax3 = axs[1, 0]  # Ext1

ax1.plot(n_steps, raff_average)
ax1.set_title("Raffinate")
ax1.set_xlabel("Step number")
ax1.set_ylabel("c [g / L]")
ax1.set_ylim(0, 0.8)

ax2.plot(n_steps, ext2_average)
ax2.set_title("Extract 2")
ax2.set_xlabel("Step number")
ax2.set_ylabel("c [g / L]")
ax2.set_ylim(0, 0.8)

ax3.plot(n_steps, ext1_average)
ax3.set_title("Extract 1")
ax3.set_xlabel("Step number")
ax3.set_ylabel("c [g / L]")
ax3.set_ylim(0, 0.8)
plt.show()


# %% [markdown]
# ### Axial concentrations
#
# The internal concentration profiles of each component in every zone of the SMB can be plotted for a desired time point within the seperation process. The following graph depicts the axial concentrations after the CSS is reached at a point **before a switching time** (200.99 steps) and **after a switch** (200.01 steps). This can be archieved by using the `axial_conc.plot_at_time` function already implemented in the `CarouselSolutionBulk` class (see {ref}`four-zone-binary). To compare the results to [Fig. 10 a,b](https://www.sciencedirect.com/science/article/pii/S002196731101363X?via%3Dihub#fig0050), the concentrations are converted from mM to g/L and plotted manually in the following code.

# %%
# Axial concentrations
from CADETProcess.modelBuilder.carouselBuilder import CarouselSolutionBulk
axial_conc = CarouselSolutionBulk(builder, simulation_results)

plotting_time = [200.01*builder.switch_time, 200.99*builder.switch_time,] 
n_cols = axial_conc.builder.n_columns
x = axial_conc.axial_coordinates

# Conversion of axial concentration plot from mM to g/L
for t in plotting_time:
    t_i = np.where(t <= axial_conc.time)[0][0]
    
    y_min_data = 0
    y_max_data = 0
    zone_counter = 0
    column_counter = 0
    _lines = []
    
    fig, axs = plt.subplots(
    ncols=n_cols,
    figsize=(n_cols*4, 6),
    gridspec_kw=dict(wspace=0.0, hspace=0.0),
    sharey='row')

    for position, ax in enumerate(axs):
        col_index = axial_conc.builder.column_indices_at_time(t, position) 
        y_data = axial_conc.solution[f'column_{col_index}'].bulk.solution[t_i, :]
        y = np.multiply(y_data, molar_mass) * 1e-3
        y_min_data = min(y_min_data, min(0, np.min(y)))
        y_max_data = max(y_max_data, 1.1*np.max(y))
        if position == 0:
            ax.set_ylabel("c [g/L]")
        if position == 2:
            ax.set_xlabel(f'axial coordinates at {t/builder.switch_time} steps')
        l = ax.plot(x, y)
      
        _lines.append(l)
    
        zone = axial_conc.builder.zones[zone_counter]
        ax.set_title(f'{zone.name}')
        plt.tight_layout()
        
        zone_counter += 1
        column_counter = 0
