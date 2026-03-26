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
# (Example_SMB)=
# # Simulated Moving Bed Chromatography (SMB)

# %% [markdown]
# The following case studies are reproductions of part of the research results published in:
# * **Model-based control of a simulated moving bed chromatographic process for the separation of fructose and glucose** (Karsten-Ulrich Klatt, Felix Hanisch, Guido Dünnebier, Journal of Process Control 2002; 12(2):203-219.) <br> https://doi.org/10.1016/S0959-1524(01)00005-1.
# * **Efficient numerical simulation of simulated moving bed chromatography with a single-column solver** (Qiao-Le He, Samuel Leweke, Eric von Lieres, Computers & Chemical Engineering 2018; 111:183-198. doi:10.1016/j.compchemeng.2017.12.022.) <br>
# https://www.sciencedirect.com/science/article/pii/S0098135417304520
# * **Improving performance of a five-zone simulated moving bed chromatography for ternary separation by simultaneous use of partial-feeding and partial-closing of the product port in charge of collecting the intermediate-affinity solute molecules** (Sungyong Mun, Journal of Chromatography A 2011;  1218(44):8060-8074) <br> https://doi.org/10.1016/j.chroma.2011.09.015.

# %% [markdown]
# ## Theoretical background
#
# Continuous processes generally have many benefits over batch processes, like higher efficiency and throughput. However, a truly continuous process in chromatography is not practically feasable. <br>
# This so called **true moving bed chromatography** would entail the solid phase of the chromatography column (the bed) moving in the opposite direction to the mobile phase. This would induce the local separation of components in the feed solution based on their column binding behaviour. Those components could then be retrieved by different outlet streams located upstream and downstream of the feed inlet.
#
# ```{figure} ./figures/true_moving_bed_chr.jpg
# :width: 600px
# <div style="text-align: center">
#
# True moving bed chromatography schematic based on [Youtube Video by orochemtech](https://www.youtube.com/watch?v=xhhJxb48tgc) 
# <div>

# %% [markdown]
# **Simulated Moving Bed Chromatography** is a method to approach such a continuous process in practice. This is realized by having multiple chromatography columns connected to each other in a carousel. By periodically switching the inlet and outlet valves connected to the columns (***"column switching"***) the movement of the solid phase can be mimicked. 
#
# A standard SMB system contains four different external units connected to the columns:
# 1. **Feed** (Inlet): Component mixture
# 2. **Raffinate** (Outlet): Faster eluting component (elutes before feed plug flow)
# 3. **Extract** (Outlet): Slower eluting component (elutes after feed plug flow), interacts more strongly with the column solid phase<br>
# There can be multiple extract or raffinate outlets in a SMB system depending on the number of components to be separated.
# 4. **Desorbent / Eluent / Solvent** (Inlet): Solution to elute extract from column before raffinate plug flow enters the column again to prevent mixing of the separated components
#
# The position of each column relative to these external units determines their specific **“zone”** in the SMB system. Based on thier external units, every zone has a different function in the seperation process and a different flow rate within the columns. The most basic SMB system for binary separation is made up of four zones with one chromatography column in each zone.

# %% [markdown]
# ```{figure} ./figures/four_zone.jpg
# :width: 600px
# :align: center
#
# Schematic of four-zone SMB chromatography for binary separations. Column positions are periodically switched in opposite direction of liquid flow.
# [Fig. 1, He et al.](https://www.sciencedirect.com/science/article/pii/S0098135417304520#fig0001)

# %% [markdown]
# ## Content
# The following two case studies examine a binary separation using a **four-zone SMB** system and a ternary separation using a **five-zone SMB** system. 
#
# ```{toctree}
# :maxdepth: 1
# :hidden:
#
# four_zone_binary
# five_zone_ternary
# ```
