#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 10:14:48 2022

@author: fabian
"""

import pypsa
import matplotlib.pyplot as plt
import cartopy.crs as ccrs


n = pypsa.Network("/home/fabian/papers/co2-network/workflow/subworkflows/pypsa-eur-sec/results/your-run-name/postnetworks/elec_s_5_lv1.0__Co2L0-3H-T-H-B-I-A-solar+p3-dist1_2050.nc")

stats = n.statistics()

# %% Plot networks

kinds = list(filter(lambda s: "pipeline" in s, n.links.carrier.unique()))
kinds.append("DC")

fig, axes = plt.subplots(len(kinds), subplot_kw={"projection", ccrs.EqualEarth()})

for kind, ax in zip(kinds, axes):

    links = n.

    n.plot(ax=ax, )