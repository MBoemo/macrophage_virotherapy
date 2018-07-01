# macrophage_virotherapy
Code supporting Boemo, Byrne 2018

This document contains instructions for how to generate each figure.

## Figure 1
- Open MOL_Figure1_and_Figure2.m
- Change l_inf from 0.2 to 0 (to turn off macrophages)
- Run MOL_Figure1_and_Figure2.m
- Run plotter.m to get heatmaps of phase spatial distributions

## Figure 2
- Open MOL_Figure1_and_Figure2.m
- Ensure that l_inf is set to 0.2
- Run MOL_Figure1_and_Figure2.m
- Run plotter.m to get heatmaps of phase spatial distributions
- Compare this radius curve with the radius curve from Figure 1 to make Figure 2 (upper)

## Figure 3 (and Appendix B figure)
- Run parameter_sensitivity.m

## Figure 4
- No treatment curve: radius curve from Figure 1
- Open MOL_virus.m
- For macrophage-delivered virotherapy alone: set mac_turn_on = 200 and rad_start = 501 and run MOL_virus.m
- For radiation alone: set mac_turn_on = 501 and rad_start = 150 and run MOL_virus.m
- For radiation and non-engineered macrophages: set mac_turn_on=200 and rad_start=150.  Then set D_phi=p_phi_max=r_rep=r_phi=k_phi=c_phi=d_phi=0 and run MOL_virus.m.
- For radiation and macrophage-delivered virotherapy: set mac_turn_on = 200 and rad_start = 150 and run MOL_virus.m.
- To generate the phase plot: set mac_turn_on = 200 and rad_start = 150.  Run MOL_virus.m.  Then run plotter_virus.m.

## Figure 5
- No treatment curve is the radius from Figure 1
- Open MOL_virus.m
- set mac_turn_on = 150 and rad_start = 150 and run MOL_virus.m
- set mac_turn_on = 200 and rad_start = 150 and run MOL_virus.m
- set mac_turn_on = 250 and rad_start = 150 and run MOL_virus.m
- set mac_turn_on = 350 and rad_start = 150 and run MOL_virus.m

## Figure 6
- open MOL_virus.m and set maxT=3000
- set mac_turn_on = 200 and rad_start = 150 and run MOL_virus.m
- set mac_turn_on = 350 and rad_start = 150 and run MOL_virus.m
- for phase plot, mac_turn_on = 350 and rad_start = 150 and run MOL_virus.m.  Then run plotter_virus.m.

## Figure 7
- open MOL_virus_twiceRad.m
- ensure that rad_start=150 and mac_turn_on=200 and maxT=1500
- No treatment is radius curve from Figure 1
- set rad_start2 = 1501 and run MOL_virus_twiceRad.m to get the single radiation dose curve
- set rad_start2 = 200 and run MOL_virus_twiceRad.m
- set rad_start2 = 250 and run MOL_virus_twiceRad.m
- set rad_start2 = 300 and run MOL_virus_twiceRad.m
- set rad_start2 = 400 and run MOL_virus_twiceRad.m






