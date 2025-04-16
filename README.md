# New_Dark_Cosmos
This code computes basics tabulated functions that are read by MPGRAFIC to generate IC for cosmological simulations and the simulation code Ramses.
The input parameters can be read from a namelist file containing the value of the reduced hubble constant h, the matter density (baryon+CDM), the CPL equation of state parameters (w_0,w_a) and a string variable specifying whether you want to include or not the radiation (RAD / NO_RAD).
The code generate two outputs, one for MPGRAFIC and one for Ramses.

MPGRAFIC output: a, a^2 E(a), D_+(a), f(a), D_+(2)(a)

where a is the scale factor, E(a) the Hubble function, D_+ the linear growth factor, f(a)=dlnD_+/dlna the linear grwoth rate and D_+^(2)(a) the linear growth factor to second order for LPT initial conditions

RAMSES output:  a, a^2 E(a), Superconf_time(z)-Superconf_time(0.d0), Lookback_time(0.d0)-Lookback_time(z), Conformal_time(z)-Conformal_time(0.d0)

