# LIF-MatLab
# Change Log
All notable changes to this project will be documented in this file.

## 1.3 - 2017-01-05
### Changed
- Minor changes here and there


## 1.2 - 2016-11-29
### Changed
- The F function (a subfunction) in delta_synaptic_strength.m is changed. F_plus and F_minus functions are combined into one subfunction F.
- The new F function in delta_synaptic_strength.m now offers two types: 'Bi1998' and 'Gutig2003'.
- Change in synaptic strength (delta_ss, synaptic strength difference vector) is now the product of two functions: F and K. The first one, F, is weight-dependent. The second one, K, is time-dependent and it is referred to as the "temporal filter".
- In LIF_network.m, when a subset of all synaptic strengths (weights) are plotted at every progress update, the histogram setting is changed to show 'countdensity' rather than 'probability'. This setting choice shows the evolution of the weight distribution more clearly.

## 1.1 - 2016-11-27
### Added
- A new test (type=5) is added to test_LIF_network.m
- Synaptic strength lower and upper bound parameter declarations are moved to the main file LIF_network.m


## 1.0 - 2016-11-24
### Added
- This CHANGELOG file will be updated as my MatLab for the LIF network model is updated.

