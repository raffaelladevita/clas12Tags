
# CLAS12 TARGETS

The engineering models are kept on 
the [engineering wiki](https://wiki.jlab.org/Hall-B/engineering/hallb_eng_wiki/index.php/Main_Page).


## Windows


## Geometry

To build the geometry, run the following command:

./targets.pl config.dat

---

# Experiments:


## RGF: Sept 22 2022

Author: Yu-Chun Hung

Experiment Description: Measure the structure function F2 of the neutron by 
tagging the spectator proton in the RTPC detector. Deuterium, hydrogen, and helium gas are used at 5.6 atm pressure.

### Variations: bonusD2, bounsH2, bonusHe

### Files used: 

### How to use it in the gcard
1. Set up the target with bonusD2 variations:
<detector name="experiments/clas12/targets/target" factory="TEXT" variation="bonusD2"/>
2. Shift the target position by 30 mm
<detector name="target">    <position x="0*cm"  y="0*cm"  z="-3.00*cm"  />  </detector>


## RGM: Sept 20 2022

Author: Justin Estee

Experiment Description: Short Range Correlations & Electrons for Neutrinos experiment 


### Variations: 
Liquid Targets (Standard 5cm liquid cell)
lH2 - liquid Hydrogen
lD2 - liquid Deuterium
lHe - liquid Helium

Solid Targets
RGM_2_C - single foil carbon target in beamline
RGM_2_Sn - single foil tin target in beamline

RGM_8_C_L - 4-foil carbon target in beamline. Large carbon target option.
RGM_8_Sn_L - 4-foil tin target in beamline. Large tin target option.

RGM_8_C_S - 4-foil carbon target in beamline. Small carbon target option.
RGM_8_Sn_S - 4-foil tin target in beamline. Small tin target option.

RGM_Ca - Calcium target cell

Scattering chamber variations
cad - standard scattering chamber and liquid cell target wall entrance and exit material
cadrgm - standard scattering chamber without liquid cell wall and upstream window, only downstream window and pipe included (use with solid targets)
### Files used


### How to use it in the gcard
Liquid Targets
Use "cad" variation with all liquid targets
<detector name="experiments/clas12/targets/cad/"   factory="CAD"/>
<detector name="experiments/clas12/targets/target" factory="TEXT" variation="lD2"/>

Solid Targets
Use "cadrgm" variation with all solid targets
<detector name="experiments/clas12/targets/cadrgm/"   factory="CAD"/>
<detector name="experiments/clas12/targets/target" factory="TEXT" variation="RGM_2_C"/>

