# Kinematically-driven 2d oceanic subduction

*This section was contributed by Anne Glerum and Yingying Li.*

This subduction example is based on the benchmark effort of Quinquis et al.,
of which initial results were published in {cite:t}`quinquis:2014`. In four 
increasingly complex cases we will go from isoviscous materials without any
temperature effects to a fully thermo-mechanical, nonlinear, strain-weakened
visco-plastic, externally-driven model of oceanic subduction. The setup is 
outlined in {numref}`fig:QQ_setup`. The models are run for 15 My
and slab tip depth, trench location, RMS velocity and temperature, and viscous
dissipation are monitored. In addition, we discuss the effects of the element
size of the subduction interface and crustal layers, viscosity averaging and
the solver tolerance.

```{figure-md} fig:QQ_setup
<img src="setup_Quinquis2014.*" style="width:80.0%" />

 Case 4 model setup. Copied from {cite:t}`quinquis:2014`.
```

## Case 1: Simple rheology

The Case 1 model setup considers seven materials (compositional fields) apart
from the background sublithospheric mantle (see {numref}`fig:QQ_case1_setup`), including (1) the Bulk Oceanic Composition (BOC), (2) Serpentinized HarzBurgite 
(SHB) and (3)"thermal" layer of the overiding plate, (4) the BOC, (5) SHB and (6) 
"thermal" layer of the subducting plate, and (7) the weak seed.

The geometry of these compositions is implemented as follows:

```{literalinclude} Case1_compositions.prm

```

No differences in material properties exist, except for density and viscosity,
so we use the multicomponent material model. To keep the boundaries between
fields as sharp as possible in terms of viscosity, we use the maximum
composition to determine the viscosity in each evaluation point (note that
this can be harder on the solver):

```{literalinclude} Case1_materialmodel.prm

```

Temperature effects are ignored. Subduction is driven by prescribed in- and
outflow through the right boundary (with a gradual transition of the flow
direction), all other boundaries are free slip. The volume of material that
flows in is balanced by the volume that is prescribed to flow out (this is
important as the model is incompressible). A weak crust along the plate
interface and the subducting lithosphere facilitates subduction. Through the
function plugin, we prescribe the in- and outflow:

```{literalinclude} Case1_velocity.prm

```

To follow the slab on its descent into the mantle, we use adaptive mesh
refinement based on viscosity in combination with the minimum refinement
strategy to ensure sufficient resolution in the crust and weak zone that allow
the slab to detach from the surface:

```{literalinclude} Case1_meshrefinement.prm

```

To monitor the model evolution, several diagnostic quantities are tracked over
time: the depth of the tip of the slab, the position of the trench, the RMS
velocity of the slab and the whole model domain, and the viscous
dissipation in the slab and total model domain. The computation of trench position
is through a plugin called 'trench_location'. It needs to be compiled as a shared library
which can be called from the input file.

```{literalinclude} Case1_postprocessing.prm

```


```{figure-md} fig:QQ_case1_setup
<img src="Case1_t0.*" style="width:60.0%" />

 Case 1 density, viscosity and velocity at time zero.
```

We run the Case 1 setup for 15 My of model time. The diagnostic quantities in
{numref}`fig:QQ_case1_diagnostics` show three stages of model evolution: first trench advance (top right plot), then free subduction (increasing slab RMS velocity), and
after about 13 My interaction between the slab and bottom boundary at 660 km
depth, which slows down the slab. The slab then curves inward along the bottom
boundary. This can also be seen in {numref}`fig:QQ_case1_results`.

```{figure-md} fig:QQ_case1_diagnostics
<img src="Case1_diagnostics.*" />

 Case 1 diagnostic quantities of ASPECT
```


```{figure-md} fig:QQ_case1_results
<img src="Case1_visc.*" style="width:60.0%" />

 Case 1 viscosity snapshots at 8 and 15 My.
```
## Case 2a: Model with thermal evolution 

Case 2a builds on Case 1 by including a heterogeneous initial temperature field 
that combines the plate cooling model for the temperature in the plates with an 
adiabatic temperature gradient below the plates. For the plate cooling model, 
we assume that the overriding plate is 40 My old and the subducting plate 70 My. 
The temperature function is too complex to comfortably write out as a function
expression in the input file, so it is put in a new initial temperature plugin
called ’subduction plate cooling’. This plugin needs to be compiled as a shared
library that can be called from the input file:

```{literalinclude} Case2a_temperatures.prm

```
The initial temperature distribution can be seen in  {numref}`fig:QQ_case2a_temperature_t0`

```{figure-md} fig:QQ_case2a_temperature_t0
<img src="Case2a_temp_t0.*" />

 Case 2a temperature distribution at time 0
```

The thermal expansivity is kept zero in Case 2a. Therefore the thermal evolution is 
not coupled to the mechanical evolution. A high conductivity (k = 183.33 Wm-1K-1) 
is defined for the upper mantle to enforce the mantle adiabat in the slow convection
system.

```{literalinclude} Case2a_materialmodel.prm

```

We can monitor the temperature evolution, however, by tracking the
root-mean-square (RMS) temperature and the isotherm depth using 
two plugins called "composition_trms_statistics" and "isotherm depth".
To this end, we add the postprocessors ’temperature statistics’ 
and "isotherm depth" to the prm file. The RMS temperature is tracked 
over the subducting plate (comprising BOC_SP, SHB_SP, thermal_SP) 
and the whole domain over time. And the isotherm temperature to be 
tracked is prescribed as 800 C.

```{literalinclude} Case2a_postprocessing.prm

```

{numref}`fig:QQ_case2a_diagnostics` shows the RMS temperature, 
800 C isotherm depth as well as the other diagnostics. The 800 C isotherm 
keeps going deeper, while in the end the descend speed becomes slower. 
The RMS temperature of the whole domain keeps decreasing due to the 
input of cold subducting plate into model domain. The slab RMS temperature 
cools at the initial stage, then it goes up due to being warmed up by the 
mantle. Notice that these other diagnostics are indeed the same as for Case 1.

```{figure-md} fig:QQ_case2a_diagnostics
<img src="Case2a_diagnostics.*" />

 ASPECT Case 2a and case 1, Elefant case 2a diagnostic quantities
```

## Case 2b: Model with a temperature-dependent density

Case2b builds on case 2a.  Compared to Case2a, the thermal expansivity in Case 2b
is set to 2.5e-5 K<sup>-1</sup> instead of 0. Therefore, temperature feeds into density 
and therefore stokes equations. In addition, the benchmark uses different reference temperatures for each composition, but only one value of reference temperature is 
allowed to be given to various compositional fields in ASPECT. Therefore, we adapted reference densities for different compositional fields. All the other parameters are 
identical with case2a:

```{literalinclude} kinematically_driven_subduction_2d_case2b.prm

```
The density evolution can be seen in {numref}`fig:QQ_case2b_density_evolution`. When the slab warms up, the density contrast between the slab and the mantle decreases.

```{figure-md} fig:QQ_case2b_density_evolution
<img src="Case2b_density_evolution.*" />
 Case2b density evolution
```
Therefore, the  subduction in this case is much slower than that in case1 and 2b (see {numref}`fig:QQ_case2b_diagnostics`). The subducting plate did not reach 
the bottom till 15 Myr, i.e., the third stage in case 1 and 2b is missing in this 
case. 

```{figure-md} fig:QQ_case2b_diagnostics
<img src="Case2b_diagnostics.*" />

 ASPECT Case2b diagnostic quantities
```

The root mean square velocity over the whole domain increases (see {numref}`fig:QQ_case2b_diagnostics`) .This is due to higher velocities in the 
mantle in the left side of the mantle where a second convection cell forms (see {numref}`fig:QQ_case2b_velocity_evolution`) 

```{figure-md} fig:QQ_case2b_velocity_evolution
<img src="Case2b_velocity_evolution.*" />

 Case2b velocity evolution
```