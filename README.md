# IPM synRM: an open motor design

A **36-slot, 8-pole permanent-magnet-assisted synchronous reluctance motor** (PMaSynRM): two N45 magnets per pole in a V, M-19 steel, a 90 mm rotor in a 128 mm stator. It is designed, simulated and optimized with open tools, and documented step by step, so that each result can be checked and reused. The motor is built to learn from, not for a specific application.

![The starting design and the optimized design, side by side: geometry and magnetic field](docs/figures/optimization/reference-vs-optimized.png)

*The starting design (left) and the optimized design (right), at the same scale and the same copper loss, each at its own best current angle: the geometry with the changed dimensions, and the magnetic field. FEMM 5.0, 2D, fine mesh.*

> **Status (November 2026).**
> - **Done, and published here:**
>   - the cross-section, optimized;
>   - its rotor, checked for stress at speed;
>   - its flux maps over the whole current plane, with cross-saturation and the iron loss;
>   - its thermal model;
>   - its copper losses, with the proximity effect, and the winding trade-offs;
>   - its torque–speed, loss and efficiency estimates ([Motor map.md](Motor%20map.md));
>   - an interactive flux page of the three designs ([open it](https://bill-gkionis.github.io/IPM-synRM/explorer/)).
> - **Next:** the MATLAB validation of the flux maps, then a flux-map-driven Simulink plant.
> - **Not yet measured:** every number here is a finite-element estimate. Nothing has been measured on hardware yet.

## The interactive flux page

[![The Flux Explorer: the optimized motor's field, currents and torque](docs/figures/explorer/flux-explorer.png)](https://bill-gkionis.github.io/IPM-synRM/explorer/)

**[Open the Flux Explorer](https://bill-gkionis.github.io/IPM-synRM/explorer/)**: the three designs, the first pass with its air band, the band removed and the optimized motor, from 11,718 FEMM solves. It shows:
- the field, with the rotor turning;
- the current dial, with the MTPA angle;
- the phase currents and their vectors;
- the torque ripple;
- the current plane, with field weakening as the speed rises;
- the 3,400 designs the optimizer scored.

One self-contained page, running in the browser. Pictures between solved points are blended, and the page says so.

## At a glance

| | |
|---|---|
| Topology | PM-assisted synchronous reluctance, interior V magnets, two per pole |
| Slots / poles | 36 / 8, double-layer distributed winding, q = 1.5, winding factor 0.945 |
| Stator | outer diameter 128 mm, bore 91.6 mm |
| Rotor | outer diameter 90 mm (air gap 0.8 mm), shaft Ø15 mm |
| Magnets | N45 sintered NdFeB, 1.57 × 9.08 × 19.14 mm |
| Steel | M-19, laminated |
| Supply | a 48–60 V DC bus (ultra-low voltage) |
| Axis convention | d = the magnet axis, so L_d < L_q |
| The optimized cross-section | outer bridge 1.4 mm, centre rib 1.0 mm, V angle 111°, slot opening 1.5 mm, tooth 3.0 mm, the same magnets. Against the starting design, at the same copper loss: **+60 % torque, and 4.8 % torque ripple instead of 18.5 %** |

**What it can do, estimated** ([Motor map.md](Motor%20map.md) has every assumption). On a 60 mm stack and a 48 V bus, with a fan on the shaft:
- **4.2 N·m continuously** at 3,000 rpm, and **8.7 N·m** for about 25 minutes from a cold start;
- **2.97 kW peak** short-term, and still 2.7 kW at 16,330 rpm, the rotor's safe speed;
- **a best efficiency of 95.2 %** (motor only).

These are 2D finite-element estimates, not measurements.

## The pipeline

```mermaid
flowchart LR
  G["Geometry<br/>drawn in Python<br/>(motorgeo)"] --> F["FEMM 5.0<br/>2D magnetostatics,<br/>parallel workers"]
  F --> O["Optimizers<br/>NSGA-II · a full grid ·<br/>a blind Bayesian check"]
  O --> V["Verification<br/>CalculiX stress · Elmer ·<br/>ROSS · SALib"]
  V --> D["Motor data<br/>flux maps · thermal model ·<br/>copper losses"]
  D --> W["Winding<br/>every candidate,<br/>with the heat"]
  W --> P["Motor map<br/>torque–speed, losses,<br/>efficiency"]
  P --> N["Next<br/>MATLAB validation ·<br/>Simulink plant · inverter"]
```

**Three principles:**
- **The geometry is generated, not drawn.** One Python recipe writes the FEMM model, a DXF for CAD, and the meshes for Gmsh and Elmer, so an optimizer can change any dimension, and every candidate is checked for manufacturability before a solver opens.
- **One heavy run at a time.** Every solver run (FEMM, CalculiX, Elmer, MATLAB) takes a slot from a queue, runs at below-normal priority under a timeout, and ends only the processes it started.
- **Every result is checked by a second route:** a known analytic answer, a different solver, or an independent optimizer that never saw the first one's results.

The code is in [`python/`](python/), with a guide to every file in [`python/README.md`](python/README.md).

## How it got here

### 1. The first pass (February–April 2026)

SolidWorks parts, a FEMM model with a rotor sweep by hand from 5° to 90°, the flux linkage extracted in MATLAB, and a first field-oriented-control plant in Simulink. The plant ran unloaded, with constant motor parameters.
- **Files:** [`cad/`](cad/), [`femm/`](femm/), [`matlab/`](matlab/), and screenshots in [`docs/figures/`](docs/figures/).

### 2. A scriptable stack (October 2026)

The cross-section was redrawn in Python (`motorgeo`), and FEMM was driven from Python on a 90° sector with periodic boundaries and an air-gap element, in parallel workers.
- The new model reproduces the hand-drawn FEMM model's areas exactly (0.000 %).
- The 90° slice matches the full 360° model within about 1.5 %.

![The cross-section drawn by motorgeo](docs/figures/geometry/motorgeo-cross-section.png)

### 3. The first study: the air band (October 2026)

The first-pass rotor had a 1.5 mm air band behind the magnets. Filled with steel, the same rotor gave **37 % more magnet flux and 31–38 % more torque** at the best current angle, with the same reluctance torque and per-cent ripple; the cogging doubled. The band sat in the magnets' return path. **It was removed.** *(FEMM 4.2.)*

![Flux and flux density with and without the air band](docs/figures/air-band/field-maps.png)

<details><summary>More from the air-band study</summary>

![Torque against current angle, with and without the band](docs/figures/air-band/torque-vs-current-angle.png)
![No-load flux linkage and back-EMF](docs/figures/air-band/no-load-flux-and-back-emf.png)

</details>

### 4. The tools, each proven on a known answer (November 2026)

- **FEMM 5.0** replaced FEMM 4.2 after a re-run of every study:
  - every mean torque and flux linkage within 1 % (fine mesh) or 2 % (coarse mesh), apart from a few cases the comparison explains;
  - the whole re-run took 8.6 minutes instead of 2 h 13 min.
- **CalculiX, Gmsh, Elmer, PrePoMax, ROSS and SALib** were each checked against a problem with a known answer, then run on the motor's own rotor in a test bench of 141 checks.
  - 116 of the 119 gating checks pass. The 3 that fail are explained: the chord error of curved edges.
- **The bench's first numbers on the band-free rotor:**
  - the bridges' peak stress: 33.45 MPa at 10,000 rpm (CalculiX);
  - the first critical speed: 22,838 rpm (ROSS, a 60 mm stack and an estimated bearing span);
  - Elmer agrees with FEMM at no load within 0.33 % on its default mesh.

<details><summary>The bench's figures</summary>

![FEMM 5.0: flux density and flux lines of the band-free rotor](docs/figures/tools/femm-flux.png)
![FEMM 5.0 against FEMM 4.2: mesh convergence](docs/figures/tools/femm-mesh-convergence.png)
![The rotor's stress in CalculiX](docs/figures/tools/calculix-bench.png)
![Elmer against FEMM](docs/figures/tools/elmer-bench.png)
![The rotor's Campbell diagram in ROSS](docs/figures/tools/ross-campbell.png)
![Morris screening of the rotor's dimensions with SALib](docs/figures/tools/salib-morris.png)
![M-19's B-H curves](docs/figures/tools/m19-bh-curves.png)

</details>

### 5. The physics, tested in FEMM (November 2026)

On the band-free rotor:
- **Torque in the d-q frame.** At 10 A the torque splits into three field-times-current terms: 0.918 + 0.461 − 0.253 = 1.127 N·m.
- **The magnets make the saliency.** With their coercivity set to zero, the unsaturated bridges short the barrier and L_d = L_q (saliency 1.00, against 1.82 with the magnets).
- **MTPA, the most torque per amp.** The best current angle climbs from 22.9° at 10 A to 43.3° at 40 A, above the constant-parameter formula's 37.5°. Fed each current's own inductances, the formula even turns the wrong way, because of **cross-saturation**. So the controller's table has to come from flux maps.

![Torque against current angle](docs/figures/physics/torque-vs-current-angle.png)

<details><summary>More</summary>

![The bridges, with the magnets on and off](docs/figures/physics/bridges-magnets-on-off.png)
![The best current angle against current](docs/figures/physics/mtpa-angle-vs-current.png)
![The MTPA path in the i_d–i_q plane](docs/figures/physics/mtpa-path.png)

</details>

### 6. The optimization (November 2026)

**The landscape came first:** sweeps of each dimension, two maps, and a Morris screening. Then each part got the optimizer its landscape called for:
- **NSGA-II for the rotor,** a full grid for the stator, NSGA-II for both together, and every combination for the winding;
- **a blind second run** with a different method (a Bayesian search on Gaussian-process surrogates) found the same design region, with a slightly better front in fewer designs. The chosen design comes from its front.

**The size of the search:** 4,156 designs, each scored with 23 FEMM solves, about 96,000 solves in all:
- the sweeps, 736;
- the screening, NSGA-II and the grid, 1,973;
- the blind run, 1,447.

**The result, on the fine mesh, at the starting design's copper loss:**

| | Starting design | Optimized |
|---|---|---|
| Mean torque | 1.129 N·m | **1.807 N·m (+60 %)** |
| Torque ripple | 18.45 % | **4.84 %** |
| Magnet flux linkage | 16.64 mWb | 24.20 mWb (+45 %, the same magnets) |

**Why:**
- the centre rib went from 3.3 to 1.0 mm, so far less of the magnets' flux is spent saturating steel;
- the bigger slot carries more current at the same loss;
- the V angle and the narrower slot opening cut the ripple.

The side-by-side picture at the top of this page shows both designs.

**Four of the six variables sit on a limit:** the bridge (1.4 mm), the rib (1.0 mm), the tooth (3.0 mm) and the slot opening (1.5 mm). Those limits are the design decisions.

![The fronts, re-scored on the fine mesh](docs/figures/optimization/fronts-fine-mesh.png)

*The torque–ripple fronts of each run (the rotor NSGA-II, the stator grid, the joint NSGA-II and its continuation, and the blind Bayesian cross-check), and their union re-scored on the fine mesh. The chosen design sits at the knee.*

<details><summary>More from the optimization</summary>

![The optimized motor, whole](docs/figures/optimization/chosen-motor.png)
![The independent Bayesian run against NSGA-II](docs/figures/optimization/independent-run.png)
![The rotor's ripple map](docs/figures/optimization/rotor-ripple-map.png)
![Morris screening of the ripple](docs/figures/optimization/morris-ripple.png)
![The rotor's front](docs/figures/optimization/rotor-front.png)
![The stator at fixed current density against fixed ampere-turns](docs/figures/optimization/stator-fixed-j-vs-at.png)
![The radius sweep](docs/figures/optimization/radius-sweep.png)
![The winding: NSGA-II against the exact front](docs/figures/optimization/winding-nsga-vs-exact.png)

</details>

### 7. The stress check (November 2026)

CalculiX, in 2D plane stress, with the magnets as loose loads on their slot walls (the upper bound):
- **The optimized rotor peaks at 39.06 MPa at 10,000 rpm,** now in the 1.0 mm centre rib. The starting design peaked at 33.45 MPa, in a bridge.
- **Safe speed:** with half the steel's yield stress at a 1.2 × overspeed, on placeholder strength values, it is safe to **16,330 rpm**. The speed range ends there.

![Stress at 10,000 rpm, starting design against optimized](docs/figures/stress/reference-vs-optimized.png)

<details><summary>More</summary>

![The peak stress against the rib](docs/figures/stress/peak-vs-rib.png)
![The peak stress against the bridge](docs/figures/stress/peak-vs-bridge.png)

</details>

### 8. The flux maps (November 2026)

The optimization scored one operating point. A drive needs the whole current plane. So the optimized motor's flux linkages ψ_d and ψ_q, its torque and its ripple were mapped over (i_d, i_q), up to 3 × the rated current and into field weakening:
- **425 operating points, 15,202 FEMM solves** on the fine mesh. The negative-i_q half was solved, not mirrored, because the winding isn't mirror-symmetric about the d axis.
- **The checks:**
  - the d-q torque 1.5 p (ψ_d i_q − ψ_q i_d), fed from the maps, matches FEMM's own torque to **0.085 %** at worst;
  - the two cross inductances agree, as reciprocity says they must (median difference 0.015 %).
- **Cross-saturation is strong.** The best current angle (MTPA) is 23.5°, 37.4° and **43.3°** at 1, 2 and 3 × rated. A constant-parameter formula says 17.3°, 26.2° and **30.9°**, and at 3 × rated it puts the torque **41 % too high**. So the controller's MTPA table and the Simulink plant both come from the maps.
- **The iron loss:**
  - **The model:** a Bertotti model fitted to M-19 data (rms **6.9 %** cross-validated), applied to every harmonic of each mesh triangle's flux-density waveform over a full electrical period.
  - **The result:** 9.0 W at the rated point at 3,000 rpm, and 129 W at 16,330 rpm (a 20 mm stack, 0.3 mm laminations).
- **For MATLAB and Simulink:** the maps, with five scripts for Simscape's FEM-Parameterized PMSM block, are in [`matlab/maps/`](matlab/maps/). They are prepared, and not yet run in MATLAB.

![Cross-saturation: ψ_d against i_d at four values of i_q, and the two cross inductances](docs/figures/maps/cross-saturation.png)

<details><summary>More from the flux maps</summary>

![ψ_d and ψ_q over the current plane](docs/figures/maps/flux-maps.png)
![The torque over the current plane, and the MTPA path from the maps against the constant-parameter formula's](docs/figures/maps/torque-and-mtpa.png)
![Apparent and incremental inductances](docs/figures/maps/inductances.png)
![The Bertotti model against the M-19 data](docs/figures/maps/m19-loss-fit.png)
![The iron loss along the MTPA path at three speeds, by region](docs/figures/maps/iron-loss.png)

</details>

### 9. The thermal model (November 2026)

**The model:** a lumped thermal network of eleven parts:
- the winding in the slots and its end turns;
- the teeth and the yoke;
- the housing, the end bells and the air in the end spaces;
- the rotor, the magnets, the shaft and the bearings.

Each loss enters at its own part. A steady state solves in a fraction of a millisecond, so the winding study and the motor map can call it at every speed.

**Calibrated first:**
- FEMM's heat flow and CalculiX matched two exact solutions within 0.13 %;
- the network was then checked against a 2D finite-element slice of the stator, and agrees **within 2 %** on the hot spot and on the winding's mean temperature;
- FEMM and CalculiX agree with each other to 0.033 %.

**The cooling** is a fan on the shaft over a finned housing. Its airflow follows the speed, and its own power grows as the cube of the speed, as a load on the shaft.

**What sheds the heat:**
- **with the shaft fan,** forced air over the housing's fins and the end bells carries 73–95 % of it;
- **in still air,** radiation carries 58 %;
- **inside the motor,** the copper in the slots passes about three quarters of its heat to the teeth.

**What limits the rating:** the N45 magnets reach their 80 °C before the winding nears its 155 °C, in every cooling, because the air gap ties the rotor to the stator's temperature. A 150 °C magnet grade would hand the limit to the winding, for 1.4–2.4 × the current.

![Where the heat goes at each cooling's continuous limit](docs/figures/thermal/heat-paths.png)

<details><summary>More from the thermal model</summary>

![The network and the thermal resistance of every path](docs/figures/thermal/network.png)
![The continuous limit against speed for each cooling](docs/figures/thermal/continuous-limit-vs-speed.png)
![The copper loss each cooling can carry, per stack length](docs/figures/thermal/sustainable-current.png)
![Heating from cold, and overloads](docs/figures/thermal/transients.png)

</details>

### 10. The copper losses and the winding (November 2026)

**The method:**
- FEMM's time-harmonic solver was first checked against the exact solutions for a round wire, over the range of strand sizes the winding uses;
- then the optimized slot was solved strand by strand (416 strands), in the machine's own field over a rotation, with the inverter's PWM ripple and the current circulating between parallel strands.

**What it found:**
- **The proximity loss is small with fine strands:** 0.42 % of the DC loss at 16,330 rpm for 0.315 mm strands, and 1.7 % for 0.63 mm. The PWM ripple adds 0.017 %.
- **The term that matters** is current circulating between a turn's parallel strands when they're spread over the slot's depth: up to 14 % at the top speed, against 0.3–0.4 % when they're bundled together.

**The winding:** every winding that fits the slot, 20,316 of them, was then enumerated against:
- the full copper loss;
- the heat, from the thermal model (step 9);
- the inverter current that the 48 V bus calls for.

The trade-off front shows the choice:
- fewer turns buy speed for current;
- each extra row of magnets raises the current.

The final pick waits for the stack length and the size of the shaft fan.

![R_ac / R_dc against frequency](docs/figures/copper-loss/rac-rdc-vs-frequency.png)

<details><summary>More on the copper losses and the winding</summary>

![One wire: the skin effect](docs/figures/copper-loss/one-wire.png)
![Two wires: the proximity effect](docs/figures/copper-loss/two-wires.png)
![The slot's strands at 20 kHz](docs/figures/copper-loss/slot-strands-20khz.png)
![The AC loss against speed](docs/figures/copper-loss/ac-loss-vs-speed.png)
![Circulating currents between parallel strands](docs/figures/copper-loss/circulating-currents.png)
![Every winding and the exact front](docs/figures/winding/exact-front.png)
![What the copper loss is made of](docs/figures/winding/copper-loss-parts.png)
![The strand arrangements](docs/figures/winding/strand-arrangements.png)

</details>

### 11. The motor map (November 2026)

The flux maps, the losses and the thermal network combine into what the motor can do across its speed range. [Motor map.md](Motor%20map.md) has every table, figure and assumption.

**On a 60 mm stack and a 48 V bus,** with the smallest of the three shaft fans studied (50 W at the top speed):
- **4.2 N·m continuously** at 3,000 rpm (1.32 kW);
- **8.7 N·m short-term,** up to 2,718 rpm, for about 25 minutes from a cold start;
- **2.97 kW peak,** and still 2.72 kW at 16,330 rpm. The short-term current equals the motor's characteristic current (84 A), so the power stays nearly flat to the top speed;
- **a best efficiency of 95.2 %** (motor only), at 3,746 rpm and 2.9 N·m.

**How it's solved and checked:**
- MTPA, field weakening and MTPV are solved at every speed from the maps;
- the energy balance closes at all 43,493 operating points;
- FEMM re-solves of deep field-weakening points reproduce the maps within 0.016 %.

**Also compared:** three fan sizes, four stack lengths, the winding's alternatives, a 150 °C magnet grade, and natural, external-fan and liquid cooling.

**A note for the inverter:**
- above 4,049 rpm the magnets' own voltage passes the 48 V bus, reaching 194 V at the top speed, so the inverter must hold a field-weakening current at speed;
- an active short circuit is mild here: 84 A steady, and a brake of 0.12 N·m at the top speed.

![Torque and power against speed with the three shaft fans](docs/figures/motor-map/envelopes-shaft-fans.png)

<details><summary>More from the motor map</summary>

![Where the loss goes along the envelope](docs/figures/motor-map/losses-fan-b.png)
![The efficiency map](docs/figures/motor-map/efficiency-fan-b.png)
![The current plane: MTPA, field weakening and MTPV](docs/figures/motor-map/current-plane.png)
![The highest speed each fan allows](docs/figures/motor-map/top-speed-per-fan.png)
![The stack table](docs/figures/motor-map/stack-table.png)
![Natural convection, an external fan and a liquid jacket](docs/figures/motor-map/envelopes-other-coolings.png)

</details>

### 12. Next

- **The MATLAB validation** of the flux maps: the scripts in [`matlab/maps/`](matlab/maps/) are prepared and not yet run in MATLAB. Then a **flux-map-driven Simulink plant** in place of the constant parameters.
- **Then:**
  - the magnets' demagnetisation and their eddy losses;
  - the rotor in 3D;
  - the protection at speed;
  - the inverter.

## The model-fidelity checklist

| Item | Status |
|---|---|
| Core saturation and cross-saturation | ✅ published (step 8, the flux maps); the MATLAB tables prepared, not yet run |
| Higher field harmonics | 🟡 measured (air-gap spectrum, back-EMF, ripple); not yet in a plant |
| Hysteresis and iron loss | ✅ published (step 8: Bertotti, per mesh element) |
| Proximity and other AC copper losses | ✅ published (step 10) |
| High-frequency loss, magnet eddy currents | ⬜ Elmer, transient and 3D |
| Torque and inertia scaling with diameter, length and air gap | ⬜ |

## Repository layout

| Path | Contents |
|---|---|
| [`python/`](python/) | The pipeline's code: geometry, FEMM driver, solver queue, studies, optimization, checks and tests ([guide](python/README.md)) |
| [Motor map.md](Motor%20map.md) | What the optimized motor can do: torque, power, losses and efficiency against speed, with every assumption |
| [`matlab/maps/`](matlab/maps/) | The flux maps for Simscape's FEM-Parameterized PMSM block, and the scripts that use them (prepared, not yet run in MATLAB) |
| [`docs/figures/`](docs/figures/) | The figures of each step |
| [`docs/explorer/`](docs/explorer/) | The interactive flux page, served by GitHub Pages ([open it](https://bill-gkionis.github.io/IPM-synRM/explorer/)) |
| [`femm/`](femm/) | The first-pass FEMM model and its rotor-sweep screenshots |
| [`matlab/`](matlab/) | The first-pass flux-linkage extraction, the Clarke and Park transforms, the Simulink plants, and `maps/` (above) |
| [`cad/`](cad/) | SolidWorks parts, assembly, drawings and equation files |

FEMM's result files (`.ans`) aren't stored: FEMM regenerates them from the models. The state before the reorganisation is tagged `v0`.

## Running the code

The pipeline was built and run on one Windows PC, with FEMM 5.0 and `uv`. The tool paths live in [`python/config.toml`](python/config.toml) and are set for that machine. A per-machine configuration is planned.

```
cd python
uv sync
uv run pytest -m "not femm"      # the tests that don't need FEMM
uv run python <script>.py        # each script's docstring says what it does and how to run it
```

## How the work is done

Each step runs as a loop:
1. set a goal;
2. decide the step;
3. test it in simulation;
4. check the numbers against a second route;
5. change the decision when the results disagree.

A decision stays a working assumption until a test backs it. The numbers are 2D finite-element estimates (FEMM 5.0 on a 20 mm slice, unless a step says otherwise), and a step says when it scales them.

## License

[Apache-2.0](LICENSE)
