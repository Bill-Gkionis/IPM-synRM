# 🧲 Custom Interior Permanent Magnet (IPM) Motor  
Electromagnetic Design • FEA Analysis • Dynamic Modeling • FOC Simulation

---

## 📌 Project Overview

This project focuses on the design and modeling of a custom Interior Permanent Magnet (IPM) synchronous motor intended for high torque density actuator applications.

The objective is to:

- Design rotor and stator electromagnetic geometry
- Perform 2D finite element analysis (FEA)
- Extract machine parameters (Ld, Lq, flux linkage)
- Develop a dynamic dq-axis motor model
- Simulate Field-Oriented Control (FOC)
- Evaluate torque-speed performance
- Structure the model for multi-objective optimization

Workflow integration:

SolidWorks → FEMM → MATLAB/Simulink (Simscape Electrical)

---

## 🧩 Electromagnetic Geometry

### Stator
- 36-slot distributed winding
- Outer diameter: 128 mm
- Inner diameter: 91.6 mm
- Slot geometry optimized for flux guidance

### Rotor
- 8-pole IPM topology
- V-shaped interior magnets
- Saliency for increased Ld/Lq ratio

---

## ⚡ Finite Element Analysis (FEMM)

- 2D magnetostatic simulation
- Material: M-19 Steel
- Magnets: N45
- ~91,000 nodes
- ~182,000 elements

### Flux Density Observations
- Peak flux density ≈ 1.3 T in stator teeth
- No severe back iron saturation under nominal excitation

---

## 📊 Parameter Extraction (In Progress)

Machine parameters extracted from FEM results:

- d-axis inductance (Ld)
- q-axis inductance (Lq)
- Flux linkage (ψf)

Integrated into dynamic performance model.

---

## 🔁 Dynamic Modeling (MATLAB / Simulink)

dq-axis formulation:

Vd = R·Id − ωe·Lq·Iq  
Vq = R·Iq + ωe·(Ld·Id + ψf)

Implemented using:
- MATLAB/Simulink
- Simscape Electrical (high-fidelity PMSM block)

---

## 📈 Performance Analysis

Simulated:
- Torque-speed characteristics
- Efficiency trade-offs
- Closed-loop FOC current response

Model structured for:
- Multi-objective optimization
- Torque density improvement
- Loss minimization

---

## 🏗 Mechanical Integration

- Rotor inertia estimation from CAD geometry
- Prepared for integration into dynamic model

---

## 🛠 Tools Used

- SolidWorks
- FEMM
- MATLAB / Simulink
- Simscape Electrical
- LTspice

---

## 🚀 Next Steps

- Automate FEM-to-MATLAB parameter extraction
- Improve model validation accuracy
- Add thermal modeling
- Extend loss modeling to include proximity effect and hysteresis losses
- Incorporate core saturation effects into the electromagnetic model
- Analyze slot harmonics and higher-order field harmonics impact on performance
- Implement optimization routine
