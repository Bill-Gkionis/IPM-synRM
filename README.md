🧲 Custom Interior Permanent Magnet (IPM) Motor

Electromagnetic Design • FEA Analysis • Dynamic Modeling • FOC Simulation

📌 Project Overview

This project focuses on the design and modeling of a custom Interior Permanent Magnet (IPM) synchronous motor intended for high torque density actuator applications.

The objective is to:

Design rotor and stator electromagnetic geometry

Perform 2D finite element analysis (FEA)

Extract machine parameters (Ld, Lq, flux linkage)

Develop a dynamic dq-axis motor model

Simulate Field-Oriented Control (FOC)

Evaluate torque-speed performance

Structure the model for multi-objective optimization

The workflow integrates:

SolidWorks → FEMM → MATLAB/Simulink (Simscape Electrical)

🧩 1. Electromagnetic Geometry Design
Stator

36-slot distributed winding configuration

Outer diameter: 128 mm

Inner diameter: 91.6 mm

Slot geometry optimized for manufacturability and flux guidance

(Insert stator CAD image here)

Rotor

8-pole IPM topology

V-shaped interior magnets

Saliency introduced to increase Ld/Lq ratio

Designed for reluctance + PM torque contribution

(Insert rotor CAD image here)

⚡ 2. Finite Element Analysis (FEMM)

2D magnetostatic simulations performed using FEMM:

Material: M-19 Steel (stator/rotor)

Neodymium permanent magnets (N45)

Air-gap field distribution analysis

Flux density mapping

Saturation evaluation

Mesh:

~91,000 nodes

~182,000 elements

Flux Density Distribution

(Insert FEM flux density plot here)

Observations:

Peak flux density ~1.3 T in stator teeth

Proper flux guidance through rotor saliency bridges

No severe saturation in back iron (under nominal excitation)

📊 3. Parameter Extraction (In Progress)

Machine parameters are being extracted from FEM results:

d-axis inductance (Ld)

q-axis inductance (Lq)

Flux linkage (ψf)

Electromagnetic torque estimation

These parameters are integrated into a dynamic motor model for performance evaluation.

🔁 4. Dynamic Modeling (MATLAB / Simulink)

The motor is modeled in dq reference frame using:

𝑉
𝑑
=
𝑅
𝐼
𝑑
−
𝜔
𝑒
𝐿
𝑞
𝐼
𝑞
V
d
	​

=RI
d
	​

−ω
e
	​

L
q
	​

I
q
	​

𝑉
𝑞
=
𝑅
𝐼
𝑞
+
𝜔
𝑒
(
𝐿
𝑑
𝐼
𝑑
+
𝜓
𝑓
)
V
q
	​

=RI
q
	​

+ω
e
	​

(L
d
	​

I
d
	​

+ψ
f
	​

)

Implemented in:

MATLAB/Simulink

Simscape Electrical (high-fidelity PMSM block)

Control Structure

Field-Oriented Control (FOC)

Closed-loop current regulation

Torque-speed simulation

(Insert Simulink model screenshot here)

📈 5. Performance Analysis

Simulated:

Torque-speed characteristics

Efficiency trade-offs

Current response under FOC

Sensitivity to parameter variation

The model is structured for future:

Multi-objective optimization

Torque density improvement

Loss minimization

🏗 6. Mechanical Integration

Rotor inertia estimation from CAD geometry

Parameter prepared for integration into dynamic model

🚀 Next Development Steps

Automate FEM-to-MATLAB parameter extraction workflow

Improve model validation accuracy

Add thermal modeling

Implement optimization routine (multi-criteria)

🛠 Tools Used

SolidWorks (geometry design)

FEMM (2D electromagnetic FEA)

MATLAB / Simulink

Simscape Electrical

LTspice (power electronics validation)
