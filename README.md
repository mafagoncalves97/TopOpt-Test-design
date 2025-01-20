# TopOpt - Test design - Documentation

### **1. Introduction**

This project comprises a set of MATLAB scripts aimed at designing heterogeneous mechanical tests using topology optimization techniques.

---

### **2. Code Structure**

#### **Main Script**
- **`main.pt`**: The main script orchestrates the entire process, including input handling and the design process. It is composed of:
  - **`inputfile_processing`**: Handles the processing and interpretation of input files, preparing parameters and settings for the optimization.
  - **`boosted_TopOpt`**: Executes the design process using topology optimization.

#### **Design Process Scripts**
These scripts and functions support the topology optimization process:

- **`FEA_analysis.m`**: Performs finite element analysis to evaluate the mechanical response of the structures.
- **`OptimizerRun.m`**: Auxiliary script to handle optimization runs.
- **`SEMDOT_analysis.m`**: Conducts sensitivity analysis for the design variables.
- **`SaveFile.m`**: Handles the saving of results and data generated during the optimization.
- **`VolumeConstraint.m`**: Imposes volume constraints to control material usage.
- **`computeB.m`**: Computes the strain-displacement matrix.
- **`computeInternalForce.m`**: Calculates internal forces within the structure.
- **`computeKE.m`**: Generates the element stiffness matrix.
- **`computeKE_ep.m`**: Computes the element stiffness matrix considering plasticity effects.
- **`computeStrainStress.m`**: Evaluates strain and stress distributions.
- **`filters.m`**: Applies filtering techniques to the design variables.
- **`init_data.m`**: Initializes data and parameters for the optimization process.

#### **Add-Ons**
These additional scripts extend the main optimization process with specific features or constraints:

- **`DamageConstraint.m`**: Implements constraints to account for potential damage within the material during optimization.
- **`BucklingConstraint.m`**: Adds constraints related to buckling behavior in the optimization process.
- **`Uderivative_dX.m`**: Computes derivatives of displacement with respect to design variables, primarily used for add-on constraints.

---

