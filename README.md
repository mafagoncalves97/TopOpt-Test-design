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

### **3. Installation and Usage**

1. **Set Up MATLAB Environment**:
   - Open MATLAB.
   - Navigate to the code folder.

2. **Prepare Input File**:
   A `.dat` file with the problem settings (details in the parameters section) is required and must be placed in the `Input` folder.

3. **Run the Main Script**:
   Execute the main script (**`main.pt`**)

---

### **4. Parameters**

These parameters need to be detailed in the input file:

- `Filename`: Problem name
  
- **Domain geometry and boundary conditions**:
  - `Domain`: x,y - number of elements in x and y
  - `Domain Type`: tot or sym, in case of total domain or only a quarter
  - `Springs`: spring_in, spring_out for the numerical stiffness values of the artificial springs in the input and output locations.
  - `BCs`: number of applied displacements (usually, 2)
  - `Fin`: load applied by the grips
    
- **Material Properties**:
  - `Material`: E, nu.
    - `E`: Young's modulus.
    - `nu`: Poisson's ratio.

- **Optimization Settings**:
  - `Vol`: Target volume fraction.
  - `Penal`: Penalization factor for intermediate densities.
  - `Radius`: Filter radius for sensitivity filtering.
  - `Constraints`: V, B and/or D, for volume, buckling, and damage constraints. 

- **Analysis settings**:
  - `Analysis_type`: n or l, for nonlinear or linear analysis.
  - `Material behavior`: l or n, for linear elastic or nonlinear material behavior (isotropic hardening and plasticity implemented).
  - `Integration_type`: f or r, for full or reduced integration.


These parameters can be adjusted within the input file .dat, which needs to be in the `Input` folder.

---

### **5. Output**

The optimization process generates:

- **Design Variables**: Optimal material distribution matrix.
- **Performance Metrics**: Objective function values and constraint violations.
- **Visualization Files**: Graphs and plots illustrating the optimization progress and results.

---

### **6. Contact**

For questions or further assistance, please contact [Mafalda Gonçalves](mafalda.goncalves@ua.pt).
