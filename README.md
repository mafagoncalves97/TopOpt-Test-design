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

### **4. Installation and Usage**

1. **Set Up MATLAB Environment**:
   - Open MATLAB.
   - Navigate to the code folder.

3. **Prepare Input File**:
   A `.dat` file with the problem settings (details in the parameters section) is required and must be placed in the `Input` folder.

4. **Run the Main Script**:
   Execute the main script (**`main.pt`**)

---

### **5. Parameters**

Key parameters influencing the optimization include:

- **Material Properties**:
  - `E`: Young's modulus.
  - `nu`: Poisson's ratio.

- **Optimization Settings**:
  - `volfrac`: Target volume fraction.
  - `penal`: Penalization factor for intermediate densities.
  - `rmin`: Filter radius for sensitivity filtering.

These parameters can be adjusted within the input file .dat, which needs to be in the `Input` folder.

---

### **6. Output**

The optimization process generates:

- **Design Variables**: Optimal material distribution matrix.
- **Performance Metrics**: Objective function values and constraint violations.
- **Visualization Files**: Graphs and plots illustrating the optimization progress and results.

---

### **7. Contact**

For questions or further assistance, please contact [Mafalda Gonçalves](mafalda.goncalves@ua.pt).
