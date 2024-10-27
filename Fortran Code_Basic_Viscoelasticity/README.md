# Viscoelastic Material Model with Prony Series (UMAT for ABAQUS/Standard)

This repository contains a UMAT subroutine for implementing a viscoelastic material model using the Prony series in ABAQUS/Standard. The code provides a time-dependent stress-strain response by accounting for viscoelastic relaxation effects through the Prony series, allowing accurate modeling of materials with history-dependent behavior.

---

## Table of Contents

1. [Overview](#overview)
2. [Material Properties and Parameters](#material-properties-and-parameters)
3. [State Variables and Internal Stress Calculation](#state-variables-and-internal-stress-calculation)
4. [Stress Calculation Procedure](#stress-calculation-procedure)
5. [Example Usage](#example-usage)
6. [Contributing](#contributing)

---

## Overview

This UMAT subroutine calculates the stress response of a viscoelastic material defined by a series of relaxation terms (Prony series). The model separates the total deformation into:
- **Volumetric (elastic) components**: Managed through bulk modulus (`K0`).
- **Deviatoric (viscoelastic) components**: Managed through the Prony series parameters (`G_i` and `TAU_i`), which govern the time-dependent relaxation response.

The subroutine utilizes state variables to track and update the internal stress state across time steps, allowing for an accurate simulation of viscoelastic behavior.

## Material Properties and Parameters

### Input Properties (PROPS array)

- **PROPS(1)** - `E0`: Instantaneous Young's Modulus
- **PROPS(2)** - `NU`: Poisson's Ratio
- **PROPS(3)** - `NPT`: Number of Prony series terms
- **PROPS(4:end)** - Alternating `G_i` (relative modulus) and `TAU_i` (relaxation time) values for each Prony term

### State Variables (STATEV array)

- The state variables store the stress-like internal variables for each Prony term. Each term has six components (to accommodate three normal and three shear directions in deviatoric stress). The state variables are used to keep track of viscoelastic history and relaxation over time.

## State Variables and Internal Stress Calculation

In the UMAT code, three key variables play an essential role in managing the time-dependent stress response:

- **`SM_OLD`**: This is the array of stress-like internal variables from the previous time increment. It represents the deviatoric stress contributions for each Prony term at the previous step, providing the "memory" of the stress state.
  
- **`SM`**: This is the array of updated stress-like internal variables for the current increment. It calculates the new viscoelastic stress contributions by adding the strain increment to the previous state stored in `SM_OLD`.

- **`SM_DOT`**: This array applies time-dependent relaxation to the updated stresses, simulating the decay of viscoelastic stresses according to the Prony series. The relaxed stress values in `SM_DOT` are stored back in `STATEV` for use in the next increment.

## Stress Calculation Procedure

The subroutine calculates the stress at each increment through the following steps:

1. **Initialize Previous State (`SM_OLD`)**: 
   - Load the previous increment’s stress state stored in `STATEV` into `SM_OLD`. This serves as the initial condition for calculating the new stress state.

2. **Calculate Current Stress Contribution (`SM`)**:
   - Using the deviatoric strain increment, the subroutine updates the stress-like internal variables. This calculation is based on both the strain increment and Prony series parameters (`G_i`, `TAU_i`), added to the previous state `SM_OLD`.

3. **Apply Relaxation (`SM_DOT`)**:
   - Each stress component in `SM` is decayed based on its relaxation time using the formula:
     \[
     \text{SM\_DOT(N, I)} = \text{EXP(-DTIME/TAU\_i(N))} \cdot \text{SM(N, I)}
     \]
   - This exponential decay simulates the natural relaxation of viscoelastic stress over the time increment `DTIME`, reducing the stress contribution as time progresses.

4. **Update Internal State (`STATEV`)**:
   - Store `SM_DOT` in `STATEV`, so that these relaxed stresses become the initial state (`SM_OLD`) in the next time increment.

5. **Compute Total Stress**:
   - The subroutine combines different components of stress:
     - **Hydrostatic Stress**: Calculated using the bulk modulus for the volumetric part.
     - **Long-term Elastic Stress**: Represents the response that remains after all relaxation occurs.
     - **Viscoelastic History Stress**: Summed from the contributions of each Prony term's internal stress state.
     - **Instantaneous Stress**: Immediate response to the deviatoric strain increment.

   - These components are combined into the total stress tensor, incorporating all time-dependent viscoelastic effects.

### Stress Calculation Equation

For normal stress components (1, 2, 3), the total stress is computed as:
\[
\text{STRESS(I)} = \text{STRESS\_HYDROSTATIC(I)} + \text{STRESS\_LONGTERM(I)} + \text{STRESS\_HISTORY(I)} + \text{STRESS\_INSTANTANEOUS(I)}
\]
For shear components (4, 5, 6), which do not include a volumetric contribution:
\[
\text{STRESS(I)} = \text{STRESS\_LONGTERM(I)} + \text{STRESS\_HISTORY(I)} + \text{STRESS\_INSTANTANEOUS(I)}
\]

## Example Usage

To use this UMAT in an ABAQUS analysis, specify the UMAT file and assign it to the material with the appropriate `PROPS` values. Ensure that `NPT`, `G_i`, and `TAU_i` values align with your Prony series terms.

Example input for ABAQUS:
```plaintext
*MATERIAL, NAME=ViscoelasticMaterial
*USER MATERIAL, CONSTANTS=10
5000.0, 0.3, 3, 0.2, 5.0, 0.1, 2.0, 0.05, 0.8
