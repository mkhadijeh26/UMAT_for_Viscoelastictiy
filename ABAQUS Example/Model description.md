# Viscoelastic Material Model

This input file provided describes the numerical solution for a multiple loading segment process of a viscoelastic material

## Model Description

The model simulates the stress response of a viscoelastic material subjected to a multi-step strain loading process. It uses a two-term Prony series to represent the material's relaxation function.

### Loading Process

The loading process consists of four steps:
1. Loading to ε = ε1
2. Holding the load for 50 seconds
3. Unloading to ε = ε2
4. Holding for another 50 seconds

## Implementation

The model is implemented in MATHCAD, with the following key components:

1. Definition of material constants (E, P1, P2, τ1, τ2)
2. Specification of loading times
3. Definition of the piecewise strain function
4. Calculation of stress for each loading step using hereditary integrals
5. Compilation of results into a piecewise stress function
6. Output of results (time, strain, and stress) to a file

## Output

The model generates three columns of data:
1. Time
2. Strain
3. Stress

These results can be used to plot stress-time and stress-strain curves, providing insights into the material's viscoelastic behavior under the specified loading conditions.

## Usage

This model can be useful for simulating and analyzing the behavior of viscoelastic materials under complex loading conditions.