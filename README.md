# Cholesky_decomp_matrix_inversion
Specialized fast, in-place matrix inversion for 5x5 symmetric positive definite 
matrices via LDL^T Cholesky decomposition.

This repo contains 2 functions that do the same thing, one is a hand unrolled 
copy of the other. See my articles on my blog [here](https://silicon42.github.io/blog/inverting_matrices_in_place) and [here](https://silicon42.github.io/blog/inverting_matrices_in_place_addenda) for a more in-depth 
explanation.

The TLDR is that since a 5x5 symmetric matrix can be represented in 
15 values, it can fit fully inside a single cache line and, depending on hardware,
may even fit entirely into registers. And using the LDL^T variant of Cholesky 
decomposition it can be quickly computed in place with minimal usage ofexpensive 
operations.

Below is a diagram of the element-wise dependencies for the Cholesky decomposition,
with the "S" column being the intermediate calculations and the "L or D" column 
being the final decomposed values, they still represent the SAME slots. The 
diagram progresses mostly right to left, bottom to top because that's just how I 
drew it and I'm not drawing it again just to fix that. Draw.io is nice but it's 
not the most convenient to use in my opinion. I made it to help myself and am
posting here so others might also benefit

![cholesky decomposition dependency graph](./cholesky%20decomp%20dependency.svg)

An arrow pointing into a box means that box depends on the box the arrow came from. Same colored arrows are required at the same time, or in other words they take part in the same calculation. All values in S have an implicit dependency on the "A" coefficients at their own position as that's what they are initialized to. Dashed arrows represent a bundle of arrows of multiple colors with the last required being the primary color. A box is ready to fulfill any outgoing dependencies that rely on it onces all of it's incoming arrows are accounted for. Boxes with multiple arrow pairs keep intermediate sum results in them and are not usable for later dependencies until all incoming arrows have been accounted for. Once a box in the "L or D" column has been written to, the corresponding box in the "S" column can no longer be used because it has been overwritten since it is the same slot, hence the need for a temporary variable.