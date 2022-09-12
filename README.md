# Hamilton-Jacobi Equations: Inverse Design
Our goal here is to study the inverse design problem associated to Hamilton Jacobi Equations (HJ)

$$ \partial_t u + H(\nabla_x u) = 0 $$

starting from an initial condition $u_0\in C^{0,1}(\mathbb{R}^n)$ and for a given superlinear convex Hamiltonian $H: \mathbb{R}^n\rightarrow\mathbb{R}$. 

More precisely, for a given target function $u_T$ and a time horizon $T>0$, we want to construct all the 
initial conditions $u_0$ such that the viscosity solution of (HJ) coincides with a target function $u_T$ at time $T$.

The study of this problem can also be motivated by considering the following question:

Given an observation of the solution to (HJ) at time $T>0$, can we construct all the possible initial data that agree we the observation at time $T$?

For this purpose, for a fixed $T>0$, we define the following nonlinear operator, which associates to any initial condition $u_0 \in C^{0,1}(R^n)$, 
the function $S_T^+ u_0 = u(T, \cdot ) \in C^{0,1}(R^n)$, where $u$ is the viscosity solution of (HJ).

The inverse design problem that we are considering is then reduced to, for $u_T$ and $T>0$ given, characterize all the initial conditions $u_0$ 
satisfying $S^+_{T} u_0= u(T, \cdot )$.

## Installation

Run one of the files Example1, Example2, Example3, or Example4. 

Example 1 requires Matlab's Statistics and Machine Learning Toolbox. 
