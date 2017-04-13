# SDRG

In the disordered spin systems, SDRG (strong-disordered renormalization group) was developed to solve for the properties of ground state and low energy excitations. The basic idea of the SDRG is to find the system’s ground state by successively eliminating degrees of freedom with the highest energy, then generating effective Hamiltonians with fewer degrees of freedom and lower energy. This technique was originally showed by Dasgupta and Ma to study the ground state and low energy behavior of the disordered Heisenberg chain. Later, Daniel Fisher solved the RG equation and showed that RG scheme gives asymptotically exact results for the ground state property of these disordered Heisenberg chains, so called infinite randomness fixed point. 
In this project, disordered spin-1/2 antiferromagnetism Heisenberg chain flows toward a "random-singlet phase", described by an infinite fixed point, for any amount of disorder. Otherwise, I'm going to study the Heisenberg-liked model with multispins coupling named Q model, do some comparison between these models and figure the difference out.

Some publication of mine:
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.174442

* The Hamiltonian of the original 1D Antiferromagetic Heisenberg Model :

$$\hat{\bf{H}}=-{J_{i}}\sum_{i}\hat{\bf{P}}_{i,i+1}$$


$\hat{\bf{P}}_{i,i+1} \equiv \frac{1}{4}-\vec{S}{i}\cdot\vec{S}{i+1}$


![equation](http://latex.codecogs.com/gif.latex?\hat{\bf{P}}_{i,i+1}\equiv\frac{1}{4}-\vec{S}_{i}\cdot\vec{S}_{i+1})
