# SDRG

In the disordered spin systems, SDRG (strong-disordered renormalization group) was developed to solve for the properties of ground state and low energy excitations. The basic idea of the SDRG is to find the system’s ground state by successively eliminating degrees of freedom with the highest energy, then generating effective Hamiltonians with fewer degrees of freedom and lower energy. This technique was originally showed by Dasgupta and Ma to study the ground state and low energy behavior of the disordered Heisenberg chain. Later, Daniel Fisher solved the RG equation of XX chain and showed that RG scheme gives asymptotically exact results for the ground state properties of these disordered Heisenberg chains, so called infinite randomness fixed point. 

In this project, disordered spin-1/2 antiferromagnetism Heisenberg chain flows toward a "random-singlet phase", described by an infinite randomness fixed point, for any amount of disorder. Otherwise, I'm going to study the Heisenberg-liked model with multispins coupling named Q model, do some comparison between these models and figure the difference out.

Some publication of mine:
https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.174442

Finally, we compared our SDRG results with Projector Monte Carlo Method. 


In this paper, we study the quantum critical behavior of strongly disordered quantum spin chains using SDRG, and compare their results to exact Quantum Monte Carlo Method results. Interestingly, we find numerically that the random singlet phase is characterized by logarithmic corrections that are missed by SDRG, which is broadly believed to yield exact results for this phase.

The asymptotic exactness of the SDRG method precisely means that universal properties do not depend on the initial disorder distribution — in general, provided the initial disorder strength is strong enough. In addition, for the Heisenberg chain, it is well known that weak disorder is a relevant perturbation to the clean fixed point, and that any initial disorder distribution should flow to the random singlet fixed point at large distances. There is no doubt in the community that the random bond Heisenberg chain is indeed in a random singlet phase for which SDRG predictions should be exact, and that the QMC results presented in this work largely confirm this.
