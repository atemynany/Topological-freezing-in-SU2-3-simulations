# Theoretical Framework: Topological Freezing and Error Analysis in Lattice QCD

## 1. Topological Freezing and Critical Slowing Down

In Lattice Quantum Chromodynamics (LQCD), topological freezing is a critical breakdown of ergodicity in Markov Chain Monte Carlo (MCMC) algorithms as the system approaches the continuum limit ($a \to 0$, $\beta \to \infty$).

While the configuration space of an $SU(N_c)$ lattice gauge theory is simply connected, the action barriers separating different topological sectors $Q$ diverge in the continuum. Consequently, algorithms based on continuous Hamiltonian dynamics (like Hybrid Monte Carlo) or local updates experience an exponential increase in the topological tunneling time.

The integrated autocorrelation time of the topological charge, $\tau_{\text{int}}(Q)$, follows the scaling law:


$$\tau_{\text{int}}(Q) \propto a^{-z}$$


where $z$ is the dynamic critical exponent. In pure gauge $SU(3)$, $z \approx 5$ has been observed.

* **$SU(3)$ Threshold:** Freezing becomes severe at $\beta \ge 6.4$ ($a \lesssim 0.05 \text{ fm}$).
* **Mitigation:** **Open Boundary Conditions (OBC)** in the temporal direction are used to allow the topological charge to flow through the boundaries, avoiding the absolute divergence of $\tau_{\text{int}}$.

## 2. Scale Setting and Interpolation

To relate the dimensionless coupling $\beta$ to physical results, we use interpolation formulas for the lattice spacing $a$.

### $SU(3)$ Scale Setting (Sommer Scale)

For pure $SU(3)$, the Sommer scale $r_0 \approx 0.5 \text{ fm}$ is used. The Necco-Sommer formula for $5.7 \le \beta \le 6.92$ is:


$$\ln\left(\frac{a}{r_0}\right) = -1.6804 - 1.7139(\beta - 6) + 0.8155(\beta - 6)^2 - 0.6667(\beta - 6)^3$$

### $SU(2)$ Scale Setting (Gradient Flow Scale)

For $SU(2)$, the gradient flow scale $t_0$ is preferred (defined where $t^2 \langle E(t) \rangle = 0.1$). For $\beta \approx 2.4 - 2.8$:


$$\ln\left(\frac{t_0}{a^2}\right) = 1.285 + 6.409(\beta - 2.600) - 0.7411(\beta - 2.600)^2$$

## 3. The Dynamical Correlation Coefficient ($C_F$)

When topology is frozen ($\tau_{exp} \to \infty$), we must determine if other physical observables are still reliable. This is done via $C_F$.

**Definition:**


$$C_F = \lim_{t\rightarrow\infty}\rho_F(t)e^{t/\tau_{exp}}$$

* **$\rho_F(t)$**: Normalized autocorrelation of observable $F$.
* **$\tau_{exp}$**: Exponential autocorrelation time of the slowest mode (topological charge).

**Application:**
If $C_F \ll 1$, the observable $F$ is **decoupled** from the topological freezing. Its statistical errors remain manageable even if the topological charge $Q$ never changes during the simulation.

## 4. Algorithmic Generation and Measurement

### Heatbath Algorithm

$SU(3)$ configurations are generated using the **Cabibbo-Marinari** decomposition, which updates the matrix via three $SU(2)$ subgroups. For each subgroup, the **Kennedy-Pendleton** algorithm samples the distribution:


$$dP(a_0) \propto \sqrt{1 - a_0^2} \exp(2\alpha a_0) da_0$$


This uses a rejection sampling method with a Gaussian envelope to maintain exactness and efficiency.

### Gradient Flow and Stability

To measure $Q$, configurations are "flowed" to remove UV noise. While standard **Wilson Flow** can be unstable for topology, over-improved actions like **DBW2** or **Iwasaki** stabilize the 1-instanton configuration because their dimension-6 operators satisfy the stability criterion $1 + 12c_1 < 0$.

---

### Primary References

1. **Schaefer, S., Sommer, R., & Virotta, F. (2011).** *Critical slowing down and error analysis in lattice QCD simulations.* Nucl. Phys. B 845. [arXiv:1009.5228]
2. **Tanizaki, Y., Tomiya, A., & Watanabe, H. (2024).** *Lattice gradient flows (de-)stabilizing topological sectors.* [arXiv:2411.14812]
3. **Lüscher, M., & Schaefer, S. (2011).** *Lattice QCD without topology barriers.* JHEP 07.
4. **Kennedy, A. D., & Pendleton, B. J. (1985).** *Improved heat bath method.* Phys. Lett. B 156.