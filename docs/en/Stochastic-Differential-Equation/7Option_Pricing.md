# Chapter 7: Option Pricing and Stopping Times

In the financial applications of Stochastic Differential Equations (SDEs), the Black-Scholes option pricing model stands as a classic milestone. In this chapter, we will utilize Itô's integral to derive the Black-Scholes Partial Differential Equation (PDE), explore the rigorous definition and properties of Stopping Times within the context of Itô integration, and finally introduce the generators of diffusion processes and the Feynman-Kac formula, establishing a profound connection between SDEs and PDEs.

---

## 1. Black-Scholes Option Pricing Model

We consider a standard European Call Option with a maturity date $T$ and a strike price $K$. At the maturity date, the payoff of the option is:

\[
V(S_T, T) = (S_T - K)^+ = \max(S_T - K, 0)
\]

### 1.1 Basic Assumptions

Before deriving the model, we establish the following assumptions based on the Black-Scholes framework:

1.  **Dynamics of the Underlying Asset**: The stock price $S_t$ follows Geometric Brownian Motion (GBM):

    \[
    \frac{dS_t}{S_t} = \mu dt + \sigma dW_t
    \]

2.  **Market Environment**: There exists a risk-free interest rate $r$ (compounded continuously), and the underlying stock **does not pay dividends**.
3.  **Transaction Friction**: There are no transaction costs or taxes, and the underlying asset is perfectly divisible and allows for short selling.
4.  **No-Arbitrage Principle**: There are no risk-free arbitrage opportunities in the market.

### 1.2 Derivation of the Black-Scholes PDE

Based on the principle of no-arbitrage pricing and the self-financing strategy, we can construct a risk-free portfolio to derive the partial differential equation that the option price $V(S_t, t)$ must satisfy.

!!! abstract "Theorem: Black-Scholes Partial Differential Equation (BS PDE)"

    Assume the option price $V(S, t)$ is a sufficiently smooth function of $S$ and $t$. Then $V$ must satisfy the following parabolic partial differential equation:

    \[
    \partial_t V + \frac{1}{2}\sigma^2 S^2 \partial_S^2 V + rS \partial_S V - rV = 0
    \]

    **Terminal Condition (for a call option)**:

    \[
    V(S, T) = (S - K)^+
    \]

??? proof "Derivation: Itô's Lemma and Delta Hedging (Infinitesimal Method) (Click to expand)"

    **Step 1: Application of Itô's Lemma**

    Given the dynamics of the underlying asset $dS_t = \mu S_t dt + \sigma S_t dW_t$, apply Itô's Lemma to the option price function $V(S_t, t)$:

    \[
    dV(S_t, t) = \partial_t V dt + \partial_S V dS_t + \frac{1}{2} \partial_S^2 V (dS_t)^2
    \]

    Substituting $(dS_t)^2 = \sigma^2 S_t^2 dt$, we obtain the differential of $V_t$:

    \[
    dV_t = \left( \partial_t V + \mu S_t \partial_S V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V \right) dt + \sigma S_t \partial_S V dW_t
    \]

    **Step 2: Constructing a Risk-Free Portfolio (Delta Hedging)**

    Construct a portfolio $\Pi_t$ consisting of 1 long position in the option and $\Delta$ short positions in the stock:

    \[
    \Pi_t = V_t - \Delta S_t
    \]

    Within an infinitesimal time interval $dt$, the change in the value of this portfolio (using the self-financing infinitesimal method) is:

    \[
    d\Pi_t = dV_t - \Delta dS_t
    \]

    Expanding $dV_t$ and $dS_t$ and substituting them in:

    \[
    d\Pi_t = \left( \partial_t V + \mu S_t \partial_S V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V - \Delta \mu S_t \right) dt + (\sigma S_t \partial_S V - \Delta \sigma S_t) dW_t
    \]

    **Step 3: Eliminating the Random Term and No-Arbitrage Pricing**

    To make the portfolio risk-free (i.e., eliminating the $dW_t$ term), we choose the stock holding $\Delta$ to be:

    \[
    \Delta = \partial_S V
    \]

    Substituting this into the equation, the random terms cancel out, yielding:

    \[
    d\Pi_t = \left( \partial_t V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V \right) dt
    \]

    According to the **no-arbitrage principle**, the return on a risk-free portfolio must equal the risk-free interest rate $r$. Thus, $d\Pi_t = r \Pi_t dt$:

    \[
    \left( \partial_t V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V \right) dt = r (V_t - S_t \partial_S V) dt
    \]

    Simplifying and rearranging the terms leads to the Black-Scholes Partial Differential Equation:

    \[
    \partial_t V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V + rS \partial_S V - rV = 0
    \]

    The derivation is complete. $\square$

---

## 2. Solving the Black-Scholes PDE

For the BS PDE, we can use a series of variable transformations to convert it into the classic Heat Equation, allowing us to find the analytical solution using the Gaussian kernel.

### 2.1 Transformation to the Heat Equation

??? proof "Solution Process: Heat Equation Transformation and Analytical Derivation (Click to expand)"

    **Step 1: Log-Transformation and Time Reversal**

    Let $\tau = T - t$ (time remaining until maturity), and $x = \ln S \Rightarrow S = e^x$.

    Define $v(x, \tau) = V(e^x, T - \tau)$. Applying the chain rule and substituting into the original PDE, the equation simplifies into a parabolic equation with constant coefficients:

    \[
    \partial_\tau v - \frac{1}{2}\sigma^2 \partial_x^2 v - \left(r - \frac{1}{2}\sigma^2\right) \partial_x v + r v = 0
    \]

    The terminal condition becomes the initial condition: $v(x, 0) = (e^x - K)^+$.

    **Step 2: Eliminating Lower-Order Derivative Terms (Variable Substitution Method)**

    To eliminate the first-order spatial derivative term $\partial_x v$ and the constant term $rv$, we perform an exponential variable substitution by setting:

    \[
    v(x, \tau) = u(x, \tau) e^{\alpha \tau + \beta x}
    \]

    Calculating partial derivatives and substituting back into the equation:

    \[
    \partial_\tau v = (\partial_\tau u + \alpha u) e^{\alpha \tau + \beta x}
    \]

    \[
    \partial_x v = (\partial_x u + \beta u) e^{\alpha \tau + \beta x}
    \]

    \[
    \partial_x^2 v = (\partial_x^2 u + 2\beta \partial_x u + \beta^2 u) e^{\alpha \tau + \beta x}
    \]

    Substituting these and rearranging, setting the coefficients of $\partial_x u$ and $u$ to zero allows us to solve for specific values of $\alpha$ and $\beta$:

    \[
    \begin{cases}
    \sigma^2 \beta + \left(r - \frac{1}{2}\sigma^2\right) = 0 \\
    \alpha - \frac{1}{2}\sigma^2 \beta^2 - \left(r - \frac{1}{2}\sigma^2\right)\beta + r = 0
    \end{cases}
    \]

    Solving these gives:

    \[
    \beta = -\frac{r}{\sigma^2} + \frac{1}{2}, \quad \alpha = -\frac{1}{2\sigma^2}\left(r - \frac{\sigma^2}{2}\right)^2 - r
    \]

    At this point, the original equation is transformed into the standard Heat Equation (Diffusion Equation):

    \[
    \partial_\tau u - \frac{1}{2}\sigma^2 \partial_x^2 u = 0
    \]

    **Step 3: Solving Using the Green's Function (Gaussian Kernel)**

    The solution to the Heat Equation can be obtained by integrating its fundamental solution (Kernel) with the initial conditions:

    \[
    u(x, \tau) = \int_{-\infty}^{+\infty} K(x-\xi, \tau) \psi(\xi) d\xi
    \]

    Where the kernel function is:

    \[
    K(x-\xi, \tau) = \frac{1}{\sigma \sqrt{2\pi \tau}} \exp\left( -\frac{(x-\xi)^2}{2\sigma^2 \tau} \right)
    \]

### 2.2 Analytical Derivation via Gaussian Kernel

??? proof "Deep Dive: From Gaussian Kernel to the Final Pricing Formula (Click to expand)"

    **Step 1: Writing the General Solution via Green's Function**

    The solution to the heat equation in full space can be obtained through the convolution of the initial condition $\psi(\xi)$ and the Gaussian kernel $K$:

    \[
    u(x, \tau) = \int_{-\infty}^{+\infty} \frac{1}{\sigma \sqrt{2\pi \tau}} \exp\left( -\frac{(x-\xi)^2}{2\sigma^2 \tau} \right) \psi(\xi) d\xi
    \]

    **Step 2: Substituting Option Initial Conditions**

    Substitute $\psi(\xi) = (e^\xi - K)^+ e^{-\beta \xi}$. Since $(e^\xi - K)^+$ is non-zero only when $\xi > \ln K$, the integration interval is narrowed:

    \[
    u(x, \tau) = \frac{1}{\sigma \sqrt{2\pi \tau}} \int_{\ln K}^{\infty} e^{-\frac{(x-\xi)^2}{2\sigma^2 \tau}} (e^\xi - K) e^{-\beta \xi} d\xi
    \]

    **Step 3: Splitting the Integral and Completing the Square**

    Split the integral into $I_1$ (containing the $e^\xi$ term) and $I_2$ (containing the $K$ term):

    1. For $I_1$, the exponential terms are combined as $-\frac{(\xi - [x + \sigma^2 \tau(1-\beta)])^2}{2\sigma^2 \tau} + \text{Const}$.
    2. For $I_2$, the exponential terms are combined as $-\frac{(\xi - [x - \sigma^2 \tau \beta])^2}{2\sigma^2 \tau} + \text{Const}$.

    **Step 4: Variable Substitution and Normalization**

    By setting $z = \frac{\xi - \text{center}}{\sigma \sqrt{\tau}}$, the integrals are transformed into the form of the standard normal cumulative distribution function $N(d)$. The lower bound $\ln K$ translates to values representing $d_1$ and $d_2$ after transformation:

    \[
    d_1 = \frac{\ln(S/K) + (r + \sigma^2/2)\tau}{\sigma\sqrt{\tau}}, \quad d_2 = d_1 - \sigma\sqrt{\tau}
    \]

    **Step 5: Returning to Original Variables**

    Finally, using $V = u e^{-(\alpha \tau + \beta x)}$ to revert to original variables, the exponential terms and constant terms outside the integral cancel out precisely, resulting in the Black-Scholes formula. $\square$

**Final Conclusion: Black-Scholes Pricing Formula**

The analytical solution for a European Call Option is:

\[
V(S,t) = S N(d_1) - K e^{-r(T-t)} N(d_2)
\]

Where:

\[
d_1 = \frac{\ln(S/K) + (r + \frac{1}{2}\sigma^2)(T-t)}{\sigma \sqrt{T-t}}
\]

\[
d_2 = d_1 - \sigma \sqrt{T-t}
\]

Here, $N(\cdot)$ denotes the cumulative distribution function of the standard normal distribution.

---

## 3. Stopping Times and Their Properties

When studying complex financial derivatives (such as American options, which allow for early exercise), the traditional fixed time $t$ is insufficient, and we must introduce the concept of **Stopping Times**.

### 3.1 Definition and Basic Properties of Stopping Times

Consider $(\Omega, \mathcal{F}, \{\mathcal{F}_t\}_{t \ge 0}, P)$ with a given filtration.

!!! abstract "Definition: Stopping Time"

    A random variable $\tau: \Omega \to [0, +\infty]$ is called a **Stopping Time** if for any given time $t \ge 0$, the event $\{\tau \le t\}$ is $\mathcal{F}_t$-measurable.

    In other words: At time $t$, given the information $\mathcal{F}_t$ available up to that point, we can definitively determine whether the process has "stopped" ($\tau \le t$).

**Common Properties of Stopping Times:**

Assume $\tau_1$ and $\tau_2$ are stopping times, then:

1.  **Equivalent Condition**: $\{\tau < t\} \in \mathcal{F}_t$.
2.  **Minima are Stopping Times**: $\min\{\tau_1, \tau_2\} = \tau_1 \wedge \tau_2$ is a stopping time.
3.  **Maxima are Stopping Times**: $\max\{\tau_1, \tau_2\} = \tau_1 \vee \tau_2$ is a stopping time.

*(Example: The time $\tau = \inf\{t > 0: S_t \ge H\}$ at which a stock price first hits a barrier $H$ is a typical stopping time.)*

### 3.2 Itô Integration with Stopping Times

We can extend the upper limit of an Itô integral to a stopping time $\tau$:

\[
\int_0^\tau G_s dW_s \triangleq \int_0^t \chi_{\{s \le \tau\}} G_s dW_s
\]

Where $\chi$ is the indicator function. Itô integrals with stopping times inherit standard martingale properties and isometry:

1.  **Zero Expectation**:

    \[
    E \left[ \int_0^\tau G_s dW_s \right] = 0
    \]

2.  **Itô Isometry**:

    \[
    E \left[ \left( \int_0^\tau G_s dW_s \right)^2 \right] = E \left[ \int_0^\tau G_s^2 ds \right]
    \]

---

## 4. Elliptic Equations and Probabilistic Representation

In the study of Brownian motion, the solution to a partial differential equation (PDE) can often be expressed as the expectation of a certain stochastic functional. We set $X_t$ as an $n$-dimensional standard Brownian motion (i.e., drift $b=0$, diffusion $B=I$), in which case its generator is $L = \frac{1}{2}\Delta$.

### 4.1 Poisson Equation and Exit Time

!!! abstract "Proposition: Expectation Representation of Exit Time"

    Consider the Poisson equation boundary value problem on a region $\Omega$:

    \[
    \begin{cases}
    -\frac{1}{2}\Delta u(x) = 1, & x \in \Omega \\
    u(x) = 0, & x \in \partial \Omega
    \end{cases}
    \]

    Then the solution to this equation can be expressed as the expected time for a Brownian motion starting from point $x$ to first exit the region $\Omega$:

    \[
    u(x) = E^x[\tau_\Omega]
    \]

??? proof "Proof: Using Dynkin's Formula (Click to expand)"

    **Step 1: Apply Dynkin’s Formula (Expectation form of Itô’s Formula)**

    For $u \in C^2(\Omega) \cap C(\bar{\Omega})$, according to Itô’s formula:

    \[
    u(X_{\tau \wedge t}) = u(X_0) + \int_0^{\tau \wedge t} \frac{1}{2}\Delta u(X_s) ds + \int_0^{\tau \wedge t} \nabla u(X_s) \cdot dW_s
    \]

    Taking the expectation $E^x$ starting from point $x$. Since the Itô integral term is a martingale, its expectation is 0:

    \[
    E^x[u(X_{\tau \wedge t})] = u(x) + E^x \left[ \int_0^{\tau \wedge t} \frac{1}{2}\Delta u(X_s) ds \right]
    \]

    **Step 2: Substitute Equation Conditions**

    Given $\frac{1}{2}\Delta u = -1$. Substituting this into the above expression:

    \[
    E^x[u(X_{\tau \wedge t})] = u(x) + E^x \left[ \int_0^{\tau \wedge t} (-1) ds \right] = u(x) - E^x[\tau \wedge t]
    \]

    **Step 3: Take the Limit as $t \to \infty$**

    According to the boundary conditions, when the Brownian motion hits the boundary, $u(X_\tau) = 0$. By the Dominated Convergence Theorem:

    \[
    0 = u(x) - E^x[\tau_\Omega] \implies u(x) = E^x[\tau_\Omega]
    \]

    The proof is complete. $\square$

---

### 4.2 Probabilistic Representation of Harmonic Functions

!!! abstract "Proposition: Probabilistic Solution to the Dirichlet Problem"

    Consider the boundary value problem for the Laplace equation (harmonic functions):

    \[
    \begin{cases}
    \Delta u(x) = 0, & x \in \Omega \\
    u(x) = g(x), & x \in \partial \Omega
    \end{cases}
    \]

    The solution $u(x)$ has the following expectation representation:

    \[
    u(x) = E^x[g(X_{\tau_\Omega})]
    \]

??? proof "Proof (Click to expand)"

    **Step 1: Use Harmonic Property**

    Since $\Delta u = 0$, we have $L u = 0$. According to Itô’s formula, the process $u(X_t)$ is a local martingale before exiting $\Omega$.

    **Step 2: Apply Optional Sampling Theorem**

    For the stopping time $\tau_\Omega$, we have:

    \[
    u(x) = E^x[u(X_0)] = E^x[u(X_{\tau_\Omega})]
    \]

    **Step 3: Combine with Boundary Conditions**

    When Brownian motion hits the boundary, $X_{\tau_\Omega} \in \partial \Omega$, at which point $u(X_{\tau_\Omega}) = g(X_{\tau_\Omega})$. Substituting this in:

    \[
    u(x) = E^x[g(X_{\tau_\Omega})]
    \]

    This indicates that the value of a harmonic function at a point is equal to the stochastic average of its boundary values. $\square$

---

### 4.3 Dirichlet Problem with Discontinuous Boundary Values

The notes mention a scenario where the boundary is divided into two parts, $\Gamma_1$ and $\Gamma_2$.

!!! abstract "Proposition: Indicator Functions and Hitting Probabilities"

    Consider the following boundary value problem:

    \[
    \begin{cases}
    \Delta u = 0, & x \in \Omega \\
    u = 1, & x \in \Gamma_1 \\
    u = 0, & x \in \Gamma_2
    \end{cases}
    \]

    Where $\partial \Omega = \Gamma_1 \cup \Gamma_2$.

??? proof "Proof: Intuitive Understanding of Probabilistic Representation (Click to expand)"

    Based on the conclusion in 4.2, the solution can be written as the expectation of boundary values. Let the boundary function $g(x) = I_{\Gamma_1}(x)$, which is 1 on $\Gamma_1$ and 0 otherwise. Then:

    \[
    u(x) = E^x[I_{\Gamma_1}(X_{\tau_\Omega})] = P^x(X_{\tau_\Omega} \in \Gamma_1)
    \]

    **Physical Meaning**: The solution $u(x)$ at that point is equal to the probability that a Brownian motion starting from that point first hits the boundary within the $\Gamma_1$ region. $\square$

---

## 5. Feynman-Kac Formula

The Feynman-Kac formula establishes a link between partial differential equations with "potential terms" or "discounting terms" and the exponential functionals of stochastic processes.

!!! abstract "Theorem: Feynman-Kac Formula"

    Consider the following boundary value problem for a second-order partial differential equation:

    \[
    \begin{cases}
    -\frac{1}{2}\Delta u(x) + q(x)u(x) = f(x), & x \in \Omega \\
    u(x) = 0, & x \in \partial \Omega
    \end{cases}
    \]

    Where $q(x) \ge 0$. The probabilistic representation of the solution $u(x)$ is:

    \[
    u(x) = E^x \left[ \int_0^{\tau_\Omega} f(X_t) \exp\left( -\int_0^t q(X_s) ds \right) dt \right]
    \]

??? proof "Proof: Using Itô Expansion of an Auxiliary Process (Click to expand)"

    **Step 1: Construct an Auxiliary Process**

    Define a process $M_t$ as follows:

    \[
    M_t = u(X_t) \exp\left( -\int_0^t q(X_s) ds \right) + \int_0^t f(X_s) \exp\left( -\int_0^s q(X_r) dr \right) ds
    \]

    **Step 2: Apply Itô’s Formula to Find the Differential**

    Let $e_t = \exp\left( -\int_0^t q(X_s) ds \right)$. Using the product rule:

    \[
    d(u(X_t) e_t) = e_t du(X_t) + u(X_t) de_t
    \]

    Where:

    \[
    du(X_t) = \nabla u \cdot dW_t + \frac{1}{2}\Delta u dt
    \]

    \[
    de_t = -q(X_t) e_t dt
    \]

    Substituting these in:

    \[
    d(u(X_t) e_t) = e_t \left( \frac{1}{2}\Delta u - q(X_t)u(X_t) \right) dt + e_t \nabla u \cdot dW_t
    \]

    **Step 3: Combine with PDE Conditions**

    Since $-\frac{1}{2}\Delta u + qu = f \implies \frac{1}{2}\Delta u - qu = -f$. Substitute this into $dM_t$:

    \[
    dM_t = d(u(X_t) e_t) + f(X_t)e_t dt = \left[ e_t(-f) + e_t f \right] dt + e_t \nabla u \cdot dW_t
    \]

    Therefore, $dM_t = e_t \nabla u \cdot dW_t$, indicating that $M_t$ is a local martingale.

    **Step 4: Take Expectation and Boundary Conditions**

    By the optional sampling theorem, $u(x) = M_0 = E^x[M_{\tau_\Omega}]$.

    When $t = \tau_\Omega$, $u(X_{\tau_\Omega}) = 0$, so the first term disappears, leaving:

    \[
    u(x) = E^x \left[ \int_0^{\tau_\Omega} f(X_s) \exp\left( -\int_0^s q(X_r) dr \right) ds \right]
    \]

    The proof is complete. $\square$