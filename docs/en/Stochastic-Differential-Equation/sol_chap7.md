---
tags:
  - Stochastic Differential Equations
  - Exercise Solutions
---

# 📝 Exercise Solutions: Option Pricing and Stopping Times

!!! abstract "About This Page"
    This page contains the core exercise solutions for Chapter 7 of the *Stochastic Differential Equations* course, focusing on stopping time theory, hitting probabilities of high-dimensional Brownian motion, and European option pricing (Black-Scholes model). The content covers rigorous proofs for stopping times, deriving Dirichlet boundary value problems using martingale properties, and the complete derivation of the Put option formula using the fundamental solution of the heat equation (Gaussian kernel).

---

## Part I: Stopping Times and Hitting Probabilities

### Exercise 1
**Problem**
Let $W$ be a one-dimensional Brownian motion. Define $\tau = \inf\{t: W(t) \in (a, b]\}$ as the first hitting time of the interval $(a, b]$. Prove that $\tau$ is a stopping time.

??? success "Solution (Click to expand)"

    To rigorously prove that $\tau$ is a stopping time, based on the logic of the assignment notes, we decompose the first entry into the interval $(a, b]$ into a combination of two fundamental hitting times.
    
    First, define two fundamental first hitting times:

    \[
    \tau_a = \inf\{t: W(t) > a\}
    \]

    \[
    \tau_b = \inf\{t: W(t) \le b\}
    \]

    **Step 1: Prove the measurability of $\tau_a$ and $\tau_b$**

    For $\tau_b$ (first entry into the closed set $(-\infty, b]$): Since the paths of Brownian motion are continuous, the first hitting time of a closed set satisfies:

    \[
    \{\tau_b \le t\} = \inf_{0 \le s \le t} W(s) \le b
    \]

    Utilizing the density of the set of rational numbers $\mathbb{Q}$, this is equivalent to taking the infimum over the rational set $\mathbb{Q}$:

    \[
    \{\tau_b \le t\} = \bigcap_{n=1}^\infty \bigcup_{s \in \mathbb{Q} \cap [0, t]} \left\{ W_s \le b + \frac{1}{n} \right\}
    \]

    Since the right-hand side is a countable intersection and union of $\mathcal{F}_t$-measurable events, $\tau_b$ is a stopping time.

    For $\tau_a$ (first entry into the open set $(a, +\infty)$): Similarly, utilizing the properties of open sets and rational approximation:

    \[
    \{\tau_a < t\} = \bigcup_{s \in \mathbb{Q} \cap [0, t)} \{ W_s > a \} \in  \mathcal{F}_t
    \]

    Combined with the right-continuity of the filtration, $\tau_a$ is also a stopping time.

    **Step 2: Combinatorial logic to prove $\tau$ is a stopping time**
    
    Now we examine the target stopping time $\tau = \inf\{t: W(t) \in (a, b]\}$. We can decompose $\{\tau \le t\}$ based on the initial state $W(0)$:
    
    1. **If $W(0) \le a$**: The process must first cross $a$ to enter $(a, b]$. The event is equivalent to: $\{W(0) \le a\} \cap \{\tau_a \le t, \inf_{\tau_a \le s \le t} W(s) \le b\}$

    2. **If $W(0) > b$**: The process must first descend to cross $b$ to enter $(a, b]$. The event is equivalent to: $\{W(0) > b\} \cap \{\tau_b \le t\}$
       
    3. **If $W(0) \in (a, b]$**: In this case, $\tau = 0 \le t$ is necessarily true. The event is equivalent to $\{W(0) \in (a, b]\}$.

    Taking the union of these three mutually exclusive and exhaustive cases, we get:
    
    \[
    \{\tau \le t\} = \left( \{W(0) \le a\} \cap \{\tau_a \le t, \inf_{\tau_a \le s \le t} W_s \le b\} \right) \cup \left( \{W(0) > b\} \cap \{\tau_b \le t\} \right) \cup \{W(0) \in (a, b]\}
    \]
    
    Since all basic events involved (including the initial state sets and the level-crossing events of $\tau_a$ and $\tau_b$) are $\mathcal{F}_t$-measurable, and the $\sigma$-algebra is closed under finite unions and intersections, the overall event $\{\tau \le t\} \in \mathcal{F}_t$.
    
    Therefore, $\tau$ satisfies the rigorous definition of a stopping time. The proof is complete.

---

### Exercise 2
**Problem**
Let $W$ be an $n$-dimensional Brownian motion ($n \ge 3$). Denote $X = W + x_0$, where $x_0$ lies within the spherical shell $R_1 < |x| < R_2$. Calculate the probability that $X$ hits the outer sphere $|x| = R_2$ before hitting the inner sphere $|x| = R_1$.

??? success "Solution (Click to expand)"

    **Step 1: Derive the Dirichlet boundary value problem using Itô's formula and martingale properties**
    
    Let $\tau$ be the first exit time of process $X$ from the spherical shell region $D = \{x: R_1 < |x| < R_2\}$. Let $u(x)$ be the probability of starting from the initial point $x$ and hitting the outer sphere first, i.e., $u(x) = P^x(|X_\tau| = R_2)$.
    
    We construct a stochastic process $M_t = u(X_{t \wedge \tau})$. Applying the multi-dimensional Itô's formula to expand $u(X_t)$:

    \[
    u(X_{t \wedge \tau}) = u(X_0) + \int_0^{t \wedge \tau} \nabla u(X_s) dW_s + \int_0^{t \wedge \tau} \frac{1}{2} \Delta u(X_s) ds
    \]

    For the expected probability to remain conserved before the stopping time, the process $M_t = u(X_{t \wedge \tau})$ must be a **martingale**.
    
    The core property of a martingale requires that the expectation of the drift term (the $ds$ term) of its differential is constantly 0. Therefore, the integrand must satisfy:

    \[
    \frac{1}{2} \Delta u(x) = 0 \implies \Delta u(x) = 0 \quad (x \in D)
    \]

    For the boundary conditions: when the process hits the outer sphere, the event definitely occurs, so the probability is 1; when it hits the inner sphere, the event does not occur, so the probability is 0. This rigorously transforms the probability problem into the following Dirichlet boundary value problem:

    \[
    \begin{cases}
    \Delta u(x) = 0, & R_1 < |x| < R_2 \\
    u(x) = 0, & |x| = R_1 \\
    u(x) = 1, & |x| = R_2
    \end{cases}
    \]

    **Step 2: Solve the partial differential equation**

    Due to the symmetry of the region, we assume the solution has a radial form $u(x) = v(r)$, where $r = |x|$. For $n \ge 3$, the fundamental solution for the Laplace operator takes the form:

    \[
    v(r) = A + B r^{2-n}
    \]

    Substitute the boundary conditions:
    At $r = R_1$:

    \[
    A + B R_1^{2-n} = 0 \implies A = -B R_1^{2-n}
    \]

    At $r = R_2$:

    \[
    A + B R_2^{2-n} = 1
    \]

    Substituting the first equation into the second, we solve for the constants $B$ and $A$:

    \[
    B = \frac{1}{R_2^{2-n} - R_1^{2-n}}, \quad A = -\frac{R_1^{2-n}}{R_2^{2-n} - R_1^{2-n}}
    \]

    Substituting these back into the original function and simplifying, we obtain the final hitting probability:

    \[
    u(x) = \frac{|x|^{2-n} - R_1^{2-n}}{R_2^{2-n} - R_1^{2-n}}
    \]

---
<br>

## Part II: European Option Pricing Practice

### Exercise 3
**Problem**
Derive the Black-Scholes formula for a European Put Option.

??? success "Solution (Click to expand)"

    **1. Establish the Black-Scholes Partial Differential Equation (PDE)**
    Let the European put option price be $P(S, t)$. Based on the no-arbitrage principle and the geometric Brownian motion assumption, the put option also satisfies the BS PDE:

    \[
    \frac{\partial P}{\partial t} + rS \frac{\partial P}{\partial S} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 P}{\partial S^2} - rP = 0
    \]

    The terminal condition is: $P(S, T) = \max(K - S, 0) = (K - S)^+$.

    **2. Variable Substitution and Equation Simplification**
    Let $x = \ln S$ (i.e., $S = e^x$), and reverse the time $\tau = T - t$. Let $V(x, \tau) = P(e^x, T - \tau)$. After substituting into the equation, to eliminate the first-order partial derivative and the constant term, we make another substitution $V(x, \tau) = u(x, \tau) e^{\alpha \tau + \beta x}$, solving for the coefficients:

    \[
    \beta = \frac{1}{2} - \frac{r}{\sigma^2}, \quad \alpha = -r - \frac{1}{2\sigma^2}\left(r - \frac{1}{2}\sigma^2\right)^2
    \]

    At this point, the equation transforms into the standard heat equation:

    \[
    \frac{\partial u}{\partial \tau} = \frac{1}{2}\sigma^2 \frac{\partial^2 u}{\partial x^2}
    \]

    **3. Integration using the Gaussian Kernel**
    The transformed initial condition is $u(y, 0) = e^{-\beta y} \max(K - e^y, 0)$. The solution to the heat equation can be given by the convolution with the Gaussian kernel:

    \[
    u(x, \tau) = \frac{1}{\sigma \sqrt{2\pi \tau}} \int_{-\infty}^{\ln K} \exp\left(-\frac{(x-y)^2}{2\sigma^2 \tau}\right) e^{-\beta y} (K - e^y) dy
    \]

    (Note that the upper limit of integration is $\ln K$, because the payoff is 0 when $e^y > K$.)
    
    We split the integral into two terms, $I_1$ (containing $K$) and $I_2$ (containing $e^y$). Apply "completing the square" to the exponential parts and use the variable substitution $z = \frac{y - \mu}{\sigma}$ to transform them into standard normal cumulative distribution functions $N(\cdot)$.

    **4. Revert to the Final Pricing Formula**
    Finally, use $P(S, t) = e^{\alpha \tau + \beta x} u(x, \tau)$ to revert the variables. The exponential terms cancel out perfectly, yielding the Black-Scholes formula for the put option:

    \[
    P(S, t) = K e^{-r(T-t)} N(-d_2) - S N(-d_1)
    \]

    Where:

    \[
    d_1 = \frac{\ln(S/K) + (r + \frac{1}{2}\sigma^2)(T-t)}{\sigma \sqrt{T-t}}
    \]

    \[
    d_2 = d_1 - \sigma \sqrt{T-t}
    \]

    The derivation is complete.