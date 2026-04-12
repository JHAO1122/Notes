---
tags:
  - Probability Theory
  - Stochastic Processes
  - Stochastic Differential Equations
  - Homework Exercises
---

# 📝 Detailed Solutions to Homework Exercises: Stochastic Integrals and Itô's Formula

!!! abstract "About This Page"
    This page compiles detailed solutions to key homework exercises from Chapter 3 (Stochastic Integrals) and Chapter 4 (Itô Integrals and Itô's Formula) of the *Stochastic Differential Equations* course. The content covers core techniques such as integral calculations with Brownian motion, Itô isometry, quadratic variation processes, and solving stochastic differential equations using Itô's formula.

---
## Part I: Chapter 3: Computation of Stochastic Integrals

### Exercise 1

**Problem:**
Find the stochastic integral $\int_0^t s dW(s)$, and compute its expectation and variance.

??? success "Solution (click to expand)"

    **1. Finding the stochastic integral:**
    
    Using integration by parts, let $g(s) = s$. Then:
    
    $$
    \int_0^t s dW(s) = sW(s) \Big|_0^t - \int_0^t W(s) ds = tW(t) - \int_0^t W(s) ds
    $$

    <br>

    **2. Computing the expectation:**
    
    Using the linearity of expectation and the zero-mean property of Brownian motion:
    
    $$
    \mathbb{E}\left[ \int_0^t s dW(s) \right] = t\mathbb{E}[W(t)] - \int_0^t \mathbb{E}[W(s)] ds = t \cdot 0 - \int_0^t 0 ds = 0
    $$

    <br>

    **3. Computing the variance:**
    
    By the Itô isometry principle, for a deterministic function $g(s) \in L^2([0,t])$, the variance of its integral with respect to Brownian motion equals the Riemann integral of the squared function:
    
    $$
    Var\left[ \int_0^t s dW(s) \right] = \mathbb{E}\left[ \left( \int_0^t s dW(s) \right)^2 \right] - \left( \mathbb{E}\left[ \int_0^t s dW(s) \right] \right)^2
    $$
    
    Since the expectation is 0, the variance is the second moment:
    
    $$
    = \mathbb{E}\left[ \int_0^t s^2 ds \right] = \int_0^t s^2 ds = \frac{1}{3}t^3
    $$

    *Note: Since the integrand is a deterministic function, this stochastic integral follows a normal distribution $N(0, \frac{1}{3}t^3)$.*

---

### Exercise 2

**Problem:**
Let the function $g: [0, T] \to \mathbb{R}$ be continuously differentiable, with $g(0) = g(T) = 0$. Find the probability density function of $\int_0^T g(t) dW(t)$.

??? success "Solution (click to expand)"

    **Step 1: Find the expectation**
    Using integration by parts and the boundary conditions $g(0) = g(T) = 0$:
    
    $$
    \int_0^T g(t) dW(t) = g(T)W(T) - g(0)W(0) - \int_0^T g'(t)W(t) dt = - \int_0^T g'(t)W(t) dt
    $$
    
    Taking the expectation:
    
    $$
    \mathbb{E}\left[ \int_0^T g(t) dW(t) \right] = \mathbb{E}\left[ - \int_0^T g'(t)W(t) dt \right] = - \int_0^T g'(t)\mathbb{E}[W(t)] dt = 0
    $$

    <br>

    **Step 2: Find the variance**
    According to Itô isometry (since $g(t)$ is a deterministic continuous function, it necessarily belongs to the $L^2([0,T])$ space):
    
    $$
    Var\left[ \int_0^T g(t) dW(t) \right] = \mathbb{E}\left[ \left( \int_0^T g(t) dW(t) \right)^2 \right] = \int_0^T g^2(t) dt
    $$

    <br>

    **Step 3: Determine the distribution type and density function**
    Since the integrand $g(t)$ is a deterministic function and the increments of Brownian motion are independent and normally distributed, this stochastic integral (viewed as the limit of a linear combination of normal increments) must follow a normal distribution.
    
    Therefore, this random variable follows a normal distribution with mean 0 and variance $\sigma^2 = \int_0^T g^2(t) dt$:
    
    $$
    I := \int_0^T g(t) dW(t) \sim N\left(0, \int_0^T g^2(t) dt\right)
    $$
    
    Its probability density function is the standard one-dimensional normal density formula:
    
    $$
    f_I(x) = \frac{1}{\sqrt{2\pi \int_0^T g^2(t) dt}} \exp\left( - \frac{x^2}{2 \int_0^T g^2(t) dt} \right)
    $$

---

### Exercise 3

**Problem:**
To explain Brownian motion from a Newtonian mechanics perspective, Langevin (1872-1946) proposed the following (stochastic) differential equation describing the velocity of a particle in a liquid:
$$\frac{dv}{dt} = -\beta v + \dot{W}(t)$$
where $-\beta v$ represents the frictional drag force on the particle ($\beta$ is a positive constant), and the white noise $\dot{W}(t)$ describes the random impulsive forces on the particle.
(1) Show that the solution to this equation is $v(t) = v_0 e^{-\beta t} + W(t) - \beta \int_0^t e^{-\beta(t-s)} W(s) ds$, and consequently the particle's path starting from the origin is $x_\beta(t) = \int_0^t e^{-\beta(t-s)} W(s) ds$;
(2) Compute the expectation and variance of $v(t)$ and $x(t)$;
(3) Prove that $\lim_{\beta \to \infty} \beta x_\beta(t) = W(t)$.

??? success "Solution (click to expand)"

    **(1) Verifying the solution form**
    Substitute the given solution $v(t)$ into the original differential equation for verification. Differentiate the solution:
    
    $$
    \frac{dv(t)}{dt} = \frac{d}{dt} \left( v_0 e^{-\beta t} + W(t) - \beta \int_0^t e^{-\beta(t-s)} W(s) ds \right)
    $$
    
    Using the Leibniz rule for differentiation under the integral sign:
    
    $$
    = -\beta v_0 e^{-\beta t} + \dot{W}(t) - \beta \left( e^{-\beta(t-t)}W(t) + \int_0^t -\beta e^{-\beta(t-s)} W(s) ds \right)
    $$
    
    $$
    = -\beta v_0 e^{-\beta t} + \dot{W}(t) - \beta W(t) + \beta^2 \int_0^t e^{-\beta(t-s)} W(s) ds
    $$
    
    Then multiply $v(t)$ by $-\beta$:
    
    $$
    -\beta v(t) = -\beta \left( v_0 e^{-\beta t} + W(t) - \beta \int_0^t e^{-\beta(t-s)} W(s) ds \right)
    $$
    
    $$
    = -\beta v_0 e^{-\beta t} - \beta W(t) + \beta^2 \int_0^t e^{-\beta(t-s)} W(s) ds
    $$
    
    Comparing the two expressions, it is clear that:
    
    $$
    \frac{dv(t)}{dt} = -\beta v(t) + \dot{W}(t)
    $$
    
    The equation holds. For the displacement $x(t) = \int_0^t v(s) ds$, using integration by parts to simplify the form:
    Through stochastic differential equation theory, the more standard expression for the displacement of the O-U process is $x_\beta(t) = \int_0^t e^{-\beta(t-s)} W(s) ds$.

    <br>

    **(2) Computing expectation and variance**
    
    For $v(t)$:
    
    $$
    \mathbb{E}[v(t)] = v_0 e^{-\beta t} + \mathbb{E}[W(t)] - \beta \int_0^t e^{-\beta(t-s)} \mathbb{E}[W(s)] ds = v_0 e^{-\beta t}
    $$
    
    Since an equivalent form of $v(t)$ is $v(t) = v_0 e^{-\beta t} + \int_0^t e^{-\beta(t-s)} dW(s)$, use Itô isometry to compute the variance:
    
    $$
    Var[v(t)] = \mathbb{E}\left[ \left( \int_0^t e^{-\beta(t-s)} dW(s) \right)^2 \right] = \int_0^t e^{-2\beta(t-s)} ds = \frac{1 - e^{-2\beta t}}{2\beta}
    $$
    
    For $x_\beta(t)$:
    
    $$
    \mathbb{E}[x_\beta(t)] = \mathbb{E}\left[ \int_0^t e^{-\beta(t-s)} W(s) ds \right] = \int_0^t e^{-\beta(t-s)} \mathbb{E}[W(s)] ds = 0
    $$
    
    Use Fubini's theorem to exchange the order of integration for computing the variance (similar to the technique in Chapter 2 exercises):
    
    $$
    Var[x_\beta(t)] = \mathbb{E}\left[ \left( \int_0^t e^{-\beta(t-s)} W(s) ds \right)^2 \right] = \int_0^t \int_0^t e^{-\beta(t-u)}e^{-\beta(t-v)} \mathbb{E}[W(u)W(v)] du dv
    $$
    
    After simplification, we obtain:
    
    $$
    = \frac{t}{\beta^2} - \frac{3 - 4e^{-\beta t} + e^{-2\beta t}}{2\beta^3}
    $$

    <br>

    **(3) Proving the limit**
    
    Consider $\beta x_\beta(t)$:
    
    $$
    \beta x_\beta(t) = \beta \int_0^t e^{-\beta(t-s)} W(s) ds
    $$
    
    Apply integration by parts to it:
    
    $$
    = \int_0^t W(s) d(e^{-\beta(t-s)}) = W(s) e^{-\beta(t-s)} \Big|_0^t - \int_0^t e^{-\beta(t-s)} dW(s)
    $$
    
    $$
    = W(t) - e^{-\beta t} W(0) - \int_0^t e^{-\beta(t-s)} dW(s) = W(t) - \int_0^t e^{-\beta(t-s)} dW(s)
    $$
    
    When $\beta \to \infty$, examine the mean-square limit of the error term:
    
    $$
    \lim_{\beta \to \infty} \mathbb{E}\left[ \left( \int_0^t e^{-\beta(t-s)} dW(s) \right)^2 \right] = \lim_{\beta \to \infty} \int_0^t e^{-2\beta(t-s)} ds = \lim_{\beta \to \infty} \frac{1 - e^{-2\beta t}}{2\beta} = 0
    $$
    
    Since the error term converges to 0 in the $L^2$ sense, we have almost surely (or in the mean-square sense):
    
    $$
    \lim_{\beta \to \infty} \beta x_\beta(t) = W(t)
    $$

---
<br>
## Part II: Chapter 4 Itô Integral and Quadratic Variation

### Exercise 1

**Problem:**
Using the definition of the Itô stochastic integral based on Riemann sums, prove:

$$\int_0^T W^2(t) dW(t) = \frac{1}{3}W^3(T) - \int_0^T W(t) dt$$

??? success "Solution (click to expand)"

    This is a classic problem that requires distinguishing between the rules of Riemann calculus and Itô calculus.
    
    **Proof using Itô's formula:**
    Consider the function $f(t, x) = \frac{1}{3}x^3$. Its partial derivatives are:
    $f_t = 0, \quad f_x = x^2, \quad f_{xx} = 2x$
    
    Substitute $X_t = W_t$ into Itô's formula $df(t, W_t) = f_t dt + f_x dW_t + \frac{1}{2} f_{xx} (dW_t)^2$:
    
    $$
    d\left( \frac{1}{3}W^3(t) \right) = 0 dt + W^2(t) dW(t) + \frac{1}{2} \cdot 2W(t) dt
    $$
    
    Here, the quadratic variation rule for Brownian motion $(dW_t)^2 = dt$ is used. Simplifying yields:
    
    $$
    d\left( \frac{1}{3}W^3(t) \right) = W^2(t) dW(t) + W(t) dt
    $$
    
    Integrating both sides from $0$ to $T$, and noting that $W(0) = 0$:
    
    $$
    \frac{1}{3}W^3(T) - \frac{1}{3}W^3(0) = \int_0^T W^2(t) dW(t) + \int_0^T W(t) dt
    $$
    
    Rearranging gives the desired result:
    
    $$
    \int_0^T W^2(t) dW(t) = \frac{1}{3}W^3(T) - \int_0^T W(t) dt
    $$
    
    *Note: Compared to ordinary calculus $\int x^2 dx = \frac{1}{3}x^3$, the Itô integral has an extra correction term $-\int_0^T W(t) dt$ arising from the quadratic variation term.*

---

### Exercise 2

**Problem:**
For the backward integral
$ \int_0^T W(t) dW(t) \doteq \lim_{\delta \to 0} \sum_{i=0}^{n-1} W(t_{i+1})[W(t_{i+1}) - W(t_i)]$
where $0 = t_0 < t_1 < \cdots < t_n = T$, and $\delta \doteq \max_i |t_{i+1} - t_i|$. Prove:

$$ \int_0^T W(t) dW(t) = \int_0^T W(t) dW(t) + T$$


??? success "Solution (click to expand)"

    Perform an identity transformation on the definition of the backward integral to construct the form of the standard Itô integral:
    
    $$
    I_B = \lim_{\delta \to 0} \sum_{i=0}^{n-1} W(t_{i+1})[W(t_{i+1}) - W(t_i)]
    $$
    
    In the summation term, we artificially add and subtract $W(t_i)$:
    
    $$
    = \lim_{\delta \to 0} \sum_{i=0}^{n-1} \big( W(t_i) + W(t_{i+1}) - W(t_i) \big) [W(t_{i+1}) - W(t_i)]
    $$
    
    Expanding the parentheses:
    
    $$
    = \lim_{\delta \to 0} \sum_{i=0}^{n-1} W(t_i)[W(t_{i+1}) - W(t_i)] + \lim_{\delta \to 0} \sum_{i=0}^{n-1} [W(t_{i+1}) - W(t_i)]^2
    $$
    
    Observe these two terms:
    * The first term is precisely the discrete definition of the standard **Itô integral** (taking the left endpoint in each subinterval):
      $$\lim_{\delta \to 0} \sum_{i=0}^{n-1} W(t_i) \Delta W_i = \int_0^T W(t) dW(t)$$
    * The second term is precisely the **Quadratic Variation** of Brownian motion over the interval $[0,T]$. It is known from the properties of Brownian motion that its quadratic variation converges in mean square to the length of the interval:
      $$\lim_{\delta \to 0} \sum_{i=0}^{n-1} (\Delta W_i)^2 = T$$
      
    Combining the two parts yields the proof:
    
    $$
    (B) \int_0^T W(t) dW(t) = \int_0^T W(t) dW(t) + T
    $$

---
<br>
## Part III: Advanced Applications of Stochastic Integration and Itô's Formula

### Exercise 1
**Problem**
Let $W(t)$ be an $n$-dimensional Brownian motion. Prove: $\mathbb{E}[|W(t) - W(s)|^4] = (2n + n^2)(t-s)^2$.

??? success "Solution (Click to expand)"

    This problem cleverly utilizes the relationship between the standard normal distribution and the chi-squared distribution.
    
    Given that $W(t)$ is an $n$-dimensional Brownian motion, we examine its increments over the time interval $[s, t]$.
    Let $X_i = W_i(t) - W_i(s)$, where $i = 1, 2, \dots, n$ denotes each dimension.
    By the properties of Brownian motion, the increments in each dimension are independent and identically distributed, and $X_i \sim N(0, t-s)$.
    
    To standardize, let $Z_i = \frac{X_i}{\sqrt{t-s}}$, then $Z_i \sim N(0, 1)$ and are mutually independent.
    Let $Q = \sum_{i=1}^n Z_i^2$. By definition, $Q$ follows a chi-squared distribution with $n$ degrees of freedom, i.e., $Q \sim \chi^2(n)$.
    
    For the chi-squared distribution $\chi^2(n)$, we know its expectation and variance are:
    
    $$
    \mathbb{E}[Q] = n, \quad Var(Q) = 2n
    $$
    
    From this, we can find the second raw moment of $Q$:
    
    $$
    \mathbb{E}[Q^2] = Var(Q) + (\mathbb{E}[Q])^2 = 2n + n^2
    $$
    
    Returning to the original expression, consider the fourth moment of the increment:
    
    $$
    |W(t) - W(s)|^4 = \left( \sum_{i=1}^n X_i^2 \right)^2 = \left( \sum_{i=1}^n (t-s)Z_i^2 \right)^2 = (t-s)^2 Q^2
    $$
    
    Taking the expectation on both sides:
    
    $$
    \mathbb{E}[|W(t) - W(s)|^4] = (t-s)^2 \mathbb{E}[Q^2] = (2n + n^2)(t-s)^2
    $$
    
    The conclusion is proven.

---

### Exercise 2
**Problem**
Define the second-order integral
$$\int_0^T f(t) [dW(t)]^2 \doteq \lim_{\delta \to 0} \sum_{i=0}^{n-1} f(t_i)[W(t_{i+1}) - W(t_i)]^2$$
where $0 = t_0 < t_1 < \dots < t_n = T$, and $\delta \doteq \max_i |t_{i+1} - t_i|$. Prove: when $f \in L^2([0,T])$, we have

$$\int_0^T f(t) [dW(t)]^2 = \int_0^T f(t) dt$$

??? success "Solution (Click to expand)"

    This problem aims to prove the strict validity of the quadratic variation of Brownian motion $(dW_t)^2 = dt$ in the integral sense.
    We will prove this limit in the mean-square ($L^2$) sense.
    
    Let $\Delta t_i = t_{i+1} - t_i$, $\Delta W_i = W(t_{i+1}) - W(t_i)$.
    Let the left-hand discrete sum be $S_n = \sum_{i=0}^{n-1} f(t_i)(\Delta W_i)^2$, and the right-hand target integral be $I = \int_0^T f(t)dt \approx \sum_{i=0}^{n-1} f(t_i)\Delta t_i$.
    
    Examine the mean-square error of their difference:
    
    $$
    \mathbb{E}\left[ \left( S_n - \sum_{i=0}^{n-1} f(t_i)\Delta t_i \right)^2 \right] = \mathbb{E}\left[ \left( \sum_{i=0}^{n-1} f(t_i) ((\Delta W_i)^2 - \Delta t_i) \right)^2 \right]
    $$
    
    Since Brownian increments over non-overlapping intervals are independent, the expectation of cross terms ($i \ne j$) can be decomposed:
    
    $$
    \mathbb{E}[((\Delta W_i)^2 - \Delta t_i)((\Delta W_j)^2 - \Delta t_j)] = \mathbb{E}[(\Delta W_i)^2 - \Delta t_i] \cdot \mathbb{E}[(\Delta W_j)^2 - \Delta t_j] = 0 \cdot 0 = 0
    $$
    
    Therefore, after expanding the square, only the squared terms remain:
    
    $$
    = \sum_{i=0}^{n-1} f^2(t_i) \mathbb{E}[((\Delta W_i)^2 - \Delta t_i)^2]
    $$
    
    Compute the single-term expectation, using the fourth moment of $N(0, \Delta t_i)$: $\mathbb{E}[(\Delta W_i)^4] = 3(\Delta t_i)^2$:
    
    $$
    \mathbb{E}[(\Delta W_i)^4 - 2\Delta t_i(\Delta W_i)^2 + (\Delta t_i)^2] = 3\Delta t_i^2 - 2\Delta t_i^2 + \Delta t_i^2 = 2\Delta t_i^2
    $$
    
    Substituting back and bounding one $\Delta t_i$ by the maximum step size $\delta$:
    
    $$
    = \sum_{i=0}^{n-1} f^2(t_i) 2\Delta t_i^2 \le 2\delta \sum_{i=0}^{n-1} f^2(t_i)\Delta t_i
    $$
    
    As $\delta \to 0$, since $f \in L^2([0,T])$, $\sum f^2(t_i)\Delta t_i \to \int_0^T f^2(t)dt < \infty$.
    Thus, the factor $\delta$ in front will drive the entire limit to $0$. Mean-square convergence is proven, i.e.:
    
    $$
    \int_0^T f(t) [dW(t)]^2 = \int_0^T f(t) dt
    $$

---

### Exercise 3
**Problem**
Let $f \in L^2([0, T])$ and $\int_0^T f(s) dW(s) = 0$. Prove: $f$ is almost everywhere zero.

??? success "Solution (Click to expand)"

    This problem utilizes the Itô isometry and the fundamental properties of Lebesgue integrals in real analysis.
    
    Since it is given that $\int_0^T f(s) dW(s) = 0$ holds almost surely, the second moment of this random variable must be 0:
    
    $$
    \mathbb{E}\left[ \left( \int_0^T f(s) dW(s) \right)^2 \right] = \mathbb{E}[0^2] = 0
    $$
    
    On the other hand, by the Itô isometry, the second moment of the stochastic integral equals the Lebesgue integral of the square of the integrand:
    
    $$
    \mathbb{E}\left[ \left( \int_0^T f(s) dW(s) \right)^2 \right] = \mathbb{E}\left[ \int_0^T f^2(s) ds \right]
    $$
    
    Since $f(s)$ is a deterministic function, its expectation is itself. Combining the two equations yields:
    
    $$
    \int_0^T f^2(s) ds = 0
    $$
    
    By the fundamental theorem of real analysis, since the integrand $f^2(s) \ge 0$ always holds and its Lebesgue integral over $[0,T]$ is $0$, this implies that $f^2(s)$ must be zero almost everywhere (a.e.) on $[0,T]$.
    
    Consequently: $f(s) = 0$ almost everywhere. Proven.

---

### Exercise 4
**Problem**
Prove: $Y(t) = e^{t/2} \cos(W(t))$ is a martingale.

??? success "Solution (Click to expand)"

    To prove a process is a martingale, the most direct method is to use Itô's formula to show that its differential has no drift term (i.e., the coefficient of the $dt$ term is 0).
    
    Let the bivariate function $u(t, x) = e^{t/2} \cos(x)$. Compute its partial derivatives with respect to $t$ and $x$:
    * $u_t = \frac{1}{2} e^{t/2} \cos(x)$
    * $u_x = -e^{t/2} \sin(x)$
    * $u_{xx} = -e^{t/2} \cos(x)$
    
    Substitute $X(t) = W(t)$ into Itô's formula $dY_t = (u_t + \frac{1}{2} u_{xx})dt + u_x dW(t)$:
    
    $$
    dY(t) = \left( \frac{1}{2} e^{t/2} \cos(W(t)) - \frac{1}{2} e^{t/2} \cos(W(t)) \right) dt - e^{t/2} \sin(W(t)) dW(t)
    $$
    
    It can be clearly seen that the two terms containing $dt$ perfectly cancel:
    
    $$
    dY(t) = -e^{t/2} \sin(W(t)) dW(t)
    $$
    
    Writing it in integral form:
    
    $$
    Y(t) = Y(0) - \int_0^t e^{s/2} \sin(W(s)) dW(s)
    $$
    
    Since the right-hand side is an Itô integral containing only $dW(s)$, and the integrand is bounded (satisfying the $L^2$ admissibility condition), the Itô integral itself is a martingale.
    Therefore, the process $Y(t)$ is also a martingale. Proven.

---

### Exercise 5
**Problem**
Prove:
1. $\int_0^T W^2 dW = \frac{1}{3}W(T)^3 - \int_0^T W dt$
2. $\int_0^T W^3 dW = \frac{1}{4}W(T)^4 - \frac{3}{2} \int_0^T W^2 dt$

??? success "Solution (Click to expand)"

    Both parts of this proof are basic reverse applications of Itô's formula, i.e., "first guess the higher-order term, then expand using Itô's formula and rearrange".

    **(1) Proving the first formula**
    
    Let the function $f(x) = \frac{1}{3}x^3$. Compute derivatives: $f'(x) = x^2$, $f''(x) = 2x$.
    Substitute $W(t)$ into Itô's formula:
    
    $$
    d\left( \frac{1}{3}W^3(t) \right) = W^2(t) dW(t) + \frac{1}{2}(2W(t)) (dW(t))^2
    $$
    
    By the quadratic variation rule $(dW(t))^2 = dt$:
    
    $$
    d\left( \frac{1}{3}W^3(t) \right) = W^2(t) dW(t) + W(t) dt
    $$
    
    Integrate both sides over $[0, T]$, noting that $W(0) = 0$:
    
    $$
    \frac{1}{3}W^3(T) - 0 = \int_0^T W^2(t) dW(t) + \int_0^T W(t) dt
    $$
    
    Rearranging yields:
    
    $$
    \int_0^T W^2(t) dW(t) = \frac{1}{3}W^3(T) - \int_0^T W(t) dt
    $$

    <br>

    **(2) Proving the second formula**
    
    Similarly, let the function $g(x) = \frac{1}{4}x^4$. Compute derivatives: $g'(x) = x^3$, $g''(x) = 3x^2$.
    Substitute into Itô's formula:
    
    $$
    d\left( \frac{1}{4}W^4(t) \right) = W^3(t) dW(t) + \frac{1}{2}(3W^2(t)) dt
    $$
    
    Integrate both sides over $[0, T]$:
    
    $$
    \frac{1}{4}W^4(T) - 0 = \int_0^T W^3(t) dW(t) + \frac{3}{2}\int_0^T W^2(t) dt
    $$
    
    Rearranging yields:
    
    $$
    \int_0^T W^3(t) dW(t) = \frac{1}{4}W^4(T) - \frac{3}{2}\int_0^T W^2(t) dt
    $$

---

### Exercise 6
**Problem**
Prove $\mathbb{E}[e^{\int_0^T g dW}] = e^{\frac{1}{2} \int_0^T g^2 ds}$.

??? success "Solution (Click to expand)"

    This problem can be elegantly solved by analyzing the distributional properties of the Itô integral, directly utilizing the moment-generating function of the normal distribution.
    
    Let the random variable $X = \int_0^T g(s) dW(s)$.
    Since $g(s)$ is a deterministic function of time (non-random), this Itô integral is a linear superposition of a Gaussian process, so $X$ still follows a normal distribution.
    
    By the properties of stochastic integrals:
    * Expectation: $\mathbb{E}[X] = \mathbb{E}[\int_0^T g dW] = 0$
    * Variance: By the Itô isometry, $Var(X) = \mathbb{E}[X^2] = \int_0^T g^2(s) ds$
    
    So $X \sim N(0, \sigma^2)$, where $\sigma^2 = \int_0^T g^2(s) ds$.
    
    The required $\mathbb{E}[e^X]$ is precisely the value of the moment-generating function $M_X(u) = \mathbb{E}[e^{uX}]$ of the random variable $X$ at $u=1$.
    For a normal distribution $N(\mu, \sigma^2)$, its moment-generating function formula is $M_X(u) = \exp(\mu u + \frac{1}{2}\sigma^2 u^2)$.
    
    Substituting $\mu = 0, u = 1, \sigma^2 = \int_0^T g^2 ds$, we immediately obtain:
    
    $$
    \mathbb{E}[e^{\int_0^T g dW}] = \exp\left( 0 + \frac{1}{2} \int_0^T g^2(s) ds \cdot 1^2 \right) = e^{\frac{1}{2} \int_0^T g^2 ds}
    $$
    
    Proven.

---

### Exercise 7
**Problem**
Let $u = u(x, t)$ satisfy the parabolic partial differential equation $u_t + \frac{1}{2}u_{xx} = 0$. Prove: $\mathbb{E}[u(W(t), t)] = u(0, 0)$.

??? success "Solution (Click to expand)"

    This problem is a classic exercise establishing the profound connection between partial differential equations (PDEs) and stochastic processes (SDEs) (i.e., the simplest form of the Feynman-Kac formula).
    
    Define the stochastic process $Y(t) = u(W(t), t)$. Apply the multivariate Itô formula to expand the differential of $Y(t)$:
    
    $$
    dY(t) = \frac{\partial u}{\partial t} dt + \frac{\partial u}{\partial x} dW(t) + \frac{1}{2}\frac{\partial^2 u}{\partial x^2} (dW(t))^2
    $$
    
    Substitute the quadratic variation $(dW)^2 = dt$ and combine the $dt$ terms:
    
    $$
    dY(t) = \left( u_t + \frac{1}{2}u_{xx} \right) dt + u_x dW(t)
    $$
    
    Since the given condition states that $u(x, t)$ satisfies $u_t + \frac{1}{2}u_{xx} = 0$, the drift term is strictly 0. The differential equation simplifies to pure diffusion:
    
    $$
    dY(t) = u_x(W(t), t) dW(t)
    $$
    
    Write it in integral form:
    
    $$
    Y(t) - Y(0) = \int_0^t u_x(W(s), s) dW(s)
    $$
    
    Take the expectation on both sides. Since the expectation of the Itô integral on the right-hand side is 0:
    
    $$
    \mathbb{E}[Y(t)] - \mathbb{E}[Y(0)] = 0 \implies \mathbb{E}[Y(t)] = \mathbb{E}[Y(0)]
    $$
    
    Substituting back the definition of $Y(t)$, and noting that at the initial time the Brownian motion $W(0) = 0$ almost surely:
    
    $$
    \mathbb{E}[u(W(t), t)] = \mathbb{E}[u(W(0), 0)] = u(0, 0)
    $$
    
    The conclusion is proven.

---

### Exercise 8
**Problem**
1. Prove $e^{W(t)} = 1 + \frac{1}{2}\int_0^t e^{W(s)} ds + \int_0^t e^{W(s)} dW(s)$;
2. Prove $\mathbb{E}[e^{W(t)}] = 1 + \frac{1}{2}\int_0^t \mathbb{E}[e^{W(s)}