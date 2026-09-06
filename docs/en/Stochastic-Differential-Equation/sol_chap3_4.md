---
tags:
  - Stochastic Differential Equations
  - Homework Exercises
---

# 📝 Detailed Solutions to Homework Exercises: Stochastic Integrals and Itô's Formula

!!! abstract "About This Page"
    This page compiles detailed solutions to key homework exercises from Chapter 3 (Stochastic Integrals) and Chapter 4 (Itô Integral and Itô's Formula) of the *Stochastic Differential Equations* course. The content covers core techniques including integral calculations with respect to Brownian motion, Itô isometry, quadratic variation processes, and solving stochastic differential equations using Itô's formula.

---

## Part I: Chapter 3 - Calculation of Stochastic Integrals

### Exercise 1

**Problem:**
Find the stochastic integral $\int_0^t s dW(s)$, and compute its expectation and variance.

??? success "Solution (Click to expand)"

    **1. Finding the stochastic integral:**
    
    Using the integration by parts formula, let $g(s) = s$, then:
    
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
    
    By the Itô Isometry principle, for a deterministic function $g(s) \in L^2([0,t])$, the variance of its integral with respect to Brownian motion equals the Riemann integral of the square of the function:
    
    $$
    Var\left[ \int_0^t s dW(s) \right] = \mathbb{E}\left[ \left( \int_0^t s dW(s) \right)^2 \right] - \left( \mathbb{E}\left[ \int_0^t s dW(s) \right] \right)^2
    $$
    
    Since the expectation is 0, the variance is the second moment:
    
    $$
    = \mathbb{E}\left[ \int_0^t s^2 ds \right] = \int_0^t s^2 ds = \frac{1}{3}t^3
    $$

    *Note: Because the integrand is a deterministic function, this stochastic integral follows a normal distribution $N(0, \frac{1}{3}t^3)$.*

---

### Exercise 2

**Problem:**
Let the function $g: [0, T] \to \mathbb{R}$ be continuously differentiable, with $g(0) = g(T) = 0$. Find the probability density function of $\int_0^T g(t) dW(t)$.

??? success "Solution (Click to expand)"

    **Step 1: Find the Expectation**
    Using integration by parts and the boundary conditions $g(0) = g(T) = 0$:
    
    $$
    \int_0^T g(t) dW(t) = g(T)W(T) - g(0)W(0) - \int_0^T g'(t)W(t) dt = - \int_0^T g'(t)W(t) dt
    $$
    
    Taking the expectation:
    
    $$
    \mathbb{E}\left[ \int_0^T g(t) dW(t) \right] = \mathbb{E}\left[ - \int_0^T g'(t)W(t) dt \right] = - \int_0^T g'(t)\mathbb{E}[W(t)] dt = 0
    $$

    <br>

    **Step 2: Find the Variance**
    According to the Itô isometry (since $g(t)$ is a deterministic continuous function, it must belong to the $L^2([0,T])$ space):
    
    $$
    Var\left[ \int_0^T g(t) dW(t) \right] = \mathbb{E}\left[ \left( \int_0^T g(t) dW(t) \right)^2 \right] = \int_0^T g^2(t) dt
    $$

    <br>

    **Step 3: Determine the Distribution Type and Density Function**
    Since the integrand $g(t)$ is a deterministic function, the increments of Brownian motion are independent and normally distributed. As the limit of linear combinations of normal increments (viewed as an infinite-dimensional linear combination), this stochastic integral must follow a normal distribution.
    
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
To explain Brownian motion from the perspective of Newtonian mechanics, Langevin (1872-1946) proposed the following (stochastic) differential equation describing the velocity of a particle in a liquid:

$$\frac{dv}{dt} = -\beta v + \dot{W}(t)$$

where $-\beta v$ represents the frictional resistance experienced by the particle's motion ($\beta$ is a positive constant), and the white noise $\dot{W}(t)$ describes the random impulsive forces on the particle.
(1) Show that the solution to this equation is $v(t) = v_0 e^{-\beta t} + W(t) - \beta \int_0^t e^{-\beta(t-s)} W(s) ds$, and thus the path of a particle starting from the origin is $x_\beta(t) = \int_0^t e^{-\beta(t-s)} W(s) ds$;
(2) Calculate the expectation and variance of $v(t), x(t)$;
(3) Prove that $\lim_{\beta \to \infty} \beta x_\beta(t) = W(t)$.

??? success "Solution (click to expand)"

    **(1) Verification of the Solution Form**
    Substitute the given solution $v(t)$ into the original differential equation for verification. Differentiate the solution:
    
    $$
    \frac{dv(t)}{dt} = \frac{d}{dt} \left( v_0 e^{-\beta t} + W(t) - \beta \int_0^t e^{-\beta(t-s)} W(s) ds \right)
    $$
    
    Using the differentiation rule for integrals with variable upper limits (Leibniz Rule):
    
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
    
    Comparing the two expressions, it is evident that:
    
    $$
    \frac{dv(t)}{dt} = -\beta v(t) + \dot{W}(t)
    $$
    
    The equation holds. For the displacement $x(t) = \int_0^t v(s) ds$, simplify the form using integration by parts:
    According to stochastic differential equation theory, the more standard expression for the displacement of the O-U process is $x_\beta(t) = \int_0^t e^{-\beta(t-s)} W(s) ds$.

    <br>

    **(2) Computing Expectation and Variance**

    For $v(t)$:

    $$
    \mathbb{E}[v(t)] = v_0 e^{-\beta t} + \mathbb{E}[W(t)] - \beta \int_0^t e^{-\beta(t-s)} \mathbb{E}[W(s)] ds = v_0 e^{-\beta t}
    $$

    Since another equivalent form of $v(t)$ is $v(t) = v_0 e^{-\beta t} + \int_0^t e^{-\beta(t-s)} dW(s)$, we compute the variance using the Itô isometry:

    $$
    Var[v(t)] = \mathbb{E}\left[ \left( \int_0^t e^{-\beta(t-s)} dW(s) \right)^2 \right] = \int_0^t e^{-2\beta(t-s)} ds = \frac{1 - e^{-2\beta t}}{2\beta}
    $$

    For $x_\beta(t)$:

    $$
    \mathbb{E}[x_\beta(t)] = \mathbb{E}\left[ \int_0^t e^{-\beta(t-s)} W(s) ds \right] = \int_0^t e^{-\beta(t-s)} \mathbb{E}[W(s)] ds = 0
    $$

    We compute the variance by exchanging the order of integration using Fubini's theorem (similar to the technique in Chapter 2 exercises):

    $$
    Var[x_\beta(t)] = \mathbb{E}\left[ \left( \int_0^t e^{-\beta(t-s)} W(s) ds \right)^2 \right] = \int_0^t \int_0^t e^{-\beta(t-u)}e^{-\beta(t-v)} \mathbb{E}[W(u)W(v)] du dv
    $$

    After simplification, we obtain:

    $$
    = \frac{t}{\beta^2} - \frac{3 - 4e^{-\beta t} + e^{-2\beta t}}{2\beta^3}
    $$

    <br>

    **(3) Proving the Limit**

    Consider $\beta x_\beta(t)$:

    $$
    \beta x_\beta(t) = \beta \int_0^t e^{-\beta(t-s)} W(s) ds
    $$

    Applying integration by parts:

    $$
    = \int_0^t W(s) d(e^{-\beta(t-s)}) = W(s) e^{-\beta(t-s)} \Big|_0^t - \int_0^t e^{-\beta(t-s)} dW(s)
    $$

    $$
    = W(t) - e^{-\beta t} W(0) - \int_0^t e^{-\beta(t-s)} dW(s) = W(t) - \int_0^t e^{-\beta(t-s)} dW(s)
    $$

    As $\beta \to \infty$, examine the mean-square limit of the error term:

    $$
    \lim_{\beta \to \infty} \mathbb{E}\left[ \left( \int_0^t e^{-\beta(t-s)} dW(s) \right)^2 \right] = \lim_{\beta \to \infty} \int_0^t e^{-2\beta(t-s)} ds = \lim_{\beta \to \infty} \frac{1 - e^{-2\beta t}}{2\beta} = 0
    $$

    Since the error term converges to 0 in the $L^2$ sense, it follows almost surely (or in the mean-square sense) that:

    $$
    \lim_{\beta \to \infty} \beta x_\beta(t) = W(t)
    $$

---

## Part II: Chapter 4 Itô Integral and Quadratic Variation

### Exercise 1

**Problem:**
Using the definition of the Itô stochastic integral based on Riemann sums, prove:

$$\int_0^T W^2(t) dW(t) = \frac{1}{3}W^3(T) - \int_0^T W(t) dt$$

??? success "Solution (click to expand)"

    This is a classic problem that requires distinguishing between the rules of Riemann calculus and Itô calculus.

    **Using Itô's formula to prove:**
    We consider the function $f(t, x) = \frac{1}{3}x^3$. Its partial derivatives are:
    $f_t = 0, \quad f_x = x^2, \quad f_{xx} = 2x$

    Substituting $X_t = W_t$ into Itô's formula $df(t, W_t) = f_t dt + f_x dW_t + \frac{1}{2} f_{xx} (dW_t)^2$:

    $$
    d\left( \frac{1}{3}W^3(t) \right) = 0 dt + W^2(t) dW(t) + \frac{1}{2} \cdot 2W(t) dt
    $$

    where the quadratic variation rule for Brownian motion $(dW_t)^2 = dt$ is used. Simplifying yields:

    $$
    d\left( \frac{1}{3}W^3(t) \right) = W^2(t) dW(t) + W(t) dt
    $$

    Integrating both sides from $0$ to $T$, and noting that $W(0) = 0$:

    $$
    \frac{1}{3}W^3(T) - \frac{1}{3}W^3(0) = \int_0^T W^2(t) dW(t) + \int_0^T W(t) dt
    $$

    Rearranging gives the conclusion:

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

    By performing an identity transformation on the definition of the backward integral, we construct the form of the standard Itô integral:

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

    Observing these two terms:
    * The first term is precisely the discrete definition of the standard **Itô integral** (taking the left endpoint in each subinterval):
    
    $$\lim_{\delta \to 0} \sum_{i=0}^{n-1} W(t_i) \Delta W_i = \int_0^T W(t) dW(t)$$

    * The second term is precisely the **Quadratic Variation** of Brownian motion on the interval $[0,T]$. From the properties of Brownian motion, it is known that its quadratic variation converges in the mean square sense to the length of the interval:$\lim_{\delta \to 0} \sum_{i=0}^{n-1} (\Delta W_i)^2 = T$

    Combining the two parts, we obtain the proof:

    $$
    (B) \int_0^T W(t) dW(t) = \int_0^T W(t) dW(t) + T
    $$

---
<br>

## Part III: Advanced Applications of Stochastic Integration and Itô's Formula

### Exercise 1
**Problem**
Let $W(t)$ be an $n$-dimensional Brownian motion. Prove: $\mathbb{E}[|W(t) - W(s)|^4] = (2n + n^2)(t-s)^2$.

??? success "Solution (click to expand)"

    This problem cleverly utilizes the relationship between the standard normal distribution and the chi-squared distribution.

    Given that $W(t)$ is an $n$-dimensional Brownian motion, we examine its increments over the time interval $[s, t]$.
    Let $X_i = W_i(t) - W_i(s)$, where $i = 1, 2, \dots, n$ denotes each dimension.
    According to the properties of Brownian motion, the increments in each dimension are independent and identically distributed, and $X_i \sim N(0, t-s)$.

    For standardization, let $Z_i = \frac{X_i}{\sqrt{t-s}}$, then $Z_i \sim N(0, 1)$, and they are mutually independent.
    Let $Q = \sum_{i=1}^n Z_i^2$, by definition, $Q$ follows a chi-squared distribution with $n$ degrees of freedom, i.e., $Q \sim \chi^2(n)$.

    For the chi-squared distribution $\chi^2(n)$, we know its expectation and variance are:

    $$
    \mathbb{E}[Q] = n, \quad Var(Q) = 2n
    $$

    From this, we can find the second raw moment of $Q$:

    $$
    \mathbb{E}[Q^2] = Var(Q) + (\mathbb{E}[Q])^2 = 2n + n^2
    $$

    Returning to the original expression, we examine the fourth moment of the increment:

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
Define the second-order integral $\int_0^T f(t) [dW(t)]^2 \doteq \lim_{\delta \to 0} \sum_{i=0}^{n-1} f(t_i)[W(t_{i+1}) - W(t_i)]^2$
where $0 = t_0 < t_1 < \dots < t_n = T$, and $\delta \doteq \max_i |t_{i+1} - t_i|$. Prove that when $f \in L^2([0,T])$, we have

$$\int_0^T f(t) [dW(t)]^2 = \int_0^T f(t) dt$$

??? success "Solution (click to expand)"

    This problem aims to prove the strict validity of the quadratic variation of Brownian motion $(dW_t)^2 = dt$ in the integral sense.
    We need to prove this limit in the mean square ($L^2$) sense.

    Let $\Delta t_i = t_{i+1} - t_i$, $\Delta W_i = W(t_{i+1}) - W(t_i)$.
    Let the left-hand discrete sum be $S_n = \sum_{i=0}^{n-1} f(t_i)(\Delta W_i)^2$, and the right-hand target integral be $I = \int_0^T f(t)dt \approx \sum_{i=0}^{n-1} f(t_i)\Delta t_i$.

    Examine the mean square error of their difference:

    $$
    \mathbb{E}\left[ \left( S_n - \sum_{i=0}^{n-1} f(t_i)\Delta t_i \right)^2 \right] = \mathbb{E}\left[ \left( \sum_{i=0}^{n-1} f(t_i) ((\Delta W_i)^2 - \Delta t_i) \right)^2 \right]
    $$

    Since Brownian increments over non-overlapping intervals are independent, the expectation of cross terms ($i \ne j$) can be factored:

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

    Substitute back into the original expression, and bound one $\Delta t_i$ by the maximum step size $\delta$:

    $$
    = \sum_{i=0}^{n-1} f^2(t_i) 2\Delta t_i^2 \le 2\delta \sum_{i=0}^{n-1} f^2(t_i)\Delta t_i
    $$

    As $\delta \to 0$, since $f \in L^2([0,T])$, $\sum f^2(t_i)\Delta t_i \to \int_0^T f^2(t)dt < \infty$.
    Therefore, the preceding factor $\delta$ will drive the entire limit to $0$. Mean square convergence is proven, i.e.:

    $$
    \int_0^T f(t) [dW(t)]^2 = \int_0^T f(t) dt
    $$

---

### Exercise 3
**Problem**
Let $f \in L^2([0, T])$ and $\int_0^T f(s) dW(s) = 0$. Prove: $f$ is almost everywhere zero.

??? success "Solution (click to expand)"

    This problem utilizes the Itô Isometry and the fundamental properties of Lebesgue integration in real analysis.
    
    Since $\int_0^T f(s) dW(s) = 0$ holds almost surely, the second moment of this random variable must be 0:
    
    $$
    \mathbb{E}\left[ \left( \int_0^T f(s) dW(s) \right)^2 \right] = \mathbb{E}[0^2] = 0
    $$
    
    On the other hand, according to the Itô Isometry, the second moment of the stochastic integral equals the Lebesgue integral of the square of the integrand:
    
    $$
    \mathbb{E}\left[ \left( \int_0^T f(s) dW(s) \right)^2 \right] = \mathbb{E}\left[ \int_0^T f^2(s) ds \right]
    $$
    
    Since $f(s)$ is a deterministic function, its expectation is itself. Combining the two equations yields:
    
    $$
    \int_0^T f^2(s) ds = 0
    $$
    
    By the fundamental theorem of real analysis, since the integrand $f^2(s) \ge 0$ always holds and its Lebesgue integral over $[0,T]$ is $0$, this implies that $f^2(s)$ must be zero almost everywhere (a.e.) on $[0,T]$.
    
    Consequently: $f(s) = 0$ holds almost everywhere. Q.E.D.

---

### Exercise 4
**Problem**
Prove: $Y(t) = e^{t/2} \cos(W(t))$ is a martingale.

??? success "Solution (click to expand)"

    To prove that a process is a martingale, the most direct method is to use Itô's formula to show that its differential term has no drift term (i.e., the coefficient of the $dt$ term is 0).

    Let the bivariate function $u(t, x) = e^{t/2} \cos(x)$, and compute its partial derivatives with respect to $t$ and $x$:

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
    Therefore, the process $Y(t)$ is also a martingale. Q.E.D.

---

### Exercise 5
**Problem**
Prove:
1. $\int_0^T W^2 dW = \frac{1}{3}W(T)^3 - \int_0^T W dt$
2. $\int_0^T W^3 dW = \frac{1}{4}W(T)^4 - \frac{3}{2} \int_0^T W^2 dt$

??? success "Solution (click to expand)"

    Both parts of the proof are basic reverse applications of Itô's formula, i.e., "first guess the higher-order term, then expand using Itô's formula and rearrange terms."

    **(1) Proving the first identity**

    Let the function $f(x) = \frac{1}{3}x^3$. Compute its derivatives: $f'(x) = x^2$, $f''(x) = 2x$.
    Substitute $W(t)$ into Itô's formula:

    $$
    d\left( \frac{1}{3}W^3(t) \right) = W^2(t) dW(t) + \frac{1}{2}(2W(t)) (dW(t))^2
    $$

    Using the quadratic variation rule $(dW(t))^2 = dt$:

    $$
    d\left( \frac{1}{3}W^3(t) \right) = W^2(t) dW(t) + W(t) dt
    $$

    Integrate both sides over $[0, T]$, noting that $W(0) = 0$:

    $$
    \frac{1}{3}W^3(T) - 0 = \int_0^T W^2(t) dW(t) + \int_0^T W(t) dt
    $$

    Rearranging gives:

    $$
    \int_0^T W^2(t) dW(t) = \frac{1}{3}W^3(T) - \int_0^T W(t) dt
    $$



    **(2) Proving the second identity**

    Similarly, let the function $g(x) = \frac{1}{4}x^4$. Compute its derivatives: $g'(x) = x^3$, $g''(x) = 3x^2$.
    Substitute into Itô's formula:

    $$
    d\left( \frac{1}{4}W^4(t) \right) = W^3(t) dW(t) + \frac{1}{2}(3W^2(t)) dt
    $$

    Integrate both sides over $[0, T]$:

    $$
    \frac{1}{4}W^4(T) - 0 = \int_0^T W^3(t) dW(t) + \frac{3}{2}\int_0^T W^2(t) dt
    $$

    Rearranging gives:

    $$
    \int_0^T W^3(t) dW(t) = \frac{1}{4}W^4(T) - \frac{3}{2}\int_0^T W^2(t) dt
    $$

---

### Exercise 6
**Problem**
Prove that $\mathbb{E}[e^{\int_0^T g dW}] = e^{\frac{1}{2} \int_0^T g^2 ds}$.

??? success "Solution (click to expand)"

    This problem can be elegantly solved by analyzing the distributional properties of the Itô integral, directly utilizing the moment generating function of the normal distribution.

    Let the random variable $X = \int_0^T g(s) dW(s)$.
    Since $g(s)$ is a deterministic function of time (non-random), this Itô integral is a linear superposition of a Gaussian process, therefore $X$ still follows a normal distribution.

    According to the properties of stochastic integrals:

    * Expectation: $\mathbb{E}[X] = \mathbb{E}[\int_0^T g dW] = 0$
    
    * Variance: By the Itô isometry, $Var(X) = \mathbb{E}[X^2] = \int_0^T g^2(s) ds$

    Thus $X \sim N(0, \sigma^2)$, where $\sigma^2 = \int_0^T g^2(s) ds$.

    The required $\mathbb{E}[e^X]$ in the original problem is precisely the value of the moment generating function $M_X(u) = \mathbb{E}[e^{uX}]$ of the random variable $X$ at $u=1$.
    For a normal distribution $N(\mu, \sigma^2)$, its moment generating function formula is $M_X(u) = \exp(\mu u + \frac{1}{2}\sigma^2 u^2)$.

    Substituting $\mu = 0, u = 1, \sigma^2 = \int_0^T g^2 ds$, we immediately obtain:

    $$
    \mathbb{E}[e^{\int_0^T g dW}] = \exp\left( 0 + \frac{1}{2} \int_0^T g^2(s) ds \cdot 1^2 \right) = e^{\frac{1}{2} \int_0^T g^2 ds}
    $$

    Q.E.D.

---

### Exercise 7
**Problem**
Let $u = u(x, t)$ satisfy the parabolic partial differential equation $u_t + \frac{1}{2}u_{xx} = 0$. Prove: $\mathbb{E}[u(W(t), t)] = u(0, 0)$.

??? success "Solution (click to expand)"

    This problem is a classic exercise establishing the profound connection between partial differential equations (PDEs) and stochastic processes (SDEs), i.e., the simplest form of the Feynman-Kac formula.

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

    Substitute back the definition of $Y(t)$, and note that at the initial time, Brownian motion $W(0) = 0$ almost surely:

    $$
    \mathbb{E}[u(W(t), t)] = \mathbb{E}[u(W(0), 0)] = u(0, 0)
    $$

    The conclusion is proven.

---

### Exercise 8
**Problem**
1. Prove that $e^{W(t)} = 1 + \frac{1}{2}\int_0^t e^{W(s)} ds + \int_0^t e^{W(s)} dW(s)$;
2. Prove that $\mathbb{E}[e^{W(t)}] = 1 + \frac{1}{2}\int_0^t \mathbb{E}[e^{W(s)}] ds$, and hence $\mathbb{E}[e^{W(t)}] = e^{t/2}$;
3. Compute $\mathbb{E}[e^{iW(t)}]$, and the variances of $e^{W(t)}, \sin W(t), \cos W(t)$.

??? success "Solution (click to expand)"

    **(1) SDE Verification**
    Let $f(x) = e^x$. Substituting $W(t)$ into Itô's formula $df(W_t) = f'(W_t)dW_t + \frac{1}{2}f''(W_t)dt$:
    
    $$
    d(e^{W(t)}) = e^{W(t)}dW(t) + \frac{1}{2}e^{W(t)}dt
    $$
    
    Integrating both sides over $[0,t]$, and using $e^{W(0)} = e^0 = 1$:
    
    $$
    e^{W(t)} - 1 = \int_0^t e^{W(s)} dW(s) + \frac{1}{2}\int_0^t e^{W(s)} ds
    $$
    
    Rearranging the terms proves the first part.

    <br>

    **(2) Solving the ODE for Expectation**
    Taking the expectation on both sides of the result from part (1). Since the Itô integral $\int_0^t e^{W(s)} dW(s)$ has expectation 0 under appropriate regularity conditions, and using Fubini's theorem to interchange the expectation and the Riemann integral:
    
    $$
    \mathbb{E}[e^{W(t)}] = 1 + \frac{1}{2}\int_0^t \mathbb{E}[e^{W(s)}] ds
    $$
    
    Let $m(t) = \mathbb{E}[e^{W(t)}]$, the above equation transforms into the integral equation $m(t) = 1 + \frac{1}{2}\int_0^t m(s) ds$.
    Differentiating it yields the initial value problem ODE:
    
    $$
    m'(t) = \frac{1}{2}m(t), \quad m(0) = 1
    $$
    
    Solving this ordinary differential equation immediately gives:
    
    $$
    m(t) = \mathbb{E}[e^{W(t)}] = e^{t/2}
    $$

    <br>

    **(3) Characteristic Function and Trigonometric Function Variance Calculation**

    *Compute $\mathbb{E}[e^{iW(t)}]$ (characteristic function):*
    Again, apply Itô's formula to the complex-valued process $Y(t) = e^{iW(t)}$:

    $$
    dY(t) = i e^{iW(t)} dW(t) + \frac{1}{2}(i)^2 e^{iW(t)} dt = i Y(t) dW(t) - \frac{1}{2} Y(t) dt
    $$

    Taking expectation and differentiating, let $m_2(t) = \mathbb{E}[Y(t)]$, we obtain the ODE: $m_2'(t) = -\frac{1}{2}m_2(t)$ with $m_2(0)=1$.
    Solving:

    $$
    \mathbb{E}[e^{iW(t)}] = e^{-t/2}
    $$

    *Compute the variance of $e^{W(t)}$:*
    By the definition of variance, $Var(e^{W(t)}) = \mathbb{E}[(e^{W(t)})^2] - (\mathbb{E}[e^{W(t)}])^2 = \mathbb{E}[e^{2W(t)}] - e^t$.
    Treating $2W(t)$ as the case with parameter $2$, from the moment generating function result we have $\mathbb{E}[e^{2W(t)}] = e^{4t/2} = e^{2t}$.

    $$
    Var(e^{W(t)}) = e^{2t} - e^t
    $$

    *Compute the variances of $\sin W(t), \cos W(t)$:*
    By Euler's formula, $\mathbb{E}[e^{iW(t)}] = \mathbb{E}[\cos W(t)] + i\mathbb{E}[\sin W(t)] = e^{-t/2}$. Comparing real and imaginary parts:

    $$
    \mathbb{E}[\cos W(t)] = e^{-t/2}, \quad \mathbb{E}[\sin W(t)] = 0
    $$

    Similarly, letting the parameter be $2i$, we have $\mathbb{E}[e^{2iW(t)}] = e^{-(2i)^2 t/(-2)} = e^{-2t}$, i.e.:

    $$
    \mathbb{E}[\cos 2W(t)] = e^{-2t}, \quad \mathbb{E}[\sin 2W(t)] = 0
    $$

    Using double-angle formulas to reduce powers:

    $$
    \mathbb{E}[\sin^2 W(t)] = \mathbb{E}\left[ \frac{1 - \cos 2W(t)}{2} \right] = \frac{1 - e^{-2t}}{2}
    $$

    $$
    \mathbb{E}[\cos^2 W(t)] = \mathbb{E}\left[ \frac{1 + \cos 2W(t)}{2} \right] = \frac{1 + e^{-2t}}{2}
    $$

    Finally, substituting into the variance formula:

    $$
    Var(\sin W(t)) = \mathbb{E}[\sin^2] - (\mathbb{E}[\sin])^2 = \frac{1 - e^{-2t}}{2} - 0 = \frac{1 - e^{-2t}}{2}
    $$

    $$
    Var(\cos W(t)) = \mathbb{E}[\cos^2] - (\mathbb{E}[\cos])^2 = \frac{1 + e^{-2t}}{2} - (e^{-t/2})^2 = \frac{1 - e^{-2t}}{2}
    $$

    *(Note: The manuscript here cleverly utilizes known results for trigonometric moments of the normal distribution.)*