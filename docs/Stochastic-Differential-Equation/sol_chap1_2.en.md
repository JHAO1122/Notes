---
tags:
  - Probability Theory
  - Stochastic Processes
  - Stochastic Differential Equations
  - Exercise Solutions
---

# 📝 Selected Exercise Solutions: Conditional Expectation & Brownian Motion Properties

!!! abstract "About This Page"
    This page contains solutions to core exercises from Chapter 1 (Independence, Conditional Expectation) and Chapter 2 (Brownian Motion and its Properties) of the *Stochastic Differential Equations* course. All solutions are presented in collapsible blocks; click to expand and view the detailed derivations.

---

## Part I: Independence and Conditional Expectation

### Exercise 1: Density Function and Change of Variable for a 1D Continuous Random Variable

**Problem:**
Let the probability density function of a random variable $X$ be $f(x) = ax(1-x)$ for $x \in (0,1)$, and $f(x)=0$ otherwise.
(1) Find the constant $a$;
(2) Let $Y = X^3$, find the probability density function of $Y$.

??? success "Solution (Click to expand)"
    
    **(1) Finding the constant $a$**
    
    By the normalization property of probability density functions, the integral over the entire space must equal 1:
    
    $$
    \int_{0}^{1} ax(1-x)dx = 1 
    $$

    Evaluating this integral:

    $$
    a \left( \frac{1}{2}x^2 - \frac{1}{3}x^3 \right) \Bigg|_{0}^{1} = 1 \implies a \left( \frac{1}{2} - \frac{1}{3} \right) = 1
    $$
    
    Solving for $a$ gives: $a = 6$.

    <br>

    **(2) Finding the probability density function of $Y=X^3$**
    
    First, we find the cumulative distribution function (CDF) of $Y$, denoted as $F_Y(y)$. Since $X \in (0,1)$, we have $Y = X^3 \in (0,1)$.
    For $y \in (0,1)$:
    
    $$
    F_Y(y) = \mathbb{P}(Y \leqslant y) = \mathbb{P}(X^3 \leqslant y) = \mathbb{P}(X \leqslant y^{1/3})
    $$
    
    Substitute the PDF of $X$ to compute the integral:

    $$
    F_Y(y) = \int_{0}^{y^{1/3}} 6x(1-x)dx = \left( 3x^2 - 2x^3 \right) \Big|_{0}^{y^{1/3}} = 3y^{2/3} - 2y
    $$
    
    Differentiating $F_Y(y)$ yields the probability density function $f_Y(y)$ of $Y$:
    
    $$
    f_Y(y) = F_Y'(y) = 3 \cdot \frac{2}{3}y^{-1/3} - 2 = 2y^{-1/3} - 2
    $$
    
    In summary, the PDF of $Y$ is:
    
    $$
    f_Y(y) = \begin{cases} 2y^{-1/3} - 2, & y \in (0,1) \\ 0, & \text{otherwise} \end{cases}
    $$

---

### Exercise 4: Properties of Conditional Expectation for Independent Random Variables

**Problem:**
Suppose $X$ and $Y$ are independent random variables, $f(x,y)$ is a bounded continuous function, and $F_X$ is the probability distribution function (CDF) of $X$. Prove that:
(1) $\mathbb{E}[f(X,Y)|Y] = \int f(x,Y)dF_X(x)$;
(2) $\mathbb{P}(X+Y \leqslant x | Y) = F_X(x-Y)$.

??? success "Solution (Click to expand)"

    **(1) Proof:**
    
    We use the standard measure-theoretic approximation method (the "Standard Machine") for a rigorous proof:
    
    **Step 1 (Indicator functions):** Let $f(x,y) = I_A(x)I_B(y)$, where $A, B$ are Borel sets.
    Since $X$ and $Y$ are independent, we use the properties of conditional expectation (taking out what is known given $Y$, and the conditional expectation of an independent variable equals its unconditional expectation):
    
    $$
    \mathbb{E}[I_A(X)I_B(Y)|Y] = I_B(Y)\mathbb{E}[I_A(X)|Y] = I_B(Y)\mathbb{E}[I_A(X)] = I_B(Y)\mathbb{P}(X \in A)
    $$
    
    On the other hand, for the integral form on the right-hand side:
    
    $$
    \int I_A(x)I_B(Y)dF_X(x) = I_B(Y) \int I_A(x)dF_X(x) = I_B(Y)\mathbb{P}(X \in A)
    $$
    
    The two sides are equal, so the proposition holds for indicator functions.
    
    **Step 2 (Simple functions):** By the linearity of expectation, the conclusion holds for simple functions (linear combinations of finitely many indicator functions).
    
    **Step 3 (Non-negative measurable functions):** For any non-negative measurable function, there exists a monotonically increasing sequence of simple functions that approximates it. By the Monotone Convergence Theorem, the conclusion holds for non-negative measurable functions.
    
    **Step 4 (Bounded continuous functions):** Any bounded continuous function can be decomposed into its positive and negative parts $f = f^+ - f^-$, both of which are Lebesgue integrable by definition. Thus, the equation holds for all bounded continuous functions. The proof is complete.

    <br>

    **(2) Proof:**
    
    Conditional probability can be directly written in the form of conditional expectation:
    
    $$
    \mathbb{P}(X+Y \leqslant x | Y) = \mathbb{E}[I_{\{X+Y \leqslant x\}} | Y]
    $$
    
    Let $g(X, Y) = I_{\{X+Y \leqslant x\}}$. Using the conclusion proven in part (1):
    
    $$
    \mathbb{E}[I_{\{X+Y \leqslant x\}} | Y] = \int I_{\{u+Y \leqslant x\}} dF_X(u) 
    $$

    Perform an equivalent transformation on the integration domain $u \leqslant x - Y$:

    $$
    \int I_{\{u \leqslant x-Y\}} dF_X(u) = \int_{-\infty}^{x-Y} dF_X(u) = F_X(x-Y)
    $$
    
    The proof is complete.

---

### Exercise 5: Calculating Conditional Expectation Based on Events

**Problem:**
Let $A, B$ be two events in a measure space $(\Omega, \mathcal{F}, \mathbb{P})$. Calculate the conditional expectation $\mathbb{E}[\chi_A|\chi_B]$. (Note: $\chi$ represents the indicator function).

??? success "Solution (Click to expand)"

    Since $\chi_B$ can only take values 0 or 1, the $\sigma$-algebra generated by it is very simple: $\sigma(\chi_B) = \{\emptyset, \Omega, B, B^c\}$.
    Therefore, the conditional expectation $\mathbb{E}[\chi_A|\chi_B]$ must be $\sigma(\chi_B)$-measurable, which means it must be constant on $B$ and $B^c$. Let it be:
    
    $$
    \mathbb{E}[\chi_A|\chi_B] = c_1 \chi_B + c_2 \chi_{B^c}
    $$
    
    According to the Radon-Nikodym derivative definition of conditional expectation, for any $\Lambda \in \sigma(\chi_B)$, the following integrals must be equal:
    
    $$
    \int_\Lambda \mathbb{E}[\chi_A|\chi_B] d\mathbb{P} = \int_\Lambda \chi_A d\mathbb{P}
    $$
    
    **Case 1: Let $\Lambda = B$**
    
    $$
    \int_B (c_1 \chi_B + c_2 \chi_{B^c}) d\mathbb{P} = \int_B \chi_A d\mathbb{P} \implies c_1 \mathbb{P}(B) = \mathbb{P}(A \cap B)
    $$
    
    If $\mathbb{P}(B) > 0$, then $c_1 = \frac{\mathbb{P}(A \cap B)}{\mathbb{P}(B)} = \mathbb{P}(A|B)$.
    
    **Case 2: Let $\Lambda = B^c$**
    
    $$
    \int_{B^c} (c_1 \chi_B + c_2 \chi_{B^c}) d\mathbb{P} = \int_{B^c} \chi_A d\mathbb{P} \implies c_2 \mathbb{P}(B^c) = \mathbb{P}(A \cap B^c)
    $$
    
    If $\mathbb{P}(B^c) > 0$, then $c_2 = \frac{\mathbb{P}(A \cap B^c)}{\mathbb{P}(B^c)} = \mathbb{P}(A|B^c)$.
    
    In conclusion, the explicit expression for the conditional expectation is:
    
    $$
    \mathbb{E}[\chi_A|\chi_B] = \mathbb{P}(A|B)\chi_B + \mathbb{P}(A|B^c)\chi_{B^c}
    $$

---

### Exercise 6: Independent $\sigma$-fields and Conditional Expectation

**Problem:**
Suppose $\mathcal{V}_1$ and $\mathcal{V}_2$ are two independent $\sigma$-fields, and $X$ is an integrable random variable. Prove that:
$\mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)|\mathcal{V}_2] = \mathbb{E}[X]$.

??? success "Solution (Click to expand)"

    **Proof:**
    
    First, by the definition of conditional expectation, the inner expectation $\mathbb{E}(X|\mathcal{V}_1)$ is necessarily a $\mathcal{V}_1$-measurable random variable.
    
    Since it is given that $\mathcal{V}_1$ and $\mathcal{V}_2$ are independent $\sigma$-fields, the $\mathcal{V}_1$-measurable random variable $\mathbb{E}(X|\mathcal{V}_1)$ and the $\sigma$-field $\mathcal{V}_2$ are naturally independent.
    
    According to the property of conditional expectation regarding independent $\sigma$-fields (if a random variable $Z$ is independent of a $\sigma$-field $\mathcal{G}$, then $\mathbb{E}[Z|\mathcal{G}] = \mathbb{E}[Z]$), we can drop the conditioning on $\mathcal{V}_2$:
    
    $$
    \mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)|\mathcal{V}_2] = \mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)]
    $$
    
    Finally, using the Law of Total Expectation (a degenerate form of the Tower Property):
    
    $$
    \mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)] = \mathbb{E}[X]
    $$
    
    Combining the two equations completes the proof. $\square$

---

### Supplementary Exercise 1: $L^p$ Continuity of Conditional Expectation

**Problem:**
Let $X$ and $\{X_n\}$ be a sequence of random variables on $(\Omega, \mathcal{F}, \mathbb{P})$. For a fixed $1 \leqslant p < \infty$, suppose $\mathbb{E}[|X_n|^p] < +\infty$. Assume $\lim_{n\to\infty}\mathbb{E}[|X_n - X|^p] = 0$. Prove that for any sub-$\sigma$-field $\mathcal{V} \subset \mathcal{F}$:
$$
\lim_{n\to\infty}\mathbb{E}\left[\big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p\right] = 0.
$$

??? success "Solution (Click to expand)"

    **Proof:**
    
    Consider the function $\Phi(x) = |x|^p$. Since $p \geqslant 1$, this is a convex function.
    Using the linearity of conditional expectation and the **Conditional Jensen's Inequality**, we merge the conditional expectations inside the absolute value, and bound it by moving the absolute value inside the conditional expectation:
    
    $$
    \big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p = \big|\mathbb{E}[X_n - X | \mathcal{V}]\big|^p \leqslant \mathbb{E}\big[|X_n - X|^p \big| \mathcal{V}\big]
    $$
    
    Taking the unconditional expectation on both sides of the inequality and applying the Law of Total Expectation:
    
    $$
    \mathbb{E}\left[\big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p\right] \leqslant \mathbb{E}\left[ \mathbb{E}\big[|X_n - X|^p \big| \mathcal{V}\big] \right] = \mathbb{E}\big[|X_n - X|^p\big]
    $$
    
    Taking the limit as $n \to \infty$ on both sides. Since it is given that $\lim_{n\to\infty}\mathbb{E}[|X_n - X|^p] = 0$, and the expectation on the left is always non-negative, by the squeeze theorem we obtain:
    
    $$
    \lim_{n\to\infty}\mathbb{E}\left[\big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p\right] = 0
    $$
    
    The proof is complete. $\square$

---

### Supplementary Exercise 2: Calculating Conditional Expectation on a Finite State Space

**Problem:**
Let $\Omega = \{1, 2, \cdots, 7, 8\}$ and $\mathcal{F} = 2^\Omega$ (the $\sigma$-field of all subsets of $\Omega$). When $i \leqslant 4$, $\mathbb{P}(\{i\}) = 1/10$; when $i > 4$, $\mathbb{P}(\{i\}) = 3/20$.
Define $X = \chi_{\{1,2,3,4\}} + 2\chi_{\{5,6,7,8\}}$ and $Y = \chi_{\{1,5\}} + 2\chi_{\{2,3,4,6,7,8\}}$.
Let $\mathcal{V}$ be the $\sigma$-field generated by $\{1,2\}$ and $\{3,4\}$, and $\mathcal{H}$ be the $\sigma$-field generated by $\{1,2,3,4\}$.
Calculate: (1) $X\mathbb{E}[Y]$; (2) $\mathbb{E}[\mathbb{E}[XY|\mathcal{V}]|\mathcal{H}]$.

??? success "Solution (Click to expand)"

    **Preparation: Clarifying the probability structure**
    
    * $\mathbb{P}(\{1,2,3,4\}) = 4 \times \frac{1}{10} = \frac{2}{5}$
    * $\mathbb{P}(\{5,6,7,8\}) = 4 \times \frac{3}{20} = \frac{3}{5}$
    * It is easy to verify that the two core atoms of $\mathcal{H} = \sigma(\{1,2,3,4\})$ are $\{1,2,3,4\}$ and $\{5,6,7,8\}$. The partition of $\mathcal{V} = \sigma(\{1,2\}, \{3,4\})$ is finer, so clearly $\mathcal{H} \subset \mathcal{V}$.

    <br>

    **(1) Calculating $X\mathbb{E}[Y]$**
    
    First, compute the unconditional expectation of $Y$:
    
    $$
    \mathbb{E}[Y] = 1 \cdot \mathbb{P}(\{1,5\}) + 2 \cdot \mathbb{P}(\{2,3,4,6,7,8\})
    $$
    
    $$
    = 1 \cdot \left(\frac{1}{10} + \frac{3}{20}\right) + 2 \cdot \left(3 \times \frac{1}{10} + 3 \times \frac{3}{20}\right) = \frac{1}{4} + 2\left(\frac{3}{4}\right) = \frac{7}{4}
    $$
    
    Multiplying it directly by $X$:
    
    $$
    X\mathbb{E}[Y] = \frac{7}{4}X = \frac{7}{4}\chi_{\{1,2,3,4\}} + \frac{14}{4}\chi_{\{5,6,7,8\}}
    $$

    <br>

    **(2) Calculating $\mathbb{E}[\mathbb{E}[XY|\mathcal{V}]|\mathcal{H}]$**
    
    Here we use the elegant Tower Property (Smoothing Property) to directly reduce the dimensionality, avoiding tedious double conditional expectation calculations.
    Since $\mathcal{H} \subset \mathcal{V}$, the smaller $\sigma$-algebra dominates on the outside:
    
    $$
    \mathbb{E}\big[\mathbb{E}[XY|\mathcal{V}]\big|\mathcal{H}\big] = \mathbb{E}[XY|\mathcal{H}]
    $$
    
    Observe the random variable $X$: it is constantly 1 on $\{1,2,3,4\}$, and constantly 2 on $\{5,6,7,8\}$, which perfectly matches the atomic partition of $\mathcal{H}$. Thus, $X$ is **$\mathcal{H}$-measurable**! As known information, it can be pulled out of the conditional expectation:
    
    $$
    \mathbb{E}[XY|\mathcal{H}] = X \mathbb{E}[Y|\mathcal{H}]
    $$
    
    Next, we calculate the conditional expectation of $Y$ on the two atoms of $\mathcal{H}$ separately:
    * On $\{1,2,3,4\}$: $\mathbb{E}[Y|\{1,2,3,4\}] = \frac{1\cdot P(\{1\}) + 2\cdot P(\{2,3,4\})}{P(\{1,2,3,4\})} = \frac{1/10 + 6/10}{4/10} = \frac{7}{4}$
    * On $\{5,6,7,8\}$: $\mathbb{E}[Y|\{5,6,7,8\}] = \frac{1\cdot P(\{5\}) + 2\cdot P(\{6,7,8\})}{P(\{5,6,7,8\})} = \frac{3/20 + 18/20}{12/20} = \frac{7}{4}$
    
    Surprisingly, $\mathbb{E}[Y|\mathcal{H}] = 7/4$ holds constantly across both partitions. Therefore:
    
    $$
    \mathbb{E}[\mathbb{E}[XY|\mathcal{V}]|\mathcal{H}] = \frac{7}{4}X
    $$

---

### Supplementary Exercise 3: Constructing a Doob Martingale

**Problem:**
Given a probability space $(\Omega, \mathcal{F}, \mathbb{P})$ and an integrable random variable $X$. Let $\{\mathcal{F}(t)\}_{t \geqslant 0}$ be a filtration of $\sigma$-fields. For $t \geqslant 0$, define $X(t) \doteq \mathbb{E}[X|\mathcal{F}(t)]$. Prove that $X(t)$ is a martingale with respect to $\mathcal{F}(t)$.

??? success "Solution (Click to expand)"

    **Proof:**
    We need to strictly verify the three core conditions of a martingale:
    
    **1. Integrability:**
    Using the contractive property of conditional expectation (a special case of Jensen's inequality for absolute values):
    
    $$
    \mathbb{E}|X(t)| = \mathbb{E}\big|\mathbb{E}[X|\mathcal{F}(t)]\big| \leqslant \mathbb{E}\big[\mathbb{E}[|X| \big| \mathcal{F}(t)]\big] = \mathbb{E}|X|
    $$
    
    Since it is given that $X$ is integrable ($\mathbb{E}|X| < \infty$), we have $\mathbb{E}|X(t)| < \infty$.
    
    **2. Adaptability:**
    By the measure-theoretic definition of conditional expectation $\mathbb{E}[X|\mathcal{F}(t)]$, it is inherently $\mathcal{F}(t)$-measurable. Thus, the process $\{X(t)\}$ is adapted to the filtration $\{\mathcal{F}(t)\}_{t \ge 0}$.
    
    **3. Martingale Property:**
    For any $0 \leqslant s \leqslant t$, since $\{\mathcal{F}(t)\}$ is a filtration (i.e., information accumulates over time), we must have $\mathcal{F}(s) \subset \mathcal{F}(t)$.
    By the Tower Property of conditional expectation:
    
    $$
    \mathbb{E}[X(t) | \mathcal{F}(s)] = \mathbb{E}\big[\mathbb{E}[X|\mathcal{F}(t)] \big| \mathcal{F}(s)\big] = \mathbb{E}[X|\mathcal{F}(s)] = X(s)
    $$
    
    Combining these three points, the process $X(t)$, formed by projecting a single integrable variable onto different filtrations, is necessarily a martingale (this is known in stochastic analysis as a Doob martingale). $\square$

---
<br>

## Part II: Brownian Motion and its Properties

### P64 Exercise 1: Translation and Scaling Invariance of Brownian Motion

**Problem:**
Let $W(t)$ be a one-dimensional Brownian motion. Prove that for any fixed $s>0$, $W(t+s)-W(s)$ is a Brownian motion; and for any positive constant $c$, $cW(t/c^2)$ is also a Brownian motion.

??? success "Solution (Click to expand)"

    **(1) Prove that $B(t) = W(t+s) - W(s)$ is a Brownian motion**
    
    We verify the four defining properties of Brownian motion one by one:
    
    1. **Initial value of 0**: $B(0) = W(0+s) - W(s) = 0$.
    2. **Independent increments**: For any $0 \le t_1 < t_2 < \cdots < t_k$, the increment sequence is:
       $$B(t_j) - B(t_{j-1}) = W(t_j+s) - W(t_{j-1}+s)$$
       Since the original process $W$ has independent increments, and the time intervals $[t_{j-1}+s, t_j+s]$ do not overlap on the time axis, these increments are mutually independent.
    3. **Stationary normal increments**: For $t > u$, the distribution of the increment is:
       $$B(t) - B(u) = W(t+s) - W(u+s) \sim N(0, (t+s) - (u+s)) = N(0, t-u)$$
    4. **Path continuity**: Since the sample paths of $W(\cdot)$ are almost surely continuous, the translated paths $B(\cdot)$ are also almost surely continuous.
    
    In conclusion, $B(t)$ is a one-dimensional Brownian motion.

    <br>

    **(2) Prove that $U(t) = cW(t/c^2)$ is a Brownian motion**
    
    Similarly, we verify the properties:
    
    1. **Initial value of 0**: $U(0) = cW(0/c^2) = cW(0) = 0$.
    2. **Independent increments**: For any $0 \le t_1 < t_2 < \cdots < t_k$, since the times $t_1/c^2 < t_2/c^2 < \cdots < t_k/c^2$ do not overlap, the increments of the original Brownian motion are mutually independent. Multiplying by a constant $c$ maintains this independence.
    3. **Stationary normal increments**: By the properties of linear transformations, the mean remains $0$, and the variance is:
       $$Var(U(t) - U(s)) = Var\left( c(W(t/c^2) - W(s/c^2)) \right) = c^2 \cdot \left( \frac{t-s}{c^2} \right) = t-s$$
       Thus, $U(t) - U(s) \sim N(0, t-s)$.
    4. **Path continuity**: Scaling does not change the continuity of the sample paths.
    
    In conclusion, $U(t) = cW(t/c^2)$ is also a Brownian motion (this is known as the scaling invariance or fractal property of Brownian motion).

---

### P64 Exercise 2: Time-Reversed Brownian Motion

**Problem:**
Let $W(t)$ be a one-dimensional Brownian motion. Let $\tilde{W}(t) = \begin{cases} tW(1/t), & t > 0, \\ 0, & t = 0. \end{cases}$
Prove that: $\tilde{W}(t) - \tilde{W}(s) \sim N(0, t-s), \forall 0 < s < t$.

??? success "Solution (Click to expand)"

    For $0 < s < t$, we identically transform the increment to create mutually independent terms:

    $$
    \tilde{W}(t) - \tilde{W}(s) = tW\left(\frac{1}{t}\right) - sW\left(\frac{1}{s}\right)
    $$

    We separate the terms by adding and subtracting $sW(1/t)$:

    $$
    = (t-s)W\left(\frac{1}{t}\right) - s\left(W\left(\frac{1}{s}\right) - W\left(\frac{1}{t}\right)\right)
    $$

    Note that here $1/t < 1/s$. Using the independent increments property of Brownian motion, the random variable $W(1/t) = W(1/t) - W(0)$ is mutually independent from the increment $W(1/s) - W(1/t)$.
    
    Since they are independent, their variances can be directly added:

    $$
    Var(\tilde{W}(t) - \tilde{W}(s)) = (t-s)^2 Var\left( W\left(\frac{1}{t}\right) \right) + s^2 Var\left( W\left(\frac{1}{s}\right) - W\left(\frac{1}{t}\right) \right)
    $$

    Substituting the variance formula $Var(W(u)) = u$:

    $$
    = (t-s)^2 \left(\frac{1}{t}\right) + s^2 \left(\frac{1}{s} - \frac{1}{t}\right)
    $$

    Finding a common denominator and simplifying:

    $$
    = \frac{t^2 - 2ts + s^2}{t} + s - \frac{s^2}{t} = \frac{t^2 - 2ts}{t} + s = t - 2s + s = t - s
    $$

    Also, by the linearity of expectation, $E[\tilde{W}(t) - \tilde{W}(s)] = 0 - 0 = 0$.
    Since it is a linear combination of jointly normal random variables, it must be normally distributed.
    
    The conclusion is proven: $\tilde{W}(t) - \tilde{W}(s) \sim N(0, t-s)$.

---

### P70 Exercise 1: Higher-Order Moments of Brownian Motion

**Problem:**
Let $W(t)$ be a Brownian motion. Prove that: $\mathbb{E}[W^{2k}(t)] = \frac{(2k)!t^k}{2^k k!}, \forall t > 0$.

??? success "Solution (Click to expand)"

    Since $W(t) \sim N(0, t)$, we use the Moment Generating Function (MGF) of the normal distribution $M_W(\lambda) = \mathbb{E}[e^{\lambda W(t)}]$:

    $$
    M_W(\lambda) = \exp\left( \frac{1}{2} \lambda^2 t \right)
    $$

    We expand both sides of the equation into Taylor Series at $\lambda = 0$.
    
    The left side expands via the linearity of expectation:
    
    $$
    \mathbb{E}[e^{\lambda W(t)}] = \sum_{m=0}^{\infty} \mathbb{E}[W^m(t)] \frac{\lambda^m}{m!}
    $$
    
    The right side expands the exponential function:
    
    $$
    \exp\left( \frac{1}{2} \lambda^2 t \right) = \sum_{k=0}^{\infty} \frac{1}{k!} \left( \frac{1}{2} \lambda^2 t \right)^k = \sum_{k=0}^{\infty} \frac{t^k}{2^k k!} \lambda^{2k}
    $$
    
    Equating the coefficients of $\lambda^{2k}$ on both sides:
    
    $$
    \frac{\mathbb{E}[W^{2k}(t)]}{(2k)!} = \frac{t^k}{2^k k!}
    $$
    
    Rearranging terms yields the proof:
    
    $$
    \mathbb{E}[W^{2k}(t)] = \frac{(2k)! t^k}{2^k k!}
    $$
    
    *(Note: For odd moments, equating coefficients shows that since there are no odd-power terms on the right side, $\mathbb{E}[W^{2k+1}(t)] = 0$)*.

---

### P70 Exercise 2: Expectation of an Exponential Brownian Increment

**Problem:**
Prove that: let $c$ be a constant and $0 < s < t$, then $\mathbb{E}[\exp(c(W(s) - W(t)))] = \exp\left(\frac{1}{2}c^2(t-s)\right)$.

??? success "Solution (Click to expand)"

    Due to the stationarity and symmetry of Brownian increments, $W(s) - W(t)$ and $W(t) - W(s)$ have the same distribution, both following $N(0, t-s)$.
    
    This is equivalent to finding the moment generating function of a normal random variable.
    Let $Z = c(W(s) - W(t))$, then $Z \sim N(0, c^2(t-s))$.
    
    According to the MGF formula for a normal distribution $N(\mu, \sigma^2)$, which is $\mathbb{E}[e^Z] = \exp(\mu + \frac{1}{2}\sigma^2)$:
    
    $$
    \mathbb{E}[\exp(c(W(s) - W(t)))] = \exp\left( 0 + \frac{1}{2} \cdot c^2(t-s) \right) = \exp\left( \frac{1}{2}c^2(t-s) \right)
    $$
    
    The proof is complete.

---

### P70 Exercise 3: Covariance of the Ornstein-Uhlenbeck (O-U) Process

**Problem:**
Let $U(t) = e^{-t}W(e^{2t})$, show that $\mathbb{E}[U(t)U(s)] = e^{-|t-s|}, \forall t, s \in \mathbb{R}$.

??? success "Solution (Click to expand)"

    Without loss of generality, let us assume $t \geqslant s$.
    
    Substitute the definition of $U(t)$ to calculate the expectation of the cross term:

    $$
    \mathbb{E}[U(t)U(s)] = \mathbb{E}\left[ e^{-t}W(e^{2t}) \cdot e^{-s}W(e^{2s}) \right] = e^{-(t+s)} \mathbb{E}\left[ W(e^{2t})W(e^{2s}) \right]
    $$

    For Brownian motion, we know its covariance function is $\mathbb{E}[W(u)W(v)] = \min(u, v)$.
    Since we assumed $t \geqslant s$, it naturally follows that $e^{2t} \geqslant e^{2s}$, thus:

    $$
    \mathbb{E}\left[ W(e^{2t})W(e^{2s}) \right] = \min(e^{2t}, e^{2s}) = e^{2s}
    $$

    Substitute this back into the original equation:

    $$
    \mathbb{E}[U(t)U(s)] = e^{-(t+s)} \cdot e^{2s} = e^{s-t}
    $$

    Since $t \geqslant s$, we have $s-t = -|t-s|$. If $s > t$, the conclusion holds symmetrically.
    Therefore, $\forall t, s \in \mathbb{R}$, $\mathbb{E}[U(t)U(s)] = e^{-|t-s|}$, completing the proof.

---

### P70 Exercise 4: Law of Large Numbers and Asymptotic Behavior of Brownian Motion at Infinity

**Problem:**
Prove that: almost surely $\lim_{m \to \infty} \frac{W(m)}{m} = 0$.

??? success "Solution (Click to expand)"

    We can construct a rigorous measure-theoretic proof using Chebyshev's Inequality combined with the Borel-Cantelli Lemma.
    
    **Step 1: Bounding the Variance**
    Consider the random sequence $Y_m = \frac{W(m)}{m}$. Since $W(m) \sim N(0, m)$, we have:

    $$
    Var\left(\frac{W(m)}{m}\right) = \frac{1}{m^2} Var(W(m)) = \frac{m}{m^2} = \frac{1}{m}
    $$

    **Step 2: Strengthening the Moment Bound (Preparing for Borel-Cantelli)**
    If we only use Chebyshev's inequality for the second moment, we obtain $P(|Y_m| \ge \epsilon) \le \frac{1}{m\epsilon^2}$. However, the harmonic series $\sum_{m=1}^{\infty} \frac{1}{m}$ diverges, so we cannot directly apply the B-C lemma.
    Therefore, we use the fourth moment (or utilize the exponential decay property of the normal distribution). For a standard normal variable $Z \sim N(0,1)$, $\mathbb{E}[Z^4] = 3$.
    
    $$
    \mathbb{E}[W(m)^4] = 3m^2 \implies \mathbb{E}\left[\left(\frac{W(m)}{m}\right)^4\right] = \frac{3m^2}{m^4} = \frac{3}{m^2}
    $$
    
    Applying Markov's inequality based on the fourth moment:
    
    $$
    \mathbb{P}\left(\left|\frac{W(m)}{m}\right| \geqslant \epsilon\right) = \mathbb{P}\left(\left(\frac{W(m)}{m}\right)^4 \geqslant \epsilon^4\right) \leqslant \frac{3}{m^2 \epsilon^4}
    $$

    **Step 3: Borel-Cantelli Lemma**
    For any fixed $\epsilon > 0$, since $\sum_{m=1}^{\infty} \frac{1}{m^2} < \infty$, the series converges absolutely:
    
    $$
    \sum_{m=1}^{\infty} \mathbb{P}\left(\left|\frac{W(m)}{m}\right| \geqslant \epsilon\right) \leqslant \sum_{m=1}^{\infty} \frac{3}{m^2 \epsilon^4} < \infty
    $$
    
    According to the first part of the Borel-Cantelli Lemma, the probability that the event $\left\{ \left|\frac{W(m)}{m}\right| \geqslant \epsilon \right\}$ occurs infinitely often is 0.
    This is equivalent to saying that almost all sample paths satisfy $\lim_{m \to \infty} \frac{W(m)}{m} = 0$. The proof is complete.

---

### P70 Exercise 5: Two Extremely Important Brownian Martingales

**Problem:**
Prove that $W(t)^2 - t$ and $\exp\left(\lambda W_t - \frac{1}{2}\lambda^2 t\right)$ $(\lambda \in \mathbb{R})$ are both martingales with respect to the history $\mathcal{F}(t)$ of $W(t)$.

??? success "Solution (Click to expand)"

    **Proof (1): $W(t)^2 - t$ is a martingale**
    
    Let $0 \leqslant s < t$. We need to prove that $\mathbb{E}[W(t)^2 - t | \mathcal{F}(s)] = W(s)^2 - s$.
    Expand the squared term using the identity $W(t) = (W(t) - W(s)) + W(s)$:

    $$
    \mathbb{E}[W(t)^2 | \mathcal{F}(s)] = \mathbb{E}[((W(t) - W(s)) + W(s))^2 | \mathcal{F}(s)]
    $$
    
    $$
    = \mathbb{E}[(W(t) - W(s))^2 | \mathcal{F}(s)] + 2\mathbb{E}[W(s)(W(t) - W(s)) | \mathcal{F}(s)] + \mathbb{E}[W(s)^2 | \mathcal{F}(s)]
    $$
    
    * Term 1: Due to independent increments of Brownian motion, $W(t)-W(s)$ is independent of $\mathcal{F}(s)$, so the conditional expectation equals the unconditional expectation $\mathbb{E}[(W(t)-W(s))^2] = t-s$.
    * Term 2: $W(s)$ is $\mathcal{F}(s)$-measurable and can be factored out. The expectation of the remaining increment is 0.
    * Term 3: $W(s)^2$ is $\mathcal{F}(s)$-measurable, so it is equal to itself.
    
    Therefore:
    
    $$
    \mathbb{E}[W(t)^2 | \mathcal{F}(s)] = (t-s) + 0 + W(s)^2
    $$
    
    Subtracting $t$:
    
    $$
    \mathbb{E}[W(t)^2 - t | \mathcal{F}(s)] = (t-s) + W(s)^2 - t = W(s)^2 - s
    $$
    
    The martingale property is proven.

    <br>

    **Proof (2): $\exp\left(\lambda W(t) - \frac{1}{2}\lambda^2 t\right)$ is a martingale (Geometric Brownian Martingale)**
    
    Again, let $0 \leqslant s < t$. Using the decomposition of exponentials:

    $$
    \mathbb{E}[\exp(\lambda W(t)) | \mathcal{F}(s)] = \mathbb{E}[\exp(\lambda(W(t) - W(s))) \cdot \exp(\lambda W(s)) | \mathcal{F}(s)]
    $$
    
    Since $\exp(\lambda W(s))$ is $\mathcal{F}(s)$-measurable, we factor it out as a constant:

    $$
    = \exp(\lambda W(s)) \cdot \mathbb{E}[\exp(\lambda(W(t) - W(s))) | \mathcal{F}(s)]
    $$
    
    Since the increment $W(t)-W(s)$ is independent of $\mathcal{F}(s)$ and follows an $N(0, t-s)$ distribution, using the MGF of the normal distribution $\mathbb{E}[e^{\lambda Z}] = e^{\lambda^2 \sigma^2 / 2}$:

    $$
    = \exp(\lambda W(s)) \cdot \exp\left( \frac{1}{2} \lambda^2 (t-s) \right)
    $$
    
    Substituting this expectation into the expression to be proven:

    $$
    \mathbb{E}\left[ \exp\left(\lambda W(t) - \frac{1}{2}\lambda^2 t\right) \Big| \mathcal{F}(s) \right] = \exp\left(-\frac{1}{2}\lambda^2 t\right) \cdot \exp(\lambda W(s)) \cdot \exp\left( \frac{1}{2} \lambda^2 (t-s) \right)
    $$
    
    Combining the exponents:
    
    $$
    -\frac{1}{2}\lambda^2 t + \frac{1}{2}\lambda^2 t - \frac{1}{2}\lambda^2 s = -\frac{1}{2}\lambda^2 s
    $$

    Finally yielding:

    $$
    = \exp\left(\lambda W(s) - \frac{1}{2}\lambda^2 s\right)
    $$

    The martingale property is proven.

---

### P70 Exercise 6: Second Moment of the Riemann Integral of a Brownian Path

**Problem:**
Let $X(t) = \int_0^t W(s)ds$. Prove that: $\mathbb{E}[X^2(t)] = \frac{t^3}{3}, \forall t > 0$.

??? success "Solution (Click to expand)"

    This is a highly representative stochastic process integral calculation. The core lies in using Fubini's Theorem to interchange the order of expectation and integration.

    First, write the squared term as a double Riemann integral:

    $$
    X^2(t) = \left( \int_0^t W(u)du \right) \left( \int_0^t W(v)dv \right) = \int_0^t \int_0^t W(u)W(v) du dv
    $$

    Take the expectation of both sides. By Fubini's Theorem (since the domain of integration is bounded and the function is absolutely integrable), we can pass the expectation inside the integrals:

    $$
    \mathbb{E}[X^2(t)] = \mathbb{E}\left[ \int_0^t \int_0^t W(u)W(v) du dv \right] = \int_0^t \int_0^t \mathbb{E}[W(u)W(v)] du dv
    $$

    Since the covariance function of Brownian motion is $\mathbb{E}[W(u)W(v)] = \min(u, v)$, substitute this into the equation:

    $$
    = \int_0^t \int_0^t \min(u, v) du dv
    $$

    Because $\min(u,v)$ is not differentiable along the diagonal $u=v$, we divide the integration region $[0,t] \times [0,t]$ along the diagonal into two parts: the region where $u \leqslant v$ and the region where $u > v$:

    $$
    = \int_0^t dv \int_0^v u du + \int_0^t du \int_0^u v dv
    $$

    Since these two regions of integration are perfectly symmetric, we can simply calculate one of them and multiply by 2:

    $$
    = 2 \int_0^t \left( \int_0^v u du \right) dv = 2 \int_0^t \frac{1}{2}v^2 dv = \int_0^t v^2 dv = \frac{1}{3}v^3 \Bigg|_0^t = \frac{t^3}{3}
    $$

    The conclusion is proven: $\mathbb{E}[X^2(t)] = \frac{t^3}{3}$.