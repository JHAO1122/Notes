---
tags:
  - Stochastic Differential Equations
  - Homework Exercises
---

# 📝 Detailed Solutions to Homework Exercises: Conditional Expectation and Properties of Brownian Motion

!!! abstract "About This Page"
    This page contains detailed solutions to key homework exercises from Chapter 1 (Independence, Conditional Expectation) and Chapter 2 (Brownian Motion and Its Properties) of the *Stochastic Differential Equations* course. All solutions are presented in collapsible boxes; click to view the detailed derivation process.

---

## Part I: Independence and Conditional Expectation

### Exercise 1

**Problem:**
Let the probability density function of the random variable $X$ be $f(x) = ax(1-x), x \in (0,1)$. For all other $x$, $f$ is zero.
(1) Find the constant $a$;
(2) Let $Y = X^3$, find the probability density function of $Y$.

??? success "Solution (Click to expand)"
    
    **(1) Finding the constant $a$**
    
    By the normalization property of the probability density function, its integral over the entire space equals 1:
    
    $$
    \int_{0}^{1} ax(1-x)dx = 1 
    $$

    Solving this integral yields:

    $$
    a \left( \frac{1}{2}x^2 - \frac{1}{3}x^3 \right) \Bigg|_{0}^{1} = 1 \implies a \left( \frac{1}{2} - \frac{1}{3} \right) = 1
    $$
    
    Solving gives: $a = 6$ .

    <br>

    **(2) Finding the probability density function of $Y=X^3$**
    
    First, find the cumulative distribution function $F_Y(y)$ of $Y$. Since $X \in (0,1)$, we have $Y = X^3 \in (0,1)$.
    For $y \in (0,1)$:
    
    $$
    F_Y(y) = \mathbb{P}(Y \leqslant y) = \mathbb{P}(X^3 \leqslant y) = \mathbb{P}(X \leqslant y^{1/3})
    $$
    
    Substituting the probability density function of $X$ and computing the integral:

    $$
    F_Y(y) = \int_{0}^{y^{1/3}} 6x(1-x)dx = \left( 3x^2 - 2x^3 \right) \Big|_{0}^{y^{1/3}} = 3y^{2/3} - 2y
    $$
    
    Differentiating $F_Y(y)$ yields the probability density function $f_Y(y)$ of $Y$ :
    
    $$
    f_Y(y) = F_Y'(y) = 3 \cdot \frac{2}{3}y^{-1/3} - 2 = 2y^{-1/3} - 2
    $$
    
    In summary, the probability density function of $Y$ is:
    
    $$
    f_Y(y) = \begin{cases} 2y^{-1/3} - 2, & y \in (0,1) \\ 0, & \text{otherwise} \end{cases}
    $$

---

### Exercise 2

**Problem:**
Let $X, Y$ be independent random variables, $f(x,y)$ be a bounded continuous function, and $F_X$ be the probability distribution function of $X$. Prove:
(1) $\mathbb{E}[f(X,Y)|Y] = \int f(x,Y)dF_X(x)$;
(2) $\mathbb{P}(X+Y \leqslant x | Y) = F_X(x-Y)$.

??? success "Solution (click to expand)"

    **(1) Proof:**
    
    We adopt the rigorous proof method of the Standard Machine from measure theory :
    
    **Step 1 (Indicator functions):** Let $f(x,y) = I_A(x)I_B(y)$, where $A, B$ are Borel sets.
    Since $X, Y$ are independent, using properties of conditional expectation (given $Y$, $I_B(Y)$ can be factored out, and the conditional expectation of independent variables equals the unconditional expectation):
    
    $$
    \mathbb{E}[I_A(X)I_B(Y)|Y] = I_B(Y)\mathbb{E}[I_A(X)|Y] = I_B(Y)\mathbb{E}[I_A(X)] = I_B(Y)\mathbb{P}(X \in A)
    $$
    
    On the other hand, for the integral form on the right-hand side:
    
    $$
    \int I_A(x)I_B(Y)dF_X(x) = I_B(Y) \int I_A(x)dF_X(x) = I_B(Y)\mathbb{P}(X \in A)
    $$
    
    The two are equal, so the statement holds for indicator functions.
    
    **Step 2 (Simple functions):** By the linearity of expectation, the conclusion holds for simple functions, which are finite linear combinations of indicator functions.
    
    **Step 3 (Non-negative measurable functions):** For any non-negative measurable function, there exists a monotone increasing sequence of simple functions approximating it. By the Monotone Convergence Theorem, the conclusion holds for non-negative measurable functions.
    
    **Step 4 (Bounded continuous functions):** Any bounded continuous function can be decomposed into its positive and negative parts $f = f^+ - f^-$, and by definition it is Lebesgue integrable. Therefore, the original formula holds for all bounded continuous functions. QED.

    <br>

    **(2) Proof:**
    
    The conditional probability can be directly written in the form of conditional expectation:
    
    $$
    \mathbb{P}(X+Y \leqslant x | Y) = \mathbb{E}[I_{\{X+Y \leqslant x\}} | Y]
    $$
    
    Let $g(X, Y) = I_{\{X+Y \leqslant x\}}$. Using the result proven in part (1):
    
    $$
    \mathbb{E}[I_{\{X+Y \leqslant x\}} | Y] = \int I_{\{u+Y \leqslant x\}} dF_X(u)
    $$

    Transform the integration domain equivalently to $u \leqslant x - Y$:

    $$
    \int I_{\{u \leqslant x-Y\}} dF_X(u) = \int_{-\infty}^{x-Y} dF_X(u) = F_X(x-Y)
    $$
    
    QED.

---

### Exercise 3

**Problem:**
Let $A, B$ be two events in the measure space $(\Omega, \mathcal{F}, \mathbb{P})$. Compute the conditional expectation $\mathbb{E}[\chi_A|\chi_B]$. (Note: $\chi$ denotes the indicator function.)

??? success "Solution (click to expand)"

    Since $\chi_B$ can only take the values 0 or 1, the $\sigma$-algebra it generates is very simple: $\sigma(\chi_B) = \{\emptyset, \Omega, B, B^c\}$.
    Therefore, the conditional expectation $\mathbb{E}[\chi_A|\chi_B]$ must be $\sigma(\chi_B)$-measurable, which means it must be constant on $B$ and $B^c$. Let us set:
    
    $$
    \mathbb{E}[\chi_A|\chi_B] = c_1 \chi_B + c_2 \chi_{B^c}
    $$
    
    According to the Radon-Nikodym derivative definition of conditional expectation, for any $\Lambda \in \sigma(\chi_B)$, the integral equality must hold:
    
    $$
    \int_\Lambda \mathbb{E}[\chi_A|\chi_B] d\mathbb{P} = \int_\Lambda \chi_A d\mathbb{P}
    $$
    
    **Case 1: Take $\Lambda = B$**
    
    $$
    \int_B (c_1 \chi_B + c_2 \chi_{B^c}) d\mathbb{P} = \int_B \chi_A d\mathbb{P} \implies c_1 \mathbb{P}(B) = \mathbb{P}(A \cap B)
    $$
    
    If $\mathbb{P}(B) > 0$, then $c_1 = \frac{\mathbb{P}(A \cap B)}{\mathbb{P}(B)} = \mathbb{P}(A|B)$.
    
    **Case 2: Take $\Lambda = B^c$**
    
    $$
    \int_{B^c} (c_1 \chi_B + c_2 \chi_{B^c}) d\mathbb{P} = \int_{B^c} \chi_A d\mathbb{P} \implies c_2 \mathbb{P}(B^c) = \mathbb{P}(A \cap B^c)
    $$
    
    If $\mathbb{P}(B^c) > 0$, then $c_2 = \frac{\mathbb{P}(A \cap B^c)}{\mathbb{P}(B^c)} = \mathbb{P}(A|B^c)$.
    
    In summary, the explicit expression for the conditional expectation is:
    
    $$
    \mathbb{E}[\chi_A|\chi_B] = \mathbb{P}(A|B)\chi_B + \mathbb{P}(A|B^c)\chi_{B^c}
    $$

---

### Exercise 4

**Problem:**
Let $\mathcal{V}_1$ and $\mathcal{V}_2$ be two independent $\sigma$-fields, and let $X$ be an integrable random variable. Prove:
$\mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)|\mathcal{V}_2] = \mathbb{E}[X]$.

??? success "Solution (click to expand)"

    **Proof:**
    
    First, by the definition of conditional expectation, the inner expectation $\mathbb{E}(X|\mathcal{V}_1)$ must be a $\mathcal{V}_1$-measurable random variable.
    
    Since the problem states that $\mathcal{V}_1$ and $\mathcal{V}_2$ are independent $\sigma$-fields, the $\mathcal{V}_1$-measurable random variable $\mathbb{E}(X|\mathcal{V}_1)$ is naturally independent of the $\sigma$-field $\mathcal{V}_2$.
    
    According to the property of conditional expectation with respect to an independent $\sigma$-field (if a random variable $Z$ is independent of a $\sigma$-field $\mathcal{G}$, then $\mathbb{E}[Z|\mathcal{G}] = \mathbb{E}[Z]$), we can remove the conditioning on $\mathcal{V}_2$:
    
    $$
    \mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)|\mathcal{V}_2] = \mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)]
    $$
    
    Finally, by the law of total expectation (i.e., the degenerate form of the Tower Property):
    
    $$
    \mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)] = \mathbb{E}[X]
    $$
    
    Combining the two equations, the original proposition is proved. $\square$

---

### Exercise 5

**Problem:**
Let $X$ and $\{X_n\}$ be a sequence of random variables on $(\Omega, \mathcal{F}, \mathbb{P})$, and for a fixed $1 \leqslant p < \infty$, $\mathbb{E}[|X_n|^p] < +\infty$. Assume $\lim_{n\to\infty}\mathbb{E}[|X_n - X|^p] = 0$. Prove that for any sub-$\sigma$-field $\mathcal{V} \subset \mathcal{F}$ of $\mathcal{F}$, the following holds:

$$
\lim_{n\to\infty}\mathbb{E}\left[\big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p\right] = 0.
$$

??? success "Solution (click to expand)"

    **Proof:**

    Consider the function $\Phi(x) = |x|^p$. Since $p \geqslant 1$, this is a convex function.
    Using the linearity of conditional expectation and the **conditional Jensen's inequality**, we combine the conditional expectations inside the absolute value and move the absolute value inside the expectation:

    $$
    \big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p = \big|\mathbb{E}[X_n - X | \mathcal{V}]\big|^p \leqslant \mathbb{E}\big[|X_n - X|^p \big| \mathcal{V}\big]
    $$

    Taking the unconditional expectation on both sides of the above inequality and applying the law of total expectation:

    $$
    \mathbb{E}\left[\big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p\right] \leqslant \mathbb{E}\left[ \mathbb{E}\big[|X_n - X|^p \big| \mathcal{V}\big] \right] = \mathbb{E}\big[|X_n - X|^p\big]
    $$

    Taking the limit as $n \to \infty$ on both sides of the inequality. Since it is given that $\lim_{n\to\infty}\mathbb{E}[|X_n - X|^p] = 0$, and the left-hand side expectation is always non-negative, by the squeeze theorem:

    $$
    \lim_{n\to\infty}\mathbb{E}\left[\big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p\right] = 0
    $$

    The original proposition is proven. $\square$

---

### Exercise 6

**Problem:**
Let $\Omega = \{1, 2, \cdots, 7, 8\}$, $\mathcal{F} = 2^\Omega$ (i.e., the $\sigma$-algebra consisting of all subsets of $\Omega$). For $i \leqslant 4$, $\mathbb{P}(\{i\}) = 1/10$, and for $i > 4$, $\mathbb{P}(\{i\}) = 3/20$.
Define $X = \chi_{\{1,2,3,4\}} + 2\chi_{\{5,6,7,8\}}$, $Y = \chi_{\{1,5\}} + 2\chi_{\{2,3,4,6,7,8\}}$.
$\mathcal{V}$ is the $\sigma$-algebra generated by $\{1,2\}$ and $\{3,4\}$, $\mathcal{H}$ is the $\sigma$-algebra generated by $\{1,2,3,4\}$.
Compute: (1) $X\mathbb{E}[Y]$; (2) $\mathbb{E}[\mathbb{E}[XY|\mathcal{V}]|\mathcal{H}]$.

??? success "Solution (click to expand)"

    **Computational Preparation: Clarifying the Probability Structure**

    *   $\mathbb{P}(\{1,2,3,4\}) = 4 \times \frac{1}{10} = \frac{2}{5}$
    *   $\mathbb{P}(\{5,6,7,8\}) = 4 \times \frac{3}{20} = \frac{3}{5}$
    *   It is easy to verify that the two core atoms of $\mathcal{H} = \sigma(\{1,2,3,4\})$ are precisely $\{1,2,3,4\}$ and $\{5,6,7,8\}$. The partition of $\mathcal{V} = \sigma(\{1,2\}, \{3,4\})$ is finer, and it is clear that $\mathcal{H} \subset \mathcal{V}$.


    **(1) Computing $X\mathbb{E}[Y]$**

    First, compute the unconditional expectation of $Y$:

    $$
    \mathbb{E}[Y] = 1 \cdot \mathbb{P}(\{1,5\}) + 2 \cdot \mathbb{P}(\{2,3,4,6,7,8\})
    $$

    $$
    = 1 \cdot \left(\frac{1}{10} + \frac{3}{20}\right) + 2 \cdot \left(3 \times \frac{1}{10} + 3 \times \frac{3}{20}\right) = \frac{1}{4} + 2\left(\frac{3}{4}\right) = \frac{7}{4}
    $$

    Multiply this directly by $X$:

    $$
    X\mathbb{E}[Y] = \frac{7}{4}X = \frac{7}{4}\chi_{\{1,2,3,4\}} + \frac{14}{4}\chi_{\{5,6,7,8\}}
    $$


    **(2) Computing $\mathbb{E}[\mathbb{E}[XY|\mathcal{V}]|\mathcal{H}]$**

    Here we utilize the extremely elegant smoothing property (Tower Property) to directly reduce dimensions, avoiding the cumbersome calculation of double conditional expectations.
    Since $\mathcal{H} \subset \mathcal{V}$, the smaller $\sigma$-algebra plays the decisive role on the outside:

    $$
    \mathbb{E}\big[\mathbb{E}[XY|\mathcal{V}]\big|\mathcal{H}\big] = \mathbb{E}[XY|\mathcal{H}]
    $$

    Observing the random variable $X$, it is constantly 1 on $\{1,2,3,4\}$ and constantly 2 on $\{5,6,7,8\}$, which exactly corresponds to the atomic partition of $\mathcal{H}$. Therefore, $X$ is **$\mathcal{H}$-measurable**! As known information, it can be factored out of the conditional expectation:

    $$
    \mathbb{E}[XY|\mathcal{H}] = X \mathbb{E}[Y|\mathcal{H}]
    $$

    Next, we compute the conditional expectation of $Y$ on the two atoms of $\mathcal{H}$:
    
    * On $\{1,2,3,4\}$: $\mathbb{E}[Y|\{1,2,3,4\}] = \frac{1\cdot P(\{1\}) + 2\cdot P(\{2,3,4\})}{P(\{1,2,3,4\})} = \frac{1/10 + 6/10}{4/10} = \frac{7}{4}$
    
    * On $\{5,6,7,8\}$: $\mathbb{E}[Y|\{5,6,7,8\}] = \frac{1\cdot P(\{5\}) + 2\cdot P(\{6,7,8\})}{P(\{5,6,7,8\})} = \frac{3/20 + 18/20}{12/20} = \frac{7}{4}$

    It is surprisingly found that $\mathbb{E}[Y|\mathcal{H}] = 7/4$ holds constant on both partitions. Hence:

    $$
    \mathbb{E}[\mathbb{E}[XY|\mathcal{V}]|\mathcal{H}] = \frac{7}{4}X
    $$

---

### Exercise 7

**Problem:**
Given a probability space $(\Omega, \mathcal{F}, \mathbb{P})$ and an integrable random variable $X$, let $\{\mathcal{F}(t)\}_{t \geqslant 0}$ be a filtration. For $t \geqslant 0$, define $X(t) \doteq \mathbb{E}[X|\mathcal{F}(t)]$. Prove that $X(t)$ is a martingale with respect to $\mathcal{F}(t)$.

??? success "Solution (click to expand)"

    **Proof:**
    We need to strictly verify the three core conditions of a martingale:
    
    **1. Integrability:**
    Using the contraction property of conditional expectation (the absolute value special case of Jensen's inequality):
    
    $$
    \mathbb{E}|X(t)| = \mathbb{E}\big|\mathbb{E}[X|\mathcal{F}(t)]\big| \leqslant \mathbb{E}\big[\mathbb{E}[|X| \big| \mathcal{F}(t)]\big] = \mathbb{E}|X|
    $$
    
    Since the problem states that $X$ is integrable ($\mathbb{E}|X| < \infty$), we have $\mathbb{E}|X(t)| < \infty$.
    
    **2. Adaptability:**
    According to the measure-theoretic definition of conditional expectation $\mathbb{E}[X|\mathcal{F}(t)]$, it is naturally $\mathcal{F}(t)$-measurable. Therefore, the process $\{X(t)\}$ is adapted to the filtration $\{\mathcal{F}(t)\}_{t \ge 0}$.
    
    **3. Martingale Property:**
    For any $0 \leqslant s \leqslant t$, since $\{\mathcal{F}(t)\}$ is a filtration, meaning information accumulates, we must have $\mathcal{F}(s) \subset \mathcal{F}(t)$.
    Using the smoothing property (Tower Property) of conditional expectation:
    
    $$
    \mathbb{E}[X(t) | \mathcal{F}(s)] = \mathbb{E}\big[\mathbb{E}[X|\mathcal{F}(t)] \big| \mathcal{F}(s)\big] = \mathbb{E}[X|\mathcal{F}(s)] = X(s)
    $$
    
    In summary, based on these three points, the process $X(t)$ constructed by projecting a single integrable variable onto different information filtrations must be a martingale (this is known as a Doob martingale in stochastic analysis). $\square$

---
<br>

## Part II: Brownian Motion and Its Properties

### Exercise 1

**Problem:**
Let $W(t)$ be a one-dimensional Brownian motion. Prove that for any fixed $s>0$, $W(t+s)-W(s)$ is a Brownian motion; and for any positive constant $c$, $cW(t/c^2)$ is also a Brownian motion.

??? success "Solution (click to expand)"

    **(1) Prove that $B(t) = W(t+s) - W(s)$ is a Brownian motion**

    We verify the four defining properties of Brownian motion one by one:

    1. **Zero initial value**: $B(0) = W(0+s) - W(s) = 0$.

    2. **Independent increments**: For any $0 \le t_1 < t_2 < \cdots < t_k$, the sequence of increments is:

    $$B(t_j) - B(t_{j-1}) = W(t_j+s) - W(t_{j-1}+s)$$

    Since the original process $W$ has the independent increments property, and the time intervals $[t_{j-1}+s, t_j+s]$ are non-overlapping on the time axis, these increments are mutually independent.

    3. **Stationary normal increments**: For $t > u$, the increment distribution is:

    $$B(t) - B(u) = W(t+s) - W(u+s) \sim N(0, (t+s) - (u+s)) = N(0, t-u)$$

    4. **Path continuity**: Because the sample paths of $W(\cdot)$ are almost surely continuous, the shifted paths $B(\cdot)$ are also almost surely continuous.

    In conclusion, $B(t)$ is a one-dimensional Brownian motion.


    **(2) Prove that $U(t) = cW(t/c^2)$ is a Brownian motion**

    Similarly, we verify the properties:

    1. **Zero initial value**: $U(0) = cW(0/c^2) = cW(0) = 0$.

    2. **Independent increments**: For any $0 \le t_1 < t_2 < \cdots < t_k$, since $t_1/c^2 < t_2/c^2 < \cdots < t_k/c^2$ are non-overlapping, the increments of the original Brownian motion are independent. Multiplying by the constant $c$ preserves independence.

    3. **Stationary normal increments**: Due to the properties of linear transformation, the mean remains $0$, and the variance is:

    $$Var(U(t) - U(s)) = Var\left( c(W(t/c^2) - W(s/c^2)) \right) = c^2 \cdot \left( \frac{t-s}{c^2} \right) = t-s$$
   
    Hence, $U(t) - U(s) \sim N(0, t-s)$.

    4. **Path continuity**: Scaling does not affect the continuity of sample paths.

    In conclusion, $U(t) = cW(t/c^2)$ is also a Brownian motion (this is known as the scaling invariance/fractal property of Brownian motion).

---

### Exercise 2

**Problem:**
Let $W(t)$ be a one-dimensional Brownian motion, and define $\tilde{W}(t) = \begin{cases} tW(1/t), & t > 0, \\ 0, & t = 0. \end{cases}$
Prove: $\tilde{W}(t) - \tilde{W}(s) \sim N(0, t-s), \forall 0 < s < t$.

??? success "Solution (click to expand)"

    For $0 < s < t$, we perform an identity transformation on the increment to arrange mutually independent terms:

    $$
    \tilde{W}(t) - \tilde{W}(s) = tW\left(\frac{1}{t}\right) - sW\left(\frac{1}{s}\right)
    $$

    We separate by adding and subtracting $sW(1/t)$:

    $$
    = (t-s)W\left(\frac{1}{t}\right) - s\left(W\left(\frac{1}{s}\right) - W\left(\frac{1}{t}\right)\right)
    $$

    Note that here $1/t < 1/s$. Using the independent increments property of Brownian motion, the random variable $W(1/t) = W(1/t) - W(0)$ is independent of the increment $W(1/s) - W(1/t)$.
    
    Due to independence, the variances can be directly added:

    $$
    Var(\tilde{W}(t) - \tilde{W}(s)) = (t-s)^2 Var\left( W\left(\frac{1}{t}\right) \right) + s^2 Var\left( W\left(\frac{1}{s}\right) - W\left(\frac{1}{t}\right) \right)
    $$

    Substituting the variance formula $Var(W(u)) = u$:

    $$
    = (t-s)^2 \left(\frac{1}{t}\right) + s^2 \left(\frac{1}{s} - \frac{1}{t}\right)
    $$

    Simplifying by common denominator:

    $$
    = \frac{t^2 - 2ts + s^2}{t} + s - \frac{s^2}{t} = \frac{t^2 - 2ts}{t} + s = t - 2s + s = t - s
    $$

    Furthermore, by the linearity of expectation, $E[\tilde{W}(t) - \tilde{W}(s)] = 0 - 0 = 0$.
    Since it is a linear combination of jointly normal random variables, it must follow a normal distribution.
    
    The conclusion is proven: $\tilde{W}(t) - \tilde{W}(s) \sim N(0, t-s)$.

---

### Exercise 3

**Problem:**
Let $W(t)$ be a Brownian motion. Prove: $\mathbb{E}[W^{2k}(t)] = \frac{(2k)!t^k}{2^k k!}, \forall t > 0$.

??? success "Solution (click to expand)"

    Since $W(t) \sim N(0, t)$, using the moment generating function (MGF) of the normal distribution $M_W(\lambda) = \mathbb{E}[e^{\lambda W(t)}]$:

    $$
    M_W(\lambda) = \exp\left( \frac{1}{2} \lambda^2 t \right)
    $$

    We expand both sides of the equation as Taylor series at $\lambda = 0$:

    The left-hand side is expanded using the linearity of expectation:

    $$
    \mathbb{E}[e^{\lambda W(t)}] = \sum_{m=0}^{\infty} \mathbb{E}[W^m(t)] \frac{\lambda^m}{m!}
    $$

    The right-hand side expands the exponential function:

    $$
    \exp\left( \frac{1}{2} \lambda^2 t \right) = \sum_{k=0}^{\infty} \frac{1}{k!} \left( \frac{1}{2} \lambda^2 t \right)^k = \sum_{k=0}^{\infty} \frac{t^k}{2^k k!} \lambda^{2k}
    $$

    Comparing the coefficients of $\lambda^{2k}$ on both sides:

    $$
    \frac{\mathbb{E}[W^{2k}(t)]}{(2k)!} = \frac{t^k}{2^k k!}
    $$

    Rearranging yields the proof:

    $$
    \mathbb{E}[W^{2k}(t)] = \frac{(2k)! t^k}{2^k k!}
    $$

    *(Note: For odd-order moments, since the right-hand side has no odd-power terms, comparing coefficients gives $\mathbb{E}[W^{2k+1}(t)] = 0$.)*

---

### Exercise 4

**Problem:**
Prove: Let $c$ be a constant, $0 < s < t$, then $\mathbb{E}[\exp(c(W(s) - W(t)))] = \exp\left(\frac{1}{2}c^2(t-s)\right)$.

??? success "Solution (click to expand)"

    Due to the stationarity and symmetry of Brownian motion increments, $W(s) - W(t)$ and $W(t) - W(s)$ have the same distribution, both following $N(0, t-s)$.
    
    This is equivalent to finding the moment generating function of a normal random variable.
    Let $Z = c(W(s) - W(t))$, then $Z \sim N(0, c^2(t-s))$.
    
    According to the moment generating function formula for a normal distribution $N(\mu, \sigma^2)$, $\mathbb{E}[e^Z] = \exp(\mu + \frac{1}{2}\sigma^2)$:
    
    $$
    \mathbb{E}[\exp(c(W(s) - W(t)))] = \exp\left( 0 + \frac{1}{2} \cdot c^2(t-s) \right) = \exp\left( \frac{1}{2}c^2(t-s) \right)
    $$
    
    Q.E.D.

---

### Exercise 5

**Problem:**
Let $U(t) = e^{-t}W(e^{2t})$, then $\mathbb{E}[U(t)U(s)] = e^{-|t-s|}, \forall t, s \in \mathbb{R}$.

??? success "Solution (click to expand)"

    Without loss of generality, we assume $t \geqslant s$.

    Substituting the definition of $U(t)$ to compute the expectation of the cross term:

    $$
    \mathbb{E}[U(t)U(s)] = \mathbb{E}\left[ e^{-t}W(e^{2t}) \cdot e^{-s}W(e^{2s}) \right] = e^{-(t+s)} \mathbb{E}\left[ W(e^{2t})W(e^{2s}) \right]
    $$

    For Brownian motion, we know its covariance function is $\mathbb{E}[W(u)W(v)] = \min(u, v)$.
    Since we assumed $t \geqslant s$, naturally $e^{2t} \geqslant e^{2s}$, therefore:

    $$
    \mathbb{E}\left[ W(e^{2t})W(e^{2s}) \right] = \min(e^{2t}, e^{2s}) = e^{2s}
    $$

    Substituting this back into the original expression:

    $$
    \mathbb{E}[U(t)U(s)] = e^{-(t+s)} \cdot e^{2s} = e^{s-t}
    $$

    Since $t \geqslant s$, $s-t = -|t-s|$. If $s > t$, the conclusion holds symmetrically.
    Thus $\forall t, s \in \mathbb{R}$, $\mathbb{E}[U(t)U(s)] = e^{-|t-s|}$, which completes the proof.

---

### Exercise 6

**Problem:**
Prove that almost surely $\lim_{m \to \infty} \frac{W(m)}{m} = 0$.

??? success "Solution (click to expand)"

    We can use Chebyshev's Inequality combined with the Borel-Cantelli lemma to give a rigorous measure-theoretic proof.
    
    **Step 1: Bounding the variance**
    Consider the random sequence $Y_m = \frac{W(m)}{m}$. Since $W(m) \sim N(0, m)$, we have:

    $$
    Var\left(\frac{W(m)}{m}\right) = \frac{1}{m^2} Var(W(m)) = \frac{m}{m^2} = \frac{1}{m}
    $$

    **Step 2: Strengthening the Moment Bound (Paving the Way for Borel-Cantelli)**

    If we only use Chebyshev's inequality with the second moment, we get $P(|Y_m| \ge \epsilon) \le \frac{1}{m\epsilon^2}$. However, $\sum_{m=1}^{\infty} \frac{1}{m}$ is the harmonic series, which does not converge, preventing the direct use of the B-C lemma.
    Therefore, we use the fourth moment (or exploit the exponential decay property of the normal distribution). For a standard normal variable $Z \sim N(0,1)$, $\mathbb{E}[Z^4] = 3$.
    
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
    
    By the first part of the Borel-Cantelli lemma, the probability that the event $\left\{ \left|\frac{W(m)}{m}\right| \geqslant \epsilon \right\}$ occurs infinitely often is 0.
    This is equivalent to: almost every sample path satisfies $\lim_{m \to \infty} \frac{W(m)}{m} = 0$, completing the proof.

---

### Exercise 7

**Problem:**
Prove that $W(t)^2 - t$ and $\exp\left(\lambda W_t - \frac{1}{2}\lambda^2 t\right)$ $(\lambda \in \mathbb{R})$ are both martingales with respect to the history $\mathcal{F}(t)$ of $W(t)$.

??? success "Solution (click to expand)"

    **Proof (1): $W(t)^2 - t$ is a Martingale**
    
    Let $0 \leqslant s < t$. We need to prove $\mathbb{E}[W(t)^2 - t | \mathcal{F}(s)] = W(s)^2 - s$.
    Using the identity transformation $W(t) = (W(t) - W(s)) + W(s)$ to expand the square term:

    $$
    \mathbb{E}[W(t)^2 | \mathcal{F}(s)] = \mathbb{E}[((W(t) - W(s)) + W(s))^2 | \mathcal{F}(s)]
    $$

    $$
    = \mathbb{E}[(W(t) - W(s))^2 | \mathcal{F}(s)] + 2\mathbb{E}[W(s)(W(t) - W(s)) | \mathcal{F}(s)] + \mathbb{E}[W(s)^2 | \mathcal{F}(s)]
    $$

    * First term: Due to the independent increments of Brownian motion, $W(t)-W(s)$ is independent of $\mathcal{F}(s)$, so the conditional expectation equals the unconditional expectation $\mathbb{E}[(W(t)-W(s))^2] = t-s$.
    * Second term: $W(s)$ is $\mathcal{F}(s)$-measurable and can be factored out. The remaining increment expectation is 0.
    * Third term: $W(s)^2$ is $\mathcal{F}(s)$-measurable, so it equals itself.

    Therefore:

    $$
    \mathbb{E}[W(t)^2 | \mathcal{F}(s)] = (t-s) + 0 + W(s)^2
    $$

    Subtract $t$:

    $$
    \mathbb{E}[W(t)^2 - t | \mathcal{F}(s)] = (t-s) + W(s)^2 - t = W(s)^2 - s
    $$

    The martingale property is proved.

    <br>

    **Proof (2): $\exp\left(\lambda W(t) - \frac{1}{2}\lambda^2 t\right)$ is a martingale (geometric Brownian martingale)**

    Again, let $0 \leqslant s < t$. Using the exponential split:

    $$
    \mathbb{E}[\exp(\lambda W(t)) | \mathcal{F}(s)] = \mathbb{E}[\exp(\lambda(W(t) - W(s))) \cdot \exp(\lambda W(s)) | \mathcal{F}(s)]
    $$

    Since $\exp(\lambda W(s))$ is $\mathcal{F}(s)$-measurable, factor it out as a constant:

    $$
    = \exp(\lambda W(s)) \cdot \mathbb{E}[\exp(\lambda(W(t) - W(s))) | \mathcal{F}(s)]
    $$

    Because the increment $W(t)-W(s)$ is independent of $\mathcal{F}(s)$ and follows a $N(0, t-s)$ distribution, using the moment generating function of the normal distribution $\mathbb{E}[e^{\lambda Z}] = e^{\lambda^2 \sigma^2 / 2}$:

    $$
    = \exp(\lambda W(s)) \cdot \exp\left( \frac{1}{2} \lambda^2 (t-s) \right)
    $$

    Substituting the above expectation into the expression to be proved:

    $$
    \mathbb{E}\left[ \exp\left(\lambda W(t) - \frac{1}{2}\lambda^2 t\right) \Big| \mathcal{F}(s) \right] = \exp\left(-\frac{1}{2}\lambda^2 t\right) \cdot \exp(\lambda W(s)) \cdot \exp\left( \frac{1}{2} \lambda^2 (t-s) \right)
    $$
    
    Combine the exponents:
    
    $$
    -\frac{1}{2}\lambda^2 t + \frac{1}{2}\lambda^2 t - \frac{1}{2}\lambda^2 s = -\frac{1}{2}\lambda^2 s
    $$

    Finally, we obtain:

    $$
    = \exp\left(\lambda W(s) - \frac{1}{2}\lambda^2 s\right)
    $$

    The martingale property is proven.

---

### Exercise 8

**Problem:**
Set $X(t) = \int_0^t W(s)ds$. Prove: $\mathbb{E}[X^2(t)] = \frac{t^3}{3}, \forall t > 0$.

??? success "Solution (click to expand)"

    This is a highly representative stochastic process integral calculation, with the core being the use of Fubini's theorem to interchange the order of expectation and integration.

    First, write the squared term as a double Riemann integral:

    $$
    X^2(t) = \left( \int_0^t W(u)du \right) \left( \int_0^t W(v)dv \right) = \int_0^t \int_0^t W(u)W(v) du dv
    $$

    Take the expectation on both sides, and use Fubini's theorem (since the integration domain is bounded and the function is absolutely integrable) to move the expectation inside the integrals:

    $$
    \mathbb{E}[X^2(t)] = \mathbb{E}\left[ \int_0^t \int_0^t W(u)W(v) du dv \right] = \int_0^t \int_0^t \mathbb{E}[W(u)W(v)] du dv
    $$

    Since the covariance function of Brownian motion is $\mathbb{E}[W(u)W(v)] = \min(u, v)$, substitute it into the above equation:

    $$
    = \int_0^t \int_0^t \min(u, v) du dv
    $$

    Because $\min(u,v)$ is not differentiable on the diagonal $u=v$, we partition the integration region $[0,t] \times [0,t]$ along the diagonal into two parts: the region $u \leqslant v$ and the region $u > v$:

    $$
    = \int_0^t dv \int_0^v u du + \int_0^t du \int_0^u v dv
    $$

    Since these two regions are completely symmetric, we can compute one and multiply by 2:

    $$
    = 2 \int_0^t \left( \int_0^v u du \right) dv = 2 \int_0^t \frac{1}{2}v^2 dv = \int_0^t v^2 dv = \frac{1}{3}v^3 \Bigg|_0^t = \frac{t^3}{3}
    $$

    The conclusion is proven: $\mathbb{E}[X^2(t)] = \frac{t^3}{3}$.