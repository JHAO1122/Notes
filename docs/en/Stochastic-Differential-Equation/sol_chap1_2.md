---
tags:
  - Probability Theory
  - Stochastic Processes
  - Stochastic Differential Equations
  - Exercises
---

# 📝 Detailed Solutions to Homework Problems: Conditional Expectation and Brownian Motion

!!! abstract "About This Page"
    This page contains detailed solutions to key homework problems from Chapter 1 (Independence, Conditional Expectation) and Chapter 2 (Brownian Motion and Its Properties) of the *Stochastic Differential Equations* course. All solutions are presented in collapsible boxes; click to view the detailed derivations.

---

## Part I: Independence and Conditional Expectation

### Exercise 1

**Problem:**
Let the probability density function of the random variable $X$ be $f(x) = ax(1-x), x \in (0,1)$. For all other $x$, $f$ is zero.
(1) Find the constant $a$;
(2) Let $Y = X^3$, find the probability density function of $Y$.

??? success "Solution (click to expand)"
    
    **(1) Finding the constant $a$**
    
    By the normalization property of the probability density function, its integral over the entire space must equal 1:
    
    $$
    \int_{0}^{1} ax(1-x)dx = 1 
    $$

    Solving this integral gives:

    $$
    a \left( \frac{1}{2}x^2 - \frac{1}{3}x^3 \right) \Bigg|_{0}^{1} = 1 \implies a \left( \frac{1}{2} - \frac{1}{3} \right) = 1
    $$
    
    Solving yields: $a = 6$.

    <br>

    **(2) Finding the probability density function of $Y=X^3$**
    
    First, find the cumulative distribution function $F_Y(y)$ of $Y$. Since $X \in (0,1)$, we have $Y = X^3 \in (0,1)$.
    For $y \in (0,1)$:
    
    $$
    F_Y(y) = P(Y \le y) = P(X^3 \le y) = P(X \le y^{1/3})
    $$
    
    Substituting the probability density function of $X$ and computing the integral:

    $$
    F_Y(y) = \int_{0}^{y^{1/3}} 6x(1-x)dx = \left( 3x^2 - 2x^3 \right) \Big|_{0}^{y^{1/3}} = 3y^{2/3} - 2y
    $$
    
    Differentiating $F_Y(y)$ gives the probability density function $f_Y(y)$ of $Y$:
    
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
(1) $E[f(X,Y)|Y] = \int f(x,Y)dF_X(x)$;
(2) $P(X+Y \le x | Y) = F_X(x-Y)$.

??? success "Solution (click to expand)"

    **(1) Proof:**
    
    We use the "standard machine" approach from measure theory for a rigorous proof:
    
    **Step 1 (Indicator functions):** Let $f(x,y) = I_A(x)I_B(y)$, where $A, B$ are Borel sets.
    Since $X, Y$ are independent, using properties of conditional expectation:
    
    $$
    E[I_A(X)I_B(Y)|Y] = I_B(Y)E[I_A(X)|Y] = I_B(Y)E[I_A(X)] = I_B(Y)P(X \in A)
    $$
    
    On the other hand, for the integral form on the right-hand side:
    
    $$
    \int I_A(x)I_B(Y)dF_X(x) = I_B(Y) \int I_A(x)dF_X(x) = I_B(Y)P(X \in A)
    $$
    
    The two are equal, so the statement holds for indicator functions.
    
    **Step 2 (Simple functions):** By the linearity of expectation, the conclusion holds for simple functions.
    
    **Step 3 (Non-negative & Bounded functions):** By the Monotone Convergence Theorem and the decomposition $f = f^+ - f^-$, the original equation holds for all bounded continuous functions. Q.E.D.

    <br>

    **(2) Proof:**
    
    The conditional probability can be written as a conditional expectation:
    
    $$
    P(X+Y \le x | Y) = E[I_{\{X+Y \le x\}} | Y]
    $$
    
    Let $g(X, Y) = I_{\{X+Y \le x\}}$. Using the result proven in part (1):
    
    $$
    E[I_{\{X+Y \le x\}} | Y] = \int I_{\{u+Y \le x\}} dF_X(u) 
    $$

    Transform the integration domain as $u \le x - Y$:

    $$
    \int I_{\{u \le x-Y\}} dF_X(u) = \int_{-\infty}^{x-Y} dF_X(u) = F_X(x-Y)
    $$
    
    Q.E.D.

---

### Exercise 3

**Problem:**
Let $A, B$ be two events. Compute the conditional expectation $E[\chi_A|\chi_B]$.

??? success "Solution (click to expand)"

    Since $\chi_B$ can only take values 0 or 1, the $\sigma$-algebra $\sigma(\chi_B)$ consists of $\{\emptyset, \Omega, B, B^c\}$.
    The conditional expectation $E[\chi_A|\chi_B]$ must be $\sigma(\chi_B)$-measurable, so it takes the form:
    
    $$
    E[\chi_A|\chi_B] = c_1 \chi_B + c_2 \chi_{B^c}
    $$
    
    By the definition of conditional expectation:
    
    1. **On set $B$**: $\int_B E[\chi_A|\chi_B] dP = \int_B \chi_A dP = P(A \cap B)$. 
       Thus $c_1 P(B) = P(A \cap B) \implies c_1 = P(A|B)$.
    2. **On set $B^c$**: $\int_{B^c} E[\chi_A|\chi_B] dP = \int_{B^c} \chi_A dP = P(A \cap B^c)$. 
       Thus $c_2 P(B^c) = P(A \cap B^c) \implies c_2 = P(A|B^c)$.
    
    Therefore:
    
    $$
    E[\chi_A|\chi_B] = P(A|B)\chi_B + P(A|B^c)\chi_{B^c}
    $$

---

### Exercise 4

**Problem:** Let $\mathcal{V}_1, \mathcal{V}_2$ be two independent $\sigma$-fields. Prove: $E[E(X|\mathcal{V}_1)|\mathcal{V}_2] = E[X]$.

??? success "Solution (click to expand)"

    **Proof:** 1. $E(X|\mathcal{V}_1)$ is $\mathcal{V}_1$-measurable.
    2. Since $\mathcal{V}_1$ and $\mathcal{V}_2$ are independent, the random variable $E(X|\mathcal{V}_1)$ is independent of $\mathcal{V}_2$.
    3. For any random variable $Z$ independent of $\mathcal{G}$, we have $E[Z|\mathcal{G}] = E[Z]$.
    
    Applying this:
    
    $$
    E[E(X|\mathcal{V}_1)|\mathcal{V}_2] = E[E(X|\mathcal{V}_1)]
    $$

    By the Law of Total Expectation, $E[E(X|\mathcal{V}_1)] = E[X]$. Thus:
    
    $$
    E[E(X|\mathcal{V}_1)|\mathcal{V}_2] = E[X]
    $$
    $\square$

---

### Exercise 5

**Problem:**
If $E[|X_n - X|^p] \to 0$ as $n \to \infty$ for $p \ge 1$, prove:
$\lim_{n\to\infty} E[|E[X_n|\mathcal{V}] - E[X|\mathcal{V}]|^p] = 0$.

??? success "Solution (click to expand)"

    **Proof:**
    
    Using the linearity of conditional expectation and the **conditional Jensen's inequality**:
    
    $$
    |E[X_n|\mathcal{V}] - E[X|\mathcal{V}]|^p = |E[X_n - X | \mathcal{V}]|^p \le E[|X_n - X|^p | \mathcal{V}]
    $$
    
    Taking expectations on both sides:
    
    $$
    E[|E[X_n|\mathcal{V}] - E[X|\mathcal{V}]|^p] \le E[E[|X_n - X|^p | \mathcal{V}]] = E[|X_n - X|^p]
    $$
    
    Since $E[|X_n - X|^p] \to 0$, the left-hand side must also converge to 0. $\square$

---

### Exercise 6

**Problem:**
$\Omega = \{1, \dots, 8\}$, $P(\{i\}) = 1/10$ for $i \le 4$, $P(\{i\}) = 3/20$ for $i > 4$.
$X = \chi_{\{1,2,3,4\}} + 2\chi_{\{5,6,7,8\}}$, $Y = \chi_{\{1,5\}} + 2\chi_{\{2,3,4,6,7,8\}}$.
Compute: (1) $X E[Y]$; (2) $E[E[XY|\mathcal{V}]|\mathcal{H}]$.

??? success "Solution (click to expand)"

    **Calculation:**
    
    (1) $E[Y] = 1 \cdot P(\{1,5\}) + 2 \cdot P(\{2,3,4,6,7,8\}) = (0.1 + 0.15) + 2(0.3 + 0.45) = 0.25 + 1.5 = 1.75 = 7/4$.
    So $X E[Y] = \frac{7}{4} X$.

    (2) Since $\mathcal{H} \subset \mathcal{V}$, by the Tower Property: $E[E[XY|\mathcal{V}]|\mathcal{H}] = E[XY|\mathcal{H}]$.
    Since $X$ is $\mathcal{H}$-measurable, $E[XY|\mathcal{H}] = X E[Y|\mathcal{H}]$.
    On $\{1,2,3,4\}$, $E[Y|\mathcal{H}] = \frac{0.1 + 2(0.3)}{0.4} = \frac{0.7}{0.4} = 7/4$.
    On $\{5,6,7,8\}$, $E[Y|\mathcal{H}] = \frac{0.15 + 2(0.45)}{0.6} = \frac{1.05}{0.6} = 7/4$.
    Thus $E[E[XY|\mathcal{V}]|\mathcal{H}] = \frac{7}{4} X$.

---

### Exercise 7

**Problem:**
Prove $X(t) = E[X|\mathcal{F}(t)]$ is a martingale.

??? success "Solution (click to expand)"

    **Proof:**
    1. **Integrability**: $E|X(t)| = E|E[X|\mathcal{F}(t)]| \le E[E[|X||\mathcal{F}(t)]] = E|X| < \infty$.
    2. **Adaptability**: $X(t)$ is $\mathcal{F}(t)$-measurable by definition.
    3. **Martingale Property**: For $s < t$, $E[X(t)|\mathcal{F}(s)] = E[E[X|\mathcal{F}(t)]|\mathcal{F}(s)] = E[X|\mathcal{F}(s)] = X(s)$.
    $\square$

---

## Part II: Brownian Motion and Its Properties

### Exercise 1

**Problem:**
Prove $W(t+s)-W(s)$ and $cW(t/c^2)$ are Brownian motions.

??? success "Solution (click to expand)"

    (1) **For $B(t) = W(t+s) - W(s)$**:
    $B(0)=0$. Increments are independent and stationary because $W$ has independent stationary increments. $B(t)-B(u) \sim N(0, (t+s)-(u+s)) = N(0, t-u)$. Paths are continuous because $W$ is continuous.
    
    (2) **For $U(t) = cW(t/c^2)$**:
    $U(0)=0$. $Var(U(t)) = c^2 Var(W(t/c^2)) = c^2 (t/c^2) = t$. Increments are normal with mean 0. Independence follows from $W$. 

---

### Exercise 2

**Problem:**
Prove $tW(1/t) - sW(1/s) \sim N(0, t-s)$.

??? success "Solution (click to expand)"

    $Var(tW(1/t) - sW(1/s)) = t^2(1/t) + s^2(1/s) - 2ts E[W(1/t)W(1/s)]$.
    Since $1/t < 1/s$, $E[W(1/t)W(1/s)] = 1/t$.
    $Var = t + s - 2ts(1/t) = t + s - 2s = t - s$.
    Mean is clearly 0. Since it's a Gaussian process, the increment is $N(0, t-s)$.

---

### Exercise 3

**Problem:**
Prove $E[W^{2k}(t)] = \frac{(2k)!t^k}{2^k k!}$.

??? success "Solution (click to expand)"

    Use the MGF of $N(0, t)$: $M(\lambda) = \exp(\frac{1}{2}\lambda^2 t) = \sum \frac{(\frac{1}{2}\lambda^2 t)^k}{k!} = \sum \frac{t^k}{2^k k!} \lambda^{2k}$.
    Also $M(\lambda) = \sum E[W^m] \frac{\lambda^m}{m!}$.
    Equating coefficients of $\lambda^{2k}$: $\frac{E[W^{2k}]}{(2k)!} = \frac{t^k}{2^k k!} \implies E[W^{2k}] = \frac{(2k)!t^k}{2^k k!}$.

---

### Exercise 4

**Problem:**
Prove $E[\exp(c(W(s) - W(t)))] = \exp(\frac{1}{2}c^2(t-s))$.

??? success "Solution (click to expand)"

    $W(s) - W(t) \sim N(0, t-s)$.
    Let $Z \sim N(0, \sigma^2)$, then $E[e^{cZ}] = e^{\frac{1}{2}c^2\sigma^2}$.
    Here $\sigma^2 = t-s$, so $E = \exp(\frac{1}{2}c^2(t-s))$.

---

### Exercise 5

**Problem:**
$U(t) = e^{-t}W(e^{2t})$. Prove $E[U(t)U(s)] = e^{-|t-s|}$.

??? success "Solution (click to expand)"

    Assume $t \ge s$. $E[U(t)U(s)] = e^{-(t+s)} E[W(e^{2t})W(e^{2s})] = e^{-(t+s)} \min(e^{2t}, e^{2s}) = e^{-(t+s)} e^{2s} = e^{s-t} = e^{-|t-s|}$.

---

### Exercise 6

**Problem:**
Prove $\lim_{m \to \infty} \frac{W(m)}{m} = 0$ a.s.

??? success "Solution (click to expand)"

    $E[(W(m)/m)^4] = 3m^2/m^4 = 3/m^2$.
    $P(|W(m)/m| > \epsilon) \le \frac{E[(W(m)/m)^4]}{\epsilon^4} = \frac{3}{\epsilon^4 m^2}$.
    Since $\sum \frac{3}{\epsilon^4 m^2} < \infty$, by Borel-Cantelli lemma, $P(|W(m)/m| > \epsilon \text{ i.o.}) = 0$.
    Thus $\lim \frac{W(m)}{m} = 0$ a.s.

---

### Exercise 7

**Problem:**
Prove $W(t)^2 - t$ and $\exp(\lambda W_t - \frac{1}{2}\lambda^2 t)$ are martingales.

??? success "Solution (click to expand)"

    (1) $E[W(t)^2 - t | \mathcal{F}(s)] = E[(W(t)-W(s)+W(s))^2 - t | \mathcal{F}(s)] = E[(W(t)-W(s))^2] + W(s)^2 - t = (t-s) + W(s)^2 - t = W(s)^2 - s$.
    
    (2) $E[\exp(\lambda W_t - \frac{1}{2}\lambda^2 t) | \mathcal{F}(s)] = \exp(-\frac{1}{2}\lambda^2 t) \exp(\lambda W_s) E[\exp(\lambda(W_t-W_s))] = \exp(-\frac{1}{2}\lambda^2 t) \exp(\lambda W_s) \exp(\frac{1}{2}\lambda^2(t-s)) = \exp(\lambda W_s - \frac{1}{2}\lambda^2 s)$.

---

### Exercise 8

**Problem:**
$X(t) = \int_0^t W(s)ds$. Prove $E[X^2(t)] = t^3/3$.

??? success "Solution (click to expand)"

    $E[X^2(t)] = E[\int_0^t W(u)du \int_0^t W(v)dv] = \int_0^t \int_0^t E[W(u)W(v)] du dv$.
    $= \int_0^t \int_0^t \min(u, v) du dv = 2 \int_0^t (\int_0^v u du) dv = 2 \int_0^t \frac{1}{2}v^2 dv = \int_0^t v^2 dv = t^3/3$.