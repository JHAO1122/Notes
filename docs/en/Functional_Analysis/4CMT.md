# Chapter 4: Fixed Point Theorems and Their Applications

In the theories of differential equations, integral equations, and various other types of equations, subjects such as the existence and uniqueness of solutions, as well as the convergence of approximate solutions, are of paramount importance. To discuss these issues, we often transform the problem into finding a **fixed point** of a specific mapping.

!!! info "Definition 4.1 (Fixed Point)"

    Let $T$ be a mapping defined on a space $X$. If there exists a point $x \in X$ such that:

    $$
    T(x) = x
    $$

    then $x$ is called a **fixed point** of the mapping $T$.

---

## 1. The Fundamental Fixed Point Theorem: Contraction Mapping Principle

### 1.1 Contraction Mapping and the Fixed Point Theorem

!!! success "Theorem 4.1 (Contraction Mapping Principle)"

    Let $X$ be a **complete metric space** with metric $\rho$. Let $T$ be a mapping from $X$ into itself. If for any $x, y \in X$, the inequality:

    $$
    \rho(Tx, Ty) \le \theta \rho(x, y)
    $$

    holds, where $\theta$ is a constant satisfying $0 \le \theta < 1$.
    
    Then, $T$ has a **unique** fixed point in $X$, i.e., there exists a unique $\overline{x}$ such that $T\overline{x} = \overline{x}$. Furthermore, this fixed point can be obtained through the iteration method.

A mapping $T$ satisfying the above conditions is called a **contraction mapping**.

??? proof "Proof of the Contraction Mapping Principle"

    **1. Construct an iteration sequence and prove it is a Cauchy sequence:**

    Starting from any point $x_0 \in X$, we construct a sequence $\{x_n\}$ by iteration:

    $$
    x_1 = Tx_0, \quad x_2 = Tx_1, \quad \dots, \quad x_{n+1} = Tx_n, \quad \dots
    $$

    First, we estimate the distance between adjacent terms. From the definition of a contraction mapping:

    $$
    \rho(x_1, x_2) = \rho(Tx_0, Tx_1) \le \theta \rho(x_0, x_1)
    $$

    $$
    \rho(x_2, x_3) = \rho(Tx_1, Tx_2) \le \theta \rho(x_1, x_2) \le \theta^2 \rho(x_0, x_1)
    $$

    By mathematical induction, we obtain the general recurrence relation:

    $$
    \rho(x_{n-1}, x_n) \le \theta^{n-1} \rho(x_0, x_1)
    $$

    Consequently:

    $$
    \rho(x_n, x_{n+1}) = \rho(Tx_{n-1}, Tx_n) \le \theta \rho(x_{n-1}, x_n) \le \theta^n \rho(x_0, x_1)
    $$

    For any natural number $p \ge 1$, using the triangle inequality and the distance estimates above:

    $$
    \rho(x_n, x_{n+p}) \le \rho(x_n, x_{n+1}) + \rho(x_{n+1}, x_{n+2}) + \dots + \rho(x_{n+p-1}, x_{n+p})
    $$

    $$
    \le (\theta^n + \theta^{n+1} + \dots + \theta^{n+p-1}) \rho(x_0, x_1)
    $$

    Using the sum formula for a geometric series and noting that $\theta < 1$:

    $$
    \le \frac{\theta^n}{1 - \theta} \rho(x_0, x_1)
    $$

    Since $0 \le \theta < 1$, as $n \rightarrow \infty$, $\theta^n \rightarrow 0$. Thus, for any given $\epsilon > 0$, when $n$ is sufficiently large, $\rho(x_n, x_{n+p}) < \epsilon$ for all $p \ge 1$.
    
    This shows that $\{x_n\}$ is a **Cauchy sequence**. Since the space $X$ is **complete**, $\{x_n\}$ must converge to some point $\overline{x} \in X$.

    **2. Prove the limit point is a fixed point:**

    Since a contraction mapping satisfies the Lipschitz condition (with constant $\theta < 1$), it is necessarily a **continuous mapping**.
    
    Taking the limit as $n \rightarrow \infty$ on both sides of $x_{n+1} = Tx_n$, and utilizing continuity:

    $$
    \lim_{n \rightarrow \infty} x_{n+1} = \lim_{n \rightarrow \infty} Tx_n = T(\lim_{n \rightarrow \infty} x_n)
    $$

    Yielding:

    $$
    \overline{x} = T\overline{x}
    $$

    This proves that $\overline{x}$ is a fixed point of the mapping $T$.

    **3. Prove the uniqueness of the fixed point:**

    Suppose there exists another fixed point $\overline{y} \in X$ such that $\overline{y} = T\overline{y}$. We calculate the distance between $\overline{x}$ and $\overline{y}$:

    $$
    \rho(\overline{x}, \overline{y}) = \rho(T\overline{x}, T\overline{y}) \le \theta \rho(\overline{x}, \overline{y})
    $$

    Which implies $(1 - \theta)\rho(\overline{x}, \overline{y}) \le 0$.
    
    Given $\theta < 1$ and the non-negativity of the metric, it follows that $\rho(\overline{x}, \overline{y}) = 0$, hence $\overline{x} = \overline{y}$. The proof is complete. $\square$

### 1.2 Error Estimation of the Iterative Sequence

As seen in the proof, to obtain the fixed point, we can start from any point $x_0 \in X$ and establish the iterative sequence $x_n = T^n x_0$.

In practice, the exact fixed point is often difficult to find, so $x_n$ is used as an approximation. Based on the estimates in the proof, the error between $x_n$ and the true fixed point $\overline{x}$ can be estimated as (letting $p \rightarrow \infty$):

$$
\rho(x_n, \overline{x}) \le \frac{\theta^n}{1 - \theta} \rho(x_0, x_1)
$$

---

## 2. Applications of the Contraction Mapping Principle

### 2.1 Existence and Uniqueness of Solutions for Initial Value Problems

Consider the first-order ordinary differential equation initial value problem (Picard's Existence and Uniqueness Theorem):

$$
\frac{dy}{dx} = f(x, y), \quad y|_{x=x_0} = y_0
$$

Assume $f(x, y)$ is continuous on a certain domain and satisfies the **Lipschitz condition** with respect to $y$:

$$
|f(x, y) - f(x, y')| \le K|y - y'|
$$

where $K > 0$ is a constant.

From calculus, this initial value problem is equivalent to the following Volterra integral equation:

$$
y(x) = y_0 + \int_{x_0}^x f(t, y(t)) dt
$$

??? proof "Proof using the Contraction Mapping Principle (Click to expand)"

    Choose a sufficiently small $\delta > 0$ such that $K\delta < 1$.
    
    We define a mapping $T$ on the space of continuous functions $C[x_0 - \delta, x_0 + \delta]$:

    $$
    (Ty)(x) = y_0 + \int_{x_0}^x f(t, y(t)) dt, \quad x \in [x_0 - \delta, x_0 + \delta]
    $$

    For any two continuous functions $y_1(x), y_2(x) \in C[x_0 - \delta, x_0 + \delta]$, calculate the distance between their images under $T$:

    $$
    \rho(Ty_1, Ty_2) = \max_{x \in [x_0 - \delta, x_0 + \delta]} \left| \int_{x_0}^x [f(t, y_1(t)) - f(t, y_2(t))] dt \right|
    $$

    Using the Lipschitz condition:

    $$
    \le \max_{x} \left| \int_{x_0}^x K |y_1(t) - y_2(t)| dt \right|
    $$

    $$
    \le K \delta \max_{t \in [x_0 - \delta, x_0 + \delta]} |y_1(t) - y_2(t)| = K\delta \rho(y_1, y_2)
    $$

    Since we chose $K\delta < 1$, $T$ is a **contraction mapping**. By the Contraction Mapping Principle, there exists a unique continuous function $y_0(x)$ in $C[x_0 - \delta, x_0 + \delta]$ satisfying the integral equation, meaning the initial value problem has a unique solution in this interval.
    
    Subsequently, this solution can be extended to the entire real axis. $\square$

### 2.2 Existence and Uniqueness of Solutions for Linear Integral Equations

Consider the Fredholm linear integral equation of the second kind:

$$
x(t) = f(t) + \lambda \int_a^b K(t, s)x(s) ds
$$

where $f \in L^2[a, b]$ is a given function and $\lambda$ is a parameter. The kernel $K(t, s)$ is a measurable function defined on the square $a \le t \le b$, $a \le s \le b$, satisfying the square-integrable condition:

$$
\int_a^b \int_a^b |K(t, s)|^2 dt ds < \infty
$$

??? proof "Proof: Existence of a unique solution for sufficiently small $|\lambda|$"

    We define a mapping $T$ on the complete metric space $L^2[a, b]$:

    $$
    (Tx)(t) = f(t) + \lambda \int_a^b K(t, s)x(s) ds
    $$

    First, it can be shown via the Cauchy-Schwarz inequality that $Tx \in L^2[a, b]$, meaning $T$ maps $L^2[a, b]$ into itself:

    $$
    \int_a^b \left| \int_a^b K(t,s)x(s)ds \right|^2 dt \le \int_a^b \left[ \int_a^b |K(t,s)|^2 ds \int_a^b x(s)^2 ds \right] dt < \infty
    $$

    Next, we choose $|\lambda|$ to be sufficiently small such that:

    $$
    \theta = |\lambda| \left[ \int_a^b \int_a^b |K(t, s)|^2 ds dt \right]^{\frac{1}{2}} < 1
    $$

    Calculating the distance between images of any $y_1, y_2 \in L^2[a,b]$ (using Cauchy-Schwarz):

    $$
    \rho(Ty_1, Ty_2) = |\lambda| \left[ \int_a^b \left| \int_a^b K(t, s)(y_1(s) - y_2(s)) ds \right|^2 dt \right]^{\frac{1}{2}}
    $$

    $$
    \le |\lambda| \left[ \int_a^b \int_a^b |K(t, s)|^2 ds dt \right]^{\frac{1}{2}} \left[ \int_a^b |y_1(s) - y_2(s)|^2 ds \right]^{\frac{1}{2}}
    $$

    $$
    = \theta \rho(y_1, y_2)
    $$

    Since $\theta < 1$, $T$ is a contraction mapping. By the principle, the equation has a unique solution $y \in L^2[a, b]$ for a sufficiently small $|\lambda|$. $\square$

---

## 3. Generalized Contraction Mapping Principle

In some cases, the mapping $T$ itself may not satisfy the strict contraction condition, but one of its powers (iterations) does (as seen later with Volterra equations). In such instances, we use the generalized version.

Let $T^2x = T(Tx)$, and generally define $T^nx = T(T^{n-1}x)$ for any natural number $n$.

!!! success "Theorem 4.2 (Generalized Contraction Mapping Principle)"

    Let $X$ be a complete metric space with metric $\rho$. Let $T$ be a mapping from $X$ into itself. If there exists a natural number $n_0$ such that for any $x, y \in X$, the inequality:

    $$
    \rho(T^{n_0}x, T^{n_0}y) \le \theta \rho(x, y)
    $$

    holds, where $\theta$ satisfies $0 \le \theta < 1$.
    
    Then, $T$ has a **unique** fixed point in $X$.

??? proof "Proof of the Generalized Contraction Mapping Principle"

    **1. Existence:**
    
    By the given condition, $T^{n_0}$ satisfies the standard Contraction Mapping Principle. Since $X$ is complete, $T^{n_0}$ has a **unique fixed point** $\overline{x} \in X$:

    $$
    T^{n_0}\overline{x} = \overline{x}
    $$

    We need to show that $\overline{x}$ is also a fixed point of $T$. Apply $T$ to both sides of the equation:

    $$
    T(T^{n_0}\overline{x}) = T\overline{x}
    $$

    Using the associativity of operators:

    $$
    T^{n_0}(T\overline{x}) = T^{n_0+1}(\overline{x}) = T(T^{n_0}\overline{x}) = T\overline{x}
    $$

    This shows that the element $T\overline{x}$ is also a fixed point of $T^{n_0}$.
    
    However, the fixed point of $T^{n_0}$ is **unique**, so we must have:

    $$
    T\overline{x} = \overline{x}
    $$

    Thus, $\overline{x}$ is indeed a fixed point of $T$.

    **2. Uniqueness:**
    
    Suppose $T$ has another fixed point $\overline{y}$ such that $T\overline{y} = \overline{y}$. Then iterating $n_0$ times yields:

    $$
    T^{n_0}\overline{y} = T^{n_0-1}(T\overline{y}) = T^{n_0-1}\overline{y} = \dots = \overline{y}
    $$

    This implies $\overline{y}$ is also a fixed point of $T^{n_0}$. By the uniqueness of the fixed point for $T^{n_0}$, it follows that $\overline{y} = \overline{x}$. $\square$

---

## 4. Application: Existence of Solutions for Volterra Integral Equations

For Volterra integral equations, the upper limit of integration is a variable $t$ rather than a constant $b$.

!!! success "Theorem 4.3"

    Let the kernel $K(t, s)$ be a continuous function defined on the triangular region $a \le t \le b$, $a \le s \le t$.
    
    Then the Volterra integral equation:

    $$
    x(t) = f(t) + \lambda \int_a^t K(t, s)x(s) ds
    $$

    has a unique solution $x \in C[a, b]$ for any given $f \in C[a, b]$ and **any constant $\lambda \ne 0$**.

??? proof "Proof: Using the Generalized Contraction Mapping Principle (Click to expand)"

    We define mapping $T$ on the complete space $C[a, b]$:

    $$
    (Tx)(t) = f(t) + \lambda \int_a^t K(t, s)x(s) ds
    $$

    Let $M = \max_{a \le t \le b, a \le s \le t} |K(t, s)|$.

    For any $y_1, y_2 \in C[a, b]$, calculate the difference after one mapping:

    $$
    |(Ty_1)(t) - (Ty_2)(t)| = |\lambda| \left| \int_a^t K(t, s)(y_1(s) - y_2(s)) ds \right|
    $$

    $$
    \le |\lambda| M \int_a^t |y_1(s) - y_2(s)| ds
    $$

    $$
    \le |\lambda| M (t - a) \max_{s \in [a, b]} |y_1(s) - y_2(s)| = |\lambda| M (t - a) \rho(y_1, y_2)
    $$

    Note that the bound depends on $t$. If $\lambda$ or the interval is large, $|\lambda| M (t - a)$ might not be less than 1, so $T$ is not necessarily a contraction.
    
    However, we can iterate. By induction:

    $$
    |(T^n y_1)(t) - (T^n y_2)(t)| \le \frac{(|\lambda| M (t - a))^n}{n!} \rho(y_1, y_2)
    $$

    *Induction Step: Assume it holds for $n$. For $n+1$:*

    $$
    |(T^{n+1} y_1)(t) - (T^{n+1} y_2)(t)| = |\lambda| \left| \int_a^t K(t, s) (T^n y_1(s) - T^n y_2(s)) ds \right|
    $$

    $$
    \le |\lambda| M \int_a^t \frac{(|\lambda| M (s - a))^n}{n!} \rho(y_1, y_2) ds
    $$

    Using $\int_a^t (s - a)^n ds = \frac{(t - a)^{n+1}}{n+1}$:

    $$
    \le \frac{(|\lambda| M (t - a))^{n+1}}{(n+1)!} \rho(y_1, y_2)
    $$

    *Induction complete.*

    Since $(t-a) \le (b-a)$ for all $t$, we take the supremum:

    $$
    \rho(T^n y_1, T^n y_2) \le \frac{(|\lambda| M (b - a))^n}{n!} \rho(y_1, y_2)
    $$

    For any constant $C = |\lambda| M (b - a)$, as $n \rightarrow \infty$:

    $$
    \lim_{n \rightarrow \infty} \frac{C^n}{n!} = 0
    $$

    Thus, one can always find a sufficiently large $n_0$ such that:

    $$
    \theta = \frac{(|\lambda| M (b - a))^{n_0}}{n_0!} < 1
    $$

    This means $T^{n_0}$ is a **contraction mapping**. By the generalized principle, a unique solution exists for any $\lambda$. $\square$