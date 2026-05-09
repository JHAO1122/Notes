# 5.1 $L^p$ Spaces, Inequalities, and Completeness

After establishing a complete theory of Lebesgue integration, we will explore function spaces composed of integrable functions—the $L^p$ spaces. They are not only the crystallization of real analysis but also the core research objects of functional analysis. In this section, we will focus on several extremely important inequalities and their structural inclusion relations within $L^p$ spaces.

---

## 1. Chebyshev's Inequality

Before exploring $L^p$ spaces, we first introduce the most fundamental tail probability control inequality in measure theory—Chebyshev's inequality (also known as the Markov-Chebyshev inequality). It establishes a connection between the $L^p$ integral of a function and the measure of its level sets.

!!! info "Theorem 5.1.1 (Chebyshev's Inequality)"

    Let $f \in L^p(X, \mu)$ (where $0 < p < \infty$), and for any $\alpha > 0$, define the level set $E_\alpha = \{x \in X : |f(x)| > \alpha\}$. Then we have:

    \[
    \mu(\{x \in X : |f(x)| > \alpha\}) \le \left( \frac{\|f\|_p}{\alpha} \right)^p = \frac{1}{\alpha^p} \int_X |f|^p \, d\mu
    \]

??? proof "Proof of Chebyshev's Inequality"

    Using the characteristic function and the monotonicity of the integral, the proof is very direct.
    
    Clearly, on the set $E_\alpha$, we have $|f(x)| > \alpha$, and thus $|f(x)|^p > \alpha^p$.
    
    Hence, we can bound the integral from below:

    \[
    \int_X |f|^p \, d\mu \ge \int_{E_\alpha} |f|^p \, d\mu \ge \int_{E_\alpha} \alpha^p \, d\mu = \alpha^p \mu(E_\alpha)
    \]

    Dividing both sides by $\alpha^p$ yields the conclusion. This shows that if a function is in an $L^p$ space, the measure of the region where it takes large values must decay very rapidly.

---

## 2. Inclusion and Interpolation Relations of $L^p$ Spaces

There are deep connections between $L^p$ spaces corresponding to different indices $p$. Their inclusion relations largely depend on the finiteness of the underlying measure space.

### 2.1 Inclusion Relations under Finite Measure

If the measure space is finite, lower-order integrability is guaranteed by higher-order integrability.

!!! success "Theorem 5.1.2 (Inclusion of $L^p$ for Finite Measures)"

    If $\mu(A) < \infty$, and $0 < p < q \le \infty$, then:

    \[
    L^q(A) \subset L^p(A)
    \]

    and there exists a constant $C$ such that $\|f\|_p \le C \|f\|_q$.

    *(Hint: This can be proven directly through Hölder's inequality by treating the constant 1 in $\int |f|^p \cdot 1 \, d\mu$ as the conjugate function.)*

### 2.2 Algebraic Decomposition and Interpolation of Spaces

When the space has infinite measure (such as $\mathbb{R}^n$), there is no absolute inclusion relation between $L^p$ and $L^q$. However, under the condition $0 < p < q < r \le \infty$, the following elegant algebraic relationship exists.

!!! success "Theorem 5.1.3 (Decomposition of $L^q$ Functions)"

    For $0 < p < q < r \le \infty$, any function in the intermediate space $L^q$ can be decomposed into two parts belonging to $L^p$ and $L^r$, respectively, namely:

    \[
    L^q \subset L^p + L^r
    \]

??? proof "Proof of the Decomposition Theorem"

    Take any $f \in L^q$. We truncate the function using a height threshold of 1.
    
    Let $E = \{x \in X : |f(x)| \ge 1\}$. Define the decomposition $f = g + h$, where:

    * **Peak part**: $g = f \cdot \chi_E$

    * **Flat part**: $h = f \cdot (1 - \chi_E)$

    **Proving $g \in L^p$**:
    Since $|f| \ge 1$ on $E$ and $p < q$, we have $|f|^p \le |f|^q$.

    \[
    \int_X |g|^p \, d\mu = \int_{E} |f|^p \, d\mu \le \int_{E} |f|^q \, d\mu \le \int_X |f|^q \, d\mu < \infty
    \]

    Thus $g \in L^p$.

    **Proving $h \in L^r$**:
    Since $|f| < 1$ on $E^c$ and $q < r$, we have $|f|^r \le |f|^q$.

    \[
    \int_X |h|^r \, d\mu = \int_{E^c} |f|^r \, d\mu \le \int_{E^c} |f|^q \, d\mu \le \int_X |f|^q \, d\mu < \infty
    \]

    Thus $h \in L^r$. Hence we have shown $f = g + h \in L^p + L^r$.

!!! success "Theorem 5.1.4 (Interpolation Property of $L^p$ Spaces / Absolute Convexity)"

    For $0 < p < q < r \le \infty$, the intersection of the two extreme spaces must contain the intermediate space:

    \[
    L^p \cap L^r \subset L^q
    \]

    More quantitatively, there exists $\lambda \in (0, 1)$ satisfying $\frac{1}{q} = \frac{\lambda}{p} + \frac{1-\lambda}{r}$ (i.e., $\frac{\lambda q}{p} + \frac{(1-\lambda)q}{r} = 1$). For any $f \in L^p \cap L^r$, we can obtain the **Log-convexity bound** through Hölder's inequality:

    \[
    \|f\|_q \le \|f\|_p^\lambda \|f\|_r^{1-\lambda}
    \]

---

## 3. Integral Operators and Schur's Test

In functional analysis, we often study integral operators defined by kernel functions. Schur's test provides a very practical sufficient condition for such operators to be bounded on $L^p$ spaces.

!!! info "Definition 5.1.2 (Integral Operator)"

    Let $(X, \mathcal{M}, \mu)$ and $(Y, \mathcal{N}, \nu)$ be finite or $\sigma$-finite measure spaces. Let $k: X \times Y \to \mathbb{C}$ be a product measurable function. Define the integral operator $T$ as:

    \[
    Tf(x) = \int_Y k(x, y) f(y) \, d\nu(y)
    \]

!!! success "Theorem 5.1.5 (Schur's Test / Generalized Young's Inequality)"

    Suppose there exists a constant $C > 0$ such that the kernel function $k$ satisfies the following two conditions for consistent row/column integral bounds:

    * For almost every $y \in Y$:

    \[
    \int_X |k(x, y)| \, d\mu(x) \le C
    \]

    * For almost every $x \in X$:

    \[
    \int_Y |k(x, y)| \, d\nu(y) \le C
    \]

    Then for any $1 \le p \le \infty$, the operator $T$ is a bounded linear operator from $L^p(\nu)$ to $L^p(\mu)$, and the operator norm satisfies:

    \[
    \|Tf\|_{L^p(\mu)} \le C \|f\|_{L^p(\nu)}
    \]

??? proof "Core of the Proof for Schur's Test (Based on Hölder's Inequality)"

    We use a clever split of $|k(x, y)|$ (for the case $1 < p < \infty$, where the conjugate index $p'$ satisfies $1/p + 1/p' = 1$).

    Split $|k(x, y)|$ in the integral into $|k(x,y)|^{1/p'} \cdot |k(x,y)|^{1/p}$:

    \[
    |Tf(x)| \le \int_Y |k(x, y)|^{1/p'} \cdot \left(|k(x, y)|^{1/p} |f(y)|\right) d\nu(y)
    \]

    Apply Hölder's inequality with respect to the measure $\nu$ to these two terms:

    \[
    |Tf(x)|^p \le \left( \int_Y |k(x, y)| \, d\nu(y) \right)^{p/p'} \left( \int_Y |k(x, y)| |f(y)|^p \, d\nu(y) \right)
    \]

    Utilizing the second condition $\int_Y |k(x, y)| d\nu(y) \le C$, substitute and simplify (noting that $p/p' = p-1$):

    \[
    |Tf(x)|^p \le C^{p-1} \int_Y |k(x, y)| |f(y)|^p \, d\nu(y)
    \]

    Now integrate both sides with respect to $x$ over $X$, and use Tonelli's theorem to swap the order of integration:

    \[
    \int_X |Tf(x)|^p \, d\mu(x) \le C^{p-1} \int_Y |f(y)|^p \left( \int_X |k(x, y)| \, d\mu(x) \right) d\nu(y)
    \]

    Utilizing the first condition $\int_X |k(x, y)| d\mu(x) \le C$:

    \[
    \|Tf\|_{L^p(\mu)}^p \le C^{p-1} \cdot C \int_Y |f(y)|^p \, d\nu(y) = C^p \|f\|_{L^p(\nu)}^p
    \]

    Taking the $p$-th root of both sides proves $\|Tf\|_p \le C \|f\|_p$.

---

## 4. Generalized Minkowski Integral Inequality

The classical Minkowski inequality establishes that $L^p$ spaces satisfy the triangle inequality ($\|f + g\|_p \le \|f\|_p + \|g\|_p$), thus making them normed vector spaces.

In product measure spaces, this inequality can be generalized to an integral form over a continuous parameter, stating that "the $L^p$ norm of an integral is less than or equal to the integral of the $L^p$ norm."

!!! success "Theorem 5.1.6 (Generalized Minkowski Integral Inequality)"

    Let $f(x, y)$ be a non-negative product measurable function on $X \times Y$, and $1 \le p < \infty$. Then we have:

    \[
    \left[ \int_X \left( \int_Y f(x, y) \, d\nu(y) \right)^p d\mu(x) \right]^{1/p} \le \int_Y \left( \int_X f(x, y)^p \, d\mu(x) \right)^{1/p} d\nu(y)
    \]

    Expressed more intuitively in the language of operator norms: Suppose $y \mapsto \|f(\cdot, y)\|_{L^p(\mu)}$ is integrable on $Y$. Then the function $F(x) = \int_Y f(x, y) d\nu(y)$ is in $L^p(\mu)$ and satisfies:

    \[
    \left\| \int_Y f(\cdot, y) \, d\nu(y) \right\|_{L^p(\mu)} \le \int_Y \|f(\cdot, y)\|_{L^p(\mu)} \, d\nu(y)
    \]

This profound inequality is an indispensable tool when dealing with partial differential equations, the theory of operator semigroups, and norm estimates for convolution operators.