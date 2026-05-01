# Chapter 7: Fundamentals of Inner Product Spaces and Hilbert Spaces

In normed linear spaces, we introduced the concepts of "length" (norm) and "distance". However, functional analysis still lacks a very important concept from analytic geometry—**angle (orthogonality)**. To extend concepts like orthogonality and orthogonal projection to more general abstract spaces, we first need to generalize the concept of "inner product" and induce a norm through it.

---

## 1. Inner Product and Inner Product Spaces

!!! info "Definition 7.1 (Inner Product and Inner Product Space)"

    Let $X$ be a linear space over the real (or complex) number field $K$. Then any pair of elements $x, y$ in $X$ constantly corresponds to a number in $K$, denoted as $(x, y)$, satisfying the following conditions:

    * (i) Linearity with respect to the first argument: $(\alpha x, y) = \alpha (x, y)$;

    * (ii) Additivity with respect to the first argument: $(x+y, z) = (x, z) + (y, z)$, where $z \in X$;

    * (iii) Conjugate symmetry: When $K$ is the real number field, $(x, y) = (y, x)$; when $K$ is the complex number field, $(x, y) = \overline{(y, x)}$;

    * (iv) Non-negativity: $(x, x) \ge 0$, and $(x, x) = 0$ if and only if $x = \theta$.

    The space $X$ satisfying the above conditions is called a **real (or complex) inner product space**, simply referred to as an inner product space. $(x, y)$ is called the inner product of elements $x$ and $y$.

The following basic facts can be directly derived from the definition of the inner product:

* When $K$ is the real number field, $(x, y)$ is linear with respect to both $x$ and $y$:

$$
(\alpha_1 x_1 + \alpha_2 x_2, y) = \alpha_1(x_1, y) + \alpha_2(x_2, y)
$$

$$
(x, \alpha_1 y_1 + \alpha_2 y_2) = \alpha_1(x, y_1) + \alpha_2(x, y_2)
$$

* When $K$ is the complex number field, $(x, y)$ is linear with respect to $x$, but **antilinear (or conjugate linear)** with respect to the second argument $y$:

$$
(x, \alpha_1 y_1 + \alpha_2 y_2) = \overline{\alpha_1}(x, y_1) + \overline{\alpha_2}(x, y_2)
$$

* When either $x$ or $y$ is equal to the zero element $\theta$, the inner product is $0$. For example, let $y = \theta$, then:

$$
(x, \theta) = (x, 0\theta) = 0(x, \theta) = 0
$$

---

## 2. Schwarz Inequality and Norm Derivation

In an inner product space, for any $x \in X$, we define the norm $\|x\| = (x, x)^{\frac{1}{2}}$. To prove that the function defined this way truly satisfies the triangle inequality of a norm, we first need to prove the famous Schwarz inequality.

!!! success "Theorem 7.1 (Schwarz Inequality)"

    In an inner product space, for any $x, y \in X$, we have:

    $$
    |(x, y)| \le \|x\| \|y\|
    $$

??? proof "Proof of the Schwarz Inequality"

    Let $X$ be a complex inner product space. When $y = \theta$, both sides are 0, and the inequality obviously holds.

    Now suppose $y \neq \theta$. For any complex number $\lambda$, by the non-negativity of the inner product, we have:

    $$
    (x + \lambda y, x + \lambda y) \ge 0
    $$

    Expanding this expression yields:

    $$
    0 \le (x, x) + \overline{\lambda}(x, y) + \lambda(y, x) + |\lambda|^2 (y, y)
    $$

    $$
    = (x, x) + \lambda(y, x) + \overline{\lambda} [ (x, y) + \lambda(y, y) ]
    $$

    We cleverly choose $\lambda = -\frac{(x, y)}{(y, y)}$ so that the terms inside the square brackets become zero. Substituting this into the above equation gives:

    $$
    0 \le (x, x) - \frac{(x, y)(y, x)}{(y, y)}
    $$

    Which is $0 \le (x, x) - \frac{|(x, y)|^2}{(y, y)}$. Rearranging the terms, we get:

    $$
    |(x, y)|^2 \le (x, x)(y, y) = \|x\|^2 \|y\|^2
    $$

    Taking the square root proves the Schwarz inequality. $\square$

Using the Schwarz inequality, we can easily prove the triangle inequality: for any $x, y \in X$,

$$
\|x + y\|^2 = (x + y, x + y) = (x, x) + (x, y) + (y, x) + (y, y)
$$

$$
\le \|x\|^2 + 2|(x, y)| + \|y\|^2 \le \|x\|^2 + 2\|x\|\|y\| + \|y\|^2 = (\|x\| + \|y\|)^2
$$

Taking the square root of both sides gives $\|x + y\| \le \|x\| + \|y\|$. Therefore, under this norm, $X$ is a normed linear space. We say that $\|\cdot\|$ is the norm induced by the inner product of $X$.

---

## 3. Continuity of the Inner Product and Polarization Identity

* **Continuity of the Inner Product**: The inner product $(x, y)$ is a continuous function of two variables with respect to $x$ and $y$.

??? proof "Proof: The Inner Product is a Continuous Function"

    Let $x_n, y_n$ be sequences in $X$ that converge in norm to $x, y \in X$, respectively. By the triangle inequality and the Schwarz inequality:

    $$
    |(x_n, y_n) - (x, y)| \le |(x_n, y_n) - (x, y_n)| + |(x, y_n) - (x, y)|
    $$

    $$
    = |(x_n - x, y_n)| + |(x, y_n - y)| \le \|x_n - x\| \|y_n\| + \|x\| \|y_n - y\|
    $$

    Since $y_n$ converges in norm, the sequence $\{\|y_n\|\}$ is bounded. As $n \to \infty$, $\|x_n - x\| \to 0$ and $\|y_n - y\| \to 0$, so the entire expression approaches 0. $\square$

* **Polarization Identity**: There is a polarization relationship between the inner product and the induced norm.

    When $K$ is the real number field:

    $$
    (x, y) = \frac{1}{4} (\|x + y\|^2 - \|x - y\|^2)
    $$

    When $K$ is the complex number field:

    $$
    (x, y) = \frac{1}{4} (\|x + y\|^2 - \|x - y\|^2 + i\|x + iy\|^2 - i\|x - iy\|^2)
    $$

??? proof "Expansion of the Proof for the Polarization Identity"

    First, calculate the real part:

    $$
    \|x + y\|^2 - \|x - y\|^2 = (x + y, x + y) - (x - y, x - y)
    $$

    $$
    = [(x, x) + (x, y) + (y, x) + (y, y)] - [(x, x) - (x, y) - (y, x) + (y, y)]
    $$

    $$
    = 2(x, y) + 2(y, x) = 4 \text{Re}(x, y)
    $$

    For a complex inner product space, we also need the imaginary part. Note that $\text{Re}(x, iy) = \text{Re}\{ \overline{i}(x, y) \} = \text{Re}\{ -i(x,y) \} = \text{Im}(x, y)$, so:

    $$
    i(\|x + iy\|^2 - \|x - iy\|^2) = i(4\text{Re}(x, iy)) = 4i\text{Im}(x, y)
    $$

    Adding the real part and the imaginary part and dividing by 4 yields the complex polarization identity. $\square$

---

## 4. Hilbert Spaces and Classical Examples

!!! info "Definition 7.2 (Hilbert Space)"

    If an inner product space is **complete** as a normed linear space, it is called a **Hilbert space**. If it is not complete, the inner product space is called a **pre-Hilbert space** (or quasi-Hilbert space).

### 4.1 Finite-Dimensional Unitary Space $\mathbb{C}^n$

For any two elements $x = (x_1, \dots, x_n)$ and $y = (y_1, \dots, y_n)$ in $\mathbb{C}^n$, define the inner product:

$$
(x, y) = \sum_{j=1}^n x_j \overline{y_j}
$$

It satisfies all conditions of an inner product. According to common algebraic terminology, $\mathbb{C}^n$ is called a unitary space, and it is a complete Hilbert space.

### 4.2 Sequence Space $l^2$

For any two elements $x = \{x_1, x_2, \dots\}, y = \{y_1, y_2, \dots\}$ in $l^2$. By the Cauchy-Schwarz inequality, the series converges absolutely:

$$
\sum_{j=1}^\infty |x_j y_j| \le \left( \sum_{j=1}^\infty |x_j|^2 \right)^{\frac{1}{2}} \left( \sum_{j=1}^\infty |y_j|^2 \right)^{\frac{1}{2}} = \|x\| \|y\| < \infty
$$

Thus, we can define the inner product:

$$
(x, y) = \sum_{j=1}^\infty x_j \overline{y_j}
$$

Under this inner product, $l^2$ is a complete, separable Hilbert space.

### 4.3 Function Space $L^2(F)$ and Weighted Space $L^2([a,b], \omega)$

For $f, g \in L^2(F)$, define the inner product:

$$
(f, g) = \int_F f(t) \overline{g(t)} dt
$$

For the weighted space, $\omega$ is a positive measurable function on $[a, b]$, and the square-integrable condition is $\|f\| = \left( \int_{[a,b]} |f(t)|^2 \omega(t) dt \right)^{1/2} < \infty$. Its inner product is defined as:

$$
(f, g)_\omega = \int_{[a, b]} f(t) \overline{g(t)} \omega(t) dt
$$

Both of these spaces are complete and separable Hilbert spaces.

---

## 5. Characteristics of Inner Product Spaces: Parallelogram Law

Conversely, given a normed linear space, what condition must its norm satisfy to guarantee that it is induced by some inner product? The answer is the "Parallelogram Law".

!!! success "Theorem 7.2 (Parallelogram Law)"

    Let $X$ be an inner product space. Then the norm induced by the inner product of $X$ must satisfy the **parallelogram law**: for any $x, y \in X$,

    $$
    \|x + y\|^2 + \|x - y\|^2 = 2(\|x\|^2 + \|y\|^2)
    $$

    **Conversely (Sufficiency)**: Let $X$ be a normed linear space. If the norm of $X$ satisfies the above equation, then an inner product $(x, y)$ can be defined in $X$ to make $X$ an inner product space, and the norm of $X$ is exactly the one induced by this inner product.

??? proof "Proof of the Theorem (Necessity and Sufficiency)"

    **1. Necessity (If it is an inner product space, it must satisfy the formula):**

    $$
    \|x + y\|^2 + \|x - y\|^2 = (x + y, x + y) + (x - y, x - y)
    $$

    $$
    = [(x, x) + (x, y) + (y, x) + (y, y)] + [(x, x) - (x, y) - (y, x) + (y, y)]
    $$

    $$
    = 2(x, x) + 2(y, y) = 2(\|x\|^2 + \|y\|^2)
    $$

    **2. Sufficiency (If the formula is satisfied, an inner product can be defined):**

    Inspired by the polarization identity, we naturally "define" the inner product as follows (taking the complex space as an example):

    $$
    (x, y) = \frac{1}{4} (\|x + y\|^2 - \|x - y\|^2 + i\|x + iy\|^2 - i\|x - iy\|^2)
    $$

    Now we need to verify that $(x,y)$ defined this way satisfies all axioms of an inner product.
    The prerequisite condition (parallelogram law) is equivalent to:

    $$
    \|a\|^2 + \|b\|^2 = 2\left\|\frac{a+b}{2}\right\|^2 + 2\left\|\frac{a-b}{2}\right\|^2
    $$

    **Proving additivity $(x+y, z) = (x, z) + (y, z)$**:
    For any $x, y, z \in X$, using the expansion of the definition:

    $$
    4(x, z) + 4(y, z) = \|x + z\|^2 - \|x - z\|^2 + i\|x + iz\|^2 - i\|x - iz\|^2
    $$

    $$
    + \|y + z\|^2 - \|y - z\|^2 + i\|y + iz\|^2 - i\|y - iz\|^2
    $$

    Grouping the terms containing $z$ and $-z$ respectively, and applying the parallelogram law:

    $$
    = 2\left\|\frac{x+y}{2} + z\right\|^2 - 2\left\|\frac{x+y}{2} - z\right\|^2 + 2i\left\|\frac{x+y}{2} + iz\right\|^2 - 2i\left\|\frac{x+y}{2} - iz\right\|^2
    $$

    This exactly equals $8\left(\frac{x+y}{2}, z\right)$.
    Therefore, $(x, z) + (y, z) = 2\left(\frac{x+y}{2}, z\right)$.

    By definition, it is obvious that $(\theta, y) = 0$. Letting $y = \theta$ in the above equation gives:

    $$
    (x, z) + (\theta, z) = 2\left(\frac{x}{2}, z\right) \implies (x, z) = 2\left(\frac{x}{2}, z\right)
    $$

    Substituting this back into the original equation yields strict additivity:

    $$
    (x, z) + (y, z) = (x + y, z)
    $$

    **Proving homogeneity $(\alpha x, y) = \alpha(x, y)$**:
    Let the function be $f(\alpha) = (\alpha x, y)$, where $\alpha$ is a real number. From the additivity proved above:

    $$
    f(\alpha + \beta) = ((\alpha + \beta)x, y) = (\alpha x + \beta x, y) = (\alpha x, y) + (\beta x, y) = f(\alpha) + f(\beta)
    $$

    Moreover, as $\alpha \to \beta$ (i.e., $\alpha - \beta \to 0$), it is easy to prove that $f(\alpha)$ is a continuous function by the continuity of the norm. According to the **auxiliary lemma** below, a continuous function satisfying $f(\alpha+\beta) = f(\alpha)+f(\beta)$ must satisfy $f(\alpha) = \alpha f(1)$.
    Thus, for real number $\alpha$, $(\alpha x, y) = \alpha(x, y)$ holds.

    For the complex case, we also need to verify $(ix, y) = i(x, y)$:

    $$
    4(ix, y) = \|ix + y\|^2 - \|ix - y\|^2 + i\|ix + iy\|^2 - i\|ix - iy\|^2
    $$

    $$
    = \|x - iy\|^2 - \|x + iy\|^2 + i\|x + y\|^2 - i\|x - y\|^2
    $$

    $$
    = i(\|x + y\|^2 - \|x - y\|^2 + i\|x + iy\|^2 - i\|x - iy\|^2) = 4i(x, y)
    $$

    Therefore, for any complex number $\alpha$, the homogeneity $(\alpha x, y) = \alpha(x, y)$ holds.
    Finally, from the definition, it is easy to directly verify conjugate symmetry $(x, y) = \overline{(y, x)}$ and non-negativity $(x, x) = \|x\|^2 \ge 0$. This completes the proof of sufficiency. $\square$

??? proof "Auxiliary Lemma: Linearity of the real function $f(\alpha+\beta) = f(\alpha)+f(\beta)$"

    **Lemma Content**: Let $f$ be a continuous real-valued function defined on $\mathbb{R}$, and for any real numbers $\alpha, \beta$, $f(\alpha+\beta) = f(\alpha)+f(\beta)$, then for any real number $\alpha$, $f(\alpha) = \alpha f(1)$.

    **Proof**: By assumption, for any natural number $n$ and real number $\alpha$, repeatedly applying additivity gives:

    $$
    f(n\alpha) = nf(\alpha)
    $$

    Taking $\alpha = \frac{1}{n}$, we get $f(1) = nf(\frac{1}{n})$, meaning $f(\frac{1}{n}) = \frac{1}{n}f(1)$.
    Thus, for any positive rational number $\frac{m}{n}$, we have:

    $$
    f\left(\frac{m}{n}\right) = m f\left(\frac{1}{n}\right) = \frac{m}{n}f(1)
    $$

    Also, because $f(0) = f(0+0) = 2f(0)$, we have $f(0) = 0$.
    From $f(\alpha) + f(-\alpha) = f(0) = 0$, we know $f(-\alpha) = -f(\alpha)$.
    Therefore, for any rational number $q = \frac{m}{n}$ (positive or negative), $f(q) = qf(1)$ holds.
    Since $f$ is continuous on the real axis and rational numbers are dense in real numbers, the conclusion holds for any real number $\alpha$. $\square$

### 5.1 Normed Spaces Failing the Parallelogram Law (Counterexamples)

Not all norms of normed linear spaces can be induced by an inner product.

* **Continuous Function Space $C[0, 1]$**:
    The norm is defined as $\|x\| = \max_{t \in [0, 1]} |x(t)|$.
    Let $f(t) = t, g(t) = 1 - t$, then $\|f\| = 1, \|g\| = 1$.
    Right side: $2(\|f\|^2 + \|g\|^2) = 4$.
    Left side: $f+g = 1 \implies \|f+g\|^2 = 1$; $f-g = 2t-1 \implies \max |2t-1| = 1 \implies \|f-g\|^2 = 1$.
    Therefore, $\|f+g\|^2 + \|f-g\|^2 = 2 \neq 4$. It does not satisfy the law.

* **Space $l^p$ (when $p \neq 2$)**:
    Take $x = e_1 = (1, 0, 0, \dots), y = e_2 = (0, 1, 0, \dots)$.
    Right side: $2(\|x\|_p^2 + \|y\|_p^2) = 2(1 + 1) = 4$.
    Left side: $\|x+y\|_p^2 + \|x-y\|_p^2 = \|(1, 1, 0, \dots)\|_p^2 + \|(1, -1, 0, \dots)\|_p^2 = (1^p+1^p)^{2/p} + (1^p+|-1|^p)^{2/p} = 2^{2/p} + 2^{2/p} = 2 \cdot 2^{2/p}$.
    Obviously, $2 \cdot 2^{2/p} = 4$ holds if and only if $p = 2$.

* **Space $L^p[0, 1]$ (when $p \neq 2$)**:
    Take the characteristic functions $f = \chi_{[0, 1/2]}, g = \chi_{[1/2, 1]}$.
    Using integration to calculate the norm, similar to the derivation for $l^p$, the left side gives $2$, and the right side gives $4 \cdot 2^{-2/p}$. The equation still only holds when $p=2$.

---

## 6. Orthogonality and Orthogonal Complement

The most important application of the inner product is introducing the concept of "orthogonality" (perpendicularity).

!!! info "Definition 7.3 (Orthogonality and Orthogonal Complement)"

    * Let $X$ be an inner product space, $x, y \in X$. If $(x, y) = 0$, then $x$ and $y$ are said to be **orthogonal**, denoted as $x \perp y$.

    * Let $M$ be a subset of $X$, and $x \in X$. If $x$ is orthogonal to any element in $M$, then $x$ is said to be orthogonal to $M$, denoted as $x \perp M$.

    * Let $N$ also be a subset of $X$. If for any $x \in M$ and any $y \in N$, we have $x \perp y$, then $M$ and $N$ are said to be orthogonal, denoted as $M \perp N$.

    * The set consisting of all elements in $X$ that are orthogonal to $M$ is called the **orthogonal complement** of $M$, denoted as $M^\perp$.

Based on the definition of orthogonality, the inner product space satisfies classical geometric theorems:

!!! success "Theorem 7.3 (Pythagorean Theorem and Properties in Inner Product Spaces)"

    * (i) Let the elements $x_1, \dots, x_n$ in $X$ be mutually orthogonal, and $x = x_1 + \dots + x_n$, then:

    $$
    \|x\|^2 = \|x_1\|^2 + \dots + \|x_n\|^2
    $$

    * (ii) If an element $x$ in $X$ is orthogonal to a dense subset $L$ in $X$, then $x = \theta$.

    * (iii) For any subset $M \subset X$, its orthogonal complement $M^\perp$ is a **closed subspace** of $X$.

??? proof "Proof of the Properties"

    **Proof of (i) Pythagorean Theorem:**
    When $n=2$, since $(x, y) = (y, x) = 0$:

    $$
    \|x + y\|^2 = (x + y, x + y) = (x, x) + (x, y) + (y, x) + (y, y) = \|x\|^2 + \|y\|^2
    $$

    For the general case, due to mutual orthogonality, $(x_1 + \dots + x_{n-1}, x_n) = 0$. By mathematical induction:

    $$
    \|x_1 + \dots + x_n\|^2 = \|x_1 + \dots + x_{n-1}\|^2 + \|x_n\|^2 = \dots = \|x_1\|^2 + \dots + \|x_n\|^2
    $$

    **Proof of (ii):**
    Since $L$ is dense in $X$, there exists a sequence $\{x_n\} \subset L$ converging to $x$. By the continuity of the inner product and $x \perp L$:

    $$
    (x, x) = \lim_{n \to \infty} (x, x_n) = 0
    $$

    According to the positive definiteness of the inner product, we have $x = \theta$.

    **Proof of (iii):**
    If $x, y \in M^\perp$, then for any given $z \in M$, we have $(\alpha x + \beta y, z) = \alpha(x, z) + \beta(y, z) = 0$, hence $M^\perp$ forms a linear subspace.
    Now suppose $x$ belongs to the closure of $M^\perp$, then there exists a sequence $\{x_n\}$ in $M^\perp$ converging to $x$. Using the continuity of the inner product, for any $z \in M$:

    $$
    (x, z) = \lim_{n \to \infty} (x_n, z) = 0
    $$

    So the limit point $x \in M^\perp$. Therefore, $M^\perp$ must be a closed set. $\square$

---

## 7. Element of Best Approximation

We will use the existence and uniqueness of the element of best approximation to pave the way for the "orthogonal decomposition" in the next section. For this purpose, we first introduce the concepts of the element of best approximation and convex sets.

!!! info "Definition 7.4 (Element of Best Approximation and Closed Convex Set)"

    * **Element of Best Approximation**: Let $M$ be a subset of an inner product space, and $x \in X$ be a given element. If there exists an element $y \in M$ such that
    $\|x - y\| = \inf_{z \in M} \|x - z\|$
    then $y$ is called the element of best approximation for $x$ in $M$.

    * **Convex Set**: Let $M$