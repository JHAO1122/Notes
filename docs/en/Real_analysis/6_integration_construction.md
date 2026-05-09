# 3.1 Construction of Abstract Integrals and Simple Functions

In Riemann integration, the integral is constructed by partitioning the domain into intervals (the vertical slicing method). In Lebesgue's abstract integration theory, the integral is constructed by partitioning the range (the horizontal slicing method). The core bridge in this construction is the **simple function**, which serves as the foundation linking measure and integration.

---

## 1. Characteristic Functions and Simple Functions

!!! info "Definition 3.1.1 (Characteristic Function and Simple Function)"

    * **Characteristic function**: For a subset \(E \in \mathcal{M}\) of a measure space \(X\), define its characteristic function as

    \[
    \chi_E(x) = \begin{cases} 1, & x \in E \\ 0, & x \notin E \end{cases}
    \]

    * **Simple function**: A measurable function \(\phi: X \to \mathbb{C}\) is called a simple function if its range is a finite set. Any simple function can be expressed as a finite linear combination of characteristic functions:

    \[
    \phi(x) = \sum_{j=1}^n a_j \chi_{E_j}(x)
    \]

    where \(a_j \in \mathbb{C}\) and the sets \(E_j = \{x \in X : \phi(x) = a_j\} \in \mathcal{M}\) form a partition of \(X\) into pairwise disjoint measurable sets. This is called the **standard representation** of a simple function.

The importance of simple functions lies in the fact that every measurable function can be perfectly approximated by simple functions.

!!! success "Theorem 3.1.1 (Simple Function Approximation Theorem)"

    Let \((X, \mathcal{M})\) be a measurable space.

    * For any non‑negative measurable function \(f: X \to [0, +\infty]\), there exists a sequence of non‑negative simple functions \(\{\phi_n\}_{n=1}^\infty\) such that they increase monotonically and converge pointwise to \(f\):

    \[
    0 \le \phi_1 \le \phi_2 \le \dots \le f, \quad \lim_{n \to \infty} \phi_n(x) = f(x)
    \]

    * If \(f\) is a general complex‑valued measurable function, there exists a sequence of complex simple functions \(\{\psi_n\}\) such that \(0 \le |\psi_1| \le |\psi_2| \le \dots \le |f|\) and \(\psi_n \to f\).

---

## 2. Integral of Simple Functions

We begin by defining “area” or “integral” from the simplest building blocks.

!!! info "Definition 3.1.2 (Integral of a Non‑negative Simple Function)"

    Let \((X, \mathcal{M}, \mu)\) be a measure space and let \(\phi = \sum_{j=1}^n a_j \chi_{E_j}\) be the standard representation of a non‑negative simple function (\(a_j \ge 0\)). The integral of \(\phi\) with respect to \(\mu\) is defined as

    \[
    \int_X \phi \, d\mu = \sum_{j=1}^n a_j \mu(E_j)
    \]

    **Note**: In measure theory we adopt the convention \(0 \cdot \infty = 0\). This means that even if a set has infinite measure, as long as the function takes the value 0 on that set, its contribution to the integral is 0. For a measurable subset \(A \subset X\), the local integral is defined by \(\int_A \phi \, d\mu = \int_X \phi \cdot \chi_A \, d\mu\).

The integral of simple functions possesses a perfect algebraic structure.

!!! success "Theorem 3.1.2 (Properties of the Integral of Simple Functions)"

    Let \(\phi, \psi\) be non‑negative simple functions and let \(c \ge 0\) be a constant. Then:

    * **Linearity**: \(\int c\phi = c \int \phi\) and \(\int (\phi + \psi) = \int \phi + \int \psi\).

    * **Monotonicity**: If \(\phi \le \psi\), then \(\int \phi \le \int \psi\).

    * **Set‑function property**: The mapping \(A \mapsto \nu(A) = \int_A \phi \, d\mu\) defines a new measure on \(\mathcal{M}\).

---

## 3. Integral of Non‑negative Measurable Functions

Using suprema, we extend the integral from simple functions to arbitrary non‑negative measurable functions.

!!! info "Definition 3.1.3 (Integral of a Non‑negative Measurable Function)"

    Let \(f: X \to [0, +\infty]\) be a non‑negative measurable function (denoted \(f \in L^+\)). The integral of \(f\) is defined as the supremum of the integrals of all non‑negative simple functions that lie below \(f\):

    \[
    \int_X f \, d\mu = \sup \left\{ \int_X \phi \, d\mu : 0 \le \phi \le f,\ \phi \text{ is a simple function} \right\}
    \]

For integrals of non‑negative measurable functions, we have the following core property involving “almost everywhere” (a.e.).

!!! success "Proposition 3.1.1 (Almost Everywhere Properties of Non‑negative Integrals)"

    Let \(f \in L^+\). Then:

    * \(\int f \, d\mu = 0\) if and only if \(f = 0\) almost everywhere (\(f = 0\) a.e.).

    * If \(\int f \, d\mu < \infty\), then \(f < \infty\) almost everywhere (i.e., an integrable function is necessarily finite almost everywhere).

??? proof "Proof: \(\int f = 0 \iff f = 0 \text{ a.e.}\)"

    **Sufficiency (\(\Leftarrow\))**: If \(f = 0\) a.e., then for any simple function \(\phi = \sum a_j \chi_{E_j}\) satisfying \(0 \le \phi \le f\), whenever \(a_j > 0\) we must have \(\mu(E_j) = 0\). Hence \(\int \phi = 0\), and by the supremum definition we obtain \(\int f = 0\).

    **Necessity (\(\Rightarrow\))**: Assume \(\int f = 0\). Define the sets \(E_n = \{x : f(x) \ge \frac{1}{n}\}\). Clearly \(\frac{1}{n} \chi_{E_n} \le f\).
  
    By monotonicity of the integral:
  
    \[
    \int f \ge \int \frac{1}{n} \chi_{E_n} = \frac{1}{n} \mu(E_n)
    \]
  
    Since \(\int f = 0\), we must have \(\mu(E_n) = 0\).
    Now \(\{x : f(x) > 0\} = \bigcup_{n=1}^\infty E_n\), and by subadditivity of the measure this set has measure zero, i.e., \(f = 0\) a.e.

---

## 4. General Integrable Functions and \(L^1\) Space

With the integral of non‑negative functions in hand, we can define the integral of general real‑valued or complex‑valued functions by splitting into positive and negative parts.

!!! info "Definition 3.1.4 (Real‑valued Integrable Functions and \(L^1\) Space)"

    Let \(f: X \to \mathbb{R}\) be a real‑valued measurable function. Decompose it into its positive and negative parts: \(f = f^+ - f^-\) (so that \(|f| = f^+ + f^-\)).

    If both the positive and negative parts have finite integrals, i.e., \(\int f^+ < \infty\) and \(\int f^- < \infty\), then \(f\) is said to be **Lebesgue integrable** on \(X\). The space of all such functions is denoted by \(L^1(\mu)\).

    In this case, the integral of \(f\) is defined as

    \[
    \int_X f \, d\mu = \int_X f^+ \, d\mu - \int_X f^- \, d\mu
    \]

    Clearly, \(f \in L^1(\mu)\) if and only if \(\int |f| \, d\mu < \infty\).

For a complex‑valued function \(f(x) = u(x) + i v(x)\), if \(|f| \in L^1(\mu)\), then its integral is defined by \(\int f = \int u + i \int v\).

!!! success "Theorem 3.1.3 (Absolute Inequality for Integrals)"

    If \(f \in L^1(\mu)\), then the following absolute inequality holds:

    \[
    \left| \int_X f \, d\mu \right| \le \int_X |f| \, d\mu
    \]

??? proof "Proof of the Absolute Inequality (Complex Case)"

    Let the value of the integral be the complex number \(z = \int f d\mu\). If \(z = 0\), the inequality is trivially true.
    If \(z \ne 0\), write \(z = r e^{i\theta}\) with \(r > 0\), and set \(\alpha = e^{-i\theta}\). Then \(\alpha z = |\int f d\mu|\).
    Since \(\alpha z\) is real, we have

    \[
    \left| \int f \right| = \alpha \int f = \int (\alpha f) = \text{Re} \left( \int \alpha f \right) = \int \text{Re}(\alpha f)
    \]

    Because for any complex number \(\text{Re}(\alpha f) \le |\alpha f| = |f|\), the monotonicity of the integral yields

    \[
    \int \text{Re}(\alpha f) \le \int |f|
    \]

    Combining these steps proves the absolute inequality.

!!! success "Proposition 3.1.2 (Invariance of the Integral under a.e. Equivalence)"

    Let \(f, g\) be measurable functions. If \(f = g\) almost everywhere (\(f = g\) a.e.) and \(f \in L^1(\mu)\), then \(g \in L^1(\mu)\) and they have the same integral:

    \[
    \int_X f \, d\mu = \int_X g \, d\mu
    \]

This property shows that the Lebesgue integral is completely “immune” to changes on sets of measure zero. In functional analysis, functions that are equal almost everywhere are usually identified as the same element in the \(L^1\) space.