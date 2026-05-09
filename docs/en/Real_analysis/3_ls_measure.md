# 1.3 Lebesgue-Stieltjes Measure and Lebesgue Measure

In the previous two sections, we established the abstract extension theory from outer measures to measure spaces. In this section, we apply this theory to the set of real numbers $\mathbb{R}$ to construct an extremely important class of Borel measures, the Lebesgue-Stieltjes measures, through monotonically increasing functions, and derive the foundation of analysis—the Lebesgue measure.

---

## 1. Distribution Functions and Borel Measures

A measure defined on the Borel $\sigma$-algebra $\mathcal{B}_{\mathbb{R}}$ of $\mathbb{R}$ is called a **Borel measure**. If this measure is finite on any bounded interval, it corresponds one-to-one with a special class of real functions—distribution functions.

!!! info "Definition 1.3.1 (Distribution Function)"

    Let $F: \mathbb{R} \to \mathbb{R}$ be a function. $F$ is called a distribution function if it satisfies the following properties:

    * **Monotonicity**: $F$ is monotonically increasing, i.e., if $x \le y$, then $F(x) \le F(y)$.

    * **Right-continuity**: For any $x_0 \in \mathbb{R}$, $F(x_0^+) = \lim_{x \to x_0^+} F(x) = F(x_0)$.

    Through the distribution function, we can define the measure value of a half-open interval as:

    \[
    \mu((a, b]) = F(b) - F(a)
    \]

---

## 2. Construction of Lebesgue-Stieltjes Measure

The first step of construction is defining a pre-measure on an algebra. Let $\mathcal{A}$ be the family consisting of finite disjoint unions of half-open intervals of the form $(a, b]$ (including $(a, \infty)$ and $(-\infty, b]$). It is easy to prove that $\mathcal{A}$ is an algebra.

!!! success "Theorem 1.3.1 (Establishment of Pre-measure)"

    Given a distribution function $F$, define the set function $\mu_0$ on the algebra $\mathcal{A}$:

    \[
    \mu_0\left( \bigcup_{j=1}^n (a_j, b_j] \right) = \sum_{j=1}^n [F(b_j) - F(a_j)]
    \]

    Then $\mu_0$ is a pre-measure on $\mathcal{A}$.

??? proof "Key Points of the Proof of Theorem 1.3.1"

    The core of proving a pre-measure lies in verifying **countable additivity**. Let $I = (a, b]$ and $I = \bigcup_{j=1}^\infty I_j$, where $I_j = (a_j, b_j]$ are disjoint intervals.

    **1. Finite Additivity**:

    Using the monotonicity of $F$, it is easy to prove that $\mu_0$ is additive for finite unions. For a union of $n$ terms, by organizing the endpoints, we obtain $\mu_0(\cup_{j=1}^n I_j) = \sum_{j=1}^n \mu_0(I_j)$.

    **2. Proof of Countable Additivity (Finite Covering Principle)**:

    * **Inequality I**: Since finite unions are contained in countable unions, monotonicity implies $\mu_0(I) \ge \sum_{j=1}^n \mu_0(I_j)$. Letting $n \to \infty$ yields $\mu_0(I) \ge \sum_{j=1}^\infty \mu_0(I_j)$.

    * **Inequality II**: Take any $\epsilon > 0$. Since $F$ is right-continuous, for $I = (a, b]$, there exists $\delta > 0$ such that $F(a+\delta) - F(a) < \epsilon$, thus the compact set $K = [a+\delta, b] \subset (a, b]$.

    * For each $I_j = (a_j, b_j]$, there exists $\delta_j > 0$ such that $F(b_j+\delta_j) - F(b_j) < \epsilon / 2^j$, thus the open interval $O_j = (a_j, b_j+\delta_j) \supset (a_j, b_j]$.

    * Since the compact set $K$ is covered by the family of open intervals $\{O_j\}$, by the **Heine-Borel Theorem**, there exist finitely many intervals $O_{j_1}, \dots, O_{j_N}$ covering $K$. Using the property of finite additivity and letting $\epsilon \to 0$, we eventually prove $\mu_0(I) \le \sum_{j=1}^\infty \mu_0(I_j)$.

!!! success "Theorem 1.3.2 (Lebesgue-Stieltjes Measure Extension)"

    The measure obtained by restricting the outer measure $\mu_F^*$, induced by the pre-measure $\mu_0$, to the family of Carathéodory measurable sets $\mathcal{M}_F$ is called the **Lebesgue-Stieltjes measure**. By the uniqueness theorem of extension, if $F$ is fixed, the restriction of this measure to the Borel $\sigma$-algebra is unique.

---

## 3. Regularity and Approximation Theorems

Lebesgue-Stieltjes measures possess good topological approximation properties, which allow us to study complex measurable sets using open or compact sets with simple structures.

!!! success "Theorem 1.3.3 (Regularity Regularity)"

    For any $E \in \mathcal{M}_F$, its measure satisfies:

    * **Outer Regularity**: $\mu_F(E) = \inf \{ \mu_F(U) : E \subset U, U \text{ is an open set} \}$.

    * **Inner Regularity**: $\mu_F(E) = \sup \{ \mu_F(K) : K \subset E, K \text{ is a compact set} \}$.

By corollary, any measurable set $E$ has the following structural descriptions:

* **$G_\delta$ Approximation**: There exists a $G_\delta$ set $V$ (countable intersection of open sets) such that $E \subset V$ and $\mu_F(V \setminus E) = 0$.

* **$F_\sigma$ Approximation**: There exists an $F_\sigma$ set $H$ (countable union of closed sets) such that $H \subset E$ and $\mu_F(E \setminus H) = 0$.

---

## 4. Lebesgue Measure (Lebesgue Measure)

When the distribution function is the identity function, i.e., $F(x) = x$, the resulting measure is the classical **Lebesgue measure**.

!!! info "Definition 1.3.2 (Lebesgue Measure)"

    The measure space $(\mathbb{R}, \mathcal{L}, m)$ induced by $F(x) = x$ is called the Lebesgue measure space. For any half-open interval $(a, b]$, its measure is exactly equal to its geometric length:

    \[
    m((a, b]) = b - a
    \]

### Invariance Properties of Lebesgue Measure

The Lebesgue measure is compatible with the addition and scalar multiplication structures of real numbers, exhibiting excellent symmetry.

!!! success "Theorem 1.3.4 (Translation and Scaling Invariance)"

    Let $E$ be a Lebesgue measurable set, for any $s, r \in \mathbb{R}$:

    * **Translation Invariance**: $m(E+s) = m(E)$, where $E+s = \{x+s : x \in E\}$.

    * **Scaling Property**: $m(rE) = |r|m(E)$, where $rE = \{rx : x \in E\}$.

---

## 5. Cantor Set and Cantor Function

The Cantor set is the most famous counterexample in measure theory for discussing "size" and "cardinality."

!!! info "Definition 1.3.3 (Cantor Ternary Set)"

    Starting from the closed interval $[0, 1]$, the first step is to remove the middle open interval $(\frac{1}{3}, \frac{2}{3})$; the second step is to remove the middle third of the remaining two segments, and so on indefinitely. The resulting set of points $C$ is called the **Cantor set**.

!!! success "Properties of the Cantor Set"

    * **Zero Measure**: The measure of the Cantor set is 0. Because the total length of the removed intervals is:

    \[
    \sum_{n=1}^\infty \frac{2^{n-1}}{3^n} = \frac{1/3}{1-2/3} = 1
    \]

    * **Uncountability**: Although the measure is 0, the cardinality of the Cantor set is the same as that of $[0, 1]$, i.e., the cardinality is $c$.

    * **Topological Properties**: $C$ is a compact set, nowhere dense, and has no isolated points (a perfect set).

!!! info "Cantor Function (Cantor Function / Devil's Staircase)"

    Based on the construction of the Cantor set, a continuous monotonic function $f: [0, 1] \to [0, 1]$ can be defined:

    * $f(x)$ is constant on the complement of the Cantor set (the removed intervals).

    * $f(0) = 0, f(1) = 1$.

    * $f'(x) = 0$ holds almost everywhere, yet the function grows continuously. This illustrates that "the derivative is 0 almost everywhere" does not necessarily mean the function is constant.