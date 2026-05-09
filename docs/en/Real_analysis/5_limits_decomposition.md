# 2.2 Limits and Decompositions of Sequences of Measurable Functions

## 5. Complex Measurable Functions and Extended Real-Valued Functions

For complex-valued functions $f: X \to \mathbb{C}$, we utilize the product space topology of $\mathbb{R}^2$ for the complex plane, namely $\mathcal{B}_{\mathbb{C}} = \mathcal{B}_{\mathbb{R}} \otimes \mathcal{B}_{\mathbb{R}}$.

!!! success "Proposition 2.1.4 (Complex Measurability)"

    A complex-valued function $f$ is measurable if and only if its real part $\text{Re } f$ and imaginary part $\text{Im } f$ are both real measurable functions.

In measure theory, we often allow functions to take values of $\pm \infty$. Define the extended real number set $\bar{\mathbb{R}} = \mathbb{R} \cup \{-\infty, +\infty\}$, whose Borel $\sigma$-algebra is defined as:

$$
\mathcal{B}_{\bar{\mathbb{R}}} = \{E \subset \bar{\mathbb{R}} : E \cap \mathbb{R} \in \mathcal{B}_{\mathbb{R}}\}
$$

若 $f$ 为广义实可测函数，我们通常约定 $0 \cdot \infty = 0$ 且避免无意义的 $\infty - \infty$。

---

## 6. Algebraic Operations and Limits of Sequences

Measurable functions possess perfect closure under various operations and limit-taking procedures.

!!! success "Theorem 2.1.3 (Closure under Algebraic Operations)"

    Let $f, g: X \to \mathbb{C}$ be measurable functions, then $f+g$ and $f \cdot g$ are both measurable functions.

??? proof "Proof of Theorem 2.1.3"

    We define the mapping $F(x) = (f(x), g(x)): X \to \mathbb{C} \times \mathbb{C}$.
    Since $f$ and $g$ are both measurable, according to the measurability of product spaces (Proposition 2.1.3), the joint mapping $F$ is measurable.
    Furthermore, since the complex addition mapping $S: \mathbb{C} \times \mathbb{C} \to \mathbb{C}$, i.e., $S(z, w) = z + w$, is a continuous mapping, it is therefore Borel measurable.
    By the composition property, $f+g = S \circ F$ must be measurable. The measurability of the product $f \cdot g$ can be proved analogously.

!!! success "Theorem 2.1.4 (Supremum/Infimum and Limits of Function Sequences)"

    Let $\{f_j\}_{j=1}^\infty$ be a sequence of extended real-valued measurable functions defined on $X$, then the following functions are also measurable:

    * $g_1(x) = \sup_{j} f_j(x)$
    * $g_2(x) = \inf_{j} f_j(x)$
    * $g_3(x) = \limsup_{j \to \infty} f_j(x)$
    * $g_4(x) = \liminf_{j \to \infty} f_j(x)$

??? proof "Proof of Theorem 2.1.4"

    We only prove the measurability of the supremum $\sup_j f_j$; the rest can be derived from algebraic relations.
    For any $a \in \mathbb{R}$, observe the set $\{x : \sup_j f_j(x) > a\}$.
    The supremum is greater than $a$ if and only if there exists at least one term in the sequence greater than $a$. Therefore:

    $$
    \{x : \sup_j f_j(x) > a\} = \bigcup_{j=1}^\infty \{x : f_j(x) > a\}
    $$

    Since each $f_j$ is measurable, every set $\{x : f_j(x) > a\}$ on the right side of the equation is in $\mathcal{M}$. Because $\mathcal{M}$ is closed under countable unions, the supremum function $g_1$ satisfies the equivalent condition for measurability (Proposition 2.1.2(b)), and is thus measurable.

    Remaining proof:
    $\inf_j f_j = -\sup_j (-f_j)$, which follows from algebraic operations and the measurability of the $\sup$.
    Since $\limsup_{j \to \infty} f_j = \inf_{k \ge 1} (\sup_{j \ge k} f_j)$, the measurability is naturally preserved after one countable $\sup$ and one countable $\inf$ operation. The same applies to $\liminf$.

!!! info "Corollary (Positive/Negative Parts and Polar Decomposition)"

    For any real measurable function $f$, we can decompose it into positive and negative parts:

    * **Positive part**: $f^+ = \max(f, 0)$
    * **Negative part**: $f^- = \max(-f, 0)$
    * At this point, we have $f = f^+ - f^-$ and $|f| = f^+ + f^-$. These two functions are also measurable.

    For complex measurable functions, we can perform polar decomposition (sign decomposition):
    
    $$
    f = (\text{sgn } f) |f|
    $$
    
    where $\text{sgn } z = z/|z|$ when $z \ne 0$, and 0 otherwise. Both $\text{sgn } f$ and $|f|$ are also measurable functions.