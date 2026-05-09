# 1.2 Outer Measure and Carathéodory Extension

In the process of establishing measure theory, a core challenge is how to construct a measure on a sufficiently large $\sigma$-algebra. Carathéodory extension theory provides a standard method: first define an "outer measure" covering all subsets, then filter out those well-behaved sets through an ingenious measurability criterion, thereby forming a $\sigma$-algebra.

---

## 1. Outer Measure

An outer measure is a set function defined on all subsets of the universal set $X$. Unlike a measure, it only needs to satisfy subadditivity rather than countable additivity.

!!! info "Definition 1.2.1 (Outer Measure)"

    Let $X$ be a non-empty set. A set function $\mu^*: 2^X \to [0, +\infty]$ is called an **outer measure**, if it satisfies the following properties:

    * **Empty Set Property**: $\mu^*(\emptyset) = 0$.

    * **Monotonicity**: If $A \subset B \subset X$, then $\mu^*(A) \le \mu^*(B)$.

    * **Countable Subadditivity**: For any sequence of subsets $\{A_j\}_{j=1}^\infty \subset 2^X$, we have:

    \[
    \mu^*\left( \bigcup_{j=1}^\infty A_j \right) \le \sum_{j=1}^\infty \mu^*(A_j)
    \]

### General Construction Method for Outer Measures

We can generate an outer measure through an elementary function $\rho$ on a simple family of sets $\Sigma$. This method is crucial in constructing the Lebesgue measure.

!!! success "Theorem 1.2.1 (Generation of Outer Measure)"

    Let $\Sigma \subset 2^X$ be a family of subsets, satisfying $\emptyset \in \Sigma$. Let $\rho: \Sigma \to [0, +\infty]$ satisfy $\rho(\emptyset) = 0$. For any $A \subset X$, define:

    \[
    \mu^*(A) = \inf \left\{ \sum_{j=1}^\infty \rho(E_j) : A \subset \bigcup_{j=1}^\infty E_j, E_j \in \Sigma \right\}
    \]

    Then $\mu^*$ is an outer measure on $X$.

---

## 2. Carathéodory Measurability Criterion

Since outer measures generally do not satisfy additivity on arbitrary sets, Carathéodory proposed a criterion to define which sets are "measurable."

!!! info "Definition 1.2.2 (Carathéodory Measurability)"

    Let $\mu^*$ be an outer measure on $X$. A set $A \subset X$ is said to be **$\mu^*$-measurable**, if for any testing set $E \subset X$, it satisfies:

    \[
    \mu^*(E) = \mu^*(E \cap A) + \mu^*(E \cap A^c)
    \]

    By convention, the collection of all $\mu^*$-measurable sets is denoted by $\mathcal{M}$.



**Intuitive Understanding**: A set $A$ is measurable if it is able to split any set $E$ in the space "neatly," such that the sum of the measures of the two parts of $E$ inside and outside of $A$ equals the measure of the whole $E$. To verify measurability, one only needs to prove the inequality $\mu^*(E) \ge \mu^*(E \cap A) + \mu^*(E \cap A^c)$, as the other half of the inequality holds automatically by subadditivity.

---

## 3. Carathéodory Extension Theorem

This is the cornerstone theorem of measure theory, revealing that the family of measurable sets filtered from an outer measure possesses an excellent algebraic structure.

!!! success "Theorem 1.2.2 (Carathéodory Theorem)"

    Let $\mu^*$ be an outer measure on $X$, and $\mathcal{M}$ be the family consisting of all $\mu^*$-measurable sets. Then:

    * **Structure**: $\mathcal{M}$ is a $\sigma$-algebra.

    * **Measurability**: $\mu^*$ restricted on $\mathcal{M}$ is a measure (satisfies countable additivity).

    * **Completeness**: The measure space $(X, \mathcal{M}, \mu^*)$ is complete, i.e., all subsets of $\mu^*$-null sets belong to $\mathcal{M}$.

??? proof "Outline of Proof for Carathéodory Theorem"

    **1. Prove that $\mathcal{M}$ is an algebra**:
    Since $A$ and $A^c$ are symmetric in the definition, it is obvious that $A \in \mathcal{M} \implies A^c \in \mathcal{M}$. For closure under union, one can pick a testing set $E$, use the measurability of $A$ for the first split, then use the measurability of $B$ to split the remaining parts, and combine terms via subadditivity.

    **2. Prove finite additivity**:
    Let $A, B \in \mathcal{M}$ and $A \cap B = \emptyset$. Utilizing the splitting effect of $A$ on $E \cap (A \cup B)$, we obtain $\mu^*(E \cap (A \cup B)) = \mu^*(E \cap A) + \mu^*(E \cap B)$.

    **3. Extend to countable additivity**:
    Let $\{A_j\}$ be a sequence of disjoint sets in $\mathcal{M}$. Let $B_n = \bigcup_{j=1}^n A_j$ and $B = \bigcup_{j=1}^\infty A_j$.
    By induction, it follows that $\mu^*(E \cap B_n) = \sum_{j=1}^n \mu^*(E \cap A_j)$.
    Since $\mu^*(E) = \mu^*(E \cap B_n) + \mu^*(E \cap B_n^c) \ge \sum_{j=1}^n \mu^*(E \cap A_j) + \mu^*(E \cap B^c)$.
    Let $n \to \infty$, and utilizing the subadditivity $\mu^*(E \cap B) \le \sum \mu^*(E \cap A_j)$, the squeeze theorem confirms $\mu^*(E) = \mu^*(E \cap B) + \mu^*(E \cap B^c)$ and the measure satisfies countable additivity.

---

## 4. Pre-measure and Uniqueness of Extension

In practical applications, measures usually start from a smaller family of sets (such as an algebra $\mathcal{A}$).

!!! info "Definition 1.2.3 (Pre-measure)"

    Let $\mathcal{A}$ be an algebra on $X$. A set function $\mu_0: \mathcal{A} \to [0, +\infty]$ is called a **pre-measure**, if it satisfies $\mu_0(\emptyset) = 0$, and for any disjoint sequence of sets $\{A_j\}$ in $\mathcal{A}$, as long as $\bigcup A_j \in \mathcal{A}$, we have:

    \[
    \mu_0\left(\bigcup_{j=1}^\infty A_j\right) = \sum_{j=1}^\infty \mu_0(A_j)
    \]

!!! success "Theorem 1.2.3 (Uniqueness of Extension)"

    Let $\mu_0$ be a pre-measure on algebra $\mathcal{A}$, and $\mu^*$ be the outer measure induced by $\mu_0$.

    * **Existence**: The restriction of $\mu^*$ on the $\sigma$-algebra $\sigma(\mathcal{A})$ generated by $\mathcal{A}$ is a measure, and this measure coincides with $\mu_0$ on $\mathcal{A}$.

    * **Uniqueness**: If $\mu_0$ is **$\sigma$-finite** (i.e., $X$ can be expressed as a countable union of sets in $\mathcal{A}$ with finite measure), then the extension of this measure on $\sigma(\mathcal{A})$ is unique.

If the pre-measure does not satisfy $\sigma$-finiteness, the measure after extension may not be unique. This conclusion emphasizes the importance of $\sigma$-finiteness in the study of measure theory.