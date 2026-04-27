# Chapter 3: Pre-compact Sets, Compact Sets, and Their Topological Properties

In the real number field $\mathbb{R}^n$, completeness has a very intuitive equivalent proposition (the Bolzano-Weierstrass Theorem): **An infinite bounded set has at least one limit point**. However, in a general metric space (even a complete metric space), this property does not necessarily hold. This indicates that the topological structure of general metric spaces is far more complex than that of finite-dimensional Euclidean spaces.

To recover excellent properties similar to those of "bounded closed intervals" in abstract spaces, we introduce an extremely core concept in functional analysis—**Compactness**.

---

## 1. Basic Concepts of Pre-compact and Compact Sets

First, we need to clarify the definition of boundedness in general metric spaces.

!!! info "Definition 3.1 (Bounded Set)"

    A subset $A$ of a metric space $X$ is called a **Bounded Set** if $A$ is contained within some open or closed ball in $X$.
    
    *(Note: Any convergent sequence and any Cauchy sequence in a metric space are bounded.)*

!!! info "Definition 3.2 (Pre-compact Set and Compact Set)"

    Let $A$ be a subset of a metric space $X$.
    
    * **Pre-compact Set**: If **every sequence** in $A$ contains a subsequence that converges to some point in $X$, then $A$ is called a pre-compact set.
    * **Compact Set**: If **every sequence** in $A$ contains a subsequence that converges to some point **in $A$**, then $A$ is called a compact set.
    
    If the space $X$ itself is a compact set, then $X$ is called a compact metric space.

**Basic properties of compact sets:**

1. Any finite set is a compact set;

2. Any subset of a pre-compact set is a pre-compact set;

3. Any closed subset of a compact set is a compact set.

!!! warning "Counterexample: A bounded closed set is not necessarily pre-compact"

    In infinite-dimensional spaces, a bounded set does not necessarily have a convergent subsequence.
    Consider the space $L^2[-\pi, \pi]$ and define the sequence of functions:

    \[
    e_k(t) = \frac{1}{\sqrt{\pi}} \sin(kt), \quad k = 1, 2, \dots
    \]

    This sequence is bounded (the norm of each element is 1). But calculating the squared distance between any two distinct elements:

    \[
    \rho(e_j, e_k)^2 = \frac{1}{\pi} \int_{-\pi}^{\pi} |\sin(jt) - \sin(kt)|^2 dt
    \]

    Expanding and utilizing the orthogonality of trigonometric functions:

    \[
    = \frac{1}{\pi} \int_{-\pi}^{\pi} (\sin^2(jt) - 2\sin(jt)\sin(kt) + \sin^2(kt)) dt = 1 - 0 + 1 = 2
    \]

    Therefore, when $j \ne k$, $\rho(e_j, e_k) = \sqrt{2}$.
    This shows that the distance between any two points in the sequence is exactly $\sqrt{2}$. It is absolutely impossible for a Cauchy subsequence to exist, and naturally, no convergent subsequence can exist. Thus, this bounded sequence is not a pre-compact set.

---

## 2. $\epsilon$-Nets and Totally Bounded Sets

To characterize pre-compactness more precisely, we introduce the concept of a "net", which is equivalent to using a finite number of points to "cover" the entire set.

!!! info "Definition 3.3 ($\epsilon$-Net)"

    Let $X$ be a metric space, $A, B \subset X$. Let $\epsilon > 0$ be a given positive number.
    If for any point $x$ in $A$, there exists a point $y$ in $B$ such that $\rho(x, y) < \epsilon$, then $B$ is called an **$\epsilon$-net** of $A$.

    *(In other words, the union of open balls centered at all points in $B$ with radius $\epsilon$ contains $A$. Note: $B$ does not necessarily need to be contained in $A$.)*

!!! success "Definition 3.4 (Totally Bounded Set)"

    Let $A$ be a subset of a metric space $X$. If for any given $\epsilon > 0$, there always exists a **finite $\epsilon$-net** for $A$, then $A$ is called a **totally bounded set**.

Totally bounded sets have the following intuitive properties:

* Any finite set is totally bounded;

* A subset of a totally bounded set is also totally bounded;

* Let $A$ be a totally bounded set. For any given $\epsilon > 0$, we can always take a **finite subset of $A$ itself** as an $\epsilon$-net of $A$.

??? proof "Proof: The net of a totally bounded set can be taken from within itself"

    Since $A$ is totally bounded, for a given $\epsilon > 0$, there exists a finite $\epsilon/2$-net $\{x_1, \dots, x_{n_0}\}$ in $X$.
    Without loss of generality, assume that for each $k$, the ball intersects $A$ (if there is no intersection, the ball can simply be discarded):

    \[
    A \cap S\left(x_k, \frac{\epsilon}{2}\right) \ne \emptyset
    \]

    We arbitrarily pick a point $y_k \in A \cap S(x_k, \epsilon/2)$ from each intersection.
    By the triangle inequality, this finite set of points $\{y_1, \dots, y_{n_0}\} \subset A$ serves as a finite $\epsilon$-net for $A$. $\square$

!!! success "Theorem 3.1"

    A totally bounded set must be bounded and separable.

??? proof "Proof of Theorem 3.1 (Click to expand)"

    **1. Proof of boundedness:**
    Suppose $A \subset X$ is totally bounded. Thus, $A$ has a finite 1-net $B_1 = \{x_1, \dots, x_{n_0}\}$.
    For any given $x \in A$, there exists $x_k \in B_1$ such that $\rho(x, x_k) < 1$. By the triangle inequality, for a fixed point $x_{n_0}$:

    \[
    \rho(x, x_{n_0}) \le \rho(x, x_k) + \rho(x_k, x_{n_0}) \le 1 + \max_{1 \le k \le n_0} \rho(x_k, x_{n_0}) < \infty
    \]

    Since the right side is the maximum of a finite number of constants, it shows that $A$ is contained within a ball of finite radius, hence $A$ is bounded.

    **2. Proof of separability:**
    For each natural number $n$, $A$ has a finite $1/n$-net, denoted as $B_n \subset A$.
    Let $B = \bigcup_{n=1}^\infty B_n$. Since the union of countably many finite sets is still an at most countable set, $B$ is a countable set.
    For any $x \in A$ and any $\epsilon > 0$, choose $k$ such that $1/k < \epsilon$. Since $B_k$ is a $1/k$-net, there must exist $x_k \in B_k \subset B$ such that $\rho(x, x_k) < 1/k < \epsilon$. This shows that $B$ is dense in $A$, hence $A$ is separable. $\square$

---

## 3. Pre-compactness and Compactness in Complete Spaces

Pre-compactness and total boundedness have a one-way implication in general metric spaces, but they are entirely equivalent in complete spaces.

!!! success "Theorem 3.2 (Equivalence of Pre-compactness and Total Boundedness)"

    Let $X$ be a metric space.
    1. If $A \subset X$ is pre-compact, then $A$ is totally bounded.
    2. If $X$ is a **complete** metric space, then when $A$ is totally bounded, $A$ must be pre-compact.

??? proof "Proof of Theorem 3.2 (Click to expand)"

    **(1) Pre-compact $\implies$ Totally bounded (Proof by contradiction):**
    Assume $A$ is pre-compact but not totally bounded. Then there must exist some $\epsilon_0 > 0$ such that $A$ does not have a finite $\epsilon_0$-net.
    Arbitrarily pick $x_1 \in A$; since $\{x_1\}$ cannot be an $\epsilon_0$-net, there must exist $x_2 \in A$ such that $\rho(x_1, x_2) > \epsilon_0$.
    Similarly, $\{x_1, x_2\}$ is not an $\epsilon_0$-net, so there must exist $x_3 \in A$ such that $\rho(x_3, x_j) > \epsilon_0$ for $j=1,2$.
    Continuing this process infinitely, we construct a sequence $\{x_n\} \subset A$ satisfying, for any $m \ne n$:

    \[
    \rho(x_m, x_n) > \epsilon_0
    \]

    Obviously, this sequence does not have any Cauchy subsequence, and naturally no convergent subsequence, which contradicts the pre-compactness of $A$! Thus, $A$ is totally bounded.

    **(2) In a complete space, Totally bounded $\implies$ Pre-compact:**
    Suppose $X$ is complete and $A$ is totally bounded. Take any sequence $\{x_n\}$ in $A$; we want to extract a convergent subsequence from it.
    (If the sequence has only finitely many distinct points, then it must have a constant convergent subsequence; thus, we assume there are infinitely many distinct points, denoted as set $B_0$).

    Because $A$ is totally bounded, its subset $B_0$ is also totally bounded. Cover $B_0$ with a finite number of open balls of radius $1/2$. By the Pigeonhole Principle, at least one ball contains infinitely many elements of $B_0$. Denote the set formed by these elements as $B_1$. Clearly, the diameter of $B_1$ does not exceed 1.
    Next, for the totally bounded set $B_1$, cover it with a finite number of open balls of radius $1/4$, and again we can select a subset $B_2$ containing infinitely many elements, whose diameter does not exceed $1/2$.
    Continuing this process, we obtain a series of nested sets:

    \[
    B_1 \supseteq B_2 \supseteq \dots \supseteq B_k \supseteq \dots
    \]

    Satisfying that each $B_k$ contains infinitely many points from the sequence, and $\text{diam}(B_k) \le 1/2^{k-1}$.
    Now construct the subsequence: take $x_{n_1} \in B_1$; pick $x_{n_2}$ with a larger index ($n_2 > n_1$) from $B_2$; and so on, taking $x_{n_k} \in B_k$ ($n_k > n_{k-1}$).
    For any $m > k$, since $x_{n_m}, x_{n_k} \in B_k$, we have:

    \[
    \rho(x_{n_k}, x_{n_m}) \le \frac{1}{2^{k-1}} \xrightarrow{k \to \infty} 0
    \]

    This indicates that the extracted subsequence $\{x_{n_k}\}$ is a Cauchy sequence. Because the space $X$ is **complete**, this Cauchy sequence must converge in $X$. Therefore, $A$ is pre-compact. $\square$

**Corollaries:**

* A compact set must be **bounded and separable**.

* In a complete metric space, $A$ is pre-compact $\iff$ for any given $\epsilon>0$, $A$ has a pre-compact $\epsilon$-net.

* **A compact set in a metric space is itself a complete metric space** (because its Cauchy sequences must have convergent subsequences, and the limits are within the set itself, thus they must converge).

!!! success "Theorem 3.3 (Nested Compact Set Theorem)"

    Let $K_n$ be a sequence of **non-empty compact sets** in a metric space $X$, satisfying the nested relationship:

    \[
    K_1 \supseteq K_2 \supseteq \cdots \supseteq K_n \supseteq \cdots
    \]

    Then their intersection $\bigcap_{n=1}^\infty K_n$ is non-empty.

??? proof "Proof of the Nested Compact Set Theorem"

    Arbitrarily pick a point $x_n$ from each compact set $K_n$ to form the sequence $\{x_n\}$.
    Since all $x_n$ are contained within the compact set $K_1$, according to the definition of compactness, the sequence $\{x_n\}$ must have a subsequence $\{x_{n_k}\}$ that converges to some point $x_0$ in $K_1$.
    
    For any given integer $n$, when $n_k \ge n$, we have $x_{n_k} \in K_{n_k} \subseteq K_n$.
    That is to say, starting from a certain term, the subsequence entirely falls within $K_n$. Since a compact set must be a closed set, the limit point must also belong to this closed set, meaning $x_0 \in K_n$.
    Because $n$ was chosen arbitrarily, $x_0 \in \bigcap_{n=1}^\infty K_n$. The intersection is non-empty. $\square$

---

## 4. Topological Characterization of Compact Sets: Finite Covers and Finite Intersections

In calculus, a bounded closed interval on the real number axis enjoys the famous Heine-Borel Finite Cover Theorem. In a general metric space, this property not only holds true, but it is also a **perfect equivalent characterization of compact sets**.

!!! info "Definition 3.5 (Open Cover and Finite Subcover)"

    Let $X$ be a metric space, $A \subset X$. Let $\{G_c\}_{c \in J}$ be a family consisting of some open sets in $X$.
    If $A \subseteq \bigcup_{c \in J} G_c$, then $\{G_c\}_{c \in J}$ is called an **open cover** of $A$.
    If the index set $J$ is a finite set, it is called a **finite open cover** of $A$.

*(Note: A bounded closed set in a general metric space does not necessarily satisfy the finite cover theorem. In the previous counterexample in $L^2[-\pi,\pi]$, the open balls $\{S(e_k, 1/2)\}$ centered at the sequence elements are mutually disjoint, making it impossible to select finitely many to cover the closed set.)*

!!! success "Theorem 3.4 (Borel Finite Cover Theorem: Criterion for Compactness)"

    A necessary and sufficient condition for a subset $A$ of a metric space $X$ to be a compact set is: from **any open cover** of $A$, a **finite subcover** can always be selected.

??? proof "Proof of the Finite Cover Theorem (Click to expand)"

    **1. Necessity (Compact set $\implies$ Can select a finite subcover):**
    Suppose $A$ is a compact set and $\{G_c\}_{c \in J}$ is an open cover of $A$.
    We first prove that there exists a so-called "Lebesgue number" $\epsilon_0 > 0$, such that for all $x \in A$, the open ball $S(x, \epsilon_0)$ is strictly contained within some single open set $G_c$.
    Proof by contradiction: If this were not the case, for every natural number $n$, there would exist $x_n \in A$ such that the open ball $S(x_n, 1/2^n)$ is not contained in any $G_c$. Because $A$ is compact, the sequence $\{x_n\}$ has a subsequence $\{x_{n_k}\}$ converging to some point $x_0$ in $A$. Because $\{G_c\}$ covers $A$, there must exist some $G_{c_0}$ such that $x_0 \in G_{c_0}$. As it is an open set, there exists $r > 0$ such that $S(x_0, r) \subset G_{c_0}$.
    When $k$ is sufficiently large, $\rho(x_{n_k}, x_0) < r/2$ and $1/2^{n_k} < r/2$. By the triangle inequality, we then have $S(x_{n_k}, 1/2^{n_k}) \subset S(x_0, r) \subset G_{c_0}$, which contradicts the previous assumption (not being contained in any $G_c$)! Hence, $\epsilon_0$ exists.

    Because $A$ is compact, it must be totally bounded, thus there exist finitely many open balls $\{S(x_j, \epsilon_0)\}_{j=1}^l$ that cover $A$.
    According to the property of $\epsilon_0$, each open ball $S(x_j, \epsilon_0)$ is contained by some open set $G_{c_j}$. Extracting these $l$ corresponding open sets, the family they form $\{G_{c_j}\}_{j=1}^l$ clearly also covers $A$, which is the finite subcover we are looking for.

    **2. Sufficiency (Every open cover has a finite subcover $\implies$ Compact set):**
    Assume the condition holds, and let $\{x_n\}$ be a sequence in $A$. We need to prove it has a convergent subsequence.
    Proof by contradiction: If $\{x_n\}$ has no convergent subsequence in $A$, then for every point $y$ in $A$, there exists a neighborhood $S(y, \delta_y)$ such that only finitely many terms of the sequence fall within this neighborhood (otherwise $y$ would be a limit point, and a convergent subsequence could be extracted).
    Obviously, this family of open balls $\{S(y, \delta_y)\}_{y \in A}$ forms an open cover of $A$.
    According to the theorem's assumption, we can select a finite number of open balls $\{S(y_j, \delta_{y_j})\}_{j=1}^l$ from it to cover $A$.
    However, each selected open ball contains only finitely many terms from the sequence $\{x_n\}$, so their finite union (i.e., the entire set $A$) can also only contain finitely many terms of the sequence! This contradicts the fact that $\{x_n\}$ is an infinite sequence. Thus, a convergent subsequence must exist, and $A$ is a compact set. $\square$

The finite cover property of compact sets can be transformed into a property of closed sets through dual conversion (De Morgan's Laws):

!!! info "Definition 3.6 (Finite Intersection Property)"

    Let $\{F_c\}_{c \in J}$ be a family of sets in a metric space $X$. If **any finite subfamily** of it has a non-empty intersection, then this family of sets is said to have the **Finite Intersection Property**.

!!! success "Theorem 3.5"

    A necessary and sufficient condition for a **closed subset $A$** in a metric space to be a compact set is: every **family of closed subsets** in $A$ that has the finite intersection property has a non-empty intersection.

??? proof "Proof of the Finite Intersection Property (Click to expand)"

    **1. Necessity (Compact set $\implies$ Non-empty intersection):**
    Proof by contradiction: Suppose $A$ is a compact set and $\{F_c\}_{c \in J}$ is a family of closed subsets possessing the finite intersection property, but the intersection of all of them is empty, i.e., $\bigcap_c F_c = \emptyset$.
    Let open sets $G_c = X \setminus F_c$. By De Morgan's Laws:

    \[
    \bigcup_c G_c = \bigcup_c (X \setminus F_c) = X \setminus \bigcap_c F_c = X \setminus \emptyset = X
    \]

    This indicates that $\{G_c\}$ forms an open cover of the entire space $X$, and naturally covers $A$ as well. Due to compactness, there exists a finite subcover $\{G_{c_j}\}_{j=1}^l$ covering $A$.
    Taking the complement again:

    \[
    \bigcap_{j=1}^l F_{c_j} = \bigcap_{j=1}^l (X \setminus G_{c_j}) = X \setminus \bigcup_{j=1}^l G_{c_j} \subseteq X \setminus A
    \]

    However, all $F_c$ are subsets of $A$, so their intersection must be contained within $A$. That is, $\bigcap_{j=1}^l F_{c_j} \subseteq A \cap (X \setminus A) = \emptyset$. This contradicts the premise "possesses the finite intersection property" (the intersection of finitely many is non-empty)!

    **2. Sufficiency (Having non-empty intersection $\implies$ Compact set):**
    Suppose any family of closed subsets possessing the finite intersection property has a non-empty intersection. We want to prove that any open cover of $A$ has a finite subcover.
    Let $\{G_c\}_{c \in J}$ be an open cover of $A$. Define closed subsets $F_c = A \setminus G_c$. Since:

    \[
    \bigcap_{c \in J} F_c = \bigcap_{c \in J} (A \setminus G_c) = A \setminus \bigcup_{c \in J} G_c = \emptyset
    \]

    We know that the family of closed subsets $\{F_c\}_{c \in J}$ **does not** possess the finite intersection property (otherwise the entire intersection would be non-empty according to the premise).
    Since it does not have the finite intersection property, it means there must exist a finite subfamily $\{F_{c_j}\}_{j=1}^l$ such that their intersection is empty:

    \[
    \bigcap_{j=1}^l F_{c_j} = \emptyset
    \]

    Taking the complement of the equation to revert to open sets:

    \[
    \bigcup_{j=1}^l G_{c_j} \supseteq \bigcup_{j=1}^l (A \setminus F_{c_j}) = A \setminus \bigcap_{j=1}^l F_{c_j} = A \setminus \emptyset = A
    \]

    This finds the finite subcover, proving that $A$ is a compact set. $\square$

---

## 5. Continuous Mappings on Compact Sets

Compact sets have excellent "preservation properties" under continuous mappings. Several of the most famous theorems regarding continuous functions on closed intervals in calculus hold true for compact sets in general metric spaces.

!!! success "Theorem 3.6 (The image of a compact set remains a compact set)"

    Let $X, Y$ both be metric spaces, and let closed subset $A \subset X$ be a compact set. If $T$ is a continuous mapping from $X$ to $Y$, then the image of $A$, $T(A)$, is a compact set in $Y$.

??? proof "Proof (Click to expand)"

    Let $\{y_n\}$ be a sequence in $T(A)$, then there must exist a corresponding preimage sequence $\{x_n\}$ in $A$ such that $y_n = T(x_n)$.
    Because $A$ is a compact set, the sequence $\{x_n\}$ must have a subsequence $\{x_{n_k}\}$ that converges to some point $x_0 \in A$ within $A$.
    Also, because the mapping $T$ is continuous, according to the sequential criterion for continuity, limits and mappings can be interchanged:

    \[
    \lim_{k \rightarrow \infty} y_{n_k} = \lim_{k \rightarrow \infty} T(x_{n_k}) = T(x_0)
    \]

    Since $x_0 \in A$, obviously $T(x_0) \in T(A)$.
    This demonstrates that any sequence in $T(A)$ has a subsequence converging within $T(A)$, therefore $T(A)$ is a compact set. $\square$

A continuous mapping on a compact set is not only pointwise continuous but also **Uniformly Continuous**.

!!! success "Theorem 3.7 (Uniform Continuity Theorem)"

    Let $X, Y$ be metric spaces, and compact set $A \subset X$. If $T: A \rightarrow Y$ is continuous, then $T$ is **uniformly continuous** on $A$.

??? proof "Proof of the Uniform Continuity Theorem (Click to expand)"

    **Proof by contradiction**.
    The definition of uniform continuity is: for any given $\epsilon > 0$, there exists a $\delta > 0$ depending only on $\epsilon$, such that when $\rho(x, y) < \delta$, we have $\rho_1(T(x), T(y)) < \epsilon$.
    
    Assume $T$ is not uniformly continuous, which means: there exists a fixed $\epsilon_0 > 0$ such that no matter how small $\delta$ is taken, one can always find points extremely close to each other whose mapped distance is strictly greater than $\epsilon_0$.
    Let $\delta_n = 1/n$, then there exist sequences $\{x_n\}$ and $\{y_n\}$ contained in $A$ such that:

    \[
    \rho(x_n, y_n) < \frac{1}{n} \quad \text{but} \quad \rho_1(T(x_n), T(y_n)) > \epsilon_0
    \]

    Because $A$ is a compact set, the sequence $\{x_n\}$ has a convergent subsequence $\{x_{n_k}\}$, converging to some point $x_0$ within $A$.
    Since $\rho(x_{n_k}, y_{n_k}) < 1/n_k \rightarrow 0$, by the triangle inequality:

    \[
    \rho(y_{n_k}, x_0) \le \rho(y_{n_k}, x_{n_k}) + \rho(x_{n_k}, x_0) \xrightarrow{k \to \infty} 0
    \]

    This indicates that the corresponding subsequence $\{y_{n_k}\}$ also converges to the exact same limit point $x_0$.
    Since $x_{n_k} \rightarrow x_0$ and $y_{n_k} \rightarrow x_0$, utilizing the continuity of $T$ at $x_0$, the image points of both should converge to $T(x_0)$.
    Therefore, by the triangle inequality:

    \[
    \rho_1(T(x_{n_k}), T(y_{n_k})) \le \rho_1(T(x_{n_k}), T(x_0)) + \rho_1(T(y_{n_k}), T(x_0)) \xrightarrow{k \to \infty} 0
    \]

    This requires that when $k$ is sufficiently large, $\rho_1(T(x_{n_k}), T(y_{n_k})) \to 0$. But this creates an irreconcilable contradiction with our initial assumption that $\rho_1(T(x_{n_k}), T(y_{n_k})) > \epsilon_0 > 0$!
    Thus, $T$ must be uniformly continuous. $\square$

Setting the target space $Y$ to be the real number field $\mathbb{R}$, we immediately obtain the classic Extreme Value Theorem.

!!! success "Corollary (Extreme Value Theorem for Continuous Functionals)"

    Let $X$ be a metric space, and $A$ be a compact set in $X$. Let $f$ be a continuous functional (real-valued function) defined on $A$, then $f$ is **bounded**, and it **must achieve its supremum and infimum** on $A$.

??? proof "Derivation of the Extreme Value Theorem"

    According to Theorem 3.6, $f(A)$ is a compact set in the real number field $\mathbb{R}$. In the real field, a compact set is equivalent to a **bounded closed set**.
    Because $f(A)$ is bounded, the functional $f$ is bounded, meaning there exists a finite supremum and infimum.
    Because $f(A)$ is a closed set, it must contain all of its limit points, and naturally contains its supremum and infimum. This means there exist points in $A$ that allow the functional to attain exactly these bounds. $\square$

---

## Summary: Logical Connections of Concepts in Metric Spaces

In general metric spaces, completeness and compactness possess the following elegant web of dualities and equivalences:

* **Core Equivalence Chain**:
  Compact set $\iff$ Complete + Totally bounded $\iff$ Complete + Pre-compact

* **Duality of Topological Characterizations**:
  Every open cover has a finite subcover $\iff$ Every family of closed subsets with the finite intersection property must have a non-empty intersection.