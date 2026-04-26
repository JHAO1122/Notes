# Chapter 2: Topology and Completeness in Metric Spaces

In Chapter 1, we axiomatically generalized the concept of "distance". With distance defined, we can naturally discuss "closeness" and "limits" in abstract spaces. In this chapter, we will further borrow intuitive geometric concepts from calculus to establish a **topological structure** (open sets, closed sets, density) in metric spaces, and explore one of the most core concepts in functional analysis—**Completeness**. Finally, we will prove the famous Baire Category Theorem, which is the cornerstone for studying linear operator theory.

---

## 1. Basic Topological Concepts in Metric Spaces

In Euclidean space, the definitions of limits and continuity rely on "neighborhoods". In general metric spaces, neighborhoods are defined by "balls".

!!! info "Definition 2.1 (Balls and Neighborhoods)"

    Let $(X, \rho)$ be a metric space, $x_0 \in X$, $r > 0$.

    * **Open Ball**: The set $\{x \in X \mid \rho(x, x_0) < r\}$ is called the open ball centered at $x_0$ with radius $r$, denoted by $S(x_0, r)$.
    * **Closed Ball**: The set $\{x \in X \mid \rho(x, x_0) \le r\}$ is called the closed ball, denoted by $\overline{S}(x_0, r)$.
    
    The open ball $S(x_0, r)$ centered at $x_0$ is also called a **(spherical) neighborhood** of $x_0$.

### 1.1 Open and Closed Sets

!!! info "Definition 2.2 (Interior Points and Open Sets)"

    * **Interior Point**: Let $G \subset X$, $x \in G$. If there exists a neighborhood $S(x, r) \subset G$ of $x$, then $x$ is called an interior point of $G$.
    * **Open Set**: If every point in the set $G$ is an interior point of $G$, then $G$ is called an **open set**. The empty set $\emptyset$ is also defined as an open set.

Open sets have the following basic properties:
1. The whole space $X$ and the empty set $\emptyset$ are both open sets;
2. The union of any number of open sets is an open set;
3. The intersection of **finitely many** open sets is an open set.

!!! info "Definition 2.3 (Limit Points, Closure, and Closed Sets)"

    * **Limit Point**: Let $A \subset X$. If for any given $\epsilon > 0$, the neighborhood $S(x_0, \epsilon)$ contains a point of $A$ distinct from $x_0$ (i.e., $S(x_0, \epsilon) \cap (A \setminus \{x_0\}) \ne \emptyset$), then $x_0$ is called a limit point of $A$.
    * **Isolated Point**: If $x_0 \in A$ but is not a limit point, it is called an isolated point.
    * **Derived Set and Closure**: The set of all limit points of $A$ is called the derived set $A'$. The union $\overline{A} = A \cup A'$ is called the **closure** of $A$.
    * **Closed Set**: If $A = \overline{A}$ (i.e., $A$ contains all its limit points), then $A$ is called a **closed set**.

Open and closed sets are dual concepts: **A set $A$ is closed if and only if its complement $A^c$ is open.** Closed sets satisfy: the intersection of any number of closed sets is a closed set; the union of finitely many closed sets is a closed set.

---

## 2. Density and Separability

In mathematical analysis, we often use a simple set to "approximate" a complex set (e.g., using rational numbers to approximate real numbers, or polynomials to approximate continuous functions). This introduces the concept of density.

!!! info "Definition 2.4 (Density)"

    Let $A, B \subset X$. If $\overline{B} \supset A$, then **$B$ is said to be dense in $A$**.

    Density has the following equivalent propositions:
    (i) For any given $x \in A$ and $r > 0$, there exists a point $y \in B$ such that $\rho(x, y) < r$.
    (ii) For any given $r > 0$, the union of all open balls centered at every point in $B$ with radius $r$ contains $A$.
    (iii) For any given $x \in A$, there always exists a sequence $\{x_n\}$ in $B$ such that $x_n \rightarrow x$.

!!! info "Definition 2.5 (Separable Space)"

    If there exists an **at most countable** subset in the metric space $X$ that is **dense** in $X$, then $X$ is said to be **separable**.

* **Example 1:** The real number field $\mathbb{R}$ is separable because the set of rational numbers $\mathbb{Q}$ is countable and dense.
* **Example 2:** The continuous function space $C[a,b]$ is separable. By the Weierstrass Approximation Theorem, the set of polynomials with rational coefficients is dense in $C[a,b]$, and this set is countable.

Not all spaces are separable. A classic non-separable space in functional analysis is the essentially bounded function space $L^\infty[a,b]$.

??? proof "Proof: $L^\infty[a,b]$ is a non-separable metric space"

    We use proof by contradiction.
    Consider the interval $[a,b]$, for any $s \in (a, b]$, define the characteristic function $f_s = \chi_{[a, s]}$.
    Obviously, when $s \ne s'$, these two functions differ on an interval with a measure greater than 0 (one is 1, the other is 0). Therefore, their essential supremum distance in $L^\infty[a,b]$ is:

    \[
    \rho(f_s, f_{s'}) = 1 \quad (s \ne s')
    \]

    Let $A = \{f_s \mid s \in (a,b]\}$ be the set consisting of all such functions. Since the real interval is uncountable, the set $A$ is also **uncountable**.

    Assume that $L^\infty[a,b]$ is separable, then it has a countable dense subset $M_0 = \{g_1, g_2, \dots\}$.
    Taking each element $g_n \in M_0$ as the center and $1/3$ as the radius, we construct open balls $S(g_n, 1/3)$. Since $M_0$ is dense, the union of these open balls covers the entire $L^\infty[a,b]$, and naturally covers the set $A$.

    Because there are only countably many open balls in the family $\{S(g_n, 1/3)\}$, while $A$ has uncountably many elements, by the Pigeonhole Principle, there must exist at least one open ball that contains two distinct elements of $A$, say $f_s$ and $f_{s'}$.
    By the triangle inequality:

    \[
    \rho(f_s, f_{s'}) \le \rho(f_s, g_n) + \rho(g_n, f_{s'}) < \frac{1}{3} + \frac{1}{3} = \frac{2}{3}
    \]

    This contradicts the previously established $\rho(f_s, f_{s'}) = 1$!
    Therefore, $L^\infty[a,b]$ must be non-separable. $\square$

---

## 3. Completeness

The existence of limits is the core of calculus. However, in actual computations, we can often only determine if a sequence is getting "closer" to each other, without knowing if that "limit point" exists within the current space. This is the motivation for introducing Cauchy sequences.

!!! info "Definition 2.6 (Fundamental Sequence / Cauchy Sequence)"

    A sequence $\{x_n\}$ in a metric space $X$ is called a **fundamental sequence (or Cauchy sequence)** if for any given $\epsilon > 0$, there exists $N > 0$ such that for all $m, n > N$, we have:

    \[
    \rho(x_n, x_m) < \epsilon
    \]

It is easy to see that **any convergent sequence must be a fundamental sequence** (derived directly from the triangle inequality). Conversely, whether a fundamental sequence necessarily converges depends on whether the space itself has "holes". For example, a fundamental sequence in the set of rational numbers $\mathbb{Q}$ might converge to an irrational number (like $\sqrt{2}$), in which case the limit is not in $\mathbb{Q}$.

!!! success "Definition 2.7 (Complete Metric Space)"

    If **every fundamental sequence** in a metric space $X$ necessarily converges to a point in $X$, then $X$ is said to be a **complete metric space**.

### 3.1 Completeness of Classical Spaces

Euclidean space $\mathbb{R}^n$ is complete. In functional analysis, we are more concerned with the completeness of function spaces.

??? proof "Proof: The continuous function space $C[a,b]$ is complete"

    Let $\{x_n\}$ be a Cauchy sequence in $C[a,b]$. That is, for any given $\epsilon > 0$, there exists $N$ such that when $m, n > N$, for all $t \in [a,b]$, the following holds:

    \[
    |x_n(t) - x_m(t)| \le \max_{t \in [a,b]} |x_n(t) - x_m(t)| = \rho(x_n, x_m) < \epsilon
    \]

    This shows that for every fixed point $t_0$ on the interval, the sequence of real numbers $\{x_n(t_0)\}$ is a Cauchy sequence in the real field $\mathbb{R}$. By the completeness of real numbers, the sequence must converge to some real number, denoted as $x(t_0)$. Let $t$ traverse $[a,b]$, and we obtain a limit function $x(t)$.

    Next, we need to prove two points: (1) $x_n \to x$ converges under the metric $\rho$ (i.e., uniform convergence); (2) the limit function $x(t) \in C[a,b]$.

    **(1) Uniform Convergence**:
    In $|x_n(t) - x_m(t)| < \epsilon$, let $m \rightarrow \infty$. Since the absolute value function is continuous, we get that when $n > N$, for all $t \in [a,b]$:

    \[
    |x_n(t) - x(t)| \le \epsilon
    \]

    That is, $\rho(x_n, x) \le \epsilon$, which means the sequence of functions $\{x_n(t)\}$ uniformly converges to $x(t)$ on $[a,b]$.

    **(2) Continuity of the Limit Function**:
    According to a classical theorem in mathematical analysis, the uniform limit of a sequence of continuous functions is necessarily a continuous function. Therefore, $x(t)$ is continuous on $[a,b]$, meaning $x \in C[a,b]$.

    In summary, the Cauchy sequence $\{x_n\}$ converges to $x \in C[a,b]$ in $C[a,b]$, thus $C[a,b]$ is complete. $\square$

### 3.2 Completeness of the Measurable Function Space $S$

The space $S$ consists of all almost everywhere finite Lebesgue measurable functions on a measurable set $E$ (assuming measure $m(E) < \infty$). Its metric is defined as:

\[
\rho(x, y) = \int_E \frac{|x(t) - y(t)|}{1 + |x(t) - y(t)|} dt
\]

Through classical bounding methods in measure theory, we know that convergence of a sequence under this metric is equivalent to **Convergence in Measure**.

??? proof "Proof: Completeness of the measurable function space $S$ (Click to expand)"

    **Step 1: Transform into a Cauchy sequence in measure**
    Let $\{x_n\}$ be a Cauchy sequence in $S$. That is, for any given $\epsilon > 0$, as $n, m \to \infty$, $\rho(x_n, x_m) \to 0$. 
    Since convergence under this integral metric is equivalent to convergence in measure, $\{x_n\}$ is also a **Cauchy sequence in measure**. That is, for any given $\sigma > 0$:
    
    \[
    m(\{t \in E \mid |x_n(t) - x_m(t)| \ge \sigma\}) \to 0 \quad (n, m \to \infty)
    \]

    **Step 2: Use Riesz's Theorem to extract an almost everywhere convergent subsequence**
    According to **Riesz's Theorem** in real analysis (a Cauchy sequence in measure must have an almost everywhere convergent subsequence), we can extract a subsequence $\{x_{n_k}\}$ from $\{x_n\}$ such that it converges almost everywhere on $E$ to some measurable function $x(t)$.
    That is, as $k \to \infty$, $x_{n_k}(t) \to x(t)$ a.e. 
    Obviously, by the closedness of measurable functions, the limit function $x \in S$.

    **Step 3: Use the Lebesgue Dominated Convergence Theorem (LDCT) to prove subsequence convergence in metric**
    Consider the integrand in the definition of the metric:
    
    \[
    f_k(t) = \frac{|x_{n_k}(t) - x(t)|}{1 + |x_{n_k}(t) - x(t)|}
    \]
    
    Since $x_{n_k}(t) \to x(t)$ a.e., we have $f_k(t) \to 0$ a.e. 
    At the same time, note that $|f_k(t)| \le 1$ holds for all $t$, and the constant function $1$ is Lebesgue integrable on a set of finite measure $E$.
    According to the **Lebesgue Dominated Convergence Theorem (LDCT)**, the limit and integral can be interchanged:
    
    \[
    \lim_{k \to \infty} \rho(x_{n_k}, x) = \lim_{k \to \infty} \int_E f_k(t) dt = \int_E \lim_{k \to \infty} f_k(t) dt = \int_E 0 \, dt = 0
    \]
    
    This shows that the subsequence $\{x_{n_k}\}$ converges to $x$ in the sense of the metric $\rho$.

    **Step 4: Prove that the original Cauchy sequence also converges to this limit**
    It is known that $\{x_n\}$ is a Cauchy sequence in a metric space, and it has a subsequence converging to $x$.
    Using the triangle inequality:
    
    \[
    \rho(x_n, x) \le \rho(x_n, x_{n_k}) + \rho(x_{n_k}, x)
    \]
    
    For any given $\epsilon > 0$:
    First term: Since $\{x_n\}$ is a Cauchy sequence, when $n, n_k$ are sufficiently large, $\rho(x_n, x_{n_k}) < \epsilon/2$;
    Second term: Since the subsequence converges, when $k$ is sufficiently large, $\rho(x_{n_k}, x) < \epsilon/2$.
    
    Therefore, as $n \to \infty$, $\rho(x_n, x) \to 0$. That is, the original sequence $\{x_n\}$ also converges in distance to the limit function $x$ in space $S$.
    Thus, it is proved that space $S$ satisfies the condition that Cauchy sequences must converge, and is therefore complete. $\square$

### 3.3 Completeness of the Sequence Space $s$

The space $s$ consists of all sequences, with its metric defined as:

\[
\rho(x, y) = \sum_{i=1}^\infty \frac{1}{2^i} \frac{|\xi_i - \eta_i|}{1 + |\xi_i - \eta_i|}
\]

??? proof "Proof: Completeness of space $s$ (Click to expand)"

    **Step 1: Term-by-term convergence**
    Let $\{x_n\}$ be a Cauchy sequence in $s$. For a fixed $k$, when $\rho(x_n, x_m) \to 0$:
    
    \[
    \frac{1}{2^k} \frac{|\xi_k^{(n)} - \xi_k^{(m)}|}{1 + |\xi_k^{(n)} - \xi_k^{(m)}|} \le \rho(x_n, x_m) < \epsilon
    \]
    
    From this, we deduce that $|\xi_k^{(n)} - \xi_k^{(m)}| \to 0$, hence for each $k$, there exists $\xi_k^{(n)} \to \xi_k$. Let $x = (\xi_1, \xi_2, \dots)$.

    **Step 2: Prove convergence in metric**
    For any given $\epsilon > 0$, choose a sufficiently large $K$ such that $\sum_{i=K+1}^\infty \frac{1}{2^i} < \frac{\epsilon}{2}$. 
    For the first $K$ terms, since $\xi_i^{(n)} \to \xi_i$, there exists $N$ such that when $n > N$:
    
    \[
    \sum_{i=1}^K \frac{1}{2^i} \frac{|\xi_i^{(n)} - \xi_i|}{1 + |\xi_i^{(n)} - \xi_i|} < \frac{\epsilon}{2}
    \]
    
    Thus, when $n > N$:
    
    \[
    \rho(x_n, x) = \sum_{i=1}^K (\dots) + \sum_{i=K+1}^\infty (\dots) < \frac{\epsilon}{2} + \frac{\epsilon}{2} = \epsilon
    \]
    
    Therefore, $x_n \to x$. Note that $s$ contains all sequences, so $x$ must be in $s$. $\square$

---

## 4. Nested Closed Sphere Theorem and Baire Category Theorem

This is the most profound and practically valuable topological property in complete metric spaces, perfectly generalizing the "Nested Interval Theorem" from the real number field to infinite-dimensional abstract spaces.

### 4.1 Nested Closed Sphere Theorem

!!! success "Theorem 2.1 (Nested Closed Sphere Theorem)"

    Let $X$ be a complete metric space. Let $K_n = \overline{S}(x_n, r_n)$ be a sequence of closed balls in $X$ satisfying the nested condition:

    \[
    K_1 \supseteq K_2 \supseteq \cdots \supseteq K_n \supseteq \cdots
    \]

    If the radius of the closed balls $r_n \rightarrow 0$ (as $n \to \infty$), then there **exists a unique** point $x_0$ in $X$ that belongs to all the closed balls, i.e., $\bigcap_{n=1}^\infty K_n = \{x_0\}$.

??? proof "Proof of the Nested Closed Sphere Theorem (Click to expand)"

    **1. Existence:**
    
    Consider the sequence of centers of the closed balls $\{x_n\}$. 
    Suppose $m > n$, due to the nested relationship $K_m \subseteq K_n$, the center $x_m$ of the $m$-th ball must be inside the $n$-th ball $K_n$. Therefore:

    \[
    \rho(x_n, x_m) \le r_n
    \]

    It is known that as $n \rightarrow \infty$, $r_n \rightarrow 0$. This indicates that the sequence $\{x_n\}$ is a fundamental sequence (Cauchy sequence). 
    Because the space $X$ is **complete**, $\{x_n\}$ must converge to some point $x_0$ in $X$.

    Next, we prove that $x_0 \in K_n$ holds for all $n$. 
    Pick any fixed integer $n_0$. When $n \ge n_0$, all $x_n$ are contained in the closed ball $K_{n_0}$. Since $K_{n_0}$ is a closed set, the limit of any sequence contained within it must also be contained within it, hence $x_0 \in K_{n_0}$. By the arbitrariness of $n_0$, we know $x_0 \in \bigcap_{n=1}^\infty K_n$.

    **2. Uniqueness:**
    
    Assume there exists another point $y_0 \ne x_0$ that also belongs to all the closed balls. By the triangle inequality, for any $n$:

    \[
    \rho(x_0, y_0) \le \rho(x_0, x_n) + \rho(x_n, y_0) \le r_n + r_n = 2r_n
    \]

    Letting $n \rightarrow \infty$, the right side tends to 0, yielding $\rho(x_0, y_0) = 0$, meaning $x_0 = y_0$. A contradiction. Therefore, the intersection point is unique. $\square$

*(Note: The converse of this theorem also holds, meaning if any sequence of nested closed balls with radii tending to zero has a non-empty intersection, then the space must be complete.)*

### 4.2 Sets of First/Second Category and Baire Category Theorem

To characterize the "thickness" of sets, we introduce the concept of categories:

* **Nowhere Dense Set**: If a subset $A$ is **not dense in any non-empty open set** of $X$, then $A$ is called nowhere dense. Equivalently, for any open ball $S_1$, there must exist a sub-open ball $S_2 \subset S_1$ such that $S_2$ does not intersect with $A$ ($S_2 \cap A = \emptyset$).
* **First Category Set (Meager)**: If $A$ can be represented as the union of **at most countably many** nowhere dense sets, then $A$ is said to be of the first category. It represents a topologically "thin" set.
* **Second Category Set (Non-meager)**: Any set that is not of the first category is called a second category set. It represents a topologically "fat" set.

In functional analysis, the conclusion we care about most is: a complete space can never be pieced together from nowhere dense sets.

!!! success "Theorem 2.2 (Baire Category Theorem)"

    **Any complete metric space $X$ must be a set of the second category.**
    
    *(In other words, a complete space cannot be represented as the union of countably many nowhere dense sets. Or put differently, the intersection of countably many dense open sets remains dense.)*

??? proof "Proof of the Baire Category Theorem (Click to expand)"

    We use **proof by contradiction**. Assume the complete metric space $X$ is a set of the first category, then there exist countably many nowhere dense sets $\{F_n\}$ such that:

    \[
    X = \bigcup_{n=1}^\infty F_n
    \]

    Take an arbitrary non-empty closed ball $\overline{S}(x_0, r_0)$ ($r_0 > 0$).
    
    1. Since $F_1$ is a nowhere dense set, it must not be dense in the interior of $\overline{S}(x_0, r_0)$. Thus, there exists a closed ball $\overline{S}(x_1, r_1) \subset \overline{S}(x_0, r_0)$ such that $\overline{S}(x_1, r_1) \cap F_1 = \emptyset$. To ensure the radii converge to 0 later, we can enforce $r_1 < 1$.
    2. For this new ball $\overline{S}(x_1, r_1)$, since $F_2$ is also a nowhere dense set, similarly there exists a closed ball $\overline{S}(x_2, r_2) \subset \overline{S}(x_1, r_1)$ such that $\overline{S}(x_2, r_2) \cap F_2 = \emptyset$, and we set $r_2 < 1/2$.
    3. Continuing this process by induction, we can construct an infinitely nested sequence of closed balls:

    \[
    \overline{S}(x_1, r_1) \supseteq \overline{S}(x_2, r_2) \supseteq \cdots \supseteq \overline{S}(x_n, r_n) \supseteq \cdots
    \]

    Satisfying two key conditions:
    (i) $r_n < 1/n$, so as $n \rightarrow \infty$, $r_n \rightarrow 0$.
    (ii) The $n$-th closed ball $\overline{S}(x_n, r_n)$ contains absolutely no points from $F_n$.

    Because the space $X$ is **complete**, according to the recently proven **Nested Closed Sphere Theorem**, the intersection of these nested closed balls is non-empty, meaning there must exist a point $y_0 \in X$ satisfying:

    \[
    y_0 \in \bigcap_{n=1}^\infty \overline{S}(x_n, r_n)
    \]

    Since $y_0$ exists in space $X$, and $X = \bigcup_{n=1}^\infty F_n$, $y_0$ must belong to some $F_k$. 
    But this presents a contradiction! According to our construction process, since $y_0$ is inside all closed balls, it is naturally in $\overline{S}(x_k, r_k)$. However, we guaranteed during construction that $\overline{S}(x_k, r_k) \cap F_k = \emptyset$, so it is absolutely impossible for $y_0$ to belong to $F_k$.

    This irreconcilable contradiction proves that the initial assumption was incorrect. Therefore, a complete metric space can never be a set of the first category. $\square$

**Significance of the Baire Category Theorem:** This theorem is the foundation of the "Three Fundamental Theorems" in functional analysis (the Uniform Boundedness Principle, the Open Mapping Theorem, and the Closed Graph Theorem). It tells us that in a complete space, if you require a point to satisfy countably many extremely strict conditions (where each condition excludes a dense open set), then the points satisfying all these conditions simultaneously not only exist, but are also dense!