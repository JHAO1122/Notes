# Homework 02 (Homework 03 & 04)

This section contains exercises related to the third and fourth homework assignments of the Point-Set Topology course.

---

## Problem 1

Let $S$ be a non-empty set.

(i) Consider the discrete metric $d: S \times S \to \mathbb{R}_{\ge 0}$, defined as: $d(x, y)=0$ when $x=y$; and $d(x, y)=1$ when $x \neq y$. Prove that the topology $\mathcal{T}_d$ induced by this metric is exactly the power set $\mathcal{P}(S)$, i.e., every subset of $S$ is an open set in the metric space $(S, d)$.

(ii) Is it always possible to find a metric on $S$ such that the topology it induces is the trivial topology (i.e., $\{ \emptyset, S \}$)?

??? success "Solution"

    **(i) Prove that the discrete metric induces the discrete topology**

    To prove $\mathcal{T}_d = \mathcal{P}(S)$, it is sufficient to prove that every singleton set $\{x\}$ in $S$ is an open set.

    For any $x \in S$, consider the open ball centered at $x$ with radius $\epsilon = 1/2$:

    $$
    B(x, 1/2) = \{ y \in S \mid d(y, x) < 1/2 \}
    $$

    According to the definition of the discrete metric, when $y \neq x$, $d(y, x) = 1 > 1/2$; only when $y = x$, $d(y, x) = 0 < 1/2$.
    Therefore, $B(x, 1/2) = \{x\}$.

    Since an open ball in a metric space must be an open set, $\{x\}$ is an open set.
    Because any subset $A$ of $S$ can be expressed as the union of singleton sets:

    $$
    A = \bigcup_{x \in A} \{x\}
    $$

    and the topology is closed under arbitrary union operations, $A$ must be an open set. Thus, $\mathcal{T}_d = \mathcal{P}(S)$.

    **(ii) Discussion on the metrizability of the trivial topology**

    **Conclusion**: It is possible when $|S| \le 1$; it does not exist when $|S| \ge 2$.

    **Reasoning**:
    Assume $S$ contains at least two distinct points $x$ and $y$. If there exists a metric $d$ inducing the trivial topology, then $d(x, y) = \epsilon > 0$.
    Consider the open ball $B(x, \epsilon) = \{ z \in S \mid d(z, x) < \epsilon \}$.

    Since $d(x, x) = 0 < \epsilon$, then $x \in B(x, \epsilon)$.
    Since $d(y, x) = \epsilon \not< \epsilon$, then $y \notin B(x, \epsilon)$.

    This means $B(x, \epsilon)$ is neither the empty set (it contains $x$) nor the entire space $S$ (it does not contain $y$).
    However, in the trivial topology, the only open sets are $\emptyset$ and $S$. This creates a contradiction.
    Therefore, if $|S| \ge 2$, there is no metric that induces the trivial topology.

---

## Problem 2

Prove that $\mathcal{B} := \{ (a, b) \mid a < b \text{ and } a, b \in \mathbb{Q} \}$ is a basis for the standard topology on $\mathbb{R}$ (induced by the metric $d(x, y) = |x - y|$). Also prove that $\text{Card}(\mathcal{B}) = \aleph_0$.

??? success "Solution"

    **1. Prove that $\mathcal{B}$ is a basis for the standard topology**

    Let $\mathcal{T}$ be the standard topology on $\mathbb{R}$. To prove that $\mathcal{B}$ is a basis, we need to prove that for any open set $U \in \mathcal{T}$ and any point $x \in U$, there exists $B \in \mathcal{B}$ such that $x \in B \subseteq U$.

    Given that $U$ is an open set, there exists $\epsilon > 0$ such that $(x - \epsilon, x + \epsilon) \subseteq U$.
    According to the **Density** of rational numbers in real numbers:
    In the interval $(x - \epsilon, x)$, there must exist a rational number $a$;
    In the interval $(x, x + \epsilon)$, there must exist a rational number $b$.

    Thus we obtain $a < x < b$, and $a, b \in \mathbb{Q}$.
    Clearly $(a, b) \in \mathcal{B}$.
    And because $x - \epsilon < a$ and $b < x + \epsilon$, we have:

    $$
    x \in (a, b) \subseteq (x - \epsilon, x + \epsilon) \subseteq U
    $$

    Proof complete.

    **2. Calculate the cardinality of the basis**

    Each element of the set $\mathcal{B}$ is uniquely determined by an ordered pair of rational numbers $(a, b)$ satisfying $a < b$.
    Therefore, there exists an injection $\Psi: \mathcal{B} \to \mathbb{Q} \times \mathbb{Q}$.
    Since $\mathbb{Q}$ is a countable set, its Cartesian product $\mathbb{Q} \times \mathbb{Q}$ is also a countable set (i.e., its cardinality is $\aleph_0$).
    Because $\mathcal{B}$ is an infinite set (e.g., containing all $(0, n), n \in \mathbb{Q}^+$), and its cardinality is not greater than $\aleph_0$, by the Cantor-Bernstein Theorem:

    $$
    \text{Card}(\mathcal{B}) = \aleph_0
    $$

---

## Problem 3

For $n \in \mathbb{N}_{\ge 1}$, consider a finite set $X_n := \{1, 2, \dots, n\}$ of size $n$.
How many topologies are there on $X_n$? Furthermore, consider the set $\Omega$ of all topologies on $X_n$. Define $\mathcal{T}_i < \mathcal{T}_j$ as $\mathcal{T}_i \subsetneq \mathcal{T}_j$. Under what conditions does $(\Omega, <)$ form a linear order relation?

??? success "Solution"

    **1. Number of Topologies**

    There is no simple general formula for the number of topologies on a finite set, usually denoted as $T(n)$. For the first few $n$, the values are:

    * $n=1: T(1) = 1$ (only the trivial topology)
    * $n=2: T(2) = 4$ (trivial, discrete, and two Sierpinski topologies)
    * $n=3: T(3) = 29$

    **2. Discussion on the Linear Order Relation**

    **Conclusion**: $(\Omega, <)$ forms a linear order relation if and only if $n=1$.

    **Reasoning**:
    A linear order requires that for any two elements $\mathcal{T}_i, \mathcal{T}_j$ in $\Omega$, it must hold that $\mathcal{T}_i \subseteq \mathcal{T}_j$ or $\mathcal{T}_j \subseteq \mathcal{T}_i$.

    When $n \ge 2$, we can construct two incomparable topologies. Suppose $X_n$ contains distinct elements $a$ and $b$.
    Define topology $\mathcal{T}_1 = \{ \emptyset, X_n, \{a\} \}$;
    Define topology $\mathcal{T}_2 = \{ \emptyset, X_n, \{b\} \}$.

    Clearly $\{a\} \notin \mathcal{T}_2$ and $\{b\} \notin \mathcal{T}_1$.
    Thus $\mathcal{T}_1 \not\subseteq \mathcal{T}_2$ and $\mathcal{T}_2 \not\subseteq \mathcal{T}_1$.
    Therefore, when $n \ge 2$, $(\Omega, <)$ is not a linear order.

---

## Problem 4

Let $S$ be a set.

(i) Prove: If $\{ \mathcal{T}_i \}_{i \in I}$ is a family of topologies on $S$, then their intersection $\bigcap_{i \in I} \mathcal{T}_i$ is also a topology on $S$.
(ii) Give an example of two topologies $\mathcal{T}_1$ and $\mathcal{T}_2$ such that $\mathcal{T}_1 \cup \mathcal{T}_2$ is not a topology. Find the unique smallest topology containing $\mathcal{T}_1 \cup \mathcal{T}_2$.

??? success "Solution"

    **(i) Prove that the intersection is a topology**

    Denote $\mathcal{T} = \bigcap_{i \in I} \mathcal{T}_i$.

    1. **Contains the empty set and the whole space**: Since for all $i \in I$, $\emptyset, S \in \mathcal{T}_i$ holds, then $\emptyset, S \in \mathcal{T}$.
    2. **Closed under arbitrary unions**: Let $\{ U_\alpha \}_{\alpha \in A} \subseteq \mathcal{T}$. Then for each $i \in I$, $\{ U_\alpha \} \subseteq \mathcal{T}_i$. Since $\mathcal{T}_i$ is a topology, $\bigcup U_\alpha \in \mathcal{T}_i$. Therefore $\bigcup U_\alpha \in \bigcap \mathcal{T}_i = \mathcal{T}$.
    3. **Closed under finite intersections**: Let $U_1, \dots, U_k \in \mathcal{T}$. Then for each $i \in I$, $U_1, \dots, U_k \in \mathcal{T}_i$, so $\bigcap_{j=1}^k U_j \in \mathcal{T}_i$. Therefore, the finite intersection belongs to $\mathcal{T}$.

    **(ii) Counterexample where the union is not a topology**

    Let $S = \{a, b, c\}$.
    Take $\mathcal{T}_1 = \{ \emptyset, S, \{a\} \}$ and $\mathcal{T}_2 = \{ \emptyset, S, \{b\} \}$.
    Then $\mathcal{T}_1 \cup \mathcal{T}_2 = \{ \emptyset, S, \{a\}, \{b\} \}$.

    However, $\{a\} \cup \{b\} = \{a, b\}$, but $\{a, b\} \notin \mathcal{T}_1 \cup \mathcal{T}_2$.
    This violates the axiom of closure under union operations, so the union is not a topology.

    **Smallest Topology**:
    The unique smallest topology containing $\mathcal{T}_1 \cup \mathcal{T}_2$ is the topology generated by this union as a **subbasis**. In this example, the topology is:

    $$
    \mathcal{T}_{min} = \{ \emptyset, S, \{a\}, \{b\}, \{a, b\} \}
    $$

---

## Problem 5

A subset of $\mathbb{R}^2$ is called **radially open** if it contains an open line segment in every direction at each of its points. Let $\mathcal{T}_{ro}$ denote the family of all radially open sets.

Prove that $\mathcal{T}_{ro}$ is a topology on $\mathbb{R}^2$. What is its relationship with the standard topology $\mathcal{T}_{std}$?

??? success "Solution"

    **1. Prove it is a topology**

    The key lies in proving closure under finite intersections.
    Let $A, B \in \mathcal{T}_{ro}$. For any point $x \in A \cap B$:

    * Since $x \in A$, there exists an open line segment $I_A(\theta) \subseteq A$ in direction $\theta$;
    * Since $x \in B$, there exists an open line segment $I_B(\theta) \subseteq B$ in direction $\theta$.

    Take the intersection of these two segments $I_A \cap I_B$. Since they share the same center and are on the same ray, their intersection is still an open line segment containing $x$, and $I_A \cap I_B \subseteq A \cap B$.
    Therefore $A \cap B$ is also radially open. Combining with the obvious inclusion of $\emptyset, S$ and the properties of union operations, $\mathcal{T}_{ro}$ is a topology.

    **2. Relationship with the standard topology**

    **Conclusion**: $\mathcal{T}_{std} \subsetneq \mathcal{T}_{ro}$.

    * **Inclusion relation**: Points in a standard open set contain an open ball $B(x, \epsilon)$. This ball clearly contains an open line segment of length $2\epsilon$ in any direction. Thus, a standard open set must be a radially open set.
    * **Strictness (Counterexample)**: There exist radially open sets that are not standard open sets.
      Construct a set: $U = \{ (r\cos\theta, r\sin\theta) \mid 0 \le r < f(\theta) \}$, where $f(\theta)$ is a positive function that approaches 0 rapidly as the polar angle approaches certain values.
      For example, if $f(\theta)$ is made extremely "sharp" near the origin, the origin contains segments in every direction but cannot contain any full disk. Such a set is radially open at the origin, but the origin is not an interior point under the standard topology.

---

## Problem 6

Let $X = \mathbb{R}^2$ endowed with the standard topology. Consider its subset (Topologist's Sine Curve):

$$
A = \{ (x, \sin(1/x)) : x \in (0, 1] \} \subset X
$$

Calculate the following four sets: $A^o$, $LP(A)$, $\overline{A}$, and $\overline{A} \cap \overline{\mathbb{R}^2 - A}$.

??? success "Solution"

    **1. Interior $A^o$**

    Since $A$ is the graph of the continuous function $y = \sin(1/x)$ on the interval $(0, 1]$, it is a curve in $\mathbb{R}^2$.
    For any point $P$ in $A$, any of its open neighborhoods (open disks) will contain points not belonging to the curve (i.e., points whose y-coordinates do not satisfy $y = \sin(1/x)$).
    Therefore, $A$ does not contain any non-empty open balls.
    Thus:

    $$
    A^o = \emptyset
    $$

    **2. Limit Point Set $LP(A)$ and Closure $\overline{A}$**

    As $x \to 0^+$, the function $\sin(1/x)$ oscillates infinitely many times between $-1$ and $1$.
    Therefore, every point in the segment $\{0\} \times [-1, 1]$ on the y-axis is a limit point of $A$.
    The points of set $A$ itself are obviously also its limit points.
    Therefore:

    $$
    LP(A) = A \cup ( \{0\} \times [-1, 1] )
    $$

    Since the closure $\overline{A} = A \cup LP(A)$, we can conclude:

    $$
    \overline{A} = A \cup ( \{0\} \times [-1, 1] )
    $$

    **3. Intersection $\overline{A} \cap \overline{\mathbb{R}^2 - A}$**

    Since $A^o = \emptyset$, every point in the whole space belongs either to the interior of $\mathbb{R}^2 - A$ or to the boundary of $A$.
    Since $\overline{\mathbb{R}^2 - A} = \mathbb{R}^2 - A^o = \mathbb{R}^2 - \emptyset = \mathbb{R}^2$.
    Therefore:

    $$
    \overline{A} \cap \overline{\mathbb{R}^2 - A} = \overline{A} \cap \mathbb{R}^2 = \overline{A}
    $$

---

## Problem 7

Let $X$ be a set, and $f: \mathcal{P}(X) \to \mathcal{P}(X)$ be a function satisfying the following properties (for any $A, B \subseteq X$):
(i) $f(\emptyset) = \emptyset$;
(ii) $A \subseteq f(A)$;
(iii) $f(A) = f(f(A))$;
(iv) $f(A \cup B) = f(A) \cup f(B)$.

Prove:
(1) The set $\mathcal{T}_f := \{ X - A : A \in \mathcal{P}(X) \text{ and } f(A) = A \}$ is a topology on $X$;
(2) The closure $\overline{A}$ of any subset $A \subseteq X$ with respect to the topology $\mathcal{T}_f$ is exactly $f(A)$.

??? success "Solution"

    **(1) Prove that $\mathcal{T}_f$ is a topology**

    * **Empty set and Whole space**:
      From property (i), $f(\emptyset) = \emptyset$, thus $X - \emptyset = X \in \mathcal{T}_f$.
      From property (ii), $X \subseteq f(X)$, while $f(X) \subseteq X$ obviously holds, so $f(X) = X$.
      Therefore $X - X = \emptyset \in \mathcal{T}_f$.

    * **Arbitrary Union Operation**:
      Let $\{ U_\alpha \}_{\alpha \in I}$ be a family of open sets in $\mathcal{T}_f$, where $U_\alpha = X - A_\alpha$ and $f(A_\alpha) = A_\alpha$.
      Their union $\bigcup U_\alpha = X - \bigcap A_\alpha$.
      From property (iv), it can be deduced that $f$ is monotonic: if $A \subseteq B$, then $f(A) \subseteq f(B)$.
      Since $\bigcap A_\alpha \subseteq A_\alpha$, then $f(\bigcap A_\alpha) \subseteq f(A_\alpha) = A_\alpha$.
      Thus $f(\bigcap A_\alpha) \subseteq \bigcap A_\alpha$.
      Combining with the reverse inclusion in property (ii), we get $f(\bigcap A_\alpha) = \bigcap A_\alpha$.
      Therefore $\bigcup U_\alpha \in \mathcal{T}_f$.

    * **Finite Intersection Operation**:
      Let $U_1, U_2 \in \mathcal{T}_f$, then $U_1 \cap U_2 = X - (A_1 \cup A_2)$.
      From property (iv), $f(A_1 \cup A_2) = f(A_1) \cup f(A_2) = A_1 \cup A_2$.
      Therefore $U_1 \cap U_2 \in \mathcal{T}_f$.

    **(2) Prove the closure $\overline{A} = f(A)$**

    According to the definition of closure, $\overline{A}$ is the intersection of all closed sets containing $A$. Under this topology, a closed set $C$ satisfies $f(C) = C$.

    * First, from property (iii), $f(f(A)) = f(A)$, so $f(A)$ itself is a closed set. Also, from property (ii), $A \subseteq f(A)$.
      Thus $f(A)$ is a closed set containing $A$, so $\overline{A} \subseteq f(A)$.

    * Second, let $C$ be any closed set containing $A$ (i.e., $A \subseteq C$ and $f(C) = C$).
      From the monotonicity of $f$, we have $f(A) \subseteq f(C) = C$.
      This means $f(A)$ is contained in every closed set containing $A$, therefore $f(A) \subseteq \bigcap C = \overline{A}$.

    In summary, $\overline{A} = f(A)$.

---

## Problem 8

Consider the function $f: [0, 1) \to S^1$, defined as $f(t) = (\cos(2\pi t), \sin(2\pi t))$.
Prove that $f$ is not a homeomorphism with respect to the standard subspace topology. Is it possible to find a homeomorphism between $[0, 1)$ and $S^1$?

??? success "Solution"

    **1. Prove $f$ is not a homeomorphism**

    Although $f$ is a continuous bijection, its inverse $f^{-1}$ is discontinuous at point $(1, 0)$.
    We can provide a more rigorous proof from the perspective of topological invariants:

    * **Method 1: Compactness**
      Under the standard Euclidean topology, $S^1$ is a bounded closed set in $\mathbb{R}^2$, and therefore it is compact.
      However, $[0, 1)$ is not a closed set in $\mathbb{R}$, and therefore it is not compact.
      Since a homeomorphism must preserve compactness, $f$ cannot be a homeomorphism.

    * **Method 2: Cut point (Connectedness)**
      If $f$ is a homeomorphism, then the connectedness of the space after removing a point should remain unchanged.
      For $S^1$, removing any point $p$ leaves the remaining space $S^1 - \{p\}$ still connected (homeomorphic to an open interval).
      For $[0, 1)$, if we remove point $t = 1/2$, the remaining space $[0, 1/2) \cup (1/2, 1)$ is disconnected.
      Therefore, these two spaces do not have the same local connected structure, so $f$ is not a homeomorphism.

    **2. Conclusion**

    **Impossible**.
    Based on the reasons of compactness and cut points above, no homeomorphism exists between $[0, 1)$ and $S^1$.

---

## Problem 9

A subset $A \subseteq X$ is called **regularly open** if $A = (\overline{A})^o$.

(i) Give an example of an open set in $\mathbb{R}$ that is not regularly open. Can you characterize regularly open sets in $\mathbb{R}$?
(ii) Prove that for any $A \subseteq X$, the set $(\overline{A})^o$ is a regularly open set.

??? success "Solution"

    **(i) Example and Characterization**

    * **Counterexample**: Take $A = (0, 1) \cup (1, 2)$.
      Then $\overline{A} = [0, 2]$, and $(\overline{A})^o = (0, 2)$.
      Clearly $A \neq (\overline{A})^o$, so $A$ is not a regularly open set. The reason is that $A$ is missing a point $\{1\}$ from the interior of the closure.

    * **Characterization**: Regularly open sets in $\mathbb{R}$ are those open sets "without extra gaps". That is, they cannot be obtained by removing some isolated points or boundary points from a larger open set.

    **(ii) Prove $(\overline{A})^o$ is regularly open**

    Let $U = $(\overline{A})^o$. We want to prove $U = (\overline{U})^o$.

    * According to the definition $U \subseteq \overline{U}$, since $U$ is an open set, $U = U^o \subseteq (\overline{U})^o$.

    * On the other hand, since $U = (\overline{A})^o \subseteq \overline{A}$, and the closure operation is monotonic and idempotent:
      It follows that $\overline{U} \subseteq \overline{(\overline{A})} = \overline{A}$.
      Taking the interior of both sides gives $(\overline{U})^o \subseteq (\overline{A})^o = U$.

    In summary, $U = (\overline{U})^o$, i.e., $(\overline{A})^o$ is a regularly open set.

---

## Problem 10

Let $A$ be a subset of a topological space $X$. Prove the equation $\overline{A} = A \cup LP(A)$.

??? success "Solution"

    **1. Prove $A \cup LP(A) \subseteq \overline{A}$**

    * It is obvious that $A \subseteq \overline{A}$.
    
    * Let $x \in LP(A)$. According to the definition of limit point, for any open neighborhood $U$ of $x$, $U \cap (A - \{x\}) \neq \emptyset$.
      This means $U \cap A \neq \emptyset$ holds for all open neighborhoods.
      According to the characterization of closure by neighborhoods, this is exactly the necessary and sufficient condition for $x \in \overline{A}$.
      Thus $LP(A) \subseteq \overline{A}$.
    
    Therefore $A \cup LP(A) \subseteq \overline{A}$.

    **2. Prove $\overline{A} \subseteq A \cup LP(A)$**

    Let $x \in \overline{A}$.

    * If $x \in A$, then the conclusion obviously holds.
    
    * If $x \notin A$, since $x \in \overline{A}$, for any open neighborhood $U$ of $x$, $U \cap A \neq \emptyset$.
      Because $x \notin A$, the points in the set $U \cap A$ must not be $x$.
      That is, $U \cap (A - \{x\}) \neq \emptyset$.
      This satisfies the definition of a limit point, so $x \in LP(A)$.
    
    Therefore $x \in A \cup LP(A)$.

    In summary, $\overline{A} = A \cup LP(A)$.