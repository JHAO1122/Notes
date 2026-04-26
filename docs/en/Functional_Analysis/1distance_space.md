# Chapter 1: Metric Spaces

In mathematical analysis, we are familiar with limits and continuity in the real number system $\mathbb{R}$ and Euclidean space $\mathbb{R}^n$. The essence of these concepts relies on the notion of "distance." The first step in functional analysis is to abstract this intuitive concept and generalize it to sets containing more complex elements, such as functions and sequences. This chapter introduces the basic axioms of metric spaces and several specific spaces and classical inequalities that are core to functional analysis.

---

## 1. Definition of Metric Spaces

!!! info "Definition 1.1 (Metric Space)"

    Let $X$ be a non-empty set. If for any two elements $x, y \in X$, there is a corresponding certain real number, denoted by $\rho(x, y)$, which satisfies the following three conditions (Metric Axioms):

    **1. Non-negativity and Positivity:**
    
    \[
    \rho(x, y) \ge 0
    \]
    
    and $\rho(x, y) = 0$ if and only if $x = y$.

    **2. Symmetry:**
    
    \[
    \rho(x, y) = \rho(y, x)
    \]

    **3. Triangle Inequality:**
    
    \[
    \rho(x, y) \le \rho(x, z) + \rho(z, y), \quad \forall z \in X
    \]

    Then $\rho$ is called a **metric** (or distance) on the set $X$, and $X$ is called a **metric space** with respect to $\rho$, denoted as $(X, \rho)$. The elements in a metric space are usually called **points**.

**Subspaces and Discrete Spaces:**

* **Subspace**: If $A$ is a non-empty subset of $X$, then $A$ itself forms a metric space under the same metric $\rho$, called a subspace of $X$.
* **Discrete Metric Space**: On any non-empty set $X$, a trivial metric can be defined: $\rho(x, y) = 0$ if $x = y$, and $\rho(x, y) = 1$ if $x \ne y$. This constitutes a discrete metric space.

---

## 2. Classical Examples in Functional Analysis

### 2.1 Euclidean Spaces $\mathbb{R}^n$ and $\mathbb{C}^n$

For the $n$-dimensional Euclidean space $\mathbb{R}^n$, for any two points $x = (x_1, \dots, x_n)$ and $y = (y_1, \dots, y_n)$, the standard distance is defined as:

\[
\rho(x, y) = \left( \sum_{j=1}^n (x_j - y_j)^2 \right)^{\frac{1}{2}}
\]

The proof of the triangle inequality for this metric relies directly on the **Cauchy Inequality** in classical algebra.
Complex spaces $\mathbb{C}^n$ can define a similar Euclidean metric, or the maximum metric $\rho_\infty(x, y) = \max_j |x_j - y_j|$.

### 2.2 Continuous Function Space $C[a,b]$

The set of all real (or complex) continuous functions on the closed interval $[a,b]$. For any $x, y \in C[a,b]$, the distance is defined as the maximum uniform deviation:

\[
\rho(x, y) = \max_{t \in [a,b]} |x(t) - y(t)|
\]

??? proof "Proof: $C[a,b]$ satisfies the triangle inequality"

    Let $x, y, z \in C[a,b]$. For any fixed $t \in [a,b]$, according to the triangle inequality for absolute values of real numbers:

    \[
    |x(t) - y(t)| \le |x(t) - z(t)| + |z(t) - y(t)|
    \]

    Bound each term on the right by its maximum value over the interval:

    \[
    |x(t) - y(t)| \le \max_{t \in [a,b]} |x(t) - z(t)| + \max_{t \in [a,b]} |z(t) - y(t)| = \rho(x, z) + \rho(z, y)
    \]

    Since the right-hand side is a constant independent of $t$ and holds for all $t \in [a,b]$, it holds for the maximum of the left-hand side:

    \[
    \max_{t \in [a,b]} |x(t) - y(t)| \le \rho(x, z) + \rho(z, y)
    \]

    Thus, $\rho(x, y) \le \rho(x, z) + \rho(z, y)$. $\square$

### 2.3 Sequence Spaces $l^p$ and $l^\infty$

* **Space $l^p$ ($1 \le p < \infty$)**: Consists of all sequences $x = \{x_1, x_2, \dots\}$ satisfying $\sum_{j=1}^\infty |x_j|^p < \infty$. Its distance is defined as:

    \[
    \rho(x, y) = \left( \sum_{j=1}^\infty |x_j - y_j|^p \right)^{\frac{1}{p}}
    \]

* **Space $l^\infty$**: Consists of all **bounded** sequences, with the distance defined as the supremum:

    \[
    \rho(x, y) = \sup_{j} |x_j - y_j|
    \]

To prove that the $l^p$ space satisfies the triangle inequality, we must introduce two of the most fundamental inequalities in functional analysis: the Hölder inequality and the Minkowski inequality.

---

## 3. Classical Inequalities and Proofs (Hölder & Minkowski)

!!! success "Theorem 3.1 (Hölder's Inequality)"

    Let $p > 1, q > 1$ satisfy the conjugate exponent relationship $\frac{1}{p} + \frac{1}{q} = 1$.
    For any sequences $a = \{a_1, a_2, \dots\} \in l^p$ and $b = \{b_1, b_2, \dots\} \in l^q$, we have:

    \[
    \sum_{j=1}^\infty |a_j b_j| \le \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} \cdot \left( \sum_{j=1}^\infty |b_j|^q \right)^{\frac{1}{q}}
    \]

??? proof "Proof of Hölder's Inequality (Click to expand)"

    **Step 1: Introduce Young's Inequality**
    
    Consider the function $f(x) = x^\alpha$ (where $0 < \alpha < 1, x \ge 0$). Since $f''(x) = \alpha(\alpha-1)x^{\alpha-2} < 0$, this function is concave. Its tangent line equation at point $(1,1)$ is $y = \alpha(x-1) + 1 = \alpha x + \beta$ (where $\beta = 1 - \alpha$).
    By concavity, the function image lies below the tangent line:

    \[
    x^\alpha \le \alpha x + \beta
    \]

    Let $x = \frac{u}{v} \; (u, v > 0)$, substitute into the above and multiply by $v$ to obtain **Young's Inequality**:

    \[
    u^\alpha v^\beta \le \alpha u + \beta v \quad (\alpha + \beta = 1)
    \]

    **Step 2: Construct Unitized Variables**
    
    Let $\alpha = \frac{1}{p}, \beta = \frac{1}{q}$. For any given $j$, let:

    \[
    u = \frac{|a_j|^p}{\sum_{j=1}^\infty |a_j|^p}, \quad v = \frac{|b_j|^q}{\sum_{j=1}^\infty |b_j|^q}
    \]

    **Step 3: Apply Young's Inequality and Sum**
    
    Substituting into Young's Inequality, we get:

    \[
    \frac{|a_j b_j|}{ \left( \sum |a_j|^p \right)^{\frac{1}{p}} \left( \sum |b_j|^q \right)^{\frac{1}{q}} } \le \frac{1}{p} \frac{|a_j|^p}{\sum |a_j|^p} + \frac{1}{q} \frac{|b_j|^q}{\sum |b_j|^q}
    \]

    Summing both sides over $j$ from $1$ to $\infty$:

    \[
    \frac{\sum_{j=1}^\infty |a_j b_j|}{ \left( \sum |a_j|^p \right)^{\frac{1}{p}} \left( \sum |b_j|^q \right)^{\frac{1}{q}} } \le \frac{1}{p} \cdot 1 + \frac{1}{q} \cdot 1 = 1
    \]

    Rearranging yields Hölder's Inequality. $\square$

!!! success "Theorem 3.2 (Minkowski's Inequality)"

    For $p \ge 1$, for any two elements $a$ and $b$ in $l^p$, then $a+b \in l^p$, and:

    \[
    \left( \sum_{j=1}^\infty |a_j + b_j|^p \right)^{\frac{1}{p}} \le \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} + \left( \sum_{j=1}^\infty |b_j|^p \right)^{\frac{1}{p}}
    \]

    *Note: This is essentially the core proof that $l^p$ spaces satisfy the triangle inequality.*

??? proof "Proof of Minkowski's Inequality (Click to expand)"

    **Step 1: Prove Closedness ($a+b \in l^p$)**
    
    Using the elementary inequality $|a_j + b_j|^p \le (2 \max\{|a_j|, |b_j|\})^p \le 2^p (|a_j|^p + |b_j|^p)$, summing this shows that as long as $a, b \in l^p$, the sum sequence must also converge absolutely, i.e., $a+b \in l^p$.

    **Step 2: Decomposition and Hölder's Inequality**
    
    When $p > 1$ (the case $p=1$ follows directly from the absolute value triangle inequality), we use the conjugate relationship $(p-1)q = p$ for algebraic manipulation:

    \[
    \sum_{j=1}^\infty |a_j + b_j|^p = \sum_{j=1}^\infty |a_j + b_j| \cdot |a_j + b_j|^{p-1} \le \sum_{j=1}^\infty (|a_j| + |b_j|) |a_j + b_j|^{p-1}
    \]

    After expansion, apply **Hölder's Inequality** to these two terms separately:

    \[
    \sum_{j=1}^\infty |a_j| |a_j + b_j|^{p-1} \le \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} \left( \sum_{j=1}^\infty |a_j + b_j|^{(p-1)q} \right)^{\frac{1}{q}}
    \]

    Since $(p-1)q = p$, the second factor on the right is $\left( \sum |a_j + b_j|^p \right)^{1/q}$. Do the same for the $|b_j|$ term.

    **Step 3: Combine Like Terms**
    
    Adding the two terms and factoring:

    \[
    \sum_{j=1}^\infty |a_j + b_j|^p \le \left[ \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} + \left( \sum_{j=1}^\infty |b_j|^p \right)^{\frac{1}{p}} \right] \left( \sum_{j=1}^\infty |a_j + b_j|^p \right)^{\frac{1}{q}}
    \]

    Dividing both sides by $\left( \sum |a_j + b_j|^p \right)^{\frac{1}{q}}$, and noting $1 - \frac{1}{q} = \frac{1}{p}$:

    \[
    \left( \sum_{j=1}^\infty |a_j + b_j|^p \right)^{\frac{1}{p}} \le \left( \sum_{j=1}^\infty |a_j|^p \right)^{\frac{1}{p}} + \left( \sum_{j=1}^\infty |b_j|^p \right)^{\frac{1}{p}}
    \]
    
    Proof complete. $\square$

*(Note: For measurable function spaces $L^p(F)$, replace the summation $\sum$ above with the Lebesgue integral $\int_F$ over the measurable set $F$ to obtain the corresponding integral forms of Hölder and Minkowski inequalities.)*

---

## 4. Convergence in Metric Spaces

After establishing the metric standard, limit operations can be naturally defined.

!!! info "Definition 4.1 (Convergence of a Sequence)"

    Let $\{x_n\}$ be a sequence in a metric space $(X, \rho)$. If there exists a point $y \in X$ such that as $n \rightarrow \infty$:

    \[
    \rho(x_n, y) \rightarrow 0
    \]

    Then $\{x_n\}$ is said to converge to $y$, denoted as $x_n \rightarrow y$ or $\lim_{n \rightarrow \infty} x_n = y$. $y$ is called the limit of the sequence.

### 4.1 Basic Properties of Convergent Sequences

For a convergent sequence $\{x_n\}$ in a metric space, it naturally inherits several excellent properties of limits on the real number line:

1. **Uniqueness of Limits**: If $x_n \rightarrow y$ and $x_n \rightarrow y'$, then it must be that $y = y'$.
2. **Boundedness**: For any given point $y_0 \in X$, the real sequence $\rho(x_n, y_0)$ is bounded.
3. **Subsequence Inheritance**: Any subsequence of $\{x_n\}$ also converges to the same limit. Conversely, if every subsequence converges, then the original sequence must converge.

??? proof "Proof: Uniqueness and Boundedness of Limits"

    **1. Proof of Uniqueness:**
    
    Using the triangle inequality:

    \[
    0 \le \rho(y, y') \le \rho(y, x_n) + \rho(x_n, y') = \rho(x_n, y) + \rho(x_n, y')
    \]

    Since $x_n \rightarrow y$ and $x_n \rightarrow y'$, the right-hand side tends to 0 as $n \rightarrow \infty$. By the squeeze theorem, $\rho(y, y') = 0$, and according to the first metric axiom, $y = y'$.

    **2. Proof of Boundedness:**
    
    Let $x_n \rightarrow y$, meaning $\rho(x_n, y) \rightarrow 0$. By the definition of convergence, there must exist $M > 0$ such that $\rho(x_n, y) < M$ for all $n$.
    Take any $y_0 \in X$, by the triangle inequality:

    \[
    \rho(x_n, y_0) \le \rho(x_n, y) + \rho(y, y_0) < M + \rho(y, y_0)
    \]

    The right-hand side is a constant independent of $n$, so the real sequence $\rho(x_n, y_0)$ is bounded. $\square$

### 4.2 Analysis of Convergence Equivalence in Different Spaces

The subtlety of functional analysis lies in: **The elements of a space may be the same, but different metrics lead to entirely different concepts of "convergence."**

* **In finite-dimensional Euclidean space $\mathbb{R}^n$**: Convergence of a sequence under the Euclidean metric is **necessary and sufficient** for each coordinate (component) to converge separately in the real field. Here, metric space convergence reduces to coordinate-wise convergence.

* **In the space of continuous functions $C[a,b]$**:
  * If assigned the maximum metric $\rho(x, y) = \max_{t \in [a,b]} |x(t) - y(t)|$, sequence convergence is equivalent to **Uniform Convergence** of the function sequence on $[a,b]$.
  * **Counter-example**: If assigned the mean-square error metric $\rho_1(x, y) = \left( \int_a^b |x(t) - y(t)|^2 dt \right)^{1/2}$ in $C[a,b]$. Consider the sequence $x_n(t) = t^n$ on $[0,1]$.
    Clearly $\rho_1(x_n, 0) \rightarrow 0$ (converges to 0 under integral metric), but as a function sequence, it does not converge uniformly to 0 on $[0,1]$ (it is always 1 at $t=1$).

This shows that the specific definition of a metric profoundly affects the topological structure and limit behavior within a space, which is also an important motivation for our subsequent study of topological properties and completeness.