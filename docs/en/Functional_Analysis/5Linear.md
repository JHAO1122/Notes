# Chapter 5: Fundamentals of Normed Linear Spaces

Normed linear spaces are the fundamental objects of study in functional analysis. They possess both the algebraic structure of linear algebra (linear operations) and a topological structure (distance and convergence) established through the introduction of a "norm." This chapter begins with linear spaces, transitions to normed linear spaces, and introduces several classical function spaces.

---

## 1. Linear Space

!!! info "Definition 5.1 (Linear Space)"

    Let $E$ be a non-empty set and $K$ be the field of real numbers $\mathbb{R}$ or complex numbers $\mathbb{C}$. $E$ is called a real (or complex) linear space if it satisfies the following properties:

    **(i) Additive Group Properties**: For any $x, y, z \in E$, there corresponds a unique sum $x+y \in E$ such that:

    * (a) Commutativity: $x+y = y+x$;

    * (b) Associativity: $(x+y)+z = x+(y+z)$;

    * (c) There exists a zero element $\theta$ such that $x+\theta = x$ for every $x \in E$;

    * (d) For any $x \in E$, there exists an additive inverse element $-x$ such that $x+(-x) = \theta$.

    **(ii) Scalar Multiplication Properties**: For any $x, y \in E$ and any numbers $\alpha, \beta \in K$, there corresponds a unique product $\alpha x \in E$ such that:

    * (a) $\alpha(\beta x) = (\alpha\beta)x$;

    * (b) $1 \cdot x = x$;

    * (c) $(\alpha+\beta)x = \alpha x + \beta x$;

    * (d) $\alpha(x+y) = \alpha x + \alpha y$.

Elements in a linear space are also called **vectors**. The addition of elements and the multiplication of numbers with elements are collectively referred to as **linear operations**. For example, $C[0,1]$ and $L^{2}[0,1]$ are linear spaces.

---

## 2. Properties of Linear Spaces

The following properties can be derived from the definition of a linear space:

* (i) $0 \cdot x = \theta$.

* (ii) $(-1)x = -x$.

* (iii) $\alpha \cdot \theta = \theta$.

* (iv) $\alpha x = \theta$ if and only if $\alpha = 0$ or $x = \theta$.

??? proof "Proof Details"

    **Regarding (i):**
    Since $2(0 \cdot x) = (2 \cdot 0)x = 0 \cdot x$, subtracting $0 \cdot x$ from both sides yields:

    $$
    0 \cdot x = 2(0 \cdot x) - 0 \cdot x = \theta
    $$

    **Regarding (ii):**
    Since $(-1)x + x = (-1 + 1)x = 0 \cdot x = \theta$, by the uniqueness of the inverse element, $(-1)x = -x$.

    **Regarding (iii):**

    $$
    \alpha \cdot \theta = \alpha(\theta + (-\theta)) = \alpha \theta + \alpha(-\theta) = \alpha \theta - \alpha \theta = \theta
    $$

    **Regarding (iv):**
    Sufficiency is obvious. For necessity: if $\alpha x = \theta$ and $\alpha \neq 0$, then:

    $$
    x = 1 \cdot x = \left( \frac{1}{\alpha} \cdot \alpha \right) x = \frac{1}{\alpha} (\alpha x) = \frac{1}{\alpha} \theta = \theta
    $$

---

## 3. Subspace, Span, and Isomorphism

### 3.1 Subspace

!!! info "Definition 5.2 (Subspace)"

    Let $E_0$ be a subset of a linear space $E$. $E_0$ is called a **linear subspace** of $E$ if it is closed under the linear operations of $E$ (i.e., for any $x, y \in E_0$ and $\alpha \in K$, then $x+y \in E_0$ and $\alpha x \in E_0$).

* A subspace of $E$ that is not equal to $E$ itself is called a **proper subspace**.

* **Linear Independence**: A set $A \subset E$ is said to be linearly independent if any finite number of elements in $A$ are linearly independent.

### 3.2 Subspace Spanned by a Subset (span L)

!!! info "Definition 5.3 (span L)"

    Let $L$ be a non-empty subset of a linear space $E$. The set consisting of all possible finite linear combinations of elements in $L$ is called the **subspace spanned by $L$**, denoted as $\text{span } L$:

    $$
    \text{span } L = \left\{ \sum_{k=1}^{n} c_k x_k \mid c_k \in K, x_k \in L, n \in \mathbb{N} \right\}
    $$

* $\text{span } L$ is the **smallest subspace** containing $L$, and it is the intersection of all subspaces containing $L$.

* For example, in $L^2[-\pi, \pi]$, the subset $L = \{ \frac{1}{\sqrt{\pi}} \sin kt, \frac{1}{\sqrt{\pi}} \cos kt \}_{k=0}^{\infty}$ spans the subspace of all trigonometric polynomials.

### 3.3 Isomorphism of Linear Spaces

* Let $E$ and $E'$ be linear spaces. $E$ is said to be **isomorphic** to $E'$ if there exists a bijective mapping $T$ from $E$ to $E'$ such that $T(x+y) = Tx + Ty$ and $T(ax) = aTx$.

---

## 4. Direct Sum

!!! info "Definition 5.4 (Direct Sum)"

    Let $L_1, \dots, L_n$ be subspaces of $E$. If any element $x \in E$ can be **uniquely** expressed as $x = x_1 + \dots + x_n$ ($x_k \in L_k$), then $E$ is called the direct sum of these subspaces, denoted by $E = L_1 \oplus \dots \oplus L_n$.

* **Property**: If $E$ is the direct sum of $L_1, \dots, L_n$, then any non-zero elements $x_1, \dots, x_n$ chosen from each respective subspace must be linearly independent.

* **Construction**: Consider the set $E$ of all ordered tuples $x = (x_1, \dots, x_n)$ where $x_k \in L_k$. By defining component-wise addition and scalar multiplication, this set forms a direct sum of the spaces $L_k$.

---

## 5. Normed Linear Space

To introduce topology into a linear space, the structure must be integrated with the linear operations.

!!! info "Definition 5.5 (Norm)"

    Let $E$ be a linear space. If for every element $x \in E$, there corresponds a real number $\|x\|$ such that:

    * (i) Non-negativity: $\|x\| \ge 0$, and $\|x\| = 0 \iff x = \theta$;

    * (ii) Homogeneity: $\|\alpha x\| = |\alpha| \|x\|$;

    * (iii) Triangle Inequality: $\|x+y\| \le \|x\| + \|y\|$.

    Then $E$ is called a **normed linear space**, and $\|x\|$ is called the **norm** of the vector $x$.

### 5.1 Distance and Strong Convergence

A normed linear space naturally induces a distance: $\rho(x, y) = \|x - y\|$.

* **Convergence in Norm (Strong Convergence)**: A sequence $\{x_n\}$ converges to $x$ if $\|x_n - x\| \rightarrow 0$.

### 5.2 Continuity of Linear Operations

* (i) The norm $\|x\|$ is a continuous functional: $| \|x_n\| - \|x\| | \le \|x_n - x\|$. If $x_n$ converges in norm, then $\{x_n\}$ is bounded.

* (ii) Continuity of Addition: If $x_n \to x$ and $y_n \to y$, then $x_n + y_n \to x+y$.

    $$
    \|x_n + y_n - (x+y)\| \le \|x_n - x\| + \|y_n - y\| \rightarrow 0
    $$

* (iii) Continuity of Scalar Multiplication: If $\alpha_n \to \alpha$ and $x_n \to x$, then $\alpha_n x_n \to \alpha x$.

??? proof "Proof: Continuity of Scalar Multiplication"

    $$
    \|\alpha_n x_n - \alpha x\| = \|\alpha_n x_n - \alpha_n x + \alpha_n x - \alpha x\|
    $$

    $$
    \le |\alpha_n| \|x_n - x\| + |\alpha_n - \alpha| \|x\| \rightarrow 0
    $$

---

## 6. Classical Examples of Normed Linear Spaces

* **1. $\mathbb{R}^n$ and $\mathbb{C}^n$**:
    For $x = (\xi_1, \dots, \xi_n)$, the Euclidean norm is defined as:

    $$
    \|x\| = \sqrt{\sum_{j=1}^n |\xi_j|^2}
    $$

* **2. Continuous Function Space $C[a, b]$**:
    The maximum norm is defined as:

    $$
    \|x\| = \max_{t \in [a, b]} |x(t)|
    $$

* **3. Sequence Space $l^p$ ($1 \le p < \infty$)**:
    Consists of sequences satisfying $\sum_{j=1}^\infty |\xi_j|^p < \infty$. The norm is:

    $$
    \|x\|_p = \left( \sum_{j=1}^\infty |\xi_j|^p \right)^{\frac{1}{p}}
    $$

    The triangle inequality $\|x+y\|_p \le \|x\|_p + \|y\|_p$ is precisely the Minkowski inequality.

* **4. Integrable Function Space $L^p(F)$ ($1 \le p < \infty$)**:
    For measurable functions on $F \subset \mathbb{R}$, the norm is defined as:

    $$
    \|f\|_p = \left( \int_F |f(t)|^p dt \right)^{\frac{1}{p}}
    $$

* **5. Space of Continuously Differentiable Functions $C^k[a, b]$**:
    Functions possessing continuous derivatives up to order $k$. The norm is:

    $$
    \|x\| = \sum_{j=0}^{k} \max_{t \in [a, b]} |x^{(j)}(t)|
    $$