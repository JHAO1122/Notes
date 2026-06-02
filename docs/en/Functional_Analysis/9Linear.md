# Chapter 9: Linear Operators

In this chapter, we delve into the core objects of functional analysis—linear operators. This part covers the basic concepts and properties of bounded linear operators, as well as the completeness structure of spaces of linear operators.

---

## 1. Basic Concepts and Properties of Bounded Linear Operators

### 1.1 Definition of Linear Operators

!!! info "Definition 9.1 (Linear Operator)"

    Let $E$ and $E_1$ both be real (or complex) linear spaces, and let $T$ be a mapping from some subspace $D$ of $E$ into the linear space $E_1$.

    If for any $x, y \in D$, we have:

    \[
    T(x+y) = T(x) + T(y)
    \]

    then $T$ is called **additive**.

    If for any real (or complex) number $a$ and any $x \in D$, we have:

    \[
    T(ax) = aT(x)
    \]

    then $T$ is called **homogeneous**.

    A mapping that satisfies both additivity and homogeneity is called a **linear mapping** or **linear operator**.

* $D$ is called the domain of $T$, often denoted by $D_T$.

* $T(D)$ is called the range of $T$, often denoted by $T(D_T)$.

* The set of elements in $D$ such that $Tx = \theta$ is called the **null space** of $T$, often denoted by $N_T$ or similar symbols.

* Let $E_1$ be the real (or complex) number field, then $T$ becomes an operator from $D$ into the real (or complex) number field. In this case, $T$ is called a **functional**. If $T$ is also linear, then $T$ is called a **linear functional**.

---

### 1.2 Continuous and Bounded Linear Operators

!!! info "Definition 9.2 (Continuous and Bounded Linear Operators)"

    Let $E$ and $E_1$ both be real (or complex) normed linear spaces, and let $T$ be a linear operator from some subspace $D$ of $E$ into the linear space $E_1$.

    * **Continuous linear operator**: If $T$ is continuous on its domain, then $T$ is called a continuous linear operator.

    * **Bounded linear operator**: If $T$ maps every bounded set in $D$ to a bounded set in $E_1$, then $T$ is called a bounded linear operator.

    * **Unbounded linear operator**: If there exists a bounded set $A$ in $D$ such that $T(A)$ is an unbounded set in $E_1$, then $T$ is called an unbounded linear operator.

* The operator that maps every element $x$ in the normed linear space $E$ to itself is called the **identity operator** on $E$, often denoted by $I$.

* The operator that maps every element $x$ in $E$ to the zero element $\theta$ is called the **zero operator**.

* The identity operator and the zero operator are both bounded linear operators and continuous linear operators.

---

### 1.3 Typical Example of an Unbounded Operator

Consider the differential operator $T = \frac{d}{dt}$ on the space of continuous functions $C[a,b]$.

Take the set $C^1[a,b]$, consisting of continuously differentiable functions on $[a,b]$, as the domain $D$ of $T$. Then $T$ is a linear operator defined on $C^1[a,b]$ and taking values in $C[a,b]$.

We can prove that $T$ is an **unbounded operator**:

Take the sequence of functions:

\[
x_n(t) = \sin(nt)
\]

Then $x_n \in C^1[0,1]$ and $\|x_n\| = 1$, but:

\[
\|Tx_n\| = \left\| \frac{d}{dt} \sin(nt) \right\| = n\|\cos(nt)\| = n
\]

Hence $T$ maps the unit sphere in $C[0,1]$ to an unbounded set in $C[0,1]$. Therefore $T$ is unbounded.

---

### 1.4 Equivalence Between Boundedness and Continuity of Linear Operators

First, consider an example of a continuous linear functional. The integral of continuous functions:

\[
f(x) = \int_a^b x(t) dt
\]

is a bounded linear functional, and also a continuous linear functional, defined on the space of continuous functions $C[a,b]$. Because:

\[
|f(x)| \le (b-a)\|x\|
\]

!!! success "Theorem 9.1"

    Let $E$ and $E_1$ both be real normed linear spaces, and let $T$ be a continuous additive operator from a subspace $D$ of $E$ into $E_1$. Then $T$ satisfies homogeneity. Therefore $T$ is a continuous linear operator.

??? proof "Proof of Theorem 9.1 (Click to Expand)"

    For any given $x \in D$, set $f(\alpha) = T(\alpha x)$. Then $f$ is continuous.
  
    If $\alpha_n \rightarrow \alpha_0$, then $\alpha_n x \rightarrow \alpha_0 x$. Hence $T(\alpha_n x) \rightarrow T(\alpha_0 x)$, i.e., $f(\alpha_n) \rightarrow f(\alpha_0)$.

    For any two real numbers $\alpha, \beta$:

    \[
    f(\alpha+\beta) = T((\alpha+\beta)x) = T(\alpha x + \beta x) = T(\alpha x) + T(\beta x) = f(\alpha) + f(\beta)
    \]

    By a lemma from real-variable functions, it follows that for any real number $f(\alpha) = \alpha f(1)$, i.e., $T(\alpha x) = \alpha T(x)$. This completes the proof. $\square$

* **Corollary**: Let $E$ and $E_1$ both be complex normed linear spaces, and let $T$ be a continuous additive operator from a subspace $D$ of $E$ into $E_1$, satisfying $T(ix) = iT(x)$. Then $T$ is homogeneous.

    Proof: Let $\alpha = a + ib$,
  
    \[
    T(\alpha x) = T(ax+ibx) = T(ax) + T(ibx) = aT(x) + bT(ix) = aT(x) + biT(x) = \alpha T(x)
    \]

!!! success "Theorem 9.2 (Equivalent Definition of Bounded Linear Operators)"

    Let $E$ and $E_1$ both be normed linear spaces, and let $T$ be a linear operator from a subspace $D$ of $E$ into $E_1$. Then $T$ is bounded if and only if there exists $M > 0$ such that for all $x \in D$,

    \[
    \|Tx\| \le M\|x\|
    \]

??? proof "Proof of Theorem 9.2 (Click to Expand)"

    **Sufficiency**: Let $A \subset D$ be any bounded set. Then there exists a number $K$ such that $\|x\| \le K$ for all $x \in A$. Hence for $x \in A$,

    \[
    \|Tx\| \le M\|x\| \le MK
    \]

    Therefore $T(A)$ is a bounded set in $E_1$, i.e., $T$ is bounded.

    **Necessity**: Take the unit sphere in $D$: $S := \{x : \|x\| = 1, x \in D\}$.
    Since $S$ is bounded, $T(S)$ is also bounded. Thus there exists a positive number $M > 0$ such that $\|Tx\| \le M$ for $x \in S$.
  
    Let $x$ be any nonzero element in $D$. Then $\frac{x}{\|x\|} \in S$. Hence:

    \[
    \left\| T\left(\frac{x}{\|x\|}\right) \right\| \le M
    \]

    By homogeneity of $T$,

    \[
    \|Tx\| \le M\|x\|
    \]

    When $x = \theta$, the inequality holds trivially. This completes the proof. $\square$

!!! success "Theorem 9.3 (Equivalence of Boundedness and Continuity)"

    Let $E$ and $E_1$ both be normed linear spaces, and let $T$ be a linear operator from a subspace $D$ of $E$ into $E_1$. Then the following properties are equivalent:

    (i) $T$ is continuous;

    (ii) $T$ is continuous at the origin $\theta$;

    (iii) $T$ is bounded.

??? proof "Proof of Theorem 9.3 (Click to Expand)"

    **(i) $\Rightarrow$ (ii)**: Trivial.

    **(ii) $\Rightarrow$ (iii)**: Since $T$ is continuous at the origin, for $\epsilon=1$, there exists $\delta>0$ such that for any $x \in D$ with $\|x\| \le \delta$, we have:

    \[
    \|Tx\| = \|Tx - T\theta\| \le 1
    \]

    Let $x$ be any nonzero element in $D$. Then $\left\| \frac{\delta x}{\|x\|} \right\| = \delta \le \delta$. Hence:

    \[
    \left\| T\left(\frac{\delta x}{\|x\|}\right) \right\| \le 1
    \]

    Extracting the coefficient by homogeneity yields:

    \[
    \|Tx\| \le \frac{1}{\delta}\|x\|
    \]

    This shows the operator is bounded.

    **(iii) $\Rightarrow$ (i)**: Let $x_n, x \in D$, and $x_n \rightarrow x$. Since $T$ is bounded,

    \[
    \|Tx_n - Tx\| = \|T(x_n - x)\| \le M\|x_n - x\| \rightarrow 0
    \]

    i.e., $Tx_n \rightarrow Tx$, proving continuity.

    For linear operators, boundedness, continuity, and continuity at the origin are all equivalent. These conditions are also equivalent to continuity of $T$ at any given point $x_0$ in $D$. $\square$

---

### 1.5 Norm of a Bounded Linear Operator

!!! info "Definition 9.3 (Operator Norm)"

    Let $E$ and $E_1$ both be normed linear spaces, and let $T$ be a bounded linear operator from a subspace $D$ of $E$ into $E_1$. The infimum of all positive numbers $M$ such that $\|Tx\| \le M\|x\|$ holds for all $x \in D$ is called the **norm** of $T$, denoted by $\|T\|$.

Since $M$ is an upper bound for the set $\{\frac{\|Tx\|}{\|x\|} : x \in D, x \ne \theta\}$, the norm $\|T\|$ of the operator $T$, being the infimum of all such upper bounds $M$, is also the least upper bound (i.e., supremum) of the above set. Thus:

\[
\|T\| = \sup_{x \in D, x \ne \theta} \frac{\|Tx\|}{\|x\|}
\]

From this we easily derive the following:

* For all $x \in D$, $\|Tx\| \le \|T\|\|x\|$.

* Other equivalent computation forms for the norm:

\[
\|T\| = \sup_{x \in D, x \ne \theta} \left\| T\left(\frac{x}{\|x\|}\right) \right\| = \sup_{x \in D, \|x\|=1} \|Tx\| = \sup_{x \in D, \|x\| \le 1} \|Tx\|
\]

The norm of an operator is an important quantity characterizing bounded linear operators. In general, finding the norm of a bounded linear operator is very difficult. Below, we illustrate through several examples how to estimate or exactly compute its norm.

---

### 1.6 Examples of Operator Norm Computation

**Example 1: Linear Operator Represented by a Matrix**

Let $(a_{ij})$ be a given $n \times n$ real square matrix. The equation:

\[
\eta_i = \sum_{j=1}^n a_{ij} \xi_j
\]

defines an operator $T : \mathbb{R}^n \to \mathbb{R}^n$, $Tx = y$. It maps the element $x = (\xi_1, \dots, \xi_n)$ to the element $y = (\eta_1, \dots, \eta_n)$.

Take another vector $x' = (\xi_1', \dots, \xi_n')$ in $\mathbb{R}^n$. It is easy to verify that $T(x+x') = Tx + Tx'$ and $T(ax) = aTx$, so $T$ is a linear operator.

We estimate its boundedness. By the Cauchy inequality,

\[
\|Tx\|^2 = \sum_{i=1}^n \eta_i^2 = \sum_{i=1}^n \left(\sum_{j=1}^n a_{ij} \xi_j\right)^2 \le \sum_{i=1}^n \left( \sum_{j=1}^n a_{ij}^2 \sum_{j=1}^n \xi_j^2 \right) \le \sum_{i=1}^n \sum_{j=1}^n a_{ij}^2 \cdot \|x\|^2
\]

Thus:

\[
\|T\| \le \left( \sum_{i,j=1}^n a_{ij}^2 \right)^{\frac{1}{2}}
\]

Hence the operator is continuous and bounded.

**Example 2: Integral Operator of Fourier Type**

$C(\mathbb{R})$: the set of bounded continuous functions defined on $\mathbb{R}$, with norm $\|y\| = \sup_{t \in \mathbb{R}} |y(t)|$. Then $C(\mathbb{R})$ is a Banach space.

Let $x \in L(\mathbb{R})$ (i.e., the space of integrable functions). Define $y = Tx$ by

\[
y(t) = \int_{-\infty}^\infty e^{-ist} x(s) ds
\]

$T$ is a linear operator defined on $L(\mathbb{R})$ with range contained in $C(\mathbb{R})$. By the norm definition:

\[
|Tx(t)| = |y(t)| \le \int_{-\infty}^\infty |e^{-ist} x(s)| ds = \int_{-\infty}^\infty |x(s)| ds = \|x\|
\]

Hence $\|T\| \le 1$. $T$ is bounded, hence continuous.

**Example 3: Lagrange Interpolation Polynomial**

Let $x \in C[a,b]$. Choose $n$ points $a \le t_1 < t_2 < \dots < t_n \le b$ in $[a,b]$, and construct the polynomials:

\[
l_k(t) = \frac{(t-t_1)\dots(t-t_{k-1})(t-t_{k+1})\dots(t-t_n)}{(t_k-t_1)\dots(t_k-t_{k-1})(t_k-t_{k+1})\dots(t_k-t_n)}
\]

Define the operator $y = L_n x = \sum_{k=1}^n x(t_k) l_k(t)$.

Then:

\[
\|L_n x\| = \max_{t \in [a,b]} \left| \sum_{k=1}^n x(t_k) l_k(t) \right| \le \max_{t \in [a,b]} \sum_{k=1}^n |l_k(t)| \max_{t \in [a,b]} |x(t)| = \left( \max_{t \in [a,b]} \sum_{k=1}^n |l_k(t)| \right) \|x\|
\]

Let $\gamma := \max_{t \in [a,b]} \sum_{k=1}^n |l_k(t)|$. Then $\|L_n\| \le \gamma$.

On the other hand, since $\sum_{k=1}^n |l_k(t)|$ is continuous on $[a,b]$, there exists $t_0 \in [a,b]$ such that:

\[
\gamma = \sum_{k=1}^n |l_k(t_0)|
\]

Take $x_0 \in C[a,b]$ satisfying $\|x_0\| = 1$ and $x_0(t_k) = \text{sgn}(l_k(t_0))$ at each node. Then:

\[
\|L_n\| \ge \|L_n(x_0)\| \ge |L_n(x_0)(t_0)| = \left| \sum_{k=1}^n l_k(t_0) x_0(t_k) \right| = \left| \sum_{k=1}^n l_k(t_0) \text{sgn}(l_k(t_0)) \right| = \sum_{k=1}^n |l_k(t_0)| = \gamma
\]

Thus we have proven the exact norm $\|L_n\| = \gamma$.

**Example 4: Integral Operator with Kernel**

Let $K(t,s)$ be a continuous function defined on $a \le t \le b$, $a \le s \le b$. Define the integral operator on $C[a,b]$ by:

\[
y(t) = Tx(t) = \int_a^b K(t,s) x(s) ds
\]

Then $T$ is a bounded linear operator from $C[a,b]$ into itself, and the norm satisfies:

\[
\|T\| = \gamma, \quad \gamma := \max_{t \in [a,b]} \int_a^b |K(t,s)| ds
\]

*Proof of upper bound*:

\[
\|Tx\| = \max_{t \in [a,b]} \left| \int_a^b K(t,s) x(s) ds \right| \le \max_{s \in [a,b]} |x(s)| \max_{t \in [a,b]} \int_a^b |K(t,s)| ds = \gamma \|x\|
\]

i.e., $\|T\| \le \gamma$.

*Proof of lower bound*:
Since $\int_a^b |K(t,s)| ds$ is continuous, there exists $t_0$ attaining the maximum $\gamma$. Take $\|\varphi\| = 1$,

\[
\|T\| \ge \|T\varphi\| = \max_{t \in [a,b]} \left| \int_a^b K(t,s) \varphi(s) ds \right| \ge \left| \int_a^b K(t_0,s) \varphi(s) ds \right|
\]

Let $E_0 = \{s : K(t_0, s) \ge 0\}$. Construct the function:

\[
\varphi_n(s) = \frac{1 - n d(s, E_0)}{1 + n d(s, E_0)}
\]

where $d(s, E_0)$ is the distance from $s$ to $E_0$. Then $\varphi_n$ is continuous on $[a,b]$, and $|\varphi_n| \le 1$.
Furthermore, $\varphi_n|_{E_0} \equiv 1$, and $\varphi_n|_{E_0^c} \rightarrow -1$.

By the Lebesgue Dominated Convergence Theorem, as $n \rightarrow \infty$:

\[
T\varphi_n(t_0) = \int_a^b K(t_0, s) \varphi_n(s) ds \rightarrow \int_a^b |K(t_0, s)| ds = \gamma
\]

Hence:

\[
\gamma = \lim_{n \rightarrow \infty} T\varphi_n(t_0) \le \lim_{n \rightarrow \infty} \|T\varphi_n\| \le \lim_{n \rightarrow \infty} \|T\|\|\varphi_n\| = \|T\|
\]

Thus $\|T\| \ge \gamma$.

---

## 2. Spaces of Linear Operators

### 2.1 The Space $\mathbb{B}(E, E_1)$

We denote by $\mathbb{B}(E, E_1)$ the set of **all bounded linear operators** from a normed linear space $E$ to a normed linear space $E_1$.

!!! success "Theorem 9.4"

    Let $E$ and $E_1$ both be normed linear spaces. Define linear operations in $\mathbb{B}(E, E_1)$ as follows:

    * Addition: $(T + T')x = Tx + T'x$
    * Scalar multiplication: $(\alpha T)x = \alpha Tx$

    where $T, T' \in \mathbb{B}(E, E_1)$ and $\alpha$ is a scalar. Then $\mathbb{B}(E, E_1)$ is a linear space under the above linear operations. If we take the operator norm as its norm, then $\mathbb{B}(E, E_1)$ is a normed linear space.

??? proof "Proof of Theorem 9.4 (Click to Expand)"

    First, by definition it is a linear space. We need to verify that the operator norm satisfies the three axioms of a norm:

    * **Non-negativity**: $\|T\| = \sup_{\|x\|=1} \|Tx\| \ge 0$. If $T = \theta$ (zero operator), clearly $\|T\| = 0$. Conversely, if $\|T\| = 0$, then $Tx = 0$ for all $x \in E$, so $T = \theta$.

    * **Homogeneity**:
  
    \[
    \|\alpha T\| = \sup_{\|x\|=1} \|\alpha Tx\| = |\alpha| \sup_{\|x\|=1} \|Tx\| = |\alpha|\|T\|
    \]

    * **Triangle inequality**:
  
    \[
    \|T + T'\| = \sup_{\|x\|=1} \|(T + T')x\| = \sup_{\|x\|=1} \|Tx + T'x\| \le \sup_{\|x\|=1} (\|Tx\| + \|T'x\|) \le \sup_{\|x\|=1} \|Tx\| + \sup_{\|x\|=1} \|T'x\| \le \|T\| + \|T'\|
    \]

    Therefore, it is a normed linear space. $\square$

*Note: We denote $\mathbb{B}(E, E)$ simply by $\mathbb{B}(E)$, representing the set of all bounded linear operators from $E$ into itself.*

---

### 2.2 Convergence in Operator Norm (Uniform Operator Topology)

Let $T, T_n \in \mathbb{B}(E, E_1)$ ($n=1,2,3,\dots$). If $T_n$ converges to $T$ in the norm of $\mathbb{B}(E, E_1)$, i.e.,

\[
\|T - T_n\| \rightarrow 0
\]

then we say that $T_n$ converges to $T$ **in operator norm**, or $T_n$ converges to $T$ in the **uniform operator topology**.
The reason for the term "convergence in the uniform operator topology" is given by the following theorem:

!!! success "Theorem 9.5"

    Let $T, T_n \in \mathbb{B}(E, E_1)$. $T_n$ converges to $T$ in operator norm if and only if $\{T_n\}$ converges **uniformly** to $T$ on **every bounded set** in $E$.

??? proof "Proof of Theorem 9.5 (Click to Expand)"

    **Necessity**: Let $A \subset E$ be a bounded set. Then there exists $K > 0$ such that $\|x\| \le K$ for $x \in A$. Hence:

    \[
    \|Tx - T_nx\| \le \|T - T_n\|\|x\| \le K\|T - T_n\|
    \]

    Since $\|T - T_n\| \rightarrow 0$, for any $\epsilon > 0$ there exists $N > 0$ such that for $n > N$, $\|T - T_n\| < \epsilon/K$. Substituting into the above gives:

    \[
    \|Tx - T_nx\| \le \epsilon
    \]

    This holds uniformly for all $x \in A$. Hence $\{T_n\}$ converges uniformly to $T$ on $A$.

    **Sufficiency**: Suppose $\{T_n\}$ converges uniformly to $T$ on every bounded set in $E$. Take the unit sphere in $E$: $S = \{x : \|x\| = 1, x \in E\}$.
    By assumption, for any $\epsilon > 0$ there exists $N > 0$ such that for $n > N$, the inequality:

    \[
    \|Tx - T_nx\| \le \epsilon
    \]

    holds uniformly for all $x \in S$. Then:

    \[
    \|T - T_n\| = \sup_{\|x\|=1} \|Tx - T_nx\| \le \epsilon
    \]

    Hence $\{T_n\}$ converges to $T$ in operator norm. $\square$

---

### 2.3 Completeness of $\mathbb{B}(E, E_1)$

The completeness of the operator space is determined by the completeness of the target space.

!!! success "Theorem 9.6"

    Let $E_1$ be a Banach space. Then the operator space $\mathbb{B}(E, E_1)$ is also a Banach space.

??? proof "Proof of Theorem 9.6 (Click to Expand)"

    Let $\{T_n\}$ be a Cauchy sequence in $\mathbb{B}(E, E_1)$. Then for any $\epsilon > 0$, there exists $N > 0$ such that for $m, n > N$:

    \[
    \|T_m - T_n\| < \epsilon
    \]

    Take any $x \in E$ and consider the sequence $\{T_n x\}$:

    \[
    \|T_mx - T_nx\| \le \|T_m - T_n\|\|x\| \le \epsilon\|x\|
    \]

    Thus, for fixed $x$, $\{T_n x\}$ is a Cauchy sequence in $E_1$. Since $E_1$ is complete, $\{T_n x\}$ converges to some element $y$ in $E_1$. Hence we can define the pointwise limit operator $T$ by:

    \[
    Tx = \lim_{n \rightarrow \infty} T_nx = y
    \]

    Now we need to prove that $T$ belongs to $\mathbb{B}(E, E_1)$ and that $T_n \to T$ in norm.

    **1. Linearity**:

    \[
    T(ax+by) = \lim_{n \rightarrow \infty} T_n(ax+by) = \lim_{n \rightarrow \infty} (aT_nx + bT_ny) = aTx + bTy
    \]

    **2. Boundedness and convergence in operator norm**:
    In the earlier Cauchy condition $\|T_mx - T_nx\| \le \epsilon\|x\|$, let $m \rightarrow \infty$. By continuity of the norm, we obtain for all $n > N$:

    \[
    \|Tx - T_nx\| \le \epsilon\|x\|
    \]

    This shows that $T - T_n$ is bounded, i.e., $T - T_n \in \mathbb{B}(E, E_1)$. Since $T_n$ is bounded, $T = (T - T_n) + T_n$ is necessarily a bounded linear operator, so $T \in \mathbb{B}(E, E_1)$.

    Moreover, this also implies that for $n > N$:

    \[
    \|T - T_n\| \le \epsilon
    \]

    Hence $\{T_n\}$ converges to $T$ in operator norm. Thus every fundamental sequence in $\mathbb{B}(E, E_1)$ has a limit, so it is a complete Banach space. $\square$

---

### 2.4 Example of Non-Convergence in Operator Norm

Although pointwise convergence can define a limit operator, pointwise convergence does not imply convergence in operator norm.

Define the **truncation shift operator** $T_n$ in $l^p$:

\[
T_n x = x_n
\]

where the input vector is $x = (\xi_1, \xi_2, \dots, \xi_n, \dots) \in l^p$, and the truncated output vector is $x_n = (\xi_n, \xi_{n+1}, \dots) \in l^p$ (according to the context of the notes, i.e., truncating the first $n-1$ components).

$T_n$ is a bounded linear operator, and $\|T_nx\|_p \le \|x\|_p$.
For any given $x \in l^p$, since the series $\sum |\xi_k|^p$ converges, its tail tends to 0, so $\|x_n\|_p \rightarrow 0$. Hence pointwise:

\[
\lim_{n \rightarrow \infty} T_n x = \theta
\]

(i.e., the operator sequence converges pointwise to the zero operator).

However, examine the operator norm $\|T_n\|$. Take the specific unit vector:

\[
y_n = (0, \dots, 0, 1, 0, \dots) \quad \text{(the } n\text{-th position is 1)}
\]

Then $T_n y_n = (1, 0, \dots)$. Hence:

\[
\|T_n\| \ge \|T_n y_n\|_p = 1
\]

So $\|T_n\| = 1$. This shows that $\|T_n - 0\| = 1 \not\rightarrow 0$.

Therefore, the sequence $\{T_n\}$ **does not converge** to the zero operator **in operator norm**.