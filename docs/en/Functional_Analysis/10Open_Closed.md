# Chapter 10: Open Mapping Theorem and Closed Graph Theorem

This part continues from the previous chapter, first introducing strong convergence of operators and algebraic operations, and then focusing on two of the three fundamental theorems of functional analysis: the **Open Mapping Theorem** and the **Closed Graph Theorem**.

---

## 1. Strong Convergence and Algebraic Operations of Operators

### 1.1 Strong Convergence

In addition to convergence in operator norm (uniform operator topology), there is another common mode of convergence for sequences of operators, namely pointwise convergence at each point.

!!! info "Definition 10.1 (Strong Convergence)"

    Let $T, T_n \in \mathbb{B}(E, E_1)$. If for every $x \in E$, we have:

    \[
    \lim_{n \rightarrow \infty} \|Tx - T_nx\| = 0
    \]

    then $\{T_n\}$ is said to converge **strongly** to $T$, or $\{T_n\}$ converges to $T$ in the **strong operator topology**. It is denoted as:

    \[
    \lim_{n \rightarrow \infty} T_n = T \text{ (strongly)}
    \]

* **Relationship between strong convergence and uniform convergence**: Any sequence of operators that converges in operator norm to some operator necessarily converges strongly to the same operator. Because:

    \[
    \|Tx - T_nx\| \le \|T - T_n\|\|x\| \rightarrow 0
    \]

    The converse does not hold. For example, the truncation shift operators $T_n$ on $l^p$ mentioned at the end of the previous chapter do not converge in operator norm to the zero operator, but they converge strongly to the zero operator, because for any $x \in l^p$:

    \[
    \|T_n x\|_p = \|x_n\|_p \rightarrow 0
    \]

---

### 1.2 Multiplication of Operators

!!! info "Definition 10.2 (Multiplication of Operators)"

    Let $E, E_1, E_2$ be normed linear spaces, $T_1 \in \mathbb{B}(E, E_1)$, $T_2 \in \mathbb{B}(E_1, E_2)$. Define the product operator $T_2 T_1$ of $T_1$ and $T_2$ as follows:

    \[
    (T_2 T_1)x = T_2(T_1 x)
    \]

Multiplication of operators satisfies the following algebraic properties:

* **Associativity**: $T_3(T_2 T_1) = (T_3 T_2)T_1$, where $T_3 \in \mathbb{B}(E_2, E_3)$.
  *Proof*: Take any $x \in E$, then
  $(T_3(T_2 T_1))x = T_3((T_2 T_1)x) = T_3(T_2(T_1 x)) = (T_3 T_2)(T_1 x) = ((T_3 T_2)T_1)x$
* **Scalar multiplication associativity**: $(\alpha T_2)T_1 = \alpha(T_2 T_1)$ and $T_2(\alpha T_1) = \alpha(T_2 T_1)$.
* **Distributivity**: $T_2(T_1 + T_1') = T_2 T_1 + T_2 T_1'$ and $(T_2 + T_2')T_1 = T_2 T_1 + T_2' T_1$.
* **Boundedness**: $T_2 T_1 \in \mathbb{B}(E, E_2)$, and $\|T_2 T_1\| \le \|T_2\|\|T_1\|$.
  *Proof*: Take any $x \in E$, then
  $\|(T_2 T_1)x\| = \|T_2(T_1 x)\| \le \|T_2\|\|T_1 x\| \le \|T_2\|\|T_1\|\|x\|$

If $T_1, T_2 \in \mathbb{B}(E)$, then both $T_2 T_1$ and $T_1 T_2$ are well-defined. However, multiplication is generally **non-commutative**, i.e., $T_2 T_1 \ne T_1 T_2$. If they are equal, then $T_1, T_2$ are said to commute.
(Note: $T^n$ denotes the product of $n$ copies of $T$, and $T^0$ denotes the identity operator $I$.)

**Example of non-commutativity**:
Let $E = C[0,1]$, and consider the following two bounded linear operators:

\[
(T_1 x)(t) = \int_0^t x(s) ds, \quad (T_2 x)(t) = t x(t)
\]

($T_1$ is the Volterra integral operator, $T_2$ is the multiplication operator). Compute their products:

\[
(T_2 T_1 x)(t) = t \int_0^t x(s) ds, \quad (T_1 T_2 x)(t) = \int_0^t s x(s) ds
\]

Take $x_0 \equiv 1$, then:

\[
(T_2 T_1 x_0)(t) = t^2 \ne \frac{t^2}{2} = (T_1 T_2 x_0)(t)
\]

Hence $T_2 T_1 \ne T_1 T_2$.

---

### 1.3 Normed Algebras and Banach Algebras

!!! info "Definition 10.3 (Normed Algebra)"

    Let $\mathscr{B}$ be a non-empty set, with algebraic operations (linear operations and multiplication) defined on $\mathscr{B}$, where multiplication is associative and distributive over addition. Additionally, a norm is defined, such that for any elements $T_1, T_2$ in $\mathscr{B}$, we have:

    \[
    \|T_2 T_1\| \le \|T_2\|\|T_1\|
    \]

    Then $\mathscr{B}$ is called a **normed algebra**. If it is also complete, it is called a **Banach algebra**.

For a normed linear space $E$, $\mathbb{B}(E)$ forms a normed algebra; if $E$ is complete (i.e., a Banach space), then $\mathbb{B}(E)$ is a Banach algebra.

---

## 2. Banach Open Mapping Theorem and Inverse Operator Theorem

### 2.1 Banach Open Mapping Theorem

!!! success "Theorem 10.1 (Banach Open Mapping Theorem)"

    Let a bounded linear operator $T$ map a Banach space $E$ onto a set $F$ of the second category in a Banach space $E_1$. Then the following conclusions hold:

    (i) $F = E_1$, i.e., the range of the operator $T$ is the whole space $E_1$;

    (ii) There exists a positive number $K$ such that for every $y \in E_1$, there exists $x \in E$ with $Tx = y$ and $\|x\| \le K\|Tx\|$.

??? proof "Proof of Theorem 10.1 (Click to Expand)"

    The proof is divided into four steps:

    **(i) Prove that the image is dense in some closed ball**:
    Let $O_n = \{x \in E : \|x\| \le n\}$, $M_n = T(O_n)$, i.e., $M_n$ is the image of $O_n$. Since $E = \bigcup_{n=1}^\infty O_n$, we have $F = \bigcup_{n=1}^\infty M_n$.
    Because $F$ is a set of the second category, by the Baire category theorem, there exists $n_0$ such that $M_{n_0}$ is not nowhere dense. Hence there exists a closed ball in $E_1$:

    \[
    Q(y_0, r_0) := \{y \in E_1 : \|y - y_0\| \le r_0\}
    \]

    such that $M_{n_0} = T(O_{n_0})$ is dense in $Q(y_0, r_0)$.

    **(ii) Prove that $M_1 = T(O_1)$ is dense in the ball $Q_{\delta_0}$ centered at the origin**:
    Let $\delta_0 = \frac{r_0}{n_0}$; define the ball centered at the origin in $E_1$: $Q_{\delta_0} := Q(\theta, \delta_0) = \{y \in E_1 : \|y\| \le \delta_0\}$.
    Take any $y \in Q_{\delta_0}$, then $y_0 \pm n_0 y \in Q(y_0, r_0)$. Hence there exist sequences $\{x_k\}$ and $\{x_k'\}$ in $O_{n_0}$ such that:

    \[
    Tx_k \rightarrow y_0 + n_0 y, \quad Tx_k' \rightarrow y_0 - n_0 y
    \]

    Then $T(x_k - x_k') \rightarrow 2n_0 y$, i.e.,

    \[
    M_1 \ni T\left(\frac{x_k - x_k'}{2n_0}\right) \rightarrow y
    \]

    because $\frac{x_k - x_k'}{2n_0} \in O_1$ (norm does not exceed $\frac{2n_0}{2n_0} = 1$). This shows that $M_1$ is dense in $Q_{\delta_0}$.

    **(iii) Prove that $M_1 \supseteq Q_{\delta_0/2}$**:
    Denote $O_{1/2^n} = \{x \in E : \|x\| \le 1/2^n\}$, $Q_{\delta_0/2^n} = \{y \in E_1 : \|y\| \le \delta_0/2^n\}$.
    By scaling from step (ii), $T(O_{1/2^n})$ is dense in $Q_{\delta_0/2^n}$.
    Take any $y \in Q_{\delta_0/2}$. Since $T(O_{1/2})$ is dense, there exists $x_1 \in O_{1/2}$ such that:

    \[
    \|y - Tx_1\| \le \frac{\delta_0}{2^2}
    \]

    i.e., $y - Tx_1 \in Q_{\delta_0/2^2}$. Since $T(O_{1/2^2})$ is dense in this ball, there exists $x_2 \in O_{1/2^2}$ such that:

    \[
    \|y - Tx_1 - Tx_2\| \le \frac{\delta_0}{2^3}
    \]

    By induction, we can prove that there exists a sequence $x_n \in O_{1/2^n}$ such that:

    \[
    \|y - T(x_1 + x_2 + \dots + x_n)\| \le \frac{\delta_0}{2^{n+1}}
    \]

    Since $E$ is a Banach space, the series $x = \sum_{n=1}^\infty x_n$ converges in $E$:

    \[
    \|x\| \le \sum_{n=1}^\infty \|x_n\| \le \sum_{n=1}^\infty \frac{1}{2^n} = 1
    \]

    Hence $x \in O_1$, and by continuity of the operator, $Tx = y$. Therefore $M_1 = T(O_1) \supseteq Q_{\delta_0/2}$.

    **(iv) Prove the conclusions of the theorem**:
    Take any $y \in E_1, y \ne \theta$. Let $y' = \frac{\delta_0}{2} \frac{y}{\|y\|} \in Q_{\delta_0/2}$.
    By (iii), there exists $x' \in O_1$ such that $Tx' = y'$. Let $x = \frac{2\|y\|}{\delta_0} x'$, then:

    \[
    Tx = y
    \]

    and

    \[
    \|x\| = \frac{2\|y\|}{\delta_0} \|x'\| \le \frac{2\|y\|}{\delta_0} = \frac{2}{\delta_0} \|Tx\|
    \]

    Thus we have proven $F = E_1$, and taking $K = \frac{2}{\delta_0}$ satisfies condition (ii) of the theorem. $\square$

---

### 2.2 Corollary of the Open Mapping Theorem

!!! success "Corollary 10.2 (Open Mapping)"

    Suppose a bounded linear operator $T$ satisfies the conditions of the above theorem. Then $T$ maps every **open set** in $E$ to an **open set** in $E_1$.

??? proof "Proof of Corollary 10.2 (Click to Expand)"

    From step (iii) of the proof of the Banach Open Mapping Theorem, we have $T(O_1) \supseteq Q_{\delta_0/2}$. By linear scaling, for any $n$, we have:

    \[
    T(O_{1/2^n}) \supseteq Q_{\delta_0/2^{n+1}}
    \]

    Now let $G \subset E$ be an open set. Take any $y \in T(G)$, then there exists $x \in G$ such that $Tx = y$. Since $G$ is open, $x$ is an interior point of $G$, so there exists $n$ such that $x + O_{1/2^n} \subset G$. Then:

    \[
    y + T(O_{1/2^n}) \subset T(G)
    \]

    Note that $Q_{\delta_0/2^{n+1}} \subseteq T(O_{1/2^n})$, therefore:

    \[
    y + Q_{\delta_0/2^{n+1}} \subset T(x + O_{1/2^n}) \subset T(G)
    \]

    This shows that $y$ is an interior point of $T(G)$, hence $T(G)$ is open. $\square$

*Note: Another corollary is that the range of a bounded linear operator $T: E \rightarrow E_1$ is either the whole $E_1$ or a set of the first category in $E_1$; one of the two must hold.*

---

### 2.3 Inverse Operator Theorem

If $T$ is bijective (injective and surjective), then the inverse mapping $T^{-1}$ of $T$ exists, called the **inverse operator** of $T$. The following theorem provides an important condition for the boundedness of the inverse operator.

!!! success "Theorem 10.3 (Inverse Operator Theorem)"

    Let a bounded linear operator $T$ map a Banach space $E$ onto a set of the second category in a Banach space $E_1$, and let $T$ be injective. Then $T$ possesses a **bounded** inverse operator.

??? proof "Proof of Theorem 10.3 (Click to Expand)"

    By the Banach Open Mapping Theorem, since the image is a set of the second category, the range of $T$ must be the whole space $E_1$.
    By assumption, $T$ is injective, so $T$ is bijective, and its inverse operator $T^{-1}$ exists.

    Write $x = T^{-1}y$. From conclusion (ii) of the Open Mapping Theorem, there exists $K$ such that:

    \[
    \|x\| \le K\|Tx\|
    \]

    Substituting $y$, we get:

    \[
    \|T^{-1}y\| = \|x\| \le K\|Tx\| = K\|y\|
    \]

    This directly proves that $T^{-1}$ is bounded. $\square$

---

## 3. Norm Equivalence

In normed spaces, studying the equivalence between different norms is crucial for understanding the topological structure of the space.

!!! info "Definition 10.4 (Norm Equivalence)"

    Let $E$ be a linear space, and $\|\cdot\|_1$ and $\|\cdot\|_2$ be two norms defined on $E$. If there exist positive numbers $K_1, K_2$ such that the inequality:

    \[
    K_1\|x\|_1 \le \|x\|_2 \le K_2\|x\|_1
    \]

    holds for all $x \in E$, then the norms $\|\cdot\|_1$ and $\|\cdot\|_2$ are said to be **equivalent**. Equivalent norms induce completely isomorphic topological structures.

Using the Inverse Operator Theorem, we can obtain a very strong conclusion on Banach spaces:

!!! success "Corollary 10.4"

    Let $(E, \|\cdot\|_1)$ and $(E, \|\cdot\|_2)$ both be Banach spaces. If there exists a positive number $K$ such that for all $x \in E$:

    \[
    \|x\|_2 \le K\|x\|_1
    \]

    then the two norms are equivalent, and hence the topologies are isomorphic.

??? proof "Proof of Corollary 10.4 (Click to Expand)"

    Let $I$ be the identity operator on $E$, which can be viewed as an operator from the Banach space $(E, \|\cdot\|_1)$ onto the Banach space $(E, \|\cdot\|_2)$. Clearly $I$ is bijective.
  
    From the inequality $\|Ix\|_2 = \|x\|_2 \le K\|x\|_1$, we see that the operator $I$ is bounded (continuous).

    By the **Inverse Operator Theorem**, the inverse operator $I^{-1}$ of $I$ is also bounded. That is, there exists $K' > 0$ such that:

    \[
    \|x\|_1 = \|I^{-1}x\|_1 \le K'\|x\|_2
    \]

    Combining the two inequalities, the norms $\|\cdot\|_1$ and $\|\cdot\|_2$ are equivalent. $\square$

---

## 4. Graph of an Operator and the Closed Graph Theorem

### 4.1 Graph of an Operator and Closed Operators

In mathematical analysis, the graph of a unary function $y = f(x)$ is a curve in the plane consisting of points $(x, f(x))$. For general linear operators, we can also introduce the graph.

Let $E, E_1$ be normed linear spaces. Form the direct sum $E \oplus E_1$, and define the norm:

\[
\|(x, y)\| = \|x\| + \|y\|
\]

It is easy to verify that $E \oplus E_1$ is a normed linear space with this norm (if both spaces are complete, then it is also a Banach space).

!!! info "Definition 10.5 (Closed Operator)"

    Let $T$ be a linear operator defined on a subspace $D$ of $E$ with range contained in $E_1$. The set of all elements of the form $(x, Tx)$ ($x \in D$) in $E \oplus E_1$ is called the **graph** of $T$, denoted by $G_T$.

    If the graph $G_T$ of $T$ is a **closed subspace** of $E \oplus E_1$, then $T$ is called a **closed linear operator** or **closed operator**.

!!! success "Property 10.5 (Equivalent Characterization of Closed Operators)"

    $T$ is a closed operator if and only if for any sequence $\{x_n\} \in D$, if $\{x_n\}$ and $\{Tx_n\}$ converge in $E$ and $E_1$ respectively to $x$ and $y$, then necessarily $x \in D$ and $Tx = y$.

??? proof "Proof of Property 10.5 (Click to Expand)"

    **Sufficiency**: Take any $(x,y) \in \overline{G_T}$. Then there exists $\{x_n\} \in D$ such that $(x_n, Tx_n) \rightarrow (x, y)$, i.e.,

    \[
    \|(x_n, Tx_n) - (x,y)\| = \|x_n - x\| + \|Tx_n - y\| \rightarrow 0
    \]

    Hence $x_n \rightarrow x$ and $Tx_n \rightarrow y$. By assumption, we have $x \in D$ and $Tx = y$, so $(x,y) = (x, Tx) \in G_T$. Therefore $G_T = \overline{G_T}$, and $T$ is closed.

    **Necessity**: Let $\{x_n\} \in D$, and suppose $x_n \rightarrow x$, $Tx_n \rightarrow y$. Then $\|(x_n, Tx_n) - (x,y)\| \rightarrow 0$. Since $G_T$ is closed, the limit must belong to $G_T$, i.e., $(x, y) \in G_T$. This implies $x \in D$ and $Tx = y$. $\square$

**Example: The differential operator is an unbounded closed operator**

The domain of the differential operator $T = \frac{d}{dt}$ is the set $C^1[a,b]$ of functions in $C[a,b]$ that have continuous derivatives.
If $\{x_n\} \in C^1[a,b]$, and in $C[a,b]$ we have $x_n \rightarrow x$ and $Tx_n \rightarrow y$. That is, the sequence of functions $x_n(t)$ and their derivatives $x_n'(t)$ converge **uniformly** to $x(t)$ and $y(t)$, respectively.

By a classical theorem in mathematical analysis, the limit function $x(t)$ must have a continuous derivative, and $x'(t) = y(t)$.
Hence $x \in C^1[a,b]$, and $Tx = y$. Therefore $T$ is a closed operator. We have already proven earlier that it is unbounded.

---

### 4.2 Closed Graph Theorem

A closed operator is not necessarily bounded, but if its domain is the whole Banach space, then it must be bounded. This leads to the **Closed Graph Theorem** in functional analysis.

!!! success "Theorem 10.6 (Closed Graph Theorem)"

    Let $T$ be a linear operator mapping a Banach space $E$ into a Banach space $E_1$ (note that the domain is the whole space $E$). Then $T$ is **bounded** if and only if $T$ is a **closed operator**.

??? proof "Proof of Theorem 10.6 (Click to Expand)"

    **Necessity**: Suppose $T$ is bounded. Take any $(x,y) \in \overline{G_T}$. Then there exists $\{x_n\} \in E$ such that $(x_n, Tx_n) \rightarrow (x, y)$, i.e.,

    \[
    x_n \rightarrow x, \quad Tx_n \rightarrow y
    \]

    Since $T$ is bounded (hence continuous), we have $Tx_n \rightarrow Tx$. By uniqueness of the limit, we must have $y = Tx$. Thus $(x,y) \in G_T$. This shows that $G_T$ is closed, i.e., $T$ is a closed operator.

    **Sufficiency**: Suppose $T$ is a closed operator. Then $G_T$ is a closed subspace of $E \oplus E_1$. Since $E$ and $E_1$ are Banach spaces, $E \oplus E_1$ is also a Banach space, and its closed subspace $G_T$ is likewise a Banach space.

    Define an operator $\tilde{T}$ from $G_T$ onto $E$ by projection:

    \[
    \tilde{T}(x, Tx) = x
    \]

    Clearly $\tilde{T}$ is linear. Since the domain is the whole space $E$, $\tilde{T}$ is surjective.
    If $\tilde{T}(x, Tx) = \theta$, then $x = \theta$, so $Tx = \theta$. That is, $(x, Tx) = (\theta, \theta)$, hence $\tilde{T}$ is injective. In summary, $\tilde{T}$ is bijective.

    Now check its boundedness:

    \[
    \|\tilde{T}(x, Tx)\| = \|x\| \le \|x\| + \|Tx\| = \|(x, Tx)\|
    \]

    This shows that $\tilde{T}$ is a bounded linear bijection. By the **Inverse Operator Theorem**, $\tilde{T}$ has a bounded inverse operator $\tilde{T}^{-1}$.

    For any $x \in E$, since $(x, Tx) = \tilde{T}^{-1}(x)$, we have:

    \[
    \|(x, Tx)\| \le \|\tilde{T}^{-1}\|\|x\|
    \]

    Because $\|Tx\| \le \|x\| + \|Tx\| = \|(x, Tx)\|$, we further obtain:

    \[
    \|Tx\| \le \|\tilde{T}^{-1}\|\|x\|
    \]

    This proves that the original operator $T$ is bounded. $\square$

**Significance of the Closed Graph Theorem**: It transforms the problem of determining whether a linear operator defined on a Banach space is bounded into the problem of determining whether the operator is closed. In practice, verifying the condition $x_n \rightarrow x$ and $Tx_n \rightarrow y \implies Tx = y$ is often much more convenient than directly proving boundedness. However, it must be noted that the premise is that the domain of the operator must be the complete whole space.