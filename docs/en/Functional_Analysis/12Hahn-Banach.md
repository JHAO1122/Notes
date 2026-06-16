# Chapter 12: Extension of Bounded Linear Functionals and the Uniform Boundedness Principle

This chapter focuses on the extension theory of bounded linear functionals, providing a detailed proof of one of the most fundamental cornerstone theorems in functional analysis—the **Hahn–Banach Extension Theorem** (including both the real and complex linear space versions)—and delving into a series of extremely important corollaries in normed linear spaces.

---

## 1. The Concept of Extension of Linear Functionals

Before discussing the extension of bounded linear functionals, a fundamental question confronts us: In a general normed linear space that contains only nonzero elements, does there always exist a nonzero bounded linear functional?

To answer this fundamental question about the “existence” of nonzero bounded linear functionals, we must thoroughly study the extension theory of bounded linear functionals.

!!! info "Definition 12.1 (Extension)"

    Let $E$ be a linear space, and let $f_1$ and $f_2$ be linear functionals defined on subspaces $G_1$ and $G_2$ of $E$, respectively. If they satisfy the following two conditions:

    (i) $G_1 \subset G_2$;

    (ii) For every $x \in G_1$, we have $f_1(x) = f_2(x)$;

    then $f_2$ is called an **extension** of $f_1$ to $G_2$.

---

## 2. Subadditive and Positively Homogeneous Functionals

To establish the mathematical tools for controlling the boundary conditions in the Hahn–Banach Theorem, we first introduce two special types of functionals.

!!! info "Definition 12.2 (Subadditive and Positively Homogeneous)"

    Let $E$ be a linear space, and let $p$ be a real‑valued functional defined on $E$.

    * If for all $x, y \in E$ we always have:
  
        \[
        p(x + y) \le p(x) + p(y)
        \]
  
        then $p$ is called **subadditive**.

    * If for every real number $a \ge 0$ and every $x \in E$ we always have:
  
        \[
        p(ax) = a p(x)
        \]
  
        then $p$ is called **positively homogeneous**.

---

## 3. The Hahn–Banach Extension Theorem (Real Space Version)

This is the most classical core algebraic version of the Hahn–Banach Theorem. It states that, under the control of a subadditive, positively homogeneous functional, a linear functional defined on a real linear subspace can always be extended to the whole space.

!!! success "Theorem 12.1 (Hahn–Banach Extension Theorem – Real Space Version)"

    Let $G$ be a subspace of a **real linear space** $E$, let $f$ be a real linear functional defined on $G$, and let $p$ be a subadditive, positively homogeneous functional defined on $E$.

    Suppose $f$ satisfies the domination condition: for every $x \in G$, $f(x) \le p(x)$.

    Then there exists a real linear functional $F$ defined on the whole space $E$ such that:

    (i) $F(x) = f(x)$ for all $x \in G$, i.e., $F$ is an extension of $f$ to $E$;

    (ii) $F(x) \le p(x)$ for all $x \in E$.

??? proof "Proof of Theorem 12.1 (Click to expand)"

    The proof uses a step‑by‑step one‑dimensional extension method, followed by transfinite induction (Zorn’s Lemma).

    **Step 1: Prove the possibility of a one‑dimensional extension**

    * If $G = E$, we simply take $F = f$ and the theorem is obviously true.

    * Assume $G \neq E$; then we can arbitrarily pick an element $x_0 \in E \setminus G$.

    * Let $G_1$ be the extended subspace spanned by the subspace $G$ together with the element $x_0$:
  
        \[
        G_1 := \{ \alpha x_0 + x \mid \alpha \in \mathbb{R}, x \in G \}
        \]
  
    * For any $x', x'' \in G$, using the linearity of $f$ and the subadditivity of $p$, we have:
  
        \[
        f(x'') - f(x') = f(x'' - x') \le p(x'' - x') = p((x'' + x_0) + (-x' - x_0)) \le p(x'' + x_0) + p(-x' - x_0)
        \]
  
    * Rearranging terms, we obtain:
  
        \[
        -p(-x' - x_0) - f(x') \le p(x'' + x_0) - f(x'')
        \]
  
    * The left‑hand side depends only on $x'$, the right‑hand side only on $x''$, and the inequality holds for all $x', x'' \in G$. Therefore the supremum of the left‑hand side does not exceed the infimum of the right‑hand side:
  
        \[
        c' = \sup_{x' \in G} \{ -p(-x' - x_0) - f(x') \} \le c'' = \inf_{x'' \in G} \{ p(x'' + x_0) - f(x'') \}
        \]
  
    * Choose an arbitrary real number $c$ in the closed interval $[c', c'']$ so that $c' \le c \le c''$. Now define a new functional $\tilde{f}$ on the extended space $G_1$:
  
        \[
        \tilde{f}(\alpha x_0 + x) := \alpha c + f(x)
        \]
  
    * Because every element in $G_1$ can be uniquely represented as $\alpha x_0 + x$, the above expression uniquely defines a linear functional $\tilde{f}$.

    **Step 2: Verify that $\tilde{f}$ satisfies the domination condition $\tilde{f}(\alpha x_0 + x) \le p(\alpha x_0 + x)$**

    * 1. If $\alpha = 0$, then $\tilde{f}(x) = f(x) \le p(x)$, and the domination condition obviously holds.

    * 2. If $\alpha > 0$, from the choice $c \le c''$ we obtain:
  
        \[
        c \le \inf_{x'' \in G} \{ p(x'' + x_0) - f(x'') \} \le p\left( \frac{x}{\alpha} + x_0 \right) - f\left( \frac{x}{\alpha} \right)
        \]
  
        Hence:
  
        \[
        c + f\left( \frac{x}{\alpha} \right) \le p\left( \frac{x}{\alpha} + x_0 \right)
        \]
  
        Multiplying both sides by the positive number $\alpha$ and using the linearity of $f$ and the positive homogeneity of $p$, we obtain:
  
        \[
        \alpha c + f(x) \le p(x + \alpha x_0)
        \]
  
    * 3. If $\alpha < 0$, from the choice $c' \le c$ we obtain:
  
        \[
        -p\left( -\frac{x}{\alpha} - x_0 \right) - f\left( -\frac{x}{\alpha} \right) \le c
        \]
  
        After changing signs and rearranging, we get:
  
        \[
        -c - f\left( \frac{x}{\alpha} \right) \le p\left( -\frac{x}{\alpha} - x_0 \right)
        \]
  
        Multiplying both sides by the positive number $-\alpha$, again using positive homogeneity, we obtain:
  
        \[
        \alpha c + f(x) \le p(x + \alpha x_0)
        \]
  
    * In summary, we have strictly proved that the linear functional can be successfully extended by one dimension, and the domination condition is preserved on the new space $G_1$.

    **Step 3: Complete the extension to the whole space using Zorn’s Lemma**

    * Let $\mathscr{F}$ denote the set of all algebraic extensions $F$ of $f$ that satisfy the domination condition $F(x) \le p(x)$ for all $x$ in the domain $D_F$ of $F$. Since Step 2 has shown that a valid extension to $G_1$ exists, the set $\mathscr{F}$ is **non‑empty**.

    * Define a partial order on $\mathscr{F}$: for any $F_1, F_2 \in \mathscr{F}$, if $F_2$ is an extension of $F_1$ (i.e., $D_{F_1} \subset D_{F_2}$ and they agree on $D_{F_1}$), then we say $F_2$ succeeds $F_1$ and write $F_1 \prec F_2$. With this partial order, $\mathscr{F}$ becomes a partially ordered set.

    * Let $\mathscr{M}$ be an arbitrary totally ordered subset of $\mathscr{F}$. Define a common domain:
  
        \[
        D = \bigcup_{F \in \mathscr{M}} D_F
        \]
  
    * Define a global functional $\varphi$ on $D$ as follows: for any $x \in D$, $x$ belongs to some $D_F$ (with $F \in \mathscr{M}$), set:
  
        \[
        \varphi(x) = F(x)
        \]
  
    * Because $\mathscr{M}$ is totally ordered, this definition is well‑defined. It is easy to verify that $\varphi$ is a linear functional with domain $D$ and satisfies $\varphi(x) \le p(x)$ for all $x \in D$.

    * Clearly $\varphi$ is an upper bound for every element in $\mathscr{M}$, so $\varphi \in \mathscr{F}$. The partially ordered set $\mathscr{F}$ satisfies the conditions of Zorn’s Lemma; therefore $\mathscr{F}$ must possess a maximal element. Denote this maximal element by $F_0$.

    * The domain $D_{F_0}$ of the maximal element $F_0$ must be the whole space $E$.

    * If $D_{F_0} \neq E$, then we could, by the method of Steps 1 and 2, extend $F_0$ by one further dimension to some $E_1$, obtaining $F_1$ that still satisfies the domination condition. This would give $F_1 \in \mathscr{F}$ with $F_0 \prec F_1$ ($F_0 \neq F_1$), contradicting the maximality of $F_0$.

    * Hence we must have $D_{F_0} = E$. Taking $F = F_0$ completes the proof of the real linear extension to the whole space. $\square$

---

## 4. The Hahn–Banach Extension Theorem (Complex Space Version)

Using the fact that a complex functional is completely determined by its real part, we can perfectly extend the theorem to complex linear spaces (also known as the Bohnenblust–Sobczyk Theorem).

!!! success "Theorem 12.2 (Hahn–Banach Extension Theorem – Complex Space Version / General Norm Version)"

    Let $G$ be a subspace of a complex (or real) linear space $E$, let $p$ be a seminorm on $E$ (satisfying subadditivity and absolute homogeneity: $p(\alpha x) = |\alpha|\,p(x)$). Let $f$ be a bounded linear functional defined on $G$ that satisfies the domination condition:
  
    \[
    |f(x)| \le p(x) \quad (\forall x \in G)
    \]
  
    Then there exists a linear functional $F$ defined on the whole space $E$ such that:

    (i) $F$ is an extension of $f$ to $E$, i.e., $F(x) = f(x)$ for all $x \in G$;

    (ii) The domination condition is preserved for every $x \in E$:
  
        \[
        |F(x)| \le p(x)
        \]

??? proof "Proof of Theorem 12.2 – Complex Space Version (Click to expand)"

    We focus on the details for complex linear spaces.

    **Step 1: Decompose the real part**

    * Let $f$ be a complex linear functional. Decompose it into its real and imaginary parts: $f(x) = \varphi_0(x) + i \psi_0(x)$, where $\varphi_0, \psi_0$ are real‑valued functionals.

    * Using complex linearity $f(ix) = i f(x)$, we expand:
  
        \[
        \varphi_0(ix) + i \psi_0(ix) = i (\varphi_0(x) + i \psi_0(x)) = -\psi_0(x) + i \varphi_0(x)
        \]
  
    * Comparing real and imaginary parts gives $\psi_0(x) = -\varphi_0(ix)$. Hence a complex linear functional is completely determined by its real part:
  
        \[
        f(x) = \varphi_0(x) - i \varphi_0(ix)
        \]

    **Step 2: Extend the real part**

    * Clearly $\varphi_0$ is a real linear functional on the real linear space $G$ and satisfies $\varphi_0(x) = \operatorname{Re} f(x) \le |f(x)| \le p(x)$.

    * By the already proved real version of the Hahn–Banach Theorem (Theorem 12.1), we can extend $\varphi_0$ to a real linear functional $\varphi$ on the whole space $E$, with the domination condition on $E$:
  
        \[
        \varphi(x) \le p(x)
        \]

    **Step 3: Reconstruct the complex functional and verify linearity**

    * Following the form of Step 1, we reconstruct a complex functional $F$ on the whole space $E$:
  
        \[
        F(x) = \varphi(x) - i \varphi(ix)
        \]
  
    * We verify that $F$ satisfies the properties of a complex linear operator. For any $x \in E$:
  
        \[
        F(ix) = \varphi(ix) - i \varphi(-x) = \varphi(ix) + i \varphi(x) = i [ \varphi(x) - i \varphi(ix) ] = i F(x)
        \]
  
    * Together with the real linearity of $\varphi$, it is easy to verify that $F$ is complex linear. When $x \in G$, because $\varphi$ extends $\varphi_0$, we clearly have $F(x) = f(x)$.

    **Step 4: Verify the absolute value domination condition**

    * Take any $x \in E$. Write the complex number $F(x)$ in polar form (argument $\alpha$): $F(x) = |F(x)| e^{i\alpha}$.

    * Then clearly $e^{-i\alpha} F(x) = |F(x)|$ is a non‑negative real number. Using the complex linearity of $F$ and the definition, we have:
  
        \[
        |F(x)| = e^{-i\alpha} F(x) = F(e^{-i\alpha}x) = \operatorname{Re} F(e^{-i\alpha}x) = \varphi(e^{-i\alpha}x)
        \]
  
    * Now apply the domination condition for $\varphi$ from the real version, together with the absolute homogeneity of the seminorm $p$, to obtain:
  
        \[
        |F(x)| \le p(e^{-i\alpha}x) = |e^{-i\alpha}| \, p(x) = p(x)
        \]
  
    * The absolute value domination condition is preserved on the whole space, completing the proof of the complex version. $\square$

---

## 5. Important Corollaries in Normed Linear Spaces

In a normed linear space, if we take $p(x) = \|f\|_G \cdot \|x\|$, the Hahn–Banach Theorem directly becomes the **norm‑preserving extension theorem** for bounded linear functionals. From this we can derive a series of extremely important geometric and topological properties.

### Corollary 12.3 (Norm‑preserving Extension)

Let $G$ be a subspace of a normed linear space $E$, and let $f$ be a bounded linear functional defined on $G$. Then there exists a bounded linear functional $F$ defined on the whole space $E$ such that $F$ is an extension of $f$ and their operator norms are exactly equal:

\[
\|F\|_E = \|f\|_G
\]

---

### Corollary 12.4 (Existence of Nonzero Functionals)

Let $E$ be a normed linear space, and let $x_0 \in E$ be any nonzero element ($x_0 \neq \theta$). Then there exists a bounded linear functional $f_0$ in the dual space $E^*$ that satisfies the following properties:

(i) $\|f_0\| = 1$;

(ii) $f_0(x_0) = \|x_0\|$.

This perfectly answers the question about the existence of nonzero bounded linear functionals posed at the beginning of this chapter.

---

### Corollary 12.5 (Dual Representation of the Norm of an Element)

For any element $x$ in a normed linear space $E$, its norm can be equivalently expressed as the supremum of the values of functionals running over the unit ball of the dual space:

\[
\|x\| = \sup_{f \in E^*,\, \|f\|=1} |f(x)|
\]

---

### Corollary 12.6 (Separation Theorem for a Point and a Subspace)

Let $G$ be a subspace of a normed linear space $E$, and suppose the distance from a point $x_0 \in E$ to the subspace $G$ is positive:

\[
d = \inf_{x \in G} \|x_0 - x\| > 0
\]

Then there exists a bounded linear functional $f$ in the dual space $E^*$ that satisfies the following separation properties:

(i) $\|f\| = 1$;

(ii) $f(x) = 0$ for all $x \in G$ (i.e., $f \in G^\perp$);

(iii) $f(x_0) = d$.

## 6. The Uniform Boundedness Principle and Its Applications

In functional analysis, besides the Open Mapping Theorem and the Closed Graph Theorem, there is another extremely fundamental theorem that deeply reveals the intrinsic relationship between pointwise boundedness and uniform boundedness of families of operators. This theorem is usually called the **Uniform Boundedness Principle**, and in classical literature it is also often called the **Resonance Theorem**.

### 6.1 The Uniform Boundedness Principle (Resonance Theorem)

!!! info "Theorem 12.3 (Uniform Boundedness Principle / Resonance Theorem)"

    Let $X$ be a Banach space and $Y$ a normed linear space.
  
    Let $\{T_\alpha\}_{\alpha \in I}$ be a family of bounded linear operators from $X$ to $Y$ (where $I$ is an arbitrary index set).
  
    If for every fixed $x \in X$ the image set of this operator family at $x$ is bounded in $Y$, i.e., there exists a constant $C_x > 0$ such that:
  
    \[
    \sup_{\alpha \in I} \|T_\alpha x\| \le C_x < \infty
    \]
  
    Then the operator family must be uniformly bounded in the operator norm sense. That is, there exists a uniform constant $M > 0$ such that for all $\alpha \in I$:
  
    \[
    \sup_{\alpha \in I} \|T_\alpha\| \le M < \infty
    \]

??? proof "Proof of the Uniform Boundedness Principle (Click to expand)"

    The proof of this theorem relies heavily on the **Baire Category Theorem** for complete metric spaces.

    For each positive integer $n$, define the set:
  
    \[
    X_n = \left\{ x \in X \;\middle|\; \|T_\alpha x\| \le n,\quad \forall \alpha \in I \right\}
    \]
  
    First, we show that $X_n$ is closed. Since each operator $T_\alpha$ is continuous (bounded linear operator), the function $x \mapsto \|T_\alpha x\|$ is also continuous. Therefore the set $\{x \in X \mid \|T_\alpha x\| \le n\}$ is closed.
  
    And $X_n$ is precisely the intersection of these closed sets over all $\alpha \in I$:
  
    \[
    X_n = \bigcap_{\alpha \in I} \left\{ x \in X \;\middle|\; \|T_\alpha x\| \le n \right\}
    \]
  
    An arbitrary intersection of closed sets is still closed, hence $X_n$ must be a closed set in $X$.
  
    By the hypothesis of the theorem, for any element $x$ in $X$, the set $\{ \|T_\alpha x\| \}_{\alpha \in I}$ is bounded. This implies that there exists some sufficiently large integer $n$ such that $\|T_\alpha x\| \le n$ for all $\alpha \in I$.
  
    This shows that every point of $X$ must belong to some $X_n$, i.e.:
  
    \[
    X = \bigcup_{n=1}^\infty X_n
    \]
  
    Since $X$ is a Banach space, it is a complete metric space. By the **Baire Category Theorem**, a complete space can never be represented as a countable union of nowhere dense sets.
  
    Therefore, among the sets $\{X_n\}_{n=1}^\infty$, at least one closed set $X_{n_0}$ is not nowhere dense. In other words, $X_{n_0}$ must have non‑empty interior.
  
    Since $X_{n_0}$ has non‑empty interior, there must exist a point $x_0 \in X_{n_0}$ and a radius $r > 0$ such that the closed ball centered at $x_0$ with radius $r$ is entirely contained in $X_{n_0}$:
  
    \[
    \overline{B}(x_0, r) = \{ x \in X \mid \|x - x_0\| \le r \} \subseteq X_{n_0}
    \]
  
    For any point $z \in \overline{B}(x_0, r)$, by the definition of $X_{n_0}$ we have, for all $\alpha \in I$:
  
    \[
    \|T_\alpha z\| \le n_0
    \]
  
    Now, take any element $x \in X$ with $\|x\| \le 1$.
  
    Construct the element $z = x_0 + r x$. Clearly $\|z - x_0\| = \|r x\| = r \|x\| \le r$, so $z \in \overline{B}(x_0, r) \subseteq X_{n_0}$.
  
    Since $x_0 \in X_{n_0}$ and $z \in X_{n_0}$, we have, for any $\alpha \in I$, $\|T_\alpha x_0\| \le n_0$ and $\|T_\alpha z\| \le n_0$.
  
    Using the linearity of $T_\alpha$ and the triangle inequality for norms, we can estimate $T_\alpha x$:
  
    \[
    \|T_\alpha(r x)\| = \|T_\alpha(z - x_0)\| = \|T_\alpha z - T_\alpha x_0\| \le \|T_\alpha z\| + \|T_\alpha x_0\| \le n_0 + n_0 = 2n_0
    \]
  
    Because $\|T_\alpha(r x)\| = r \|T_\alpha x\|$, we obtain:
  
    \[
    r \|T_\alpha x\| \le 2n_0
    \]
  
    Hence:
  
    \[
    \|T_\alpha x\| \le \frac{2n_0}{r}
    \]
  
    This inequality holds for all $x$ with $\|x\| \le 1$ and for all $\alpha \in I$. By the definition of the operator norm:
  
    \[
    \|T_\alpha\| = \sup_{\|x\| \le 1} \|T_\alpha x\| \le \frac{2n_0}{r}
    \]
  
    Taking $M = \frac{2n_0}{r}$, we obtain $\sup_{\alpha \in I} \|T_\alpha\| \le M < \infty$. This completes the proof of the Uniform Boundedness Principle. $\square$

### 6.2 Important Applications of the Uniform Boundedness Principle

The Uniform Boundedness Principle reveals the bridge between "local (pointwise) properties" and "global (uniform) properties", and it has extremely wide applications in functional analysis.

#### 6.2.1 Strong Convergence of Operator Sequences (Banach–Steinhaus Theorem)

!!! info "Theorem 12.4 (Banach–Steinhaus Theorem)"

    Let $X$ be a Banach space and $Y$ a normed linear space.
  
    Let $\{T_n\}_{n=1}^\infty$ be a sequence of bounded linear operators from $X$ to $Y$. If for every point $x \in X$, the sequence $\{T_n x\}$ converges in $Y$ (i.e., the operator sequence converges strongly), then:
  
    (i) There exists a uniform bound $M > 0$ such that $\sup_{n} \|T_n\| \le M$.
  
    (ii) Define the limit operator $T(x) = \lim_{n \to \infty} T_n x$; then $T$ must also be a **bounded linear operator**, i.e., $T \in \mathcal{B}(X, Y)$.

??? proof "Proof (Click to expand)"

    * Because for each $x \in X$ the sequence $\{T_n x\}$ converges, a convergent sequence is necessarily bounded. Hence there exists a constant $C_x$ such that $\sup_n \|T_n x\| \le C_x$.
  
    * $X$ is a Banach space, which exactly satisfies the conditions of the Uniform Boundedness Principle. Applying the Uniform Boundedness Principle directly yields $\sup_n \|T_n\| \le M < \infty$, which is conclusion (i).

    * For conclusion (ii), first, by the linearity of limits, the limit operator $T$ obviously satisfies additivity and homogeneity; therefore $T$ is a linear operator.

    * Next we prove the boundedness of $T$. For any $x \in X$, by the continuity of the norm, we have:
  
    \[
    \|T x\| = \left\| \lim_{n \to \infty} T_n x \right\| = \lim_{n \to \infty} \|T_n x\|
    \]
  
    * Since $\|T_n x\| \le \|T_n\| \|x\| \le M \|x\|$ for all $n$, taking the limit naturally yields:
  
    \[
    \|T x\| \le M \|x\|
    \]
  
    * This perfectly proves that the operator $T$ is also bounded, and its norm satisfies $\|T\| \le M$. $\square$

#### 6.2.2 Equivalence of Weakly Bounded Sets and Strongly Bounded Sets

In a normed linear space, the boundedness of a subset can be viewed from two different topological perspectives.

!!! info "Theorem 12.5"

    Let $E$ be a normed linear space and $A$ a subset of $E$.
  
    If $A$ is **weakly bounded** (i.e., for every bounded linear functional $f$ in the dual space $E^*$, the set $\{f(x) \mid x \in A\}$ is bounded), then $A$ must be **strongly bounded** (i.e., there exists a constant $M$ such that $\|x\| \le M$ for all $x \in A$).

??? proof "Proof (Click to expand)"

    * To apply the Uniform Boundedness Principle, we reverse the roles, viewing elements of the space as functionals. Using the canonical embedding, for each $x \in A$ we define a functional $x^{**}$ in the bidual space $E^{**}$:
  
    \[
    x^{**}(f) = f(x), \quad \forall f \in E^*
    \]
  
    * By the hypothesis of the theorem, for any fixed $f \in E^*$, the set of numbers $\{ f(x) \mid x \in A \}$ is bounded. This is equivalent to saying: for every $f \in E^*$, there exists a constant $C_f$ such that:
  
    \[
    \sup_{x \in A} |x^{**}(f)| \le C_f < \infty
    \]
  
    * We note that the dual space $E^*$ is always a complete **Banach space**, and $\{x^{**}\}_{x \in A}$ is a family of bounded linear operators defined on $E^*$.
  
    * This perfectly satisfies the conditions of the Uniform Boundedness Principle. Applying the theorem, this family of functionals must be uniformly bounded in the operator norm sense, i.e., there exists a constant $M$ such that:
  
    \[
    \sup_{x \in A} \|x^{**}\| \le M
    \]
  
    * Recalling the property of the canonical embedding, we have already proved $\|x^{**}\| = \|x\|$. Therefore:
  
    \[
    \sup_{x \in A} \|x\| \le M
    \]
  
    * This shows that the subset $A$ is bounded in the norm of the original space, i.e., $A$ is strongly bounded. $\square$

#### 6.2.3 Divergence of Fourier Series of Continuous Functions

This is one of the most famous and striking applications of the Uniform Boundedness Principle in classical analysis: it provides a non‑constructive proof that there exists a continuous function whose Fourier series diverges at some point.

!!! info "Theorem 12.6"

    In the space $C_{2\pi}$ of continuous functions with period $2\pi$, there must exist a continuous function $f$ such that its Fourier series diverges (does not converge) at the point $t=0$.

??? proof "Proof Outline (Click to expand)"

    * Let $X = C_{2\pi}$ be the Banach space of all continuous functions with period $2\pi$, equipped with the supremum norm $\|f\|_\infty = \max_{t} |f(t)|$.
  
    * Let $Y = \mathbb{R}$. For any $f \in C_{2\pi}$, the value at $t=0$ of the $n$-th partial sum of its Fourier series can be expressed by the Dirichlet integral:
  
    \[
    S_n(f, 0) = \frac{1}{2\pi} \int_{-\pi}^{\pi} f(t) D_n(t) dt
    \]
  
    * where $D_n(t) = \frac{\sin((n+1/2)t)}{\sin(t/2)}$ is the Dirichlet kernel.
  
    * Define a sequence of bounded linear functionals $T_n : C_{2\pi} \to \mathbb{R}$ by:
  
    \[
    T_n(f) = S_n(f, 0)
    \]
  
    * By a precise calculation, one can prove that the operator norm of $T_n$ equals the $L^1$ norm of the Dirichlet kernel:
  
    \[
    \|T_n\| = \frac{1}{2\pi} \int_{-\pi}^{\pi} |D_n(t)| dt
    \]
  
    * A classical result in mathematical analysis shows that the integral of the Dirichlet kernel is unbounded (its growth is logarithmic), i.e., as $n \to \infty$:
  
    \[
    \|T_n\| \to \infty
    \]
  
    * We argue by contradiction. Suppose that for every continuous function $f \in C_{2\pi}$, its Fourier series converges at $t=0$. Then for each given $f$, the sequence $\{T_n(f)\}$ would be convergent, and hence necessarily bounded:
  
    \[
    \sup_n |T_n(f)| < \infty \quad (\forall f \in C_{2\pi})
    \]
  
    * Since $X$ is a Banach space and the pointwise boundedness condition holds, by the Uniform Boundedness Principle there must exist a uniform constant $M$ such that:
  
    \[
    \sup_n \|T_n\| \le M < \infty
    \]
  
    * This contradicts the previously established fact that $\|T_n\| \to \infty$ in an irreconcilable way!
  
    * The contradiction shows that the original assumption must be false. Therefore, there exists at least one continuous function $f \in C_{2\pi}$ for which the sequence $\{S_n(f, 0)\}$ is unbounded, and hence diverges at $t=0$. $\square$
