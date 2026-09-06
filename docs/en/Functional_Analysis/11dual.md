# Chapter 11: Dual Spaces

This chapter mainly discusses the dual spaces of normed linear spaces, and provides detailed proofs of the Riesz Representation Theorem for Hilbert spaces and the dual spaces of $L^p$ spaces. Subsequently, it introduces and proves the extremely important Hahn-Banach Theorem (the extension theorem for bounded linear functionals) and the Uniform Boundedness Principle in functional analysis, along with their related applications.

---

## 1. Dual Spaces and Bidual Spaces

### 1.1 Definition of the Dual Space

!!! info "Definition 11.1 (Dual Space)"

    The set of all **bounded linear functionals** on a normed linear space $E$ forms a normed linear space. This space is called the **dual space** (or conjugate space) of $E$, denoted by $E^*$.

    Since the field $\mathbb{K}$ of real (or complex) numbers is complete and $E^* = \mathbb{B}(E, \mathbb{K})$, and because the codomain $\mathbb{K}$ is complete, regardless of whether $E$ itself is complete, $E^*$ as a space of operators is always complete. That is, $E^*$ is always a **Banach space**.

    **Definition of linear operations in $E^*$**: Let $f_1, f_2 \in E^*$. For any $x \in E$:

    \[
    (f_1 + f_2)(x) = f_1(x) + f_2(x)
    \]

    \[
    (\alpha f)(x) = \alpha f(x)
    \]

    **Definition of the norm**:

    \[
    \|f\| = \sup_{x \in E, x \neq \theta} \frac{|f(x)|}{\|x\|} = \sup_{\|x\|=1} |f(x)|
    \]

### 1.2 Bidual Spaces and Reflexive Spaces

Since $E^*$ itself is a Banach space, $E^*$ can also define its own dual space, i.e., $(E^*)^*$, which is called the **bidual space** of $E$ and is denoted by $E^{**}$. By analogy, we can successively define the third dual space $E^{***}$, and so on.

Let $x \in E$, $f \in E^*$. In the original definition of the dual space, we usually regard the functional $f$ as a given mapping, and let the variable $x$ run through the entire space $E$.
Conversely, if we fix an element $x$ in the linear space and let the functional $f$ run through the entire dual space $E^*$, then the expression $f(x)$ becomes a functional defined on $E^*$, which we denote by $x^{**}$. That is, for any $f \in E^*$:

\[
x^{**}(f) = f(x)
\]

!!! success "Theorem 11.1 (Isometric Isomorphism of the Canonical Mapping)"

    The mapping $E \rightarrow E^{**}: x \mapsto x^{**}$ has the following properties:
  
    (i) The mapping is **linear**, i.e., if $x_1 \mapsto x_1^{**}$, $x_2 \mapsto x_2^{**}$, then for any real (or complex) numbers $\alpha, \beta$:
  
    \[
    \alpha x_1 + \beta x_2 \mapsto \alpha x_1^{**} + \beta x_2^{**}
    \]
  
    (ii) The mapping is **isometric**, hence it is necessarily **injective**.

??? proof "Proof of Theorem 11.1 (Click to expand)"

    **Step 1: Show that $x^{**}$ is a bounded linear functional (i.e., $x^{**} \in E^{**}$)**

    * Linearity: From the definitions of addition and scalar multiplication of functionals in the dual space, $(f_1 + f_2)(x) = f_1(x) + f_2(x)$ and $(\alpha f)(x) = \alpha f(x)$, we have:
  
    \[
    x^{**}(f_1 + f_2) = (f_1 + f_2)(x) = f_1(x) + f_2(x) = x^{**}(f_1) + x^{**}(f_2)
    \]
  
    \[
    x^{**}(\alpha f) = (\alpha f)(x) = \alpha f(x) = \alpha x^{**}(f)
    \]

    * Boundedness: By the definition of the norm of a bounded linear functional:
  
    \[
    |x^{**}(f)| = |f(x)| \le \|x\| \|f\|
    \]
  
    This shows that when $x$ is fixed, this operator is bounded, and its operator norm satisfies $\|x^{**}\| \le \|x\|$.

    **Step 2: Prove property (i) the mapping is linear**

    * For any $f \in E^*$, expanding directly gives:
  
    \[
    (\alpha x_1^{**} + \beta x_2^{**})(f) = \alpha x_1^{**}(f) + \beta x_2^{**}(f) = \alpha f(x_1) + \beta f(x_2)
    \]
  
    \[
    = f(\alpha x_1) + f(\beta x_2) = f(\alpha x_1 + \beta x_2) = (\alpha x_1 + \beta x_2)^{**}(f)
    \]
  
    * Since this holds for all $f \in E^*$, by the equality of functionals we obtain:
  
    \[
    \alpha x_1^{**} + \beta x_2^{**} = (\alpha x_1 + \beta x_2)^{**}
    \]

    **Step 3: Prove property (ii) the mapping is isometric**

    * In Step 1 we already established the inequality in one direction:
  
    \[
    \|x^{**}\| \le \|x\|
    \]
  
    * On the other hand, by an important corollary of the Hahn-Banach Theorem, for a given non-zero element $x \in E$, there exists a bounded linear functional $f_0 \in E^*$ such that:
  
    \[
    \|f_0\| = 1 \quad \text{and} \quad f_0(x) = \|x\|
    \]
  
    * Then, substituting this special functional $f_0$ into the definition of $x^{**}$, we have:
  
    \[
    \|x^{**}\| = \sup_{\|f\|=1} |x^{**}(f)| \ge |x^{**}(f_0)| = |f_0(x)| = \|x\|
    \]
  
    * Combining the two inequalities, we obtain:
  
    \[
    \|x^{**}\| = \|x\|
    \]
  
    Since this mapping perfectly preserves the norm of space elements, it is necessarily injective. $\square$

!!! info "Definition 11.2 (Reflexive Space)"

    The mapping $x \mapsto x^{**}$ is an isometric embedding from $E$ into $E^{**}$. In functional analysis, it is usually called the **canonical mapping** from $E$ to $E^{**}$.
  
    * The canonical mapping is always isometric and isomorphic, but it is not necessarily surjective. When it is surjective (i.e., the image of the canonical mapping fills the entire bidual space $E^{**}$), we call the space $E$ a **reflexive space**.
  
    * We can regard $E$ and its image under the canonical mapping as completely identical. Therefore, up to isometric isomorphism, $E$ can be viewed as a subspace of $E^{**}$. From this perspective, if $E = E^{**}$, the space is called reflexive.

---

## 2. Bounded Linear Functionals on Hilbert Spaces

In a Hilbert space $H$, due to the presence of the inner product structure, bounded linear functionals have a very intuitive and elegant algebraic representation, given by the Riesz Representation Theorem.

!!! success "Theorem 11.2 (Riesz Representation Theorem)"

    Let $f$ be a bounded linear functional on a Hilbert space $H$. There exists a **unique** $u \in H$ such that:

    \[
    f(x) = (x, u)
    \]

    and the isometric norm relation $\|f\| = \|u\|$ holds.

    Conversely, for any given element $u \in H$, the equality $f(x) = (x, u)$ uniquely determines a bounded linear functional on $H$, and its operator norm equals $\|u\|$.

??? proof "Proof of the Riesz Representation Theorem (Click to expand)"

    **1. Proof of Existence:**

    * Let $f = 0$ be the zero functional. Then obviously we only need to take $u = \theta$, and the equality and norm relation hold trivially. Below we assume $f \neq 0$.

    * Let the null space of the functional be:
  
    \[
    L := \{x \in H \mid f(x) = 0\}
    \]
  
    * Since the functional $f$ is linear and continuous, its kernel $L$ is a proper closed subspace of $H$. By the orthogonal decomposition theorem of Hilbert spaces, there exists a non-zero element in its orthogonal complement, i.e., there exists a non-zero $x_0 \in L^\perp$. For convenience, we normalize it: let $\|x_0\| = 1$. Then the entire space can be decomposed as:
  
    \[
    H = L \oplus \mathbb{K}\{x_0\}
    \]
  
    * Let the constant $\alpha = f(x_0)$. For any element $x \in H$, we construct the following element:
  
    \[
    x - \frac{f(x)}{\alpha} x_0
    \]
  
    * Applying the functional $f$ to it, using its linearity:
  
    \[
    f\left(x - \frac{f(x)}{\alpha} x_0\right) = f(x) - \frac{f(x)}{\alpha} f(x_0) = f(x) - \frac{f(x)}{\alpha} \cdot \alpha = 0
    \]
  
    * This means the element $x - \frac{f(x)}{\alpha} x_0$ belongs to the null space $L$. Hence, any $x$ can be uniquely decomposed as the sum of two parts:
  
    \[
    x = \left(x - \frac{f(x)}{\alpha} x_0\right) + \frac{f(x)}{\alpha} x_0
    \]
  
    * Consider taking the inner product of both sides with the previously chosen special element $x_0$. Since the first part belongs to $L$ and $x_0 \in L^\perp$, their inner product is zero:
  
    \[
    (x, x_0) = \left(x - \frac{f(x)}{\alpha} x_0, x_0\right) + \left(\frac{f(x)}{\alpha} x_0, x_0\right) = 0 + \frac{f(x)}{\alpha} \|x_0\|^2 = \frac{f(x)}{\alpha}
    \]
  
    * From this we solve for the expression of the functional value:
  
    \[
    f(x) = \alpha (x, x_0)
    \]
  
    * By the conjugate linearity of the inner product with respect to the second argument, we can absorb the constant into the inner product:
  
    \[
    f(x) = (x, \overline{\alpha} x_0)
    \]
  
    * Now, simply set $u = \overline{\alpha} x_0 = \overline{f(x_0)} x_0$, and we have found an element $u \in H$ that represents the functional $f$. Existence is proved.

    **2. Proof of Uniqueness:**

    * Suppose there exists another element $u' \in H$ such that for all $x \in H$, we have the representation $f(x) = (x, u')$.

    * Then we have:
  
    \[
    (x, u) = (x, u') \implies (x, u - u') = 0
    \]
  
    * This must hold for every $x \in H$. This means the vector $u - u'$ is orthogonal to the entire space $H$.

    * Choose the specific $x = u - u'$ and substitute it in, obtaining $\|u - u'\|^2 = 0$. By the positivity of the inner product, we must have $u - u' = \theta$, i.e., $u = u'$. Uniqueness is proved.

    **3. Proof of the Isometric Norm Relation $\|f\| = \|u\|$:**

    * On one hand, by the Cauchy-Schwarz inequality, for any $x \in H$:
  
    \[
    |f(x)| = |(x, u)| \le \|x\| \|u\|
    \]
  
    By the definition of the supremum norm of a bounded linear functional, we directly obtain $\|f\| \le \|u\|$.

    * On the other hand, choose the specific element $x = u$ in the space and substitute it in:
  
    \[
    |f(u)| = |(u, u)| = \|u\|^2
    \]
  
    If $u \neq \theta$, then by the definition of the norm:
  
    \[
    \|f\| = \sup_{x \neq \theta} \frac{|f(x)|}{\|x\|} \ge \frac{|f(u)|}{\|u\|} = \frac{\|u\|^2}{\|u\|} = \|u\|
    \]
  
    * Combining both inequalities, we establish the strict isometric relation between the operator norm and the vector norm:
  
    \[
    \|f\| = \|u\|
    \]

    * Conversely, if an element $u \in H$ is given, using the properties of the inner product and the Cauchy-Schwarz inequality, it is easy to verify that $f(x) = (x, u)$ does satisfy boundedness and linearity. The proof of the norm equality is exactly as above. $\square$

### 2.1 Self-Duality of Hilbert Spaces

According to the Riesz Representation Theorem, we can establish a mapping from a Hilbert space $H$ to its dual space $H^*$: $u \mapsto f = (\cdot, u)$.
Suppose another element $v \mapsto g = (\cdot, v)$. Then, by the algebraic properties of the inner product, clearly:

\[
u + v \mapsto f + g = (\cdot, u) + (\cdot, v) = (\cdot, u + v)
\]

For any given real or complex scalar $\alpha$, since the inner product is conjugate linear in the second argument, we have:

\[
\alpha u \mapsto \overline{\alpha} f
\]

This is because if we map the element $\alpha u$ to a functional $h$, then for any $x \in H$, its functional value expands as:

\[
h(x) = (x, \alpha u) = \overline{\alpha}(x, u) = \overline{\alpha}f(x)
\]

Similar to the basic definitions of inner product spaces:

* When $H$ is a real space, the mapping $u \mapsto f$ is completely linear;

* When $H$ is a complex space, the mapping $u \mapsto f$ satisfies the complex conjugate rule, which we call **conjugate linear** (or antilinear).

From the isometric norm relation $\|f\| = \|u\|$ established in the Riesz Representation Theorem, this mapping constitutes an isometric isomorphism (or isometric conjugate isomorphism) from the space onto its dual space.

We can naturally induce an inner product structure on the dual space $H^*$: for any functionals $f$ (originating from $u$) and $g$ (originating from $v$) obtained via the mapping, define their inner product as:

\[
(f, g)_{H^*} := \overline{(u, v)}_H = (v, u)_H
\]

Under this definition of the inner product, $H^*$ itself perfectly forms a Hilbert space. Thus, the two Hilbert spaces $H$ and $H^*$ can be perfectly identified through this isometric isomorphism. We regard them as the same space, i.e., $\mathcal{H} = \mathcal{H}^*$. This extremely important property is called the **self-duality** or **self-conjugacy** of Hilbert spaces.

---

## 3. Dual Spaces of $L^p[a,b]$ Spaces

Besides the complete inner product space $L^2[a,b]$, functional analysis can also give a very definite and concrete algebraic representation theorem for the dual spaces of general Lebesgue integrable function spaces $L^p[a,b]$.

!!! success "Theorem 11.3"

    Let $1 < p < \infty$. Then the dual space of $L^p[a,b]$ can be isometrically identified with its conjugate exponent space:
  
    \[
    (L^p[a,b])^* = L^q[a,b]
    \]
  
    Here $p$ and $q$ satisfy the classic conjugate exponent relation: $\frac{1}{p} + \frac{1}{q} = 1$.

??? proof "Proof of Theorem 11.3 (Click to expand)"

    **Step 1: Construct the operator mapping $T$ from $L^q[a,b]$ to $(L^p[a,b])^*$**

    * Take any function $y \in L^q[a,b]$. For any given function $x \in L^p[a,b]$, define a functional by the integral expression:
  
    \[
    f(x) = \int_a^b x(t)y(t) dt
    \]

    * By Hölder's inequality from measure theory, the absolute value of the integral satisfies:
  
    \[
    \int_a^b |x(t)y(t)| dt \le \left(\int_a^b |x(t)|^p dt\right)^{\frac{1}{p}} \left(\int_a^b |y(t)|^q dt\right)^{\frac{1}{q}} = \|x\|_p \|y\|_q
    \]
  
    This shows that for all $x \in L^p[a,b]$, this integral is bounded and well-defined.

    * From the basic properties of integrals, the functional $f$ is linear on $L^p$. Moreover, the above inequality shows $|f(x)| \le \|y\|_q \|x\|_p$, which further implies that the functional is bounded, and its operator norm satisfies the inequality in one direction:
  
    \[
    \|f\| \le \|y\|_q
    \]

    * We construct the mapping $T: L^q[a,b] \longrightarrow (L^p[a,b])^*$ defined by $y \mapsto f$. Clearly, the operator $T$ is linear and, by the above, bounded. Thus $T$ is a bounded linear operator from $L^q$ to the dual of $L^p$.

    **Step 2: Prove that the mapping $T$ is injective by contradiction**

    * Suppose the mapping is not injective. Then there exists a non-zero element $y \neq \theta$ such that its image functional is the zero functional, i.e., $f = Ty = \theta$.

    * Since $y$ is non-zero in the space, there must exist some sufficiently small positive number $N > 0$ such that the set $e_N = \{t \in [a,b] \mid |y(t)| \ge N\}$ has Lebesgue measure strictly positive, i.e., $m(e_N) > 0$.

    * For this set, we construct a test function as follows:
  
    \[
    \chi_N(t) = \begin{cases} \text{sgn } y(t), & t \in e_N \\ 0, & t \notin e_N \end{cases}
    \]
  
    Clearly, $\chi_N(t)$, being a bounded function with support of finite measure, belongs to $L^p[a,b]$.

    * Substitute this test function into the integral expression of the functional $f$:
  
    \[
    f(\chi_N) = \int_a^b \chi_N(t)y(t) dt = \int_{e_N} \text{sgn}(y(t)) \cdot y(t) dt = \int_{e_N} |y(t)| dt \ge N \cdot m(e_N) > 0
    \]
  
    This shows that the functional $f$ is not the zero functional ($f \neq \theta$), which contradicts the initial assumption. Therefore, the mapping $T$ must be injective.

    **Step 3: Prove that the mapping $T$ is surjective (absolute continuity representation of functionals)**

    * Let $f$ be an arbitrary given bounded linear functional on $L^p[a,b]$. For any point $s \in [a,b]$, let $\chi_s$ denote the characteristic function of the subinterval $[a,s]$. Define an auxiliary function:
  
    \[
    g(s) = f(\chi_s)
    \]
  
    We will rigorously prove that $g$ is an **absolutely continuous function** on $[a,b]$.

    * Let $\delta_j = [s_j, t_j]$ ($j=1,2,\dots,n$) be a collection of closed subintervals of $[a,b]$ with disjoint interiors. Define the signs $\epsilon_j = \text{sgn}[g(t_j) - g(s_j)]$.

    * Using the linearity and boundedness of the functional $f$, we derive:
  
    \[
    \sum_{j=1}^n |g(t_j) - g(s_j)| = \sum_{j=1}^n \epsilon_j [g(t_j) - g(s_j)] = \sum_{j=1}^n \epsilon_j [f(\chi_{t_j}) - f(\chi_{s_j})]
    \]
  
    \[
    = f\left( \sum_{j=1}^n \epsilon_j [\chi_{t_j} - \chi_{s_j}] \right) \le \|f\| \cdot \left\| \sum_{j=1}^n \epsilon_j [\chi_{t_j} - \chi_{s_j}] \right\|_p
    \]
  
    \[
    = \|f\| \cdot \left( \int_a^b \left| \sum_{j=1}^n \epsilon_j [\chi_{t_j} - \chi_{s_j}] \right|^p dt \right)^{\frac{1}{p}}
    \]

    * Note that these intervals $\delta_j$ are disjoint. Hence the integral splits completely into a sum of integrals over each individual interval:
  
    \[
    = \|f\| \cdot \left( \sum_{j=1}^n \int_{\delta_j} |\epsilon_j|^p dt \right)^{\frac{1}{p}} = \|f\| \cdot \left( \sum_{j=1}^n m(\delta_j) \right)^{\frac{1}{p}}
    \]

    * From this final estimate, given any $\epsilon > 0$, we simply choose the control variable $\delta = \left( \frac{\epsilon}{\|f\|} \right)^p$. Then, whenever the total length of the intervals satisfies $\sum_{j=1}^n m(\delta_j) < \delta$, we strictly have:
  
    \[
    \sum_{j=1}^n |g(t_j) - g(s_j)| < \epsilon
    \]
  
    This perfectly matches the definition of an absolutely continuous function. Therefore, $g$ must be absolutely continuous.

    * By the classical calculus theorem for absolutely continuous functions, $g$ is differentiable almost everywhere, and its derivative is integrable. Define $y(s) = g'(s)$. Then $y \in L[a,b]$. Together with the boundary condition $g(a) = f(\chi_a) = f(\theta) = 0$, via the fundamental theorem of calculus we can write:
  
    \[
    f(\chi_s) = g(s) = \int_a^s y(t) dt = \int_a^b \chi_s(t)y(t) dt
    \]

    * This proves that the integral representation holds for all characteristic functions of intervals. Since any step function can be written as a finite linear combination of such characteristic functions, by the linearity of $f$, this integral representation also holds strictly for all step functions.

    * Next, let $x$ be any bounded measurable function defined on $[a,b]$. By density in real variable function theory, we can choose a uniformly bounded sequence of step functions $\{x_n\}$ that converges almost everywhere to $x$.

    * By the Bounded Convergence Theorem in measure theory, we have:
  
    \[
    \lim_{n \to \infty} \int_a^b x_n(t)y(t) dt = \int_a^b x(t)y(t) dt
    \]
  
    Moreover, this sequence also converges strongly in the $L^p$ norm:
  
    \[
    \|x_n - x\|_p = \left( \int_a^b |x_n(t) - x(t)|^p dt \right)^{\frac{1}{p}} \longrightarrow 0
    \]
  
    * Taking the limit $n \rightarrow \infty$ on both sides of the relation $f(x_n) = \int_a^b x_n(t)y(t) dt$, by the bounded continuity of the functional $f$, we successfully extend the conclusion to all bounded measurable functions:
  
    \[
    f(x) = \int_a^b x(t)y(t) dt
    \]

    **Step 4: Prove that the generated kernel function $y$ belongs to $L^q[a,b]$ and satisfies the norm estimate**

    * To remove the restriction of bounded measurability, we need a quantitative estimate of the norm of $y$. Introduce a sequence of bounded measurable functions on $[a,b]$ with a truncation property:
  
    \[
    y_N(t) = \begin{cases} |y(t)|^{q-1} \text{sgn } y(t), & |y(t)| \le N \\ 0, & |y(t)| > N \end{cases}
    \]
  
    Here $N$ is an arbitrary positive number. Denote the corresponding truncated set as $E_N = \{t \in [a,b] \mid |y(t)| \le N\}$.

    * Substitute this bounded measurable function $y_N$ as the argument into the functional expression established in Step 3:
  
    \[
    f(y_N) = \int_a^b y_N(t)y(t) dt = \int_{E_N} |y(t)|^{q-1} \text{sgn}(y(t)) \cdot y(t) dt = \int_{E_N} |y(t)|^q dt
    \]

    * On the other hand, using the boundedness of the functional $f$, we have $f(y_N) \le \|f\| \|y_N\|_p$, and thus:
  
    \[
    f(y_N) \le \|f\| \cdot \left( \int_a^b |y_N(t)|^p dt \right)^{\frac{1}{p}} = \|f\| \cdot \left( \int_{E_N} \left( |y(t)|^{q-1} \right)^p dt \right)^{\frac{1}{p}}
    \]
  
    * Note that the conjugate exponent relation gives $(q-1)p = q$. Substituting this into the power of the integral above, we obtain:
  
    \[
    f(y_N) \le \|f\| \cdot \left( \int_{E_N} |y(t)|^q dt \right)^{\frac{1}{p}}
    \]

    * Combining the two inequalities for $f(y_N)$ above and below, we get:
  
    \[
    \int_{E_N} |y(t)|^q dt \le \|f\| \cdot \left( \int_{E_N} |y(t)|^q dt \right)^{\frac{1}{p}}
    \]
  
    * Transposing the term on the right side (equivalently dividing both sides by the $1/p$ power of the integral), since $1 - \frac{1}{p} = \frac{1}{q}$, we successfully derive the key boundedness relation:
  
    \[
    \left( \int_{E_N} |y(t)|^q dt \right)^{\frac{1}{q}} \le \|f\|
    \]

    * Now let the truncation parameter go to infinity ($N \rightarrow \infty$). By the Monotone Convergence Theorem, the local integral on the left side transitions perfectly to the integral over the whole space, establishing:
  
    \[
    \|y\|_q = \left( \int_a^b |y(t)|^q dt \right)^{\frac{1}{q}} \le \|f\|
    \]
  
    This not only rigorously proves that the kernel function $y \in L^q[a,b]$, but also gives a lower bound estimate for the operator norm.

    **Step 5: Final extension to the whole space $L^p[a,b]$**

    * For any general function $x \in L^p[a,b]$, since the set of bounded measurable functions is dense in $L^p$, we can always find a sequence of bounded measurable functions $\{x_n\}$ such that in the $L^p$ norm sense:
  
    \[
    \int_a^b |x_n(t) - x(t)|^p dt \longrightarrow 0
    \]

    * For each bounded measurable function $x_n$, the relation $f(x_n) = \int_a^b x_n(t)y(t) dt$ is already strictly established. Now let $n \rightarrow \infty$ on both sides.
  
    * The left side converges to $f(x)$ by continuity of the functional; the right side also converges to $\int_a^b x(t)y(t) dt$ by Hölder's inequality, since $y \in L^q$ is already proven.

    * Hence, we have shown that for every element in the entire space $L^p[a,b]$, the following integral representation holds:
  
    \[
    f(x) = \int_a^b x(t)y(t) dt
    \]
  
    Moreover, combining the upper bound from Step 1 and the lower bound from Step 4 gives $\|f\| = \|y\|_q$. This shows that the mapping $T$ is indeed a surjective isometric isomorphism. $\square$