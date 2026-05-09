# 3.2 The Three Fundamental Convergence Theorems and Convergence in Measure

In Riemann integration, exchanging limits and integrals usually requires the extremely harsh condition of "uniform convergence." In the framework of Lebesgue integration, since integration is built on measure, we only need very weak conditions to achieve the exchange of integration and limit. The cornerstone of this grand theory is the famous "Three Fundamental Convergence Theorems."

---

## 1. Monotone Convergence Theorem (MCT)

The Monotone Convergence Theorem (MCT) is the first cornerstone in establishing integration theory. It states that for a sequence of non-negative increasing functions, the integral of the limit equals the limit of the integrals.

!!! success "Theorem 3.2.1 (Monotone Convergence Theorem MCT)"

    Let $\{f_n\}_{n=1}^\infty \subset L^+$ be a sequence of non-negative measurable functions. If it is monotonically increasing almost everywhere and converges pointwise to $f$:

    \[
    0 \le f_1 \le f_2 \le \dots \le f_n \le \dots \text{ a.e.} \quad \text{and} \quad \lim_{n \to \infty} f_n(x) = f(x) \text{ a.e.}
    \]

    Then $f \in L^+$, and we have:

    \[
    \int_X f \, d\mu = \lim_{n \to \infty} \int_X f_n \, d\mu
    \]

??? proof "Proof of Monotone Convergence Theorem"

    First, since $f_n \le f$, the monotonicity of the integral clearly gives $\int f_n \le \int f$, thus $\lim_{n \to \infty} \int f_n \le \int f$.

    Next, we need to prove the reverse inequality. Take any non-negative simple function $\phi$ satisfying $0 \le \phi \le f$, and any constant $\alpha \in (0, 1)$.
    Define a sequence of measurable sets:

    \[
    E_n = \{x \in X : f_n(x) \ge \alpha \phi(x)\}
    \]

    Since $f_n$ is monotonically increasing, the sequence of sets is increasing, i.e., $E_1 \subset E_2 \subset \dots$.
    Also, because $f_n \to f \ge \phi > \alpha \phi$ (when $\phi > 0$), we have $\bigcup_{n=1}^\infty E_n = X$.

    According to the properties of the integral, we have:

    \[
    \int_X f_n \, d\mu \ge \int_{E_n} f_n \, d\mu \ge \alpha \int_{E_n} \phi \, d\mu
    \]

    Letting $n \to \infty$, and using the continuity from below of the measure (and the integral of simple functions):

    \[
    \lim_{n \to \infty} \int_X f_n \, d\mu \ge \alpha \int_X \phi \, d\mu
    \]

    Since $\alpha \in (0, 1)$ is arbitrary, letting $\alpha \to 1$ gives $\lim \int f_n \ge \int \phi$. Taking the supremum over all simple functions satisfying $0 \le \phi \le f$, we get:

    \[
    \lim_{n \to \infty} \int_X f_n \, d\mu \ge \int_X f \, d\mu
    \]

    Combining both inequalities, the theorem is proved.

Using MCT, we can immediately derive the term-by-term summation property for the integral of a series.

!!! info "Corollary (Term-by-term Integration of Non-negative Series)"

    Let $\{f_n\}_{n=1}^\infty \subset L^+$ be a sequence of non-negative measurable functions, then:

    \[
    \int_X \left( \sum_{n=1}^\infty f_n \right) d\mu = \sum_{n=1}^\infty \int_X f_n \, d\mu
    \]

---

## 2. Fatou's Lemma

For a sequence of non-negative functions that is not necessarily monotonic, we cannot directly obtain an equality, but we can obtain a highly valuable inequality. It serves as a bridge connecting MCT and DCT.

!!! success "Theorem 3.2.2 (Fatou's Lemma)"

    Let $\{f_n\}_{n=1}^\infty \subset L^+$ be a sequence of non-negative measurable functions, then:

    \[
    \int_X \left( \liminf_{n \to \infty} f_n \right) d\mu \le \liminf_{n \to \infty} \int_X f_n \, d\mu
    \]

??? proof "Proof of Fatou's Lemma"

    Define a sequence of functions $g_n(x) = \inf_{k \ge n} f_k(x)$.

    Clearly $\{g_n\}$ is also a sequence of non-negative measurable functions, and since the range for the infimum shrinks as $n$ increases, $\{g_n\}$ is monotonically increasing, i.e., $g_1 \le g_2 \le \dots$.
    By the definition of limit inferior:

    \[
    \lim_{n \to \infty} g_n(x) = \liminf_{n \to \infty} f_n(x)
    \]

    Since $g_n \le f_k$ for all $k \ge n$, and specifically $g_n \le f_n$, we have:

    \[
    \int_X g_n \, d\mu \le \int_X f_n \, d\mu
    \]

    Taking the limit inferior on both sides:

    \[
    \liminf_{n \to \infty} \int_X g_n \, d\mu \le \liminf_{n \to \infty} \int_X f_n \, d\mu
    \]

    On the other hand, since $g_n$ is monotonically increasing, apply the Monotone Convergence Theorem (MCT) to $\{g_n\}$:

    \[
    \liminf_{n \to \infty} \int_X g_n \, d\mu = \lim_{n \to \infty} \int_X g_n \, d\mu = \int_X \left( \lim_{n \to \infty} g_n \right) d\mu = \int_X \left( \liminf_{n \to \infty} f_n \right) d\mu
    \]

    Combining the above results yields Fatou's Lemma.

---

## 3. Dominated Convergence Theorem (DCT)

This is the most frequently used convergence theorem in real analysis. It relaxes the requirement for monotonicity at the cost of requiring an integrable "dominating function."

!!! success "Theorem 3.2.3 (Lebesgue Dominated Convergence Theorem DCT)"

    Let $\{f_n\}$ be a sequence of measurable functions converging almost everywhere to $f$:

    \[
    \lim_{n \to \infty} f_n(x) = f(x) \quad \text{a.e.}
    \]

    If there exists a non-negative **integrable function** $g \in L^1(\mu)$ such that for all $n$, almost everywhere we have:

    \[
    |f_n(x)| \le g(x) \quad \text{a.e.}
    \]

    Then $f \in L^1(\mu)$, and:

    \[
    \lim_{n \to \infty} \int_X f_n \, d\mu = \int_X f \, d\mu \quad \text{and} \quad \lim_{n \to \infty} \int_X |f_n - f| \, d\mu = 0
    \]

??? proof "Proof of Dominated Convergence Theorem"

    Since $|f_n| \le g$ and $f_n \to f$, taking the limit gives $|f| \le g$ a.e. Since $g \in L^1$, we have $f \in L^1$.
    
    We construct two sequences of non-negative functions: $g + f_n \ge 0$ and $g - f_n \ge 0$.

    Apply Fatou's Lemma to $g + f_n$:

    \[
    \int_X (g + f) \, d\mu \le \liminf_{n \to \infty} \int_X (g + f_n) \, d\mu = \int_X g \, d\mu + \liminf_{n \to \infty} \int_X f_n \, d\mu
    \]

    Since $g$ is integrable, we can cancel $\int g$ on both sides to get:

    \[
    \int_X f \, d\mu \le \liminf_{n \to \infty} \int_X f_n \, d\mu
    \]

    Similarly, apply Fatou's Lemma to $g - f_n$:

    \[
    \int_X (g - f) \, d\mu \le \liminf_{n \to \infty} \int_X (g - f_n) \, d\mu = \int_X g \, d\mu - \limsup_{n \to \infty} \int_X f_n \, d\mu
    \]

    Again, cancelling $\int g$ and simplifying yields:

    \[
    \limsup_{n \to \infty} \int_X f_n \, d\mu \le \int_X f \, d\mu
    \]

    Combining these two inequalities, we must have $\liminf = \limsup = \lim$, thereby proving $\lim_{n \to \infty} \int f_n = \int f$.
    Replacing $f_n$ with $|f_n - f|$ (which is dominated by $2g$) in the above proof, we can similarly prove $\lim \int |f_n - f| = 0$.

---

## 4. Convergence in Measure and Classical Relationships

In addition to almost everywhere convergence (a.e.) and $L^1$ convergence, there is another type of convergence in measure theory that is extremely important in probability theory—convergence in measure (corresponding to convergence in probability).

!!! info "Definition 3.2.1 (Convergence in Measure)"

    Let $f_n, f$ be measurable functions. If for any $\epsilon > 0$, the measure of the set satisfying the following condition tends to 0:

    \[
    \lim_{n \to \infty} \mu(\{x \in X : |f_n(x) - f(x)| \ge \epsilon\}) = 0
    \]

    Then $f_n$ is said to **converge in measure** to $f$, denoted as $f_n \xrightarrow{\mu} f$.

Similarly, if $\mu(\{x : |f_n(x) - f_m(x)| \ge \epsilon\}) \to 0$ as $n, m \to \infty$, then $\{f_n\}$ is called a **Cauchy sequence in measure**.

There is a subtle relationship between convergence in measure and almost everywhere convergence. The famous Riesz Theorem points out the relationship between these two types of convergence via subsequences.

!!! success "Theorem 3.2.4 (Riesz Theorem)"

    * If $\{f_n\}$ is a Cauchy sequence in measure, then it must converge in measure to some measurable function $f$.

    * If $f_n \xrightarrow{\mu} f$, then there must exist a subsequence $\{f_{n_k}\}$ of $\{f_n\}$ that converges to $f$ almost everywhere:

    \[
    f_{n_k}(x) \to f(x) \quad \text{a.e. as } k \to \infty
    \]

On a finite measure space, almost everywhere convergence can even derive a property stronger than convergence in measure, known as **almost uniform convergence**.

!!! success "Theorem 3.2.5 (Egoroff's Theorem)"

    Let $\mu(X) < \infty$. If $\{f_n\}$ converges to $f$ almost everywhere on $X$, then for any $\epsilon > 0$, there exists a measurable subset $E \subset X$ such that:

    * $\mu(E) < \epsilon$ (The measure of the excluded "bad set" is minimal)

    * $f_n$ converges to $f$ **uniformly** on the remaining set $E^c = X \setminus E$.

??? proof "Core of the Proof for Egoroff's Theorem"

    Construct sets $E_n(k) = \bigcup_{m=n}^\infty \{x : |f_m(x) - f(x)| \ge \frac{1}{k}\}$.
    Since $f_n \to f$ a.e., for a fixed $k$, as $n \to \infty$, the sequence of sets $E_n(k)$ is monotonically decreasing and the measure of its intersection is 0.
    Since the measure of the entire space is finite, utilizing the continuity from above of the measure, we must have $\lim_{n \to \infty} \mu(E_n(k)) = 0$.
    
    For each $k$, we can select an $N_k$ large enough such that $\mu(E_{N_k}(k)) < \epsilon \cdot 2^{-k}$.
    Let the bad set be $E = \bigcup_{k=1}^\infty E_{N_k}(k)$. By subadditivity, $\mu(E) < \sum \epsilon 2^{-k} = \epsilon$.
    On $E^c$, for any $k$, as long as $n > N_k$, we have $|f_n - f| < 1/k$, which is exactly the definition of uniform convergence.