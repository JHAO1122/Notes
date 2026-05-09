# 4.2 Absolute Continuity and the Radon-Nikodym Theorem

In the previous section, we explored the internal structure of a single signed measure (the Hahn-Jordan decomposition). In this section, we shift our perspective to the **relationship between two different measures** on the same measurable space. If "mutually singular" characterizes the extreme case where two measures are completely unrelated, then "absolute continuity" describes the other extreme: one measure is entirely governed by another.

The pinnacle of exploring this relationship is a core pillar of modern analysis and probability theory—the **Lebesgue-Radon-Nikodym Theorem**.

---

## 1. Absolute Continuity of Measures

In calculus, absolutely continuous functions possess favorable properties. In measure theory, absolute continuity reflects the "heritability of null sets."

!!! info "Definition 4.2.1 (Absolute Continuity)"

    Let $\mu$ be a non-negative measure on a measurable space $(X, \mathcal{M})$, and $\nu$ be a signed measure (or complex measure) on the same space.
    
    If for any measurable set $E \in \mathcal{M}$ satisfying $\mu(E) = 0$, we have:

    \[
    \nu(E) = 0
    \]

    then $\nu$ is said to be **absolutely continuous** with respect to $\mu$, denoted as $\nu \ll \mu$.

* **Note**: Since the value of the total variation measure $|\nu|$ on $E$ is determined by the values of $\nu$ on subsets of $E$, $\nu \ll \mu$ if and only if $|\nu| \ll \mu$.

Absolute continuity has an extremely intuitive $\epsilon-\delta$ characterization, which explains the origin of the term "continuity": as long as the $\mu$-measure of a set is sufficiently small, its $\nu$-measure can be made arbitrarily small.

!!! success "Theorem 4.2.1 ($\epsilon-\delta$ Characterization of Absolute Continuity)"

    Let $\nu$ be a **finite** signed measure (i.e., $|\nu|(X) < \infty$). Then $\nu \ll \mu$ if and only if: for any $\epsilon > 0$, there exists $\delta > 0$ such that for any measurable set $E \in \mathcal{M}$ satisfying $\mu(E) < \delta$, we have:

    \[
    |\nu(E)| < \epsilon
    \]

??? proof "Proof of Theorem 4.2.1"

    **Sufficiency ($\Leftarrow$)**: If $\mu(E) = 0$, then for any $\delta > 0$, the condition $\mu(E) < \delta$ is satisfied. From the condition, it follows that for any $\epsilon > 0$, $|\nu(E)| < \epsilon$; therefore, we must have $\nu(E) = 0$, i.e., $\nu \ll \mu$.

    **Necessity ($\Rightarrow$)**: We use proof by contradiction. Suppose the conclusion does not hold; then there exists some $\epsilon > 0$ such that for every $n \in \mathbb{N}$ (taking $\delta = 1/2^n$), there exists a corresponding set $E_n \in \mathcal{M}$ satisfying:

    \[
    \mu(E_n) < \frac{1}{2^n} \quad \text{but} \quad |\nu(E_n)| \ge \epsilon
    \]

    Define the "limit superior" of the sequence of sets:

    \[
    F_k = \bigcup_{n=k}^\infty E_n, \quad F = \bigcap_{k=1}^\infty F_k
    \]

    Due to the subadditivity of the measure $\mu$:

    \[
    \mu(F_k) \le \sum_{n=k}^\infty \mu(E_n) < \sum_{n=k}^\infty \frac{1}{2^n} = \frac{1}{2^{k-1}}
    \]

    Clearly, as $k \to \infty$, $\mu(F_k) \to 0$. By continuity from above, we have $\mu(F) = 0$.
    
    Since it is known that $\nu \ll \mu$, we must have $|\nu|(F) = 0$.
    
    However, since $\nu$ is a finite measure, utilizing continuity from above:

    \[
    |\nu|(F) = \lim_{k \to \infty} |\nu|(F_k) \ge \limsup_{k \to \infty} |\nu(E_k)| \ge \epsilon > 0
    \]

    This leads to a severe contradiction! Thus, the assumption is false and the original proposition is proven.

---

## 2. Lebesgue-Radon-Nikodym Theorem

In many application scenarios, given an integral measure $d\nu = f d\mu$, we can easily determine that $\nu \ll \mu$. Conversely, if it is known that $\nu \ll \mu$, can we find a function $f$ such that $d\nu = f d\mu$?

The Lebesgue-Radon-Nikodym Theorem answers this question perfectly and provides a more general decomposition structure.

!!! success "Theorem 4.2.2 (Lebesgue-Radon-Nikodym Theorem)"

    Let $\mu$ and $\nu$ be two **$\sigma$-finite measures** on a measurable space $(X, \mathcal{M})$ ($\nu$ can be a signed measure).
    
    Then, there exists a **unique** pair of measures $\lambda$ and $\rho$ such that:

    \[
    \nu = \lambda + \rho
    \]

    and satisfying the following two core conditions:

    * **Singularity**: $\lambda \perp \mu$ (i.e., there exist complementary sets that separate them completely).

    * **Absolute Continuity**: $\rho \ll \mu$. Furthermore, there exists a $\mu$-almost everywhere unique function $f \in L^1(\mu)$ (if $\nu$ does not take negative infinity, this is an extended locally integrable function) such that for all measurable sets $E \in \mathcal{M}$:

    \[
    \rho(E) = \int_E f \, d\mu
    \]

* **Note**: This decomposition is called the **Lebesgue decomposition** of $\nu$ with respect to $\mu$. $\lambda$ is called the singular part, and $\rho$ is called the absolutely continuous part. The function $f$ serves exactly as the "density" between the two.

??? proof "Core Constructive Proof of the Lebesgue-Radon-Nikodym Theorem (for the finite non-negative case)"

    For simplicity, assume $\mu$ and $\nu$ are finite non-negative measures. We seek $f$ by constructing the largest "under-function."

    Define the following family of functions:

    \[
    \mathcal{F} = \left\{ g : X \to [0, +\infty] \text{ is a measurable function } \mid \int_E g \, d\mu \le \nu(E), \forall E \in \mathcal{M} \right\}
    \]

    Clearly, the constant function $g = 0$ is in $\mathcal{F}$, so $\mathcal{F}$ is non-empty. Moreover, if $g_1, g_2 \in \mathcal{F}$, it can be proven that $\max(g_1, g_2) \in \mathcal{F}$.

    Let $M = \sup_{g \in \mathcal{F}} \int_X g \, d\mu$. Since $\nu(X) < \infty$, we have $M < \infty$.
    We can obtain a sequence $\{g_n\} \subset \mathcal{F}$ such that $\int g_n d\mu \to M$.
    Let $f_n = \max(g_1, \dots, g_n)$, then $\{f_n\}$ is a monotonically increasing sequence of functions, and $f_n \in \mathcal{F}$.
    
    According to the Monotone Convergence Theorem (MCT), let $f = \lim_{n \to \infty} f_n = \sup_n f_n$, then $f \in \mathcal{F}$ and:

    \[
    \int_X f \, d\mu = M
    \]

    Now, we define the remaining measure part $\lambda$ as:

    \[
    \lambda(E) = \nu(E) - \int_E f \, d\mu
    \]

    From $f \in \mathcal{F}$, it follows that $\lambda$ is a non-negative measure. We now only need to prove $\lambda \perp \mu$.
    If $\lambda$ were not singular to $\mu$, then according to measure theory techniques, we could necessarily find a tiny perturbation $\epsilon \mu \le \lambda$ on some set of positive measure, such that $f + \epsilon \chi_A$ would still belong to $\mathcal{F}$.
    However, this would result in $\int (f + \epsilon \chi_A) d\mu = M + \epsilon \mu(A) > M$, contradicting the fact that $M$ is the supremum!
    
    Therefore, $\lambda$ must be singular with respect to $\mu$, proving the non-negative case of the theorem. Signed measures or $\sigma$-finite measures can be reduced to this case through the Hahn decomposition and space partitioning.

---

## 3. Radon-Nikodym Theorem and Derivatives

As a direct corollary of the Lebesgue decomposition theorem, when it is known in advance that $\nu \ll \mu$, the singular part $\lambda$ collapses to 0.

!!! success "Theorem 4.2.3 (Radon-Nikodym Theorem)"

    Let $\mu, \nu$ be $\sigma$-finite measures, and $\nu \ll \mu$.
    
    Then there must exist a $\mu$-almost everywhere unique function $f$ such that for any measurable set $E \in \mathcal{M}$:

    \[
    \nu(E) = \int_E f \, d\mu
    \]

This miraculous function $f$ is called the **Radon-Nikodym derivative** of $\nu$ with respect to $\mu$, usually denoted by the fractional symbol:

\[
f = \frac{d\nu}{d\mu}
\]

The Radon-Nikodym derivative perfectly simulates the derivative chain rule from calculus within integral operations.

!!! success "Proposition 4.2.1 (Properties of the Radon-Nikodym Derivative)"

    Assume all measures are $\sigma$-finite measures.

    * **Integration Rule**: If $g \in L^1(\nu)$, then $g \frac{d\nu}{d\mu} \in L^1(\mu)$, and there is an intuitive change-of-variable formula:

    \[
    \int_X g \, d\nu = \int_X g \frac{d\nu}{d\mu} \, d\mu
    \]

    * **Chain Rule**: If $\nu \ll \mu$ and $\mu \ll \lambda$, then $\nu \ll \lambda$, and almost everywhere we have:

    \[
    \frac{d\nu}{d\lambda} = \frac{d\nu}{d\mu} \frac{d\mu}{d\lambda}
    \]

    * **Reciprocal Rule**: If $\mu \ll \nu$ and $\nu \ll \mu$ (i.e., they are **equivalent**, possessing identical null sets), then almost everywhere we have:

    \[
    \frac{d\nu}{d\mu} = \left( \frac{d\mu}{d\nu} \right)^{-1}
    \]

---

## 4. Complex Measures and Polar Decomposition

The values of a signed measure are restricted to real numbers. If we allow measures to take complex values, we obtain **complex measures**, which are very important in functional analysis.

!!! info "Definition 4.2.2 (Complex Measure)"

    Let $(X, \mathcal{M})$ be a measurable space. A set function $\nu: \mathcal{M} \to \mathbb{C}$ is called a **complex measure** if for any sequence of disjoint measurable sets $\{E_j\}$ it satisfies:

    \[
    \nu\left(\bigcup_{j=1}^\infty E_j\right) = \sum_{j=1}^\infty \nu(E_j)
    \]

    Note: To ensure the sum of a complex series is independent of the order of terms, the series must be **absolutely convergent**. This means complex measures **naturally take only finite values** and cannot take $\infty$.

Any complex measure $\nu$ can be uniquely decomposed into real and imaginary parts: $\nu = \nu_r + i \nu_i$, where $\nu_r$ and $\nu_i$ are both finite signed measures.

!!! info "Definition 4.2.3 (Total Variation)"

    For a complex measure $\nu$, define its total variation measure $|\nu|$ as:

    \[
    |\nu|(E) = \sup \sum_{j=1}^n |\nu(E_j)|
    \]

    where the supremum is taken over all finite measurable partitions $\{E_j\}_{j=1}^n$ of the set $E$. $|\nu|$ is also a finite non-negative measure.

With the support of the Radon-Nikodym theorem, we can perform an exquisite "polar decomposition" of a complex measure, much like representing a complex number in polar form $z = r e^{i\theta}$.

!!! success "Theorem 4.2.4 (Polar Decomposition of Complex Measures)"

    Let $\nu$ be a complex measure. Since from the definition of total variation it is clear that $\nu \ll |\nu|$, by the Radon-Nikodym theorem, there must exist a complex measurable function $h = \frac{d\nu}{d|\nu|}$ such that:

    \[
    d\nu = h \, d|\nu|
    \]

    More importantly, this function $h$ satisfies, $|\nu|$-almost everywhere:

    \[
    |h(x)| = 1
    \]

Using polar decomposition, integration with respect to a complex measure can be strictly defined as a conversion to integration with respect to its total variation measure (a non-negative measure). If $f \in L^1(|\nu|)$, we define:

\[
\int_X f \, d\nu = \int_X f \cdot h \, d|\nu|
\]

From this definition, we can very directly derive the inequality bound for the absolute value of a complex integral:

\[
\left| \int_X f \, d\nu \right| \le \int_X |f| \, d|\nu|
\]