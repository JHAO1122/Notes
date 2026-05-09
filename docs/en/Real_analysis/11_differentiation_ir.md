# 4.3 Differentiation on the Real Line and the Fundamental Theorem of Calculus

In classical mathematical analysis, the Fundamental Theorem of Calculus (Newton–Leibniz formula) reveals the inverse relationship between differentiation and integration. Within the framework of Lebesgue integration and measure theory, not only can integration be extended to abstract measure spaces, but the concept of "derivative" can also be generalized to the derivative of a measure with respect to another measure (i.e., the Radon–Nikodym derivative). In this section, we return to the real line $\mathbb{R}$ and $\mathbb{R}^n$, exploring the differential properties of locally integrable functions, functions of bounded variation, absolutely continuous functions, and finally present the Fundamental Theorem of Calculus in the Lebesgue sense.

---

## 1. Differentiation of Measures and the Hardy–Littlewood Maximal Function

In $\mathbb{R}^n$, we attempt to define the "derivative" of a measure by examining the average density of the measure in small balls near a point. Let $\nu$ be a complex measure or signed measure on $\mathbb{R}^n$ and let $m$ be Lebesgue measure. Define the derivative at a point $x$ as

\[
F(x) := \lim_{r \to 0} \frac{\nu(B(r, x))}{m(B(r, x))}
\]

If $d\nu = f\,dm$, i.e., $\nu$ is absolutely continuous with respect to $m$ with density $f$, then the above limit can be written as the average of $f$ over the ball $B(r, x)$:

\[
\frac{\nu(B(r, x))}{m(B(r, x))} = \frac{\int_{B(r, x)} f(y) \, dm(y)}{\int_{B(r, x)} 1 \, dm(y)}
\]

We naturally hope that as $r \to 0$, this average converges to $f(x)$ almost everywhere. To prove this, we need to introduce a powerful tool in analysis—the Hardy–Littlewood maximal function.

!!! info "Definition 4.3.1 (Locally Integrable and Averaging Operator)"

    * **Locally integrable space $L_{\text{loc}}^1(\mathbb{R}^n)$**: A function $f: \mathbb{R}^n \to \mathbb{C}$ is said to belong to $L_{\text{loc}}^1(\mathbb{R}^n)$ if it is Lebesgue integrable on every compact (or bounded) set, i.e., for every measurable set $K$ of finite measure we have $\int_K |f|\,dm < \infty$.

    * **Averaging operator $A_r f(x)$**: For $f \in L_{\text{loc}}^1(\mathbb{R}^n)$, define its average over the open ball $B(r, x)$ of radius $r$ centered at $x$ as

    \[
    A_r f(x) = \frac{1}{m(B(r, x))} \int_{B(r, x)} f(y) \, dy
    \]

!!! info "Definition 4.3.2 (Hardy–Littlewood Maximal Function)"

    For $f \in L_{\text{loc}}^1(\mathbb{R}^n)$, define the maximal function as the supremum of absolute averages over all radii $r > 0$:

    \[
    Hf(x) = \sup_{r > 0} A_r |f|(x) = \sup_{r > 0} \frac{1}{m(B(r, x))} \int_{B(r, x)} |f(y)| \, dy
    \]

The maximal function $Hf(x)$ describes the most extreme local oscillation of $f$ near $x$. The most important property of the maximal function is that it satisfies a **weak type (1,1) inequality**, i.e., the tail probability is controlled by the $L^1$ norm of the original function.

!!! success "Theorem 4.3.1 (Maximal Theorem)"

    If $f \in L^1(\mathbb{R}^n)$, then for any $\alpha > 0$, the Lebesgue measure of the level set $E_\alpha = \{x : Hf(x) > \alpha\}$ satisfies

    \[
    m(\{x : Hf(x) > \alpha\}) \le \frac{3^n}{\alpha} \int_{\mathbb{R}^n} |f| \, dm
    \]

??? proof "Proof Outline of the Maximal Theorem (Covering Lemma)"

    Take any $x \in E_\alpha$. By the definition of supremum, there exists a radius $r_x > 0$ such that the average over that ball exceeds $\alpha$:

    \[
    \frac{1}{m(B(r_x, x))} \int_{B(r_x, x)} |f| \, dm > \alpha \quad \Longrightarrow \quad m(B(r_x, x)) < \frac{1}{\alpha} \int_{B(r_x, x)} |f| \, dm
    \]

    All these balls $\{B(r_x, x)\}_{x \in E_\alpha}$ form an open cover of $E_\alpha$.
  
    Using the **Wiener/Vitali covering lemma**, we can select a **disjoint** finite or countable subcollection $\{B_j\}_{j=1}^\infty$ such that after expanding the radii of these selected balls by a factor of 3, the expanded balls cover the core part of the original cover, even every bounded subset. By the geometric property of the covering lemma, the measure magnification factor is $3^n$:

    \[
    m(E_\alpha) \le m\left(\bigcup_{j} 3 B_j\right) \le 3^n \sum_j m(B_j)
    \]

    Since the balls $B_j$ are disjoint, we can estimate the sum:

    \[
    m(E_\alpha) \le 3^n \sum_j \frac{1}{\alpha} \int_{B_j} |f| \, dm \le \frac{3^n}{\alpha} \int_{\bigcup B_j} |f| \, dm \le \frac{3^n}{\alpha} \int_{\mathbb{R}^n} |f| \, dm
    \]

    The theorem is proved.

Using the maximal theorem, we can prove one of the most beautiful theorems in real analysis—the Lebesgue differentiation theorem.

!!! success "Theorem 4.3.2 (Lebesgue Differentiation Theorem)"

    If $f \in L_{\text{loc}}^1(\mathbb{R}^n)$, then for almost every $x \in \mathbb{R}^n$ (i.e., except on a Lebesgue null set), we have

    \[
    \lim_{r \to 0} A_r f(x) = f(x)
    \]

    More strongly, almost every $x$ is actually a **Lebesgue point** of $f$, satisfying

    \[
    \lim_{r \to 0} \frac{1}{m(B(r, x))} \int_{B(r, x)} |f(y) - f(x)| \, dy = 0
    \]

---

## 2. Functions of Bounded Variation (BV) and Borel Measures

On the one-dimensional real line $\mathbb{R}$, signed (or complex) measures can be described by cumulative distribution functions. However, while monotone functions correspond to non‑negative measures, signed measures correspond to **functions of bounded variation**.

!!! info "Definition 4.3.3 (Function of Bounded Variation, BV)"

    Let $F: \mathbb{R} \to \mathbb{C}$. If for any finite strictly increasing partition points $x_0 < x_1 < \dots < x_n$ in $\mathbb{R}$, the sum of absolute differences has a uniform upper bound, i.e., the total variation

    \[
    T_F(x) = \sup \sum_{j=1}^n |F(x_j) - F(x_{j-1})| < \infty
    \]

    then $F$ is called a **function of bounded variation**, denoted $F \in BV$.

Any function of bounded variation can be decomposed as the difference of two monotone increasing functions (Jordan decomposition). To establish a one‑to‑one correspondence with Borel measures, we need to normalize BV functions, analogous to requiring measures to be right‑continuous and vanishing at $-\infty$.

!!! info "Definition 4.3.4 (Normalized Bounded Variation Space NBV)"

    The space of **normalized bounded variation functions (NBV)** consists of all BV functions $F$ satisfying the following two conditions:

    * $F$ is right‑continuous, i.e., $F(x^+) = F(x)$.

    * $F$ vanishes at $-\infty$, i.e., $F(-\infty) = \lim_{x \to -\infty} F(x) = 0$.

!!! success "Theorem 4.3.3 (Correspondence between NBV and Borel Measures)"

    * If $\nu$ is a finite complex Borel measure on $\mathbb{R}$, define $F(x) = \nu((-\infty, x])$. Then $F \in NBV$.

    * Conversely, if $F \in NBV$, then there exists a **unique** finite complex Borel measure $\mu_F$ such that for any interval, $\mu_F((a, b]) = F(b) - F(a)$.

    * In this case, the total variation measure corresponds perfectly to the total variation function: $|\mu_F| = \mu_{T_F}$.

Since any function of bounded variation can be decomposed into monotone functions, and monotone functions are differentiable almost everywhere by Lebesgue's result, we have the following corollary:

!!! success "Proposition 4.3.1 (Differentiability of NBV Functions)"

    * If $F \in NBV$, then $F$ is differentiable almost everywhere and its derivative $F' \in L^1(m)$.

    * If $\mu_F \perp m$ (i.e., $\mu_F$ is singular with respect to Lebesgue measure), then $F'(x) = 0$ a.e.

---

## 3. Absolutely Continuous Functions (AC) and the Fundamental Theorem of Calculus

Although NBV functions are differentiable almost everywhere, the Newton–Leibniz formula does not hold for all of them. To make the Fundamental Theorem of Calculus strictly valid, we need a class of functions that is stronger than continuous functions, where "small local fluctuations cannot accumulate into a large jump" — **absolutely continuous functions**.

!!! info "Definition 4.3.5 (Absolutely Continuous Function, AC)"

    A function $F: [a, b] \to \mathbb{C}$ is said to be **absolutely continuous** on $[a, b]$ if for any given $\epsilon > 0$, there exists $\delta > 0$ such that for any finite collection of **pairwise disjoint open intervals** $\{(a_j, b_j)\}_{j=1}^n$ in $[a, b]$, if their total length satisfies

    \[
    \sum_{j=1}^n (b_j - a_j) < \delta
    \]

    then the corresponding total fluctuation of the function values satisfies

    \[
    \sum_{j=1}^n |F(b_j) - F(a_j)| < \epsilon
    \]

Clearly, absolute continuity (AC) implies uniform continuity, hence continuity. Moreover, at the measure level it is equivalent to absolute continuity of the associated measure:

!!! success "Proposition 4.3.2"

    If $F \in NBV$, then $F \in AC$ if and only if the measure induced by $F$ is absolutely continuous with respect to Lebesgue measure, i.e., $\mu_F \ll m$.

Based on this relationship, we can finally state the Fundamental Theorem of Calculus in real analysis:

!!! success "Theorem 4.3.4 (The Fundamental Theorem of Calculus)"

    Let $F: [a, b] \to \mathbb{C}$ be a function. Then the following statements are equivalent:

    * **1. AC property**: $F$ is absolutely continuous on $[a, b]$ ($F \in AC$).

    * **2. Integral representation**: There exists a locally integrable function $f \in L^1([a, b])$ such that for all $x \in [a, b]$,

    \[
    F(x) - F(a) = \int_a^x f(t) \, dt
    \]

    * **3. Almost everywhere differentiable and derivative integrable recovers the function**: $F$ is differentiable almost everywhere, $F' \in L^1([a, b])$, and the Newton–Leibniz formula holds exactly:

    \[
    F(x) - F(a) = \int_a^x F'(t) \, dt
    \]

---

## 4. Lebesgue Decomposition of Measures and Singular Measures

Combining the Radon–Nikodym theorem and the Lebesgue decomposition theorem, any finite Borel measure $\nu$ on $\mathbb{R}$ can be completely decomposed into three mutually singular pure parts:

\[
\nu = \nu_{ac} + \nu_{cs} + \nu_{d}
\]

* **1. Discrete measure $\nu_d$**:
  It is concentrated on a finite or countable set of points, of the form $\nu_d = \sum c_j \delta_{x_j}$, where $\delta_{x_j}$ is the Dirac point mass. Its corresponding distribution function is a pure "step function", jumping only at countably many points, and its derivative is zero almost everywhere.

* **2. Absolutely continuous measure $\nu_{ac}$**:
  It is dominated by Lebesgue measure $m$, i.e., $\nu_{ac} \ll m$. There exists a density function $f \in L^1$ such that $d\nu_{ac} = f\,dm$. Its corresponding distribution function is absolutely continuous (AC).

* **3. Continuous singular measure $\nu_{cs}$**:
  This is the most counterintuitive existence in measure theory. It is a singular measure ($\nu_{cs} \perp m$), meaning it is entirely concentrated on a set of Lebesgue measure zero; yet it is also a continuous measure, i.e., it has no point masses (no jumps) at any single point.

  **Classic example**: The Cantor measure built on the Cantor ternary set; the corresponding Cantor function (Devil's staircase) is continuous, monotone increasing, and has derivative zero almost everywhere (it is flat outside the Cantor set, which is a null set), yet it rises from 0 to 1. This demonstrates that not every continuous and almost everywhere differentiable function can be recovered by integrating its derivative!