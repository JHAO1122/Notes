# Chapter 6: Banach Spaces and Quotient Spaces

In Chapter 5, we introduced the concept of normed linear spaces. In order to thoroughly discuss limit operations and calculus within these spaces, we require the spaces to possess a "no-hole" property. This chapter will focus on the completeness of normed linear spaces, namely Banach spaces, and introduce the topological constructions of quotient spaces and direct sum spaces.

---

## 1. Banach Space

!!! info "Definition 6.1 (Banach Space)"

    A complete normed linear space is called a **Banach space**.
    That is, in a normed linear space, if every Cauchy sequence (fundamental sequence) converges in norm to a point within that space, then the space is said to be complete.

As proven in Chapter 1, finite-dimensional spaces $\mathbb{R}^n$, $\mathbb{C}^n$, and the continuous function space $C[a, b]$ are all complete, hence they are all Banach spaces.

---

## 2. Completeness Proofs for Typical Spaces

### 2.1 Completeness of the Continuously Differentiable Space $C^k[a, b]$

$C^k[a, b]$ is the set consisting of all functions that possess continuous derivatives up to order $k$ on the interval $[a, b]$. Its norm is defined as $\|x\| = \sum_{j=0}^{k} \max_{t \in [a, b]} |x^{(j)}(t)|$.

??? proof "Proof: $C^k[a, b]$ is a Banach Space"

    * Let $\{x_n\} \subset C^k[a, b]$ be a Cauchy sequence.

    * Then for any given $\epsilon > 0$, there exists a positive integer $N > 0$ such that when $m, n > N$, $||x_n - x_m|| < \epsilon$. That is:

    $$
    \sum_{j=0}^{k} \max_{t \in [a, b]} |x_n^{(j)}(t) - x_m^{(j)}(t)| < \epsilon
    $$

    * Therefore, for each $(0 \le j \le k)$, the inequality:

    $$
    |x_n^{(j)}(t) - x_m^{(j)}(t)| < \epsilon
    $$

    * holds uniformly with respect to $t \in [a, b]$.

    * From mathematical analysis, it is known that $\{x_n^{(j)}\}$ converges uniformly to a continuous function $y_j(t)$, and $y_{j+1}(t)$ is also the derivative of $y_j(t)$ $(0 \le j \le k-1)$.

    * Therefore, $y_0(t)$ has continuous derivatives up to order $k$, and $y_j(t) = y_0^{(j)}(t)$.

    * When $n > N$:

    $$
    \sum_{j=0}^{k} \max_{t \in [a, b]} |x_n^{(j)}(t) - y_0^{(j)}(t)| \le \epsilon
    $$

    * Thus, $\{x_n\}$ converges to $y_0(t)$ under the norm of $C^k[a, b]$, hence $C^k[a, b]$ is complete, and is therefore a Banach space. $\square$

### 2.2 Completeness of the Integrable Function Space $L^p(F)$

??? proof "Proof: $L^p(F)$ is a Banach Space"

    * Let $\{f_n\} \subset L^p(F)$ be a Cauchy sequence. Then for any given $\epsilon > 0$, there exists a positive integer $m_k > 0$ such that when $m, n > m_k$:

    $$
    ||f_n - f_m||_p < \frac{1}{2^k}, \quad k = 1, 2, \dots
    $$

    * Choose $n_k \ge m_k$ such that $n_1 < n_2 < \dots < n_k < \dots$, then:

    $$
    ||f_{n_k} - f_{n_{k+1}}||_p < \frac{1}{2^k}, \quad k = 1, 2, \dots
    $$

    * Therefore:

    $$
    \sum_{k=1}^{\infty} ||f_{n_{k+1}} - f_{n_k}||_p \le \sum_{k=1}^{\infty} \frac{1}{2^k} < \infty
    $$

    * However, $1 \in L^q(F)$. By Hölder's inequality, we obtain:

    $$
    \int_F |f_{n_{k+1}}(t) - f_{n_k}(t)| dt \le ||f_{n_{k+1}} - f_{n_k}||_p m(F)^{\frac{1}{q}}
    $$

    * Thus the series:

    $$
    \sum_{k=1}^{\infty} \int_F |f_{n_{k+1}}(t) - f_{n_k}(t)| dt \quad \text{converges.}
    $$

    * By Levi's Theorem, the series $\sum_{k=1}^{\infty} |f_{n_{k+1}}(t) - f_{n_k}(t)|$ converges almost everywhere on $F$.

    * Furthermore, it follows that:

    $$
    f_{n_k}(t) = f_{n_1}(t) + \sum_{j=1}^{k-1} (f_{n_{j+1}}(t) - f_{n_j}(t))
    $$

    * converges almost everywhere to a measurable function $f(t)$. Next, we prove that $f \in L^p(F)$ and that the norm converges:

    * Since $\{f_n\}$ is a Cauchy sequence in $L^p(F)$, for any given $\epsilon > 0$, there exists a positive integer $N > 0$ such that when $m, n > N$, $||f_n - f_m||_p < \epsilon$.

    * Choose a sufficiently large $k_0$ such that $n_{k_0} > N$. Then when $k \ge k_0$ and $n \ge N$, we have:

    $$
    \int_F |f_n(t) - f_{n_k}(t)|^p dt = ||f_n - f_{n_k}||_p^p < \epsilon^p
    $$

    * As $k \rightarrow \infty$, the sequence of functions on $F$ satisfies almost everywhere (a.e.):

    $$
    |f_n(t) - f_{n_k}(t)|^p \rightarrow |f_n(t) - f(t)|^p
    $$

    * By Fatou's Lemma, we know that $|f_n(t) - f(t)|^p$ is integrable, and:

    $$
    \int_F |f_n(t) - f(t)|^p dt \le \liminf_{k \rightarrow \infty} \int_F |f_n(t) - f_{n_k}(t)|^p dt \le \epsilon^p
    $$

    * This shows that $f - f_n \in L^p(F)$, and when $n \ge N$:

    $$
    ||f_n - f||_p \le \epsilon
    $$

    * Since $L^p(F)$ is a linear space, and $f_n \in L^p(F)$, then $f = (f - f_n) + f_n \in L^p(F)$.

    * In conclusion, $f_n \rightarrow f$ holds in $L^p(F)$, and $L^p(F)$ is a Banach space. $\square$

---

## 3. Existence of Incomplete Normed Linear Spaces

Not all normed linear spaces are complete. Completeness depends not only on the elements of the space itself but also on the norm assigned to it.

* Consider the space $C[a, b] \subset L^2[a, b]$, and choose $a < c < b$.

* Define the step function $f_0 = \chi_{[a, c]} - \chi_{(c, b]}$.

* Since the continuous function space $C[a, b]$ is dense in $L^2[a, b]$, there must exist a sequence of continuous functions $\{f_n\}$ in $C[a, b]$ such that:

    $$
    ||f_n - f_0||_2 \rightarrow 0
    $$

* However, the limit function $f_0$ is discontinuous at point $c$, thus $f_0 \notin C[a, b]$. Therefore, $C[a, b]$ is incomplete under the $L^2$ norm $|| \cdot ||_2$.

---

## 4. Quotient Space

### 4.1 Definition of Equivalence Classes and Quotient Spaces

* Let $E$ be a given linear space, and $L$ be a subspace of $E$. We can define an equivalence relation on $E$: $x \sim y$, if $x - y \in L$.

* This is an equivalence relation satisfying the following three properties:

* (i) Reflexivity: $x \sim x$;

* (ii) Symmetry: if $x \sim y$, then $y \sim x$;

* (iii) Transitivity: if $x \sim y$ and $y \sim z$, then $x \sim z$.

* By properties (i)-(iii), equivalence classes can be defined in $E$: for any two elements $x, y$ in $E$, if $x - y \in L$, then $x$ and $y$ are grouped into the same equivalence class. These are denoted by symbols such as $\xi, \eta$.

* It is easy to prove that any element in $E$ must belong to exactly one equivalence class, and any two distinct equivalence classes share no common elements. Thus, space $E$ is partitioned into the union of several equivalence classes.

* For a given equivalence class $\xi$, pick any $x \in \xi$, then the set:

    $$
    x + L := \{x + y \mid y \in L\}
    $$

    is exactly $\xi$. Therefore, equivalence classes can also be denoted as $x + L$, $y + L$, etc.

### 4.2 Linear Operations in Quotient Spaces

* Let $\hat{E}$ (or $E/L$) denote the set composed of all equivalence classes in $E$. Now we define linear operations in $\hat{E}$: let $\xi, \eta \in \hat{E}$, pick any $x \in \xi$ and $y \in \eta$, we stipulate:

    $$
    \xi + \eta := x + y + L
    $$

    $$
    \alpha \xi := \alpha x + L
    $$

* We need to prove that these operations are well-defined: that is, $x + y + L$ is independent of the choice of $x \in \xi, y \in \eta$, and the scalar multiplication $\alpha x + L$ is independent of the choice of $x \in \xi$.

* Proof: Pick another $x' \in \xi, y' \in \eta$, we have:

    $$
    x' + y' + L = x + y + (x' - x) + (y' - y) + L = x + y + L
    $$

* After defining linear operations in $\hat{E}$, it can be proven that $\hat{E}$ is a linear space under such linear operations, which is called the quotient space of $E$ with respect to $L$, and it is easily seen that its zero element is the subspace $L$.

### 4.3 Example: $L^p(F)$ is a Quotient Space

* In $L^p(F)$, we stipulate that any two functions that are equal almost everywhere represent the same element.

* Define the subset $\tilde{L} := \{g \mid g \text{ is equal to } 0 \text{ almost everywhere on } F\}$.

* Then the elements in $L^p(F)$ are actually equivalence classes of the form $f + \tilde{L}$, where $f$ satisfies $\int_F |f(t)|^p dt < \infty$.

* Therefore, $L^p(F)$ is actually a quotient space of the linear space composed of all $p$-integrable functions on $F$ with respect to the subspace $\tilde{L}$.

---

## 5. Norming of Quotient Spaces and Direct Sum Spaces

### 5.1 $E/L$ Normed Linear Space

* Let $E$ be a normed linear space, and $L$ be a subspace of $E$. $L$ is clearly a normed linear space under the norm of $E$. If $L$ is closed, then $L$ is called a **closed subspace** of $E$.

* When $L$ is a closed subspace of $E$, a norm can be introduced in the quotient space $E/L$. For any $\xi \in E/L$, let:

    $$
    ||\xi|| = \inf_{x \in \xi} ||x||
    $$

??? proof "Proof: This satisfies all conditions for a norm"

    * If $||\xi|| = 0$, let $\xi = x_0 + L$, where $x_0 \in \xi$. Then:

    $$
    0 = \inf_{x \in \xi} ||x|| = \inf_{y \in L} ||x_0 + y||
    $$

    * For any positive integer $n$, there exists $y_n \in L$ such that:

    $$
    ||x_0 + y_n|| < \frac{1}{n}
    $$

    * Then $x_0 + y_n \rightarrow 0$, which means $y_n \rightarrow -x_0$. Because $L$ is a closed set, the limit point $-x_0 \in L$. This implies that $x_0 \in L$, that is, $\xi = L$ (the zero element of the quotient space).

    * Proof of the triangle inequality:

    $$
    ||\xi + \eta|| = \inf_{x \in \xi, y \in \eta} ||x + y|| \le ||x|| + ||y|| \le ||\xi|| + \epsilon + ||\eta|| + \epsilon
    $$

    * Due to the arbitrariness of $\epsilon$, $E/L$ is a normed linear space under this norm. $\square$

* It can be further proven that: when $E$ is a Banach space, and $L$ is a closed subspace of $E$, then $E/L$ is also a Banach space.

### 5.2 Norms on Direct Sum Spaces

* Let $L_1, \dots, L_n$ be normed linear spaces, and $E$ be the direct sum of $L_1, \dots, L_n$, that is, $E = L_1 \oplus \dots \oplus L_n$. A norm in $E$ is defined as follows:

    $$
    ||x|| = ||x_1|| + \dots + ||x_n||
    $$

* Where $x = (x_1, \dots, x_n) \in E$, and $x_k \in L_k$ $(k=1, 2, \dots, n)$. $E$ is a normed linear space under this norm.

* If $||x|| = ||x_1|| + \dots + ||x_n|| = 0$, then it must be that $||x_1|| = \dots = ||x_n|| = 0$, meaning each $x_k$ is the zero element of $L_k$.

* Other equivalent norms can also be defined in $E$:

    $$
    ||x|| = \max \{||x_1||, \dots, ||x_n||\}
    $$

* Or the 2-norm:

    $$
    ||x|| = \sqrt{||x_1||^2 + \dots + ||x_n||^2}
    $$

* They are all norms on $E$, thus $E$ is also a normed linear space under these norms.

* If the subspaces $L_1, \dots, L_n$ are all Banach spaces, $E$ must also be a Banach space under these equivalent norms.

---

## 6. Chapter Summary

* **Linear Space**: The linear operations in a linear space belong to the realm of algebra; in a linear space, there is neither distance nor norm, thus there is no convergence, nor open sets, closed sets, denseness, separability, etc.

* **Normed Linear Space**: A norm can be introduced in certain linear spaces to make them normed linear spaces. When we introduce distance into a linear space, we should require the linear operations to be continuous with respect to the introduced distance.

* **Topological Linear Space**: A general topology can be introduced in certain linear spaces to obtain a topological linear space. It should be required that the linear operations are continuous with respect to the introduced topology.