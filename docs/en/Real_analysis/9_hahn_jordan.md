# 4.1 Signed Measures and Hahn-Jordan Decomposition

In previous chapters, measures were always defined as non-negative set functions (representing length, area, volume, or probability). However, in many physical and mathematical scenarios (such as charge distributions which can be positive or negative, or the difference between two integrals), we need to allow "measures" to take negative values. This leads to the concept of **Signed Measures**, also known as **variable-sign measures**.

---

## 1. Basic Concepts of Signed Measures

!!! info "Definition 4.1.1 (Signed Measure)"

    Let $(X, \mathcal{M})$ be a measurable space. A set function $\nu: \mathcal{M} \to [-\infty, +\infty]$ is called a **signed measure** if it satisfies the following three conditions:

    * **Zero for the empty set**:
    
    \[
    \nu(\emptyset) = 0
    \]

    * **Takes at most one infinity**: $\nu$ can take at most one of the values $+\infty$ or $-\infty$ (to avoid meaningless operations like $\infty - \infty$).

    * **Countable additivity**: For any sequence of disjoint sets $\{E_j\}_{j=1}^\infty$ in $\mathcal{M}$, we have:
    
    \[
    \nu\left(\bigcup_{j=1}^\infty E_j\right) = \sum_{j=1}^\infty \nu(E_j)
    \]
    
    Furthermore, if the measure of the union $\nu\left(\bigcup_{j=1}^\infty E_j\right)$ is finite, then the above series must be **absolutely convergent**.

* **Note 1**: Ordinary measures (non-negative measures) are naturally a special case of signed measures.
* **Note 2**: The most typical example of a signed measure is one induced by an integral. Let $\mu$ be a general measure and $f$ be a variable-sign $\mu$-integrable function (or at least the integral of its positive or negative part is finite). Define:

\[
\nu(E) = \int_E f \, d\mu
\]

Then $\nu$ is a signed measure.

---

## 2. Positive, Negative, and Null Sets

For signed measures, since the measure of a set can be canceled out by positive and negative components, simply knowing $\nu(E) > 0$ does not imply that $E$ is "positive" internally (it might just mean the positive part is larger than the negative part). Therefore, we need to strictly define what a "purely" positive set is.

!!! info "Definition 4.1.2 (Positive, Negative, and Null Sets)"

    Let $\nu$ be a signed measure on $(X, \mathcal{M})$, and $E \in \mathcal{M}$.

    * **Positive Set**: If for **every** measurable subset $F \subset E$ ($F \in \mathcal{M}$), we have $\nu(F) \ge 0$, then $E$ is called a positive set for $\nu$.

    * **Negative Set**: If for every measurable subset $F \subset E$, we have $\nu(F) \le 0$, then $E$ is called a negative set for $\nu$.

    * **Null Set**: If for every measurable subset $F \subset E$, we have $\nu(F) = 0$, then $E$ is called a null set for $\nu$.

Obviously, a null set is both a positive set and a negative set. At the same time, it is important to distinguish between a "null set" and a "set with measure zero": a set with measure zero might contain subsets with positive and negative measures that cancel each other out, but all subsets of a null set must strictly have measure 0.

!!! success "Lemma 4.1.1"

    The union of any countable collection of positive sets is still a positive set.

??? proof "Proof of Lemma 4.1.1"

    Let $\{P_j\}_{j=1}^\infty$ be a sequence of positive sets. We wish to prove that $P = \bigcup_{j=1}^\infty P_j$ is also a positive set.
    
    First, we orthogonalize (disjointify) this sequence of sets. Define:

    \[
    E_n = P_n \setminus \bigcup_{j=1}^{n-1} P_j
    \]

    Clearly $E_n \subset P_n$. Since any measurable subset of a positive set remains a positive set, all $E_n$ are positive sets, they are pairwise disjoint, and $\bigcup_{n=1}^\infty E_n = P$.

    Now take any $F \subset P$ with $F \in \mathcal{M}$. We can decompose $F$ as $F = \bigcup_{n=1}^\infty (F \cap E_n)$.
    Since $F \cap E_n \subset E_n$ and $E_n$ is a positive set, $\nu(F \cap E_n) \ge 0$.

    Using the countable additivity of signed measures:

    \[
    \nu(F) = \sum_{n=1}^\infty \nu(F \cap E_n) \ge 0
    \]

    This proves that $P$ is a positive set.

---

## 3. Hahn Decomposition Theorem

The Hahn Decomposition Theorem points to an extremely beautiful geometric picture: for any signed measure, we can always "cleanly" cut the entire space into a completely positive half and a completely negative half, with no ambiguous zones.

!!! success "Theorem 4.1.1 (Hahn Decomposition Theorem)"

    Let $\nu$ be a signed measure on the measurable space $(X, \mathcal{M})$.
    
    Then there must exist a positive set $P$ and a negative set $N$ satisfying:

    * $P \cup N = X$
    
    * $P \cap N = \emptyset$

    This pair of decomposition $(P, N)$ is called the **Hahn decomposition** of space $X$ with respect to $\nu$. Furthermore, this decomposition is unique up to a $\nu$-null set.

??? proof "Proof of Hahn Decomposition Theorem (Core Idea Outline)"

    Without loss of generality, we assume the signed measure $\nu$ does not take the value $+\infty$ (if it does not take $-\infty$, simply prove for $-\nu$). This means $\nu$ is bounded from above on all measurable sets.

    **Step 1: Extract the largest positive set**
    
    Define $m = \sup\{\nu(A) \mid A \text{ is a positive set}\}$. Since the empty set is a positive set and $\nu(\emptyset) = 0$, we have $m \ge 0$.
    By the definition of supremum, there exists a sequence of positive sets $\{P_j\}$ such that $\lim_{j \to \infty} \nu(P_j) = m$.
    Let $P = \bigcup_{j=1}^\infty P_j$. By Lemma 4.1.1, $P$ is also a positive set. Therefore $\nu(P) \le m$.
    Since $P_j \subset P$ and $P$ is a positive set, $\nu(P) = \nu(P_j) + \nu(P \setminus P_j) \ge \nu(P_j)$. Letting $j \to \infty$, we obtain $\nu(P) \ge m$.
    In summary, $\nu(P) = m$. Since $\nu$ does not take $+\infty$, $m < \infty$.

    **Step 2: Prove the remaining part is a negative set**
    
    Let $N = X \setminus P$. We claim that $N$ is a negative set.
    Suppose for contradiction that $N$ is not a negative set, then $N$ must contain a subset $E$ with measure strictly greater than 0. However, $E$ might not be a positive set (it might contain negative components).
    Through a process of continuously removing negative subsets (an exhaustion method), we can extract a true positive set $A \subset E$ from $E$ with $\nu(A) > 0$.
    
    Now, consider the set $P \cup A$. Since $P$ is a positive set, $A$ is a positive set, and they are disjoint, $P \cup A$ is a new positive set.
    Calculate its measure:

    \[
    \nu(P \cup A) = \nu(P) + \nu(A) = m + \nu(A) > m
    \]

    This contradicts the fact that $m$ is the supremum of the measures of all positive sets!
    Therefore, $N$ cannot contain any subset with positive measure; $N$ must be a negative set.

    **Step 3: Uniqueness**
    
    Suppose there exists another Hahn decomposition $(P', N')$. Consider $P \setminus P'$.
    Since $P \setminus P' \subset P$, it is a positive set with measure $\ge 0$.
    At the same time, $P \setminus P' \subset N'$ (since the space is split complementarily), so it is also a negative set with measure $\le 0$.
    This implies $P \setminus P'$ must be a null set. Similarly, all other differences between the parts can be proved to be null sets.

---

## 4. Jordan Decomposition Theorem and Total Variation

Using the Hahn decomposition, we can naturally split a signed measure into the difference of two independent, non-negative ordinary measures.

First, we need to define the "independence" between measures, namely **Mutually Singular**.

!!! info "Definition 4.1.3 (Mutually Singular Measures)"

    Let $\mu_1, \mu_2$ be two non-negative measures on the same measurable space. If there exist disjoint sets $E, F \in \mathcal{M}$ such that:

    * $E \cup F = X$

    * $\mu_1(F) = 0$ and $\mu_2(E) = 0$

    Then the measures $\mu_1$ and $\mu_2$ are said to be **mutually singular**, denoted as $\mu_1 \perp \mu_2$. Intuitively, this means the two measures "live" on disjoint territories.

Based on the concept of mutual singularity, we introduce the Jordan decomposition:

!!! success "Theorem 4.1.2 (Jordan Decomposition Theorem)"

    Any signed measure $\nu$ can be **uniquely** decomposed into the difference of two mutually singular non-negative measures, namely:

    \[
    \nu = \nu^+ - \nu^-
    \]

    And $\nu^+ \perp \nu^-$. This is called the **Jordan decomposition** of $\nu$. $\nu^+$ is called the positive variation, and $\nu^-$ is called the negative variation.

??? proof "Construction and Proof of Jordan Decomposition Theorem"

    Let $(P, N)$ be the Hahn decomposition of $X$ with respect to $\nu$.
    For any measurable set $E \in \mathcal{M}$, we define:

    \[
    \nu^+(E) = \nu(E \cap P)
    \]

    \[
    \nu^-(E) = -\nu(E \cap N)
    \]

    Clearly:
    1. Since $P$ is a positive set and $N$ is a negative set, $\nu^+$ and $\nu^-$ are both non-negative measures.
    2. For any $E \in \mathcal{M}$, we have $\nu(E) = \nu(E \cap P) + \nu(E \cap N) = \nu^+(E) - \nu^-(E)$.
    3. Taking $E=P, F=N$, obviously $P \cup N = X$. Since $\nu^-(P) = -\nu(P \cap N) = \nu(\emptyset) = 0$, and similarly $\nu^+(N) = 0$, this proves $\nu^+ \perp \nu^-$.
    
    Regarding uniqueness, any decomposition satisfying $\nu = \lambda^+ - \lambda^-$ and $\lambda^+ \perp \lambda^-$ must have support sets that constitute a Hahn decomposition of the original signed measure. From the uniqueness of the Hahn decomposition (up to a null set), the uniqueness of the Jordan decomposition measures follows.

!!! info "Definition 4.1.4 (Total Variation)"

    For a signed measure $\nu = \nu^+ - \nu^-$, we define its **total variation** as:

    \[
    |\nu| = \nu^+ + \nu^-
    \]

    $|\nu|$ is a standard non-negative measure. If $|\nu|(X) < \infty$, then $\nu$ is called a finite signed measure.
    For any set $E$, the total variation measure satisfies $|\nu(E)| \le |\nu|(E)$.

### Integration of Signed Measures

For a signed measure $\nu = \nu^+ - \nu^-$, since $\nu^+$ and $\nu^-$ are both ordinary non-negative measures, we can naturally define integration with respect to a signed measure.

Let $f$ be a measurable function. If $\int |f| d|\nu| < \infty$ (i.e., $f \in L^1(|\nu|)$), then the integral of $f$ with respect to $\nu$ is defined as:

\[
\int_X f \, d\nu = \int_X f \, d\nu^+ - \int_X f \, d\nu^-
\]

Thus, it can be seen that the space of integrable functions is entirely determined by the total variation measure, i.e., $L^1(\nu) = L^1(|\nu|)$.