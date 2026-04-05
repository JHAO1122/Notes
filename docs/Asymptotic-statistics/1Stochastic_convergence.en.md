# Chapter 1: Data Distributions and Stochastic Convergence

In statistical inference, finding the exact distribution $F_{\hat{\theta}_n}(x)$ of a statistic $\hat{\theta}_n$ under finite samples is often extremely difficult. However, by letting the sample size $n \to \infty$, we can utilize asymptotic theory to greatly simplify the problem and obtain exceptionally high-quality approximate distributions. This not only helps us construct approximate confidence intervals and hypothesis tests but also theoretically evaluates the asymptotic efficiency of different inference methods.

---

## 1. Basic Definitions of Stochastic Convergence

Let $\{X_n\}$ be a sequence of $\mathbb{R}^p$-valued random vectors defined on the same probability space $(\Omega, \mathcal{A}, P)$, and let $d(x,y)$ be the Euclidean distance.

!!! abstract "Definition: Four Types of Stochastic Convergence"

    **1. Almost-Sure Convergence ($X_n \xrightarrow{a.s.} X$)**:
    
    $$
    P(\lim_{n \to \infty} d(X_n, X) = 0) = 1
    $$
    
    *(Intuitive understanding: 100% certain and 100% accurate).*

    **2. Convergence in Probability ($X_n \xrightarrow{P} X$)**:
    
    For any given $\epsilon > 0$:
    
    $$
    \lim_{n \to \infty} P(d(X_n, X) < \epsilon) = 1
    $$
    
    *(Intuitive understanding: 100% certain, but not necessarily absolutely accurate).*

    **3. Convergence in $r$-th Mean ($X_n \xrightarrow{L^r} X$)**:
    
    $$
    \lim_{n \to \infty} E[d(X_n, X)^r] = 0
    $$

    **4. Convergence in Distribution / Weak Convergence ($X_n \xrightarrow{d} X$)**:
    
    Let $F_n$ and $F$ be the cumulative distribution functions (CDFs) of $X_n$ and $X$, respectively. If for all continuity points $x \in \mathcal{C}_F$ of $F$:
    
    $$
    \lim_{n \to \infty} F_n(x) = F(x)
    $$
    
    Then $X_n$ is said to converge in distribution to $X$.
    
    > **Note**: Convergence in distribution is the weakest form of convergence, and it is also the most central convergence in statistical inference. It **does not require** $X_n$ and $X$ to be defined on the same probability space. Moreover, the set of discontinuity points $\mathcal{C}_F^c$ is at most countable.

---

## 2. Polya's Theorem and Asymptotic Normality

In calculus, a continuous function on a closed interval is uniformly continuous. This can be generalized to cumulative distribution functions over the entire space:

!!! info "Lemma 1.2"
    
    If $F$ is a continuous distribution function on $\mathbb{R}$, then $F$ is uniformly continuous on $\mathbb{R}$.

Based on this, we can obtain an extremely elegant conclusion that strengthens convergence in distribution:

!!! success "Theorem 1.3: Polya's Theorem"

    Suppose $X_n \xrightarrow{d} X$, and the cumulative distribution function (CDF) $F(x)$ of the limiting random variable $X$ is **continuous**. 
    Then, this pointwise convergence automatically upgrades to **uniform convergence**:
    
    $$
    \sup_{x \in \mathbb{R}} |F_n(x) - F(x)| \rightarrow 0 \quad \text{as } n \rightarrow \infty
    $$
    
    ??? proof "Detailed Proof of Polya's Theorem (Click to expand)"
        
        We need to prove that for any $\epsilon > 0$, there exists an $N$ such that for all $n > N$, $\sup_{x \in \mathbb{R}} |F_n(x) - F(x)| < \epsilon$.

        **1. Constructing a Finite Partition**
        
        Since $F(x)$ is a continuous distribution function, its range is $[0, 1]$. For a given $\epsilon > 0$, we can find a finite number of points $-\infty = x_0 < x_1 < x_2 < \dots < x_K = \infty$ such that:
        
        $$
        F(x_i) - F(x_{i-1}) < \frac{\epsilon}{2}, \quad \forall i = 1, \dots, K
        $$
        
        (Note: We define $F(x_0) = 0$ and $F(x_K) = 1$).

        **2. Utilizing Pointwise Convergence**
        
        Since $X_n \xrightarrow{d} X$, for each of the finite grid points $x_i$ mentioned above (where $1 \le i \le K-1$), by the definition of convergence in distribution, as $n \to \infty$:
        
        $$
        F_n(x_i) \rightarrow F(x_i)
        $$
        
        Because there are only finitely many grid points, there must exist an $N$ such that for all $n > N$ and for all $i=1, \dots, K-1$:
        
        $$
        |F_n(x_i) - F(x_i)| < \frac{\epsilon}{2}
        $$

        **3. The Sandwich Argument via Monotonicity**
        
        For any arbitrary point $x \in \mathbb{R}$, it must fall within some interval $[x_{i-1}, x_i]$. Utilizing the non-decreasing property of $F_n$ and $F$:
        
        * **Upper Bound:**
        
        $$
        F_n(x) - F(x) \le F_n(x_i) - F(x_{i-1}) = [F_n(x_i) - F(x_i)] + [F(x_i) - F(x_{i-1})]
        $$
        
        When $n > N$, substituting the results from the previous two steps yields:
        
        $$
        F_n(x) - F(x) < \frac{\epsilon}{2} + \frac{\epsilon}{2} = \epsilon
        $$
        
        * **Lower Bound:**
        
        $$
        F_n(x) - F(x) \ge F_n(x_{i-1}) - F(x_i) = [F_n(x_{i-1}) - F(x_{i-1})] - [F(x_i) - F(x_{i-1})]
        $$
        
        When $n > N$, similarly we obtain:
        
        $$
        F_n(x) - F(x) > -\frac{\epsilon}{2} - \frac{\epsilon}{2} = -\epsilon
        $$

        **4. Conclusion**
        
        Combining the upper and lower bounds, for all $x \in \mathbb{R}$, as long as $n > N$, we have:
        
        $$
        |F_n(x) - F(x)| < \epsilon
        $$
        
        This proves that:
        
        $$
        \sup_{x \in \mathbb{R}} |F_n(x) - F(x)| \rightarrow 0 \quad (n \rightarrow \infty)
        $$
        
        $\square$

The most common type of convergence in statistics is convergence to a normal distribution:

!!! abstract "Definition 1.4 & 1.5: Asymptotic Normality (AN)"

    **1. Univariate Asymptotic Normality**:
    
    A sequence $\{X_n\}$ is said to have an asymptotic normal distribution with "mean" $\mu_n$ and "variance" $\sigma_n^2 > 0$, denoted as $X_n \sim AN(\mu_n, \sigma_n^2)$, if for sufficiently large $n$, it satisfies:
    
    $$
    \frac{X_n - \mu_n}{\sigma_n} \xrightarrow{d} N(0, 1)
    $$
    
    *(Note: Here $\mu_n$ and $\sigma_n^2$ are not necessarily the true mean and variance of $X_n$; sometimes the true moments of $X_n$ may not even exist!)*
    
    **2. Multivariate Asymptotic Normality**:
    
    For a sequence of random vectors $X_n$, it is said to follow $AN(\mu_n, \Sigma_n)$ if, for any non-zero vector $a \in \mathbb{R}^p$, its 1D projection satisfies:
    
    $$
    a^\top X_n \sim AN(a^\top \mu_n, a^\top \Sigma_n a)
    $$

---

## 3. The Cornerstone of Stochastic Convergence: Portmanteau Lemma

Besides the definition via CDF, convergence in distribution can also be equivalently described through expectations, open sets/closed sets, and other topological approaches. This constitutes the most core tool in asymptotic theory.

!!! info "Theorem 1.6: Portmanteau Lemma"

    For any random vectors $X_n$ and $X$, the following statements are **completely equivalent**:

    **(i)** $X_n \xrightarrow{d} X$;
    
    **(ii)** For any **bounded continuous function** $f \in C_B$, $E[f(X_n)] \to E[f(X)]$;
    
    **(iii)** For any **bounded Lipschitz continuous function** $f \in C_{B, Lip}$, $E[f(X_n)] \to E[f(X)]$;
    
    **(iv)** For any **non-negative continuous function** $f$, $\liminf E[f(X_n)] \ge E[f(X)]$;
    
    **(v)** For any **open set** $G$, $\liminf P(X_n \in G) \ge P(X \in G)$;
    
    **(vi)** For any **closed set** $F$, $\limsup P(X_n \in F) \le P(X \in F)$;
    
    **(vii)** For any Borel set $B$ with **boundary measure zero** (i.e., $P(X \in \partial B) = 0$), $P(X_n \in B) \to P(X \in B)$.

    ??? proof "Rigorous Derivation of the Core Steps in the Portmanteau Lemma (Click to expand)"

        **Proof (i) $\Rightarrow$ (ii)**:
        
        Without loss of generality, assume $\sup |f(x)| \le 1$. For any $\epsilon > 0$, choose a sufficiently large rectangular region $I$ such that the tail probability $P(X \in I^c) < \epsilon$.
        Partition $I$ into a finite number of non-overlapping small rectangles $I = \cup_{j=1}^K I_j$, and pick a representative point $x_j$ within each small rectangle. Construct a simple step function:
        
        $$
        f_\epsilon(x) = \sum_{j=1}^K f(x_j) \mathbb{I}(x \in I_j)
        $$
        
        Since $f$ is uniformly continuous on the compact set $I$, a sufficiently fine partition ensures that $|f(x) - f_\epsilon(x)| < \epsilon$ for all $x \in I$. By construction, we also have $\sup |f_\epsilon| \le \sup |f| \le 1$.
        
        Using indicator functions, we **decompose** the expectation over the entire space based on whether $X_n$ falls into the region $I$:
        
        $$
        \begin{aligned}
        |E[f(X_n)] - E[f_\epsilon(X_n)]| &\le E\left[ |f(X_n) - f_\epsilon(X_n)| \right] \\
        &= E\left[ |f(X_n) - f_\epsilon(X_n)| \cdot \mathbb{I}(X_n \in I) \right] + E\left[ |f(X_n) - f_\epsilon(X_n)| \cdot \mathbb{I}(X_n \in I^c) \right]
        \end{aligned}
        $$
        
        For the first term, since the error is bounded by $\epsilon$ within $I$:
        
        $$
        E\left[ |f(X_n) - f_\epsilon(X_n)| \cdot \mathbb{I}(X_n \in I) \right] < \epsilon \cdot P(X_n \in I) \le \epsilon
        $$
        
        For the second term, using $\sup |f| \le 1$ and $\sup |f_\epsilon| \le 1$, we have $|f - f_\epsilon| \le 2$:
        
        $$
        E\left[ |f(X_n) - f_\epsilon(X_n)| \cdot \mathbb{I}(X_n \in I^c) \right] \le 2 \cdot P(X_n \in I^c)
        $$
        
        Combining both parts yields the upper bound for the expectation error:
        
        $$
        |E[f(X_n)] - E[f_\epsilon(X_n)]| \le \epsilon + 2P(X_n \in I^c)
        $$
        
        For the simple function part:
        
        $$
        |Ef_\epsilon(X_n) - Ef_\epsilon(X)| \le \sum_{j=1}^K |f(x_j)| |P(X_n \in I_j) - P(X \in I_j)| \to 0
        $$
        
        Since $K$ is finite and we can construct the boundaries of $I_j$ to be continuity sets of $X$, combining these three components proves that $E[f(X_n)] \to E[f(X)]$.

        **Proof (iii) $\Rightarrow$ (v)**:
        
        For any open set $G$, we construct a sequence of non-negative Lipschitz functions to approximate its indicator function: let $f_m(x) = (m \cdot d(x, G^c)) \wedge 1$.
        As $m \to \infty$, $f_m \uparrow \mathbb{I}_G$. For a fixed $m$:
        
        $$
        \liminf_{n \to \infty} P(X_n \in G) \ge \liminf_{n \to \infty} E[f_m(X_n)] = E[f_m(X)]
        $$
        
        By the Monotone Convergence Theorem, as we let $m \to \infty$, the right side monotonically increases to $P(X \in G)$.

        **Proof (v) $\Leftrightarrow$ (vi)**:
        
        Utilizing the complementary relationship between open and closed sets (De Morgan's Laws), simply take the complement to reverse the direction of the inequality.

        **Proof (v) + (vi) $\Rightarrow$ (vii)**:
        
        Let $B^\circ$ be the interior of $B$ and $\overline{B}$ be its closure. Using the previous two properties:
        
        $$
        P(X \in B^\circ) \le \liminf P(X_n \in B) \le \limsup P(X_n \in B) \le P(X \in \overline{B})
        $$
        
        Since it is given that the boundary measure is zero, i.e., $P(X \in \partial B) = 0$, we have $P(X \in B^\circ) = P(X \in \overline{B})$. This squeezes the limit in the middle, proving that it exists and equals $P(X \in B)$. $\square$

        **Proof (vii) $\Rightarrow$ (i)**:
        
        For any real number $x$, construct a left-infinite closed interval $B = (-\infty, x]$. The boundary of this set is simply the singleton $\partial B = \{x\}$.
        If $x$ is a continuity point of the cumulative distribution function $F(x) = P(X \le x)$, then the probability measure at this point is zero, meaning $P(X \in \partial B) = P(X = x) = 0$.
        Since the boundary measure of this set is zero, condition (vii) implies:
        
        $$
        P(X_n \le x) = P(X_n \in B) \to P(X \in B) = P(X \le x)
        $$
        
        Because this equality holds for all continuity points $x$ of $F(x)$, it precisely matches the strict definition of convergence in distribution $X_n \xrightarrow{d} X$. $\square$


!!! success "Supplement to Theorem 1.6: Lévy's Continuity Theorem"

    Besides the equivalent topological conditions given by the Portmanteau Lemma, convergence in distribution has another extremely important and computationally valuable equivalent characterization: pointwise convergence of the **Characteristic Function**.
    
    Let $\{X_n\}$ and $X$ be random vectors in $\mathbb{R}^d$, and let $\phi_{X_n}(t)$ and $\phi_X(t)$ be their respective characteristic functions (defined as $\phi_X(t) = E[e^{i t^\top X}]$). Then:
    
    $$
    X_n \xrightarrow{d} X \iff \phi_{X_n}(t) \to \phi_X(t), \quad \forall t \in \mathbb{R}^d
    $$
    
    > *(Note: The convergence of characteristic functions is the most commonly used tool when proving asymptotic distributions such as the Central Limit Theorem (CLT)!)*

---

## 4. Continuous Mapping Theorem (CMT)

If a sequence of random variables is convergent, does the convergence property still hold when they are transformed by a "sufficiently good" function mapping? The mapping theorem provides an affirmative answer.

!!! success "Theorem 1.7: Continuous Mapping Theorem (Mapping Theorem)"

    Let the function $g: \mathbb{R}^k \to \mathbb{R}^m$ be continuous on the set of continuity points $\mathcal{C}_g$, and satisfy $P(X \in \mathcal{C}_g) = 1$ (i.e., $X$ almost surely falls on the continuity points of $g$).
    Then, the mapping operator $g(\cdot)$ perfectly inherits and transfers the following three types of convergence:
    
    1. If $X_n \xrightarrow{a.s.} X$, then $g(X_n) \xrightarrow{a.s.} g(X)$
    2. If $X_n \xrightarrow{P} X$, then $g(X_n) \xrightarrow{P} g(X)$
    3. If $X_n \xrightarrow{d} X$, then $g(X_n) \xrightarrow{d} g(X)$

    ??? proof "Rigorous Proof of the Mapping Theorem (Click to expand)"

        We focus here on proving the case for **convergence in distribution** $X_n \xrightarrow{d} X$. We will utilize the remarkably clever closed-set property (vi) of the Portmanteau Lemma.
        
        For any closed set $F$, consider its preimage $g^{-1}(F) = \{x : g(x) \in F\}$.
        Since $g$ is not necessarily continuous everywhere, we need to analyze the structure of the closure of the preimage $\overline{g^{-1}(F)}$:
        
        $$
        g^{-1}(F) \subset \overline{g^{-1}(F)} \subset g^{-1}(F) \cup \mathcal{C}_g^c
        $$
        
        *(Explanation: If a limit point $x$ is a continuity point, i.e., $x \in \mathcal{C}_g$, then for a sequence $x_m \to x$, we must have $g(x_m) \to g(x)$. Since $F$ is a closed set, naturally $g(x) \in F$, hence $x \in g^{-1}(F)$).*
        
        Applying Portmanteau Lemma (vi) to this set:
        
        $$
        \limsup P(g(X_n) \in F) \le \limsup P\left(X_n \in \overline{g^{-1}(F)}\right)
        $$
        
        By the lemma's property, the above is less than or equal to the probability of the limit on the closure:
        
        $$
        \le P\left(X \in \overline{g^{-1}(F)}\right)
        $$
        
        Decompose this into continuous and non-continuous point parts:
        
        $$
        \le P(X \in g^{-1}(F)) + P(X \notin \mathcal{C}_g)
        $$
        
        According to the premise of the theorem, $P(X \notin \mathcal{C}_g) = 0$. Therefore:
        
        $$
        \limsup P(g(X_n) \in F) \le P(g(X) \in F)
        $$
        
        Using the reverse deduction of the Portmanteau Lemma (vi) $\Rightarrow$ (i), it immediately proves that $g(X_n) \xrightarrow{d} g(X)$. $\square$

        **Proof of mapping property for Convergence in Probability (ii)**:
        
        We need to prove: for any given $\epsilon > 0$, $P(|g(X_n) - g(X)| > \epsilon) \to 0$.
        
        Fix any $\epsilon > 0$. For any $\delta > 0$, we define a "bad set" $B_\delta$, which contains all points $x$ where the function values may experience severe abrupt changes:
        
        $$
        B_\delta = \{x : \exists y \text{ such that } |x - y| < \delta \text{ but } |g(x) - g(y)| > \epsilon\}
        $$
        
        Now examine the event $\{|g(X_n) - g(X)| > \epsilon\}$. If this event occurs, and the limit variable $X$ happens to *not* be in the "bad set" $B_\delta$ (i.e., $X \notin B_\delta$), then it must be because $|X_n - X| \ge \delta$.
        Using the law of total probability for bounding, we get:
        
        $$
        P(|g(X_n) - g(X)| > \epsilon) \le P(X \in B_\delta) + P(|X_n - X| \ge \delta)
        $$
        
        Next, we take the limit for both terms on the right side:
        
        **First term**: Since $g$ is continuous on $\mathcal{C}_g$, as $\delta \downarrow 0$, the intersection of the set $B_\delta$ and the continuity set $\mathcal{C}_g$ must be empty. Moreover, because we know $P(X \in \mathcal{C}_g) = 1$, as $\delta \to 0$, $P(X \in B_\delta) \to 0$.
        
        **Second term**: For any fixed $\delta > 0$, since it is known that $X_n \xrightarrow{P} X$, as $n \to \infty$, $P(|X_n - X| \ge \delta) \to 0$.
        
        Combining both terms, by letting $n \to \infty$ and then letting $\delta \downarrow 0$, we successfully prove that the original probability goes to 0, which means $g(X_n) \xrightarrow{P} g(X)$. $\square$

!!! tip "Classic Application Examples of the Mapping Theorem (Applications)"

    The Mapping Theorem is a "divine weapon" when deriving the asymptotic distributions of complex statistics:

    1. **Derivation of the Chi-Square Distribution**:
       If a 1D sequence $X_n \xrightarrow{d} X \sim N(0,1)$, applying the continuous mapping $g(x) = x^2$ immediately yields $X_n^2 \xrightarrow{d} \chi_1^2$.
       
    2. **Derivation of the Cauchy Distribution**:
       If a 2D sequence $(X_n, Y_n)^\top \xrightarrow{d} N_2(0, I_2)$, taking the mapping $g(x,y) = x/y$ (which is discontinuous at $y=0$, but under the standard normal distribution $P(Y=0)=0$, satisfying the almost sure continuity condition), then $X_n/Y_n \xrightarrow{d} \text{Cauchy}$.
       
    3. **Convergence in Probability of the Sample Variance**:
       By the Law of Large Numbers, $(\overline{X}, \frac{1}{n}\sum X_i^2)^\top \xrightarrow{P} (\mu, \mu_2)^\top$. Taking the continuous function $g(x,y) = y - x^2$, we directly obtain the sample variance $S_n^2 = g(\overline{X}, \frac{1}{n}\sum X_i^2) \xrightarrow{P} \mu_2 - \mu^2 = \sigma^2$.
       
    4. **Affine Transformation of Multivariate Normal**:
       If $X_n \xrightarrow{d} N_p(\mu, \Sigma)$, for any constant matrix $C \in \mathbb{R}^{m \times p}$, we have $C X_n \xrightarrow{d} N_m(C\mu, C\Sigma C^\top)$.

## 5. Relations of Stochastic Convergence

There are strictly defined strong and weak implications among the four types of stochastic convergence. In practical applications, we often need to utilize these relations to transition between convergence in probability and convergence in distribution.

!!! success "Theorem 1.8: Relations of Stochastic Convergence"

    Let $X_n$, $X$, and $Y_n$ be random vectors. Then the following corollaries hold:

    1. **Almost-Sure $\Rightarrow$ in Probability**: If $X_n \xrightarrow{a.s.} X$, then $X_n \xrightarrow{P} X$.
    2. **in Probability $\Rightarrow$ in Distribution**: If $X_n \xrightarrow{P} X$, then $X_n \xrightarrow{d} X$.
    3. **Equivalence for Convergence to a Constant**: $X_n \xrightarrow{P} c$ (where $c$ is a constant) if and only if $X_n \xrightarrow{d} c$.
    4. **Distance Convergence Transfer (Convergence Lemma)**: If $X_n \xrightarrow{d} X$ and $d(X_n, Y_n) \xrightarrow{P} 0$, then $Y_n \xrightarrow{d} X$.
    5. **Joint Distribution Convergence (1)**: If $X_n \xrightarrow{d} X$ and $Y_n \xrightarrow{P} c$ (constant), then $(X_n, Y_n) \xrightarrow{d} (X, c)$.
    6. **Joint Distribution Convergence (2)**: If $X_n \xrightarrow{P} X$ and $Y_n \xrightarrow{P} Y$, then $(X_n, Y_n) \xrightarrow{P} (X, Y)$.

    > **Important Note**: Marginal convergence in probability implies joint convergence in probability. However, **marginal convergence in distribution generally does NOT imply joint convergence in distribution** (unless Copulas are used or one converges to a constant, as in property 5).

    ??? proof "Rigorous Proof of Core Properties (Click to expand)"

        **(1) Almost-Sure $\Rightarrow$ in Probability**:
        Define the sequence of events $A_n = \cup_{m \ge n} \{||X_m - X|| > \epsilon\}$. This sequence of sets is monotonically decreasing.
        If $X_n(\omega) \to X(\omega)$ for all $\omega \in \Omega$, then $A_n$ decreases to the empty set as $n \to \infty$.
        If $X_n \xrightarrow{a.s.} X$, by the continuity of probability:
        
        $$
        P(||X_n - X|| > \epsilon) \le P(A_n) \to 0
        $$

        **(4) Distance Convergence Transfer (Foundation for proofs)**:
        We prove this using Portmanteau Lemma (iii). For any bounded Lipschitz continuous function $f \in C_{B, Lip}$, let its Lipschitz constant be $L$ and $\sup |f| \le M$.
        Consider the absolute difference of the expectations, truncated by introducing an indicator function:
        
        $$
        \begin{aligned}
        |Ef(X_n) - Ef(Y_n)| &\le E\left[ |f(X_n) - f(Y_n)| \cdot \mathbb{I}(||X_n - Y_n|| \le \epsilon) \right] \\
        &\quad + E\left[ |f(X_n) - f(Y_n)| \cdot \mathbb{I}(||X_n - Y_n|| > \epsilon) \right]
        \end{aligned}
        $$
        
        For the first part, we bound it using the Lipschitz property as $L \epsilon \cdot P(||X_n - Y_n|| \le \epsilon) \le L \epsilon$;
        For the second part, we bound it using boundedness as $2M \cdot P(||X_n - Y_n|| > \epsilon)$.
        Therefore:
        
        $$
        |Ef(X_n) - Ef(Y_n)| \le L\epsilon + 2M \cdot P(||X_n - Y_n|| > \epsilon)
        $$
        
        Since $d(X_n, Y_n) \xrightarrow{P} 0$, the second term goes to 0 as $n \to \infty$. Because $\epsilon$ is arbitrarily small, $Ef(X_n) - Ef(Y_n) \to 0$.
        Combined with the premise $X_n \xrightarrow{d} X$ (i.e., $Ef(X_n) \to Ef(X)$), we must have $Ef(Y_n) \to Ef(X)$. By the Portmanteau Lemma again, it is proven that $Y_n \xrightarrow{d} X$.

        **(2) in Probability $\Rightarrow$ in Distribution**:
        Split $X_n$ into $X_n = X + (X_n - X)$.
        Since trivially $X \xrightarrow{d} X$, and we know $X_n - X \xrightarrow{P} 0$, we can directly apply the newly proven property (4) (letting $X_n$ in property 4 play the role of $X$, and $Y_n$ play the role of $X_n$), immediately yielding $X_n \xrightarrow{d} X$.

        **(3) Equivalence for Convergence to a Constant**:
        Sufficiency is guaranteed by property (2). To prove necessity ($X_n \xrightarrow{d} c \Rightarrow X_n \xrightarrow{P} c$):
        The event $\{||X_n - c|| \ge \epsilon\}$ is equivalent to $X_n \in B(c, \epsilon)^c$. This is a closed set.
        By Portmanteau Lemma (vi):
        
        $$
        \limsup P(||X_n - c|| \ge \epsilon) \le P(c \in B(c, \epsilon)^c) = 0
        $$
        
        Thus the limit is 0, proving convergence in probability.

        **(5) $(X_n, Y_n) \xrightarrow{d} (X, c)$**:
        Decompose the joint variables as $(X_n, Y_n) = (X_n, c) + (0, Y_n - c)$.
        Since $Y_n \xrightarrow{P} c$, the error term $(0, Y_n - c) \xrightarrow{P} (0,0)$.
        By property (4), we only need to prove the main term $(X_n, c) \xrightarrow{d} (X, c)$. For any bivariate continuous bounded function $f(x,y)$, fixing $y=c$ makes the marginal function $f_m(x) = f(x,c)$ also continuous and bounded.
        Because $X_n \xrightarrow{d} X$, we have $E[f(X_n, c)] = E[f_m(X_n)] \to E[f_m(X)] = E[f(X, c)]$. Proven!
        
        **(6) $(X_n, Y_n) \xrightarrow{P} (X, Y)$**:
        By the triangle inequality $||(X_n, Y_n) - (X, Y)|| \le ||X_n - X|| + ||Y_n - Y||$.
        
        $$
        P(||(X_n, Y_n) - (X, Y)|| > \epsilon) \le P(||X_n - X|| > \epsilon/2) + P(||Y_n - Y|| > \epsilon/2) \to 0
        $$
        
        Proven. $\square$

---

## 6. Slutsky's Theorem

As a direct corollary of Theorem 1.8 and the Continuous Mapping Theorem (CMT), Slutsky's Theorem provides highly convenient criteria for algebraic operations between statistics. It is the cornerstone for constructing asymptotic distributions like the $t$-statistic and Wald statistic.

!!! success "Theorem 1.9: Slutsky's Theorem (Slutsky, 1925)"

    Let $X_n, X, Y_n$ be random vectors or scalar random variables.
    If $X_n \xrightarrow{d} X$ and $Y_n \xrightarrow{P} c$ (constant), then:

    1. **Addition Rule**: $X_n + Y_n \xrightarrow{d} X + c$
    2. **Multiplication Rule**: $Y_n X_n \xrightarrow{d} cX$
    3. **Division Rule**: $Y_n^{-1} X_n \xrightarrow{d} c^{-1}X$ (provided $c \neq 0$ and is invertible for matrices)

    ??? proof "Extremely Concise Proof (Click to expand)"

        According to property (5) of Theorem 1.8, given $X_n \xrightarrow{d} X$ and $Y_n \xrightarrow{P} c$, we can directly deduce their joint distribution convergence:
        
        $$
        (X_n, Y_n) \xrightarrow{d} (X, c)
        $$
        
        Next, construct the bivariate continuous mapping functions $g(x,y) = x+y$, $g(x,y) = yx$, and $g(x,y) = y^{-1}x$ respectively.
        Applying the Continuous Mapping Theorem (CMT) directly $g(X_n, Y_n) \xrightarrow{d} g(X, c)$, the three rules of Slutsky's Theorem are immediately proven. $\square$

!!! tip "Classic Application: Asymptotic Normality of the $t$-Statistic"

    Let $Y_1, \dots, Y_n$ be i.i.d. samples satisfying $E[Y_1] = 0, E[Y_1^2] = \sigma^2$.
    We want to derive the limit distribution of the $t$-statistic $t_n := \frac{\sqrt{n}\overline{Y}}{S_n}$, where $S_n^2$ is the sample variance.

    **Step 1: Analyze the Numerator**
    By the Central Limit Theorem (CLT), the asymptotic distribution of the mean is:
    
    $$
    \sqrt{n}\overline{Y} \xrightarrow{d} N(0, \sigma^2)
    $$
    
    **Step 2: Analyze the Denominator**
    By the Law of Large Numbers, the sample variance converges in probability to the population variance: $S_n^2 \xrightarrow{P} \sigma^2$. By the mapping theorem ($g(x)=\sqrt{x}$), we get $S_n \xrightarrow{P} \sigma$.

    **Step 3: Apply Slutsky's Theorem**
    Treat the numerator as $X_n$ and the denominator as $Y_n^{-1}$:
    
    $$
    \frac{\sqrt{n}\overline{Y}}{S_n} \xrightarrow{d} \frac{1}{\sigma} N(0, \sigma^2) \stackrel{d}{=} N(0, 1)
    $$
    
    This theoretically proves that under large samples, the $t$-test can approximately use the critical values of the standard normal distribution.

---

## 7. Tightness & Stochastic Boundedness

When studying whether a sequence of random variables converges, a core question is: "Will they escape to infinity?" This leads to the concept of tightness.

!!! abstract "Definition 1.10: Stochastically Bounded / Tight"

    A sequence $\{X_n\}$ is said to be **stochastically bounded** (or **tight**) if, for any given $\epsilon > 0$, there exists a finite constant $M_\epsilon > 0$ such that for all $n$:
    
    $$
    \sup_n P(||X_n|| > M_\epsilon) < \epsilon
    $$
    
    In asymptotic statistics, we typically denote this property with the big-O notation: **$X_n = O_p(1)$**.

A single random variable $X$ itself is inherently tight (since its distribution function satisfies $F(\infty)=1, F(-\infty)=0$). Furthermore, any set of **finite** random variables is also tight. What truly requires vigilance is the escape of an infinite sequence $\{X_n\}$.

!!! info "Theorem 1.11: Prohorov's Theorem"

    There is an extremely profound topological connection between convergence in distribution and tightness:

    1. If $X_n \xrightarrow{d} X$, then $\{X_n\}$ is **necessarily tight** (i.e., convergence implies boundedness).
    2. If $\{X_n\}$ is tight, then it **must have a subsequence converging in distribution** $\{X_{n_i}\}$ such that $X_{n_i} \xrightarrow{d} X$ as $n_i \to \infty$. (This can be viewed as the probability space generalization of the Bolzano-Weierstrass theorem in real analysis, which states that "every bounded sequence has a convergent subsequence").

    ??? proof "Rigorous Proof of Theorem 1.11 (1) (Click to expand)"

        Since the limit random variable $X$ is a single random variable, it is tight. For any $\epsilon > 0$, we can find a sufficiently large $M_\epsilon$ that satisfies $P(||X|| = M_\epsilon) = 0$ (avoiding discontinuity points), such that:
        
        $$
        P(||X|| \ge M_\epsilon) < \epsilon
        $$
        
        Consider the closed set $F = \{x : ||x|| \ge M_\epsilon\}$. By Portmanteau Lemma (vi) and $X_n \xrightarrow{d} X$:
        
        $$
        \limsup_{n \to \infty} P(||X_n|| \ge M_\epsilon) \le P(||X|| \ge M_\epsilon) < \epsilon
        $$
        
        Since the limit superior is strictly less than $\epsilon$, there must exist a positive integer $N$ such that for all $n \ge N$:
        
        $$
        P(||X_n|| \ge M_\epsilon) \le 2\epsilon
        $$
        
        For the preceding finite $N-1$ random variables $\{X_1, \dots, X_{N-1}\}$, since a finite set is necessarily tight, we can appropriately enlarge $M_\epsilon$ such that $P(||X_n|| > M_\epsilon) < \epsilon$ holds for all $n \in \mathbb{N}_+$. Thus, the sequence is globally tight. $\square$

---

## 8. Stochastic Order ($o_p$ and $O_p$)

To perform calculus-like Taylor expansions in probability theory, we need a notation system to describe "infinitesimal in probability" and "of the same order in probability".

!!! abstract "Definition 1.12: Stochastic small $o_p$ and large $O_p$"

    Let $\{a_n\}$ be a sequence of constants. For a sequence of random vectors $\{X_n\}$:

    * **Large $O_p$ (Stochastically of order $a_n$)**: If $\frac{X_n}{a_n} = O_p(1)$, meaning $\frac{X_n}{a_n}$ is tight (stochastically bounded), it is denoted as $X_n = O_p(a_n)$.
    * **Small $o_p$ (Stochastically of smaller order $a_n$)**: If $\frac{X_n}{a_n} \xrightarrow{P} 0$, meaning it decays faster than $a_n$, it is denoted as $X_n = o_p(a_n)$.
    
    > *Special case*: $X_n = o_p(1)$ means $X_n \xrightarrow{P} 0$.

When determining the order of a statistic, Chebyshev's Inequality is the most commonly used tool. For instance, if $E[T_n] = \mu_n$ and $Var(T_n) = \sigma_n^2$, then it is guaranteed that $T_n - \mu_n = O_p(\sigma_n)$.

**Rules of Calculus for $o_p$ and $O_p$:**
When performing asymptotic expansions, we can manipulate these stochastic symbols just like deterministic limits (reading from left to right):

* $o_p(1) + o_p(1) = o_p(1)$
* $O_p(1) + o_p(1) = O_p(1)$
* $O_p(1) \cdot o_p(1) = o_p(1)$
* $(1 + o_p(1))^{-1} = O_p(1)$
* $O_p(1) + O_p(1) = O_p(1)$
* $o_p(O_p(1)) = O_p(o_p(1)) = o_p(1)$

---

## 9. Lemma of Stochastic Plug-in

This is the core prerequisite lemma we use in statistical inference to derive the Delta Method and the asymptotic normality of Maximum Likelihood Estimators (MLE).

!!! success "Lemma 1.13: Lemma of Stochastic Plug-in"

    Let $R: \mathbb{R}^k \to \mathbb{R}$ be a real function satisfying $R(0)=0$. Let $\{X_n\}$ be a sequence of random vectors with probability limit 0 (i.e., $X_n \xrightarrow{P} 0$). Then for any $p > 0$:

    1. If the deterministic limit $R(h) = o(||h||^p)$ as $h \to 0$, then $R(X_n) = o_p(||X_n||^p)$.
    2. If the deterministic limit $R(h) = O(||h||^p)$ as $h \to 0$, then $R(X_n) = O_p(||X_n||^p)$.

    > *Note: This lemma is extremely powerful. It allows us to seamlessly convert the Taylor expansion remainder $o(|x|)$ from deterministic calculus directly into $o_p(|X_n|)$ under random variables, without worrying whether the function $R$ is continuous elsewhere.*

    ??? proof "Proof of $o_p$ and $O_p$ Truncation (Click to expand)"

        **Proof for (1) small $o_p$ case**:
        Define the auxiliary function:
        
        $$
        g(h) = \begin{cases} \frac{R(h)}{||h||^p}, & h \neq 0 \\ 0, & h = 0 \end{cases}
        $$
        
        Then we can write $R(X_n) = ||X_n||^p g(X_n)$. To prove $R(X_n) = o_p(||X_n||^p)$, we only need to show $g(X_n) \xrightarrow{P} 0$.
        Due to the premise $R(h) = o(||h||^p)$, as $h \to 0$, $g(h) \to 0 = g(0)$. That is, $g(h)$ is continuous at 0.
        Since it is known that $X_n \xrightarrow{P} 0$, applying the Continuous Mapping Theorem (CMT), we immediately get $g(X_n) \xrightarrow{P} g(0) = 0$. Proven!

        **Proof for (2) large $O_p$ case**:
        Similarly use the $g(h)$ from above. Since $R(h) = O(||h||^p)$, $g(h)$ is bounded near $h=0$.
        This means there exist $M > 0$ and $\delta > 0$ such that when $||h|| < \delta$, $|g(h)| \le M$.
        Now examine the event $\{|g(X_n)| > M\}$. If this event occurs, it indicates $X_n$ must have escaped the safe region of radius $\delta$. Therefore, the following set inclusion holds:
        
        $$
        \{\omega : |g(X_n(\omega))| > M\} \subset \{\omega : ||X_n(\omega)|| > \delta\}
        $$
        
        Taking the probability measure:
        
        $$
        P(|g(X_n)| > M) \le P(||X_n|| > \delta)
        $$
        
        Because $X_n \xrightarrow{P} 0$, for any given $\epsilon > 0$, when $n$ is sufficiently large, the right side $P(||X_n|| > \delta) < \epsilon$.
        Thus:
        
        $$
        P(|g(X_n)| > M) < \epsilon, \quad \forall n \ge N
        $$
        
        This perfectly satisfies the definition of being stochastically bounded $O_p(1)$, meaning $g(X_n) = O_p(1)$, and thereby $R(X_n) = O_p(||X_n||^p)$. $\square$

!!! tip "Comprehensive Application Example: Asymptotic Distribution of Sample Variance"

    Suppose $X_1, \dots, X_n$ are i.i.d. $\sim F(\mu, \sigma^2)$ and the fourth moment exists $E[X^4] < \infty$. Derive the limit distribution of the biased variance $S_n^2$.
    
    First, expand the variance:
    
    $$
    S_n^2 = \frac{1}{n}\sum_{i=1}^n (X_i - \overline{X})^2 = \frac{1}{n}\sum_{i=1}^n (X_i - \mu)^2 - (\overline{X} - \mu)^2
    $$
    
    For the first part, since the fourth moment exists, using the standard Central Limit Theorem:
    
    $$
    \sqrt{n}\left( \frac{1}{n}\sum_{i=1}^n (X_i - \mu)^2 - \sigma^2 \right) \xrightarrow{d} N(0, Var((X_i - \mu)^2)) := N(0, \nu^2)
    $$
    
    For the second part, by CLT we know $\overline{X} - \mu = O_p(n^{-1/2})$. Using the Lemma of Stochastic Plug-in (setting $p=2$), squaring it yields:
    
    $$
    (\overline{X} - \mu)^2 = O_p(n^{-1}) = o_p(n^{-1/2})
    $$
    
    Therefore, when multiplied by $\sqrt{n}$:
    
    $$
    \sqrt{n}(\overline{X} - \mu)^2 = o_p(1)
    $$
    
    Applying Slutsky's Theorem, the higher-order error term vanishes in the limit:
    
    $$
    \sqrt{n}(S_n^2 - \sigma^2) = \sqrt{n}\left( \frac{1}{n}\sum_{i=1}^n (X_i - \mu)^2 - \sigma^2 \right) - \sqrt{n}(\overline{X} - \mu)^2 \xrightarrow{d} N(0, \nu^2) - 0 = N(0, \nu^2)
    $$