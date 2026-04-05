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

    Assume $X_n \xrightarrow{d} X$, and the distribution function $F(x)$ of the limit random variable $X$ is **continuous**.
    Then, this pointwise convergence automatically upgrades to **uniform convergence**:
    
    $$
    \sup_x |F_n(x) - F(x)| \rightarrow 0 \quad \text{as } n \rightarrow \infty
    $$
    
    ??? proof "Rigorous Proof of Polya's Theorem (Click to expand)"
        
        This proof uses the classic **Covering Method**.
        
        For any $\epsilon > 0$, since $F(x)$ is a valid continuous CDF, $\lim_{x \to -\infty} F(x) = 0$ and $\lim_{x \to \infty} F(x) = 1$. We can choose a sufficiently large $M > 0$ such that the tail probabilities are extremely small:
        
        $$
        F(-M) < \frac{\epsilon}{2} \quad \text{and} \quad 1 - F(M) < \frac{\epsilon}{2}
        $$
        
        Because $F$ is continuous on the compact interval $[-M, M]$, it is uniformly continuous. Thus, we can partition $[-M, M]$ into a finite number of grid points $-M = x_0 < x_1 < \dots < x_K = M$ such that for all $i = 1, \dots, K$:
        
        $$
        F(x_i) - F(x_{i-1}) < \frac{\epsilon}{2}
        $$
        
        For any $x \in [-M, M]$, there exists an $i$ such that $x_{i-1} \le x \le x_i$. Due to the monotonicity of CDFs, we can bound the difference:
        
        $$
        F_n(x) - F(x) \le F_n(x_i) - F(x_{i-1}) = F_n(x_i) - F(x_i) + F(x_i) - F(x_{i-1}) \le |F_n(x_i) - F(x_i)| + \frac{\epsilon}{2}
        $$
        
        Similarly, for the lower bound:
        
        $$
        F_n(x) - F(x) \ge F_n(x_{i-1}) - F(x_{i}) \ge -|F_n(x_{i-1}) - F(x_{i-1})| - \frac{\epsilon}{2}
        $$
        
        Therefore, the maximum error on the interval $[-M, M]$ is bounded by the grid points:
        
        $$
        \sup_{x \in [-M, M]} |F_n(x) - F(x)| \le \max_{0 \le i \le K} |F_n(x_i) - F(x_i)| + \frac{\epsilon}{2}
        $$
        
        Because $X_n \xrightarrow{d} X$ and $F$ has no points of discontinuity, for these **finite** $K+1$ grid points, we must have $\max_{i} |F_n(x_i) - F(x_i)| \to 0$ as $n \to \infty$. 
        
        The tail regions $(-\infty, -M)$ and $(M, \infty)$ can be bounded similarly using the limits of $F(-M)$ and $1-F(M)$. Thus, by taking the limit superior:
        
        $$
        \limsup_{n \to \infty} \sup_x |F_n(x) - F(x)| \le \epsilon
        $$
        
        Since $\epsilon$ is arbitrary, the uniform convergence over the entire space is proven. $\square$

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
        
        By making the partition sufficiently fine, we can guarantee that $|f(x) - f_\epsilon(x)| < \epsilon$ for all $x \in I$.
        Using the triangle inequality to decompose the expectation error:
        
        $$
        |Ef(X_n) - Ef_\epsilon(X_n)| \le \epsilon + 2P(X_n \in I^c)
        $$
        
        Similarly, for the limit variable:
        
        $$
        |Ef(X) - Ef_\epsilon(X)| \le \epsilon + 2P(X \in I^c) < 3\epsilon
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