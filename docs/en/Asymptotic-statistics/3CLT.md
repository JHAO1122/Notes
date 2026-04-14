# Chapter 3: Central Limit Theorem (Part I)

Unlike the classical central limit theorem (which requires independent and identically distributed, i.i.d., variables), this chapter explores more general central limit theorems, specifically the case where random variables are **independent but not identically distributed (i.n.i.d.)**.

---

## 1. Double Arrays & A Lemma on Complex Limits

When dealing with sums of random variables from different distributions, we often represent them in the form of a "double array" (or triangular array).

!!! abstract "Definition 3.1: Double Array of Independent Random Vectors"

    For each $n \ge 1$, let $\{X_{n1}, X_{n2}, \dots, X_{nk_n}\}$ be a set of random vectors defined on a probability space $(\Omega_n, \mathcal{F}_n, P_n)$, such that for a given $n$, $X_{n1}, \dots, X_{nk_n}$ are mutually independent, and $k_n \to \infty$ as $n \to \infty$.
    Then $\{X_{nj} : 1 \le j \le k_n\}_{n \ge 1}$ is called a **double array of independent random vectors**.

    **Common notations in this chapter**:
    
    * Expectation: $\alpha_{nj} = E(X_{nj})$, row total expectation $\alpha_n = \sum_{j=1}^{k_n} E(X_{nj}) = \sum_{j=1}^{k_n} \alpha_{nj}$
    * Partial sum: $S_n = \sum_{j=1}^{k_n} X_{nj}$
    * Variance: $\sigma_{nj}^2 = Var(X_{nj})$, row total variance $\sigma_n^2 = \sum_{j=1}^{k_n} \sigma_{nj}^2$

To handle products of characteristic functions, we need to introduce a lemma on the limit of products of complex sequences.

!!! info "Lemma 3.2: Limit of Products of Complex Sequences"

    Let $\{\theta_{nj} : 1 \le j \le k_n\}_{n \ge 1}$ be a double array of complex numbers satisfying, as $n \to \infty$:
    
    **(i)** Uniform convergence to 0: $\max_{1 \le j \le k_n} |\theta_{nj}| \to 0$
    **(ii)** Uniformly bounded absolute sum: $\sum_{j=1}^{k_n} |\theta_{nj}| \le M < \infty$ (where $M$ is independent of $n$)
    **(iii)** Sum convergence: $\sum_{j=1}^{k_n} \theta_{nj} \to \theta$ (where $\theta$ is a finite complex number)
    
    Then their product converges to the exponential function:
    
    $$
    \prod_{j=1}^{k_n} (1 + \theta_{nj}) \rightarrow e^\theta
    $$
    
    > *(Note: This generalizes the classical calculus result $\lim_{n \to \infty} (1 + \theta/n)^n = e^\theta$, which corresponds to substituting $\theta_{nj} \equiv \theta/n$ into this lemma.)*

    ??? proof "Detailed Proof of Lemma 3.2 (Click to expand)"

        For a non-zero complex number $z$, the principal value of the complex logarithm is defined as $Log~z = \log|z| + i Arg~z$, where $Arg~z \in [-\pi, \pi]$.
        When $|z| < 1$, the complex logarithm has the following Taylor series expansion:
        
        $$
        \log(1+z) = \sum_{m=1}^\infty (-1)^{m-1} \frac{z^m}{m}
        $$
        
        By condition (i), there exists $n_0$ such that for all $n > n_0$, $\max_{1 \le j \le k_n} |\theta_{nj}| \le 1/2$.
        Then $|\theta_{nj}| < 1$ and $1+\theta_{nj} \neq 0$. Consider the truncation error between its logarithmic expansion and the linear term:
        
        $$
        |\log(1+\theta_{nj}) - \theta_{nj}| = \left| \sum_{m \ge 2} (-1)^{m-1} \frac{\theta_{nj}^m}{m} \right| \le \sum_{m \ge 2} \frac{|\theta_{nj}|^m}{m}
        $$
        
        Extracting the quadratic term and bounding the remainder using a geometric series:
        
        $$
        \le \frac{|\theta_{nj}|^2}{2} \sum_{m=2}^\infty \left(\frac{1}{2}\right)^{m-2} = |\theta_{nj}|^2 < 1
        $$
        
        Since the absolute error is bounded by $|\theta_{nj}|^2$, we can write:
        
        $$
        \log(1+\theta_{nj}) = \theta_{nj} + \Lambda_{nj}|\theta_{nj}|^2, \quad \text{where } |\Lambda_{nj}| < 1
        $$
        
        Summing over a row:
        
        $$
        \sum_{j=1}^{k_n} \log(1+\theta_{nj}) = \sum_{j=1}^{k_n} \theta_{nj} + \sum_{j=1}^{k_n} \Lambda_{nj}|\theta_{nj}|^2
        $$
        
        Using conditions (i) and (ii), estimate the total error term:
        
        $$
        \left| \sum_{j=1}^{k_n} \Lambda_{nj}|\theta_{nj}|^2 \right| \le \max_{1 \le j \le k_n} |\theta_{nj}| \sum_{j=1}^{k_n} |\theta_{nj}| \le \left( \max_{1 \le j \le k_n} |\theta_{nj}| \right) \cdot M \rightarrow 0
        $$
        
        Combined with condition (iii) $\sum_{j=1}^{k_n} \theta_{nj} \to \theta$, we obtain:
        
        $$
        \sum_{j=1}^{k_n} \log(1+\theta_{nj}) \rightarrow \theta
        $$
        
        Taking the complex exponential on both sides completes the proof. $\square$

---

## 2. Liapounov's Central Limit Theorem (Liapounov's CLT)

If the sequence of random variables possesses moments higher than the second order, we can provide a very easily verifiable sufficient condition.

!!! success "Theorem 3.3: Liapounov's CLT (Liapounov's Theorem)"

    For the double array $\{X_{nj} : 1 \le j \le k_n\}_{n \ge 1}$, define the sum of their third-order central absolute moments as $\Gamma_n = \sum_{j=1}^{k_n} E|X_{nj} - \alpha_{nj}|^3$, assuming this value is finite for each $n$.
    If the **Liapounov's Condition** is satisfied:
    
    $$
    \frac{\Gamma_n}{\sigma_n^3} = \frac{1}{\sigma_n^3} \sum_{j=1}^{k_n} E|X_{nj} - \alpha_{nj}|^3 \rightarrow 0 \quad \text{as } n \to \infty
    $$
    
    then the standardized sum converges in distribution to the standard normal:
    
    $$
    \frac{S_n - \alpha_n}{\sigma_n} \rightarrow N(0,1)
    $$
    
    > *(Note: The third-order moment can be relaxed to a $2+\delta$ order moment, where $\delta > 0$.)*

    ??? proof "Rigorous Proof of Theorem 3.3 (click to expand)"

        Let $\gamma_{nj} = E|X_{nj} - \alpha_{nj}|^3$. By Liapounov's moment inequality:
        
        $$
        \sigma_{nj} = (E|X_{nj} - \alpha_{nj}|^2)^{1/2} \le (E|X_{nj} - \alpha_{nj}|^3)^{1/3} = \gamma_{nj}^{1/3}
        $$
        
        So $\sigma_{nj}^3 \le \gamma_{nj}$. Consequently:
        
        $$
        \max_{1 \le j \le k_n} \sigma_{nj}^3 \le \max_{1 \le j \le k_n} \gamma_{nj} \le \sum_{j=1}^{k_n} \gamma_{nj} = \Gamma_n
        $$
        
        Let $\phi_{nj}(t)$ be the characteristic function of the standardized variable $(X_{nj} - \alpha_{nj})/\sigma_n$. Since $\gamma_{nj}$ is finite, the characteristic function can be expanded to the third order via Taylor series:
        
        $$
        \phi_{nj}(t) = 1 - \frac{\sigma_{nj}^2 t^2}{2\sigma_n^2} + \frac{\Lambda_{nj}}{6} \frac{\gamma_{nj} t^3}{\sigma_n^3}, \quad \text{where } |\Lambda_{nj}| < 1
        $$
        
        To apply Lemma 3.2, we set $\theta_{nj} = \phi_{nj}(t) - 1$ and verify its three conditions:
        
        **Verification of (i)**:
        
        $$
        \max_j |\phi_{nj}(t) - 1| \le \frac{t^2}{2\sigma_n^2} \max_j \sigma_{nj}^2 + \frac{|t|^3}{6\sigma_n^3} \max_j \gamma_{nj}
        $$
        
        Since $\sigma_{nj}^2 = (\sigma_{nj}^3)^{2/3} \le (\max_j \sigma_{nj}^3)^{2/3}$, we have:
        
        $$
        \frac{\max \sigma_{nj}^2}{\sigma_n^2} \le \left( \frac{\max \sigma_{nj}^3}{\sigma_n^3} \right)^{2/3} \le \left( \frac{\Gamma_n}{\sigma_n^3} \right)^{2/3} \rightarrow 0
        $$
        
        Also $\max_j \gamma_{nj} / \sigma_n^3 \le \Gamma_n / \sigma_n^3 \rightarrow 0$. Therefore $\max_j |\theta_{nj}| \to 0$, condition (i) holds.
        
        **Verification of (ii)**:
        
        $$
        \sum_{j=1}^{k_n} |\phi_{nj}(t) - 1| \le \frac{t^2}{2\sigma_n^2} \sum_{j=1}^{k_n} \sigma_{nj}^2 + \frac{|t|^3}{6} \frac{\Gamma_n}{\sigma_n^3} = \frac{t^2}{2} + \frac{|t|^3}{6} \frac{\Gamma_n}{\sigma_n^3} \le M(t)
        $$
        
        This is a bounded quantity, condition (ii) holds.
        
        **Verification of (iii)**:
        
        Since the sum of error terms satisfies:
        
        $$
        \left| \sum_{j=1}^{k_n} \frac{\Lambda_{nj} \gamma_{nj}}{\sigma_n^3} \right| \le \frac{\Gamma_n}{\sigma_n^3} \rightarrow 0
        $$
        
        Therefore, the limit of the sum of characteristic function offsets is:
        
        $$
        \sum_{j=1}^{k_n} (\phi_{nj}(t) - 1) = -\frac{t^2}{2} + t^3 \sum_{j=1}^{k_n} \frac{\Lambda_{nj} \gamma_{nj}}{\sigma_n^3} \rightarrow -\frac{t^2}{2}
        $$
        
        Condition (iii) holds.
        
        In summary, according to Lemma 3.2, the characteristic function of the sum of independent standardized variables satisfies:
        
        $$
        \prod_{j=1}^{k_n} \phi_{nj}(t) = \prod_{j=1}^{k_n} (1 + (\phi_{nj}(t) - 1)) \rightarrow e^{-t^2/2}
        $$
        
        By the Lévy-Cramér continuity theorem, the result is proven. $\square$

For a general sequence with a single subscript, this corollary also applies:

!!! info "Corollary 3.4 (Single Sequence)"

    Let $\{X_n\}_{n \ge 1}$ be a sequence of independent random vectors. Let $\alpha_j = E(X_j)$, $\sigma_j^2 = Var(X_j)$ and $\gamma_j = E|X_j - \alpha_j|^3 < \infty$.
    Let $P_n = \sum_{j=1}^n \gamma_j$. If $P_n / \sigma_n^3 \rightarrow 0$, then:
    
    $$
    \frac{S_n - \sum_{j=1}^n \alpha_j}{\sigma_n} \rightarrow N(0,1)
    $$

---

## 3. Lindeberg's Telescoping Method

The CLT can also be proved directly using analytical techniques (a precursor to the Stein method) without using characteristic functions.

Assume $\alpha_j = 0$. Introduce a sequence of auxiliary random variables $Y_1, \dots, Y_n$ following a standard normal distribution, such that they are independent of the $X_j$, and $Y_j \sim N(0, \sigma_j^2)$ matches the first two moments.
Let $Y_0 = \sum_{i=1}^n Y_i / \sigma_n \sim N(0,1)$.

Our goal is to prove that for all bounded test functions $f \in C_B^\infty$ with bounded derivatives of all orders:

$$
E\left[f\left(\frac{\sum_{i=1}^n X_i}{\sigma_n}\right)\right] - E[f(Y_0)] \rightarrow 0
$$

We rely on the following theorem:

!!! success "Theorem 3.5 (Chung 6.1.6)"

    Let $\{\mu_n\}$ be a sequence of probability measures. If for all infinitely differentiable test functions $f \in C_B^\infty$ with bounded derivatives of all orders, we have:

    $$
    \int f(x) \mu_n(dx) \rightarrow \int f(x) \mu(dx)
    $$

    where $\mu$ is a probability measure, then $\mu_n$ converges weakly to $\mu$.

According to Theorem 3.5 (Chung 6.1.6), if the above expectation converges for all test functions, then the probability measures converge weakly.

!!! success "Core of the Proof Based on Telescoping Expansion"

    Construct the mixed partial sum sequence $Z_j$:
    $Z_j = Y_1 + \dots + Y_{j-1} + X_{j+1} + \dots + X_n$ (for $2 \le j \le n-1$).
    The boundaries are $Z_1 = X_2 + \dots + X_n$ and $Z_n = Y_1 + \dots + Y_{n-1}$.

    We write the total difference as a telescoping sum of term-by-term replacements:

    $$
    E\left[f\left(\frac{\sum X_i}{\sigma_n}\right)\right] - E\left[f\left(\frac{\sum Y_i}{\sigma_n}\right)\right] = \sum_{i=1}^n \left[ E\left[f\left(\frac{Z_i + X_i}{\sigma_n}\right)\right] - E\left[f\left(\frac{Z_i + Y_i}{\sigma_n}\right)\right] \right]
    $$

    Perform a third-order Taylor expansion of $f$ at $Z_i/\sigma_n$:

    $$
    f\left(\frac{Z_i + X_i}{\sigma_n}\right) = f\left(\frac{Z_i}{\sigma_n}\right) + f'\left(\frac{Z_i}{\sigma_n}\right)\frac{X_i}{\sigma_n} + \frac{1}{2}f''\left(\frac{Z_i}{\sigma_n}\right)\frac{X_i^2}{\sigma_n^2} + \theta_i^{(1)}\frac{X_i^3}{3!\sigma_n^3}
    $$

    Since $X_i$ and $Y_i$ have identical first two moments, the first and second-order terms cancel perfectly upon taking expectations.
    Only the third-order error term remains, bounded by the upper bound $M$ of the third derivative:

    $$
    \le \sum_{i=1}^n \frac{M}{3!\sigma_n^3} (\gamma_i + \sqrt{\frac{8}{\pi}}\sigma_i^3) \le \frac{M_1}{6}\frac{\Gamma_n}{\sigma_n^3} \rightarrow 0
    $$

    This method provides the foundation for the Stein method and Gaussian approximation for high-dimensional random vectors.

!!! Info "Corollary 3.6 Preliminary Theorem for Truncation Method"

    If there exists a double array $|X_{nj}/\sigma_n| \le M_{nj}$ a.s. and $\lim_{n \to \infty} \max_j M_{nj} = 0$, then the standardized sum converges to a normal distribution.

---

## 4. Null Arrays

To explore the sufficient and necessary conditions for the CLT to hold, we need to exclude pathological cases where a single variable dominates the overall variance.

!!! abstract "Definition 3.7: Null Array"

    A double array is called a **null array** if for any $\epsilon > 0$,
    
    $$
    \lim_{n \to \infty} \max_{1 \le j \le k_n} P(|X_{nj} - \alpha_{nj}| > \epsilon \sigma_n) = 0
    $$
    
    This is equivalent to saying each component $(X_{nj} - \alpha_j)/\sigma_n$ converges in probability to 0 uniformly in $j$ as $n \to \infty$.

Using characteristic functions, we can give a very convenient equivalent characterization:

!!! info "Proposition 3.8: Equivalent Form via Characteristic Functions"

    A double array $\{X_{nj}\}$ is a null array if and only if for all $t \in \mathbb{R}$:
    
    $$
    \lim_{n \to \infty} \max_{1 \le j \le k_n} |\phi_{nj}(t) - 1| = 0
    $$
    
    and this convergence is uniform on any finite interval $[-K, K]$.

    ??? proof "Proof of Equivalence (click to expand)"

        **($\Rightarrow$ direction)**: Without loss of generality, assume $\alpha_j = 0$.
        Decompose the expectation into parts inside and outside the threshold $\epsilon\sigma_n$:
        
        $$
        |\phi_{nj}(t) - 1| \le E\left[ |e^{itX_{nj}/\sigma_n} - 1| \mathbb{I}(|X_{nj}| > \epsilon\sigma_n) \right] + E\left[ |e^{itX_{nj}/\sigma_n} - 1| \mathbb{I}(|X_{nj}| \le \epsilon\sigma_n) \right]
        $$
        
        Using the inequality $|e^{itu} - 1| \le |tu|$ and the bound 2 for the modulus of the complex exponential:
        
        $$
        \le 2 P(|X_{nj}| > \epsilon\sigma_n) + |t| \epsilon
        $$
        
        Since we are on the bounded closed set $[-K, K]$:
        
        $$
        \sup_{|t|\le K} \max_j |\phi_{nj}(t) - 1| \le 2 \max_j P(|X_{nj}| > \epsilon\sigma_n) + K\epsilon
        $$
        
        As $n \to \infty$, the first term tends to 0 by the definition of a null array, and the second term can be made arbitrarily small because $\epsilon$ is arbitrary, thus proving uniform convergence.
        
        **($\Leftarrow$ direction)**:
        Using the classical characteristic function inequality:
        
        $$
        P\left(\left|\frac{X_{nj}}{\sigma_n}\right| > \frac{2}{\delta}\right) \le \frac{1}{\delta} \int_{|t|\le\delta} (1 - \phi_{nj}(t)) dt \le \frac{1}{\delta} \int_{|t|\le\delta} |1 - \phi_{nj}(t)| dt
        $$
        
        Taking the max on both sides:
        
        $$
        \max_j P\left(\left|\frac{X_{nj}}{\sigma_n}\right| > \frac{2}{\delta}\right) \le \max_j \frac{1}{\delta} \int_{|t|\le\delta} |1 - \phi_{nj}(t)| dt
        $$
        
        By the Bounded Convergence Theorem (BCT) for integrals and the given condition, the limit of the above expression tends to 0. The proposition is proved $\square$.

---

## 5. Lindeberg-Feller Central Limit Theorem (Lindeberg-Feller CLT)

When we only know that the second moment exists, while the third moment may not exist, the Lyapunov condition fails. At this point, the **Lindeberg Condition (LC)** becomes the most precise sufficient condition for the CLT of independent variables to hold.

!!! abstract "Definition 3.9: Lindeberg Condition (LC)"

    For a double array $\{X_{nj}\}$, if for any $\epsilon > 0$, its truncated variance ratio satisfies:
    
    $$
    \lim_{n \to \infty} \frac{1}{\sigma_n^2} \sum_{j=1}^{k_n} E \left[ (X_{nj} - \alpha_{nj})^2 \mathbb{I}(|X_{nj} - \alpha_{nj}| > \epsilon \sigma_n) \right] = 0
    $$
    
    then the array is said to satisfy the Lindeberg condition.

This, combined with the following lemma, leads to the famous theorem in the asymptotic distribution theory of independent variables:

!!! info "Lemma 3.10 (Diagonal Construction Method)"

    Let $u(m, n)$ be a function defined on positive integers $m$ and $n$, such that for each fixed $m$, we have:
    
    $$
    \lim_{n \to \infty} u(m, n) = 0
    $$
    
    Then there exists a monotonically increasing sequence $\{m_n\}$ tending to infinity ($m_n \to \infty$), such that:
    
    $$
    \lim_{n \to \infty} u(m_n, n) = 0
    $$

    ??? proof "Proof of Lemma 3.10 (click to expand)"

        Since for each fixed $m$, $\lim_{n \to \infty} u(m, n) = 0$, by the definition of limit, there must exist an index $n_m$ such that for all $n \ge n_m$, we have:
        
        $$
        u(m, n_m) \le \frac{1}{m}
        $$
        
        In this way, we can construct a strictly monotonic increasing sequence $\{n_m\}_{m \ge 1}$ that tends to infinity.
        
        For any $n$ satisfying $n_m \le n < n_{m+1}$, we set $m_n \equiv m$.
        
        Thus, when $n \ge n_m$, by construction we have:
        
        $$
        u(m_n, n) = u(m, n) \le \frac{1}{m}
        $$
        
        Since as $n \to \infty$, the index $m$ corresponding to $n_m$ tends to infinity, $m_n$ also monotonically increases to infinity.
        By the squeeze theorem from the above inequality, we have:
        
        $$
        \lim_{n \to \infty} u(m_n, n) = 0
        $$
        
        Proof complete. $\square$

!!! success "Theorem 3.11: Lindeberg-Feller CLT"

    Assume $Var(X_{nj}) = \sigma_{nj}^2 < \infty$, $S_n = \sum_{j=1}^{k_n} X_{nj}$. Then the following two sets of propositions are equivalent:
    
    * **(i)** $\frac{S_n - E S_n}{\sigma_n} \rightarrow N(0,1)$ and **(ii)** the double array is a null array.
    * **$\Longleftrightarrow$** The double array satisfies the **Lindeberg Condition (LC)**.

    ??? proof "Core Proof Derivation of the Theorem (click to expand)"

        **(1) Sufficiency Proof: LC $\Rightarrow$ CLT and Null Array**
        
        Assume $E(X_{nj})=0$ and $\sigma_n^2 = 1$. Define truncated random variables based on a truncation point $\eta \in (0, 1)$:
        $X_{nj}' = X_{nj}$ (if $|X_{nj}| < \eta$), otherwise 0.
        
        Compute the expectation and variance after truncation:
        
        $$
        |E[X_{nj}']| = \left| \int_{|x|<\eta} x dF_{nj}(x) \right| = \left| -\int_{|x|\ge\eta} x dF_{nj}(x) \right| \le \frac{1}{\eta} \int_{|x|\ge\eta} x^2 dF_{nj}(x)
        $$
        
        Summing over $j$, due to the LC condition, the total expectation tends to 0. Meanwhile, the truncated variance $\sigma_n'^2 \to 1 = \sigma_n^2$.
        According to Lemma 3.10 (diagonal construction principle), we can select a monotonically increasing sequence $m_n \to \infty$ and set $\eta_n = m_n^{-1} \to 0$.
        Using $\eta_n$ as the truncation threshold ensures $|X_{nj}'| \le \eta_n := M_{nj}$. Since $\max M_{nj} = \eta_n \to 0$, by the earlier Corollary 3.6 (CLT for bounded variables), we have $(S_n' - ES_n')/\sigma_n' \to N(0,1)$, hence $S_n' \to N(0,1)$.
        
        Finally, evaluate the difference between the untruncated sum and the truncated sum:
        
        $$
        P(S_n \neq S_n') \le \sum_{j=1}^{k_n} P(|X_{nj}| \ge \eta_n) \le \sum_{j=1}^{k_n} \frac{1}{\eta_n^2} \int_{|x| > \eta_n} x^2 dF_{nj}(x) \rightarrow 0 \quad \text{(guaranteed by LC)}
        $$
        
        By Slutsky's theorem, since the error is $o_p(1)$, we finally obtain $S_n \to N(0,1)$.
        
        **(2) Necessity Proof: CLT + Null Array $\Rightarrow$ LC**
        
        From $S_n \xrightarrow{d} N(0,1)$, we know the logarithm of the product of characteristic functions:
        
        $$
        \lim_{n \to \infty} \sum_{j=1}^{k_n} \log \phi_{nj}(t) = -\frac{t^2}{2}
        $$
        
        The null array guarantees $\max_j |\phi_{nj}(t) - 1| \to 0$. Using the equivalence between $\log \phi_{nj}(t)$ and $\phi_{nj}(t) - 1$:
        
        $$
        \lim_{n \to \infty} \sum_{j=1}^{k_n} (\phi_{nj}(t) - 1) = -\frac{t^2}{2}
        $$
        
        Extracting its real part:
        
        $$
        \lim_{n \to \infty} \sum_j \int_{-\infty}^{\infty} (1 - \cos tx) dF_{nj}(x) = \frac{t^2}{2}
        $$
        
        Split the integral into two parts: $|x| \le \eta$ and $|x| > \eta$, and use the inequality $0 \le 1 - \cos tx \le (tx)^2/2$:
        
        $$
        \frac{t^2}{2} - \sum_j \int_{|x|\le\eta} (1 - \cos tx) dF_{nj}(x) \ge \frac{t^2}{2} \sum_j \int_{|x| \ge \eta} x^2 dF_{nj}(x)
        $$
        
        By bounding the outer part of the integral controlled by Chebyshev's inequality ($\le 2/\eta^2 + \epsilon$), taking the limit as $t \to \infty$, we finally force $\sum_{j=1}^{k_n} E[X_{nj}^2 \mathbb{I}(|X_{nj}| > \eta)] \to 0$, which is precisely the Lindeberg condition LC. $\square$

---

## 6. Applications & Further Conditions

!!! tip "Application Example: Ordinary Least Squares Regression (OLS Regression)"

    Consider the classical linear regression model $y_j = x_j \beta + \epsilon_j$, where the error terms $\epsilon_j \sim (0, \sigma_\epsilon^2)$ and are i.i.d.
    The design matrix satisfies $\max_{1 \le j \le n} \frac{|x_j|}{a_n} \to 0$, where $a_n^2 = \sum_{j=1}^n x_j^2$.
    The OLS estimator is $\hat{\beta}_{LS} = \sum x_j y_j / a_n^2$.
    
    Construct the standardized double array $X_{nj} = \frac{x_j \epsilon_j}{\sqrt{\sum x_j^2}}$.
    We bound the truncated integral in the Lindeberg condition: let $m_n = \max_j |x_j/a_n|$,
    
    $$
    \frac{1}{\sigma_n^2 a_n^2} \sum_{j=1}^n x_j^2 E\left[ \epsilon_j^2 \mathbb{I}(|x_j \epsilon_j / a_n| > \delta \sigma_n) \right] \le \frac{1}{\sigma_\epsilon^2} E\left[ \epsilon_j^2 \mathbb{I}(|\epsilon_j| > \delta \sigma_\epsilon m_n^{-1}) \right]
    $$
    
    Since the $\epsilon_j$ are identically distributed and have finite second moments, this expectation tends to 0 as $m_n \to 0$. By Theorem 3.11, we immediately obtain:
    
    $$
    a_n(\hat{\beta}_{LS} - \beta) \xrightarrow{d} N(0, \sigma_\epsilon^2)
    $$

When verifying the LC condition, besides using bounded truncation, we can also utilize the existence of higher-order moments, which is an **extension of the Lyapunov-type condition**.

!!! info "Proposition 3.12: Sufficient Criterion for the Lindeberg Condition"

    For a double array $\{X_{nj}\}$, if there exists some real number $\nu > 2$ such that:
    
    $$
    \sum_{j=1}^{k_n} E|X_{nj} - \mu_{nj}|^\nu = o(\sigma_n^\nu)
    $$
    
    then the array necessarily satisfies the Lindeberg condition.

    ??? proof "Proof Derivation (click to expand)"

        In the truncation region $|t - \mu_{nj}| > \epsilon \sigma_n$, we bound the second-moment integral:
        
        $$
        E\left[ (X_{nj} - \mu_{nj})^2 \mathbb{I}(|X_{nj} - \mu_{nj}| > \epsilon \sigma_n) \right] = \int_{|t - \mu_{nj}| > \epsilon \sigma_n} (t - \mu_{nj})^2 dF_{nj}(t)
        $$
        
        By forcibly introducing the $\nu$-th power and extracting a constant factor:
        
        $$
        \le (\epsilon \sigma_n)^{2-\nu} \int_{|t - \mu_{nj}| > \epsilon \sigma_n} (t - \mu_{nj})^\nu dF_{nj}(t) \le (\epsilon \sigma_n)^{2-\nu} E|X_{nj} - \mu_{nj}|^\nu
        $$
        
        Summing over $j$ and dividing by $\sigma_n^2$:
        
        $$
        \frac{1}{\sigma_n^2} \sum_{j=1}^{k_n} E\left[ \dots \right] \le \epsilon^{2-\nu} \frac{\sum_{j=1}^{k_n} E|X_{nj} - \mu_{nj}|^\nu}{\sigma_n^\nu} \rightarrow 0
        $$
        
        The Lindeberg condition is thus proved. $\square$