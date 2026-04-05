# Chapter 3: Central Limit Theorem

Unlike the classical Central Limit Theorem (which requires independent and identically distributed, i.i.d. variables), in this chapter we will explore a more general Central Limit Theorem, where the random variables are **independent but not identically distributed (i.n.i.d.)**.

---

## 1. Double Arrays & A Lemma

When dealing with the sum of random variables with different distributions, we usually represent them in the form of a "double array" (or triangular array).

!!! abstract "Definition 3.1: Double Array of Independent Random Vectors"

    For each $n \ge 1$, let $\{X_{n1}, X_{n2}, \dots, X_{nk_n}\}$ be a collection of random vectors defined on a probability space $(\Omega_n, \mathcal{F}_n, P_n)$, such that for a given $n$, $X_{n1}, \dots, X_{nk_n}$ are mutually independent, and $k_n \to \infty$ as $n \to \infty$.
    Then $\{X_{nj} : 1 \le j \le k_n\}_{n \ge 1}$ is called a **double array of independent random vectors**.

    **Common notations in this chapter**:
    
    * Expectation: $\alpha_{nj} = E(X_{nj})$, row total expectation $\alpha_n = \sum_{j=1}^{k_n} E(X_{nj}) = \sum_{j=1}^{k_n} \alpha_{nj}$
    * Partial sum: $S_n = \sum_{j=1}^{k_n} X_{nj}$
    * Variance: $\sigma_{nj}^2 = Var(X_{nj})$, row total variance $\sigma_n^2 = \sum_{j=1}^{k_n} \sigma_{nj}^2$

To handle the product of characteristic functions, we need to introduce a limit lemma for products of complex sequences.

!!! info "Lemma 3.2: Limit of Products of Complex Sequences"

    Let $\{\theta_{nj} : 1 \le j \le k_n\}_{n \ge 1}$ be a double array of complex numbers satisfying as $n \to \infty$:
    
    **(i)** Uniformly tends to 0: $\max_{1 \le j \le k_n} |\theta_{nj}| \to 0$
    **(ii)** Absolute sum is uniformly bounded: $\sum_{j=1}^{k_n} |\theta_{nj}| \le M < \infty$ (where $M$ is independent of $n$)
    **(iii)** Sum converges: $\sum_{j=1}^{k_n} \theta_{nj} \to \theta$ (where $\theta$ is a finite complex number)
    
    Then its product converges to the exponential function:
    
    $$
    \prod_{j=1}^{k_n} (1 + \theta_{nj}) \rightarrow e^\theta
    $$
    
    > *(Note: This generalizes the classic calculus conclusion $\lim_{n \to \infty} (1 + \theta/n)^n = e^\theta$ by substituting $\theta_{nj} \equiv \theta/n$ into this lemma.)*

    ??? proof "Detailed Proof of Lemma 3.2 (Click to expand)"

        For a non-zero complex number $z$, the principal value of the complex logarithm is defined as $Log~z = \log|z| + i Arg~z$, where $Arg~z \in [-\pi, \pi]$.
        When $|z| < 1$, the complex logarithm has the following Taylor series expansion:
        
        $$
        \log(1+z) = \sum_{m=1}^\infty (-1)^{m-1} \frac{z^m}{m}
        $$
        
        From condition (i), there exists $n_0$ such that for all $n > n_0$, $\max_{1 \le j \le k_n} |\theta_{nj}| \le 1/2$.
        At this time $|\theta_{nj}| < 1$ and $1+\theta_{nj} \neq 0$. Consider its logarithmic expansion and the truncation error of the linear term:
        
        $$
        |\log(1+\theta_{nj}) - \theta_{nj}| = \left| \sum_{m \ge 2} (-1)^{m-1} \frac{\theta_{nj}^m}{m} \right| \le \sum_{m \ge 2} \frac{|\theta_{nj}|^m}{m}
        $$
        
        Extract the quadratic term, and use the geometric series sum to bound the remaining terms:
        
        $$
        \le \frac{|\theta_{nj}|^2}{2} \sum_{m=2}^\infty \left(\frac{1}{2}\right)^{m-2} = |\theta_{nj}|^2 < 1
        $$
        
        Since the absolute value of the error is bounded by $|\theta_{nj}|^2$, we can write it as:
        
        $$
        \log(1+\theta_{nj}) = \theta_{nj} + \Lambda_{nj}|\theta_{nj}|^2, \quad \text{where } |\Lambda_{nj}| < 1
        $$
        
        Summing over a row:
        
        $$
        \sum_{j=1}^{k_n} \log(1+\theta_{nj}) = \sum_{j=1}^{k_n} \theta_{nj} + \sum_{j=1}^{k_n} \Lambda_{nj}|\theta_{nj}|^2
        $$
        
        Using conditions (i) and (ii), estimate the sum of the error terms:
        
        $$
        \left| \sum_{j=1}^{k_n} \Lambda_{nj}|\theta_{nj}|^2 \right| \le \max_{1 \le j \le k_n} |\theta_{nj}| \sum_{j=1}^{k_n} |\theta_{nj}| \le \left( \max_{1 \le j \le k_n} |\theta_{nj}| \right) \cdot M \rightarrow 0
        $$
        
        Combined with condition (iii) $\sum_{j=1}^{k_n} \theta_{nj} \to \theta$, we have:
        
        $$
        \sum_{j=1}^{k_n} \log(1+\theta_{nj}) \rightarrow \theta
        $$
        
        Taking the complex exponential mapping on both sides, the lemma is proven. $\square$

---

## 2. Liapounov's CLT

If the sequence of random variables possesses moments higher than the second order, we can provide a very easily verifiable sufficient condition.

!!! success "Theorem 3.3: Liapounov's Theorem"

    For a double array $\{X_{nj} : 1 \le j \le k_n\}_{n \ge 1}$, define the sum of its third-order absolute central moments as $\Gamma_n = \sum_{j=1}^{k_n} E|X_{nj} - \alpha_{nj}|^3$, assuming this value is finite for each $n$.
    If **Liapounov's Condition** is satisfied:
    
    $$
    \frac{\Gamma_n}{\sigma_n^3} = \frac{1}{\sigma_n^3} \sum_{j=1}^{k_n} E|X_{nj} - \alpha_{nj}|^3 \rightarrow 0 \quad \text{as } n \to \infty
    $$
    
    Then the standardized sum will converge in distribution to the standard normal:
    
    $$
    \frac{S_n - \alpha_n}{\sigma_n} \rightarrow N(0,1)
    $$
    
    > *(Note: The third-order moment can be relaxed to the $(2+\delta)$-th order moment, where $\delta > 0$.)*

    ??? proof "Rigorous Proof of Theorem 3.3 (Click to expand)"

        Let $\gamma_{nj} = E|X_{nj} - \alpha_{nj}|^3$. From Liapounov's moment inequality we have:
        
        $$
        \sigma_{nj} = (E|X_{nj} - \alpha_{nj}|^2)^{1/2} \le (E|X_{nj} - \alpha_{nj}|^3)^{1/3} = \gamma_{nj}^{1/3}
        $$
        
        So $\sigma_{nj}^3 \le \gamma_{nj}$. Furthermore:
        
        $$
        \max_{1 \le j \le k_n} \sigma_{nj}^3 \le \max_{1 \le j \le k_n} \gamma_{nj} \le \sum_{j=1}^{k_n} \gamma_{nj} = \Gamma_n
        $$
        
        Let $\phi_{nj}(t)$ be the characteristic function of the standardized variable $(X_{nj} - \alpha_{nj})/\sigma_n$. Since $\gamma_{nj}$ is finite, the characteristic function can be expanded to the third-order Taylor series:
        
        $$
        \phi_{nj}(t) = 1 - \frac{\sigma_{nj}^2 t^2}{2\sigma_n^2} + \frac{\Lambda_{nj}}{6} \frac{\gamma_{nj} t^3}{\sigma_n^3}, \quad \text{where } |\Lambda_{nj}| < 1
        $$
        
        To apply Lemma 3.2, we let $\theta_{nj} = \phi_{nj}(t) - 1$, and verify its three conditions:
        
        **Verify (i)**:
        
        $$
        \max_j |\phi_{nj}(t) - 1| \le \frac{t^2}{2\sigma_n^2} \max_j \sigma_{nj}^2 + \frac{|t|^3}{6\sigma_n^3} \max_j \gamma_{nj}
        $$
        
        Since $\sigma_{nj}^2 = (\sigma_{nj}^3)^{2/3} \le (\max_j \sigma_{nj}^3)^{2/3}$, we have:
        
        $$
        \frac{\max \sigma_{nj}^2}{\sigma_n^2} \le \left( \frac{\max \sigma_{nj}^3}{\sigma_n^3} \right)^{2/3} \le \left( \frac{\Gamma_n}{\sigma_n^3} \right)^{2/3} \rightarrow 0
        $$
        
        Also $\max_j \gamma_{nj} / \sigma_n^3 \le \Gamma_n / \sigma_n^3 \rightarrow 0$. Thus $\max_j |\theta_{nj}| \to 0$, condition (i) holds.
        
        **Verify (ii)**:
        
        $$
        \sum_{j=1}^{k_n} |\phi_{nj}(t) - 1| \le \frac{t^2}{2\sigma_n^2} \sum_{j=1}^{k_n} \sigma_{nj}^2 + \frac{|t|^3}{6} \frac{\Gamma_n}{\sigma_n^3} = \frac{t^2}{2} + \frac{|t|^3}{6} \frac{\Gamma_n}{\sigma_n^3} \le M(t)
        $$
        
        This is a bounded quantity, condition (ii) holds.
        
        **Verify (iii)**:
        
        Since the sum of the error terms satisfies:
        
        $$
        \left| \sum_{j=1}^{k_n} \frac{\Lambda_{nj} \gamma_{nj}}{\sigma_n^3} \right| \le \frac{\Gamma_n}{\sigma_n^3} \rightarrow 0
        $$
        
        Therefore, the limit of the sum of the characteristic function offsets is:
        
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

Without using characteristic functions, directly applying analytical techniques (the predecessor of Stein's method) can also prove the CLT.

Assume $\alpha_j = 0$. Introduce an auxiliary sequence of random variables $Y_1, \dots, Y_n$ following standard normal distributions, such that they are independent of $X_j$, and $Y_j \sim N(0, \sigma_j^2)$ matching the first two moments.
Let $Y_0 = \sum_{i=1}^n Y_i / \sigma_n \sim N(0,1)$.

Our goal is to prove that for all bounded test functions $f \in C_B^\infty$ with bounded derivatives of all orders:

$$
E\left[f\left(\frac{\sum_{i=1}^n X_i}{\sigma_n}\right)\right] - E[f(Y_0)] \rightarrow 0
$$

According to Theorem 3.5 (Chung 6.1.6), if the above expectation converges for all test functions, then the probability measure converges weakly.

!!! success "Core of the Proof Based on Telescoping Expansion"

    Construct the mixed partial sum sequence $Z_j$:
    $Z_j = Y_1 + \dots + Y_{j-1} + X_{j+1} + \dots + X_n$ (for $2 \le j \le n-1$).
    The boundaries are $Z_1 = X_2 + \dots + X_n$ and $Z_n = Y_1 + \dots + Y_{n-1}$.

    We write the total difference as a telescoping sum of term-by-term replacements:
    
    $$
    E\left[f\left(\frac{\sum X_i}{\sigma_n}\right)\right] - E\left[f\left(\frac{\sum Y_i}{\sigma_n}\right)\right] = \sum_{i=1}^n \left[ E\left[f\left(\frac{Z_i + X_i}{\sigma_n}\right)\right] - E\left[f\left(\frac{Z_i + Y_i}{\sigma_n}\right)\right] \right]
    $$
    
    Perform a third-order Taylor expansion on $f$ at $Z_i/\sigma_n$:
    
    $$
    f\left(\frac{Z_i + X_i}{\sigma_n}\right) = f\left(\frac{Z_i}{\sigma_n}\right) + f'\left(\frac{Z_i}{\sigma_n}\right)\frac{X_i}{\sigma_n} + \frac{1}{2}f''\left(\frac{Z_i}{\sigma_n}\right)\frac{X_i^2}{\sigma_n^2} + \theta_i^{(1)}\frac{X_i^3}{3!\sigma_n^3}
    $$
    
    Since the first two moments of $X_i$ and $Y_i$ are completely identical, the first-order and second-order terms perfectly cancel out after taking expectations.
    Only the third-order error term remains, bounded by the supremum $M$ of the third derivative:
    
    $$
    \le \sum_{i=1}^n \frac{M}{3!\sigma_n^3} \left(\gamma_i + \sqrt{\frac{8}{\pi}}\sigma_i^3\right) \le \frac{M_1}{6}\frac{\Gamma_n}{\sigma_n^3} \rightarrow 0
    $$
    
    This method provides the cornerstone for Stein's method and Gaussian approximation of high-dimensional random vectors.

Corollary 3.6 Truncation Method Preparation Theorem:
If there exists a double array $|X_{nj}/\sigma_n| \le M_{nj}$ a.s. and $\lim_{n \to \infty} \max_j M_{nj} = 0$, then the standardized sum converges to a normal distribution.

---

## 4. Null Arrays

To explore the necessary and sufficient conditions for the CLT to hold, we need to rule out the pathological case where a single variable dominates the total variance.

!!! abstract "Definition 3.7: Null Array"

    If a double array satisfies: for any $\epsilon > 0$,
    
    $$
    \lim_{n \to \infty} \max_{1 \le j \le k_n} P(|X_{nj} - \alpha_{nj}| > \epsilon \sigma_n) = 0
    $$
    
    Then this double array is called a **Null array**.
    This is equivalent to saying that each component $(X_{nj} - \alpha_j)/\sigma_n$ uniformly degenerates to 0 in probability with respect to $j$ as $n \to \infty$.

Using characteristic functions, we can provide a very computationally convenient equivalent characterization:

!!! info "Proposition 3.8: Equivalent Characteristic Function Form for Null Arrays"

    A double array $\{X_{nj}\}$ is a null array if and only if for all $t \in \mathbb{R}$:
    
    $$
    \lim_{n \to \infty} \max_{1 \le j \le k_n} |\phi_{nj}(t) - 1| = 0
    $$
    
    And this convergence is uniform on any finite interval $[-K, K]$.

    ??? proof "Proof of Equivalence (Click to expand)"

        **($\Rightarrow$ Direction)**: Without loss of generality, assume $\alpha_j = 0$.
        Decompose the expectation into two parts, inside and outside the threshold $\epsilon\sigma_n$:
        
        $$
        |\phi_{nj}(t) - 1| \le E\left[ |e^{itX_{nj}/\sigma_n} - 1| \mathbb{I}(|X_{nj}| > \epsilon\sigma_n) \right] + E\left[ |e^{itX_{nj}/\sigma_n} - 1| \mathbb{I}(|X_{nj}| \le \epsilon\sigma_n) \right]
        $$
        
        Using the inequality $|e^{itu} - 1| \le |tu|$ and the bound 2 for the modulus of complex exponentials:
        
        $$
        \le 2 P(|X_{nj}| > \epsilon\sigma_n) + |t| \epsilon
        $$
        
        Because it is on a bounded closed set $[-K, K]$:
        
        $$
        \sup_{|t|\le K} \max_j |\phi_{nj}(t) - 1| \le 2 \max_j P(|X_{nj}| > \epsilon\sigma_n) + K\epsilon
        $$
        
        As $n \to \infty$, the first term goes to 0 by the definition of a null array, and the second term can be arbitrarily small since $\epsilon$ is arbitrary. Thus, uniform convergence is proven.
        
        **($\Leftarrow$ Direction)**:
        Using the classic characteristic function inequality:
        
        $$
        P\left(\left|\frac{X_{nj}}{\sigma_n}\right| > \frac{2}{\delta}\right) \le \frac{1}{\delta} \int_{|t|\le\delta} (1 - \phi_{nj}(t)) dt \le \frac{1}{\delta} \int_{|t|\le\delta} |1 - \phi_{nj}(t)| dt
        $$
        
        Taking max on both sides:
        
        $$
        \max_j P\left(\left|\frac{X_{nj}}{\sigma_n}\right| > \frac{2}{\delta}\right) \le \max_j \frac{1}{\delta} \int_{|t|\le\delta} |1 - \phi_{nj}(t)| dt
        $$
        
        By the Bounded Convergence Theorem (BCT) for integrals and the given conditions, the limit of the above expression tends to 0. The proposition is proven. $\square$

---

## 5. Lindeberg-Feller CLT

When we only know that the second moment exists, and the third moment may not, Liapounov's condition fails. At this time, the **Lindeberg Condition** becomes the most precise necessary and sufficient condition for the CLT of independent variables to hold.

!!! abstract "Definition 3.9: Lindeberg Condition (LC)"

    For a double array $\{X_{nj}\}$, if for any $\epsilon > 0$, its truncated variance ratio satisfies:
    
    $$
    \lim_{n \to \infty} \frac{1}{\sigma_n^2} \sum_{j=1}^{k_n} E \left[ (X_{nj} - \alpha_{nj})^2 \mathbb{I}(|X_{nj} - \alpha_{nj}| > \epsilon \sigma_n) \right] = 0
    $$
    
    Then the array is said to satisfy the Lindeberg condition.

This leads to the pinnacle of the asymptotic distribution theory for independent variables:

!!! success "Theorem 3.11: Lindeberg-Feller CLT"

    Assume $Var(X_{nj}) = \sigma_{nj}^2 < \infty$, $S_n = \sum_{j=1}^{k_n} X_{nj}$. Then the following two sets of propositions are equivalent:
    
    * **(i)** $\frac{S_n - E S_n}{\sigma_n} \rightarrow N(0,1)$ and **(ii)** the double array is a null array.
    * **$\Longleftrightarrow$** The double array satisfies the **Lindeberg Condition (LC)**.

    ??? proof "Core Proof Derivation of the Theorem (Click to expand)"

        **(1) Sufficiency proof: LC $\Rightarrow$ CLT and Null Array**
        
        Assume $E(X_{nj})=0$ and $\sigma_n^2 = 1$. Define a truncated random variable based on the truncation point $\eta \in (0, 1)$:
        $X_{nj}' = X_{nj}$ (if $|X_{nj}| < \eta$), otherwise 0.
        
        Calculate the truncated expectation and variance:
        
        $$
        |E[X_{nj}']| = \left| \int_{|x|<\eta} x dF_{nj}(x) \right| = \left| -\int_{|x|\ge\eta} x dF_{nj}(x) \right| \le \frac{1}{\eta} \int_{|x|\ge\eta} x^2 dF_{nj}(x)
        $$
        
        After summation, due to the LC condition, the total expectation tends to 0. Meanwhile, the truncated variance $\sigma_n'^2 \to 1 = \sigma_n^2$.
        According to Lemma 3.10 (diagonal construction rule), we can select a monotonically increasing sequence $m_n \to \infty$, and let $\eta_n = m_n^{-1} \to 0$.
        Using $\eta_n$ as the truncation threshold ensures that $|X_{nj}'| \le \eta_n := M_{nj}$. Since $\max M_{nj} = \eta_n \to 0$, according to the previous Corollary 3.6 (CLT for bounded variables), we have $(S_n' - ES_n')/\sigma_n' \to N(0,1)$, hence $S_n' \to N(0,1)$.
        
        Finally, evaluate the difference between the untruncated sum and the truncated sum:
        
        $$
        P(S_n \neq S_n') \le \sum_{j=1}^{k_n} P(|X_{nj}| \ge \eta_n) \le \sum_{j=1}^{k_n} \frac{1}{\eta_n^2} \int_{|x| > \eta_n} x^2 dF_{nj}(x) \rightarrow 0 \quad \text{(Guaranteed by LC)}
        $$
        
        By Slutsky's Theorem, since the error is $o_p(1)$, we finally get $S_n \to N(0,1)$.
        
        **(2) Necessity proof: CLT + Null Array $\Rightarrow$ LC**
        
        From $S_n \xrightarrow{d} N(0,1)$, we know the logarithm of the product of characteristic functions:
        
        $$
        \lim_{n \to \infty} \sum_{j=1}^{k_n} \log \phi_{nj}(t) = -\frac{t^2}{2}
        $$
        
        The null array guarantees $\max_j |\phi_{nj}(t) - 1| \to 0$. Using the equivalence relationship between $\log \phi_{nj}(t)$ and $\phi_{nj}(t) - 1$:
        
        $$
        \lim_{n \to \infty} \sum_{j=1}^{k_n} (\phi_{nj}(t) - 1) = -\frac{t^2}{2}
        $$
        
        Extracting its real part:
        
        $$
        \lim_{n \to \infty} \sum_j \int_{-\infty}^{\infty} (1 - \cos tx) dF_{nj}(x) = \frac{t^2}{2}
        $$
        
        Split the integral into $|x| \le \eta$ and $|x| > \eta$, and use the inequality $0 \le 1 - \cos tx \le (tx)^2/2$:
        
        $$
        \frac{t^2}{2} - \sum_j \int_{|x|\le\eta} (1 - \cos tx) dF_{nj}(x) \ge \frac{t^2}{2} \sum_j \int_{|x| \ge \eta} x^2 dF_{nj}(x)
        $$
        
        By bounding the outer part of the integral with Chebyshev's inequality ($\le 2/\eta^2 + \epsilon$), letting $t \to \infty$ to take the limit, we finally force out $\sum_{j=1}^{k_n} E[X_{nj}^2 \mathbb{I}(|X_{nj}| > \eta)] \to 0$, which is exactly the Lindeberg condition LC. $\square$

---

## 6. Applications & Further Conditions

!!! tip "Application Example: Ordinary Least Squares (OLS) Regression"

    Consider the classic linear regression model $y_j = x_j \beta + \epsilon_j$, where error terms $\epsilon_j \sim (0, \sigma_\epsilon^2)$ and are i.i.d..
    The design matrix satisfies $\max_{1 \le j \le n} \frac{|x_j|}{a_n} \to 0$, where $a_n^2 = \sum_{j=1}^n x_j^2$.
    The OLS estimator is $\hat{\beta}_{LS} = \sum_{j=1}^n x_j y_j / a_n^2$.
    
    Construct the standardized double array $X_{nj} = \frac{x_j \epsilon_j}{\sqrt{\sum x_j^2}}$.
    We bound the truncated integral in the Lindeberg condition: let $m_n = \max_j |x_j/a_n|$,
    
    $$
    \frac{1}{\sigma_n^2 a_n^2} \sum_{j=1}^n x_j^2 E\left[ \epsilon_j^2 \mathbb{I}(|x_j \epsilon_j / a_n| > \delta \sigma_n) \right] \le \frac{1}{\sigma_\epsilon^2} E\left[ \epsilon_j^2 \mathbb{I}(|\epsilon_j| > \delta \sigma_\epsilon m_n^{-1}) \right]
    $$
    
    Since $\epsilon_j$ are identically distributed with a finite second moment, as $m_n \to 0$, this expectation tends to 0. By Theorem 3.11 we immediately obtain:
    
    $$
    a_n(\hat{\beta}_{LS} - \beta) \xrightarrow{d} N(0, \sigma_\epsilon^2)
    $$

When verifying the LC condition, besides using bounded truncation, we can also utilize the existence of higher-order moments. This is the **extension of the Liapounov type condition**.

!!! info "Proposition 3.12: Sufficient Criterion for Lindeberg Condition"

    For a double array $\{X_{nj}\}$, if there exists some real number $\nu > 2$, satisfying:
    
    $$
    \sum_{j=1}^{k_n} E|X_{nj} - \mu_{nj}|^\nu = o(\sigma_n^\nu)
    $$
    
    Then, this array necessarily satisfies the Lindeberg condition.

    ??? proof "Proof Derivation (Click to expand)"

        In the truncated region $|t - \mu_{nj}| > \epsilon \sigma_n$, we bound its second moment integral:
        
        $$
        E\left[ (X_{nj} - \mu_{nj})^2 \mathbb{I}(|X_{nj} - \mu_{nj}| > \epsilon \sigma_n) \right] = \int_{|t - \mu_{nj}| > \epsilon \sigma_n} (t - \mu_{nj})^2 dF_{nj}(t)
        $$
        
        By forcing a $\nu$-th power and extracting the constant term:
        
        $$
        \le (\epsilon \sigma_n)^{2-\nu} \int_{|t - \mu_{nj}| > \epsilon \sigma_n} (t - \mu_{nj})^\nu dF_{nj}(t) \le (\epsilon \sigma_n)^{2-\nu} E|X_{nj} - \mu_{nj}|^\nu
        $$
        
        Summing and dividing by $\sigma_n^2$:
        
        $$
        \frac{1}{\sigma_n^2} \sum_{j=1}^{k_n} E\left[ \dots \right] \le \epsilon^{2-\nu} \frac{\sum_{j=1}^{k_n} E|X_{nj} - \mu_{nj}|^\nu}{\sigma_n^\nu} \rightarrow 0
        $$
        
        The Lindeberg condition is proven. $\square$