# Chapter 5: Weakly Dependent Data I

In classical limit theory, we usually assume that data is independent and identically distributed (i.i.d.). However, in practical applications (e.g., daily, weekly, or yearly sampled financial time series), data often exhibits serial correlation. This chapter explores stationary time series, various **Mixing Coefficients** that characterize data dependence, and a series of crucial covariance inequalities under weakly dependent conditions.

Let $Z_1, \dots, Z_n \in \mathbb{R}^d$ be equally sampled time series data, where $d$ is the dimension of the observation vector.

---

## 1. Stationarity and Time Series Models

To perform statistical inference on dependent data, we generally require the data to satisfy certain stationarity conditions.

!!! info "Definition of Stationarity"

    * **Strictly Stationary**: For any integers $l$ and $m$, the sequence of random vectors $(Z_{i_1}, \dots, Z_{i_m})^T$ and the shifted sequence $(Z_{i_1+l}, \dots, Z_{i_m+l})^T$ have the exact same joint distribution (i.e., **strong shift invariance**).
    * **Weak Stationary / Second Order Stationary**: The first and second moments of the process depend only on the time interval and do not change with time, that is:
      
    \[
    E(Z_i) = E(Z_{i+1}), \quad Var(Z_i) = Var(Z_{i+l}), \quad Cov(Z_i, Z_j) = Cov(Z_{i+l}, Z_{j+l})
    \]
      
      (i.e., **weak shift invariance**).

    *Note: In practice, common methods to transform non-stationary time series into stationary ones include differencing, taking the square root, or logarithmic transformations.*

### 1.1 Linear Time Series Models

!!! info "Definition 4.1: Autoregressive Moving Average Model ARMA(p, q)"

    A sequence $\{Z_t\}_{t \in \mathbb{Z}}$ is said to be an **ARMA(p, q)** model if $\{Z_t\}$ is weakly stationary and for any $t$ satisfies:
    
    \[
    Z_t - \theta_1 Z_{t-1} - \dots - \theta_p Z_{t-p} = \epsilon_t - \eta_1 \epsilon_{t-1} - \dots - \eta_q \epsilon_{t-q}
    \]

    Where $\{\epsilon_t\}$ is an independent white noise process, denoted as $WN(0, \sigma^2)$.

Let $\theta(z)$ and $\eta(z)$ be polynomials of degree $p$ and $q$ respectively, defined as:

\[
\theta(z) = 1 - \theta_1 z - \dots - \theta_p z^p, \quad \eta(z) = 1 - \eta_1 z - \dots - \eta_q z^q
\]

Introduce the backward shift operator $B$, such that $B^j Z_t = Z_{t-j}$.

!!! info "Definition 4.2: Causality"

    An ARMA(p, q) process $\{Z_t\}$ is said to be **causal** if there exists an absolutely summable sequence of coefficients $\{\psi_j\}_{j=0}^\infty$ (i.e., $\sum_{j=0}^\infty |\psi_j| < \infty$), such that for all $t \in \mathbb{Z}$:
    
    \[
    Z_t = \sum_{j=0}^\infty \psi_j \epsilon_{t-j}
    \]

    *Note: A causal ARMA(p, q) model is essentially an $MA(\infty)$ process.*

!!! success "Theorem 4.3 (Causality Condition for ARMA Processes)"

    Let $\{Z_t\}$ be an ARMA(p, q) process, represented as $\theta(B)Z_t = \eta(B)\epsilon_t$. Assuming the polynomials $\theta(z)$ and $\eta(z)$ have no common roots, then $\{Z_t\}$ is causal **if and only if**:
    
    \[
    \theta(z) \ne 0, \quad \forall z \in \mathbb{C}, \ |z| \le 1
    \]

    That is, the roots of the characteristic polynomial all lie strictly outside the unit circle. In this case, the coefficients $\{\psi_j\}$ are determined by the expansion $\psi(z) = \eta(z) / \theta(z)$.

If we allow the summation index to extend to negative infinity, we can define the generalized **Linear Process**:

\[
Z_t = \sum_{j=-\infty}^\infty \psi_j \epsilon_{t-j}, \quad \epsilon_t \sim i.i.d.\ F(0, \sigma^2)
\]

*(Under certain conditions, the ARMA process is a special case of the linear process.)*

### 1.2 ARCH(p) Model: Conditional Heteroskedasticity

The **Auto-Regression Conditionally Heteroskedasticity (ARCH)** model is a profoundly important non-linear model in financial time series:

\[
Z_t = m(\vec{Z}_{t,p}) + \sigma(\vec{Z}_{t,p})\epsilon_t
\]

Where $\vec{Z}_{t,p} = (Z_{t-1}, \dots, Z_{t-p})^T$.
It generalizes the classical AR(p) model ($Z_t = \theta_0 + \theta_1 Z_{t-1} + \dots + \theta_p Z_{t-p} + \epsilon_t$) in two aspects:
1. It generalizes the linear conditional mean to a **non-linear conditional mean function** $m(\cdot)$.
2. It generalizes the constant conditional variance to a **state-dependent conditional variance function** $\sigma(\cdot)$.

*Note: ARCH models and general non-linear time series models are not necessarily always stationary. However, under certain parameter constraints, they can be guaranteed to be "asymptotically stationary" — meaning they tend towards stationarity after a period of "pre-burning".*

---

## 2. Mixing Coefficients: Measures of Dependence

For weakly dependent sequences, the influence of past events on the future should gradually decay to zero as the time gap lengthens. We use Mixing Coefficients to precisely measure this decay rate.

Let $\{Z_1, \dots, Z_t, \dots\}$ be a strictly stationary process. For positive integers $l \le m$, define $\mathcal{F}_l^m$ as the $\sigma$-algebra generated by the random variables $\{Z_i\}_{i=l}^m$. For a time gap $k \ge 1$, we define the following mixing coefficients:

1. **$\alpha$-mixing coefficient (or strong mixing)**:
   
    \[
    \alpha(k) = \sup_{B \in \mathcal{F}_{-\infty}^t, C \in \mathcal{F}_{t+k}^\infty} |P(B \cap C) - P(B)P(C)|
    \]

2. **$\beta$-mixing coefficient (absolute regularity)**:
   
    \[
    \beta(k) = E \left[ \sup_{C \in \mathcal{F}_{t+k}^\infty} |P(C) - P(C | \mathcal{F}_{-\infty}^t)| \right]
    \]

3. **$\phi$-mixing coefficient (uniform mixing)**:
   
    \[
    \phi(k) = \sup_{B \in \mathcal{F}_{-\infty}^t, C \in \mathcal{F}_{t+k}^\infty} |P(C) - P(C | B)|
    \]

4. **$\rho$-mixing coefficient (maximal correlation)**:
   
    \[
    \rho(k) = \sup_{X \in L^2(\mathcal{F}_{-\infty}^t), Y \in L^2(\mathcal{F}_{t+k}^\infty)} |Corr(X, Y)| = \sup_{X, Y} \left| \frac{Cov(X, Y)}{\sqrt{Var(X)Var(Y)}} \right|
    \]
   
   *(Here, $L^2(\mathcal{F})$ denotes the set of all random variables that are measurable with respect to the $\sigma$-algebra $\mathcal{F}$ and have finite second moments, $EX^2 < \infty$.)*

### 2.1 Relationships and Properties of Mixing Coefficients

There exist strict inequality relationships among the various mixing coefficients:

\[
2\alpha(k) \le \beta(k) \le \phi(k), \quad 4\alpha(k) \le \rho(k) \le 2\phi^{1/2}(k)
\]

**Definition:**
If $\lim_{k \to \infty} \alpha(k) = 0$, the process $\{Z_t\}$ is said to be **$\alpha$-mixing**. We can similarly define $\phi$-mixing, $\rho$-mixing, etc. These mixing properties measure the dependence between the past ($\mathcal{F}_{-\infty}^t$) and the future ($\mathcal{F}_{t+k}^\infty$). As the time gap $k \to \infty$, the mixing coefficients approaching zero implies that the system exhibits **asymptotic independence**.

From the above inequalities, the following implications hold:

\[
\phi\text{-mixing} \implies \beta\text{-mixing}, \quad \rho\text{-mixing} \implies \alpha\text{-mixing}
\]

*Note: $\alpha$-mixing is the **weakest condition** among all the mixing conditions above (the easiest to satisfy), but ironically, it was historically named "Strong Mixing".*

### 2.2 When are Linear and Markov Processes Mixing?

* **Linear Processes**: For a causal linear process $Z_t = \sum \psi_j \epsilon_{t-j}$, Gorodetskii showed that under certain conditions, it possesses the $\alpha$-mixing property. Pham and Tran further proved that if the coefficients decay exponentially as $j \to \infty$, i.e., $\psi_j = O(r^j)$ (where $0 < r < 1$), then the process is **geometric $\alpha$-mixing**, meaning there exist constants $C$ and $\rho \in [0, 1)$ such that $\alpha(k) \le C \rho^k$.
* **Markov Processes / ARCH(p) Models**: For processes of the Markov form $Y_i = m(X_i) + \sigma(X_i)\epsilon_i$, Masry and Tjosheim provided specific conditions under which the process is ergodic and geometrically $\alpha$-mixing.

---

## 3. Core Inequalities: Covariance Bounds for Mixing Sequences

To study the limit theory of dependent variables (e.g., variance convergence in limit theorems), we need to control the covariance of random variables separated by $k$ time steps. The following four lemmas form the cornerstone of the theory for weakly dependent data.

### 3.1 Billingsley's Inequality (Covariance Bound for Bounded Variables)

!!! success "Lemma 4.4 (Billingsley's Inequality)"

    Suppose $\{Z_i\}$ is $\alpha$-mixing (**stationarity is NOT required**), and the random variables $X \in \mathcal{F}_{-\infty}^t$ and $Y \in \mathcal{F}_{t+k}^\infty$ are uniformly bounded, i.e., $|X| \le C_1$ and $|Y| \le C_2$. Then their covariance is bounded by:
    
    \[
    |Cov(X, Y)| \le 4C_1 C_2 \alpha(k)
    \]

??? proof "Proof of Billingsley's Inequality (Click to expand)"

    Considering the definition of covariance, we can expand it using the tower property of conditional expectation:

    \[
    |Cov(X, Y)| = |E(XY) - E(X)E(Y)| = |E[X \{ E(Y | \mathcal{F}_{-\infty}^t) - EY \}]|
    \]

    Since $|X| \le C_1$ is $\mathcal{F}_{-\infty}^t$-measurable, we have:

    \[
    |Cov(X, Y)| \le C_1 E|E(Y | \mathcal{F}_{-\infty}^t) - EY| = C_1 E[\xi (E(Y | \mathcal{F}_{-\infty}^t) - EY)]
    \]

    Here, we introduce the sign function $\xi = \text{sgn}(E(Y | \mathcal{F}_{-\infty}^t) - EY)$. Since the conditional expectation itself is $\mathcal{F}_{-\infty}^t$-measurable, the directional indicator variable $\xi$ is also $\mathcal{F}_{-\infty}^t$-measurable.

    Reverting this to the form of unconditional expectation:

    \[
    |Cov(X, Y)| \le C_1 |E(\xi Y) - E\xi EY| = C_1 |Cov(\xi, Y)|  \quad \dots (1)
    \]

    For $|Cov(\xi, Y)|$, we can apply the exact same trick to $Y$. Define $\eta = \text{sgn}(E(\xi | \mathcal{F}_{t+k}^\infty) - E\xi) \in \mathcal{F}_{t+k}^\infty$. Since $|Y| \le C_2$, we obtain:

    \[
    |Cov(\xi, Y)| \le C_2 |E(\xi \eta) - E\xi E\eta| \quad \dots (2)
    \]

    Now we need to bound $|E(\xi \eta) - E\xi E\eta|$. Note that $\xi, \eta$ can only take values in $\{1, -1\}$.
    Define the set $A = \{\xi = 1\}$, so $A^c = \{\xi = -1\}$; similarly, $B = \{\eta = 1\}$ and $B^c = \{\eta = -1\}$.
    Clearly, $A, A^c \in \mathcal{F}_{-\infty}^t$ and $B, B^c \in \mathcal{F}_{t+k}^\infty$.

    Expanding the expectation:

    \[
    \begin{aligned}
    |E(\xi \eta) - E\xi E\eta| &= |P(AB) + P(A^c B^c) - P(A^c B) - P(A B^c) - (P(A) - P(A^c))(P(B) - P(B^c))| \\
    &= |4 (P(A \cap B) - P(A)P(B))| \\
    &\le 4 \alpha(k) \quad \dots (3)
    \end{aligned}
    \]

    *(The second step is obtained by substituting $P(A^c) = 1 - P(A)$ and $P(B^c) = 1 - P(B)$ and simplifying. The final step directly applies the definition of $\alpha$-mixing.)*

    Combining equations (1), (2), and (3) completes the proof of Billingsley's Inequality:

    \[
    |Cov(X, Y)| \le 4C_1 C_2 \alpha(k)
    \]

    $\square$

---

### 3.2 Truncation Technique and $L^p$ Norm Bounds

!!! success "Lemma 4.5 (Case where $X$ has higher-order moments and $Y$ is bounded)"

    Suppose $\{Z_i\}$ is $\alpha$-mixing (stationarity is not required), $X \in \mathcal{F}_{-\infty}^t$ and there exists $p > 1$ such that $E|X|^p < \infty$; simultaneously, $Y \in \mathcal{F}_{t+k}^\infty$ and is bounded $|Y| \le C$. Then:
    
    \[
    |Cov(X, Y)| \le 6C \|X\|_p \alpha^{1/q}(k)
    \]

    Where $q$ is the conjugate index of $p$ (satisfying $\frac{1}{p} + \frac{1}{q} = 1$), and $\|X\|_p = (E|X|^p)^{1/p}$.

??? proof "Proof of Lemma 4.5: Truncation Method (Click to expand)"

    To utilize Lemma 4.4 (which requires bounded variables), we introduce a truncation constant $M > 0$ for the unbounded $X$.
    Let $X_M = X \mathbb{I}(|X| \le M)$ and the tail part $X'_M = X - X_M = X \mathbb{I}(|X| > M)$.
    By the bilinearity of covariance and the triangle inequality:

    \[
    |Cov(X, Y)| = |Cov(X_M, Y) + Cov(X'_M, Y)| \le |Cov(X_M, Y)| + |Cov(X'_M, Y)|
    \]

    For the bounded first term, we directly apply Lemma 4.4 (here $X_M$ is bounded by $M$, and $Y$ is bounded by $C$):

    \[
    |Cov(X_M, Y)| \le 4CM \alpha(k)
    \]

    For the second tail term, we first control its expectation. For $X'_M$, since the region of integration is $|X| > M$, the ratio $|x|/M > 1$:

    \[
    E|X'_M| = \int_{|x| > M} |x| dF(x) \le \int_{|x| > M} |x| \left(\frac{|x|}{M}\right)^{p-1} dF(x)
    \]

    Simplifying this gives the bound for the absolute expectation of the tail:

    \[
    E|X'_M| \le M^{-p+1} E|X|^p \quad \dots (4)
    \]

    Based on this, we can control the covariance of the tail:

    \[
    \begin{aligned}
    |Cov(X'_M, Y)| &= |E(X'_M Y) - EX'_M EY| \\
    &\le E|X'_M Y| + C E|X'_M| \\
    &\le 2C E|X'_M| \\
    &\le 2C M^{-p+1} E|X|^p \quad \dots (5)
    \end{aligned}
    \]

    Combining both terms, the upper bound for the total covariance is: $4CM\alpha(k) + 2C M^{-p+1} E|X|^p$.
    To minimize this upper bound, we cleverly choose the truncation point $M$ to be:

    \[
    M = \|X\|_p \{ \alpha(k) \}^{-1/p}
    \]

    Substituting this in yields:

    \[
    |Cov(X, Y)| \le 4C \|X\|_p \alpha(k)^{1 - 1/p} + 2C \|X\|_p^p \left(\|X\|_p \alpha(k)^{-1/p}\right)^{1-p}
    \]

    Since $1 - 1/p = 1/q$, the above expression simplifies to:

    \[
    |Cov(X, Y)| \le 4C \|X\|_p \alpha^{1/q}(k) + 2C \|X\|_p \alpha^{1/q}(k) = 6C \|X\|_p \alpha^{1/q}(k)
    \]

    The proof is complete. $\square$

---

### 3.3 Rio's Inequality and Quantile Functions

!!! success "Lemma 4.6 (Rio's Inequality)"

    Let $X$ and $Y$ be two integrable real-valued random variables satisfying the limit $\lim_{c \to \infty} E\{|X| \mathbb{I}(|X| > c)\} = 0$. Define $Q_X(u) = \inf\{t: P(|X| > t) \le u\}$ as the **upper quantile function** of $|X|$.
    
    If $Q_X Q_Y$ is integrable over the interval $(0, 1)$, then:
    
    \[
    |Cov(X, Y)| \le 2 \int_0^{2\alpha} Q_X(u) Q_Y(u) du
    \]
    
    Where $\alpha = \alpha(\sigma(X), \sigma(Y))$ is the $\alpha$-mixing coefficient between the $\sigma$-algebras generated by $X$ and $Y$.

*(The proof of this inequality involves deep knowledge of non-parametric statistics and stochastic processes, see Bosq, D. (1998) for details. It serves as a crucial foundation for deriving stronger inequalities subsequently.)*

---

### 3.4 Davydov Inequality (Generalized $L^q-L^r$ Covariance Bound)

Davydov's Inequality is the most commonly used covariance inequality. It removes the restriction that variables must be bounded, requiring only that they possess appropriate higher-order moments.

!!! success "Lemma 4.7 (Davydov Inequality)"

    Let $X$ and $Y$ be two real-valued random variables such that $X \in L^q(\mathcal{F}_{-\infty}^t)$ and $Y \in L^r(\mathcal{F}_{t+k}^\infty)$. Where $q > 1$, $r > 1$, and there exists a $p$ satisfying the Hölder-style relationship:
    
    \[
    \frac{1}{q} + \frac{1}{r} = 1 - \frac{1}{p} \quad (\text{where } p > 1)
    \]
    
    Then we have:
    
    \[
    |Cov(X, Y)| \le 2p (2\alpha(k))^{1/p} \|X\|_q \|Y\|_r
    \]

??? proof "Proof of Davydov Inequality (Click to expand)"

    The core of the proof is to utilize Rio's Inequality (Lemma 4.6) by discussing three cases based on the finiteness of the parameters.

    **(i) Assume both $q$ and $r$ are finite:**
    According to Markov's Inequality, for any $u \in (0, 1]$, we have:
    
    \[
    P\left( |X| > \frac{\|X\|_q}{u^{1/q}} \right) \le \frac{E|X|^q}{(\|X\|_q / u^{1/q})^q} = u
    \]
    
    Based on the definition of the upper quantile function, this directly implies that for all $0 < u \le 1$:
    
    \[
    Q_X(u) \le \frac{\|X\|_q}{u^{1/q}}
    \]
    
    Symmetrically, for $Y$ we also have $Q_Y(u) \le \frac{\|Y\|_r}{u^{1/r}}$.
    Substituting these two terms into the integral of Rio's Inequality:
    
    \[
    \begin{aligned}
    |Cov(X,Y)| &\le 2 \int_0^{2\alpha(k)} \frac{\|X\|_q}{u^{1/q}} \frac{\|Y\|_r}{u^{1/r}} du \\
    &= 2 \|X\|_q \|Y\|_r \int_0^{2\alpha(k)} u^{-\frac{1}{q} - \frac{1}{r}} du
    \end{aligned}
    \]
    
    Since $\frac{1}{q} + \frac{1}{r} = 1 - \frac{1}{p}$, the exponent part inside the integral is $u^{\frac{1}{p} - 1}$. Integrating this power function:
    
    \[
    \int_0^{2\alpha(k)} u^{\frac{1}{p} - 1} du = \left[ p \cdot u^{1/p} \right]_0^{2\alpha(k)} = p (2\alpha(k))^{1/p}
    \]
    
    Substituting this back into the original expression yields:
    
    \[
    |Cov(X,Y)| \le 2p (2\alpha(k))^{1/p} \|X\|_q \|Y\|_r
    \]

    **(ii) Assume $r = +\infty$, and $q$ is finite:**
    According to the parameter relationship, $\frac{1}{q} + \frac{1}{p} = 1$.
    For $Y$ in the $L^\infty$ space, its upper quantile function has a natural hard upper bound $Q_Y(u) \le Q_Y(0) = \|Y\|_\infty = \sup |Y|$.
    *(Because the smallest $t$ that makes $P(|Y| > t) = 0$ is exactly $\|Y\|_\infty$)*.
    
    Applying Rio's Inequality again:
    
    \[
    \begin{aligned}
    |Cov(X,Y)| &\le 2 \int_0^{2\alpha(k)} \frac{\|X\|_q}{u^{1/q}} \|Y\|_\infty du \\
    &= 2 \|X\|_q \|Y\|_\infty \int_0^{2\alpha(k)} u^{\frac{1}{p} - 1} du \\
    &= 2p (2\alpha(k))^{1/p} \|X\|_q \|Y\|_\infty
    \end{aligned}
    \]

    *(Note: This result is highly similar to the previously derived Lemma 4.5, but the constant factor differs slightly.)*

    **(iii) Assume $r = +\infty$ and $q = +\infty$:**
    In this case, both variables are essentially bounded. Substituting into the parameter relationship must yield $1/p = 1 \implies p = 1$.
    Here, we apply the hard bounds for both quantile functions: $Q_X(u) \le \|X\|_\infty$ and $Q_Y(u) \le \|Y\|_\infty$.
    Substituting into Rio's Inequality:
    
    \[
    |Cov(X,Y)| \le 2 \int_0^{2\alpha(k)} \|X\|_\infty \|Y\|_\infty du = 4 \|X\|_\infty \|Y\|_\infty \alpha(k)
    \]
    
    This completely degenerates exactly into the form of Billingsley's Inequality (Lemma 4.4).
    
    Summarizing all cases, the lemma is proven. $\square$