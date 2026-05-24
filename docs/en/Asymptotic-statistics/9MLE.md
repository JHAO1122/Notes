# Chapter 9: Maximum Likelihood Estimation (II)

In the previous chapter, we explored the consistency and asymptotic normality of the maximum likelihood estimator (MLE). In this chapter, we will start from the famous Wilks' theorem and delve into the asymptotic properties of the likelihood ratio test. Next, we will discuss how to use the Bartlett correction to improve the accuracy of the likelihood ratio test, and introduce counterexamples when parameters are not identifiable. Finally, we will transcend the limitations of parametric models and introduce the core tool of modern nonparametric statistics—**Empirical Likelihood**—and demonstrate its elegant mathematical structure and the nonparametric version of Wilks' theorem.

---

## 1. Wilks' Theorem for Likelihood Ratio

Consider the composite hypothesis testing problem $H_0: \theta \in \Theta_0$ VS $H_1: \theta \in \Theta_1$. We need to clarify the mathematical structure of the null hypothesis space $\Theta_0$.

### 1.1 Structural Assumption for $\Theta_0$ (Condition S1)

Assume there exists a set $A \subset \mathbb{R}^r$ with interior points (i.e., an $r$-dimensional set), and $p$ functions $g_i$ defined on $\mathbb{R}^r$, such that $\theta_i = g_i(\phi_1, ..., \phi_r)$ for $i=1,...,p$.
This establishes a $1-1$ correspondence between $\Theta_0$ and $A$. Moreover, each $g_i$ is twice continuously differentiable at any interior point of $A$. The true parameter $\theta_0 \in \Theta_0$ under $H_0$ corresponds to an interior point $\phi_0 \in A$.

This condition implies that the family of distributions under the null hypothesis can be equivalently parameterized as:

\[
\{ \tilde{f}(x, \phi) = f(x, g_1(\phi_1, ..., \phi_r), ..., g_p(\phi_1, ..., \phi_r)) : \phi = (\phi_1, ..., \phi_r) \in A \}
\]

**Equivalent statement:**
$\Theta_0$ can also be equivalently defined by $p-r$ constraints:

\[
R_i(\theta) = 0, \quad 1 \le i \le p-r
\]

* **Example of simple hypothesis**: For $H_0: \theta = \theta_0$, we have $\Theta_0 = \{\theta_0\}$, and the constraint functions can be taken as $R_i(\theta) = \theta_i - \theta_{i0}$ ($1 \le i \le p$).
* **Example of composite hypothesis**: If $p=3$ and $H_0: \theta_1 = \theta_{10}$ (only the first parameter is restricted), then the constraint is $R_1(\theta) = \theta_1 - \theta_{10}$.

### 1.2 Wilks' Theorem and Three Test Statistics

Let $\hat{\theta}_n$ be the global maximum likelihood estimator, $\hat{\phi}_n$ be the restricted MLE under the family $\{f(x, \phi), \phi \in A\}$, and the corresponding restricted parameter estimate is $\tilde{\theta}_n = g(\hat{\phi}_n)$.

Define the log-likelihood ratio as:

\[
LR_n = -2 \log \frac{L_n(\tilde{\theta}_n)}{L_n(\hat{\theta}_n)}
\]

!!! success "Theorem 9.1 (Wilks' Theorem, Wilks 1938)"

    Under Condition S1 and the regularity conditions for the asymptotic normality of the global MLE, under the null hypothesis $H_0$ we have:

    \[
    LR_n \xrightarrow{d} \chi^2_{p-r}
    \]

Under the null hypothesis, the following three major test statistics are **asymptotically equivalent**:

1.  **Likelihood Ratio Statistic (LR statistic)**: $LR_n$

2.  **Wald Test Statistic (based on Delta method)**:

    \[
    \mathcal{W}_n = R(\hat{\theta}_n)^T \left( \frac{\partial}{\partial \theta} R(\hat{\theta}_n) \left( -\frac{\partial^2 \log L_n(\hat{\theta}_n)}{\partial \theta \partial \theta^T} \right)^{-1} \frac{\partial}{\partial \theta} R(\hat{\theta}_n)^T \right)^{-1} R(\hat{\theta}_n)
    \]

3.  **Rao Score Test Statistic (Lagrange Multiplier method)**:

    \[
    R_n = \frac{\partial \log L_n(\tilde{\theta}_n)}{\partial \theta^T} \left( -\frac{\partial^2 \log L_n(\tilde{\theta}_n)}{\partial \theta \partial \theta^T} \right)^{-1} \frac{\partial \log L_n(\tilde{\theta}_n)}{\partial \theta}
    \]


??? proof "Proof Sketch of Wilks' Theorem (click to expand)"

    Under $H_0$, via a second-order Taylor expansion of the log-likelihood function, we obtain:

    \[
    2(\log L_n(\hat{\theta}_n) - \log L_n(\theta_0)) = \frac{1}{n} \left(\frac{\partial \log L_n(\theta_0)}{\partial \theta}\right)^T I_1(\theta_0)^{-1} \left(\frac{\partial \log L_n(\theta_0)}{\partial \theta}\right) + o_p(1)
    \]

    Similarly, using the parameterization $\phi$ under $H_0$, the expansion of the restricted log-likelihood is:

    \[
    2(\log L_n(g(\hat{\phi}_n)) - \log L_n(g(\phi_0))) = \frac{1}{n} \left(\frac{\partial \log L_n(g(\phi_0))}{\partial \phi}\right)^T \tilde{I}_1(\phi_0)^{-1} \left(\frac{\partial \log L_n(g(\phi_0))}{\partial \phi}\right) + o_p(1)
    \]

    Using the chain rule:

    \[
    \frac{\partial \log L_n(g(\phi_0))}{\partial \phi} = \frac{\partial g(\phi_0)^T}{\partial \phi} \frac{\partial \log L_n(\theta_0)}{\partial \theta}
    \]

    Combining the above results, the likelihood ratio statistic $LR_n$ can be written as:

    \[
    LR_n = \frac{1}{n} \left( \frac{\partial \log L_n(\theta_0)}{\partial \theta} \right)^T \left( I_1(\theta_0)^{-1} - \frac{\partial g(\phi_0)}{\partial \phi^T} \tilde{I}_1(\phi_0)^{-1} \frac{\partial g(\phi_0)^T}{\partial \phi} \right) \left( \frac{\partial \log L_n(\theta_0)}{\partial \theta} \right) + o_p(1)
    \]

    By the Central Limit Theorem (CLT), the score function satisfies:

    \[
    \frac{1}{\sqrt{n}} I_1(\theta_0)^{-1/2} \frac{\partial \log L_n(\theta_0)}{\partial \theta} \xrightarrow{d} Z \sim N(0, I_{p \times p})
    \]

    Denote $A = I_1(\theta_0)$, $C = \tilde{I}_1(\phi_0)$, $D = \frac{\partial g(\phi_0)^T}{\partial \phi}$, and $B = A^{-1} - D^T C^{-1} D$. It can be verified that $C = D A D^T$.
    Therefore $LR_n \Rightarrow Z^T A^{1/2} B A^{1/2} Z$.

    We need to show that the matrix $A^{1/2} B A^{1/2}$ is **idempotent**:

    \[
    (A^{1/2} B A^{1/2})^2 = A^{1/2} B A B A^{1/2}
    \]
    \[
    = A^{1/2} (A^{-1} - D^T C^{-1} D) A (A^{-1} - D^T C^{-1} D) A^{1/2}
    \]
    \[
    = (I_p - A^{1/2} D^T C^{-1} D A^{1/2})(I_p - A^{1/2} D^T C^{-1} D A^{1/2})
    \]
    \[
    = I_p - 2 A^{1/2} D^T C^{-1} D A^{1/2} + A^{1/2} D^T C^{-1} D A D^T C^{-1} D A^{1/2}
    \]

    Since $D A D^T = C$, the last term simplifies to $A^{1/2} D^T C^{-1} C C^{-1} D A^{1/2} = A^{1/2} D^T C^{-1} D A^{1/2}$. Substituting back yields:

    \[
    = I_p - A^{1/2} D^T C^{-1} D A^{1/2} = A^{1/2} B A^{1/2}
    \]

    Since this matrix is idempotent, its quadratic form follows a chi-squared distribution with degrees of freedom equal to its trace, i.e., $Tr(A^{1/2} B A^{1/2}) = p - r$. Therefore $LR_n \xrightarrow{d} \chi^2_{p-r}$.

---

## 2. Nuisance Parameters, Confidence Regions, and Counterexamples

### 2.1 Nuisance Parameter Case

Suppose $\theta = (\beta, \gamma)$, where $\beta$ is the parameter of interest, and $\gamma$ is a **nuisance parameter**.
Test $H_0: \beta = \beta_0$ VS $H_1: \beta \ne \beta_0$.

The likelihood ratio statistic is:

\[
LR_n(\beta_0) = 2 \{ l_n(\hat{\beta}_n, \hat{\gamma}_n) - l_n(\beta_0, \tilde{\gamma}_n) \}
\]

where $\hat{\theta}_n = (\hat{\beta}_n, \hat{\gamma}_n)$ is the unrestricted global MLE, and $\tilde{\gamma}_n = \arg\sup_{\gamma} l_n(\beta_0, \gamma)$ is the restricted MLE.
Assume the parameter space satisfies $\Theta = \Theta_\beta \oplus \Theta_\gamma$, and $dim(\Theta_\gamma) = r$. By Wilks' theorem, we have:

\[
LR_n(\beta_0) \xrightarrow{d} \chi^2_{p-r}
\]

The "local linearity" structure of the hypothesis test is crucial for the chi-squared approximation.

### 2.2 Likelihood Ratio Confidence Regions

Using the likelihood ratio statistic, we can construct a $(1-\alpha)$-level confidence region for $\beta$:

\[
I_{1-\alpha}(\beta) = \{ \beta | LR_n(\beta) \le \chi^2_{p-r, 1-\alpha} \}
\]

This confidence region has a **natural shape and direction**, completely determined by the data and model automatically, without subjective choice of boundaries (unlike the symmetric ellipsoids constructed using Wald tests). Although "plotting" this confidence region may require computational or graphical tools, modern computing power has made it quite accessible.

### 2.3 Counterexample Where Wilks' Theorem Fails

Wilks' theorem relies on strict regularity conditions. If these are not satisfied, anomalous behavior occurs.

**Counterexample**: Let $X_1, \dots, X_n \stackrel{i.i.d}{\sim} \text{Uniform}(0, \theta)$, $\theta > 0$.
Here, the support of the density function depends on the parameter $\theta$ (violating a basic regularity condition of the theorem).

We want to test $H_0: \theta = \theta_0$ VS $H_1: \theta \ne \theta_0$.
The unrestricted MLE is the maximum order statistic $\hat{\theta}_n = X_{(n)}$.
The likelihood ratio statistic is:

\[
LR_n = -2 \log \frac{\theta_0^{-n}}{X_{(n)}^{-n}} = -2n \log \frac{X_{(n)}}{\theta_0}
\]

It can be shown that the asymptotic distribution of $LR_n$ is some deformation of $\chi^2_2$ or related to an exponential distribution, **not** the standard $\chi^2_1$.

---

## 3. Bartlett Correction for Likelihood Ratio

Bartlett (1937) discovered a correction method that improves the accuracy of the asymptotic approximation of the likelihood ratio test.

Under $H_0$, $LR_n \rightarrow \chi^2_{p-r}$. Suppose the mean of $LR_n$ has the following asymptotic expansion:

\[
E(LR_n) = (p-r) + \frac{b_1}{n} + o(n^{-1}) = (p-r)\left\{ 1 + \frac{b}{n} + o(n^{-1}) \right\}
\]

Here $b$ or $(1 + b/n)$ is called the **Bartlett factor**. Let $\hat{b}_n$ be a consistent estimator of $b$ ($\hat{b}_n \rightarrow b$). By Slutsky's theorem, the corrected statistic:

\[
LR_{n, bc} := \frac{LR_n}{1 + b/n} \xrightarrow{d} \chi^2_{p-r}
\]
\[
\hat{LR}_{n, bc} := \frac{LR_n}{1 + \hat{b}/n} \xrightarrow{d} \chi^2_{p-r}
\]

**Bartlett's Conjecture and Edgeworth Expansion**:
Wilks' theorem only guarantees $P(LR_n \le x) = P(\chi^2_{p-r} \le x) + o(1)$.
If the data tails are sufficiently smooth and additional moment conditions (Cramér's condition) hold, via Edgeworth expansion it can be shown (Barndorff-Nielsen and Cox, 1984):

\[
P(LR_n(\theta_0) \le x) = P(\chi^2_{p-r} \le x) - \beta x g_p(x) n^{-1} + O(n^{-2})
\]

where $g_p$ is the probability density function of $\chi^2_{p-r}$.
After the Bartlett correction, the convergence rate improves substantially:

\[
P\left(\frac{LR_n}{1 + \beta n^{-1}} \le x\right) = P(\chi^2_{p-r} \le x) + O(n^{-2})
\]
\[
P\left(\frac{LR_n}{1 + \hat{\beta} n^{-1}} \le x\right) = P(\chi^2_{p-r} \le x) + O(n^{-3/2})
\]

**Significance of the Bartlett Correction**:
Using the corrected statistic to construct a confidence region $I_{\alpha, bc}$, the coverage error $P(\theta \in I_{\alpha, bc}) - (1-\alpha) = O(n^{-3/2})$, which has a smaller higher-order error compared to the uncorrected $O(n^{-1})$. It provides a highly accurate approximate significance level for the likelihood ratio test.

---

## 4. From Parametric to Nonparametric: Empirical Likelihood

In standard parametric MLE, we assume the data come from some known distribution family $F_\theta$. However, if the model is misspecified, the MLE can be severely biased and inconsistent.

**Question**: Can we construct a likelihood function without assuming a specific parametric distribution family?
**Intuition**: Let the data "speak for themselves". The empirical cumulative distribution function (ECDF):

\[
F_n(x) = \frac{1}{n} \sum_{i=1}^n I(X_i \le x)
\]

is a natural, consistent nonparametric estimator of the true CDF $F(x)$. Can we define a likelihood function that is maximized by $F_n$?

!!! info "Definition 9.1 (Nonparametric Likelihood of a Distribution)"

    Let $X_1, ..., X_n$ be i.i.d. samples from a completely unknown distribution $F$. Given the data, the nonparametric likelihood of a distribution $F$ is defined as:

    \[
    L(F) = \prod_{i=1}^n dF(X_i)
    \]

    where $dF(X_i) = P_F(X = X_i)$ is the probability mass assigned to the observed point $X_i$ by the distribution $F$.

To maximize $L(F)$, the distribution $F$ must place all of its probability mass strictly on the observed data points.
Let $p_i = dF(X_i)$. We need to satisfy $p_i \ge 0$ and $\sum_{i=1}^n p_i = 1$. The likelihood function simplifies to $L(p) = \prod_{i=1}^n p_i$.

!!! success "Theorem 9.2 (ECDF Maximizes Nonparametric Likelihood)"

    The empirical distribution $F_n(x)$ maximizes the nonparametric likelihood $L(F)$ among all possible distributions.

??? proof "Proof: Using the Arithmetic-Geometric Mean Inequality (AM-GM)"

    We aim to maximize $\prod_{i=1}^n p_i$ subject to $p_i \ge 0$ and $\sum_{i=1}^n p_i = 1$.
    By the AM-GM inequality:

    \[
    \left( \prod_{i=1}^n p_i \right)^{1/n} \le \frac{1}{n} \sum_{i=1}^n p_i = \frac{1}{n}
    \]

    Equality holds if and only if $p_1 = p_2 = \dots = p_n = \frac{1}{n}$.
    This means that, without any constraints, the nonparametric likelihood is maximized when each observation is assigned mass $1/n$, which exactly recovers the empirical distribution $F_n$. The maximum likelihood value is $(1/n)^n$.

!!! info "Definition 9.2 (Empirical Likelihood Ratio)"

    Introduced by Art Owen (1988, 1990), the empirical likelihood ratio of a candidate distribution $F$ is defined as the ratio of its likelihood to the unrestricted maximum likelihood:

    \[
    R(F) = \frac{L(F)}{L(F_n)} = \frac{\prod_{i=1}^n p_i}{\prod_{i=1}^n (1/n)} = \prod_{i=1}^n (n p_i)
    \]

    Note: For all valid distributions $F$, $R(F) \le 1$. This ratio measures the "plausibility" of the candidate distribution $F$ relative to the unrestricted empirical distribution $F_n$.

---

## 5. Empirical Likelihood for the Mean and Optimization

Suppose we are interested in a parameter vector, for example, the mean $\mu = E_F(X) = \int x dF(x)$.
The **profile empirical likelihood ratio** $\mathcal{R}(\mu)$ for $\mu$ is defined as the maximum of $R(F)$ over all distributions $F$ satisfying the parameter constraint:

\[
\mathcal{R}(\mu) = \sup \left\{ \prod_{i=1}^n (n p_i) \Bigg| p_i \ge 0, \sum_{i=1}^n p_i = 1, \sum_{i=1}^n p_i X_i = \mu \right\}
\]

This transforms the statistical inference into a **constrained optimization problem**. We need to find a set of weights $p_i$ that are closest in likelihood sense to $1/n$, while forcing their expectation to equal the candidate $\mu$.

### Solving via Lagrange Multipliers

We aim to maximize the log-likelihood ratio:

\[
\log \mathcal{R}(\mu) = \sum_{i=1}^n \log(n p_i)
\]

subject to $\sum p_i = 1$ and $\sum p_i(X_i - \mu) = 0$. Introduce Lagrange multipliers $\gamma$ and $n\lambda$ (the scaling by $n$ for $\lambda$ is purely for algebraic convenience later), and construct the Lagrangian:

\[
\mathcal{L}(p, \gamma, \lambda) = \sum_{i=1}^n \log(n p_i) - \gamma \left( \sum_{i=1}^n p_i - 1 \right) - n\lambda^T \left( \sum_{i=1}^n p_i(X_i - \mu) \right)
\]

Take the partial derivative with respect to $p_i$ and set it to zero:

\[
\frac{\partial \mathcal{L}}{\partial p_i} = \frac{1}{p_i} - \gamma - n\lambda^T (X_i - \mu) = 0
\]

Multiply the entire equation by $p_i$ and sum over all $i$:

\[
\sum_{i=1}^n 1 - \gamma \sum_{i=1}^n p_i - n\lambda^T \sum_{i=1}^n p_i(X_i - \mu) = 0
\]

Using the constraints $\sum p_i = 1$ and $\sum p_i(X_i - \mu) = 0$, the equation elegantly collapses to:

\[
n - \gamma(1) - n\lambda^T(0) = 0 \implies \gamma = n
\]

Substituting $\gamma = n$ back into the first-order condition and solving for $p_i$, we obtain the optimal probability weights:

\[
p_i = \frac{1}{n} \frac{1}{1 + \lambda^T(X_i - \mu)}
\]

The Lagrange multiplier $\lambda$ is implicitly determined by the remaining moment constraint equation. Substituting $p_i$ yields the **score equation** for $\lambda$:

\[
\frac{1}{n} \sum_{i=1}^n \frac{X_i - \mu}{1 + \lambda^T(X_i - \mu)} = 0
\]

Substituting the optimal weights $p_i$ back into the objective function, we obtain the **Profile Empirical Log-Likelihood Ratio Statistic**:

\[
l_E(\mu) = -2 \log \mathcal{R}(\mu) = -2 \sum_{i=1}^n \log(n p_i) = -2 \sum_{i=1}^n \log \left( \frac{1}{1 + \lambda^T(X_i - \mu)} \right) = 2 \sum_{i=1}^n \log \{ 1 + \lambda^T(X_i - \mu) \}
\]

Amazingly, this function behaves very similarly to the traditional parametric likelihood ratio statistic in its asymptotic behavior.

---

## 6. Wilks' Theorem for Empirical Likelihood

!!! success "Theorem 9.3 (Wilks' Theorem for Empirical Likelihood, Owen 1988, 1990)"

    Let $X_1, ..., X_n$ be i.i.d. random vectors in $\mathbb{R}^d$ with true mean $\mu_0$ and a finite and positive definite covariance matrix $\Sigma$.
    Then, as $n \rightarrow \infty$:

    \[
    l_E(\mu_0) = -2 \log \mathcal{R}(\mu_0) \xrightarrow{d} \chi^2_d
    \]

**The extraordinary significance of this theorem**: We obtain a nonparametric pivotal quantity. To construct confidence intervals for $\mu$, we only need to rely on the quantiles of $\chi^2$, without assuming any parametric distribution, **and without explicitly estimating the covariance matrix $\Sigma$ at all**!

??? proof "Proof of Wilks' Theorem for Empirical Likelihood (for one-dimensional $d=1$, click to expand)"

    For algebraic simplicity, assume dimension $d=1$. Let $Z_i = X_i - \mu_0$. By definition, $E(Z_i) = 0$ and $E(Z_i^2) = \sigma^2$.

    The equation determining $\lambda$ is $g(\lambda) = \frac{1}{n} \sum_{i=1}^n \frac{Z_i}{1 + \lambda Z_i} = 0$. Applying the geometric series expansion $\frac{1}{1+x} = 1 - x + \frac{x^2}{1+x}$:

    \[
    0 = \frac{1}{n} \sum_{i=1}^n Z_i \left( 1 - \lambda Z_i + \frac{\lambda^2 Z_i^2}{1 + \lambda Z_i} \right) = \bar{Z} - \lambda \left( \frac{1}{n} \sum_{i=1}^n Z_i^2 \right) + \frac{1}{n} \sum_{i=1}^n \frac{\lambda^2 Z_i^3}{1 + \lambda Z_i}
    \]

    Since the variance is finite, we have $\max_i |Z_i| = o_p(n^{1/2})$, and the error term can be controlled. Solving for $\lambda$ gives its asymptotic representation:

    \[
    \lambda = \frac{\bar{Z}}{S_Z^2} + o_p(n^{-1/2}) = O_p(n^{-1/2})
    \]

    where $S_Z^2 = \frac{1}{n}\sum Z_i^2 \xrightarrow{p} \sigma^2$.

    Next, we use a Taylor series expansion of the empirical log-likelihood ratio statistic $l_E(\mu_0)$: $\log(1+x) = x - x^2/2 + x^3/3 - \dots$
    Since $\lambda = O_p(n^{-1/2})$ and $\max |Z_i| = o_p(n^{1/2})$, we have $\max |\lambda Z_i| \rightarrow 0$, allowing a safe Taylor expansion:

    \[
    l_E(\mu_0) = 2 \sum_{i=1}^n \log(1 + \lambda Z_i) = 2 \sum_{i=1}^n \left( \lambda Z_i - \frac{1}{2}\lambda^2 Z_i^2 + \eta_i \right)
    \]

    The remainder $2 \sum \eta_i \le C |\lambda|^3 \sum |Z_i|^3 = o_p(1)$. Therefore:

    \[
    l_E(\mu_0) = 2\lambda \sum_{i=1}^n Z_i - \lambda^2 \sum_{i=1}^n Z_i^2 + o_p(1)
    \]

    Substituting the leading term $\lambda \approx \frac{\bar{Z}}{\sigma^2}$ and $\sum Z_i^2 \approx n\sigma^2$ into the expansion:

    \[
    l_E(\mu_0) \approx 2 \left(\frac{\bar{Z}}{\sigma^2}\right)(n\bar{Z}) - \left(\frac{\bar{Z}}{\sigma^2}\right)^2 (n\sigma^2) = \frac{2n\bar{Z}^2}{\sigma^2} - \frac{n\bar{Z}^2}{\sigma^2} = \frac{n\bar{Z}^2}{\sigma^2}
    \]

    Through algebra we have successfully simplified it to:

    \[
    l_E(\mu_0) \approx \left( \frac{\sqrt{n}\bar{Z}}{\sigma} \right)^2 = \left( \frac{\sqrt{n}(\bar{X} - \mu_0)}{\sigma} \right)^2
    \]

    By the classical Central Limit Theorem (CLT), it is known that $\frac{\sqrt{n}(\bar{X} - \mu_0)}{\sigma} \xrightarrow{d} N(0,1)$.
    By the continuous mapping theorem, the square of a standard normal random variable converges to a chi-squared distribution with 1 degree of freedom:

    \[
    -2 \log \mathcal{R}(\mu_0) \xrightarrow{d} \chi^2_1
    \]

### Key Advantages of Empirical Likelihood

1.  **Data-Driven Shape**: Unlike Wald-type intervals (e.g., $\bar{X} \pm 1.96 \frac{\sigma}{\sqrt{n}}$), the confidence regions of empirical likelihood are entirely determined by the data. They need not be symmetric and can perfectly capture the skewness of the underlying distribution.

2.  **Range Preserving**: The EL interval is strictly contained within the convex hull of the data. The estimates never logically violate natural boundaries (e.g., negative probability values will not occur).

3.  **Internal Studentization**: The scaling (studentization) of the variance is implicitly accomplished within empirical likelihood. We do not need to explicitly plug in an estimator of the covariance $\Sigma$.

4.  **Transformation Respecting**: Empirical likelihood regions are invariant under smooth reparameterizations of the parameter.