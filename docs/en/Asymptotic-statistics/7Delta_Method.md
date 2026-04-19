# Chapter 7: Delta Method

In previous chapters, we studied the convergence of sequences of random variables (such as the sample mean) themselves. This chapter will explore a more central question in statistics: **Is convergence in distribution preserved under smooth transformations?** This is the famous **Delta Method** and its wide applications in various statistical inferences (e.g., confidence intervals, hypothesis testing, variance stabilization).

Suppose we have a sequence of estimators $\{T_n\}_{n \ge 1}$ (taking values in $\mathbb{R}^k$) for a parameter $\theta \in \mathbb{R}^k$.

* For a parameter function of interest $\phi(\theta)$ (where $\phi: \mathbb{R}^k \rightarrow \mathbb{R}^m$), according to the Continuous Mapping Theorem, if $T_n \xrightarrow{p} \theta$ and $\phi$ is continuous at $\theta$, then $\phi(T_n) \xrightarrow{p} \phi(\theta)$.

* But in statistics, a more interesting and practical question is: if it is known that $\sqrt{n}(T_n - \theta) \Rightarrow T$, does the non-linearly transformed sequence $\sqrt{n}(\phi(T_n) - \phi(\theta))$ also converge to a definite distribution?

---

## 1. The Derivative of Vector-valued Functions and the First-Order Delta Method

Recalling multivariable calculus, if a function $\phi(\cdot)$ is differentiable at $\theta$, it means there exists a linear mapping (matrix) $\phi'_{\theta}: \mathbb{R}^k \mapsto \mathbb{R}^m$ such that:

\[
\phi(\theta + h) - \phi(\theta) = \phi'(\theta)h + R(h)
\]

where the remainder term satisfies $R(h) = o(\|h\|)$ as $h \rightarrow 0$.

The specific form of this derivative mapping (Jacobian Matrix) is:

\[
\phi'_{\theta} = \begin{pmatrix} 
\frac{\partial \phi_1}{\partial \theta_1}(\theta) & \cdots & \frac{\partial \phi_1}{\partial \theta_k}(\theta) \\ 
\vdots & \ddots & \vdots \\ 
\frac{\partial \phi_m}{\partial \theta_1}(\theta) & \cdots & \frac{\partial \phi_m}{\partial \theta_k}(\theta) 
\end{pmatrix}_{m \times k}
\]

*(Note: If $m=1, k>1$, this derivative mapping is exactly the **Gradient** of the function.)*

### 1.1 First Order Delta Method

!!! success "Theorem 5.1 (First Order Delta Method)"

    If $\phi$ is differentiable at $\theta$, and the derivative matrix $\phi'(\theta) \ne 0$. Suppose there is a deterministic divergent sequence $\{r_n\}$ (typically $r_n = \sqrt{n}$) satisfying $r_n \rightarrow \infty$, and $r_n(T_n - \theta) \Rightarrow T$, then:

    (i) $r_n(\phi(T_n) - \phi(\theta)) - \phi'(\theta)(r_n(T_n - \theta)) \xrightarrow{p} 0$

    (ii) $r_n(\phi(T_n) - \phi(\theta)) \Rightarrow \phi'(\theta)T$

??? proof "Proof of Theorem 5.1 (Click to expand)"

    **Proof of (i):**
    
    It is known that $r_n(T_n - \theta) \Rightarrow T$. Since $r_n \rightarrow \infty$, by Stochastic Boundedness ($O_p(1)$), we must have:

    \[
    T_n - \theta \xrightarrow{p} 0
    \]

    Using the Taylor expansion of $\phi$ at $\theta$ (differentiability), for a sufficiently small $h = T_n - \theta$:

    \[
    \phi(T_n) - \phi(\theta) - \phi'(\theta)(T_n - \theta) = o(\|T_n - \theta\|)
    \]

    Multiplying both sides by $r_n$:

    \[
    r_n [ \phi(T_n) - \phi(\theta) - \phi'(\theta)(T_n - \theta) ] = r_n \cdot o_p(\|T_n - \theta\|)
    \]

    Rewriting the right side as:

    \[
    o_p(1) \cdot r_n \|T_n - \theta\|
    \]

    Since $r_n(T_n - \theta) = O_p(1)$, we have $o_p(1) \cdot O_p(1) = o_p(1)$. This proves conclusion (i).

    **Proof of (ii):**
    
    Rearranging the conclusion of (i) yields:

    \[
    r_n(\phi(T_n) - \phi(\theta)) = \phi'(\theta) r_n(T_n - \theta) + o_p(1)
    \]

    Since $r_n(T_n - \theta) \Rightarrow T$, according to the Continuous Mapping Theorem, $\phi'(\theta) r_n(T_n - \theta) \Rightarrow \phi'(\theta)T$.
    Finally, applying **Slutsky's Theorem** (adding a term that converges in probability to 0 does not change the convergence in distribution), we directly obtain:

    \[
    r_n(\phi(T_n) - \phi(\theta)) \Rightarrow \phi'(\theta)T
    \]

    The proof is complete. $\square$

**Classic Application: Normal Delta Method**

If the estimator satisfies asymptotic normality: $\sqrt{n}(T_n - \theta) \xrightarrow{d} N(0, \sigma^2(\theta))$.
For any scalar function $g: \mathbb{R} \rightarrow \mathbb{R}$ that is differentiable at $\theta$ with derivative $g'(\theta) \ne 0$, we have:

\[
\sqrt{n}[g(T_n) - g(\theta)] \xrightarrow{d} N(0, [g'(\theta)]^2 \sigma^2(\theta))
\]

---

## 2. High Order Delta Method

The First Order Delta Method relies heavily on $\phi'(\theta) \ne 0$. If we encounter a degenerate case where $\phi'(\theta) = 0$ but $\phi''(\theta) \ne 0$, the first-order method fails (resulting in a degenerate point mass distribution). In this case, we need to introduce a higher-order Taylor expansion.

Expanding to the second-order term:

\[
\phi(T_n) = \phi(\theta) + \frac{1}{2}\phi''(\theta)(T_n - \theta)^2 + \cdots
\]

Multiplying by $n$ (note it is $n$ here instead of $\sqrt{n}$, due to the presence of the squared term):

\[
n(\phi(T_n) - \phi(\theta)) = \frac{1}{2}\phi''(\theta)[\sqrt{n}(T_n - \theta)]^2 \Rightarrow \frac{1}{2}\phi''(\theta)T^2
\]

!!! success "Theorem 5.2 (High Order Delta Method)"

    Suppose a univariate function $\phi$ is $m$ times differentiable at $\theta$, and satisfies $\phi^{(m)}(\theta) \ne 0$ but all preceding lower-order derivatives are zero (i.e., $\phi^{(j)}(\theta) = 0, \forall j < m$). If $r_n(T_n - \theta) \Rightarrow T$, then:

    \[
    \frac{r_n^m (\phi(T_n) - \phi(\theta))}{\frac{1}{m!} \phi^{(m)}(\theta)} \Rightarrow T^m
    \]

### 2.1 High Order Delta Method Application Example

Suppose $X_1, \dots, X_n$ is an i.i.d. sequence with mean $\mu$ and known variance $\sigma^2$. We want to test the null hypothesis $H_0: \mu = 0$.
Under the null hypothesis, the statistic $n\bar{X}_n^2 / \sigma^2 \rightarrow [N(0,1)]^2 = \chi_1^2$.

Now consider the limiting behavior of the random variable $\cos(\bar{X}_n)$:

* **If we force the use of the First Order Delta Method**: Since the derivative of the function $g(x) = \cos(x)$ at $x=0$ is $g'(0) = -\sin(0) = 0$, the normalization term $\sqrt{n}$ will lead to:
  
    \[
    \sqrt{n}(\cos(\bar{X}_n) - 1) \xrightarrow{p} 0
    \]
  
    This provides no useful distributional information, indicating that $\sqrt{n}$ is not the correct convergence rate.

* **Using the Second Order Delta Method**: Because at $x=0$, the second derivative $\cos''(0) = -\cos(0) = -1 \ne 0$. Expanding yields:
  
    \[
    \cos(\bar{X}_n) - \cos(0) = (\bar{X}_n - 0) \cdot 0 + \frac{1}{2}(\bar{X}_n - 0)^2 \cdot (-1) + o_p(\bar{X}_n^2)
    \]
  
    Multiplying by $-2n$: $-2n(\cos(\bar{X}_n) - 1) = n\bar{X}_n^2 + o_p(1) \Rightarrow \chi_1^2 \cdot \sigma^2$
  
    This gives the correct non-degenerate limiting distribution.

---

## 3. Asymptotic Normality and Classic Applications of the Delta Method

### 3.1 Limit Distribution of Sample Variance and Standard Deviation

Let $X_1, \dots, X_n \sim i.i.d. F$, with finite 4th moments. Let the population central moments be $\alpha_i = E(X_1^i)$, and sample moments be $m_{ni} = n^{-1}\sum X_j^i$.
The sample variance can be written as a function of two sample moments:

\[
S_n = n^{-1}\sum_{i=1}^n (X_i - \bar{X})^2 = m_{n2} - m_{n1}^2 = \phi(m_{n1}, m_{n2})
\]

where the non-linear transformation function is $\phi(x_1, x_2) = x_2 - x_1^2$. Its gradient vector is:

\[
\phi'(\alpha_1, \alpha_2) = (-2\alpha_1, 1)
\]

By the Multivariate Central Limit Theorem (Multivariate CLT):

\[
\sqrt{n} \left[ \begin{pmatrix} m_{n1} \\ m_{n2} \end{pmatrix} - \begin{pmatrix} \alpha_1 \\ \alpha_2 \end{pmatrix} \right] \xrightarrow{d} N\left( 0, Var \begin{pmatrix} X_1 \\ X_1^2 \end{pmatrix} \right)
\]

Applying the multivariate first-order Delta method, the limiting distribution of the sample variance is:

\[
\sqrt{n}(S_n - \sigma^2) \xrightarrow{d} N\left( 0, (-2\alpha_1, 1) Var \begin{pmatrix} X_1 \\ X_1^2 \end{pmatrix} \begin{pmatrix} -2\alpha_1 \\ 1 \end{pmatrix} \right)
\]

By expanding the quadratic form, it can be cleverly simplified into the form of central moments: $E(X_1 - \alpha_1)^4 - [E(X_1 - \alpha_1)^2]^2 = c_4 - c_2^2$ (i.e., the fourth central moment minus the square of the variance).
Therefore:

\[
\sqrt{n}(S_n - \sigma^2) \xrightarrow{d} N(0, c_4 - c_2^2)
\]

**Corollary: Unbiased Sample Variance and Sample Standard Deviation**

For the unbiased variance $S_{n-1} = \frac{n}{n-1}S_n$, since the constant factor converges to 1 in the limit, and the difference term $\sqrt{n}(\frac{n}{n-1} - 1)S_n = o_p(1)$, it has the **same limiting distribution**.

For the **sample standard deviation** $S_n^{1/2}$, applying the univariate Delta method, we take $\phi(x) = \sqrt{x}$, and its derivative is $\phi'(x) = \frac{1}{2}x^{-1/2}$. Substituting the value at $\sigma^2$:

\[
\sqrt{n}(S_n^{1/2} - \sigma) \xrightarrow{d} N\left( 0, \frac{c_4 - c_2^2}{4\sigma^2} \right)
\]

### 3.2 More Examples of Common Transformations

Suppose the base sequence $X_n$ satisfies Asymptotic Normality ($AN$): $X_n \sim AN(\mu, \sigma_n^2)$ and $\sigma_n \rightarrow 0$.

* (i) $X_n^2 \sim AN(\mu^2, 4\mu^2 \sigma_n^2)$ (requires $\mu \ne 0$)

* (ii) $\frac{1}{X_n} \sim AN(\mu^{-1}, \frac{\sigma_n^2}{\mu^4})$ (requires $\mu \ne 0$)

* (iii) $e^{X_n} \sim AN(e^\mu, e^{2\mu} \sigma_n^2)$ (for any $\mu$)

* (iv) $\log|X_n| \sim AN(\log|\mu|, \mu^{-2} \sigma_n^2)$ (requires $\mu \ne 0$). If $\mu = 0$ and $\sigma_n = 1/\sqrt{n}$, then by the Continuous Mapping Theorem, the limiting distribution is related to $\log|N(0,1)|$.

**Weighted $\chi^2$ Distribution for Multidimensional Quadratic Forms:**

Let $X_1, \dots, X_n \sim i.i.d.\ F$ in the space $\mathbb{R}^p$, with mean $\mu$ and covariance $\Sigma$. Consider the target $\hat{\theta} = \bar{X}^T \bar{X}$.

* **If $\mu \ne 0$**: Applying the first-order Delta method, $\phi'(\mu) = 2\mu^T$, we have $\sqrt{n}(\bar{X}^T\bar{X} - \mu^T\mu) \xrightarrow{d} N(0, 4\mu^T \Sigma \mu)$.

* **If $\mu = 0$**: The first derivative is 0 (because $\mu^T \Sigma \mu = 0$). In this case, higher-order mapping is needed. Since $\sqrt{n}\bar{X} \xrightarrow{d} N_p(0, \Sigma)$, we have:
  
    \[
    n\bar{X}^T\bar{X} \Rightarrow N_p^T(0, \Sigma) N_p(0, \Sigma) \stackrel{d}{=} Z^T \Sigma^{1/2} \Sigma^{1/2} Z = Z^T \Sigma Z
    \]
  
    where $Z \sim N_p(0, I_p)$. Through eigenvalue decomposition $\Sigma = U^T \text{diag}(\lambda_1, \dots, \lambda_p) U$, the above expression is equivalent to the linear combination: $\sum_{i=1}^p \lambda_i \chi_{1i}^2$
  
    This is a **Weighted $\chi^2$ distribution**.

---

## 4. Asymptotic Theory in Hypothesis Testing

### 4.1 $\chi^2$ Test for Variance and the Effect of Excessive Kurtosis

Let $X_1, \dots, X_n \sim i.i.d.\ F$, with $EX_1^4 < \infty$. We want to test $H_0: \sigma^2 \le 1$ VS $H_1: \sigma^2 > 1$.
Under the normality assumption, the test statistic is $nS_n$, and the rejection region is $nS_n > \chi^2_{n-1, \alpha}$. The size of the test is exactly $\alpha$.

However, if the data distribution $F$ is **not a normal distribution**, and there is an **Excessive Kurtosis** $\kappa = \frac{E(X-\mu)^4}{\sigma^4} - 3 \ne 0$, the situation changes fundamentally.

It is known that for a chi-square distribution formed by the sum of standard normal variables, when $n$ is large:

\[
\frac{\chi^2_{n-1} - (n-1)}{\sqrt{2(n-1)}} \xrightarrow{d} N(0, 1)
\]

And the asymptotic distribution of the true sample variance (known from Section 3.1) is:

\[
\sqrt{n}\left( \frac{S_n}{\sigma^2} - 1 \right) \xrightarrow{d} N(0, \kappa + 2) \ne N(0, 2)
\]

**The Actual Size of the Test (Type I Error Rate):**

Using the asymptotic expansion of the chi-square critical value $\chi^2_{n-1, \alpha} \approx (n-1) + Z_\alpha \sqrt{2(n-1)}$, when the true variance lies on the boundary $\sigma^2 = 1$:

\[
P_{\sigma^2=1}(nS_n > \chi^2_{n-1, \alpha}) \approx P\left( \sqrt{n}(S_n - 1) > \frac{Z_\alpha \sqrt{2n}}{\sqrt{n}} \right) \rightarrow P(N(0, \kappa+2) > \sqrt{2}Z_\alpha)
\]

After standardization:

\[
= 1 - \Phi\left( \frac{\sqrt{2} Z_\alpha}{\sqrt{\kappa + 2}} \right)
\]

* **Conclusion**: For a distribution with heavy tails ($\kappa > 0$), the true Size will be strictly greater than the nominal $\alpha$. This explains why the traditional chi-square test for variance produces too many false positives under non-normal data.

### 4.2 Multinomial Vectors and the Pearson $\chi^2$ Statistic

Consider a multinomial distribution $(n_1, \dots, n_K) \sim Multinomial(n; p_1, \dots, p_K)$.
Define the standardized frequency vector $X_n = \sqrt{n}(\frac{n_1}{n} - p_1, \dots, \frac{n_K}{n} - p_K)^T \xrightarrow{d} N(0, \Sigma)$.
where the elements of the covariance matrix $\Sigma$ are $\sigma_{ii} = p_i(1-p_i)$ and $\sigma_{ij} = -p_i p_j$.

The Pearson $\chi^2$ statistic for Goodness-of-fit can be written as a quadratic form:

\[
T_n = \sum_{i=1}^K \frac{(n_i - np_i)^2}{np_i} = X_n^T C X_n
\]

where $C = \text{diag}(p_1^{-1}, \dots, p_K^{-1})$.

By the mapping theorem, the limiting distribution is the quadratic form $Z^T \Sigma^{1/2} C \Sigma^{1/2} Z$.
It can be proven that the matrix $A = \Sigma^{1/2} C \Sigma^{1/2}$ is an **Idempotent matrix ($A^2 = A$)**.
A quadratic form with an idempotent matrix follows a chi-square distribution, with degrees of freedom equal to its Trace:

\[
\text{tr}(\Sigma^{1/2} C \Sigma^{1/2}) = \text{tr}(C \Sigma) = \sum_{i=1}^K p_i^{-1} p_i(1-p_i) = K - 1
\]

Therefore, the Pearson statistic $T_n \Rightarrow \chi^2_{K-1}$.

### 4.3 Wald Test

For a multi-dimensional hypothesis test $H_0: \mu = \mu_0$ VS $H_1: \mu \ne \mu_0$, the commonly used Wald statistic is:

\[
W_n = n(\bar{X} - \mu_0)^T S_n^{-1} (\bar{X} - \mu_0)
\]

By the Law of Large Numbers, the sample covariance matrix $S_n \xrightarrow{p} \Sigma$, thus $S_n^{-1} \xrightarrow{p} \Sigma^{-1}$.
Through the Plug-in method and asymptotic expansion:

\[
W_n = \sqrt{n}(\bar{X} - \mu_0)^T \Sigma^{-1} \sqrt{n}(\bar{X} - \mu_0) + \sqrt{n}(\bar{X} - \mu_0)^T (S_n^{-1} - \Sigma^{-1}) \sqrt{n}(\bar{X} - \mu_0)
\]

Since $\sqrt{n}(\bar{X} - \mu_0) = O_p(1)$ and $S_n^{-1} - \Sigma^{-1} = o_p(1)$, the second term is $o_p(1)$.
The first term is exactly a standard multivariate normal quadratic form, therefore:

\[
W_n \xrightarrow{d} \chi^2_p
\]

---

## 5. Variance Stabilizing Transform (VST)

When constructing confidence intervals using asymptotic normality:

\[
T_n \pm Z_{1-\alpha/2} \frac{\sigma(\hat{\theta})}{\sqrt{n}}
\]

We find that the width of the interval fluctuates dramatically with the change of the unknown parameter $\theta$ (reflected in $\sigma(\theta)$).
The purpose of the **Variance Stabilizing Transform (VST)** is to find a smooth transformation $\phi(\cdot)$ such that the limiting variance of the transformed statistic no longer depends on the parameter $\theta$:

\[
\sqrt{n}(\phi(T_n) - \phi(\theta)) \xrightarrow{d} N(0, c^2)
\]

where $c > 0$ is a constant.

From the first-order Delta method, it is known that the transformed variance is $(\phi'(\theta))^2 \sigma^2(\theta)$. Setting this equal to a constant $c^2$:

\[
\phi'(\theta) \sigma(\theta) = c \implies \phi'(\theta) = \frac{c}{\sigma(\theta)}
\]

Integrating both sides, we obtain the core construction formula of VST:

\[
\phi(\theta) = \int \frac{d\theta}{\sigma(\theta)}
\]

### 5.1 VST Application: Tukey's Hanging Rootgram

In Nonparametric Kernel Density Estimation (KDE):

\[
\hat{f}_{nh}(x) = \frac{1}{nh} \sum_{i=1}^n K\left( \frac{x - X_i}{h} \right)
\]

It is known that under appropriate bandwidth conditions:

\[
\sqrt{nh}(\hat{f}_{nh}(x) - f(x)) \Rightarrow N(0, f(x))
\]

That is, the variance of the original estimator is proportional to the density function itself $f(x)$. To stabilize the variance for drawing uniform error bands, we apply VST. Here the variance term is $\sigma^2(f) = f$. Substituting into the construction formula:

\[
\phi(f) = \int \frac{df}{\sqrt{f}} = f^{1/2}
\]

(Ignoring constants of integration and multiples). Therefore, we take the square root of the density estimator ("Root-gram"):

\[
\hat{f}_{nh}^{1/2}(x) \sim AN\left( f^{1/2}(x), \frac{1}{4nh} \right)
\]

At this point, the asymptotic variance only relates to the sample size and bandwidth, perfectly eliminating the dependence on the density value $f(x)$.

---

## 6. Uniform Integrability and Asymptotic Approximation of Moments

The Delta method can be used not only to study the convergence of distributions but also to approximate the expectation and variance of estimators. But this requires a bridge connecting "convergence in distribution" and "convergence of moments" — **Uniform Integrability**.

!!! info "Definition 5.3 (Asymptotically Uniformly Integrable, u.i.)"

    A sequence $\{Y_n\}_{n \ge 0}$ is called asymptotically uniformly integrable if it satisfies:

    \[
    \lim_{M \rightarrow \infty} \limsup_{n \rightarrow \infty} E[|Y_n| \mathbb{I}_{\{|Y_n| > M\}}] = 0
    \]

Uniform integrability is the key to ensuring the validity of taking limits under the expectation.

!!! success "Theorem 5.4"

    Let $f: \mathbb{R}^k \rightarrow \mathbb{R}$ be measurable and continuous at every point in a set $C$. If $X_n \xrightarrow{d} X$ and $X$ takes its values in $C$. Then:
    
    \[
    E[f(X_n)] \rightarrow E[f(X)] \quad \text{if and only if the sequence } f(X_n) \text{ is asymptotically u.i.}
    \]

**Taylor Approximation of Moments (Moment Approximation)**

If we want to use a second-order Taylor expansion to approximate $E[\phi(T_n)]$ and $Var(\phi(T_n))$:

\[
\phi(T_n) = \phi(\theta) + \phi'(\theta)(T_n - \theta) + \frac{1}{2}\phi''(\theta)(T_n - \theta)^2 + \cdots
\]

After taking the expectation and variance, we *expect* to obtain:

* $E[\phi(T_n)] \approx \phi(\theta) + \phi'(\theta)\text{Bias}(T_n) + \frac{1}{2}\phi''(\theta)\text{MSE}(T_n)$

* $Var(\phi(T_n)) \approx [\phi'(\theta)]^T Var(T_n) [\phi'(\theta)]$

**The Strict Prerequisite for Validity**: For the above approximate equalities to strictly hold, we must ensure that the expectation of the remainder term converges. This requires the random sequence $\phi(T_n) - \phi(\theta)$ to be **uniformly integrable (u.i.)**. Generally, if the base deviation $T_n - \theta$ is uniformly integrable, and the function $\phi$ satisfies the Lipschitz continuity condition, then the transformed sequence is also uniformly integrable, thus ensuring the asymptotic validity of the moment approximation.