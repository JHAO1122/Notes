# Chapter 8: Maximum Likelihood Estimation (I)

Maximum Likelihood Estimator (MLE) is one of the most central parameter estimation methods in modern statistics. In this chapter, we will systematically explore the asymptotic theory of maximum likelihood estimation, including Fisher information, Cramér-Rao lower bound, consistency of MLE, and its asymptotic normality. Finally, we will introduce likelihood-based hypothesis testing methods, particularly the Likelihood Ratio Test.

Assume $X = \{X_1, ..., X_n\}$ are independent and identically distributed (i.i.d.) samples, whose distribution $F_\theta$ belongs to a parametric family $\mathcal{F} = \{F_\theta : \theta = (\theta_1, ..., \theta_k)^T \in \Theta\}$, and suppose the distribution $F_\theta$ has probability density function $f_\theta(x)$. The **likelihood function** of the sample $X$ is defined as:

\[
L(\theta; X) = \prod_{i=1}^n f_\theta(X_i)
\]

The **maximum likelihood estimator (MLE)** $\hat{\theta}$ is the point that maximizes the likelihood function in the parameter space:

\[
\hat{\theta} = \arg \max_{\theta \in \Theta} \log L(\theta; X)
\]

Typically, the MLE can be obtained by solving the **score equations**:

\[
\frac{\partial \log L(\theta; X)}{\partial \theta_j} \Bigg|_{\theta=\hat{\theta}} = 0, \quad j=1, 2, \dots, k
\]

The variance of the score function plays a decisive role in the asymptotic normality of the MLE.

---

## 1. Fisher Information

!!! info "Definition 8.1 (Fisher Information Regularity FI Regularity)"

    Suppose the family of distributions $\mathcal{P} = \{P_\theta, \theta \in \Theta\}$ is dominated by a $\sigma$-finite measure $\mu$. If there exists an open neighborhood $\Theta_\theta$ of $\theta$ such that the following conditions hold, then $\mathcal{P}$ is said to be **FI regular** at $\theta \in \Theta$:

    (i) $f_\theta(x) := \frac{dP_\theta(x)}{d\mu} > 0$ for all $x$ and $\theta \in \Theta_\theta$.

    (ii) For every $x$, $f_\theta(x)$ is differentiable at $\theta$.

    (iii) The integral $\int f_\theta(x) \mu(dx)$ can be differentiated under the integral sign with respect to $\theta$, i.e.:

    \[
    \int \frac{d}{d\theta'} f_{\theta'}(x) \Big|_{\theta'=\theta} \mu(dx) = 0
    \]

!!! info "Definition 8.2 (Fisher Information)"

    If the model $\mathcal{P} = \{P_\theta, \theta \in \Theta\}$ is FI regular, then:

    \[
    I_1(\theta) = E_\theta \left[ \left( \frac{d}{d\theta'} \log f_{\theta'}(X) \Big|_{\theta'=\theta} \right)^2 \right]
    \]

    is called the **Fisher Information** based on a single observation $X$ at $\theta$.

**Important properties of Fisher information:**

1.  **Variance representation**: From the property of differentiating under the integral sign in the definition, the expectation of the score function is zero, i.e., $E_\theta[\frac{d}{d\theta'} \log f_{\theta'}(X)|_{\theta'=\theta}] = 0$. Thus Fisher information is the variance of the score function:

    \[
    I_1(\theta) = var\left( \frac{d}{d\theta'} \log f_{\theta'}(X) \Big|_{\theta'=\theta} \right)
    \]

2.  **Multivariate case**: If the parameter space $\Theta \subset \mathbb{R}^K$ ($K > 1$), where $\theta = (\theta_1, \dots, \theta_K)^T$, then the score function is a $K$-dimensional vector $\nabla_\theta \log f_{\theta'}(x)$. In this case, Fisher information is a matrix (FI matrix):

    \[
    I(\theta) = E_\theta \left[ \left( \nabla_\theta \log f_{\theta'}(X) \right) \left( \nabla_\theta \log f_{\theta'}(X) \right)^T \Big|_{\theta'=\theta} \right]
    \]

3.  **Second derivative form**: If $\mathcal{P}$ is FI regular, and for every $x$, $f_\theta(x)$ is **twice differentiable** at $\theta$, and the identity $1 = \int f_\theta(x)\mu(dx)$ can be differentiated twice under the integral sign (i.e., $\int \frac{d^2}{d\theta'^2} f_{\theta'}(x)|_{\theta'=\theta} \mu(dx) = 0$), then:

    \[
    I(\theta) = -E_\theta \left[ \frac{d^2}{d\theta'^2} \log f_{\theta'}(X) \Big|_{\theta'=\theta} \right]
    \]

---

## 2. Cramér-Rao Lower Bound and Bhattacharya Inequality

### 2.1 Cramér-Rao (C-R) Lower Bound

!!! success "Theorem 8.1 (Cramér-Rao Lower Bound)"

    Let $(\mathbb{X}, \mathcal{X}, \mathcal{P}=\{P_\theta, \theta \in \Theta\})$ be the probability space of a random variable $X$, $\mathcal{P}$ be absolutely continuous with respect to a $\sigma$-finite measure $\mu$, with density $f_\theta(x) = \frac{dP_\theta}{d\mu}$. Assume the following conditions hold:
  
    (i) $\Theta \subset \mathbb{R}$ is an open set.
    (ii) $A = \text{support}(f_\theta)$ does not depend on $\theta$.
    (iii) $\forall \theta \in \Theta$, the partial derivative $\frac{df_\theta(x)}{d\theta}$ exists.
    (iv) For any $\theta \in \Theta$, $E_\theta[\frac{\partial}{\partial\theta}\log f_\theta(X)] = \int \frac{\partial f_\theta(x)}{\partial\theta}\mu(dx) = 0$.
    (v) $I(\theta) > 0$ for all $\theta \in \Theta$.
    (vi) The target function $g: \Theta \rightarrow \mathbb{R}$ has derivative $\frac{dg(\theta)}{d\theta}$. And there exists an estimator $\hat{g}: \mathbb{X} \rightarrow \Theta$ such that $\hat{g}(X)$ is an **unbiased estimator** of $g(\theta)$.
    (vii) Integration and differentiation are interchangeable: $\frac{d}{d\theta}\int \hat{g}(x)f_\theta(x)\mu(dx) = \int \hat{g}(x)\frac{df_\theta(x)}{d\theta}\mu(dx)$.

    Then, the variance of the unbiased estimator $\hat{g}(X)$ satisfies the Cramér-Rao lower bound:

    \[
    var_\theta(\hat{g}(X)) \ge \frac{[g'(\theta)]^2}{I_n(\theta)}
    \]

    For the multivariate case, we have $var_\theta(\hat{g}(X)) \ge [g'(\theta)]^T I_n^{-1}(\theta) [g'(\theta)]$.

In the above theorem, conditions (iv) and (vii) (regarding the interchange of differentiation and integration) are the most restrictive. The following lemma provides a set of sufficient conditions (application of the dominated convergence theorem) that guarantee their validity.

!!! info "Lemma 8.1 (Sufficient Conditions for Interchangeability)"

    In addition to conditions (i)-(iii) of the previous theorem, assume there exists an envelope function $G: \mathbb{X} \times \Theta \rightarrow \mathbb{R}$ such that:
  
    (a) For every $\theta \in \Theta$, $G(x, \theta)$ is $\mathcal{X}$-measurable.
    (b) For every $\theta \in \Theta$, $E_\theta G^2(X, \theta) < \infty$.
    (c) For every $\theta \in \Theta$, there exists $\epsilon_\theta > 0$ such that when $|\theta - \theta'| < \epsilon_\theta$ and $x \in A$ we have:
  
    \[
    \left| \frac{df_{\theta'}(x)}{d\theta'} \right| \le G(x, \theta)f_\theta(x)
    \]

    Then condition (iv) must hold; and for any unbiased estimator $\hat{g}(X)$ with finite variance ($E_\theta(\hat{g}(X))^2 < \infty$), condition (vii) also holds.

??? proof "Proof of Lemma 8.1 (click to expand)"

    **Proof uses the Mean Value Theorem (MVT) and Dominated Convergence Theorem (DCT):**

    For any $\theta \in \Theta$ and $\theta' \in \Theta$ satisfying $|\theta - \theta'| < \epsilon_\theta$:
    Since $\int_{\mathcal{X}} f_\theta(x)\mu(dx) = \int_{\mathcal{X}} f_{\theta'}(x)\mu(dx) = 1$, we have:

    \[
    \int_{\mathcal{X}} \frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'} \mu(dx) = 0
    \]

    By the Mean Value Theorem (MVT) and condition (c), there exists $\tilde{\theta}$ between $\theta$ and $\theta'$ such that:

    \[
    \left| \frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'} \right| = \left| \frac{df_{\tilde{\theta}}(x)}{d\tilde{\theta}} \right| \le G(x, \theta)f_\theta(x)
    \]

    Note that the dominating function is $\mu$-integrable (using Cauchy-Schwarz):

    \[
    \int_{\mathcal{X}} G(x, \theta)f_\theta(x)\mu(dx) = E_\theta G(X, \theta) \le E_\theta^{1/2}[G^2(X, \theta)] < \infty
    \]

    Therefore, using the Dominated Convergence Theorem (DCT) as $\theta' \rightarrow \theta$:

    \[
    \int_{\mathcal{X}} \frac{df_\theta(x)}{d\theta}\mu(dx) = \lim_{\theta' \rightarrow \theta} \int_{\mathcal{X}} \frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'} \mu(dx) = 0
    \]

    This proves condition (iv).

    On the other hand, let $\hat{g}(X)$ be an unbiased estimator of $g(\theta)$ with $E_\theta \hat{g}^2(X) < \infty$. Consider the increment:

    \[
    \int_{\mathcal{X}} \hat{g}(x)\frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'}\mu(dx) = \frac{g(\theta) - g(\theta')}{\theta - \theta'}
    \]

    Using the same envelope bound:

    \[
    \left| \hat{g}(x)\frac{f_\theta(x) - f_{\theta'}(x)}{\theta - \theta'} \right| \le |\hat{g}(x)|G(x, \theta)f_\theta(x)
    \]

    Verify the integrability of the new dominating function:

    \[
    \int_{\mathcal{X}} |\hat{g}(x)|G(x, \theta)f_\theta(x)\mu(dx) = E_\theta[|\hat{g}(X)|G(X, \theta)] \le [E_\theta \hat{g}^2(X) \cdot E_\theta G^2(X, \theta)]^{1/2} < \infty
    \]

    Applying DCT again, as $\theta' \rightarrow \theta$, we obtain condition (vii). The proof is complete.

### 2.2 Bhattacharya Inequality

The C-R lower bound is sometimes "too low", insufficient to provide a tighter variance lower bound. The Bhattacharya inequality, by introducing higher-order derivatives, is a natural generalization of the C-R lower bound.

!!! success "Theorem 8.2 (Bhattacharya Inequality)"

    In addition to conditions (i) and (ii) of the C-R lower bound theorem, if we strengthen other conditions to:
  
    (iii)* For $i=1, \dots, K$, the $i$-th derivative of $f_\theta(x)$ exists and $\int_{\mathcal{X}} \frac{\partial^i f_\theta(x)}{\partial \theta^i}\mu(dx) = 0$.
    (iv)* For $i=1, \dots, K$, the variance of the higher-order logarithmic derivatives is finite: $\int_{\mathcal{X}} \frac{1}{f_\theta(x)} \left( \frac{\partial^i f_\theta(x)}{\partial \theta^i} \right)^2 \mu(dx) < \infty$.
    (v)* $\hat{g}(X)$ is an unbiased estimator of $g(\theta)$ with finite variance, and for every $i=1, \dots, K$, the interchange of integration and differentiation holds:
  
    \[
    g^{(i)}(\theta) = \frac{\partial^i}{\partial\theta^i}g(\theta) = \int_{\mathcal{X}} \hat{g}(x)\frac{\partial^i f_\theta(x)}{\partial\theta^i}\mu(dx)
    \]

    Then the variance of the estimator satisfies:

    \[
    var_\theta(\hat{g}(X)) \ge \tilde{g}^T(\theta) V^{-1}(\theta) \tilde{g}(\theta)
    \]

    where $V(\theta) = (V_{ij}(\theta))$ is a matrix with elements $V_{ij}(\theta) = E_\theta[\frac{1}{f^2_\theta(X)} \frac{\partial^i f_\theta(X)}{\partial\theta^i} \frac{\partial^j f_\theta(X)}{\partial\theta^j}]$, and the derivative vector $\tilde{g}(\theta) = (g'(\theta), \dots, g^{(K)}(\theta))^T$.

??? proof "Proof of Bhattacharya Inequality (click to expand)"

    Define the vector $S = S_\theta(x) = (S_\theta^{(1)}(x), \dots, S_\theta^{(K)}(x))^T$, where the $i$-th component is:

    \[
    S_\theta^{(i)}(x) = \frac{1}{f_\theta(x)} \frac{\partial^i f_\theta(x)}{\partial\theta^i}
    \]

    * From condition (iii)*, $E_\theta[S] = 0$.
    * From condition (iv)* and the definition of $V(\theta)$, $var_\theta(S) = V(\theta)$.
    * From condition (v)*, the covariance $cov_\theta(\hat{g}(X), S_\theta^{(i)}(X)) = E_\theta[\hat{g}(X) S_\theta^{(i)}(X)] - 0 = g^{(i)}(\theta)$.

    Consider the block covariance matrix $A$ of the combined vector $(\hat{g}, S^T)^T$:

    \[
    A := var_\theta \begin{pmatrix} \hat{g} \\ S \end{pmatrix} = \begin{pmatrix} var_\theta(\hat{g}(X)) & \tilde{g}^T(\theta) \\ \tilde{g}(\theta) & V(\theta) \end{pmatrix}
    \]

    Since the covariance matrix $A$ is positive semidefinite, its determinant $|A| \ge 0$. Using the formula for the determinant of a block matrix:

    \[
    |A| = |V(\theta)| \left[ var_\theta(\hat{g}(X)) - \tilde{g}^T(\theta) V^{-1}(\theta) \tilde{g}(\theta) \right]
    \]

    Since $|V(\theta)| > 0$, we must have:

    \[
    var_\theta(\hat{g}(X)) - \tilde{g}^T(\theta) V^{-1}(\theta) \tilde{g}(\theta) \ge 0
    \]

    This completes the proof. Clearly when $K=1$, this reduces to the Cramér-Rao lower bound.

---

## 3. Kullback-Leibler Divergence and Identifiability

To prove the consistency of the maximum likelihood estimator, we need a tool to measure the "distance" between two probability distributions, and to define the identifiability of parameters.

!!! info "Definition 8.3 (Kullback-Leibler Divergence)"

    The Kullback-Leibler (KL) divergence from probability measure $P_\theta$ to $P_\eta$ is defined as:

    \[
    D_{KL}(P_\eta || P_\theta) = E_\eta \left[ \log \frac{p_\eta(X)}{p_\theta(X)} \right], \quad X \sim P_\eta
    \]

    where $p_\theta$ and $p_\eta$ are the density functions of $P_\theta$ and $P_\eta$, respectively.

* KL divergence **is not a true metric**, because in general it is asymmetric: $D_{KL}(P || Q) \ne D_{KL}(Q || P)$.
* By the concavity of the logarithm (Jensen's inequality), we always have $D_{KL}(P || Q) \ge 0$, with equality if and only if $P = Q$ (under model identifiability).

!!! info "Definition 8.4 (Identifiability)"

    If for any $\theta_1 \ne \theta_2$ ($\theta_1, \theta_2 \in \Theta$), we have:

    \[
    \mu(x : f_{\theta_1}(x) \ne f_{\theta_2}(x)) > 0
    \]

    (i.e., the probability distributions given by two different parameters are not almost everywhere equal under the base measure $\mu$), then the parametric family $\mathbb{P}_\Theta := \{f_\theta(x) : \theta \in \Theta\}$ is said to be **identifiable**.

Identifiability is a necessary prerequisite for the consistency of MLE: if the parameter is not identifiable, a consistent estimator cannot exist at all.

!!! success "Lemma 8.2 (Minimizing KL Distance)"

    Let $\mathbb{P}_\Theta$ be an identifiable parametric family. If $E_{\theta_0} \log f_{\theta_0}(X) < \infty$, then the objective function $M(\theta) := E_{\theta_0} [\log f_\theta(X) / f_{\theta_0}(X)]$ attains its maximum **uniquely** at the true parameter $\theta_0$, i.e.:

    \[
    E_{\theta_0} \log f_\theta(X) \le E_{\theta_0} \log f_{\theta_0}(X) < \infty
    \]

    **Proof sketch:** For $\theta \in \Theta$, since $-\log(t)$ is strictly convex, applying Jensen's inequality:
  
    \[
    E_{\theta_0} \log \frac{f_\theta}{f_{\theta_0}}(X) \le \log E_{\theta_0} \left[ \frac{f_\theta}{f_{\theta_0}}(X) \right] = \log \int \frac{f_\theta(x)}{f_{\theta_0}(x)} f_{\theta_0}(x) dx = \log \int f_\theta(x) dx = \log 1 = 0
    \]
  
    By the identifiability condition, equality holds if and only if $\theta = \theta_0$. Therefore, the expected log-likelihood is maximized at the true parameter.

---

## 4. Consistency of MLE

!!! success "Theorem 8.3 (Existence of a Consistent Root of the Likelihood Equation)"

    Let $X_1, \dots, X_n$ be i.i.d. $\sim P_\theta$ and satisfy the Cramér-Rao conditions. Suppose there exists an open neighborhood $\Theta_\theta$ of $\theta$ such that:
  
    (i) The support $A := \{x | f_\theta(x) > 0\}$ does not depend on $\theta$.
    (ii) For every $x \in A$, $f_{\theta}(x)$ is differentiable with respect to any $\theta' \in \Theta_\theta$.
    (iii) The expectation $E_\theta \log f_{\theta'}(X)$ exists and is finite for all $\theta' \in \Theta_\theta$.
    (iv) The parametric family is identifiable.
  
    Then, for any $\epsilon > 0, \delta > 0$, there exists $m_{\epsilon, \delta} > 0$ such that when $n > m_{\epsilon, \delta}$:

    \[
    P_\theta\left\{ \text{the likelihood equation } \frac{d}{d\theta'} \sum_{i=1}^n \log f_{\theta'}(X_i) = 0 \text{ has a root in } (\theta-\epsilon, \theta+\epsilon) \right\} \ge 1 - \delta
    \]

??? proof "Proof of Theorem 8.3 (click to expand)"

    Without loss of generality, assume $\epsilon$ is small enough so that $[\theta-\epsilon, \theta+\epsilon] \subset \Theta_\theta$. By the Weak Law of Large Numbers (WLLN) and condition (iii):

    \[
    \frac{1}{n}\sum_{i=1}^n \log \frac{f_{\theta\pm\epsilon}(X_i)}{f_\theta(X_i)} \xrightarrow{P_\theta} E_\theta \log \frac{f_{\theta\pm\epsilon}(X)}{f_\theta(X)} := -\eta_{\theta\pm\epsilon} < 0
    \]

    (The last inequality is guaranteed by Lemma 8.2, since the expected log-likelihood is maximized at the true parameter.)
    Therefore, for any $\delta > 0, \xi > 0$, there exists $m = m_{\epsilon, \delta}$ such that for all $n > m$:

    \[
    P_\theta \left\{ \left| \frac{1}{n}\sum_{i=1}^n \log \frac{f_{\theta\pm\epsilon}(X_i)}{f_\theta(X_i)} + \eta_{\theta\pm\epsilon} \right| < \xi \right\} \ge 1 - \frac{\delta}{2}
    \]

    By choosing $0 < \xi < \min\{\eta_{\theta-\epsilon}, \eta_{\theta+\epsilon}\}$ sufficiently small, the above bounds imply that for $n > m$, with high probability the log-likelihood at the boundary points is lower than at the center:

    \[
    P_\theta(A) := P_\theta \left( \frac{1}{n}\sum_{i=1}^n \log f_\theta(X_i) > \frac{1}{n}\sum_{i=1}^n \log f_{\theta+\epsilon}(X_i) \right) \ge 1 - \frac{\delta}{2}
    \]

    \[
    P_\theta(B) := P_\theta \left( \frac{1}{n}\sum_{i=1}^n \log f_\theta(X_i) > \frac{1}{n}\sum_{i=1}^n \log f_{\theta-\epsilon}(X_i) \right) \ge 1 - \frac{\delta}{2}
    \]

    Since $P(AB) = P(A) - P(AB^C) \ge P(A) - P(B^C) \ge 1 - \frac{\delta}{2} - \frac{\delta}{2} = 1 - \delta$, we have:

    \[
    P_\theta \Big( l_n(\theta-\epsilon) < l_n(\theta) \text{ and } l_n(\theta+\epsilon) < l_n(\theta) \Big) \ge 1 - \delta
    \]

    Since the function values at the endpoints are smaller than at the interior point, and the log-likelihood function $l_n(\theta')$ is continuously differentiable, there must exist a local maximum point inside $(\theta-\epsilon, \theta+\epsilon)$. That is, the probability that a root of the derivative exists is at least $1-\delta$.

**Note**: The above theorem only guarantees the existence of "a consistent root" near the true parameter, but **does not guarantee** that this root is the global maximum likelihood estimator (MLE) we are looking for.

!!! success "Theorem 8.4 (Consistency of MLE)"

    Under the conditions of Theorem 8.3, define $\hat{\theta}_n$ as the root of the likelihood equation (when the equation has exactly one root). If:

    \[
    \lim_{n \rightarrow \infty} P_\theta(\text{the likelihood equation has a single root}) = 1
    \]

    Then:

    \[
    \hat{\theta}_n \xrightarrow{P_\theta} \theta
    \]

The proof is straightforward: when $n$ is large enough, the probability that a root exists in $(\theta-\epsilon, \theta+\epsilon)$ converges to 1 (Theorem 8.3), and simultaneously the probability that there is a single root in the whole space also converges to 1. The intersection of these two events has probability converging to 1, meaning that this unique root must lie in the tiny interval $(\theta-\epsilon, \theta+\epsilon)$.

**Counterexample and Improvement: Cauchy MLE**

For the Cauchy distribution $f_\theta(x) = \frac{1}{\pi\{1+(x-\theta)^2\}}$, the log-likelihood equation often has **multiple local extrema** (Reeds, 1985 proved that the number of roots is asymptotically $2 \times \text{Poisson}(1/\pi) + 1$). The global maximum point can usually be well-separated from other extremum points.
In practice, to avoid falling into wrong local extrema, we often use the **One-step update method** (Newton-Raphson):

1. First, compute a robust initial estimator $\hat{\theta}_0$ with $\sqrt{n}$-consistency (e.g., the sample median).

2. Perform one-step Newton-Raphson update: $\hat{\theta}_n := \hat{\theta}_0 - (\mathbb{P}_n \ddot{m}(\hat{\theta}_0))^{-1} \mathbb{P}_n \dot{m}(\hat{\theta}_0)$

This ensures consistency while improving statistical efficiency.

---

## 5. Asymptotic Normality of MLE

!!! success "Theorem 8.5 (Asymptotic Normality of MLE)"

    Let $X_1, \dots, X_n$ be i.i.d. $\sim P_{\theta_0}$, $\Theta \subset \mathbb{R}$, and in an open neighborhood $\Theta_0$ of the true parameter $\theta_0$ satisfy:
  
    (i) $f_{\theta'}(x) > 0, \forall x, \theta' \in \Theta_0$.
    (ii) $\forall x$, $f_{\theta'}(x)$ is 3 times differentiable in $\Theta_0$.
    (iii) There exists a nonnegative function $M(x)$ with $E_{\theta_0}M(X) < \infty$ such that for any $\theta' \in \Theta_0$, $\left| \frac{d^3}{d\theta'^3}\log f_{\theta'}(x) \right| \le M(x)$.
    (iv) For $l=1,2$, the integral equation $\int f_{\theta'}(x)\mu(dx)=1$ can be differentiated twice under the integral sign.
    (v) $0 < I_1(\theta') < \infty$.
    (vi) $\lim_{n\rightarrow \infty} P_\theta(\hat{\theta}_n \text{ is a root of the likelihood equation}) = 1$ and $E_{\theta} |\log f_{\theta'}(X)| < \infty$.
    (vii) The MLE satisfies consistency $\hat{\theta}_n \rightarrow \theta_0$ and identifiability.

    Then, the asymptotic distribution of the MLE is:

    \[
    \sqrt{n}(\hat{\theta}_n - \theta_0) \xrightarrow{d} N(0, I_1^{-1}(\theta_0))
    \]

??? proof "Proof of Theorem 8.5 (Taylor expansion method, click to expand)"

    Let $l_n(\theta') = \frac{1}{n} \sum_{i=1}^n \log f_{\theta'}(X_i)$. Since $\hat{\theta}_n$ is a root of the score equation, perform a second-order Taylor expansion of $l'_n(\hat{\theta}_n)$ at $\theta_0$:

    \[
    0 = l'_n(\hat{\theta}_n) = l'_n(\theta_0) + l''_n(\theta_0)(\hat{\theta}_n - \theta_0) + \frac{1}{2}l'''_n(\theta_1)(\hat{\theta}_n - \theta_0)^2
    \]

    where $\theta_1$ lies between $\hat{\theta}_n$ and $\theta_0$. By the Weak Law of Large Numbers (WLLN) and the definition of Fisher information:

    \[
    l''_n(\theta_0) = \frac{1}{n} \sum_{i=1}^n \frac{d^2}{d\theta^2}\log f_\theta(X_i) \xrightarrow{P_\theta} -I_1(\theta_0) \in (-\infty, 0)
    \]

    That is, we can write $l''_n(\theta_0) = -I_1(\theta_0) + o_{P_\theta}(1)$.
  
    For the third-order remainder term, by condition (iii):

    \[
    |l'''_n(\theta_1)| = \left| \frac{1}{n}\sum_{i=1}^n \frac{d^3}{d\theta^3}\log f_\theta(X_i) \Big|_{\theta=\theta_1} \right| \le \frac{1}{n}\sum_{i=1}^n M(X_i) \xrightarrow{P_\theta} E_{\theta_0} M(X) < \infty
    \]

    Therefore, the sequence $\{l'''_n(\theta_1)\}$ is bounded in probability, i.e., tight, denoted as $O_{P_\theta}(1)$.
    It is known that $\hat{\theta}_n \xrightarrow{P_\theta} \theta_0$, i.e., $\hat{\theta}_n - \theta_0 = o_{P_\theta}(1)$. Thus the entire third-order remainder term is:

    \[
    (\hat{\theta}_n - \theta_0)^2 l'''_n(\theta_1) = o_{P_\theta}(\hat{\theta}_n - \theta_0)
    \]

    Substituting these asymptotic terms back into the expansion:

    \[
    0 = l'_n(\theta_0) + (-I_1(\theta_0) + o_{P_\theta}(1))(\hat{\theta}_n - \theta_0) + o_{P_\theta}(\hat{\theta}_n - \theta_0)
    \]

    Multiply both sides by $\sqrt{n}$ and rearrange to solve for $\sqrt{n}(\hat{\theta}_n - \theta_0)$:

    \[
    \sqrt{n}(\hat{\theta}_n - \theta_0) = -I_1^{-1}(\theta_0)\sqrt{n}l'_n(\theta_0) + o_{P_\theta}(\sqrt{n}(\hat{\theta}_n - \theta_0))
    \]

    By the Central Limit Theorem (CLT), the standardized score function converges to a normal distribution:

    \[
    \sqrt{n}l'_n(\theta_0) = \frac{1}{\sqrt{n}}\sum_{i=1}^n \frac{\partial \log f_\theta(X_i)}{\partial\theta} \xrightarrow{d} N(0, I_1(\theta_0))
    \]

    Finally, applying Slutsky's theorem (multiplication by a constant):

    \[
    -I_1^{-1}(\theta_0)\sqrt{n}l'_n(\theta_0) \xrightarrow{d} N(0, I_1^{-1}(\theta_0))
    \]

    Hence, $\sqrt{n}(\hat{\theta}_n - \theta_0) \xrightarrow{d} N(0, I_1^{-1}(\theta_0))$.

Moreover, for **compact convex parameter spaces**, we also have the following theorem (Bijma & Jonker & Van der Vaart, 2017):
If the parameter space $\Theta$ is compact and convex, the model is identifiable; the mapping $\vartheta \mapsto \log p_\vartheta(x)$ is continuously differentiable, and the derivative $|l_\vartheta(x)| \le L(x)$ has a square-integrable envelope $E_\theta L^2(X_1) < \infty$; and the true parameter $\theta_0$ is an interior point of $\Theta$, the Fisher information matrix $I(\vartheta)$ is continuous and positive definite, then:

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) \sim N(0, I^{-1}(\theta_0))
\]

---

## 6. Likelihood-based Statistical Inference and Likelihood Ratio Test

From the asymptotic normality of $\hat{\theta}_n$, we can construct confidence intervals or perform hypothesis tests for $\theta$. Using Slutsky's theorem and the continuous mapping theorem, if $I_1(\theta)$ is continuous:

\[
\sqrt{n} I_1^{1/2}(\hat{\theta}_n) (\hat{\theta}_n - \theta) \xrightarrow{d} N_p(0, I_p)
\]

For the multidimensional hypothesis test $H_0: \theta = \theta_0$ VS $H_1: \theta \ne \theta_0$, the commonly used **Wald test statistic** is:

\[
\mathcal{W}_n \hat{=} n(\hat{\theta}_n - \theta_0)^T I_1(\theta_0) (\hat{\theta}_n - \theta_0) \rightarrow \chi^p_2
\]

If $W_n > \chi^2_{p, 1-\alpha}$, we reject the null hypothesis.

### Likelihood Ratio Test (LRT)

According to the Neyman-Pearson lemma, for simple null versus simple alternative, the likelihood ratio test is a uniformly most powerful test (UMP). We consider more general composite hypothesis tests:

\[
H_0: \theta \in \Theta_0 \quad \text{VS} \quad H_1: \theta \in \Theta_1
\]

where the full space $\Theta = \Theta_0 \cup \Theta_1$ and $\Theta_0 \cap \Theta_1 = \emptyset$.
Let $\hat{\theta}_n$ be the global MLE over $\Theta$, and $\hat{\theta}_{n,0}$ be the restricted MLE under the constraint $\Theta_0$.

**Log Likelihood Ratio Statistic** is defined as:

\[
LR_n = -2 \log \frac{L_n(\hat{\theta}_{n,0})}{L_n(\hat{\theta}_n)} = 2 \left\{ l_n(\hat{\theta}_n) - l_n(\hat{\theta}_{n,0}) \right\}
\]

The general decision rule is: when $LR_n$ is too large, reject the null hypothesis $H_0$.

**Viewing LR as an Asymptotic Distance Measure**

If the log-likelihood function is sufficiently smooth (e.g., 3 times continuously differentiable) and $\tilde{\theta}_n$ is $\sqrt{n}$-consistent, by the first-order optimality condition (zero derivative) at the global MLE $\hat{\theta}_n$ and Taylor expansion, we have:

\[
\log L_n(\tilde{\theta}_n) - \log L_n(\hat{\theta}_n) = 0 + \frac{1}{2}(\tilde{\theta}_n - \hat{\theta}_n)^T \frac{\partial^2 \log L_n(\hat{\theta}_n)}{\partial\theta\partial\theta^T} (\tilde{\theta}_n - \hat{\theta}_n) + o_p(1)
\]

Hence, the likelihood ratio $LR_n$ can essentially be viewed as an **asymptotic distance between the two MLEs under the null and the full space** (a quadratic form weighted by Fisher information):

\[
LR_n = (\tilde{\theta}_n - \hat{\theta}_n)^T \left( -\frac{\partial^2 \log L_n(\theta_0)}{\partial\theta\partial\theta^T} \right) (\tilde{\theta}_n - \hat{\theta}_n) + o_p(1)
\]