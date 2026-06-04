# Chapter 14: Bootstrap

---

## 1. From Wald Inference to Bootstrap

In previous chapters, we have obtained asymptotic normality for many estimators. For example:

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) \rightarrow N(0, \Sigma)
\]

For regular M-estimators, we have:

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) \rightarrow N(0, A^{-1}BA^{-1})
\]

For regular maximum likelihood estimators (MLE), we have:

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) \rightarrow N(0, I(\theta_0)^{-1})
\]

**Raising the question:** Since we already have asymptotic normality, why do we still need the Bootstrap?

### 1.1 First-order approximation

Let $T_n$ be a standardized or studentized statistic:

\[
T_n = \frac{\hat{\theta}_n - \theta_0}{\hat{se}_n}
\]

The usual Wald approximation is:

\[
P(T_n \le x) \approx \Phi(x)
\]

This is called a **first-order approximation**. In many regular problems, the error order is:

\[
P(T_n \le x) = \Phi(x) + O(n^{-1/2})
\]

* This approximation successfully captures the center and variance of the limiting distribution.

* However, it typically ignores skewness, bias, and curvature effects.

* These ignored effects are usually of order $O(n^{-1/2})$.

### 1.2 Edgeworth Expansion

To obtain higher-order approximations, we can use an Edgeworth expansion:

\[
P(T_n \le x) = \Phi(x) + n^{-1/2} p_1(x) \phi(x) + n^{-1} p_2(x) \phi(x) + \dots + O(n^{-k/2})
\]

* Here $p_1(x)$ is a polynomial related to the population **skewness**.

* $p_2(x)$ is a polynomial related to the population **kurtosis**.

* However, these polynomials $p_1, p_2$ depend on unknown population moments and are difficult to compute analytically in practice.

---

## 2. Core Idea of Bootstrap

The Bootstrap (resampling method) proposed by Efron (1979) provides an extremely ingenious solution.

**Core idea:** Replace the unknown population distribution $P$ with the empirical distribution $\mathbb{P}_n$.
Sample from $\mathbb{P}_n$ to mimic the process of sampling from the true distribution $P$.

What is most remarkable is that this extremely simple data-driven operation **automatically captures higher-order terms!**

### 2.1 Correspondence between the Real World and the Bootstrap World

* **The Real World**:
  Unknown population distribution $P \implies$ draw sample $X_1, \dots, X_n \implies$ compute estimator $\hat{\theta}_n$

* **The Bootstrap World**:
  Empirical distribution $\mathbb{P}_n \implies$ draw Bootstrap sample $X_1^*, \dots, X_n^* \implies$ compute Bootstrap estimator $\hat{\theta}_n^*$

The cornerstone of Bootstrap inference is: **the distribution of $\hat{\theta}_n^*$ relative to $\hat{\theta}_n$ in the Bootstrap world mimics extremely well the distribution of $\hat{\theta}_n$ relative to $\theta_0$ in the real world.**

### 2.2 Two Main Uses of Bootstrap

* **1. Estimating distribution and standard errors**: When the analytic form of the estimator's variance is too complex or difficult to derive, Bootstrap provides a computation-based nonparametric estimation method.

* **2. Refining approximations**: Used to replace the first-order normal approximation, automatically implementing Edgeworth corrections, thus obtaining confidence intervals with smaller error bounds (e.g., $O(n^{-1})$ instead of $O(n^{-1/2})$).

### 2.3 Notations

* **Original data**: $X_1, \dots, X_n \sim P$

* **Empirical measure**: $\mathbb{P}_n = \frac{1}{n} \sum_{i=1}^n \delta_{X_i}$

* **Bootstrap sample**: $X_1^*, \dots, X_n^*$ are drawn independently **with replacement** from $\mathbb{P}_n$. That is, $X_i^* \sim \mathbb{P}_n$.

* **Bootstrap empirical measure**: $\mathbb{P}_n^* = \frac{1}{n} \sum_{i=1}^n \delta_{X_i^*}$

* **Estimator**: $\hat{\theta}_n = T(\mathbb{P}_n)$

* **Bootstrap estimator**: $\hat{\theta}_n^* = T(\mathbb{P}_n^*)$

---

## 3. Consistency of Bootstrap

!!! info "Definition 10.1 (Consistency of Bootstrap)"

    If, given the original data, the conditional distribution of $\sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n)$ converges weakly to the same limit as the unconditional distribution of $\sqrt{n}(\hat{\theta}_n - \theta_0)$, then the Bootstrap is called **consistent**.

Specifically, define the true cumulative distribution function as:

\[
R_n(x, P) = P(\sqrt{n}(\hat{\theta}_n - \theta_0) \le x)
\]

The Bootstrap estimated distribution function (conditional distribution) is:

\[
R_n(x, \mathbb{P}_n) = P^*(\sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) \le x \mid X_1, \dots, X_n)
\]

Then Bootstrap consistency means:

\[
\sup_x |R_n(x, \mathbb{P}_n) - R_n(x, P)| \xrightarrow{P} 0 \quad (\text{or almost surely a.s.})
\]

### 3.1 Example: Bootstrap for Sample Mean

!!! example "Example 10.2 (Sample Mean)"

    Let $X_i \sim P$, with mean $\mu$ and variance $\sigma^2$. The target parameter $\theta_0 = \mu$, estimator $\hat{\theta}_n = \bar{X}_n$.
  
    By the Central Limit Theorem (CLT), in the real world:

    \[
    \sqrt{n}(\bar{X}_n - \mu) \xrightarrow{d} N(0, \sigma^2)
    \]

    Now examine the Bootstrap world: Draw $X_1^*, \dots, X_n^* \sim \mathbb{P}_n$.
  
    The Bootstrap estimator is $\hat{\theta}_n^* = \bar{X}_n^* = \frac{1}{n}\sum_{i=1}^n X_i^*$.

    Given the original data $(X_1, \dots, X_n)$, compute conditional expectation and conditional variance:

    * **Conditional expectation**:
  
    \[
    E^*(\bar{X}_n^*) = E^*(X_1^*) = \frac{1}{n} \sum_{i=1}^n X_i = \bar{X}_n
    \]

    * **Conditional variance**:
  
    \[
    Var^*(\sqrt{n}\bar{X}_n^*) = Var^*(X_1^*) = \frac{1}{n} \sum_{i=1}^n (X_i - \bar{X}_n)^2 = \hat{\sigma}_n^2
    \]

### 3.2 Bootstrap CLT for Sample Mean

Since $X_1^*, \dots, X_n^*$ are i.i.d. conditional on the original data, we can apply the Lindeberg-Feller CLT or Berry-Esseen theorem to the Bootstrap sample.

Also, because the sample variance is consistent by the law of large numbers: $\hat{\sigma}_n^2 \xrightarrow{p} \sigma^2$, so:

\[
\sup_x \left| P^*(\sqrt{n}(\bar{X}_n^* - \bar{X}_n) \le x) - \Phi(x/\hat{\sigma}_n) \right| \xrightarrow{P} 0
\]

At the same time, since $\hat{\sigma}_n \xrightarrow{p} \sigma$ and the normal cumulative distribution function $\Phi$ is continuous:

\[
\Phi(x/\hat{\sigma}_n) \xrightarrow{P} \Phi(x/\sigma)
\]

Thus we conclude that the Bootstrap for the sample mean is consistent.

**Remark: Why not directly use $N(0, \hat{\sigma}_n^2)$?**
For the simple sample mean, using the normal approximation directly is indeed sufficient. However, for complex estimators (such as sample median, general M-estimators, etc.), the analytic form of the asymptotic variance $\Sigma$ may be very complex or difficult to estimate accurately. Bootstrap cleverly bypasses the step of explicitly computing the asymptotic variance.

---

## 4. Example Where Bootstrap Fails: Maximum

The Bootstrap is **not** consistent in all cases. For extremal statistics or parameters on the boundary, the standard Bootstrap fails.

!!! warning "Example 10.3 (Bootstrap Failure: The Maximum)"

    Suppose $X_i \sim U(0, \theta)$ are i.i.d. samples from a uniform distribution.
  
    We want to estimate the parameter $\theta$. The natural choice is the maximum statistic:
  
    \[
    \hat{\theta}_n = X_{(n)} = \max_{1 \le i \le n}(X_i)
    \]

    In the real world, the maximum has the following extreme value limiting distribution:

    \[
    n(\theta - X_{(n)}) \xrightarrow{d} \text{Exp}(1/\theta)
    \]

    where the limiting distribution is an exponential distribution.

### 4.1 Behavior of the Maximum in the Bootstrap World

Now let's see what happens in the Bootstrap world.
The Bootstrap estimator is $\hat{\theta}_n^* = X_{(n)}^* = \max_{1 \le i \le n}(X_i^*)$.

Compute the probability that the Bootstrap maximum equals the original sample maximum:

\[
P^*(X_{(n)}^* = X_{(n)}) = P^*(\text{at least one } X_i^* \text{ equals } X_{(n)})
\]

This equals $1$ minus the probability that none of the Bootstrap samples equals $X_{(n)}$. Since the probability of not drawing $X_{(n)}$ in one draw is $1 - 1/n$:

\[
P^*(X_{(n)}^* = X_{(n)}) = 1 - \left(1 - \frac{1}{n}\right)^n
\]

As $n \rightarrow \infty$, we have the very famous limit:

\[
1 - \left(1 - \frac{1}{n}\right)^n \rightarrow 1 - e^{-1} \approx 0.632
\]

* **Real world**: Since the population is continuous, $P(X_{(n)} = \theta) = 0$.

* **Bootstrap world**: The Bootstrap distribution has a huge **point mass** of approximately $0.632$ at the point $0$!

Clearly, a Bootstrap distribution containing a huge point mass can never converge to a continuous exponential distribution. Therefore, **Bootstrap fails here**.

### 4.2 Why Does It Fail?

The fundamental reason is the conflict between discreteness and continuity:
The empirical distribution $\mathbb{P}_n$ is essentially a **discrete distribution**, while the true population distribution $P$ is **continuous** at the boundary.
Extremal statistics are extremely sensitive to such local discrete properties of the empirical distribution.

*(A common approach to solve such problems is to use "m-out-of-n bootstrap" or subsampling, where the resampling size $m$ is much smaller than the original sample size $n$.)*

---

## 5. Bootstrap Delta Method

If a statistic $T_n$ is Bootstrap-consistent, is $g(T_n)$ after a nonlinear smooth transformation still consistent? The answer is yes.

!!! success "Theorem 10.4 (Bootstrap Delta Method)"

    Suppose in the real world $\sqrt{n}(T_n - \theta) \xrightarrow{d} T$.
  
    Let the function $g$ be differentiable at $\theta$ with derivative $g'(\theta)$.
  
    If $T_n$ is Bootstrap-consistent, then the Bootstrap is also consistent for $g(T_n)$. That is, conditional on the data, we have:

    \[
    \sqrt{n}(g(T_n^*) - g(T_n)) \xrightarrow{d} g'(\theta) T \quad \text{(under } P^* \text{ measure)}
    \]

??? proof "Proof Outline of Theorem 10.4 (Click to Expand)"

    Use the Mean Value Theorem or Taylor expansion:

    \[
    g(T_n^*) - g(T_n) = g'(\tilde{T}_n^*)(T_n^* - T_n)
    \]

    where $\tilde{T}_n^*$ lies between $T_n^*$ and $T_n$.

    Since $T_n$ is Bootstrap-consistent, under the conditional measure $\sqrt{n}(T_n^* - T_n) = O_{P^*}(1)$, i.e., it is bounded in $P^*$ probability.
  
    At the same time, since in the real world $T_n \xrightarrow{p} \theta$, this implies $T_n^* \xrightarrow{p^*} \theta$.

    Because $\tilde{T}_n^*$ is sandwiched between them, $\tilde{T}_n^* \xrightarrow{p^*} \theta$.

    Since the derivative is continuous, we have $g'(\tilde{T}_n^*) \xrightarrow{p^*} g'(\theta)$.

    Finally, applying Slutsky's theorem to the product, since one part converges in conditional probability to a constant matrix and the other part converges in conditional distribution to $T$, we obtain:

    \[
    \sqrt{n}(g(T_n^*) - g(T_n)) \xrightarrow{d} g'(\theta) T
    \]

    This completes the proof. $\square$

---

## 6. Generalized Delta Method and Functional Derivatives

For parametric models, parameters are typically just functions of moments, and the classical Delta method suffices.
However, in nonparametric problems, the parameter of interest is often a functional of the entire distribution:

\[
\theta = \phi(P)
\]

The corresponding plug-in estimator is:

\[
\hat{\theta}_n = \phi(\mathbb{P}_n)
\]

To apply the Delta method to such functionals, we need the concept of a **generalized derivative for functionals**.

### 6.1 Gateaux Derivative

The Gateaux derivative is analogous to the **directional derivative** in multivariate calculus.

!!! info "Definition (Gateaux Derivative)"

    If the limit:

    \[
    \lim_{t \downarrow 0} \frac{\phi((1-t)P + tH) - \phi(P)}{t} = \int \phi'_P d(H-P)
    \]

    exists and is linear in the direction $(H-P)$, then the functional $\phi$ is said to be **Gateaux differentiable** at $P$ in the direction $H-P$.

### 6.2 Frechet and Hadamard Differentiability

Gateaux differentiability only guarantees convergence along a single fixed direction. In functional analysis this is often insufficient; some form of uniform convergence is required.

* **Frechet differentiability**: requires the limit to be uniform over **all directions**. This is often too strong.

* **Hadamard differentiability**: also called **compact differentiability**. It requires the limit to be uniform over **compact sets of directions**.

!!! info "Definition (Hadamard Differentiability)"

    Let $\phi: \mathcal{P} \to \mathbb{R}$. If for any sequence of directions $H_t$ such that $H_t \to H$ as $t \to 0$, we have:

    \[
    \frac{\phi(P + t H_t) - \phi(P)}{t} \to \phi'_P(H)
    \]

    then $\phi$ is said to be Hadamard differentiable at $P$.

**Why does statistics favor Hadamard differentiability?**
Because, by Prohorov's theorem, weak convergence of empirical processes is inherently related to compact sets in infinite-dimensional spaces. Hadamard differentiability provides exactly the requisite uniformity for applying the Continuous Mapping Theorem and the nonparametric Delta method.

## 7. Asymptotic Linear Estimators and Bootstrap

Many widely used estimators in practice can be expressed as (or approximated by) averages of i.i.d. random variables. This leads to the concept of asymptotic linear estimators.

!!! info "Definition 10.5 (Asymptotically linear estimator)"

    An estimator $\hat{\theta}_n$ is said to be an **asymptotically linear estimator** of the parameter $\theta(P_0)$ if it can be expanded as:

    \[
    \hat{\theta}_n - \theta_0 = \frac{1}{n} \sum_{i=1}^n \dot{\phi}_{P_0}(X_i) + o_P(n^{-1/2})
    \]

    where the function $\dot{\phi}_{P_0}(x)$ is called the **influence function** and satisfies:

    \[
    E_{P_0} \dot{\phi}_{P_0}(X) = 0, \quad E_{P_0} \|\dot{\phi}_{P_0}(X)\|^2 < \infty
    \]

### 7.1 Z-estimators as Asymptotically Linear Estimators

Z-estimators are an excellent example of asymptotically linear estimators. Recalling the theorem on asymptotic normality of Z-estimators, if the estimator satisfies the estimating equation $\Psi_n(\hat{\theta}_n) = 0$, then under certain regularity conditions we have:

\[
\sqrt{n}(\hat{\theta}_n - \theta_0) = \frac{1}{\sqrt{n}} \sum_{i=1}^n (-V_{\theta_0}^{-1} \psi_{\theta_0}(X_i)) + o_P(1)
\]

Comparing with the definition, we see that the Z-estimator is asymptotically linear with **influence function**:

\[
\dot{\phi}_{P_0}(x) = -V_{\theta_0}^{-1} \psi_{\theta_0}(x)
\]

---

## 8. Bootstrap Consistency of Asymptotically Linear Estimators

!!! success "Theorem 10.5 (Bootstrap Consistency of Asymptotically Linear Estimators)"

    Suppose $\hat{\theta}_n$ is an asymptotically linear estimator of $\theta(P_0)$ with influence function $\dot{\phi}_{P_0}(x)$. Let $\hat{\theta}_n^*$ be the bootstrap estimator obtained by applying the same estimation procedure to the bootstrap sample $\mathbb{P}_n^*$.
  
    If $E_{P_0} \|\dot{\phi}_{P_0}(X)\|^2 < \infty$, and conditional on the original sample, under the bootstrap probability measure $P^*$ we have the following asymptotic linear expansion:

    \[
    \sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) = \frac{1}{\sqrt{n}} \sum_{i=1}^n (\dot{\phi}_{P_0}(X_i^*) - \mathbb{P}_n \dot{\phi}_{P_0}) + o_{P^*}(1)
    \]

    then the bootstrap is consistent.

??? proof "Proof of Theorem 10.5 (Click to Expand)"

    **Step 1: Convert to an i.i.d. mean sequence**

    Let $Y_i = \dot{\phi}_{P_0}(X_i)$. By the law of large numbers and the given conditions:
  
    \[
    \bar{Y}_n = \mathbb{P}_n \dot{\phi}_{P_0} \xrightarrow{p} E_{P_0} \dot{\phi}_{P_0}(X) = 0
    \]
  
    In the bootstrap sampling, $X_1^*, \dots, X_n^*$ are i.i.d. draws with replacement from the empirical distribution $\mathbb{P}_n$. Define the corresponding bootstrap random variables:
  
    \[
    Y_i^* = \dot{\phi}_{P_0}(X_i^*)
    \]
  
    Given the original data $X_1, \dots, X_n$, the variables $Y_1^*, \dots, Y_n^*$ are i.i.d., and their expectation and variance under the bootstrap measure are:
  
    * **Conditional expectation**:
  
    \[
    E^*(Y_1^*) = \bar{Y}_n = \mathbb{P}_n \dot{\phi}_{P_0}
    \]
  
    * **Conditional variance**:
  
    \[
    Var^*(Y_1^*) = \frac{1}{n} \sum_{i=1}^n (Y_i - \bar{Y}_n)^2 = \frac{1}{n} \sum_{i=1}^n Y_i^2 - (\bar{Y}_n)^2
    \]

    Since $E_{P_0} Y_1^2 < \infty$ and $E_{P_0} Y_1 = 0$, by the law of large numbers:
  
    \[
    Var^*(Y_1^*) \xrightarrow{p} E_{P_0} Y_1^2 \quad (\text{as } n \rightarrow \infty)
    \]

    **Step 2: Apply the Lindeberg-Feller central limit theorem**

    We need to verify that, conditional on the data, the sequence $\frac{1}{\sqrt{n}} \sum_{i=1}^n (Y_i^* - \bar{Y}_n)$ satisfies the Lindeberg condition. Specifically, for any $\epsilon > 0$, we need to show:
  
    \[
    E^* \left[ (Y_1^* - \bar{Y}_n)^2 I(|Y_1^* - \bar{Y}_n| > \epsilon \sqrt{n}) \right] \xrightarrow{p} 0
    \]
  
    Expand the expectation under the bootstrap measure:
  
    \[
    = \frac{1}{n} \sum_{i=1}^n (Y_i - \bar{Y}_n)^2 I(|Y_i - \bar{Y}_n| > \epsilon \sqrt{n})
    \]
  
    Using the inequality $(a-b)^2 \le 2a^2 + 2b^2$ and noting that when $|Y_i - \bar{Y}_n| > \epsilon \sqrt{n}$, we must have either $|Y_i| > \frac{\epsilon}{2}\sqrt{n}$ or $|\bar{Y}_n| > \frac{\epsilon}{2}\sqrt{n}$, we obtain:
  
    \[
    \le \frac{1}{n} \sum_{i=1}^n 2(Y_i^2 + \bar{Y}_n^2) \left[ I\left(|Y_i| > \frac{\epsilon}{2}\sqrt{n}\right) + I\left(|\bar{Y}_n| > \frac{\epsilon}{2}\sqrt{n}\right) \right]
    \]
  
    Since $E_{P_0} Y_1^2 < \infty$, by the dominated convergence theorem:
  
    \[
    E_{P_0} \left[ Y_1^2 I\left(|Y_1| > \frac{\epsilon}{2}\sqrt{n}\right) \right] \rightarrow 0
    \]
  
    And we know $\bar{Y}_n \xrightarrow{p} 0$. Therefore, the above bound converges to 0 in probability. This verifies the Lindeberg condition.

    **Step 3: Conclude**
  
    By the Lindeberg-Feller central limit theorem, we obtain:
  
    \[
    \frac{1}{\sqrt{n}} \sum_{i=1}^n (Y_i^* - \bar{Y}_n) \xrightarrow{d} N(0, E_{P_0} Y_1^2) \quad \text{(under } P^* \text{)}
    \]
  
    Combined with the given condition:
  
    \[
    \sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) = \frac{1}{\sqrt{n}} \sum_{i=1}^n (Y_i^* - \bar{Y}_n) + o_{P^*}(1)
    \]
  
    By Slutsky’s theorem, the limiting distribution of the bootstrap replicates that of $\sqrt{n}(\hat{\theta}_n - \theta_0)$ in the real world, so the bootstrap is consistent. $\square$

---

## 9. Bootstrap for Z-estimators

For a Z-estimator, the original estimator $\hat{\theta}_n$ satisfies an estimating equation:

\[
\Psi_n(\theta) = \mathbb{P}_n \psi_\theta = 0 \quad (\text{or } = o_P(n^{-1/2}))
\]

In the bootstrap world, the corresponding bootstrap Z-estimator $\hat{\theta}_n^*$ is defined to satisfy the estimating equation based on the empirical measure $\mathbb{P}_n^*$:

\[
\Psi_n^*(\theta) = \mathbb{P}_n^* \psi_\theta = 0 \quad (\text{or } = o_{P^*}(n^{-1/2}))
\]

!!! success "Theorem 10.6 (Bootstrap Consistency of Z-estimators)"

    Assume all standard conditions for asymptotic normality of Z-estimators hold. Further assume the bootstrap estimator $\hat{\theta}_n^*$ satisfies:
  
    \[
    \sqrt{n}\mathbb{P}_n^* \psi_{\hat{\theta}_n^*} = o_{P^*}(1)
    \]
  
    and is consistent in conditional probability: $\hat{\theta}_n^* \xrightarrow{P} \theta_0$.
  
    Then, under the bootstrap measure, analogous to the original Taylor expansion, we have the following asymptotic linear representation:
  
    \[
    \sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) = -V_{\theta_0}^{-1} \sqrt{n}(\mathbb{P}_n^* - \mathbb{P}_n)\psi_{\theta_0} + o_{P^*}(1)
    \]

**Remark on Theorem 10.6:**

* We need to show that $\hat{\theta}_n^*$ is consistent for $\theta_0$ in $P^*$-probability.
* In practical computational implementation, we do not necessarily have to solve this optimization problem exactly, which is often non-convex. To avoid computational cost, we often start from $\hat{\theta}_n$ and perform a **one-step** Newton-Raphson iteration; the resulting one-step bootstrap estimator already guarantees asymptotic consistency.

---

## 10. Bootstrap for M-estimators

For an M-estimator, the original estimator $\hat{\theta}_n$ maximizes the random criterion function $\mathbb{P}_n m_\theta$.
Similarly, the bootstrap M-estimator $\hat{\theta}_n^*$ maximizes the bootstrap criterion function $\mathbb{P}_n^* m_\theta$.

In the theory of M-estimators, we expand the objective function around the maximum as a quadratic (similar to a first-order Taylor expansion of the derivative for Z-estimators).

!!! success "Theorem 10.7 (Bootstrap Consistency of M-estimators)"

    Assume all standard conditions for asymptotic normality of M-estimators hold, and $\hat{\theta}_n^*$ is consistent in conditional probability: $\hat{\theta}_n^* \xrightarrow{P} \theta_0$.
  
    Then, under the bootstrap measure, we have the following expansion:
  
    \[
    \sqrt{n}(\hat{\theta}_n^* - \hat{\theta}_n) = -V_{\theta_0}^{-1} \sqrt{n}(\mathbb{P}_n^* - \mathbb{P}_n)\dot{m}_{\theta_0} + o_{P^*}(1)
    \]
  
    This shows that M-estimators are bootstrap-consistent and can be expressed in an asymptotically linear form.

---

## 11. Estimating Standard Errors and Bias with Bootstrap

Through the bootstrap, we can directly estimate the standard error and bias of a statistic using computation and resampling, without deriving complex analytical asymptotic expressions.

### 11.1 Estimating Standard Errors

The steps to compute the bootstrap standard error are as follows:

1. Draw $B$ independent bootstrap samples with replacement from the original sample.
2. For each bootstrap sample, compute the corresponding estimator, denoted $\hat{\theta}_{n,1}^*, \dots, \hat{\theta}_{n,B}^*$.
3. Compute the sample standard deviation of these $B$ estimators as the bootstrap estimate of the standard error $\widehat{se}_B$:

\[
\widehat{se}_B = \sqrt{ \frac{1}{B-1} \sum_{b=1}^B (\hat{\theta}_{n,b}^* - \bar{\hat{\theta}}_n^*)^2 }
\]

where $\bar{\hat{\theta}}_n^* = \frac{1}{B} \sum_{b=1}^B \hat{\theta}_{n,b}^*$.

### 11.2 Estimating Bias

The bias of an estimator $\hat{\theta}_n$ is defined as:

\[
Bias = E(\hat{\theta}_n) - \theta
\]

By replacing the population distribution with the empirical distribution in the bootstrap world, we can approximate $E(\hat{\theta}_n)$ by $\bar{\hat{\theta}}_n^*$ and $\theta$ by $\hat{\theta}_n$. Therefore, the bootstrap bias estimator is:

\[
\widehat{Bias} = E^*(\hat{\theta}_n^*) - \hat{\theta}_n \approx \bar{\hat{\theta}}_n^* - \hat{\theta}_n
\]

*(Note: In practice, if $\widehat{Bias}$ is small, it is usually not advisable to apply a bias correction directly, as it may increase the overall mean squared error (MSE); correction is only recommended when the bias is significantly larger than the standard error.)*

---

## 12. Bootstrap Confidence Intervals

Bootstrap is often used to construct confidence intervals. Three of the most common methods are introduced below.

### 12.1 Normal Confidence Interval

The simplest method assumes the bootstrap distribution is approximately normal, and constructs the interval using the estimated standard error:

\[
C_n = \left[ \hat{\theta}_n - z_{1-\alpha/2} \widehat{se}, \quad \hat{\theta}_n + z_{1-\alpha/2} \widehat{se} \right]
\]

* **Advantages**: Simple to compute.
* **Disadvantages**: It is not strictly nonparametric; it imposes assumptions of symmetry and normality.

---

### 12.2 Basic Bootstrap Confidence Interval (Pivotal Method)

This method is based on constructing a **pivotal quantity**.

Define the root $R_n = \hat{\theta}_n - \theta$. We approximate the distribution of the real-world root $R_n$ by the distribution of the bootstrap root $R_n^* = \hat{\theta}_n^* - \hat{\theta}_n$.

Let $q_{\alpha/2}^*$ and $q_{1-\alpha/2}^*$ be the $\alpha/2$ and $1-\alpha/2$ quantiles of $R_n^*$, respectively. Since we approximate the true distribution by the bootstrap distribution, we have:

\[
P(q_{\alpha/2}^* \le \hat{\theta}_n - \theta \le q_{1-\alpha/2}^*) \approx 1 - \alpha
\]

Rearranging the inequality for $\theta$, we obtain the confidence interval for $\theta$:

\[
\left[ \hat{\theta}_n - q_{1-\alpha/2}^*, \quad \hat{\theta}_n - q_{\alpha/2}^* \right]
\]

Furthermore, if we let $\theta_p^*$ denote the $p$-quantile of the bootstrap estimator $\hat{\theta}_n^*$ itself, then $q_p^* = \theta_p^* - \hat{\theta}_n$. Substituting into the formula, the basic bootstrap confidence interval can be rewritten in a "pivoted" form:

\[
C_n = \left[ 2\hat{\theta}_n - \theta_{1-\alpha/2}^*, \quad 2\hat{\theta}_n - \theta_{\alpha/2}^* \right]
\]

---

### 12.3 Percentile Confidence Interval

The percentile interval is the most intuitive and most commonly used method in practice. It directly takes the lower and upper quantiles of the bootstrap distribution $\hat{\theta}_n^*$ as the endpoints of the interval:

\[
C_n = \left[ \theta_{\alpha/2}^*, \quad \theta_{1-\alpha/2}^* \right]
\]

**Why is the percentile method valid? (Justification)**

The percentile method has an excellent property: **monotone transformation invariance**.
Suppose there exists an unknown monotone increasing transformation $g$ such that the transformed estimator is exactly normally distributed, i.e.,

\[
g(\hat{\theta}_n) - g(\theta) \sim N(0, c^2)
\]

Since $g$ is monotone, the relative positions of quantiles are preserved under the transformation (quantiles of the transformation equal the transformation of the quantiles):

\[
g(\theta_\alpha^*) = (g(\theta))_\alpha^*
\]

If we construct a confidence interval on the transformed scale for $g(\theta)$, the interval is exactly $[ (g(\theta))_{\alpha/2}^*, (g(\theta))_{1-\alpha/2}^* ]$. Transforming back to the original parameter $\theta$ via $g^{-1}$, we obtain exactly:

\[
\left[ g^{-1}((g(\theta))_{\alpha/2}^*), \quad g^{-1}((g(\theta))_{1-\alpha/2}^*) \right] = \left[ \theta_{\alpha/2}^*, \quad \theta_{1-\alpha/2}^* \right]
\]

This means that **as long as there exists some monotone transformation that normalises the distribution, even if we do not know what this transformation $g$ is, the percentile interval implicitly performs this normalisation and thus yields correct interval bounds**.

## 13. Studentized Bootstrap (Bootstrap-t)

In the standard normal approximation, we typically use a standardized or studentized statistic (i.e., a pivot):

\[
T_n = \frac{\hat{\theta}_n - \theta_0}{\widehat{se}}
\]

In the Bootstrap world, we construct the corresponding studentized statistic:

\[
T_n^* = \frac{\hat{\theta}_n^* - \hat{\theta}_n}{\widehat{se}^*}
\]

* Note that here $\widehat{se}^*$ is the standard error estimate computed from the **current Bootstrap sample**, not from the original sample.

### 13.1 Studentized Bootstrap Confidence Interval

Let $t_{\alpha/2}^*$ and $t_{1-\alpha/2}^*$ be the $\alpha/2$ and $1-\alpha/2$ quantiles of the Bootstrap statistic $T_n^*$, respectively.
We approximate the true distribution using the Bootstrap distribution, i.e., we assume that $P(t_{\alpha/2}^* \le T_n \le t_{1-\alpha/2}^*) \approx 1-\alpha$.

Expanding this and solving for $\theta_0$, we obtain the studentized Bootstrap confidence interval:

\[
C_n = \left[ \hat{\theta}_n - t_{1-\alpha/2}^* \widehat{se}, \quad \hat{\theta}_n - t_{\alpha/2}^* \widehat{se} \right]
\]

### 13.2 Why Studentize? (Second-order accuracy)

* The Basic Bootstrap Interval and the Percentile Interval typically have a coverage error of $O(n^{-1/2})$. This is the same order as the interval based on the normal approximation (first-order accurate).

* **The studentized Bootstrap is second-order accurate**. Its coverage error is $O(n^{-1})$.

**Reason: Edgeworth Expansion**

The studentized statistic $T_n$ is an **asymptotic pivot**; its limiting distribution $N(0,1)$ does not depend on any unknown population parameters.
When we resample the pivot in the Bootstrap, the Bootstrap can automatically estimate and correct the **skewness correction term** of order $n^{-1/2}$ in the Edgeworth expansion, thereby eliminating the first-order error and reducing the final approximation error to $O(n^{-1})$.

---

## 14. Double Bootstrap (Iterated Bootstrap)

In the studentized Bootstrap, we need to compute $T_n^* = (\hat{\theta}_n^* - \hat{\theta}_n) / \widehat{se}^*$. However, in some complex models, the standard error $\widehat{se}^*$ may not have an analytic expression or may be very difficult to compute directly for each Bootstrap sample.

**Solution: Nest another Bootstrap inside the Bootstrap!**

**Specific algorithm steps:**

1. Draw a first-level Bootstrap sample $\mathbb{P}_n^*$ with replacement from the original empirical distribution $\mathbb{P}_n$. Compute the estimator $\hat{\theta}_n^*$ using this sample.

2. From the current first-level Bootstrap sample $\mathbb{P}_n^*$, draw a second-level (inner) Bootstrap sample $\mathbb{P}_n^{**}$ with replacement again, and compute $\hat{\theta}_n^{**}$.

3. Repeat step 2 $M$ times, and compute the standard error $\widehat{se}^*$ as the empirical standard deviation of these $M$ values $\hat{\theta}_n^{**}$.

4. Compute the current Bootstrap pivot $T_n^* = \frac{\hat{\theta}_n^* - \hat{\theta}_n}{\widehat{se}^*}$.

5. Repeat steps 1 through 4 $B$ times, obtaining the full empirical distribution of $T_n^*$, from which quantiles are extracted to construct the confidence interval.

* **Disadvantage**: Extremely high computational cost. If the outer loop requires $B$ resamples and the inner loop requires $M$ resamples, a total of $B \times M$ model fits are needed.

---

## 15. Bootstrap for Regression

Consider the standard linear or nonlinear regression model:

\[
Y_i = x_i^T \beta + \epsilon_i
\]

In regression analysis, how we resample depends on the degree of assumption we are willing to make about the model. There are generally two main approaches:

### 15.1 Method 1: Residual Bootstrap (Model-based Bootstrap)

If we strongly believe in the correctness of the model structure and assume that the errors $\epsilon_i$ are i.i.d. and independent of the covariates $x_i$ (homoscedasticity).

**Steps:**

1. Fit the model to the original data, obtaining the parameter estimate $\hat{\beta}$ and the corresponding residuals $\hat{\epsilon}_i = Y_i - x_i^T \hat{\beta}$.

2. Center the residuals (recommended): $\tilde{\epsilon}_i = \hat{\epsilon}_i - \bar{\hat{\epsilon}}$.

3. Randomly draw $n$ Bootstrap residuals $\epsilon_1^*, \dots, \epsilon_n^*$ with replacement from the set of centered residuals $\{\tilde{\epsilon}_1, \dots, \tilde{\epsilon}_n\}$.

4. Keep the covariates $x_i$ unchanged, and generate new Bootstrap response variables:

\[
Y_i^* = x_i^T \hat{\beta} + \epsilon_i^*
\]

5. Refit the model using the generated $(x_i, Y_i^*)$ to obtain the Bootstrap estimator $\hat{\beta}^*$.

* **Characteristics**: Fully utilizes the parametric structure of the model; extremely efficient when the model is correctly specified and the errors satisfy the homoscedasticity assumption.

### 15.2 Method 2: Pairs Bootstrap (Nonparametric Bootstrap)

If we are uncertain about homoscedasticity (i.e., heteroscedasticity may be present) or doubt the correctness of the linear model, we can use a more nonparametric approach.

**Steps:**

1. Treat each observation as a data pair $Z_i = (x_i, Y_i)$.

2. Randomly draw $n$ data pairs with replacement directly from the dataset $\{Z_1, \dots, Z_n\}$, denoted as $Z_1^*, \dots, Z_n^*$ (where $Z_i^* = (x_i^*, Y_i^*)$).

3. Refit the model to this new set of pairs to obtain the Bootstrap estimator $\hat{\beta}^*$.

* **Characteristics**: Highly robust to heteroscedasticity and model misspecification. This is because the resampling preserves the joint distribution structure between $x_i$ and $Y_i$.

---

## 16. When does Bootstrap fail?

Although Bootstrap is extremely powerful, it is not a panacea. In the following typical scenarios, the standard Bootstrap fails (i.e., consistency does not hold):

### 16.1 1. Heavy tails / Infinite variance

* If the population distribution has infinite variance (e.g., $E[X^2] = \infty$), then the sample mean does not satisfy the classical normal limiting distribution (it may converge in distribution to some $\alpha$-stable distribution).

* In this case, the conditional distribution of the Bootstrap sample mean does not converge weakly to a deterministic non-random measure; it remains random in the limit.

### 16.2 2. Boundary parameters

* The **maximum statistic** we discussed earlier is a typical example (e.g., $X_i \sim U(0, \theta)$ with estimator $X_{(n)}$).

* When the parameter lies on the boundary of the parameter space, the discreteness of the empirical distribution $\mathbb{P}_n$ introduces a large probability "point mass" at the extreme value, making it impossible to mimic the true continuous limiting distribution.

### 16.3 3. Non-differentiable functionals

* If the parameter of interest is a non-smooth or non-differentiable functional of the empirical measure, the Bootstrap may fail.

* **Example 1: Hodges estimator** (super-efficient estimator). Such estimators exhibit a non-differentiable "shrinkage" behavior at a point (e.g., 0), so the Bootstrap cannot capture its local abrupt nature.

* **Example 2: Sample median when there are ties**. Although the sample median is Bootstrap-consistent for continuous distributions, if the true distribution is not sufficiently smooth around the population median (e.g., a distribution with a discontinuity or high discreteness), the performance of the Bootstrap median can be very unstable.

---

## 17. Remedies for Failure: Subsampling and m-out-of-n Bootstrap

When faced with situations where the standard Bootstrap fails, we usually adopt the following two alternative strategies:

### 17.1 Subsampling

* **Mechanism**: Draw **without replacement** a subsample of size $m$ from the original sample of size $n$.

* **Condition**: Requires that the subsample size $m \rightarrow \infty$ and $m/n \rightarrow 0$.

* **Advantage**: The conditions for consistency are very mild. It only requires that the original statistic converges to some non-degenerate limiting distribution under the true measure, without requiring the functional to be differentiable or not on the boundary.

### 17.2 m-out-of-n Bootstrap

* **Mechanism**: Draw **with replacement** a Bootstrap sample of size $m$ from the original sample of size $n$.

* **Condition**: Also requires $m \rightarrow \infty$ and $m/n \rightarrow 0$.

* **Advantage**: By reducing the resampling size $m$, it breaks the excessive "ties" and clustering at extreme values caused by resampling $n$ times, thereby restoring consistency in boundary parameter and some heavy-tailed cases.

---

## 18. Final Summary

1. **First-order inference**: Asymptotic normality theory provides the primary first-order inference tool (Asymptotic normality gives first-order inference).

2. **Finite-sample error**: The Edgeworth expansion theoretically explains the leading finite-sample error sources (Edgeworth expansion explains the leading finite-sample error).

3. **Automatic correction**: The power of Bootstrap lies in its ability to automatically estimate the higher-order correction terms in the Edgeworth expansion using the data (Bootstrap often estimates Edgeworth correction terms automatically).

4. **Second-order accuracy**: **Studentized Bootstrap** is the main and recommended route to obtain second-order accurate confidence intervals (Studentized bootstrap is the main route to second-order accurate intervals).

5. **Theoretical consistency**: For M-estimators, Z-estimators, and maximum likelihood estimators (MLE), the consistency of the Bootstrap can be derived by the same linearization argument used to prove asymptotic normality.

6. **Limitations**: Bootstrap is an extremely powerful tool, but it is **not universally valid** (Bootstrap is powerful, but it is not universally valid). Pay attention to the prerequisites before using it.

!!! success "One-line summary"

    The essence of Bootstrap is to imitate the real world with a "parallel universe": the real fluctuation $\hat{\theta}_n - \theta_0$ is perfectly mimicked by $\hat{\theta}_n^* - \hat{\theta}_n$ in the Bootstrap world; and the invisible god—the population distribution $P$, is replaced by the data in our hands—the empirical distribution $\mathbb{P}_n$.