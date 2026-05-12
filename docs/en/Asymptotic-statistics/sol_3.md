---
tags:
  - Large Sample Theory
  - Statistics
  - Homework Exercise
---

# 📈 Detailed Solutions to Homework Exercises: Large Sample Theory Assignment 3

## Exercise 1: Second-Order Delta Method

**Problem:**
Let $X_1, ..., X_n \stackrel{i.i.d}{\sim} Bernoulli(\theta)$.

(a) When $\theta = 1/2$, find the asymptotic distribution of $\hat{\theta}_n = \overline{X}_n$.

(b) The Kullback-Leibler (KL) divergence from $P_\theta$ to $P_\eta$ is defined as:

\[
D_{KL}(P_\eta || P_\theta) := -\mathbb{E}_\eta \log \frac{p_\theta}{p_\eta}(X) \quad , X \sim P_\eta
\]

where $p_\theta, p_\eta$ are the probability mass (density) functions of $P_\theta, P_\eta$ respectively. Prove that:

\[
n D_{KL}(P_\theta || P_{\hat{\theta}_n}) \xrightarrow{d} \frac{1}{2}\chi_1^2
\]

??? success "Solution (click to expand)"

    **(a) Asymptotic distribution**

    Since $X_1, ..., X_n$ are i.i.d. Bernoulli random variables, we have:

    \[
    E[X_1] = \theta \quad \text{and} \quad Var(X_1) = \theta(1-\theta)
    \]

    The estimator is the sample mean $\hat{\theta}_n = \overline{X}_n = \frac{1}{n}\sum_{i=1}^n X_i$. By the classical Central Limit Theorem (CLT),

    \[
    \sqrt{n}(\overline{X}_n - \theta) \xrightarrow{d} \mathcal{N}(0, \theta(1-\theta))
    \]

    Substituting $\theta = 1/2$, the mean is $1/2$ and the variance is $(1/2)(1 - 1/2) = 1/4$. Therefore,

    \[
    \sqrt{n} \left(\overline{X}_n - \frac{1}{2}\right) \xrightarrow{d} \mathcal{N} \left(0, \frac{1}{4}\right)
    \]

    Equivalently, we can write it in standard normal form:

    \[
    \sqrt{n}(2\overline{X}_n - 1) \xrightarrow{d} \mathcal{N}(0, 1)
    \]

    **(b) Asymptotic distribution of the KL divergence**

    The KL divergence from $P_\theta$ to $P_{\hat{\theta}_n}$ is defined as:

    \[
    D_{KL}(P_\theta || P_{\hat{\theta}_n}) = \mathbb{E}_\theta \left[ \log \frac{p_\theta(X)}{p_{\hat{\theta}_n}(X)} \right]
    \]

    For the Bernoulli distribution, the probability mass function is $p_t(x) = t^x (1-t)^{1-x}$. To apply the Delta method, we treat the KL divergence as a function $g(t)$ of $t = \hat{\theta}_n$:

    \[
    g(t) = D_{KL}(P_\theta || P_t) = \mathbb{E}_\theta \left[ X \log \frac{\theta}{t} + (1-X) \log \frac{1-\theta}{1-t} \right]
    \]

    Since the expectation is taken under the true $\theta$ (i.e., $\mathbb{E}_\theta[X] = \theta$), this becomes:

    \[
    g(t) = \theta \log \left(\frac{\theta}{t}\right) + (1-\theta) \log \left(\frac{1-\theta}{1-t}\right)
    \]

    Our goal is to find the asymptotic distribution of $n g(\hat{\theta}_n) = n(g(\hat{\theta}_n) - g(\theta))$ (since $g(\theta) = D_{KL}(P_\theta || P_\theta) = 0$).

    First, compute the first and second derivatives of $g(t)$:

    \[
    g'(t) = -\frac{\theta}{t} + \frac{1-\theta}{1-t}
    \]

    Evaluated at the true parameter $t = \theta$:

    \[
    g'(\theta) = -\frac{\theta}{\theta} + \frac{1-\theta}{1-\theta} = -1 + 1 = 0
    \]

    Since the first derivative is zero, we need to use the **second-order Delta method**. Compute the second derivative:

    \[
    g''(t) = \frac{\theta}{t^2} + \frac{1-\theta}{(1-t)^2}
    \]

    Evaluated at $t = \theta$:

    \[
    g''(\theta) = \frac{\theta}{\theta^2} + \frac{1-\theta}{(1-\theta)^2} = \frac{1}{\theta} + \frac{1}{1-\theta} = \frac{1}{\theta(1-\theta)}
    \]

    From the CLT, $\sqrt{n}(\hat{\theta}_n - \theta) \xrightarrow{d} \mathcal{N}(0, \sigma^2)$ with asymptotic variance $\sigma^2 = \theta(1-\theta)$.

    According to the second-order Delta method, since $g'(\theta) = 0$ and $g''(\theta) \neq 0$, we have:

    \[
    n(g(\hat{\theta}_n) - g(\theta)) \xrightarrow{d} \frac{1}{2} g''(\theta) \sigma^2 \chi_1^2
    \]

    Substituting $\sigma^2$ and $g''(\theta)$:

    \[
    g''(\theta)\sigma^2 = \left( \frac{1}{\theta(1-\theta)} \right) \cdot \theta(1-\theta) = 1
    \]

    This simplifies perfectly to 1. Hence,

    \[
    n D_{KL}(P_\theta || P_{\hat{\theta}_n}) = n(g(\hat{\theta}_n) - 0) \xrightarrow{d} \frac{1}{2} \chi_1^2
    \]

---

## Exercise 2: Hodges (1951) Superefficiency Estimator

**Problem:**
Suppose $X_1, X_2, \dots, X_n$ are i.i.d. observations from $\mathcal{N}(\theta, 1)$ with $\theta \in \mathbb{R}$. Then $\overline{X}_n = n^{-1}\sum_{i=1}^n X_i$ is the MLE of $\theta$. Show that when $\theta \neq 0$,

\[
P(|\overline{X}_n| \le n^{-1/4}) \rightarrow 0
\]

??? success "Solution (click to expand)"

    **Proof:**

    Since $X_1, X_2, \dots, X_n \stackrel{i.i.d}{\sim} \mathcal{N}(\theta, 1)$, the exact distribution of the sample mean is:

    \[
    \overline{X}_n \sim \mathcal{N} \left( \theta, \frac{1}{n} \right)
    \]

    Standardize $\overline{X}_n$ to a standard normal random variable $Z$:

    \[
    Z = \sqrt{n}(\overline{X}_n - \theta) \sim \mathcal{N}(0, 1)
    \]

    Now rewrite the target probability using $Z$:

    \[
    P(|\overline{X}_n| \le n^{-1/4}) = P(-n^{-1/4} \le \overline{X}_n \le n^{-1/4})
    \]

    \[
    = P(-n^{-1/4} - \theta \le \overline{X}_n - \theta \le n^{-1/4} - \theta)
    \]

    \[
    = P(-n^{1/4} - \sqrt{n}\theta \le \sqrt{n}(\overline{X}_n - \theta) \le n^{1/4} - \sqrt{n}\theta)
    \]

    \[
    = P(-n^{1/4} - \sqrt{n}\theta \le Z \le n^{1/4} - \sqrt{n}\theta)
    \]

    To show that this probability converges to 0 when $\theta \neq 0$, consider two cases:

    **Case 1: $\theta > 0$**

    Examine the upper bound of $Z$:

    \[
    n^{1/4} - \sqrt{n}\theta = \sqrt{n}(n^{-1/4} - \theta)
    \]

    As $n \rightarrow \infty$, $n^{-1/4} \rightarrow 0$. Since $\theta > 0$, $(n^{-1/4} - \theta) \rightarrow -\theta < 0$. Therefore,

    \[
    \lim_{n \rightarrow \infty} (n^{1/4} - \sqrt{n}\theta) = -\infty
    \]

    The probability is bounded above by the CDF of the standard normal at the upper bound:

    \[
    P(-n^{1/4} - \sqrt{n}\theta \le Z \le n^{1/4} - \sqrt{n}\theta) \le P(Z \le n^{1/4} - \sqrt{n}\theta)
    \]

    As $n \rightarrow \infty$, $P(Z \le n^{1/4} - \sqrt{n}\theta) \rightarrow \Phi(-\infty) = 0$.

    **Case 2: $\theta < 0$**

    Examine the lower bound of $Z$:

    \[
    -n^{1/4} - \sqrt{n}\theta = -\sqrt{n}(n^{-1/4} + \theta)
    \]

    As $n \rightarrow \infty$, $n^{-1/4} \rightarrow 0$. Since $\theta < 0$, $(n^{-1/4} + \theta) \rightarrow \theta < 0$. Therefore,

    \[
    \lim_{n \rightarrow \infty} -\sqrt{n}(n^{-1/4} + \theta) = +\infty
    \]

    The probability is bounded above by the right tail probability of the standard normal:

    \[
    P(-n^{1/4} - \sqrt{n}\theta \le Z \le n^{1/4} - \sqrt{n}\theta) \le P(Z \ge -n^{1/4} - \sqrt{n}\theta)
    \]

    As $n \rightarrow \infty$, $P(Z \ge -n^{1/4} - \sqrt{n}\theta) \rightarrow 1 - \Phi(+\infty) = 0$.

    In both cases where $\theta \neq 0$, the probability vanishes. Hence, as $n \rightarrow \infty$,

    \[
    P(|\overline{X}_n| \le n^{-1/4}) \rightarrow 0
    \]

---

## Exercise 3: Irregular Maximum Likelihood Estimator

**Problem:**
Let $X_1, X_2, \dots, X_n$ be a random sample from the uniform distribution $U(0, \theta)$ ($\theta > 0$).

(a) Obtain the MLE $\hat{\theta}_n$ of $\theta$ and show that $\sqrt{n}(\hat{\theta}_n - \theta) \xrightarrow{p} 0$.

(b) Find an appropriate scaling factor such that the asymptotic distribution is non-degenerate, and find that asymptotic distribution.

??? success "Solution (click to expand)"

    **(a) Derive the MLE and prove convergence**

    The likelihood function is:

    \[
    L(\theta) = \prod_{i=1}^n \frac{1}{\theta} I(0 \le X_i \le \theta) = \frac{1}{\theta^n} I(\max_{1 \le i \le n} X_i \le \theta)
    \]

    where $I(\cdot)$ is the indicator function. To maximize $L(\theta)$, we need to minimize $\theta^n$ subject to $\theta \ge \max_{1 \le i \le n} X_i$.

    Since $\theta^n$ is strictly increasing in $\theta$ for $\theta > 0$, the value that maximizes the likelihood is the smallest possible $\theta$ satisfying the constraint, i.e., the maximum order statistic:

    \[
    \hat{\theta}_n = X_{(n)} = \max_{1 \le i \le n} X_i
    \]

    To show $\sqrt{n}(\hat{\theta}_n - \theta) \xrightarrow{p} 0$, analyze the bias $\theta - \hat{\theta}_n$. Clearly $\hat{\theta}_n \le \theta$ almost surely, so $\theta - \hat{\theta}_n \ge 0$. Use Markov's inequality.

    First compute the expectation of $\hat{\theta}_n$. The PDF of $X_{(n)}$ for $0 < x < \theta$ is:

    \[
    f_{(n)}(x) = \frac{n x^{n-1}}{\theta^n}
    \]

    Thus the expectation is:

    \[
    E[\hat{\theta}_n] = \int_0^\theta x \cdot \frac{n x^{n-1}}{\theta^n} dx = \frac{n}{\theta^n} \left[ \frac{x^{n+1}}{n+1} \right]_0^\theta = \frac{n}{n+1}\theta
    \]

    The expected bias is:

    \[
    E[\theta - \hat{\theta}_n] = \theta - \frac{n}{n+1}\theta = \frac{\theta}{n+1}
    \]

    For any $\epsilon > 0$, by Markov's inequality:

    \[
    P(|\sqrt{n}(\hat{\theta}_n - \theta)| \ge \epsilon) = P(\sqrt{n}(\theta - \hat{\theta}_n) \ge \epsilon) \le \frac{E[\sqrt{n}(\theta - \hat{\theta}_n)]}{\epsilon} = \frac{\sqrt{n}\theta}{\epsilon(n+1)}
    \]

    As $n \rightarrow \infty$, $\frac{\sqrt{n}}{n+1} \rightarrow 0$. The probability tends to 0, which means:

    \[
    \sqrt{n}(\hat{\theta}_n - \theta) \xrightarrow{p} 0
    \]

    **(b) Appropriate scaling and non-degenerate asymptotic distribution**

    From (a) we have $\sqrt{n}(\theta - \hat{\theta}_n) \xrightarrow{p} 0$, so $\hat{\theta}_n$ converges to $\theta$ at a rate faster than $\sqrt{n}$. The appropriate scaling factor here is $n$.

    We study the asymptotic distribution of the non‑negative random variable $Y_n = n(\theta - \hat{\theta}_n)$ by deriving the limit of its CDF. For any $x > 0$:

    \[
    P(n(\theta - \hat{\theta}_n) \le x) = P \left( \theta - \hat{\theta}_n \le \frac{x}{n} \right)
    \]

    \[
    = P \left( \hat{\theta}_n \ge \theta - \frac{x}{n} \right) = 1 - P \left( \hat{\theta}_n < \theta - \frac{x}{n} \right)
    \]

    Since $\hat{\theta}_n = \max_{1 \le i \le n} X_i$, the sample maximum is strictly less than some value iff every individual observation is strictly less than that value. Because $X_i$ are i.i.d. $U(0, \theta)$,

    \[
    P \left( \hat{\theta}_n < \theta - \frac{x}{n} \right) = \prod_{i=1}^n P \left( X_i < \theta - \frac{x}{n} \right) = \left( \frac{\theta - x/n}{\theta} \right)^n = \left( 1 - \frac{x/\theta}{n} \right)^n
    \]

    Taking the limit as $n \rightarrow \infty$ using the standard limit definition of the exponential:

    \[
    \lim_{n \rightarrow \infty} \left( 1 - \frac{x/\theta}{n} \right)^n = e^{-x/\theta}
    \]

    Therefore, the asymptotic CDF of $n(\theta - \hat{\theta}_n)$ is:

    \[
    \lim_{n \rightarrow \infty} P(n(\theta - \hat{\theta}_n) \le x) = 1 - e^{-x/\theta} \quad (\text{for } x > 0)
    \]

    This is exactly the CDF of an **exponential distribution** with mean $\theta$ (rate parameter $1/\theta$).

    Conclusion, the non-degenerate asymptotic distribution is:

    \[
    n(\theta - \hat{\theta}_n) \xrightarrow{d} Exp(mean = \theta)
    \]

## Exercise 4: Different Convergence Rates of MLEs

**Problem:**
Suppose $X_1, ..., X_n \stackrel{i.i.d}{\sim} Exp(\mu, \sigma)$ with density $f(x) = \frac{1}{\sigma}e^{-\frac{x-\mu}{\sigma}}$, where $x \ge \mu$. Let $(\hat{\mu}, \hat{\sigma})$ be the MLE of $(\mu, \sigma)$.
Find the asymptotic distributions of $\hat{\mu}$ and $\hat{\sigma}$ respectively. (Marginal distributions, not joint.)

??? success "Solution (click to expand)"

    **Step 1: Derive the MLEs of $\mu$ and $\sigma$**

    The likelihood function is:

    \[
    L(\mu, \sigma) = \prod_{i=1}^n \frac{1}{\sigma}e^{-\frac{X_i-\mu}{\sigma}}I(X_i \ge \mu) = \frac{1}{\sigma^n} \exp \left\{ -\frac{\sum_{i=1}^n X_i - n\mu}{\sigma} \right\} I(\mu \le \min_{1 \le i \le n} X_i)
    \]

    Under the constraint $\mu \le X_{(1)}$, the log-likelihood is:

    \[
    \ln L(\mu, \sigma) = -n \ln \sigma - \frac{\sum_{i=1}^n X_i - n\mu}{\sigma}
    \]

    First maximize with respect to $\mu$. Note that:

    \[
    \frac{\partial \ln L}{\partial \mu} = \frac{n}{\sigma} > 0
    \]

    Thus the log-likelihood is strictly increasing in $\mu$. To maximize it, we must make $\mu$ as large as possible subject to $\mu \le X_{(1)}$.
    Hence, the MLE of $\mu$ is the sample minimum:

    \[
    \hat{\mu} = X_{(1)} = \min_{1 \le i \le n} X_i
    \]

    Next, differentiate with respect to $\sigma$ and set to zero:

    \[
    \frac{\partial \ln L}{\partial \sigma} = -\frac{n}{\sigma} + \frac{\sum_{i=1}^n X_i - n\mu}{\sigma^2} = 0
    \]

    Solving gives $\hat{\sigma} = \overline{X}_n - \mu$. Substituting $\hat{\mu} = X_{(1)}$ yields the MLE of $\sigma$:

    \[
    \hat{\sigma} = \overline{X}_n - X_{(1)}
    \]

    **Step 2: Asymptotic distribution of $\hat{\mu}$**

    We use the CDF method for $\hat{\mu} = X_{(1)}$. The appropriate scaling factor is $n$.
    For any $x > 0$:

    \[
    P(n(\hat{\mu} - \mu) > x) = P \left( \hat{\mu} > \mu + \frac{x}{n} \right)
    \]

    \[
    = P \left( \min_{1 \le i \le n} X_i > \mu + \frac{x}{n} \right) = \prod_{i=1}^n P \left( X_i > \mu + \frac{x}{n} \right)
    \]

    Since $X_i \sim Exp(\mu, \sigma)$, for $t \ge \mu$, $P(X_i > t) = e^{-(t-\mu)/\sigma}$. Therefore:

    \[
    P \left( X_i > \mu + \frac{x}{n} \right) = \exp \left\{ -\frac{\mu + x/n - \mu}{\sigma} \right\} = e^{-\frac{x}{n\sigma}}
    \]

    Substituting back into the product:

    \[
    P(n(\hat{\mu} - \mu) > x) = \left( e^{-\frac{x}{n\sigma}} \right)^n = e^{-x/\sigma}
    \]

    This means the limiting CDF is $1 - e^{-x/\sigma}$, which is the CDF of an **exponential distribution** with mean $\sigma$ (rate $1/\sigma$). Hence:

    \[
    n(\hat{\mu} - \mu) \xrightarrow{d} Exp(mean = \sigma)
    \]

    **Step 3: Asymptotic distribution of $\hat{\sigma}$**

    We need the asymptotic distribution of $\hat{\sigma} = \overline{X}_n - X_{(1)}$.

    First, the moments of $X_1$ are: $E[X_1] = \mu + \sigma$ and $Var(X_1) = \sigma^2$. By the standard CLT for $\overline{X}_n$:

    \[
    \sqrt{n}(\overline{X}_n - (\mu + \sigma)) \xrightarrow{d} \mathcal{N}(0, \sigma^2)
    \]

    Rewrite the standardized estimator $\sqrt{n}(\hat{\sigma} - \sigma)$:

    \[
    \sqrt{n}(\hat{\sigma} - \sigma) = \sqrt{n}(\overline{X}_n - X_{(1)} - \sigma)
    \]

    \[
    = \sqrt{n}(\overline{X}_n - \mu - \sigma) - \sqrt{n}(X_{(1)} - \mu)
    \]

    From Step 2, $n(X_{(1)} - \mu) \xrightarrow{d} Exp(\sigma)$, which implies $n(X_{(1)} - \mu) = O_p(1)$.
    Consequently,

    \[
    \sqrt{n}(X_{(1)} - \mu) = \frac{1}{\sqrt{n}} \cdot n(X_{(1)} - \mu) = o(1) \cdot O_p(1) = o_p(1)
    \]

    That is, $\sqrt{n}(X_{(1)} - \mu) \xrightarrow{p} 0$.

    By Slutsky’s theorem, since the second term converges in probability to 0, the asymptotic distribution is entirely determined by the first term. Hence:

    \[
    \sqrt{n}(\hat{\sigma} - \sigma) \xrightarrow{d} \mathcal{N}(0, \sigma^2)
    \]

---

## Exercise 5: Inconsistent MLE and the Neyman-Scott Problem

**Problem:**
Let $X_{ij}$ ($i=1,...,n$, $j=1,...,k$) be independent with $X_{ij} \sim \mathcal{N}(\mu_i, \sigma^2)$.
This is essentially a balanced one‑way ANOVA design; assume $k$ is fixed and $n \rightarrow \infty$, so that the number of parameters grows with the sample size.
Show that the MLE $\hat{\sigma}^2$ of the common variance satisfies:

\[
\hat{\sigma}^2 \xrightarrow{p} \left(1 - \frac{1}{k}\right)\sigma^2
\]

??? success "Solution (click to expand)"

    **Step 1: Derive the MLEs of $\mu_i$ and $\sigma^2$**

    The likelihood function is:

    \[
    L(\mu_1, ..., \mu_n, \sigma^2) = \prod_{i=1}^n \prod_{j=1}^k \frac{1}{\sqrt{2\pi\sigma^2}} \exp \left( -\frac{(X_{ij} - \mu_i)^2}{2\sigma^2} \right)
    \]

    Taking logarithms gives:

    \[
    \ln L = -\frac{nk}{2} \ln(2\pi) - \frac{nk}{2} \ln(\sigma^2) - \frac{1}{2\sigma^2} \sum_{i=1}^n \sum_{j=1}^k (X_{ij} - \mu_i)^2
    \]

    To find the MLE of $\mu_i$, take the derivative and set it to zero:

    \[
    \frac{\partial \ln L}{\partial \mu_i} = \frac{1}{\sigma^2} \sum_{j=1}^k (X_{ij} - \mu_i) = 0 \implies \hat{\mu}_i = \frac{1}{k} \sum_{j=1}^k X_{ij} = \overline{X}_{i.}
    \]

    Next, differentiate with respect to $\sigma^2$ and set to zero:

    \[
    \frac{\partial \ln L}{\partial (\sigma^2)} = -\frac{nk}{2\sigma^2} + \frac{1}{2(\sigma^2)^2} \sum_{i=1}^n \sum_{j=1}^k (X_{ij} - \mu_i)^2 = 0
    \]

    Substituting $\hat{\mu}_i = \overline{X}_{i.}$ gives the MLE of the common variance:

    \[
    \hat{\sigma}^2 = \frac{1}{nk} \sum_{i=1}^n \sum_{j=1}^k (X_{ij} - \overline{X}_{i.})^2
    \]

    **Step 2: Prove convergence in probability**

    Let $Y_i = \sum_{j=1}^k (X_{ij} - \overline{X}_{i.})^2$. For each group $i$, because $X_{ij} \sim \mathcal{N}(\mu_i, \sigma^2)$ for $j=1,...,k$, the scaled sum of squares follows a chi‑square distribution:

    \[
    \frac{Y_i}{\sigma^2} = \frac{1}{\sigma^2} \sum_{j=1}^k (X_{ij} - \overline{X}_{i.})^2 \sim \chi_{k-1}^2
    \]

    Hence $Y_i \sim \sigma^2 \chi_{k-1}^2$. Since samples across groups ($i=1,...,n$) are independent, the sequence $Y_1, Y_2, ..., Y_n$ consists of i.i.d. random variables.

    The expectation of $Y_i$ is $E[Y_i] = \sigma^2(k-1)$.

    When $n \rightarrow \infty$ (while $k$ is fixed), apply the Weak Law of Large Numbers (WLLN) to $\{Y_i\}_{i=1}^n$:

    \[
    \frac{1}{n} \sum_{i=1}^n Y_i \xrightarrow{p} E[Y_1] = (k-1)\sigma^2
    \]

    Now rewrite the estimator $\hat{\sigma}^2$ in terms of this sequence:

    \[
    \hat{\sigma}^2 = \frac{1}{nk} \sum_{i=1}^n Y_i = \frac{1}{k} \left( \frac{1}{n} \sum_{i=1}^n Y_i \right)
    \]

    By the continuous mapping theorem (or basic properties of convergence in probability), multiplying by the constant $1/k$ yields:

    \[
    \hat{\sigma}^2 \xrightarrow{p} \frac{1}{k} \cdot (k-1)\sigma^2 = \left( 1 - \frac{1}{k} \right)\sigma^2
    \]

    This completes the proof. This result illustrates the famous **Neyman-Scott problem**: the MLE of the common variance is inconsistent because the number of nuisance parameters ($\mu_1, ..., \mu_n$) grows at the same rate as the sample size $n$.

---

## Exercise 6: Consistency and Asymptotic Normality of Sample Quantiles

**Problem:**
(a) Define the check function $\rho_\tau(y) = y(\tau - I(y < 0))$ as the loss function, and let the quantile function of a distribution $F$ be $F^{-1}(p) := \inf\{x \in \mathbb{R} : F(x) \ge p\}$. Prove that:

\[
\theta_0 = \arg \min_\theta E\rho_\tau(X-\theta) = F^{-1}(\tau)
\]

(b) Verify the consistency of the sample quantile defined by the following Z‑estimator:

\[
\Psi_n(\theta) := \frac{1}{n} \sum_{i=1}^n \left( (1-\tau)I(X_i < \theta) - \tau I(X_i > \theta) \right) = 0
\]

Hint: Use Theorem 5.7 or Lemma 5.10 from van der Vaart [1] (add regularity conditions if necessary).

(c) (van der Vaart [1] Ex. 5.11) Assuming consistency, use Theorem 5.23 to prove the asymptotic normality of the sample quantile.

??? success "Solution (click to expand)"

    **(a) Show that the minimizer of the expected loss is the population quantile**

    Assume a regularity condition: the random variable $X$ has a continuous density $f(x)$ and distribution function $F(x)$.
    Define the objective function $M(\theta) = E[\rho_\tau(X-\theta)]$. Using the definition $\rho_\tau(y) = y(\tau - I(y < 0))$, write the expectation as an integral:

    \[
    M(\theta) = \int_{-\infty}^{\infty} (x-\theta)(\tau - I(x-\theta < 0)) f(x) dx
    \]

    Split the integral:

    \[
    M(\theta) = \int_{-\infty}^{\theta} (x-\theta)(\tau - 1) f(x) dx + \int_{\theta}^{\infty} (x-\theta)\tau f(x) dx
    \]

    Differentiate with respect to $\theta$ using Leibniz's rule. Note that when $x=\theta$, the term $(x-\theta)$ vanishes, so boundary terms disappear:

    \[
    M'(\theta) = \int_{-\infty}^{\theta} \frac{\partial}{\partial \theta} [(x-\theta)(\tau-1)] f(x) dx + \int_{\theta}^{\infty} \frac{\partial}{\partial \theta} [(x-\theta)\tau] f(x) dx
    \]

    \[
    = \int_{-\infty}^{\theta} (1-\tau) f(x) dx + \int_{\theta}^{\infty} (-\tau) f(x) dx
    \]

    Computing the integrals:

    \[
    = (1-\tau)F(\theta) - \tau(1-F(\theta)) = F(\theta) - \tau F(\theta) - \tau + \tau F(\theta) = F(\theta) - \tau
    \]

    Set the derivative to zero to find the minimizer:

    \[
    M'(\theta_0) = 0 \implies F(\theta_0) = \tau
    \]

    By the definition of the generalized inverse distribution function (quantile function) $F^{-1}(\tau) = \inf\{x \in \mathbb{R} : F(x) \ge \tau\}$, we conclude:

    \[
    \theta_0 = F^{-1}(\tau)
    \]

    **(b) Verify consistency of the sample quantile**

    We will use **Lemma 5.10** to prove consistency of $\hat{\theta}_n$.

    !!! info "Lemma 5.10 (van der Vaart)"

        Let $\Psi_n(\theta)$ be random functions and $\Psi(\theta)$ a fixed function such that for each $\theta$, $\Psi_n(\theta)$ converges in probability to $\Psi(\theta)$.
        Suppose the map $\theta \mapsto \Psi_n(\theta)$ is **nondecreasing** and $\Psi_n(\hat{\theta}_n) = o_p(1)$.
        If $\theta_0$ satisfies that for every $\epsilon > 0$, $\Psi(\theta_0 - \epsilon) < 0 < \Psi(\theta_0 + \epsilon)$, then $\hat{\theta}_n \xrightarrow{p} \theta_0$.

    We verify the conditions of this lemma one by one:

    * **Condition 1: Monotonicity of $\Psi_n(\theta)$**

        As $\theta$ increases, the indicator $I(X_i < \theta)$ is nondecreasing, while $I(X_i > \theta)$ is nonincreasing (so $-\tau I(X_i > \theta)$ is also nondecreasing). Since $(1-\tau) > 0$ and $\tau > 0$, the map $\theta \mapsto \Psi_n(\theta)$ is a positive linear combination of nondecreasing functions, hence it is **nondecreasing**.

    * **Condition 2: Convergence in probability to a fixed function $\Psi(\theta)$**

        By the Weak Law of Large Numbers, for any fixed $\theta$, $\Psi_n(\theta)$ converges in probability to its expectation.
        Assume $X$ is continuous (so $P(X=\theta)=0$):

        \[
        \Psi(\theta) = E[\Psi_n(\theta)] = (1-\tau)P(X < \theta) - \tau P(X > \theta) = (1-\tau)F(\theta) - \tau(1-F(\theta)) = F(\theta) - \tau
        \]

        Thus $\Psi_n(\theta) \xrightarrow{p} \Psi(\theta)$ holds.

    * **Condition 3: Unique zero crossing (add a regularity condition)**

        The lemma requires that for any $\epsilon > 0$, $\Psi(\theta_0 - \epsilon) < 0 < \Psi(\theta_0 + \epsilon)$.
        For this we add a regularity condition: **assume the distribution function $F(x)$ is strictly increasing in a neighborhood of the true parameter $\theta_0 = F^{-1}(\tau)$**.
        Under this condition, for any $\epsilon > 0$,

        \[
        F(\theta_0 - \epsilon) < F(\theta_0) = \tau < F(\theta_0 + \epsilon)
        \]

        Subtracting $\tau$ gives exactly:

        \[
        \Psi(\theta_0 - \epsilon) < 0 < \Psi(\theta_0 + \epsilon)
        \]

        Since $\theta \mapsto \Psi_n(\theta)$ is nondecreasing, the sample quantile satisfies $\Psi_n(\hat{\theta}_n) = 0$ (or $o_p(1)$), and the true parameter $\theta_0$ uniquely crosses the zero of the limiting function $\Psi(\theta)$, all conditions of Lemma 5.10 are satisfied. Hence the sample quantile is consistent:

        \[
        \hat{\theta}_n \xrightarrow{p} \theta_0
        \]

    **(c) Asymptotic normality of the sample quantile**

    To prove asymptotic normality, we apply **Theorem 5.23**. To match the notation of the textbook, we denote the quantile as $p$ (which corresponds to $\tau$ in the previous part).

    !!! info "Theorem 5.23 (Asymptotic normality of M‑estimators)"

        Suppose the map $\theta \mapsto m_\theta(x)$ is differentiable at $\theta_0$ for $P$-almost every $x$, with derivative $\dot{m}_{\theta_0}(x)$. Moreover, assume a Lipschitz condition:
        $|m_{\theta_1}(x) - m_{\theta_2}(x)| \le \dot{m}(x)|\theta_1 - \theta_2|$, where the envelope satisfies $P\dot{m}^2 < \infty$.
        Also assume that $\theta \mapsto P m_\theta$ admits a second‑order Taylor expansion at the maximum $\theta_0$ with a non‑singular symmetric second derivative matrix $V_{\theta_0}$.
        If $\hat{\theta}_n \xrightarrow{p} \theta_0$ and approximately maximizes the empirical objective function, then:

        $\sqrt{n}(\hat{\theta}_n - \theta_0) = -V_{\theta_0}^{-1} \frac{1}{\sqrt{n}} \sum_{i=1}^n \dot{m}_{\theta_0}(X_i) + o_p(1)$.

    The sample $p$-quantile $\hat{\theta}_n = \mathbb{F}_n^{-1}(p)$ minimizes the empirical objective function $\frac{1}{n}\sum \rho_p(X_i - \theta)$. To fit the "maximization" framework of Theorem 5.23, define the objective function as:

    \[
    m_\theta(x) = -\rho_p(x-\theta) = -(x-\theta)(p - I(x < \theta))
    \]

    We verify the conditions of Theorem 5.23 at the true parameter $\theta_0 = F^{-1}(p)$:

    * **Condition 1: $P$-almost everywhere differentiability**

        The map $\theta \mapsto m_\theta(x)$ is differentiable everywhere except at $\theta = x$. Because the distribution $F$ is differentiable at $\theta_0$ with density $f(\theta_0) > 0$, the probability of an observation exactly equal to $\theta_0$ is zero, i.e., $P(X=\theta_0)=0$. Hence $m_\theta(x)$ is differentiable at $\theta_0$ for $P$-almost every $x$, with derivative:

        \[
        \dot{m}_{\theta_0}(x) = \frac{\partial}{\partial \theta} [-(x-\theta)p + (x-\theta)I(x < \theta)] \Big|_{\theta=\theta_0} = p - I(x \le \theta_0)
        \]

    * **Condition 2: Lipschitz continuity and square‑integrable envelope**

        The check function $\rho_p(u)$ is Lipschitz continuous with constant bounded by $\max(p, 1-p) \le 1$.
        Therefore, for any $\theta_1, \theta_2$:

        \[
        |m_{\theta_1}(x) - m_{\theta_2}(x)| = |\rho_p(x-\theta_1) - \rho_p(x-\theta_2)| \le 1 \cdot |\theta_1 - \theta_2|
        \]

        We can choose the envelope $\dot{m}(x) = 1$. The condition $P\dot{m}^2 = E[1^2] = 1 < \infty$ is clearly satisfied.

    * **Condition 3: Second‑order Taylor expansion and non‑singular $V_{\theta_0}$**

        The expected objective function is $Pm_\theta = -E[\rho_p(X-\theta)]$. From part (a), its first derivative is:

        \[
        \frac{d}{d\theta} (Pm_\theta) = p - F(\theta)
        \]

        Differentiating again, we obtain the second derivative at $\theta_0$:

        \[
        V_{\theta_0} = \frac{d}{d\theta} (p - F(\theta)) \Big|_{\theta=\theta_0} = -f(\theta_0)
        \]

        The problem assumes that $F$ has a positive derivative at the population $p$-quantile $F^{-1}(p)$, i.e., $f(\theta_0) > 0$. Therefore, the scalar $V_{\theta_0} = -f(\theta_0)$ is strictly negative, satisfying the "non‑singular" condition.

    * **Applying Theorem 5.23**

        We have already established consistency $\hat{\theta}_n \xrightarrow{p} \theta_0$ in part (b). All conditions are met. Applying Theorem 5.23:

        \[
        \sqrt{n}(\hat{\theta}_n - \theta_0) = -V_{\theta_0}^{-1} \frac{1}{\sqrt{n}} \sum_{i=1}^n \dot{m}_{\theta_0}(X_i) + o_p(1)
        \]

        Substituting $V_{\theta_0}$ and $\dot{m}_{\theta_0}(X_i)$:

        \[
        \sqrt{n}(\hat{\theta}_n - \theta_0) = \frac{1}{f(\theta_0)} \frac{1}{\sqrt{n}} \sum_{i=1}^n (p - I(X_i \le \theta_0)) + o_p(1)
        \]

        Let $Y_i = p - I(X_i \le \theta_0)$. The sequence $\{Y_i\}_{i=1}^n$ is i.i.d. with:
        $E[Y_i] = p - P(X_i \le \theta_0) = p - F(\theta_0) = p - p = 0$
        $Var(Y_i) = Var(I(X_i \le \theta_0)) = p(1-p)$

        By the classical CLT:

        \[
        \frac{1}{\sqrt{n}} \sum_{i=1}^n Y_i \xrightarrow{d} \mathcal{N}(0, p(1-p))
        \]

        Finally, by Slutsky’s theorem, multiplying by the constant scaling $1/f(\theta_0)$ gives the asymptotic normal distribution:

        \[
        \sqrt{n}(\hat{\theta}_n - \theta_0) \xrightarrow{d} \mathcal{N}\left(0, \frac{p(1-p)}{f(\theta_0)^2}\right)
        \]