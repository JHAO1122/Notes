# 🎲 Probability Theory

This module covers core inequalities in probability theory, convergence theory of random variables, and large sample limit theorems.

## 1. Core Probability Inequalities

!!! info "Markov's Inequality"
    For a non-negative random variable $X \ge 0$ and any $a > 0$:
    
    $$
    P(X \ge a) \le \frac{E[X]}{a}
    $$

!!! info "Chebyshev's Inequality"
    For any random variable $X$ and $a > 0$:

    $$
    P(|X - E[X]| \ge a) \le \frac{Var(X)}{a^2}
    $$

!!! success "Kolmogorov's Maximal Inequality"
    If $X_1, \dots, X_n$ are independent and identically distributed (i.i.d.) with mean 0 and finite variance:
    
    $$
    P\left(\max_{1\le k\le n} \left|\sum_{i=1}^k X_i\right| \ge \lambda\right) \le \frac{1}{\lambda^2} \sum_{i=1}^n Var(X_i)
    $$

??? note "Borel-Cantelli Lemma"
    Let the event sequence $\{A_n\}$ occurring infinitely often (i.o.) be denoted as $\limsup A_n$:
    
    1. **First Lemma**: If $\sum_{n=1}^\infty P(A_n) < \infty$, then $P(A_n \text{ i.o.}) = 0$.
    2. **Second Lemma**: If $\{A_n\}$ are **mutually independent**, and $\sum_{n=1}^\infty P(A_n) = \infty$, then $P(A_n \text{ i.o.}) = 1$.

---

## 2. Convergence of Random Variables

!!! abstract "Definitions of Four Modes of Convergence"
    1. **Convergence in Probability ($X_n \xrightarrow{P} X$)**: For any $\epsilon > 0$, $\lim_{n \to \infty} P(|X_n - X| > \epsilon) = 0$.
    2. **Almost Sure Convergence ($X_n \xrightarrow{a.s.} X$)**: $P(\lim_{n \to \infty} X_n = X) = 1$.
    3. **Convergence in $L^p$ ($X_n \xrightarrow{L^p} X$)**: $\lim_{n \to \infty} E[|X_n - X|^p] = 0$.
    4. **Convergence in Distribution ($X_n \xrightarrow{d} X$)**:
        * **Definition**: For all continuity points $x$ of the limit distribution function $F(x)$, $\lim_{n \to \infty} F_n(x) = F(x)$.
        * **Lévy's Continuity Theorem**: $X_n \xrightarrow{d} X$ if and only if its characteristic function $\phi_{X_n}(t) \to \phi_X(t)$ holds for every $t$, and $\phi_X(t)$ is continuous at $t=0$.

??? tip "Relationships Between Modes of Convergence (Important!)"
    * **Implication Chains**: $L^p \implies P \implies d$ and $a.s. \implies P \implies d$.
    * **Reverse Implications (Conditional)**:
        * **$P \to L^p$**: If $X_n \xrightarrow{P} X$ and $\{|X_n|^p\}$ is **Uniformly Integrable (UI)**, then $X_n \xrightarrow{L^p} X$.
        * **$P \to a.s.$**: If $X_n \xrightarrow{P} X$, then there necessarily exists a **subsequence** $\{n_k\}$ such that $X_{n_k} \xrightarrow{a.s.} X$.
        * **$d \to P$**: If the limit distribution $X=c$ is a degenerate constant, then convergence in distribution implies convergence in probability.
    * **Slutsky's Theorem**: If $X_n \xrightarrow{d} X$ and $Y_n \xrightarrow{P} c$ (a constant), then $X_n Y_n \xrightarrow{d} cX$ and $X_n + Y_n \xrightarrow{d} X + c$.

---

## 3. Law of Large Numbers and Central Limit Theorem 

!!! info "Khinchin's Weak Law of Large Numbers (WLLN)"
    Let $X_1, X_2, \dots, X_n$ be a sequence of **independent and identically distributed (i.i.d.)** random variables with finite expectation $E[X_i] = \mu$ (i.e., $E[|X_1|] < \infty$).
    Then the sample mean **converges in probability** to $\mu$:

    $$
    \bar{X}_n = \frac{1}{n}\sum_{i=1}^n X_i \xrightarrow{P} \mu \quad (n \to \infty)
    $$

!!! info "Kolmogorov's Strong Law of Large Numbers (SLLN)"
    If $X_i$ are i.i.d., then $E[|X_1|] < \infty$ is a necessary and sufficient condition for $\bar{X}_n \xrightarrow{a.s.} E[X]$.

??? success "Central Limit Theorem (CLT) and Conditions"
    For a sequence of independent but not necessarily identically distributed random variables $X_1, X_2, \dots$, let $E[X_i]=\mu_i$, $Var(X_i)=\sigma_i^2$, and $S_n^2 = \sum_{i=1}^n \sigma_i^2$.
    The standardized sum is $Z_n = \frac{1}{S_n} \sum_{i=1}^n (X_i - \mu_i)$.
    
    * **Lindeberg's Condition (Necessary and Sufficient)**: For any $\epsilon > 0$, it satisfies:

    $$
    \lim_{n \to \infty} \frac{1}{S_n^2} \sum_{i=1}^n E[(X_i - \mu_i)^2 I(|X_i - \mu_i| > \epsilon S_n)] = 0
    $$
    
    *(Intuition: Ensures that no single random variable's extreme tail values dominate the total variance)*

    * **Lyapunov's Condition (Sufficient)**: If there exists a $\delta > 0$ such that:

    $$
    \lim_{n \to \infty} \frac{1}{S_n^{2+\delta}} \sum_{i=1}^n E[|X_i - \mu_i|^{2+\delta}] = 0
    $$
    
    *(Note: Lyapunov's condition is often easier to verify because if a higher-order moment exists and satisfies this limit, simple integral bounding directly implies Lindeberg's condition holds, thus $Z_n \xrightarrow{d} N(0,1)$)*

---

## 4. Quick Reference Table for Common Distributions 

| Distribution Name | Notation | Expectation $E[X]$ | Variance $Var(X)$ | PDF/PMF | Characteristic Function $\phi_X(t)$ | Moment Generating Function $M_X(t)$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Normal Distribution** | $N(\mu, \sigma^2)$ | $\mu$ | $\sigma^2$ | $\frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{(x-\mu)^2}{2\sigma^2}}$ | $e^{it\mu - \frac{1}{2}\sigma^2 t^2}$ | $e^{t\mu + \frac{1}{2}\sigma^2 t^2}$ |
| **Exponential Distribution** | $Exp(\lambda)$ | $\frac{1}{\lambda}$ | $\frac{1}{\lambda^2}$ | $\lambda e^{-\lambda x} \quad (x \ge 0)$ | $\frac{\lambda}{\lambda - it}$ | $\frac{\lambda}{\lambda - t} \quad (t < \lambda)$ |
| **Uniform Distribution** | $U(a, b)$ | $\frac{a+b}{2}$ | $\frac{(b-a)^2}{12}$ | $\frac{1}{b-a} \quad (a \le x \le b)$ | $\frac{e^{itb} - e^{ita}}{it(b-a)}$ | $\frac{e^{tb} - e^{ta}}{t(b-a)}$ |
| **Poisson Distribution** | $Pois(\lambda)$ | $\lambda$ | $\lambda$ | $e^{-\lambda} \frac{\lambda^k}{k!} \quad (k \in \mathbb{N})$ | $e^{\lambda(e^{it} - 1)}$ | $e^{\lambda(e^t - 1)}$ |
| **Binomial Distribution** | $Bin(n, p)$ | $np$ | $np(1-p)$ | $\binom{n}{k} p^k (1-p)^{n-k}$ | $(1 - p + pe^{it})^n$ | $(1 - p + pe^t)^n$ |

---

## V. Change of Variables & PDF Transformation

!!! info "Univariate Transformation"
    Let $X$ have PDF $f_X(x)$, and let $Y = g(X)$. If $g(x)$ is a strictly monotonic and differentiable function with inverse $x = g^{-1}(y)$, the PDF of $Y$ is:
    
    $$
    f_Y(y) = f_X(g^{-1}(y)) \cdot \left| \frac{d}{dy} g^{-1}(y) \right|
    $$

!!! success "Multivariate Transformation & the Jacobian"
    Let $\mathbf{X} = (X_1, \dots, X_n)^T$ have a joint PDF $f_{\mathbf{X}}(\mathbf{x})$, and let $\mathbf{Y} = \mathbf{g}(\mathbf{X})$ be a one-to-one mapping.
    Let the inverse be $\mathbf{x} = \mathbf{h}(\mathbf{y})$ where $\mathbf{h} = \mathbf{g}^{-1}$. The joint PDF of $\mathbf{Y}$ is:
    
    $$
    f_{\mathbf{Y}}(\mathbf{y}) = f_{\mathbf{X}}(\mathbf{h}(\mathbf{y})) \cdot | \det(J_{\mathbf{h}}(\mathbf{y})) |
    $$
    
    Where $J_{\mathbf{h}}(\mathbf{y})$ is the **Jacobian Matrix**, defined as:
    
    $$
    J_{\mathbf{h}}(\mathbf{y}) = \frac{\partial (x_1, \dots, x_n)}{\partial (y_1, \dots, y_n)} = \begin{pmatrix} 
    \frac{\partial x_1}{\partial y_1} & \dots & \frac{\partial x_1}{\partial y_n} \\ 
    \vdots & \ddots & \vdots \\ 
    \frac{\partial x_n}{\partial y_1} & \dots & \frac{\partial x_n}{\partial y_n} 
    \end{pmatrix}
    $$