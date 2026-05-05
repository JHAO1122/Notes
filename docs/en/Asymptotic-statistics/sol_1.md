---
tags:
  - Large Sample Theory
  - Statistics
  - Homework Solutions
---

# 📈 Exercise Solutions: Large Sample Theory Assignment 1

!!! abstract "About this page"
    This page contains solutions for the core exercises of the first unit of the "Large Sample Theory" course (Stochastic Convergence, Asymptotic Normality). All solutions are organized based on standard course answers and use collapsible callouts to display detailed proofs and construction processes.

---

## Exercise 1: Stochastic Convergence

**Problem:**
(a) Find an example of a sequence of random variables such that $X_{n} \xrightarrow{d} 0$ but $E[X_{n}] \rightarrow \infty$.
(b) Find an example of sequences such that $X_{n} \rightarrow X$ and $Y_{n} \rightarrow Y$, but the joint sequence $(X_{n}, Y_{n})$ does not converge in law.

??? success "Solution (Click to expand)"

    **(a) Constructing a sequence with diverging expectation but converging in distribution to 0**
    
    We construct the following discrete random variable $X_n$:
    
    $$
    P(X_{n}=n^{2}) = \frac{1}{n}, \quad P(X_{n}=0) = 1 - \frac{1}{n}
    $$

    * **Convergence in Probability**: For any $\epsilon > 0$, when $n$ is sufficiently large:
        $P(|X_{n}| > \epsilon) = P(X_{n} = n^{2}) = \frac{1}{n} \xrightarrow{n \rightarrow \infty} 0$.
        Since convergence in probability implies convergence in distribution, $X_{n} \xrightarrow{d} 0$.
    * **Limit of Expectation**:
        $E[X_{n}] = n^{2} \cdot \frac{1}{n} + 0 \cdot (1 - \frac{1}{n}) = n \xrightarrow{n \rightarrow \infty} \infty$.
    
    Proved.

    <br>

    **(b) Constructing an example where the joint distribution does not converge**
    
    Let $X \sim \mathcal{N}(0,1)$. Define sequences $X_{n}$ and $Y_{n}$ as follows:
    
    $$
    X_{n} = X \quad (\forall n), \quad Y_{n} = \begin{cases} X, & \text{if } n \text{ is odd} \\ -X, & \text{if } n \text{ is even} \end{cases}
    $$
    
    * **Marginal Convergence**: Clearly $X_{n} \xrightarrow{d} \mathcal{N}(0,1)$. Due to the symmetry of the standard normal distribution, $-X \sim \mathcal{N}(0,1)$, thus each component of $Y_{n}$ also follows $\mathcal{N}(0,1)$, meaning $Y_{n} \xrightarrow{d} \mathcal{N}(0,1)$.
    * **Non-convergence of Joint Distribution**: The joint distribution $(X_{n}, Y_{n})$ alternates between $(X, X)$ and $(X, -X)$. Since these two joint distributions are clearly different, the joint sequence does not converge.

---

## Exercise 2: Sums of Convergent Sequences

**Problem:**
Does the convergence $X_{n} \rightarrow X$ and $Y_{n} \rightarrow Y$ always imply that $X_{n} + Y_{n} \rightarrow X + Y$? Here "$\rightarrow$" can represent convergence in distribution, convergence in probability, or almost sure convergence. Please prove the claim or find a counter-example for at least one of these scenarios.

??? success "Solution (Click to expand)"

    **Conclusion: The claim does not hold for Convergence in Distribution.**

    **Counter-example Construction:**
    Let $X$ and $Y$ be independent standard normal random variables, $X, Y \sim i.i.d. \mathcal{N}(0,1)$.
    Define the sequences:
    
    $$
    X_{n} = X \quad \text{and} \quad Y_{n} = -X \quad (\forall n)
    $$
    
    * **Marginal Convergence**: Clearly $X_{n} \xrightarrow{d} X$. Since $-X \sim \mathcal{N}(0,1)$, $Y_{n} \xrightarrow{d} Y$.
    * **Convergence of the Sum**: Under this construction, $X_{n} + Y_{n} \equiv 0$, therefore $X_{n} + Y_{n} \xrightarrow{d} 0$.
    * **Distribution of the Target Limit**: However, $X + Y$ is the sum of two independent normal variables, following a $\mathcal{N}(0, 2)$ distribution, which is clearly not a degenerate distribution at 0.
    
    Since the distribution at 0 is not equal to $\mathcal{N}(0, 2)$, the claim fails for convergence in distribution.

---

## Exercise 3: Convergence Properties and Cauchy Distribution

**Problem:**
(a) Let $X_{1}, ..., X_{n}$ be a sequence of random variables. Prove that $X_{n} \xrightarrow{p} 0$ if and only if $E\left(\frac{|X_{n}|}{1+|X_{n}|}\right) \rightarrow 0$.
(b) Let $X_{1}, ..., X_{n}$ be i.i.d. Cauchy $(0,1)$ random variables. Prove that $\bar{X}_{n} \rightarrow Cauchy(0,1)$.
(c) Let $X_{n}, 1 \le n \le \infty$ be integer-valued random variables. Prove that $X_{n} \xrightarrow{d} X_{\infty}$ if and only if $P(X_{n}=m) \rightarrow P(X_{\infty}=m)$ for all $m$.

??? success "Solution (Click to expand)"

    **(a) Proof of equivalent characterization of convergence in probability**
    Let $f(x) = \frac{|x|}{1+|x|}$. Note that $f$ is a continuous, bounded function ($0 \le f(x) < 1$) that is strictly increasing on $[0, \infty)$.
    
    * **Necessity ($\Rightarrow$)**: If $X_{n} \xrightarrow{p} 0$, then by the Continuous Mapping Theorem, $f(X_{n}) \xrightarrow{p} 0$. Since $f(X_{n})$ is bounded by 1, applying the Bounded Convergence Theorem (BCT) yields $E[f(X_{n})] \rightarrow 0$.
    * **Sufficiency ($\Leftarrow$)**: For any $\epsilon > 0$, by Markov's Inequality:

        $$P(|X_{n}| > \epsilon) = P(f(X_{n}) > f(\epsilon)) \le \frac{E[f(X_{n})]}{f(\epsilon)}$$
        
        Since $f(\epsilon) > 0$ is a constant and $E[f(X_{n})] \rightarrow 0$, then $P(|X_{n}| > \epsilon) \rightarrow 0$, meaning $X_{n} \xrightarrow{p} 0$.

    <br>

    **(b) Properties of the sample mean of Cauchy distributions**
    The characteristic function of $Cauchy(0,1)$ is $\varphi_{X}(t) = \exp(-|t|)$.
    By independence, the characteristic function of the sample mean $\bar{X}_{n} = \frac{1}{n} \sum X_{k}$ is:
    
    $$
    \varphi_{\bar{X}_{n}}(t) = \prod_{k=1}^{n} \varphi_{X}\left(\frac{t}{n}\right) = \left(\exp\left(-\left|\frac{t}{n}\right|\right)\right)^{n} = \exp(-|t|)
    $$
    
    Since the characteristic function corresponds uniquely to the distribution and $\varphi_{\bar{X}_{n}}(t) = \varphi_{X}(t)$ for all $n$, $\bar{X}_{n}$ follows $Cauchy(0,1)$ for all $n$. By Lévy's Continuity Theorem, $\bar{X}_{n} \xrightarrow{d} Cauchy(0,1)$.

    <br>

    **(c) Convergence in distribution for integer-valued random variables**
    
    * **Necessity ($\Rightarrow$)**: Define $f_{m}(x) = \max(0, 1 - |x - m|)$, which is a bounded continuous function. Since $X_n$ is integer-valued, $E[f_{m}(X_{n})] = P(X_{n}=m)$.
        By the Portmanteau Lemma, $X_{n} \xrightarrow{d} X_{\infty}$ implies $E[f(X_n)] \rightarrow E[f(X)]$ for all bounded continuous functions $f$, thus $P(X_{n}=m) \rightarrow P(X_{\infty}=m)$.
    * **Sufficiency ($\Leftarrow$)**: The characteristic function is $\varphi_{X_{n}}(t) = \sum_{m \in \mathbb{Z}} e^{itm} P(X_{n}=m)$. Since the series is bounded by 1 and $\sum P(X_{\infty}=m)=1$, we apply the Dominated Convergence Theorem (DCT) for series:

        $$\lim_{n \rightarrow \infty} \varphi_{X_{n}}(t) = \sum_{m \in \mathbb{Z}} e^{itm} \lim_{n \rightarrow \infty} P(X_{n}=m) = \varphi_{X_{\infty}}(t)$$
        
    By Lévy's Continuity Theorem, the result holds.

---

## Exercise 4: Asymptotic Normality

**Problem:**
If $X_{n}$ is $AN(\mu_{n}, \sigma_{n}^{2})$ and $Y_{n} \sim \mathcal{N}(\mu_{n}, \sigma_{n}^{2})$, prove:
(a) $\sup_{t \in \mathbb{R}} |P(X_{n} \le t) - P(Y_{n} \le t)| \rightarrow 0$ as $n \rightarrow \infty$.
(b) $X_{n}$ is $AN(\bar{\mu}_{n}, \bar{\sigma}_{n}^{2})$ if and only if $\frac{\bar{\sigma}_{n}}{\sigma_{n}} \rightarrow 1$ and $\frac{\mu_{n} - \bar{\mu}_{n}}{\sigma_{n}} \rightarrow 0$.
(c) $a_{n}X_{n} + b_{n}$ is $AN(\mu_{n}, \sigma_{n}^{2})$ if and only if $a_{n} \rightarrow 1$ and $\frac{\mu_{n}(a_{n}-1) + b_{n}}{\sigma_{n}} \rightarrow 0$.

??? success "Solution (Click to expand)"

    **(a) Uniform convergence of Cumulative Distribution Functions**
    Let $Z_{n} = \frac{X_{n} - \mu_{n}}{\sigma_{n}}$. By the definition of asymptotic normality, $Z_{n} \xrightarrow{d} \mathcal{N}(0,1)$.
    Let $\Phi(t)$ be the CDF of the standard normal distribution.
    
    $$
    P(X_{n} \le t) = P\left(Z_{n} \le \frac{t - \mu_{n}}{\sigma_{n}}\right), \quad P(Y_{n} \le t) = \Phi\left(\frac{t - \mu_{n}}{\sigma_{n}}\right)
    $$
    
    Since $\Phi(\cdot)$ is continuous everywhere, **Pólya's Theorem** guarantees uniform convergence of the CDF:
    
    $$\sup_{z \in \mathbb{R}} |P\left(Z_{n} \le z\right) - \Phi(z)| \rightarrow 0$$
    
    Substituting $z = \frac{t - \mu_{n}}{\sigma_{n}}$ completes the proof.

    <br>

    **(b) Conditions for parameter substitution**
    Define $W_{n} = \frac{X_{n} - \bar{\mu}_{n}}{\bar{\sigma}_{n}} = \frac{\sigma_{n}}{\bar{\sigma}_{n}} Z_{n} + \frac{\mu_{n} - \bar{\mu}_{n}}{\bar{\sigma}_{n}}$.
    
    * **Sufficiency ($\Leftarrow$)**: If the limit conditions hold, by Slutsky's Theorem, $W_{n} \xrightarrow{d} 1 \cdot Z + 0 \sim \mathcal{N}(0,1)$.
    * **Necessity ($\Rightarrow$)**: If $W_{n} \xrightarrow{d} \mathcal{N}(0,1)$ and we know $Z_{n} \xrightarrow{d} \mathcal{N}(0,1)$, according to the **Convergence of Types Theorem**, the scale multiplier ratio must tend to 1 and the location shift ratio must tend to 0.

    <br>

    **(c) Asymptotic normality of linear transformations**
    Let $T_{n} = \frac{(a_{n}X_{n} + b_{n}) - \mu_{n}}{\sigma_{n}} = a_{n} Z_{n} + \frac{\mu_{n}(a_{n}-1) + b_{n}}{\sigma_{n}}$.
    Similarly, using Slutsky's Theorem and the Convergence of Types Theorem, $T_{n} \xrightarrow{d} \mathcal{N}(0,1)$ if and only if $a_{n} \rightarrow 1$ and the shift term $\frac{\mu_{n}(a_{n}-1) + b_{n}}{\sigma_{n}} \rightarrow 0$.

---

## Exercise 5: Portmanteau Lemma

**Problem:**
Prove the equivalence of (ii) and (iv) in the Portmanteau Lemma:
(ii) $E[f(X_{n})] \rightarrow E[f(X)]$ for all bounded continuous functions $f$.
(iv) $\liminf_{n \rightarrow \infty} E[f(X_{n})] \ge E[f(X)]$ for all non-negative continuous functions $f$.

??? success "Solution (Click to expand)"

    **Proof:**

    **1. Prove (ii) $\Rightarrow$ (iv):**
    Let $f \ge 0$ be any non-negative continuous function.
    For any constant $M > 0$, define the truncation function $f_{M}(x) = \min(f(x), M)$.
    Clearly $f_{M}(x)$ is bounded and continuous. Since $f(X_{n}) \ge f_{M}(X_{n})$, by the monotonicity of expectation:
    
    $$E[f(X_{n})] \ge E[f_{M}(X_{n})]$$
    
    Taking the limit inferior on both sides and applying condition (ii):
    
    $$\liminf_{n \rightarrow \infty} E[f(X_{n})] \ge \liminf_{n \rightarrow \infty} E[f_{M}(X_{n})] = E[f_{M}(X)]$$
    
    As $M \rightarrow \infty$, $f_{M}(X) \uparrow f(X)$. By the Monotone Convergence Theorem (MCT):
    
    $$E[f_{M}(X)] \rightarrow E[f(X)]$$
    
    Therefore, $\liminf_{n \rightarrow \infty} E[f(X_{n})] \ge E[f(X)]$.

    <br>

    **2. Prove (iv) $\Rightarrow$ (ii):**
    Let $f$ be a bounded continuous function such that $|f(x)| \le M$ ($M > 0$).
    Then $f(x) + M \ge 0$ and $M - f(x) \ge 0$ are both non-negative continuous functions.
    
    * Apply (iv) to $f(x) + M$:
    
        $$\liminf_{n \rightarrow \infty} E[f(X_{n}) + M] \ge E[f(X) + M] \implies \liminf_{n \rightarrow \infty} E[f(X_{n})] \ge E[f(X)] \quad (*)$$
    
    * Apply (iv) to $M - f(x)$:
    
        $$\liminf_{n \rightarrow \infty} E[M - f(X_{n})] \ge E[M - f(X)] \implies M - \limsup_{n \rightarrow \infty} E[f(X_{n})] \ge M - E[f(X)]$$
        
        This is equivalent to:
    
        $$\limsup_{n \rightarrow \infty} E[f(X_{n})] \le E[f(X)] \quad (**)$$
    
    Combining $(*)$ and $(**)$, we get:
    
    $$\limsup_{n \rightarrow \infty} E[f(X_{n})] \le E[f(X)] \le \liminf_{n \rightarrow \infty} E[f(X_{n})]$$
    
    This shows the limit exists and $\lim_{n \rightarrow \infty} E[f(X_{n})] = E[f(X)]$.

---

## Exercise 6: Basic rules of stochastic order

**Problem:**
Prove the following basic rules regarding stochastic order:
(a) $O_{p}(1)o_{p}(1) = o_{p}(1)$
(b) $(1 + o_{p}(1))^{-1} = O_{p}(1)$
(c) $o_{p}(O_{p}(1)) = o_{p}(1)$

??? success "Solution (Click to expand)"

    **(a) Prove $O_{p}(1)o_{p}(1) = o_{p}(1)$**
    Let $X_{n} = O_{p}(1)$ and $Y_{n} = o_{p}(1)$. We want to prove that for any $\epsilon > 0$, $P(|X_{n}Y_{n}| > \epsilon) \rightarrow 0$.
    For any fixed constant $M > 0$, the following inequality holds:
    
    $$P(|X_{n}Y_{n}| > \epsilon) \le P(|X_{n}| > M) + P(|Y_{n}| > \frac{\epsilon}{M})$$
    
    Given $\eta > 0$:
    1. Since $X_{n} = O_{p}(1)$, there exists $M > 0$ and $N_{1}$ such that for $n > N_{1}$, $P(|X_{n}| > M) < \frac{\eta}{2}$.
    2. Fix this $M$. Since $Y_{n} = o_{p}(1)$, for $\frac{\epsilon}{M} > 0$, there exists $N_{2}$ such that for $n > N_{2}$, $P(|Y_{n}| > \frac{\epsilon}{M}) < \frac{\eta}{2}$.
    
    Let $N = \max(N_{1}, N_{2})$, then for $n > N$:
    
    $$P(|X_{n}Y_{n}| > \epsilon) < \frac{\eta}{2} + \frac{\eta}{2} = \eta$$
    
    Thus $X_{n}Y_{n} = o_{p}(1)$.

    <br>

    **(b) Prove $(1 + o_{p}(1))^{-1} = O_{p}(1)$**
    Let $X_{n} = o_{p}(1)$. We need to prove that for any $\epsilon > 0$, there exists $M > 0$ and $N$ such that for $n > N$, $P(|(1 + X_{n})^{-1}| > M) < \epsilon$.
    Consider $M > 1$, the event $\{|(1 + X_{n})^{-1}| > M\}$ is equivalent to $\{|1 + X_{n}| < \frac{1}{M}\}$.
    By the reverse triangle inequality $|1 + X_{n}| \ge 1 - |X_{n}|$, if $|1 + X_{n}| < \frac{1}{M}$, then $1 - |X_{n}| < \frac{1}{M}$, meaning $|X_{n}| > 1 - \frac{1}{M}$. Therefore:
    
    $$P\left(\left|\frac{1}{1 + X_{n}}\right| > M\right) \le P\left(|X_{n}| > 1 - \frac{1}{M}\right)$$
    
    Since $M > 1$, then $1 - \frac{1}{M} > 0$. Because $X_{n} = o_{p}(1)$, the probability tends to 0 as $n \rightarrow \infty$. For a given $\epsilon$, we can always find a sufficiently large $n$ to satisfy the requirement, thus the conclusion holds.

    <br>

    **(c) Prove $o_{p}(O_{p}(1)) = o_{p}(1)$**
    Let $X_{n} = O_{p}(1)$ and $W_{n} = o_{p}(X_{n})$.
    By the definition of the little $o_p$ notation, $W_{n}$ can be expressed as $W_{n} = Y_{n}X_{n}$ where $Y_{n} = o_{p}(1)$.
    Thus we have:
    
    $$o_{p}(X_{n}) = Y_{n}X_{n} = o_{p}(1) \cdot O_{p}(1)$$
    
    According to the result proved in part (a) of this problem, $o_{p}(1)O_{p}(1) = o_{p}(1)$.
    Therefore $W_{n} = o_{p}(1)$.

---

## Exercise 7: Characteristic functions

**Problem:**
Let $\phi$ be a characteristic function (ch.f.) on $\mathcal{R}^{k}$:
(a) Prove that $|\phi| \le 1$ and is uniformly continuous on $\mathcal{R}^{k}$.
(b) Find an example of two random variables $X$ and $Y$ such that $X, Y$ are not independent, but their characteristic functions satisfy the multiplication formula $\phi_{X}(t)\phi_{Y}(t) = \phi_{X+Y}(t)$ for all $t \in \mathcal{R}$.

??? success "Solution (Click to expand)"

    **(a) Proof of $|\phi| \le 1$ and uniform continuity**
    
    * **Boundedness**:
       By definition $\phi(t) = E[\exp(it^{T}X)]$. Using the properties of expectation:
       
       $|\phi(t)| = |E[\exp(it^{T}X)]| \le E[|\exp(it^{T}X)|] = E[1] = 1$
    

    * **Uniform Continuity**:
       For any $t, h \in \mathcal{R}^{k}$, consider the difference:
       
       $|\phi(t + h) - \phi(t)| = |E[\exp(i(t + h)^{T}X) - \exp(it^{T}X)]|$
       
       $\le E[|\exp(it^{T}X)(\exp(ih^{T}X) - 1)|] = E[|\exp(ih^{T}X) - 1|]$
       
       Notice that the upper bound $E[|\exp(ih^{T}X) - 1|]$ is independent of $t$.
       Let the random variable $Z_{h} = |\exp(ih^{T}X) - 1|$. Clearly $Z_{h} \le 2$ and as $h \rightarrow 0$, $Z_{h} \xrightarrow{a.s.} 0$.
       By the Lebesgue Dominated Convergence Theorem (LDCT):
       
       $\lim_{h \rightarrow 0} E[|\exp(ih^{T}X) - 1|] = E[0] = 0$
       
       Since this convergence is uniform for all $t$, $\phi(t)$ is uniformly continuous on $\mathcal{R}^{k}$.

    <br>

    **(b) Constructing an example that is not independent but satisfies the ch.f. multiplication formula**
    
    Let $X \sim Cauchy(0,1)$ and let $Y = X$.
    Clearly $X$ and $Y$ are perfectly dependent, hence not independent.
    
    * The ch.f. of the standard Cauchy distribution is $\phi_{X}(t) = \exp(-|t|)$.
    * Since $Y = X$, then $\phi_{Y}(t) = \exp(-|t|)$ as well.
    * The product of their ch.f.s is:
        
        $$\phi_{X}(t)\phi_{Y}(t) = \exp(-|t|) \cdot \exp(-|t|) = \exp(-2|t|)$$
    
    Consider their sum $X + Y = 2X$. According to the properties of ch.f.s:
    
    $$\phi_{X+Y}(t) = E[\exp(it(2X))] = E[\exp(i(2t)X)] = \phi_{X}(2t)$$
    
    Substituting the form of the Cauchy ch.f.:
    
    $$\phi_{X}(2t) = \exp(-|2t|) = \exp(-2|t|)$$
    
    We observe that $\phi_{X}(t)\phi_{Y}(t) = \exp(-2|t|) = \phi_{X+Y}(t)$ holds for all $t \in \mathcal{R}$.
    This counter-example proves that the multiplicative property of characteristic functions is not a sufficient condition for the independence of random variables.