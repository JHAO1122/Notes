---
tags:
  - Large Sample Theory
  - Statistics
  - Homework Solutions
---

# 📈 Detailed Solutions to Homework: Large Sample Theory Assignment 2

## Problem 1: Moments and expansion of CFs

**Problem:**
If $E|X|^{r}<\infty$, show that the characteristic function $\phi_{X}(t)$ has the following expansion:

\[
\phi_{X}(t)=\sum_{j=0}^{r}\frac{(it)^{j}}{j!}EX^{j}+o(|t|^{r})
\]

??? success "Solution (click to expand)"

    **Lemma: Bound on the remainder of the complex exponential Taylor expansion**
    
    Let $R_{r}(x)=e^{ix}-\sum_{j=0}^{r}\frac{(ix)^{j}}{j!}$ be the remainder of the $r$-th order Taylor expansion of $e^{ix}$. For any real $x$ and integer $r\ge0$, we have:

    \[
    |R_{r}(x)|\le \min \left( \frac{2|x|^{r}}{r!}, \frac{|x|^{r+1}}{(r+1)!} \right)
    \]

    **Proof of Lemma:**
    We use mathematical induction on $r$. First, establish the recursive integral relation for the remainder:

    \[
    \int_{0}^{x}iR_{r-1}(u)du=\int_{0}^{x}i \left( e^{iu}-\sum_{j=0}^{r-1}\frac{(iu)^{j}}{j!} \right) du = \left[ e^{iu}-\sum_{j=0}^{r-1}\frac{(iu)^{j+1}}{(j+1)!} \right]_{0}^{x} = e^{ix}-\sum_{k=0}^{r}\frac{(ix)^{k}}{k!} = R_{r}(x)
    \]

    Hence we have the recurrence $R_{r}(x)=i\int_{0}^{x}R_{r-1}(u)du$.

    * **Base case ($r=0$)**:
        $R_{0}(x)=e^{ix}-1$. On one hand, by the triangle inequality: $|R_{0}(x)|=|e^{ix}-1|\le|e^{ix}|+1=2$.
        On the other hand, using the integral form: $|R_{0}(x)|=|i\int_{0}^{x}e^{iu}du|\le\int_{0}^{|x|}|e^{iu}|du=|x|$.
        Combining gives $|R_{0}(x)|\le \min(2,|x|)$, and the statement holds.

    * **Inductive step**:
        Assume the inequality holds for $r-1$, i.e., $|R_{r-1}(x)|\le\frac{2|x|^{r-1}}{(r-1)!}$ and $|R_{r-1}(x)|\le\frac{|x|^{r}}{r!}$.
        Using the integral relation $R_{r}(x)=i\int_{0}^{x}R_{r-1}(u)du$, taking absolute values:

        \[
        |R_{r}(x)|\le\int_{0}^{|x|}|R_{r-1}(u)|du
        \]

        1. Integrating the first part of the inductive hypothesis: $|R_{r}(x)|\le\int_{0}^{|x|}\frac{2u^{r-1}}{(r-1)!}du = \frac{2|x|^{r}}{r!}$.
        2. Integrating the second part of the inductive hypothesis: $|R_{r}(x)|\le\int_{0}^{|x|}\frac{u^{r}}{r!}du = \frac{|x|^{r+1}}{(r+1)!}$.
        Hence, the lemma holds for all $r\ge0$.

    <br>

    **Proof of Theorem:**
    By Taylor's theorem, for any real $x$:

    \[
    e^{ix}=\sum_{j=0}^{r}\frac{(ix)^{j}}{j!}+R_{r}(x)
    \]

    Substituting $x=tX$ and taking expectation, we obtain the characteristic function:

    \[
    \phi_{X}(t)=E[e^{itX}]=\sum_{j=0}^{r}\frac{(it)^{j}}{j!}E[X^{j}]+E[R_{r}(tX)]
    \]

    We need to show $E[R_{r}(tX)]=o(|t|^{r})$, i.e., $\frac{1}{|t|^{r}}E[R_{r}(tX)]\rightarrow0$ as $t\rightarrow0$.
    Define the random variable $Z_{t}=\frac{R_{r}(tX)}{|t|^{r}}$. Using the bound from the lemma:

    \[
    |Z_{t}|=\frac{|R_{r}(tX)|}{|t|^{r}}\le \min \left( \frac{2|X|^{r}}{r!}, \frac{|t|\cdot|X|^{r+1}}{(r+1)!} \right)
    \]

    * **Almost sure convergence**: For fixed $X$, as $t\rightarrow0$, $\frac{|t|\cdot|X|^{r+1}}{(r+1)!}\rightarrow0$, so $Z_{t} \xrightarrow{a.s.} 0$.
    * **Dominating function**: For all $t$, we have $|Z_{t}|\le\frac{2}{r!}|X|^{r}$. Since $E|X|^{r}<\infty$, this bound is integrable.
    * **Applying the Dominated Convergence Theorem (DCT)**:

    \[
    \lim_{t\rightarrow0}E[Z_{t}]=E[\lim_{t\rightarrow0}Z_{t}]=E[0]=0
    \]

    Hence $E[R_{r}(tX)]=o(|t|^{r})$. The conclusion is proved.

---

## Problem 2: Corollary of Liapounov CLT

**Problem:**
Let $\{X_{nj}:1\le j\le k_n\}_{n\ge1}$ be a triangular array of independent random variables, and denote $\sigma_{n}^{2}=Var(\sum_{j=1}^{k_n}X_{nj})$. Suppose there exist constants such that $|X_{nj}/\sigma_{n}|\le M_{nj}$ a.e., and $\lim_{n\rightarrow\infty} \max_{1\le j \le k_n} M_{nj}=0$. Show that $S_{n}=\sum_{j=1}^{k_{n}}X_{nj}$ satisfies:

\[
\frac{S_{n}-E(S_{n})}{\sigma_{n}} \xrightarrow{d} \mathcal{N}(0,1)
\]

??? success "Solution (click to expand)"

    **Proof:**
    Let $M_{n}=\max_{1\le j\le k_{n}}M_{nj}$. It is given that $M_{n}\rightarrow0$ as $n\rightarrow\infty$. By the condition, we have:

    \[
    |X_{nj}|\le M_{nj}\sigma_{n}\le M_{n}\sigma_{n} \quad \text{a.e.}
    \]

    To apply the Liapounov central limit theorem, we verify the Liapounov condition for $\delta=1$, i.e., show that the sum of third absolute central moments tends to 0:

    \[
    \lim_{n\rightarrow\infty} \frac{1}{\sigma_{n}^{3}}\sum_{j=1}^{k_{n}}E|X_{nj}-E[X_{nj}]|^{3} = 0
    \]

    1.  **Deviation bound**: By the triangle inequality,
        $|X_{nj}-E[X_{nj}]|\le|X_{nj}|+E|X_{nj}|\le M_{n}\sigma_{n}+M_{n}\sigma_{n}=2M_{n}\sigma_{n}$.

    2.  **Third moment bound**:
        $E|X_{nj}-E[X_{nj}]|^{3}=E[|X_{nj}-E[X_{nj}]|\cdot|X_{nj}-E[X_{nj}]|^{2}] \le E[(2M_{n}\sigma_{n})\cdot|X_{nj}-E[X_{nj}]|^{2}] = 2M_{n}\sigma_{n}Var(X_{nj})$.

    3.  **Summing and simplifying**:
        Since the random variables within each row are independent, $\sum_{j=1}^{k_{n}}Var(X_{nj})=\sigma_{n}^{2}$. Then:

        \[
        \frac{1}{\sigma_{n}^{3}}\sum_{j=1}^{k_{n}}E|X_{nj}-E[X_{nj}]|^{3} \le \frac{2M_{n}\sigma_{n}}{\sigma_{n}^{3}} \sum_{j=1}^{k_{n}}Var(X_{nj}) = \frac{2M_{n}\sigma_{n}}{\sigma_{n}^{3}} \cdot \sigma_{n}^{2} = 2M_{n}
        \]

    Since $M_{n}\rightarrow0$, the Liapounov condition is satisfied. By the Liapounov CLT, the conclusion follows.

---

## Problem 3: Null array

**Problem:**
Let $\{X_{nj}:1\le j\le k_{n}\}_{n\ge1}$ be a triangular array, with $E(X_{nj})=\alpha_{nj}$. Show the following implications: $(d)\Rightarrow(c)\Rightarrow(b)\Rightarrow(a)$:

(a) $\lim_{n\rightarrow\infty}P(|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})=0$ for each $j$.

(b) $\lim_{n\rightarrow\infty} \max_{j} P(|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})=0$.

(c) $\lim_{n\rightarrow\infty}P(\max_{j} |X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})=0$.

(d) $\lim_{n\rightarrow\infty} \sum_{j=1}^{k_{n}}P(|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})=0$.

??? success "Solution (click to expand)"

    **Proof:**

    * **$(d)\Rightarrow(c)$**:
        The event that the maximum exceeds the threshold is the union of the events for each component. By **Boole's inequality (subadditivity of probability)**:

        \[
        P(\max_{1\le j\le k_{n}}|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n}) = P \left( \bigcup_{j=1}^{k_{n}} \{|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n}\} \right) \le \sum_{j=1}^{k_{n}} P(|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})
        \]

        If (d) holds, the right-hand side tends to 0, so (c) must hold.

    * **$(c)\Rightarrow(b)$**:
        For any fixed $j_{0}$, we have:
        $\{|X_{nj_{0}}-\alpha_{nj_{0}}|>\epsilon\sigma_{n}\} \subseteq \{\max_{1\le j\le k_{n}}|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n}\}$.
        By monotonicity of probability:

        \[
        P(|X_{nj_{0}}-\alpha_{nj_{0}}|>\epsilon\sigma_{n}) \le P(\max_{1\le j\le k_{n}}|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})
        \]

        Taking the maximum over the left-hand side preserves the inequality. Hence (c) implies (b).

    * **$(b)\Rightarrow(a)$**:
        This implication is obvious. The probability for each component is bounded by the maximum probability. If the maximum tends to 0, then each individual probability must also tend to 0.

---

## Problem 4: Sufficiency of Lindeberg's condition

**Problem:**
Let $X_{n}\sim\mathcal{N}(0,2^{-n})$ be a sequence of independent random variables. Show that the sequence satisfies the central limit theorem, but **Lindeberg's condition fails**.

??? success "Solution (click to expand)"

    **Proof:**

    **Step 1: Compute the variance and the limiting distribution**
    
    The variance of the sum of independent random variables is the sum of variances:

    \[
    \sigma_{n}^{2} = Var \left( \sum_{j=1}^{n}X_{j} \right) = \sum_{j=1}^{n}2^{-j} = \frac{1/2(1-(1/2)^{n})}{1-1/2} = 1-(1/2)^{n}
    \]

    As $n\rightarrow\infty$, $\sigma_{n}^{2}\rightarrow1$. Compute the characteristic function of the standardized sum $T_{n} = S_{n}/\sigma_{n}$:

    \[
    \phi_{T_{n}}(t) = \prod_{j=1}^{n}\phi_{X_{j}}(t/\sigma_{n}) = \prod_{j=1}^{n} \exp \left( -\frac{1}{2}\cdot2^{-j}\cdot \left( \frac{t}{\sigma_{n}} \right) ^{2} \right) = \exp \left( -\frac{t^{2}}{2\sigma_{n}^{2}} \sum_{j=1}^{n}2^{-j} \right) = \exp \left( -\frac{t^{2}}{2} \right)
    \]

    This is exactly the characteristic function of $\mathcal{N}(0,1)$. By Lévy's continuity theorem, $S_{n}/\sigma_{n} \xrightarrow{d} \mathcal{N}(0,1)$, i.e., the central limit theorem holds.

    **Step 2: Verify the failure of Lindeberg's condition**
    
    The Lindeberg quantity is defined as $L_{n}(\epsilon) = \frac{1}{\sigma_{n}^{2}}\sum_{j=1}^{n}E[X_{j}^{2}I(|X_{j}|>\epsilon\sigma_{n})]$. For any $n\ge1$:

    \[
    L_{n}(\epsilon) \ge \frac{1}{\sigma_{n}^{2}}E[X_{1}^{2}I(|X_{1}|>\epsilon\sigma_{n})]
    \]

    As $n\rightarrow\infty$, $\sigma_{n}\rightarrow1$. Therefore:

    \[
    \liminf_{n\rightarrow\infty} L_{n}(\epsilon) \ge E[X_{1}^{2}I(|X_{1}|>\epsilon)]
    \]

    Since $X_{1}\sim\mathcal{N}(0,1/2)$, for any fixed $\epsilon>0$, the above expectation is a positive constant $C>0$.
    Hence $L_{n}(\epsilon)$ does not converge to 0, and Lindeberg's condition does not hold. This shows that Lindeberg's condition is not necessary for the central limit theorem to hold.

## Problem 5: M-dependence

**Problem:**
Let $X_{1}, X_{2}, \dots$ be independent and identically distributed random variables with mean $\mu$, variance $\sigma^{2}$, and finite fourth moment.

(a) Find the joint asymptotic distribution of $\overline{X}_{n}$ and $\overline{Z}_{n} = (1/n)\sum_{i=1}^{n} X_{i}X_{i+1}$.

(b) Obtain the asymptotic distribution of the sample autocorrelation coefficient $r_{n} = \frac{\overline{Z}_{n} - \overline{X}_{n}^{2}}{(1/n)\sum X_{i}^{2} - \overline{X}_{n}^{2}}$.

??? success "Solution (click to expand)"

    **(a) Proof of the joint asymptotic distribution**

    * **Construct a linear combination**:
        Define $Y_{i} = aX_{i} + bX_{i}X_{i+1}$. Since $\{X_{i}\}$ are i.i.d., $\{Y_{i}\}$ is an $m=1$ dependent sequence.
        Let $W_{n} = a\overline{X}_{n} + b\overline{Z}_{n} = \frac{1}{n}\sum_{i=1}^{n} Y_{i}$. By the central limit theorem for $m$-dependent sequences:

    \[
    \sqrt{n}(W_{n} - E[Y_{1}]) \xrightarrow{d} \mathcal{N}(0, V)
    \]

    \[
    V = Var(Y_{1}) + 2Cov(Y_{1}, Y_{2})
    \]

    * **Compute $Var(Y_{1})$**:
        Introduce centered variables $q_{i} = X_{i} - \mu$. Then $E[q_{i}] = 0$, $E[q_{i}^{2}] = \sigma^{2}$.
        Compute the deviation of $Y_{1}$:

    \[
    Y_{1} - E[Y_{1}] = a(X_{1} - \mu) + b(X_{1}X_{2} - \mu^{2})
    \]

    \[
    = aq_{1} + b[(q_{1} + \mu)(q_{2} + \mu) - \mu^{2}]
    \]

    \[
    = (a + b\mu)q_{1} + b\mu q_{2} + bq_{1}q_{2}
    \]

    Using independence and moments of $q_{i}$:

    \[
    Var(Y_{1}) = E[(a + b\mu)^{2}q_{1}^{2} + (b\mu)^{2}q_{2}^{2} + b^{2}q_{1}^{2}q_{2}^{2}]
    \]

    \[
    = (a + b\mu)^{2}\sigma^{2} + b^{2}\mu^{2}\sigma^{2} + b^{2}\sigma^{4}
    \]

    * **Compute $Cov(Y_{1}, Y_{2})$**:
        Similarly, $Y_{2} - E[Y_{2}] = (a + b\mu)q_{2} + b\mu q_{3} + bq_{2}q_{3}$.
        Compute the expectation of the product; only terms involving $q_{2}^{2}$ are non-zero:

    \[
    Cov(Y_{1}, Y_{2}) = E[(b\mu q_{2}) \cdot ((a + b\mu)q_{2})] = ab\mu\sigma^{2} + b^{2}\mu^{2}\sigma^{2}
    \]

    * **Combine to obtain $V$**:

    \[
    V = (a^{2} + 2ab\mu + b^{2}\mu^{2})\sigma^{2} + b^{2}\mu^{2}\sigma^{2} + b^{2}\sigma^{4} + 2ab\mu\sigma^{2} + 2b^{2}\mu^{2}\sigma^{2}
    \]

    \[
    = a^{2}\sigma^{2} + 4ab\mu\sigma^{2} + b^{2}(4\mu^2\sigma^2 + \sigma^4)
    \]

    * **Conclusion**:
        By the Cramer-Wold device, the joint distribution corresponds to the covariance matrix $\Sigma$:

    \[
    \sqrt{n} \begin{pmatrix} \overline{X}_{n} - \mu \\ \overline{Z}_{n} - \mu^{2} \end{pmatrix} \xrightarrow{d} \mathcal{N} \left( 0, \begin{pmatrix} \sigma^{2} & 2\mu\sigma^{2} \\ 2\mu\sigma^{2} & 4\mu^{2}\sigma^{2} + \sigma^{4} \end{pmatrix} \right)
    \]

    **(b) Proof of the asymptotic distribution of $r_{n}$**

    * **Denominator analysis**:
        By the law of large numbers, $\frac{1}{n}\sum X_{i}^{2} \xrightarrow{p} E[X^{2}] = \mu^{2} + \sigma^{2}$ and $\overline{X}_{n} \xrightarrow{p} \mu$.
        Hence the denominator $\frac{1}{n}\sum X_{i}^{2} - \overline{X}_{n}^{2} \xrightarrow{p} \sigma^{2}$.

    * **Numerator analysis (first-order Delta method)**:
        Let $N_{n} = \overline{Z}_{n} - \overline{X}_{n}^{2}$. Define the function $f(A, C) = C - A^{2}$.
        The gradient at $(\mu, \mu^{2})$ is:

    \[
    \nabla f(\mu, \mu^{2}) = \begin{pmatrix} -2\mu \\ 1 \end{pmatrix}
    \]

    The asymptotic variance of the numerator is $\nabla f^{T} \Sigma \nabla f$:

    \[
    \begin{pmatrix} -2\mu & 1 \end{pmatrix} \begin{pmatrix} \sigma^{2} & 2\mu\sigma^{2} \\ 2\mu\sigma^{2} & 4\mu^{2}\sigma^{2} + \sigma^{4} \end{pmatrix} \begin{pmatrix} -2\mu \\ 1 \end{pmatrix} = \sigma^{4}
    \]

    * **Final result**:
        Combining with Slutsky's theorem:

    \[
    \sqrt{n}r_{n} = \frac{1}{\sigma^{2} + o_{p}(1)} \sqrt{n}(\overline{Z}_{n} - \overline{X}_{n}^{2}) \xrightarrow{d} \frac{1}{\sigma^{2}} \mathcal{N}(0, \sigma^{4}) \sim \mathcal{N}(0, 1)
    \]

---

## Problem 6: First-Order Delta Method

**Problem:**
(a) Find the joint limiting distribution of $(\sqrt{n}(\overline{X} - \mu), \sqrt{n}(S^{2} - \sigma^{2}))$, and discuss the condition for asymptotic independence.

(b) Let $X_{1}, \dots, X_{n}$ be a Poisson sample with mean $\theta$. Find the variance-stabilizing transformation of the sample mean and construct a confidence interval.

??? success "Solution (click to expand)"

    **(a) Joint distribution and independence**

    * **Multivariate central limit theorem**:
        Let $M_{n}^{(2)} = \frac{1}{n}\sum_{i=1}^{n}(X_{i} - \mu)^{2}$. The elements of the asymptotic covariance matrix $\Sigma^{*}$ are:

    \[
    \Sigma_{11}^{*} = Var(X_{1}) = \sigma^{2}
    \]

    \[
    \Sigma_{22}^{*} = Var((X_{1} - \mu)^{2}) = E[(X_{1} - \mu)^{4}] - \sigma^{4}
    \]

    \[
    \Sigma_{12}^{*} = Cov(X_{1} - \mu, (X_{1} - \mu)^{2}) = E[(X_{1} - \mu)^{3}]
    \]

    * **Asymptotic equivalence**:
        Since $S_{n}^{2} = \frac{n}{n-1} [M_{n}^{(2)} - (\overline{X}_{n} - \mu)^{2}]$, the term $\sqrt{n}(\overline{X}_{n} - \mu)^{2} = o_{p}(1)$.
        Therefore $\sqrt{n}(S_{n}^{2} - \sigma^{2})$ is asymptotically equivalent to $\sqrt{n}(M_{n}^{(2)} - \sigma^{2})$.

    * **Conclusion**:
        The joint distribution is:

    \[
    \sqrt{n} \begin{pmatrix} \overline{X}_{n} - \mu \\ S_{n}^{2} - \sigma^{2} \end{pmatrix} \xrightarrow{d} \mathcal{N} \left( 0, \begin{pmatrix} \sigma^{2} & E[(X_{1} - \mu)^{3}] \\ E[(X_{1} - \mu)^{3}] & E[(X_{1} - \mu)^{4}] - \sigma^{4} \end{pmatrix} \right)
    \]

    *Condition for asymptotic independence*: The covariance term $E[(X_{1} - \mu)^{3}] = 0$, i.e., the distribution has zero third central moment (skewness), e.g., symmetric distributions.

    **(b) Variance-stabilizing transformation and confidence interval**

    * **Derivation of the variance-stabilizing transformation**:
        By the CLT, $\sqrt{n}(\overline{X}_{n} - \theta) \xrightarrow{d} \mathcal{N}(0, \theta)$.
        Let the transformation be $g(x)$. By the Delta method, the asymptotic variance is $[g'(\theta)]^{2}\theta$.
        Set $[g'(\theta)]^{2}\theta = 1$, then $g'(\theta) = \theta^{-1/2}$.
        Integrating gives $g(\theta) = 2\sqrt{\theta}$.

    * **Property after transformation**:

    \[
    \sqrt{n}(2\sqrt{\overline{X}_{n}} - 2\sqrt{\theta}) \xrightarrow{d} \mathcal{N}(0, 1)
    \]

    * **Construct the confidence interval**:
        For a $(1-\alpha)$ confidence level:

    \[
    P \left( -z_{\alpha/2} \le \sqrt{n}(2\sqrt{\overline{X}_{n}} - 2\sqrt{\theta}) \le z_{\alpha/2} \right) \approx 1 - \alpha
    \]

    Rearranging:

    \[
    \sqrt{\overline{X}_{n}} - \frac{z_{\alpha/2}}{2\sqrt{n}} \le \sqrt{\theta} \le \sqrt{\overline{X}_{n}} + \frac{z_{\alpha/2}}{2\sqrt{n}}
    \]

    Squaring both sides gives the confidence interval for $\theta$:

    \[
    \theta \in \left[ \left( \sqrt{\overline{X}_{n}} - \frac{z_{\alpha/2}}{2\sqrt{n}} \right)^{2}, \left( \sqrt{\overline{X}_{n}} + \frac{z_{\alpha/2}}{2\sqrt{n}} \right)^{2} \right]
    \]