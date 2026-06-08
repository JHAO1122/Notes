# 📈 After-Class Exercise Detailed Explanation: Large Sample Theory Assignment 4

## 1. One-step estimator

**Problem:**
Let $\{P_{\theta}\}_{\theta\in\Theta}$ be a model family, where $\Theta\subset\mathbb{R}^{d}$ is an open set, and $X_{1},...,X_{n} \stackrel{i.i.d}{\sim} P_{\theta_{0}}$. The log-likelihood function $l_{\theta}(x)=\log p_{\theta}(x)$ satisfies the given regularity conditions. Let $\hat{\theta}_{n}$ be an $n$-consistent estimator. The one-step estimator $\delta_{n}$ satisfies the following equation:

\[
\nabla L_{n}(\hat{\theta}_{n})+\nabla^{2}L_{n}(\hat{\theta}_{n})(\delta_{n}-\hat{\theta}_{n})=0
\]

(a) Assume the model family $\{P_{\theta}\}_{\theta\in\mathbb{R}}$ is the Cauchy family, with probability density function $p_{\theta}(x)= \frac{1}{\pi(1+(x-\theta)^{2})}$. Let $X_{1},...,X_{n} \stackrel{i.i.d}{\sim} P_{\theta}$ and define $\hat{\theta}_{n}=\text{Median}(X_{1},...,X_{n})$. Find the limiting distribution of $\sqrt{n}(\hat{\theta}_{n}-\theta)$.

(b) Let $\delta_{n}$ be the one-step estimator for this Cauchy family. What is its asymptotic distribution? In (a), is $\delta_{n}$ more efficient than $\hat{\theta}_{n}$?

??? success "Solution (click to expand)"

    **(a) Find the asymptotic distribution of the sample median**

    The sample median $\hat{\theta}_{n}$ corresponds to the sample quantile at $p=1/2$. According to the established asymptotic theory for sample quantiles (e.g., Theorem 5.23 in van der Vaart), if the cumulative distribution function $F$ has a continuous and strictly positive derivative (density) $f$ at the true population quantile $\theta$, then the sample median is asymptotically normal:

    \[
    \sqrt{n}(\hat{\theta}_{n}-\theta)\xrightarrow{d}\mathcal{N}\left(0,\frac{p(1-p)}{f(\theta)^{2}}\right)
    \]

    For $p=1/2,$ the numerator of the variance is $1/4.$ We now compute the Cauchy density at the true median $\theta$:

    \[
    f(\theta)=p_{\theta}(\theta)=\frac{1}{\pi(1+(\theta-\theta)^{2})}=\frac{1}{\pi}
    \]

    Substituting $f(\theta)=1/\pi$ into the asymptotic variance formula:

    \[
    V=\frac{1/4}{(1/\pi)^{2}}=\frac{\pi^{2}}{4}
    \]

    Therefore, the limiting distribution of the sample median is:

    \[
    \sqrt{n}(\hat{\theta}_{n}-\theta)\xrightarrow{d}\mathcal{N}\left(0,\frac{\pi^{2}}{4}\right)
    \]

    **(b) Find the asymptotic distribution of the one-step estimator and compare efficiency**

    According to the definition of the one-step estimator, we can rewrite it as:

    \[
    \sqrt{n}(\delta_{n}-\theta)=\sqrt{n}(\hat{\theta}_{n}-\theta)-\left[\frac{1}{n}\nabla^{2}L_{n}(\hat{\theta}_{n})\right]^{-1}\frac{1}{\sqrt{n}}\nabla L_{n}(\hat{\theta}_{n})
    \]

    We use Taylor's theorem to expand the score function $\frac{1}{\sqrt{n}}\nabla L_{n}(\hat{\theta}_{n})$ around the true parameter $\theta$:

    \[
    \frac{1}{\sqrt{n}}\nabla L_{n}(\hat{\theta}_{n})=\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)+\frac{1}{n}\nabla^{2}L_{n}(\tilde{\theta})\sqrt{n}(\hat{\theta}_{n}-\theta)
    \]

    where $\tilde{\theta}$ lies between $\hat{\theta}_{n}$ and $\theta$. Since $\hat{\theta}_{n}$ is $\sqrt{n}$-consistent, $\tilde{\theta}\xrightarrow{P}\theta$. By the Weak Law of Large Numbers and the Lipschitz continuity of the second derivative, we have:

    \[
    -\frac{1}{n}\nabla^{2}L_{n}(\hat{\theta}_{n})\xrightarrow{p}I(\theta) \quad \text{and} \quad -\frac{1}{n}\nabla^{2}L_{n}(\tilde{\theta})\xrightarrow{p}I(\theta)
    \]

    where $I(\theta)=-E[\nabla^{2}l_{\theta}(X)]$ is the Fisher information. Substituting the Taylor expansion back into our estimator, we obtain:

    \[
    \sqrt{n}(\delta_{n}-\theta)=\sqrt{n}(\hat{\theta}_{n}-\theta)+I(\theta)^{-1}\left(\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)-I(\theta)\sqrt{n}(\hat{\theta}_{n}-\theta)\right)+o_{p}(1)
    \]

    \[
    =\sqrt{n}(\hat{\theta}_{n}-\theta)+I(\theta)^{-1}\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)-\sqrt{n}(\hat{\theta}_{n}-\theta)+o_{p}(1)
    \]

    \[
    =I(\theta)^{-1}\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)+o_{p}(1)
    \]

    By the Central Limit Theorem, the standardized score function $\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)\xrightarrow{d}\mathcal{N}(0,I(\theta))$. By Slutsky's theorem, the one-step estimator achieves the efficiency of the maximum likelihood estimator:

    \[
    \sqrt{n}(\delta_{n}-\theta)\xrightarrow{d}\mathcal{N}(0,I(\theta)^{-1})
    \]

    Now, we compute the Fisher information $I(\theta)$ for the Cauchy distribution. The log-density is:

    \[
    l_{\theta}(x)=-\log\pi-\log(1+(x-\theta)^{2})
    \]

    The first derivative (score function) is:

    \[
    \nabla l_{\theta}(x)=\frac{2(x-\theta)}{1+(x-\theta)^{2}}
    \]

    The Fisher information is the expectation of the squared score:

    \[
    I(\theta)=E[(\nabla l_{\theta}(X))^{2}]=\int_{-\infty}^{\infty}\left(\frac{2(x-\theta)}{1+(x-\theta)^{2}}\right)^{2}\frac{1}{\pi(1+(x-\theta)^{2})}dx
    \]

    Let $y=x-\theta$. Since the integrand is an even function:

    \[
    I(\theta)=\frac{4}{\pi}\int_{-\infty}^{\infty}\frac{y^{2}}{(1+y^{2})^{3}}dy=\frac{8}{\pi}\int_{0}^{\infty}\frac{y^{2}}{(1+y^{2})^{3}}dy
    \]

    We can solve this via the trigonometric substitution $y=\tan z$, which gives $dy=\sec^{2}z dz$ and $1+y^{2}=\sec^{2}z$. The interval changes from $[0,\infty)$ to $[0,\pi/2)$:

    \[
    I(\theta)=\frac{8}{\pi}\int_{0}^{\pi/2}\frac{\tan^{2}z}{\sec^{6}z}\sec^{2}z dz
    \]

    \[
    =\frac{8}{\pi}\int_{0}^{\pi/2}\sin^{2}z \cos^{2}z dz
    \]

    \[
    =\frac{8}{\pi}\int_{0}^{\pi/2}\frac{1}{4}\sin^{2}(2z)dz
    \]

    \[
    =\frac{2}{\pi}\int_{0}^{\pi/2}\sin^{2}(2z)dz
    \]

    Let $u=2z$, so $du=2dz$, and the interval becomes $[0, \pi]$:

    \[
    I(\theta)=\frac{1}{\pi}\int_{0}^{\pi}\sin^{2}u du=\frac{1}{\pi}\cdot\frac{\pi}{2}=\frac{1}{2}
    \]

    Therefore, the asymptotic variance of the one-step estimator is $I(\theta)^{-1}=2$.

    \[
    \sqrt{n}(\delta_{n}-\theta)\xrightarrow{d}\mathcal{N}(0,2)
    \]

    **Efficiency comparison**: From part (a), the asymptotic variance of the sample median $\hat{\theta}_{n}$ is $\frac{\pi^{2}}{4}\approx2.467$. Since $2<\frac{\pi^{2}}{4}$, the asymptotic variance of the one-step estimator $\delta_{n}$ is strictly smaller than that of the initial median estimator $\hat{\theta}_{n}$. Therefore, $\delta_{n}$ is more efficient.

---

## 2. Unbiased estimator and U-statistic

**Problem:**
Let $X_{1},...,X_{n}$ be i.i.d. random variables, with finite $\mu=\mathbb{E}(X_{1})$ and $\overline{\mu}=\mathbb{E}(X_{1}^{-1})$.
Find a U-statistic that is an unbiased estimator of $\overline{\mu}\mu$, and derive its variance and asymptotic distribution.

??? success "Solution (click to expand)"

    **Step 1: Construct a symmetric kernel and the U-statistic**

    We want to estimate $\theta=\mu\overline{\mu}$. Consider a pair of independent observations $(X_{1},X_{2})$. Since they are independent, we have:

    \[
    \mathbb{E}[X_{1}X_{2}^{-1}]=\mathbb{E}[X_{1}]\mathbb{E}[X_{2}^{-1}]=\mu\overline{\mu}
    \]

    Thus, $g(x_{1},x_{2})=x_{1}x_{2}^{-1}$ is an unbiased estimator of $\theta$, but it is asymmetric. We construct a symmetric kernel of degree $m=2$ by averaging over permutations:

    \[
    h(X_{1},X_{2})=\frac{1}{2}(g(X_{1},X_{2})+g(X_{2},X_{1}))=\frac{1}{2}\left(\frac{X_{1}}{X_{2}}+\frac{X_{2}}{X_{1}}\right)
    \]

    Clearly, $\mathbb{E}[h(X_{1},X_{2})]=\frac{1}{2}(\mu\overline{\mu}+\mu\overline{\mu})=\mu\overline{\mu}$. The corresponding U-statistic is:

    \[
    U_{n}=\binom{n}{2}^{-1}\sum_{1\le i<j\le n}h(X_{i},X_{j})=\binom{n}{2}^{-1}\sum_{1\le i<j\le n}\frac{1}{2}\left(\frac{X_{i}}{X_{j}}+\frac{X_{j}}{X_{i}}\right)
    \]

    **Step 2: Derive the variance of $U_{n}$**

    According to Hoeffding's theorem, the variance of a U-statistic of degree $m=2$ is given by:

    \[
    Var(U_{n})=\binom{n}{2}^{-1}\sum_{k=1}^{2}\binom{2}{k}\binom{n-2}{2-k}\zeta_{k}=\frac{2}{n(n-1)}(2(n-2)\zeta_{1}+\zeta_{2})
    \]

    First, we compute $h_{1}(x_{1})=\mathbb{E}[h(X_{1},X_{2})|X_{1}=x_{1}]$:

    \[
    h_{1}(x_{1})=\mathbb{E}\left[\frac{1}{2}\left(\frac{x_{1}}{X_{2}}+\frac{X_{2}}{x_{1}}\right)\right]=\frac{1}{2}\overline{\mu}x_{1}+\frac{1}{2}\frac{\mu}{x_{1}}
    \]

    Now we compute $\zeta_{1}=Var(h_{1}(X_{1}))$. Let $\sigma^{2}=Var(X_{1})$ and $\overline{\sigma}^{2}=Var(X_{1}^{-1})$. Using the property $Var(aX+bY)=a^{2}Var(X)+b^{2}Var(Y)+2abCov(X,Y)$:

    \[
    \zeta_{1}=Var\left(\frac{1}{2}\overline{\mu}X_{1}+\frac{1}{2}\mu X_{1}^{-1}\right)
    \]

    \[
    =\frac{1}{4}\overline{\mu}^{2}Var(X_{1})+\frac{1}{4}\mu^{2}Var(X_{1}^{-1})+2\left(\frac{1}{2}\overline{\mu}\right)\left(\frac{1}{2}\mu\right)Cov(X_{1},X_{1}^{-1})
    \]

    \[
    =\frac{1}{4}(\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(\mathbb{E}[X_{1}X_{1}^{-1}]-\mathbb{E}[X_{1}]\mathbb{E}[X_{1}^{-1}]))
    \]

    \[
    =\frac{1}{4}(\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(1-\mu\overline{\mu}))
    \]

    Next, we compute $\zeta_{2}=Var(h(X_{1},X_{2}))=\mathbb{E}[h^{2}(X_{1},X_{2})]-(\mu\overline{\mu})^{2}$:

    \[
    \mathbb{E}[h^{2}(X_{1},X_{2})]=\mathbb{E}\left[\frac{1}{4}\left(\frac{X_{1}^{2}}{X_{2}^{2}}+\frac{X_{2}^{2}}{X_{1}^{2}}+2\frac{X_{1}X_{2}}{X_{2}X_{1}}\right)\right]
    \]

    \[
    =\frac{1}{4}(\mathbb{E}[X_{1}^{2}]\mathbb{E}[X_{2}^{-2}]+\mathbb{E}[X_{2}^{2}]\mathbb{E}[X_{1}^{-2}]+2)
    \]

    \[
    =\frac{1}{2}\mathbb{E}[X_{1}^{2}]\mathbb{E}[X_{1}^{-2}]+\frac{1}{2}
    \]

    Recall the relations between second moments and variances: $\mathbb{E}[X_{1}^{2}]=\sigma^{2}+\mu^{2}$ and $\mathbb{E}[X_{1}^{-2}]=\overline{\sigma}^{2}+\overline{\mu}^{2}$. Substituting these into the equation:

    \[
    \zeta_{2}=\frac{1}{2}(\sigma^{2}+\mu^{2})(\overline{\sigma}^{2}+\overline{\mu}^{2})+\frac{1}{2}-\mu^{2}\overline{\mu}^{2}
    \]

    \[
    =\frac{1}{2}(\sigma^{2}\overline{\sigma}^{2}+\sigma^{2}\overline{\mu}^{2}+\mu^{2}\overline{\sigma}^{2}+\mu^{2}\overline{\mu}^{2})+\frac{1}{2}-\mu^{2}\overline{\mu}^{2}
    \]

    \[
    =\frac{1}{2}(\sigma^{2}\overline{\sigma}^{2}+\sigma^{2}\overline{\mu}^{2}+\mu^{2}\overline{\sigma}^{2}-\mu^{2}\overline{\mu}^{2}+1)
    \]

    Substituting $\zeta_{1}$ and $\zeta_{2}$ back into the formula gives the exact variance:

    \[
    Var(U_{n})=\frac{4(n-2)}{n(n-1)}\zeta_{1}+\frac{2}{n(n-1)}\zeta_{2}
    \]

    \[
    =\frac{n-2}{n(n-1)}(\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(1-\mu\overline{\mu}))+\frac{1}{n(n-1)}(\sigma^{2}\overline{\sigma}^{2}+\sigma^{2}\overline{\mu}^{2}+\mu^{2}\overline{\sigma}^{2}-\mu^{2}\overline{\mu}^{2}+1)
    \]

    **Step 3: Derive the asymptotic distribution**

    According to the projection theorem for U-statistics (Hoeffding's CLT), since $\mathbb{E}[h^{2}]<\infty$ and $\zeta_{1}>0$, the U-statistic is asymptotically normal, driven by the first-order projection:

    \[
    \sqrt{n}(U_{n}-\mu\overline{\mu})\rightarrow\mathcal{N}(0,m^{2}\zeta_{1})
    \]

    where the kernel degree is $m=2$. Therefore, the asymptotic variance is $4\zeta_{1}$:

    \[
    4\zeta_{1}=\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(1-\mu\overline{\mu})
    \]

    Thus, the asymptotic distribution is:

    \[
    \sqrt{n}(U_{n}-\mu\overline{\mu})\rightarrow\mathcal{N}(0,\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(1-\mu\overline{\mu}))
    \]

---

## 3. Hoeffding's theorem

**Problem:**
Let $X_{1},...,X_{n} \stackrel{i.i.d}{\sim} P\in\mathcal{P}$ and $\mathbb{E}_{P}h^{2}(X_{1},...,X_{m})<\infty$, then:

\[
Var_{P}(U_{n})=\binom{n}{m}^{-1}\sum_{k=1}^{m}\binom{m}{k}\binom{n-m}{m-k}\xi_{k}
\]

where $\xi_{k}=Var_{P}(h_{k}(X_{1},...,X_{k}))$.

(a) Prove $\frac{m^{2}}{n}\xi_{1}\le Var_{P}(U_{n})\le\frac{m}{n}\xi_{m}$

(b) Prove $(n+1)Var_{P}(U_{n+1})\le nVar_{P}(U_{n})$

??? success "Solution (click to expand)"

    **(a) Prove the upper and lower bounds**

    According to the Hoeffding H-decomposition (ANOVA decomposition) for U-statistics, we can rewrite the variance of $U_{n}$ in terms of the projection variances $\eta_{j}$ (the variances of the $j$-th order orthogonal projection kernels):

    \[
    Var_{P}(U_{n})=\sum_{j=1}^{m}\binom{m}{j}^{2}\binom{n}{j}^{-1}\eta_{j}
    \]

    By definition, for all $j=1,...,m$, $\eta_{j}\ge0$. The relationship between $\xi_{k}$ and $\eta_{j}$ is given by:

    \[
    \xi_{k}=\sum_{j=1}^{k}\binom{k}{j}\eta_{j}
    \]

    In particular, for $k=1$, we have $\xi_{1}=\binom{1}{1}\eta_{1}=\eta_{1}$. For $k=m$, we have $\xi_{m}=\sum_{j=1}^{m}\binom{m}{j}\eta_{j}$.

    **Lower bound:** Since all terms in the variance sum are non-negative $(\eta_{j}\ge0)$, the total variance is bounded below by its first term $(j=1)$:

    \[
    Var_{P}(U_{n})\ge\binom{m}{1}^{2}\binom{n}{1}^{-1}\eta_{1}=\frac{m^{2}}{n}\eta_{1}=\frac{m^{2}}{n}\xi_{1}
    \]

    **Upper bound:** We analyze the combinatorial coefficients in the variance sum. We can rewrite them as:

    \[
    \binom{m}{j}^{2}\binom{n}{j}^{-1}=\binom{m}{j}\left[\binom{m}{j}\binom{n}{j}^{-1}\right]
    \]

    Expanding the term in brackets:

    \[
    \binom{m}{j}\binom{n}{j}^{-1}=\frac{m!}{j!(m-j)!}\frac{j!(n-j)!}{n!}=\frac{m(m-1)...(m-j+1)}{n(n-1)...(n-j+1)}=\prod_{i=0}^{j-1}\frac{m-i}{n-i}
    \]

    Since $m\le n$, for any $i\ge0$, we have $nn-in\le mn-im$, which implies $\frac{m-i}{n-i}\le\frac{m}{n}$. Therefore, the product is bounded by:

    \[
    \prod_{i=0}^{j-1}\frac{m-i}{n-i}\le\frac{m}{n}\cdot1...1=\frac{m}{n}
    \]

    Substituting this inequality back into the variance formula, we obtain:

    \[
    Var_{P}(U_{n})=\sum_{j=1}^{m}\binom{m}{j}\left[\binom{m}{j}\binom{n}{j}^{-1}\right]\eta_{j}\le\sum_{j=1}^{m}\binom{m}{j}\left(\frac{m}{n}\right)\eta_{j}=\frac{m}{n}\sum_{j=1}^{m}\binom{m}{j}\eta_{j}
    \]

    Since $\xi_{m}=\sum_{j=1}^{m}\binom{m}{j}\eta_{j}$, we conclude:

    \[
    Var_{P}(U_{n})\le\frac{m}{n}\xi_{m}
    \]

    **(b) Prove the sample size recursion inequality**

    Using the variance formula from the H-decomposition, we write the scaled variances for sample sizes $n$ and $n+1$:

    \[
    nVar_{P}(U_{n})=\sum_{j=1}^{m}\binom{m}{j}^{2}\left[n\binom{n}{j}^{-1}\right]\eta_{j}
    \]

    \[
    (n+1)Var_{P}(U_{n+1})=\sum_{j=1}^{m}\binom{m}{j}^{2}\left[(n+1)\binom{n+1}{j}^{-1}\right]\eta_{j}
    \]

    To prove the inequality, it suffices to show that the coefficient of each $\eta_{j}$ is non-increasing when the sample size increases from $n$ to $n+1$. We compare the terms in brackets. For size $n+1$:

    \[
    (n+1)\binom{n+1}{j}^{-1}=(n+1)\frac{j!(n+1-j)!}{(n+1)!}=\frac{j!(n-j+1)!}{n!}
    \]

    We can factor out $(n-j+1)$ from the numerator:

    \[
    \frac{j!(n-j+1)!}{n!}=(n-j+1)\frac{j!(n-j)!}{n!}=(n-j+1)\binom{n}{j}^{-1}
    \]

    Since for all terms in the sum $j\ge1$, we clearly have $n-j+1\le n$. Therefore:

    \[
    (n+1)\binom{n+1}{j}^{-1}=(n-j+1)\binom{n}{j}^{-1}\le n\binom{n}{j}^{-1}
    \]

    Because $\binom{m}{j}^{2}\eta_{j}\ge0$ for all $j$, multiplying this non-negative constant to both sides of the inequality preserves the direction. Summing over all $j=1,...,m$ yields:

    \[
    \sum_{j=1}^{m}\binom{m}{j}^{2}\left[(n+1)\binom{n+1}{j}^{-1}\right]\eta_{j}\le\sum_{j=1}^{m}\binom{m}{j}^{2}\left[n\binom{n}{j}^{-1}\right]\eta_{j}
    \]

    This is precisely:

    \[
    (n+1)Var_{P}(U_{n+1})\le nVar_{P}(U_{n})
    \]

## Exercise 4: U-statistics and V-statistics

**Problem:**
Let $X_1, X_2, \dots$ be a sequence of i.i.d. random variables, and let $h$ be a symmetric kernel function satisfying $0 < E[h^2(X_1, X_2)] < \infty$ and $E[h^2(X_1, X_1)] < \infty$. Define $V_n$ as:

\[
V_n := \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n h(X_i, X_j)
\]

(a) Use the central limit theorem for U-statistics to prove that $V_n$ is asymptotically normal.

(b) Can $V_n$ be written as a U-statistic? Why?

(c) Now assume $E X_1 = 0$ and $E X_1^4 < \infty$. Take the kernel $h(x, y) = xy$. Find the asymptotic relative efficiency (ARE) of $V_n$ relative to the U-statistic $U_n := \frac{2}{n(n-1)} \sum_{1 \le i < j \le n} X_i X_j$, i.e., $\frac{Avar(U_n)}{Avar(V_n)}$.

??? success "Solution (click to expand)"

    **(a) Prove the asymptotic normality of $V_n$**

    Let $\theta = E[h(X_1, X_2)]$ be the expectation of the kernel for independent observations. The corresponding U-statistic is:

    \[
    U_n = \frac{2}{n(n-1)} \sum_{1 \le i < j \le n} h(X_i, X_j)
    \]

    We can decompose the V-statistic into the U-statistic part (off-diagonal terms) and the diagonal part:

    \[
    V_n = \frac{1}{n^2} \sum_{i \ne j} h(X_i, X_j) + \frac{1}{n^2} \sum_{i=1}^n h(X_i, X_i)
    \]

    \[
    = \frac{n(n-1)}{n^2} U_n + \frac{1}{n^2} \sum_{i=1}^n h(X_i, X_i)
    \]

    \[
    = \left( 1 - \frac{1}{n} \right) U_n + \frac{1}{n^2} \sum_{i=1}^n h(X_i, X_i)
    \]

    To analyze its asymptotic distribution, we scale the centered V-statistic by $\sqrt{n}$:

    \[
    \sqrt{n}(V_n - \theta) = \sqrt{n} \left( \left( 1 - \frac{1}{n} \right) U_n + \frac{1}{n^2} \sum_{i=1}^n h(X_i, X_i) - \theta \right)
    \]

    \[
    = \left( 1 - \frac{1}{n} \right) \sqrt{n}(U_n - \theta) - \frac{1}{\sqrt{n}} \theta + \frac{1}{n\sqrt{n}} \sum_{i=1}^n h(X_i, X_i)
    \]

    Let $\mu_D = E[h(X_1, X_1)]$. By the Weak Law of Large Numbers (WLLN), $\frac{1}{n} \sum_{i=1}^n h(X_i, X_i) \xrightarrow{p} \mu_D$.
    Therefore, the last term in the expression is $\frac{1}{n} (\frac{1}{n} \sum h(X_i, X_i)) = O_p(n^{-1})$, which is $o_p(1)$.
    Meanwhile, the scaling factor $(1 - 1/n) \rightarrow 1$.
  
    By Slutsky's theorem, the scaled difference of the V-statistic converges to that of the U-statistic:

    \[
    \sqrt{n}(V_n - \theta) = \sqrt{n}(U_n - \theta) + o_p(1)
    \]

    According to the central limit theorem for U-statistics, $\sqrt{n}(U_n - \theta) \xrightarrow{d} \mathcal{N}(0, 4\zeta_1)$. Hence $V_n$ has the same asymptotic normal distribution:

    \[
    \sqrt{n}(V_n - \theta) \xrightarrow{d} \mathcal{N}(0, 4\zeta_1)
    \]

    **(b) Can $V_n$ be written as a U-statistic?**

    No. By definition, a U-statistic $U_n$ of size $n$ must be an unbiased estimator of its population parameter $\theta$, i.e., $E[U_n] = \theta$ for all $n \ge m$.
    Compute the expectation of $V_n$:

    \[
    E[V_n] = \left( 1 - \frac{1}{n} \right) \theta + \frac{1}{n} \mu_D = \theta + \frac{\mu_D - \theta}{n}
    \]

    Since $E[V_n]$ explicitly depends on the sample size $n$ (unless $\theta = \mu_D$, which is not generally true), $V_n$ is a biased estimator and cannot be a U-statistic.

    **(c) Compute the asymptotic relative efficiency (ARE)**

    Since $E X_1 = 0$ and $h(x, y) = xy$, we compute $\theta$ and the projection variance $\zeta_1$:

    \[
    \theta = E[X_1 X_2] = E[X_1] E[X_2] = 0 \cdot 0 = 0
    \]

    The projection function is:

    \[
    h_1(x_1) = E[h(x_1, X_2)] = x_1 E[X_2] = 0
    \]

    Because $h_1(x_1) \equiv 0$, we have $\zeta_1 = Var(h_1(X_1)) = 0$. This means $U_n$ is a degenerate U-statistic, and the standard $\sqrt{n}$ scaling yields a limit of 0.
    We need to use a scaling factor of $n$.

    **Analysis of $n V_n$**:
  
    \[
    V_n = \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n X_i X_j = \left( \frac{1}{n} \sum_{i=1}^n X_i \right)^2 = \bar{X}_n^2
    \]

    By the CLT, $\sqrt{n}\bar{X}_n \xrightarrow{d} \mathcal{N}(0, \sigma^2)$. By the Continuous Mapping Theorem (CMT):

    \[
    n V_n = (\sqrt{n} \bar{X}_n)^2 \xrightarrow{d} \sigma^2 \chi_1^2
    \]

    Its asymptotic variance is $Var(\sigma^2 \chi_1^2) = \sigma^4 Var(\chi_1^2) = \sigma^4 \cdot 2 = 2\sigma^4$.

    **Analysis of $n U_n$**:

    \[
    n U_n = \frac{n}{n(n-1)} \sum_{i \ne j} X_i X_j = \frac{n}{n-1} \left[ \left( \sum X_i \right)^2 - \sum X_i^2 \right] \cdot \frac{1}{n}
    \]

    \[
    = \frac{1}{n-1} \left[ n (\sqrt{n} \bar{X}_n)^2 - \sum X_i^2 \right] \xrightarrow{d} \sigma^2 \chi_1^2 - \sigma^2 = \sigma^2(\chi_1^2 - 1)
    \]

    Its asymptotic variance is $Var(\sigma^2(\chi_1^2 - 1)) = \sigma^4 Var(\chi_1^2) = 2\sigma^4$.

    Hence, the ARE is:

    \[
    ARE(U_n, V_n) = \frac{Avar(n U_n)}{Avar(n V_n)} = \frac{2\sigma^4}{2\sigma^4} = 1
    \]

---

## Exercise 5: Divide and Conquered U-statistics for big data

**Problem:**
Let $U_N$ be a U-statistic with kernel $h(x_1, \dots, x_m)$. Assume $\theta = E[h(X_1, \dots, X_m)]$ and $E[|h|] < \infty$. Define $h_1(x_1) = E[h(X_1, X_2, \dots, X_m) | X_1 = x_1]$ and assume $\zeta_1 = Var(h_1(X_1)) > 0$. Suppose the total sample size is $N = nK$, and the data are partitioned into $K$ non-overlapping subsets $D_k$ (each of size $n$). Let $U_{kn}$ be the U-statistic based on the $k$-th subset $D_k$. Define $AU_N = K^{-1} \sum_{k=1}^K U_{kn}$. Prove: if $K = o(N)$, then $\sqrt{N}(AU_N - \theta) \xrightarrow{d} \mathcal{N}(0, m^2\zeta_1)$.

??? success "Solution (click to expand)"

    **Proof:**

    We need to show that $\sqrt{N}(AU_N - \theta) \xrightarrow{d} \mathcal{N}(0, m^2\zeta_1)$. We rewrite the scaled averaged U-statistic as a sum of independent random variables.
    Since $N = nK$, we have:

    \[
    \sqrt{N}(AU_N - \theta) = \sqrt{nK} \left( \frac{1}{K} \sum_{k=1}^K (U_{kn} - \theta) \right) = \sum_{k=1}^K \sqrt{\frac{n}{K}} (U_{kn} - \theta)
    \]

    For a given total sample size $N$ (hence $n$ and $K$ fixed), define $Z_{N,k} = \sqrt{\frac{n}{K}} (U_{kn} - \theta)$. Because the $K$ subsets are non-overlapping, $\{Z_{N,k}\}_{k=1}^K$ forms a triangular array of row-wise independent random variables. We will use the Lindeberg-Feller central limit theorem.

    **Step 1: Convergence of the row sum variance**

    Let $S_N = \sum_{k=1}^K Z_{N,k}$. Since the $Z_{N,k}$ are i.i.d. for fixed $N$:

    \[
    Var(S_N) = \sum_{k=1}^K Var(Z_{N,k}) = K \cdot Var\left( \sqrt{\frac{n}{K}} U_{kn} \right) = n Var(U_{kn})
    \]

    From standard U-statistic theory, for a sample of size $n$, $n Var(U_{kn}) \rightarrow m^2 \zeta_1$ (as $n \rightarrow \infty$). Because $K = o(N)$, i.e., $K/N \rightarrow 0$, this implies $n = N/K \rightarrow \infty$. Therefore:

    \[
    Var(S_N) \rightarrow m^2 \zeta_1
    \]

    **Step 2: Verification of the Lindeberg condition**

    We need to show that for any $\epsilon > 0$:

    \[
    \lim_{N \rightarrow \infty} \sum_{k=1}^K E[Z_{N,k}^2 I(|Z_{N,k}| > \epsilon)] = 0
    \]

    Since the $Z_{N,k}$ are identically distributed, the sum simplifies to:

    \[
    K \cdot E[Z_{N,1}^2 I(|Z_{N,1}| > \epsilon)] = K \cdot E\left[ \frac{n}{K} (U_{1n} - \theta)^2 I\left(\sqrt{\frac{n}{K}}|U_{1n} - \theta| > \epsilon\right) \right]
    \]

    \[
    = E[n(U_{1n} - \theta)^2 I(| \sqrt{n}(U_{1n} - \theta) | > \epsilon \sqrt{K})]
    \]

    Let $W_n = \sqrt{n}(U_{1n} - \theta)$. The above expectation becomes $E[W_n^2 I(|W_n| > \epsilon \sqrt{K})]$.
    By the CLT for regular U-statistics, $W_n \xrightarrow{d} \mathcal{N}(0, m^2\zeta_1)$.
    Moreover, we know $E[W_n^2] = n Var(U_{1n}) \rightarrow m^2 \zeta_1$.
    Convergence in distribution together with convergence of second moments implies that $\{W_n^2\}$ is **uniformly integrable (UI)**.
    For any uniformly integrable sequence, if the truncation threshold tends to infinity, the truncated expectation tends to 0. Since $K \rightarrow \infty$, the threshold $\epsilon \sqrt{K} \rightarrow \infty$. Hence:

    \[
    \lim_{N \rightarrow \infty} E[W_n^2 I(|W_n| > \epsilon \sqrt{K})] = 0
    \]

    This satisfies the Lindeberg condition.

    **Step 3: Conclusion**

    Since the Lindeberg condition holds and the variance of the row sums converges to $m^2\zeta_1$, by the Lindeberg-Feller central limit theorem:

    \[
    \sum_{k=1}^K Z_{N,k} = \sqrt{N}(AU_N - \theta) \xrightarrow{d} \mathcal{N}(0, m^2 \zeta_1)
    \]

    The proof is complete. $\square$