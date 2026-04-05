# Chapter 2: Characteristic Functions

In asymptotic statistical theory, we urgently need powerful tools to derive convergence in distribution (weak convergence). The **Characteristic Function (cf)** is such a "divine weapon," providing a way to examine probability distributions from a **Frequency Domain** perspective and characterizing a distribution perfectly and uniquely.

---

## 1. Definition and Basic Properties of Characteristic Functions

!!! abstract "Definition 2.1: Characteristic Function"

    For any random vector $X$ with distribution function $F$, its **characteristic function (cf)** is defined as:

    $$
    \phi_X(t) = E[e^{itX}] = \int e^{itx} dF(x)
    $$

    Using Euler's formula, it can be expanded into real and imaginary parts:

    $$
    \phi_X(t) = E[\cos tX] + iE[\sin tX], \quad \text{for any } t \in \mathbb{R}
    $$

    > **Note**: Compared to the Moment Generating Function (MGF) $M_X(t) = E[e^{tX}]$, which may not exist for certain distributions (such as the Cauchy distribution), the characteristic function **always exists** for any probability distribution because $|e^{itX}| = 1$.



Characteristic functions inherit many extremely elegant mathematical properties, which play a decisive role in the subsequent derivation of limit theorems.

!!! info "Basic Properties of Characteristic Functions (Properties of CFs)"

    For a univariate random variable $X$, its characteristic function $\phi_X(t)$ satisfies:

    * **(i) Boundedness**: $|\phi_X(t)| \le \phi_X(0) = 1$.

    * **(ii) Conjugate Symmetry**: $\overline{\phi_X(t)} = \phi_X(-t)$.

    * **(iii) Uniform Continuity**: $\phi_X(t)$ is uniformly continuous on $\mathbb{R}$.

    * **(iv) Operational Closure**: $\overline{\phi_X}$, $|\phi_X|^2$, and $Re(\phi_X)$ correspond to the characteristic functions of $-X$, $X-Y$ (where $X, Y$ are i.i.d. following $F$), and the mixture distribution $(F_X + F_{-X})/2$ respectively.

    * **(v) Lattice Distribution Criterion**: If there exists $t_0 \neq 0$ s.t. $|\phi_X(t_0)| = 1$, then there must exist $a \in \mathbb{R}$ and $a \neq 0$ s.t. $P(X \in \{a + jh : j \in \mathbb{Z}\}) = 1$. That is, $X$ is a lattice random vector.

    * **(vi) Riemann-Lebesgue Lemma**: If $F$ is absolutely continuous (i.e., a density function exists), then $\lim_{|t|\to\infty} |\phi_X(t)| = 0$.

    * **(vii) Uniqueness and Fourier Inversion**: Two random variables are equal in distribution $X \stackrel{d}{=} Y$ if and only if $\phi_X(t) = \phi_Y(t)$ for all $t$. If $\phi_X$ is absolutely integrable (i.e., $\phi_X \in \mathcal{L}^1(\mathbb{R})$), then $F$ has a continuous density function, which can be obtained by the inverse transform:
    
    $$
    f(x) = \frac{1}{2\pi} \int e^{-itx} \phi_X(t) dt
    $$

    ??? proof "Supplementary Proof of the First Four Basic Properties (Click to expand)"

        **Proof (i)**:
        Using the absolute value inequality for integrals:
        
        $$
        |\phi_X(t)| = \left| E[e^{itX}] \right| \le E[|e^{itX}|] = E[1] = 1 = \phi_X(0)
        $$

        **Proof (ii)**:
        By the properties of complex conjugates:
        
        $$
        \overline{\phi_X(t)} = \overline{E[\cos tX + i\sin tX]} = E[\cos tX - i\sin tX] = E[\cos(-tX) + i\sin(-tX)] = \phi_X(-t)
        $$

        **Proof (iii)**:
        For any $t$ and increment $h$:
        
        $$
        |\phi_X(t+h) - \phi_X(t)| = \left| E[e^{i(t+h)X} - e^{itX}] \right| \le E\left[ |e^{itX}| \cdot |e^{ihX} - 1| \right] = E[|e^{ihX} - 1|]
        $$
        
        Since $|e^{ihX} - 1| \le 2$ (bounded) and $e^{ihX} - 1 \to 0$ as $h \to 0$, by the Dominated Convergence Theorem (DCT), the expectation of the expression tends to 0. Since this limit is independent of $t$, it is uniformly continuous.

        **Proof (iv) for the $|\phi_X|^2$ property**:
        Let $X, Y$ be i.i.d. Then the characteristic function of $X-Y$ is:
        
        $$
        \phi_{X-Y}(t) = E[e^{it(X-Y)}] = E[e^{itX}] E[e^{-itY}] = \phi_X(t) \phi_Y(-t)
        $$
        
        Since $X, Y$ are identically distributed, $\phi_Y(-t) = \phi_X(-t) = \overline{\phi_X(t)}$, therefore:
        
        $$
        \phi_{X-Y}(t) = \phi_X(t) \overline{\phi_X(t)} = |\phi_X(t)|^2
        $$

---

## 2. Multivariate Characteristic Functions

The above concepts can be naturally extended to high-dimensional spaces.

!!! abstract "Definition: Multivariate Characteristic Function"

    Let $X$ be a $p$-dimensional random vector. Its characteristic function is defined as:

    $$
    \phi_X(t) = E[e^{it^\top X}] = \int_{\mathbb{R}^p} e^{it^\top x} dF_X(x), \quad \text{for any } t \in \mathbb{R}^p
    $$

Multivariate characteristic functions perfectly inherit the properties of univariate ones and add properties related to matrix calculus:

* **Affine Transformation**: For a scalar $b \neq 0$, $\phi_{X/b}(t) = \phi_X(t/b)$. For a constant vector $c$, $\phi_{X+c}(t) = \exp\{it^\top c\} \phi_X(t)$.

* **Independence and Summation**: If $X$ and $Y$ are independent, then $\phi_{X+Y}(t) = \phi_X(t)\phi_Y(t)$.

* **Relationship between Moments and Derivatives**:
  * If $E\|X\| < \infty$, then the gradient $\nabla \phi_X(t)$ exists and is continuous, and $\nabla \phi_X(0) = i\mu$ (where $\mu = EX$).
  * If $E\|X\|^2 < \infty$, then the Hessian matrix $\nabla^2 \phi_X(t)$ exists and is continuous, and $\nabla^2 \phi_X(0) = -E[XX^\top]$.

* **Multivariate Normal Distribution Special Case**: 

If $X \sim N_d(\mu, \Sigma)$, 
its characteristic function is an extremely elegant quadratic exponential form:

$$
\phi_X(t) = \exp\left\{ it^\top \mu - \frac{1}{2} t^\top \Sigma t \right\}
$$

---

## 3. Lévy Continuity Theorem and Limit Applications

The most powerful application of characteristic functions is that they convert the **convergence of probability measures (Weak Convergence)** into the **pointwise convergence of complex-valued functions**.

!!! success "Theorem 2.2: Lévy-Cramér Theorem (Lévy's Continuity Theorem)"

    Let $\{X_n\}$ and $X$ be random vectors in $\mathbb{R}^d$. Then:

    $$
    X_n \xrightarrow{d} X \iff \phi_{X_n}(t) \to \phi_X(t), \quad \forall t \in \mathbb{R}^d
    $$

    ??? proof "Proof Sketch Based on Portmanteau Lemma (Click to expand)"

        **$\Rightarrow$ Direction**:
        Since the complex exponential function $e^{it^\top x} = \cos(t^\top x) + i\sin(t^\top x)$ is **bounded and continuous**.
        Directly applying Portmanteau Lemma (ii): $Ef(X_n) \to Ef(X)$ holds for any bounded continuous function $f \in C_B$, thus the characteristic function must converge pointwise.
        
        **$\Leftarrow$ Direction**:
        This is the difficult part of the theorem. The core idea is to first use the continuity of the characteristic function near the origin to prove that the sequence $\{X_n\}$ is **Tight**. By Prohorov's Theorem, a tight sequence must have a convergent subsequence. Then, using the uniqueness theorem for characteristic functions, prove that the limit distribution of all convergent subsequences must be the same as the distribution of $X$, thereby concluding that the entire sequence converges in distribution to $X$.

With this theorem, proving the Weak Law of Large Numbers (WLLN) and the Central Limit Theorem (CLT) becomes pure algebraic expansion.

!!! tip "Application 1: Central Limit Theorem for Poisson Distribution"

    Suppose $X_1, \dots, X_n$ are independent and identically distributed as $Poisson(\lambda)$.
    We know the characteristic function of $X_j$ is $\phi_X(t) = \exp\{\lambda(e^{it}-1)\}$.
    Let $\overline{X} = n^{-1}\sum X_i$. We examine the characteristic function of the standardized statistic $\frac{\overline{X} - \lambda}{\sqrt{\lambda/n}}$:

    ??? proof "Derivation Process (Click to expand)"

        Using affine transformation and independence properties:
        
        $$
        \phi_{\frac{\overline{X} - \lambda}{\sqrt{\lambda/n}}}(t) = \exp\{-it\sqrt{n\lambda}\} \cdot \phi_{\overline{X}}\left(\frac{t}{\sqrt{\lambda/n}}\right) = \exp\{-it\sqrt{n\lambda}\} \cdot \phi_X^n\left(\frac{t}{\sqrt{n\lambda}}\right)
        $$
        
        Substituting the Poisson characteristic function:
        
        $$
        = \exp\{-it\sqrt{n\lambda}\} \cdot \exp\left\{ n\lambda \left( e^{\frac{it}{\sqrt{n\lambda}}} - 1 \right) \right\}
        $$
        
        Perform a Taylor expansion of the internal exponential function $e^x = 1 + x + x^2/2 + o(x^2)$:
        
        $$
        = \exp\left\{ -it\sqrt{n\lambda} + n\lambda \left( \frac{it}{\sqrt{n\lambda}} + \frac{i^2 t^2}{2n\lambda} + o\left(\frac{1}{n\lambda}\right) \right) \right\}
        $$
        
        Expand and cancel the first-order terms:
        
        $$
        = \exp\left\{ -it\sqrt{n\lambda} + it\sqrt{n\lambda} - \frac{t^2}{2} + o(1) \right\} = \exp\left\{ -t^2/2 + o(1) \right\}
        $$
        
        As $n \to \infty$, this characteristic function converges to $e^{-t^2/2}$, which is the characteristic function of the standard normal distribution $N(0,1)$.
        Therefore, by Lévy's Continuity Theorem:
        
        $$
        \frac{\overline{X} - \lambda}{\sqrt{\lambda/n}} \xrightarrow{d} N(0, 1)
        $$

!!! tip "Application 2: Weak Law of Large Numbers (WLLN)"

    Let $Y_1, \dots, Y_n$ be i.i.d. random variables, and $\phi_Y(t)$ be differentiable at $t=0$, with derivative $i\mu = \phi'(0)$ (this is equivalent to the existence of a finite first-order moment). Then the sample mean $\overline{Y} \xrightarrow{P} \mu$.

    ??? proof "Derivation Process (Click to expand)"

        Since $\phi(0)=1$ and $\phi'(0)$ exists, there is a Taylor expansion at $t \to 0$:
        
        $$
        \phi_Y(t) = 1 + t\phi'(0) + o(t)
        $$
        
        Examine the characteristic function of the sample mean $\overline{Y}$:
        
        $$
        \phi_{\overline{Y}}(t) = \phi_Y^n\left(\frac{t}{n}\right) = \left( 1 + \frac{t}{n}\phi'(0) + o\left(\frac{t}{n}\right) \right)^n
        $$
        
        Substituting $\phi'(0) = i\mu$:
        
        $$
        = \left( 1 + \frac{it\mu}{n} + o\left(\frac{1}{n}\right) \right)^n
        $$
        
        Using the limit formula from calculus $\lim_{n \to \infty} (1 + x/n)^n = e^x$：
        
        $$
        \lim_{n \to \infty} \phi_{\overline{Y}}(t) = e^{it\mu}
        $$
        
        This is the characteristic function of a degenerate distribution (the constant $\mu$). Therefore $\overline{Y} \xrightarrow{d} \mu$. Since convergence in distribution to a constant is equivalent to convergence in probability, it is proven that $\overline{Y} \xrightarrow{P} \mu$.

---

## 4. Moments and Taylor Expansion of Characteristic Functions

As seen in the previous section, the core of asymptotic theory lies in the **Taylor Expansion of characteristic functions**. This is directly linked to the moments of random variables.

If the $r$-th moment of the random variable $X$ exists, then $\phi_X(t)$ is $r$-th order differentiable, and:

$$
\phi_X^{(r)}(t) = \int (ix)^r e^{itx} dF(x) = E[(iX)^r e^{itX}]
$$

This results in the derivative value at the origin directly giving the moment about the origin: $\phi_X^{(r)}(0) = i^r E[X^r]$.

!!! info "Theorem 2.3: Expansion of Characteristic Functions"

    If $E|X|^r < \infty$, then its characteristic function can be expanded as:

    $$
    \phi_X(t) = \sum_{j=0}^r \frac{(it)^j}{j!} E[X^j] + o(|t|^r)
    $$

    > **Note (The Moment Problem)**:
    > The characteristic function determines all moments of $X$. However, conversely, can a sequence of all moments $\{m_r := E[X^r]\}_{r=1}^\infty$ uniquely determine the distribution of $X$?
    > This is called the **Moment Problem**. The answer is: **No**. Only when **Carleman's Condition** is satisfied can the distribution be uniquely determined:
    > 
    > $$
    > \sum_{r=1}^\infty m_{2r}^{-\frac{1}{2r}} = +\infty
    > $$

Using high-order Taylor expansion, we can also prove the Central Limit Theorem in a general case extremely concisely:

??? proof "Proof of the General Central Limit Theorem (CLT) (Click to expand)"

    Suppose $X_1, \dots, X_n$ are i.i.d., with mean $\mu = E[X]$, and variance $\sigma^2 = E[X^2] < \infty$.
    Let the centered variable be $Y = X - \mu$, then $E[Y]=0, E[Y^2]=\sigma^2$.
    Its characteristic function expanded to the second order is:
    
    $$
    \phi_{X-\mu}(t) = 1 + \frac{1}{2}(it)^2 \sigma^2 + o(t^2) = 1 - \frac{t^2 \sigma^2}{2} + o(t^2)
    $$
    
    For the standardized sum $Z_n = \frac{n\overline{X} - n\mu}{\sqrt{n\sigma^2}}$, its characteristic function is:
    
    $$
    \phi_{Z_n}(t) = \phi_{X-\mu}^n\left(\frac{t}{\sigma\sqrt{n}}\right) = \left( 1 - \frac{1}{2}\left(\frac{t}{\sigma\sqrt{n}}\right)^2 \sigma^2 + o\left(\frac{t^2}{\sigma^2 n}\right) \right)^n
    $$
    
    After simplification:
    
    $$
    = \left( 1 - \frac{t^2}{2n} + o\left(\frac{1}{n}\right) \right)^n \xrightarrow{n \to \infty} e^{-t^2/2}
    $$
    
    From Lévy's Continuity Theorem, $Z_n \xrightarrow{d} N(0,1)$ is proven.

---

## 5. Cumulants and Edgeworth Expansion

If we continue to expand the characteristic function to higher orders $(r > 2)$, for example, to the fourth order:

$$
\phi_{\frac{\overline{X}-\mu}{\sqrt{\sigma^2/n}}}(t) = \left( 1 - \frac{1}{2}\frac{t^2}{n} - \frac{1}{6}\frac{it^3}{n^{3/2}} \left(\frac{m_3}{\sigma}\right)^3 + \frac{1}{24}\frac{t^4}{n^2} \left(\frac{m_4}{\sigma}\right)^4 + \dots \right)^n
$$

This leads to extremely complex algebraic expressions. To simplify this expansion for the sum of $n$ i.i.d. variables, we introduce **Cumulants (Semi-Invariants)**.

!!! abstract "Definition 2.4: Cumulant Generating Function"

    We do not expand $\phi_X(t)$ itself, but its logarithm $K_X(t) = \log \phi_X(t)$ in a Taylor series. The coefficients $\kappa_j$ of the expansion are the **Cumulants**:

    $$
    K_X(t) := \log \phi_X(t) = \sum_{j \ge 1} \frac{(it)^j}{j!} \kappa_j = \log \left\{ 1 + \sum_{j \ge 1} \frac{1}{j!} m_j (it)^j \right\}
    $$

    Using the series expansion $\log(1+x) = x - x^2/2 + x^3/3 - \dots$ to match the coefficients, we can obtain the transformation relationship between moments and cumulants (setting $\kappa_1 = m_1 = EX$):

    * $\kappa_2 = m_2 - m_1^2 = E(X - EX)^2 =: c_2$ (i.e., Variance)
    * $\kappa_3 = m_3 - 3m_1 m_2 + 2m_1^3 = E(X - EX)^3 =: c_3$
    * $\kappa_4 = m_4 - 4m_1 m_3 - 3m_2^2 + 12m_1^2 m_2 - 6m_1^4 = c_4 - 3c_2^2$

    > **Note**: Higher-order $(j > 3)$ cumulants are different from central moments. For standardized variables $Y_i = (X_i - \mu)/\sigma$, $\kappa_1=0, \kappa_2=1$. $\kappa_3$ is called **Skewness**, and $\kappa_4$ is called **Kurtosis**.

The great advantage of expanding $\log \phi(t)$ is that when independent variables are added, **cumulants are directly linear and additive**.

### Edgeworth Expansion

Through cumulants, we can write the characteristic function of the standardized sum $S_n = \frac{\overline{X}-\mu}{\sqrt{\sigma^2/n}}$ as:

$$
\phi_{S_n}(t) = \phi_Y^n\left(\frac{t}{\sqrt{n}}\right) = \exp\left\{ -\frac{t^2}{2} + \sum_{j \ge 3} \kappa_j \frac{(it)^j}{j!} n^{-\frac{j}{2}+1} \right\}
$$

Expanding the exponent terms in powers of $n^{-1/2}$:

$$
= e^{-t^2/2} \left\{ 1 + \sum_{j \ge 1} n^{-\frac{j}{2}} r_j(it) \right\}
$$

Where $r_j(\cdot)$ is a polynomial with real coefficients, with a maximum degree of $3j$ (for example, $r_1(u) = \frac{1}{6}\kappa_3 u^3$).

!!! success "Higher-Order Asymptotic Approximation: Edgeworth Expansion"

    Using the idea of Fourier inversion, since the characteristic function can be written as the product of the aforementioned polynomial and the normal characteristic function, then the **Cumulative Distribution Function $P(S_n \le x)$** must also be written as a modified form of the standard normal CDF $\Phi(x)$:

    $$
    P(S_n \le x) = \Phi(x) + n^{-\frac{1}{2}} R_1(x) + n^{-1} R_2(x) + \dots
    $$

    This is called the **Edgeworth Expansion**. It provides a more precise convergence rate and finite-sample correction than the simple CLT.

    ??? proof "Calculation of the Correction Term $R_j(x)$ and Hermite Polynomials (Click to expand)"

        To solve for $R_j(x)$, we need to find a function such that its Fourier-Stieltjes transform is exactly equal to $e^{-t^2/2} r_j(it)$:
        
        $$
        e^{-t^2/2} r_j(it) = \int e^{itx} dR_j(x)
        $$

        We use the properties of the standard normal distribution and repeated integration by parts:
        
        $$
        e^{-t^2/2} = (-it)^{-j} \int e^{itx} d\Phi^{(j)}(x)
        $$

        By treating $r_j(it)$ as a differential operator $r_j(-D)$ acting on $\Phi(x)$ where $D = d/dx$:
        
        $$
        \int e^{itx} d\{r_j(-D)\Phi(x)\} = r_j(it) e^{-t^2/2}
        $$

        This implies:
        
        $$
        R_j(x) = r_j(-D)\Phi(x)
        $$

        The derivatives of the normal distribution are exactly generated by the famous **Hermite Polynomials $He_{j}(x)$**:
        
        $$
        (-D)^j \Phi(x) = -He_{j-1}(x) e^{-t^2/2} \cdot \frac{1}{\sqrt{2\pi}}
        $$

        Thus, $R_j(x)$ can be expressed precisely by the standard normal density function and its Hermite polynomials. This is an extremely fundamental tool in high-order asymptotic theory (such as Bootstrap theory).