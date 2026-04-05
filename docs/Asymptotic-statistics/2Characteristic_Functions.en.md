# Chapter 2: Characteristic Functions

In asymptotic statistical theory, we urgently need powerful tools to derive convergence in distribution (weak convergence). [cite_start]The **Characteristic Function (cf)** is such a "divine weapon," providing a way to examine probability distributions from a **Frequency Domain** perspective and characterizing a distribution perfectly and uniquely[cite: 697, 704].

---

## 1. Definition and Basic Properties of Characteristic Functions

!!! abstract "Definition 2.1: Characteristic Function"

    [cite_start]For any random vector $X$ with distribution function $F$, its **characteristic function (cf)** is defined as[cite: 699]:

    $$
    \phi_X(t) = E[e^{itX}] = \int e^{itx} dF(x)
    $$

    Using Euler's formula, it can be expanded into real and imaginary parts[cite: 701]:

    $$
    \phi_X(t) = E[\cos tX] + iE[\sin tX], \quad \text{for any } t \in \mathbb{R}
    $$

    > **Note**: Compared to the Moment Generating Function (MGF) $M_X(t) = E[e^{tX}]$, which may not exist for certain distributions (such as the Cauchy distribution), the characteristic function **always exists** for any probability distribution because $|e^{itX}| = 1$[cite: 702, 703].

Characteristic functions inherit many extremely elegant mathematical properties, which play a decisive role in the subsequent derivation of limit theorems.

!!! info "Basic Properties of Characteristic Functions (Properties of CFs)"

    For a univariate random variable $X$, its characteristic function $\phi_X(t)$ satisfies[cite: 710]:

    * **(i) Boundedness**: $|\phi_X(t)| [cite_start]\le \phi_X(0) = 1$[cite: 711].

    * [cite_start]**(ii) Conjugate Symmetry**: $\overline{\phi_X(t)} = \phi_X(-t)$[cite: 712].

    * [cite_start]**(iii) Uniform Continuity**: $\phi_X(t)$ is uniformly continuous on $\mathbb{R}$[cite: 713].

    * [cite_start]**(iv) Operational Closure**: $\overline{\phi_X}$, $|\phi_X|^2$, and $Re(\phi_X)$ correspond to the characteristic functions of $-X$, $X-Y$ (where $X, Y$ are i.i.d. following $F$), and the mixture distribution $(F_X + F_{-X})/2$ respectively[cite: 714].

    * **(v) Lattice Distribution Criterion**: If there exists $t_0 \neq 0$ s.t. $|\phi_X(t_0)| = 1$, then there must exist $a \in \mathbb{R}$ and $a \neq 0$ s.t. $P(X \in \{a + jh : j \in \mathbb{Z}\}) = 1$. [cite_start]That is, $X$ is a lattice random vector[cite: 716, 717].

    * **(vi) Riemann-Lebesgue Lemma**: If $F$ is absolutely continuous (i.e., a density function exists), then $\lim_{|t|\to\infty} |\phi_X(t)| [cite_start]= 0$[cite: 725].

    * [cite_start]**(vii) Uniqueness and Fourier Inversion**: Two random variables are equal in distribution $X \stackrel{d}{=} Y$ if and only if $\phi_X(t) = \phi_Y(t)$ for all $t$[cite: 726]. [cite_start]If $\phi_X$ is absolutely integrable (i.e., $\phi_X \in \mathcal{L}^1(\mathbb{R})$), then $F$ has a continuous density function, which can be obtained by the inverse transform[cite: 728]:
    
    $$
    f(x) = \frac{1}{2\pi} \int e^{-itx} \phi_X(t) dt
    $$ [cite: 729]

    ??? proof "Supplementary Proof of the First Four Basic Properties (Click to expand)"

        **Proof (i)**:
        Using the absolute value inequality for integrals[cite: 711]:
        
        $$
        |\phi_X(t)| = \left| E[e^{itX}] \right| \le E[|e^{itX}|] = E[1] = 1 = \phi_X(0)
        $$

        **Proof (ii)**:
        By the properties of complex conjugates[cite: 712]:
        
        $$
        \overline{\phi_X(t)} = \overline{E[\cos tX + i\sin tX]} = E[\cos tX - i\sin tX] = E[\cos(-tX) + i\sin(-tX)] = \phi_X(-t)
        $$

        **Proof (iii)**:
        For any $t$ and increment $h$[cite: 713]:
        
        $$
        |\phi_X(t+h) - \phi_X(t)| = \left| E[e^{i(t+h)X} - e^{itX}] \right| \le E\left[ |e^{itX}| \cdot |e^{ihX} - 1| \right] = E[|e^{ihX} - 1|]
        $$
        
        Since $|e^{ihX} - 1| \le 2$ (bounded) and $e^{ihX} - 1 \to 0$ as $h \to 0$, by the Dominated Convergence Theorem (DCT), the expectation of the expression tends to 0. Since this limit is independent of $t$, it is uniformly continuous.

        **Proof (iv) for the $|\phi_X|^2$ property**:
        Let $X, Y$ be i.i.d. Then the characteristic function of $X-Y$ is[cite: 714]:
        
        $$
        \phi_{X-Y}(t) = E[e^{it(X-Y)}] = E[e^{itX}] E[e^{-itY}] = \phi_X(t) \phi_Y(-t)
        $$
        
        Since $X, Y$ are identically distributed, $\phi_Y(-t) = \phi_X(-t) = \overline{\phi_X(t)}$, therefore:
        
        $$
        \phi_{X-Y}(t) = \phi_X(t) \overline{\phi_X(t)} = |\phi_X(t)|^2
        $$

---

## 2. Multivariate Characteristic Functions

[cite_start]The above concepts can be naturally extended to high-dimensional spaces[cite: 734].

!!! abstract "Definition: Multivariate Characteristic Function"

    Let $X$ be a $p$-dimensional random vector. [cite_start]Its characteristic function is defined as[cite: 735]:

    $$
    \phi_X(t) = E[e^{it^\top X}] = \int_{\mathbb{R}^p} e^{it^\top x} dF_X(x), \quad \text{for any } t \in \mathbb{R}^p
    $$ [cite: 736]

[cite_start]Multivariate characteristic functions perfectly inherit the properties of univariate ones and add properties related to matrix calculus[cite: 739, 746]:

* [cite_start]**Affine Transformation**: For a scalar $b \neq 0$, $\phi_{X/b}(t) = \phi_X(t/b)$[cite: 752]. [cite_start]For a constant vector $c$, $\phi_{X+c}(t) = \exp\{it^\top c\} \phi_X(t)$[cite: 753].

* [cite_start]**Independence and Summation**: If $X$ and $Y$ are independent, then $\phi_{X+Y}(t) = \phi_X(t)\phi_Y(t)$[cite: 754].

* **Relationship between Moments and Derivatives**:
  * If $E\|X\| [cite_start]< \infty$, then the gradient $\nabla \phi_X(t)$ exists and is continuous, and $\nabla \phi_X(0) = i\mu$ (where $\mu = EX$)[cite: 755].
  * [cite_start]If $E\|X\|^2 < \infty$, then the Hessian matrix $\nabla^2 \phi_X(t)$ exists and is continuous, and $\nabla^2 \phi_X(0) = -E[XX^\top]$[cite: 756].

* [cite_start]**Multivariate Normal Distribution Special Case**: If $X \sim N_d(\mu, \Sigma)$, its characteristic function is an extremely elegant quadratic exponential form[cite: 756]:

  $$
  \phi_X(t) = \exp\left\{ it^\top \mu - \frac{1}{2} t^\top \Sigma t \right\}
  $$

---

## 3. Lévy Continuity Theorem and Limit Applications

[cite_start]The most powerful application of characteristic functions is that they convert the **convergence of probability measures (Weak Convergence)** into the **pointwise convergence of complex-valued functions**[cite: 762, 772].

!!! success "Theorem 2.2: Lévy-Cramér Theorem (Lévy's Continuity Theorem)"

    Let $\{X_n\}$ and $X$ be random vectors in $\mathbb{R}^d$. Then[cite: 763]:

    $$
    X_n \xrightarrow{d} X \iff \phi_{X_n}(t) \to \phi_X(t), \quad \forall t \in \mathbb{R}^d
    $$ [cite: 764]

    ??? proof "Proof Sketch Based on Portmanteau Lemma (Click to expand)"

        **$\Rightarrow$ Direction**:
        Since the complex exponential function $e^{it^\top x} = \cos(t^\top x) + i\sin(t^\top x)$ is **bounded and continuous**[cite: 768].
        Directly applying Portmanteau Lemma (ii): $Ef(X_n) \to Ef(X)$ holds for any bounded continuous function $f \in C_B$, thus the characteristic function must converge pointwise[cite: 768].
        
        **$\Leftarrow$ Direction**:
        This is the difficult part of the theorem. The core idea is to first use the continuity of the characteristic function near the origin to prove that the sequence $\{X_n\}$ is **Tight**. By Prohorov's Theorem, a tight sequence must have a convergent subsequence. Then, using the uniqueness theorem for characteristic functions, prove that the limit distribution of all convergent subsequences must be the same as the distribution of $X$, thereby concluding that the entire sequence converges in distribution to $X$[cite: 770].

With this theorem, proving the Weak Law of Large Numbers (WLLN) and the Central Limit Theorem (CLT) becomes pure algebraic expansion.

!!! tip "Application 1: Central Limit Theorem for Poisson Distribution"

    Suppose $X_1, \dots, X_n$ are independent and identically distributed as $Poisson(\lambda)$[cite: 778].
    We know the characteristic function of $X_j$ is $\phi_X(t) = \exp\{\lambda(e^{it}-1)\}$[cite: 779].
    Let $\overline{X} = n^{-1}\sum X_i$. We examine the characteristic function of the standardized statistic $\frac{\overline{X} - \lambda}{\sqrt{\lambda/n}}$[cite: 780]:

    ??? proof "Derivation Process (Click to expand)"

        Using affine transformation and independence properties[cite: 781]:
        
        $$
        \phi_{\frac{\overline{X} - \lambda}{\sqrt{\lambda/n}}}(t) = \exp\{-it\sqrt{n\lambda}\} \cdot \phi_{\overline{X}}\left(\frac{t}{\sqrt{\lambda/n}}\right) = \exp\{-it\sqrt{n\lambda}\} \cdot \phi_X^n\left(\frac{t}{\sqrt{n\lambda}}\right)
        $$
        
        Substituting the Poisson characteristic function[cite: 781]:
        
        $$
        = \exp\{-it\sqrt{n\lambda}\} \cdot \exp\left\{ n\lambda \left( e^{\frac{it}{\sqrt{n\lambda}}} - 1 \right) \right\}
        $$
        
        Perform a Taylor expansion of the internal exponential function $e^x = 1 + x + x^2/2 + o(x^2)$[cite: 781]:
        
        $$
        = \exp\left\{ -it\sqrt{n\lambda} + n\lambda \left( \frac{it}{\sqrt{n\lambda}} + \frac{i^2 t^2}{2n\lambda} + o\left(\frac{1}{n\lambda}\right) \right) \right\}
        $$
        
        Expand and cancel the first-order terms[cite: 781]:
        
        $$
        = \exp\left\{ -it\sqrt{n\lambda} + it\sqrt{n\lambda} - \frac{t^2}{2} + o(1) \right\} = \exp\left\{ -t^2/2 + o(1) \right\}
        $$
        
        As $n \to \infty$, this characteristic function converges to $e^{-t^2/2}$, which is the characteristic function of the standard normal distribution $N(0,1)$[cite: 781].
        Therefore, by Lévy's Continuity Theorem[cite: 782]:
        
        $$
        \frac{\overline{X} - \lambda}{\sqrt{\lambda/n}} \xrightarrow{d} N(0, 1)
        $$ [cite: 783]

!!! tip "Application 2: Weak Law of Large Numbers (WLLN)"

    Let $Y_1, \dots, Y_n$ be i.i.d. random variables, and $\phi_Y(t)$ be differentiable at $t=0$, with derivative $i\mu = \phi'(0)$ (this is equivalent to the existence of a finite first-order moment)[cite: 789]. Then the sample mean $\overline{Y} \xrightarrow{P} \mu$[cite: 790].

    ??? proof "Derivation Process (Click to expand)"

        Since $\phi(0)=1$ and $\phi'(0)$ exists, there is a Taylor expansion at $t \to 0$[cite: 791, 792]:
        
        $$
        \phi_Y(t) = 1 + t\phi'(0) + o(t)
        $$
        
        Examine the characteristic function of the sample mean $\overline{Y}$[cite: 792]:
        
        $$
        \phi_{\overline{Y}}(t) = \phi_Y^n\left(\frac{t}{n}\right) = \left( 1 + \frac{t}{n}\phi'(0) + o\left(\frac{t}{n}\right) \right)^n
        $$
        
        Substituting $\phi'(0) = i\mu$[cite: 792]:
        
        $$
        = \left( 1 + \frac{it\mu}{n} + o\left(\frac{1}{n}\right) \right)^n
        $$
        
        Using the limit formula from calculus $\lim_{n \to \infty} (1 + x/n)^n = e^x$[cite: 792]:
        
        $$
        \lim_{n \to \infty} \phi_{\overline{Y}}(t) = e^{it\mu}
        $$
        
        This is the characteristic function of a degenerate distribution (the constant $\mu$). Therefore $\overline{Y} \xrightarrow{d} \mu$. Since convergence in distribution to a constant is equivalent to convergence in probability, it is proven that $\overline{Y} \xrightarrow{P} \mu$[cite: 793].

---

## 4. Moments and Taylor Expansion of Characteristic Functions

As seen in the previous section, the core of asymptotic theory lies in the **Taylor Expansion of characteristic functions**. [cite_start]This is directly linked to the moments of random variables[cite: 798].

[cite_start]If the $r$-th moment of the random variable $X$ exists, then $\phi_X(t)$ is $r$-th order differentiable, and[cite: 799]:

$$
\phi_X^{(r)}(t) = \int (ix)^r e^{itx} dF(x) = E[(iX)^r e^{itX}]
[cite_start]$$ [cite: 800]

[cite_start]This results in the derivative value at the origin directly giving the moment about the origin: $\phi_X^{(r)}(0) = i^r E[X^r]$[cite: 801].

!!! info "Theorem 2.3: Expansion of Characteristic Functions"

    [cite_start]If $E|X|^r < \infty$, then its characteristic function can be expanded as[cite: 804]:

    $$
    \phi_X(t) = \sum_{j=0}^r \frac{(it)^j}{j!} E[X^j] + o(|t|^r)
    $$ [cite: 805]

    > **Note (The Moment Problem)**:
    > The characteristic function determines all moments of $X$. However, conversely, can a sequence of all moments $\{m_r := E[X^r]\}_{r=1}^\infty$ uniquely determine the distribution of $X$? [cite: 829]
    > This is called the **Moment Problem**. The answer is: **No**. Only when **Carleman's Condition** is satisfied can the distribution be uniquely determined[cite: 829, 830]:
    > 
    > $$
    > \sum_{r=1}^\infty m_{2r}^{-\frac{1}{2r}} = +\infty
    > $$ [cite: 831]

[cite_start]Using high-order Taylor expansion, we can also prove the Central Limit Theorem in a general case extremely concisely[cite: 810, 811]:

??? proof "Proof of the General Central Limit Theorem (CLT) (Click to expand)"

    Suppose $X_1, \dots, X_n$ are i.i.d., with mean $\mu = E[X]$, and variance $\sigma^2 = E[X^2] < \infty$[cite: 811].
    Let the centered variable be $Y = X - \mu$, then $E[Y]=0, E[Y^2]=\sigma^2$.
    Its characteristic function expanded to the second order is[cite: 813]:
    
    $$
    \phi_{X-\mu}(t) = 1 + \frac{1}{2}(it)^2 \sigma^2 + o(t^2) = 1 - \frac{t^2 \sigma^2}{2} + o(t^2)
    $$
    
    For the standardized sum $Z_n = \frac{n\overline{X} - n\mu}{\sqrt{n\sigma^2}}$, its characteristic function is[cite: 814]:
    
    $$
    \phi_{Z_n}(t) = \phi_{X-\mu}^n\left(\frac{t}{\sigma\sqrt{n}}\right) = \left( 1 - \frac{1}{2}\left(\frac{t}{\sigma\sqrt{n}}\right)^2 \sigma^2 + o\left(\frac{t^2}{\sigma^2 n}\right) \right)^n
    $$ [cite: 815]
    
    After simplification[cite: 815]:
    
    $$
    = \left( 1 - \frac{t^2}{2n} + o\left(\frac{1}{n}\right) \right)^n \xrightarrow{n \to \infty} e^{-t^2/2}
    $$
    
    From Lévy's Continuity Theorem, $Z_n \xrightarrow{d} N(0,1)$ is proven[cite: 816].

---

## 5. Cumulants and Edgeworth Expansion

[cite_start]If we continue to expand the characteristic function to higher orders $(r > 2)$, for example, to the fourth order[cite: 824, 825]:

$$
\phi_{\frac{\overline{X}-\mu}{\sqrt{\sigma^2/n}}}(t) = \left( 1 - \frac{1}{2}\frac{t^2}{n} - \frac{1}{6}\frac{it^3}{n^{3/2}} \left(\frac{m_3}{\sigma}\right)^3 + \frac{1}{24}\frac{t^4}{n^2} \left(\frac{m_4}{\sigma}\right)^4 + \dots \right)^n
[cite_start]$$ [cite: 826]

This leads to extremely complex algebraic expressions. [cite_start]To simplify this expansion for the sum of $n$ i.i.d. variables, we introduce **Cumulants (Semi-Invariants)**[cite: 837].

!!! abstract "Definition 2.4: Cumulant Generating Function"

    We do not expand $\phi_X(t)$ itself, but its logarithm $K_X(t) = \log \phi_X(t)$ in a Taylor series. [cite_start]The coefficients $\kappa_j$ of the expansion are the **Cumulants**[cite: 839]:

    $$
    K_X(t) := \log \phi_X(t) = \sum_{j \ge 1} \frac{(it)^j}{j!} \kappa_j = \log \left\{ 1 + \sum_{j \ge 1} \frac{1}{j!} m_j (it)^j \right\}
    $$ [cite: 840]

    Using the series expansion $\log(1+x) = x - x^2/2 + x^3/3 - \dots$ to match the coefficients, we can obtain the transformation relationship between moments and cumulants (setting $\kappa_1 = m_1 = EX$)[cite: 840, 841]:

    * [cite_start]$\kappa_2 = m_2 - m_1^2 = E(X - EX)^2 =: c_2$ (i.e., Variance) [cite: 842]
    * [cite_start]$\kappa_3 = m_3 - 3m_1 m_2 + 2m_1^3 = E(X - EX)^3 =: c_3$ [cite: 842]
    * [cite_start]$\kappa_4 = m_4 - 4m_1 m_3 - 3m_2^2 + 12m_1^2 m_2 - 6m_1^4 = c_4 - 3c_2^2$ [cite: 843]

    > [cite_start]**Note**: Higher-order $(j > 3)$ cumulants are different from central moments[cite: 844]. [cite_start]For standardized variables $Y_i = (X_i - \mu)/\sigma$, $\kappa_1=0, \kappa_2=1$[cite: 852]. [cite_start]$\kappa_3$ is called **Skewness**, and $\kappa_4$ is called **Kurtosis**[cite: 853].

[cite_start]The great advantage of expanding $\log \phi(t)$ is that when independent variables are added, **cumulants are directly linear and additive**[cite: 856].

### Edgeworth Expansion

[cite_start]Through cumulants, we can write the characteristic function of the standardized sum $S_n = \frac{\overline{X}-\mu}{\sqrt{\sigma^2/n}}$ as[cite: 875]:

$$
\phi_{S_n}(t) = \phi_Y^n\left(\frac{t}{\sqrt{n}}\right) = \exp\left\{ -\frac{t^2}{2} + \sum_{j \ge 3} \kappa_j \frac{(it)^j}{j!} n^{-\frac{j}{2}+1} \right\}
[cite_start]$$ [cite: 867]

[cite_start]Expanding the exponent terms in powers of $n^{-1/2}$[cite: 867]:

$$
= e^{-t^2/2} \left\{ 1 + \sum_{j \ge 1} n^{-\frac{j}{2}} r_j(it) \right\}
[cite_start]$$ [cite: 867]

[cite_start]Where $r_j(\cdot)$ is a polynomial with real coefficients, with a maximum degree of $3j$ (for example, $r_1(u) = \frac{1}{6}\kappa_3 u^3$)[cite: 868, 869].

!!! success "Higher-Order Asymptotic Approximation: Edgeworth Expansion"

    [cite_start]Using the idea of Fourier inversion, since the characteristic function can be written as the product of the aforementioned polynomial and the normal characteristic function [cite: 877][cite_start], then the **Cumulative Distribution Function $P(S_n \le x)$** must also be written as a modified form of the standard normal CDF $\Phi(x)$[cite: 878, 879]:

    $$
    P(S_n \le x) = \Phi(x) + n^{-\frac{1}{2}} R_1(x) + n^{-1} R_2(x) + \dots
    [cite_start]$$ [cite: 880]

    [cite_start]This is called the **Edgeworth Expansion**[cite: 879]. [cite_start]It provides a more precise convergence rate and finite-sample correction than the simple CLT[cite: 883].

    ??? proof "Calculation of the Correction Term $R_j(x)$ and Hermite Polynomials (Click to expand)"

        [cite_start]To solve for $R_j(x)$, we need to find a function such that its Fourier-Stieltjes transform is exactly equal to $e^{-t^2/2} r_j(it)$[cite: 881]:
        
        $$
        e^{-t^2/2} r_j(it) = \int e^{itx} dR_j(x)
        [cite_start]$$ [cite: 882]

        [cite_start]We use the properties of the standard normal distribution and repeated integration by parts[cite: 889]:
        
        $$
        e^{-t^2/2} = (-it)^{-j} \int e^{itx} d\Phi^{(j)}(x)
        [cite_start]$$ [cite: 890]

        Define the differential operator $D = d/dx$. [cite_start]This inspires us to replace $r_j(it)$ with the differential operator polynomial $r_j(-D)$ and apply it to $\Phi(x)$[cite: 891]:
        
        $$
        \int e^{itx} d\{r_j(-D)\Phi(x)\} = r_j(it) e^{-t^2/2}
        [cite_start]$$ [cite: 893]

        [cite_start]This means[cite: 894]:
        
        $$
        R_j(x) = r_j(-D)\Phi(x)
        $$

        [cite_start]The derivatives of the normal distribution are exactly generated by the famous **Hermite Polynomials $He_{j}(x)$**[cite: 895, 897]:
        
        $$
        (-D)^j \Phi(x) = -He_{j-1}(x) e^{-t^2/2} \cdot \frac{1}{\sqrt{2\pi}}
        $$

        Thus, $R_j(x)$ can be expressed precisely by the standard normal density function and its Hermite polynomials. This is an extremely fundamental tool in high-order asymptotic theory (such as Bootstrap theory).