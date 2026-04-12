# Chapter 2: Brownian Motion

After establishing the measure-theoretic foundations of conditional expectations and martingale theory, we formally introduce a classic continuous-time stochastic process—**Brownian Motion (BM)**, also known as the Wiener Process. It serves as the foundation for constructing the integral theory of stochastic differential equations (SDEs), such as Itô integration.

## 1. Basic Definitions and Properties

!!! abstract "Definition: Standard Brownian Motion"

    Let $(\Omega, \mathcal{F}, P)$ be a probability space, on which a real-valued stochastic process $W = \{W(t), t \ge 0\}$ is defined. If $W$ satisfies the following four conditions, it is called a **Standard Brownian Motion**:

    **1. Initial Zero Point**: $P(W(0) = 0) = 1$ a.s.
    
    **2. Independent Increments**: For any time partition $0 \le t_1 < t_2 < \dots < t_n$, the sequence of increments:

    $$
    W(t_1), W(t_2) - W(t_1), \dots, W(t_n) - W(t_{n-1})
    $$

    are mutually independent.
    
    **3. Stationary Gaussian Increments**: For any $0 \le s < t$, the increment follows a normal distribution with mean 0 and variance equal to the time difference:

    $$
    W(t) - W(s) \sim \mathcal{N}(0, t-s)
    $$
    
    **4. Continuous Paths**: Almost all sample paths $t \mapsto W(t, \omega)$ are continuous (i.e., a.s. continuous).

From the above definition, we can immediately derive the low-order moments and covariance structure of Brownian motion, which is key to deriving the properties of white noise later.

!!! info "Basic Properties: Moments and Covariance Structure"

    **1. Mean and Variance**:
    Since $W(t) = W(t) - W(0) \sim \mathcal{N}(0, t)$, we obviously have:

    $$
    E[W(t)] = 0, \quad Var(W(t)) = t
    $$

    **2. Autocovariance Function**:
    For any $s, t \ge 0$, we have:

    $$
    R(s, t) = Cov(W(s), W(t)) = E[W(s)W(t)] = s \wedge t = \min(s, t)
    $$

    ??? proof "Derivation of the Autocovariance (click to expand)"

        Without loss of generality, assume $0 \le s \le t$. We decompose $W(t)$ into an incremental form:

        $$
        E[W(s)W(t)] = E\big[W(s) \big( W(t) - W(s) + W(s) \big)\big]
        $$

        Expanding gives:

        $$
        = E[W(s)(W(t) - W(s))] + E[W(s)^2]
        $$

        Due to the independent increments property of Brownian motion, $W(t) - W(s)$ and $W(s) = W(s) - W(0)$ are independent. Furthermore, since the increments have mean 0, the first term is 0:

        $$
        = E[W(s)] E[W(t) - W(s)] + Var(W(s)) = 0 + s = s
        $$

        Since we assumed $s \le t$, the general case can be written as $s \wedge t$. $\square$

!!! tip "Joint Probability Density and Transition Density of Brownian Motion"

    The values $(W(t_1), \dots, W(t_n))$ of Brownian motion at time points $0 < t_1 < t_2 < \dots < t_n$ follow a multivariate normal distribution.
    Due to the Markov property, its Joint Probability Density (Joint PDF) can be elegantly expressed as a product of **Transition Density Functions**:

    Define the Gaussian transition density (from spatial point $y$ to $x$ over time $t$):

    $$
    g(x, t | y) = \frac{1}{\sqrt{2\pi t}} e^{-\frac{|x-y|^2}{2t}}
    $$

    Then the joint probability density is:

    $$
    p(x_1, t_1; x_2, t_2; \dots; x_n, t_n) = g(x_1, t_1 | 0) g(x_2, t_2 - t_1 | x_1) \dots g(x_n, t_n - t_{n-1} | x_{n-1})
    $$

---

## 2. Introducing White Noise and the Prototype of SDE

When transitioning from an ordinary differential equation (ODE) $\frac{dX(t)}{dt} = b(X(t), t)$ to a stochastic differential equation (SDE), we need to add a noise term $\xi(t)$.
The ideal white noise $\xi(t)$ should be completely uncorrelated at different times, meaning its covariance exhibits the property of the Dirac $\delta$ function: $E[\xi(s)\xi(t)] = \delta(s-t)$.
Mathematically, this so-called white noise is precisely the **"formal derivative" of Brownian motion** $\dot{W}(t)$.

!!! info "Theorem: The Covariance Limit of the Increment Quotient of Brownian Motion is the $\delta$ Function"

    Consider the difference quotient process of Brownian motion $\xi_h(t) = \frac{W(t+h) - W(t)}{h}$ ($h > 0$).
    As $h \to 0$, its autocovariance function converges to the Dirac $\delta$ function in the sense of generalized functions (distributions).

    ??? proof "Derivation of the Limit (Click to Expand)"

        We compute the covariance $E[\xi_h(s) \xi_h(t)]$ of the difference quotient at different times $s, t$.
        Using the covariance property of Brownian motion $E[W(u)W(v)] = u \wedge v$:

        $$
        E\left[ \frac{W(s+h)-W(s)}{h} \frac{W(t+h)-W(t)}{h} \right] = \frac{1}{h^2} \Big( (s+h \wedge t+h) - (s \wedge t+h) - (s+h \wedge t) + (s \wedge t) \Big)
        $$

        Assuming $s \le t$, we analyze the non-zero region of the above expression:
        1. When the time difference $|t-s| \ge h$, the four terms cancel each other out, resulting in $0$. This shows that as long as the time interval is greater than $h$, the difference quotient process is uncorrelated.
        2. When the time difference $|t-s| < h$, there is an overlapping interval. The calculation yields the covariance as:

        $$
        \varphi_h(t-s) = \frac{1}{h^2} (h - |t-s|)
        $$

        This is an isosceles triangular function with a base width of $2h$ and a height of $\frac{1}{h}$.
        Clearly, the integral $\int_{-\infty}^\infty \varphi_h(x) dx = 1$.
        As $h \to 0$, this function tends to 0 at non-zero points, tends to infinity at 0, and its integral remains constant at 1. This is precisely the definition of the Dirac $\delta$ function:

        $$
        \lim_{h \to 0} E[\xi_h(s) \xi_h(t)] = \delta(t-s)
        $$

        This explains why SDEs are usually written in differential form $dX(t) = b(X,t)dt + \sigma(X,t)dW(t)$, because the true derivative of $W(t)$ does not exist and can only be treated as a generalized function. $\square$

---

## 3. Multidimensional Brownian Motion and Core Properties

In quantitative finance and multi-particle systems, we often encounter high-dimensional stochastic phenomena.

!!! abstract "Definition: Multidimensional Brownian Motion (n-dimensional BM)"

    An $n$-dimensional Brownian motion is defined as a vector process $W(t) = (W^1(t), W^2(t), \dots, W^n(t))^T$, where:
    1. Each component $W^k(t)$ is a standard one-dimensional Brownian motion.
    2. The components are mutually independent: for any $k \neq l$, the $\sigma$-algebras $\sigma(W^k(t), t \ge 0)$ and $\sigma(W^l(t), t \ge 0)$ are independent.
    
    The covariance structure between its components is:

    $$
    E[W^k(t) W^l(s)] = (t \wedge s) \delta_{kl}
    $$

    *(Here $\delta_{kl}$ is the Kronecker delta, not the Dirac $\delta$ function)*

The classic nature of Brownian motion lies in its combination of **Martingale** and **Markov Process** characteristics.

!!! success "Theorem: Brownian Motion is a Continuous Martingale"

    Let $\{\mathcal{F}_t\}_{t \ge 0}$ be the natural filtration generated by the Brownian motion, $\mathcal{F}_t = \sigma(W(u), 0 \le u \le t)$. Then $W(t)$ is a martingale.

    ??? proof "Proof of Martingale Property (click to expand)"

        For any $s \le t$, we need to prove $E[W(t) | \mathcal{F}_s] = W(s)$.
        Construct independence via increment decomposition:

        $$
        E[W(t) | \mathcal{F}_s] = E[W(t) - W(s) + W(s) | \mathcal{F}_s]
        $$

        By the linearity of conditional expectation:

        $$
        = E[W(t) - W(s) | \mathcal{F}_s] + E[W(s) | \mathcal{F}_s]
        $$

        Because Brownian motion has independent increments, the future increment $W(t) - W(s)$ is completely independent of the historical information flow $\mathcal{F}_s$, so independence implies irrelevance (equals the unconditional expectation);
        And $W(s)$ itself is $\mathcal{F}_s$-measurable, so it is known as a constant (directly taken out):

        $$
        = E[W(t) - W(s)] + W(s) = 0 + W(s) = W(s)
        $$

        Q.E.D. $\square$

!!! success "Theorem: Brownian Motion is a Markov Process"

    Brownian motion satisfies the Markov property: i.e., "the future depends only on the present, not on the past".
    For any Borel set $B \in \mathcal{B}(\mathbb{R}^n)$ and $s \le t$:

    $$
    P(W(t) \in B | \mathcal{F}_s) = P(W(t) \in B | W(s)) \quad a.s.
    $$

    ??? proof "Rigorous Measure-Theoretic Proof of Markov Property (click to expand)"

        Using the indicator function $\chi_B$ (or denoted as $I_B$), we can write the probability as a conditional expectation:

        $$
        P(W(t) \in B | \mathcal{F}_s) = E[\chi_B(W(t)) | \mathcal{F}_s]
        $$

        Introduce the increment, and define the function $f(x, y) = \chi_B(x + y)$. Decompose $W(t)$ into $W(s)$ and $W(t) - W(s)$:

        $$
        E[\chi_B(W(t)) | \mathcal{F}_s] = E[f(W(s), W(t) - W(s)) | \mathcal{F}_s]
        $$

        Here we apply a core lemma in measure theory concerning independence and conditional expectation (the freezing lemma/substitution theorem):
        Because $W(s)$ is $\mathcal{F}_s$-measurable, and $W(t) - W(s)$ is independent of $\mathcal{F}_s$, we can treat $W(s)$ as a "frozen" constant $x$, take the unconditional expectation of the other part, and then substitute $W(s)$ back:

        $$
        = E[f(x, W(t) - W(s))]\Big|_{x = W(s)}
        $$

        Let the increment $Z = W(t) - W(s) \sim \mathcal{N}(0, t-s)$, with density function $g(z)$. The above expression equals:

        $$
        = \int_{\mathbb{R}^n} \chi_B(x + z) g(z) dz \Bigg|_{x = W(s)}
        $$

        Change variables by letting $y = x + z$, then $dz = dy$:

        $$
        = \int_B g(y - x) dy \Bigg|_{x = W(s)} = \int_B \frac{1}{\sqrt{2\pi(t-s)}} e^{-\frac{|y - W(s)|^2}{2(t-s)}} dy
        $$

        The result of this integral clearly depends only on the value of the random variable $W(s)$, and not on any information in $\mathcal{F}_s$ prior to time $s$.
        By the definition of conditional expectation, this is exactly equal to $E[\chi_B(W(t)) | W(s)]$, i.e., $P(W(t) \in B | W(s))$. $\square$

---

## 4. Kolmogorov Continuity Theorem and Path Properties

The definition of Brownian motion directly assumes "continuous paths." However, mathematically we need to ask: given a consistent set of finite-dimensional distributions, does there **necessarily exist** a modification with continuous paths? This requires the extremely powerful Kolmogorov continuity theorem.

To quantify the "roughness" of continuity, we introduce Hölder spaces.

!!! abstract "Definition: Hölder Continuous Space $C^\gamma$"

    A function $f(t)$ is said to be Hölder continuous with exponent $\gamma$ if there exists a constant $K > 0$ such that:

    $$
    |f(t) - f(s)| \le K |t - s|^\gamma
    $$

    It is denoted as $f \in C^\gamma$.
    (Note: When $\gamma = 1$, this is Lipschitz continuity. The paths of Brownian motion are extremely rough; they are nowhere differentiable, thus not Lipschitz continuous).

!!! info "Theorem: Kolmogorov Continuity Theorem"

    Let $X(t)$ be a stochastic process defined on an interval $[0, T]$. If there exist constants $\alpha > 0, \beta > 0, C > 0$ such that for all $s, t \in [0, T]$:

    $$
    E[|X(t) - X(s)|^\beta] \le C |t - s|^{1+\alpha}
    $$

    then $X(t)$ has a modification with continuous sample paths. Furthermore, for any $\gamma \in \left(0, \frac{\alpha}{\beta}\right)$, the sample paths of this modification are almost surely (a.s.) locally $\gamma$-Hölder continuous.

    ??? proof "Core Proof Idea (Based on the Borel-Cantelli Lemma) (Click to expand)"

        The rigorous proof of the theorem is quite involved, but its core mechanism is elegant.
        
        **1. Examine dyadic rational points:** Consider the dyadic grid points on the interval: $t_i = \frac{i}{2^n}$.
        Define the event $A_n$ where the increment between adjacent points is too large:

        $$
        A_n = \left\{ \omega : \max_{0 \le i < 2^n} \left| X\left(\frac{i+1}{2^n}\right) - X\left(\frac{i}{2^n}\right) \right| > \left(\frac{1}{2^n}\right)^\gamma \right\}
        $$
        
        **2. Bound using Chebyshev/Markov inequality:**

        $$
        P(A_n) \le \sum_{i=0}^{2^n-1} P\left( |X_{i+1} - X_i| > 2^{-n\gamma} \right) \le \sum_{i=0}^{2^n-1} \frac{E[|X_{i+1} - X_i|^\beta]}{2^{-n\gamma\beta}}
        $$

        Substituting the condition given in the theorem, $E[|X_{i+1} - X_i|^\beta] \le C (2^{-n})^{1+\alpha}$, yields:

        $$
        P(A_n) \le 2^n \cdot C 2^{-n(1+\alpha)} \cdot 2^{n\gamma\beta} = C 2^{-n(\alpha - \gamma\beta)}
        $$
        
        **3. Apply the Borel-Cantelli Lemma:**
        Since we chose $\gamma < \frac{\alpha}{\beta}$, we have $\alpha - \gamma\beta > 0$. This forms a geometric series with ratio less than 1, therefore:

        $$
        \sum_{n=1}^\infty P(A_n) < \infty
        $$

        By the Borel-Cantelli Lemma (the first lemma), $P(\limsup A_n) = 0$. That is, almost surely, there exists $N(\omega)$ such that for all $n \ge N(\omega)$, the increments on the dyadic grid are tightly controlled within $2^{-n\gamma}$. Through uniform convergence, this guarantees the Hölder continuity of the limiting process. $\square$

We then apply this theorem to the continuity analysis of Brownian motion, obtaining its important $C^{\frac{1}{2}-}$ property.

!!! success "Conclusion: The Hölder Roughness of Brownian Motion Paths is $C^{\frac{1}{2}-}$"

    The sample paths of Brownian motion $W(t)$ are almost everywhere $\gamma$-Hölder continuous for any $\gamma \in \left(0, \frac{1}{2}\right)$.
    (That is, it approaches $1/2$-Hölder continuity arbitrarily closely, but does not achieve $1/2$).

    ??? proof "Derivation (Using Gaussian Moments) (Click to expand)"

        Since the increment of Brownian motion $W(t) - W(s) \sim \mathcal{N}(0, |t-s|)$, we know the explicit formula for higher even moments of the normal distribution.
        For even powers $2m$ ($m \in \mathbb{N}$):

        $$
        E[|W(t) - W(s)|^{2m}] = (2m-1)!! \big( Var(W(t)-W(s)) \big)^m = (2m-1)!! |t - s|^m
        $$

        This perfectly matches the condition of the Kolmogorov continuity theorem. Let:
        - $\beta = 2m$
        - $1 + \alpha = m \implies \alpha = m - 1$
        
        According to the theorem, the Hölder exponent $\gamma$ of the paths must satisfy:

        $$
        \gamma < \frac{\alpha}{\beta} = \frac{m - 1}{2m} = \frac{1}{2} - \frac{1}{2m}
        $$

        Since this holds for all positive integers $m$, we can let $m \to \infty$, where $\frac{1}{2m} \to 0$.
        Therefore, for any $\gamma$ strictly less than $\frac{1}{2}$, Brownian motion is $\gamma$-Hölder continuous. Consequently, we arrive at the conclusion that **Brownian motion paths are $C^{\frac{1}{2}-}$ continuous**. $\square$