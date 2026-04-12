# Chapter 2: Brownian Motion

After establishing the measure-theoretic foundation of conditional expectation and martingale theory, we officially introduce a classic continuous-time stochastic process—**Brownian Motion (BM)**, also known as the Wiener Process. It serves as the foundation for constructing the integration theory of Stochastic Differential Equations (SDEs) (such as Itô Calculus).

## 1. Basic Definition and Properties

!!! abstract "Definition: Standard Brownian Motion"

    Let $(\Omega, \mathcal{F}, P)$ be a probability space on which a real-valued stochastic process $W = \{W(t), t \ge 0\}$ is defined. If $W$ satisfies the following four conditions, it is called a **Standard Brownian Motion**:

    * **1. Initial Zero**: $P(W(0) = 0) = 1$ a.s.

    * **2. Independent Increments**: For any time partition $0 \le t_1 < t_2 < \dots < t_n$, the sequence of increments:

    $$
    W(t_1), W(t_2) - W(t_1), \dots, W(t_n) - W(t_{n-1})
    $$

    are mutually independent.

    * **3. Stationary Gaussian Increments**: For any $0 \le s < t$, the increment follows a normal distribution with mean 0 and variance equal to the time difference:

    $$
    W(t) - W(s) \sim \mathcal{N}(0, t-s)
    $$

    * **4. Continuous Paths**: Almost all sample paths $t \mapsto W(t, \omega)$ are continuous (i.e., continuous a.s.).

From the above definition, we can immediately derive the low-order moments and covariance structure of Brownian motion, which are crucial for subsequently deriving the properties of white noise.

!!! info "Basic Properties: Moments and Covariance Structure"

    * **1. Mean and Variance**:

    Since $W(t) = W(t) - W(0) \sim \mathcal{N}(0, t)$, it is obvious that:

    $$
    E[W(t)] = 0, \quad Var(W(t)) = t
    $$

    * **2. Autocovariance**:

    For any $s, t \ge 0$, we have:

    $$
    R(s, t) = Cov(W(s), W(t)) = E[W(s)W(t)] = s \wedge t = \min(s, t)
    $$

    ??? proof "Derivation of Autocovariance (Click to expand)"

        Without loss of generality, assume $0 \le s \le t$. We split $W(t)$ into incremental forms:

        $$
        E[W(s)W(t)] = E\big[W(s) \big( W(t) - W(s) + W(s) \big)\big]
        $$

        Expanding this yields:

        $$
        = E[W(s)(W(t) - W(s))] + E[W(s)^2]
        $$

        Since Brownian motion has independent increments, $W(t) - W(s)$ and $W(s) = W(s) - W(0)$ are mutually independent. Moreover, since the mean of the increment is 0, the first term vanishes:

        $$
        = E[W(s)] E[W(t) - W(s)] + Var(W(s)) = 0 + s = s
        $$

        Because we assumed $s \le t$, the general case can be written as $s \wedge t$. $\square$

!!! tip "Joint PDF and Transition Density of Brownian Motion"

    The values of Brownian motion at time points $0 < t_1 < t_2 < \dots < t_n$, denoted as $(W(t_1), \dots, W(t_n))$, follow a multivariate normal distribution.
    Due to the Markov property, its Joint Probability Density Function (Joint PDF) can be elegantly expressed as a product of **Transition Density Functions**:

    Define the Gaussian transition density (transitioning from spatial point $y$ to $x$ over time $t$):

    $$
    g(x, t | y) = \frac{1}{\sqrt{2\pi t}} e^{-\frac{|x-y|^2}{2t}}
    $$

    Then the joint probability density is:

    $$
    p(x_1, t_1; x_2, t_2; \dots; x_n, t_n) = g(x_1, t_1 | 0) g(x_2, t_2 - t_1 | x_1) \dots g(x_n, t_n - t_{n-1} | x_{n-1})
    $$

---

## 2. Introducing White Noise and the SDE Prototype

When transitioning from Ordinary Differential Equations (ODEs) $\frac{dX(t)}{dt} = b(X(t), t)$ to Stochastic Differential Equations (SDEs), we need to add a noise term $\xi(t)$.
An ideal white noise $\xi(t)$ should be completely uncorrelated at different times, meaning its covariance exhibits the property of a Dirac $\delta$ function: $E[\xi(s)\xi(t)] = \delta(s-t)$.
Mathematically, this so-called white noise is precisely the **"formal derivative" of Brownian motion**, denoted as $\dot{W}(t)$.

!!! info "Theorem: The Covariance Limit of the Difference Quotient of BM is the $\delta$ Function"
    
    Consider the difference quotient process of Brownian motion $\xi_h(t) = \frac{W(t+h) - W(t)}{h}$ ($h > 0$).
    As $h \to 0$, its autocovariance function converges to the Dirac $\delta$ function in the sense of generalized functions (distributions).

    ??? proof "Limit Derivation Process (Click to expand)"

        We calculate the covariance $E[\xi_h(s) \xi_h(t)]$ of the difference quotient at different times $s, t$.
        Using the covariance property of Brownian motion $E[W(u)W(v)] = u \wedge v$:
        
        $$
        E\left[ \frac{W(s+h)-W(s)}{h} \frac{W(t+h)-W(t)}{h} \right] = \frac{1}{h^2} \Big( (s+h \wedge t+h) - (s \wedge t+h) - (s+h \wedge t) + (s \wedge t) \Big)
        $$

        Assuming $s \le t$, we analyze the non-zero regions of the above equation:

        * 1. When the time difference $|t-s| \ge h$, the four terms above cancel each other out, resulting in $0$. This indicates that as long as the time interval is greater than $h$, the difference quotient process is uncorrelated.

        * 2. When the time difference $|t-s| < h$, there is an overlapping interval. Calculation shows the covariance is:

        $$
        \varphi_h(t-s) = \frac{1}{h^2} (h - |t-s|)
        $$

        This is an isosceles triangle function with a base width of $2h$ and a height of $\frac{1}{h}$.
        Obviously, the integral $\int_{-\infty}^\infty \varphi_h(x) dx = 1$.
        As $h \to 0$, this function tends to 0 at non-zero points and approaches infinity at 0, while its integral remains constantly 1. This perfectly aligns with the definition of the Dirac $\delta$ function:

        $$
        \lim_{h \to 0} E[\xi_h(s) \xi_h(t)] = \delta(t-s)
        $$

        This explains why SDEs are usually written in differential form $dX(t) = b(X,t)dt + \sigma(X,t)dW(t)$, because the true derivative of $W(t)$ does not exist and can only be treated as a generalized function. $\square$

---

## 3. Multi-dimensional Brownian Motion and Core Properties

In quantitative finance and multi-particle systems, we often face high-dimensional stochastic phenomena.

!!! abstract "Definition: Multi-dimensional Brownian Motion (n-dimensional BM)"

    An $n$-dimensional Brownian motion is defined as a vector process $W(t) = (W^1(t), W^2(t), \dots, W^n(t))^T$, where:

    * 1. Each component $W^k(t)$ is a standard one-dimensional Brownian motion.

    * 2. The components are mutually independent: i.e., for any $k \neq l$, the $\sigma$-algebra $\sigma(W^k(t), t \ge 0)$ is independent of $\sigma(W^l(t), t \ge 0)$.
    
    The covariance structure between its components is:

    $$
    E[W^k(t) W^l(s)] = (t \wedge s) \delta_{kl}
    $$

    *(Here $\delta_{kl}$ is the Kronecker delta, not the Dirac $\delta$ function)*

The classical beauty of Brownian motion lies in the fact that it possesses the characteristics of both a **Martingale** and a **Markov Process**.

!!! success "Theorem: Brownian Motion is a Continuous Martingale"

    Let $\{\mathcal{F}_t\}_{t \ge 0}$ be the natural filtration generated by Brownian motion, $\mathcal{F}_t = \sigma(W(u), 0 \le u \le t)$. Then $W(t)$ is a martingale.

    ??? proof "Proof of Martingale Property (Click to expand)"

        For any $s \le t$, we need to prove $E[W(t) | \mathcal{F}_s] = W(s)$.
        Construct independence by splitting the increments:

        $$
        E[W(t) | \mathcal{F}_s] = E[W(t) - W(s) + W(s) | \mathcal{F}_s]
        $$

        By the linearity property of conditional expectation:

        $$
        = E[W(t) - W(s) | \mathcal{F}_s] + E[W(s) | \mathcal{F}_s]
        $$

        Because Brownian motion has independent increments, the future increment $W(t) - W(s)$ is completely independent of the historical filtration $\mathcal{F}_s$, thus independence implies irrelevance (equals the unconditional expectation);
        Meanwhile, $W(s)$ itself is $\mathcal{F}_s$-measurable, so taking out what is known (pull it out directly):

        $$
        = E[W(t) - W(s)] + W(s) = 0 + W(s) = W(s)
        $$

        The proof is complete. $\square$

!!! success "Theorem: Brownian Motion is a Markov Process"

    Brownian motion satisfies the Markov property: "the future depends only on the present and is independent of the past".
    For any Borel set $B \in \mathcal{B}(\mathbb{R}^n)$ and $s \le t$:

    $$
    P(W(t) \in B | \mathcal{F}_s) = P(W(t) \in B | W(s)) \quad a.s.
    $$

    ??? proof "Strict Measure-Theoretic Proof of Markov Property (Click to expand)"

        Using the indicator function $\chi_B$ (or denoted as $I_B$), we can write the probability as a conditional expectation:

        $$
        P(W(t) \in B | \mathcal{F}_s) = E[\chi_B(W(t)) | \mathcal{F}_s]
        $$

        Introducing increments, let the function be $f(x, y) = \chi_B(x + y)$. Split $W(t)$ into $W(s)$ and $W(t) - W(s)$:

        $$
        E[\chi_B(W(t)) | \mathcal{F}_s] = E[f(W(s), W(t) - W(s)) | \mathcal{F}_s]
        $$

        Here we apply the core lemma of independence and conditional expectation in measure theory (the Freezing Lemma / Substitution Theorem):
        Since $W(s)$ is $\mathcal{F}_s$-measurable, and $W(t) - W(s)$ is independent of $\mathcal{F}_s$, we can treat $W(s)$ as a "frozen" constant $x$, take the unconditional expectation with respect to the other part, and then substitute $W(s)$ back in:

        $$
        = E[f(x, W(t) - W(s))]\Big|_{x = W(s)}
        $$

        Let the increment be $Z = W(t) - W(s) \sim \mathcal{N}(0, t-s)$, with its density function $g(z)$. The above equation equals:

        $$
        = \int_{\mathbb{R}^n} \chi_B(x + z) g(z) dz \Bigg|_{x = W(s)}
        $$

        Using integration by substitution, let $y = x + z$, then $dz = dy$:

        $$
        = \int_B g(y - x) dy \Bigg|_{x = W(s)} = \int_B \frac{1}{\sqrt{2\pi(t-s)}} e^{-\frac{|y - W(s)|^2}{2(t-s)}} dy
        $$

        This integration result obviously depends only on the value of the random variable $W(s)$, and does not depend on any information in $\mathcal{F}_s$ prior to time $s$.
        According to the definition of conditional expectation, this is exactly equal to $E[\chi_B(W(t)) | W(s)]$, which is $P(W(t) \in B | W(s))$. $\square$

---

## 4. Kolmogorov Continuity Theorem and Path Properties

The definition of Brownian motion directly assumes "continuous paths". But mathematically we must ask: given consistent finite-dimensional distributions, does there **necessarily exist** a modified version with continuous paths? This requires the extremely powerful Kolmogorov Continuity Theorem.

To quantify the "roughness" of continuity, we introduce Hölder spaces.

!!! abstract "Definition: Hölder Continuous Space $C^\gamma$"

    If a function $f(t)$ satisfies the condition that there exists a constant $K > 0$ such that:

    $$
    |f(t) - f(s)| \le K |t - s|^\gamma
    $$

    Then it is said to have Hölder continuity with exponent $\gamma$. Denoted as $f \in C^\gamma$.
    (Note: When $\gamma = 1$, it is Lipschitz continuous. The path of Brownian motion is extremely rough; it is nowhere differentiable, hence it cannot achieve Lipschitz continuity).

!!! info "Theorem: Kolmogorov Continuity Theorem"

    Let $X(t)$ be a stochastic process defined on the interval $[0, T]$. If there exist constants $\alpha > 0, \beta > 0, C > 0$ such that for all $s, t \in [0, T]$:

    $$
    E[|X(t) - X(s)|^\beta] \le C |t - s|^{1+\alpha}
    $$

    Then there exists a modification of $X(t)$ that has continuous sample paths. Furthermore, for any $\gamma \in \left(0, \frac{\alpha}{\beta}\right)$, the sample paths of this modification are almost surely (a.s.) locally $\gamma$-Hölder continuous.

    ??? proof "Core Proof Idea (Based on Borel-Cantelli Lemma) (Click to expand)"

        The rigorous proof of the theorem is quite tedious, but its core mechanism is very elegant.
        
        * **1. Consider dyadic rational points**: We take dyadic grid points $t_i = \frac{i}{2^n}$ on the interval.
        Consider the probability event $A_n$ where the increment between two adjacent points is too large:

        $$
        A_n = \left\{ \omega : \max_{0 \le i < 2^n} \left| X\left(\frac{i+1}{2^n}\right) - X\left(\frac{i}{2^n}\right) \right| > \left(\frac{1}{2^n}\right)^\gamma \right\}
        $$
        
        * **2. Use Chebyshev/Markov Inequality to bound**:

        $$
        P(A_n) \le \sum_{i=0}^{2^n-1} P\left( |X_{i+1} - X_i| > 2^{-n\gamma} \right) \le \sum_{i=0}^{2^n-1} \frac{E[|X_{i+1} - X_i|^\beta]}{2^{-n\gamma\beta}}
        $$

        Substitute the condition given by the theorem $E[|X_{i+1} - X_i|^\beta] \le C (2^{-n})^{1+\alpha}$, we obtain:

        $$
        P(A_n) \le 2^n \cdot C 2^{-n(1+\alpha)} \cdot 2^{n\gamma\beta} = C 2^{-n(\alpha - \gamma\beta)}
        $$
        
        * **3. Apply the Borel-Cantelli Lemma**:
        Since we chose $\gamma < \frac{\alpha}{\beta}$, we have $\alpha - \gamma\beta > 0$. This is a geometric series with a common ratio less than 1, therefore:

        $$
        \sum_{n=1}^\infty P(A_n) < \infty
        $$

        By the Borel-Cantelli Lemma (First Lemma), we know $P(\limsup A_n) = 0$. That is, almost surely, there exists $N(\omega)$ such that for all $n \ge N(\omega)$, the increments on the dyadic grid are tightly bounded within $2^{-n\gamma}$. Through uniform convergence, this guarantees the Hölder continuity of the limit process. $\square$

Later, applying this theorem to the continuity analysis of Brownian motion yields its important $C^{\frac{1}{2}-}$ property.

!!! success "Conclusion: The Hölder Roughness of Brownian Motion Paths is $C^{\frac{1}{2}-}$"

    The sample path of Brownian motion $W(t)$ is almost everywhere $\gamma$-Hölder continuous, for any $\gamma \in \left(0, \frac{1}{2}\right)$.
    (i.e., it infinitely approaches $1/2$-Hölder continuity, but cannot reach $1/2$).

    ??? proof "Derivation Process (Using Gaussian Moments) (Click to expand)"

        Since the increment of Brownian motion is $W(t) - W(s) \sim \mathcal{N}(0, |t-s|)$, we know that the high-order even moments of the normal distribution have clear formulas.
        For an even power $2m$ ($m \in \mathbb{N}$):

        $$
        E[|W(t) - W(s)|^{2m}] = (2m-1)!! \big( Var(W(t)-W(s)) \big)^m = (2m-1)!! |t - s|^m
        $$

        This perfectly matches the conditions of the Kolmogorov Continuity Theorem. We set:

        * $\beta = 2m$
        
        * $1 + \alpha = m \implies \alpha = m - 1$
        
        According to the theorem, the Hölder exponent $\gamma$ of the path must satisfy:

        $$
        \gamma < \frac{\alpha}{\beta} = \frac{m - 1}{2m} = \frac{1}{2} - \frac{1}{2m}
        $$

        Since this property holds for all positive integers $m$, we can let $m \to \infty$, at which point $\frac{1}{2m} \to 0$.
        Therefore, for any $\gamma$ strictly less than $\frac{1}{2}$, Brownian motion is $\gamma$-Hölder continuous. Furthermore, this leads to the conclusion that **Brownian motion paths are $C^{\frac{1}{2}-}$ continuous**. $\square$