# Chapter 4: Itô Integral and Stochastic Differential Equations

In the previous chapter, we extended integration to the $L^2$ space. Now, we formally define the Itô integral for **stochastic processes** and establish a framework distinct from classical calculus—**Itô Calculus**.

## 1. Itô Integral for Adapted Processes and Isometry

For a stochastic process $G(t, \omega)$, due to the roughness of Brownian motion $W(t)$ paths, we cannot use pointwise-defined Riemann-Stieltjes integrals. The definition must be based on the left endpoints of time partitions (not knowing the future).

!!! abstract "Definition: The Itô Integral"

    Let $\{W(t)\}_{t \ge 0}$ be a standard Brownian motion and $\{\mathcal{F}(t)\}$ be its filtration.
    Assume the stochastic process $G(t) \in L^2(0, T)$, and $G(t)$ is an **Adapted Process** (i.e., its value at time $t$ depends only on the historical information $\mathcal{F}(t)$).

    For a partition $P = \{0 = t_0 < t_1 < \dots < t_m = T\}$ of the interval $[0, T]$, define its Riemann sum as:

    $$
    R_m = \sum_{k=0}^{m-1} G(t_k) \big( W(t_{k+1}) - W(t_k) \big)
    $$

    When the mesh size $|P| \to 0$, if this sum converges in the $L^2(\Omega, P)$ sense, its limit is defined as the **Itô integral** of $G(t)$ with respect to Brownian motion:

    $$
    \int_0^T G(t) dW(t) = \lim_{|P| \to 0} \sum_{k=0}^{m-1} G(t_k) \big( W(t_{k+1}) - W(t_k) \big)
    $$

Using the method of approximating by simple functions (step processes) to define this integral yields the following properties:

!!! info "Theorem 1: Core Properties of the Itô Integral"

    Let $G(t), H(t) \in L^2(0, T)$ be adapted processes, and $a, b \in \mathbb{R}$ be constants.

    **(1) Linearity**:

    $$
    \int_0^T (aG(t) + bH(t)) dW(t) = a \int_0^T G(t) dW(t) + b \int_0^T H(t) dW(t) \quad a.s.
    $$

    **(2) Zero Mean**:

    $$
    E\left[ \int_0^T G(t) dW(t) \right] = 0
    $$

    **(3) Itô Isometry**:

    $$
    E\left[ \left( \int_0^T G(t) dW(t) \right)^2 \right] = \int_0^T E[G(t)^2] dt
    $$

    ??? proof "Rigorous Proof of Zero Mean and Itô Isometry (Click to Expand)"

        We prove the properties using simple functions (Step Functions). Let $G_k = G(t_k)$.

        **1. Proving Zero Mean**:

        $$
        E\left[ \sum_{k} G_k \big(W(t_{k+1}) - W(t_k)\big) \right] = \sum_k E\Big[ G_k \big(W(t_{k+1}) - W(t_k)\big) \Big]
        $$

        Since $G(t)$ is adapted, $G_k$ is $\mathcal{F}(t_k)$-measurable; and Brownian motion has independent increments, so the increment $\Delta W_k = W(t_{k+1}) - W(t_k)$ is **completely independent** of $\mathcal{F}(t_k)$. Using the property of expectations for independent products (or conditional expectation):

        $$
        E[G_k \Delta W_k] = E[G_k] \cdot E[\Delta W_k] = E[G_k] \cdot 0 = 0
        $$

        Therefore, the expectation of the sum is 0.

        **2. Proving Itô Isometry**:
        Expand the expectation of the squared integral (using Fubini's theorem to interchange expectation and summation):

        $$
        E\left[ \left(\sum_k G_k \Delta W_k \right)^2 \right] = E\left[ \sum_{k} \sum_{j} G_k G_j \Delta W_k \Delta W_j \right]
        $$

        Split the double sum into three parts: $k > j$, $k < j$, and $k = j$.

        **Analyzing cross-terms ($k \neq j$)**: Without loss of generality, assume $k > j$. Then the time points satisfy $t_j < t_{j+1} \le t_k < t_{k+1}$.
        Among the four random variables $G_k, G_j, \Delta W_j, \Delta W_k$, the first three belong entirely to the historical information $\mathcal{F}(t_k)$, while the last increment $\Delta W_k$ occurs after $t_k$ and is independent of $\mathcal{F}(t_k)$.
        Using the Tower Property, first take conditional expectation with respect to $\mathcal{F}(t_k)$:

        $$
        E\Big[ G_k G_j \Delta W_j \Delta W_k \Big] = E\Big[ E\big[ G_k G_j \Delta W_j \Delta W_k \mid \mathcal{F}(t_k) \big] \Big]
        $$

        Since $G_k, G_j, \Delta W_j$ are known, they can be factored out:

        $$
        = E\Big[ G_k G_j \Delta W_j \underbrace{ E[\Delta W_k \mid \mathcal{F}(t_k)] }_{= 0} \Big] = 0
        $$

        Therefore, the expectation of all cross-terms is 0.

        **Analyzing diagonal terms ($k = j$)**:
        Only the squared terms on the diagonal remain:

        $$
        \sum_{k} E\Big[ G_k^2 (\Delta W_k)^2 \Big]
        $$

        Again using conditional expectation, factor out the square of $\Delta W_k$:

        $$
        = \sum_{k} E\Big[ E\big[ G_k^2 (\Delta W_k)^2 \mid \mathcal{F}(t_k) \big] \Big] = \sum_k E\Big[ G_k^2 E\big[ (\Delta W_k)^2 \mid \mathcal{F}(t_k) \big] \Big]
        $$

        Since the increments are independent and have variance $\Delta t_k = t_{k+1} - t_k$, we have $E[(\Delta W_k)^2 | \mathcal{F}(t_k)] = \Delta t_k$. Substituting yields:

        $$
        = \sum_k E[G_k^2] (t_{k+1} - t_k)
        $$

        When $|P| \to 0$, this Riemann sum directly converges to the Riemann integral $\int_0^T E[G(t)^2] dt$. The isometry is proved! $\square$

        > Approximation by simple functions in the $L^2(0, T)$ space suffices.
---

## 2. Indefinite Integrals and Continuous Martingale Properties

If we replace the upper limit of integration $T$ with a variable $t$, we obtain the **Indefinite Integral** of the stochastic process. This essentially defines a stochastic process.

!!! abstract "Definition: Itô Indefinite Integral"

    Let $G \in L^2(0, T)$. Its indefinite integral is defined as the stochastic process $I(t)$:

    $$
    I(t) = \int_0^t G(s) dW(s), \quad 0 \le t \le T
    $$

    The initial condition is obviously $I(0) = 0$.

!!! success "Theorem: The Itô Integral is a Continuous Square-Integrable Martingale"

    The indefinite integral process $\{I(t)\}_{t \ge 0}$ defined by the Itô integral possesses exceptionally elegant mathematical properties:
    It not only has **continuous sample paths** almost surely, but is also a **Martingale** with respect to the natural filtration.

    ??? proof "Proof of Martingale Property (Click to expand)"

        For any $0 \le s \le t \le T$, we need to prove $E[I(t) \mid \mathcal{F}(s)] = I(s)$ a.s.

        Split the interval $[0, t]$ at the point $s$:

        $$
        I(t) = \int_0^s G(\tau) dW(\tau) + \int_s^t G(\tau) dW(\tau) = I(s) + \int_s^t G(\tau) dW(\tau)
        $$

        Taking conditional expectation with respect to $\mathcal{F}(s)$ on both sides:

        $$
        E[I(t) \mid \mathcal{F}(s)] = E\left[ I(s) + \int_s^t G(\tau) dW(\tau) \bigg| \mathcal{F}(s) \right]
        $$

        Since the integration domain of $I(s)$ lies within $[0, s]$, it is clearly $\mathcal{F}(s)$-measurable, hence treated as a known constant (pulled out directly).
        For the latter part, use the conditional version of the property that the expectation of an Itô integral is 0:

        $$
        = I(s) + E\left[ \int_s^t G(\tau) dW(\tau) \bigg| \mathcal{F}(s) \right] = I(s) + 0 = I(s)
        $$

        Therefore, $I(t)$ is a martingale.

        *(Note: A rigorous proof of continuity requires the use of Doob's maximal inequality and the Borel-Cantelli lemma to construct uniform convergence of an $L^2$ Cauchy sequence. The manuscript mentions this, and the idea is very similar to the construction of Brownian motion. The intricate analytical details are omitted here.)* $\square$

---

## 3. Itô Processes and the Product Rule (Integration by Parts)

With the integral defined, the next step is to study its differential form.

!!! abstract "Definition: Itô Process and SDE"

    Let $F(t) \in L^1(0,T)$ and $G(t) \in L^2(0,T)$ be adapted processes. Define the stochastic process $X(t)$ as:

    $$
    X(t) = X(0) + \int_0^t F(s) ds + \int_0^t G(s) dW(s)
    $$

    This is often written in the intuitive form of a **stochastic differential equation (SDE)**:

    $$
    dX(t) = F(t) dt + G(t) dW(t)
    $$

    Here, $F(t)dt$ is called the **drift term**, representing the deterministic trend, and $G(t)dW(t)$ is called the **diffusion term**, representing the stochastic perturbation.

In classical Leibniz calculus, $d(X_1 X_2) = X_1 dX_2 + X_2 dX_1$. However, in stochastic analysis, because the quadratic variation of Brownian motion is non-zero, we must introduce an additional **quadratic correction term**.

!!! info "Theorem: Itô Product Rule"

    Consider two Itô processes $dX_i = F_i dt + G_i dW \quad (i=1,2)$.
    Then their product $X_1(t)X_2(t)$ is also an Itô process and satisfies:

    $$
    d(X_1 X_2) = X_1 dX_2 + X_2 dX_1 + dX_1 dX_2
    $$

    To expand $dX_1 dX_2$, we need to apply the famous **Itô Multiplication Table**:

    | $\times$ | $dt$ | $dW$ |
    | :---: | :---: | :---: |
    | **$dt$** | 0 | 0 |
    | **$dW$** | 0 | $dt$ |

    > *(Key intuition: Since $dW \sim \sqrt{dt}$, we have $(dW)^2 = dt$. Terms involving $dt$ raised to a power greater than one vanish as higher-order infinitesimals in the limit.)*

    **Direct expansion and application**:

    $$
    dX_1 dX_2 = (F_1 dt + G_1 dW)(F_2 dt + G_2 dW) = G_1 G_2 (dW)^2 = G_1 G_2 dt
    $$

    Therefore, the fully expanded form is:

    $$
    d(X_1 X_2) = (X_1 F_2 + X_2 F_1 + G_1 G_2) dt + (X_1 G_2 + X_2 G_1) dW
    $$

    > **Example 1:** $d(W^2) = W dW + W dW + (dW)^2 = 2W dW + dt$
    > (From this, the integral form follows immediately: $\int_0^t W dW = \frac{1}{2}W^2(t) - \frac{1}{2}t$, which matches the Riemann sum limit from the previous chapter!)
    >
    > **Example 2:** $d(tW) = t dW + W dt + dt \cdot dW = t dW + W dt$

---

## 4. Itô's Formula in Stochastic Calculus

Itô's Formula is the "Chain Rule" of stochastic calculus, which can be used to solve equations for various non-random phenomena.

!!! success "Theorem: Itô's Formula"

    Let $X(t)$ be an Itô process, $dX(t) = F dt + G dW$.
    Let $U(x, t): \mathbb{R} \times [0, T] \to \mathbb{R}$ be a deterministic function that is twice continuously differentiable ($C^{2,1}$).
    Define a new stochastic process $Y(t) = U(X(t), t)$. Then the stochastic differential of $Y(t)$ satisfies the truncated form of the Taylor expansion up to the second order:

    $$
    dU(X(t), t) = \frac{\partial U}{\partial t} dt + \frac{\partial U}{\partial x} dX + \frac{1}{2} \frac{\partial^2 U}{\partial x^2} (dX)^2
    $$

    Substituting the expression for $dX$ and the Itô multiplication table $(dX)^2 = G^2 dt$ and rearranging, we obtain the standard **Itô's Formula**:

    $$
    dU(X(t), t) = \left( \frac{\partial U}{\partial t} + \frac{\partial U}{\partial x} F + \frac{1}{2} \frac{\partial^2 U}{\partial x^2} G^2 \right) dt + \frac{\partial U}{\partial x} G dW
    $$

    ??? proof "Example Application of Itô's Formula to Polynomials (Click to Expand)"

        Derive and verify for the polynomial $f(x) = x^m$,
        using Itô's Formula:

        $$
        df(X) = f'(X) dX + \frac{1}{2} f''(X) (dX)^2
        $$

        Substituting the derivatives $f'(x) = m x^{m-1}$ and $f''(x) = m(m-1) x^{m-2}$:

        $$
        d(X^m) = m X^{m-1} dX + \frac{1}{2} m(m-1) X^{m-2} G^2 dt
        $$

        This is entirely consistent with the conclusion obtained by stepwise recursion using the product rule $d(x \cdot x^{k-1})$ from earlier, confirming the self-consistency of the operator algebra system in stochastic calculus.

---

## 5. Fokker-Planck Equation

With Itô's formula, we can not only study the microscopic single stochastic trajectory but also directly derive the deterministic law for the macroscopic **probability density function**. This is the Fokker-Planck equation, also known as the Kolmogorov forward equation (KFE).

!!! abstract "Derivation: Fokker-Planck Equation"

    Consider a general diffusion process (SDE):

    $$
    dX(t) = b(X, t) dt + \sigma(X, t) dW(t)
    $$

    Let the probability density function of the system being in state $x$ at time $t$ be $p(x, t)$.
    We take an arbitrary bounded, twice-differentiable **test function** $\phi(x)$ that vanishes at infinity.
    Expand $d\phi(X(t))$ using Itô's formula:

    $$
    d\phi(X(t)) = \phi'(X) dX + \frac{1}{2} \phi''(X) (dX)^2 = \phi'(X) \big[b(X, t) dt + \sigma(X, t) dW\big] + \frac{1}{2} \phi''(X) \sigma^2(X, t) dt
    $$

    **Take the expectation** on both sides (the diffusion term $\int \sigma dW$ has martingale property, so its expectation vanishes):

    $$
    E[\phi(X(t))] = E[\phi(X(0))] + \int_0^t E\left[ \phi'(X) b(X, s) + \frac{1}{2} \phi''(X) \sigma^2(X, s) \right] ds
    $$

    Using the density function $p(x, t)$, express the expectation as an integral over space $x$:

    $$
    \int_{\mathbb{R}} \phi(x) p(x, t) dx = \int_{\mathbb{R}} \phi(x) p(x, 0) dx + \int_0^t \int_{\mathbb{R}} \left( \phi'(x) b(x, s) + \frac{1}{2} \phi''(x) \sigma^2(x, s) \right) p(x, s) dx ds
    $$

    Differentiate both sides with respect to time $t$:

    $$
    \int_{\mathbb{R}} \phi(x) \frac{\partial p}{\partial t} dx = \int_{\mathbb{R}} \left( \phi'(x) b(x, t) + \frac{1}{2} \phi''(x) \sigma^2(x, t) \right) p(x, t) dx
    $$

    **The key step: Integration by parts twice.**
    To transfer the derivative operators from the test function $\phi$ to the density function $p$, we perform integration by parts on the right-hand side (assuming $\phi$ and its derivatives decay to 0 at the boundaries):
    
    First term:

    $$
    \int_{\mathbb{R}} \phi'(x) \big(b(x, t) p(x, t)\big) dx = -\int_{\mathbb{R}} \phi(x) \frac{\partial}{\partial x}\big(b(x, t) p(x, t)\big) dx
    $$

    Second term (integration by parts twice consecutively):

    $$
    \frac{1}{2} \int_{\mathbb{R}} \phi''(x) \big(\sigma^2(x, t) p(x, t)\big) dx = \frac{1}{2} \int_{\mathbb{R}} \phi(x) \frac{\partial^2}{\partial x^2}\big(\sigma^2(x, t) p(x, t)\big) dx
    $$

    After combining, since the equation holds for **any** test function $\phi(x)$, the integrands themselves must be equal. This yields the **Fokker-Planck equation**:

    $$
    \frac{\partial p}{\partial t} = - \frac{\partial}{\partial x} \big( b(x, t) p(x, t) \big) + \frac{1}{2} \frac{\partial^2}{\partial x^2} \big( \sigma^2(x, t) p(x, t) \big)
    $$