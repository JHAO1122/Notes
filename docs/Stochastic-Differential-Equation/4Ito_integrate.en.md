# Chapter 4: Itô Integral & Stochastic Differential

In the previous chapter, we extended the integral to the $L^2$ space. Now, we formally define the Itô integral for **stochastic processes** and establish a mathematical system distinct from classical calculus—**Itô Calculus**.

## 1. The Itô Integral of Adapted Processes and Itô Isometry

For a stochastic process $G(t, \omega)$, due to the rough paths of Brownian motion $W(t)$, we cannot employ the pointwise-defined Riemann-Stieltjes integral. Instead, we must define it based on the left endpoints of a time grid (without knowledge of the future).

!!! abstract "Definition: The Itô Integral"

    Let $\{W(t)\}_{t \ge 0}$ be a standard Brownian motion, and $\{\mathcal{F}(t)\}$ be its filtration.
    Assume the stochastic process $G(t) \in L^2(0, T)$ and $G(t)$ is an **Adapted Process** (i.e., its value at time $t$ depends only on the historical information $\mathcal{F}(t)$).

    For a partition $P = \{0 = t_0 < t_1 < \dots < t_m = T\}$ of the interval $[0, T]$, define its Riemann sum as:

    $$
    R_m = \sum_{k=0}^{m-1} G(t_k) \big( W(t_{k+1}) - W(t_k) \big)
    $$

    As the mesh size $|P| \to 0$, if this sum converges in the sense of $L^2(\Omega, P)$, its limit is defined as the **Itô integral** of $G(t)$ with respect to Brownian motion:

    $$
    \int_0^T G(t) dW(t) = \lim_{|P| \to 0} \sum_{k=0}^{m-1} G(t_k) \big( W(t_{k+1}) - W(t_k) \big)
    $$

After defining this integral using the approximation method of simple functions (step processes), it possesses the following properties:

!!! info "Theorem 1: Core Properties of the Itô Integral"

    Let $G(t), H(t) \in L^2(0, T)$ both be adapted processes, and $a, b \in \mathbb{R}$ be constants.

    **(1) Linearity**:

    $$
    \int_0^T (aG(t) + bH(t)) dW(t) = a \int_0^T G(t) dW(t) + b \int_0^T H(t) dW(t) \quad a.s.
    $$

    **(2) Zero Expectation**:

    $$
    E\left[ \int_0^T G(t) dW(t) \right] = 0
    $$

    **(3) Itô Isometry**:

    $$
    E\left[ \left( \int_0^T G(t) dW(t) \right)^2 \right] = \int_0^T E[G(t)^2] dt
    $$

    ??? proof "Rigorous Proof of Zero Expectation and Itô Isometry (Click to expand)"

        We use simple functions (step functions) to prove these properties. Let $G_k = G(t_k)$.

        **1. Proof of Zero Expectation**:

        $$
        E\left[ \sum_{k} G_k \big(W(t_{k+1}) - W(t_k)\big) \right] = \sum_k E\Big[ G_k \big(W(t_{k+1}) - W(t_k)\big) \Big]
        $$

        Since $G(t)$ is an adapted process, $G_k$ is measurable with respect to $\mathcal{F}(t_k)$; meanwhile, Brownian motion has independent increments, so the increment $\Delta W_k = W(t_{k+1}) - W(t_k)$ is **completely independent** of $\mathcal{F}(t_k)$. Using the independent multiplication property of expectations (or conditional expectation):

        $$
        E[G_k \Delta W_k] = E[G_k] \cdot E[\Delta W_k] = E[G_k] \cdot 0 = 0
        $$

        Therefore, the expectation of the sum is 0.

        **2. Proof of Itô Isometry**:
        
        Expanding the expectation of the squared integral (using Fubini's Theorem to swap expectation and summation):

        $$
        E\left[ \left(\sum_k G_k \Delta W_k \right)^2 \right] = E\left[ \sum_{k} \sum_{j} G_k G_j \Delta W_k \Delta W_j \right]
        $$

        Split the double summation into three parts: $k > j$, $k < j$, and $k = j$.
        
        **Analyzing the cross terms ($k \neq j$)**: Without loss of generality, assume $k > j$. The time points satisfy $t_j < t_{j+1} \le t_k < t_{k+1}$.
        Among the four random variables $G_k, G_j, \Delta W_j, \Delta W_k$, the first three belong entirely to the historical information $\mathcal{F}(t_k)$, while the last increment $\Delta W_k$ occurs after $t_k$ and is independent of $\mathcal{F}(t_k)$.
        Using the Tower Property, we first take the conditional expectation with respect to $\mathcal{F}(t_k)$:

        $$
        E\Big[ G_k G_j \Delta W_j \Delta W_k \Big] = E\Big[ E\big[ G_k G_j \Delta W_j \Delta W_k \mid \mathcal{F}(t_k) \big] \Big]
        $$

        Since $G_k, G_j, \Delta W_j$ are known, they can be pulled out:

        $$
        = E\Big[ G_k G_j \Delta W_j \underbrace{ E[\Delta W_k \mid \mathcal{F}(t_k)] }_{= 0} \Big] = 0
        $$

        Thus, the expectations of all cross terms are exactly 0.

        **Analyzing the diagonal terms ($k = j$)**:
        
        Only the squared terms on the diagonal remain:

        $$
        \sum_{k} E\Big[ G_k^2 (\Delta W_k)^2 \Big]
        $$

        Again, using conditional expectation, we extract the squared term of $\Delta W_k$:

        $$
        = \sum_{k} E\Big[ E\big[ G_k^2 (\Delta W_k)^2 \mid \mathcal{F}(t_k) \big] \Big] = \sum_k E\Big[ G_k^2 E\big[ (\Delta W_k)^2 \mid \mathcal{F}(t_k) \big] \Big]
        $$

        Since the increments are independent and the variance is $\Delta t_k = t_{k+1} - t_k$, we have $E[(\Delta W_k)^2 | \mathcal{F}(t_k)] = \Delta t_k$. Substituting this yields:

        $$
        = \sum_k E[G_k^2] (t_{k+1} - t_k)
        $$

        As $|P| \to 0$, this Riemann sum converges directly to the Riemann integral $\int_0^T E[G(t)^2] dt$. Itô Isometry is proven! $\square$

        > Approximating to the $L^2(0, T)$ space with simple functions completes the argument.

---

## 2. Indefinite Integrals and Continuous Martingale Properties

If we replace the upper limit of the integral $T$ with a variable $t$, we obtain the **Indefinite Integral** of a stochastic process. This fundamentally defines a new stochastic process.

!!! abstract "Definition: Itô Indefinite Integral"

    Let $G \in L^2(0, T)$, define its indefinite integral as the stochastic process $I(t)$:

    $$
    I(t) = \int_0^t G(s) dW(s), \quad 0 \le t \le T
    $$

    Obviously, the initial condition is $I(0) = 0$.

!!! success "Theorem: The Itô Integral is a Continuous Square-Integrable Martingale"

    The indefinite integral process $\{I(t)\}_{t \ge 0}$ defined by the Itô integral has perfectly elegant mathematical properties:
    Not only does it almost surely possess **continuous sample paths**, but it is also a **Martingale** with respect to the natural filtration.

    ??? proof "Proof of Martingale Property (Click to expand)"

        For any $0 \le s \le t \le T$, we need to prove $E[I(t) \mid \mathcal{F}(s)] = I(s)$ a.s.

        Split the interval $[0, t]$ at point $s$:

        $$
        I(t) = \int_0^s G(\tau) dW(\tau) + \int_s^t G(\tau) dW(\tau) = I(s) + \int_s^t G(\tau) dW(\tau)
        $$

        Take the conditional expectation with respect to $\mathcal{F}(s)$ on both sides:

        $$
        E[I(t) \mid \mathcal{F}(s)] = E\left[ I(s) + \int_s^t G(\tau) dW(\tau) \bigg| \mathcal{F}(s) \right]
        $$

        Since the integration domain of $I(s)$ is within $[0, s]$, it is clearly $\mathcal{F}(s)$-measurable, so it is "known and therefore a constant" (pulled out directly).
        For the second part, using the conditional version of the zero expectation property of the Itô integral:

        $$
        = I(s) + E\left[ \int_s^t G(\tau) dW(\tau) \bigg| \mathcal{F}(s) \right] = I(s) + 0 = I(s)
        $$

        Therefore, $I(t)$ is a martingale.

        *(Note: The rigorous proof of continuity requires using Doob's maximal inequality and the Borel-Cantelli lemma to construct the uniform convergence of an $L^2$ Cauchy sequence. This is mentioned in the manuscript, and the idea is extremely similar to the construction of Brownian motion; the tedious analytic details are omitted here.)* $\square$

---

## 3. Itô Process and the Product Rule (Integration by Parts)

Having established the integral, the next step is to study its differential form.

!!! abstract "Definition: Itô Process and SDE"

    Let $F(t) \in L^1(0,T)$ and $G(t) \in L^2(0,T)$ both be adapted processes. Define the stochastic process $X(t)$:

    $$
    X(t) = X(0) + \int_0^t F(s) ds + \int_0^t G(s) dW(s)
    $$

    The above equation is typically written in the intuitive form of a **Stochastic Differential Equation (SDE)**:

    $$
    dX(t) = F(t) dt + G(t) dW(t)
    $$

    Here, $F(t)dt$ is called the **Drift** term, representing the deterministic trend; $G(t)dW(t)$ is called the **Diffusion** term, representing the stochastic perturbation.

In classical Leibniz calculus, $d(X_1 X_2) = X_1 dX_2 + X_2 dX_1$. However, in stochastic analysis, because the quadratic variation of Brownian motion is not 0, we must introduce an additional **quadratic correction term**.

!!! info "Theorem: Itô Product Rule"

    Assume we have two Itô processes $dX_i = F_i dt + G_i dW \quad (i=1,2)$.
    Then their product $X_1(t)X_2(t)$ is also an Itô process and satisfies:

    $$
    d(X_1 X_2) = X_1 dX_2 + X_2 dX_1 + dX_1 dX_2
    $$

    Here, we need to apply the following famous **Itô Multiplication Table** to expand $dX_1 dX_2$:

    | $\times$ | $dt$ | $dW$ |
    | :---: | :---: | :---: |
    | **$dt$** | 0 | 0 |
    | **$dW$** | 0 | $dt$ |

    > *(Core intuition: Because $dW \sim \sqrt{dt}$, we have $(dW)^2 = dt$. Terms containing powers of $dt$ higher than one become higher-order infinitesimals in the limit, meaning they are 0).*

    **Direct application of expansion**:

    $$
    dX_1 dX_2 = (F_1 dt + G_1 dW)(F_2 dt + G_2 dW) = G_1 G_2 (dW)^2 = G_1 G_2 dt
    $$

    Therefore, the fully expanded form is:

    $$
    d(X_1 X_2) = (X_1 F_2 + X_2 F_1 + G_1 G_2) dt + (X_1 G_2 + X_2 G_1) dW
    $$

    > **Example 1:** $d(W^2) = W dW + W dW + (dW)^2 = 2W dW + dt$
    > 
    > (From this, the integral form can be immediately derived: $\int_0^t W dW = \frac{1}{2}W^2(t) - \frac{1}{2}t$, which perfectly aligns with the Riemann sum limit from the previous chapter!)
    >
    > **Example 2:** $d(tW) = t dW + W dt + dt \cdot dW = t dW + W dt$

---

## 4. Itô's Formula in Stochastic Differential Calculus

Itô's formula is the "Chain Rule" of the entire stochastic calculus and can be used to solve equations for various stochastic phenomena.

!!! success "Theorem: Itô's Formula"

    Let $X(t)$ be an Itô process, $dX(t) = F dt + G dW$.
    Let $U(x, t): \mathbb{R} \times [0, T] \to \mathbb{R}$ be a twice continuously differentiable ($C^{2,1}$) deterministic function.
    Define a new stochastic process $Y(t) = U(X(t), t)$. Then the stochastic differential of $Y(t)$ satisfies the truncated Taylor expansion up to the second order:

    $$
    dU(X(t), t) = \frac{\partial U}{\partial t} dt + \frac{\partial U}{\partial x} dX + \frac{1}{2} \frac{\partial^2 U}{\partial x^2} (dX)^2
    $$

    Substituting the expression for $dX$ and the Itô multiplication table rule $(dX)^2 = G^2 dt$ and rearranging yields the standard **Itô's Formula**:

    $$
    dU(X(t), t) = \left( \frac{\partial U}{\partial t} + \frac{\partial U}{\partial x} F + \frac{1}{2} \frac{\partial^2 U}{\partial x^2} G^2 \right) dt + \frac{\partial U}{\partial x} G dW
    $$

    ??? proof "Example Application of Itô's Formula to Polynomials (Click to expand)"

        We derive and verify this for the polynomial $f(x) = x^m$.
        By Itô's formula:

        $$
        df(X) = f'(X) dX + \frac{1}{2} f''(X) (dX)^2
        $$

        Substitute the derivatives $f'(x) = m x^{m-1}$ and $f''(x) = m(m-1) x^{m-2}$:

        $$
        d(X^m) = m X^{m-1} dX + \frac{1}{2} m(m-1) X^{m-2} G^2 dt
        $$

        This is perfectly consistent with the conclusion derived step-by-step using the product rule $d(x \cdot x^{k-1})$ earlier, confirming the self-consistency of the algebraic system of stochastic calculus operators.

---

## 5. The Fokker-Planck Equation

With Itô's formula, we can study not only microscopic individual stochastic paths but also directly derive the **deterministic laws governing macroscopic probability density functions**. This is the Fokker-Planck equation, also known as the Kolmogorov Forward Equation (KFE).

!!! abstract "Derivation: Fokker-Planck Equation"

    Consider a general diffusion process (SDE):

    $$
    dX(t) = b(X, t) dt + \sigma(X, t) dW(t)
    $$

    Let $p(x, t)$ be the probability density function of the system being in state $x$ at time $t$.
    We take an arbitrary bounded, twice-differentiable **test function** $\phi(x)$ that vanishes at infinity.
    Expand $d\phi(X(t))$ according to Itô's formula:

    $$
    d\phi(X(t)) = \phi'(X) dX + \frac{1}{2} \phi''(X) (dX)^2 = \phi'(X) \big[b(X, t) dt + \sigma(X, t) dW\big] + \frac{1}{2} \phi''(X) \sigma^2(X, t) dt
    $$

    Take the **expectation** on both sides of the equation (the diffusion term $\int \sigma dW$ possesses the martingale property, so its expectation vanishes):

    $$
    E[\phi(X(t))] = E[\phi(X(0))] + \int_0^t E\left[ \phi'(X) b(X, s) + \frac{1}{2} \phi''(X) \sigma^2(X, s) \right] ds
    $$

    Using the density function $p(x, t)$, write the expectation in the form of a spatial integral over $x$:

    $$
    \int_{\mathbb{R}} \phi(x) p(x, t) dx = \int_{\mathbb{R}} \phi(x) p(x, 0) dx + \int_0^t \int_{\mathbb{R}} \left( \phi'(x) b(x, s) + \frac{1}{2} \phi''(x) \sigma^2(x, s) \right) p(x, s) dx ds
    $$

    Take the derivative with respect to time $t$ on both sides:

    $$
    \int_{\mathbb{R}} \phi(x) \frac{\partial p}{\partial t} dx = \int_{\mathbb{R}} \left( \phi'(x) b(x, t) + \frac{1}{2} \phi''(x) \sigma^2(x, t) \right) p(x, t) dx
    $$

    **The Key: Integration by Parts Twice.**
    
    To peel the derivative operator off the test function $\phi$ and transfer it to the density function $p$, we perform integration by parts on the right side (assuming $\phi$ and its derivatives decay to 0 at the boundaries):
    
    First term:

    $$
    \int_{\mathbb{R}} \phi'(x) \big(b(x, t) p(x, t)\big) dx = -\int_{\mathbb{R}} \phi(x) \frac{\partial}{\partial x}\big(b(x, t) p(x, t)\big) dx
    $$

    Second term (applying continuous integration by parts twice):

    $$
    \frac{1}{2} \int_{\mathbb{R}} \phi''(x) \big(\sigma^2(x, t) p(x, t)\big) dx = \frac{1}{2} \int_{\mathbb{R}} \phi(x) \frac{\partial^2}{\partial x^2}\big(\sigma^2(x, t) p(x, t)\big) dx
    $$

    After combining them, since the equation holds for **any** test function $\phi(x)$, the integrands themselves must be equal. This yields the **Fokker-Planck Equation**:

    $$
    \frac{\partial p}{\partial t} = - \frac{\partial}{\partial x} \big( b(x, t) p(x, t) \big) + \frac{1}{2} \frac{\partial^2}{\partial x^2} \big( \sigma^2(x, t) p(x, t) \big)
    $$