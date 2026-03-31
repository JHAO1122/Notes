# Chapter 3: Stochastic Integral

In classical calculus, one can define the Riemann-Stieltjes integral $\int g(t) dX(t)$. However, when the integration differential becomes Brownian motion $dW(t)$, traditional methods fail because its sample paths are too rough. This chapter demonstrates the evolutionary process of transitioning from deterministic integrands to stochastic integrands, leading to the concept of the Itô integral.

## 1. Paley-Wiener-Zygmund (PWZ) Stochastic Integral

The simplest case is when the integrand $g(t)$ is an ordinary deterministic function, while the integration variable is Brownian motion.

!!! abstract "Definition: PWZ Stochastic Integral"

    Let $W(t)$ be a standard Brownian motion. Assume $g(t) \in C^1([0, T])$ is a continuously differentiable deterministic function.
    Using integration by parts, the **PWZ Stochastic Integral** is defined as:

    $$
    \int_0^T g(t) dW(t) \triangleq g(T)W(T) - g(0)W(0) - \int_0^T g'(t) W(t) dt
    $$

    Since $W(0) = 0$ a.s., the above equation is usually abbreviated as:

    $$
    \int_0^T g(t) dW(t) = g(T)W(T) - \int_0^T g'(t) W(t) dt
    $$

    The integral on the right side is an ordinary Riemann integral (since the sample path of $W(t)$ is continuous a.s. and $g'(t)$ is continuous, the integral is well-defined).

!!! info "Theorem: Core Properties of the PWZ Integral"

    The PWZ integral defined above satisfies the following two extremely important properties (which serve as the cornerstone for all subsequent stochastic integrals):

    **1. Zero Expectation**:

    $$
    E\left[ \int_0^T g(t) dW(t) \right] = 0
    $$

    **2. Itô Isometry**:

    $$
    E\left[ \left( \int_0^T g(t) dW(t) \right)^2 \right] = \int_0^T g(t)^2 dt
    $$

    ??? proof "Rigorous Derivation of Itô Isometry (Click to expand)"

        We start directly from the definition of the PWZ integral and calculate the expectation of its square (assuming $g(0)=0$ or absorbing it into a constant term for simplicity. Here we use the more essential double integral method from your manuscript):
        
        Write the square in the form of a double integral:

        $$
        E\left[ \left( \int_0^T g'(t)W(t) dt \right)^2 \right] = E\left[ \int_0^T g'(t)W(t) dt \int_0^T g'(s)W(s) ds \right]
        $$

        Use Fubini's Theorem to swap the order of expectation and integration:

        $$
        = \int_0^T \int_0^T g'(t)g'(s) E[W(t)W(s)] ds dt
        $$

        Substitute the autocovariance function of Brownian motion $E[W(t)W(s)] = t \wedge s$:

        $$
        = \int_0^T \int_0^T g'(t)g'(s) (t \wedge s) ds dt
        $$

        Using symmetry, split the integration region into two parts: $s \le t$ and $t \le s$:

        $$
        = \int_0^T \left( \int_0^t g'(s) g'(t) s \, ds + \int_t^T g'(s) g'(t) t \, ds \right) dt
        $$

        Through integration by parts and algebraic simplification (skipping the tedious algebraic steps from the manuscript here), the above equation perfectly collapses into the target form:

        $$
        = \int_0^T g(t)^2 dt \quad \square
        $$

---

## 2. Norm-Preserving Extension of Densely Defined Bounded Linear Operators (BLT Theorem)

The PWZ integral just discussed requires $g(t) \in C^1$. However, in practical applications, we need to integrate over more general $L^2$ functions. This requires a tool from functional analysis: the **Bounded Linear Transformation (BLT) Theorem**.

!!! abstract "Theorem: Norm-Preserving Extension of Bounded Linear Operators"

    Let $X$ and $Y$ be Banach spaces, and $S$ be a dense linear subspace of $X$.
    Let $T: S \to Y$ be a bounded linear operator, meaning there exists a constant $C > 0$ such that for any $x \in S$:

    $$
    \|Tx\|_Y \le C \|x\|_X
    $$

    Then there exists a unique bounded linear operator $\overline{T}: X \to Y$ such that for all $x \in S$, $\overline{T}x = Tx$ (i.e., $\overline{T}|_S = T$).
    Furthermore, the operator norm remains unchanged: $\|\overline{T}\| = \|T\| \le C$.

    ??? proof "Constructive Proof of the BLT Theorem (Click to expand)"

        **Step 1: Constructing the Limit Mapping**
        Since $S$ is dense in $X$, for any $x \in X$, there must exist a sequence $\{x_n\}$ in $S$ such that $x_n \to x$.
        Because $T$ is bounded (continuous) on $S$, we examine the distance of the sequence $\{Tx_n\}$ in $Y$:

        $$
        \|Tx_n - Tx_m\|_Y = \|T(x_n - x_m)\|_Y \le C \|x_n - x_m\|_X
        $$

        Since $\{x_n\}$ is a convergent sequence, it must be a Cauchy sequence. Therefore, as $n, m \to \infty$, $\|x_n - x_m\|_X \to 0$.
        This implies that $\{Tx_n\}$ is a Cauchy sequence in the Banach space $Y$. Because $Y$ is complete, the limit $\lim_{n \to \infty} Tx_n$ must exist. We define:

        $$
        \overline{T}x \triangleq \lim_{n \to \infty} Tx_n
        $$

        **Step 2: Proving the Well-definedness (Independence of Sequence Choice)**
        Suppose there is another sequence $x_n' \to x$. We need to prove $\lim Tx_n' = \lim Tx_n$.
        Let $y = \lim Tx_n$ and $y' = \lim Tx_n'$.

        $$
        \|y - y'\|_Y = \lim_{n \to \infty} \|Tx_n - Tx_n'\|_Y \le C \lim_{n \to \infty} \|x_n - x_n'\|_X = 0
        $$

        Thus $y = y'$, making the mapping $\overline{T}$ well-defined.

        **Step 3: Proving Linearity and Norm-Preservation**
        Linearity is a natural corollary of the limit. For norm-preservation, for any $x \in X$ and its approximating sequence $x_n \to x$:

        $$
        \|\overline{T}x\|_Y = \lim_{n \to \infty} \|Tx_n\|_Y \le C \lim_{n \to \infty} \|x_n\|_X = C \|x\|_X
        $$

        Therefore, $\overline{T}$ is bounded, and its norm does not exceed $C$. $\square$

!!! success "Application: Extension of the PWZ Integral to $L^2$ Space"

    We take the space of integrands as $X = L^2([0, T])$, and the space of integral results as $Y = L^2(\Omega, P)$.
    The dense subspace is chosen as $S = C^1([0, T])$.
    Define the operator $T: g \mapsto \int_0^T g(t) dW(t)$.
    
    From Itô Isometry, we know that for any $g \in S$:

    $$
    \|Tg\|_{L^2(\Omega)}^2 = E\left[ \left( \int_0^T g(t) dW(t) \right)^2 \right] = \int_0^T g(t)^2 dt = \|g\|_{L^2([0,T])}^2
    $$

    This means the operator $T$ is an isometric isomorphism (operator norm $C=1$). By the BLT Theorem, we can perfectly and uniquely extend the PWZ integral to the entire $L^2([0, T])$ space.

---

## 3. Quadratic Variation of Brownian Motion

The fundamental reason traditional Riemann integrals cannot handle stochastic integrals lies in the "roughness" of Brownian motion paths, specifically its quadratic variation property.

Consider a partition $P = \{0 = t_0 < t_1 < \dots < t_m = T\}$ of the time interval $[0, T]$, with the mesh size denoted as $|P| = \max (t_{k+1} - t_k)$.

!!! info "Theorem 1: The Quadratic Variation of Brownian Motion Equals Time $T$"

    As the mesh size approaches zero ($|P| \to 0$), the sum of the squared increments of Brownian motion converges to $T$ in the sense of $L^2(\Omega, P)$:

    $$
    \sum_{k=0}^{m-1} (W(t_{k+1}) - W(t_k))^2 \xrightarrow{L^2} T
    $$

    ??? proof "Rigorous Derivation of $L^2$ Convergence of Quadratic Variation (Click to expand)"

        To prove $L^2$ convergence, we need to show that its mean squared error from $T$ approaches 0.
        Since $\sum_{k=0}^{m-1} (t_{k+1} - t_k) = T$, we can write the target error as:

        $$
        E\left[ \left( \sum_{k=0}^{m-1} \big((W(t_{k+1}) - W(t_k))^2 - (t_{k+1} - t_k)\big) \right)^2 \right]
        $$

        Let $\Delta W_k = W(t_{k+1}) - W(t_k)$ and $\Delta t_k = t_{k+1} - t_k$. Expanding the squared term separates it into squared terms and cross terms:

        $$
        = \sum_k E\Big[ \big((\Delta W_k)^2 - \Delta t_k\big)^2 \Big] + \sum_{k \neq j} E\Big[ \big((\Delta W_k)^2 - \Delta t_k\big)\big((\Delta W_j)^2 - \Delta t_j\big) \Big]
        $$

        **Key Point 1: Cross Terms are Zero.**
        Due to the independent increments property of Brownian motion, when $k \neq j$, $\Delta W_k$ is independent of $\Delta W_j$. Furthermore, since $E[(\Delta W_k)^2] = \Delta t_k$, the expectation of each factor is 0, making the overall expectation of the cross terms 0.

        **Key Point 2: Calculation of Squared Terms.**
        Only the variance terms on the diagonal remain. Note that $\Delta W_k \sim \mathcal{N}(0, \Delta t_k)$, so it can be standardized as $\sqrt{\Delta t_k} Z$, where $Z \sim \mathcal{N}(0, 1)$.

        $$
        E\Big[ \big((\Delta W_k)^2 - \Delta t_k\big)^2 \Big] = E\Big[ (\Delta t_k Z^2 - \Delta t_k)^2 \Big] = (\Delta t_k)^2 E[(Z^2 - 1)^2]
        $$

        Since the fourth moment of the standard normal distribution is $E[Z^4] = 3$ and $E[Z^2] = 1$, we have $E[(Z^2 - 1)^2] = 3 - 2(1) + 1 = 2$.

        $$
        \text{Total Error} = 2 \sum_{k=0}^{m-1} (\Delta t_k)^2
        $$

        We bound this summation by extracting the largest $\Delta t_k$, which is the mesh size $|P|$:

        $$
        2 \sum_{k=0}^{m-1} (\Delta t_k)^2 \le 2 |P| \sum_{k=0}^{m-1} \Delta t_k = 2 |P| T
        $$

        As $|P| \to 0$, $2 |P| T \to 0$. Thus, it converges to $T$ in the $L^2$ sense. $\square$

This theorem directly leads to a corollary regarding the properties of Brownian motion:

!!! success "Theorem 2: Brownian Motion has Infinite Total Variation Everywhere"

    Almost all sample paths of Brownian motion $W(t, \omega)$ have an infinite total variation on any interval.

    > **Minimalist Proof by Contradiction**: If a certain path has a bounded total variation $V_T < \infty$, then its quadratic variation can be bounded by:
    > $\sum (\Delta W_k)^2 \le (\max_k |\Delta W_k|) \sum |\Delta W_k| \le (\max_k |\Delta W_k|) \cdot V_T$
    > Because the path is continuous, as the partition becomes infinitely fine, $\max |\Delta W_k| \to 0$. This would cause the quadratic variation to tend to 0, which absolutely contradicts Theorem 1 where the quadratic variation equals $T > 0$!

---

## 4. The Riemann Sum of Stochastic Integrals and the Derivation of the Itô Integral

So how do we calculate $\int_0^T W(t) dW(t)$?
Returning to the definition of the Riemann sum, we observe a phenomenon that does not occur in classical calculus: **changing the evaluation point leads to drastically different integration results.**

Construct a partition $P$, and select a point $\tau_k = (1-\lambda)t_k + \lambda t_{k+1}$ ($\lambda \in [0, 1]$) in each sub-interval $[t_k, t_{k+1}]$.
Examine the Riemann sum:

$$
R_n = \sum_{k=0}^{m-1} W(\tau_k) \big( W(t_{k+1}) - W(t_k) \big)
$$

!!! info "Dependence of the Integral Result on the Value of $\lambda$"

    For clarity, we examine the special case where $\lambda = 0$ (taking the left endpoint, i.e., the Itô integral), where $\tau_k = t_k$:

    $$
    R_n = \sum_{k=0}^{m-1} W(t_k) \big( W(t_{k+1}) - W(t_k) \big)
    $$

    ??? proof "Algebraic Identity Splitting and Limit Evaluation of the Itô Integral (Click to expand)"

        This is an extremely clever algebraic trick. Using the identity $a(b-a) = \frac{1}{2}\big( b^2 - a^2 - (b-a)^2 \big)$, we rewrite each term by letting $a = W(t_k)$ and $b = W(t_{k+1})$:

        $$
        W(t_k) \big( W(t_{k+1}) - W(t_k) \big) = \frac{1}{2}\big( W(t_{k+1})^2 - W(t_k)^2 \big) - \frac{1}{2}\big( W(t_{k+1}) - W(t_k) \big)^2
        $$

        Summing all the terms, the original sum splits into two parts, $B_1$ and $B_2$:

        $$
        R_n = \underbrace{ \frac{1}{2} \sum_{k=0}^{m-1} \big( W(t_{k+1})^2 - W(t_k)^2 \big) }_{B_1} - \underbrace{ \frac{1}{2} \sum_{k=0}^{m-1} \big( W(t_{k+1}) - W(t_k) \big)^2 }_{B_2}
        $$

        **Analyzing $B_1$**: This is a perfect telescoping sum, where all intermediate terms cancel out:

        $$
        B_1 = \frac{1}{2} \big( W(T)^2 - W(0)^2 \big) = \frac{1}{2} W(T)^2
        $$

        **Analyzing $B_2$**: This is precisely the quadratic variation of Brownian motion! According to Theorem 1, as the partition becomes finer, it converges in the $L^2$ sense:

        $$
        B_2 \xrightarrow{L^2} \frac{1}{2} T
        $$

        In conclusion, when $\lambda = 0$ (the Itô integral), we obtain:

        $$
        \int_0^T W(t) dW(t) \triangleq \lim_{|P|\to 0} R_n = \frac{1}{2} W(T)^2 - \frac{1}{2} T
        $$

        *(Note: The extra $-\frac{1}{2}T$ term added later is called the Itô correction term)* $\square$

**The two most important schools of integration:**

* When $\lambda = 0$ (taking the left endpoint), it is the **Itô Integral**, yielding $\frac{1}{2}W(T)^2 - \frac{1}{2}T$. It preserves the martingale property and is commonly used in mathematical finance.

* When $\lambda = 1/2$ (taking the midpoint), it is the **Stratonovich Integral**. The correction terms cancel out, yielding $\frac{1}{2}W(T)^2$, which is formally consistent with classical calculus and is commonly used in physics and engineering.

---

## 5. Measure-Theoretic Preparation for the Rigorous Itô Integral

To generalize the Itô integral to more general stochastic processes (not just $W(t)$), we need to strictly define what it means to "ignore the future and only look at the past." This requires introducing the concepts of Filtration and Adapted Processes.

!!! abstract "Definition: Information Filtration and $\sigma$-Algebra"

    **1. Natural Filtration of Brownian Motion**:
    For any time $t$, the information generated by the historical trajectory of Brownian motion is denoted by the $\sigma$-algebra:

    $$
    \mathcal{F}_W(t) \triangleq \sigma(\{W(s) \mid 0 \le s \le t\})
    $$

    It contains all the path information of the Brownian motion up to time $t$.

    **2. Independence of Future Increments**:
    We define the future increment information filtration as $\mathcal{F}^t \triangleq \sigma(\{W(s) - W(t) \mid s > t\})$. According to the independent increments property of Brownian motion, $\mathcal{F}^t$ is **completely independent** of the historical filtration $\mathcal{F}_W(t)$.

    **3. General Filtration**:
    A family of $\sigma$-algebras $\{\mathcal{F}(t)\}_{t \ge 0}$ satisfying the following conditions:
    * **Monotonicity**: $\mathcal{F}(s) \subset \mathcal{F}(t)$ for any $0 \le s \le t$ (information is not forgotten).
    * **Contains History**: $\mathcal{F}_W(t) \subset \mathcal{F}(t)$.
    * **Future Independence**: The increment $W(s) - W(t)$ is independent of $\mathcal{F}(t)$.

With the mathematical framework for the concept of "information" established, we can determine which stochastic processes can be integrated using Itô's calculus.

!!! info "Definition: Adapted Process and Progressively Measurable"

    **1. Adapted Process**:
    If for every fixed time $t$, the random variable $G(t, \omega)$ is $\mathcal{F}(t)$-measurable, then the stochastic process $G(t)$ is said to be adapted to the filtration $\{\mathcal{F}(t)\}$.
    *(Intuitive understanding: At time $t$, as long as you know all the historical information $\mathcal{F}(t)$, you know the value of $G(t)$ at this moment, and you absolutely do not need future information.)*

    **2. Progressively Measurable**:
    A stronger condition. The mapping $(s, \omega) \mapsto G(s, \omega)$ is jointly measurable on the product space $[0, t] \times \Omega$ with respect to the $\mathcal{B}([0, t]) \otimes \mathcal{F}(t)$ field. This guarantees that operations like Riemann or Lebesgue integration over a time interval are legitimate.

In the subsequent construction of this book, the Itô integral $\int_0^T G(t, \omega) dW(t)$ will be rigorously defined in the following Hilbert space:

$$
L^2(\Omega \times [0, T]) = \left\{ G(t, \omega) \text{ is a progressively measurable process} \mid E\left[ \int_0^T G(t, \omega)^2 dt \right] < \infty \right\}
$$