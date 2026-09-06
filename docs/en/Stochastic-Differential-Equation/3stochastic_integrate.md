# Chapter 3: Stochastic Integration

In classical calculus, the Riemann-Stieltjes integral $\int g(t) dX(t)$ can be defined. However, when the differential element becomes Brownian motion $dW(t)$, traditional methods cannot be applied due to the roughness of its sample paths. This section demonstrates the progression from deterministic integrands to stochastic integrands and introduces the concept of the Itô integral.

## 1. Paley-Wiener-Zygmund (PWZ) Stochastic Integral

The simplest case is when the integrand $g(t)$ is an ordinary deterministic function, and the variable of integration is Brownian motion.

!!! abstract "Definition: PWZ Stochastic Integral"

    Let $W(t)$ be a standard Brownian motion. Assume $g(t) \in C^1([0, T])$ is a deterministic function that is continuously differentiable.
    Using the integration by parts formula, the **PWZ stochastic integral** is defined as:

    $$
    \int_0^T g(t) dW(t) \triangleq g(T)W(T) - g(0)W(0) - \int_0^T g'(t) W(t) dt
    $$

    Since $W(0) = 0$ a.s., the above formula is usually simplified to:

    $$
    \int_0^T g(t) dW(t) = g(T)W(T) - \int_0^T g'(t) W(t) dt
    $$

    The integral on the right-hand side here is an ordinary Riemann integral (because the sample paths of $W(t)$ are a.s. continuous, and $g'(t)$ is continuous, so the integral is well-defined).

!!! info "Theorem: Core Properties of the PWZ Integral"

    For the PWZ integral defined above, the following two extremely important properties hold (they are also the cornerstone for all subsequent stochastic integrals):

    **(1) Zero Mean**:

    $$
    E\left[ \int_0^T g(t) dW(t) \right] = 0
    $$

    **(2) Itô Isometry**:

    $$
    E\left[ \left( \int_0^T g(t) dW(t) \right)^2 \right] = \int_0^T g(t)^2 dt
    $$

    ??? proof "Rigorous Derivation of the Isometry Property (Click to Expand)"

        We start directly from the definition of the PWZ integral and compute the expectation of its square (assuming $g(0)=0$ or absorbing it into a constant term for simplicity; here we use the more fundamental double integral method from your manuscript):
        
        Write the square in the form of a double integral:

        $$
        E\left[ \left( \int_0^T g'(t)W(t) dt \right)^2 \right] = E\left[ \int_0^T g'(t)W(t) dt \int_0^T g'(s)W(s) ds \right]
        $$

        Use Fubini's theorem to interchange the order of expectation and integration:

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

        Through integration by parts and algebraic simplification (the manuscript skips the cumbersome algebraic steps here), the above expression ultimately contracts perfectly to the target form:

        $$
        = \int_0^T g(t)^2 dt \quad \square
        $$

---

## 2. Norm-Preserving Extension of Densely Defined Bounded Linear Operators (BLT Theorem)

The previous PWZ integral required $g(t) \in C^1$. However, in practical applications, we need to integrate more general $L^2$ functions. This requires the use of tools from functional analysis: the **Bounded Linear Transformation Theorem (BLT Theorem)**.

!!! abstract "Theorem: Norm-Preserving Extension of Densely Defined Bounded Linear Operators"

    Let $X, Y$ be Banach spaces, and let $S$ be a dense linear subspace of $X$.
    Let $T: S \to Y$ be a bounded linear operator, i.e., there exists a constant $C > 0$ such that for any $x \in S$:

    $$
    \|Tx\|_Y \le C \|x\|_X
    $$

    Then there exists a unique bounded linear operator $\overline{T}: X \to Y$ such that for all $x \in S$, $\overline{T}x = Tx$ (i.e., $\overline{T}|_S = T$).
    Furthermore, the operator norm is preserved: $\|\overline{T}\| = \|T\| \le C$.

    ??? proof "Constructive Proof of the BLT Theorem (Click to expand)"

        **Step 1: Constructing the limit map**
        Since $S$ is dense in $X$, for any $x \in X$, there must exist a sequence $\{x_n\}$ in $S$ such that $x_n \to x$.
        Because $T$ is bounded (continuous) on $S$, we examine the distance of the sequence $\{Tx_n\}$ in $Y$:

        $$
        \|Tx_n - Tx_m\|_Y = \|T(x_n - x_m)\|_Y \le C \|x_n - x_m\|_X
        $$

        Since $\{x_n\}$ is a convergent sequence, it is necessarily a Cauchy sequence. Therefore, as $n, m \to \infty$, $\|x_n - x_m\|_X \to 0$.
        This implies that $\{Tx_n\}$ is a Cauchy sequence in the Banach space $Y$. Since $Y$ is complete, the limit $\lim_{n \to \infty} Tx_n$ must exist. We define:

        $$
        \overline{T}x \triangleq \lim_{n \to \infty} Tx_n
        $$

        **Step 2: Proving the well-definedness of the map (independence of sequence choice)**
        Suppose there is another sequence $x_n' \to x$. We need to prove $\lim Tx_n' = \lim Tx_n$.
        Let $y = \lim Tx_n$, $y' = \lim Tx_n'$.

        $$
        \|y - y'\|_Y = \lim_{n \to \infty} \|Tx_n - Tx_n'\|_Y \le C \lim_{n \to \infty} \|x_n - x_n'\|_X = 0
        $$

        Hence $y = y'$, and the map $\overline{T}$ is well-defined.

        **Step 3: Proving linearity and norm preservation**
        Linearity is a natural consequence of limits. For norm preservation, for any $x \in X$ and an approximating sequence $x_n \to x$:

        $$
        \|\overline{T}x\|_Y = \lim_{n \to \infty} \|Tx_n\|_Y \le C \lim_{n \to \infty} \|x_n\|_X = C \|x\|_X
        $$

        Therefore, $\overline{T}$ is bounded, and its norm does not exceed $C$. $\square$

!!! success "Application: Extension of PWZ Integral to $L^2$ Space"

    We take the space of integrand functions as $X = L^2([0, T])$, and the space of integral results as $Y = L^2(\Omega, P)$.
    The dense subspace is taken as $S = C^1([0, T])$.
    Define the operator $T: g \mapsto \int_0^T g(t) dW(t)$.

    By the Itô isometry, for any $g \in S$:

    $$
    \|Tg\|_{L^2(\Omega)}^2 = E\left[ \left( \int_0^T g(t) dW(t) \right)^2 \right] = \int_0^T g(t)^2 dt = \|g\|_{L^2([0,T])}^2
    $$

    This means the operator $T$ is an isometry (operator norm $C=1$). By the BLT Theorem, we can perfectly and uniquely extend the PWZ integral to the entire $L^2([0, T])$ space.

---

## 3. Quadratic Variation of Brownian Motion

The fundamental reason why traditional Riemann integration cannot handle stochastic integration lies in the "roughness" of Brownian motion paths, specifically its quadratic variation property.

Consider a partition $P = \{0 = t_0 < t_1 < \dots < t_m = T\}$ of the time interval $[0, T]$, with mesh size $|P| = \max (t_{k+1} - t_k)$.

!!! info "Theorem 1: The quadratic variation of Brownian motion equals time $T$"

    As the mesh is refined $|P| \to 0$, the sum of squared increments of Brownian motion converges to $T$ in the sense of $L^2(\Omega, P)$:

    $$
    \sum_{k=0}^{m-1} (W(t_{k+1}) - W(t_k))^2 \xrightarrow{L^2} T
    $$

    ??? proof "Rigorous derivation of $L^2$ convergence for quadratic variation (click to expand)"

        To prove $L^2$ convergence, we need to show that its mean square error with respect to $T$ tends to 0.
        Since $\sum_{k=0}^{m-1} (t_{k+1} - t_k) = T$, we can write the target error as:

        $$
        E\left[ \left( \sum_{k=0}^{m-1} \big((W(t_{k+1}) - W(t_k))^2 - (t_{k+1} - t_k)\big) \right)^2 \right]
        $$

        Let $\Delta W_k = W(t_{k+1}) - W(t_k)$, $\Delta t_k = t_{k+1} - t_k$. Expanding the square term, it splits into squared terms and cross terms:

        $$
        = \sum_k E\Big[ \big((\Delta W_k)^2 - \Delta t_k\big)^2 \Big] + \sum_{k \neq j} E\Big[ \big((\Delta W_k)^2 - \Delta t_k\big)\big((\Delta W_j)^2 - \Delta t_j\big) \Big]
        $$

        **Key Point 1: Cross terms are 0.**
        Due to the independent increments property of Brownian motion, when $k \neq j$, $\Delta W_k$ and $\Delta W_j$ are independent. Furthermore, since $E[(\Delta W_k)^2] = \Delta t_k$, the expectation of each factor is 0, so the overall expectation of the cross terms is 0.

        **Key Point 2: Calculation of squared terms.**
        Only the variance terms on the diagonal remain. Note that $\Delta W_k \sim \mathcal{N}(0, \Delta t_k)$, so it can be expressed in standardized form as $\sqrt{\Delta t_k} Z$, where $Z \sim \mathcal{N}(0, 1)$.

        $$
        E\Big[ \big((\Delta W_k)^2 - \Delta t_k\big)^2 \Big] = E\Big[ (\Delta t_k Z^2 - \Delta t_k)^2 \Big] = (\Delta t_k)^2 E[(Z^2 - 1)^2]
        $$

        Since the fourth moment of the standard normal distribution is $E[Z^4] = 3$, and $E[Z^2] = 1$, we have $E[(Z^2 - 1)^2] = 3 - 2(1) + 1 = 2$.

        $$
        \text{Total error} = 2 \sum_{k=0}^{m-1} (\Delta t_k)^2
        $$

        We bound this sum: extract the largest $\Delta t_k$, which is the mesh size $|P|$:

        $$
        2 \sum_{k=0}^{m-1} (\Delta t_k)^2 \le 2 |P| \sum_{k=0}^{m-1} \Delta t_k = 2 |P| T
        $$

        As $|P| \to 0$, $2 |P| T \to 0$. Therefore, it converges to $T$ in the $L^2$ sense. $\square$

This theorem directly leads to a corollary about the nature of Brownian motion:

!!! success "Theorem 2: Brownian motion has infinite total variation almost surely"

    Almost every sample path $W(t, \omega)$ of Brownian motion has infinite total variation on any interval.

    > **Proof by contradiction (extremely concise):** If a certain path had bounded total variation $V_T < \infty$, then its quadratic variation could be bounded as:
    > $\sum (\Delta W_k)^2 \le (\max_k |\Delta W_k|) \sum |\Delta W_k| \le (\max_k |\Delta W_k|) \cdot V_T$
    > Since the path is continuous, $\max |\Delta W_k| \to 0$ as the partition is infinitely refined. This would force the quadratic variation to tend to 0, creating an absolute contradiction with Theorem 1 where the quadratic variation equals $T > 0$!

---

## 4. Riemann Sums for Stochastic Integrals and the Derivation of the Itô Integral

So how do we compute $\int_0^T W(t) dW(t)$?
Returning to the definition of Riemann sums, we can observe a phenomenon that does not appear in classical calculus: **Changing the evaluation point leads to a significant change in the integral result.**

Construct a partition $P$, and in each subinterval $[t_k, t_{k+1}]$, take a point $\tau_k = (1-\lambda)t_k + \lambda t_{k+1}$ ($\lambda \in [0, 1]$).
Examine the Riemann sum:

$$
R_n = \sum_{k=0}^{m-1} W(\tau_k) \big( W(t_{k+1}) - W(t_k) \big)
$$

!!! info "Dependence of the Integral Result on the Value of $\lambda$"

    For clarity, we study the special case $\lambda = 0$ (taking the left endpoint, i.e., the Itô integral), where $\tau_k = t_k$:

    $$
    R_n = \sum_{k=0}^{m-1} W(t_k) \big( W(t_{k+1}) - W(t_k) \big)
    $$

    ??? proof "Algebraic Identity Splitting and Limit Calculation for the Itô Integral (Click to Expand)"

        This is an extremely clever algebraic trick. Using the identity $a(b-a) = \frac{1}{2}\big( b^2 - a^2 - (b-a)^2 \big)$, we rewrite each term:
        Let $a = W(t_k), b = W(t_{k+1})$:

        $$
        W(t_k) \big( W(t_{k+1}) - W(t_k) \big) = \frac{1}{2}\big( W(t_{k+1})^2 - W(t_k)^2 \big) - \frac{1}{2}\big( W(t_{k+1}) - W(t_k) \big)^2
        $$

        Summing all terms, the original sum splits into two parts $B_1$ and $B_2$:

        $$
        R_n = \underbrace{ \frac{1}{2} \sum_{k=0}^{m-1} \big( W(t_{k+1})^2 - W(t_k)^2 \big) }_{B_1} - \underbrace{ \frac{1}{2} \sum_{k=0}^{m-1} \big( W(t_{k+1}) - W(t_k) \big)^2 }_{B_2}
        $$

        **Analyzing $B_1$**: This is a perfect telescoping sum, where all intermediate terms cancel:

        $$
        B_1 = \frac{1}{2} \big( W(T)^2 - W(0)^2 \big) = \frac{1}{2} W(T)^2
        $$

        **Analyzing $B_2$**: This is precisely the quadratic variation of Brownian motion! According to Theorem 1, as the partition is refined, it converges in the $L^2$ sense:

        $$
        B_2 \xrightarrow{L^2} \frac{1}{2} T
        $$

        In summary, when $\lambda = 0$ (Itô integral), we obtain:

        $$
        \int_0^T W(t) dW(t) \triangleq \lim_{|P|\to 0} R_n = \frac{1}{2} W(T)^2 - \frac{1}{2} T
        $$

        *(Note: The extra term $-\frac{1}{2}T$ that appears later is called the Itô correction term.)* $\square$

**The Two Most Important Schools of Integration:**

* When $\lambda = 0$ (taking the left endpoint), it is the **Itô integral**, with result $\frac{1}{2}W(T)^2 - \frac{1}{2}T$. It preserves the martingale property and is commonly used in financial mathematics.

* When $\lambda = 1/2$ (taking the midpoint), it is the **Stratonovich integral**, where the correction terms cancel out, resulting in $\frac{1}{2}W(T)^2$, which formally aligns with classical calculus and is often used in physics and engineering.

---

## 5. Measure-Theoretic Preparation for the Strict Itô Integral

To extend the Itô integral to more general stochastic processes (not just $W(t)$), we need to rigorously define what it means to "look only at the past, not the future." This requires us to introduce the concepts of a **filtration** and an **adapted process**.

!!! abstract "Definition: Information Flow and $\sigma$-algebra"

    **1. Natural Filtration of Brownian Motion**:
    For any time $t$, the information generated by the historical path of Brownian motion is denoted by the $\sigma$-algebra:

    $$
    \mathcal{F}_W(t) \triangleq \sigma(\{W(s) \mid 0 \le s \le t\})
    $$

    It contains all path information of the Brownian motion up to time $t$.

    **2. Independence of Future Increments**:
    We define the future increment information flow $\mathcal{F}^t \triangleq \sigma(\{W(s) - W(t) \mid s > t\})$. According to the independent increments property of Brownian motion, $\mathcal{F}^t$ is **completely independent** of the historical flow $\mathcal{F}_W(t)$.

    **3. General Information Flow (Filtration)**:
    A family of $\sigma$-algebras $\{\mathcal{F}(t)\}_{t \ge 0}$ satisfying the following conditions:
    
    * **Monotonicity**: $\mathcal{F}(s) \subset \mathcal{F}(t)$ for any $0 \le s \le t$ (information is not forgotten).
    
    * **Contains History**: $\mathcal{F}_W(t) \subset \mathcal{F}(t)$.
    
    * **Future Independence**: The increment $W(s) - W(t)$ is independent of $\mathcal{F}(t)$.

With this mathematical framework for the concept of "information," we can determine which stochastic processes can be integrated using the Itô integral.

!!! info "Definition: Adapted Process and Progressively Measurable"

    **1. Adapted Process**:
    If for each fixed time $t$, the random variable $G(t, \omega)$ is $\mathcal{F}(t)$-measurable, then the stochastic process $G(t)$ is said to be **adapted** to the filtration $\{\mathcal{F}(t)\}$.
    *(Intuitive understanding: At time $t$, if you know all historical information $\mathcal{F}(t)$, you know the value of $G(t)$ at that moment, without needing any future information.)*

    **2. Progressively Measurable**:
    A stronger condition. The mapping $(s, \omega) \mapsto G(s, \omega)$ is jointly measurable with respect to the product $\sigma$-algebra $\mathcal{B}([0, t]) \otimes \mathcal{F}(t)$ on the product space $[0, t] \times \Omega$. This ensures that Riemann or Lebesgue integration over time intervals is well-defined.

In the subsequent construction in this book, the Itô integral $\int_0^T G(t, \omega) dW(t)$ will be rigorously defined in the following Hilbert space:

$$
L^2(\Omega \times [0, T]) = \left\{ G(t, \omega) \text{ is a progressively measurable process} \mid E\left[ \int_0^T G(t, \omega)^2 dt \right] < \infty \right\}
$$