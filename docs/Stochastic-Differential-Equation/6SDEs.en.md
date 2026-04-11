# Chapter 6: Existence and Uniqueness Theorem

In the theory of Ordinary Differential Equations (ODEs), we use the Picard iteration method combined with the Lipschitz condition to prove the existence and uniqueness of the solution to the initial value problem $y' = f(x,y), y(x_0)=y_0$.

For Stochastic Differential Equations (SDEs), we face a complex situation involving integration with respect to Brownian motion:

$$
dX_t = b(X_t, t)dt + \sigma(X_t, t)dW_t, \quad X_0 = x_0
$$

To guarantee that this equation has a unique strong solution, we need to impose similar restrictions on the drift coefficient $b$ and the diffusion coefficient $\sigma$. This section will detail this core theorem and its rigorous proof.

---

## 1. Existence and Uniqueness Theorem

!!! abstract "Theorem: Existence and Uniqueness of SDE Strong Solutions"

    Let $T > 0$. Suppose the functions $b(x,t): \mathbb{R}^n \times [0,T] \to \mathbb{R}^n$ and $\sigma(x,t): \mathbb{R}^n \times [0,T] \to \mathbb{R}^{n \times m}$ satisfy the following two fundamental conditions:

    **1. Global Lipschitz Condition**:
    There exists a constant $L > 0$ such that for all $x, y \in \mathbb{R}^n$ and $t \in [0,T]$:

    $$
    |b(x,t) - b(y,t)| + |\sigma(x,t) - \sigma(y,t)| \le L|x - y|
    $$

    **2. Linear Growth Condition**:
    There exists a constant $L > 0$ such that for all $x \in \mathbb{R}^n$ and $t \in [0,T]$:

    $$
    |b(x,t)|^2 + |\sigma(x,t)|^2 \le L^2(1 + |x|^2)
    $$

    **Initial Condition**:
    Let $\xi$ be a random variable independent of the Brownian motion $W_t$, satisfying a finite second moment $E|\xi|^2 < \infty$.

    **Conclusion**:
    Then, the initial value problem $X_t = \xi + \int_0^t b(X_s, s)ds + \int_0^t \sigma(X_s, s)dW_s$ has a **unique, square-integrable, continuous strong solution** $X_t$ on the interval $[0,T]$, which satisfies:

    $$
    E\left[ \sup_{0 \le t \le T} |X_t|^2 \right] < \infty
    $$

---

## 2. Proof of the Theorem

The proof is divided into two main parts: **Uniqueness** (using Gronwall's inequality) and **Existence** (using Picard iteration and the Borel-Cantelli lemma).

### 2.1 Proof of Uniqueness

??? proof "Proof of Uniqueness: Application of Gronwall's Inequality (Click to expand)"

    Assume there exist two solutions $X_t$ and $\tilde{X}_t$ satisfying the same initial condition $X_0 = \tilde{X}_0 = \xi$.

    Let the error process be $Y_t = X_t - \tilde{X}_t$, then $Y_0 = 0$. According to the integral form of the equation, we have:

    $$
    Y_t = \int_0^t [b(X_s, s) - b(\tilde{X}_s, s)]ds + \int_0^t [\sigma(X_s, s) - \sigma(\tilde{X}_s, s)]dW_s
    $$

    Using the basic inequality $(a+b)^2 \le 2a^2 + 2b^2$, taking the squared absolute value on both sides and taking the expectation:

    $$
    E|Y_t|^2 \le 2E \left| \int_0^t [b(X_s, s) - b(\tilde{X}_s, s)]ds \right|^2 + 2E \left| \int_0^t [\sigma(X_s, s) - \sigma(\tilde{X}_s, s)]dW_s \right|^2
    $$

    For the first term (Riemann integral), using the Cauchy-Schwarz inequality:

    $$
    E \left| \int_0^t [b(X_s, s) - b(\tilde{X}_s, s)]ds \right|^2 \le t E \int_0^t |b(X_s, s) - b(\tilde{X}_s, s)|^2 ds
    $$

    For the second term (Itô integral), using the Itô Isometry:

    $$
    E \left| \int_0^t [\sigma(X_s, s) - \sigma(\tilde{X}_s, s)]dW_s \right|^2 = E \int_0^t |\sigma(X_s, s) - \sigma(\tilde{X}_s, s)|^2 ds
    $$

    Substituting the Lipschitz conditions $|b(X) - b(\tilde{X})|^2 \le L^2 |X - \tilde{X}|^2 = L^2 |Y_s|^2$ and $|\sigma(X) - \sigma(\tilde{X})|^2 \le L^2 |Y_s|^2$ into the above equation:

    $$
    E|Y_t|^2 \le 2t L^2 \int_0^t E|Y_s|^2 ds + 2L^2 \int_0^t E|Y_s|^2 ds \le 2L^2(T+1) \int_0^t E|Y_s|^2 ds
    $$

    Letting the constant $C = 2L^2(T+1)$, we obtain:

    $$
    E|Y_t|^2 \le C \int_0^t E|Y_s|^2 ds
    $$

    Since $E|Y_t|^2$ is non-negative, **it follows directly from Gronwall's inequality** that for all $t \in [0,T]$:

    $$
    E|Y_t|^2 = 0
    $$

    This implies that for any given $t$, $X_t = \tilde{X}_t$ almost surely. Combined with the continuity of the sample paths, we conclude that these two solutions are indistinguishable on the entire interval $[0,T]$. Uniqueness is proven. $\square$

---

### 2.2 Proof of Existence (Picard Iteration Method)

The core idea of existence is to construct an approximating sequence. The logical chain of derivation in this part is relatively long and requires a step-by-step approach.

??? proof "Proof of Existence: Picard Iteration Method (Click to expand)"

    **Step 1: Construct the Picard iteration sequence**

    Define the initial approximation as $X^{(0)}_t = \xi$.
    
    For $k \ge 0$, define recursively:

    $$
    X^{(k+1)}_t = \xi + \int_0^t b(X^{(k)}_s, s)ds + \int_0^t \sigma(X^{(k)}_s, s)dW_s
    $$

    **Step 2: Prove the mean-square boundedness of the iteration sequence**

    We need to ensure that the second moment $E|X^{(k)}_t|^2$ of each element in the sequence is finite. This can be proven via mathematical induction and the linear growth condition:

    $$
    E|X^{(k)}_t|^2 \le C(1 + E|\xi|^2) e^{Ct}
    $$
    
    *(The derivation here is similar to the Cauchy-Schwarz inequality and Itô isometry used in the uniqueness proof, utilizing the linear growth condition. In the notes, it was abbreviated as $E|X^k_t|^2 \le Ce^{ct}$; the repetitive steps are omitted here).*

    **Step 3: Estimate the decay rate of the difference between adjacent iteration terms**

    Define the mean-square error between adjacent terms as $d_k(t) = E|X^{(k+1)}_t - X^{(k)}_t|^2$.
    
    Using the exact same bounding techniques as in the uniqueness proof (Cauchy-Schwarz + Itô Isometry + Lipschitz condition), we can obtain the recurrence relation:

    $$
    d_k(t) \le 2L^2(T+1) \int_0^t d_{k-1}(s) ds = M \int_0^t d_{k-1}(s) ds
    $$

    where $M = 2L^2(T+1)$.

    For $k=0$:
    
    $$
    d_0(t) = E|X^{(1)}_t - X^{(0)}_t|^2 \le M t (1 + E|\xi|^2) = C_0 t
    $$

    Recursively integrating backward:

    $$
    d_k(t) \le M \int_0^t d_{k-1}(s) ds \le M^k C_0 \int_0^t \frac{s^{k-1}}{(k-1)!} ds = C_0 \frac{M^k t^k}{k!}
    $$

    This shows that the expectation of the error decays extremely fast (at a factorial rate).

    **Step 4: Use Doob's maximal inequality to bound the supremum of the path**

    To prove uniform convergence, we need to examine the maximum error of the path over the entire interval. We use Doob's Martingale Inequality and the Burkholder-Davis-Gundy (BDG) inequality to handle the maximum of the Itô integral term:

    $$
    E\left[ \sup_{0 \le s \le T} |X^{(k+1)}_s - X^{(k)}_s|^2 \right] \le C_1 \int_0^T E|X^{(k)}_s - X^{(k-1)}_s|^2 ds \le C_2 \frac{(MT)^k}{k!}
    $$

    According to Markov's Inequality:

    $$
    P\left( \sup_{0 \le t \le T} |X^{(k+1)}_t - X^{(k)}_t| > 2^{-k} \right) \le \frac{E[\sup_{0 \le t \le T} |X^{(k+1)}_t - X^{(k)}_t|^2]}{(2^{-k})^2} \le 2^{2k} C_2 \frac{(MT)^k}{k!}
    $$

    **Step 5: Borel-Cantelli Lemma and uniform convergence**

    Since the factorial grows much faster than the exponential, the series converges:

    $$
    \sum_{k=0}^{\infty} 2^{2k} C_2 \frac{(MT)^k}{k!} < \infty
    $$

    By the **Borel-Cantelli Lemma**, the probability that the event $\{\sup_{0 \le t \le T} |X^{(k+1)}_t - X^{(k)}_t| > 2^{-k}\}$ occurs infinitely often is 0.

    This implies that for almost all sample paths $\omega$, there exists a $K(\omega)$ such that for $k > K(\omega)$, $\sup_{0 \le t \le T} |X^{(k+1)}_t - X^{(k)}_t| \le 2^{-k}$.
    
    Since $\sum 2^{-k}$ converges, this guarantees that the series $\sum (X^{(k+1)}_t - X^{(k)}_t)$ converges absolutely and uniformly. Therefore, the limit process exists:

    $$
    X_t = X^{(0)}_t + \sum_{k=0}^{\infty} (X^{(k+1)}_t - X^{(k)}_t) = \lim_{k \to \infty} X^{(k)}_t
    $$

    And this convergence is uniform over $t \in [0,T]$, ensuring that the limit process $X_t$ is continuous.

    **Step 6: Consistency - Verifying the limit is a solution**

    Finally, we need to verify that this almost surely uniformly convergent limit $X_t$ indeed satisfies the original equation.
    
    Taking the limit $k \to \infty$ on both sides of the Picard iteration formula. Since $b$ and $\sigma$ satisfy Lipschitz continuity, combined with the Dominated Convergence Theorem and the isometry property of the Itô integral, the limit can pass through the integral sign:

    $$
    \int_0^t \sigma(X^{(k)}_s) dW_s \xrightarrow{L^2} \int_0^t \sigma(X_s) dW_s
    $$

    This verifies that $X_t = \xi + \int_0^t b(X_s, s)ds + \int_0^t \sigma(X_s, s)dW_s$. Existence is proven. $\square$

---

## 3. Higher-Order Moments and Initial Value Dependence of the Solution

In the final part of the notes, bounds on higher-order moments of the solution and its continuous dependence on the initial value are mentioned. These are very important corollaries in SDE theory.

!!! info "Properties: Higher-Order Moments and Dependence"

    **1. Bounds on Higher-Order Moments**
    
    If the initial condition $\xi$ satisfies $E|\xi|^{2p} < \infty$ ($p \ge 1$), then under the global Lipschitz and linear growth conditions, the strong solution $X_t$ satisfies:

    $$
    E\left[ \sup_{0 \le t \le T} |X_t|^{2p} \right] \le C(1 + E|\xi|^{2p}) e^{CT}
    $$
    
    *(The proof similarly utilizes Itô's formula, the BDG inequality, and Gronwall's inequality).*

    **2. Continuous Dependence on the Initial Value**
    
    Let $X^x_t$ and $X^y_t$ be strong solutions with initial values $X_0 = x$ and $X_0 = y$ respectively, then the difference between them is bounded by the difference in their initial values:

    $$
    E\left[ \sup_{0 \le t \le T} |X^x_t - X^y_t|^2 \right] \le C |x - y|^2 e^{CT}
    $$
    
    This indicates that the solution to an SDE is Lipschitz continuously dependent on its initial state, which plays a crucial role in studying the connection between Markov semigroups and Partial Differential Equations (PDEs).