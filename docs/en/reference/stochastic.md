# 🌊 Stochastic Processes

This module focuses on content related to stochastic processes.

## 1. Algebraic Properties of Conditional Expectation

!!! abstract "Core Computational Properties"
    The conditional expectation $E[X \mid \mathcal{F}]$ itself is a random variable measurable with respect to $\mathcal{F}$.
    
    1. **Pulling out knowns**: If $Y$ is $\mathcal{F}$-measurable, then $E[XY \mid \mathcal{F}] = Y E[X \mid \mathcal{F}]$.
    2. **Tower Property**: If $\mathcal{F}_1 \subset \mathcal{F}_2$, then $E[E[X \mid \mathcal{F}_2] \mid \mathcal{F}_1] = E[X \mid \mathcal{F}_1]$.
    3. **Jensen's Inequality**: For a convex function $\phi$, we have $\phi(E[X \mid \mathcal{F}]) \le E[\phi(X) \mid \mathcal{F}]$. **This inequality is often used to prove that a convex function of a martingale is a submartingale**.

---

## 2. Martingales and Stopping Times

!!! info "Definition of Martingale"
    An integrable sequence $X_n$ adapted to the filtration $\mathcal{F}_n$, satisfying:

    \[ 
    E[X_{n+1} \mid \mathcal{F}_n] = X_n 
    \]

    (If it is $\ge$, it is a Submartingale; if it is $\le$, it is a Supermartingale).

!!! info "Stopping Time"
    A random variable $T$ taking values in $\mathbb{N} \cup \{\infty\}$ is called a stopping time if, for any $n$, the event $\{T \le n\}$ is completely determined by $\mathcal{F}_n$ (i.e., whether it occurs depends only on the information available up to that time).
    
    * **Example**: The first hitting time $\inf\{n: S_n = 0\}$ is a stopping time.

??? success "Optional Stopping Theorem"
    Under certain conditions, the expectation of a martingale at a stopping time $T$ is equal to its initial expectation: $E[X_T] = E[X_0]$. It holds if **any one** of the following conditions is met:
    
    1. **Bounded time**: The stopping time $T$ is almost surely bounded (i.e., there exists a constant $N$ such that $P(T \le N) = 1$).
    2. **Uniformly bounded process**: For all $n \le T$, there exists a constant $M$ such that $|X_n| \le M$ almost surely.
    3. **Finite expected time and bounded increments**: $E[T] < \infty$, and the single-step increments of the martingale are bounded ($|X_m - X_n| \le C$).

---

## 3. Core Limit Theorems of Martingales

!!! abstract "Doob's Maximal Inequality"
    For a non-negative submartingale $X_n$ (or the absolute value of a martingale), the tail probability of its running maximum is bounded by its final expectation:

    \[ 
    \lambda P(\max_{1\le k\le n} X_k \ge \lambda) \le E[X_n I(\max_{1\le k\le n} X_k \ge \lambda)] \le E[X_n] 
    \]

!!! abstract "Doob's Martingale Convergence Theorem"
    If $X_n$ is a submartingale and the supremum of its positive part's expectation is bounded (i.e., $\sup_n E[X_n^+] < \infty$), then there exists an integrable random variable $X$ such that:

    \[ 
    X_n \xrightarrow{a.s.} X 
    \]

    * Note: This only guarantees almost sure convergence, not $L^1$ convergence. If $E[X_n] \to E[X]$ is required, the condition of **Uniform Integrability (UI)** must be added.

??? note "Doob-Meyer Decomposition Theorem"
    Any submartingale $X_n$ can be uniquely decomposed into the sum of a martingale $M_n$ and a non-decreasing predictable sequence $A_n$:

    \[ 
    X_n = M_n + A_n 
    \]
    
    Here, $A_n$ being predictable means that $A_n$ is completely known at time $n-1$ (i.e., it is $\mathcal{F}_{n-1}$-measurable).