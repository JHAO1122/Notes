# Chapter 1: Foundations of Conditional Expectation and Martingales

Before diving into stochastic differential equations, we need to move past how basic probability theory defines conditional expectation. Instead, we will rigorously define it using measure theory. This is the foundation for our study of martingale theory and stochastic analysis.

## 1. Measure-Theoretic Definition of Conditional Expectation

In classical probability theory, the conditional probability of event $A$ given event $B$ is defined as:
$P(A|B) \triangleq \frac{P(AB)}{P(B)}.$

The corresponding conditional expectation is:
$E[X|B] \triangleq \frac{1}{P(B)}\int_B X dP.$

But what if we are given a continuous random variable $Y$? The probability of the event $\{Y=y\}$ is 0, which makes the classical definition fail.

To fix this, we need to redefine it at the $\sigma$-algebra level using the Radon-Nikodym theorem.

!!! abstract "Definition: Conditional Expectation"

    Let $X$ be an integrable random variable (meaning $E[|X|] < \infty$) on a probability space $(\Omega, \mathcal{F}, P)$. Let $\mathcal{G}$ be a sub-$\sigma$-algebra of $\mathcal{F}$ (often written as $\mathcal{G} = \sigma(Y)$).
    We call the random variable $Z$ the **conditional expectation** of $X$ given $\mathcal{G}$, written as $E[X|\mathcal{G}] = Z$, if it meets these two conditions:

    **1. Measurability**: $Z$ is a $\mathcal{G}$-measurable random variable.

    **2. Partial Averaging (Integral Matching)**: For any $A \in \mathcal{G}$, we have:

    $$
    \int_A X dP = \int_A Z dP
    $$

    The Radon-Nikodym theorem guarantees that such a $Z$ exists and is unique almost surely (a.s.).

## 2. Core Properties

Let $\mathcal{G}$ and $\mathcal{H}$ be sub-$\sigma$-algebras.

!!! info "Basic Properties"

    **1. Linearity**:

    $$
    E[aX + bY | \mathcal{G}] = aE[X|\mathcal{G}] + bE[Y|\mathcal{G}] \quad a.s.
    $$

    **2. Expectation Preservation**:

    $$
    E[E[X|\mathcal{G}]] = E[X]
    $$

    > (You can prove this by setting $A = \Omega$ in the definition).

    **3. Knowns act as Constants (Taking out what is known)**:

    If $X$ is $\mathcal{G}$-measurable, then:

    $$
    E[XY|\mathcal{G}] = X E[Y|\mathcal{G}] \quad a.s.
    $$

    **4. Independence means Irrelevance**:

    If $X$ is independent of $\mathcal{G}$, then:

    $$
    E[X|\mathcal{G}] = E[X] \quad a.s.
    $$

    **5. Tower Property**:

    If $\mathcal{H} \subset \mathcal{G} \subset \mathcal{F}$, then:

    $$
    E\big[ E[X|\mathcal{G}] \big| \mathcal{H} \big] = E[X|\mathcal{H}] \quad a.s.
    $$

    > *(Intuition: The smaller information set dictates the result.)*

!!! tip "Geometric Interpretation: Orthogonal Projection in $L^2$ Space"

    Conditional expectation has a beautiful geometric interpretation.

    Consider the square-integrable space $L^2(\Omega, \mathcal{F}, P)$. This is a Hilbert space with the inner product defined as:

    $$
    (X,Y) = E[XY]
    $$

    Because $\mathcal{G} \subset \mathcal{F}$, the subspace $L^2(\Omega, \mathcal{G}, P)$ is also a closed subspace.

    From this perspective, **the conditional expectation $E[X|\mathcal{G}]$ is simply the orthogonal projection of $X$ onto the subspace $L^2(\Omega, \mathcal{G}, P)$**.

    It minimizes the mean squared error $E[(X - Z)^2]$ among all $\mathcal{G}$-measurable random variables $Z$:

    $$
    \min_{Z \in \mathcal{G}} E[(X - Z)^2] = E\big[(X - E[X|\mathcal{G}])^2\big]
    $$

## 3. Jensen's Inequality

In advanced probability, Jensen's inequality is a sharp tool for handling convex functions and limit theorems.

!!! success "Theorem: Conditional Jensen's Inequality"

    Let $\Phi: \mathbb{R} \rightarrow \mathbb{R}$ be a convex function, and $E[|\Phi(X)|] < \infty$. Then:

    $$
    \Phi(E[X|\mathcal{G}]) \le E[\Phi(X)|\mathcal{G}] \quad a.s.
    $$

    **Proof Idea**:
    We use the convex function property $\Phi(\sum a_i b_i) \le \sum b_i \Phi(a_i)$ (where $\sum b_i = 1$). Since conditional expectation acts like a weighted average, we can build simple functions using indicator functions to approximate it. Then we take the limit to generalize it. This is a standard proof method in measure theory: going from indicator functions to simple functions, and then to measurable functions.

## 4. Martingales and Doob's Inequalities

To study random phenomena that evolve over time (like Brownian motion), we use filtrations and martingales.

!!! abstract "Definition: Filtration and Martingale"

    **1. Filtration**: A family of increasing $\sigma$-algebras $\{\mathcal{F}_t\}_{t \ge 0}$. If $s \le t$, then $\mathcal{F}_s \subset \mathcal{F}_t$. It represents the history of information up to time $t$.

    **2. Martingale**: Let a process $\{X_t\}$ be adapted to the filtration $\mathcal{F}_t$, and $E[|X_t|] < \infty$. If for any $s \le t$:

    $$
    E[X_t | \mathcal{F}_s] = X_s \quad a.s.
    $$

    Then we call $\{X_t\}$ a **martingale**. (This represents a fair game, where the expected future value is just the current value).

    **3. Sub/Super-martingales**: If $E[X_t | \mathcal{F}_s] \ge X_s$, it is a **Submartingale**. If $E[X_t | \mathcal{F}_s] \le X_s$, it is a **Supermartingale**.

In stochastic analysis, we often need to bound the maximum value of a process along its path. Doob's inequality is the tool for this.

!!! info "Theorem: Doob's Maximal Inequality"

    **(1) Sublinear Bound**: Let $\{X_n\}_{n=1}^N$ be a non-negative submartingale. For any $\lambda > 0$:

    $$
    P\left(\max_{1 \le k \le n} X_k \ge \lambda\right) \le \frac{1}{\lambda} E[X_n I_{\{\max X_k \ge \lambda\}}] \le \frac{1}{\lambda} E[X_n]
    $$

    > *(It looks like Chebyshev's inequality, but it controls the extreme value of the entire path.)*

    **(2) $L^p$ Maximal Inequality**: Let $p > 1$ and $\frac{1}{p} + \frac{1}{q} = 1$. If $\{X_n\}$ is a non-negative submartingale and $X_n \in L^p$:

    $$
    E\left[\max_{1 \le k \le n} X_k^p\right] \le \left(\frac{p}{p-1}\right)^p E[X_n^p]
    $$

    **Key Steps for the $L^p$ Bound Derivation**:

    The trick is to use Fubini's theorem to swap the order of integration, then apply Hölder's inequality.
    Let $X^* = \max_{1 \le k \le n} X_k$. Using the identity:

    $$
    E[(X^*)^p] = \int_0^\infty p \lambda^{p-1} P(X^* \ge \lambda) d\lambda
    $$

    Substitute the sublinear bound, then use Fubini's theorem to swap the integral over $\lambda$ into the expectation:

    $$
    \le \int_{\Omega} X_n \left( \int_0^{X^*} p \lambda^{p-2} d\lambda \right) dP = \frac{p}{p-1} E[X_n (X^*)^{p-1}]
    $$

    Finally, apply Hölder's inequality $\int f g \le (\int f^p)^{1/p} (\int g^q)^{1/q}$ to split $X_n$ and $(X^*)^{p-1}$. Done.