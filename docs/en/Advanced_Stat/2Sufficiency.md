# Chapter 2: Sufficiency, Minimal Sufficiency, and Exponential Families

This chapter studies the problem of exact data reduction in statistical inference. We first define sufficient statistics via a common conditional distribution kernel, then prove the Fisher–Neyman factorization theorem, and illustrate how to use the factorization criterion through Bernoulli, Poisson, normal, uniform, and Gamma-type models. Next, we introduce exponential families and the geometric structure of their natural parameter spaces, discuss minimal sufficiency, and prove that the natural statistic of a full affine rank exponential family is a minimal sufficient statistic.

The core idea of this chapter is: **If all dependence of the sample on the parameter is transmitted through a statistic $T$, then after conditioning on $T$, the remaining randomness in the sample is independent of the parameter.**

---

## 1. Why Sufficiency Is an Exact Data Reduction Principle

Let

\[
\mathcal{P}=\{P_\theta:\theta\in\Theta\}
\]

be a statistical model on a measurable sample space $(\mathcal{X},\mathcal{B})$, and let $T:\mathcal{X}\to\mathcal{T}$ be a statistic. The basic question is: Does observing $T(X)$ retain all the information about $\theta$ contained in observing the full sample $X$?

!!! info "Definition 1.1 (Sufficiency via a Common Kernel)"

    Assume all spaces involved are standard Borel spaces, so that regular conditional distributions exist. If there exists a single Markov kernel $K(t,B)$ independent of $\theta$ such that for every $B\in\mathcal{B}$ and every $\theta$,

    \[
    P_\theta(X\in B\mid T)=K\{T(X),B\}
    \qquad P_\theta\text{-a.s.},
    \]

    then the statistic $T$ is called **sufficient** for the model $\mathcal{P}$.

    Thus, once $T$ is known, the remaining randomness in $X$ no longer depends on the parameter. The phrase “single kernel” is important: under each $P_\theta$, conditional probabilities are defined only almost surely, and sufficiency requires that we can choose a version of the conditional distribution that is consistent across the entire model.

!!! example "Example 1 (Bernoulli Sample)"

    Let

    \[
    X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
    \qquad
    T=\sum_{i=1}^nX_i.
    \]

    If $x\in\{0,1\}^n$ and $\sum_i x_i=t$, then

    \[
    P_p(X=x\mid T=t)
    =\frac{p^t(1-p)^{n-t}}{\binom{n}{t}p^t(1-p)^{n-t}}
    =\binom{n}{t}^{-1}.
    \]

    The conditional distribution is uniform over all binary sequences containing exactly $t$ ones, and hence does not depend on $p$. The statistic $T$ discards the order in which successes occur, but discards no information about $p$.

In the above example, directly computing the conditional distribution is very intuitive, but it is inconvenient for continuous data and almost infeasible for samples with complex structure. The Fisher–Neyman factorization theorem will provide a practically useful criterion.

---

## 2. The Domination Lemma: Why One Reference Measure Is Enough

Assume every $P_\theta$ is dominated by the same $\sigma$-finite measure $\mu$. Proving the Fisher–Neyman factorization theorem requires the following result, commonly called the **Halmos–Savage lemma**.

!!! info "Lemma (Countable Domination / Halmos–Savage)"

    For a family of probability measures $\mathcal{P}$, the following three statements are equivalent:

    1. There exists a $\sigma$-finite measure $\mu$ such that $P_\theta\ll\mu$ for every $\theta$;
    2. There exist $P_{\theta_1},P_{\theta_2},\ldots\in\mathcal{P}$ such that for every $A\in\mathcal{B}$,

    \[
    \bigl[P_\theta(A)=0\text{ for all }\theta\bigr]
    \quad\Longleftrightarrow\quad
    \bigl[P_{\theta_j}(A)=0\text{ for all }j\bigr];
    \]

    3. There exists a probability measure of the form

    \[
    \lambda=\sum_{j=1}^{\infty}w_jP_{\theta_j},
    \qquad
    w_j>0,
    \qquad
    \sum_{j=1}^{\infty}w_j=1
    \]

    such that $P_\theta\ll\lambda$ for every $\theta$.

The nontrivial direction is from statement 1 to statement 3; its full proof is given in Section 10. The key point of this result is that even if the parameter space $\Theta$ is uncountable, a mixture of countably many distributions in the model can still capture all null sets involved in the entire statistical experiment.

**Two simple implication directions:** If statement 2 holds, one may take $w_j=2^{-j}$ in the mixture. Then $\lambda(A)=0$ if and only if $P_{\theta_j}(A)=0$ for all selected $j$, and by statement 2 this implies $P_\theta(A)=0$ for all $\theta$, so statement 2 implies statement 3. Statement 3 implies statement 1 directly by taking $\mu=\lambda$.

---

## 3. The Fisher–Neyman Factorization Theorem

Let

\[
p_\theta=\frac{dP_\theta}{d\mu}.
\]

We allow the support of the distributions to depend on the parameter; zeros in the density $p_\theta$ are also part of the factorization.

!!! success "Theorem 1 (Fisher–Neyman Factorization Theorem)"

    Under the above assumptions, $T$ is sufficient for $\mathcal{P}$ if and only if there exist nonnegative measurable functions $g_\theta$ and $h$ such that

    \[
    p_\theta(x)=g_\theta\{T(x)\}h(x)
    \qquad \mu\text{-a.e.},
    \]

    where $h$ does not depend on $\theta$.

??? proof "Proof of Theorem 1 (click to expand)"

    Let

    \[
    \lambda=\sum_jw_jP_{\theta_j}
    \]

    be the countable mixture given by the domination lemma. Then for every $\theta$,

    \[
    P_\theta\ll\lambda\ll\mu.
    \]

    **Sufficiency implies factorization.** Let $K(T,B)$ be the common kernel in Definition 1.1. Since under each selected $P_{\theta_j}$ the same kernel is a conditional distribution, after taking the mixture it is also a conditional distribution under $\lambda$:

    \[
    K(T,B)=E_\lambda(\mathbf{1}_B\mid T)
    \qquad \lambda\text{-a.s.}
    \]

    Write

    \[
    f_\theta=\frac{dP_\theta}{d\lambda}.
    \]

    For any $B\in\mathcal{B}$,

    \[
    \begin{aligned}
    P_\theta(B)
    &=\int K(T,B)\,dP_\theta\\
    &=\int E_\lambda(\mathbf{1}_B\mid T)f_\theta\,d\lambda\\
    &=\int E_\lambda(\mathbf{1}_B\mid T)E_\lambda(f_\theta\mid T)\,d\lambda\\
    &=\int_B E_\lambda(f_\theta\mid T)\,d\lambda.
    \end{aligned}
    \]

    By uniqueness of the Radon–Nikodym derivative,

    \[
    f_\theta=E_\lambda(f_\theta\mid T)
    \qquad \lambda\text{-a.s.}
    \]

    Hence $f_\theta$ is $\sigma(T)$-measurable. By the Doob–Dynkin lemma, there exists a measurable function $g_\theta$ such that

    \[
    f_\theta(x)=g_\theta\{T(x)\}.
    \]

    If we set $h=d\lambda/d\mu$, then by the chain rule for Radon–Nikodym derivatives,

    \[
    p_\theta(x)
    =\frac{dP_\theta}{d\lambda}(x)
    \frac{d\lambda}{d\mu}(x)
    =g_\theta\{T(x)\}h(x),
    \]

    which gives the required factorization.

    **Factorization implies sufficiency.** Suppose

    \[
    p_\theta(x)=g_\theta\{T(x)\}h(x).
    \]

    For the same mixture $\lambda$, we have

    \[
    p_\lambda(x)
    =\sum_jw_jp_{\theta_j}(x)
    =h(x)q\{T(x)\},
    \qquad
    q(t)=\sum_jw_jg_{\theta_j}(t).
    \]

    Since every $P_\theta\ll\lambda$, on $\{q>0\}$,

    \[
    \frac{dP_\theta}{d\lambda}(x)
    =r_\theta\{T(x)\},
    \qquad
    r_\theta(t)=\frac{g_\theta(t)}{q(t)};
    \]

    define this ratio arbitrarily elsewhere. Let $K_\lambda(t,B)$ be a regular conditional distribution of $X$ given $T=t$ under $\lambda$. If $C\in\sigma(T)$, then $\mathbf{1}_Cr_\theta(T)$ is $\sigma(T)$-measurable, and

    \[
    \begin{aligned}
    \int_CK_\lambda(T,B)\,dP_\theta
    &=\int_CK_\lambda(T,B)r_\theta(T)\,d\lambda\\
    &=\int\mathbf{1}_C\mathbf{1}_Br_\theta(T)\,d\lambda\\
    &=P_\theta(B\cap C).
    \end{aligned}
    \]

    Therefore, $K_\lambda(T,B)$ is also a conditional probability of event $B$ under every $P_\theta$. It is a common kernel independent of the parameter, so $T$ is sufficient. $\square$

!!! note "Interpretation of the Proof"

    In the forward proof, sufficiency forces every likelihood ratio $dP_\theta/d\lambda$ to be a function of $T$. In the reverse proof, these likelihood ratios only reweight the distribution of $T$, without changing the conditional distribution of $X$ given $T$. This is the rigorous meaning of “all parameter dependence is transmitted through $T$.”

---

## 4. Using the Factorization Theorem: Four Examples

### 4.1 Bernoulli and Poisson Samples

If

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
\]

then

\[
p_p(x)=p^{\sum_i x_i}(1-p)^{n-\sum_i x_i},
\]

so $\sum_iX_i$ is a sufficient statistic.

Similarly, if

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Poisson}(\lambda),
\]

then

\[
p_\lambda(x)
=e^{-n\lambda}\lambda^{\sum_i x_i}
\prod_{i=1}^n\frac{1}{x_i!},
\]

so $\sum_iX_i$ is sufficient for $\lambda$.

### 4.2 Mean and Variance of the Normal Distribution

If

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2),
\]

then the joint density can be written as

\[
p_{\mu,\sigma^2}(x)
=(2\pi\sigma^2)^{-n/2}
\exp\left\{-\frac{1}{2\sigma^2}
\left(\sum_i x_i^2-2\mu\sum_i x_i+n\mu^2\right)\right\}.
\]

Therefore

\[
S(X)=\left(\sum_iX_i,\sum_iX_i^2\right)
\]

is sufficient for $(\mu,\sigma^2)$. Also, because

\[
\sum_iX_i=n\overline{X},
\qquad
\sum_iX_i^2=(n-1)S^2+n\overline{X}^{,2},
\]

$S(X)$ is equivalent to $(\overline{X},S^2)$ when $n\geq2$.

### 4.3 Endpoint of the Uniform Distribution: Parameter-Dependent Support

If

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}(0,\theta),
\]

then

\[
p_\theta(x)
=\theta^{-n}\mathbf{1}\{0<x_{(1)}\}
\mathbf{1}\{x_{(n)}<\theta\}.
\]

Therefore $X_{(n)}$ is a sufficient statistic. The indicator function involving $\theta$ belongs to $g_\theta\{X_{(n)}\}$; the factorization theorem does not require the model to have common support. However, because the support of this model changes with $\theta$, subsequent regular likelihood theory must handle it separately.

### 4.4 A Gamma-Type Family

Assume

\[
f_{a,b}(x)=H(a,b)x^ae^{-bx^c}\mathbf{1}\{x>0\},
\qquad
a>-1,
\qquad
b>0,
\]

where $c$ is known. For an iid sample,

\[
\prod_{j=1}^nf_{a,b}(x_j)
=H(a,b)^n
\exp\left\{a\sum_j\log x_j-b\sum_jx_j^c\right\}
\prod_j\mathbf{1}\{x_j>0\}.
\]

Therefore

\[
\left(\sum_j\log X_j,\sum_jX_j^c\right)
\]

is sufficient for $(a,b)$.

---

## 5. Exponential Families

!!! info "Definition 5.1 (Exponential Family)"

    If a dominated statistical model can be written as

    \[
    p_\vartheta(x)
    =\exp\{\eta(\vartheta)^TS(x)-B(\vartheta)\}h(x),
    \]

    then it is called an **exponential family** with a $k$-dimensional statistic. Here $S=(S_1,\ldots,S_k)^T$ is the natural statistic, and $\eta(\vartheta)$ is the natural parameter. By Theorem 1, $S$ is a sufficient statistic.

The original parameter $\vartheta$ and the natural parameter $\eta$ need not be the same. This distinction is especially important for the normal family, because

\[
\eta_1=\frac{\mu}{\sigma^2},
\qquad
\eta_2=-\frac{1}{2\sigma^2}.
\]

!!! info "Definition 5.2 (Canonical Form and Natural Parameter Space)"

    The **canonical form** of an exponential family is

    \[
    p_\eta(x)=\exp\{\eta^TS(x)-A(\eta)\}h(x),
    \qquad
    \eta\in\mathcal{H},
    \]

    where the natural parameter space is

    \[
    \mathcal{H}
    =\left\{\eta\in\mathbb{R}^k:
    0<\int e^{\eta^TS(x)}h(x)\,d\mu(x)<\infty\right\},
    \]

    and the log-normalizing function is

    \[
    A(\eta)
    =\log\int e^{\eta^TS(x)}h(x)\,d\mu(x).
    \]

For iid observations, the canonical form becomes

\[
p_\eta(x_1,\ldots,x_n)
=\exp\left\{\eta^T\sum_{i=1}^nS(x_i)-nA(\eta)\right\}
\prod_{i=1}^nh(x_i).
\]

Therefore, the dimension of the sufficient statistic $\sum_iS(X_i)$ does not grow with the sample size $n$.

!!! success "Proposition (Geometry and Moments of the Log-Normalizing Function)"

    The natural parameter space $\mathcal{H}$ is convex, and the function $A$ is convex. At interior points where differentiation and integration may be interchanged,

    \[
    \nabla A(\eta)=E_\eta S(X),
    \qquad
    \nabla^2A(\eta)=\operatorname{Cov}_\eta\{S(X)\}.
    \]

    Hence, in a direction $a$ for which $a^TS(X)$ is not constant, $A$ is strictly convex.

??? proof "Proof of the Proposition (click to expand)"

    For any $\eta,\zeta\in\mathcal{H}$ and $0<t<1$, by Hölder's inequality,

    \[
    \int e^{\{(1-t)\eta+t\zeta\}^TS}h\,d\mu
    \leq
    \left(\int e^{\eta^TS}h\,d\mu\right)^{1-t}
    \left(\int e^{\zeta^TS}h\,d\mu\right)^t.
    \]

    Therefore $(1-t)\eta+t\zeta\in\mathcal{H}$. Taking logarithms on both sides gives

    \[
    A\{(1-t)\eta+t\zeta\}
    \leq(1-t)A(\eta)+tA(\zeta),
    \]

    so $\mathcal{H}$ and $A$ have the stated convexity properties.

    Differentiating the normalizing identity

    \[
    e^{A(\eta)}=\int e^{\eta^TS}h\,d\mu
    \]

    gives

    \[
    \partial_jA(\eta)
    =e^{-A(\eta)}\int S_j(x)e^{\eta^TS(x)}h(x)\,d\mu(x)
    =E_\eta S_j(X).
    \]

    Differentiating once more gives

    \[
    \partial_{j\ell}^2A
    =E_\eta(S_jS_\ell)-E_\eta S_jE_\eta S_\ell,
    \]

    proving the matrix identity. Finally,

    \[
    a^T\nabla^2A(\eta)a
    =\operatorname{Var}_\eta\{a^TS(X)\},
    \]

    which is positive exactly when $a^TS(X)$ is not constant. $\square$

| Model | Natural statistic $S(x)$ | Natural parameter $\eta$ | Natural parameter space |
| --- | --- | --- | --- |
| $\operatorname{Bernoulli}(p)$ | $x$ | $\log\{p/(1-p)\}$ | $\mathbb{R}$ |
| $\operatorname{Poisson}(\lambda)$ | $x$ | $\log\lambda$ | $\mathbb{R}$ |
| $\operatorname{Exponential}(\lambda)$ | $x$ | $-\lambda$ | $(-\infty,0)$ |
| $N(\mu,\sigma^2)$ | $(x,x^2)$ | $(\mu/\sigma^2,-1/(2\sigma^2))$ | $\mathbb{R}\times(-\infty,0)$ |
| $\operatorname{Gamma}(\alpha,\beta)$ | $(\log x,x)$ | $(\alpha-1,-\beta)$ | $(-1,\infty)\times(-\infty,0)$ |

---

## 6. Minimal Sufficiency

The full sample itself is always sufficient. Minimal sufficiency seeks the coarsest sufficient summary of the data, in which transformations that are one-to-one with each other are regarded as equivalent.

!!! info "Definition 6.1 (Minimal Sufficient Statistic)"

    A sufficient statistic $T$ is called **minimal sufficient** if for every sufficient statistic $U$, there exists a measurable function $H$ such that

    \[
    T=H(U)
    \qquad P_\theta\text{-a.s. for every }\theta.
    \]

Therefore, any two minimal sufficient statistics are functions of each other; after removing null sets common to the entire model, they encode the same partition of the sample space.

### 6.1 The Likelihood Proportionality Criterion

This subsection assumes that all densities have common support $\mathcal{S}$.

!!! success "Theorem 2 (Likelihood Ratio Criterion)"

    Assume $T$ is a sufficient statistic and that for any $x,y\in\mathcal{S}$,

    \[
    T(x)=T(y)
    \quad\Longleftrightarrow\quad
    \frac{p_\theta(x)}{p_\theta(y)}
    \text{ is independent of }\theta.
    \]

    Then $T$ is a minimal sufficient statistic.

??? proof "Proof of Theorem 2 (click to expand)"

    Let $U$ be an arbitrary sufficient statistic. By the factorization theorem, we can write

    \[
    p_\theta(x)=a_\theta\{U(x)\}b(x).
    \]

    If $U(x)=U(y)$, then

    \[
    \frac{p_\theta(x)}{p_\theta(y)}
    =\frac{b(x)}{b(y)},
    \]

    which is independent of $\theta$. By the reverse implication in the criterion, $T(x)=T(y)$. Therefore, $T$ is constant on every fiber of $U$. By the measurable factorization lemma, there exists a measurable function $H$ such that $T=H(U)$. Since $U$ is an arbitrary sufficient statistic, $T$ is minimal sufficient. $\square$

!!! note "Why Sufficiency Is Assumed Separately"

    In the standard Euclidean models used in this course, textbooks often give a slightly stronger version: under routine measurability conditions, the likelihood ratio condition alone can also prove sufficiency by selecting a likelihood representative from each fiber. Here, sufficiency in all applications has already been established by the factorization theorem. Writing the sufficiency assumption explicitly avoids hiding a measurable selection step inside the likelihood ratio argument.

!!! note "What the Criterion Says"

    Two samples belong to the same minimal sufficiency equivalence class if and only if their likelihood functions have the same shape as functions of the parameter; they differ at most by a multiplicative constant depending only on the data.

!!! example "Example 2 (Normal Family)"

    For the normal model with both mean and variance unknown,

    \[
    \log\frac{p_{\mu,\sigma^2}(x)}{p_{\mu,\sigma^2}(y)}
    =-\frac{\sum_i x_i^2-\sum_i y_i^2}{2\sigma^2}
    +\frac{\mu}{\sigma^2}\left(\sum_i x_i-\sum_i y_i\right).
    \]

    This expression is independent of $(\mu,\sigma^2)$ if and only if both differences in the expression are zero. Therefore,

    \[
    \left(\sum_iX_i,\sum_iX_i^2\right)
    \]

    is a minimal sufficient statistic.

!!! example "Example 3 (Endpoint of the Uniform Distribution)"

    Although this model does not have common support, the same likelihood contour idea remains valid. On the positive quadrant, two likelihood functions from Section 4.3 are proportional for all $\theta$ if and only if the maxima of the two samples are equal. If the maxima are different, one can choose $\theta$ strictly between them so that one likelihood is zero and the other is positive. Therefore, $X_{(n)}$ is a minimal sufficient statistic.

    This conclusion uses the more general cross-product version of the likelihood ratio criterion that allows likelihoods to take the value zero.

---

## 7. Why Full Affine Rank Implies Minimal Sufficiency

We prove the exponential family result by reducing the model to a finite subfamily. The following two lemmas isolate the key steps in the argument.

!!! info "Lemma 3 (Likelihood Ratio Statistic for a Finite Family)"

    Let $P_0,\ldots,P_k$ have common support, with densities $p_0,\ldots,p_k$, respectively. Then

    \[
    R(X)=\left(\frac{p_1(X)}{p_0(X)},\ldots,
    \frac{p_k(X)}{p_0(X)}\right)
    \]

    is a minimal sufficient statistic for this finite family of distributions.

??? proof "Proof of Lemma 3 (click to expand)"

    The factorization

    \[
    p_j(x)=R_j(x)p_0(x),
    \qquad
    j=1,\ldots,k,
    \qquad
    R_0\equiv1
    \]

    proves that $R$ is sufficient.

    If $U$ is an arbitrary sufficient statistic, write

    \[
    p_j(x)=a_j\{U(x)\}b(x).
    \]

    On the common support,

    \[
    \frac{p_j(x)}{p_0(x)}
    =\frac{a_j\{U(x)\}}{a_0\{U(x)\}},
    \]

    so every component of $R$ is a function of $U$. Therefore, for every sufficient statistic $U$, we have $R=H(U)$, and hence $R$ is minimal sufficient. $\square$

!!! info "Lemma 4 (Minimal Sufficiency on a Subfamily)"

    Assume all distributions in $\mathcal{P}$ have common support. If $T$ is sufficient for $\mathcal{P}$ and minimal sufficient for some subfamily $\mathcal{P}_0\subset\mathcal{P}$, then $T$ is minimal sufficient for $\mathcal{P}$.

??? proof "Proof of Lemma 4 (click to expand)"

    Let $U$ be sufficient for $\mathcal{P}$; then it is also sufficient for $\mathcal{P}_0$. The minimal sufficiency of $T$ on $\mathcal{P}_0$ gives

    \[
    T=H(U)
    \]

    almost everywhere on the common support. Since every distribution in $\mathcal{P}$ has this support, the same equality holds almost surely under every $P\in\mathcal{P}$. Therefore, $T$ is a function of every sufficient statistic for the full family. $\square$

!!! info "Definition 7.1 (Full Affine Rank and Curved Families)"

    For a canonical exponential family, let the parameter set be $H\subseteq\mathcal{H}$. We say that $H$ has **full affine rank** if the affine hull of $H$ is $\mathbb{R}^k$, equivalently, if $H$ contains $k+1$ affinely independent points.

    If $H$ contains a nonempty open set in $\mathbb{R}^k$, the exponential family is called **full-rank**. Full rank is sufficient for full affine rank, but not necessary. If a family has no redundancy in its displayed coordinates while the interior of the parameter set is empty, it is often called a **curved family**.

!!! success "Theorem 5 (Minimal Sufficiency in Exponential Families)"

    Consider the exponential family

    \[
    p_\eta(x)=\exp\{\eta^TS(x)-A(\eta)\}h(x),
    \qquad
    \eta\in H\subseteq\mathcal{H},
    \]

    and assume its common support is $\{h>0\}$. If $H$ contains $k+1$ affinely independent points, then $S$ is a minimal sufficient statistic. In particular, the conclusion holds for every full-rank exponential family.

??? proof "Proof of Theorem 5 (click to expand)"

    The factorization theorem first shows that $S$ is sufficient for the entire family. Choose affinely independent points

    \[
    \eta^{(0)},\eta^{(1)},\ldots,\eta^{(k)},
    \]

    and consider the corresponding finite subfamily. For $j=1,\ldots,k$,

    \[
    \log\frac{p_{\eta^{(j)}}(x)}{p_{\eta^{(0)}}(x)}
    =\{\eta^{(j)}-\eta^{(0)}\}^TS(x)
    -\{A(\eta^{(j)})-A(\eta^{(0)})\}.
    \]

    Let $D$ be a $k\times k$ matrix whose $j$th row is $\{\eta^{(j)}-\eta^{(0)}\}^T$. Affine independence implies that $D$ is nonsingular. If $r(x)$ is the vector of the above $k$ log-likelihood ratios and $a$ is the vector of differences of log-normalizing functions, then

    \[
    r(x)=DS(x)-a,
    \qquad
    S(x)=D^{-1}\{r(x)+a\}.
    \]

    Therefore, $S$ is in one-to-one correspondence with the likelihood ratio statistic in Lemma 3. Hence $S$ is minimal sufficient for the selected finite subfamily, and by Lemma 4 its minimal sufficiency extends to the entire family. $\square$

!!! success "Corollary 6 (IID Observations)"

    For an iid sample from the above exponential family, if $H$ has full affine rank, then

    \[
    \sum_{i=1}^nS(X_i)
    \]

    is a minimal sufficient statistic.

??? proof "Proof of Corollary 6 (click to expand)"

    The joint density of an iid sample is still a canonical exponential family, with natural statistic $\sum_iS(X_i)$ and the same natural parameter set $H$. A direct application of Theorem 5 gives the result. $\square$

!!! note "Why the Open Set Condition Is Stronger Than Needed"

    An open subset of $\mathbb{R}^k$ necessarily contains $k+1$ affinely independent points, so the full-rank condition is convenient to state. But the proof actually uses only that the parameter set has a full affine hull. Some curved parameter sets also have a full affine hull, and hence still yield the same minimal sufficient statistic.

---

## 8. Full-Rank, Curved, and Redundant Normal Submodels

The following three examples distinguish several concepts that are often confused.

### Case 1: The Full Normal Family

In the full family $N(\mu,\sigma^2)$, the natural parameter space is

\[
\mathbb{R}\times(-\infty,0),
\]

which is open. Therefore

\[
\left(\sum_iX_i,\sum_iX_i^2\right)
\]

is a minimal sufficient statistic.

### Case 2: The Curved Submodel $N(\sigma,\sigma^2)$

In the submodel $N(\sigma,\sigma^2)$ with $\sigma>0$, the natural parameters satisfy

\[
\eta_1=\frac{1}{\sigma},
\qquad
\eta_2=-\frac{1}{2\sigma^2}
=-\frac{1}{2}\eta_1^2.
\]

This parameter set is a curve, so the model is not full-rank. However, this curve contains three noncollinear points, so it has full affine hull in $\mathbb{R}^2$; Theorem 5 still applies.

### Case 3: The Redundant Representation $N(\sigma^2,\sigma^2)$

In the submodel $N(\sigma^2,\sigma^2)$, the coefficient of $x$ in the density is identically $1$. This term can be absorbed into the baseline density, leaving only $x^2$ as the statistic depending on the parameter. Therefore, the two-component representation is redundant, and

\[
\sum_iX_i^2
\]

is sufficient and minimal sufficient.

This distinction is very important: **the dimension of the original parameter, the number of displayed natural statistics, and the affine dimension of the natural parameter set need not be the same.**

---

## 9. What Sufficiency Does and Does Not Guarantee

**Sufficiency does not automatically imply minimal sufficiency.** If $T$ is sufficient, then for any statistic $W$, $(T,W)$ is also sufficient. Minimal sufficiency is responsible for removing such unnecessary additional information.

**Completeness is a separate concept.** Completeness studies whether, if $E_\theta a(T)=0$ for every $\theta$, it necessarily follows that $a(T)=0$ almost everywhere. This is a different property and will be developed in Chapter 3. Full-rank exponential families often have both minimal sufficiency and completeness, but their proofs and uses are not the same.

**Sufficiency depends on the model.** In an iid Bernoulli model, $\sum_iX_i$ is sufficient for the Bernoulli parameter. But if correlation, measurement error, or additional nuisance parameters are introduced, it may no longer be sufficient.

**Decision-theoretic preview.** The common kernel in Definition 1.1 allows us to condition on $T$ and average over the parameter-independent random variation in $X$. Chapter 3 will turn this observation into the Rao–Blackwell improvement, and after adding completeness, into the Lehmann–Scheffé theorem. The relevant proofs are left to the next chapter to avoid repetition.

---

## 10. Detailed Proof of the Countable Domination Lemma

It suffices to prove that statement 1 implies statement 3 in the Halmos–Savage lemma of Section 2.

Start from a $\sigma$-finite dominating measure $\mu$. After discarding null sets, choose a measurable partition

\[
\mathcal{X}=\bigcup_{m\geq1}A_m,
\qquad
0<\mu(A_m)<\infty,
\]

and define

\[
\nu(B)=\sum_{m=1}^{\infty}2^{-m}
\frac{\mu(B\cap A_m)}{\mu(A_m)}.
\]

Then $\nu$ is a probability measure and has exactly the same null sets as $\mu$. Therefore, every $P_\theta\ll\nu$.

Let $\mathcal{G}$ be the set of all countable mixtures:

\[
Q=\sum_{j=1}^{\infty}c_jP_{\theta_j},
\qquad
c_j>0,
\qquad
\sum_jc_j=1.
\]

For $Q\in\mathcal{G}$, write

\[
q=\frac{dQ}{d\nu},
\qquad
S_Q=\{x:q(x)>0\},
\]

and let

\[
a=\sup_{Q\in\mathcal{G}}\nu(S_Q)\leq1.
\]

Choose $Q_m\in\mathcal{G}$ such that

\[
\nu(S_{Q_m})>a-\frac{1}{m},
\]

and define

\[
\lambda=\sum_{m=1}^{\infty}2^{-m}Q_m.
\]

Expanding the double series as a single series shows that $\lambda$ itself is also a countable mixture of distributions in $\mathcal{P}$. If $q_m=dQ_m/d\nu$, then

\[
\frac{d\lambda}{d\nu}=\sum_m2^{-m}q_m,
\qquad
S_\lambda=\bigcup_mS_{Q_m}
\quad \nu\text{-a.e.}
\]

Because $\lambda\in\mathcal{G}$, the definition of $a$ gives $\nu(S_\lambda)\leq a$; the choice of $Q_m$ gives the reverse inequality. Hence

\[
\nu(S_\lambda)=a.
\]

Now fix an arbitrary $P\in\mathcal{P}$, let $p=dP/d\nu$, and suppose the set

\[
B=\{p>0\}\setminus S_\lambda
\]

has positive $\nu$-measure. Then $P(B)>0$, and the mixture

\[
Q=\frac{1}{2}(\lambda+P)\in\mathcal{G}
\]

has support $S_\lambda\cup\{p>0\}$, whose $\nu$-measure is strictly greater than $a$, contradicting the definition of $a$. Therefore

\[
p=0
\qquad \nu\text{-a.e. on }S_\lambda^c.
\]

If $\lambda(C)=0$, then since $d\lambda/d\nu$ is positive on $S_\lambda$, we must have

\[
\nu(C\cap S_\lambda)=0.
\]

Combining this with $p=0$ almost everywhere on $S_\lambda^c$, we obtain

\[
P(C)
=\int_{C\cap S_\lambda}p\,d\nu
+\int_{C\cap S_\lambda^c}p\,d\nu
=0.
\]

Thus $P\ll\lambda$. Since $P$ is arbitrary, $\lambda$ dominates the entire family, proving statement 3. Finally, listing the model distributions appearing in the expansion of the mixture $\lambda$ gives the countable subfamily in statement 2. $\square$