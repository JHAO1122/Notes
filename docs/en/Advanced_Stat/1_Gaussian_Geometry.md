# Chapter 1: Gaussian Geometry and Classical Sampling Distributions

This chapter starts from the projection geometry of Gaussian vectors and derives in a unified way classical sampling distributions such as the normal, chi-square, $t$, and $F$ distributions. The core idea is: **linear contrasts produce normal distributions, squared lengths of projected vectors produce chi-square distributions, studentization by an independently estimated scale produces $t$ distributions, and the ratio of two independent mean squares produces $F$ distributions.**

When the projected mean is nonzero, the corresponding chi-square, $t$, or $F$ distribution becomes a noncentral distribution. Hence, noncentral distributions naturally describe fixed alternative hypotheses and are used in power analysis and sample size calculation.

---

## 1. A One-Page Overview of the Statistical Problem

The statistical model is a family of distributions

\[
\mathcal{P}=\{P_\theta:\theta\in\Theta\},
\]

that describes the distributions the observed data $X$ may follow. The parameter $\theta$ is fixed but unknown, while the sample $X$ is random before it is observed. A **statistic** is a measurable function $T(X)$ of $X$, which must not contain unknown parameters.

The same statistic can serve different statistical inference tasks:

* **estimator** $\widehat{g}(X)$ used to approximate the target $g(\theta)$;
* **test** $\phi(X)\in[0,1]$ measures evidence against the null hypothesis, where $E_\theta\phi(X)$ equals the rejection probability;
* **confidence set** $C(X)$ has repeated-sampling coverage

\[
P_\theta\{\theta\in C(X)\}\geq 1-\alpha.
\]

In each case, statistical inference requires us to know the distribution of the statistic under $P_\theta$, which is the **sampling distribution** of the statistic.

!!! example "Example 1 (Running model for this chapter)"

    If $X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2)$, then $\theta=(\mu,\sigma^2)$. The statistics $\overline{X}$ and

    \[
    S^2=\frac{1}{n-1}\sum_{i=1}^n(X_i-\overline{X})^2
    \]

    estimate $\mu$ and $\sigma^2$, respectively.

    Exact inference about $\mu$ is feasible because the standardized $\overline{X}$ is normally distributed, the residual sum of squares is chi-square distributed, and the two are independent. The $t$ distribution is exactly the distribution of their studentized ratio.

---

## 2. Why These Distributions Keep Appearing

The normal, chi-square, $t$, and $F$ distributions are not a collection of distributions to be memorized in isolation. In the Gaussian model, they correspond respectively to four statistical operations:

\[
\begin{aligned}
\text{take a linear contrast} &\longrightarrow \text{normal distribution},\\
\text{take the squared length of a projected vector} &\longrightarrow \text{chi-square distribution},\\
\frac{\text{normal contrast}}{\text{independently estimated scale}} &\longrightarrow t\text{ distribution},\\
\text{compare two independent mean squares} &\longrightarrow F\text{ distribution}.
\end{aligned}
\]

| Distribution | Structural Form | Typical Exact Statistic | Main Use |
| --- | --- | --- | --- |
| Normal distribution | Gaussian linear contrast | Sample mean or regression coefficient under normal errors | Estimation and testing with known variance |
| Chi-square distribution | Squared norm of a central Gaussian projection | Residual sum of squares divided by $\sigma^2$ | Variance inference |
| Noncentral chi-square distribution | Squared norm retaining a nonzero mean component | Regression sum of squares under an alternative hypothesis | Signal strength and power |
| $t$ distribution | Central normal variable divided by an independent chi-square scale | Studentized mean or coefficient under the null hypothesis | Inference with unknown variance |
| Noncentral $t$ distribution | Shifted normal variable divided by an independent chi-square scale | Studentized mean or coefficient under an alternative hypothesis | Power and sample size |
| $F$ distribution | Ratio of two independent central mean squares | ANOVA or nested-model statistic under the null hypothesis | Testing multiple constraints |
| Noncentral $F$ distribution | Noncentral numerator divided by central residual mean square | ANOVA or nested-model statistic under an alternative hypothesis | Power and experimental design |

!!! warning "Exact versus Approximate Distributions"

    The results in this chapter are exact in finite samples under a normal model with covariance $\sigma^2I$. For nonnormal samples, the central limit theorem may make linear statistics approximately normal, but the exact independence between the estimated mean and the estimated variance usually no longer holds; in that case robust or asymptotic methods may be more appropriate.

    Therefore, we should always also ask: **What is the statistic? Which assumptions guarantee its reference distribution?**

---

## 3. Gaussian Vectors and Projection Geometry

Let $X\sim N_n(\mu,\Sigma)$. For a fixed matrix $A$ and a fixed vector $b$,

\[
AX+b\sim N_m(A\mu+b,A\Sigma A^T).
\]

The closure of the Gaussian distribution under linear transformations is the source of all exact results in this chapter.

!!! success "Theorem 1 (Zero Covariance Implies Independence in a Gaussian Vector)"

    Suppose $(U^T,V^T)^T$ is jointly Gaussian. Then

    \[
    U\perp V \quad\Longleftrightarrow\quad \operatorname{Cov}(U,V)=0.
    \]

??? proof "Proof of Theorem 1 (click to expand)"

    When second moments exist, independence always implies zero cross-covariance. We now prove the converse.

    Write the joint mean as $(\mu_U,\mu_V)$, and assume the cross-covariance is zero. For vectors $s$ and $t$, the joint characteristic function is

    \[
    \begin{aligned}
    E\exp\{i(s^TU+t^TV)\}
    &=\exp\left\{i(s^T\mu_U+t^T\mu_V)
    -\frac{1}{2}
    \begin{pmatrix}s\\t\end{pmatrix}^{T}
    \begin{pmatrix}\Sigma_U&0\\0&\Sigma_V\end{pmatrix}
    \begin{pmatrix}s\\t\end{pmatrix}\right\}\\
    &=\exp\left\{is^T\mu_U-\frac{1}{2}s^T\Sigma_Us\right\}
    \exp\left\{it^T\mu_V-\frac{1}{2}t^T\Sigma_Vt\right\}.
    \end{aligned}
    \]

    The last expression is the product of the two marginal characteristic functions. By uniqueness of characteristic functions, $U$ and $V$ are independent. $\square$

!!! success "Theorem 2 (Equivalent Characterization of Orthogonal Projections)"

    For an $n\times n$ real matrix $P$, the following three conditions are equivalent:

    1. $P=P^T=P^2$;
    2. there exists a subspace $\mathcal{S}\subseteq\mathbb{R}^n$ such that for every $x$, $Px$ is the orthogonal projection of $x$ onto $\mathcal{S}$;
    3. there exists an $n\times r$ matrix $U$ satisfying $U^TU=I_r$ such that $P=UU^T$.

    Under these conditions,

    \[
    \mathcal{S}=\operatorname{col}(P)=\operatorname{col}(U),
    \qquad
    \operatorname{rank}(P)=\operatorname{tr}(P)=r,
    \]

    and $I-P$ is the orthogonal projection onto $\mathcal{S}^{\perp}$.

??? proof "Proof of Theorem 2 (click to expand)"

    First assume condition 1 holds. If $y=Px\in\operatorname{col}(P)$, then

    \[
    Py=P^2x=Px=y.
    \]

    For any $z=Pw\in\operatorname{col}(P)$,

    \[
    z^T(x-Px)=w^TP^T(I-P)x=w^TP(I-P)x=0.
    \]

    Hence $Px\in\operatorname{col}(P)$, and $x-Px\perp\operatorname{col}(P)$, so condition 2 holds.

    Now assume condition 2 holds. Choose an orthonormal basis $u_1,\ldots,u_r$ of $\mathcal{S}$, and let $U=(u_1,\ldots,u_r)$. By the orthogonal projection formula,

    \[
    Px=\sum_{j=1}^ru_ju_j^Tx=UU^Tx,
    \]

    so condition 3 holds.

    Finally, if condition 3 holds, then

    \[
    P^T=UU^T=P,
    \qquad
    P^2=U(U^TU)U^T=P,
    \]

    so condition 1 holds.

    Every eigenvalue of an idempotent matrix satisfies $\lambda^2=\lambda$, so eigenvalues can only be $0$ or $1$. Since $P$ is symmetric, $P$ can be orthogonally diagonalized, and both its rank and trace equal the number of eigenvalues equal to $1$. The conclusion about $I-P$ follows directly. $\square$

!!! info "Lemma 3 (Nested Projections)"

    Let $P_0$ and $P_1$ be orthogonal projections. If $\operatorname{col}(P_0)\subseteq\operatorname{col}(P_1)$, then

    \[
    P_1P_0=P_0P_1=P_0.
    \]

    Hence $P_1-P_0$ is the orthogonal projection onto $\operatorname{col}(P_1)\cap\operatorname{col}(P_0)^\perp$, and it is orthogonal to both $P_0$ and $I-P_1$.

??? proof "Proof of Lemma 3 (click to expand)"

    For any $x$, $P_0x\in\operatorname{col}(P_0)\subseteq\operatorname{col}(P_1)$, so $P_1P_0x=P_0x$. Taking transposes on both sides gives $P_0P_1=P_0$. Thus

    \[
    (P_1-P_0)^2=P_1-P_1P_0-P_0P_1+P_0=P_1-P_0.
    \]

    Also, since $P_1-P_0$ is symmetric, it is a projection. Its image space is contained in $\operatorname{col}(P_1)$, and it annihilates $\operatorname{col}(P_0)$. At the same time,

    \[
    \operatorname{rank}(P_1-P_0)
    =\operatorname{tr}(P_1-P_0)
    =\operatorname{rank}(P_1)-\operatorname{rank}(P_0),
    \]

    so its image space is exactly $\operatorname{col}(P_1)\cap\operatorname{col}(P_0)^\perp$. Finally,

    \[
    P_0(P_1-P_0)=0,
    \qquad
    (I-P_1)(P_1-P_0)=0,
    \]

    which gives the required orthogonality. $\square$

---

## 4. Central and Noncentral Chi-Square Distributions

!!! info "Definition 4.1 (Noncentral Chi-Square Distribution)"

    Let $Z_1,\ldots,Z_r$ be mutually independent, with $Z_j\sim N(a_j,1)$. Define

    \[
    Q=\sum_{j=1}^rZ_j^2\sim\chi_r^2(\lambda),
    \qquad
    \lambda=\sum_{j=1}^ra_j^2.
    \]

    The number $\lambda\geq0$ is called the **noncentrality parameter**. When $\lambda=0$, it reduces to the central chi-square distribution $\chi_r^2$.

!!! note "Notation Convention"

    Some lecture notes write the noncentrality parameter as $\delta^2$ and denote the same distribution by $\chi_r^2(\delta^2)$. This chapter uses the common convention $\lambda=\delta^2$, with the standard terminology “noncentral.”

!!! success "Theorem 4 (Gaussian Quadratic Forms)"

    Let $X\sim N_n(\mu,\sigma^2I_n)$, where $\sigma^2>0$, and let $P$ be an orthogonal projection of rank $r$. Then

    \[
    \frac{X^TPX}{\sigma^2}\sim\chi_r^2(\lambda),
    \qquad
    \lambda=\frac{\mu^TP\mu}{\sigma^2}
    =\frac{\lVert P\mu\rVert^2}{\sigma^2}.
    \]

    In particular, this distribution is central chi-square if and only if $P\mu=0$.

??? proof "Proof of Theorem 4 (click to expand)"

    By the spectral theorem, there exists an orthogonal matrix $O=(O_1,O_0)$ such that

    \[
    O^TPO=
    \begin{pmatrix}
    I_r&0\\
    0&0
    \end{pmatrix},
    \]

    where the columns of $O_1$ form an orthonormal basis of $\operatorname{col}(P)$. Let $Z=O^TX/\sigma$. Then

    \[
    Z\sim N_n(O^T\mu/\sigma,I_n),
    \]

    and the coordinates of $Z$ are mutually independent. Therefore

    \[
    \frac{X^TPX}{\sigma^2}
    =Z^T
    \begin{pmatrix}
    I_r&0\\
    0&0
    \end{pmatrix}Z
    =\sum_{j=1}^rZ_j^2.
    \]

    By definition, it follows a noncentral chi-square distribution with noncentrality parameter

    \[
    \lambda
    =\sum_{j=1}^r\left(\frac{o_j^T\mu}{\sigma}\right)^2
    =\frac{\mu^TO_1O_1^T\mu}{\sigma^2}
    =\frac{\mu^TP\mu}{\sigma^2}.
    \]

    Since $P=P^T=P^2$, we have $\mu^TP\mu=\lVert P\mu\rVert^2$. This value is zero if and only if $P\mu=0$. $\square$

!!! success "Corollary 5 (Mean and Variance)"

    If $Q\sim\chi_r^2(\lambda)$, then

    \[
    EQ=r+\lambda,
    \qquad
    \operatorname{Var}(Q)=2(r+2\lambda).
    \]

??? proof "Proof of Corollary 5 (click to expand)"

    Write $Q=\sum_{j=1}^rZ_j^2$, where the mutually independent $Z_j\sim N(a_j,1)$ and $\sum_ja_j^2=\lambda$. Because

    \[
    EZ_j^2=1+a_j^2,
    \qquad
    EZ_j^4=3+6a_j^2+a_j^4,
    \]

    so

    \[
    \operatorname{Var}(Z_j^2)
    =EZ_j^4-(EZ_j^2)^2
    =2+4a_j^2.
    \]

    Summing using independence gives the result. $\square$

---

## 5. Orthogonal Gaussian Components and Cochran's Theorem

!!! success "Theorem 6 (Independence of Orthogonal Projections)"

    Let $X\sim N_n(\mu,\sigma^2I_n)$. If $P$ and $Q$ are orthogonal projections satisfying $PQ=0$, then $PX$ and $QX$ are independent. Consequently, $X^TPX$ and $X^TQX$ are also independent.

??? proof "Proof of Theorem 6 (click to expand)"

    Because $P$ and $Q$ are both symmetric,

    \[
    QP=(PQ)^T=0.
    \]

    The vector pair $(PX,QX)$ is jointly Gaussian, and

    \[
    \operatorname{Cov}(PX,QX)
    =P(\sigma^2I_n)Q^T
    =\sigma^2PQ
    =0.
    \]

    By Theorem 1, $PX\perp QX$. Also, because

    \[
    X^TPX=\lVert PX\rVert^2,
    \qquad
    X^TQX=\lVert QX\rVert^2,
    \]

    the two quadratic forms are functions of two independent vectors, respectively, so they are also independent. $\square$

!!! success "Theorem 7 (Cochran's Theorem: Projection Form)"

    Let $P_1,\ldots,P_k$ be pairwise orthogonal projections, that is, $P_iP_j=0$ when $i\neq j$. If $X\sim N_n(\mu,\sigma^2I_n)$, then

    \[
    Q_j=\frac{X^TP_jX}{\sigma^2},
    \qquad j=1,\ldots,k,
    \]

    are mutually independent, and

    \[
    Q_j\sim\chi_{r_j}^2(\lambda_j),
    \qquad
    r_j=\operatorname{rank}(P_j),
    \qquad
    \lambda_j=\frac{\mu^TP_j\mu}{\sigma^2}.
    \]

    If, moreover, $\sum_jP_j=I_n$, then

    \[
    \frac{\lVert X\rVert^2}{\sigma^2}=\sum_jQ_j,
    \qquad
    \sum_jr_j=n,
    \qquad
    \sum_j\lambda_j=\frac{\lVert\mu\rVert^2}{\sigma^2}.
    \]

??? proof "Proof of Theorem 7 (click to expand)"

    Theorem 4 gives the marginal distribution of each $Q_j$. Stack the projected vectors as

    \[
    W=(P_1X,\ldots,P_kX).
    \]

    It is a jointly Gaussian vector, and when $i\neq j$, its $(i,j)$th cross-covariance block is

    \[
    \sigma^2P_iP_j=0.
    \]

    Hence the covariance matrix of $W$ is block diagonal. Its characteristic function factors into the product of the marginal characteristic functions, so the projected vectors are mutually independent, and their squared norms are also mutually independent.

    If $\sum_jP_j=I_n$, then

    \[
    \lVert X\rVert^2
    =X^T\left(\sum_jP_j\right)X
    =\sum_jX^TP_jX.
    \]

    Taking traces gives $\sum_jr_j=n$, and

    \[
    \sum_j\lambda_j
    =\frac{\mu^T(\sum_jP_j)\mu}{\sigma^2}
    =\frac{\lVert\mu\rVert^2}{\sigma^2}.
    \]

    This completes the proof. $\square$

!!! success "Theorem 8 (Cochran's Theorem: Rank Form)"

    Let $A_1,\ldots,A_k$ be symmetric positive semidefinite matrices satisfying

    \[
    \sum_jA_j=I_n.
    \]

    Write $r_j=\operatorname{rank}(A_j)$. If $\sum_jr_j=n$, then each $A_j$ is an orthogonal projection, and these projections are pairwise orthogonal. Therefore, for $X\sim N_n(\mu,\sigma^2I_n)$,

    \[
    \frac{X^TA_jX}{\sigma^2}
    \sim
    \chi_{r_j}^2\left(\frac{\mu^TA_j\mu}{\sigma^2}\right),
    \]

    and these quadratic forms are mutually independent.

??? proof "Proof of Theorem 8 (click to expand)"

    Because $A_j\succeq0$ and

    \[
    I_n-A_j=\sum_{\ell\neq j}A_\ell\succeq0,
    \]

    every eigenvalue of $A_j$ lies in $[0,1]$. Thus

    \[
    \operatorname{tr}(A_j)\leq\operatorname{rank}(A_j)=r_j.
    \]

    But

    \[
    n=\operatorname{tr}(I_n)
    =\sum_j\operatorname{tr}(A_j)
    \leq\sum_jr_j
    =n.
    \]

    Therefore all inequalities must be equalities. For each $j$, the $r_j$ nonzero eigenvalues of $A_j$ lie in $(0,1]$ and sum to $r_j$, so they all equal $1$, and hence $A_j^2=A_j$.

    Fix $i\neq j$, and take $x\in\operatorname{col}(A_j)$. Since $A_j$ is already a projection, $A_jx=x$, so

    \[
    0=x^T(I_n-A_j)x
    =\sum_{\ell\neq j}x^TA_\ell x.
    \]

    Each term is nonnegative, so $x^TA_ix=0$. Also, because $A_i$ is a projection, $x^TA_ix=\lVert A_ix\rVert^2$, so $A_ix=0$. This holds for all $x\in\operatorname{col}(A_j)$, so $A_iA_j=0$. Finally, Theorem 7 yields the distribution and independence conclusions. $\square$

!!! note "Why the Rank Condition Matters"

    The rank condition forces these quadratic forms to use exactly $n$ mutually perpendicular Gaussian coordinates: they neither reuse any direction nor omit any direction.

---

## 6. Normal Samples: Chi-Square and Student Distributions

Let $X=(X_1,\ldots,X_n)^T$, where $X_i\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2)$, and define

\[
P_1=\frac{1}{n}\mathbf{1}\mathbf{1}^T,
\qquad
M=I_n-P_1.
\]

$P_1$ and $M$ project onto $\operatorname{span}(\mathbf{1})$ and its orthogonal complement, respectively, with ranks $1$ and $n-1$.

!!! success "Theorem 9 (Mean–Variance Decomposition for Normal Samples)"

    For a normal sample,

    \[
    \overline{X}\sim N\left(\mu,\frac{\sigma^2}{n}\right),
    \qquad
    \overline{X}\perp S^2,
    \qquad
    \frac{(n-1)S^2}{\sigma^2}\sim\chi_{n-1}^2.
    \]

??? proof "Proof of Theorem 9 (click to expand)"

    Because $\overline{X}=n^{-1}\mathbf{1}^TX$ is a linear Gaussian statistic,

    \[
    \overline{X}
    \sim
    N\left(n^{-1}\mathbf{1}^T(\mu\mathbf{1}),
    n^{-2}\mathbf{1}^T(\sigma^2I_n)\mathbf{1}\right)
    =N\left(\mu,\frac{\sigma^2}{n}\right).
    \]

    In addition,

    \[
    P_1X=\overline{X}\mathbf{1},
    \qquad
    MX=X-\overline{X}\mathbf{1},
    \qquad
    X^TMX=\sum_{i=1}^n(X_i-\overline{X})^2.
    \]

    Because $P_1M=0$, Theorem 6 shows that $P_1X$ and $MX$ are independent, so $\overline{X}$ and $S^2$ are independent. Finally, $M$ has rank $n-1$, and $M(\mu\mathbf{1})=0$, so by Theorem 4,

    \[
    \frac{X^TMX}{\sigma^2}
    =\frac{(n-1)S^2}{\sigma^2}
    \sim\chi_{n-1}^2.
    \]

    This completes the proof. $\square$

!!! info "Definition 6.1 (Central and Noncentral Student Distributions)"

    Let $Z\sim N(\delta,1)$ and $V\sim\chi_\nu^2$, and let the two be independent. Then

    \[
    T=\frac{Z}{\sqrt{V/\nu}}\sim t_\nu(\delta).
    \]

    When $\delta=0$, it is the central Student distribution $t_\nu$.

!!! success "Corollary 10 (One-Sample Statistic under the Null and Alternative Hypotheses)"

    For any reference value $\mu_0$,

    \[
    T=\frac{\sqrt{n}(\overline{X}-\mu_0)}{S}
    \sim t_{n-1}(\delta),
    \qquad
    \delta=\frac{\sqrt{n}(\mu-\mu_0)}{\sigma}.
    \]

    Under $H_0:\mu=\mu_0$, $T\sim t_{n-1}$; under a fixed alternative hypothesis, $T$ follows a noncentral $t$ distribution.

??? proof "Proof of Corollary 10 (click to expand)"

    Let

    \[
    Z=\frac{\sqrt{n}(\overline{X}-\mu_0)}{\sigma},
    \qquad
    V=\frac{(n-1)S^2}{\sigma^2}.
    \]

    Theorem 9 gives

    \[
    Z\sim N\left(\frac{\sqrt{n}(\mu-\mu_0)}{\sigma},1\right),
    \qquad
    V\sim\chi_{n-1}^2,
    \qquad
    Z\perp V.
    \]

    Because

    \[
    \frac{Z}{\sqrt{V/(n-1)}}
    =\frac{\sqrt{n}(\overline{X}-\mu_0)}{S},
    \]

    the conclusion follows directly from the definition. $\square$

---

## 7. Central and Noncentral $F$ Distributions

!!! info "Definition 7.1 (Central and Singly Noncentral $F$ Distributions)"

    Let $U\sim\chi_r^2(\lambda)$ and $V\sim\chi_s^2$, and let the two be independent. Then

    \[
    F=\frac{U/r}{V/s}\sim F_{r,s}(\lambda).
    \]

    When $\lambda=0$, it is the central $F$ distribution $F_{r,s}$. In the standard singly noncentral definition, the denominator is a central chi-square variable; if both numerator and denominator are noncentral chi-square variables, the ratio follows a doubly noncentral $F$ distribution, which this chapter does not need to use.

### 7.1 Variances of Two Independent Samples

Let

\[
X_1,\ldots,X_m\overset{\mathrm{iid}}{\sim}N(\mu_X,\sigma_X^2),
\qquad
Y_1,\ldots,Y_n\overset{\mathrm{iid}}{\sim}N(\mu_Y,\sigma_Y^2),
\]

and let the two samples be independent.

!!! success "Proposition 11 (Variance Ratio)"

    We have

    \[
    \frac{S_X^2/\sigma_X^2}{S_Y^2/\sigma_Y^2}
    \sim F_{m-1,n-1},
    \qquad
    \frac{S_X^2}{S_Y^2}
    \sim\frac{\sigma_X^2}{\sigma_Y^2}F_{m-1,n-1}.
    \]

    Under $H_0:\sigma_X^2=\sigma_Y^2$, the unstandardized variance ratio follows a central $F_{m-1,n-1}$ distribution.

??? proof "Proof of Proposition 11 (click to expand)"

    By Theorem 9, the following two variables are independent:

    \[
    U=\frac{(m-1)S_X^2}{\sigma_X^2}\sim\chi_{m-1}^2,
    \qquad
    V=\frac{(n-1)S_Y^2}{\sigma_Y^2}\sim\chi_{n-1}^2.
    \]

    Therefore

    \[
    \frac{U/(m-1)}{V/(n-1)}
    =\frac{S_X^2/\sigma_X^2}{S_Y^2/\sigma_Y^2}
    \sim F_{m-1,n-1}.
    \]

    This completes the proof. $\square$

!!! warning "An Important Distinction"

    When $\sigma_X^2\neq\sigma_Y^2$, the unstandardized variance ratio follows a **scaled central $F$ distribution**, not a noncentral $F$ distribution. Noncentrality comes from a nonzero Gaussian mean in a squared projection; different variances only change the scale.

### 7.2 Nested Normal Linear Models

This is the main reason for introducing the noncentral $F$ distribution. Let $Y\sim N_n(\mu,\sigma^2I_n)$, and consider nested design spaces

\[
\mathcal{S}_R\subseteq\mathcal{S}_F,
\qquad
\dim(\mathcal{S}_R)=p_R,
\qquad
\dim(\mathcal{S}_F)=p_F,
\]

with projection matrices $H_R,H_F$, respectively, and let $q=p_F-p_R$.

!!! success "Theorem 12 ($F$ Statistic for Nested Models)"

    Assume $\mu\in\mathcal{S}_F$ and $n>p_F$. Then

    \[
    F=
    \frac{Y^T(H_F-H_R)Y/q}
    {Y^T(I_n-H_F)Y/(n-p_F)}
    \sim F_{q,n-p_F}(\lambda),
    \]

    where

    \[
    \lambda
    =\frac{\mu^T(H_F-H_R)\mu}{\sigma^2}
    =\frac{\lVert(H_F-H_R)\mu\rVert^2}{\sigma^2}.
    \]

    For testing $H_0:\mu\in\mathcal{S}_R$ against $H_1:\mu\in\mathcal{S}_F\setminus\mathcal{S}_R$, under the null hypothesis the statistic follows a central $F_{q,n-p_F}$ distribution; under a fixed alternative hypothesis, $\lambda>0$.

??? proof "Proof of Theorem 12 (click to expand)"

    By Lemma 3,

    \[
    P_N=H_F-H_R,
    \qquad
    P_D=I_n-H_F
    \]

    are orthogonal projections of ranks $q$ and $n-p_F$, respectively, and $P_NP_D=0$. By Theorem 4 and Theorem 6, the following two variables are independent:

    \[
    U=\frac{Y^TP_NY}{\sigma^2}\sim\chi_q^2(\lambda),
    \qquad
    V=\frac{Y^TP_DY}{\sigma^2}\sim\chi_{n-p_F}^2(\lambda_D).
    \]

    Because $\mu\in\mathcal{S}_F$, we have $P_D\mu=0$, so $\lambda_D=0$. Therefore

    \[
    \frac{U/q}{V/(n-p_F)}\sim F_{q,n-p_F}(\lambda).
    \]

    If $\mu\in\mathcal{S}_R$, then $H_R\mu=H_F\mu=\mu$, so $\lambda=0$. Conversely, within $\mathcal{S}_F$ there is the orthogonal decomposition

    \[
    \mu=H_R\mu+(H_F-H_R)\mu.
    \]

    If $\mu\notin\mathcal{S}_R$, then the second component is nonzero, so $\lambda>0$. $\square$

!!! example "Example 2 (Testing a Set of Regression Coefficients)"

    In the full model

    \[
    Y=X_1\beta_1+X_2\beta_2+\varepsilon,
    \qquad
    \varepsilon\sim N_n(0,\sigma^2I_n)
    \]

    test $H_0:\beta_2=0$. Suppose $X_2$ adds $q$ linearly independent directions outside $\operatorname{col}(X_1)$, and take

    \[
    \mathcal{S}_R=\operatorname{col}(X_1),
    \qquad
    \mathcal{S}_F=\operatorname{col}(X_1,X_2).
    \]

    The numerator of Theorem 12 is the increase in fitted sum of squares after adding $X_2$, and the denominator is the residual mean square of the full model. Under the alternative hypothesis, the component of $X_2\beta_2$ that cannot be explained by $X_1$ makes $\lambda>0$. Therefore

    \[
    P_{\beta_2}\{F>F_{q,n-p_F;1-\alpha}\}
    \]

    is the exact power; conversely, solving the equation for this probability with respect to $n$ allows sample size calculation.

**Geometric summary:** The data admit the orthogonal decomposition

\[
Y=H_RY+(H_F-H_R)Y+(I_n-H_F)Y,
\]

where the three terms correspond in order to the **reduced-model fit, extra fit, and full-model residual**. Cochran's theorem converts the squared lengths of these mutually perpendicular Gaussian components into independent chi-square variables, and the $F$ statistic compares the extra fitted mean square with the residual mean square.

---

## Appendix A: Moment Generating Function and Poisson Mixture Representation

!!! info "Lemma 13"

    If $Z\sim N(a,1)$, then for $t<1/2$,

    \[
    Ee^{tZ^2}
    =(1-2t)^{-1/2}
    \exp\left\{\frac{a^2t}{1-2t}\right\}.
    \]

??? proof "Proof of Lemma 13 (click to expand)"

    By completing the square,

    \[
    \begin{aligned}
    Ee^{tZ^2}
    &=\frac{e^{-a^2/2}}{\sqrt{2\pi}}
    \int_{\mathbb{R}}
    \exp\left\{-\frac{1-2t}{2}z^2+az\right\}\,dz\\
    &=\frac{e^{-a^2/2}}{\sqrt{2\pi}}
    \exp\left\{\frac{a^2}{2(1-2t)}\right\}
    \int_{\mathbb{R}}
    \exp\left\{-\frac{1-2t}{2}
    \left(z-\frac{a}{1-2t}\right)^2\right\}\,dz\\
    &=(1-2t)^{-1/2}
    \exp\left\{\frac{a^2t}{1-2t}\right\}.
    \end{aligned}
    \]

    This completes the proof. $\square$

!!! success "Proposition 14 (Moment Generating Function, Additivity, and Poisson Mixture)"

    If $Q\sim\chi_r^2(\lambda)$, then

    \[
    M_Q(t)
    =(1-2t)^{-r/2}
    \exp\left\{\frac{\lambda t}{1-2t}\right\}.
    \]

    Therefore, when mutually independent noncentral chi-square variables are added, their degrees of freedom and noncentrality parameters add separately.

    In addition, if

    \[
    K\sim\operatorname{Poisson}(\lambda/2),
    \qquad
    Q\mid K\sim\chi_{r+2K}^2,
    \]

    then the marginal distribution of $Q$ is $\chi_r^2(\lambda)$.

??? proof "Proof of Proposition 14 (click to expand)"

    Applying Lemma 13 to each independent coordinate and multiplying the moment generating functions gives the moment generating function formula above. The product form of the moment generating function also directly proves additivity.

    For the Poisson mixture representation,

    \[
    \begin{aligned}
    Ee^{tQ}
    &=(1-2t)^{-r/2}E(1-2t)^{-K}\\
    &=(1-2t)^{-r/2}
    \exp\left\{\frac{\lambda}{2}\left[(1-2t)^{-1}-1\right]\right\}\\
    &=(1-2t)^{-r/2}
    \exp\left\{\frac{\lambda t}{1-2t}\right\}.
    \end{aligned}
    \]

    By uniqueness of the moment generating function, $Q\sim\chi_r^2(\lambda)$. $\square$

---

## Appendix B: Converse for Gaussian Quadratic Forms

!!! success "Theorem 15 (Necessity of Being a Projection)"

    Let $Z\sim N_n(0,I_n)$, and let $A$ be a real symmetric matrix. If

    \[
    Z^TAZ\sim\chi_r^2,
    \]

    then $A$ must be an orthogonal projection of rank $r$.

??? proof "Proof of Theorem 15 (click to expand)"

    Diagonalize $A$ as

    \[
    A=O\operatorname{diag}(a_1,\ldots,a_n)O^T.
    \]

    Let $W=O^TZ$. Then $W\sim N_n(0,I_n)$, and

    \[
    Z^TAZ=\sum_{j=1}^na_jW_j^2.
    \]

    There cannot be any negative $a_j$. For example, if $a_1<0$, then the event “$W_1^2$ is sufficiently large while the other squared terms remain bounded” has positive probability and makes the quadratic form negative, which contradicts the nonnegativity of a chi-square variable.

    For $t$ near zero, equality of moment generating functions gives

    \[
    \prod_{j=1}^n(1-2a_jt)^{-1/2}
    =(1-2t)^{-r/2}.
    \]

    Taking reciprocals and then squaring both sides gives the polynomial identity

    \[
    \prod_{j=1}^n(1-2a_jt)=(1-2t)^r.
    \]

    Every positive $a_j$ produces a root $(2a_j)^{-1}$ on the left; the right side has only the root $1/2$, with multiplicity $r$. Therefore exactly $r$ eigenvalues equal $1$, and the remaining eigenvalues equal $0$. Hence $A^2=A$, and together with symmetry this implies that $A$ is an orthogonal projection of rank $r$. $\square$

!!! note "Why Symmetry Is Natural"

    Every quadratic form depends only on the symmetric part of the matrix, because

    \[
    x^TAx=x^T\left(\frac{A+A^T}{2}\right)x.
    \]

---

## Appendix C: Density Derivations for Central Distributions

### C.1 Chi-Square Distribution

If $Z\sim N(0,1)$ and $W=Z^2$, using the two inverse-function branches $z=\pm\sqrt{w}$, we obtain

\[
f_W(w)
=\frac{\phi(\sqrt{w})+\phi(-\sqrt{w})}{2\sqrt{w}}
=\frac{w^{-1/2}e^{-w/2}}{2^{1/2}\Gamma(1/2)},
\qquad w>0.
\]

Therefore $Z^2$ follows a Gamma distribution with shape parameter $1/2$ and scale parameter $2$. Adding independent Gamma variables with the same scale parameter gives

\[
f_{\chi_\nu^2}(v)
=\frac{v^{\nu/2-1}e^{-v/2}}{2^{\nu/2}\Gamma(\nu/2)},
\qquad v>0.
\]

### C.2 Student Distribution

Let $Z\sim N(0,1)$ and $V\sim\chi_\nu^2$, and let them be independent. For

\[
T=\frac{Z}{\sqrt{V/\nu}},
\]

make the change of variables $z=t\sqrt{v/\nu}$, whose Jacobian is $\sqrt{v/\nu}$. Then

\[
\begin{aligned}
f_T(t)
&=\int_0^\infty
\phi\left(t\sqrt{v/\nu}\right)
f_{\chi_\nu^2}(v)
\sqrt{\frac{v}{\nu}}\,dv\\
&=\frac{\Gamma((\nu+1)/2)}{\sqrt{\nu\pi}\,\Gamma(\nu/2)}
\left(1+\frac{t^2}{\nu}\right)^{-(\nu+1)/2},
\qquad t\in\mathbb{R}.
\end{aligned}
\]

The last step follows from the Gamma integral.

### C.3 $F$ Distribution

Let $U\sim\chi_r^2$ and $V\sim\chi_s^2$, and let them be independent. For

\[
F=\frac{U/r}{V/s},
\]

make the change of variables $u=(r/s)fv$ while retaining $v$, whose Jacobian is $(r/s)v$. Using the Gamma integral gives

\[
f_F(f)
=\frac{\Gamma((r+s)/2)}{\Gamma(r/2)\Gamma(s/2)}
\left(\frac{r}{s}\right)^{r/2}
f^{r/2-1}
\left(1+\frac{r}{s}f\right)^{-(r+s)/2},
\qquad f>0.
\]

---

## Appendix D: Proof Checklist and Distribution Selection Checklist

1. Is the vector **exactly Gaussian**, or only asymptotically normal?
2. Is the covariance matrix $\sigma^2I$? If not, can the vector be whitened?
3. Is the quadratic-form matrix symmetric and idempotent? Its degrees of freedom equal the rank of the matrix.
4. Compute

\[
\lambda=\frac{\mu^TP\mu}{\sigma^2};
\]

do not incorrectly set it to zero when away from the null hypothesis.
5. For Gaussian projections, check $P_iP_j=0$ before asserting independence.
6. Write the $F$ statistic as $(U/r)/(V/s)$, and confirm that in the usual noncentral $F$ distribution the denominator is a central chi-square variable.
7. Clarify whether the distribution is used for null calibration, confidence coverage, power, or sample size calculation.