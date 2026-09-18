# 第二章：充分性、最小充分性与指数族

本章研究统计推断中的精确数据降维问题。我们首先用共同条件分布核定义充分统计量，随后证明 Fisher–Neyman 分解定理，并通过 Bernoulli、Poisson、正态、均匀和 Gamma 型模型说明如何使用分解准则。接着，我们介绍指数族及其自然参数空间的几何结构，讨论最小充分性，并证明满仿射秩指数族的自然统计量是最小充分统计量。

本章的核心思想是：**如果样本对参数的全部依赖都通过统计量 $T$ 传递，那么在给定 $T$ 后，样本中剩余的随机性便与参数无关。**

---

## 1. 为什么充分性是一条精确的数据降维原理

设

\[
\mathcal{P}=\{P_\theta:\theta\in\Theta\}
\]

是可测样本空间 $(\mathcal{X},\mathcal{B})$ 上的统计模型，并设 $T:\mathcal{X}\to\mathcal{T}$ 是一个统计量。基本问题是：观测 $T(X)$ 是否保留了观测完整样本 $X$ 所包含的关于 $\theta$ 的全部信息？

!!! info "定义 1.1（通过共同核定义充分性）"

    假设所涉及的空间均为标准 Borel 空间，因此正则条件分布存在。如果存在一个不依赖于 $\theta$ 的单一 Markov 核 $K(t,B)$，使得对每个 $B\in\mathcal{B}$ 和每个 $\theta$，都有

    \[
    P_\theta(X\in B\mid T)=K\{T(X),B\}
    \qquad P_\theta\text{-a.s.},
    \]

    则称统计量 $T$ 对模型 $\mathcal{P}$ 是**充分的（sufficient）**。

    因而，一旦 $T$ 已知，$X$ 中剩余的随机性便不再依赖于参数。“单一核”这几个字非常重要：在每个 $P_\theta$ 下，条件概率只在几乎处处意义下定义，而充分性要求我们能够选择在整个模型中一致的条件分布版本。

!!! example "例 1（Bernoulli 样本）"

    设

    \[
    X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
    \qquad
    T=\sum_{i=1}^nX_i.
    \]

    若 $x\in\{0,1\}^n$ 且 $\sum_i x_i=t$，则

    \[
    P_p(X=x\mid T=t)
    =\frac{p^t(1-p)^{n-t}}{\binom{n}{t}p^t(1-p)^{n-t}}
    =\binom{n}{t}^{-1}.
    \]

    条件分布在所有恰好含有 $t$ 个 $1$ 的二元序列上均匀分布，因此与 $p$ 无关。统计量 $T$ 舍弃了成功出现的次序，却没有舍弃任何关于 $p$ 的信息。

在上例中，直接计算条件分布非常直观，但对于连续数据并不方便，对于结构复杂的样本则几乎不可行。Fisher–Neyman 分解定理将给出实际可用的判别准则。

---

## 2. 控制引理：为什么一个参考测度就足够

假设每个 $P_\theta$ 都被同一个 $\sigma$-有限测度 $\mu$ 控制。证明 Fisher–Neyman 分解定理需要下面的结果，它通常称为 **Halmos–Savage 引理**。

!!! info "引理（可数控制 / Halmos–Savage）"

    对一个概率测度族 $\mathcal{P}$，下列三个命题等价：

    1. 存在一个 $\sigma$-有限测度 $\mu$，使得对每个 $\theta$ 都有 $P_\theta\ll\mu$；
    2. 存在 $P_{\theta_1},P_{\theta_2},\ldots\in\mathcal{P}$，使得对每个 $A\in\mathcal{B}$，

    \[
    \bigl[P_\theta(A)=0\text{ 对所有 }\theta\bigr]
    \quad\Longleftrightarrow\quad
    \bigl[P_{\theta_j}(A)=0\text{ 对所有 }j\bigr];
    \]

    3. 存在形如

    \[
    \lambda=\sum_{j=1}^{\infty}w_jP_{\theta_j},
    \qquad
    w_j>0,
    \qquad
    \sum_{j=1}^{\infty}w_j=1
    \]

    的概率测度，使得对每个 $\theta$ 都有 $P_\theta\ll\lambda$。

非平凡的方向是由命题 1 推出命题 3，其完整证明见第 10 节。这个结果的要点是：即使参数空间 $\Theta$ 不可数，模型中可数多个分布的混合仍能捕捉整个统计试验涉及的全部零测集。

**两个简单的推论方向：** 若命题 2 成立，可以在混合分布中取 $w_j=2^{-j}$。此时 $\lambda(A)=0$ 当且仅当所有选出的 $P_{\theta_j}(A)=0$，由命题 2 可知所有 $P_\theta(A)=0$，所以命题 2 推出命题 3。命题 3 推出命题 1 则是直接的，只需令 $\mu=\lambda$。

---

## 3. Fisher–Neyman 分解定理

令

\[
p_\theta=\frac{dP_\theta}{d\mu}.
\]

我们允许分布的支撑依赖于参数；密度 $p_\theta$ 中的零值也是分解的一部分。

!!! success "定理 1（Fisher–Neyman 分解定理）"

    在上述假设下，$T$ 对 $\mathcal{P}$ 充分，当且仅当存在非负可测函数 $g_\theta$ 和 $h$，使得

    \[
    p_\theta(x)=g_\theta\{T(x)\}h(x)
    \qquad \mu\text{-a.e.},
    \]

    其中 $h$ 与 $\theta$ 无关。

??? proof "定理 1 的证明（点击展开）"

    令

    \[
    \lambda=\sum_jw_jP_{\theta_j}
    \]

    为控制引理给出的可数混合分布。于是对每个 $\theta$，都有

    \[
    P_\theta\ll\lambda\ll\mu.
    \]

    **充分性推出分解。** 设 $K(T,B)$ 是定义 1.1 中的共同核。因为在每个选出的 $P_{\theta_j}$ 下，同一个核都是条件分布，所以取混合后，它也是 $\lambda$ 下的条件分布：

    \[
    K(T,B)=E_\lambda(\mathbf{1}_B\mid T)
    \qquad \lambda\text{-a.s.}
    \]

    记

    \[
    f_\theta=\frac{dP_\theta}{d\lambda}.
    \]

    对任意 $B\in\mathcal{B}$，

    \[
    \begin{aligned}
    P_\theta(B)
    &=\int K(T,B)\,dP_\theta\\
    &=\int E_\lambda(\mathbf{1}_B\mid T)f_\theta\,d\lambda\\
    &=\int E_\lambda(\mathbf{1}_B\mid T)E_\lambda(f_\theta\mid T)\,d\lambda\\
    &=\int_B E_\lambda(f_\theta\mid T)\,d\lambda.
    \end{aligned}
    \]

    由 Radon–Nikodym 导数的唯一性，

    \[
    f_\theta=E_\lambda(f_\theta\mid T)
    \qquad \lambda\text{-a.s.}
    \]

    因而 $f_\theta$ 是 $\sigma(T)$-可测的。由 Doob–Dynkin 引理，存在可测函数 $g_\theta$，使得

    \[
    f_\theta(x)=g_\theta\{T(x)\}.
    \]

    若令 $h=d\lambda/d\mu$，则由 Radon–Nikodym 导数的链式法则，

    \[
    p_\theta(x)
    =\frac{dP_\theta}{d\lambda}(x)
    \frac{d\lambda}{d\mu}(x)
    =g_\theta\{T(x)\}h(x),
    \]

    从而得到所需分解。

    **分解推出充分性。** 假设

    \[
    p_\theta(x)=g_\theta\{T(x)\}h(x).
    \]

    对同一个混合分布 $\lambda$，有

    \[
    p_\lambda(x)
    =\sum_jw_jp_{\theta_j}(x)
    =h(x)q\{T(x)\},
    \qquad
    q(t)=\sum_jw_jg_{\theta_j}(t).
    \]

    因为每个 $P_\theta\ll\lambda$，所以在 $\{q>0\}$ 上，

    \[
    \frac{dP_\theta}{d\lambda}(x)
    =r_\theta\{T(x)\},
    \qquad
    r_\theta(t)=\frac{g_\theta(t)}{q(t)};
    \]

    在其他位置任意定义该比值。令 $K_\lambda(t,B)$ 是 $\lambda$ 下给定 $T=t$ 时 $X$ 的一个正则条件分布。若 $C\in\sigma(T)$，则 $\mathbf{1}_Cr_\theta(T)$ 是 $\sigma(T)$-可测的，并且

    \[
    \begin{aligned}
    \int_CK_\lambda(T,B)\,dP_\theta
    &=\int_CK_\lambda(T,B)r_\theta(T)\,d\lambda\\
    &=\int\mathbf{1}_C\mathbf{1}_Br_\theta(T)\,d\lambda\\
    &=P_\theta(B\cap C).
    \end{aligned}
    \]

    因此，$K_\lambda(T,B)$ 在每个 $P_\theta$ 下也都是事件 $B$ 的条件概率。它是一个共同且与参数无关的核，所以 $T$ 充分。$\square$

!!! note "证明的解释"

    在正向证明中，充分性迫使每个似然比 $dP_\theta/d\lambda$ 都成为 $T$ 的函数。在反向证明中，这些似然比只重新加权 $T$ 的分布，而不改变给定 $T$ 后 $X$ 的条件分布。这就是“所有参数依赖都通过 $T$ 传递”的严格含义。

---

## 4. 使用分解定理：四个例子

### 4.1 Bernoulli 与 Poisson 样本

若

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Bernoulli}(p),
\]

则

\[
p_p(x)=p^{\sum_i x_i}(1-p)^{n-\sum_i x_i},
\]

所以 $\sum_iX_i$ 是充分统计量。

类似地，若

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Poisson}(\lambda),
\]

则

\[
p_\lambda(x)
=e^{-n\lambda}\lambda^{\sum_i x_i}
\prod_{i=1}^n\frac{1}{x_i!},
\]

因此 $\sum_iX_i$ 对 $\lambda$ 充分。

### 4.2 正态分布的均值与方差

若

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2),
\]

则联合密度可以写成

\[
p_{\mu,\sigma^2}(x)
=(2\pi\sigma^2)^{-n/2}
\exp\left\{-\frac{1}{2\sigma^2}
\left(\sum_i x_i^2-2\mu\sum_i x_i+n\mu^2\right)\right\}.
\]

因此

\[
S(X)=\left(\sum_iX_i,\sum_iX_i^2\right)
\]

对 $(\mu,\sigma^2)$ 充分。又因为

\[
\sum_iX_i=n\overline{X},
\qquad
\sum_iX_i^2=(n-1)S^2+n\overline{X}^{,2},
\]

所以当 $n\geq2$ 时，$S(X)$ 与 $(\overline{X},S^2)$ 等价。

### 4.3 均匀分布的端点：依赖参数的支撑

若

\[
X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}\operatorname{Uniform}(0,\theta),
\]

则

\[
p_\theta(x)
=\theta^{-n}\mathbf{1}\{0<x_{(1)}\}
\mathbf{1}\{x_{(n)}<\theta\}.
\]

因此 $X_{(n)}$ 是充分统计量。包含 $\theta$ 的示性函数属于 $g_\theta\{X_{(n)}\}$；分解定理并不要求模型具有共同支撑。不过，由于该模型的支撑随 $\theta$ 改变，后续的正则似然理论必须单独处理它。

### 4.4 一个 Gamma 型分布族

假设

\[
f_{a,b}(x)=H(a,b)x^ae^{-bx^c}\mathbf{1}\{x>0\},
\qquad
a>-1,
\qquad
b>0,
\]

其中 $c$ 已知。对于独立同分布样本，

\[
\prod_{j=1}^nf_{a,b}(x_j)
=H(a,b)^n
\exp\left\{a\sum_j\log x_j-b\sum_jx_j^c\right\}
\prod_j\mathbf{1}\{x_j>0\}.
\]

因此

\[
\left(\sum_j\log X_j,\sum_jX_j^c\right)
\]

对 $(a,b)$ 充分。

---

## 5. 指数族

!!! info "定义 5.1（指数族）"

    如果一个被控制的统计模型可以写成

    \[
    p_\vartheta(x)
    =\exp\{\eta(\vartheta)^TS(x)-B(\vartheta)\}h(x),
    \]

    则称它为一个具有 $k$ 维统计量的**指数族（exponential family）**。其中 $S=(S_1,\ldots,S_k)^T$ 是自然统计量，$\eta(\vartheta)$ 是自然参数。由定理 1，$S$ 是充分统计量。

原参数 $\vartheta$ 与自然参数 $\eta$ 不一定相同。这个区别在正态分布族中尤其重要，因为

\[
\eta_1=\frac{\mu}{\sigma^2},
\qquad
\eta_2=-\frac{1}{2\sigma^2}.
\]

!!! info "定义 5.2（规范形式与自然参数空间）"

    指数族的**规范形式（canonical form）** 为

    \[
    p_\eta(x)=\exp\{\eta^TS(x)-A(\eta)\}h(x),
    \qquad
    \eta\in\mathcal{H},
    \]

    其中自然参数空间为

    \[
    \mathcal{H}
    =\left\{\eta\in\mathbb{R}^k:
    0<\int e^{\eta^TS(x)}h(x)\,d\mu(x)<\infty\right\},
    \]

    对数规范化函数为

    \[
    A(\eta)
    =\log\int e^{\eta^TS(x)}h(x)\,d\mu(x).
    \]

对于独立同分布观测，规范形式变为

\[
p_\eta(x_1,\ldots,x_n)
=\exp\left\{\eta^T\sum_{i=1}^nS(x_i)-nA(\eta)\right\}
\prod_{i=1}^nh(x_i).
\]

因此，充分统计量 $\sum_iS(X_i)$ 的维数不会随着样本量 $n$ 增长。

!!! success "命题（对数规范化函数的几何与矩）"

    自然参数空间 $\mathcal{H}$ 是凸集，函数 $A$ 是凸函数。在可以交换微分与积分次序的内点处，

    \[
    \nabla A(\eta)=E_\eta S(X),
    \qquad
    \nabla^2A(\eta)=\operatorname{Cov}_\eta\{S(X)\}.
    \]

    因而，在 $a^TS(X)$ 非常数的方向 $a$ 上，$A$ 严格凸。

??? proof "命题的证明（点击展开）"

    对任意 $\eta,\zeta\in\mathcal{H}$ 和 $0<t<1$，由 Hölder 不等式，

    \[
    \int e^{\{(1-t)\eta+t\zeta\}^TS}h\,d\mu
    \leq
    \left(\int e^{\eta^TS}h\,d\mu\right)^{1-t}
    \left(\int e^{\zeta^TS}h\,d\mu\right)^t.
    \]

    所以 $(1-t)\eta+t\zeta\in\mathcal{H}$。两边取对数可得

    \[
    A\{(1-t)\eta+t\zeta\}
    \leq(1-t)A(\eta)+tA(\zeta),
    \]

    从而 $\mathcal{H}$ 和 $A$ 分别具有所述凸性。

    对规范化恒等式

    \[
    e^{A(\eta)}=\int e^{\eta^TS}h\,d\mu
    \]

    求导，得到

    \[
    \partial_jA(\eta)
    =e^{-A(\eta)}\int S_j(x)e^{\eta^TS(x)}h(x)\,d\mu(x)
    =E_\eta S_j(X).
    \]

    再求一次导数，得到

    \[
    \partial_{j\ell}^2A
    =E_\eta(S_jS_\ell)-E_\eta S_jE_\eta S_\ell,
    \]

    从而证明矩阵恒等式。最后，

    \[
    a^T\nabla^2A(\eta)a
    =\operatorname{Var}_\eta\{a^TS(X)\},
    \]

    它恰好在 $a^TS(X)$ 非常数时为正。$\square$

| 模型 | 自然统计量 $S(x)$ | 自然参数 $\eta$ | 自然参数空间 |
| --- | --- | --- | --- |
| $\operatorname{Bernoulli}(p)$ | $x$ | $\log\{p/(1-p)\}$ | $\mathbb{R}$ |
| $\operatorname{Poisson}(\lambda)$ | $x$ | $\log\lambda$ | $\mathbb{R}$ |
| $\operatorname{Exponential}(\lambda)$ | $x$ | $-\lambda$ | $(-\infty,0)$ |
| $N(\mu,\sigma^2)$ | $(x,x^2)$ | $(\mu/\sigma^2,-1/(2\sigma^2))$ | $\mathbb{R}\times(-\infty,0)$ |
| $\operatorname{Gamma}(\alpha,\beta)$ | $(\log x,x)$ | $(\alpha-1,-\beta)$ | $(-1,\infty)\times(-\infty,0)$ |

---

## 6. 最小充分性

完整样本本身总是充分的。最小充分性所寻找的是最粗糙的充分数据摘要，其中彼此一一对应的变换被视为等价。

!!! info "定义 6.1（最小充分统计量）"

    如果充分统计量 $T$ 满足：对每个充分统计量 $U$，都存在可测函数 $H$，使得

    \[
    T=H(U)
    \qquad P_\theta\text{-a.s. 对每个 }\theta,
    \]

    则称 $T$ 是**最小充分的（minimal sufficient）**。

因此，任意两个最小充分统计量互为彼此的函数；除去整个模型共同的零测集后，它们在样本空间上编码相同的划分。

### 6.1 似然成比例判据

本小节假设所有密度具有共同支撑 $\mathcal{S}$。

!!! success "定理 2（似然比判据）"

    假设 $T$ 是充分统计量，并且对任意 $x,y\in\mathcal{S}$，满足

    \[
    T(x)=T(y)
    \quad\Longleftrightarrow\quad
    \frac{p_\theta(x)}{p_\theta(y)}
    \text{ 与 }\theta\text{ 无关}.
    \]

    则 $T$ 是最小充分统计量。

??? proof "定理 2 的证明（点击展开）"

    设 $U$ 是任意充分统计量。由分解定理，可以写成

    \[
    p_\theta(x)=a_\theta\{U(x)\}b(x).
    \]

    若 $U(x)=U(y)$，则

    \[
    \frac{p_\theta(x)}{p_\theta(y)}
    =\frac{b(x)}{b(y)},
    \]

    该比值与 $\theta$ 无关。由判据中的反向蕴含可得 $T(x)=T(y)$。因此，$T$ 在 $U$ 的每个纤维上都是常数。由可测分解引理，存在可测函数 $H$ 使得 $T=H(U)$。因为 $U$ 是任意充分统计量，所以 $T$ 最小充分。$\square$

!!! note "为什么要单独假设充分性"

    在本课程使用的标准 Euclidean 模型中，教材常给出一个略强的版本：在常规可测性条件下，仅凭似然比条件也可以通过从每个纤维中选取一个似然代表来证明充分性。这里，所有应用中的充分性已经由分解定理建立。显式写出充分性条件，可以避免把一个可测选择步骤隐藏在似然比论证之中。

!!! note "该判据表达了什么"

    两个样本属于同一个最小充分性等价类，当且仅当它们的似然函数作为参数函数具有相同形状；二者最多相差一个只依赖数据的乘法常数。

!!! example "例 2（正态分布族）"

    对均值和方差均未知的正态模型，

    \[
    \log\frac{p_{\mu,\sigma^2}(x)}{p_{\mu,\sigma^2}(y)}
    =-\frac{\sum_i x_i^2-\sum_i y_i^2}{2\sigma^2}
    +\frac{\mu}{\sigma^2}\left(\sum_i x_i-\sum_i y_i\right).
    \]

    上式与 $(\mu,\sigma^2)$ 无关，当且仅当式中两个差都为零。因此，

    \[
    \left(\sum_iX_i,\sum_iX_i^2\right)
    \]

    是最小充分统计量。

!!! example "例 3（均匀分布的端点）"

    尽管该模型不具有共同支撑，相同的似然轮廓思想仍然有效。在正象限上，来自第 4.3 节的两个似然函数对所有 $\theta$ 成比例，当且仅当两个样本的最大值相同。如果最大值不同，可以选择严格位于二者之间的 $\theta$，使一个似然为零而另一个为正。因此，$X_{(n)}$ 是最小充分统计量。

    该结论使用了允许似然取零值的似然比判据之更一般的交叉乘积版本。

---

## 7. 为什么满仿射秩推出最小充分性

我们通过把模型约化到一个有限子族来证明指数族的结论。下面两个引理将论证中的关键步骤分离出来。

!!! info "引理 3（有限族的似然比统计量）"

    设 $P_0,\ldots,P_k$ 具有共同支撑，其密度分别为 $p_0,\ldots,p_k$。则

    \[
    R(X)=\left(\frac{p_1(X)}{p_0(X)},\ldots,
    \frac{p_k(X)}{p_0(X)}\right)
    \]

    是该有限分布族的最小充分统计量。

??? proof "引理 3 的证明（点击展开）"

    分解

    \[
    p_j(x)=R_j(x)p_0(x),
    \qquad
    j=1,\ldots,k,
    \qquad
    R_0\equiv1
    \]

    证明了 $R$ 的充分性。

    若 $U$ 是任意充分统计量，写成

    \[
    p_j(x)=a_j\{U(x)\}b(x).
    \]

    在共同支撑上，

    \[
    \frac{p_j(x)}{p_0(x)}
    =\frac{a_j\{U(x)\}}{a_0\{U(x)\}},
    \]

    因而 $R$ 的每个分量都是 $U$ 的函数。因此，对每个充分统计量 $U$，都有 $R=H(U)$，从而 $R$ 最小充分。$\square$

!!! info "引理 4（子族上的最小充分性）"

    假设 $\mathcal{P}$ 中所有分布具有共同支撑。如果 $T$ 对 $\mathcal{P}$ 充分，并且对某个子族 $\mathcal{P}_0\subset\mathcal{P}$ 最小充分，则 $T$ 对 $\mathcal{P}$ 最小充分。

??? proof "引理 4 的证明（点击展开）"

    设 $U$ 对 $\mathcal{P}$ 充分，则它对 $\mathcal{P}_0$ 也充分。$T$ 在 $\mathcal{P}_0$ 上的最小充分性给出

    \[
    T=H(U)
    \]

    在共同支撑上几乎处处成立。因为 $\mathcal{P}$ 中的每个分布都具有该支撑，所以同一个等式在每个 $P\in\mathcal{P}$ 下几乎处处成立。因此，$T$ 是完整分布族的每个充分统计量的函数。$\square$

!!! info "定义 7.1（满仿射秩与曲线族）"

    对规范指数族，设参数集合 $H\subseteq\mathcal{H}$。若 $H$ 的仿射包为 $\mathbb{R}^k$，等价地，若 $H$ 包含 $k+1$ 个仿射无关点，则称 $H$ 具有**满仿射秩（full affine rank）**。

    若 $H$ 包含 $\mathbb{R}^k$ 中的非空开集，则称该指数族是**满秩的（full-rank）**。满秩是满仿射秩的充分条件，但不是必要条件。如果一个分布族的展示坐标中没有冗余，而参数集合的内部为空，则通常称其为**曲线族（curved family）**。

!!! success "定理 5（指数族中的最小充分性）"

    考虑指数族

    \[
    p_\eta(x)=\exp\{\eta^TS(x)-A(\eta)\}h(x),
    \qquad
    \eta\in H\subseteq\mathcal{H},
    \]

    并假设其共同支撑为 $\{h>0\}$。如果 $H$ 包含 $k+1$ 个仿射无关点，则 $S$ 是最小充分统计量。特别地，该结论对每个满秩指数族都成立。

??? proof "定理 5 的证明（点击展开）"

    分解定理首先表明 $S$ 对整个分布族充分。选取仿射无关点

    \[
    \eta^{(0)},\eta^{(1)},\ldots,\eta^{(k)},
    \]

    并考虑与之对应的有限子族。对 $j=1,\ldots,k$，有

    \[
    \log\frac{p_{\eta^{(j)}}(x)}{p_{\eta^{(0)}}(x)}
    =\{\eta^{(j)}-\eta^{(0)}\}^TS(x)
    -\{A(\eta^{(j)})-A(\eta^{(0)})\}.
    \]

    令 $D$ 是一个 $k\times k$ 矩阵，其第 $j$ 行是 $\{\eta^{(j)}-\eta^{(0)}\}^T$。仿射无关性意味着 $D$ 非奇异。若 $r(x)$ 是上述 $k$ 个对数似然比组成的向量，$a$ 是对数规范化函数差组成的向量，则

    \[
    r(x)=DS(x)-a,
    \qquad
    S(x)=D^{-1}\{r(x)+a\}.
    \]

    因此，$S$ 与引理 3 中的似然比统计量互为一一函数。于是 $S$ 对选出的有限子族最小充分，再由引理 4，其最小充分性可以提升到整个分布族。$\square$

!!! success "推论 6（独立同分布观测）"

    对来自上述指数族的独立同分布样本，如果 $H$ 具有满仿射秩，则

    \[
    \sum_{i=1}^nS(X_i)
    \]

    是最小充分统计量。

??? proof "推论 6 的证明（点击展开）"

    独立同分布样本的联合密度仍是规范指数族，其自然统计量为 $\sum_iS(X_i)$，自然参数集合仍为 $H$。直接应用定理 5 即可。$\square$

!!! note "为什么开集条件强于证明所需条件"

    $\mathbb{R}^k$ 的开子集一定包含 $k+1$ 个仿射无关点，因此满秩条件表述起来很方便。但证明真正使用的只是参数集合具有满仿射包。有些曲线参数集合也具有满仿射包，因此仍然产生相同的最小充分统计量。

---

## 8. 满秩、曲线与冗余的正态子模型

下面三个例子区分了几个经常被混淆的概念。

### 情形 1：完整正态分布族

在完整分布族 $N(\mu,\sigma^2)$ 中，自然参数空间为

\[
\mathbb{R}\times(-\infty,0),
\]

它是开集。因此

\[
\left(\sum_iX_i,\sum_iX_i^2\right)
\]

是最小充分统计量。

### 情形 2：曲线子模型 $N(\sigma,\sigma^2)$

在 $\sigma>0$ 的子模型 $N(\sigma,\sigma^2)$ 中，自然参数满足

\[
\eta_1=\frac{1}{\sigma},
\qquad
\eta_2=-\frac{1}{2\sigma^2}
=-\frac{1}{2}\eta_1^2.
\]

该参数集合是一条曲线，所以模型不是满秩的。但是，这条曲线包含三个不共线的点，因此它在 $\mathbb{R}^2$ 中具有满仿射包；定理 5 仍然适用。

### 情形 3：冗余表示 $N(\sigma^2,\sigma^2)$

在子模型 $N(\sigma^2,\sigma^2)$ 中，密度中 $x$ 的系数恒为 $1$。该项可以被吸收到基准密度中，只留下 $x^2$ 作为依赖参数的统计量。因此，二分量表示是冗余的，而

\[
\sum_iX_i^2
\]

是充分且最小充分的。

这个区别非常重要：**原参数的维数、展示出来的自然统计量个数以及自然参数集合的仿射维数不一定相同。**

---

## 9. 充分性能够保证什么，又不能保证什么

**充分性不会自动推出最小充分性。** 如果 $T$ 是充分统计量，那么对任意统计量 $W$，$(T,W)$ 也充分。最小充分性负责去除这种不必要的附加信息。

**完备性是另一个独立概念。** 完备性研究的是：如果对每个 $\theta$ 都有 $E_\theta a(T)=0$，是否必然推出 $a(T)=0$ 几乎处处成立。这是一种不同的性质，将在第三章展开。满秩指数族经常同时具有最小充分性与完备性，但二者的证明和用途并不相同。

**充分性依赖于模型。** 在独立同分布 Bernoulli 模型中，$\sum_iX_i$ 对 Bernoulli 参数充分。但如果引入相关性、测量误差或额外的干扰参数，它可能不再充分。

**决策论预告。** 定义 1.1 中的共同核使我们能够在给定 $T$ 后，对 $X$ 中与参数无关的随机变化取条件平均。第三章将把这一观察转化为 Rao–Blackwell 改进，并在加入完备性后得到 Lehmann–Scheffé 定理。相关证明留到下一章，以避免重复。

---

## 10. 可数控制引理的详细证明

只需证明第 2 节 Halmos–Savage 引理中的命题 1 推出命题 3。

从一个 $\sigma$-有限控制测度 $\mu$ 出发。舍弃零测集后，选择一个可测划分

\[
\mathcal{X}=\bigcup_{m\geq1}A_m,
\qquad
0<\mu(A_m)<\infty,
\]

并定义

\[
\nu(B)=\sum_{m=1}^{\infty}2^{-m}
\frac{\mu(B\cap A_m)}{\mu(A_m)}.
\]

则 $\nu$ 是概率测度，并且与 $\mu$ 具有完全相同的零测集。因此，每个 $P_\theta\ll\nu$。

令 $\mathcal{G}$ 为所有可数混合分布组成的集合：

\[
Q=\sum_{j=1}^{\infty}c_jP_{\theta_j},
\qquad
c_j>0,
\qquad
\sum_jc_j=1.
\]

对 $Q\in\mathcal{G}$，记

\[
q=\frac{dQ}{d\nu},
\qquad
S_Q=\{x:q(x)>0\},
\]

并令

\[
a=\sup_{Q\in\mathcal{G}}\nu(S_Q)\leq1.
\]

选取 $Q_m\in\mathcal{G}$，使得

\[
\nu(S_{Q_m})>a-\frac{1}{m},
\]

并定义

\[
\lambda=\sum_{m=1}^{\infty}2^{-m}Q_m.
\]

把二重级数展开为单个级数可知，$\lambda$ 本身也是 $\mathcal{P}$ 中分布的一个可数混合。若 $q_m=dQ_m/d\nu$，则

\[
\frac{d\lambda}{d\nu}=\sum_m2^{-m}q_m,
\qquad
S_\lambda=\bigcup_mS_{Q_m}
\quad \nu\text{-a.e.}
\]

因为 $\lambda\in\mathcal{G}$，由 $a$ 的定义可知 $\nu(S_\lambda)\leq a$；而 $Q_m$ 的选择给出相反方向的不等式。因此

\[
\nu(S_\lambda)=a.
\]

现在固定任意 $P\in\mathcal{P}$，令 $p=dP/d\nu$，并假设集合

\[
B=\{p>0\}\setminus S_\lambda
\]

具有正的 $\nu$-测度。于是 $P(B)>0$，并且混合分布

\[
Q=\frac{1}{2}(\lambda+P)\in\mathcal{G}
\]

的支撑为 $S_\lambda\cup\{p>0\}$，其 $\nu$-测度严格大于 $a$，这与 $a$ 的定义矛盾。因此

\[
p=0
\qquad \nu\text{-a.e. on }S_\lambda^c.
\]

若 $\lambda(C)=0$，由于 $d\lambda/d\nu$ 在 $S_\lambda$ 上为正，必有

\[
\nu(C\cap S_\lambda)=0.
\]

结合 $p$ 在 $S_\lambda^c$ 上几乎处处为零，得到

\[
P(C)
=\int_{C\cap S_\lambda}p\,d\nu
+\int_{C\cap S_\lambda^c}p\,d\nu
=0.
\]

所以 $P\ll\lambda$。由于 $P$ 是任意的，$\lambda$ 控制整个分布族，从而证明命题 3。最后，列出展开混合分布 $\lambda$ 时出现的模型分布，便得到命题 2 中的可数子族。$\square$
