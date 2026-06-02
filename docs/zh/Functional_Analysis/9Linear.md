# 第九章：线性算子

本章我们将深入探讨泛函分析中的核心对象——线性算子。这部分内容将涵盖有界线性算子的基本概念、性质，以及线性算子空间的完备性结构。

---

## 1. 有界线性算子的基本概念与性质

### 1.1 线性算子的定义

!!! info "定义 9.1 (线性算子)"

    设 $E$ 及 $E_1$ 都是实（或复）的线性空间，$T$ 是由 $E$ 的某个子空间 $D$ 到线性空间 $E_1$ 中的映射。

    如果对任意的 $x, y \in D$，有：

    \[
    T(x+y) = T(x) + T(y)
    \]

    则称 $T$ 是**可加的**。

    若对任意的实（或复）数 $a$ 及任意的 $x \in D$，有：

    \[
    T(ax) = aT(x)
    \]

    则称 $T$ 是**齐次的**。

    既满足可加性又满足齐次性的映射称为**线性映射**或**线性算子**。

* $D$ 称为 $T$ 的定义域，通常记为 $D_T$。

* $T(D)$ 称为 $T$ 的值域，通常记为 $T(D_T)$。

* $D$ 中使 $Tx = \theta$ 的元素的集合称为 $T$ 的**零空间**，常用 $N_T$ 或类似记号表示。

* 设 $E_1$ 是实（或复）数域，于是 $T$ 成为由 $D$ 到实（或复）数域的算子。这时称 $T$ 为**泛函**。如果 $T$ 还是线性的，则称 $T$ 为**线性泛函**。

---

### 1.2 连续与有界线性算子

!!! info "定义 9.2 (连续与有界线性算子)"

    设 $E$ 及 $E_1$ 都是实（或复）的赋范线性空间，$T$ 是由 $E$ 的某个子空间 $D$ 到线性空间 $E_1$ 中的线性算子。

    * **连续线性算子**：如果 $T$ 在定义域上是连续的，则称 $T$ 为连续线性算子。

    * **有界线性算子**：如果 $T$ 将 $D$ 中任一有界集映成 $E_1$ 中的有界集，则称 $T$ 是有界线性算子。

    * **无界线性算子**：如果存在 $D$ 中的有界集 $A$ 使得 $T(A)$ 是 $E_1$ 中的无界集，则称 $T$ 是无界线性算子。

* 将赋范线性空间 $E$ 中的每个元素 $x$ 映成 $x$ 自身的算子称为 $E$ 上的**单位算子**，常以 $I$ 表示。

* 将 $E$ 中的每个元素 $x$ 映成零元素 $\theta$ 的算子称为**零算子**。

* 单位算子与零算子既是有界线性算子也是连续线性算子。

---

### 1.3 无界算子的典型示例

在连续函数空间 $C[a,b]$ 中讨论微分算子 $T = \frac{d}{dt}$。

将 $[a,b]$ 上连续可微函数构成的集 $C^1[a,b]$ 作为 $T$ 的定义域 $D$。则 $T$ 是定义在 $C^1[a,b]$ 上，并且在 $C[a,b]$ 中取值的线性算子。

我们可以证明 $T$ 是一个**无界算子**：

取函数序列：

\[
x_n(t) = \sin(nt)
\]

则 $x_n \in C^1[0,1]$ 且 $\|x_n\| = 1$，但是：

\[
\|Tx_n\| = \left\| \frac{d}{dt} \sin(nt) \right\| = n\|\cos(nt)\| = n
\]

故 $T$ 将 $C[0,1]$ 中的单位球面映成 $C[0,1]$ 中的无界集。因此 $T$ 是无界的。

---

### 1.4 线性算子有界性与连续性的等价关系

首先看一个连续线性泛函的例子。连续函数的积分：

\[
f(x) = \int_a^b x(t) dt
\]

是定义在连续函数空间 $C[a,b]$ 上的一个有界线性泛函，也是连续线性泛函。因为：

\[
|f(x)| \le (b-a)\|x\|
\]

!!! success "定理 9.1"

    设 $E$ 及 $E_1$ 都是实赋范线性空间，$T$ 是由 $E$ 的子空间 $D$ 到 $E_1$ 的连续可加算子，则 $T$ 满足齐次性。因此 $T$ 是连续线性算子。

??? proof "定理 9.1 的证明（点击展开）"

    对任给的 $x \in D$，令 $f(\alpha) = T(\alpha x)$，则 $f$ 连续。
    
    若 $\alpha_n \rightarrow \alpha_0$，则 $\alpha_n x \rightarrow \alpha_0 x$。所以 $T(\alpha_n x) \rightarrow T(\alpha_0 x)$，即 $f(\alpha_n) \rightarrow f(\alpha_0)$。

    对任意两个实数 $\alpha, \beta$：

    \[
    f(\alpha+\beta) = T((\alpha+\beta)x) = T(\alpha x + \beta x) = T(\alpha x) + T(\beta x) = f(\alpha) + f(\beta)
    \]

    由实变函数引理可知，对任何实数有 $f(\alpha) = \alpha f(1)$，即 $T(\alpha x) = \alpha T(x)$。证明完毕。 $\square$

* **推论**：设 $E$ 及 $E_1$ 都是复赋范线性空间，$T$ 是由 $E$ 的子空间 $D$ 到 $E_1$ 的连续可加算子，且满足 $T(ix) = iT(x)$，则 $T$ 满足齐次性。

    证明如下：令 $\alpha = a + ib$，
    
    \[
    T(\alpha x) = T(ax+ibx) = T(ax) + T(ibx) = aT(x) + bT(ix) = aT(x) + biT(x) = \alpha T(x)
    \]

!!! success "定理 9.2 (线性算子有界的等价定义)"

    设 $E$ 及 $E_1$ 都是赋范线性空间，$T$ 是由 $E$ 的子空间 $D$ 到 $E_1$ 的线性算子。则 $T$ 有界的充分必要条件是存在 $M > 0$ 使得对一切 $x \in D$，有：

    \[
    \|Tx\| \le M\|x\|
    \]

??? proof "定理 9.2 的证明（点击展开）"

    **充分性**：设 $A \subset D$ 为任一有界集。则存在数 $K$，使得对一切 $x \in A$ 有 $\|x\| \le K$。故当 $x \in A$ 时，

    \[
    \|Tx\| \le M\|x\| \le MK
    \]

    因此 $T(A)$ 是 $E_1$ 中的有界集，即 $T$ 有界。

    **必要性**：在 $D$ 中取单位球面 $S := \{x : \|x\| = 1, x \in D\}$。
    因 $S$ 有界，故 $T(S)$ 也有界。于是存在正数 $M > 0$ 使得当 $x \in S$ 时，$\|Tx\| \le M$。
    
    设 $x$ 是 $D$ 中任一非零元，则 $\frac{x}{\|x\|} \in S$。所以：

    \[
    \left\| T\left(\frac{x}{\|x\|}\right) \right\| \le M
    \]

    由 $T$ 的齐次性可知，

    \[
    \|Tx\| \le M\|x\|
    \]

    当 $x = \theta$ 时，不等式显然成立。证明完毕。 $\square$

!!! success "定理 9.3 (有界性与连续性等价)"

    设 $E$ 及 $E_1$ 都是赋范线性空间，$T$ 是由 $E$ 的子空间 $D$ 到 $E_1$ 的线性算子。则下列性质等价：

    (i) $T$ 连续；

    (ii) $T$ 在原点 $\theta$ 处连续；

    (iii) $T$ 有界。

??? proof "定理 9.3 的证明（点击展开）"

    **(i) $\Rightarrow$ (ii)**：显然成立。

    **(ii) $\Rightarrow$ (iii)**：因 $T$ 在原点连续，故对 $\epsilon=1$，存在 $\delta>0$，使得对任给的 $x \in D$ 当 $\|x\| \le \delta$ 时，有：

    \[
    \|Tx\| = \|Tx - T\theta\| \le 1
    \]

    设 $x$ 是 $D$ 中任一非零元，则 $\left\| \frac{\delta x}{\|x\|} \right\| = \delta \le \delta$。所以：

    \[
    \left\| T\left(\frac{\delta x}{\|x\|}\right) \right\| \le 1
    \]

    由齐次性提出系数，即得：

    \[
    \|Tx\| \le \frac{1}{\delta}\|x\|
    \]

    说明算子是有界的。

    **(iii) $\Rightarrow$ (i)**：设 $x_n, x \in D$，且 $x_n \rightarrow x$。因 $T$ 有界，故：

    \[
    \|Tx_n - Tx\| = \|T(x_n - x)\| \le M\|x_n - x\| \rightarrow 0
    \]

    即 $Tx_n \rightarrow Tx$，这证明了连续性。

    对线性算子来说，有界性、连续性以及在原点的连续性均相互等价。这些条件也与 $T$ 在 $D$ 中任一给定点 $x_0$ 处的连续性等价。 $\square$

---

### 1.5 有界线性算子的范数

!!! info "定义 9.3 (算子范数)"

    设 $E$ 及 $E_1$ 都是赋范线性空间，$T$ 是由 $E$ 的子空间 $D$ 到 $E_1$ 的有界线性算子。使得 $\|Tx\| \le M\|x\|$ 对一切 $x \in D$ 都成立的正数 $M$ 的下确界称为 $T$ 的**范数**，记为 $\|T\|$。

因 $M$ 是集合 $\{\frac{\|Tx\|}{\|x\|} : x \in D, x \ne \theta\}$ 的一个上界，因此算子 $T$ 的范数 $\|T\|$ 作为所有上界 $M$ 的下确界，也是上述集合的最小上界（即上确界）。于是有：

\[
\|T\| = \sup_{x \in D, x \ne \theta} \frac{\|Tx\|}{\|x\|}
\]

由此容易导出下列结论：

* 对一切 $x \in D$，有 $\|Tx\| \le \|T\|\|x\|$。

* 范数的其他等价计算形式：

\[
\|T\| = \sup_{x \in D, x \ne \theta} \left\| T\left(\frac{x}{\|x\|}\right) \right\| = \sup_{x \in D, \|x\|=1} \|Tx\| = \sup_{x \in D, \|x\| \le 1} \|Tx\|
\]

算子的范数是刻画有界线性算子的一个重要的量。一般情形下求有界线性算子的范数是很困难的，下面我们通过几个例子说明如何估计或精确求出其范数。

---

### 1.6 算子范数的计算示例

**示例 1：矩阵表示的线性算子**

设 $(a_{ij})$ 为一给定的 $n \times n$ 实方阵，由等式：

\[
\eta_i = \sum_{j=1}^n a_{ij} \xi_j
\]

定义了一个由 $\mathbb{R}^n$ 到 $\mathbb{R}^n$ 的算子 $T : Tx = y$。它将元素 $x = (\xi_1, \dots, \xi_n)$ 映成元素 $y = (\eta_1, \dots, \eta_n)$。

在 $\mathbb{R}^n$ 中取另一个向量 $x' = (\xi_1', \dots, \xi_n')$，易证 $T(x+x') = Tx + Tx'$ 及 $T(ax) = aTx$，因此 $T$ 是线性算子。

我们估计其有界性。由 Cauchy 不等式，有：

\[
\|Tx\|^2 = \sum_{i=1}^n \eta_i^2 = \sum_{i=1}^n \left(\sum_{j=1}^n a_{ij} \xi_j\right)^2 \le \sum_{i=1}^n \left( \sum_{j=1}^n a_{ij}^2 \sum_{j=1}^n \xi_j^2 \right) \le \sum_{i=1}^n \sum_{j=1}^n a_{ij}^2 \cdot \|x\|^2
\]

由此可得：

\[
\|T\| \le \left( \sum_{i,j=1}^n a_{ij}^2 \right)^{\frac{1}{2}}
\]

因此算子连续且有界。

**示例 2：Fourier 形式的积分算子**

$C(\mathbb{R})$：定义在 $\mathbb{R}$ 上有界连续函数构成的集，定义其范数为 $\|y\| = \sup_{t \in \mathbb{R}} |y(t)|$。则 $C(\mathbb{R})$ 是一个 Banach 空间。

设 $x \in L(\mathbb{R})$（即可积函数空间）。令 $y = Tx$，且

\[
y(t) = \int_{-\infty}^\infty e^{-ist} x(s) ds
\]

$T$ 是定义在 $L(\mathbb{R})$ 上而值域包含在 $C(\mathbb{R})$ 中的线性算子。由范数定义：

\[
|Tx(t)| = |y(t)| \le \int_{-\infty}^\infty |e^{-ist} x(s)| ds = \int_{-\infty}^\infty |x(s)| ds = \|x\|
\]

所以 $\|T\| \le 1$。$T$ 有界，因此连续。

**示例 3：Lagrange 插值多项式**

设 $x \in C[a,b]$，在 $[a,b]$ 中任取 $n$ 个点 $a \le t_1 < t_2 < \dots < t_n \le b$，作多项式：

\[
l_k(t) = \frac{(t-t_1)\dots(t-t_{k-1})(t-t_{k+1})\dots(t-t_n)}{(t_k-t_1)\dots(t_k-t_{k-1})(t_k-t_{k+1})\dots(t_k-t_n)}
\]

定义算子 $y = L_n x = \sum_{k=1}^n x(t_k) l_k(t)$。

则有：

\[
\|L_n x\| = \max_{t \in [a,b]} \left| \sum_{k=1}^n x(t_k) l_k(t) \right| \le \max_{t \in [a,b]} \sum_{k=1}^n |l_k(t)| \max_{t \in [a,b]} |x(t)| = \left( \max_{t \in [a,b]} \sum_{k=1}^n |l_k(t)| \right) \|x\|
\]

令 $\gamma := \max_{t \in [a,b]} \sum_{k=1}^n |l_k(t)|$，则有 $\|L_n\| \le \gamma$。

另一方面，由于 $\sum_{k=1}^n |l_k(t)|$ 在 $[a,b]$ 上连续，故存在 $t_0 \in [a,b]$ 使得：

\[
\gamma = \sum_{k=1}^n |l_k(t_0)|
\]

取 $x_0 \in C[a,b]$ 满足 $\|x_0\| = 1$ 且在各节点处 $x_0(t_k) = \text{sgn}(l_k(t_0))$，于是：

\[
\|L_n\| \ge \|L_n(x_0)\| \ge |L_n(x_0)(t_0)| = \left| \sum_{k=1}^n l_k(t_0) x_0(t_k) \right| = \left| \sum_{k=1}^n l_k(t_0) \text{sgn}(l_k(t_0)) \right| = \sum_{k=1}^n |l_k(t_0)| = \gamma
\]

从而证明了精确范数 $\|L_n\| = \gamma$。

**示例 4：含核积分算子**

设 $K(t,s)$ 是定义在 $a \le t \le b$, $a \le s \le b$ 上的连续函数。在空间 $C[a,b]$ 上定义如下的积分算子：

\[
y(t) = Tx(t) = \int_a^b K(t,s) x(s) ds
\]

则 $T$ 为 $C[a,b]$ 到其自身的有界线性算子，且范数满足：

\[
\|T\| = \gamma, \quad \gamma := \max_{t \in [a,b]} \int_a^b |K(t,s)| ds
\]

*上界的证明*：

\[
\|Tx\| = \max_{t \in [a,b]} \left| \int_a^b K(t,s) x(s) ds \right| \le \max_{s \in [a,b]} |x(s)| \max_{t \in [a,b]} \int_a^b |K(t,s)| ds = \gamma \|x\|
\]

即 $\|T\| \le \gamma$。

*下界的证明*：
由于 $\int_a^b |K(t,s)| ds$ 连续，存在 $t_0$ 使其达到最大值 $\gamma$。取 $\|\varphi\| = 1$，

\[
\|T\| \ge \|T\varphi\| = \max_{t \in [a,b]} \left| \int_a^b K(t,s) \varphi(s) ds \right| \ge \left| \int_a^b K(t_0,s) \varphi(s) ds \right|
\]

令 $E_0 = \{s : K(t_0, s) \ge 0\}$。作函数：

\[
\varphi_n(s) = \frac{1 - n d(s, E_0)}{1 + n d(s, E_0)}
\]

其中 $d(s, E_0)$ 为 $s$ 与 $E_0$ 的距离。则 $\varphi_n$ 于 $[a,b]$ 上连续，且 $|\varphi_n| \le 1$。
更进一步，$\varphi_n|_{E_0} \equiv 1$，且 $\varphi_n|_{E_0^c} \rightarrow -1$。

由 Lebesgue 控制收敛定理，当 $n \rightarrow \infty$ 时：

\[
T\varphi_n(t_0) = \int_a^b K(t_0, s) \varphi_n(s) ds \rightarrow \int_a^b |K(t_0, s)| ds = \gamma
\]

所以：

\[
\gamma = \lim_{n \rightarrow \infty} T\varphi_n(t_0) \le \lim_{n \rightarrow \infty} \|T\varphi_n\| \le \lim_{n \rightarrow \infty} \|T\|\|\varphi_n\| = \|T\|
\]

即得 $\|T\| \ge \gamma$。

---

## 2. 线性算子空间

### 2.1 $\mathbb{B}(E, E_1)$ 空间

我们将赋范线性空间 $E$ 到赋范线性空间 $E_1$ 的**有界线性算子全体**构成的集合记为 $\mathbb{B}(E, E_1)$。

!!! success "定理 9.4"

    设 $E$ 及 $E_1$ 都是赋范线性空间，在 $\mathbb{B}(E, E_1)$ 中定义线性运算如下：

    * 加法：$(T + T')x = Tx + T'x$
    * 数乘：$(\alpha T)x = \alpha Tx$

    其中 $T, T' \in \mathbb{B}(E, E_1)$，$\alpha$ 为数。那么 $\mathbb{B}(E, E_1)$ 按照上述线性运算是一个线性空间。若以算子的范数作为其范数，则 $\mathbb{B}(E, E_1)$ 是一个赋范线性空间。

??? proof "定理 9.4 的证明（点击展开）"

    首先，按定义可知它是线性空间。需要验证算子范数满足范数的三大公理：

    * **非负性**： $\|T\| = \sup_{\|x\|=1} \|Tx\| \ge 0$。若 $T = \theta$ (零算子)，显然有 $\|T\| = 0$。反之，若 $\|T\| = 0$，则对一切 $x \in E$ 有 $Tx = 0$，故 $T = \theta$。

    * **齐次性**：
    
    \[
    \|\alpha T\| = \sup_{\|x\|=1} \|\alpha Tx\| = |\alpha| \sup_{\|x\|=1} \|Tx\| = |\alpha|\|T\|
    \]

    * **三角不等式**：
    
    \[
    \|T + T'\| = \sup_{\|x\|=1} \|(T + T')x\| = \sup_{\|x\|=1} \|Tx + T'x\| \le \sup_{\|x\|=1} (\|Tx\| + \|T'x\|) \le \sup_{\|x\|=1} \|Tx\| + \sup_{\|x\|=1} \|T'x\| \le \|T\| + \|T'\|
    \]

    因此，这是一个赋范线性空间。 $\square$

*注：我们将 $\mathbb{B}(E, E)$ 简记为 $\mathbb{B}(E)$，代表 $E$ 到其自身的有界线性算子全体。*

---

### 2.2 按算子范数收敛 (一致算子拓扑)

设 $T, T_n \in \mathbb{B}(E, E_1)$ ($n=1,2,3,\dots$)。若 $T_n$ 按 $\mathbb{B}(E, E_1)$ 中的范数收敛于 $T$，即：

\[
\|T - T_n\| \rightarrow 0
\]

则称 $T_n$ **按算子范数收敛**于 $T$，或称 $T_n$ 按**一致算子拓扑**收敛于 $T$。
使用“按一致算子拓扑收敛”这一名称的原因在于下面的定理：

!!! success "定理 9.5"

    设 $T, T_n \in \mathbb{B}(E, E_1)$。$T_n$ 按算子范数收敛于 $T$ 的充分必要条件是：$\{T_n\}$ 在 $E$ 中的**任一有界集上一致收敛**于 $T$。

??? proof "定理 9.5 的证明（点击展开）"

    **必要性**：设 $A \subset E$ 为有界集。则存在正数 $K > 0$，使得当 $x \in A$ 时，$\|x\| \le K$。故：

    \[
    \|Tx - T_nx\| \le \|T - T_n\|\|x\| \le K\|T - T_n\|
    \]

    因为 $\|T - T_n\| \rightarrow 0$，则任给 $\epsilon > 0$，存在 $N > 0$，使得当 $n > N$ 时，$\|T - T_n\| < \epsilon/K$。代入上式得：

    \[
    \|Tx - T_nx\| \le \epsilon
    \]

    这对于所有的 $x \in A$ 是一致成立的。故 $\{T_n\}$ 在 $A$ 上一致收敛于 $T$。

    **充分性**：设 $\{T_n\}$ 在 $E$ 中的任一有界集上一致收敛于 $T$。取 $E$ 中的单位球面 $S = \{x : \|x\| = 1, x \in E\}$。
    根据假定，对任给的 $\epsilon > 0$，存在 $N > 0$，使得当 $n > N$ 时，不等式：

    \[
    \|Tx - T_nx\| \le \epsilon
    \]

    对于所有的 $x \in S$ 一致地成立。于是：

    \[
    \|T - T_n\| = \sup_{\|x\|=1} \|Tx - T_nx\| \le \epsilon
    \]

    故 $\{T_n\}$ 按算子范数收敛于 $T$。 $\square$

---

### 2.3 $\mathbb{B}(E, E_1)$ 空间的完备性

算子空间本身的完备性由目标空间的完备性决定。

!!! success "定理 9.6"

    设 $E_1$ 是 Banach 空间，则算子空间 $\mathbb{B}(E, E_1)$ 也是 Banach 空间。

??? proof "定理 9.6 的证明（点击展开）"

    设 $\{T_n\}$ 是 $\mathbb{B}(E, E_1)$ 中的 Cauchy 列。则对任给的 $\epsilon > 0$，存在 $N > 0$，使得当 $m, n > N$ 时：

    \[
    \|T_m - T_n\| < \epsilon
    \]

    任取 $x \in E$，考察序列 $\{T_n x\}$：

    \[
    \|T_mx - T_nx\| \le \|T_m - T_n\|\|x\| \le \epsilon\|x\|
    \]

    故对于固定的 $x$，序列 $\{T_n x\}$ 是 $E_1$ 中的 Cauchy 列。因 $E_1$ 是完备的，故 $\{T_n x\}$ 在 $E_1$ 中收敛于某一元素 $y$。于是可以定义点态极限算子 $T$：

    \[
    Tx = \lim_{n \rightarrow \infty} T_nx = y
    \]

    接下来我们需要证明 $T$ 属于 $\mathbb{B}(E, E_1)$，且 $T_n \to T$ 是依范数收敛的。

    **1. 线性性**：

    \[
    T(ax+by) = \lim_{n \rightarrow \infty} T_n(ax+by) = \lim_{n \rightarrow \infty} (aT_nx + bT_ny) = aTx + bTy
    \]

    **2. 有界性及算子范数收敛**：
    在前面的 Cauchy 列条件 $\|T_mx - T_nx\| \le \epsilon\|x\|$ 中，令 $m \rightarrow \infty$，由范数的连续性可得，对于所有的 $n > N$：

    \[
    \|Tx - T_nx\| \le \epsilon\|x\|
    \]

    这说明 $T - T_n$ 是有界算子，即 $T - T_n \in \mathbb{B}(E, E_1)$。既然 $T_n$ 是有界的，从而 $T = (T - T_n) + T_n$ 必然也是有界线性算子，故 $T \in \mathbb{B}(E, E_1)$。

    同时这也表明，当 $n > N$ 时：

    \[
    \|T - T_n\| \le \epsilon
    \]

    故 $\{T_n\}$ 按算子范数收敛于 $T$。由此可知，$\mathbb{B}(E, E_1)$ 中任一基本点列必有极限，因此它是完备的 Banach 空间。 $\square$

---

### 2.4 不按算子范数收敛的例子

虽然点态收敛可以定义一个极限算子，但点态收敛并不蕴含算子范数收敛。

在 $l^p$ 空间中定义**截断位移算子** $T_n$：

\[
T_n x = x_n
\]

其中输入向量为 $x = (\xi_1, \xi_2, \dots, \xi_n, \dots) \in l^p$，截断输出向量为 $x_n = (\xi_n, \xi_{n+1}, \dots) \in l^p$（这里符号根据讲义上下文，即截断前 $n-1$ 项）。

$T_n$ 是有界线性算子，且 $\|T_nx\|_p \le \|x\|_p$。
对任意给定的 $x \in l^p$，由于级数 $\sum |\xi_k|^p$ 收敛，其尾部必然趋于 0，有 $\|x_n\|_p \rightarrow 0$。因此在逐点意义下：

\[
\lim_{n \rightarrow \infty} T_n x = \theta
\]

（即算子列点态收敛于零算子）。

然而，考察算子范数 $\|T_n\|$。取特定单位向量：

\[
y_n = (0, \dots, 0, 1, 0, \dots) \quad \text{（第 } n \text{ 个位置为 1）}
\]

则 $T_n y_n = (1, 0, \dots)$。故：

\[
\|T_n\| \ge \|T_n y_n\|_p = 1
\]

所以 $\|T_n\| = 1$。这说明 $\|T_n - 0\| = 1 \not\rightarrow 0$。

因此，序列 $\{T_n\}$ **不按算子范数收敛**于零算子。