#  中心极限定理（二）：$m$-相依序列与稳定分布

本节将突破传统中心极限定理中“独立同分布 (i.i.d.)”的限制，首先探讨具有局部依赖性的 **$m$-相依序列 (m-dependent sequence)** 的中心极限定理；随后通过 **Cramér-Wold 定理** 将一维结论推广到多维；最后深入研究极限分布的广义形式——**稳定分布 (Stable Distributions)** 以及它们的**吸引域 (Domain of Attraction, DA)**。

---

## 1. $m$-相依随机变量序列的 CLT

在实际应用中（如时间序列分析），数据往往存在序列相关性。我们首先考虑一种最简单的相依结构：局部相依。

!!! info "定义 3.13：$m$-相依序列 ($m$-dependent sequence)"

    一个随机变量序列 $\{X_n\}_{n \ge 1}$ 被称为 **$m$-相依的 ($m$-dependent)**，如果存在一个正整数 $m$，使得对于任意的 $n \ge 1$ 和 $j \ge m$，随机变量 $X_{n+j}$ 独立于由前 $n$ 个变量生成的 $\sigma$-代数 $\mathcal{F}_n = \sigma\{X_i, 1 \le i \le n\}$。
    
    *例子：阶数为 $q$ 的滑动平均模型 MA($q$) 是一个 $q+1$ 相依序列。*

!!! success "定理 3.14 ($m$-相依序列的 CLT)"

    设 $\{X_n\}_{n \ge 1}$ 是一个 $m$-相依序列。假设这些随机变量是**一致有界的 (uniformly bounded)**（即存在常数 $M$ 使得 $\sup_n |X_n| \le M$）。
    记 $S_n = \sum_{i=1}^n X_i$，$\sigma_n^2 = Var(S_n)$。如果满足方差增长条件：
    
    $$
    \frac{\sigma_n}{m n^{1/3}} \rightarrow \infty \quad \text{当 } n \rightarrow \infty
    $$
    
    并且 $m = o(n^{1/3})$，那么有：
    
    $$
    \frac{S_n - E(S_n)}{\sigma_n} \xrightarrow{d} N(0,1)
    $$

    *(注：通过引入 Lindeberg 条件，可以去掉“一致有界”的假设，详见 S. Janson (2021))*。

??? proof "定理 3.14 的证明：分块技术 (Blocking Technique)（点击展开）"

    **核心思想：** 将整个序列划分为“大块”和“小块”交替的形式。由于 $m$-相依性，只要大块之间的间隔（小块的长度）大于等于 $m$，大块之间就是相互独立的。
    
    不失一般性，假设 $E(X_j) = 0$。由于序列一致有界，存在 $M$ 使得 $\sup_n |X_n| \le M$。
    
    **第一步：构造大块与小块**
    
    令大块的长度为 $k = [n^{1/3}]$，小块的长度为 $m$。
    则总的块数为 $p_n = [\frac{n}{k+m}] = O(n^{2/3})$。记 $B_j = j(k+m)$。
    
    我们构造：
    
    * **大块 (Large blocks):** $Y_j = X_{B_{j-1}+1} + \cdots + X_{B_{j-1}+k}$ （共 $p_n$ 块）

    * **小块 (Small blocks):** $Z_j = X_{B_{j-1}+k+1} + \cdots + X_{B_j}$ （共 $p_n$ 块）
    
    * **剩余块 (Residual block):** $R_p = X_{B_{p_n}+1} + \cdots + X_n$
    
    由于当 $n$ 足够大时 $k \gg m$，且大块 $Y_j$ 之间的间隔为 $m$，因此序列 $\{Y_j\}_{j=1}^{p_n}$ 相互独立。同理，$\{Z_j\}_{j=1}^{p_n}$ 也相互独立。
    
    我们将总和分解为三部分：
    
    $$
    S_n = \sum_{j=1}^{p_n} Y_j + \sum_{j=1}^{p_n} Z_j + R_p := S_n' + S_n'' + S_n'''
    $$
    
    **第二步：控制小块和剩余块的方差**
    
    由于 $\sup_j |X_j| \le M$，我们有协方差界 $|E(X_j X_l)| \le M^2$。
    对于剩余块 $S_n'''$：
    
    $$
    Var(S_n''') = E[(S_n''')^2] = \left| \sum_{j,l} E(X_{\dots} X_{\dots}) \right| \le (n - p_n(k+m))^2 M^2 \le (k+m)^2 M^2
    $$
    
    因此，依概率有界：
    
    $$
    S_n''' = O_p(\sqrt{Var(S_n''')}) = O_p(k+m) = O_p(n^{1/3})
    $$
    
    同理，对于小块之和 $S_n''$：
    
    $$
    E[Z_j^2] = E\left[ \left(\sum X_{\dots}\right)^2 \right] \le m^2 M^2
    $$
    
    由于 $Z_j$ 相互独立，有 $Var(S_n'') = \sum Var(Z_j) \le p_n m^2 M^2$。从而：
    
    $$
    S_n'' = O_p(p_n^{1/2} m) = O_p(n^{1/3} m)
    $$
    
    **第三步：证明小块和剩余块可以忽略**
    
    利用题设条件 $\sigma_n / (m n^{1/3}) \rightarrow \infty$，我们有：
    
    $$
    \frac{S_n''}{\sigma_n} = \frac{S_n''}{m n^{1/3}} \times \frac{m n^{1/3}}{\sigma_n} = O_p(1) \cdot o(1) = o_p(1)
    $$
    
    同理，因为 $k = O(n^{1/3})$，有 $S_n''' / \sigma_n = o_p(1)$。
    
    因此，总和的标准化形式可以写为：
    
    $$
    \frac{S_n}{\sigma_n} = \frac{S_n'}{\sigma_n} + o_p(1) = \frac{\sigma_n'}{\sigma_n} \frac{S_n'}{\sigma_n'} + o_p(1)
    $$
    
    其中 $\sigma_n'^2 = Var(S_n')$。接下来只需要证明 $\sigma_n'^2 / \sigma_n^2 \rightarrow 1$ 且 $S_n' / \sigma_n' \xrightarrow{d} N(0,1)$。
    
    **第四步：方差渐近等价**
    
    展开 $S_n$ 的方差：
    
    $$
    E(S_n^2) = E(S_n'^2) + E(S_n''^2) + E(S_n'''^2) + 2E(S_n' S_n'') + \dots
    $$
    
    其中交叉项由于 $m$-相依性，大部分协方差为 0：
    
    $$
    E(S_n' S_n'') = \sum_{j,l=1}^{p_n} Cov(Y_j, Z_l) = \sum_{j=1}^{p_n} [Cov(Y_j, Z_j) + Cov(Y_j, Z_{j-1})] \le 2p_n (mM)^2
    $$
    
    综合各项误差的阶数，可以得到：
    
    $$
    \left| 1 - \frac{\sigma_n'^2}{\sigma_n^2} \right| = O\left( \frac{m^2 n^{2/3}}{\sigma_n^2} \right) \rightarrow 0
    $$
    
    故 $\sigma_n'^2 / \sigma_n^2 \rightarrow 1$。
    
    **第五步：大块的中心极限定理**
    
    由于大块 $Y_j$ 之间是相互独立的，我们可以对 $\{Y_j\}$ 验证 Lindeberg 条件。
    由于 $|Y_j| \le kM = O(n^{1/3}) = o(\sigma_n')$，对于任意 $\eta > 0$，当 $n$ 足够大时，指示函数 $\mathbb{I}(|Y_j| \ge \eta \sigma_n')$ 将恒为 0：
    
    $$
    \frac{1}{\sigma_n'^2} \sum_{j=1}^{p_n} E\left[ Y_j^2 \mathbb{I}(|Y_j| \ge \eta \sigma_n') \right] \rightarrow 0
    $$
    
    因此 Lindeberg 条件成立。由 Lindeberg-Feller CLT，我们得到 $S_n' / \sigma_n' \xrightarrow{d} N(0,1)$。
    结合 Slutsky 定理，最终结论得证。 $\square$

---

## 2. 多维中心极限定理与 Cramér-Wold 定理

为了将一维的中心极限定理推广到多维随机向量，我们借助 Cramér-Wold 定理。它的核心思想是：**多维随机向量的弱收敛，等价于其在任意一维方向上投影的弱收敛。**

!!! info "定理 3.15：Cramér-Wold 定理"

    设 $X_n$ 是 $\mathbb{R}^d$ 中的随机向量序列，$X$ 是随机向量。则 $X_n$ 依分布收敛于 $X$ 当且仅当对于任意的线性组合方向 $a \in \mathbb{R}^d$，都有：
    
    $$
    X_n \xrightarrow{d} X \iff a^T X_n \xrightarrow{d} a^T X, \quad \forall a \in \mathbb{R}^d
    $$

??? proof "Cramér-Wold 定理的证明（点击展开）"

    "$\implies$"：由连续映射定理 (Continuous Mapping Theorem) 显然成立，因为内积函数 $g(x) = a^T x$ 是连续函数。
    
    "$\impliedby$"：利用特征函数 (Characteristic Function)。设 $X_n = (X_{n1}, \dots, X_{nd})^T$。
    任取 $c = (c_1, \dots, c_d)^T \in \mathbb{R}^d$。已知条件意味着：
    
    $$
    c^T X_n = c_1 X_{n1} + \dots + c_d X_{nd} \xrightarrow{d} c_1 X_1 + \dots + c_d X_d = c^T X
    $$
    
    根据 Lévy 连续性定理，一维随机变量的弱收敛意味着其特征函数逐点收敛。
    对于 $c^T X_n$，其在参数 $t$ 处的特征函数为：
    
    $$
    \phi_{c^T X_n}(t) = E\left[ e^{it(c_1 X_{n1} + \dots + c_d X_{nd})} \right]
    $$
    
    特别地，取 $t=1$，则有：
    
    $$
    \lim_{n \rightarrow \infty} E\left[ e^{i(c_1 X_{n1} + \dots + c_d X_{nd})} \right] = E\left[ e^{i(c_1 X_1 + \dots + c_d X_d)} \right]
    $$
    
    这恰好就是多维随机向量 $X_n$ 在 $c$ 处的联合特征函数 $\phi_n(c_1, \dots, c_d)$。即：
    
    $$
    \lim_{n \rightarrow \infty} \phi_n(c) = \phi_X(c), \quad \forall c \in \mathbb{R}^d
    $$
    
    联合特征函数处处收敛，由多维 Lévy 连续性定理，即证 $X_n \xrightarrow{d} X$。 $\square$

利用 Cramér-Wold 定理，多维 i.i.d. 序列的 CLT 就变得非常直接：

!!! success "定理 3.16 (多元中心极限定理 Multivariate CLT)"

    设 $X_1, X_2, \dots$ 是 i.i.d. 的 $d$ 维随机向量，具有均值向量 $\mu$ 和有限的协方差矩阵 $\Sigma$。记 $\overline{X}_n = \frac{1}{n} \sum_{i=1}^n X_i$，则：
    
    $$
    \sqrt{n}(\overline{X}_n - \mu) \xrightarrow{d} N_d(0, \Sigma)
    $$

---

## 3. 稳定分布 (Stable Distributions)

对于独立同分布且方差有限的序列，其标准化的和收敛于正态分布。一个自然的问题是：**如果我们放宽方差有限的条件，标准化的和还能收敛到哪些非退化的分布？**
（想想 Cauchy 分布的例子）。

这就引出了稳定分布的概念。

!!! info "定义 3.17：稳定分布 (Stable Distribution)"

    一个分布 $F$ 被称为是**稳定的 (Stable)**，如果对于独立且服从 $F$ 的随机变量 $X_1, X_2$，以及任意非负常数 $c_1, c_2$，都存在常数 $a(c_1, c_2)$ 和 $b(c_1, c_2) > 0$，使得：
    
    $$
    c_1 X_1 + c_2 X_2 \stackrel{d}{=} b(c_1, c_2) X + a(c_1, c_2)
    $$
    
    其中 $X \sim F$ 且独立于 $X_1, X_2$。

**定理 3.18 (稳定分布的极限性质)：**
非退化的稳定分布族，恰好等价于所有可能的、通过适当中心化和标准化后的 i.i.d. 随机变量之和的**非退化极限分布族**。

### 3.1 稳定分布的谱表示 (Spectral Representation)

由于除了极少数情况（正态分布、Cauchy 分布、Lévy 分布），稳定分布没有闭式解析的概率密度函数，我们通常通过特征函数来刻画它们。

!!! success "定理 3.19 (稳定分布的特征函数)"

    一个稳定分布的特征函数具有如下形式：
    
    $$
    \phi_X(t) = E(e^{itX}) = \exp\left\{ i\gamma t - c|t|^\alpha (1 - i\beta \text{sgn}(t) z(t, \alpha)) \right\}
    $$
    
    其中参数范围为：位置参数 $\gamma \in \mathbb{R}$，尺度参数 $c > 0$，特征指数 $\alpha \in (0, 2]$，偏度参数 $\beta \in [-1, 1]$。并且有：
    
    $$
    z(t, \alpha) = \begin{cases} \tan\left(\frac{\pi \alpha}{2}\right), & \text{if } \alpha \ne 1 \\ -\frac{2}{\pi} \ln|t|, & \text{if } \alpha = 1 \end{cases}
    $$

!!! info "定义 3.20 ($\alpha$-稳定分布)"
    公式中的参数 $\alpha$ 被称为**特征指数 (characteristic exponent)**。对应的分布被称为 **$\alpha$-稳定分布 ($\alpha$-stable)**。

**关于参数的 Remark：**

* $\alpha = 2$ 对应正态分布 $N(\gamma, 2c)$。

* $\alpha = 1, \beta = 0$ 对应对称 Cauchy 分布。

* $\beta$ 描述了分布的偏度。当 $\beta = 0$ 时，特征函数是实值的，说明分布是对称的。

* 稳定分布在统计推断中存在困难，因为很难模拟，且极大似然估计难以直接写出，通常需要借助**经验特征函数 (Empirical Characteristic Functions)** 进行参数估计。

---

## 4. 吸引域 (Domain of Attraction, DA)

给定一个 $\alpha$-稳定分布 $G_\alpha$，我们需要知道什么样的原始分布 $F$，在经过怎样的标准化 $b_n > 0$ 和中心化 $a_n$ 后，其样本和的极限会落入 $G_\alpha$？

!!! info "定义 3.21：吸引域 (Domain of Attraction)"

    如果对于 i.i.d. 且服从分布 $F$ 的随机变量之和 $S_n = \sum_{i=1}^n X_i$，存在常数 $a_n \in \mathbb{R}$ 和 $b_n > 0$，使得：
    
    $$
    b_n^{-1} (S_n - a_n) \xrightarrow{d} G_\alpha
    $$
    
    则称分布 $F$ 属于 $G_\alpha$ 的**吸引域 (Domain of Attraction)**，记为 $X \in DA(G_\alpha)$ 或 $X \in DA(\alpha)$。

### 4.1 吸引域的特征刻画

为了刻画吸引域，需要用到**缓变函数 (Slowly varying function)** 的概念：若函数 $L$ 满足对所有 $t > 0$ 都有 $\lim_{x \to \infty} L(tx)/L(x) = 1$，则称 $L$ 为缓变函数。

!!! success "定理 3.22 与 推论 3.23 (吸引域的刻画)"

    1. **正态分布的吸引域 $DA(2)$:** $F \in DA(2)$ 当且仅当 $L(x) = \int_{|y|<x} y^2 dF(y)$ 是缓变函数。
       这等价于尾部概率满足：
       
    $$
    P(|X| > x) = o\left( x^{-2} \int_{|y|<x} y^2 dF(y) \right) \quad \text{当 } x \rightarrow \infty
    $$
       
       特别地，**所有二阶矩有限 ($E(X^2) < \infty$) 的分布都属于正态分布的吸引域**。
       
    2. **$\alpha < 2$ 稳定分布的吸引域 $DA(\alpha)$:**
       $F \in DA(\alpha)$ 当且仅当它的左右尾部具有 Pareto 衰减性质：
       
    $$
    F(-x) = \frac{c_1 + o(1)}{x^\alpha} L(x), \quad 1 - F(x) = \frac{c_2 + o(1)}{x^\alpha} L(x) \quad \text{当 } x \rightarrow \infty
    $$
       
       其中 $c_1, c_2 \ge 0$ 且 $c_1 + c_2 > 0$。

**推论 3.24 (矩的性质)：**
如果 $X \in DA(\alpha)$，那么：

* 对于 $\delta < \alpha$，有 $E(|X|^\delta) < \infty$。

* 对于 $\delta > \alpha$ (且 $\alpha < 2$)，有 $E(|X|^\delta) = \infty$。
（这也意味着如果 $\alpha < 2$，则方差必然无穷大；如果 $\alpha \le 1$，则均值也无穷大）。

---

### 4.2 标准化与中心化常数的选取

**命题 3.25 (标准化常数 $b_n$):**
对于 $F \in DA(\alpha)$，标准化常数 $b_n$ 可以选为如下方程的唯一解：

$$
G(b_n) + K(b_n) = n^{-1}, \quad n \ge 1
$$

其中 $G(x) = P(|X| > x)$，且 $K(x) = x^{-2} \int_{|y| \le x} y^2 dF(y)$。
如果 $E(X^2) < \infty$，那么 $b_n \sim \sigma \sqrt{n}$。如果 $\alpha < 2$，$b_n$ 的形式通常为 $n^{1/\alpha} L(n)$。

**命题 3.26 (中心化常数 $a_n$):**
中心化常数可以选为：

$$
a_n = n \int_{|y| \le b_n} y dF(y)
$$

---

## 5. 广义中心极限定理与正规吸引域 (DNA)

综合以上性质，我们得到最一般形式的中心极限定理：

!!! success "定理 3.27 (一般中心极限定理 General CLT)"

    假设 $F \in DA(\alpha)$，其中 $\alpha \in (0, 2]$。
    
    1. 如果 $E(X^2) < \infty$，那么：
       
    $$
    \frac{S_n - n\mu}{\sigma \sqrt{n}} \xrightarrow{d} N(0,1)
    $$
       
    2. 如果 $E(X^2) = \infty$ 且 $\alpha = 2$；或者 $\alpha < 2$，那么：
       
    $$
    \frac{S_n - a_n}{n^{1/\alpha} L_1(n)} \xrightarrow{d} G_\alpha
    $$
       
       其中 $G_\alpha$ 是某个 $\alpha$-稳定分布，$L_1$ 是适当的缓变函数。

!!! info "定义 3.28 与 推论 3.29：正规吸引域 (Domain of Normal Attraction, DNA)"

    注意到一般 CLT 的分母中包含缓变函数 $L_1(n)$。如果在极限中，标准化常数可以直接取**纯幂次形式** $b_n = c n^{1/\alpha}$（即 $L_1(n)$ 退化为常数），我们就称该分布属于 **正规吸引域 (DNA)**，记为 $F \in DNA(G_\alpha)$。
    
    * $F \in DNA(2)$ 当且仅当 $E(X^2) < \infty$。
    * 当 $\alpha < 2$ 时，$F \in DNA(\alpha)$ 当且仅当尾部严格服从幂律衰减（即没有额外的缓变函数干扰）：
      
    $$
    F(-x) \sim c_1 x^{-\alpha}, \quad 1 - F(x) \sim c_2 x^{-\alpha}
    $$

    特别地，每一个 $\alpha$-稳定分布都属于它自己的正规吸引域。