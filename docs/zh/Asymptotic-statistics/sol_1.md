---
tags:
  - 大样本理论
  - 统计学
  - 课后习题
---

# 📈 课后习题精解：大样本理论 Assignment 1

!!! abstract "关于本页"
    本页面收录了《大样本理论》课程第一单元（随机收敛、渐近正态性）的核心习题解答。所有解答均整理自课程提供的标准答案，并采用折叠框形式展示详细的证明与构造过程。

---

## 习题 1: Stochastic Convergence

**题目：**
(a) 寻找一个随机变量序列的例子，使得 $X_{n} \xrightarrow{d} 0$ 但 $E[X_{n}] \rightarrow \infty$。
(b) 寻找两个随机变量序列的例子，使得 $X_{n} \rightarrow X$ 且 $Y_{n} \rightarrow Y$，但其联合序列 $(X_{n}, Y_{n})$ 并不依分布收敛。

??? success "解答 (点击展开)"

    **(a) 构造期望发散但依分布收敛到 0 的序列**
    
    我们构造如下离散型随机变量 $X_n$：
    
    $$
    P(X_{n}=n^{2}) = \frac{1}{n}, \quad P(X_{n}=0) = 1 - \frac{1}{n}
    $$

    * **依概率收敛性**：对于任意 $\epsilon > 0$，当 $n$ 足够大时：
        $P(|X_{n}| > \epsilon) = P(X_{n} = n^{2}) = \frac{1}{n} \xrightarrow{n \rightarrow \infty} 0$。
        由于依概率收敛蕴含依分布收敛，故 $X_{n} \xrightarrow{d} 0$。
    * **期望的极限**：
        $E[X_{n}] = n^{2} \cdot \frac{1}{n} + 0 \cdot (1 - \frac{1}{n}) = n \xrightarrow{n \rightarrow \infty} \infty$。
    
    得证。

    <br>

    **(b) 构造联合分布不收敛的例子**
    
    设 $X \sim \mathcal{N}(0,1)$。定义序列 $X_{n}$ 和 $Y_{n}$ 如下：
    
    $$
    X_{n} = X \quad (\forall n), \quad Y_{n} = \begin{cases} X, & \text{if } n \text{ is odd} \\ -X, & \text{if } n \text{ is even} \end{cases}
    $$
    
    * **边缘分布收敛性**：显然 $X_{n} \xrightarrow{d} \mathcal{N}(0,1)$。由于标准正态分布具有对称性，$-X \sim \mathcal{N}(0,1)$，因此 $Y_{n}$ 的每个分量也服从 $\mathcal{N}(0,1)$，故 $Y_{n} \xrightarrow{d} \mathcal{N}(0,1)$。
    * **联合分布的非收敛性**：联合分布 $(X_{n}, Y_{n})$ 在 $(X, X)$ 和 $(X, -X)$ 之间交替切换。由于这两个联合分布显然不同，该序列的联合分布不收敛。

---

## 习题 2: Sums of Convergent Sequences

**题目：**
收敛性 $X_{n} \rightarrow X$ 和 $Y_{n} \rightarrow Y$ 是否总是蕴含 $X_{n} + Y_{n} \rightarrow X + Y$？此处 "$\rightarrow$" 可以代表依分布收敛、依概率收敛或几乎处处收敛。请针对其中至少一种情形证明结论或给出反例。

??? success "解答 (点击展开)"

    **结论：对于依分布收敛（Convergence in Distribution），该结论不成立。**

    **反例构造：**
    设 $X$ 和 $Y$ 是相互独立的标准正态随机变量，$X, Y \sim i.i.d. \mathcal{N}(0,1)$。
    定义序列：
    
    $$
    X_{n} = X \quad \text{和} \quad Y_{n} = -X \quad (\forall n)
    $$
    
    * **边缘分布收敛**：显然 $X_{n} \xrightarrow{d} X$。由于 $-X \sim \mathcal{N}(0,1)$，故 $Y_{n} \xrightarrow{d} Y$。
    * **和的收敛性**：在此构造下，$X_{n} + Y_{n} \equiv 0$，因此 $X_{n} + Y_{n} \xrightarrow{d} 0$。
    * **目标极限的分布**：然而，$X + Y$ 是两个独立正态变量之和，服从 $\mathcal{N}(0, 2)$ 分布，这显然不是在 0 处的退化分布。
    
    由于 $0$ 分布不等于 $\mathcal{N}(0, 2)$，原命题在依分布收敛下不成立。

---

## 习题 3: Convergence Properties and Cauchy Distribution

**题目：**
(a) 设 $X_{1}, ..., X_{n}$ 为一随机变量序列。证明 $X_{n} \xrightarrow{p} 0$ 当且仅当 $E\left(\frac{|X_{n}|}{1+|X_{n}|}\right) \rightarrow 0$。
(b) 设 $X_{1}, ..., X_{n}$ 是独立同分布的柯西分布 $Cauchy(0,1)$ 随机变量。证明 $\bar{X}_{n} \rightarrow Cauchy(0,1)$。
(c) 设 $X_{n}, 1 \le n \le \infty$ 为取整数值的随机变量。证明 $X_{n} \xrightarrow{d} X_{\infty}$ 当且仅当对于所有 $m$，$P(X_{n}=m) \rightarrow P(X_{\infty}=m)$ 成立。

??? success "解答 (点击展开)"

    **(a) 依概率收敛的等价特征证明**
    令 $f(x) = \frac{|x|}{1+|x|}$。注意 $f$ 是连续且在 $[0, \infty)$ 上严格递增的有界函数（$0 \le f(x) < 1$）。
    
    * **必要性 ($\Rightarrow$)**：若 $X_{n} \xrightarrow{p} 0$，根据连续映射定理，$f(X_{n}) \xrightarrow{p} 0$。由于 $f(X_{n})$ 被 1 有界控制，由有界收敛定理（Bounded Convergence Theorem）可得 $E[f(X_{n})] \rightarrow 0$。
    * **充分性 ($\Leftarrow$)**：对于任何 $\epsilon > 0$，由马尔可夫不等式（Markov's Inequality）：

        $$P(|X_{n}| > \epsilon) = P(f(X_{n}) > f(\epsilon)) \le \frac{E[f(X_{n})]}{f(\epsilon)}$$
        
        因为 $f(\epsilon) > 0$ 为常数且 $E[f(X_{n})] \rightarrow 0$，故 $P(|X_{n}| > \epsilon) \rightarrow 0$，即 $X_{n} \xrightarrow{p} 0$。

    <br>

    **(b) 柯西分布样本均值的性质**
    $Cauchy(0,1)$ 的特征函数为 $\varphi_{X}(t) = \exp(-|t|)$。
    利用独立性，样本均值 $\bar{X}_{n} = \frac{1}{n} \sum X_{k}$ 的特征函数为：
    
    $$
    \varphi_{\bar{X}_{n}}(t) = \prod_{k=1}^{n} \varphi_{X}\left(\frac{t}{n}\right) = \left(\exp\left(-\left|\frac{t}{n}\right|\right)\right)^{n} = \exp(-|t|)
    $$
    
    由于特征函数与分布一一对应，且 $\varphi_{\bar{X}_{n}}(t) = \varphi_{X}(t)$ 对所有 $n$ 成立，故 $\bar{X}_{n}$ 始终服从 $Cauchy(0,1)$。根据 Levy 连续性定理，$\bar{X}_{n} \xrightarrow{d} Cauchy(0,1)$。

    <br>

    **(c) 整数值随机变量的依分布收敛**
    
    * **必要性 ($\Rightarrow$)**：定义 $f_{m}(x) = \max(0, 1 - |x - m|)$。这是一个有界连续函数。由于 $X_n$ 是整数值，有 $E[f_{m}(X_{n})] = P(X_{n}=m)$。
        根据 Portmanteau 引理，$X_{n} \xrightarrow{d} X_{\infty}$ 蕴含对所有有界连续函数 $f$ 都有 $E[f(X_n)] \rightarrow E[f(X)]$，故 $P(X_{n}=m) \rightarrow P(X_{\infty}=m)$。
    * **充分性 ($\Leftarrow$)**：特征函数为 $\varphi_{X_{n}}(t) = \sum_{m \in \mathbb{Z}} e^{itm} P(X_{n}=m)$。由于级数被 1 控制且 $\sum P(X_{\infty}=m)=1$，根据级数的受控收敛定理（DCT）：

        $$\lim_{n \rightarrow \infty} \varphi_{X_{n}}(t) = \sum_{m \in \mathbb{Z}} e^{itm} \lim_{n \rightarrow \infty} P(X_{n}=m) = \varphi_{X_{\infty}}(t)$$
        
    由 Levy 连续性定理，结论成立。

---

## 习题 4: Asymptotic Normality

**题目：**
若 $X_{n}$ 为 $AN(\mu_{n}, \sigma_{n}^{2})$ 且 $Y_{n} \sim \mathcal{N}(\mu_{n}, \sigma_{n}^{2})$，证明：
(a) 当 $n \rightarrow \infty$ 时，$\sup_{t \in \mathbb{R}} |P(X_{n} \le t) - P(Y_{n} \le t)| \rightarrow 0$。
(b) $X_{n}$ 是 $AN(\bar{\mu}_{n}, \bar{\sigma}_{n}^{2})$ 当且仅当 $\frac{\bar{\sigma}_{n}}{\sigma_{n}} \rightarrow 1$ 且 $\frac{\mu_{n} - \bar{\mu}_{n}}{\sigma_{n}} \rightarrow 0$。
(c) $a_{n}X_{n} + b_{n}$ 是 $AN(\mu_{n}, \sigma_{n}^{2})$ 当且仅当 $a_{n} \rightarrow 1$ 且 $\frac{\mu_{n}(a_{n}-1) + b_{n}}{\sigma_{n}} \rightarrow 0$。

??? success "解答 (点击展开)"

    **(a) 累积分布函数的均匀收敛**
    令 $Z_{n} = \frac{X_{n} - \mu_{n}}{\sigma_{n}}$。根据渐近正态性定义，$Z_{n} \xrightarrow{d} \mathcal{N}(0,1)$。
    令 $\Phi(t)$ 为标准正态分布的 CDF。
    
    $$
    P(X_{n} \le t) = P\left(Z_{n} \le \frac{t - \mu_{n}}{\sigma_{n}}\right), \quad P(Y_{n} \le t) = \Phi\left(\frac{t - \mu_{n}}{\sigma_{n}}\right)
    $$
    
    由于 $\Phi(\cdot)$ 在全空间连续，根据 **Pólya 定理**，依分布收敛蕴含 CDF 的均匀收敛：
    
    $$\sup_{z \in \mathbb{R}} |P(Z_{n} \le z) - \Phi(z)| \rightarrow 0$$
    
    代入 $z = \frac{t - \mu_{n}}{\sigma_{n}}$ 即得证。

    <br>

    **(b) 参数替换的等价条件**
    定义 $W_{n} = \frac{X_{n} - \bar{\mu}_{n}}{\bar{\sigma}_{n}} = \frac{\sigma_{n}}{\bar{\sigma}_{n}} Z_{n} + \frac{\mu_{n} - \bar{\mu}_{n}}{\bar{\sigma}_{n}}$。
    
    * **充分性 ($\Leftarrow$)**：若极限条件成立，由 Slutsky 定理，$W_{n} \xrightarrow{d} 1 \cdot Z + 0 \sim \mathcal{N}(0,1)$。
    * **必要性 ($\Rightarrow$)**：若 $W_{n} \xrightarrow{d} \mathcal{N}(0,1)$，已知 $Z_{n} \xrightarrow{d} \mathcal{N}(0,1)$。根据 **收敛类型定理（Convergence of Types Theorem）**，尺度参数之比必须趋于 1，位置偏移之比必须趋于 0。

    <br>

    **(c) 线性变换的渐近正态性**
    令 $T_{n} = \frac{(a_{n}X_{n} + b_{n}) - \mu_{n}}{\sigma_{n}} = a_{n} Z_{n} + \frac{\mu_{n}(a_{n}-1) + b_{n}}{\sigma_{n}}$。
    同理利用 Slutsky 定理与收敛类型定理，$T_{n} \xrightarrow{d} \mathcal{N}(0,1)$ 当且仅当 $a_{n} \rightarrow 1$ 且偏移项 $\frac{\mu_{n}(a_{n}-1) + b_{n}}{\sigma_{n}} \rightarrow 0$。

---

## 习题 5: Portmanteau Lemma

**题目：**
证明 Portmanteau 引理中 (ii) 与 (iv) 的等价性：
(ii) 对于所有有界连续函数 $f$，$E[f(X_{n})] \rightarrow E[f(X)]$。
(iv) 对于所有非负连续函数 $f$，$\liminf_{n \rightarrow \infty} E[f(X_{n})] \ge E[f(X)]$。

??? success "解答 (点击展开)"

    **证明：**

    **1. 证明 (ii) $\Rightarrow$ (iv)：**
    设 $f \ge 0$ 为任一非负连续函数。
    对于任何常数 $M > 0$，定义截断函数 $f_{M}(x) = \min(f(x), M)$。
    显然 $f_{M}(x)$ 是有界且连续的。由于 $f(X_{n}) \ge f_{M}(X_{n})$，由期望的单调性可知：
    
    $$E[f(X_{n})] \ge E[f_{M}(X_{n})]$$
    
    对两边取极限下确界，并应用条件 (ii)：
    
    $$\liminf_{n \rightarrow \infty} E[f(X_{n})] \ge \liminf_{n \rightarrow \infty} E[f_{M}(X_{n})] = E[f_{M}(X)]$$
    
    当 $M \rightarrow \infty$ 时，$f_{M}(X) \uparrow f(X)$。根据单调收敛定理（Monotone Convergence Theorem）：
    
    $$E[f_{M}(X)] \rightarrow E[f(X)]$$
    
    因此，$\liminf_{n \rightarrow \infty} E[f(X_{n})] \ge E[f(X)]$。

    <br>

    **2. 证明 (iv) $\Rightarrow$ (ii)：**
    设 $f$ 为一个有界连续函数，且满足 $|f(x)| \le M$（$M > 0$）。
    那么 $f(x) + M \ge 0$ 且 $M - f(x) \ge 0$ 都是非负连续函数。
    
    * 对 $f(x) + M$ 应用 (iv)：
    
        $$\liminf_{n \rightarrow \infty} E[f(X_{n}) + M] \ge E[f(X) + M] \implies \liminf_{n \rightarrow \infty} E[f(X_{n})] \ge E[f(X)] \quad (*)$$
    
    * 对 $M - f(x)$ 应用 (iv)：
    
        $$\liminf_{n \rightarrow \infty} E[M - f(X_{n})] \ge E[M - f(X)] \implies M - \limsup_{n \rightarrow \infty} E[f(X_{n})] \ge M - E[f(X)]$$
        
        这等价于：
    
        $$\limsup_{n \rightarrow \infty} E[f(X_{n})] \le E[f(X)] \quad (**)$$
    
    结合 $(*)$ 和 $(**)$，我们得到：
    
    $$\limsup_{n \rightarrow \infty} E[f(X_{n})] \le E[f(X)] \le \liminf_{n \rightarrow \infty} E[f(X_{n})]$$
    
    由此说明极限存在且 $\lim_{n \rightarrow \infty} E[f(X_{n})] = E[f(X)]$。

---

## 习题 6: Basic rules of stochastic order

**题目：**
证明下列关于随机阶（stochastic order）的基本规则：
(a) $O_{p}(1)o_{p}(1) = o_{p}(1)$
(b) $(1 + o_{p}(1))^{-1} = O_{p}(1)$
(c) $o_{p}(O_{p}(1)) = o_{p}(1)$

??? success "解答 (点击展开)"

    **(a) 证明 $O_{p}(1)o_{p}(1) = o_{p}(1)$**
    设 $X_{n} = O_{p}(1)$ 且 $Y_{n} = o_{p}(1)$。我们要证明对于任意 $\epsilon > 0$，$P(|X_{n}Y_{n}| > \epsilon) \rightarrow 0$。
    对于任何固定常数 $M > 0$，有以下不等式成立：
    
    $$P(|X_{n}Y_{n}| > \epsilon) \le P(|X_{n}| > M) + P(|Y_{n}| > \frac{\epsilon}{M})$$
    
    给定 $\eta > 0$：
    1. 因为 $X_{n} = O_{p}(1)$，存在 $M > 0$ 和 $N_{1}$ 使得当 $n > N_{1}$ 时，$P(|X_{n}| > M) < \frac{\eta}{2}$。
    2. 固定此 $M$。因为 $Y_{n} = o_{p}(1)$，对于 $\frac{\epsilon}{M} > 0$，存在 $N_{2}$ 使得当 $n > N_{2}$ 时，$P(|Y_{n}| > \frac{\epsilon}{M}) < \frac{\eta}{2}$。
    
    令 $N = \max(N_{1}, N_{2})$，则当 $n > N$ 时：
    
    $$P(|X_{n}Y_{n}| > \epsilon) < \frac{\eta}{2} + \frac{\eta}{2} = \eta$$
    
    故 $X_{n}Y_{n} = o_{p}(1)$。

    <br>

    **(b) 证明 $(1 + o_{p}(1))^{-1} = O_{p}(1)$**
    设 $X_{n} = o_{p}(1)$。需证明对于任意 $\epsilon > 0$，存在 $M > 0$ 和 $N$ 使得当 $n > N$ 时，$P(|(1 + X_{n})^{-1}| > M) < \epsilon$。
    考虑 $M > 1$，事件 $\{|(1 + X_{n})^{-1}| > M\}$ 等价于 $\{|1 + X_{n}| < \frac{1}{M}\}$。
    根据反向三角不等式 $|1 + X_{n}| \ge 1 - |X_{n}|$，若 $|1 + X_{n}| < \frac{1}{M}$，则必有 $1 - |X_{n}| < \frac{1}{M}$，即 $|X_{n}| > 1 - \frac{1}{M}$。因此：
    
    $$P\left(\left|\frac{1}{1 + X_{n}}\right| > M\right) \le P\left(|X_{n}| > 1 - \frac{1}{M}\right)$$
    
    由于 $M > 1$，则 $1 - \frac{1}{M} > 0$。因为 $X_{n} = o_{p}(1)$，当 $n \rightarrow \infty$ 时该概率趋于 0。对于给定的 $\epsilon$，总能找到足够大的 $n$ 满足要求，故结论成立。

    <br>

    **(c) 证明 $o_{p}(O_{p}(1)) = o_{p}(1)$**
    设 $X_{n} = O_{p}(1)$ 且 $W_{n} = o_{p}(X_{n})$。
    根据小 $o_p$ 符号定义，$W_{n}$ 可以表示为 $W_{n} = Y_{n}X_{n}$，其中 $Y_{n} = o_{p}(1)$。
    于是我们有：
    
    $$o_{p}(X_{n}) = Y_{n}X_{n} = o_{p}(1) \cdot O_{p}(1)$$
    
    根据本题 (a) 部分已证的结论，$o_{p}(1)O_{p}(1) = o_{p}(1)$。
    因此 $W_{n} = o_{p}(1)$。

---

## 习题 7: Characteristic functions

**题目：**
设 $\phi$ 是 $\mathcal{R}^{k}$ 上的一个特征函数 (ch.f.)：
(a) 证明 $|\phi| \le 1$ 且在 $\mathcal{R}^{k}$ 上一致连续。
(b) 寻找两个随机变量 $X$ 和 $Y$ 的例子，使得 $X, Y$ 不独立，但它们的特征函数满足乘法公式 $\phi_{X}(t)\phi_{Y}(t) = \phi_{X+Y}(t)$ 对所有 $t \in \mathcal{R}$ 成立。

??? success "解答 (点击展开)"

    **(a) 证明 $|\phi| \le 1$ 及其一致连续性**
    
    * **有界性**：
       根据定义 $\phi(t) = E[\exp(it^{T}X)]$。利用期望的性质：
       
       $|\phi(t)| = |E[\exp(it^{T}X)]| \le E[|\exp(it^{T}X)|] = E[1] = 1$
    

    * **一致连续性**：
       对于任意 $t, h \in \mathcal{R}^{k}$，考虑差值：
       
       $|\phi(t + h) - \phi(t)| = |E[\exp(i(t + h)^{T}X) - \exp(it^{T}X)]|$
       
       $\le E[|\exp(it^{T}X)(\exp(ih^{T}X) - 1)|] = E[|\exp(ih^{T}X) - 1|]$
       
       注意到上界 $E[|\exp(ih^{T}X) - 1|]$ 与 $t$ 无关。
       设随机变量 $Z_{h} = |\exp(ih^{T}X) - 1|$。显然 $Z_{h} \le 2$ 且当 $h \rightarrow 0$ 时，$Z_{h} \xrightarrow{a.s.} 0$。
       根据勒贝格控制收敛定理（LDCT）：
       
       $\lim_{h \rightarrow 0} E[|\exp(ih^{T}X) - 1|] = E[0] = 0$
       
       由于此趋近过程对所有 $t$ 是均匀的，故 $\phi(t)$ 在 $\mathcal{R}^{k}$ 上一致连续。

    <br>

    **(b) 构造不独立但特征函数满足乘法公式的例子**
    
    设 $X \sim Cauchy(0,1)$，并令 $Y = X$。
    显然 $X$ 和 $Y$ 具有完全的相关性，因此它们不是独立的。
    
    * 标准柯西分布的特征函数为 $\phi_{X}(t) = \exp(-|t|)$。
    * 因为 $Y = X$，故 $\phi_{Y}(t) = \exp(-|t|)$。
    * 两者特征函数的乘积为：
        
        $$\phi_{X}(t)\phi_{Y}(t) = \exp(-|t|) \cdot \exp(-|t|) = \exp(-2|t|)$$
    
    考虑它们的和 $X + Y = 2X$。根据特征函数的性质：
    
    $$\phi_{X+Y}(t) = E[\exp(it(2X))] = E[\exp(i(2t)X)] = \phi_{X}(2t)$$
    
    代入柯西分布的特征函数形式：
    
    $$\phi_{X}(2t) = \exp(-|2t|) = \exp(-2|t|)$$
    
    我们观察到 $\phi_{X}(t)\phi_{Y}(t) = \exp(-2|t|) = \phi_{X+Y}(t)$ 对所有 $t \in \mathcal{R}$ 均成立。
    此反例证明了特征函数的乘法性质并不是判断随机变量独立的充分条件。