# 🎲 Probability Theory

本模块涵盖了概率论中的核心不等式、随机变量的收敛性理论以及大样本极限定理。

## 一、 核心概率不等式 (Inequalities)

!!! info "Markov 不等式"
    对于非负随机变量 $X \ge 0$ 及任意 $a > 0$：
    
    $$
    P(X \ge a) \le \frac{E[X]}{a}
    $$

!!! info "Chebyshev 不等式"
    对于任意随机变量 $X$ 及 $a > 0$：

    $$
    P(|X - E[X]| \ge a) \le \frac{Var(X)}{a^2}
    $$

!!! success "Kolmogorov 极大值不等式"
    若 $X_1, \dots, X_n$ 独立同分布且均值为 0，方差有限：
    
    $$
    P\left(\max_{1\le k\le n} \left|\sum_{i=1}^k X_i\right| \ge \lambda\right) \le \frac{1}{\lambda^2} \sum_{i=1}^n Var(X_i)
    $$

??? note "Borel-Cantelli 引理"
    记事件序列 $\{A_n\}$ 无限次发生（infinitely often, i.o.）为 $\limsup A_n$：
    
    1. **第一引理**：若 $\sum_{n=1}^\infty P(A_n) < \infty$，则 $P(A_n \text{ i.o.}) = 0$。
    2. **第二引理**：若 $\{A_n\}$ **相互独立**，且 $\sum_{n=1}^\infty P(A_n) = \infty$，则 $P(A_n \text{ i.o.}) = 1$。

---

## 二、 随机变量的收敛性 (Modes of Convergence)

!!! abstract "四种收敛定义"
    1. **依概率收敛 ($X_n \xrightarrow{P} X$)**：对于任意 $\epsilon > 0$，有 $\lim_{n \to \infty} P(|X_n - X| > \epsilon) = 0$。
    2. **几乎处处收敛 ($X_n \xrightarrow{a.s.} X$)**：$P(\lim_{n \to \infty} X_n = X) = 1$。
    3. **$L^p$ 收敛 ($X_n \xrightarrow{L^p} X$)**：$\lim_{n \to \infty} E[|X_n - X|^p] = 0$。
    4. **依分布收敛 ($X_n \xrightarrow{d} X$)**：
        * **定义**：对于极限分布函数 $F(x)$ 的所有连续点 $x$，满足 $\lim_{n \to \infty} F_n(x) = F(x)$。
        * **Lévy 连续性定理**：$X_n \xrightarrow{d} X$ 当且仅当其特征函数 $\phi_{X_n}(t) \to \phi_X(t)$ 对每个 $t$ 成立，且 $\phi_X(t)$ 在 $t=0$ 处连续。

??? tip "收敛性之间的转换关系 (Important!)"
    * **蕴含链**：$L^p \implies P \implies d$ 以及 $a.s. \implies P \implies d$。
    * **反向转换补充**：
        * **$P \to L^p$**：若 $X_n \xrightarrow{P} X$ 且 $\{|X_n|^p\}$ 是**一致可积 (Uniformly Integrable, UI)** 的，则有 $X_n \xrightarrow{L^p} X$。
        * **$P \to a.s.$**：若 $X_n \xrightarrow{P} X$，则必定存在**子序列** $\{n_k\}$ 使得 $X_{n_k} \xrightarrow{a.s.} X$。
        * **$d \to P$**：若极限分布 $X=c$ 是一个退化的常数，则依分布收敛可推回依概率收敛。
    * **Slutsky 定理**：若 $X_n \xrightarrow{d} X$ 且 $Y_n \xrightarrow{P} c$（常数），则 $X_n Y_n \xrightarrow{d} cX$，$X_n + Y_n \xrightarrow{d} X + c$。

---

## 三、 大数定律与中心极限定理 (Limit Theorems)

!!! info "辛钦大数定律 (Khinchin's WLLN)"
    设 $X_1, X_2, \dots, X_n$ 是**独立同分布 (i.i.d.)** 的随机变量序列，且其期望 $E[X_i] = \mu$ 有限（$E[|X_1|] < \infty$）。
    则样本均值**依概率收敛**于 $\mu$：

    $$
    \bar{X}_n = \frac{1}{n}\sum_{i=1}^n X_i \xrightarrow{P} \mu \quad (n \to \infty)
    $$

!!! info "Kolmogorov 强大数定律 (SLLN)"
    若 $X_i$ i.i.d.，则 $E[|X_1|] < \infty$ 是 $\bar{X}_n \xrightarrow{a.s.} E[X]$ 的充要条件。

??? success "中心极限定理 (CLT) 与成立条件"
    对于独立但不一定同分布的序列 $X_1, X_2, \dots$，记 $E[X_i]=\mu_i, Var(X_i)=\sigma_i^2$，$S_n^2 = \sum_{i=1}^n \sigma_i^2$。
    标准化后的和为 $Z_n = \frac{1}{S_n} \sum_{i=1}^n (X_i - \mu_i)$。
    
    * **Lindeberg 条件（充要条件）**：对于任意 $\epsilon > 0$，满足：

    $$
    \lim_{n \to \infty} \frac{1}{S_n^2} \sum_{i=1}^n E[(X_i - \mu_i)^2 I(|X_i - \mu_i| > \epsilon S_n)] = 0
    $$
    
    *(直觉：保证了序列中没有任何单个随机变量的尾部极值在总方差中占主导地位)*

    * **Lyapunov 条件（充分条件）**：若存在 $\delta > 0$，使得：

    $$
    \lim_{n \to \infty} \frac{1}{S_n^{2+\delta}} \sum_{i=1}^n E[|X_i - \mu_i|^{2+\delta}] = 0
    $$
    
    *(说明：Lyapunov 条件往往更容易验证，因为只要高阶矩存在并满足该极限，通过简单的积分放缩就能直接推导出 Lindeberg 条件成立，从而 $Z_n \xrightarrow{d} N(0,1)$)*

---

## 四、 常见概率分布速查表 (Common Distributions)

| 分布名称 | 记号 | 期望 $E[X]$ | 方差 $Var(X)$ | 概率密度/质量函数 (PDF/PMF) | 特征函数 $\phi_X(t)$ | 矩母函数 $M_X(t)$ |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **正态分布** | $N(\mu, \sigma^2)$ | $\mu$ | $\sigma^2$ | $\frac{1}{\sqrt{2\pi}\sigma} e^{-\frac{(x-\mu)^2}{2\sigma^2}}$ | $e^{it\mu - \frac{1}{2}\sigma^2 t^2}$ | $e^{t\mu + \frac{1}{2}\sigma^2 t^2}$ |
| **指数分布** | $Exp(\lambda)$ | $\frac{1}{\lambda}$ | $\frac{1}{\lambda^2}$ | $\lambda e^{-\lambda x} \quad (x \ge 0)$ | $\frac{\lambda}{\lambda - it}$ | $\frac{\lambda}{\lambda - t} \quad (t < \lambda)$ |
| **均匀分布** | $U(a, b)$ | $\frac{a+b}{2}$ | $\frac{(b-a)^2}{12}$ | $\frac{1}{b-a} \quad (a \le x \le b)$ | $\frac{e^{itb} - e^{ita}}{it(b-a)}$ | $\frac{e^{tb} - e^{ta}}{t(b-a)}$ |
| **泊松分布** | $Pois(\lambda)$ | $\lambda$ | $\lambda$ | $e^{-\lambda} \frac{\lambda^k}{k!} \quad (k \in \mathbb{N})$ | $e^{\lambda(e^{it} - 1)}$ | $e^{\lambda(e^t - 1)}$ |
| **二项分布** | $Bin(n, p)$ | $np$ | $np(1-p)$ | $\binom{n}{k} p^k (1-p)^{n-k}$ | $(1 - p + pe^{it})^n$ | $(1 - p + pe^t)^n$ |

---

## 五、 变量代换与 PDF 转换 (Variable Transformations)

!!! info "一元连续型随机变量代换 (Univariate)"
    设 $X$ 的 PDF 为 $f_X(x)$，令 $Y = g(X)$。若 $g(x)$ 是严格单调且可微的函数，其反函数记为 $x = g^{-1}(y)$，则 $Y$ 的 PDF 为：
    
    $$
    f_Y(y) = f_X(g^{-1}(y)) \cdot \left| \frac{d}{dy} g^{-1}(y) \right|
    $$

!!! success "多元变量代换与雅可比矩阵 (Multivariate & Jacobian)"
    设 $\mathbf{X} = (X_1, \dots, X_n)^T$ 的联合 PDF 为 $f_{\mathbf{X}}(\mathbf{x})$，令 $\mathbf{Y} = \mathbf{g}(\mathbf{X})$ 为 $n$ 维一一映射。
    记反函数为 $\mathbf{x} = \mathbf{h}(\mathbf{y})$，其中 $\mathbf{h} = \mathbf{g}^{-1}$。则 $\mathbf{Y}$ 的联合 PDF 为：
    
    $$
    f_{\mathbf{Y}}(\mathbf{y}) = f_{\mathbf{X}}(\mathbf{h}(\mathbf{y})) \cdot | \det(J_{\mathbf{h}}(\mathbf{y})) |
    $$
    
    其中 $J_{\mathbf{h}}(\mathbf{y})$ 是 **Jacobian 矩阵**，其定义为：
    
    $$
    J_{\mathbf{h}}(\mathbf{y}) = \frac{\partial (x_1, \dots, x_n)}{\partial (y_1, \dots, y_n)} = \begin{pmatrix} 
    \frac{\partial x_1}{\partial y_1} & \dots & \frac{\partial x_1}{\partial y_n} \\ 
    \vdots & \ddots & \vdots \\ 
    \frac{\partial x_n}{\partial y_1} & \dots & \frac{\partial x_n}{\partial y_n} 
    \end{pmatrix}
    $$
