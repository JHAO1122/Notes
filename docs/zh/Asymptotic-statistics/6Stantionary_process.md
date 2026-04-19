# 第六章：弱相依数据（二）(Weakly Dependent Data II)

## 4. 弱相依平稳过程与中心极限定理

通过前面的混合系数和协方差不等式，我们现在可以正式定义时间序列中的“弱相依”与“长记忆”现象，并建立混合过程的中心极限定理 (CLT)。

### 4.1 弱相依与长记忆的定义

假设 $\{X_i\}$ 是一个具有有限二阶矩的弱平稳过程。记其自协方差函数为 $\gamma(j) = Cov(X_i, X_{i+j})$。

!!! info "定义 4.8：弱相依与长记忆 (Weakly Dependent & Long Memory)"

    如果自协方差绝对可和：

    \[
    \sum_{k=0}^{\infty} |\gamma(k)| < \infty
    \]

    则称该过程是**弱相依的 (weakly dependent)**（或短记忆过程）。
    
    如果自协方差绝对和发散：

    \[
    \sum_{k=0}^{\infty} |\gamma(k)| = \infty
    \]

    则称该过程是**长记忆过程 (long memory process)**。

**混合系数与弱相依的关系：**

设 $\alpha(k)$ 为由 $\{X_i\}_{i \in \mathbb{Z}}$ 生成的 $\sigma$-代数上的强混合（$\alpha$-混合）系数。
根据 Davydov 不等式（引理 4.7，取 $r=q$），如果存在 $q > 2$ 使得 $E|X_i|^q < \infty$，且对于 $p = \frac{q}{q-2}$ 满足级数收敛条件 $\sum_{k=0}^{\infty} \alpha^{1/p}(k) < \infty$，那么：

\[
\sum_{k=0}^{\infty} |\gamma(k)| \le 2p \|X\|_q^2 \sum_{k=0}^{\infty} \alpha^{1/p}(k) < \infty
\]

这说明满足该混合衰减速度的过程必然是**弱相依（短记忆）**的。
特别地，如果过程是**几何强混合 (Geometric Strong Mixing, GSM)** 的，即 $\alpha(k) \le C\rho^k$ (其中 $\rho \in (0,1)$)，那么：

\[
\sum_{k=0}^{\infty} \alpha^{1/p}(k) \le C \sum_{k=0}^{\infty} \rho^{k/p} = \frac{C}{1 - \rho^{1/p}} < \infty
\]

这自动保证了弱相依性。一般而言，为了保证级数收敛，我们只需要混合系数满足多项式衰减 $\alpha(k) \sim k^{-p(1+\eta)}$（$\eta > 0$）即可。

---

### 4.2 强混合过程的渐近方差与 CLT

在应用中心极限定理时，样本均值的方差极限是一个核心量。

!!! success "引理 4.9 (强混合过程的渐近方差)"

    设 $\{X_t\}_{t \in \mathbb{Z}}$ 是一个零均值、实值的弱平稳过程。假设存在 $r > 2$，使得：

    \[
    \sup_{t \in \mathbb{Z}} E|X_t|^r < \infty, \quad \sum_{k \ge 1} \alpha(k)^{1 - \frac{2}{r}} < +\infty
    \]

    那么，级数 $\sum_{k \in \mathbb{Z}} \gamma(k)$ 绝对收敛，且收敛到一个非负的常数 $\sigma^2$。并且部分和 $S_n = \sum_{t=1}^n X_t$ 的方差满足：

    \[
    \lim_{n \rightarrow \infty} n Var\left(\frac{S_n}{n}\right) = \sigma^2
    \]

??? proof "引理 4.9 的证明（点击展开）"

    首先，我们利用 Davydov 不等式（引理 4.7）来研究级数 $\sum_{k \in \mathbb{Z}} \gamma(k)$ 的绝对收敛性。
    在引理 4.7 中取 $q = r$，根据参数关系式 $\frac{1}{q} + \frac{1}{r} = 1 - \frac{1}{p}$，解得：

    \[
    \frac{1}{p} = 1 - \frac{2}{r} \implies p = \frac{r}{r-2}
    \]

    代入 Davydov 不等式，得到协方差 $\gamma(k) = Cov(X_0, X_k)$ 的界：

    \[
    |\gamma(k)| \le \frac{2r}{r-2} (E|X_0|^r)^{2/r} (2\alpha(k))^{1 - 2/r}
    \]

    由于题设条件已知 $\sum_{k \ge 1} \alpha(k)^{1 - 2/r} < +\infty$，因此由比较判别法，级数 $\sum_{k \in \mathbb{Z}} |\gamma(k)|$ 绝对收敛。

    接下来考察标准化和的方差。由于 $\{X_t\}$ 是弱平稳的，我们可以展开其方差：

    \[
    n Var\left(\frac{S_n}{n}\right) = n^{-1} \sum_{0 \le s, t \le n-1} Cov(X_s, X_t) = \sum_{k=-(n-1)}^{n-1} \left( 1 - \frac{|k|}{n} \right) \gamma(k)
    \]

    因为 $\gamma(k)$ 绝对可和，根据控制收敛定理（或 Kronecker 引理），当 $n \rightarrow \infty$ 时，权重 $\left( 1 - \frac{|k|}{n} \right) \rightarrow 1$。因此：

    \[
    \lim_{n \rightarrow \infty} n Var\left(\frac{S_n}{n}\right) = \sum_{k=-\infty}^{\infty} \gamma(k) = \sigma^2 \ge 0
    \]

    方差的非负性保证了极限 $\sigma^2 \ge 0$。得证。 $\square$

有了渐近方差的保证，我们可以直接陈述 $\alpha$-混合过程的中心极限定理：

!!! success "定理 4.13 ($\alpha$-混合过程的 CLT)"

    设 $\{X_t\}_{t \in \mathbb{Z}}$ 为零均值、实值的严格平稳过程。假设存在 $r > 2$ 和 $\beta > 0$ 使得：

    \[
    E|X_t|^r < \infty, \quad \alpha(k) \le a k^{-\beta}
    \]

    其中常数 $a > 0$ 且衰减阶数 $\beta > r / (r - 2)$。如果长期方差 $\sigma^2 = \sum_{k=-\infty}^{\infty} \gamma(k) > 0$，那么我们有：

    \[
    \frac{S_n}{\sigma \sqrt{n}} \xrightarrow{d} N(0, 1)
    \]

---

### 4.3 耦合方法与指数不等式 (Coupling Method & Exponential Inequalities)

为了证明复杂的混合极限定理，统计学家发明了**耦合方法 (Coupling Method)**。其核心思想是：将相依的平稳随机序列，构造性地替换为具有相同分布的**独立**序列，从而利用独立序列的成熟结论。

!!! success "引理 4.10 (Bradley's Lemma / 耦合引理)"

    设 $(X, Y)$ 是取值于 $\mathbb{R}^d \times \mathbb{R}$ 的随机向量，且对于某个 $p \in [1, \infty)$ 有 $Y \in L^p(P)$。
    设 $c$ 为实数使得 $\|Y+c\|_p > 0$，且 $\xi \in (0, \|Y+c\|_p]$。那么，存在一个辅助随机变量 $Y^*$ 满足：
    
    1. $P_{Y^*} = P_Y$（即 $Y^*$ 与 $Y$ 同分布），并且 $Y^*$ 独立于 $X$。
    2. 它们之间的距离受混合系数控制：

    \[
    P(|Y - Y^*| > \xi) \le 11 \left( \xi^{-1} \|Y+c\|_p \right)^{\frac{p}{2p+1}} \{ \alpha(\sigma(X), \sigma(Y)) \}^{\frac{2p}{2p+1}}
    \]

借助 Bradley 引理，我们可以将经典的独立序列指数不等式（如 Hoeffding 和 Bernstein 不等式）推广到混合序列中。

**回顾：独立序列的经典不等式 (定理 4.11)**
设 $X_1, \dots, X_n$ 为独立零均值随机变量，$S_n = \sum X_i$。

* **Hoeffding 不等式**：若 $a_i \le X_i \le b_i$，则 $P(|S_n| \ge t) \le 2 \exp\left\{ -\frac{2t^2}{\sum (b_i - a_i)^2} \right\}$。

* **Bernstein 不等式**：若满足 Cramér 条件 $E|X_i|^p \le c^{p-2} p! EX_i^2 < \infty$，则 $P(|S_n| \ge t) \le 2 \exp\left\{ -\frac{2t^2}{4\sum EX_i^2 + 2ct} \right\}$。

**推广：混合序列的指数不等式 (定理 4.12, Bosq 1998)**
对于零均值实值过程 $(X_t)$，若其一致有界 $\sup_t \|X_t\|_\infty \le b$。则对于整数 $q \in [1, n/2]$ 和 $\epsilon > 0$：

\[
P(|S_n| \ge n\epsilon) \le 4 \exp\left(-\frac{\epsilon^2}{8b^2} q\right) + 22 \left(1 + \frac{4b}{\epsilon}\right)^{1/2} q \alpha\left(\left[\frac{n}{2q}\right]\right)
\]

*(注：该不等式的第一项类似独立序列的指数衰减，第二项则是由于相依性引入的由 $\alpha$-混合系数控制的惩罚项。)*

---

## 5. 长期协方差 $\sigma^2$ 的谱估计方法 (Spectral Method)

前面我们在 CLT 中定义了长期协方差 $\sigma^2 = \sum_{k=-\infty}^{\infty} \gamma(k)$。在实际数据中，我们需要对它进行一致估计。频域方法（谱分析）为此提供了一个极为优雅的视角。

!!! info "定义与定理 4.14：谱密度函数 (Spectral Density Function)"

    定义时间序列 $\{X_t\}$ 的**谱密度函数** $f(\lambda)$ 为自协方差函数 $\gamma(k)$ 的傅里叶变换 (Fourier Transform)：
    对于频率 $\lambda \in (-\pi, \pi)$：

    \[
    f(\lambda) = \frac{1}{2\pi} \sum_{k=-\infty}^{\infty} \gamma(k) \exp(-ik\lambda)
    \]

    如果 $\sum |\gamma(k)| < \infty$，则过程存在谱密度 $f(\lambda)$。特别地，在频率 $\lambda = 0$ 处：

    \[
    \sum_{k=-\infty}^{\infty} \gamma(k) = 2\pi f(0)
    \]

**核心思想：** 估计长期方差 $\sigma^2$ 的问题，等价于估计谱密度在零频率处的值 $2\pi f(0)$。

### 5.1 周期图 (Periodograms)

为了估计谱密度，我们引入周期图的概念。给定样本 $\{X_1, \dots, X_n\}$，在其傅里叶频率 $\omega_j = 2\pi j / n \in [-\pi, \pi]$ 处，周期图定义为：

\[
l_n(\omega_j) = \frac{1}{n} \left| \sum_{t=1}^n X_t e^{-it\omega_j} \right|^2, \quad j \in T = \{0, \pm 1, \dots, \pm [n/2]\}
\]

它可以等价地展开为样本自协方差函数 $\hat{\gamma}(k)$ 的傅里叶变换：

\[
\begin{cases} 
l_n(0) = n|\overline{X}|^2 \\ 
l_n(\omega_j) = \sum_{|k|<n} \hat{\gamma}(k) e^{-ik\omega_j} & \text{if } \omega_j \ne 0 
\end{cases}
\]

其中 $\hat{\gamma}(k) = n^{-1} \sum_{t=1}^{n-|k|} (X_t - \overline{X})(X_{t+|k|} - \overline{X})$。

为了在连续频率域上分析，我们定义**扩展周期图 (Extended Periodogram)** $I_n(\omega)$：对于任意 $\omega \in [-\pi, \pi]$，将 $I_n(\omega)$ 定义为离 $\omega$ 最近的傅里叶频率点处的 $l_n(\omega_k)$ 值（即阶梯函数）。

!!! success "命题 4.16 (周期图的期望性质)"

    如果 $\{X_t\}$ 是均值为 $\mu$ 且自协方差绝对可和的平稳序列，那么：

    \[
    E(I_n(0)) - n\mu^2 \rightarrow 2\pi f(0)
    \]

    对于非零频率 $\omega \ne 0$：

    \[
    E(l_n(\omega)) \rightarrow 2\pi f(\omega)
    \]

    特别地，如果真实均值 $\mu = 0$，则 $E(I_n(\omega))$ 在 $[-\pi, \pi]$ 上一致收敛于 $2\pi f(\omega)$。

??? proof "命题 4.16 的证明（点击展开）"

    首先，对于零频率 $\omega = 0$：

    \[
    E(l_n(0)) - n\mu^2 = n E(\overline{X}^2) - n\mu^2 = n Var(\overline{X}) = n Var(S_n/n)
    \]

    由引理 4.9 的结论，当 $n \rightarrow \infty$ 时，上式收敛于 $\sum_{k=-\infty}^{\infty} \gamma(k) = 2\pi f(0)$。

    现在考察 $\omega \in (0, \pi]$。利用 $l_n(\omega)$ 的等价表示展开其期望：

    \[
    E(l_n(\omega)) = \sum_{|k|<n} \left( 1 - \frac{|k|}{n} \right) \gamma(k) e^{-ik g(n,\omega)}
    \]

    其中 $g(n,\omega)$ 是离 $\omega$ 最近的傅里叶频率。
    因为自协方差 $\gamma(\cdot)$ 是绝对可和的，序列 $\sum_{|k|<n} \left( 1 - \frac{|k|}{n} \right) \gamma(k) e^{-ik\lambda}$ 会一致收敛到其傅里叶变换 $2\pi f(\lambda)$。
    
    又因为当 $n \to \infty$ 时，$g(n,\omega) \rightarrow \omega$，故有：

    \[
    E(l_n(\omega)) \rightarrow 2\pi f(\omega)
    \]

    如果 $\mu=0$，利用 $f$ 在闭区间 $[-\pi, \pi]$ 上的均匀连续性，可得一致收敛结论。得证。 $\square$

### 5.2 线性过程周期图的渐近分布

!!! success "定理 4.17"

    设 $\{X_t\}$ 为线性过程 $X_t = \sum \psi_j \epsilon_{t-j}$ 且 $\sum |\psi_j| < \infty$，$\epsilon_t \sim i.i.d.\ F(0, \sigma^2)$。
    如果谱密度 $f(\lambda) > 0$，那么对于 $m$ 个不同的频率 $0 < \lambda_1 < \dots < \lambda_m < \pi$，随机向量：

    \[
    (I_n(\lambda_1), \dots, I_n(\lambda_m))^T
    \]

    依分布收敛于一个由**独立的指数分布 (Exponential Distribution)** 随机变量组成的向量，其中第 $i$ 个分量的均值为 $2\pi f(\lambda_i)$。

*(注：该定理不仅给出了极限分布，还指出了不同频率上的周期图在渐近意义上是不相关的，这为后续的频域回归奠定了基础。)*

---

## 6. 谱密度的非参数核平滑估计

定理 4.17 揭示了一个严重的问题：**单一的周期图 $I_n(\omega)$ 并不是谱密度 $2\pi f(\omega)$ 的一致估计量**。
（因为它的极限是一个方差不为 0 的指数分布，而不是退化到常数）。且当 $\mu \ne 0$ 时，$I_n(0)$ 甚至是有偏的。

**解决方案：** 既然不同频率 $\omega_j$ 处的周期图是渐近独立的，且它们的期望围绕着真实谱密度 $f(\omega)$ 波动，我们可以将相邻频率的周期图进行**局部加权平均 (Locally weighted averaging)**，这就是**非参数核回归 (Nonparametric Kernel Regression)** 的核心思想。

### 6.1 构建对数周期图回归模型

我们可以将周期图的关系写成一个乘性模型：

\[
I_n(\omega_j) = 2\pi f(\omega_j) e_j + R_j \quad \text{for } j \in T \setminus \{0\}
\]

其中 $\{e_j\}$ 是独立的 $Exp(1)$ 随机变量，$\{R_j\}$ 是高阶可忽略项。
对两边取自然对数，转化为标准的加性非参数回归模型：

\[
\log\left(\frac{l_n(\omega_j)}{2\pi}\right) = \log(f(\omega_j)) + \log(e_j)
\]

已知指数分布对数的性质：$E(\log(e_j)) = -0.57721$ (欧拉常数) 且 $Var(\log(e_j)) = \pi^2/6$。
通过中心化，令：

* $\eta_j = \log(e_j) + 0.57721$ （零均值，方差 $\pi^2/6$ 的 i.i.d. 误差项）

* $W_j = \log\left(\frac{l_n(\omega_j)}{2\pi}\right) + 0.57721$ （新的响应变量）

* $m(\omega) = \log(f(\omega))$ （未知的目标平滑函数）

于是我们得到了一个固定设计的非参数回归模型：

\[
W_j = m(\omega_j) + \eta_j, \quad j \in T \setminus \{0\}
\]

### 6.2 Nadaraya-Watson (NW) 核估计量

现在我们的目标是通过核平滑估计 $m(0) = \log(f(0))$。
给定核函数 $K(\cdot)$ 和平滑带宽 (bandwidth) $b$，$m(\omega)$ 的 NW 估计量定义为：

\[
\hat{m}_b(\omega) = \frac{\sum_{j \in T} K\left(\frac{\omega - \omega_j}{b}\right) W_j}{\sum_{j \in T} K\left(\frac{\omega - \omega_j}{b}\right)}
\]

其中带宽需要满足渐近条件：当 $n \rightarrow \infty$ 时，$b \rightarrow 0$ 且 $nb \rightarrow \infty$。

最后，通过指数逆变换，我们就能得到长期协方差（零频率谱密度）的一致估计：

\[
\hat{f}(0) = \exp(\hat{m}_b(0))
\]

有了 $\hat{f}(0)$，我们就成功得到了 CLT 标准化中所需的长期方差 $\sigma^2 = 2\pi \hat{f}(0)$，从而闭环了整个弱相依数据的渐近推断体系。