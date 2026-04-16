# 第七章：期权定价与停时

在随机微分方程的金融应用中，Black-Scholes 期权定价模型是最为经典的里程碑。本章我们将利用 Itô 积分推导 Black-Scholes 偏微分方程，探讨停时（Stopping Time）的严格定义及其在 Itô 积分中的性质，最后引入扩散过程的生成元与 Feynman-Kac 公式，建立 SDE 与 PDE（偏微分方程）之间的深刻联系。

---

## 1. Black-Scholes 期权定价模型

我们考虑一个标准的欧式看涨期权（European Call Option），到期日为 $T$，执行价格为 $K$。在到期日，期权的收益（Payoff）为：

\[
V(S_T, T) = (S_T - K)^+ = \max(S_T - K, 0)
\]

### 1.1 基本假设

在推导模型之前，我们基于 Black-Scholes 框架设定以下假设：

1. **标的资产动态**：股票价格 $S_t$ 服从几何布朗运动：

1. **标的资产动态**：股票价格 $S_t$ 服从几何布朗运动：

    \[
    \frac{dS_t}{S_t} = \mu dt + \sigma dW_t
    \]

2. **市场环境**：存在无风险利率 $r$（以连续复利计），且标的股票**不分红（无股息）**。
3. **交易摩擦**：不存在交易费用和税费，且可以无限分割和卖空。
4. **无套利原则**：市场中不存在无风险套利机会（No Arbitrage）。

### 1.2 Black-Scholes 偏微分方程的推导

基于无套利定价原理和自融资策略（Self-financing），我们可以构造一个无风险的投资组合，从而推导出期权价格 $V(S_t, t)$ 必须满足的偏微分方程。

!!! abstract "定理：Black-Scholes 偏微分方程 (BS PDE)"

    设期权价格 $V(S,t)$ 是关于 $S$ 和 $t$ 足够光滑的函数，则 $V$ 必定满足如下的抛物型偏微分方程：

    $$
    \partial_t V + \frac{1}{2}\sigma^2 S^2 \partial_S^2 V + rS \partial_S V - rV = 0
    $$

    **终端条件（对于看涨期权）**：
    
    $$
    V(S,T) = (S - K)^+
    $$

??? proof "推导证明：Itô 引理与 Delta 对冲（微元法）（点击展开）"

    **步骤 1：应用 Itô 引理**
    
    已知标的资产动态为 $dS_t = \mu S_t dt + \sigma S_t dW_t$。对期权价格函数 $V(S_t, t)$ 应用 Itô 引理：

    $$
    dV(S_t, t) = \partial_t V dt + \partial_S V dS_t + \frac{1}{2} \partial_S^2 V (dS_t)^2
    $$

    代入 $(dS_t)^2 = \sigma^2 S_t^2 dt$，我们得到 $V_t$ 的微分为：

    $$
    dV_t = \left( \partial_t V + \mu S_t \partial_S V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V \right) dt + \sigma S_t \partial_S V dW_t
    $$

    **步骤 2：构造无风险投资组合 (Delta Hedging)**
    
    构造一个投资组合 $\Pi_t$，包含 1 份期权多头和 $\Delta$ 份股票空头：

    $$
    \Pi_t = V_t - \Delta S_t
    $$

    在微小时间 $dt$ 内，该投资组合的价值变化（自融资微元法）为：

    $$
    d\Pi_t = dV_t - \Delta dS_t
    $$

    将 $dV_t$ 和 $dS_t$ 展开并代入：

    $$
    d\Pi_t = \left( \partial_t V + \mu S_t \partial_S V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V - \Delta \mu S_t \right) dt + (\sigma S_t \partial_S V - \Delta \sigma S_t) dW_t
    $$

    **步骤 3：消除随机项与无套利定价**
    
    为了让投资组合无风险（即消除 $dW_t$ 项），我们选择股票持有量 $\Delta$ 为：

    $$
    \Delta = \partial_S V
    $$

    代入上式，随机项完全抵消，得到：

    $$
    d\Pi_t = \left( \partial_t V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V \right) dt
    $$

    根据**无套利原则**，一个无风险投资组合的收益率必须等于无风险利率 $r$。即 $d\Pi_t = r \Pi_t dt$：

    $$
    \left( \partial_t V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V \right) dt = r (V_t - S_t \partial_S V) dt
    $$

    化简并移项，即得 Black-Scholes 偏微分方程：

    $$
    \partial_t V + \frac{1}{2} \sigma^2 S_t^2 \partial_S^2 V + rS \partial_S V - rV = 0
    $$
    
    推导完毕。 $\square$

---

## 2. Black-Scholes 偏微分方程的求解

对于 BS PDE，我们可以通过一系列变量代换，将其转化为经典的热传导方程（Heat Equation），从而利用高斯核求出解析解。

### 2.1 转化为热传导方程

??? proof "求解过程：热传导方程转换与解析解推导（点击展开）"

    **步骤 1：对数变换与时间反转**
    
    令 $\tau = T - t$（距离到期日的时间），且 $x = \ln S \Rightarrow S = e^x$。
    
    设 $v(x, \tau) = V(e^x, T - \tau)$。根据链式法则求导并代入原 PDE，经过化简后方程变为常系数抛物型方程：

    $$
    \partial_\tau v - \frac{1}{2}\sigma^2 \partial_x^2 v - \left(r - \frac{1}{2}\sigma^2\right) \partial_x v + r v = 0
    $$

    终端条件变为了初始条件：$v(x, 0) = (e^x - K)^+$。

    **步骤 2：消除低阶导数项（变量代换法）**
    
    为了消去一阶空间导数项 $\partial_x v$ 和常数项 $rv$，我们进行指数形式的变量代换，令：

    $$
    v(x, \tau) = u(x, \tau) e^{\alpha \tau + \beta x}
    $$

    求偏导数并代入方程：
    
    $$
    \partial_\tau v = (\partial_\tau u + \alpha u) e^{\alpha \tau + \beta x}
    $$
    
    $$
    \partial_x v = (\partial_x u + \beta u) e^{\alpha \tau + \beta x}
    $$
    
    $$
    \partial_x^2 v = (\partial_x^2 u + 2\beta \partial_x u + \beta^2 u) e^{\alpha \tau + \beta x}
    $$

    将其代入并整理，令 $\partial_x u$ 和 $u$ 的系数全为 0，我们可以解出特定的 $\alpha$ 和 $\beta$：

    $$
    \begin{cases}
    \sigma^2 \beta + \left(r - \frac{1}{2}\sigma^2\right) = 0 \\
    \alpha - \frac{1}{2}\sigma^2 \beta^2 - \left(r - \frac{1}{2}\sigma^2\right)\beta + r = 0
    \end{cases}
    $$

    解得：
    
    $$
    \beta = -\frac{r}{\sigma^2} + \frac{1}{2}, \quad \alpha = -\frac{1}{2\sigma^2}\left(r - \frac{\sigma^2}{2}\right)^2 - r
    $$

    此时，原方程彻底化为了标准的热传导方程（扩散方程）：

    $$
    \partial_\tau u - \frac{1}{2}\sigma^2 \partial_x^2 u = 0
    $$

    **步骤 3：利用 Green 函数（高斯核）求解**
    
    热传导方程的解可以通过积分其基本解（Kernel）与初始条件获得：

    $$
    u(x, \tau) = \int_{-\infty}^{+\infty} K(x-\xi, \tau) \psi(\xi) d\xi
    $$

    其中核函数为：

    $$
    K(x-\xi, \tau) = \frac{1}{\sigma \sqrt{2\pi \tau}} \exp\left( -\frac{(x-\xi)^2}{2\sigma^2 \tau} \right)
    $$


### 2.2 利用 Gaussian 核进行解析推导

??? proof "深度推导：如何从 Gaussian 核算到最终定价公式（点击展开）"

    **步骤 1：利用 Green 函数写出通用解**
    热传导方程在全空间的解可以通过初始条件 $\psi(\xi)$ 与 Gaussian 核 $K$ 的卷积得到：

    \[
    u(x, \tau) = \int_{-\infty}^{+\infty} \frac{1}{\sigma \sqrt{2\pi \tau}} \exp\left( -\frac{(x-\xi)^2}{2\sigma^2 \tau} \right) \psi(\xi) d\xi
    \]

    **步骤 2：代入期权初始条件**
    代入 $\psi(\xi) = (e^\xi - K)^+ e^{-\beta \xi}$。由于 $(e^\xi - K)^+$ 仅在 $\xi > \ln K$ 时非零，积分区间缩短：

    \[
    u(x, \tau) = \frac{1}{\sigma \sqrt{2\pi \tau}} \int_{\ln K}^{\infty} e^{-\frac{(x-\xi)^2}{2\sigma^2 \tau}} (e^\xi - K) e^{-\beta \xi} d\xi
    \]

    **步骤 3：拆分积分与配方**
    将积分拆分为 $I_1$（含 $e^\xi$ 项）和 $I_2$（含 $K$ 项）：
    1. 对于 $I_1$，指数项合并为 $-\frac{(\xi - [x + \sigma^2 \tau(1-\beta)])^2}{2\sigma^2 \tau} + \text{Const}$。
    2. 对于 $I_2$，指数项合并为 $-\frac{(\xi - [x - \sigma^2 \tau \beta])^2}{2\sigma^2 \tau} + \text{Const}$。

    **步骤 4：变量替换归一化**
    通过令 $z = \frac{\xi - \text{中心}}{\sigma \sqrt{\tau}}$，将积分转化为标准正态分布累积函数 $N(d)$ 的形式。
    其中下界 $\ln K$ 转换后的负值即为 $d_1$ 和 $d_2$：

    \[
    d_1 = \frac{\ln(S/K) + (r + \sigma^2/2)\tau}{\sigma\sqrt{\tau}}, \quad d_2 = d_1 - \sigma\sqrt{\tau}
    \]

    **步骤 5：回到原变量**
    最后利用 $V = u e^{-(\alpha \tau + \beta x)}$ 还原，指数项与积分外的常数项精准抵消，最终得到 Black-Scholes 公式。 $\square$

**最终结论：Black-Scholes 定价公式**

欧式看涨期权的解析解为：

\[
V(S,t) = S N(d_1) - K e^{-r(T-t)} N(d_2)
\]

其中：

\[
d_1 = \frac{\ln(S/K) + (r + \frac{1}{2}\sigma^2)(T-t)}{\sigma \sqrt{T-t}}
\]

\[
d_2 = d_1 - \sigma \sqrt{T-t}
\]

这里 $N(\cdot)$ 是标准正态分布的累积分布函数。

---

## 3. 停时 (Stopping Times) 及其性质

在研究复杂的金融衍生品（如美式期权，允许提前行权）时，传统的固定时间 $t$ 就不够用了，我们需要引入**停时**的概念。

### 3.1 停时的定义与基本性质

设 $(\Omega, \mathcal{F}, \{\mathcal{F}_t\}_{t \ge 0}, P)$ 带有给定的信息流（Filtration）。

!!! abstract "定义：停时"

    随机变量 $\tau: \Omega \to [0, +\infty]$ 被称为一个**停时 (Stopping Time)**，如果对于任意给定的时间 $t \ge 0$，事件 $\{\tau \le t\}$ 都是 $\mathcal{F}_t$ 可测的。

    即：在时间 $t$，通过截止到 $t$ 的信息 $\mathcal{F}_t$，我们一定能明确知道“是否已经停止” ($\tau \le t$)。

**停时的常见性质：**
假设 $\tau_1, \tau_2$ 都是停时，则：

1. **等价条件**：$\{\tau < t\} \in \mathcal{F}_t$。
2. **极小值也是停时**：$\min\{\tau_1, \tau_2\} = \tau_1 \wedge \tau_2$。
3. **极大值也是停时**：$\max\{\tau_1, \tau_2\} = \tau_1 \vee \tau_2$。

*(例子：股票价格首次击穿某个屏障 $H$ 的时间 $\tau = \inf\{t > 0: S_t \ge H\}$ 就是一个典型的停时。)*

### 3.2 带停时的 Itô 积分

我们可以将 Itô 积分的积分上限推广到停时 $\tau$：

\[
\int_0^\tau G_s dW_s \triangleq \int_0^t \chi_{\{s \le \tau\}} G_s dW_s
\]

这里 $\chi$ 是示性函数。带停时的 Itô 积分继承了标准的鞅性质和等距同构：

1. **期望为零**：
   
    \[
    E \left[ \int_0^\tau G_s dW_s \right] = 0
    \]

2. **Itô 等距**：
   
    \[
    E \left[ \left( \int_0^\tau G_s dW_s \right)^2 \right] = E \left[ \int_0^\tau G_s^2 ds \right]
    \]

---

## 4. 椭圆型方程与随机表示 (Elliptic Equations & Probabilistic Representation)

在研究布朗运动时，偏微分方程（PDE）的解往往可以表示为某种随机泛函的期望。我们设定 $X_t$ 为 $n$ 维标准布朗运动（即漂移项 $b=0$，扩散项 $B=I$），此时其生成元为 $L = \frac{1}{2}\Delta$。

### 4.1 泊松方程 (Poisson Equation) 与退出时间

!!! abstract "命题：退出时间的期望表示"

    考虑区域 $\Omega$ 上的泊松方程边值问题：

    $$
    \begin{cases}
    -\frac{1}{2}\Delta u(x) = 1, & x \in \Omega \\
    u(x) = 0, & x \in \partial \Omega
    \end{cases}
    $$

    则该方程的解可以表示为布朗运动从 $x$ 点出发首次离开区域 $\Omega$ 的期望时间：

    $$
    u(x) = E^x[\tau_\Omega]
    $$

??? proof "证明：利用 Dynkin 公式（点击展开）"

    **步骤 1：应用 Dynkin 公式（Itô 公式的期望形式）**

    对于 $u \in C^2(\Omega) \cap C(\bar{\Omega})$，根据 Itô 公式，有：

    $$
    u(X_{\tau \wedge t}) = u(X_0) + \int_0^{\tau \wedge t} \frac{1}{2}\Delta u(X_s) ds + \int_0^{\tau \wedge t} \nabla u(X_s) \cdot dW_s
    $$

    对两边取从 $x$ 点出发的期望 $E^x$。由于 Itô 积分项是鞅，其期望为 0：

    $$
    E^x[u(X_{\tau \wedge t})] = u(x) + E^x \left[ \int_0^{\tau \wedge t} \frac{1}{2}\Delta u(X_s) ds \right]
    $$

    **步骤 2：代入方程条件**

    已知 $\frac{1}{2}\Delta u = -1$。代入上式得：

    $$
    E^x[u(X_{\tau \wedge t})] = u(x) + E^x \left[ \int_0^{\tau \wedge t} (-1) ds \right] = u(x) - E^x[\tau \wedge t]
    $$

    **步骤 3：取极限 $t \to \infty$**

    根据边界条件，当布朗运动到达边界时 $u(X_\tau) = 0$。由控制收敛定理：

    $$
    0 = u(x) - E^x[\tau_\Omega] \implies u(x) = E^x[\tau_\Omega]
    $$

    推导完毕。 $\square$

---

### 4.2 调和函数 (Harmonic Functions) 的概率表示

!!! abstract "命题：狄里克雷问题的概率解"

    考虑拉普拉斯方程（调和函数）的边值问题：

    $$
    \begin{cases}
    \Delta u(x) = 0, & x \in \Omega \\
    u(x) = g(x), & x \in \partial \Omega
    \end{cases}
    $$

    则解 $u(x)$ 具有如下期望表示：

    $$
    u(x) = E^x[g(X_{\tau_\Omega})]
    $$

??? proof "证明（点击展开）"

    **步骤 1：利用调和性质**

    因为 $\Delta u = 0$，所以 $L u = 0$。根据 Itô 公式，过程 $u(X_t)$ 在离开 $\Omega$ 之前是一个局部鞅。

    **步骤 2：应用可选采样定理**

    对于停时 $\tau_\Omega$，我们有：

    $$
    u(x) = E^x[u(X_0)] = E^x[u(X_{\tau_\Omega})]
    $$

    **步骤 3：结合边界条件**

    当布朗运动到达边界时，$X_{\tau_\Omega} \in \partial \Omega$，此时 $u(X_{\tau_\Omega}) = g(X_{\tau_\Omega})$。代入得：

    $$
    u(x) = E^x[g(X_{\tau_\Omega})]
    $$

    这说明调和函数在一点的值等于其在边界上值的随机平均。 $\square$

---

### 4.3 边值不连续的狄里克雷问题

笔记中提到了将边界划分为两部分 $\Gamma_1$ 和 $\Gamma_2$ 的情况。

!!! abstract "命题：特征函数与击中概率"

    考虑如下边值问题：

    $$
    \begin{cases}
    \Delta u = 0, & x \in \Omega \\
    u = 1, & x \in \Gamma_1 \\
    u = 0, & x \in \Gamma_2
    \end{cases}
    $$

    其中 $\partial \Omega = \Gamma_1 \cup \Gamma_2$。

??? proof "证明：概率表示的直观理解（点击展开）"

    根据 4.2 的结论，解可以写为边界值的期望。令边界函数 $g(x) = I_{\Gamma_1}(x)$，即在 $\Gamma_1$ 上为 1，其余为 0。则：

    $$
    u(x) = E^x[I_{\Gamma_1}(X_{\tau_\Omega})] = P^x(X_{\tau_\Omega} \in \Gamma_1)
    $$

    **物理意义**：在该点处的解 $u(x)$ 等于布朗运动从该点出发，首次穿过边界时落在 $\Gamma_1$ 区域的概率。 $\square$

---

## 5. Feynman-Kac 公式 (Feynman-Kac Formula)

Feynman-Kac 公式建立了带有“势能项”或“折扣项”的偏微分方程与随机过程指数泛函之间的联系。

!!! abstract "定理：Feynman-Kac 公式"

    考虑如下二阶偏微分方程的边值问题：

    $$
    \begin{cases}
    -\frac{1}{2}\Delta u(x) + q(x)u(x) = f(x), & x \in \Omega \\
    u(x) = 0, & x \in \partial \Omega
    \end{cases}
    $$

    其中 $q(x) \ge 0$。则解 $u(x)$ 的随机表示为：

    $$
    u(x) = E^x \left[ \int_0^{\tau_\Omega} f(X_t) \exp\left( -\int_0^t q(X_s) ds \right) dt \right]
    $$

??? proof "证明：利用辅助过程的 Itô 展开（点击展开）"

    **步骤 1：构造辅助过程**

    定义过程 $M_t$ 如下：

    $$
    M_t = u(X_t) \exp\left( -\int_0^t q(X_s) ds \right) + \int_0^t f(X_s) \exp\left( -\int_0^s q(X_r) dr \right) ds
    $$

    **步骤 2：应用 Itô 公式求导**

    令 $e_t = \exp\left( -\int_0^t q(X_s) ds \right)$。利用乘积法则：

    $$
    d(u(X_t) e_t) = e_t du(X_t) + u(X_t) de_t
    $$

    其中：
    
    $$
    du(X_t) = \nabla u \cdot dW_t + \frac{1}{2}\Delta u dt
    $$
    
    $$
    de_t = -q(X_t) e_t dt
    $$

    代入得：

    $$
    d(u(X_t) e_t) = e_t \left( \frac{1}{2}\Delta u - q(X_t)u(X_t) \right) dt + e_t \nabla u \cdot dW_t
    $$

    **步骤 3：结合 PDE 条件**

    由于 $-\frac{1}{2}\Delta u + qu = f \implies \frac{1}{2}\Delta u - qu = -f$。代入 $dM_t$：

    $$
    dM_t = d(u(X_t) e_t) + f(X_t)e_t dt = \left[ e_t(-f) + e_t f \right] dt + e_t \nabla u \cdot dW_t
    $$

    因此，$dM_t = e_t \nabla u \cdot dW_t$，说明 $M_t$ 是一个局部鞅。

    **步骤 4：取期望与边界条件**

    由可选采样定理，$u(x) = M_0 = E^x[M_{\tau_\Omega}]$。
    
    当 $t = \tau_\Omega$ 时，$u(X_{\tau_\Omega}) = 0$，因此第一项消失，剩下：

    $$
    u(x) = E^x \left[ \int_0^{\tau_\Omega} f(X_s) \exp\left( -\int_0^s q(X_r) dr \right) ds \right]
    $$

    证明完毕。 $\square$