# 第五章：多元 Itô 积分与随机微分方程

前面的章节已经建立了一维布朗运动下的 Itô 积分理论。然而在金融数学与统计物理中，系统往往受到多个独立噪声源的影响。本章我们将 Itô 积分推广到高维空间，建立多元 Itô 公式，并正式引入刻画连续时间随机动态系统的工具——**随机微分方程 (Stochastic Differential Equations, SDE)**。

## 1. 多元 Itô 积分与等距同构

我们首先定义多维的布朗运动及其对应的信息流。

!!! abstract "定义：多元布朗运动与适应矩阵过程"

    设 $W(t) = (W^1(t), W^2(t), \dots, W^m(t))^T$ 是一个 $m$ 维标准布朗运动，它的每一个分量 $W^j(t)$ 都是独立的一维标准布朗运动。
    
    定义其产生的信息流为 $\mathcal{F}(t) = \sigma\{W(s) \mid 0 \le s \le t\}$。对于任意 $t < s$，未来的增量 $W(s) - W(t)$ 与 $\mathcal{F}(t)$ 完全独立。

    设 $G(t, \omega) = (G^{ij}(t, \omega))_{n \times m}$ 是一个 $n \times m$ 阶的随机矩阵过程。如果 $G(t)$ 的每一个分量 $G^{ij}(t)$ 都是 $\mathcal{F}(t)$-适应的，且满足平方可积条件：
    
    $$
    E\left[ \int_0^T \sum_{i=1}^n \sum_{j=1}^m (G^{ij}(t))^2 dt \right] < \infty
    $$
    
    则称矩阵过程 $G \in L_{n \times m}^2(0, T)$。

对于这样的矩阵过程，多元 Itô 积分被自然地定义为向量形式。

!!! info "定义与定理：多元 Itô 积分"

    对于 $G \in L_{n \times m}^2(0, T)$，我们定义多元 Itô 积分 $\int_0^T G dW$ 为一个 $n$ 维列向量，其第 $i$ 个分量定义为：
    
    $$
    \left( \int_0^T G dW \right)^i = \sum_{j=1}^m \int_0^T G^{ij} dW^j
    $$
    
    在这个定义下，多元 Itô 积分依然保持完美的鞅性质和等距性质：
    
    **(1) 期望为零**：$E\left[ \int_0^T G dW \right] = 0$ （即每个分量的期望均为 0）。
    
    **(2) 多元 Itô 等距同构 (Multivariate Itô Isometry)**：定义矩阵的范数为 $|G| = \sqrt{\sum_{i,j} (G^{ij})^2}$，则有：
    
    $$
    E\left[ \left| \int_0^T G dW \right|^2 \right] = E\left[ \int_0^T |G(t)|^2 dt \right]
    $$

    ??? proof "多元 Itô 等距同构的严格证明（点击展开）"

        由向量范数的定义，积分模长的平方可以写为各个分量平方的和：
        
        $$
        E\left[ \left| \int_0^T G dW \right|^2 \right] = E\left[ \sum_{i=1}^n \left( \sum_{j=1}^m \int_0^T G^{ij} dW^j \right)^2 \right]
        $$
        
        将内部的平方项展开为双重求和形式：
        
        $$
        = \sum_{i=1}^n E\left[ \sum_{j=1}^m \int_0^T G^{ij} dW^j \sum_{k=1}^m \int_0^T G^{ik} dW^k \right]
        $$
        
        由于求和与期望是线性算子，我们交换顺序，重点考察交叉项：
        
        $$
        = \sum_{i=1}^n \sum_{j=1}^m \sum_{k=1}^m E\left[ \int_0^T G^{ij} dW^j \int_0^T G^{ik} dW^k \right]
        $$
        
        **关键点：布朗运动分量间的独立性**。
        当 $j \neq k$ 时，$W^j$ 和 $W^k$ 是相互独立的布朗运动。利用一维随机积分的极化恒等式和独立增量性质，这部分的交叉期望为 $0$。
        
        因此，非零项仅存在于对角线 $j = k$ 处：
        
        $$
        = \sum_{i=1}^n \sum_{j=1}^m E\left[ \left( \int_0^T G^{ij} dW^j \right)^2 \right]
        $$
        
        对每一项应用一维的 Itô 等距同构 $E[(\int G^{ij} dW^j)^2] = E[\int (G^{ij})^2 dt]$：
        
        $$
        = \sum_{i=1}^n \sum_{j=1}^m E\left[ \int_0^T (G^{ij})^2 dt \right] = E\left[ \int_0^T \sum_{i,j} (G^{ij})^2 dt \right] = E\left[ \int_0^T |G(t)|^2 dt \right]
        $$
        
        证明完毕。$\square$

---

## 2. 多元 Itô 公式与乘法表

在进入随机微分方程的求解前，我们必须掌握多维情况下的“链式法则”。

!!! abstract "定义：多元 Itô 过程 (Multivariate Itô Process)"

    一个 $n$ 维 Itô 过程 $X(t) = (X^1(t), \dots, X^n(t))^T$ 的微元形式定义为：
    
    $$
    dX(t) = F(t) dt + G(t) dW(t)
    $$
    
    写成分量形式即为：
    
    $$
    dX^i(t) = F^i(t) dt + \sum_{j=1}^m G^{ij}(t) dW^j(t), \quad i = 1, \dots, n
    $$
    
    其中漂移项 $F \in L_n^1(0, T)$，扩散项 $G \in L_{n \times m}^2(0, T)$。

处理多元随机微分的核心在于推广版的 **Itô 乘法表**。由于不同的布朗运动分量是独立的，它们的“二次交叉变差”为 0。

!!! info "定理：多元 Itô 乘法法则"

    对于时间微元 $dt$ 和不同维度的布朗运动微元 $dW^j, dW^k$，有如下乘法规则：
    
    * $(dt)^2 = 0$
    * $dt \cdot dW^j = 0$
    * $dW^j \cdot dW^k = \delta_{jk} dt = \begin{cases} dt, & j = k \\ 0, & j \neq k \end{cases}$

    利用上述规则，我们可以直接得出两个一维 Itô 过程 $X_1, X_2$ 的**乘积法则 (Product Rule)**：
    
    $$
    d(X_1 X_2) = X_1 dX_2 + X_2 dX_1 + dX_1 dX_2
    $$

    ??? proof "乘积法则中交叉项 $dX_1 dX_2$ 的推导计算（点击展开）"

        设两个过程分别为 $dX_1 = F_1 dt + \sum_j G_1^j dW^j$，且 $dX_2 = F_2 dt + \sum_k G_2^k dW^k$。
        展开它们的乘积：
        
        $$
        dX_1 dX_2 = \left( F_1 dt + \sum_j G_1^j dW^j \right) \left( F_2 dt + \sum_k G_2^k dW^k \right)
        $$
        
        根据多元乘法表，所有含有 $dt$ 的二次项（如 $(dt)^2, dt \cdot dW$）均视为高阶无穷小，直接消去：
        
        $$
        = \sum_{j=1}^m \sum_{k=1}^m G_1^j G_2^k dW^j dW^k
        $$
        
        由于 $dW^j dW^k = \delta_{jk} dt$，只有当 $j=k$ 时非零，和式坍缩为一重求和：
        
        $$
        = \sum_{j=1}^m G_1^j G_2^j dt
        $$
        
        将其写成向量内积的形式，即为 $(G_1^T G_2) dt$。这就是二次变差带来的修正项。$\square$

基于多元乘法表，我们可以将泰勒展开推广，得到现代随机分析最重要的基石。

!!! success "定理：多元 Itô 公式 (Multivariate Itô's Formula)"

    设 $X(t)$ 为一个 $n$ 维 Itô 过程，$U(t, x_1, \dots, x_n)$ 为 $\mathbb{R}^+ \times \mathbb{R}^n \to \mathbb{R}$ 的具有 $C^{1,2}$ 连续偏导数的多元纯量函数。
    则复合随机过程 $U(t, X(t))$ 依然是一个 Itô 过程，且其微分为：
    
    $$
    dU(t, X(t)) = \frac{\partial U}{\partial t} dt + \sum_{i=1}^n \frac{\partial U}{\partial x^i} dX^i + \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \frac{\partial^2 U}{\partial x^i \partial x^j} dX^i dX^j
    $$
    
    如果将 $dX^i$ 展开并代入 $dW^j dW^k = \delta_{jk} dt$，上述公式可写为更紧凑的矩阵迹 (Trace) 形式：
    
    $$
    dU = \left( U_t + \nabla U^T F + \frac{1}{2} \text{Tr}(G^T D_x^2 U G) \right) dt + (\nabla U^T G) dW
    $$
    
    其中 $\nabla U$ 为空间梯度向量，$D_x^2 U$ 为 Hessian 矩阵。

---

## 3. 线性随机微分方程 (LSDE) 求解实战

我们现在正式考虑形如 $dX(t) = b(t, X(t)) dt + \sigma(t, X(t)) dW(t)$ 的随机微分方程。由于非线性 SDE 的显式解往往不存在，本节我们重点剖析一类在金融学中具有统治地位的方程：**线性随机微分方程 (Linear SDE)**。

### 3.1 指数鞅与 Hermite 多项式

考虑最简单的齐次线性方程，其漂移项为 0：

$$
\begin{cases} dY(t) = \lambda Y(t) dW(t) \\ Y(0) = 1 \end{cases}
$$

这代表系统的随机波动与系统当前的状态成正比。

!!! tip "方程的求解方法 (Itô 公式反向构造)"

    我们猜想解为某种指数形式 $Y(t) = e^{f(t, W(t))}$。为了抵消 Itô 公式展开时产生的 $\frac{1}{2} (dW)^2 = \frac{1}{2}dt$ 修正项，我们在指数中预先加入“补偿项”。
    
    定义函数 $U(x, t) = e^{\lambda x - \frac{\lambda^2}{2}t}$，将 $x$ 替换为 $W(t)$。应用一维 Itô 公式：
    
    $$
    dU(W(t), t) = \frac{\partial U}{\partial t} dt + \frac{\partial U}{\partial x} dW + \frac{1}{2} \frac{\partial^2 U}{\partial x^2} (dW)^2
    $$
    
    代入偏导数 $U_t = -\frac{\lambda^2}{2} U$，$U_x = \lambda U$，$U_{xx} = \lambda^2 U$：
    
    $$
    dY(t) = -\frac{\lambda^2}{2} Y(t) dt + \lambda Y(t) dW(t) + \frac{1}{2} \lambda^2 Y(t) dt = \lambda Y(t) dW(t)
    $$
    
    这完美匹配了目标方程。因此其解为 **指数鞅 (Exponential Martingale)**：
    
    $$
    Y(t) = e^{\lambda W(t) - \frac{\lambda^2}{2}t}
    $$

这个解在概率论中有着极深的内涵，它的泰勒展开与正交多项式直接挂钩。

??? note "拓展：Hermite 多项式与随机积分的递推关系"
    
    回顾概率论中的 Hermite 多项式母函数生成公式：
    
    $$
    e^{\lambda x - \frac{\lambda^2}{2}t} = \sum_{n=0}^{\infty} \frac{\lambda^n}{n!} H_n(x, t) = \sum_{n=0}^{\infty} \lambda^n h_n(x, t)
    $$
    
    我们将解 $Y(t)$ 的积分形式与级数展开形式进行对比。已知 $Y(t) = 1 + \lambda \int_0^t Y(s) dW(s)$，代入展开式：
    
    $$
    \sum_{n=0}^{\infty} \lambda^n h_n(W(t), t) = 1 + \lambda \int_0^t \sum_{n=0}^{\infty} \lambda^n h_n(W(s), s) dW(s)
    $$
    
    比较方程两边同次幂 $\lambda^{n+1}$ 的系数，我们得到了一个极其优美的关于布朗运动的积分递推公式：
    
    $$
    h_{n+1}(W(t), t) = \int_0^t h_n(W(s), s) dW(s)
    $$
    
    这个公式展示了如何通过反复对布朗运动进行随机积分，自然地生成出高阶的 Hermite 多项式（由于常数项 $h_0 = 1$，则 $h_1 = W(t)$，$h_2 = \frac{1}{2}(W(t)^2 - t)$，与上一章黎曼和极限完全一致）。

### 3.2 黑尔斯-斯科尔斯 (Black-Scholes) 模型

如果我们在上面的方程中加入线性的时间漂移，就得到了金融学中最著名的资产定价几何布朗运动模型 (GBM)。

$$
\frac{dS(t)}{S(t)} = \mu dt + \sigma dW(t)
$$

其中 $\mu$ 称为预期收益率（漂移项），$\sigma$ 称为资产波动率。

!!! success "几何布朗运动的精确解"

    我们将方程改写为 $dS(t) = \mu S(t) dt + \sigma S(t) dW(t)$。
    启发于前面的解法，我们考虑对数变换。设 $U(x) = \ln x$，根据 Itô 公式对 $\ln S(t)$ 求微分：
    
    $$
    d(\ln S(t)) = U'(S(t)) dS(t) + \frac{1}{2} U''(S(t)) (dS(t))^2
    $$
    
    代入导数 $U'(x) = 1/x$ 和 $U''(x) = -1/x^2$，以及 $(dS)^2 = \sigma^2 S^2 dt$：
    
    $$
    d(\ln S(t)) = \frac{1}{S(t)} (\mu S(t) dt + \sigma S(t) dW(t)) - \frac{1}{2} \frac{1}{S(t)^2} (\sigma^2 S(t)^2 dt)
    $$
    
    化简得到：
    
    $$
    d(\ln S(t)) = \left(\mu - \frac{\sigma^2}{2}\right) dt + \sigma dW(t)
    $$
    
    此时右侧不再依赖于 $S(t)$，变成了一个普通的积分！对方程两端从 $0$ 到 $t$ 积分：
    
    $$
    \ln S(t) - \ln S(0) = \left(\mu - \frac{\sigma^2}{2}\right) t + \sigma W(t)
    $$
    
    两边取指数，得到终极解：
    
    $$
    S(t) = S_0 e^{\left(\mu - \frac{\sigma^2}{2}\right)t + \sigma W(t)}
    $$

    这个结果告诉我们，虽然对数收益率的期望表现为线性的 $\mu t$，但由于方差的惩罚，实际资产价格的长期增长率被削弱为了 $\mu - \frac{\sigma^2}{2}$，这就是“波动性拖累 (Volatility Drag)”的数学本质。

---

## 4. 布朗桥 (Brownian Bridge) 

最后我们探讨一个具有强约束边界条件的 SDE 案例。

在统计学中，我们常常需要研究被“钉在”两个固定点之间的随机游走（例如 Kolmogorov-Smirnov 检验中的极限过程）。这就是布朗桥过程。

!!! abstract "定义：标准布朗桥方程"

    考虑以下具有奇异漂移项的随机微分方程，初始条件 $B(0) = 0$，时间限定在 $t \in [0, 1)$：
    
    $$
    dB(t) = -\frac{B(t)}{1-t} dt + dW(t)
    $$
    
    直观上理解：随着时间 $t$ 逼近 $1$，分母 $(1-t) \to 0$，漂移项会产生一股极其巨大的“拉力”，强行将 $B(t)$ 的轨迹拉回 $0$。

我们可以用经典 ODE 中的**积分因子法**配合 Itô 乘积法则来求解此 SDE。

!!! info "布朗桥的显式解"

    我们将含有 $B$ 的项移到方程左边：
    
    $$
    dB(t) + \frac{B(t)}{1-t} dt = dW(t)
    $$
    
    方程两边同乘积分因子 $\frac{1}{1-t}$：
    
    $$
    \frac{1}{1-t} dB(t) + \frac{B(t)}{(1-t)^2} dt = \frac{dW(t)}{1-t}
    $$
    
    观察左侧，它恰好是联合微元 $d\left(\frac{B(t)}{1-t}\right)$ 的展开（因为根据普通乘积求导法则，这里没有两个随机过程相乘，不存在二次修正项）：
    
    $$
    d\left(\frac{B(t)}{1-t}\right) = \frac{dW(t)}{1-t}
    $$
    
    对两边从 $0$ 到 $t$ 积分（注意 $B(0) = 0$）：
    
    $$
    \frac{B(t)}{1-t} - 0 = \int_0^t \frac{dW(s)}{1-s}
    $$
    
    整理即得布朗桥的显式积分解：
    
    $$
    B(t) = (1-t) \int_0^t \frac{dW(s)}{1-s}, \quad 0 \le t < 1
    $$

    ??? proof "端点 $t=1$ 处的连续性与零值证明框架（点击展开）"

        上述积分在 $t \to 1$ 时分母会出现奇异性。为了证明几乎所有轨道都能连续地“降落”在 $B(1) = 0$，需要借助高级概率论工具进行界定。
        
        首先分析积分部分的方差。记 $M(t) = \int_0^t \frac{dW(s)}{1-s}$。由于这是对于非随机函数的 Itô 积分，其期望为 0，方差由等距同构给出：
        
        $$
        E[M(t)^2] = \int_0^t \frac{1}{(1-s)^2} ds = \frac{1}{1-t} - 1 = \frac{t}{1-t}
        $$
        
        因此 $B(t) = (1-t) M(t)$ 的方差为 $E[B(t)^2] = (1-t)^2 \frac{t}{1-t} = t(1-t)$。显然当 $t \to 1$ 时，方差趋于 0。但这仅说明了依概率或 $L^2$ 收敛。
        
        为了证明**几乎必然 (a.s.) 连续收敛**，需采用 Borel-Cantelli 引理与极大值不等式。
        取二进网络点 $t_n = 1 - 2^{-n}$。对于区间 $[t_n, t_{n+1}]$，利用 Chebyshev 不等式和鞅的极大值不等式，估计在该区间内偏离目标值的概率界：
        
        $$
        P\left( \max_{t \in [t_n, t_{n+1}]} |B(t)| > \varepsilon \right) \le \frac{C}{\varepsilon^2} 2^{-\eta n}
        $$
        
        由于右侧构成收敛的几何级数 $\sum 2^{-\eta n} < \infty$，根据 Borel-Cantelli 第一引理，大偏差事件发生的概率为 0。因此，几乎所有的布朗桥轨道在 $t \to 1$ 时都被强制收敛到了 $0$。$\square$