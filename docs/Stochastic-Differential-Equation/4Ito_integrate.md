# 第四章：Itô 积分与随机微分(Itô Integral & Stochastic Differential)

上一章将积分推广到了 $L^2$ 空间，现在我们正式定义**随机过程**的 Itô 积分，并建立起一套区别于经典微积分的体系——**Itô 微积分 (Itô Calculus)**。

## 1. 适应过程的 Itô 积分与等距同构

对于一个随机过程 $G(t, \omega)$，由于布朗运动 $W(t)$ 的路径粗糙，我们无法采用逐点定义的黎曼-斯蒂尔杰斯积分，须基于时间网格的左端点（不知未来）定义。

!!! abstract "定义：Itô 积分 (The Itô Integral)"

    设 $\{W(t)\}_{t \ge 0}$ 为标准布朗运动，$\{\mathcal{F}(t)\}$ 为其信息流。
    假设随机过程 $G(t) \in L^2(0, T)$，且 $G(t)$ 是**适应过程 (Adapted Process)**（即在时刻 $t$ 的取值只依赖于历史信息 $\mathcal{F}(t)$）。

    对于区间 $[0, T]$ 的划分 $P = \{0 = t_0 < t_1 < \dots < t_m = T\}$，定义其 Riemann 和为：

    $$
    R_m = \sum_{k=0}^{m-1} G(t_k) \big( W(t_{k+1}) - W(t_k) \big)
    $$

    当网格模长 $|P| \to 0$ 时，如果该和式在 $L^2(\Omega, P)$ 意义下收敛，其极限即定义为 $G(t)$ 关于布朗运动的 **Itô 积分**：

    $$
    \int_0^T G(t) dW(t) = \lim_{|P| \to 0} \sum_{k=0}^{m-1} G(t_k) \big( W(t_{k+1}) - W(t_k) \big)
    $$

利用简单函数（阶梯过程）逼近的方法定义这一积分后有如下性质：

!!! info "定理 1：Itô 积分的核心性质"

    设 $G(t), H(t) \in L^2(0, T)$ 均为适应过程，$a, b \in \mathbb{R}$ 为常数。

    **(1) 线性性**：

    $$
    \int_0^T (aG(t) + bH(t)) dW(t) = a \int_0^T G(t) dW(t) + b \int_0^T H(t) dW(t) \quad a.s.
    $$

    **(2) 期望为零**：

    $$
    E\left[ \int_0^T G(t) dW(t) \right] = 0
    $$

    **(3) Itô 等距同构 (Itô Isometry)**：

    $$
    E\left[ \left( \int_0^T G(t) dW(t) \right)^2 \right] = \int_0^T E[G(t)^2] dt
    $$

    ??? proof "期望为零与 Itô 等距同构的严格证明（点击展开）"

        我们利用简单函数（Step Functions）对性质进行证明。设 $G_k = G(t_k)$。

        **1. 证明期望为零**：

        $$
        E\left[ \sum_{k} G_k \big(W(t_{k+1}) - W(t_k)\big) \right] = \sum_k E\Big[ G_k \big(W(t_{k+1}) - W(t_k)\big) \Big]
        $$

        由于 $G(t)$ 是适应过程，$G_k$ 关于 $\mathcal{F}(t_k)$ 可测；而布朗运动具有独立增量性质，增量 $\Delta W_k = W(t_{k+1}) - W(t_k)$ 与 $\mathcal{F}(t_k)$ **完全独立**。利用期望的独立乘法性质（或条件期望）：

        $$
        E[G_k \Delta W_k] = E[G_k] \cdot E[\Delta W_k] = E[G_k] \cdot 0 = 0
        $$

        因此和式的期望为 0。

        **2. 证明 Itô 等距同构**：
        展开积分平方的期望（利用 Fubini 定理交换期望与求和）：

        $$
        E\left[ \left(\sum_k G_k \Delta W_k \right)^2 \right] = E\left[ \sum_{k} \sum_{j} G_k G_j \Delta W_k \Delta W_j \right]
        $$

        将双重求和拆分为三部分：$k > j$、$k < j$ 以及 $k = j$。
        
        **分析交叉项 ($k \neq j$)**：不失一般性，假设 $k > j$。此时时间点满足 $t_j < t_{j+1} \le t_k < t_{k+1}$。
        在这四个随机变量 $G_k, G_j, \Delta W_j, \Delta W_k$ 中，前三个完全属于历史信息 $\mathcal{F}(t_k)$，而最后一个增量 $\Delta W_k$ 在 $t_k$ 之后，与 $\mathcal{F}(t_k)$ 独立。
        利用塔牌性质（Tower Property），先对 $\mathcal{F}(t_k)$ 求条件期望：

        $$
        E\Big[ G_k G_j \Delta W_j \Delta W_k \Big] = E\Big[ E\big[ G_k G_j \Delta W_j \Delta W_k \mid \mathcal{F}(t_k) \big] \Big]
        $$

        由于 $G_k, G_j, \Delta W_j$ 已知可提取：

        $$
        = E\Big[ G_k G_j \Delta W_j \underbrace{ E[\Delta W_k \mid \mathcal{F}(t_k)] }_{= 0} \Big] = 0
        $$

        因此，所有交叉项的期望全为 0。

        **分析对角项 ($k = j$)**：
        只剩下对角线上的平方项：

        $$
        \sum_{k} E\Big[ G_k^2 (\Delta W_k)^2 \Big]
        $$

        同样利用条件期望，将 $\Delta W_k$ 的平方提取出来：

        $$
        = \sum_{k} E\Big[ E\big[ G_k^2 (\Delta W_k)^2 \mid \mathcal{F}(t_k) \big] \Big] = \sum_k E\Big[ G_k^2 E\big[ (\Delta W_k)^2 \mid \mathcal{F}(t_k) \big] \Big]
        $$

        由于增量独立且方差为 $\Delta t_k = t_{k+1} - t_k$，所以 $E[(\Delta W_k)^2 | \mathcal{F}(t_k)] = \Delta t_k$。代入得：

        $$
        = \sum_k E[G_k^2] (t_{k+1} - t_k)
        $$

        当 $|P| \to 0$ 时，这个黎曼和直接收敛于黎曼积分 $\int_0^T E[G(t)^2] dt$。等距同构得证！$\square$

        > 简单函数逼近到$L^2(0, T)$空间即可。
---

## 2. 不定积分与连续鞅性质

如果我们把积分的上限 $T$ 换成一个变量 $t$，我们就得到了随机过程的**不定积分 (Indefinite Integral)**。这在本质上是在定义随机过程。

!!! abstract "定义：Itô 不定积分"

    设 $G \in L^2(0, T)$，定义其不定积分为随机过程 $I(t)$：

    $$
    I(t) = \int_0^t G(s) dW(s), \quad 0 \le t \le T
    $$

    显然初始条件为 $I(0) = 0$。

!!! success "定理：Itô 积分是连续平方可积鞅"

    由 Itô 积分定义的不定积分过程 $\{I(t)\}_{t \ge 0}$ 具有极其完美的数学性质：
    它不仅几乎必然具有**连续的样本轨道**，而且是一个关于自然信息流的**鞅 (Martingale)**。

    ??? proof "鞅性质证明（点击展开）"

        对于任意 $0 \le s \le t \le T$，我们需要证明 $E[I(t) \mid \mathcal{F}(s)] = I(s)$ a.s.。

        将区间 $[0, t]$ 在 $s$ 点拆断：

        $$
        I(t) = \int_0^s G(\tau) dW(\tau) + \int_s^t G(\tau) dW(\tau) = I(s) + \int_s^t G(\tau) dW(\tau)
        $$

        两边同时对 $\mathcal{F}(s)$ 取条件期望：

        $$
        E[I(t) \mid \mathcal{F}(s)] = E\left[ I(s) + \int_s^t G(\tau) dW(\tau) \bigg| \mathcal{F}(s) \right]
        $$

        由于 $I(s)$ 的积分域在 $[0, s]$ 内，它显然是 $\mathcal{F}(s)$-可测的，故已知即常数（直接提出）。
        对于后半部分，利用 Itô 积分期望为 0 的性质的条件版本：

        $$
        = I(s) + E\left[ \int_s^t G(\tau) dW(\tau) \bigg| \mathcal{F}(s) \right] = I(s) + 0 = I(s)
        $$

        因此，$I(t)$ 是一个鞅。

        *(注：关于连续性的严格证明需要用到 Doob 极大值不等式与 Borel-Cantelli 引理构造 $L^2$ Cauchy 列的一致收敛，手稿中有提及，思路与布朗运动构造极其相似，此处略去解析学繁冗细节)* $\square$

---

## 3. Itô 过程与乘积法则 (Integration by Parts)

有了积分，下一步就是研究其微分形式。

!!! abstract "定义：Itô 过程 (Itô Process) 与 SDE"

    设 $F(t) \in L^1(0,T)$，$G(t) \in L^2(0,T)$ 均为适应过程。定义随机过程 $X(t)$：

    $$
    X(t) = X(0) + \int_0^t F(s) ds + \int_0^t G(s) dW(s)
    $$

    上式通常被写为直观的**随机微分方程 (SDE)** 形式：

    $$
    dX(t) = F(t) dt + G(t) dW(t)
    $$

    这里 $F(t)dt$ 称为**漂移项 (Drift)**，代表确定性趋势；$G(t)dW(t)$ 称为**扩散项 (Diffusion)**，代表随机扰动。

在经典的莱布尼茨微积分中，$d(X_1 X_2) = X_1 dX_2 + X_2 dX_1$。但在随机分析中，由于布朗运动的二次变差不为 0，我们必须引入额外的**二次修正项**。

!!! info "定理：Itô 乘积法则 (Product Rule)"

    设有两个 Itô 过程 $dX_i = F_i dt + G_i dW \quad (i=1,2)$。
    那么它们的乘积 $X_1(t)X_2(t)$ 也是一个 Itô 过程，且满足：

    $$
    d(X_1 X_2) = X_1 dX_2 + X_2 dX_1 + dX_1 dX_2
    $$

    这里需要应用以下著名的**Itô 乘法表 (Itô Multiplication Table)** 展开 $dX_1 dX_2$：

    | $\times$ | $dt$ | $dW$ |
    | :---: | :---: | :---: |
    | **$dt$** | 0 | 0 |
    | **$dW$** | 0 | $dt$ |

    > *(核心直觉：因为 $dW \sim \sqrt{dt}$，所以 $(dW)^2 = dt$。而带有高于一次 $dt$ 的项在极限下皆为高阶无穷小 0)*。

    **直接展开应用**：

    $$
    dX_1 dX_2 = (F_1 dt + G_1 dW)(F_2 dt + G_2 dW) = G_1 G_2 (dW)^2 = G_1 G_2 dt
    $$

    所以展开后的全写形式为：

    $$
    d(X_1 X_2) = (X_1 F_2 + X_2 F_1 + G_1 G_2) dt + (X_1 G_2 + X_2 G_1) dW
    $$

    > **例 1：** $d(W^2) = W dW + W dW + (dW)^2 = 2W dW + dt$
    > （由此可立刻得出积分形式：$\int_0^t W dW = \frac{1}{2}W^2(t) - \frac{1}{2}t$，与上章黎曼和极限完全一致！）
    >
    > **例 2：** $d(tW) = t dW + W dt + dt \cdot dW = t dW + W dt$

---

## 4. 随机微分中的Itô 公式 (Itô's Formula)

Itô 公式是整个随机微积分的“链式法则 (Chain Rule)”，可用于求解各种非随机现象的方程。

!!! success "定理：Itô 公式"

    设 $X(t)$ 为一个 Itô 过程，$dX(t) = F dt + G dW$。
    设 $U(x, t): \mathbb{R} \times [0, T] \to \mathbb{R}$ 为一个二阶连续可导 ($C^{2,1}$) 的确定性函数。
    定义新的随机过程 $Y(t) = U(X(t), t)$。则 $Y(t)$ 的随机微分满足泰勒展开到二阶的截断形式：

    $$
    dU(X(t), t) = \frac{\partial U}{\partial t} dt + \frac{\partial U}{\partial x} dX + \frac{1}{2} \frac{\partial^2 U}{\partial x^2} (dX)^2
    $$

    将 $dX$ 的表达式和 Itô 乘法表 $(dX)^2 = G^2 dt$ 代入重整，得到标准的 **Itô 公式**：

    $$
    dU(X(t), t) = \left( \frac{\partial U}{\partial t} + \frac{\partial U}{\partial x} F + \frac{1}{2} \frac{\partial^2 U}{\partial x^2} G^2 \right) dt + \frac{\partial U}{\partial x} G dW
    $$

    ??? proof "Itô 公式多项式应用示例（点击展开）"

        针对多项式 $f(x) = x^m$ 进行推导验证，
        由 Itô 公式：

        $$
        df(X) = f'(X) dX + \frac{1}{2} f''(X) (dX)^2
        $$

        代入导数 $f'(x) = m x^{m-1}$ 和 $f''(x) = m(m-1) x^{m-2}$：

        $$
        d(X^m) = m X^{m-1} dX + \frac{1}{2} m(m-1) X^{m-2} G^2 dt
        $$

        这与使用前面乘积法则 $d(x \cdot x^{k-1})$ 逐步递推得出的结论完全一致，印证了随机微积分算子代数体系的自洽性。

---

## 5. Fokker-Planck 方程

有了 Itô 公式，我们不仅能研究微观的单条随机轨道，还能直接得出宏观的**概率密度函数的确定性规律**。这就是 Fokker-Planck 方程，也被称为 Kolmogorov 向前方程 (KFE)。

!!! abstract "推导：Fokker-Planck 方程"

    考虑一般的扩散过程 (SDE)：

    $$
    dX(t) = b(X, t) dt + \sigma(X, t) dW(t)
    $$

    设在时刻 $t$，系统处于状态 $x$ 的概率密度函数为 $p(x, t)$。
    我们取一个任意的有界且在无穷远处消失的二次可导**测试函数** $\phi(x)$。
    根据 Itô 公式展开 $d\phi(X(t))$：

    $$
    d\phi(X(t)) = \phi'(X) dX + \frac{1}{2} \phi''(X) (dX)^2 = \phi'(X) \big[b(X, t) dt + \sigma(X, t) dW\big] + \frac{1}{2} \phi''(X) \sigma^2(X, t) dt
    $$

    对等式两边同时**取期望**（此时扩散项 $\int \sigma dW$ 具有鞅性质，其期望消去）：

    $$
    E[\phi(X(t))] = E[\phi(X(0))] + \int_0^t E\left[ \phi'(X) b(X, s) + \frac{1}{2} \phi''(X) \sigma^2(X, s) \right] ds
    $$

    利用密度函数 $p(x, t)$，将期望写成关于空间 $x$ 的积分形式：

    $$
    \int_{\mathbb{R}} \phi(x) p(x, t) dx = \int_{\mathbb{R}} \phi(x) p(x, 0) dx + \int_0^t \int_{\mathbb{R}} \left( \phi'(x) b(x, s) + \frac{1}{2} \phi''(x) \sigma^2(x, s) \right) p(x, s) dx ds
    $$

    两边对时间 $t$ 求导：

    $$
    \int_{\mathbb{R}} \phi(x) \frac{\partial p}{\partial t} dx = \int_{\mathbb{R}} \left( \phi'(x) b(x, t) + \frac{1}{2} \phi''(x) \sigma^2(x, t) \right) p(x, t) dx
    $$

    **关键所在：两次分部积分 (Integration by Parts)。**
    为了将导数算子从测试函数 $\phi$ 上剥离转移到密度函数 $p$ 上，我们在右侧进行分部积分（假设边界处 $\phi$ 及其导数衰减为 0）：
    
    第一项：

    $$
    \int_{\mathbb{R}} \phi'(x) \big(b(x, t) p(x, t)\big) dx = -\int_{\mathbb{R}} \phi(x) \frac{\partial}{\partial x}\big(b(x, t) p(x, t)\big) dx
    $$

    第二项（连续分部积分两次）：

    $$
    \frac{1}{2} \int_{\mathbb{R}} \phi''(x) \big(\sigma^2(x, t) p(x, t)\big) dx = \frac{1}{2} \int_{\mathbb{R}} \phi(x) \frac{\partial^2}{\partial x^2}\big(\sigma^2(x, t) p(x, t)\big) dx
    $$

    合并后，因为方程对于**任意**的测试函数 $\phi(x)$ 都成立，所以积分解内的函数本身必然相等。这样便得到了**Fokker-Planck 方程**：

    $$
    \frac{\partial p}{\partial t} = - \frac{\partial}{\partial x} \big( b(x, t) p(x, t) \big) + \frac{1}{2} \frac{\partial^2}{\partial x^2} \big( \sigma^2(x, t) p(x, t) \big)
    $$