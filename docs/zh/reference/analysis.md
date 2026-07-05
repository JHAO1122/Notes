# ♾️ 数学分析 (Advanced Mathematical Analysis)

本模块涵盖了多元微积分与函数序列的核心拓扑与分析学定理，重点关注与实变函数、泛函分析高度关联的硬分析（Hard Analysis）工具。

---

## 一、 多元微分学与隐函数定理 (Differential Calculus)

!!! info "多元全微分与连续可微"
    设映射 $f: U \subset \mathbb{R}^n \to \mathbb{R}^m$。若在点 $x_0 \in U$ 处存在线性映射 $L: \mathbb{R}^n \to \mathbb{R}^m$ 满足：

    \[
    \lim_{h \to 0} \frac{\|f(x_0 + h) - f(x_0) - L(h)\|}{\|h\|} = 0
    \]

    则称 $f$ 在 $x_0$ 处**可微**，记 $L = Df(x_0)$（雅可比矩阵）。若 $Df(x)$ 在整个 $U$ 上存在且连续，则称 $f$ 是 **$C^1$ 连续可微**的。

!!! success "隐函数定理 (Implicit Function Theorem)"
    设 $F: \mathbb{R}^n \times \mathbb{R}^m \to \mathbb{R}^m$ 在开集 $U$ 上是 $C^1$ 的。若在某点 $(x_0, y_0) \in U$ 满足：
    
    1. $F(x_0, y_0) = 0$
    
    2. 局部关于 $y$ 的雅可比矩阵可逆，即 $\det \left( \frac{\partial F}{\partial y} (x_0, y_0) \right) \neq 0$
    
    则存在 $x_0$ 的邻域 $V \subset \mathbb{R}^n$ 与唯一的 $C^1$ 映射 $g: V \to \mathbb{R}^m$ 满足 $g(x_0) = y_0$，使得在 $V$ 上恒有：

    \[
    F(x, g(x)) = 0
    \]

    且隐函数的导数（灵敏度）可通过链式法则直接求得：

    \[
    Dg(x) = - \left( \frac{\partial F}{\partial y}(x, g(x)) \right)^{-1} \frac{\partial F}{\partial x}(x, g(x))
    \]

---

## 二、 函数序列一致收敛性与级数 (Uniform Convergence)

!!! abstract "一致收敛的定义"
    定义在区间 $I$ 上的函数列 $\{f_n\}$ **一致收敛**于极限函数 $f$，是指对任意 $\epsilon > 0$，存在一个与自变量 $x$ 无关的整数 $N$，使得当 $n > N$ 时，对**所有** $x \in I$ 恒有：

    \[
    |f_n(x) - f(x)| < \epsilon
    \]

    * **泛函分析视角的等价性**：这等价于 $\{f_n\}$ 在区间 $I$ 的上确界范数（Sup-norm）下收敛：$\lim_{n \to \infty} \|f_n - f\|_\infty = 0$。

!!! success "一致收敛的连续性保持定理（泛函 $C^k[a, b]$ 完备性的基石）"
    若函数列 $\{f_n\}$ 在区间 $[a, b]$ 上一致收敛于 $f$，且每个 $f_n$ 都是 $[a, b]$ 上的**连续函数**，则极限函数 $f$ 在 $[a, b]$ 上也**必然连续**。

??? note "拓展：泛函空间 $C^k[a, b]$ 的完备性联系"
    在泛函分析中，我们定义 Banach 空间 $C^k[a, b]$ 的范数构筑为：

    \[
    \|f\|_{C^k} = \sum_{j=0}^k \max_{t \in [a, b]} |f^{(j)}(t)|
    \]

    要证明该空间完备，即证明其中的 Cauchy 列必收敛在空间内。数学分析的**导数一致收敛定理**为其提供了底层支持：
    
    若各阶导数函数列 $\{f_n^{(j)}\}$ 在 $[a, b]$ 上分别**一致收敛**于 $g_j$，则原函数列的极限 $f$ 必满足 $f^{(j)} = g_j$，从而通过上方的连续性保持定理，确保最终的极限函数 $f \in C^k[a, b]$。

---

## 三、 极限与积分/微分次序的交换 (Interchange of Limits)

!!! info "黎曼积分号下求极限定理"
    设 $\{f_n\}$ 是 $[a, b]$ 上的黎曼可积函数列。若 $f_n \xrightarrow{\text{一致}} f$，则极限函数 $f$ 在 $[a, b]$ 上亦黎曼可积，且可以交换极限与积分号：

    \[
    \lim_{n \to \infty} \int_a^b f_n(x) dx = \int_a^b \left( \lim_{n \to \infty} f_n(x) \right) dx = \int_a^b f(x) dx
    \]

    *(注：在实变函数中，若将黎曼积分升级为 Lebesgue 积分，该定理将被条件更宽松的 Lebesgue 控制收敛定理 (LDCT) 完全取代。)*

!!! success "黎曼积分号下求导定理 (Leibniz 积分法则)"
    设 $f(x, y)$ 及其偏导数 $\frac{\partial f}{\partial y}(x, y)$ 在矩形区域 $[a, b] \times [c, d]$ 上连续。定义变参数积分 $F(y) = \int_a^b f(x, y) dx$，则 $F(y)$ 可导且满足：

    \[
    \frac{d}{dy} \int_a^b f(x, y) dx = \int_a^b \frac{\partial f}{\partial y}(x, y) dx
    \]

---

## 四、 紧性理论与阿泽拉-阿斯科利定理 (Compactness & Arzela-Ascoli)

这是数学分析中最高级、也是与泛函分析（特别是研究紧算子理论）联系最紧密的定理之一，它给出了函数空间中子集紧性的判别准则。

!!! info "一致有界与等度连续 (Equicontinuity)"
    设 $\mathcal{F}$ 是定义在紧区间 $[a, b]$ 上的一个函数族：
    
    * **一致有界 (Uniformly Bounded)**：存在常数 $M > 0$，使得对所有 $f \in \mathcal{F}$ 以及所有 $t \in [a, b]$，恒有 $|f(t)| \le M$。
    
    * **等度连续 (Equicontinuous)**：对任意 $\epsilon > 0$，存在一个**仅与 $\epsilon$ 有关**（而与具体函数 $f$ 及点 $t$ 无关）的 $\delta > 0$，只要 $|t_1 - t_2| < \delta$，就有：

    \[
    |f(t_1) - f(t_2)| < \epsilon \quad (\forall f \in \mathcal{F})
    \]

!!! success "Arzelà-Ascoli (阿泽拉-阿斯科利) 定理"
    设 $\mathcal{F} \subset C[a, b]$。则 $\mathcal{F}$ 是 $C[a, b]$ 中的**列紧集**（即 $\mathcal{F}$ 中的任意函数序列都存在一个一致收敛的子序列）的**充要条件**是：
    $\mathcal{F}$ 在 $[a, b]$ 上是**一致有界**且**等度连续**的。