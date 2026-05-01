# 第五章：赋范线性空间基础

赋范线性空间是泛函分析研究的基本对象。它不仅具有线性代数中的代数结构（线性运算），还通过引入“范数”建立了拓扑结构（距离与收敛）。本章将从线性空间出发，逐步过渡到赋范线性空间，并介绍几类经典的函数空间。

---

## 1. 线性空间 (Linear Space)

!!! info "定义 5.1 (线性空间)"

    设 $E$ 是一个非空集合，$K$ 是实数域 $\mathbb{R}$ 或复数域 $\mathbb{C}$。若 $E$ 满足以下性质，则称 $E$ 为线性空间：

    **(i) 加法群性质**：对任意 $x, y, z \in E$，对应唯一的一个和 $x+y \in E$，满足：

    * (a) 交换律：$x+y = y+x$；

    * (b) 结合律：$(x+y)+z = x+(y+z)$；

    * (c) 存在零元素 $\theta$，使得对任一 $x \in E$，有 $x+\theta = x$；

    * (d) 任何 $x \in E$，存在加法逆元素 $-x$ 使得 $x+(-x) = \theta$。

    **(ii) 数乘性质**：对任意 $x, y \in E$ 以及任何数 $\alpha, \beta \in K$，对应唯一的一个积 $\alpha x \in E$，满足：

    * (a) $\alpha(\beta x) = (\alpha\beta)x$；

    * (b) $1 \cdot x = x$；

    * (c) $(\alpha+\beta)x = \alpha x + \beta x$；

    * (d) $\alpha(x+y) = \alpha x + \alpha y$。

线性空间中的元素又称为**向量**，元素的相加以及数与元素的相乘统称为**线性运算**。例如 $C[0,1]$ 和 $L^2[0,1]$ 都是线性空间。

---

## 2. 线性空间的性质

由线性空间的定义可导出下列性质：

* (i) $0 \cdot x = \theta$。

* (ii) $(-1)x = -x$。

* (iii) $\alpha \cdot \theta = \theta$。

* (iv) $\alpha x = \theta$ 的充分必要条件是 $\alpha = 0$ 或 $x = \theta$。

??? proof "证明细节"

    **关于 (i) 的证明：**
    因为 $2(0 \cdot x) = (2 \cdot 0)x = 0 \cdot x$，两边减去 $0 \cdot x$ 得：

    $$
    0 \cdot x = 2(0 \cdot x) - 0 \cdot x = \theta
    $$

    **关于 (ii) 的证明：**
    由于 $(-1)x + x = (-1 + 1)x = 0 \cdot x = \theta$，根据逆元素的唯一性知 $(-1)x = -x$。

    **关于 (iii) 的证明：**

    $$
    \alpha \cdot \theta = \alpha(\theta + (-\theta)) = \alpha \theta + \alpha(-\theta) = \alpha \theta - \alpha \theta = \theta
    $$

    **关于 (iv) 的证明：**
    充分性显然。必要性：若 $\alpha x = \theta$ 且 $\alpha \neq 0$，则：

    $$
    x = 1 \cdot x = \left( \frac{1}{\alpha} \cdot \alpha \right) x = \frac{1}{\alpha} (\alpha x) = \frac{1}{\alpha} \theta = \theta
    $$

---

## 3. 子空间、张成与同构

### 3.1 子空间 (Subspace)

!!! info "定义 5.2 (子空间)"

    设 $E_0$ 是线性空间 $E$ 的一个子集。若 $E_0$ 对 $E$ 中的线性运算封闭（即对任意 $x, y \in E_0$ 及 $\alpha \in K$，有 $x+y \in E_0$ 且 $\alpha x \in E_0$），则称 $E_0$ 为 $E$ 的**线性子空间**。

* 不同于 $E$ 本身的子空间称为 $E$ 的**真子空间**。

* 线性无关：若集合 $A \subset E$ 中任意有限个元素均线性无关，则称 $A$ 线性无关。

### 3.2 张成的子空间 span L

!!! info "定义 5.3 (span L)"

    设 $L$ 是线性空间 $E$ 的一个非空子集，由 $L$ 中元素的所有可能有限线性组合构成的集合称为由 $L$ **张成的子空间**，记为 $\text{span } L$：

    $$
    \text{span } L = \left\{ \sum_{k=1}^{n} c_k x_k \mid c_k \in K, x_k \in L, n \in \mathbb{N} \right\}
    $$

* $\text{span } L$ 是包含 $L$ 的**最小子空间**，也是所有包含 $L$ 的子空间的交。

* 例如在 $L^2[-\pi, \pi]$ 中，子集 $L = \{ \frac{1}{\sqrt{\pi}} \sin kt, \frac{1}{\sqrt{\pi}} \cos kt \}_{k=0}^{\infty}$ 张成的子空间即为全体三角多项式。

### 3.3 线性空间的同构

* 设 $E, E'$ 都是线性空间。若存在一个从 $E$ 到 $E'$ 的双射 $T$，满足 $T(x+y) = Tx + Ty$ 且 $T(ax) = aTx$，则称 $E$ 与 $E'$ 同构。

---

## 4. 直接和 (Direct Sum)

!!! info "定义 5.4 (直接和)"

    设 $L_1, \dots, L_n$ 是 $E$ 的子空间。如果任一元素 $x \in E$ 可以**唯一**地表示成 $x = x_1 + \dots + x_n$ ($x_k \in L_k$)，则称 $E$ 是 $L_1, \dots, L_n$ 的直接和，记为 $E = L_1 \oplus \dots \oplus L_n$。

* 性质：若 $E$ 是 $L_1, \dots, L_n$ 的直接和，则从中任取的非零元素 $x_1, \dots, x_n$ 必线性无关。

* 构造：考虑有序组 $x = (x_1, \dots, x_n)$ 构成的集，定义分量加法和数乘，则该集构成 $L_k$ 的直接和。

---

## 5. 赋范线性空间 (Normed Linear Space)

在线性空间中引入拓扑，必须与线性运算结合。

!!! info "定义 5.5 (范数)"

    设 $E$ 是线性空间。如果对于 $E$ 中每个元素 $x$ 都有一个实数 $\|x\|$ 与之对应，且满足：

    * (i) 非负性：$\|x\| \ge 0$，且 $\|x\| = 0 \iff x = \theta$；

    * (ii) 齐次性：$\|\alpha x\| = |\alpha| \|x\|$；

    * (iii) 三角不等式：$\|x+y\| \le \|x\| + \|y\|$。

    则称 $E$ 为赋范线性空间，$\|x\|$ 称为范数。

### 5.1 距离与强收敛

赋范线性空间自然诱导距离 $\rho(x, y) = \|x - y\|$。

* 依范数收敛（强收敛）：点列 $\{x_n\}$ 收敛于 $x$ 是指 $\|x_n - x\| \rightarrow 0$。

### 5.2 线性运算的连续性

* (i) 范数 $\|x\|$ 是连续泛函：$| \|x_n\| - \|x\| | \le \|x_n - x\|$。若 $x_n$ 依范数收敛，则 $\{x_n\}$ 有界。

* (ii) 加法连续性：若 $x_n \to x, y_n \to y$，则 $x_n + y_n \to x+y$。

    $$
    \|x_n + y_n - (x+y)\| \le \|x_n - x\| + \|y_n - y\| \rightarrow 0
    $$

* (iii) 数乘连续性：若 $\alpha_n \to \alpha, x_n \to x$，则 $\alpha_n x_n \to \alpha x$。

??? proof "证明：数乘连续性"

    $$
    \|\alpha_n x_n - \alpha x\| = \|\alpha_n x_n - \alpha_n x + \alpha_n x - \alpha x\|
    $$

    $$
    \le |\alpha_n| \|x_n - x\| + |\alpha_n - \alpha| \|x\| \rightarrow 0
    $$

---

## 6. 经典赋范线性空间实例


* **1. $\mathbb{R}^n$ 与 $\mathbb{C}^n$**：
    对于 $x = (\xi_1, \dots, \xi_n)$，定义欧氏范数：

    $$
    \|x\| = \sqrt{\sum_{j=1}^n |\xi_j|^2}
    $$

* **2. 连续函数空间 $C[a, b]$**：
    定义极大值范数：

    $$
    \|x\| = \max_{t \in [a, b]} |x(t)|
    $$

* **3. 序列空间 $l^p$ ($1 \le p < \infty$)**：
    由满足 $\sum_{j=1}^\infty |\xi_j|^p < \infty$ 的数列构成，范数为：

    $$
    \|x\|_p = \left( \sum_{j=1}^\infty |\xi_j|^p \right)^{\frac{1}{p}}
    $$

    其三角不等式 $\|x+y\|_p \le \|x\|_p + \|y\|_p$ 即为 Minkowski 不等式。

* **4. 可积函数空间 $L^p(F)$ ($1 \le p < \infty$)**：
    设 $F \subset \mathbb{R}$ 为可测集，范数定义为：

    $$
    \|f\|_p = \left( \int_F |f(t)|^p dt \right)^{\frac{1}{p}}
    $$

* **5. 连续导数空间 $C^k[a, b]$**：
    具有直到 $k$ 阶连续导数的函数全体，范数定义为：

    $$
    \|x\| = \sum_{j=0}^{k} \max_{t \in [a, b]} |x^{(j)}(t)|
    $$