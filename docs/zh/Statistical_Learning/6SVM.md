# 第六章：支持向量机（SVM）

支持向量机（Support Vector Machine, SVM）是一种经典的监督学习算法，主要用于二分类任务。其核心思想是在特征空间中寻找一个能够最大化两类样本几何间隔的超平面。本章将从完全线性可分情况下的硬间隔最大化出发，逐步过渡到线性不可分情况下的软间隔最大化，最后引入核技巧以处理非线性分类问题，并给出完整的 Lagrange 对偶问题推导。

---

## 1. 线性可分支持向量机与硬间隔最大化

!!! info "定义 1.1 (分离超平面与分类规则)"

    给定训练数据集 $\mathcal{D}_{n}=\{(x_{i}, y_{i})\}_{i=1}^{n}$，其中 $x_{i} \in \mathbb{R}^p$，标签 $y_{i} \in \{-1, +1\}$。若数据集完全线性可分，则存在一个超平面：

    \[
    f(x) = w^{\top}x + b = 0
    \]

    能够将两类样本完全正确地分开。对应的分类决策规则为 $C(x) = \text{sign}(w^{\top}x + b)$。也就是说，当 $y_i = +1$ 时，有 $w^{\top}x_i + b > 0$；当 $y_i = -1$ 时，有 $w^{\top}x_i + b < 0$。统一写成不等式形式为：

    \[
    y_i (w^{\top}x_i + b) > 0, \quad \forall i = 1, \dots, n
    \]

### 1.1 函数间隔与几何间隔

* **函数间隔（Functional Margin）**：定义为 $\hat{\gamma}_i = y_i (w^{\top}x_i + b)$。通过成比例地放大 $w$ 和 $b$，函数间隔也会按比例放大，因此它不能客观反映点到超平面的物理距离。

* **几何间隔（Geometric Margin）**：对法向量 $w$ 施加 $L_2$ 范数约束，定义点到超平面的真几何距离为：

\[
\gamma_i = \frac{y_i (w^{\top}x_i + b)}{\|w\|_2}
\]

SVM 的目标是最大化整个训练集上的最小几何间隔 $\gamma = \min_{i} \gamma_i$。因此，原始优化问题可写为：

\[
\max_{w, b, \gamma} \gamma
\]

\[
\text{subject to } \frac{y_i (w^{\top}x_i + b)}{\|w\|_2} \ge \gamma, \quad \forall i = 1, \dots, n
\]

为了消去比例缩放的影响，令最小函数间隔 $\hat{\gamma} = 1$。此时几何间隔 $\gamma = \frac{1}{\|w\|_2}$。最大化 $\frac{1}{\|w\|_2}$ 等价于最小化 $\frac{1}{2}\|w\|_2^2$。

!!! note "定理 1.1 (硬间隔 SVM 凸二次规划问题)"

    硬间隔支持向量机的标准原始优化问题表达如下：

    \[
    \min_{w, b} \frac{1}{2} \|w\|_2^2
    \]

    \[
    \text{subject to } y_i (w^{\top}x_i + b) \ge 1, \quad \forall i = 1, \dots, n
    \]

---

## 2. 硬间隔对偶问题的拉格朗日乘子法推导

为了求解上述条件约束优化问题，并为后续引入核函数做准备，我们利用 Lagrange 对偶性（Duality）将原始问题转化为对偶问题。

!!! note "定理 2.1 (硬间隔对偶优化问题)"

    硬间隔 SVM 原始问题的对偶问题是一个关于拉格朗日乘子向量 $\alpha = (\alpha_1, \dots, \alpha_n)^{\top}$ 的最大化问题：

    \[
    \max_{\alpha \in \mathbb{R}^n} \left\{ \sum_{i=1}^n \alpha_i - \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j \right\}
    \]

    \[
    \text{subject to } \sum_{i=1}^n \alpha_i y_i = 0, \quad \alpha_i \ge 0, \quad \forall i = 1, \dots, n
    \]

??? proof "证明：硬间隔 Lagrange 对偶问题的严格推导"

    首先，引入非负的拉格朗日乘子 $\alpha_i \ge 0$，构建拉格朗日函数如下：

    \[
    L(w, b, \alpha) = \frac{1}{2} \|w\|_2^2 - \sum_{i=1}^n \alpha_i \left[ y_i (w^{\top}x_i + b) - 1 \right]
    \]

    根据对偶理论，我们需要先对原始变量 $w$ 和 $b$ 最小化 $L(w, b, \alpha)$。为此，分别求一阶偏导数并令其等于 0：

    1. 关于 $w$ 求导：

    \[
    \frac{\partial L}{\partial w} = w - \sum_{i=1}^n \alpha_i y_i x_i = 0 \implies w = \sum_{i=1}^n \alpha_i y_i x_i
    \]

    2. 关于 $b$ 求导：

    \[
    \frac{\partial L}{\partial b} = -\sum_{i=1}^n \alpha_i y_i = 0 \implies \sum_{i=1}^n \alpha_i y_i = 0
    \]

    将得到的 $w$ 表达式和 $\sum_{i=1}^n \alpha_i y_i = 0$ 约束代回拉格朗日函数中进行消元：

    \[
    L(w, b, \alpha) = \frac{1}{2} \left( \sum_{i=1}^n \alpha_i y_i x_i \right)^{\top} \left( \sum_{j=1}^n \alpha_j y_j x_j \right) - \sum_{i=1}^n \alpha_i y_i \left( \sum_{j=1}^n \alpha_j y_j x_j \right)^{\top} x_i - b \sum_{i=1}^n \alpha_i y_i + \sum_{i=1}^n \alpha_i
    \]

    由于 $b \sum_{i=1}^n \alpha_i y_i = 0$，上式可化简为：

    \[
    L(w, b, \alpha) = \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j - \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j + \sum_{i=1}^n \alpha_i
    \]

    合并同类项，即得：

    \[
    \min_{w, b} L(w, b, \alpha) = \sum_{i=1}^n \alpha_i - \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j
    \]

    最后，对乘子 $\alpha$ 求极大值，结合其自身边界与约束条件，即证对偶目标函数形式。

### 2.2 互补松弛性（Complementary Slackness）与支持向量

根据 KKT（Karush-Kuhn-Tucker）条件，最优解必须满足互补松弛条件：

\[
\alpha_i \left[ y_i (w^{\top}x_i + b) - 1 \right] = 0, \quad \forall i = 1, \dots, n
  \]

* *支持向量的定义*：由上式可知，若 $\alpha_i > 0$，则必定有 $y_i (w^{\top}x_i + b) - 1 = 0$，即这些样本点精确地落在间隔边界上。这类样本被称为**支持向量（Support Vectors）**。最终的超平面参数 $w$ 仅仅由这极少数的 $\alpha_i > 0$ 的支持向量决定，其余大部分 $\alpha_i = 0$ 的样本点对最终决策边界没有任何贡献。

---

## 3. 线性不可分支持向量机与软间隔最大化

在实际数据中，完美的线性可分条件很难满足，经常会有噪声或两类重叠。为了提高模型的容错能力，需要允许部分样本违反间隔约束。

!!! info "定义 3.1 (松弛变量与软间隔)"

    对每个样本引入一个松弛变量（Slack Variable）$\xi_i \ge 0$，将硬约束放宽为：

    \[
    y_i (w^{\top}x_i + b) \ge 1 - \xi_i, \quad \forall i = 1, \dots, n
    \]

    * 当 $0 < \xi_i \le 1$ 时，样本点被正确分类，但落在了间隔内部。
    * 当 $\xi_i > 1$ 时，样本点跨越了决策边界，被完全错分。

!!! note "定理 3.1 (软间隔 SVM 原始优化问题)"

    为了平衡“最大化几何间隔”与“减少错分样本数量”，引入惩罚参数 $C > 0$，软间隔优化目标函数表达为：

    \[
    \min_{w, b, \xi} \left\{ \frac{1}{2} \|w\|_2^2 + C \sum_{i=1}^n \xi_i \right\}
    \]

    \[
    \text{subject to } y_i (w^{\top}x_i + b) \ge 1 - \xi_i, \quad \xi_i \ge 0, \quad \forall i = 1, \dots, n
    \]

### 3.1 软间隔对偶问题与 KKT 条件

利用同样的方法构造 Lagrange 函数并求偏导数，可推导得出软间隔的对偶问题在形式上与硬间隔几乎完全相同，唯一的区别是拉格朗日乘子 $\alpha_i$ 被施加了一个上界 $C$（被称为盒约束 Box Constraint）：

\[
\max_{\alpha \in \mathbb{R}^n} \left\{ \sum_{i=1}^n \alpha_i - \frac{1}{2} \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j y_i y_j x_i^{\top}x_j \right\}
\]

\[
\text{subject to } \sum_{i=1}^n \alpha_i y_i = 0, \quad 0 \le \alpha_i \le C, \quad \forall i = 1, \dots, n
\]

其对应的 KKT 互补松弛条件变形为：

\[
\alpha_i \left[ y_i (w^{\top}x_i + b) - 1 + \xi_i \right] = 0
\]

\[
(C - \alpha_i) \xi_i = 0
\]

* 当 $\alpha_i = 0$ 时，$\xi_i = 0$，样本点被正确分类且落在边界之外。
* 当 $0 < \alpha_i < C$ 时，由第二式可得 $\xi_i = 0$，此时样本点恰好精确落在间隔边界上。
* 当 $\alpha_i = C$ 时，$\xi_i > 0$，样本点落在边界内部或被错误分类。

---

## 4. 非线性支持向量机与核技巧（Kernel Trick）

当样本在原始低维输入空间中呈现高度非线性分布时，我们可以通过一个非线性映射 $\phi(x): \mathbb{R}^p \to \mathcal{H}$ 将其映射到更高维的希尔伯特空间 $\mathcal{H}$ 中，使得数据在高维空间中变得线性可分。

### 4.1 核函数定义与 Mercer 定理

在高维空间中求解对偶问题时，我们需要计算高维向量的内积 $\phi(x_i)^{\top}\phi(x_j)$。如果先映射再求内积，计算复杂度会发生爆炸。为了解决这一问题，模型引入了**核函数（Kernel Function）**。

!!! info "定义 4.1 (核函数)"

    若存在一个从输入空间到高维空间的映射 $\phi(x)$，使得对于所有 $x, z$，满足：

    \[
    K(x, z) = \langle \phi(x), \phi(z) \rangle
    \]

    则称 $K(x, z)$ 为核函数。

常用的核函数包括：

* **多项式核函数（Polynomial Kernel）**：

\[
K(x, z) = (x^{\top}z + 1)^d
\]

* **高斯径向基核函数（RBF Kernel）**：其能够将数据映射到无穷维的希尔伯特空间。

\[
K(x, z) = \exp\left( -\gamma \|x - z\|_2^2 \right)
\]

将核函数代入对偶问题中，所有内积项 $x_i^{\top}x_j$ 均可直接替换为 $K(x_i, x_j)$。非线性支持向量机的最终决策函数变为：

\[
f(x) = \sum_{i=1}^n \alpha_i y_i K(x_i, x) + b
\]

---

## 5. 惩罚损失视角：合页损失函数（Hinge Loss）

除了从几何间隔的角度进行推导外，支持向量机还可以完美地纳入力度正则化经验风险最小化（SRM）的通用机器学习框架中。

!!! note "定理 5.1 (SVM 的 Hinge 损失等价性)"

    软间隔支持向量机的原始优化问题，在数学上等价于以下采用合页损失函数（Hinge Loss）的无约束正则化问题：

    \[
    \min_{w, b} \left\{ \sum_{i=1}^n \max\left( 0, 1 - y_i(w^{\top}x_i + b) \right) + \frac{1}{2C} \|w\|_2^2 \right\}
    \]

??? proof "证明：Hinge 损失等价性的推导"

    令原始变量中的松弛变量满足边界极限情况，根据软间隔的不等式约束：

    \[
    y_i (w^{\top}x_i + b) \ge 1 - \xi_i \implies \xi_i \ge 1 - y_i (w^{\top}x_i + b)
    \]

    同时由于原始问题要求 $\xi_i \ge 0$，为了使目标函数 $\sum_{i=1}^n \xi_i$ 达到最小，$\xi_i$ 应当恰好取到其能够满足约束的最小下界。因此有：

    \[
    \xi_i = \max\left( 0, 1 - y_i (w^{\top}x_i + b) \right)
    \]

    将其直接代入软间隔目标函数中：

    \[
    \min_{w, b} \frac{1}{2} \|w\|_2^2 + C \sum_{i=1}^n \max\left( 0, 1 - y_i (w^{\top}x_i + b) \right)
    \]

    将整个目标函数同时除以常数 $C$，并令新正则化参数 $\lambda = \frac{1}{2C}$，即可得到：

    \[
    \min_{w, b} \sum_{i=1}^n \max\left( 0, 1 - y_i(w^{\top}x_i + b) \right) + \lambda \|w\|_2^2
    \]

    这正是传统的“损失项 + $L_2$ 正则项”框架，即证等价性。