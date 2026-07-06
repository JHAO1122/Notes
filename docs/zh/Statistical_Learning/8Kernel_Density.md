# 第八章：核密度估计与非参数核回归

在前述章节中，模型往往依赖于对数据潜在分布（如多元正态分布）或功能关系（如线性模型）的参数化假设。当这些假设失效时，非参数估计方法能直接从数据中学习其结构。本章聚焦于基于局部平滑（Locally Smoothing）理论的非参数估计方法，系统探讨**核密度估计（Kernel Density Estimation, KDE）**与**纳达拉亚-沃森核回归（Nadaraya-Watson Kernel Regression）**的数学原理与统计性质。

---

## 1. 核密度估计（Kernel Density Estimation）

!!! info "定义 1.1 (直方图估计量的启发)"

    假设有独立同分布的连续随机变量样本 $X_1, X_2, \dots, X_n \sim f(\cdot)$，其中 $f(x)$ 为待估计的概率密度函数。根据微积分，在点 $x$ 处的密度可表示为小区间内概率的极限：

    \[
    f(x) = \lim_{\lambda \to 0} \frac{1}{\lambda} P\left( x - \frac{\lambda}{2} \le X \le x + \frac{\lambda}{2} \right)
    \]

    传统的直方图估计（Histogram Estimator）将样本空间划分为固定的独立箱子（Bins）。若不划分固定箱子，而是以目标点 $x$ 为中心，构建一个长度为 $\lambda$ 的滑动区间，通过统计落入该区间的样本数，可得到朴素的密度估计量：

    \[
    \hat{f}_{\text{naive}}(x) = \frac{1}{n\lambda} \sum_{i=1}^n \mathbb{I}\left( x - \frac{\lambda}{2} \le X_i \le x + \frac{\lambda}{2} \right) = \frac{1}{n\lambda} \sum_{i=1}^n \mathbb{I}\left( \left| \frac{x - X_i}{\lambda} \right| \le \frac{1}{2} \right)
    \]

### 1.1 核函数的抽象化与平滑

上述朴素估计量对应的权重函数是一个阶跃的方形窗（Uniform Kernel），它在边界处不连续且赋予所有内部点相同的权重。为了获得更平滑的密度估计，我们引入通用的**核函数（Kernel Function）** $K(t)$。

一个合格的核函数必须是一个合法的概率密度函数，即满足以下基本数学性质：

1. **非负性**：$K(t) \ge 0, \quad \forall t \in \mathbb{R}$
2. **积分归一性**：$\int_{-\infty}^{\infty} K(t) dt = 1$
3. **对称性（通常满足）**：$K(-t) = K(t)$，这也意味着其一阶矩为 0：$\int_{-\infty}^{\infty} t K(t) dt = 0$

常用的平滑核函数包括高斯核（Gaussian Kernel）：

\[
K(t) = \frac{1}{\sqrt{2\pi}} \exp\left( -\frac{t^2}{2} \right)
\]

!!! info "定义 1.2 (核密度估计量 KDE)"

    引入窗口宽度（Bandwidth） $h > 0$（对应上述的 $\lambda$），核密度估计量定义为：

    \[
    \hat{f}_h(x) = \frac{1}{nh} \sum_{i=1}^n K\left( \frac{x - X_i}{h} \right)
    \]

---

## 2. KDE 的统计性质：渐近偏置与方差推导

为了评估估计量 $\hat{f}_h(x)$ 的好坏，我们需要对其进行渐近分析。设真实的密度函数 $f(x)$ 具有二阶连续导数。令 $\sigma_K^2 = \int t^2 K(t) dt$，并简写 $R(K) = \int K(t)^2 dt$。

!!! note "定理 2.1 (KDE 的渐近期望、偏置与方差)"

    当 $n \to \infty, h \to 0$ 且 $nh \to \infty$ 时，核密度估计量 $\hat{f}_h(x)$ 的渐近偏置（Bias）与渐近方差（Variance）分别满足：

    \[
    \text{Bias}(\hat{f}_h(x)) = \frac{1}{2} h^2 \sigma_K^2 f''(x) + o(h^2)
    \]

    \[
    \text{Var}(\hat{f}_h(x)) = \frac{1}{nh} f(x) R(K) + o\left( \frac{1}{nh} \right)
    \]

??? proof "证明：KDE 偏置与方差的严格积分推导"

    **1. 期望与偏置推导：**

    利用期望的线性性质以及样本的独立同分布性质：

    \[
    \mathbb{E}[\hat{f}_h(x)] = \mathbb{E}\left[ \frac{1}{nh} \sum_{i=1}^n K\left( \frac{x - X_i}{h} \right) \right] = \frac{1}{h} \mathbb{E}\left[ K\left( \frac{x - X_1}{h} \right) \right]
    \]

    将其写为连续型随机变量的积分形式：

    \[
    \mathbb{E}[\hat{f}_h(x)] = \frac{1}{h} \int_{-\infty}^{\infty} K\left( \frac{x - u}{h} \right) f(u) du
    \]

    引入换元变量 $t = \frac{x - u}{h}$，则 $u = x - ht$，$du = -h dt$。代入积分限转化后得到：

    \[
    \mathbb{E}[\hat{f}_h(x)] = \int_{-\infty}^{\infty} K(t) f(x - ht) dt
    \]

    对真实密度项 $f(x - ht)$ 在点 $x$ 处进行二阶泰勒展开（Taylor Expansion）：

    \[
    f(x - ht) = f(x) - ht f'(x) + \frac{1}{2} h^2 t^2 f''(x) + o(h^2)
    \]

    将展开式代回积分中，并利用核函数的性质（积分归一、对称一阶矩为0）：

    \[
    \mathbb{E}[\hat{f}_h(x)] = f(x) \int K(t) dt - h f'(x) \int t K(t) dt + \frac{1}{2} h^2 f''(x) \int t^2 K(t) dt + o(h^2)
    \]

    \[
    \mathbb{E}[\hat{f}_h(x)] = f(x) \cdot 1 - 0 + \frac{1}{2} h^2 \sigma_K^2 f''(x) + o(h^2)
    \]

    因此，偏置项即为：

    \[
    \text{Bias}(\hat{f}_h(x)) = \mathbb{E}[\hat{f}_h(x)] - f(x) = \frac{1}{2} h^2 \sigma_K^2 f''(x) + o(h^2)
    \]

    **2. 方差推导：**

    由于样本独立，方差具有可加性：

    \[
    \text{Var}(\hat{f}_h(x)) = \text{Var}\left( \frac{1}{nh} \sum_{i=1}^n K\left( \frac{x - X_i}{h} \right) \right) = \frac{1}{n h^2} \text{Var}\left( K\left( \frac{x - X_1}{h} \right) \right)
    \]

    利用方差公式 $\text{Var}(Z) = \mathbb{E}[Z^2] - (\mathbb{E}[Z])^2$：

    \[
    \text{Var}\left( K\left( \frac{x - X_1}{h} \right) \right) = \int K\left( \frac{x - u}{h} \right)^2 f(u) du - \left( \int K\left( \frac{x - u}{h} \right) f(u) du \right)^2
    \]

    对于前项平方积分项，同样采用换元 $t = \frac{x - u}{h}$：

    \[
    \int K\left( \frac{x - u}{h} \right)^2 f(u) du = h \int K(t)^2 f(x - ht) dt
    \]

    对 $f(x - ht)$ 进行一阶泰勒展开 $f(x - ht) = f(x) + o(1)$，代入得：

    \[
    h \int K(t)^2 [f(x) + o(1)] dt = h f(x) R(K) + o(h)
    \]

    而对于后项 $(\mathbb{E}[\dots])^2$ 而言，其量级为 $(h f(x) + o(h))^2 = \mathcal{O}(h^2)$。当 $h \to 0$ 时，相对于前项的 $\mathcal{O}(h)$ 可以忽略不计。将其代回原式：

    \[
    \text{Var}(\hat{f}_h(x)) = \frac{1}{n h^2} \left[ h f(x) R(K) + o(h) \right] = \frac{1}{nh} f(x) R(K) + o\left( \frac{1}{nh} \right)
    \]

    即证完毕。

### 1.2 渐近均方误差与最优窗宽

结合偏置与方差，我们可以计算点 $x$ 处的渐近均方误差（Asymptotic Mean Squared Error, AMSE）：

\[
\text{AMSE}(\hat{f}_h(x)) = \frac{1}{4} h^4 (\sigma_K^2)^2 [f''(x)]^2 + \frac{1}{nh} f(x) R(K)
\]

通过求导 $\frac{d \text{AMSE}}{dh} = 0$，可以求得理论上的最优窗宽（Optimal Bandwidth）满足：

\[
h_{\text{opt}} = \left( \frac{R(K)}{n (\sigma_K^2)^2 \int [f''(x)]^2 dx} \right)^{1/5} \propto n^{-1/5}
  \]

这表明非参数估计的收敛速度为 $n^{-2/5}$，慢于传统参数化模型的 $n^{-1/2}$。

---

## 3. 非参数核回归（Kernel Regression）

现考虑回归模型 $Y = m(X) + \epsilon$，其中目标是在不预设形式的前提下估计条件期望函数（即回归函数）：

\[
m(x) = \mathbb{E}[Y \mid X = x] = \int y f(y \mid x) dy = \frac{\int y f(x, y) dy}{f_X(x)}
\]

### 3.1 纳达拉亚-沃森估计量（Nadaraya-Watson Estimator）

为了估计上述表达式的分子与分母，我们使用二维核函数对联合密度 $f(x, y)$ 进行核密度估计。假定二维核采用乘积形式 $K_h(x) K_b(y)$。

!!! note "定理 3.1 (Nadaraya-Watson 估计量形式)"

    通过对联合密度项和边际密度项分别套用 KDE 并进行积分化简，可得到著名的 N-W 核回归估计量：

    \[
    \hat{m}_{\text{NW}}(x) = \sum_{i=1}^n W_i(x) Y_i
    \]

    其中权重函数（称为 Nadaraya-Watson 权重算子）定义为：

    \[
    W_i(x) = \frac{K\left( \frac{x - X_i}{h} \right)}{\sum_{j=1}^n K\left( \frac{x - X_j}{h} \right)}
    \]

??? proof "证明：Nadaraya-Watson 估计量的数学推导"

    根据条件概率定义，分子项可写为 $g(x) = \int y f(x, y) dy$。我们使用多元 KDE 理论，将观测样本 $(X_i, Y_i)$ 的联合密度估计代入：

    \[
    \hat{f}(x, y) = \frac{1}{n h_x h_y} \sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) K_y\left( \frac{y - Y_i}{h_y} \right)
    \]

    现在对分子 $g(x)$ 进行估计，将 $\hat{f}(x, y)$ 代入积分：

    \[
    \hat{g}(x) = \int y \left[ \frac{1}{n h_x h_y} \sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) K_y\left( \frac{y - Y_i}{h_y} \right) \right] dy
    \]

    利用积分的可加性将求和号与积分号互换位置：

    \[
    \hat{g}(x) = \frac{1}{n h_x} \sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) \left[ \int \frac{y}{h_y} K_y\left( \frac{y - Y_i}{h_y} \right) dy \right]
    \]

    对于括号内部关于 $y$ 的积分，引入换元 $t = \frac{y - Y_i}{h_y}$，则 $y = Y_i + h_y t$，$dy = h_y dt$：

    \[
    \int \frac{Y_i + h_y t}{h_y} K_y(t) h_y dt = \int (Y_i + h_y t) K_y(t) dt = Y_i \int K_y(t) dt + h_y \int t K_y(t) dt
    \]

    由于核函数 $K_y$ 满足积分归一且对称（一阶矩为0），上式化简为：

    \[
    Y_i \cdot 1 + 0 = Y_i
    \]

    将此核心结果代回分子估计式 $\hat{g}(x)$ 中：

    \[
    \hat{g}(x) = \frac{1}{n h_x} \sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) Y_i
    \]

    同时，分母项即为对 $X$ 的边际核密度估计：

    \[
    \hat{f}_X(x) = \frac{1}{n h_x} \sum_{j=1}^n K_x\left( \frac{x - X_j}{h_x} \right)
    \]

    两式相除，消去常数项 $\frac{1}{n h_x}$，即证：

    \[
    \hat{m}_{\text{NW}}(x) = \frac{\hat{g}(x)}{\hat{f}_X(x)} = \frac{\sum_{i=1}^n K_x\left( \frac{x - X_i}{h_x} \right) Y_i}{\sum_{j=1}^n K_x\left( \frac{x - X_j}{h_x} \right)} = \sum_{i=1}^n W_i(x) Y_i
    \]

---

## 4. 局部线性回归（Local Linear Regression）

虽然 N-W 核回归构造直观，但它在处理**边界点**或**设计矩阵非均匀分布**的数据时，会产生严重的系统性偏置。这是因为 N-W 估计实质上是在局部加权计算一个常数项（零阶逼近）。为了修正这一缺陷，我们引入局部线性回归（Local Linear Regression）。

局部线性回归不直接计算局部均值，而是认为函数在目标点 $x_0$ 的邻域内可以由一条直线局部近似：

\[
m(x) \approx m(x_0) + m'(x_0)(x - x_0) = \alpha + \beta (x - x_0)
\]

根据该设定，我们在点 $x_0$ 处求解一个局部加权最小二乘（Locally Weighted Least Squares）问题：

\[
\min_{\alpha, \beta} \sum_{i=1}^n K\left( \frac{X_i - x_0}{h} \right) \left[ Y_i - \alpha - \beta(X_i - x_0) \right]^2
\]

求解出的 $\hat{\alpha}(x_0)$ 即为对该点处回归函数值 $m(x_0)$ 的估计量。通过矩阵求导可以证明，局部线性回归能够完美地抵消边界处的偏置效应（Boundary Bias Cancellation），使其在整体区域内具备更稳健的收敛性能。