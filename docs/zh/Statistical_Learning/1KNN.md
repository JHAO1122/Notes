# 第一章：K-近邻算法与偏置-方差权衡

在统计学习中，K-近邻算法（K-Nearest Neighbor, KNN）是一种经典且直观的非参数化方法。其核心思想是利用特征空间中距离目标点最近的 $K$ 个已知样本的标签，来对目标点进行预测。本章将详细推导 KNN 在回归与分类任务中的数学表述，并深入探讨机器学习中至关重要的偏置-方差权衡以及自由度等概念。

---

## 1. KNN 回归模型及其数学定义

!!! info "定义 1.1 (回归模型基本设定)"

    设输入变量为 $X \in \mathbb{R}^p$，输出变量为 $Y \in \mathbb{R}$。我们假定现实中的数据生成机制满足以下回归模型：

    \[
    Y = f(X) + \epsilon
    \]

    其中，$f(X)$ 为需要估计的真实回归函数，$\epsilon$ 为随机误差项，满足：

    \[
    \mathbb{E}(\epsilon) = 0, \quad \text{Var}(\epsilon) = \sigma^2
    \]

    且 $\epsilon$ 与 $X$ 独立。现收集到一个包含 $n$ 个独立同分布（i.i.d.）样本的训练数据集：

    \[
    \mathcal{D}_n = \{(x_1, y_1), (x_2, y_2), \dots, (x_n, y_n)\}
    \]

### 1.1 K-近邻回归估计量

对于给定的测试目标点 $x_0$，KNN 算法通过在训练集 $\mathcal{D}_n$ 中寻找与之几何距离最近的 $k$ 个样本来估计该点的回归函数值。

设 $N_k(x_0)$ 为包含这 $k$ 个最近邻样本特征值的集合，则目标点 $x_0$ 处的 KNN 预测值 $\hat{y}$（或记作 $\hat{f}(x_0)$）定义为：

\[
\hat{f}(x_0) = \frac{1}{k} \sum_{x_i \in N_k(x_0)} y_i
\]

---

## 2. 偏置-方差分解与预测误差

为了评估 KNN 模型的泛化能力，我们采用平方误差损失函数（Squared-Error Loss）。在目标点 $x_0$ 处的期望预测误差（Expected Prediction Error, EPE）可以被分解为三个独立的部分。

!!! note "定理 2.1 (偏置-方差分解定理)"

    设 $\hat{f}(x_0)$ 是基于训练集 $\mathcal{D}_n$ 得到的估计量，则在测试点 $x_0$ 处，新观测值 $Y_0 = f(x_0) + \epsilon_0$ 的期望误差具有如下分解式：

    \[
    \text{Err}(x_0) = \mathbb{E}_{\mathcal{D}_n, Y_0} \left[ (Y_0 - \hat{f}(x_0))^2 \right] = \sigma^2 + \left[ \text{Bias}(\hat{f}(x_0)) \right]^2 + \text{Var}(\hat{f}(x_0))
    \]

??? proof "证明：偏置-方差分解的详细推导"

    为了简化记号，将 $\hat{f}(x_0)$ 记为 $\hat{f}$，将 $f(x_0)$ 记为 $f$。首先，将 $Y_0 - \hat{f}$ 进行恒等变形并展开：

    \[
    Y_0 - \hat{f} = (f + \epsilon_0) - \hat{f} = (f - \mathbb{E}[\hat{f}]) + (\mathbb{E}[\hat{f}] - \hat{f}) + \epsilon_0
    \]

    对两边取平方项后求期望 $\mathbb{E}_{\mathcal{D}_n, Y_0}[\cdot]$：

    \[
    \mathbb{E}[(Y_0 - \hat{f})^2] = \mathbb{E}\left[ \left( (f - \mathbb{E}[\hat{f}]) + (\mathbb{E}[\hat{f}] - \hat{f}) + \epsilon_0 \right)^2 \right]
    \]

    将其展开为平方项与交叉项之和：

    \[
    \mathbb{E}[(Y_0 - \hat{f})^2] = \mathbb{E}[(f - \mathbb{E}[\hat{f}])^2] + \mathbb{E}[(\mathbb{E}[\hat{f}] - \hat{f})^2] + \mathbb{E}[\epsilon_0^2] + 2\mathbb{E}[(f - \mathbb{E}[\hat{f}])(\mathbb{E}[\hat{f}] - \hat{f})] + 2\mathbb{E}[(f - \mathbb{E}[\hat{f}])\epsilon_0] + 2\mathbb{E}[(\mathbb{E}[\hat{f}] - \hat{f})\epsilon_0]
    \]

    下面依次分析各个项：

    1. *第一项*：因为 $f - \mathbb{E}[\hat{f}]$ 对于训练集抽样是常数，故：

    \[
    \mathbb{E}[(f - \mathbb{E}[\hat{f}])^2] = (f - \mathbb{E}[\hat{f}])^2 = \left[ \text{Bias}(\hat{f}) \right]^2
    \]

    2. *第二项*：根据方差的定义，显然有：

    \[
    \mathbb{E}[(\mathbb{E}[\hat{f}] - \hat{f})^2] = \text{Var}(\hat{f})
    \]

    3. *第三项*：由于新测试误差 $\epsilon_0$ 满足 $\mathbb{E}[\epsilon_0] = 0$ 且 $\text{Var}(\epsilon_0) = \sigma^2$，故：

    \[
    \mathbb{E}[\epsilon_0^2] = \sigma^2
    \]

    4. *交叉项 1*：因为 $\mathbb{E}[\mathbb{E}[\hat{f}] - \hat{f}] = \mathbb{E}[\hat{f}] - \mathbb{E}[\hat{f}] = 0$，且 $(f - \mathbb{E}[\hat{f}])$ 为常数，故：

    \[
    \mathbb{E}[(f - \mathbb{E}[\hat{f}])(\mathbb{E}[\hat{f}] - \hat{f})] = (f - \mathbb{E}[\hat{f}]) \cdot \mathbb{E}[\mathbb{E}[\hat{f}] - \hat{f}] = 0
    \]

    5. *交叉项 2* 和 *交叉项 3*：由于测试集误差 $\epsilon_0$ 与训练集 $\mathcal{D}_n$ 相互独立，且 $\mathbb{E}[\epsilon_0] = 0$，因此：

    \[
    \mathbb{E}[(f - \mathbb{E}[\hat{f}])\epsilon_0] = (f - \mathbb{E}[\hat{f}]) \cdot \mathbb{E}[\epsilon_0] = 0
    \]

    \[
    \mathbb{E}[(\mathbb{E}[\hat{f}] - \hat{f})\epsilon_0] = \mathbb{E}[\mathbb{E}[\hat{f}] - \hat{f}] \cdot \mathbb{E}[\epsilon_0] = 0
    \]

    综上所述，所有交叉项均为 0。将各非零项相加，即证：

    \[
    \text{Err}(x_0) = \sigma^2 + \left[ \text{Bias}(\hat{f}(x_0)) \right]^2 + \text{Var}(\hat{f}(x_0))
    \]

### 2.2 三种误差的统计学含义

* **不可约误差 (Irreducible Error)**：即 $\sigma^2$。这是系统固有的随机噪声，无法通过改进算法来消除。

* **偏置的平方 (Bias$^2$)**：反映了模型的期望预测与真实函数 $f(x_0)$ 之间的系统性偏差。

* **方差 (Variance)**：反映了模型在不同训练集重复抽样下预测值的波动剧烈程度。

### 2.3 KNN 中的偏置-方差渐近性质

在 KNN 算法中，超参数 $k$ 决定了模型的复杂度：

* 当 $k=1$（1-近邻）时，由于估计量强行拟合最近的一个训练样本，随着 $n \to \infty$，最近邻点收敛到 $x_0$，使得：

\[
\text{Bias}^2(\hat{f}(x_0)) \to 0
\]

然而，由于只使用了一个观测值，其预测方差在极限下大约为：

\[
\text{Var}(\hat{f}(x_0)) = \text{Var}(y_i) = \sigma^2
\]

* 当 $k$ 增大时，由于引入了较远特征点的响应值，*模型的偏置增大*（产生欠拟合）；但由于对 $k$ 个样本取平均，*预测方差降低*，近似为：

\[
\text{Var}(\hat{f}(x_0)) \approx \frac{\sigma^2}{k}
\]

---

## 3. 模型复杂度与有效自由度

在广义加性或非参数模型中，模型的复杂度不便用简单的参数个数来衡量。为此，我们引入通用自由度的定义。

!!! info "定义 3.1 (有效自由度 Degrees of Freedom)"

    假设特征值 $X_i = x_i$ 是固定的常数（非随机）。设真实观测向量为 $Y = (Y_1, \dots, Y_n)^T$，模型拟合向量为 $\hat{Y} = (\hat{Y}_1, \dots, \hat{Y}_n)^T$。则该模型的有效自由度定义为：

    \[
    \text{df}(\hat{f}) = \frac{1}{\sigma^2} \sum_{i=1}^n \text{Cov}(\hat{Y}_i, Y_i) = \frac{1}{\sigma^2} \text{Trace}\left( \text{Cov}(\hat{Y}, Y) \right)
    \]

### 3.1 典型模型的自由度验证

* **1-近邻模型**：每个点的预测值就是其自身对应的观测值（即 $\hat{Y}_i = Y_i$），因此：

\[
\text{Cov}(\hat{Y}_i, Y_i) = \text{Cov}(Y_i, Y_i) = \sigma^2
\]

代入定义式可得自由度：

\[
\text{df} = \frac{1}{\sigma^2} \sum_{i=1}^n \sigma^2 = n
\]

* **$n$-近邻模型**：每个点的预测值均为全样本的均值（即 $\hat{Y}_i = \frac{1}{n}\sum_{j=1}^n Y_j$），由此可得：

\[
\text{Cov}(\hat{Y}_i, Y_i) = \text{Cov}\left( \frac{1}{n} \sum_{j=1}^n Y_j, Y_i \right) = \frac{1}{n} \sigma^2
\]

代入定义式求和可得自由度：

\[
\text{df} = \frac{1}{\sigma^2} \sum_{i=1}^n \left( \frac{1}{n}\sigma^2 \right) = 1
\]

* **$k$-近邻模型**：同理，每个样本对应的自由度贡献为 $\frac{1}{k}\sigma^2$，其整体有效自由度为：

\[
\text{df} = \frac{n}{k}
\]

---

## 4. KNN 分类模型与 1NN 误差上界

对于分类任务，KNN 采用多数投票制（Majority Vote）。

!!! info "定义 4.1 (KNN 分类决策)"

    设分类目标为 $c \in \mathcal{C}$。给定测试点 $x_0$，其邻域 $N_k(x_0)$ 内包含 $k$ 个样本，则 KNN 的硬分类预测结果为：

    \[
    \hat{y} = \arg\max_{c \in \mathcal{C}} \sum_{x_i \in N_k(x_0)} \mathbb{I}(y_i = c)
    \]

    其中，$\mathbb{I}(\cdot)$ 为指示函数。

### 4.1 1NN 渐近误差率上界证明

在二分类问题中（假设标签 $Y \in \{0, 1\}$），当样本量 $n \to \infty$ 时，1NN 的渐近误差率存在一个经典的、由 Bayes 错误率控制的上界。

!!! note "定理 4.1 (Cover-Hart 1NN 误差上界定理)"

    设 $P(Y=1 \mid x)$ 是定义在特征空间上的真实条件概率。记 $x_0$ 处的最佳 Bayes 错误率（最小可能错误率）为：

    \[
    P^* = \min\{P(Y=1 \mid x_0), 1 - P(Y=1 \mid x_0)\}
    \]

    若概率密度和条件概率平滑，当 $n \to \infty$ 时，1-近邻算法的期望错误率 $P_{1NN}$ 满足：

    \[
    P_{1NN} \le 2P^*(1 - P^*) \le 2P^*
    \]

??? proof "证明：1NN 渐近误差上界的推导"

    设 $x_{1NN}$ 为在训练集中距离目标点 $x_0$ 最近的样本点。当 $n \to \infty$ 时，由特征空间的稠密性，最近邻点在几何上逼近目标点：

    \[
    d(x_0, x_{1NN}) \to 0
    \]

    由条件概率函数的连续性假设，有：

    \[
    P(Y=1 \mid x_{1NN}) \to P(Y=1 \mid x_0)
    \]

    为了简写，记 $p = P(Y=1 \mid x_0)$。在测试点 $x_0$ 处，1NN 模型发生错分的情况有两种：一是真实值为 1 但邻域预测为 0，二是真实值为 0 但邻域预测为 1。因此，当 $n \to \infty$ 时，1NN 的渐近条件错误率为：

    \[
    P_{1NN} = P(Y_0=1 \mid x_0)P(Y_{1NN}=0 \mid x_{1NN}) + P(Y_0=0 \mid x_0)P(Y_{1NN}=1 \mid x_{1NN})
    \]

    代入渐近极限值，可得：

    \[
    P_{1NN} = p(1 - p) + (1 - p)p = 2p(1 - p)
    \]

    下面引入最佳 Bayes 错误率 $P^* = \min\{p, 1-p\}$。根据此定义，显然有：

    \[
    p(1 - p) = P^*(1 - P^*)
    \]

    因此，1NN 的错误率可以精确表示为：

    \[
    P_{1NN} = 2P^*(1 - P^*)
    \]

    由于 $P^* \in [0, 0.5]$，从而 $(1 - P^*) \le 1$。由此立刻得出上界：

    \[
    P_{1NN} \le 2P^*
    \]

    即证，极限情况下 1NN 的错误率不会超过最佳 Bayes 错误率的两倍。

---

## 5. 距离度量与维数灾难

### 5.1 常用距离度量数学表达

为了界定“近邻”，通常需要计算两点间的距离。设 $u, v \in \mathbb{R}^p$：

* **欧氏距离 (Minkowski $L_2$ 范数)**：

\[
d_2(u, v) = \sqrt{\sum_{j=1}^p (u_j - v_j)^2}
\]

* **马氏距离 (Mahalanobis Distance)**：能够消除尺度影响并考虑特征间相关性，设 $\Sigma$ 均为特征协方差矩阵：

\[
d_{\Sigma}(u, v) = \sqrt{(u - v)^T \Sigma^{-1} (u - v)}
\]

### 5.2 维数灾难的几何数学推导

维数灾难（Curse of Dimensionality）是指在高维空间中，数据的稀疏性会呈指数级增长。

!!! note "高维立方体邻域的边长推导"

    假设 $n$ 个样本点均匀分布在 $p$ 维单位超立方体 $[0, 1]^p$ 中。为了捕获固定比例 $r = \frac{k}{n}$ 的本地样本，我们需要构造一个体积同样占总体积比例为 $r$ 的子超立方体。

    设该子超立方体的目标边长为 $l$，则其体积为 $l^p$。由于数据均匀分布，根据体积公式有：

    \[
    l^p = r \implies l = r^{\frac{1}{p}}
    \]

    我们代入具体的数值来审视不同维度 $p$ 下的边长 $l$ 变化（设 $r = 0.1$，即试图捕获 $10\%$ 的数据）：

    1. 当 $p = 1$（一维线段）时：

    \[
    l = 0.1^1 = 0.10
    \]

    2. 当 $p = 2$（二维平面）时：

    \[
    l = 0.1^{0.5} \approx 0.32
    \]

    3. 当 $p = 10$（十维空间）时：

    \[
    l = 0.1^{0.1} \approx 0.63
    \]

    4. 当 $p = 100$（百维空间）时：

    \[
    l = 0.1^{0.01} \approx 0.955
    \]

    这表明，在 100 维空间中，为了获取最近邻的 $10\%$ 样本，必须在每个坐标轴上都覆盖多达 $95.5\%$ 的跨度。此时，所谓的“局部近邻”在几何上已经不再具有“局部”的亲近性，导致基于距离邻域的估计彻底失效。