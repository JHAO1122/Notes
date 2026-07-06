# 第二章：线性回归模型及其统计推断

线性回归（Linear Regression）是统计学习中最基础、应用最广泛的参数化方法之一。它假设输入变量与输出变量之间存在线性关系。本章将从矩阵形式的最小二乘估计出发，推导其统计性质，并深入探讨高维设定下的模型选择方法。

---

## 1. 线性回归模型的基本设定与矩阵表示

!!! info "定义 1.1 (线性回归模型)"

    观察到一个包含 $n$ 个独立同分布（i.i.d.）样本的训练数据集：

    \[
    \mathcal{D}_{n} = \{(x_{i}, y_{i})\}_{i=1}^{n}
    \]

    其中每个输入向量 $x_{i}$ 是一个 $p$ 维的特征向量：

    \[
    x_{i} = (x_{i1}, x_{i2}, \dots, x_{ip})^{\top}
    \]

    响应变量 $y_{i} \in \mathbb{R}$ 为连续型变量。我们假设潜在的真实数据生成机制满足以下多元线性模型：

    \[
    Y = \beta_0 + \beta_1 X_1 + \beta_2 X_2 + \dots + \beta_p X_p + \epsilon
    \]

    其中，$\beta_0$ 为截距项，$\beta_j$ ($j=1,\dots,p$) 为回归系数，$\epsilon$ 为随机误差项，满足经典的 Gauss-Markov 假设：

    \[
    \mathbb{E}(\epsilon) = 0, \quad \text{Var}(\epsilon) = \sigma^2
    \]

    且各个样本的误差项之间相互独立。

### 1.1 数据的矩阵表示

为了便于数学推导，通常将截距项吸收到特征向量中。我们令首列特征全为 1，即 $x_{i} = (1, x_{i1}, \dots, x_{ip})^{\top}$。定义设计矩阵（Design Matrix）$X \in \mathbb{R}^{n \times (p+1)}$ 以及观测向量 $Y \in \mathbb{R}^n$ 如下：

\[
X = \begin{pmatrix} 1 & x_{11} & \dots & x_{1p} \\ 1 & x_{21} & \dots & x_{2p} \\ \vdots & \vdots & \ddots & \vdots \\ 1 & x_{n1} & \dots & x_{np} \end{pmatrix}, \quad Y = \begin{pmatrix} y_1 \\ y_2 \\ \vdots \\ y_n \end{pmatrix}
\]

同时令回归系数向量 $\beta = (\beta_0, \beta_1, \dots, \beta_p)^{\top}$，则整个样本系统的多元线性回归模型可以紧凑地表示为矩阵形式：

\[
Y = X\beta + \epsilon
\]

---

## 2. 最小二乘估计（OLS）

为了估计未知参数 $\beta$，我们通常采用最小二乘法（Ordinary Least Squares, OLS），其核心思想是最小化残差平方和（Residual Sum of Squares, RSS）。

!!! note "定理 2.1 (OLS 估计量的显式解)"

    若设计矩阵 $X$ 为满列秩（即 $\text{Rank}(X) = p+1 < n$），则使得残差平方和最小的参数估计量 $\hat{\beta}$ 具有如下唯一的闭式解：

    \[
    \hat{\beta} = (X^{\top}X)^{-1}X^{\top}Y
    \]

??? proof "证明：基于矩阵微积分的 OLS 显式解推导"

    首先，根据定义写出残差平方和 $RSS(\beta)$ 的矩阵表达式：

    \[
    RSS(\beta) = \sum_{i=1}^n (y_i - x_i^{\top}\beta)^2 = (Y - X\beta)^{\top}(Y - X\beta)
    \]

    利用矩阵乘法将其展开：

    \[
    RSS(\beta) = Y^{\top}Y - Y^{\top}X\beta - \beta^{\top}X^{\top}Y + \beta^{\top}X^{\top}X\beta
    \]

    由于 $Y^{\top}X\beta$ 是一个标量，它的转置等于其自身，即 $(Y^{\top}X\beta)^{\top} = \beta^{\top}X^{\top}Y$。因此上式可以简化为：

    \[
    RSS(\beta) = Y^{\top}Y - 2\beta^{\top}X^{\top}Y + \beta^{\top}X^{\top}X\beta
    \]

    为了求得使 $RSS(\beta)$ 最小的 $\beta$，我们需要对向量 $\beta$ 求偏导数。根据矩阵微分公式 $\frac{\partial (\beta^{\top}A)}{\partial \beta} = A$ 以及 $\frac{\partial (\beta^{\top}B\beta)}{\partial \beta} = 2B\beta$（其中 $B$ 为对称矩阵），可得：

    \[
    \frac{\partial RSS(\beta)}{\partial \beta} = -2X^{\top}Y + 2X^{\top}X\beta
    \]

    令导数等于 0，得到著名的正规方程组（Normal Equations）：

    \[
    X^{\top}X\beta = X^{\top}Y
    \]

    由于假设 $X$ 为满列秩，矩阵 $X^{\top}X$ 必定是可逆的。在正规方程组两边同时左乘 $(X^{\top}X)^{-1}$，即证：

    \[
    \hat{\beta} = (X^{\top}X)^{-1}X^{\top}Y
    \]

---

## 3. OLS 估计量的统计性质

求得参数估计量 $\hat{\beta}$ 后，我们需要定量分析其作为随机向量的统计特性（均值与方差-协方差矩阵）。

!!! note "定理 3.1 (OLS 估计量的无偏性与方差)"

    在 Gauss-Markov 假设下，最小二乘估计量 $\hat{\beta}$ 是真实参数 $\beta$ 的无偏估计，且其协方差矩阵为：

    \[
    \text{Var}(\hat{\beta}) = \sigma^2 (X^{\top}X)^{-1}
    \]

??? proof "证明：无偏性与方差矩阵的推导"

    **1. 无偏性证明：**

    将真实模型 $Y = X\beta + \epsilon$ 代入 $\hat{\beta}$ 的表达式中：

    \[
    \hat{\beta} = (X^{\top}X)^{-1}X^{\top}(X\beta + \epsilon) = (X^{\top}X)^{-1}X^{\top}X\beta + (X^{\top}X)^{-1}X^{\top}\epsilon
    \]

    由于 $(X^{\top}X)^{-1}X^{\top}X = I$（单位矩阵），因此有：

    \[
    \hat{\beta} = \beta + (X^{\top}X)^{-1}X^{\top}\epsilon
    \]

    在给定设计矩阵 $X$（通常视为非随机常数）的条件下，对 $\hat{\beta}$ 两边取期望：

    \[
    \mathbb{E}(\hat{\beta}) = \mathbb{E}\left[ \beta + (X^{\top}X)^{-1}X^{\top}\epsilon \right] = \beta + (X^{\top}X)^{-1}X^{\top}\mathbb{E}(\epsilon)
    \]

    因为 $\mathbb{E}(\epsilon) = 0$，所以：

    \[
    \mathbb{E}(\hat{\beta}) = \beta
    \]

    无偏性得证。

    **2. 协方差矩阵证明：**

    根据线性变换的方差性质，对于常数矩阵 $A$ 和随机向量 $Z$，有 $\text{Var}(AZ) = A\text{Var}(Z)A^{\top}$。由此我们计算 $\hat{\beta}$ 的方差：

    \[
    \text{Var}(\hat{\beta}) = \text{Var}\left( \beta + (X^{\top}X)^{-1}X^{\top}\epsilon \right) = \text{Var}\left( (X^{\top}X)^{-1}X^{\top}\epsilon \right)
    \]

    将 $(X^{\top}X)^{-1}X^{\top}$ 视为常数矩阵 $A$，$\epsilon$ 视为随机向量 $Z$。因为各个样本误差独立同分布，其协方差矩阵为 $\text{Var}(\epsilon) = \sigma^2 I_n$。因此：

    \[
    \text{Var}(\hat{\beta}) = \left[ (X^{\top}X)^{-1}X^{\top} \right] \cdot \text{Var}(\epsilon) \cdot \left[ (X^{\top}X)^{-1}X^{\top} \right]^{\top}
    \]

    \[
    \text{Var}(\hat{\beta}) = (X^{\top}X)^{-1}X^{\top} \cdot (\sigma^2 I_n) \cdot X \left[ (X^{\top}X)^{-1} \right]^{\top}
    \]

    由于 $(X^{\top}X)$ 是对称矩阵，其逆矩阵也是对称的，即 $\left[ (X^{\top}X)^{-1} \right]^{\top} = (X^{\top}X)^{-1}$。化简上式：

    \[
    \text{Var}(\hat{\beta}) = \sigma^2 (X^{\top}X)^{-1}X^{\top}X(X^{\top}X)^{-1} = \sigma^2 (X^{\top}X)^{-1}
    \]

    方差矩阵得证。

### 3.1 总体方差 $\sigma^2$ 的无偏估计

由于真实的误差项方差 $\sigma^2$ 是未知的，我们常用残差平方和进行估计。定义 $n$ 维拟合向量为 $\hat{Y} = X\hat{\beta}$。则 $\sigma^2$ 的无偏估计量 $s^2$（或记为 $\hat{\sigma}^2$）定义为：

\[
s^2 = \frac{RSS(\hat{\beta})}{n - p - 1} = \frac{(Y - \hat{Y})^{\top}(Y - \hat{Y})}{n - p - 1}
\]

分母中的 $n - p - 1$ 代表自由度，修正了因为估计 $p+1$ 个回归系数而造成的自由度损失。

---

## 4. 训练误差、测试误差与模型选择准则

在机器学习与统计学中，我们必须区分模型在训练集上的表现与未见过的测试集上的表现。

* **训练误差（Training Error）**：模型在训练数据集 $\mathcal{D}_n$ 上计算得到的平均损失，随着模型复杂度的增加（即引入更多特征），训练误差总是单调递减。

* **测试误差（Testing Error）**：模型在独立抽样的全新测试点上的期望损失，呈现出先降后升的 U 型曲线。

当特征维度 $p$ 很大时，全模型容易导致严重的过拟合。为了挑选出泛化能力最佳的子模型，我们需要利用特定的统计学准则在模型复杂度和拟合优度之间做出权衡。

### 4.1 Mallows' $C_p$ 准则

对于包含 $d$ 个特征的子模型，其 $C_p$ 统计量定义为：

\[
C_p = \frac{1}{n} \left( RSS_d + 2d\hat{\sigma}^2 \right)
\]

其中 $RSS_d$ 是该包含 $d$ 个特征模型的残差平方和，$\hat{\sigma}^2$ 是包含所有特征的全模型对总体方差的估计。

### 4.2 AIC（赤池信息准则）与 BIC（贝叶斯信息准则）

在随机误差满足正态分布的假设下，基于最大似然估计的模型选择准则表达如下：

* **AIC 准则**：

\[
AIC = \frac{1}{n\hat{\sigma}^2} \left( RSS_d + 2d\hat{\sigma}^2 \right)
\]

* **BIC 准则**：BIC 相比 AIC 对高复杂度模型施加了更重的惩罚（当 $\ln n > 2$ 时）：

\[
BIC = \frac{1}{n\hat{\sigma}^2} \left( RSS_d + \ln(n)d\hat{\sigma}^2 \right)
\]

我们应当选择使得 $C_p$、AIC 或 BIC 达到最小值的特征子集。

---

## 5. 模型选择算法

当特征总数 $p$ 较大时，遍历所有 $2^p$ 个可能的特征组合在计算上是不现实的。因此我们需要采用计算高效的搜索算法。

### 5.1 最优子集选择（Best Subset Selection）

* **基本步骤**：对于每一个固定的模型大小 $k = 1, 2, \dots, p$，穷举所有 $\binom{p}{k}$ 种特征组合，找出其中令 $RSS$ 最小的那个模型。接着，利用 $C_p$、AIC 或 BIC 准则在这些最优的交叉尺寸模型中挑选出最终的全局最佳模型。

* **局限性**：由于需要检查 $2^p$ 个模型，该算法仅在 $p$ 较小（通常 $p < 50$）时可行。

### 5.2 逐步回归算法（Stepwise Regression）

逐步回归是一类贪心算法，虽然无法保证找到全局最优解，但在实际应用中速度极快且效果优良：

* *前向逐步选择（Forward Stepwise Selection）*：从不包含任何特征的零模型开始，每一步都测试将剩余特征逐个加入模型后的效果，选择让 $RSS$ 下降最显著的那个特征加入，直至满足停止准则。

* *后向逐步选择（Backward Elimination）*：从包含所有 $p$ 个特征的全模型开始，每一步都尝试剔除一个特征，选择使 $RSS$ 增加最小（即最不显著）的那个特征并将其移除，直至所有剩余特征都达到显著性要求。