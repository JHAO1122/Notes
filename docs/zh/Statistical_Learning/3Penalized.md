# 第三章：惩罚线性回归与收缩估计

当特征维度 $p$ 较大、甚至大于样本量 $n$（即 $p > n$ 高维设定）时，传统的最小二乘估计（OLS）会面临严重的过拟合问题，甚至由于设计矩阵 $X^{\top}X$ 不可逆而导致解不唯一。本章引入惩罚线性回归（Penalized Linear Regression），通过在损失函数中加入对系数的惩罚项（正则化），利用偏置-方差权衡来提高模型的预测精度与可解释性。

---

## 1. 脊回归（Ridge Regression）

!!! info "定义 1.1 (脊回归优化问题)"

    脊回归（岭回归）在残差平方和（RSS）的基础上引入了回归系数的 $L_2$ 范数惩罚项。其目标函数定义为：

    \[
    \min_{\beta \in \mathbb{R}^p} \left\{ \sum_{i=1}^n \left( y_i - \beta_0 - \sum_{j=1}^p \beta_j x_{ij} \right)^2 + \lambda \sum_{j=1}^p \beta_j^2 \right\}
    \]

    其中，$\lambda \ge 0$ 为调节复杂度的超参数（Tuning Parameter）。

    * 当 $\lambda = 0$ 时，脊回归退化为传统的 OLS 估计。
    * 当 $\lambda \to \infty$ 时，惩罚项占据主导，使得所有回归系数 $\beta_j \to 0$（不含截距项）。

### 1.1 数据的中心化与矩阵表示

在实际求解前，通常先对特征进行标准化，并对响应变量与特征进行中心化（Centering），使得中心化后的数据满足 $\sum_{i=1}^n y_i = 0$ 且 $\sum_{i=1}^n x_{ij} = 0$。此时截距项的估计量 $\hat{\beta}_0 = 0$，我们可以直接写出参数估计 $\beta = (\beta_1, \dots, \beta_p)^{\top}$ 的矩阵目标函数：

\[
S_{\lambda}(\beta) = (Y - X\beta)^{\top}(Y - X\beta) + \lambda \beta^{\top}\beta
\]

!!! note "定理 1.1 (脊回归估计量的显式解)"

    对于任意给定的 $\lambda > 0$，脊回归的参数估计量 $\hat{\beta}_{\text{ridge}}$ 具有如下唯一的闭式解：

    \[
    \hat{\beta}_{\text{ridge}} = (X^{\top}X + \lambda I_p)^{-1}X^{\top}Y
    \]

    其中 $I_p$ 为 $p \times p$ 的单位矩阵。

??? proof "证明：脊回归显式解的推导"

    首先，将矩阵目标函数 $S_{\lambda}(\beta)$ 展开：

    \[
    S_{\lambda}(\beta) = Y^{\top}Y - 2\beta^{\top}X^{\top}Y + \beta^{\top}X^{\top}X\beta + \lambda \beta^{\top}\beta
    \]

    利用矩阵微积分对参数向量 $\beta$ 求偏导数：

    \[
    \frac{\partial S_{\lambda}(\beta)}{\partial \beta} = -2X^{\top}Y + 2X^{\top}X\beta + 2\lambda \beta
    \]

    令偏导数向量等于 0，得到修正后的正规方程组：

    \[
    (X^{\top}X + \lambda I_p)\beta = X^{\top}Y
    \]

    由于当 $\lambda > 0$ 时，即使 $X^{\top}X$ 基底矩阵是奇异的（例如 $p > n$），矩阵 $(X^{\top}X + \lambda I_p)$ 必然是严格正定的，因而永远可逆。在两边同左乘逆矩阵，即证：

    \[
    \hat{\beta}_{\text{ridge}} = (X^{\top}X + \lambda I_p)^{-1}X^{\top}Y
    \]

### 1.2 脊回归的期望与方差推导

现假设真实模型为 $Y = X\beta + \epsilon$，其中 $\mathbb{E}(\epsilon) = 0$，$\text{Var}(\epsilon) = \sigma^2 I_n$。

!!! note "定理 1.2 (脊回归的统计性质)"

    脊回归估计量 $\hat{\beta}_{\text{ridge}}$ 是一个有偏估计，其期望向量与协方差矩阵分别为：

    \[
    \mathbb{E}(\hat{\beta}_{\text{ridge}}) = (X^{\top}X + \lambda I_p)^{-1}X^{\top}X\beta
    \]

    \[
    \text{Var}(\hat{\beta}_{\text{ridge}}) = \sigma^2 (X^{\top}X + \lambda I_p)^{-1}X^{\top}X(X^{\top}X + \lambda I_p)^{-1}
    \]

??? proof "证明：脊回归期望与方差的数学推导"

    **1. 期望推导：**

    代入 $Y = X\beta + \epsilon$ 进入显式解中：

    \[
    \hat{\beta}_{\text{ridge}} = (X^{\top}X + \lambda I_p)^{-1}X^{\top}(X\beta + \epsilon) = (X^{\top}X + \lambda I_p)^{-1}X^{\top}X\beta + (X^{\top}X + \lambda I_p)^{-1}X^{\top}\epsilon
    \]

    两边取条件期望，由于 $\mathbb{E}(\epsilon) = 0$：

    \[
    \mathbb{E}(\hat{\beta}_{\text{ridge}}) = (X^{\top}X + \lambda I_p)^{-1}X^{\top}X\beta
    \]

    显然，除非 $\lambda = 0$，否则 $\mathbb{E}(\hat{\beta}_{\text{ridge}}) \neq \beta$，即估计量是有偏的。

    **2. 方差矩阵推导：**

    根据线性变换的方差性质 $\text{Var}(AZ) = A\text{Var}(Z)A^{\top}$，由于 $\beta$ 为非随机常数：

    \[
    \text{Var}(\hat{\beta}_{\text{ridge}}) = \text{Var}\left( (X^{\top}X + \lambda I_p)^{-1}X^{\top}\epsilon \right)
    \]

    提取常数矩阵并结合 $\text{Var}(\epsilon) = \sigma^2 I_n$ 展开，即证：

    \[
    \text{Var}(\hat{\beta}_{\text{ridge}}) = (X^{\top}X + \lambda I_p)^{-1}X^{\top} \left( \sigma^2 I_n \right) \left[ (X^{\top}X + \lambda I_p)^{-1}X^{\top} \right]^{\top}
    \]

    \[
    \text{Var}(\hat{\beta}_{\text{ridge}}) = \sigma^2 (X^{\top}X + \lambda I_p)^{-1}X^{\top}X(X^{\top}X + \lambda I_p)^{-1}
    \]

---

## 2. Lasso 回归（Least Absolute Shrinkage and Selection Operator）

!!! info "定义 2.1 (Lasso 优化问题)"

    Lasso 回归通过引入回归系数的 $L_1$ 范数惩罚项来实现系数的收缩。中心化设计矩阵下的目标函数定义为：

    \[
    \min_{\beta \in \mathbb{R}^p} \left\{ \sum_{i=1}^n \left( y_i - \sum_{j=1}^p \beta_j x_{ij} \right)^2 + \lambda \sum_{j=1}^p |\beta_j| \right\}
    \]

### 2.1 稀疏性（Sparsity）的几何解释

由于 $L_1$ 惩罚项中绝对值函数 $|\beta_j|$ 在零点处不可导，其对应的等高约束域是一个多面体（在二维空间中是菱形），而脊回归的 $L_2$ 约束域是一个超球体。当 RSS 的椭圆形等高线与约束域首次相交时，Lasso 的交点极易落在多面体的顶点或棱边上（即某些轴上），从而迫使部分回归系数**精确等于 0**。这使得 Lasso 兼具参数收缩与特征选择（Feature Selection）的功能。

---

## 3. 正交设计矩阵下的解析解与特例分析

为了深入理解不同惩罚项对系数产生的收缩效应，我们考虑一种特例：设计矩阵是严格正交的，即满足 $X^{\top}X = I_p$。此时，传统的 OLS 估计量为 $\hat{\beta}_{\text{OLS}} = X^{\top}Y$。

!!! note "定理 3.1 (正交设计下的特殊估计量解)"

    在 $X^{\top}X = I_p$ 的正交设计下，脊回归与 Lasso 的解可以表示为 $\hat{\beta}_{\text{OLS}}$ 的显式单变量函数：

    * **脊回归的解**：

    \[
    \hat{\beta}_{j,\text{ridge}} = \frac{1}{1 + \lambda} \hat{\beta}_{j,\text{OLS}}
    \]

    * **Lasso 的解（软阈值算子 Soft Thresholding）**：

    \[
    \hat{\beta}_{j,\text{Lasso}} = \text{sign}(\hat{\beta}_{j,\text{OLS}}) \max \left( |\hat{\beta}_{j,\text{OLS}}| - \frac{\lambda}{2}, 0 \right)
    \]

??? proof "证明：正交设计下显式收缩解的推导"

    **1. 脊回归正交解：**

    将 $X^{\top}X = I_p$ 代入脊回归的解析解中：

    \[
    \hat{\beta}_{\text{ridge}} = (I_p + \lambda I_p)^{-1}X^{\top}Y = \frac{1}{1 + \lambda} X^{\top}Y
    \]

    由于 $\hat{\beta}_{\text{OLS}} = X^{\top}Y$，故对于每个分量：

    \[
    \hat{\beta}_{j,\text{ridge}} = \frac{1}{1 + \lambda} \hat{\beta}_{j,\text{OLS}}
    \]

    **2. Lasso 正交解：**

    在正交设计下，Lasso 目标函数可以拆解为 $p$ 个独立的单变量优化问题。对第 $j$ 个分量有：

    \[
    Q(\beta_j) = (\hat{\beta}_{j,\text{OLS}} - \beta_j)^2 + \lambda |\beta_j|
    \]

    我们通过对 $\beta_j$ 的正负号进行分情况讨论来求解 $\min Q(\beta_j)$：

    * 当 $\beta_j > 0$ 时，$|\beta_j| = \beta_j$，求导：

    \[
    \frac{dQ}{d\beta_j} = -2(\hat{\beta}_{j,\text{OLS}} - \beta_j) + \lambda = 0 \implies \hat{\beta}_j = \hat{\beta}_{j,\text{OLS}} - \frac{\lambda}{2}
    \]

    此解成立的必要条件是 $\hat{\beta}_{j,\text{OLS}} > \frac{\lambda}{2}$。

    * 当 $\beta_j < 0$ 时，$|\beta_j| = -\beta_j$，求导：

    \[
    \frac{dQ}{d\beta_j} = -2(\hat{\beta}_{j,\text{OLS}} - \beta_j) - \lambda = 0 \implies \hat{\beta}_j = \hat{\beta}_{j,\text{OLS}} + \frac{\lambda}{2}
    \]

    此解成立的必要条件是 $\hat{\beta}_{j,\text{OLS}} < -\frac{\lambda}{2}$。

    * 当 $-\frac{\lambda}{2} \le \hat{\beta}_{j,\text{OLS}} \le \frac{\lambda}{2}$ 时，最小值在边界点 $\beta_j = 0$ 处取得。

    综合以上三种情况，利用符号函数 $\text{sign}(\cdot)$ 和截断函数 $\max(\cdot, 0)$，即证：

    \[
    \hat{\beta}_{j,\text{Lasso}} = \text{sign}(\hat{\beta}_{j,\text{OLS}}) \max \left( |\hat{\beta}_{j,\text{OLS}}| - \frac{\lambda}{2}, 0 \right)
    \]

---

## 4. 惩罚项的扩展：$L_q$ 范数与弹性网

为了结合脊回归和 Lasso 的优势，可以对惩罚项的形式进行一般化推广。

### 4.1 广义 $L_q$ 惩罚（Bridge Penalty）

我们考虑通用形式的回归系数惩罚项：

\[
\lambda \sum_{j=1}^p |\beta_j|^q, \quad q \ge 0
\]

* 当 $q = 2$ 时，对应脊回归（$L_2$ 范数）。
* 当 $q = 1$ 时，对应 Lasso 回归（$L_1$ 范数）。
* 当 $q = 0$ 时，惩罚项变为非零系数的个数 $\sum \mathbb{I}(\beta_j \neq 0)$，对应最优子集选择。
* *几何特征*：当且仅当 $q \le 1$ 时，约束域在零点处存在尖角，模型才具备产生稀疏解（特征选择）的能力；当 $q < 1$ 时，约束域是非凸的，优化问题将变得非常难以求解。

### 4.2 弹性网（Elastic Net）

弹性网通过将 $L_1$ 与 $L_2$ 惩罚项进行线性组合，克服了 Lasso 在处理高相关特征组时的局限性（Lasso 在高度相关的特征中往往随机选择一个，而弹性网倾向于将其成组保留）。其目标函数定义为：

\[
\min_{\beta} \left\{ \sum_{i=1}^n \left( y_i - \sum_{j=1}^p \beta_j x_{ij} \right)^2 + \lambda_1 \sum_{j=1}^p |\beta_j| + \lambda_2 \sum_{j=1}^p \beta_j ^2 \right\}
\]
