# 第四章：高级惩罚线性回归与选择一致性

在上一章中，我们讨论了脊回归与 Lasso 的基本设定。虽然 Lasso 通过 $L_1$ 惩罚实现了变量选择，但从统计学理论来看，它在收缩大系数时会产生显著的系统性偏置。本章将从渐近理论出发，深入探讨惩罚估计量的“神谕性质”（Oracle Properties）、偏置修正方法（Adaptive Lasso, SCAD, MCP）以及针对特殊数据结构的惩罚项扩展。

---

## 1. 理论动机：高维稀疏设定与神谕性质

!!! info "定义 1.1 (高维稀疏模型设定)"

    假设观测数据满足多元线性模型 $Y = X\beta + \epsilon$。在特征维度 $p$ 很大甚至满足 $p > n$ 的高维设定下，我们假定真实的回归系数向量 $\beta^* = (\beta_1^*, \dots, \beta_p^*)^{\top}$ 满足**稀疏性（Sparsity）**假设。

    不失一般性，设前 $p_0$ 个变量为真正相关的特征，其余 $p - p_0$ 个变量对应的系数为 0。定义真实非零系数的指标集合为：

    \[
    \mathcal{S} = \{j : \beta_j^* \neq 0\}, \quad |\mathcal{S}| = p_0
    \]

    相应地，真实的系数向量可写为分块矩阵形式：$\beta^* = (\beta_{\mathcal{S}}^{*\top}, \beta_{\mathcal{S}^c}^{*\top})^{\top}$，其中 $\beta_{\mathcal{S}^c}^* = 0$。

### 1.1 神谕性质（Oracle Properties）

一个理想的参数化估计量 $\hat{\beta}$ 在处理高维稀疏模型时，应当如同拥有“神谕”（预知哪个变量真实有效）一般完美。统计学中将满足以下两个条件的估计量称为具有**神谕性质**：

1. **选择一致性（Sparsity/Selection Consistency）**：能够以趋于 1 的概率正确识别出真实非零变量的集合，即：

\[
\lim_{n \to \infty} P(\hat{\mathcal{S}} = \mathcal{S}) = 1
\]

2. **渐近正态性（Asymptotic Normality）**：对于非零系数部分，其估计量的收敛速度与真正了解稀疏结构时的 OLS 估计量完全相同，即：

\[
\sqrt{n}(\hat{\beta}_{\mathcal{S}} - \beta_{\mathcal{S}}^*) \xrightarrow{d} N(0, \Sigma^*)
\]

### 1.2 Lasso 的局限性与不可调和的矛盾

Lasso 估计量采用统一的惩罚力度 $\lambda |\beta_j|$。为了消除不相关变量（使得 $\hat{\beta}_{\mathcal{S}^c} = 0$），必须令 $\lambda$ 足够大；然而，大惩罚力度会对真正的大系数 $\beta_{\mathcal{S}}^*$ 产生过度的收缩，从而导致显著的系统性偏置，使其在数学上难以同时完美满足选择一致性与渐近正态性。

---

## 2. 偏置修正方法（一）：自适应 Lasso 与非负加洛特

为了克服 Lasso 的系统性偏置，首要思想是对不同的系数施加动态调整的、不平等的惩罚权重。

### 2.1 非负加洛特估计（Non-negative Garrote）

非负加洛特是最早提出的收缩与选择方法之一。它首先利用 OLS（或脊回归）得到初始无偏估计 $\hat{\beta}^{\text{init}}$，然后通过解以下约束优化问题对初始估计进行非负缩放：

\[
\min_{d_1, \dots, d_p} \left\{ \frac{1}{2n} \left\| Y - \sum_{j=1}^p d_j X_j \hat{\beta}_j^{\text{init}} \right\|^2 + \lambda \sum_{j=1}^p d_j \right\}
\]

\[
\text{subject to } d_j \ge 0, \quad \forall j = 1, \dots, p
\]

最终的系数估计值为 $\hat{\beta}_{j}^{\text{garrote}} = d_j \hat{\beta}_j^{\text{init}}$。

### 2.2 自适应 Lasso（Adaptive Lasso）

自适应 Lasso 进一步推广了这一思想，直接在 $L_1$ 惩罚项内为每个特征引入特定的自适应权重。

!!! info "定义 2.1 (自适应 Lasso 优化问题)"

    自适应 Lasso 的目标函数定义为：

    \[
    \min_{\beta} \left\{ \frac{1}{2n} \|Y - X\beta\|^2 + \lambda \sum_{j=1}^p w_j |\beta_j| \right\}
    \]

    其中的自适应权重向量 $w = (w_1, \dots, w_p)^{\top}$ 通常利用一个已有的相合估计量 $\hat{\beta}^{\text{init}}$（如 OLS 估计量或脊回归估计量）构造：

    \[
    w_j = \frac{1}{|\hat{\beta}_j^{\text{init}}|^\gamma}
    \]

    其中 $\gamma > 0$ 是一个常数（通常设为 $\gamma = 1$ 或 2）。

* *统计直觉*：如果某个特征的真实系数 $\beta_j^* \neq 0$，随着样本量增加，其初始估计 $|\hat{\beta}_j^{\text{init}}|$ 会远离 0，导致权重 $w_j$ 变小，因而其受到的惩罚非常微弱，**消除了偏置**；若真实系数 $\beta_j^* = 0$，其初始估计 $|\hat{\beta}_j^{\text{init}}| \to 0$，导致权重 $w_j \to \infty$，使其受到极强的惩罚，从而**实现精准的变量剔除**。理论证明，自适应 Lasso 完美具备神谕性质。

---

## 3. 偏置修正方法（二）：非凸惩罚算子（SCAD 与 MCP）

除了引入自适应权重外，另一种更为根本的解决方案是直接修改惩罚函数 $P_{\lambda}(\cdot)$ 的几何形态，令其对大系数的导数衰减为 0。

!!! note "定理 3.1 (神谕估计量的三项基本特性准则)"

    一个优秀的惩罚函数 $P_{\lambda}(|\beta|)$ 若想产生具备神谕性质的估计，其导数 $P_{\lambda}'(|\beta|)$ 应当同时满足以下三个数学微积分准则：

    1. **稀疏性（Sparsity）**：$\min_{t} \{t + P_{\lambda}'(t)\} > 0$，使得小系数能被直接归零。
    2. **无偏性（Unbiasedness）**：当 $|t|$ 很大时，$P_{\lambda}'(t) = 0$，保证大系数不承受额外偏置。
    3. **连续性（Continuity）**：$\arg\min_{t} \{ \frac{1}{2}(z - t)^2 + P_{\lambda}(t) \}$ 作为 $z$ 的函数必须是连续的，防止预测值产生剧烈跳跃。

### 3.1 SCAD 惩罚算子（Smoothly Clipped Absolute Deviation）

SCAD 惩罚函数的导数形式由如下分段函数定义（其中超参数 $a > 2$，通常推荐 $a = 3.7$）：

\[
P_{\lambda}'(\theta) = \lambda \mathbb{I}(\theta \le \lambda) + \frac{(a\lambda - \theta)_+}{a - 1} \mathbb{I}(\theta > \lambda)
\]

??? proof "证明：正交设计下 SCAD 显式解（阈值算子）的推导"

    在正交设计矩阵 $X^{\top}X = I$ 设定下，SCAD 估计量对初始 OLS 估计量 $z_j = \hat{\beta}_{j}^{\text{OLS}}$ 的单变量解析解可以通过对导数进行积分后分类求极值导出。

    对其分段求解可得如下硬平滑阈值算子式：

    \[
    \hat{\beta}_{j}^{\text{SCAD}} = \begin{cases} \text{sign}(z_j)(|z_j| - \lambda)_+, & \text{当 } |z_j| \le 2\lambda \\ \frac{(a-1)z_j - \text{sign}(z_j)a\lambda}{a-2}, & \text{当 } 2\lambda < |z_j| \le a\lambda \\ z_j, & \text{当 } |z_j| > a\lambda \end{cases}
    \]

    由该解析解可见，当 OLS 估计量 $|z_j| > a\lambda$ 时，估计值直接等于 $z_j$（即回归本身），**大系数偏置完全转为 0**，满足无偏性准则。

### 3.2 MCP 惩罚算子（Minimax Concave Penalty）

MCP 相比 SCAD 其凹度变化更加迅捷，其导数定义为（其中 $a > 1$）：

\[
P_{\lambda}'(\theta) = \left( \lambda - \frac{\theta}{a} \right)_+
\]

当 $\theta > a\lambda$ 时，导数直接变为 0。在正交设计下，其单变量解析解为：

\[
\hat{\beta}_{j}^{\text{MCP}} = \begin{cases} \text{sign}(z_j) \frac{(|z_j| - \lambda)_+}{1 - 1/a}, & \text{当 } |z_j| \le a\lambda \\ z_j, & \text{当 } |z_j| > a\lambda \end{cases}
\]

---

## 4. 针对特殊数据结构的惩罚项扩展

在众多实际场景中，变量之间天生存在着某种特殊的结构（如分组结构或序列结构），普通的 Lasso 会割裂这些联系。

### 4.1 分组 Lasso（Group Lasso）

!!! info "定义 4.1 (分组 Lasso 优化问题)"

    假设 $p$ 个预测变量被划分为 $K$ 个互不相交的组。令 $X^{(k)}$ 表示属于第 $k$ 组的特征子矩阵，$\beta^{(k)}$ 表示对应的系数子向量。分组 Lasso 目标函数定义为：

    \[
    \min_{\beta} \left\{ \frac{1}{2} \left\| Y - \sum_{k=1}^K X^{(k)}\beta^{(k)} \right\|^2 + \lambda \sum_{k=1}^K \sqrt{p_k} \|\beta^{(k)}\|_2 \right\}
    \]

    其中 $p_k$ 为第 $k$ 组中包含的变量个数，$\|\beta^{(k)}\|_2$ 是回归系数向量的欧氏 $L_2$ 范数。

* *统计学含义*：该惩罚项在组与组之间构成了 $L_1$ 范数（产生组间稀疏性），而在组内部构成了 $L_2$ 范数。因此，它能够迫使整个组的系数**要么全为 0（被整体剔除），要么全部不为 0（被整体保留）**。

### 4.2 融合 Lasso（Fused Lasso）

当特征存在某种天然的先后顺序（例如时序信号、基因组沿染色体排列的位点）时，我们不仅希望系数本身稀疏，还希望相邻系数之间具有连续性。其目标函数定义为：

\[
\min_{\beta} \left\{ \frac{1}{2} \sum_{i=1}^n \left( y_i - \sum_{j=1}^p \beta_j x_{ij} \right)^2 + \lambda_1 \sum_{j=1}^p |\beta_j| + \lambda_2 \sum_{j=2}^p |\beta_j - \beta_{j-1}| \right\}
\]

通过对相邻系数差值 $|\beta_j - \beta_{j-1}|$ 施加 $L_1$ 惩罚，融合 Lasso 可以迫使相邻系数产生**精确相等的阶梯状分块常数效应**。

---

## 5. 协同回归（Collaborative Regression）

在生物医学等多源数据集成分析中，我们常拥有来自不同渠道提取的特征（例如 $X$ 代表 DNA 甲基化矩阵，$Z$ 代表 RNA 表达矩阵）。协同回归通过在目标函数中加入鼓励两组预测结果对齐的惩罚项来联合建模。

其基本的三项平方损失优化目标为：

\[
\min_{\beta, \theta} \left\{ \|Y - X\beta\|_2^2 + \|Y - Z\theta\|_2^2 + \gamma \|X\beta - Z\theta\|_2^2 + \lambda_1 P(\beta) + \lambda_2 P(\theta) \right\}
\]

!!! note "定理 5.1 (协同回归的设计矩阵增强求解机制)"

    上述包含两特征源的协同回归优化问题，在忽略惩罚项时，可以通过构造如下的高维“巨型”增强设计矩阵（Design Matrix Augmentation）与增强观测向量，完美等价转化为一个标准的最小二乘问题：

    \[
    \tilde{X} = \begin{pmatrix} X & 0 \\ 0 & Z \\ \sqrt{\gamma}X & -\sqrt{\gamma}Z \end{pmatrix}, \quad \tilde{Y} = \begin{pmatrix} Y \\ Y \\ 0 \end{pmatrix}
    \]

??? proof "证明：设计矩阵增强等价性的推导"

    写出增强系统下标准残差平方和的矩阵表达式： $\|\tilde{Y} - \tilde{X} \cdot (\beta^{\top}, \theta^{\top})^{\top} \|_2^2$。

    根据矩阵的分块乘法，我们将其按行块展开：

    \[
    \tilde{Y} - \tilde{X}\begin{pmatrix} \beta \\ \theta \end{pmatrix} = \begin{pmatrix} Y \\ Y \\ 0 \end{pmatrix} - \begin{pmatrix} X\beta + 0 \\ 0 + Z\theta \\ \sqrt{\gamma}X\beta - \sqrt{\gamma}Z\theta \end{pmatrix} = \begin{pmatrix} Y - X\beta \\ Y - Z\theta \\ -\sqrt{\gamma}(X\beta - Z\theta) \end{pmatrix}
    \]

    对其求欧氏 $L_2$ 范数的平方（即各分块欧氏平方和相加）：

    \[
    \left\| \tilde{Y} - \tilde{X}\begin{pmatrix} \beta \\ \theta \end{pmatrix} \right\|_2^2 = \|Y - X\beta\|_2^2 + \|Y - Z\theta\|_2^2 + \left\| -\sqrt{\gamma}(X\beta - Z\theta) \right\|_2^2
    \]

    将常数项开方后移出范数：

    \[
    \left\| \tilde{Y} - \tilde{X}\begin{pmatrix} \beta \\ \theta \end{pmatrix} \right\|_2^2 = \|Y - X\beta\|_2^2 + \|Y - Z\theta\|_2^2 + \gamma \|X\beta - Z\theta\|_2^2
    \]

    上式与协同回归的主体损失函数完全一致。因此，我们可以直接利用这个变换后的巨型设计矩阵 $\tilde{X}$ 和响应向量 $\tilde{Y}$，套用通用的 Lasso 或脊回归算法完成求解。