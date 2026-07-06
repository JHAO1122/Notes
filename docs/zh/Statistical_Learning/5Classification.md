# 第五章：分类问题与判别分析

在前面的章节中，我们主要讨论了响应变量 $Y$ 为连续型变量的回归问题。本章将聚焦于分类问题（Classification），即响应变量 $Y$ 取值为有限个离散标签（如 $\{0, 1\}$ 或 $\{1, 2, \dots, K\}$）的情形。我们将系统地推导软分类器（Soft Classifier）的两大基石：逻辑回归（Logistic Regression）与线性/二次判别分析（LDA/QDA），并严格给出相关的数学推导与性质证明。

---

## 1. 分类问题的基本设定与决策理论

!!! info "定义 1.1 (分类器与 0-1 损失)"

    设训练数据集为 $\mathcal{D}_{n}=\{(x_{i}, y_{i})\}_{i=1}^{n}$，其中特征向量 $x_{i} \in \mathbb{R}^p$，响应变量 $y_{i} \in \{0, 1\}$（或二分类中的 $\{-1, +1\}$）。分类任务的目标是学习一个映射函数（即分类器）：

    \[
    f: \mathbb{R}^p \longrightarrow \{0, 1\}
    \]

    其表现通常由 0-1 损失函数（0-1 Loss Function）来衡量：

    \[
    L(f(X), Y) = \mathbb{I}(f(X) \neq Y) = \begin{cases} 0 & \text{若 } f(X) = Y \\ 1 & \text{若 } f(X) \neq Y \end{cases}
    \]

### 1.1 贝叶斯最优分类器（Bayes Optimal Classifier）

为了使整体的期望风险（Expected Risk）达到最小，我们需要寻找最优的分类规则。

!!! note "定理 1.1 (贝叶斯决策准则)"

    在 0-1 损失下，使期望风险 $R(f) = \mathbb{E}_{X,Y}[L(f(X), Y)]$ 最小化的决策函数为贝叶斯最优分类器 $f^*(x)$，其表达式为：

    \[
    f^*(x) = \arg\max_{c \in \{0, 1\}} P(Y = c \mid X = x)
    \]

??? proof "证明：贝叶斯决策最优性的推导"

    期望风险可以利用全期望公式，通过对特征 $X$ 施加条件期望来展开：

    \[
    R(f) = \mathbb{E}_X \left[ \mathbb{E}_{Y \mid X} [ \mathbb{I}(f(X) \neq Y) \mid X ] \right]
    \]

    为了使整体期望最小，只需对于每一个固定的 $X = x$，都使其条件期望（条件风险）达到最小：

    \[
    \mathbb{E}_{Y \mid X} [ \mathbb{I}(f(x) \neq Y) \mid X = x ] = P(f(x) \neq Y \mid X = x)
    \]

    对于二分类问题，上式可以写为：

    \[
    P(f(x) \neq Y \mid X = x) = 1 - P(Y = f(x) \mid X = x)
    \]

    显然，要使 $1 - P(Y = f(x) \mid X = x)$ 最小，等价于使 $P(Y = f(x) \mid X = x)$ 最大。因此，分类器在点 $x$ 处的预测值 $f(x)$ 应当选择使得条件概率最大化的类别：

    \[
    f^*(x) = \arg\max_{c \in \{0, 1\}} P(Y = c \mid X = x)
    \]

    即证，后验概率最大化准则（贝叶斯准则）在 0-1 损失下是最优的。

---

## 2. 逻辑回归（Logistic Regression）

逻辑回归是一种直接对后验条件概率 $P(Y=1 \mid X=x)$ 进行建模的参数化方法。

!!! info "定义 2.1 (Logistic 模型与 Logit 变换)"

    令 $\eta(x) = P(Y = 1 \mid X = x)$。为了将定义域为 $\mathbb{R}^p$ 的线性预测器映射到 $[0, 1]$ 概率区间，模型引入了 Logit 变换（对数几率）：

    \[
    \ln \left( \frac{\eta(x)}{1 - \eta(x)} \right) = \beta_0 + \beta^{\top}x = \beta_0 + \sum_{j=1}^p \beta_j X_j
    \]

    通过代数反解，可直接得到条件概率的 Sigmoid 函数形式：

    \[
    \eta(x) = \frac{\exp(\beta_0 + \beta^{\top}x)}{1 + \exp(\beta_0 + \beta^{\top}x)}
    \]

### 2.1 最大似然估计（Maximum Likelihood Estimation, MLE）

由于响应变量 $Y_i \mid x_i \sim \text{Bernoulli}(\eta(x_i))$，我们可以写出参数 $\theta = (\beta_0, \beta^{\top})^{\top}$ 的似然函数。为了表达紧凑，令样本扩展特征 $x_i = (1, x_{i1}, \dots, x_{ip})^{\top}$。

!!! note "定理 2.2 (逻辑回归的对数似然函数与梯度)"

    逻辑回归参数的最大似然估计通过最小化负对数似然函数（即交叉熵损失）实现：

    \[
    \ell(\theta) = -\sum_{i=1}^n \left[ y_i \ln \eta(x_i) + (1 - y_i) \ln (1 - \eta(x_i)) \right]
    \]

    其关于参数 $\theta$ 的梯度向量（一阶导数）为：

    \[
    \nabla \ell(\theta) = \sum_{i=1}^n (\eta(x_i) - y_i) x_i
    \]

??? proof "证明：梯度公式的详细微积分推导"

    首先，计算 Sigmoid 函数 $\eta_i = \eta(x_i) = \frac{1}{1 + e^{-\theta^{\top}x_i}}$ 对 $\theta$ 的导数。由链式法则有：

    \[
    \frac{\partial \eta_i}{\partial \theta} = \frac{e^{-\theta^{\top}x_i}}{(1 + e^{-\theta^{\top}x_i})^2} \cdot x_i = \eta_i (1 - \eta_i) x_i
    \]

    接下来，对对数似然项 $L_i = y_i \ln \eta_i + (1 - y_i) \ln (1 - \eta_i)$ 求关于 $\theta$ 的偏导数：

    \[
    \frac{\partial L_i}{\partial \theta} = \frac{y_i}{\eta_i} \frac{\partial \eta_i}{\partial \theta} - \frac{1 - y_i}{1 - \eta_i} \frac{\partial \eta_i}{\partial \theta}
    \]

    将 $\frac{\partial \eta_i}{\partial \theta} = \eta_i (1 - \eta_i) x_i$ 代入上式中：

    \[
    \frac{\partial L_i}{\partial \theta} = \left[ \frac{y_i}{\eta_i} \eta_i (1 - \eta_i) - \frac{1 - y_i}{1 - \eta_i} \eta_i (1 - \eta_i) \right] x_i
    \]

    \[
    \frac{\partial L_i}{\partial \theta} = \left[ y_i (1 - \eta_i) - (1 - y_i) \eta_i \right] x_i = (y_i - \eta_i) x_i
    \]

    因此，负对数似然函数 $\ell(\theta) = -\sum_{i=1}^n L_i$ 的梯度向量为：

    \[
    \nabla \ell(\theta) = -\sum_{i=1}^n (y_i - \eta_i) x_i = \sum_{i=1}^n (\eta(x_i) - y_i) x_i
    \]

    即证完毕。由于该目标函数关于 $\theta$ 是严格对数凹的（Log-concave），实际计算中通常采用牛顿-拉夫森法（Newton-Raphson）或梯度下降法迭代求解。

---

## 3. 生成式模型：判别分析（Discriminant Analysis）

与逻辑回归不同，判别分析属于**生成式模型（Generative Model）**。它先对每一个类别 $k$ 内特征的条件概率密度 $f_k(x) = f(X=x \mid Y=k)$ 以及类先验概率 $\pi_k = P(Y=k)$ 进行建模，随后利用贝叶斯公式反推后验概率。

根据贝叶斯公式，已知 $X=x$ 时属于第 $k$ 类的后验概率为：

\[
P(Y = k \mid X = x) = \frac{\pi_k f_k(x)}{\sum_{l=1}^K \pi_l f_l(x)}
\]

我们假设每一类中的数据均服从多元正态分布（Multivariate Normal Distribution）：

\[
f_k(x) = \frac{1}{(2\pi)^{p/2} |\Sigma_k|^{1/2}} \exp \left( -\frac{1}{2} (x - \mu_k)^{\top} \Sigma_k^{-1} (x - \mu_k) \right)
\]

### 3.1 线性判别分析（Linear Discriminant Analysis, LDA）

!!! info "定义 3.1 (LDA 基本假设)"

    线性判别分析（LDA）假设所有类别的特征协方差矩阵均相等，即：

    \[
    \Sigma_k = \Sigma, \quad \forall k = 1, 2, \dots, K
    \]

!!! note "定理 3.1 (LDA 线性决策函数)"

    在各类别协方差矩阵相等的假设下，最大化贝叶斯后验概率等价于最大化如下关于 $x$ 的线性判别函数（Linear Discriminant Function）$\delta_k(x)$：

    \[
    \delta_k(x) = x^{\top} \Sigma^{-1} \mu_k - \frac{1}{2} \mu_k^{\top} \Sigma^{-1} \mu_k + \ln \pi_k
    \]

??? proof "证明：LDA 线性判别函数的数学推导"

    要最大化 $P(Y=k \mid X=x)$，由于贝叶斯公式的分母 $\sum_{l=1}^K \pi_l f_l(x)$ 对所有类别 $k$ 都是完全相同的，因此只需要最大化分子项：

    \[
    \arg\max_k P(Y=k \mid X=x) = \arg\max_k \left[ \pi_k f_k(x) \right]
    \]

    由于对数函数是单调递增的，我们在分子项上取自然对数：

    \[
    \ln \left( \pi_k f_k(x) \right) = \ln \pi_k + \ln f_k(x)
    \]

    将多元正态分布密度函数 $f_k(x)$（其中 $\Sigma_k = \Sigma$）代入展开：

    \[
    \ln \left( \pi_k f_k(x) \right) = \ln \pi_k - \frac{p}{2}\ln(2\pi) - \frac{1}{2}\ln|\Sigma| - \frac{1}{2}(x - \mu_k)^{\top} \Sigma^{-1} (x - \mu_k)
    \]

    注意到项 $-\frac{p}{2}\ln(2\pi) - \frac{1}{2}\ln|\Sigma|$ 是与类别 $k$ 无关的常数，在进行关于 $k$ 的最大化决策时可以将其剔除。我们将剩余部分展开：

    \[
    - \frac{1}{2}(x - \mu_k)^{\top} \Sigma^{-1} (x - \mu_k) + \ln \pi_k = -\frac{1}{2} \left[ x^{\top}\Sigma^{-1}x - 2x^{\top}\Sigma^{-1}\mu_k + \mu_k^{\top}\Sigma^{-1}\mu_k \right] + \ln \pi_k
    \]

    再次注意到 $x^{\top}\Sigma^{-1}x$ 这一项也与类别 $k$ 无关，因此可以再次忽略。化简合并后，最终留下与 $k$ 相关的核心部分：

    \[
    \delta_k(x) = x^{\top} \Sigma^{-1} \mu_k - \frac{1}{2} \mu_k^{\top} \Sigma^{-1} \mu_k + \ln \pi_k
    \]

    由于上式中关于未知变量 $x$ 的最高次项为一次项，因此不同类别间的决策边界（即 $\delta_k(x) = \delta_l(x)$ 的超平面）在几何上呈现为**线性超平面**。证毕。

### 3.2 二次判别分析（Quadratic Discriminant Analysis, QDA）

!!! info "定义 3.2 (QDA 基本假设)"

    二次判别分析（QDA）放宽了协方差矩阵相等的限制，允许不同类别的群体拥有各自独立的特征协方差矩阵 $\Sigma_k$。

此时，无法在对数似然表达式中消去二次项与行列式项。保留完整的特征参数后，QDA 的二次判别函数（Quadratic Discriminant Function）表述为：

\[
\delta_k(x) = -\frac{1}{2} \ln |\Sigma_k| - \frac{1}{2} (x - \mu_k)^{\top} \Sigma_k^{-1} (x - \mu_k) + \ln \pi_k
\]

由于项 $x^{\top}\Sigma_k^{-1}x$ 无法被消去，不同类别之间的决策边界将呈现为**二次曲面**（如双曲线、椭圆或抛物面）。

---

## 4. 逻辑回归与 LDA 的深层统计学联系

如果我们在二分类问题中仔细审视 LDA 下两类的对数后验几率（Log-odds ratio），会发现一个惊人的结论。

!!! note "定理 4.1 (LDA 的隐含 Logit 形式)"

    对于 LDA 模型，类别 1 与类别 0 之间的对数后验几率关于特征 $x$ 具有精确的线性形式：

    \[
    \ln \left( \frac{P(Y=1 \mid X=x)}{P(Y=0 \mid X=x)} \right) = \alpha_0 + \alpha^{\top}x
    \]

    其中参数对应关系为：

    \[
    \alpha_0 = \ln \frac{\pi_1}{\pi_0} - \frac{1}{2}\mu_1^{\top}\Sigma^{-1}\mu_1 + \frac{1}{2}\mu_0^{\top}\Sigma^{-1}\mu_0
    \]

    \[
    \alpha = \Sigma^{-1}(\mu_1 - \mu_0)
    \]

### 4.1 核心区别：条件似然 vs. 全似然

虽然由上可见，LDA 与逻辑回归在形式上都导出了线性的后验几率边界，但两者的根本区别在于其**参数估计的基准和假设强度**：

* **逻辑回归（条件似然）**：仅仅对条件概率 $P(Y \mid X)$ 进行判别式建模。它不需要对特征 $X$ 的边缘分布做出任何正态性或独立性假设。因此，当特征数据不满足正态分布（例如包含离散型哑变量）时，逻辑回归表现得更加稳健。

* **LDA（全似然）**：基于特征服从联合正态分布这一极强的生成式假设。当该假设真实成立时，LDA 能够充分榨取特征空间的信息，其在极限情况下的参数估计效率（方差）要显著优于逻辑回归。

---

## 5. 朴素贝叶斯分类器（Naive Bayes）

当特征维度 $p$ 极大时，在有限样本下精确估计 LDA 中的全协方差矩阵 $\Sigma$ 或是 QDA 中的多个 $\Sigma_k$ 会遇到严重的“维数灾难”（参数过多导致估计方差爆炸）。朴素贝叶斯算法提出了一个极端的条件独立性假设以化简计算。

!!! info "定义 5.1 (朴素贝叶斯条件独立假设)"

    朴素贝叶斯分类器假定：在给定类别标签 $Y=k$ 的条件下，各个特征变量 $X_1, X_2, \dots, X_p$ 之间是彼此条件独立的。

基于该假设，联合条件密度函数可以级联拆解为各个单变量边际密度的乘积：

\[
f_k(x) = \prod_{j=1}^p f_{kj}(x_j)
\]

代入贝叶斯决策决策准则，朴素贝叶斯的最终判别准变为：

\[
f_{\text{NB}}(x) = \arg\max_{k} \left\{ \pi_k \prod_{j=1}^p f_{kj}(x_j) \right\}
\]

* *算法优势*：通过将多维密度估计化解为 $p$ 个独立的单维密度估计 $f_{kj}(x_j)$，朴素贝叶斯将参数估计的数量级由 $\mathcal{O}(p^2)$ 直接削减到 $\mathcal{O}(p)$。其中，对于连续特征，一维密度 $f_{kj}$ 可由一维正态分布或一维核密度估计（KDE）拟合；对于离散特征，则直接使用频率直方图进行计数估计。