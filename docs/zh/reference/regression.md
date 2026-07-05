# 📈 Regression Analysis (回归分析)

本模块涵盖了线性回归模型的核心理论，包括普通最小二乘法 (OLS)、高斯-马尔可夫定理以及模型违背经典假设时的修正。

---

## 一、 多元线性回归模型与 OLS 估计

!!! info "模型矩阵形式"
    多元线性回归模型可以表示为：

    \[
    Y = X\beta + \varepsilon
    \]

    其中 $Y$ 为 $n \times 1$ 观测向量，$X$ 为 $n \times k$ 自变量设计矩阵（包含常数项），$\beta$ 为 $k \times 1$ 未知参数向量，$\varepsilon$ 为 $n \times 1$ 随机误差项。

!!! success "普通最小二乘法 (OLS)"
    通过最大程度减少残差平方和（RSS）来估计参数：

    \[
    \min_{\beta} \text{RSS} = \hat{\varepsilon}^T \hat{\varepsilon} = (Y - X\beta)^T (Y - X\beta)
    \]

    对其求导并令导数为 0 得到**正规方程组 (Normal Equations)**：$X^TX\beta = X^TY$。若 $X$ 列满秩，则最优估计量为：

    \[
    \hat{\beta} = (X^T X)^{-1} X^T Y
    \]

---

## 二、 经典高斯-马尔可夫假设与性质

!!! info "经典线性回归模型 (CLRM) 假设"
    * **1. 线性假定**：总体模型在参数上是线性的。
    
    * **2. 严格外生性**：误差项的条件期望为 0，即 $E[\varepsilon \mid X] = 0$。
    
    * **3. 无多重共线性**：秩 $\text{rank}(X) = k < n$，即自变量之间不存在完全线性相关。
    
    * **4. 同方差与无自相关**：条件方差矩阵满足：

    \[
    \text{Var}(\varepsilon \mid X) = \sigma^2 I_n
    \]

!!! success "高斯-马尔可夫定理 (Gauss-Markov Theorem)"
    在经典线性回归模型的假设下，OLS 估计量 $\hat{\beta}$ 是所有线性无偏估计量中方差最小的，即 $\hat{\beta}$ 是 **BLUE (Best Linear Unbiased Estimator，最佳线性无偏估计量)**。

??? note "估计量的统计性质"
    * **无偏性**：$E[\hat{\beta}] = \beta$。
    
    * **协方差矩阵**：$\text{Var}(\hat{\beta}) = \sigma^2 (X^T X)^{-1}$。
    
    * **随机误差项方差 $\sigma^2$ 的无偏估计**：

    \[
    s^2 = \hat{\sigma}^2 = \frac{\text{RSS}}{n - k}
    \]

---

## 三、 拟合优度与假设检验

!!! info "平方和分解 (TSS = ESS + RSS)"
    * **总平方和 (TSS)**：$\sum (y_i - \bar{y})^2$。
    
    * **回归平方和 (ESS)**：$\sum (\hat{y}_i - \bar{y})^2$。
    
    * **残差平方和 (RSS)**：$\sum (y_i - \hat{y}_i)^2$。
    
    * **判定系数 ($R^2$)**：衡量模型解释力的指标：

    \[
    R^2 = \frac{\text{ESS}}{\text{TSS}} = 1 - \frac{\text{RSS}}{\text{TSS}}
    \]

!!! success "一般线性约束与 F 检验"
    考虑对参数 $\beta$ 施加 $m$ 个线性约束，即原假设 $H_0: R\beta = r$，其中 $R$ 为 $m \times k$ 的满行秩矩阵。
    
    * **$\text{RSS}_H$**：满足原假设 $H_0$ 的约束模型下的残差平方和（施加约束后模型的拟合能力变差，故 $\text{RSS}_H \ge \text{RSS}$）。
    
    * **$\text{RSS}$**：无约束全模型下的残差平方和。
    
    在 $H_0$ 成立且误差项满足正态分布的假定下，由于 $(\text{RSS}_H - \text{RSS})/\sigma^2 \sim \chi^2(m)$ 且 $\text{RSS}/\sigma^2 \sim \chi^2(n - \text{rank}(X))$，两者相互独立，可构造如下 **F 统计量**：

    \[
    F = \frac{(\text{RSS}_H - \text{RSS}) / m}{\text{RSS} / (n - \text{rank}(X))} \sim F(m, n - \text{rank}(X))
    \]

    *(直觉：若原假设 $H_0$ 成立，施加约束不应导致残差平方和显著增大，此时 $F$ 值较小；反之若 $F$ 超过临界值，则拒绝 $H_0$)*

??? note "常用于单系数显著性的 t 检验"
    * **t 检验**：检验单个系数 $H_0: \beta_j = 0$ 的显著性：

    \[
    t = \frac{\hat{\beta}_j}{\text{se}(\hat{\beta}_j)} \sim t(n - k)
    \]

---

## 四、 违背经典假设的情形及修正

| 违背情形 | 定义/影响 | 检验方法 | 修正方法 |
| :--- | :--- | :--- | :--- |
| **异方差 (Heteroskedasticity)** | $\text{Var}(\varepsilon_i) = \sigma_i^2 \neq \sigma^2$。导致 OLS 估计量非有效，标准误偏误。 | White 检验、BP 检验 | 加权最小二乘法 (WLS)、稳健标准误 (White SE) |
| **自相关 (Autocorrelation)** | $\text{Cov}(\varepsilon_i, \varepsilon_j) \neq 0 \ (i \neq j)$。常见于时间序列，导致 OLS 估计量非有效。 | Durbin-Watson (DW) 检验、BG 检验 | 广义最小二乘法 (GLS)、Newey-West 稳健标准误 |
| **多重共线性 (Multicollinearity)** | 自变量之间存在高度线性相关。导致变动增大，单个系数不显著。 | 方差膨胀因子 (VIF > 10) | 剔除变量、岭回归 (Ridge)、主成分分析 (PCA) |