# 📖 数理统计 (Mathematical Statistics)

本模块持续收录数理统计相关的核心基础概念。

---

## 1. 假设检验与混淆矩阵 (Hypothesis Testing)

在频率学派的 Neyman-Pearson 框架下，假设检验本质上是一个基于样本对真实参数空间进行二元决策的过程。

!!! abstract "假设检验混淆矩阵 (Confusion Matrix)"
    
    | 真实状态 \ 模型预测 | 接受 $H_0$ (Do not reject) | 拒绝 $H_0$ (Reject, Discover) |
    | :--- | :--- | :--- |
    | **$H_0$ 为真** (Null is True) | ✅ 正确推断 (True Negative, $1-\alpha$) | ❌ **第一类错误** (False Positive, $\alpha$) |
    | **$H_1$ 为真** (Alt is True) | ❌ **第二类错误** (False Negative, $\beta$) | ✅ **统计功效** (True Positive, Power: $1-\beta$) |

### 核心概念速查

!!! info "第一类错误 (Type I Error / False Positive)"
    * **定义**：原假设 $H_0$ 为真，却被错误地拒绝了。即弃真错误（假阳性）。
    * **数学表达**：$\alpha = P(\text{Reject } H_0 \mid H_0 \text{ is True}) = \frac{FP}{FP+TN}$
    * **别名**：检验的**显著性水平 (Significance Level)**。在经典统计中，我们通常首先严格控制 $\alpha$（如设定为 $0.05$）。

!!! info "第二类错误 (Type II Error / False Negative)"
    * **定义**：备择假设 $H_1$ 为真，却没有能够拒绝 $H_0$。即取伪错误（假阴性）。
    * **数学表达**：$\beta = P(\text{Accept } H_0 \mid H_1 \text{ is True}) = \frac{FN}{TP+FN}$

!!! success "统计功效 (Statistical Power)"
    * **定义**：备择假设 $H_1$ 为真时，正确拒绝 $H_0$ 的概率。即检验出真实效应的能力。
    * **数学表达**：$\text{Power} = 1 - \beta = P(\text{Reject } H_0 \mid H_1 \text{ is True}) = \frac{TP}{TP+FN}$
    * **直觉**：在给定 $\alpha$ 的前提下，我们希望寻找能够使 Power 最大化的检验方法（即 UMP 检验, Uniformly Most Powerful test）。

!!! info "p-value (p 值)"
    * **定义**：在原假设 $H_0$ 为真的条件下，观察到当前样本统计量（或比其更极端情况）的概率。
    * **避坑**：p 值**绝对不是**“原假设为真的概率”（即 $p \ne P(H_0 \mid \text{Data})$）。它反映的是**数据与原假设的不一致程度**。p 值越小，拒绝 $H_0$ 的理由越充分。

---

### 现代前沿：多重检验与 FDR 控制

在现代高维统计中，我们往往需要同时进行成千上万次检验。此时传统的 $\alpha$ 控制会彻底失效。

??? note "多重检验的灾难 (The Multiple Testing Problem)"
    假设我们独立测试了 $m$ 个完全无效的因子（即 $m$ 个 $H_0$ 均成立），设单次检验的显著性水平为 $\alpha = 0.05$。
    
    那么至少犯一次第一类错误的概率（Family-Wise Error Rate, FWER）为：
    
    \[
    \text{FWER} = P(V \ge 1) = 1 - (1 - \alpha)^m
    \]
    
    当 $m = 100$ 时，$\text{FWER} \approx 0.994$。这意味着只要测试足够多的因子，必然会挖掘出看似显著的“伪信号”。这就引入了 FDR 的概念。

??? success "FDR (False Discovery Rate) 错误发现率"
    由 Benjamini 和 Hochberg 在 1995 年提出，是高维推断的核心。
    
    * **定义**：在所有被拒绝的原假设（即所有声称有显著效应的发现 $R$）中，错误拒绝（即第一类错误 $V$）所占比例的**期望值**。
    
    \[
    \text{FDR} = E\left[ \frac{V}{\max(R, 1)} \right] = \frac{FP}{TP+FP}
    \]
    
    * **BH 过程 (Benjamini-Hochberg Procedure)**：
        1. 将 $m$ 个假设的 p 值从小到大排序：$p_{(1)} \le p_{(2)} \le \dots \le p_{(m)}$。
        2. 找到最大的整数 $k$，使得 $p_{(k)} \le \frac{k}{m} \alpha$。
        3. 拒绝前 $k$ 个原假设（即 $H_{(1)}, \dots, H_{(k)}$）。
        4. **定理**：在独立性（或正相依）假设下，该过程能够严格控制 FDR $\le \alpha$。