---
tags:
  - 组会报告
  - Knockoff
  - FDR控制
  - 统计机器学习
---

# 🎙️ 25秋季组会：Derandomized Knockoff

> **"Randomness in Model-X Knockoffs provides exact FDR control, but Derandomized Knockoffs restore selection stability without sacrificing statistical guarantees."**

!!! abstract "报告综述 (Seminar Abstract)"
    * **主题**：Derandomized Knockoff 方法与理论分析
    * **核心问题**：传统的 Model-X Knockoff 在生成过程中存在随机性，导致特征选择结果存在不稳定性。如何通过去随机化过程获得稳定的特征选择集合？
    * **理论工具**：通过引入 $e$-value 和 $e$-process，揭示了 Knockoff 方法与 e-BH 过程的内在等价性。
    * **主要结论**：通过多次生成 Knockoff 变量并聚合 $e$-value，能够在理论上严格保证 False Discovery Rate (FDR) 控制的前提下，大幅提高特征选择的稳定性。

---

## 1. Model-X 方法回顾

### 1.1 背景引入

在传统的统计与高维特征选择模型中，数据的设定往往会带来不同的处理思路：

* **Fixed-X 方法**：要求设计矩阵 $X$ 是已知的、非随机的，即将其视为一个固定的常数矩阵。
* **Model-X 方法**：现实中的大部分观测数据都是随机生成的。因此，假设 $X$ 作为数据矩阵，服从某个已知的或者未知的联合概率分布。
* 两者的核心目标相似：在保证 **FDR (False Discovery Rate)** 控制的前提下，筛选出与响应变量 $Y$ 真正相关的特征。

### 1.2 方法核心

!!! info "定义 1.1 (Model-X Knockoff 变量)"

    设 $\tilde{X} = (\tilde{X}_1, \tilde{X}_2, \dots, \tilde{X}_p)$ 是基于原始变量 $X = (X_1, X_2, \dots, X_p)$ 进行 Model-X Knockoff 变换生成的新变量族。则其必须满足以下性质：

    * **Swap 性质 (交换不变性)**：对于任意的特征子集 $S \subseteq \{1, \dots, p\}$，在交换 $S$ 中的特征与对应的 Knockoff 变量后，联合分布保持不变：

    $$
    (X, \tilde{X})_{\operatorname{swap}(S)} \stackrel{d}{=} (X, \tilde{X})
    $$

    * **条件独立性**：对于响应变量 $Y$，给定原始特征 $X$ 的条件下，Knockoff 变量 $\tilde{X}$ 应当与 $Y$ 独立：

    $$
    \tilde{X} \perp\!\!\!\perp Y \mid X
    $$

对于特征选择的统计量，通常采用**系数差统计量 (LCD, Lasso Coefficient-Difference)**：

$$
W_j = Z_j - \tilde{Z}_j
$$

其中，特征重要度通过正则化路径获得：

$$
Z_j = |\hat{b}_j(\lambda)|, \quad \tilde{Z}_j = |\hat{b}_{j+p}(\lambda)|
$$

对于 **Model-X Gauss** 设定，典型的显式构造形式为：

$$
\tilde{X} = X(I - \Sigma^{-1} S) + EC
$$

其中，要求 $X \sim N(0, \Sigma)$，且 $S = \text{diag}(s_1, s_2, \dots, s_p)$ 满足 $0 \le S \le 2\Sigma$。误差项矩阵对应的系数 $C$ 必须满足：

$$
CC^T = 2S - S\Sigma^{-1}S
$$

统计量 $W_j$ 满足重要的**硬币反转性质 (Coin-flip property)**，这是保证 FDR 控制的对称性基础。

* （补充）对于多元高斯分布的随机变量，已有成熟的 Knockoff R 包可供使用：

```
Xk <- create.gaussian(X, mu, Sigma, diag_s = diags)
```

### 2. Derandomized Knockoff 方法介绍

### 2.1 背景与动机

在 Model-X 方法中，Knockoff 变量的生成具有内在随机性。这种随机性导致每次运行算法时，最终选择的特征集合都可能发生变化，由此引发了特征选择结果的不稳定性。这种不稳定性可以从 Knockoff 变量的生成公式中直观看出（以高斯分布为例）：$$\widetilde{X} = X\bigl(I - \Sigma^{-1} S\bigr) + E C$$

相比于传统假设检验中的 $p$-value，$e$-value 具有极为优良的组合性质（例如可以直接取平均值或乘积）。因此，文章基于 $e$-value 设计了去随机化 (Derandomization) 的过程。

### 2.2 核心概念定义

!!! info "定义 2.1 (p-value)"
    * **定义**：在原假设 $H_0$ 成立的条件下，观察到检验统计量等于当前数据或更加极端的概率。
    * 在 $H_0$ 成立且检验统计量为连续型时，$p$-value 服从均匀分布 $U[0,1]$。
    * **决策准则**：如果 $p$-value 小于预设的显著性水平 $\alpha$，则拒绝原假设。

!!! info "定义 2.2 (e-value 与 e-process)"
    
    * **定义 (e-value)**：在带滤子 (filtration) 的概率空间上，如果非负随机变量 $E \ge 0$ (a.s.) 满足在 $H_0$ 成立时的期望不超过 1，即：
    
    $$
    \mathbb{E}_{H_0}[E] \le 1
    $$
    
    则称随机变量 $E$ 是一个 **e-value**。

    * **常见的 e-value 构造**：
        1. 基于 p-value 的构造：$E = \frac{1(p \le \alpha)}{\alpha}$。
        2. 似然比 (Likelihood ratio)：$E = \frac{f_1(x)}{f_0(x)}$，其中 $f_0, f_1$ 分别为 $H_0$ 和 $H_1$ 下的密度函数。

    * **e-process**：在序贯检验中，构造一系列随机变量 $E_t \ge 0$，且 $E_t$ 是适应于滤子 $\mathcal{F}_t$ 的。如果它在 $H_0$ 下构成一个**上鞅 (supermartingale)**，即：
    
    $$
    \mathbb{E}_{H_0}(E_t \mid \mathcal{F}_{t-1}) \le E_{t-1}, \quad E_0 \le 1
    $$
    
    则称其为 e-process。此时可以通过 Optional Stopping 定理来控制错误率。

基于上述 e-BH 的性质，我们现在可以将多次生成的 Knockoff e-values 进行安全聚合，而不会破坏 FDR 控制。

### 2.3 主要理论与性质

!!! success "定理 2.1 (等价性定理 THM1)"
    
    Knockoff 方法与基于上述 $e$-value 的 e-BH 方法是等价的。即它们所选择的特征集合完全一致：

    $$
    S_{kn} = S_{ebh}
    $$

??? proof "THM 1 证明思路 (点击展开)"

    基于 $e_j$ 的构造，对于任意 $j \in S_{kn}$，由于 $W_j \ge T$，其指示函数值为 1，可得：
    $e_j = \frac{p}{1 + \sum_{k \in [p]} \mathbb{1}(W_k \le -T)} \ge \frac{p}{\alpha \hat{k}}$。
    因此，如果有特征在 $S_{kn}$ 中，它也必定满足 e-BH 的截断条件，即 $j \in S_{ebh}$。反之，若 $j \notin S_{kn}$，则 $e_j = 0$，永远不会被选中。


!!! success "定理 2.2 (e-BH 过程的 FDR 控制 THM2)"

    假设任意生成的 $e_1, e_2, \dots, e_p$ 满足原假设下的期望和边界约束：

    $$
    \sum_{j \in \mathcal{H}_0} \mathbb{E}[e_j] \le p
    $$

    那么 e-BH 程序所选择的特征集合 $S_{ebh}$ 满足严格的 FDR 控制：

    $$
    \text{FDR} < \alpha
    $$

    *(注：其中 $\text{FDR} = \mathbb{E}\left[ \frac{\sum_{j \in \mathcal{H}_0} \mathbb{1}\left(j \in S_{ebh}\right)}{|\mathcal{S}_{\text{ebh}}| \vee 1} \right]$)*

??? proof "THM 2 证明思路 (点击展开)"

    对于任意原假设成立的特征 $j \in H_0$，若 $j \in S_{ebh}$，则必然满足 $e_j \ge \frac{p}{\alpha |S_{ebh}|}$。
    因此，我们可以利用指示函数的放缩：
    $\mathbb{1}\{j \in S_{ebh}\} = \mathbb{1}\left\{ e_j \ge \frac{p}{\alpha |S_{ebh}|} \right\} \le \frac{\alpha |S_{ebh}| e_j}{p}$。
    将其代入 FDR 的定义求期望，项 $|S_{ebh}|$ 被约去，最终得证 $\text{FDR} \le \alpha$。


基于上述 e-BH 的性质，我们现在可以将多次生成的 Knockoff e-values 进行安全聚合，而不会破坏 FDR 控制。

!!! success "定理 2.3 (Derandomized Knockoff 的 FDR 控制 THM3)"
    
    对于任意的 $\alpha_{\text{kn}}, \alpha_{\text{ebh}} \in (0, 1)$，以及任意给定的 Knockoff 生成次数 $M \ge 1$，由 Derandomized 算法计算得到的最终选择集 $\mathcal{S}_{\text{kn-derand}}$ 始终满足：

    $$
    \text{FDR} \le \alpha_{\text{ebh}}
    $$

    * **强大数律下的渐进性质**：当聚合次数 $M \to \infty$ 时，平均 $e$-value 会几乎处处 (a.s.) 收敛到其条件期望：

    $$
    e_j^{\text{avg}} \xrightarrow[M \to \infty]{\text{a.s.}} e_j^{\infty} := \mathbb{E}[e_j^{(1)} \mid X, Y]
    $$

### 2.4 去随机化 Knockoff 算法流程

利用上述理论支撑，算法的核心在于生成 $M$ 次不同的 Knockoff 副本，转化为对应的 $e$-value 后求平均，最后代入 e-BH 过程中进行筛选。

!!! abstract "算法：Derandomized Knockoff Procedure"

    **输入**: 设计矩阵 $X$，响应变量 $Y$；目标 e-BH FDR 控制水平 $\alpha_{ebh}$，基础 Knockoff 截断水平 $\alpha_{kn}$；迭代次数 $M$。

    **步骤 1 到 4 (生成与转换)**：
    对于 $m = 1, 2, \dots, M$ 循环执行：

    * 采样生成第 $m$ 个 Knockoff 矩阵 $\tilde{X}^{(m)}$。
    
    * 计算特征统计量：

    $$
    W^{(m)} = W([X, \tilde{X}^{(m)}], Y)
    $$

    * 根据 $\alpha_{kn}$ 计算第 $m$ 次的阈值 $T^{(m)}$：

    $$
    T^{(m)} = \inf\left\{ t > 0 : \frac{1 + \sum_{j} \mathbb{1}(W_j^{(m)} \le -t)}{\sum_{j} \mathbb{1}(W_j^{(m)} \ge t)} \le \alpha_{kn} \right\}
    $$

    * 对于所有特征 $j \in [p]$，计算单次的 $e$-value：

    $$
    e_j^{(m)} = p \cdot \frac{\mathbb{1}(W_j^{(m)} \ge T^{(m)})}{1 + \sum_{k=1}^p \mathbb{1}(W_k^{(m)} \le -T^{(m)})}
    $$

    **步骤 5 (聚合)**：
    结束循环后，计算每个特征的平均 $e$-value：

    $$
    e_j^{avg} = \frac{1}{M}\sum_{m=1}^{M}e_j^{(m)}, \quad \forall j \in [p]
    $$

    **步骤 6 (e-BH 截断)**：
    寻找截断秩 $\hat{k}$：

    $$
    \hat{k} = \max\left\{ k \in [p] : e_{(k)}^{avg} \ge \frac{p}{\alpha_{ebh} k} \right\}
    $$

    **输出**: 最终的去随机化特征选择集合：

    $$
    S_{kn-derand} := \left\{ j \in [p] : e_j^{avg} \ge \frac{p}{\alpha_{ebh} \hat{k}} \right\}
    $$

    * **参数选择推荐**：在实际应用中，原论文推导并建议使用以下参数比例，在定理 3 的保证下能够达到最优的势 (Power) 与 FDR 平衡：

    $$
    \alpha_{ebh} = 2\alpha_{kn}
    $$


## 3. 扩展：Copula-Model-X Knockoff

在实际应用中，如果数据的联合分布 $F(x_1, \dots, x_p)$ 无法直接得知为高斯分布，而是只知道各个边缘分布 $F_j(x_j)$，我们可以利用 **Copula 模型** 将问题转化为高斯 Knockoff 来处理。

**实现步骤**：

1. **边缘分布估计**：估计每个特征的经验 CDF，记为 $\hat{F}_j$。

2. **正态分位数变换 (Normal Quantile Transform)**：通过概率积分变换，将原始数据投影到 $[0, 1]$ 空间，再映射到潜变量的高斯空间：
   $$
   \hat{U}_{ij} = \hat{F}_j(X_{ij}), \quad \tilde{Z}_{ij} = \Phi^{-1}(\hat{U}_{ij})
   $$

3. **生成潜变量 Knockoff**：在转换后的高斯空间中，生成对应于潜变量 $\tilde{Z}$ 的高斯 Knockoff 变量。

4. **逆变换回原空间**：最后通过逆变换将 Knockoff 变量映射回原特征的数据尺度：
   $$
   U_{ij}^{\tilde{Z}} = \Phi(\tilde{Z}_{ij}), \quad \tilde{X}_{ij} = \hat{F}_j^{-1}(U_{ij}^{\tilde{Z}})
   $$