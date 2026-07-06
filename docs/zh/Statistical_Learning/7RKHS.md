# 第七章：再生核希尔伯特空间（RKHS）与核脊回归

在第六章支持向量机中，我们通过核技巧非正式地引入了高维特征映射。本章将从泛函分析的角度，严格建立**再生核希尔伯特空间（Reproducing Kernel Hilbert Space, RKHS）**的数学基石。RKHS 构成了一个足够灵活的无限维函数空间，同时其优良的拓扑性质又确保了其中无限维优化问题的解可以退化到有限维空间中求解，这也是核方法在统计学习中取得成功的核心理论所在。

---

## 1. 希尔伯特空间与评估泛函

在深入再生核之前，我们需要先明确函数空间的完备性与线性泛函的连续性。

!!! info "定义 1.1 (希尔伯特空间 Hilbert Space)"

    设 $\mathcal{H}$ 是一个定义了内积 $\langle \cdot, \cdot \rangle_{\mathcal{H}}$ 的线性空间。由该内积导出的范数为 $\|f\|_{\mathcal{H}} = \sqrt{\langle f, f \rangle_{\mathcal{H}}}$。若 $\mathcal{H}$ 在该范数拓扑下是完备的（即其中任意 Cauchy 序列均收敛到空间中的元素），则称 $\mathcal{H}$ 为**希尔伯特空间**。

### 1.1 狄拉克评估泛函（Dirac Evaluation Functional）

设 $\mathcal{H}$ 是定义在集合 $\mathcal{X}$ 上的实值函数空间。对于每一个固定的点 $x \in \mathcal{X}$，我们可以定义一个映射 $L_x: \mathcal{H} \to \mathbb{R}$，其作用是提取函数 $f$ 在点 $x$ 处的取值。

!!! info "定义 1.2 (评估泛函)"

    线性算子 $L_x: \mathcal{H} \to \mathbb{R}$ 称为点 $x$ 处的**狄拉克评估泛函（Evaluation Functional）**，其定义为：

    \[
    L_x(f) = f(x), \quad \forall f \in \mathcal{H}
    \]

在许多传统的希尔伯特空间（例如平方可积空间 $L_2$）中，点评估泛函甚至是没有意义或不连续的（因为改变单点的值不影响 $L_2$ 内积积分）。而 RKHS 的核心特征正是要求点评估泛函必须是连续的。

---

## 2. 再生核希尔伯特空间的定义与再生性

!!! info "定义 2.1 (再生核希尔伯特空间 RKHS)"

    一个函数希尔伯特空间 $\mathcal{H}$ 被称为**再生核希尔伯特空间（RKHS）**，当且仅当对于任意的 $x \in \mathcal{X}$，其对应的点评估泛函 $L_x$ 在 $\mathcal{H}$ 上都是**有界（连续）线性算子**。也就是说，存在一个常数 $M_x > 0$，使得：

    \[
    |L_x(f)| = |f(x)| \le M_x \|f\|_{\mathcal{H}}, \quad \forall f \in \mathcal{H}
    \]

### 2.1 再生核（Reproducing Kernel）的引入

根据泛函分析中的 **Riesz 表示定理（Riesz Representation Theorem）**，定义在希尔伯特空间上的任何连续线性泛函，都可以由该空间中一个唯一的元素通过内积的形式表现出来。

由于 RKHS 中的 $L_x$ 连续，因此必然存在一个唯一的函数（记为 $K_x \in \mathcal{H}$），满足：

\[
L_x(f) = \langle f, K_x \rangle_{\mathcal{H}}, \quad \forall f \in \mathcal{H}
    \]

由于 $K_x$ 本身是定义在 $\mathcal{X}$ 上的函数，我们可以将两个自变量统一表示。定义二元函数 $K: \mathcal{X} \times \mathcal{X} \to \mathbb{R}$ 为：

\[
K(x, z) = K_x(z)
\]

!!! note "定理 2.1 (核的再生性 Reproducing Property)"

    根据上述构造，二元函数 $K(\cdot, \cdot)$ 满足以下两项决定性的**再生性质**：

    1. 对于任意固定的 $x \in \mathcal{X}$，函数 $K(x, \cdot)$ 属于空间 $\mathcal{H}$。
    2. 对于任意的 $x \in \mathcal{X}$ 和 $f \in \mathcal{H}$，满足：

    \[
    \langle f, K(x, \cdot) \rangle_{\mathcal{H}} = f(x)
    \]

    特别地，将 $f$ 替换为 $K(z, \cdot)$ 得到：

    \[
    \langle K(x, \cdot), K(z, \cdot) \rangle_{\mathcal{H}} = K(x, z)
    \]

---

## 3. 代表定理（Representer Theorem）

在统计学习中，我们在无限维函数空间 $\mathcal{H}$ 中寻找最优函数时，面临的是一个无限维度的优化问题。代表定理奠定了核方法在计算上的可行性。

!!! note "定理 3.1 (代表定理)"

    设 $\mathcal{H}$ 是一个由再生核 $K$ 导出的 RKHS。给定训练集 $\mathcal{D}_n = \{(x_i, y_i)\}_{i=1}^n$。考虑如下包含任意经验损失函数 $L$ 以及单调递增正则化项 $\Omega(\|f\|_{\mathcal{H}})$ 的通用正则化经验风险最小化问题：

    \[
    \min_{f \in \mathcal{H}} \left\{ \frac{1}{n} \sum_{i=1}^n L(y_i, f(x_i)) + \Omega(\|f\|_{\mathcal{H}}^2) \right\}
    \]

    则该无限维优化问题的任意全局最优解 $f^*$ 必定可以精确地表达为训练样本点处核函数的有限线性组合：

    \[
    f^*(x) = \sum_{i=1}^n \alpha_i K(x_i, x)
    \]

    其中 $\alpha_i \in \mathbb{R}$ ($i=1,\dots,n$) 为有限个待求的实数标量系数。

??? proof "证明：基于正交分解的代表定理证明"

    我们可以将整个希尔伯特空间 $\mathcal{H}$ 沿着由当前观测样本点构成的子空间进行正交分解。

    定义由所有训练样本点的核映射张成的有限维线性子空间 $\mathcal{H}_0$ 如下：

    \[
    \mathcal{H}_0 = \text{span}\{K(x_1, \cdot), K(x_2, \cdot), \dots, K(x_n, \cdot)\}
    \]

    则空间中的任意函数 $f \in \mathcal{H}$ 均可唯一地分解为子空间 $\mathcal{H}_0$ 中的成分与正交补空间 $\mathcal{H}_0^{\perp}$ 中的成分之和：

    \[
    f = f_0 + f_{\perp}
    \]

    其中 $f_0 \in \mathcal{H}_0$ 且 $f_{\perp} \in \mathcal{H}_0^{\perp}$。这意味着对于任意的 $i = 1, \dots, n$，有：

    \[
    \langle f_{\perp}, K(x_i, \cdot) \rangle_{\mathcal{H}} = 0
    \]

    接下来，我们考察函数 $f$ 在任意训练数据点 $x_i$ 处的取值。根据核的再生性：

    \[
    f(x_i) = \langle f, K(x_i, \cdot) \rangle_{\mathcal{H}} = \langle f_0 + f_{\perp}, K(x_i, \cdot) \rangle_{\mathcal{H}}
    \]

    \[
    f(x_i) = \langle f_0, K(x_i, \cdot) \rangle_{\mathcal{H}} + \langle f_{\perp}, K(x_i, \cdot) \rangle_{\mathcal{H}} = \langle f_0, K(x_i, \cdot) \rangle_{\mathcal{H}} = f_0(x_i)
    \]

    由此可见，垂直分量 $f_{\perp}$ 对函数在训练集上的预测取值没有任何贡献，因此改变 $f_{\perp}$ 完全不会改变损失函数项 $\frac{1}{n} \sum_{i=1}^n L(y_i, f(x_i))$ 的大小。

    最后，我们考察正则化项。根据 Pythagorean 定理以及正交性：

    \[
    \|f\|_{\mathcal{H}}^2 = \|f_0 + f_{\perp}\|_{\mathcal{H}}^2 = \|f_0\|_{\mathcal{H}}^2 + \|f_{\perp}\|_{\mathcal{H}}^2 \ge \|f_0\|_{\mathcal{H}}^2
    \]

    由于正则化惩罚函数 $\Omega$ 是单调递增的，所以：

    \[
    \Omega(\|f\|_{\mathcal{H}}^2) \ge \Omega(\|f_0\|_{\mathcal{H}}^2)
    \]

    为了使整体目标函数最小化，我们必须强制令非必要的残差分量 $f_{\perp} = 0$。因此，最优解 $f^*$ 必定落在有限维子空间 $\mathcal{H}_0$ 中，即可以表示为基底的线性组合：

    \[
    f^*(\cdot) = f_0(\cdot) = \sum_{i=1}^n \alpha_i K(x_i, \cdot)
    \]

    代入自变量 $x$，即证：

    \[
    f^*(x) = \sum_{i=1}^n \alpha_i K(x_i, x)
    \]

---

## 4. 核脊回归（Kernel Ridge Regression）

核脊回归（Kernel Ridge Regression, KRR）是代表定理的一个经典应用，它将传统的线性脊回归通过平方损失无缝推广到非线性空间。

!!! info "定义 4.1 (核脊回归泛函优化问题)"

    给定样本，核脊回归的目标是在 RKHS 中寻找一个函数 $f$，以最小化带有 $L_2$ 正则化惩罚项的平方损失：

    \[
    \min_{f \in \mathcal{H}} \left\{ \frac{1}{n} \sum_{i=1}^n (y_i - f(x_i))^2 + \lambda \|f\|_{\mathcal{H}}^2 \right\}
    \]

### 4.1 矩阵格式的有限维转化

根据代表定理，我们代入 $f(x) = \sum_{j=1}^n \alpha_j K(x_j, x)$。定义核矩阵（Gram Matrix）$K \in \mathbb{R}^{n \times n}$ 满足 $K_{ij} = K(x_i, x_j)$，以及系数向量 $\alpha = (\alpha_1, \dots, \alpha_n)^{\top}$。

我们分别重写目标函数的两部分：

1. 拟合预测向量为 $\hat{Y} = K\alpha$。因此经验平方损失项变为：

\[
\frac{1}{n} \|Y - K\alpha\|_2^2
\]

2. 根据内积的双线性性质，计算函数范数平方：

\[
\|f\|_{\mathcal{H}}^2 = \left\langle \sum_{i=1}^n \alpha_i K(x_i, \cdot), \sum_{j=1}^n \alpha_j K(x_j, \cdot) \right\rangle_{\mathcal{H}} = \sum_{i=1}^n \sum_{j=1}^n \alpha_i \alpha_j \langle K(x_i, \cdot), K(x_j, \cdot) \rangle_{\mathcal{H}} = \alpha^{\top}K\alpha
\]

结合以上两项，核脊回归的有穷参数目标函数表示为：

\[
S(\alpha) = \frac{1}{n} (Y - K\alpha)^{\top}(Y - K\alpha) + \lambda \alpha^{\top}K\alpha
\]

!!! note "定理 4.2 (核脊回归的显式系数解)"

    若核矩阵 $K$ 是正定的，则使上述目标函数最小化的系数向量 $\alpha$ 具有如下闭式解：

    \[
    \alpha = (K + n\lambda I_n)^{-1}Y
    \]

??? proof "证明：核脊回归参数解的矩阵微积分推导"

    首先，将矩阵目标函数 $S(\alpha)$ 完整展开：

    \[
    S(\alpha) = \frac{1}{n} \left( Y^{\top}Y - 2\alpha^{\top}K Y + \alpha^{\top}K K \alpha \right) + \lambda \alpha^{\top}K\alpha
    \]

    由于 Gram 矩阵 $K$ 是对称的（$K^{\top} = K$），利用矩阵微分公式对列向量 $\alpha$ 求梯度偏导数：

    \[
    \frac{\partial S(\alpha)}{\partial \alpha} = \frac{1}{n} \left( -2KY + 2KK\alpha \right) + 2\lambda K\alpha
    \]

    令梯度向量等于 0，以求得极值点：

    \[
    -\frac{2}{n}KY + \frac{2}{n}KK\alpha + 2\lambda K\alpha = 0
    \]

    在等式两边同时左乘 $\frac{n}{2}$ 进行常数化简：

    \[
    -KY + KK\alpha + n\lambda K\alpha = 0
    \]

    整理项，将未知系数 $\alpha$ 提取到右边：

    \[
    K(K + n\lambda I_n)\alpha = KY
    \]

    由于对于任意的 $\lambda > 0$，正则化矩阵 $(K + n\lambda I_n)$ 总是严格正定可逆的，我们在两边同时左乘其逆矩阵 $(K + n\lambda I_n)^{-1}$（或在原正规方程两边消去半正定的 $K$ 算子），即证：

    \[
    \alpha = (K + n\lambda I_n)^{-1}Y
    \]

最终，对于任意一个全新输入的测试点 $x_0$，核脊回归给出的非线性预测输出值为：

\[
\hat{f}(x_0) = \sum_{i=1}^n \alpha_i K(x_i, x_0) = K(x_0, \cdot)^{\top}(K + n\lambda I_n)^{-1}Y
\]