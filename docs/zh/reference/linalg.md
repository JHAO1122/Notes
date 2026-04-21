# 🔢 线性代数 (Linear Algebra)

本模块目前聚焦于**线性代数**的相关内容。

## 一、矩阵求导

矩阵求导的方法一直是笔者在学习过程中最头疼的点。在笔者上过的回归分析、多元统计分析以及优化实用算法中都或多或少接触过矩阵求导，但是从未进行过系统地整理。因在*优化实用算法*期末考试中深陷繁琐的求和指标而痛失 15 分后，笔者痛定思痛，系统整理矩阵求导的方法，并使用更加简单的矩阵微分方法。


!!! warning "布局约定 (Layout Convention)"
    本速查表统一采用机器学习与优化领域最常用的**分母布局 (Denominator Layout)**。
    即：标量 $f$ 对矩阵 $X_{m \times n}$ 的导数 $\frac{\partial f}{\partial X}$，其维度与 $X$ 保持一致，同为 $m \times n$。

---

### （一）核心方法：矩阵微分法

传统的对元素求偏导非常容易在求和符号中迷失。更优雅的做法是使用**全微分 (Differential)**。

!!! info "微分与导数的核心联系"
    标量函数 $f(X)$ 的全微分 $df$ 总是可以表示为某个矩阵与 $dX$ 乘积的迹 (Trace)：
    
    $$
    df = \text{tr}\left( \left(\frac{\partial f}{\partial X}\right)^T dX \right)
    $$
    
    **四步法：**

    1. 求全微分 $df$。

    2. 利用迹的性质化简，凑成 $\text{tr}(M^T dX)$ 的形式。
    
    3. 扒掉 $\text{tr}$ 和 $dX$，剩下的 $M$ 就是导数 $\frac{\partial f}{\partial X}$。

??? note "矩阵微分的基本运算法则"
    * **加减法**：$d(X \pm Y) = dX \pm dY$
    * **乘法**：$d(XY) = (dX)Y + X(dY)$
    * **转置**：$d(X^T) = (dX)^T$
    * **迹**：$d(\text{tr}(X)) = \text{tr}(dX)$
    * **迹的轮换对称性**：$\text{tr}(ABC) = \text{tr}(CAB) = \text{tr}(BCA)$
    * **迹的转置不变性**：$\text{tr}(A) = \text{tr}(A^T)$

---

### （二） 迹的求导 (Trace Derivatives)

在推导 OLS (最小二乘法) 或构建损失函数时，经常需要对迹求导。

!!! success "常见迹的导数公式"
    设 $X$ 为变量矩阵，$A, B$ 为常量矩阵：
    
    1. $$\frac{\partial \text{tr}(AX)}{\partial X} = A^T$$

    2. $$\frac{\partial \text{tr}(AXB)}{\partial X} = A^T B^T$$

    3. $$\frac{\partial \text{tr}(X^TAX)}{\partial X} = (A + A^T)X$$
       
    *(特例：若 $A$ 是对称矩阵，则结果为 $2AX$)*


??? note "推导过程：$\frac{\partial \text{tr}(AX)}{\partial X}$"
    令 $f = \text{tr}(AX)$，求微分：
    
    $$
    df = d(\text{tr}(AX)) = \text{tr}(d(AX)) = \text{tr}(A dX)
    $$
    
    为了凑成 $\text{tr}(M^T dX)$ 的形式，由于 $\text{tr}(A dX) = \text{tr}((A^T)^T dX)$，所以 $M^T = A \implies M = A^T$。
    
    故：$\frac{\partial f}{\partial X} = A^T$。

??? note "推导过程：$\frac{\partial \text{tr}(X^TAX)}{\partial X}$"
    令 $f = \text{tr}(X^TAX)$，求微分（利用乘法法则 $d(UV) = dU V + U dV$）：
    
    $$
    df = \text{tr}(d(X^T) AX + X^T d(AX)) = \text{tr}((dX)^T AX) + \text{tr}(X^T A dX)
    $$
    
    利用迹的转置不变性 $\text{tr}(M) = \text{tr}(M^T)$，处理第一项：
    
    $$
    \text{tr}((dX)^T AX) = \text{tr}(((dX)^T AX)^T) = \text{tr}(X^T A^T dX)
    $$
    
    将两项合并：
    
    $$
    df = \text{tr}(X^T A^T dX) + \text{tr}(X^T A dX) = \text{tr}((X^T A^T + X^T A) dX)
    $$
    
    提取公因式并利用迹的轮换对称性凑成 $\text{tr}(M^T dX)$：
    
    $$
    df = \text{tr}( (X^T(A^T + A)) dX ) = \text{tr}( ((A^T + A)^T X)^T dX )
    $$
    
    故 $M = (A^T + A)^T X = (A + A^T)X$。

---

### （三）行列式与逆矩阵的求导 (Determinant & Inverse)

在极大似然估计 (MLE) 中，如果涉及多元正态分布的协方差矩阵 $\Sigma$，必考行列式和逆的求导。

!!! success "逆矩阵与行列式的导数公式"
    设 $X$ 为可逆方阵：
    
    1. **逆矩阵微分**： $d(X^{-1}) = -X^{-1} (dX) X^{-1}$

    2. **行列式微分 (Jacobi's Formula)**： $d|X| = |X| \text{tr}(X^{-1} dX)$
    
    3. **行列式的导数**： $\frac{\partial |X|}{\partial X} = |X| (X^{-1})^T$
    
    4. **对数行列式的导数**：$\frac{\partial \ln |X|}{\partial X} = (X^{-1})^T$

??? note "推导过程：逆矩阵的微分 $d(X^{-1})$"
    从恒等式 $XX^{-1} = I$ 出发，两边同时求微分：
    
    $$
    d(XX^{-1}) = d(I) \implies (dX)X^{-1} + X d(X^{-1}) = 0
    $$
    
    移项并左乘 $X^{-1}$：
    
    $$
    X d(X^{-1}) = -(dX)X^{-1} \implies d(X^{-1}) = -X^{-1}(dX)X^{-1}
    $$

??? note "行列式求导 $d|X|$（代数余子式求和法）"
    本推导旨在打通标量求和、伴随矩阵与矩阵迹（Trace）之间的联系。

    1. 元素级求导：$\frac{\partial |X|}{\partial x_{ij}} = C_{ij}$
    根据行列式的 **Laplace 展开**，沿第 $i$ 行展开有：

    $$|X| = \sum_{k=1}^n x_{ik} C_{ik}$$
    
    由于代数余子式 $C_{ik}$ 的计算不包含第 $i$ 行的任何元素，因此对 $x_{ij}$ 而言，$C_{ik}$ 均为常数。求导得：

    $$\frac{\partial |X|}{\partial x_{ij}} = \frac{\partial}{\partial x_{ij}} (x_{i1}C_{i1} + \dots + x_{ij}C_{ij} + \dots + x_{in}C_{in}) = C_{ij}$$
    
    **结论**：行列式对某个元素的偏导数等于该位置的代数余子式。

    2. 微分与迹的转换：$d|X| = \text{tr}(X^* dX)$
    利用全微分定义与伴随矩阵 $X^*$ 的定义（$(X^*)_{ji} = C_{ij}$）：
    
    $$
    d|X| = \sum_{i=1}^n \sum_{j=1}^n \frac{\partial |X|}{\partial x_{ij}} dx_{ij} = \sum_{i=1}^n \sum_{j=1}^n C_{ij} dx_{ij}
    $$
    
    将其写为矩阵乘积的对角线元素求和形式：
    
    $$
    d|X| = \sum_{j=1}^n \left( \sum_{i=1}^n (X^*)_{ji} (dX)_{ij} \right) = \text{tr}(X^* dX)
    $$

    3. 导出 Jacobi 公式：$d|X| = |X| \text{tr}(X^{-1} dX)$
    利用逆矩阵与伴随矩阵的经典关系 $X^* = |X| X^{-1}$，代入上式：
    
    $$d|X| = \text{tr}(|X| X^{-1} dX) = |X| \text{tr}(X^{-1} dX)$$

    这便是 **Jacobi's Formula**，它将行列式的变化率与矩阵的逆和迹联系了起来。

    4. 对数行列式的梯度 $\frac{\partial \ln |X|}{\partial X}$
    对 $f = \ln |X|$ 求全微分，并代入 Jacobi 公式：
    
    $$df = d(\ln |X|) = \frac{1}{|X|} d|X| = \text{tr}(X^{-1} dX)$$

    对比矩阵微分标准型 $df = \text{tr}\left( \left(\frac{\partial f}{\partial X}\right)^T dX \right)$：

    $$\left( \frac{\partial \ln |X|}{\partial X} \right)^T = X^{-1} \implies \frac{\partial \ln |X|}{\partial X} = (X^{-1})^T$$

    *(注：若 $X$ 为对称矩阵，则结果简化为 $X^{-1}$)*


??? note "行列式求导 $d|X|$（矩阵微分法）"
    本推导完全脱离元素级的 Laplace 展开，仅利用行列式的乘法性质与迹的线性性质。

    核心引理（单位矩阵处的微元）：
    对于单位矩阵 $I$ 的一个无穷小扰动 $dE$，其行列式满足：
    
    $$|I + dE| = 1 + \text{tr}(dE) + o(\|dE\|)$$
    
    **直觉理解**：若 $dE$ 的特征值为 $\lambda_i$，则 $|I+dE| = \prod (1+\lambda_i) = 1 + \sum \lambda_i + \dots = 1 + \text{tr}(dE)$。

    一般矩阵 $X$ 的全微分推导：
    我们要计算 $d|X| = |X + dX| - |X|$。利用行列式的乘法性质，强制提取 $X$ 构造单位矩阵形式：
    
    $$
    \begin{aligned}
    |X + dX| &= |X(I + X^{-1}dX)| \\
    &= |X| \cdot |I + X^{-1}dX|
    \end{aligned}
    $$
    
    令微元 $dE = X^{-1}dX$，代入上述引理：
    
    $$
    \begin{aligned}
    |X + dX| &= |X| \cdot (1 + \text{tr}(X^{-1}dX)) \\
    &= |X| + |X|\text{tr}(X^{-1}dX)
    \end{aligned}
    $$
    
    由微分定义 $d|X| = |X + dX| - |X|$，立即得到：
    
    $$d|X| = |X|\text{tr}(X^{-1}dX)$$


??? note "推导过程：行列式的导数 $\frac{\partial |X|}{\partial X}$"
    令 $f = |X|$。根据雅可比公式 (Jacobi's Formula)，行列式的微分为：
    
    $$
    df = |X| \text{tr}(X^{-1} dX)
    $$
    
    将上式重写为迹的内积形式：
    
    $$
    df = \text{tr}(|X| X^{-1} dX) = \text{tr}( (|X| (X^{-1})^T)^T dX )
    $$
    
    对照 $df = \text{tr}(M^T dX)$，即可得到：
    
    $$
    \frac{\partial |X|}{\partial X} = |X| (X^{-1})^T
    $$

---

### （四）特征值的求导 (Eigenvalue Derivatives)

在推导主成分分析 (PCA) 或研究协方差矩阵的扰动时，极其重要。

!!! success "特征值的导数公式"
    设方阵 $X$ 有非零特征值 $\lambda$，其对应的**右特征向量**为 $u$ ($Xu = \lambda u$)，**左特征向量**为 $v$ ($v^T X = \lambda v^T$)。
    
    假设左右特征向量已经归一化满足 $v^T u = 1$，则特征值 $\lambda$ 对矩阵 $X$ 的导数为：
    
    $$
    \frac{\partial \lambda}{\partial X} = v u^T
    $$
    
    *(注意：若 $X$ 是实对称矩阵，则左、右特征向量相等，即 $u=v$，此时导数为 $uu^T$)*

??? note "推导过程：$\frac{\partial \lambda}{\partial X}$"
    从右特征向量的定义式出发：
    
    $$
    Xu = \lambda u
    $$
    
    对两边同时求全微分：
    
    $$
    (dX)u + X(du) = (d\lambda)u + \lambda(du)
    $$
    
    在等式两边同时**左乘**左特征向量 $v^T$：
    
    $$
    v^T(dX)u + v^TX(du) = v^T(d\lambda)u + v^T\lambda(du)
    $$
    
    因为 $v^T$ 是左特征向量，满足 $v^TX = \lambda v^T$，代入上式左边第二项：
    
    $$
    v^T(dX)u + \lambda v^T(du) = d\lambda (v^Tu) + \lambda v^T(du)
    $$
    
    等式两边的 $\lambda v^T(du)$ 刚好抵消！我们得到：
    
    $$
    v^T(dX)u = d\lambda (v^Tu) \implies d\lambda = \frac{v^T(dX)u}{v^Tu}
    $$
    
    由于标量的迹等于其自身，我们可以加上 $\text{tr}$ 并利用轮换对称性：
    
    $$
    d\lambda = \text{tr}\left( \frac{v^T(dX)u}{v^Tu} \right) = \text{tr}\left( \frac{u v^T}{v^Tu} dX \right)
    $$
    
    由于 $v^T u$ 只是一个标量，若将其提取出来并与标准形式 $df = \text{tr}(M^T dX)$ 对比：
    
    $$
    M^T = \frac{u v^T}{v^Tu} \implies M = \frac{v u^T}{v^Tu}
    $$
    
    如果事先约定 $v^Tu = 1$，则最终结果简化为：
    
    $$
    \frac{\partial \lambda}{\partial X} = v u^T
    $$