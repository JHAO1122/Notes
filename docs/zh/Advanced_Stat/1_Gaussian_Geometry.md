# 第一章：高斯几何与经典抽样分布

本章从高斯向量的投影几何出发，统一推导正态分布、卡方分布、$t$ 分布和 $F$ 分布等经典抽样分布。核心思想是：**线性对比产生正态分布，投影向量的平方长度产生卡方分布，用独立估计的尺度学生化后产生 $t$ 分布，而两个独立均方之比产生 $F$ 分布。**

当投影后的均值不为零时，相应的卡方、$t$ 或 $F$ 分布会变为非中心分布。因此，非中心分布自然地描述固定备择假设，并用于功效分析与样本量计算。

---

## 1. 一页概览统计问题

统计模型是一个分布族

\[
\mathcal{P}=\{P_\theta:\theta\in\Theta\},
\]

它描述观测数据 $X$ 可能服从的分布。参数 $\theta$ 是固定但未知的，而样本 $X$ 在被观测之前是随机的。**统计量（statistic）** 是 $X$ 的可测函数 $T(X)$，其中不能包含未知参数。

同一个统计量可以服务于不同的统计推断任务：

* **估计量（estimator）** $\widehat{g}(X)$ 用于逼近目标 $g(\theta)$；
* **检验（test）** $\phi(X)\in[0,1]$ 衡量反对原假设的证据，其中 $E_\theta\phi(X)$ 等于拒绝概率；
* **置信集（confidence set）** $C(X)$ 具有重复抽样覆盖率

\[
P_\theta\{\theta\in C(X)\}\geq 1-\alpha.
\]

在每一种情形中，统计推断都要求我们知道统计量在 $P_\theta$ 下的分布，这就是统计量的**抽样分布（sampling distribution）**。

!!! example "例 1（贯穿本章的模型）"

    若 $X_1,\ldots,X_n\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2)$，则 $\theta=(\mu,\sigma^2)$。统计量 $\overline{X}$ 和

    \[
    S^2=\frac{1}{n-1}\sum_{i=1}^n(X_i-\overline{X})^2
    \]

    分别估计 $\mu$ 和 $\sigma^2$。

    关于 $\mu$ 的精确推断之所以可行，是因为标准化后的 $\overline{X}$ 服从正态分布，残差平方和服从卡方分布，并且二者相互独立。$t$ 分布恰好就是二者学生化比值的分布。

---

## 2. 为什么这些分布反复出现

正态、卡方、$t$ 和 $F$ 分布并不是一组需要孤立记忆的分布。在高斯模型中，它们分别对应四种统计操作：

\[
\begin{aligned}
\text{取线性对比} &\longrightarrow \text{正态分布},\\
\text{取投影向量的平方长度} &\longrightarrow \text{卡方分布},\\
\frac{\text{正态对比}}{\text{独立估计的尺度}} &\longrightarrow t\text{ 分布},\\
\text{比较两个独立均方} &\longrightarrow F\text{ 分布}.
\end{aligned}
\]

| 分布 | 结构形式 | 典型精确统计量 | 主要用途 |
| --- | --- | --- | --- |
| 正态分布 | 高斯线性对比 | 正态误差下的样本均值或回归系数 | 估计与已知方差检验 |
| 卡方分布 | 中心高斯投影的平方范数 | 残差平方和除以 $\sigma^2$ | 方差推断 |
| 非中心卡方分布 | 保留非零均值分量的平方范数 | 备择假设下的回归平方和 | 信号强度与功效 |
| $t$ 分布 | 中心正态变量除以独立的卡方尺度 | 原假设下学生化的均值或系数 | 未知方差推断 |
| 非中心 $t$ 分布 | 平移正态变量除以独立的卡方尺度 | 备择假设下学生化的均值或系数 | 功效与样本量 |
| $F$ 分布 | 两个独立中心均方之比 | 原假设下的 ANOVA 或嵌套模型统计量 | 多个约束的检验 |
| 非中心 $F$ 分布 | 非中心分子除以中心残差均方 | 备择假设下的 ANOVA 或嵌套模型统计量 | 功效与试验设计 |

!!! warning "精确分布与近似分布"

    本章结果在协方差为 $\sigma^2I$ 的正态模型下是有限样本精确的。对于非正态样本，中心极限定理可能使线性统计量近似正态，但估计均值与估计方差之间的精确独立性通常不再成立，此时稳健方法或渐近方法可能更合适。

    因而，我们总要同时追问：**统计量是什么？哪些假设保证了它的参考分布？**

---

## 3. 高斯向量与投影几何

设 $X\sim N_n(\mu,\Sigma)$。对于固定矩阵 $A$ 和固定向量 $b$，有

\[
AX+b\sim N_m(A\mu+b,A\Sigma A^T).
\]

高斯分布在线性变换下的封闭性，是本章所有精确结果的来源。

!!! success "定理 1（高斯向量中零协方差等价于独立）"

    假设 $(U^T,V^T)^T$ 联合高斯，则

    \[
    U\perp V \quad\Longleftrightarrow\quad \operatorname{Cov}(U,V)=0.
    \]

??? proof "定理 1 的证明（点击展开）"

    在二阶矩存在时，独立总能推出互协方差为零。下面证明反方向。

    将联合均值写为 $(\mu_U,\mu_V)$，并假设互协方差为零。对于向量 $s$ 和 $t$，联合特征函数为

    \[
    \begin{aligned}
    E\exp\{i(s^TU+t^TV)\}
    &=\exp\left\{i(s^T\mu_U+t^T\mu_V)
    -\frac{1}{2}
    \begin{pmatrix}s\\t\end{pmatrix}^{T}
    \begin{pmatrix}\Sigma_U&0\\0&\Sigma_V\end{pmatrix}
    \begin{pmatrix}s\\t\end{pmatrix}\right\}\\
    &=\exp\left\{is^T\mu_U-\frac{1}{2}s^T\Sigma_Us\right\}
    \exp\left\{it^T\mu_V-\frac{1}{2}t^T\Sigma_Vt\right\}.
    \end{aligned}
    \]

    最后一个表达式是两个边缘特征函数的乘积。由特征函数的唯一性，$U$ 与 $V$ 相互独立。$\square$

!!! success "定理 2（正交投影的等价刻画）"

    对于 $n\times n$ 实矩阵 $P$，以下三个条件等价：

    1. $P=P^T=P^2$；
    2. 存在子空间 $\mathcal{S}\subseteq\mathbb{R}^n$，使得对每个 $x$，$Px$ 都是 $x$ 在 $\mathcal{S}$ 上的正交投影；
    3. 存在满足 $U^TU=I_r$ 的 $n\times r$ 矩阵 $U$，使得 $P=UU^T$。

    在这些条件下，

    \[
    \mathcal{S}=\operatorname{col}(P)=\operatorname{col}(U),
    \qquad
    \operatorname{rank}(P)=\operatorname{tr}(P)=r,
    \]

    且 $I-P$ 是到 $\mathcal{S}^{\perp}$ 的正交投影。

??? proof "定理 2 的证明（点击展开）"

    先假设条件 1 成立。若 $y=Px\in\operatorname{col}(P)$，则

    \[
    Py=P^2x=Px=y.
    \]

    对任意 $z=Pw\in\operatorname{col}(P)$，有

    \[
    z^T(x-Px)=w^TP^T(I-P)x=w^TP(I-P)x=0.
    \]

    因而 $Px\in\operatorname{col}(P)$，且 $x-Px\perp\operatorname{col}(P)$，条件 2 成立。

    再假设条件 2 成立。选择 $\mathcal{S}$ 的一组标准正交基 $u_1,\ldots,u_r$，并令 $U=(u_1,\ldots,u_r)$。由正交投影公式，

    \[
    Px=\sum_{j=1}^ru_ju_j^Tx=UU^Tx,
    \]

    所以条件 3 成立。

    最后，若条件 3 成立，则

    \[
    P^T=UU^T=P,
    \qquad
    P^2=U(U^TU)U^T=P,
    \]

    因而条件 1 成立。

    幂等矩阵的每个特征值都满足 $\lambda^2=\lambda$，所以特征值只能是 $0$ 或 $1$。又因为 $P$ 对称，所以 $P$ 可正交对角化，秩和迹都等于特征值 $1$ 的个数。关于 $I-P$ 的结论可直接得到。$\square$

!!! info "引理 3（嵌套投影）"

    设 $P_0$ 和 $P_1$ 是正交投影。若 $\operatorname{col}(P_0)\subseteq\operatorname{col}(P_1)$，则

    \[
    P_1P_0=P_0P_1=P_0.
    \]

    因而 $P_1-P_0$ 是到 $\operatorname{col}(P_1)\cap\operatorname{col}(P_0)^\perp$ 的正交投影，并且它与 $P_0$ 和 $I-P_1$ 都正交。

??? proof "引理 3 的证明（点击展开）"

    对任意 $x$，$P_0x\in\operatorname{col}(P_0)\subseteq\operatorname{col}(P_1)$，所以 $P_1P_0x=P_0x$。两边取转置可得 $P_0P_1=P_0$。于是

    \[
    (P_1-P_0)^2=P_1-P_1P_0-P_0P_1+P_0=P_1-P_0.
    \]

    又因为 $P_1-P_0$ 对称，所以它是一个投影。其像空间包含于 $\operatorname{col}(P_1)$，并且它消去 $\operatorname{col}(P_0)$。同时，

    \[
    \operatorname{rank}(P_1-P_0)
    =\operatorname{tr}(P_1-P_0)
    =\operatorname{rank}(P_1)-\operatorname{rank}(P_0),
    \]

    因而其像空间正是 $\operatorname{col}(P_1)\cap\operatorname{col}(P_0)^\perp$。最后，

    \[
    P_0(P_1-P_0)=0,
    \qquad
    (I-P_1)(P_1-P_0)=0,
    \]

    从而得到所需的正交性。$\square$

---

## 4. 中心与非中心卡方分布

!!! info "定义 4.1（非中心卡方分布）"

    设 $Z_1,\ldots,Z_r$ 相互独立，且 $Z_j\sim N(a_j,1)$。定义

    \[
    Q=\sum_{j=1}^rZ_j^2\sim\chi_r^2(\lambda),
    \qquad
    \lambda=\sum_{j=1}^ra_j^2.
    \]

    数值 $\lambda\geq0$ 称为**非中心参数（noncentrality parameter）**。当 $\lambda=0$ 时，它退化为中心卡方分布 $\chi_r^2$。

!!! note "记号约定"

    一些讲义把非中心参数写成 $\delta^2$，并将同一分布记作 $\chi_r^2(\delta^2)$。本章采用常见约定 $\lambda=\delta^2$，标准术语为“非中心”。

!!! success "定理 4（高斯二次型）"

    设 $X\sim N_n(\mu,\sigma^2I_n)$，其中 $\sigma^2>0$，并设 $P$ 是秩为 $r$ 的正交投影，则

    \[
    \frac{X^TPX}{\sigma^2}\sim\chi_r^2(\lambda),
    \qquad
    \lambda=\frac{\mu^TP\mu}{\sigma^2}
    =\frac{\lVert P\mu\rVert^2}{\sigma^2}.
    \]

    特别地，该分布为中心卡方分布，当且仅当 $P\mu=0$。

??? proof "定理 4 的证明（点击展开）"

    由谱定理，存在正交矩阵 $O=(O_1,O_0)$，使得

    \[
    O^TPO=
    \begin{pmatrix}
    I_r&0\\
    0&0
    \end{pmatrix},
    \]

    其中 $O_1$ 的列构成 $\operatorname{col}(P)$ 的一组标准正交基。令 $Z=O^TX/\sigma$，则

    \[
    Z\sim N_n(O^T\mu/\sigma,I_n),
    \]

    且 $Z$ 的坐标相互独立。因此

    \[
    \frac{X^TPX}{\sigma^2}
    =Z^T
    \begin{pmatrix}
    I_r&0\\
    0&0
    \end{pmatrix}Z
    =\sum_{j=1}^rZ_j^2.
    \]

    根据定义，它服从非中心卡方分布，其非中心参数为

    \[
    \lambda
    =\sum_{j=1}^r\left(\frac{o_j^T\mu}{\sigma}\right)^2
    =\frac{\mu^TO_1O_1^T\mu}{\sigma^2}
    =\frac{\mu^TP\mu}{\sigma^2}.
    \]

    因为 $P=P^T=P^2$，所以 $\mu^TP\mu=\lVert P\mu\rVert^2$。该值为零当且仅当 $P\mu=0$。$\square$

!!! success "推论 5（均值与方差）"

    若 $Q\sim\chi_r^2(\lambda)$，则

    \[
    EQ=r+\lambda,
    \qquad
    \operatorname{Var}(Q)=2(r+2\lambda).
    \]

??? proof "推论 5 的证明（点击展开）"

    写成 $Q=\sum_{j=1}^rZ_j^2$，其中相互独立的 $Z_j\sim N(a_j,1)$ 且 $\sum_ja_j^2=\lambda$。因为

    \[
    EZ_j^2=1+a_j^2,
    \qquad
    EZ_j^4=3+6a_j^2+a_j^4,
    \]

    所以

    \[
    \operatorname{Var}(Z_j^2)
    =EZ_j^4-(EZ_j^2)^2
    =2+4a_j^2.
    \]

    利用独立性求和，即得结论。$\square$

---

## 5. 正交高斯分量与 Cochran 定理

!!! success "定理 6（正交投影的独立性）"

    设 $X\sim N_n(\mu,\sigma^2I_n)$。若 $P$ 和 $Q$ 是满足 $PQ=0$ 的正交投影，则 $PX$ 与 $QX$ 相互独立。因此，$X^TPX$ 与 $X^TQX$ 也相互独立。

??? proof "定理 6 的证明（点击展开）"

    因为 $P$ 和 $Q$ 都对称，所以

    \[
    QP=(PQ)^T=0.
    \]

    向量对 $(PX,QX)$ 联合高斯，并且

    \[
    \operatorname{Cov}(PX,QX)
    =P(\sigma^2I_n)Q^T
    =\sigma^2PQ
    =0.
    \]

    由定理 1，$PX\perp QX$。又因为

    \[
    X^TPX=\lVert PX\rVert^2,
    \qquad
    X^TQX=\lVert QX\rVert^2,
    \]

    两个二次型分别是两个独立向量的函数，所以也相互独立。$\square$

!!! success "定理 7（Cochran 定理：投影形式）"

    设 $P_1,\ldots,P_k$ 是两两正交的投影，即当 $i\neq j$ 时 $P_iP_j=0$。若 $X\sim N_n(\mu,\sigma^2I_n)$，则

    \[
    Q_j=\frac{X^TP_jX}{\sigma^2},
    \qquad j=1,\ldots,k,
    \]

    相互独立，并且

    \[
    Q_j\sim\chi_{r_j}^2(\lambda_j),
    \qquad
    r_j=\operatorname{rank}(P_j),
    \qquad
    \lambda_j=\frac{\mu^TP_j\mu}{\sigma^2}.
    \]

    若进一步有 $\sum_jP_j=I_n$，则

    \[
    \frac{\lVert X\rVert^2}{\sigma^2}=\sum_jQ_j,
    \qquad
    \sum_jr_j=n,
    \qquad
    \sum_j\lambda_j=\frac{\lVert\mu\rVert^2}{\sigma^2}.
    \]

??? proof "定理 7 的证明（点击展开）"

    定理 4 给出每个 $Q_j$ 的边缘分布。将投影向量堆叠为

    \[
    W=(P_1X,\ldots,P_kX).
    \]

    它是联合高斯向量，并且当 $i\neq j$ 时，其第 $(i,j)$ 个互协方差块为

    \[
    \sigma^2P_iP_j=0.
    \]

    因而 $W$ 的协方差矩阵是分块对角的。其特征函数可分解为各边缘特征函数的乘积，所以投影向量相互独立，它们的平方范数也相互独立。

    若 $\sum_jP_j=I_n$，则

    \[
    \lVert X\rVert^2
    =X^T\left(\sum_jP_j\right)X
    =\sum_jX^TP_jX.
    \]

    取迹可得 $\sum_jr_j=n$，并且

    \[
    \sum_j\lambda_j
    =\frac{\mu^T(\sum_jP_j)\mu}{\sigma^2}
    =\frac{\lVert\mu\rVert^2}{\sigma^2}.
    \]

    证明完毕。$\square$

!!! success "定理 8（Cochran 定理：秩形式）"

    设 $A_1,\ldots,A_k$ 是对称半正定矩阵，满足

    \[
    \sum_jA_j=I_n.
    \]

    记 $r_j=\operatorname{rank}(A_j)$。若 $\sum_jr_j=n$，则每个 $A_j$ 都是正交投影，并且这些投影两两正交。因此，对于 $X\sim N_n(\mu,\sigma^2I_n)$，

    \[
    \frac{X^TA_jX}{\sigma^2}
    \sim
    \chi_{r_j}^2\left(\frac{\mu^TA_j\mu}{\sigma^2}\right),
    \]

    且这些二次型相互独立。

??? proof "定理 8 的证明（点击展开）"

    因为 $A_j\succeq0$ 且

    \[
    I_n-A_j=\sum_{\ell\neq j}A_\ell\succeq0,
    \]

    所以 $A_j$ 的每个特征值都位于 $[0,1]$。于是

    \[
    \operatorname{tr}(A_j)\leq\operatorname{rank}(A_j)=r_j.
    \]

    但是

    \[
    n=\operatorname{tr}(I_n)
    =\sum_j\operatorname{tr}(A_j)
    \leq\sum_jr_j
    =n.
    \]

    因此所有不等式都必须取等号。对每个 $j$，$A_j$ 的 $r_j$ 个非零特征值位于 $(0,1]$ 且和为 $r_j$，故它们全都等于 $1$，从而 $A_j^2=A_j$。

    固定 $i\neq j$，取 $x\in\operatorname{col}(A_j)$。因为 $A_j$ 已经是投影，$A_jx=x$，所以

    \[
    0=x^T(I_n-A_j)x
    =\sum_{\ell\neq j}x^TA_\ell x.
    \]

    每一项均非负，因此 $x^TA_ix=0$。又因为 $A_i$ 是投影，$x^TA_ix=\lVert A_ix\rVert^2$，所以 $A_ix=0$。这对所有 $x\in\operatorname{col}(A_j)$ 都成立，因此 $A_iA_j=0$。最后由定理 7 得到分布与独立性结论。$\square$

!!! note "为什么秩条件很重要"

    秩条件迫使这些二次型恰好使用 $n$ 个相互垂直的高斯坐标：既不重复使用任何方向，也不遗漏任何方向。

---

## 6. 正态样本：卡方分布与 Student 分布

令 $X=(X_1,\ldots,X_n)^T$，其中 $X_i\overset{\mathrm{iid}}{\sim}N(\mu,\sigma^2)$，并定义

\[
P_1=\frac{1}{n}\mathbf{1}\mathbf{1}^T,
\qquad
M=I_n-P_1.
\]

$P_1$ 和 $M$ 分别投影到 $\operatorname{span}(\mathbf{1})$ 及其正交补，秩分别为 $1$ 和 $n-1$。

!!! success "定理 9（正态样本的均值—方差分解）"

    对于正态样本，

    \[
    \overline{X}\sim N\left(\mu,\frac{\sigma^2}{n}\right),
    \qquad
    \overline{X}\perp S^2,
    \qquad
    \frac{(n-1)S^2}{\sigma^2}\sim\chi_{n-1}^2.
    \]

??? proof "定理 9 的证明（点击展开）"

    因为 $\overline{X}=n^{-1}\mathbf{1}^TX$ 是线性高斯统计量，所以

    \[
    \overline{X}
    \sim
    N\left(n^{-1}\mathbf{1}^T(\mu\mathbf{1}),
    n^{-2}\mathbf{1}^T(\sigma^2I_n)\mathbf{1}\right)
    =N\left(\mu,\frac{\sigma^2}{n}\right).
    \]

    此外，

    \[
    P_1X=\overline{X}\mathbf{1},
    \qquad
    MX=X-\overline{X}\mathbf{1},
    \qquad
    X^TMX=\sum_{i=1}^n(X_i-\overline{X})^2.
    \]

    因为 $P_1M=0$，定理 6 表明 $P_1X$ 与 $MX$ 相互独立，所以 $\overline{X}$ 与 $S^2$ 相互独立。最后，$M$ 的秩为 $n-1$，且 $M(\mu\mathbf{1})=0$，由定理 4，

    \[
    \frac{X^TMX}{\sigma^2}
    =\frac{(n-1)S^2}{\sigma^2}
    \sim\chi_{n-1}^2.
    \]

    证明完毕。$\square$

!!! info "定义 6.1（中心与非中心 Student 分布）"

    设 $Z\sim N(\delta,1)$、$V\sim\chi_\nu^2$，且二者相互独立。则

    \[
    T=\frac{Z}{\sqrt{V/\nu}}\sim t_\nu(\delta).
    \]

    当 $\delta=0$ 时，它就是中心 Student 分布 $t_\nu$。

!!! success "推论 10（原假设与备择假设下的单样本统计量）"

    对任意参考值 $\mu_0$，

    \[
    T=\frac{\sqrt{n}(\overline{X}-\mu_0)}{S}
    \sim t_{n-1}(\delta),
    \qquad
    \delta=\frac{\sqrt{n}(\mu-\mu_0)}{\sigma}.
    \]

    在 $H_0:\mu=\mu_0$ 下，$T\sim t_{n-1}$；在固定备择假设下，$T$ 服从非中心 $t$ 分布。

??? proof "推论 10 的证明（点击展开）"

    令

    \[
    Z=\frac{\sqrt{n}(\overline{X}-\mu_0)}{\sigma},
    \qquad
    V=\frac{(n-1)S^2}{\sigma^2}.
    \]

    定理 9 给出

    \[
    Z\sim N\left(\frac{\sqrt{n}(\mu-\mu_0)}{\sigma},1\right),
    \qquad
    V\sim\chi_{n-1}^2,
    \qquad
    Z\perp V.
    \]

    因为

    \[
    \frac{Z}{\sqrt{V/(n-1)}}
    =\frac{\sqrt{n}(\overline{X}-\mu_0)}{S},
    \]

    所以结论由定义直接得到。$\square$

---

## 7. 中心与非中心 $F$ 分布

!!! info "定义 7.1（中心与单重非中心 $F$ 分布）"

    设 $U\sim\chi_r^2(\lambda)$、$V\sim\chi_s^2$，且二者相互独立。则

    \[
    F=\frac{U/r}{V/s}\sim F_{r,s}(\lambda).
    \]

    当 $\lambda=0$ 时，它就是中心 $F$ 分布 $F_{r,s}$。在标准的单重非中心定义中，分母是中心卡方变量；若分子和分母均为非中心卡方变量，则比值服从双重非中心 $F$ 分布，本章不需要使用它。

### 7.1 两个独立样本的方差

设

\[
X_1,\ldots,X_m\overset{\mathrm{iid}}{\sim}N(\mu_X,\sigma_X^2),
\qquad
Y_1,\ldots,Y_n\overset{\mathrm{iid}}{\sim}N(\mu_Y,\sigma_Y^2),
\]

且两个样本相互独立。

!!! success "命题 11（方差比）"

    有

    \[
    \frac{S_X^2/\sigma_X^2}{S_Y^2/\sigma_Y^2}
    \sim F_{m-1,n-1},
    \qquad
    \frac{S_X^2}{S_Y^2}
    \sim\frac{\sigma_X^2}{\sigma_Y^2}F_{m-1,n-1}.
    \]

    在 $H_0:\sigma_X^2=\sigma_Y^2$ 下，未经标准化的方差比服从中心 $F_{m-1,n-1}$ 分布。

??? proof "命题 11 的证明（点击展开）"

    由定理 9，以下两个变量相互独立：

    \[
    U=\frac{(m-1)S_X^2}{\sigma_X^2}\sim\chi_{m-1}^2,
    \qquad
    V=\frac{(n-1)S_Y^2}{\sigma_Y^2}\sim\chi_{n-1}^2.
    \]

    因此

    \[
    \frac{U/(m-1)}{V/(n-1)}
    =\frac{S_X^2/\sigma_X^2}{S_Y^2/\sigma_Y^2}
    \sim F_{m-1,n-1}.
    \]

    证明完毕。$\square$

!!! warning "一个重要区别"

    当 $\sigma_X^2\neq\sigma_Y^2$ 时，未经标准化的方差比服从**缩放后的中心 $F$ 分布**，而不是非中心 $F$ 分布。非中心性来自平方投影中非零的高斯均值；方差不同只会改变尺度。

### 7.2 嵌套正态线性模型

这是引入非中心 $F$ 分布的主要原因。设 $Y\sim N_n(\mu,\sigma^2I_n)$，并考虑嵌套设计空间

\[
\mathcal{S}_R\subseteq\mathcal{S}_F,
\qquad
\dim(\mathcal{S}_R)=p_R,
\qquad
\dim(\mathcal{S}_F)=p_F,
\]

其投影矩阵分别为 $H_R,H_F$，并令 $q=p_F-p_R$。

!!! success "定理 12（嵌套模型的 $F$ 统计量）"

    假设 $\mu\in\mathcal{S}_F$ 且 $n>p_F$。则

    \[
    F=
    \frac{Y^T(H_F-H_R)Y/q}
    {Y^T(I_n-H_F)Y/(n-p_F)}
    \sim F_{q,n-p_F}(\lambda),
    \]

    其中

    \[
    \lambda
    =\frac{\mu^T(H_F-H_R)\mu}{\sigma^2}
    =\frac{\lVert(H_F-H_R)\mu\rVert^2}{\sigma^2}.
    \]

    对于检验 $H_0:\mu\in\mathcal{S}_R$ 对 $H_1:\mu\in\mathcal{S}_F\setminus\mathcal{S}_R$，原假设下统计量服从中心 $F_{q,n-p_F}$ 分布；固定备择假设下则有 $\lambda>0$。

??? proof "定理 12 的证明（点击展开）"

    由引理 3，

    \[
    P_N=H_F-H_R,
    \qquad
    P_D=I_n-H_F
    \]

    是秩分别为 $q$ 和 $n-p_F$ 的正交投影，并且 $P_NP_D=0$。由定理 4 和定理 6，以下两个变量相互独立：

    \[
    U=\frac{Y^TP_NY}{\sigma^2}\sim\chi_q^2(\lambda),
    \qquad
    V=\frac{Y^TP_DY}{\sigma^2}\sim\chi_{n-p_F}^2(\lambda_D).
    \]

    因为 $\mu\in\mathcal{S}_F$，所以 $P_D\mu=0$，从而 $\lambda_D=0$。因此

    \[
    \frac{U/q}{V/(n-p_F)}\sim F_{q,n-p_F}(\lambda).
    \]

    若 $\mu\in\mathcal{S}_R$，则 $H_R\mu=H_F\mu=\mu$，所以 $\lambda=0$。反之，在 $\mathcal{S}_F$ 内有正交分解

    \[
    \mu=H_R\mu+(H_F-H_R)\mu.
    \]

    若 $\mu\notin\mathcal{S}_R$，则第二个分量非零，所以 $\lambda>0$。$\square$

!!! example "例 2（检验一组回归系数）"

    在完整模型

    \[
    Y=X_1\beta_1+X_2\beta_2+\varepsilon,
    \qquad
    \varepsilon\sim N_n(0,\sigma^2I_n)
    \]

    中，检验 $H_0:\beta_2=0$。假设 $X_2$ 在 $\operatorname{col}(X_1)$ 之外增加了 $q$ 个线性无关方向，取

    \[
    \mathcal{S}_R=\operatorname{col}(X_1),
    \qquad
    \mathcal{S}_F=\operatorname{col}(X_1,X_2).
    \]

    定理 12 的分子是加入 $X_2$ 后增加的拟合平方和，分母是完整模型的残差均方。在备择假设下，$X_2\beta_2$ 中不能被 $X_1$ 解释的分量使 $\lambda>0$。因此

    \[
    P_{\beta_2}\{F>F_{q,n-p_F;1-\alpha}\}
    \]

    是精确功效；反过来求解该概率关于 $n$ 的方程，就可以进行样本量计算。

**几何总结：** 数据可以正交分解为

\[
Y=H_RY+(H_F-H_R)Y+(I_n-H_F)Y,
\]

三项依次对应**约简模型拟合、额外拟合、完整模型残差**。Cochran 定理将这些相互垂直的高斯分量的平方长度转化为独立的卡方变量，而 $F$ 统计量比较额外拟合均方与残差均方。

---

## 附录 A：矩母函数与 Poisson 混合表示

!!! info "引理 13"

    若 $Z\sim N(a,1)$，则对 $t<1/2$，

    \[
    Ee^{tZ^2}
    =(1-2t)^{-1/2}
    \exp\left\{\frac{a^2t}{1-2t}\right\}.
    \]

??? proof "引理 13 的证明（点击展开）"

    通过配方，

    \[
    \begin{aligned}
    Ee^{tZ^2}
    &=\frac{e^{-a^2/2}}{\sqrt{2\pi}}
    \int_{\mathbb{R}}
    \exp\left\{-\frac{1-2t}{2}z^2+az\right\}\,dz\\
    &=\frac{e^{-a^2/2}}{\sqrt{2\pi}}
    \exp\left\{\frac{a^2}{2(1-2t)}\right\}
    \int_{\mathbb{R}}
    \exp\left\{-\frac{1-2t}{2}
    \left(z-\frac{a}{1-2t}\right)^2\right\}\,dz\\
    &=(1-2t)^{-1/2}
    \exp\left\{\frac{a^2t}{1-2t}\right\}.
    \end{aligned}
    \]

    证明完毕。$\square$

!!! success "命题 14（矩母函数、可加性与 Poisson 混合）"

    若 $Q\sim\chi_r^2(\lambda)$，则

    \[
    M_Q(t)
    =(1-2t)^{-r/2}
    \exp\left\{\frac{\lambda t}{1-2t}\right\}.
    \]

    因此，相互独立的非中心卡方变量相加时，其自由度和非中心参数分别相加。

    此外，若

    \[
    K\sim\operatorname{Poisson}(\lambda/2),
    \qquad
    Q\mid K\sim\chi_{r+2K}^2,
    \]

    则 $Q$ 的边缘分布为 $\chi_r^2(\lambda)$。

??? proof "命题 14 的证明（点击展开）"

    对各个独立坐标应用引理 13 并将矩母函数相乘，就得到上面的矩母函数公式。矩母函数的乘积形式也直接证明了可加性。

    对于 Poisson 混合表示，

    \[
    \begin{aligned}
    Ee^{tQ}
    &=(1-2t)^{-r/2}E(1-2t)^{-K}\\
    &=(1-2t)^{-r/2}
    \exp\left\{\frac{\lambda}{2}\left[(1-2t)^{-1}-1\right]\right\}\\
    &=(1-2t)^{-r/2}
    \exp\left\{\frac{\lambda t}{1-2t}\right\}.
    \end{aligned}
    \]

    由矩母函数的唯一性，$Q\sim\chi_r^2(\lambda)$。$\square$

---

## 附录 B：高斯二次型的逆命题

!!! success "定理 15（投影的必要性）"

    设 $Z\sim N_n(0,I_n)$，并设 $A$ 为实对称矩阵。若

    \[
    Z^TAZ\sim\chi_r^2,
    \]

    则 $A$ 必须是秩为 $r$ 的正交投影。

??? proof "定理 15 的证明（点击展开）"

    将 $A$ 对角化为

    \[
    A=O\operatorname{diag}(a_1,\ldots,a_n)O^T.
    \]

    令 $W=O^TZ$，则 $W\sim N_n(0,I_n)$，并且

    \[
    Z^TAZ=\sum_{j=1}^na_jW_j^2.
    \]

    不可能存在负的 $a_j$。例如，若 $a_1<0$，则事件“$W_1^2$ 足够大而其他平方项保持有界”具有正概率，并会使二次型取负值，这与卡方变量非负矛盾。

    对零点附近的 $t$，矩母函数相等给出

    \[
    \prod_{j=1}^n(1-2a_jt)^{-1/2}
    =(1-2t)^{-r/2}.
    \]

    两边先取倒数再平方，得到多项式恒等式

    \[
    \prod_{j=1}^n(1-2a_jt)=(1-2t)^r.
    \]

    每个正的 $a_j$ 都会在左侧产生根 $(2a_j)^{-1}$；而右侧只有根 $1/2$，重数为 $r$。因此，恰有 $r$ 个特征值等于 $1$，其余特征值等于 $0$。所以 $A^2=A$，再结合对称性可知 $A$ 是秩为 $r$ 的正交投影。$\square$

!!! note "为什么对称性是自然的"

    每个二次型只依赖于矩阵的对称部分，因为

    \[
    x^TAx=x^T\left(\frac{A+A^T}{2}\right)x.
    \]

---

## 附录 C：中心分布的密度推导

### C.1 卡方分布

若 $Z\sim N(0,1)$ 且 $W=Z^2$，利用两个反函数分支 $z=\pm\sqrt{w}$，可得

\[
f_W(w)
=\frac{\phi(\sqrt{w})+\phi(-\sqrt{w})}{2\sqrt{w}}
=\frac{w^{-1/2}e^{-w/2}}{2^{1/2}\Gamma(1/2)},
\qquad w>0.
\]

因此 $Z^2$ 服从形状参数为 $1/2$、尺度参数为 $2$ 的 Gamma 分布。将尺度参数相同的独立 Gamma 变量相加可得

\[
f_{\chi_\nu^2}(v)
=\frac{v^{\nu/2-1}e^{-v/2}}{2^{\nu/2}\Gamma(\nu/2)},
\qquad v>0.
\]

### C.2 Student 分布

设 $Z\sim N(0,1)$、$V\sim\chi_\nu^2$ 且相互独立。对于

\[
T=\frac{Z}{\sqrt{V/\nu}},
\]

作变量变换 $z=t\sqrt{v/\nu}$，其 Jacobian 为 $\sqrt{v/\nu}$。于是

\[
\begin{aligned}
f_T(t)
&=\int_0^\infty
\phi\left(t\sqrt{v/\nu}\right)
f_{\chi_\nu^2}(v)
\sqrt{\frac{v}{\nu}}\,dv\\
&=\frac{\Gamma((\nu+1)/2)}{\sqrt{\nu\pi}\,\Gamma(\nu/2)}
\left(1+\frac{t^2}{\nu}\right)^{-(\nu+1)/2},
\qquad t\in\mathbb{R}.
\end{aligned}
\]

最后一步由 Gamma 积分得到。

### C.3 $F$ 分布

设 $U\sim\chi_r^2$、$V\sim\chi_s^2$ 且相互独立。对于

\[
F=\frac{U/r}{V/s},
\]

作变量变换 $u=(r/s)fv$ 并保留 $v$，其 Jacobian 为 $(r/s)v$。利用 Gamma 积分可得

\[
f_F(f)
=\frac{\Gamma((r+s)/2)}{\Gamma(r/2)\Gamma(s/2)}
\left(\frac{r}{s}\right)^{r/2}
f^{r/2-1}
\left(1+\frac{r}{s}f\right)^{-(r+s)/2},
\qquad f>0.
\]

---

## 附录 D：证明核查与分布选择清单

1. 向量是**精确高斯**的，还是仅有渐近正态性？
2. 协方差矩阵是否为 $\sigma^2I$？如果不是，能否对向量进行白化？
3. 二次型矩阵是否对称且幂等？其自由度等于矩阵的秩。
4. 计算

\[
\lambda=\frac{\mu^TP\mu}{\sigma^2};
\]

不要在偏离原假设时把它错误地设为零。
5. 对高斯投影，在断言独立之前检查 $P_iP_j=0$。
6. 将 $F$ 统计量写成 $(U/r)/(V/s)$，并确认在通常的非中心 $F$ 分布中，分母是中心卡方变量。
7. 明确该分布用于原假设校准、置信覆盖率、功效还是样本量计算。
