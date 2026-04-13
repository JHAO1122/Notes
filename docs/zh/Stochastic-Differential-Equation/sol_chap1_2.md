---
tags:
  - 随机微分方程
  - 课后习题
---

# 📝 课后习题精解：条件期望与布朗运动性质

!!! abstract "关于本页"
    本页面收录了《随机微分方程》课程第一章（独立性、条件期望）与第二章（布朗运动及其性质）的核心课后习题解答。所有解答均采用折叠框形式，点击即可查看详细的推导过程。

---

## Part I: 独立性与条件期望

### 习题 1

**题目：**
设随机变量 $X$ 的概率密度函数是 $f(x) = ax(1-x), x \in (0,1)$。除此之外，$f$ 都取零。
(1) 求常数 $a$；
(2) 设 $Y = X^3$，求 $Y$ 的概率密度函数。

??? success "解答 (点击展开)"
    
    **(1) 求解常数 $a$**
    
    由概率密度函数的归一性可知，在全空间上的积分为 1：
    
    $$
    \int_{0}^{1} ax(1-x)dx = 1 
    $$

    解此积分得：

    $$
    a \left( \frac{1}{2}x^2 - \frac{1}{3}x^3 \right) \Bigg|_{0}^{1} = 1 \implies a \left( \frac{1}{2} - \frac{1}{3} \right) = 1
    $$
    
    解得：$a = 6$ [cite: 1519]。

    <br>

    **(2) 求解 $Y=X^3$ 的概率密度函数**
    
    首先求 $Y$ 的累积分布函数 $F_Y(y)$。因为 $X \in (0,1)$，故 $Y = X^3 \in (0,1)$。
    对于 $y \in (0,1)$：
    
    $$
    F_Y(y) = \mathbb{P}(Y \leqslant y) = \mathbb{P}(X^3 \leqslant y) = \mathbb{P}(X \leqslant y^{1/3})
    $$
    
    代入 $X$ 的概率密度函数计算积分：

    $$
    F_Y(y) = \int_{0}^{y^{1/3}} 6x(1-x)dx = \left( 3x^2 - 2x^3 \right) \Big|_{0}^{y^{1/3}} = 3y^{2/3} - 2y
    $$
    
    对 $F_Y(y)$ 求导即可得到 $Y$ 的概率密度函数 $f_Y(y)$ [cite: 1519]：
    
    $$
    f_Y(y) = F_Y'(y) = 3 \cdot \frac{2}{3}y^{-1/3} - 2 = 2y^{-1/3} - 2
    $$
    
    综上所述，$Y$ 的概率密度函数为：
    
    $$
    f_Y(y) = \begin{cases} 2y^{-1/3} - 2, & y \in (0,1) \\ 0, & \text{其他} \end{cases}
    $$

---

### 习题 2

**题目：**
设 $X, Y$ 是独立的随机变量，$f(x,y)$ 是有界的连续函数，$F_X$ 是 $X$ 的概率分布函数。证明：
(1) $\mathbb{E}[f(X,Y)|Y] = \int f(x,Y)dF_X(x)$；
(2) $\mathbb{P}(X+Y \leqslant x | Y) = F_X(x-Y)$。

??? success "解答 (点击展开)"

    **(1) 证明：**
    
    [cite_start]我们采用测度论中的标准测度逼近方法（Standard Machine）进行严谨证明 [cite: 1520]：
    
    **第一步（示性函数）：** 设 $f(x,y) = I_A(x)I_B(y)$，其中 $A, B$ 为 Borel 集。
    因为 $X, Y$ 独立，利用条件期望的性质（已知 $Y$ 时 $I_B(Y)$ 可提出，且独立变量的条件期望等于无条件期望）：
    
    $$
    \mathbb{E}[I_A(X)I_B(Y)|Y] = I_B(Y)\mathbb{E}[I_A(X)|Y] = I_B(Y)\mathbb{E}[I_A(X)] = I_B(Y)\mathbb{P}(X \in A)
    $$
    
    另一方面，对右侧的积分形式：
    
    $$
    \int I_A(x)I_B(Y)dF_X(x) = I_B(Y) \int I_A(x)dF_X(x) = I_B(Y)\mathbb{P}(X \in A)
    $$
    
    两者相等，故对示性函数成立。
    
    **第二步（简单函数）：** 由期望的线性性质，该结论对有限个示性函数线性组合的简单函数成立 [cite: 1520]。
    
    **第三步（非负可测函数）：** 对于任意非负可测函数，存在单调递增的简单函数序列逼近它，由单调收敛定理 (Monotone Convergence Theorem)，该结论对非负可测函数成立 [cite: 1520]。
    
    **第四步（有界连续函数）：** 任何有界连续函数均可拆分为正部和负部 $f = f^+ - f^-$，且依定义勒贝格可积，故原式对所有有界连续函数成立。得证 [cite: 1520]。

    <br>

    **(2) 证明：**
    
    条件概率可以直接写为条件期望的形式：
    
    $$
    \mathbb{P}(X+Y \leqslant x | Y) = \mathbb{E}[I_{\{X+Y \leqslant x\}} | Y]
    $$
    
    我们令 $g(X, Y) = I_{\{X+Y \leqslant x\}}$，利用第(1)问已证的结论：
    
    $$
    \mathbb{E}[I_{\{X+Y \leqslant x\}} | Y] = \int I_{\{u+Y \leqslant x\}} dF_X(u) 
    $$

    对积分域进行等价转换 $u \leqslant x - Y$：

    $$
    \int I_{\{u \leqslant x-Y\}} dF_X(u) = \int_{-\infty}^{x-Y} dF_X(u) = F_X(x-Y)
    $$
    
    得证 [cite: 1521]。

---

### 习题 3

**题目：**
设 $A, B$ 是测度空间 $(\Omega, \mathcal{F}, \mathbb{P})$ 中的两个事件，计算条件期望 $\mathbb{E}[\chi_A|\chi_B]$。（注：$\chi$ 为示性函数）

??? success "解答 (点击展开)"

    由于 $\chi_B$ 只能取 0 或 1，由其生成的 $\sigma$-代数非常简单：$\sigma(\chi_B) = \{\emptyset, \Omega, B, B^c\}$。
    因此，条件期望 $\mathbb{E}[\chi_A|\chi_B]$ 必须是 $\sigma(\chi_B)$-可测的，这意味着它在 $B$ 和 $B^c$ 上必须是常数。我们设为：
    
    $$
    \mathbb{E}[\chi_A|\chi_B] = c_1 \chi_B + c_2 \chi_{B^c}
    $$
    
    根据条件期望的 Radon-Nikodym 导数定义，对于任意 $\Lambda \in \sigma(\chi_B)$，必须满足积分相等：
    
    $$
    \int_\Lambda \mathbb{E}[\chi_A|\chi_B] d\mathbb{P} = \int_\Lambda \chi_A d\mathbb{P}
    $$
    
    **情况一：取 $\Lambda = B$**
    
    $$
    \int_B (c_1 \chi_B + c_2 \chi_{B^c}) d\mathbb{P} = \int_B \chi_A d\mathbb{P} \implies c_1 \mathbb{P}(B) = \mathbb{P}(A \cap B)
    $$
    
    若 $\mathbb{P}(B) > 0$，则 $c_1 = \frac{\mathbb{P}(A \cap B)}{\mathbb{P}(B)} = \mathbb{P}(A|B)$ [cite: 1523]。
    
    **情况二：取 $\Lambda = B^c$**
    
    $$
    \int_{B^c} (c_1 \chi_B + c_2 \chi_{B^c}) d\mathbb{P} = \int_{B^c} \chi_A d\mathbb{P} \implies c_2 \mathbb{P}(B^c) = \mathbb{P}(A \cap B^c)
    $$
    
    若 $\mathbb{P}(B^c) > 0$，则 $c_2 = \frac{\mathbb{P}(A \cap B^c)}{\mathbb{P}(B^c)} = \mathbb{P}(A|B^c)$ [cite: 1524]。
    
    综上所述，条件期望的显式表达为：
    
    $$
    \mathbb{E}[\chi_A|\chi_B] = \mathbb{P}(A|B)\chi_B + \mathbb{P}(A|B^c)\chi_{B^c}
    $$

---

### 习题 4

**题目：**
设 $\mathcal{V}_1$ 和 $\mathcal{V}_2$ 是两个独立的 $\sigma$-域，$X$ 是个可积的随机变量。证明：
$\mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)|\mathcal{V}_2] = \mathbb{E}[X]$。

??? success "解答 (点击展开)"

    **证明：**
    
    首先，由条件期望的定义可知，内层期望 $\mathbb{E}(X|\mathcal{V}_1)$ 必然是一个 $\mathcal{V}_1$-可测的随机变量 [cite: 1529]。
    
    由于题目已知 $\mathcal{V}_1$ 和 $\mathcal{V}_2$ 是相互独立的 $\sigma$-域，因此 $\mathcal{V}_1$-可测的随机变量 $\mathbb{E}(X|\mathcal{V}_1)$ 与 $\sigma$-域 $\mathcal{V}_2$ 自然是相互独立的。
    
    根据条件期望关于独立 $\sigma$-域的性质（如果随机变量 $Z$ 独立于 $\sigma$-域 $\mathcal{G}$，则 $\mathbb{E}[Z|\mathcal{G}] = \mathbb{E}[Z]$），我们可以去掉关于 $\mathcal{V}_2$ 的条件 [cite: 1529]：
    
    $$
    \mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)|\mathcal{V}_2] = \mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)]
    $$
    
    最后，根据条件期望的全期望公式（即 Tower Property 的退化形式）：
    
    $$
    \mathbb{E}[\mathbb{E}(X|\mathcal{V}_1)] = \mathbb{E}[X]
    $$
    
    联合两式，原命题得证 [cite: 1529]。 $\square$

---

### 习题 5

**题目：**
设 $X$ 和 $\{X_n\}$ 是 $(\Omega, \mathcal{F}, \mathbb{P})$ 上的一列随机变量，对固定的 $1 \leqslant p < \infty$，$\mathbb{E}[|X_n|^p] < +\infty$。假设 $\lim_{n\to\infty}\mathbb{E}[|X_n - X|^p] = 0$，证明：对 $\mathcal{F}$ 的任意子 $\sigma$-域 $\mathcal{V} \subset \mathcal{F}$，成立

$$
\lim_{n\to\infty}\mathbb{E}\left[\big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p\right] = 0.
$$

??? success "解答 (点击展开)"

    **证明：**
    
    考察函数 $\Phi(x) = |x|^p$。由于 $p \geqslant 1$，这是一个凸函数。
    利用条件期望的线性性质和 **条件 Jensen 不等式**，我们将绝对值内的条件期望合并，并把绝对值放缩到期望内部 [cite: 1530]：
    
    $$
    \big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p = \big|\mathbb{E}[X_n - X | \mathcal{V}]\big|^p \leqslant \mathbb{E}\big[|X_n - X|^p \big| \mathcal{V}\big]
    $$
    
    对上述不等式两边同时取无条件期望，应用全期望公式 [cite: 1530, 1531]：
    
    $$
    \mathbb{E}\left[\big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p\right] \leqslant \mathbb{E}\left[ \mathbb{E}\big[|X_n - X|^p \big| \mathcal{V}\big] \right] = \mathbb{E}\big[|X_n - X|^p\big]
    $$
    
    对不等式两边令 $n \to \infty$ 取极限。由于题目已知 $\lim_{n\to\infty}\mathbb{E}[|X_n - X|^p] = 0$，且左侧期望恒非负，由夹逼定理可得 [cite: 1532]：
    
    $$
    \lim_{n\to\infty}\mathbb{E}\left[\big|\mathbb{E}[X_n|\mathcal{V}] - \mathbb{E}[X|\mathcal{V}]\big|^p\right] = 0
    $$
    
    原命题得证。 $\square$

---

### 习题 6

**题目：**
取 $\Omega = \{1, 2, \cdots, 7, 8\}$，$\mathcal{F} = 2^\Omega$（即 $\Omega$ 的所有子集组成的 $\sigma$-域）。当 $i \leqslant 4$ 时 $\mathbb{P}(\{i\}) = 1/10$，当 $i > 4$ 时 $\mathbb{P}(\{i\}) = 3/20$。
定义 $X = \chi_{\{1,2,3,4\}} + 2\chi_{\{5,6,7,8\}}$，$Y = \chi_{\{1,5\}} + 2\chi_{\{2,3,4,6,7,8\}}$。
$\mathcal{V}$ 是由 $\{1,2\}$ 和 $\{3,4\}$ 生成的 $\sigma$-域，$\mathcal{H}$ 是 $\{1,2,3,4\}$ 生成的 $\sigma$-域。
计算：(1) $X\mathbb{E}[Y]$；(2) $\mathbb{E}[\mathbb{E}[XY|\mathcal{V}]|\mathcal{H}]$。

??? success "解答 (点击展开)"

    **计算准备：明确概率结构**
    
    * $\mathbb{P}(\{1,2,3,4\}) = 4 \times \frac{1}{10} = \frac{2}{5}$
    * $\mathbb{P}(\{5,6,7,8\}) = 4 \times \frac{3}{20} = \frac{3}{5}$
    * 容易验证，$\mathcal{H} = \sigma(\{1,2,3,4\})$ 的两个核心原子就是 $\{1,2,3,4\}$ 和 $\{5,6,7,8\}$。而 $\mathcal{V} = \sigma(\{1,2\}, \{3,4\})$ 的划分更细，显然有 $\mathcal{H} \subset \mathcal{V}$。

    <br>

    **(1) 计算 $X\mathbb{E}[Y]$**
    
    首先计算 $Y$ 的无条件期望 [cite: 1533]：
    
    $$
    \mathbb{E}[Y] = 1 \cdot \mathbb{P}(\{1,5\}) + 2 \cdot \mathbb{P}(\{2,3,4,6,7,8\})
    $$
    
    $$
    = 1 \cdot \left(\frac{1}{10} + \frac{3}{20}\right) + 2 \cdot \left(3 \times \frac{1}{10} + 3 \times \frac{3}{20}\right) = \frac{1}{4} + 2\left(\frac{3}{4}\right) = \frac{7}{4}
    $$
    
    直接将其乘上 $X$：
    
    $$
    X\mathbb{E}[Y] = \frac{7}{4}X = \frac{7}{4}\chi_{\{1,2,3,4\}} + \frac{14}{4}\chi_{\{5,6,7,8\}}
    $$

    <br>

    **(2) 计算 $\mathbb{E}[\mathbb{E}[XY|\mathcal{V}]|\mathcal{H}]$**
    
    这里我们利用极其优雅的平滑性质 (Tower Property) 直接降维，避免双重条件期望的繁琐计算 [cite: 1533]。
    因为 $\mathcal{H} \subset \mathcal{V}$，较小的 $\sigma$-代数在外部起决定性作用：
    
    $$
    \mathbb{E}\big[\mathbb{E}[XY|\mathcal{V}]\big|\mathcal{H}\big] = \mathbb{E}[XY|\mathcal{H}]
    $$
    
    观察随机变量 $X$，它在 $\{1,2,3,4\}$ 上恒为 1，在 $\{5,6,7,8\}$ 上恒为 2，这恰好完全对应 $\mathcal{H}$ 的原子划分。因此 $X$ 是 **$\mathcal{H}$-可测的**！作为已知信息，可以将其提出条件期望 [cite: 1535]：
    
    $$
    \mathbb{E}[XY|\mathcal{H}] = X \mathbb{E}[Y|\mathcal{H}]
    $$
    
    接着我们在 $\mathcal{H}$ 的两个原子上分别计算 $Y$ 的条件期望：
    * 在 $\{1,2,3,4\}$ 上：$\mathbb{E}[Y|\{1,2,3,4\}] = \frac{1\cdot P(\{1\}) + 2\cdot P(\{2,3,4\})}{P(\{1,2,3,4\})} = \frac{1/10 + 6/10}{4/10} = \frac{7}{4}$
    * 在 $\{5,6,7,8\}$ 上：$\mathbb{E}[Y|\{5,6,7,8\}] = \frac{1\cdot P(\{5\}) + 2\cdot P(\{6,7,8\})}{P(\{5,6,7,8\})} = \frac{3/20 + 18/20}{12/20} = \frac{7}{4}$
    
    惊人地发现，在两块划分上 $\mathbb{E}[Y|\mathcal{H}] = 7/4$ 恒成立。故：
    
    $$
    \mathbb{E}[\mathbb{E}[XY|\mathcal{V}]|\mathcal{H}] = \frac{7}{4}X
    $$

---

### 习题 7

**题目：**
给定概率空间 $(\Omega, \mathcal{F}, \mathbb{P})$ 及可积的随机变量 $X$，又设 $\{\mathcal{F}(t)\}_{t \geqslant 0}$ 是一个 $\sigma$-域流。对 $t \geqslant 0$，定义 $X(t) \doteq \mathbb{E}[X|\mathcal{F}(t)]$。证明：$X(t)$ 关于 $\mathcal{F}(t)$ 是一个鞅。

??? success "解答 (点击展开)"

    **证明：**
    [cite_start]我们要严格验证鞅的三大核心条件 [cite: 1546, 1547]：
    
    **1. 可积性 (Integrability)：**
    利用条件期望的收缩性（Jensen 不等式的绝对值特例）：
    
    $$
    \mathbb{E}|X(t)| = \mathbb{E}\big|\mathbb{E}[X|\mathcal{F}(t)]\big| \leqslant \mathbb{E}\big[\mathbb{E}[|X| \big| \mathcal{F}(t)]\big] = \mathbb{E}|X|
    $$
    
    由于题目给定 $X$ 是可积的（$\mathbb{E}|X| < \infty$），故 $\mathbb{E}|X(t)| < \infty$。
    
    **2. 适应性 (Adaptability)：**
    根据条件期望 $\mathbb{E}[X|\mathcal{F}(t)]$ 的测度论定义，它天然是 $\mathcal{F}(t)$-可测的。因此过程 $\{X(t)\}$ 适应于信息流 $\{\mathcal{F}(t)\}_{t \ge 0}$。
    
    **3. 鞅性质 (Martingale Property)：**
    对于任意的 $0 \leqslant s \leqslant t$，由于 $\{\mathcal{F}(t)\}$ 是 $\sigma$-域流，即信息是不断累积的，必有 $\mathcal{F}(s) \subset \mathcal{F}(t)$。
    根据条件期望的平滑性质 (Tower Property)：
    
    $$
    \mathbb{E}[X(t) | \mathcal{F}(s)] = \mathbb{E}\big[\mathbb{E}[X|\mathcal{F}(t)] \big| \mathcal{F}(s)\big] = \mathbb{E}[X|\mathcal{F}(s)] = X(s)
    $$
    
    综上三点，由单一可积变量在不同信息流上的投影所构成的过程 $X(t)$，必是一个鞅（这在随机分析中被称为 Doob 鞅）。 $\square$

---
<br>

## Part II: 布朗运动及其性质

### 习题 1

**题目：**
设 $W(t)$ 是一维 Brown 运动，证明对任意固定的 $s>0$，$W(t+s)-W(s)$ 是 Brown 运动；对任意正数 $c$，$cW(t/c^2)$ 也是 Brown 运动。

??? success "解答 (点击展开)"

    **(1) 证明 $B(t) = W(t+s) - W(s)$ 是布朗运动**
    
    我们逐一验证布朗运动的四条定义 [cite: 1551]：
    
    1. **零初值**：$B(0) = W(0+s) - W(s) = 0$。
    2. **独立增量**：对于任意 $0 \le t_1 < t_2 < \cdots < t_k$，其增量序列为：
       $$B(t_j) - B(t_{j-1}) = W(t_j+s) - W(t_{j-1}+s)$$
       由于原过程 $W$ 具有独立增量性质，且时间区间 $[t_{j-1}+s, t_j+s]$ 在时间轴上互不重叠，因此这些增量是相互独立的。
    3. **平稳正态增量**：对于 $t > u$，增量分布为：
       $$B(t) - B(u) = W(t+s) - W(u+s) \sim N(0, (t+s) - (u+s)) = N(0, t-u)$$
    4. **路径连续性**：因为 $W(\cdot)$ 的样本路径几乎必然连续，所以平移后的路径 $B(\cdot)$ 也几乎必然连续。
    
    综上，$B(t)$ 是一维布朗运动。

    <br>

    **(2) 证明 $U(t) = cW(t/c^2)$ 是布朗运动**
    
    同理进行验证 [cite: 1552, 1553]：
    
    1. **零初值**：$U(0) = cW(0/c^2) = cW(0) = 0$。
    2. **独立增量**：对于任意 $0 \le t_1 < t_2 < \cdots < t_k$，由于 $t_1/c^2 < t_2/c^2 < \cdots < t_k/c^2$ 不重叠，故原布朗运动的增量相互独立，乘上常数 $c$ 后依然独立。
    3. **平稳正态增量**：由于线性变换的性质，均值仍为 $0$，方差为：
       $$Var(U(t) - U(s)) = Var\left( c(W(t/c^2) - W(s/c^2)) \right) = c^2 \cdot \left( \frac{t-s}{c^2} \right) = t-s$$
       故 $U(t) - U(s) \sim N(0, t-s)$。
    4. **路径连续性**：尺度变换不改变样本路径的连续性。
    
    综上，$U(t) = cW(t/c^2)$ 也是布朗运动（这被称为布朗运动的标度不变性/分形性质）。

---

### 习题 2

**题目：**
设 $W(t)$ 是一维 Brown 运动，记 $\tilde{W}(t) = \begin{cases} tW(1/t), & t > 0, \\ 0, & t = 0. \end{cases}$
证明：$\tilde{W}(t) - \tilde{W}(s) \sim N(0, t-s), \forall 0 < s < t$。

??? success "解答 (点击展开)"

    对于 $0 < s < t$，我们将增量进行恒等变形，凑出相互独立的项 [cite: 1557, 1558]：

    $$
    \tilde{W}(t) - \tilde{W}(s) = tW\left(\frac{1}{t}\right) - sW\left(\frac{1}{s}\right)
    $$

    我们通过加一项减一项 $sW(1/t)$ 来分离：

    $$
    = (t-s)W\left(\frac{1}{t}\right) - s\left(W\left(\frac{1}{s}\right) - W\left(\frac{1}{t}\right)\right)
    $$

    注意，这里 $1/t < 1/s$。利用布朗运动的独立增量性质，随机变量 $W(1/t) = W(1/t) - W(0)$ 与增量 $W(1/s) - W(1/t)$ 是相互独立的。
    
    由于相互独立，方差可以直接相加 [cite: 1561]：

    $$
    Var(\tilde{W}(t) - \tilde{W}(s)) = (t-s)^2 Var\left( W\left(\frac{1}{t}\right) \right) + s^2 Var\left( W\left(\frac{1}{s}\right) - W\left(\frac{1}{t}\right) \right)
    $$

    代入方差公式 $Var(W(u)) = u$：

    $$
    = (t-s)^2 \left(\frac{1}{t}\right) + s^2 \left(\frac{1}{s} - \frac{1}{t}\right)
    $$

    通分化简：

    $$
    = \frac{t^2 - 2ts + s^2}{t} + s - \frac{s^2}{t} = \frac{t^2 - 2ts}{t} + s = t - 2s + s = t - s
    $$

    又因为期望的线性性质，$E[\tilde{W}(t) - \tilde{W}(s)] = 0 - 0 = 0$ [cite: 1562]。
    由于它是联合正态随机变量的线性组合，因此必然服从正态分布。
    
    结论得证：$\tilde{W}(t) - \tilde{W}(s) \sim N(0, t-s)$。

---

### 习题 3

**题目：**
设 $W(t)$ 是 Brown 运动，证明：$\mathbb{E}[W^{2k}(t)] = \frac{(2k)!t^k}{2^k k!}, \forall t > 0$。

??? success "解答 (点击展开)"

    由于 $W(t) \sim N(0, t)$，利用正态分布的矩母函数 (Moment Generating Function) $M_W(\lambda) = \mathbb{E}[e^{\lambda W(t)}]$ [cite: 1564]：

    $$
    M_W(\lambda) = \exp\left( \frac{1}{2} \lambda^2 t \right)
    $$

    我们将等式两边同时在 $\lambda = 0$ 处展开为泰勒级数 (Taylor Series)：
    
    左边根据期望的线性展开：
    
    $$
    \mathbb{E}[e^{\lambda W(t)}] = \sum_{m=0}^{\infty} \mathbb{E}[W^m(t)] \frac{\lambda^m}{m!}
    $$
    
    右边将指数函数展开：
    
    $$
    \exp\left( \frac{1}{2} \lambda^2 t \right) = \sum_{k=0}^{\infty} \frac{1}{k!} \left( \frac{1}{2} \lambda^2 t \right)^k = \sum_{k=0}^{\infty} \frac{t^k}{2^k k!} \lambda^{2k}
    $$
    
    对比两边 $\lambda^{2k}$ 次方的系数：
    
    $$
    \frac{\mathbb{E}[W^{2k}(t)]}{(2k)!} = \frac{t^k}{2^k k!}
    $$
    
    移项即可得证 [cite: 1571]：
    
    $$
    \mathbb{E}[W^{2k}(t)] = \frac{(2k)! t^k}{2^k k!}
    $$
    
    *(注：对于奇数阶矩，对比两边系数由于右边没有奇数次项，故 $\mathbb{E}[W^{2k+1}(t)] = 0$)*。

---

### 习题 4

**题目：**
证明：设 $c$ 是常数，$0 < s < t$，那么 $\mathbb{E}[\exp(c(W(s) - W(t)))] = \exp\left(\frac{1}{2}c^2(t-s)\right)$。

??? success "解答 (点击展开)"

    由于布朗运动的增量具有平稳性与对称性，$W(s) - W(t)$ 与 $W(t) - W(s)$ 同分布，均服从 $N(0, t-s)$ 分布。
    
    这等价于求正态随机变量的矩母函数 [cite: 1573]。
    令 $Z = c(W(s) - W(t))$，则 $Z \sim N(0, c^2(t-s))$。
    
    根据正态分布 $N(\mu, \sigma^2)$ 的矩母函数公式 $\mathbb{E}[e^Z] = \exp(\mu + \frac{1}{2}\sigma^2)$ [cite: 1573, 1574]：
    
    $$
    \mathbb{E}[\exp(c(W(s) - W(t)))] = \exp\left( 0 + \frac{1}{2} \cdot c^2(t-s) \right) = \exp\left( \frac{1}{2}c^2(t-s) \right)
    $$
    
    得证。

---

### 习题 5

**题目：**
设 $U(t) = e^{-t}W(e^{2t})$，则 $\mathbb{E}[U(t)U(s)] = e^{-|t-s|}, \forall t, s \in \mathbb{R}$。

??? success "解答 (点击展开)"

    不失一般性，我们假设 $t \geqslant s$。
    
    代入 $U(t)$ 的定义计算交叉项的期望：

    $$
    \mathbb{E}[U(t)U(s)] = \mathbb{E}\left[ e^{-t}W(e^{2t}) \cdot e^{-s}W(e^{2s}) \right] = e^{-(t+s)} \mathbb{E}\left[ W(e^{2t})W(e^{2s}) \right]
    $$

    对于布朗运动，我们已知其协方差函数为 $\mathbb{E}[W(u)W(v)] = \min(u, v)$。
    由于我们假定了 $t \geqslant s$，自然有 $e^{2t} \geqslant e^{2s}$，因此：

    $$
    \mathbb{E}\left[ W(e^{2t})W(e^{2s}) \right] = \min(e^{2t}, e^{2s}) = e^{2s}
    $$

    将其代回原式 [cite: 1576, 1577]：

    $$
    \mathbb{E}[U(t)U(s)] = e^{-(t+s)} \cdot e^{2s} = e^{s-t}
    $$

    由于 $t \geqslant s$，$s-t = -|t-s|$。若 $s > t$ 结论对称同样成立。
    故 $\forall t, s \in \mathbb{R}$，$\mathbb{E}[U(t)U(s)] = e^{-|t-s|}$，得证。

---

### 习题 6

**题目：**
证明：几乎必然成立 $\lim_{m \to \infty} \frac{W(m)}{m} = 0$。

??? success "解答 (点击展开)"

    我们可以利用切比雪夫不等式 (Chebyshev's Inequality) 结合 Borel-Cantelli 引理来进行严格测度论意义上的证明。
    
    **第一步：方差界定**
    考察随机序列 $Y_m = \frac{W(m)}{m}$。由于 $W(m) \sim N(0, m)$，有：

    $$
    Var\left(\frac{W(m)}{m}\right) = \frac{1}{m^2} Var(W(m)) = \frac{m}{m^2} = \frac{1}{m}
    $$

    **第二步：强化矩界限（为 Borel-Cantelli 铺垫）**
    如果仅用二阶矩的切比雪夫不等式，会得到 $P(|Y_m| \ge \epsilon) \le \frac{1}{m\epsilon^2}$ [cite: 1583]。但 $\sum_{m=1}^{\infty} \frac{1}{m}$ 是调和级数，不收敛，无法直接使用 B-C 引理。
    因此我们使用四阶矩（或利用正态分布的指数衰减性质）。对于标准正态变量 $Z \sim N(0,1)$，$\mathbb{E}[Z^4] = 3$。
    
    $$
    \mathbb{E}[W(m)^4] = 3m^2 \implies \mathbb{E}\left[\left(\frac{W(m)}{m}\right)^4\right] = \frac{3m^2}{m^4} = \frac{3}{m^2}
    $$
    
    应用基于四阶矩的马尔可夫不等式：
    
    $$
    \mathbb{P}\left(\left|\frac{W(m)}{m}\right| \geqslant \epsilon\right) = \mathbb{P}\left(\left(\frac{W(m)}{m}\right)^4 \geqslant \epsilon^4\right) \leqslant \frac{3}{m^2 \epsilon^4}
    $$

    **第三步：Borel-Cantelli 引理**
    对于任意固定的 $\epsilon > 0$，由于 $\sum_{m=1}^{\infty} \frac{1}{m^2} < \infty$，级数绝对收敛：
    
    $$
    \sum_{m=1}^{\infty} \mathbb{P}\left(\left|\frac{W(m)}{m}\right| \geqslant \epsilon\right) \leqslant \sum_{m=1}^{\infty} \frac{3}{m^2 \epsilon^4} < \infty
    $$
    
    根据 Borel-Cantelli 引理第一部分，事件 $\left\{ \left|\frac{W(m)}{m}\right| \geqslant \epsilon \right\}$ 发生无穷次的概率为 0。
    这等价于，几乎所有的样本路径满足 $\lim_{m \to \infty} \frac{W(m)}{m} = 0$，得证。

---

### 习题 7

**题目：**
证明 $W(t)^2 - t$ 和 $\exp\left(\lambda W_t - \frac{1}{2}\lambda^2 t\right)$ $(\lambda \in \mathbb{R})$ 关于 $W(t)$ 的历史 $\mathcal{F}(t)$ 都是鞅。

??? success "解答 (点击展开)"

    **证明 (1): $W(t)^2 - t$ 是鞅**
    
    设 $0 \leqslant s < t$。我们需要证明 $\mathbb{E}[W(t)^2 - t | \mathcal{F}(s)] = W(s)^2 - s$。
    利用恒等变形 $W(t) = (W(t) - W(s)) + W(s)$ 展开平方项：

    $$
    \mathbb{E}[W(t)^2 | \mathcal{F}(s)] = \mathbb{E}[((W(t) - W(s)) + W(s))^2 | \mathcal{F}(s)]
    $$
    
    $$
    = \mathbb{E}[(W(t) - W(s))^2 | \mathcal{F}(s)] + 2\mathbb{E}[W(s)(W(t) - W(s)) | \mathcal{F}(s)] + \mathbb{E}[W(s)^2 | \mathcal{F}(s)]
    $$
    
    * [cite_start]第一项：由于布朗运动独立增量，$W(t)-W(s)$ 独立于 $\mathcal{F}(s)$，故条件期望等于无条件期望 $\mathbb{E}[(W(t)-W(s))^2] = t-s$ [cite: 1585, 1586]。
    * 第二项：$W(s)$ 是 $\mathcal{F}(s)$-可测的，可提出。剩余增量期望为 0。
    * 第三项：$W(s)^2$ 是 $\mathcal{F}(s)$-可测的，直接等于自身。
    
    因此：
    
    $$
    \mathbb{E}[W(t)^2 | \mathcal{F}(s)] = (t-s) + 0 + W(s)^2
    $$
    
    减去 $t$：
    
    $$
    \mathbb{E}[W(t)^2 - t | \mathcal{F}(s)] = (t-s) + W(s)^2 - t = W(s)^2 - s
    $$
    
    鞅性质得证。

    <br>

    **证明 (2): $\exp\left(\lambda W(t) - \frac{1}{2}\lambda^2 t\right)$ 是鞅 (几何布朗鞅)**
    
    同样设 $0 \leqslant s < t$。利用指数的拆分：

    $$
    \mathbb{E}[\exp(\lambda W(t)) | \mathcal{F}(s)] = \mathbb{E}[\exp(\lambda(W(t) - W(s))) \cdot \exp(\lambda W(s)) | \mathcal{F}(s)]
    $$
    
    因为 $\exp(\lambda W(s))$ 是 $\mathcal{F}(s)$-可测的，将其作为常数提出：

    $$
    = \exp(\lambda W(s)) \cdot \mathbb{E}[\exp(\lambda(W(t) - W(s))) | \mathcal{F}(s)]
    $$
    
    由于增量 $W(t)-W(s)$ 独立于 $\mathcal{F}(s)$ 且服从 $N(0, t-s)$ 分布，利用正态分布矩母函数 $\mathbb{E}[e^{\lambda Z}] = e^{\lambda^2 \sigma^2 / 2}$：

    $$
    = \exp(\lambda W(s)) \cdot \exp\left( \frac{1}{2} \lambda^2 (t-s) \right)
    $$
    
    将上述期望代入要证明的式子 [cite: 1587]：

    $$
    \mathbb{E}\left[ \exp\left(\lambda W(t) - \frac{1}{2}\lambda^2 t\right) \Big| \mathcal{F}(s) \right] = \exp\left(-\frac{1}{2}\lambda^2 t\right) \cdot \exp(\lambda W(s)) \cdot \exp\left( \frac{1}{2} \lambda^2 (t-s) \right)
    $$
    
    指数合并：
    
    $$
    -\frac{1}{2}\lambda^2 t + \frac{1}{2}\lambda^2 t - \frac{1}{2}\lambda^2 s = -\frac{1}{2}\lambda^2 s
    $$

    最终得到：

    $$
    = \exp\left(\lambda W(s) - \frac{1}{2}\lambda^2 s\right)
    $$

    鞅性质得证。

---

### 习题 8

**题目：**
置 $X(t) = \int_0^t W(s)ds$。证明：$\mathbb{E}[X^2(t)] = \frac{t^3}{3}, \forall t > 0$。

??? success "解答 (点击展开)"

    [cite_start]这是一个极具代表性的随机过程积分计算，核心在于利用 Fubini 定理交换期望和积分顺序 [cite: 1588]。

    首先将平方项写成双重黎曼积分：

    $$
    X^2(t) = \left( \int_0^t W(u)du \right) \left( \int_0^t W(v)dv \right) = \int_0^t \int_0^t W(u)W(v) du dv
    $$

    对等式两边取期望，并利用 Fubini 定理（因积分域有界且函数绝对可积）将期望穿透到积分内部：

    $$
    \mathbb{E}[X^2(t)] = \mathbb{E}\left[ \int_0^t \int_0^t W(u)W(v) du dv \right] = \int_0^t \int_0^t \mathbb{E}[W(u)W(v)] du dv
    $$

    由于布朗运动的协方差函数为 $\mathbb{E}[W(u)W(v)] = \min(u, v)$，代入上式：

    $$
    = \int_0^t \int_0^t \min(u, v) du dv
    $$

    由于 $\min(u,v)$ 在对角线 $u=v$ 处不可微，我们将积分区域 $[0,t] \times [0,t]$ 沿着对角线划分为两块：$u \leqslant v$ 区域和 $u > v$ 区域 [cite: 1589, 1590]：

    $$
    = \int_0^t dv \int_0^v u du + \int_0^t du \int_0^u v dv
    $$

    由于这两块区域积分完全对称，计算其中一块并乘以 2 即可：

    $$
    = 2 \int_0^t \left( \int_0^v u du \right) dv = 2 \int_0^t \frac{1}{2}v^2 dv = \int_0^t v^2 dv = \frac{1}{3}v^3 \Bigg|_0^t = \frac{t^3}{3}
    $$

    结论得证：$\mathbb{E}[X^2(t)] = \frac{t^3}{3}$。