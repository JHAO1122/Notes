# 第二章：特征函数

在渐近统计理论中，我们迫切需要强有力的工具来推导依分布收敛（弱收敛）。特征函数 (Characteristic Function) 正是这样一种“神兵利器”，它提供了一种从**频域 (Frequency Domain)** 视角审视概率分布的方法，并能完美且唯一地刻画一个分布。

---

## 1. 特征函数的定义与基本性质

!!! abstract "定义 2.1：特征函数 (Characteristic Function)"

    对于任意具有分布函数 $F$ 的随机向量 $X$，其**特征函数 (cf)** 定义为：

    $$
    \phi_X(t) = E[e^{itX}] = \int e^{itx} dF(x)
    $$

    利用欧拉公式，它可以展开为实部与虚部：

    $$
    \phi_X(t) = E[\cos tX] + iE[\sin tX], \quad \text{for any } t \in \mathbb{R}
    $$

    > **注**：相较于矩母函数 (MGF) $M_X(t) = E[e^{tX}]$ 可能在某些分布（如柯西分布）下不存在，特征函数由于 $|e^{itX}| = 1$，对**任何**概率分布都必然存在。

特征函数继承了许多极其优美的数学性质，这些性质在后续的极限定理推导中起着决定性作用。

!!! info "特征函数的基本性质 (Properties of CFs)"

    对于一维随机变量 $X$，其特征函数 $\phi_X(t)$ 满足：

    * **(i) 有界性**：$|\phi_X(t)| \le \phi_X(0) = 1$。

    * **(ii) 共轭对称性**：$\overline{\phi_X(t)} = \phi_X(-t)$。

    * **(iii) 一致连续性**：$\phi_X(t)$ 在 $\mathbb{R}$ 上是一致连续的。

    * **(iv) 运算封闭性**：$\overline{\phi_X}$、$|\phi_X|^2$ 和 $Re(\phi_X)$ 分别对应于随机变量 $-X$、$X-Y$（其中 $X,Y$ 独立同分布于 $F$）以及混合分布 $(F_X + F_{-X})/2$ 的特征函数。

    * **(v) 格点分布判定**：如果存在 $t_0 \neq 0$ 使得 $|\phi_X(t_0)| = 1$，那么必然存在 $a \in \mathbb{R}$ 且 $a \neq 0$ 使得 $P(X \in \{a + jh : j \in \mathbb{Z}\}) = 1$。即 $X$ 是一个格点 (lattice) 随机变量。

    * **(vi) 黎曼-勒贝格引理 (Riemann-Lebesgue)**：如果 $F$ 是绝对连续的（即存在密度函数），则 $\lim_{|t|\to\infty} |\phi_X(t)| = 0$。

    * **(vii) 唯一性与傅里叶逆变换**：两个随机变量同分布 $X \stackrel{d}{=} Y$ 当且仅当 $\forall t, \phi_X(t) = \phi_Y(t)$。若 $\phi_X$ 绝对可积（即 $\phi_X \in \mathcal{L}^1(\mathbb{R})$），则 $F$ 具有连续密度函数，且可由逆变换求得：
    
    $$
    f(x) = \frac{1}{2\pi} \int e^{-itx} \phi_X(t) dt
    $$

    ??? proof "前四条基本性质的证明补充（点击展开）"

        **证明 (i)**：
        利用积分的绝对值不等式：
        
        $$
        |\phi_X(t)| = \left| E[e^{itX}] \right| \le E[|e^{itX}|] = E[1] = 1 = \phi_X(0)
        $$

        **证明 (ii)**：
        由复数共轭的性质：
        
        $$
        \overline{\phi_X(t)} = \overline{E[\cos tX + i\sin tX]} = E[\cos tX - i\sin tX] = E[\cos(-tX) + i\sin(-tX)] = \phi_X(-t)
        $$

        **证明 (iii)**：
        对于任意 $t$ 和增量 $h$：
        
        $$
        |\phi_X(t+h) - \phi_X(t)| = \left| E[e^{i(t+h)X} - e^{itX}] \right| \le E\left[ |e^{itX}| \cdot |e^{ihX} - 1| \right] = E[|e^{ihX} - 1|]
        $$
        
        由于 $|e^{ihX} - 1| \le 2$（有界），且当 $h \to 0$ 时 $e^{ihX} - 1 \to 0$。由控制收敛定理 (DCT)，上式期望趋于 0，且该极限与 $t$ 无关，故为一致连续。

        **证明 (iv) 中的 $|\phi_X|^2$ 性质**：
        设 $X, Y$ 独立同分布。则 $X-Y$ 的特征函数为：
        
        $$
        \phi_{X-Y}(t) = E[e^{it(X-Y)}] = E[e^{itX}] E[e^{-itY}] = \phi_X(t) \phi_Y(-t)
        $$
        
        因为 $X, Y$ 同分布，$\phi_Y(-t) = \phi_X(-t) = \overline{\phi_X(t)}$，故：
        
        $$
        \phi_{X-Y}(t) = \phi_X(t) \overline{\phi_X(t)} = |\phi_X(t)|^2
        $$

---

## 2. 多元特征函数 (Multivariate CF)

上述概念可以自然地推广到高维空间。

!!! abstract "定义：多元特征函数"

    设 $X$ 为 $p$ 维随机向量，其特征函数定义为：

    $$
    \phi_X(t) = E[e^{it^\top X}] = \int_{\mathbb{R}^p} e^{it^\top x} dF_X(x), \quad \text{for any } t \in \mathbb{R}^p
    $$

多元特征函数完美继承了一维特征函数的性质，并增加了矩阵相关的微积分性质：

* **仿射变换**：对于标量 $b \neq 0$，$\phi_{X/b}(t) = \phi_X(t/b)$；对于常数向量 $c$，$\phi_{X+c}(t) = \exp\{it^\top c\} \phi_X(t)$。

* **独立性和**：若 $X$ 和 $Y$ 独立，则 $\phi_{X+Y}(t) = \phi_X(t)\phi_Y(t)$。

* **矩与导数的关系**：
  * 若 $E\|X\| < \infty$，则梯度 $\nabla \phi_X(t)$ 存在且连续，且 $\nabla \phi_X(0) = i\mu$（其中 $\mu = EX$）。
  * 若 $E\|X\|^2 < \infty$，则 Hessian 矩阵 $\nabla^2 \phi_X(t)$ 存在且连续，且 $\nabla^2 \phi_X(0) = -E[XX^\top]$。

* **多元正态分布特例**：

  若 $X \sim N_d(\mu, \Sigma)$，其特征函数为极其优美的二次型指数形式：

  $$
  \phi_X(t) = \exp\left\{ it^\top \mu - \frac{1}{2} t^\top \Sigma t \right\}
  $$

---

## 3. Lévy 连续性定理与极限应用

特征函数最强大的应用在于它将**概率测度的收敛 (Weak Convergence)** 转化为了**复值函数的逐点收敛**。

!!! success "定理 2.2：Lévy-Cramér 定理 (Lévy's Continuity Theorem)"

    设 $\{X_n\}$ 和 $X$ 为 $\mathbb{R}^d$ 中的随机向量。则：

    $$
    X_n \xrightarrow{d} X \iff \phi_{X_n}(t) \to \phi_X(t), \quad \forall t \in \mathbb{R}^d
    $$

    ??? proof "基于 Portmanteau 引理的证明提示（点击展开）"

        **$\Rightarrow$ 方向**：
        由于复指数函数 $e^{it^\top x} = \cos(t^\top x) + i\sin(t^\top x)$ 是**有界且连续**的。
        直接应用第一章中的 Portmanteau 引理 (ii)：对任意有界连续函数 $f \in C_B$，$E[f(X_n)] \to E[f(X)]$ 成立，故特征函数必然逐点收敛。
        
        **$\Leftarrow$ 方向**：
        这是定理的难点。核心思路是先利用特征函数在原点附近的连续性证明序列 $\{X_n\}$ 是**胎紧的 (Tight)**。
        由 Prokhorov 定理，胎紧序列必然存在收敛子列。再利用特征函数的唯一性定理，证明所有收敛子列的极限分布都必定与 $X$ 的分布相同，从而得出整个序列依分布收敛于 $X$。

有了这个定理，证明大数定律 (WLLN) 和中心极限定理 (CLT) 就变成了纯粹的代数展开。

!!! tip "应用 1：泊松分布的中心极限定理"

    假设 $X_1, \dots, X_n$ 独立同分布于 $Poisson(\lambda)$。
    我们已知 $X_j$ 的特征函数为 $\phi_X(t) = \exp\{\lambda(e^{it}-1)\}$。
    令 $\overline{X} = n^{-1}\sum X_i$，我们来考察标准化统计量 $\frac{\overline{X} - \lambda}{\sqrt{\lambda/n}}$ 的特征函数：

    ??? proof "推导过程（点击展开）"

        利用仿射变换和独立性性质：
        
        $$
        \phi_{\frac{\overline{X} - \lambda}{\sqrt{\lambda/n}}}(t) = \exp\{-it\sqrt{n\lambda}\} \cdot \phi_{\overline{X}}\left(\frac{t}{\sqrt{\lambda/n}}\right) = \exp\{-it\sqrt{n\lambda}\} \cdot \phi_X^n\left(\frac{t}{\sqrt{n\lambda}}\right)
        $$
        
        代入 Poisson 特征函数：
        
        $$
        = \exp\{-it\sqrt{n\lambda}\} \cdot \exp\left\{ n\lambda \left( e^{\frac{it}{\sqrt{n\lambda}}} - 1 \right) \right\}
        $$
        
        对内部的指数函数进行泰勒展开 $e^x = 1 + x + x^2/2 + o(x^2)$：
        
        $$
        = \exp\left\{ -it\sqrt{n\lambda} + n\lambda \left( \frac{it}{\sqrt{n\lambda}} + \frac{i^2 t^2}{2n\lambda} + o\left(\frac{1}{n\lambda}\right) \right) \right\}
        $$
        
        展开并消去一次项：
        
        $$
        = \exp\left\{ -it\sqrt{n\lambda} + it\sqrt{n\lambda} - \frac{t^2}{2} + o(1) \right\} = \exp\left\{ -t^2/2 + o(1) \right\}
        $$
        
        当 $n \to \infty$ 时，该特征函数收敛于 $e^{-t^2/2}$，这正是标准正态分布 $N(0,1)$ 的特征函数。
        故由 Lévy 连续性定理：
        
        $$
        \frac{\overline{X} - \lambda}{\sqrt{\lambda/n}} \xrightarrow{d} N(0, 1)
        $$

!!! tip "应用 2：弱大数定律 (WLLN)"

    设 $Y_1, \dots, Y_n$ 是 i.i.d. 随机变量，且 $\phi_Y(t)$ 在 $t=0$ 处可导，导数为 $i\mu = \phi'(0)$（这等价于存在有限一阶矩）。那么样本均值 $\overline{Y} \xrightarrow{P} \mu$。

    ??? proof "推导过程（点击展开）"

        由于 $\phi(0)=1$ 且 $\phi'(0)$ 存在，在 $t \to 0$ 时有泰勒展开：
        
        $$
        \phi_Y(t) = 1 + t\phi'(0) + o(t)
        $$
        
        考察样本均值 $\overline{Y}$ 的特征函数：
        
        $$
        \phi_{\overline{Y}}(t) = \phi_Y^n\left(\frac{t}{n}\right) = \left( 1 + \frac{t}{n}\phi'(0) + o\left(\frac{t}{n}\right) \right)^n
        $$
        
        代入 $\phi'(0) = i\mu$：
        
        $$
        = \left( 1 + \frac{it\mu}{n} + o\left(\frac{1}{n}\right) \right)^n
        $$
        
        利用微积分中的极限公式 $\lim_{n \to \infty} (1 + x/n)^n = e^x$：
        
        $$
        \lim_{n \to \infty} \phi_{\overline{Y}}(t) = e^{it\mu}
        $$
        
        这是退化分布（即常数 $\mu$）的特征函数。故 $\overline{Y} \xrightarrow{d} \mu$。又因为向常数依分布收敛等价于依概率收敛，得证 $\overline{Y} \xrightarrow{P} \mu$。$\square$

---

## 4. 矩与特征函数的泰勒展开

通过上节可知，渐近理论的核心在于**特征函数的泰勒展开**。这直接与随机变量的矩相关联。

如果随机变量 $X$ 的 $r$ 阶矩存在，那么 $\phi_X(t)$ 就是 $r$ 阶可导的，且：

$$
\phi_X^{(r)}(t) = \int (ix)^r e^{itx} dF(x) = E[(iX)^r e^{itX}]
$$

这导致了在原点的导数值直接给出原点矩：$\phi_X^{(r)}(0) = i^r E[X^r]$。

!!! info "定理 2.3：特征函数的展开式"

    如果 $E|X|^r < \infty$，那么其特征函数可以展开为：

    $$
    \phi_X(t) = \sum_{j=0}^r \frac{(it)^j}{j!} E[X^j] + o(|t|^r)
    $$

    > **注记 (The Moment Problem)**：
    > 特征函数决定了 $X$ 的所有矩。但是反过来，一个序列的所有矩 $\{m_r := E[X^r]\}_{r=1}^\infty$ 能否唯一确定 $X$ 的分布呢？
    > 这被称为**矩问题 (Moment Problem)**。答案是：**不能**。只有当满足 **Carleman 条件** 时，分布才能被唯一确定：
    > 
    > $$
    > \sum_{r=1}^\infty m_{2r}^{-\frac{1}{2r}} = +\infty
    > $$

利用高阶泰勒展开，我们同样可以极其简洁地证明一般情况下的中心极限定理：

??? proof "一般分布的中心极限定理 (CLT) 证明（点击展开）"

    假设 $X_1, \dots, X_n$ i.i.d.，且均值 $\mu = E[X]$，方差 $\sigma^2 = E[X^2] < \infty$。
    令中心化变量 $Y = X - \mu$，则 $E[Y]=0, E[Y^2]=\sigma^2$。
    其特征函数展开到二阶为：
    
    $$
    \phi_{X-\mu}(t) = 1 + \frac{1}{2}(it)^2 \sigma^2 + o(t^2) = 1 - \frac{t^2 \sigma^2}{2} + o(t^2)
    $$
    
    对于标准化和 $Z_n = \frac{n\overline{X} - n\mu}{\sqrt{n\sigma^2}}$，其特征函数为：
    
    $$
    \phi_{Z_n}(t) = \phi_{X-\mu}^n\left(\frac{t}{\sigma\sqrt{n}}\right) = \left( 1 - \frac{1}{2}\left(\frac{t}{\sigma\sqrt{n}}\right)^2 \sigma^2 + o\left(\frac{t^2}{\sigma^2 n}\right) \right)^n
    $$
    
    化简后：
    
    $$
    = \left( 1 - \frac{t^2}{2n} + o\left(\frac{1}{n}\right) \right)^n \xrightarrow{n \to \infty} e^{-t^2/2}
    $$
    
    由 Lévy-Cramér 定理即得 $Z_n \xrightarrow{d} N(0,1)$。$\square$

---

## 5. 累积量 (Cumulants) 与 Edgeworth 展开

如果我们继续对特征函数进行高阶展开 $(r > 2)$，比如展开到第四阶：

$$
\phi_{\frac{\overline{X}-\mu}{\sqrt{\sigma^2/n}}}(t) = \left( 1 - \frac{t^2}{2n} - \frac{i t^3}{6 n^{3/2}} \left(\frac{m_3}{\sigma}\right)^3 + \frac{t^4}{24 n^2} \left(\frac{m_4}{\sigma}\right)^4 + \dots \right)^n
$$

这导致了极其复杂的代数表达。为了简化这种针对 $n$ 个独立同分布变量求和的展开，我们引入**累积量 (Cumulants 或 Semi-Invariants)**。

!!! abstract "定义 2.4：累积量生成函数 (Cumulant Generating Function)"

    我们不对 $\phi_X(t)$ 本身进行展开，而是对其对数 $K_X(t) = \log \phi_X(t)$ 进行泰勒展开，展开式的系数 $\kappa_j$ 即为**累积量**：

    $$
    K_X(t) := \log \phi_X(t) = \sum_{j \ge 1} \frac{(it)^j}{j!} \kappa_j = \log \left\{ 1 + \sum_{j \ge 1} \frac{1}{j!} m_j (it)^j \right\}
    $$

    利用 $\log(1+x) = x - x^2/2 + x^3/3 - \dots$ 的级数匹配系数，我们可以得到矩与累积量的转换关系（设 $\kappa_1 = m_1 = EX$）：

    * $\kappa_2 = m_2 - m_1^2 = E(X - EX)^2 =: c_2$ （即方差）
    * $\kappa_3 = m_3 - 3m_1 m_2 + 2m_1^3 = E(X - EX)^3 =: c_3$
    * $\kappa_4 = m_4 - 4m_1 m_3 - 3m_2^2 + 12m_1^2 m_2 - 6m_1^4 = c_4 - 3c_2^2$

    > **注**：高阶 $(j > 3)$ 的累积量不同于中心矩。对于标准化变量 $Y_i = (X_i - \mu)/\sigma$，其 $\kappa_1=0, \kappa_2=1$。而 $\kappa_3$ 被称为**偏度 (Skewness)**，$\kappa_4$ 被称为**峰度 (Kurtosis)**。

展开 $\log \phi(t)$ 的巨大好处在于独立变量相加时，**累积量是直接线性相加的**！

### Edgeworth 展开 (Edgeworth Expansion)

通过累积量，我们可以将标准化和 $S_n = \frac{\overline{X}-\mu}{\sqrt{\sigma^2/n}}$ 的特征函数写为：

$$
\phi_{S_n}(t) = \phi_Y^n\left(\frac{t}{\sqrt{n}}\right) = \exp\left\{ -\frac{t^2}{2} + \sum_{j \ge 3} \kappa_j \frac{(it)^j}{j!} n^{-\frac{j}{2}+1} \right\}
$$

将指数项重新按 $n^{-1/2}$ 的幂次展开：

$$
= e^{-t^2/2} \left\{ 1 + \sum_{j \ge 1} n^{-\frac{j}{2}} r_j(it) \right\}
$$

其中 $r_j(\cdot)$ 是实系数多项式，最高次数为 $3j$（例如 $r_1(u) = \frac{1}{6}\kappa_3 u^3$）。

!!! success "高阶渐近逼近：Edgeworth 展开"

    利用傅里叶逆变换思想，既然特征函数可以写成上述多项式乘以正态特征函数的形式，那么**累积分布函数 $P(S_n \le x)$ 必然也能写成标准正态 CDF $\Phi(x)$ 的修正形式**：

    $$
    P(S_n \le x) = \Phi(x) + n^{-\frac{1}{2}} R_1(x) + n^{-1} R_2(x) + \dots
    $$

    这被称为 **Edgeworth 展开**。它比单纯的 CLT 提供了更精确的收敛速率和有限样本修正。

    ??? proof "修正项 $R_j(x)$ 与 Hermite 多项式的计算（点击展开）"

        为了求解 $R_j(x)$，我们需要找到一个函数，使得它的傅里叶-斯蒂尔杰斯变换恰好等于 $e^{-t^2/2} r_j(it)$：
        
        $$
        e^{-t^2/2} r_j(it) = \int e^{itx} dR_j(x)
        $$

        我们利用标准正态分布的性质，反复使用分部积分：
        
        $$
        e^{-t^2/2} = (-it)^{-j} \int e^{itx} d\Phi^{(j)}(x)
        $$

        定义微分算子 $D = d/dx$。这启发我们将 $r_j(it)$ 替换为微分算子多项式 $r_j(-D)$，作用于 $\Phi(x)$ 上：
        
        $$
        \int e^{itx} d\{r_j(-D)\Phi(x)\} = r_j(it) e^{-t^2/2}
        $$

        这意味着：
        
        $$
        R_j(x) = r_j(-D)\Phi(x)
        $$

        而正态分布的各阶导数恰好由著名的 **Hermite 多项式 $He_{j}(x)$** 生成：
        
        $$
        (-D)^j \Phi(x) = -He_{j-1}(x) e^{-t^2/2} \cdot \frac{1}{\sqrt{2\pi}}
        $$

        因此，$R_j(x)$ 可以由标准正态密度函数及其 Hermite 多项式精确表达。这是高阶渐近理论（如 Bootstrap 理论）中极其基础的工具。$\square$