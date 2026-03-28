# 第一章：条件期望Conditional Expectation

在开始随机微分方程前，我们需要抛弃初等概率论中条件期望的定义方法，用测度论的语言严格定义条件期望。这是我们后续探索鞅理论、随机分析的基础。

> 在测度论中，$\sigma$-代数是对“信息”的数学刻画
对于集合 $\Omega$ 上的一个子集族 $\mathcal{F}$，若满足以下条件，则称其为一个 **$\sigma$-代数**：
1.  **包含全集**：$\Omega \in \mathcal{F}$。
2.  **对补集封闭**：若 $A \in \mathcal{F}$，则 $A^c \in \mathcal{F}$。
3.  **对可列并封闭**：若 $A_1, A_2, \dots \in \mathcal{F}$，则 $\bigcup_{n=1}^\infty A_n \in \mathcal{F}$。

> Radon-Nikodym 定理：若测度 $\nu$ 关于 $P$ 是 **绝对连续 (Absolutely Continuous)** 的，即满足 $\nu \ll P$（若 $P(A)=0$ 则 $\nu(A)=0$），则存在一个 $P$-几乎处处唯一的 $\mathcal{F}$-可测函数 $f: \Omega \to [0, \infty)$，使得对于任意 $A \in \mathcal{F}$，有：
> $$
> \nu(A) = \int_A f(\omega) \, dP(\omega)
> $$

> Radon-Nikodym 导数：称上述可测函数 $f$ 为 $\nu$ 关于 $P$ 的 **Radon-Nikodym 导数**（也称为 $\nu$ 的密度函数），记作：$f = \frac{d\nu}{dP}$

## 1. 条件期望的测度论定义

在古典概率论中，定义事件 $A$ 在事件 $B$ 发生的条件下的条件概率为 $P(A|B) \triangleq \frac{P(AB)}{P(B)}.$

相应的条件期望为 $E[X|B] \triangleq \frac{1}{P(B)}\int_B X dP.$

但当给定的是一个连续随机变量 $Y$ 时，事件 $\{Y=y\}$ 的概率为 0，古典定义失效。

此时我们需要借助 Radon-Nikodym 定理在 $\sigma$-代数层面重新定义。

!!! abstract "定义：条件期望 (Conditional Expectation)"

    设 $X$ 为概率空间 $(\Omega, \mathcal{F}, P)$ 上的可积随机变量（即 $E[|X|] < \infty$），$\mathcal{G}$ 是 $\mathcal{F}$ 的一个子 $\sigma$-代数（常记为 $\mathcal{G} = \sigma(Y)$）。
    我们称随机变量 $Z$ 为 $X$ 给定 $\mathcal{G}$ 下的**条件期望**，记作 $E[X|\mathcal{G}] = Z$，如果它满足以下两个条件：

    **1. 可测性**：$Z$ 是 $\mathcal{G}$-可测的随机变量。

    **2. 局部平均相等（积分匹配）**：对于任意的 $A \in \mathcal{G}$，都有：

    $$
    \int_A X dP = \int_A Z dP
    $$

    由 Radon-Nikodym 定理保证，这样的 $Z$ 存在且在 almost sure (a.s.) 的意义下唯一。

!!! note "补充：初等定义与测度论定义的过渡性质"
    基于随机变量 $Y$ 产生的 $\sigma$-代数，我们有以下三个直观的基本性质：
    1. $E[X|Y] = E[X|\sigma(Y)]$ （对随机变量求条件期望，本质是对其产生的 $\sigma$-代数求条件期望）。
    2. $E\big[E[X|\mathcal{G}]\big] = E[X]$ （保期望性质）。
    3. 若 $\mathcal{G} = \left\{\Omega，\emptyset\right\}$ ，则 $E[X|\mathcal{G}] = X$ a.s. （平凡的 $\sigma$-代数无信息）。

## 2. 核心性质

> 绝大部分定义证明或者简单变式即可。

设 $\mathcal{G}, \mathcal{H}$ 为子 $\sigma$-代数。

!!! info "基本性质"

    **1. 线性性**：

    $$
    E[aX + bY | \mathcal{G}] = aE[X|\mathcal{G}] + bE[Y|\mathcal{G}] \quad a.s.
    $$

    **2. 保期望性**：

    如果 $X$ 是 $\mathcal{G}$-可测的，那么：

    $$
    E[X|\mathcal{G}] = X \quad a.s.
    $$

    > *符合投影直观，可测意味着完全在$\mathcal{G}$张成的超平面上。*

    **3. 已知即常数（提出已知因子）**：

    如果 $X$ 是 $\mathcal{G}$-可测的，那么：

    $$
    E[XY|\mathcal{G}] = X E[Y|\mathcal{G}] \quad a.s.
    $$

    **4. 独立则无关**：

    如果 $X$ 与 $\mathcal{G}$ 独立，那么：

    $$
    E[X|\mathcal{G}] = E[X] \quad a.s.
    $$

    > *符合投影直观，独立意味着垂直$\mathcal{G}$张成的超平面上。*

    **5. 塔牌性质 (Tower Property)**：

    若 $\mathcal{H} \subset \mathcal{G} \subset \mathcal{F}$，则：

    $$
    E\big[ E[X|\mathcal{G}] \big| \mathcal{H} \big] = E\big[ E[X|\mathcal{H}] \big| \mathcal{G} \big] = E[X|\mathcal{H}] \quad a.s.
    $$

    > *取信息少的。*

    ??? proof "后三个核心性质的完整推导"

        **3. 提出已知因子证明 (Standard Machine)**：
        第一步：设 $X = \chi_B$（其中 $B \in \mathcal{G}$）。对于任意测试集 $A \in \mathcal{G}$，由于 $A \cap B \in \mathcal{G}$，有：

        $$
        \int_A \chi_B E[Y|\mathcal{G}] dP = \int_{A \cap B} E[Y|\mathcal{G}] dP = \int_{A \cap B} Y dP = \int_A \chi_B Y dP
        $$

        第二步：由线性性推广到简单函数 $X = \sum a_i \chi_{B_i}$。
        第三步：对非负可测函数 $X \ge 0, Y \ge 0$，存在非负简单函数列 $X_n \uparrow X$，由单调收敛定理极限可与积分交换。一般情况拆分正负部即得证。

        **4. 独立则无关证明**：
        常数 $E[X]$ 自然是 $\mathcal{G}$-可测的。对于任意 $A \in \mathcal{G}$，由于 $X$ 与 $\mathcal{G}$ 独立，所以 $X$ 与指示函数 $\chi_A$ 独立：

        $$
        \int_A X dP = E[X \cdot \chi_A] = E[X] \cdot E[\chi_A] = E[X] P(A) = \int_A E[X] dP
        $$

        完全符合条件期望定义。

        **5. 塔牌性质证明**：
        设 $Z = E[X|\mathcal{G}]$。因为子代数嵌套 $\mathcal{H} \subset \mathcal{G}$，对于任意 $A \in \mathcal{H}$，必然有 $A \in \mathcal{G}$。
        由 $Z$ 的定义，对该集合 $A$ 有：

        $$
        \int_A Z dP = \int_A X dP
        $$

        而由 $E[X|\mathcal{H}]$ 的定义，对同一个 $A \in \mathcal{H}$ 有：

        $$
        \int_A E[X|\mathcal{H}] dP = \int_A X dP
        $$

        两式结合即得 $\int_A Z dP = \int_A E[X|\mathcal{H}] dP$。由于等式对所有 $A \in \mathcal{H}$ 成立，且二者均 $\mathcal{H}$-可测，故 a.s. 相等。

!!! tip "几何解释：$L^2$ 空间上的正交投影"

    条件期望有一个极其优美的几何直观。

    考虑平方可积空间 $L^2(\Omega, \mathcal{F}, P)$，它是一个 Hilbert 空间，内积定义为：

    $$
    (X,Y) = E[XY]
    $$

    由于 $\mathcal{G} \subset \mathcal{F}$，子空间 $L^2(\Omega, \mathcal{G}, P)$ 也是一个闭子空间。

    在这个视角下，**条件期望 $E[X|\mathcal{G}]$ 实际上就是元素 $X$ 在子空间 $L^2(\Omega, \mathcal{G}, P)$ 上的正交投影**。

    它使得均方误差 $E[(X - Z)^2]$ 在所有 $\mathcal{G}$-可测的随机变量 $Z$ 中达到最小：

    $$
    \min_{Z \in \mathcal{G}} E[(X - Z)^2] = E\big[(X - E[X|\mathcal{G}])^2\big]
    $$

## 3. Jensen 不等式

Jesen不等式是高等概率论中非常重要的关于函数凸性的不等式。

!!! success "定理：条件 Jensen 不等式"

    设 $\Phi: \mathbb{R} \rightarrow \mathbb{R}$ 为凸函数，且 $E[|\Phi(X)|] < \infty$，则有：

    $$
    \Phi(E[X|\mathcal{G}]) \le E[\Phi(X)|\mathcal{G}] \quad a.s.
    $$

    **证明思路 (基于简单函数逼近)**：
    利用凸函数的离散性质 $\Phi(\sum a_i b_i) \le \sum b_i \Phi(a_i)$（其中 $b_i \ge 0, \sum b_i = 1$）。由于条件期望可以看作是一种加权平均，我们可以用指示函数构造简单函数来逼近，利用实变函数中从指示函数到简单函数再到可测函数的标准途径来证明。

    ??? proof "Jensen 不等式的证明"
        
        **第一步：对简单函数证明**
        假设 $X$ 是一个简单函数，可表示为：

        $$
        X = \sum_{i=1}^n a_i \chi_{B_i}
        $$

        其中 $B_i$ 互不相交且 $\bigcup_{i=1}^n B_i = \Omega$。
        取关于 $\mathcal{G}$ 的条件期望，由线性性：

        $$
        E[X|\mathcal{G}] = \sum_{i=1}^n a_i E[\chi_{B_i}|\mathcal{G}]
        $$

        令 $b_i = E[\chi_{B_i}|\mathcal{G}]$。显然 $b_i \ge 0$，且 $\sum_{i=1}^n b_i = E\big[\sum_{i=1}^n \chi_{B_i} \big| \mathcal{G}\big] = E[1|\mathcal{G}] = 1$。
        这说明 $b_i$ 构成了一组凸组合的权重。
        根据凸函数的离散性质 $\Phi(\sum a_i b_i) \le \sum b_i \Phi(a_i)$，我们有：

        $$
        \Phi(E[X|\mathcal{G}]) = \Phi\left(\sum_{i=1}^n a_i b_i\right) \le \sum_{i=1}^n b_i \Phi(a_i)
        $$

        将 $b_i$ 的定义代回：

        $$
        \sum_{i=1}^n E[\chi_{B_i}|\mathcal{G}] \Phi(a_i) = E\left[\sum_{i=1}^n \Phi(a_i) \chi_{B_i} \bigg| \mathcal{G}\right] = E[\Phi(X)|\mathcal{G}]
        $$

        故对简单函数，结论成立。

        **第二步：推广至一般可测函数**
        对于一般的可积随机变量 $X$，我们可以找到一列简单函数 $X_n$ 使得 $X_n \to X$ a.s.。
        利用极限过程（控制收敛定理或单调收敛定理），并结合 $\Phi$ 作为凸函数的连续性，两边取极限即可将不等式推广到一般情况。$\square$


## 4. 鞅 (Martingale) 与 Doob 不等式

在研究随机微分方程前，我们需要从离散时间步入连续时间，引入信息流和鞅的概念。

!!! abstract "定义：离散鞅 与 连续信息流"

    **1. 离散时间的鞅 (Discrete Martingale)**：
    设 $X_1, X_2, \dots$ 为一列随机变量，且满足 $E[|X_i|] < \infty$。如果对于任意 $j > i$，都有：

    $$
    E[X_j | X_1, \dots, X_i] = X_i \quad a.s.
    $$

    则称 $\{X_n\}$ 为一个**离散鞅**。（代表公平游戏，已知历史信息下，未来的期望等于现在）。
    若 $E[X_j | X_1, \dots, X_i] \ge X_i$，称为**下鞅 (Submartingale)**。

    **2. 连续时间的信息流 (Filtration)**：
    在连续时间 $t \in [0, \infty)$ 下，我们定义一族单调递增的 $\sigma$-代数 $\{\mathcal{F}_t\}_{t \ge 0}$，即当 $s \le t$ 时，$\mathcal{F}_s \subset \mathcal{F}_t$。它代表了随着时间 $t$ 累积的“历史信息”。（例如布朗运动自然生成的流 $\mathcal{F}_t = \sigma(B_s, 0 \le s \le t)$）。

    **3. 连续鞅**：设随机过程 $\{X_t\}$ 适应于流 $\mathcal{F}_t$，如果对于任意 $s \le t$，都有 $E[X_t | \mathcal{F}_s] = X_s \ a.s.$，则称 $\{X_t\}$ 为一个**连续鞅**。


在随机分析中，我们常常需要控制整个过程在其路径上的最大值（例如研究布朗运动的最大位移）。Doob 不等式成为解决此类问题最重要的工具。

!!! info "定理：Doob 极大值不等式 (Doob's Maximal Inequality)"

    **(1) 次线性界**：设 $\{X_n\}_{n=1}^N$ 为非负下鞅，则对于任意 $\lambda > 0$：

    $$
    P\left(\max_{1 \le k \le n} X_k \ge \lambda\right) \le \frac{1}{\lambda} E[X_n I_{\{\max X_k \ge \lambda\}}] \le \frac{1}{\lambda} E[X_n]
    $$

    > *(注意它与 Chebyshev 不等式的形式相似，但它强大在于控制的是整个路径的极大值，而不仅是单点)*

    **(2) $L^p$ 极大不等式**：设 $p > 1$ 且 $\frac{1}{p} + \frac{1}{q} = 1$。如果 $\{X_n\}$ 是非负下鞅且 $X_n \in L^p$，则：

    $$
    E\left[\max_{1 \le k \le n} X_k^p\right] \le \left(\frac{p}{p-1}\right)^p E[X_n^p]
    $$

    ??? proof "Doob 不等式完整推导（包含次线性界与 $L^p$ 边界）"

        **第一部分：次线性界证明**
        设 $\{X_n\}$ 是非负下鞅。对于任意 $\lambda > 0$，我们定义首次越界集 $A_k$（即在时刻 $k$ 第一次突破 $\lambda$）：

        $$
        A_k = \{X_1 < \lambda, \dots, X_{k-1} < \lambda, X_k \ge \lambda\}
        $$

        1. 显然 $A_k$ 互不相交，且它们的并集正是我们要估计的事件：

        $$
        A = \bigcup_{k=1}^n A_k = \left\{\max_{1 \le k \le n} X_k \ge \lambda\right\}
        $$

        2. 由于 $A_k \in \mathcal{F}_k$ 且 $X_n$ 是下鞅（即 $E[X_n | \mathcal{F}_k] \ge X_k$），我们在 $A_k$ 上积分有：

        $$
        \int_{A_k} X_n dP \ge \int_{A_k} X_k dP
        $$

        3. 根据 $A_k$ 的定义，在集合 $A_k$ 上必然有 $X_k \ge \lambda$，因此：

        $$
        \int_{A_k} X_k dP \ge \int_{A_k} \lambda dP = \lambda P(A_k)
        $$

        4. 将所有的 $k = 1, \dots, n$ 加起来。由于 $A_k$ 是不交并，且 $X_n \ge 0$：

        $$
        E[X_n] \ge \int_A X_n dP = \sum_{k=1}^n \int_{A_k} X_n dP \ge \lambda \sum_{k=1}^n P(A_k) = \lambda P(A)
        $$

        由此即得 $P(\max X_k \ge \lambda) \le \frac{1}{\lambda} \int_A X_n dP \le \frac{1}{\lambda} E[X_n]$。次线性界得证。

        ---

        **第二部分：$L^p$ 边界推导**
        设极大值变量为 $X^* = \max_{1 \le k \le n} X_k$。利用非负随机变量期望的积分恒等式：

        $$
        E[(X^*)^p] = \int_0^\infty p \lambda^{p-1} P(X^* \ge \lambda) d\lambda
        $$

        将第一部分证明的更紧的边界 $P(X^* \ge \lambda) \le \frac{1}{\lambda} \int_{\{X^* \ge \lambda\}} X_n dP$ 代入：

        $$
        E[(X^*)^p] \le \int_0^\infty p \lambda^{p-2} \left( \int_{\{X^* \ge \lambda\}} X_n dP \right) d\lambda
        $$

        **关键点 1：Fubini 定理交换积分顺序**。先对 $\lambda$ 积分（注意积分上限受制于 $\lambda \le X^*$）：

        $$
        = \int_{\Omega} X_n \left( \int_0^{X^*} p \lambda^{p-2} d\lambda \right) dP = \int_{\Omega} X_n \left( \frac{p}{p-1} (X^*)^{p-1} \right) dP
        $$

        提出常数，得到：

        $$
        E[(X^*)^p] \le \frac{p}{p-1} E\left[X_n (X^*)^{p-1}\right]
        $$

        **关键点 2：Hölder 不等式剥离**。对上式右侧的乘积应用 Hölder 不等式（共轭指数 $\frac{1}{p} + \frac{1}{q} = 1$）：

        $$
        E\left[X_n (X^*)^{p-1}\right] \le \left( E[X_n^p] \right)^{\frac{1}{p}} \left( E\left[((X^*)^{p-1})^q\right] \right)^{\frac{1}{q}}
        $$

        由于 $\frac{1}{p} + \frac{1}{q} = 1$，所以 $q = \frac{p}{p-1}$，进而在第二个因子中 $(p-1)q = p$：

        $$
        \left( E\left[((X^*)^{p-1})^q\right] \right)^{\frac{1}{q}} = \left( E[(X^*)^p] \right)^{1 - \frac{1}{p}}
        $$

        将 Hölder 的结果代回原式：

        $$
        E[(X^*)^p] \le \frac{p}{p-1} \left( E[X_n^p] \right)^{\frac{1}{p}} \left( E[(X^*)^p] \right)^{1 - \frac{1}{p}}
        $$

        假设 $E[(X^*)^p] < \infty$（严格证明中可通过截断处理），两边同时除以 $\left( E[(X^*)^p] \right)^{1 - 1/p}$：

        $$
        \left( E[(X^*)^p] \right)^{\frac{1}{p}} \le \frac{p}{p-1} \left( E[X_n^p] \right)^{\frac{1}{p}}
        $$

        两边同取 $p$ 次方，即得最终结论：

        $$
        E\left[\max_{1 \le k \le n} X_k^p\right] \le \left(\frac{p}{p-1}\right)^p E[X_n^p] \quad \square
        $$


## 5. Borel-Cantelli 引理 (Borel-Cantelli Lemmas)

在研究随机过程的轨道性质和几乎必然 (almost surely, a.s.) 收敛时，我们需要一种能将“概率的级数求和”转化为“事件发生频率”的强大工具。这就是测度论中极其著名的 Borel-Cantelli 引理。

在介绍引理前，我们需要先定义序列事件的“上极限”：

!!! abstract "定义：上极限集 (Limit Superior of Events)"

    设 $\{A_n\}_{n=1}^\infty$ 为概率空间 $(\Omega, \mathcal{F}, P)$ 中的一列事件。我们定义这列事件的**上极限 (limsup)** 为：

    $$
    \limsup_{n \to \infty} A_n = \bigcap_{n=1}^\infty \bigcup_{k=n}^\infty A_k
    $$

    **概率论直觉**：
    样本点 $\omega \in \limsup_{n \to \infty} A_n$ 当且仅当 $\omega$ 属于无穷多个 $A_k$。换句话说，上极限集代表了那些**“无限次发生 (infinitely often, 简记为 i.o.)”** 的事件的集合。
    因此，我们常将其记作 $P(A_n \text{ i.o.})$。

Borel-Cantelli 引理分为两个部分，分别给出了事件无限次发生的充分条件和必要条件。

!!! info "定理：Borel-Cantelli 第一引理 (收敛部分)"

    如果事件列 $\{A_n\}_{n=1}^\infty$ 满足其概率级数收敛，即：

    $$
    \sum_{n=1}^\infty P(A_n) < \infty
    $$

    那么这些事件无限次发生的概率为 0：

    $$
    P(\limsup_{n \to \infty} A_n) = 0
    $$

    > *(注意：第一引理**不需要**事件之间具有任何独立性假设，这是一个极其强大的普适结论！)*

    ??? proof "B-C 第一引理的严格证明（点击展开）"

        设 $B_n = \bigcup_{k=n}^\infty A_k$。显然，序列 $\{B_n\}$ 是一个单调递减的集合序列，即 $B_1 \supset B_2 \supset B_3 \dots$。
        根据概率测度的连续性（连续性由上），上极限的概率可以写为极限的概率：

        $$
        P\left( \bigcap_{n=1}^\infty B_n \right) = \lim_{n \to \infty} P(B_n) = \lim_{n \to \infty} P\left( \bigcup_{k=n}^\infty A_k \right)
        $$

        利用概率的次可加性（Boole 不等式），我们对并集进行放缩：

        $$
        P\left( \bigcup_{k=n}^\infty A_k \right) \le \sum_{k=n}^\infty P(A_k)
        $$

        由于已知整个无穷级数 $\sum_{k=1}^\infty P(A_k)$ 是收敛的，根据微积分中收敛级数的性质，其“尾部余项 (Tail sum)”在 $n \to \infty$ 时必然趋于 0：

        $$
        \lim_{n \to \infty} \sum_{k=n}^\infty P(A_k) = 0
        $$

        将上下两部分结合，由于概率非负，由夹逼定理即可得到：

        $$
        P\left( \limsup_{n \to \infty} A_n \right) = 0 \quad \square
        $$

与第一引理相对，第二引理探讨了当级数发散时的情况，但它额外要求了一个极其苛刻的条件：独立性。

!!! success "定理：Borel-Cantelli 第二引理 (发散部分)"

    如果事件列 $\{A_n\}_{n=1}^\infty$ 是**相互独立 (Mutually Independent)** 的，且其概率级数发散，即：

    $$
    \sum_{n=1}^\infty P(A_n) = \infty
    $$

    那么这些事件无限次发生的概率为 1：

    $$
    P(\limsup_{n \to \infty} A_n) = 1
    $$

    ??? proof "B-C 第二引理的严格证明（点击展开）"

        直接证明某个事件发生概率为 1 往往比较困难，我们转而证明其对立事件的概率为 0。
        利用 De Morgan 定律，上极限集的补集（即 $\liminf$ 下极限）为：

        $$
        \left( \limsup_{n \to \infty} A_n \right)^c = \left( \bigcap_{n=1}^\infty \bigcup_{k=n}^\infty A_k \right)^c = \bigcup_{n=1}^\infty \bigcap_{k=n}^\infty A_k^c
        $$

        这就意味着，我们只需要证明对于任意给定的 $n \ge 1$，都有 $P\left( \bigcap_{k=n}^\infty A_k^c \right) = 0$ 即可。
        
        对于任意的有限整数 $m > n$，由于事件序列 $\{A_k\}$ 相互独立，它们的补集 $\{A_k^c\}$ 也是相互独立的。因此交集的概率等于概率的乘积：

        $$
        P\left( \bigcap_{k=n}^m A_k^c \right) = \prod_{k=n}^m P(A_k^c) = \prod_{k=n}^m (1 - P(A_k))
        $$

        这里我们引入一个在概率论放缩中极其常用的初等不等式：对于任意 $x \ge 0$，都有 $1 - x \le e^{-x}$。
        将其应用到上述连乘中：

        $$
        \prod_{k=n}^m (1 - P(A_k)) \le \prod_{k=n}^m e^{-P(A_k)} = \exp\left( -\sum_{k=n}^m P(A_k) \right)
        $$

        现在，让 $m \to \infty$。由于已知条件给出级数 $\sum_{k=1}^\infty P(A_k) = \infty$，故对于任何固定的起始项 $n$，其部分和 $\sum_{k=n}^\infty P(A_k)$ 也必然趋于 $\infty$。
        因此：

        $$
        P\left( \bigcap_{k=n}^\infty A_k^c \right) = \lim_{m \to \infty} P\left( \bigcap_{k=n}^m A_k^c \right) \le \exp(-\infty) = 0
        $$

        由概率的非负性，该交集的概率必然为 0。
        既然对任意的 $n$ 概率都为 0，那么这些零测集的可列并（即原集合的补集）概率依然为 0：

        $$
        P\left( \bigcup_{n=1}^\infty \bigcap_{k=n}^\infty A_k^c \right) \le \sum_{n=1}^\infty P\left( \bigcap_{k=n}^\infty A_k^c \right) = 0
        $$

        故原事件发生的概率为 $1 - 0 = 1$。证明完毕。$\square$