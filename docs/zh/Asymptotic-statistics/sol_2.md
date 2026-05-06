---
tags:
  - 大样本理论
  - 统计学
  - 课后习题
---

# 📈 课后习题精解：大样本理论 Assignment 2 


## 习题 1: Moments and expansion of CFs

**题目：**
若 $E|X|^{r}<\infty$，证明特征函数 $\phi_{X}(t)$ 具有如下展开式：

\[
\phi_{X}(t)=\sum_{j=0}^{r}\frac{(it)^{j}}{j!}EX^{j}+o(|t|^{r})
\]

??? success "解答 (点击展开)"

    **引理：复指数泰勒展开余项的界限 (Bound on the remainder of the complex exponential Taylor expansion)**
    
    设 $R_{r}(x)=e^{ix}-\sum_{j=0}^{r}\frac{(ix)^{j}}{j!}$ 为 $e^{ix}$ 的 $r$ 阶泰勒展开余项。对于任何实数 $x$ 和整数 $r\ge0$，有：

    \[
    |R_{r}(x)|\le \min \left( \frac{2|x|^{r}}{r!}, \frac{|x|^{r+1}}{(r+1)!} \right)
    \]

    **引理证明：**
    我们对 $r$ 使用数学归纳法。首先建立余项的递归积分关系：

    \[
    \int_{0}^{x}iR_{r-1}(u)du=\int_{0}^{x}i \left( e^{iu}-\sum_{j=0}^{r-1}\frac{(iu)^{j}}{j!} \right) du = \left[ e^{iu}-\sum_{j=0}^{r-1}\frac{(iu)^{j+1}}{(j+1)!} \right]_{0}^{x} = e^{ix}-\sum_{k=0}^{r}\frac{(ix)^{k}}{k!} = R_{r}(x)
    \]

    故得递归关系 $R_{r}(x)=i\int_{0}^{x}R_{r-1}(u)du$。

    * **基础步骤 ($r=0$)**：
        $R_{0}(x)=e^{ix}-1$。一方面，由三角不等式：$|R_{0}(x)|=|e^{ix}-1|\le|e^{ix}|+1=2$。
        另一方面，由积分形式：$|R_{0}(x)|=|i\int_{0}^{x}e^{iu}du|\le\int_{0}^{|x|}|e^{iu}|du=|x|$。
        结合得 $|R_{0}(x)|\le \min(2,|x|)$，命题成立。

    * **归纳步骤**：
        假设不等式对 $r-1$ 成立，即 $|R_{r-1}(x)|\le\frac{2|x|^{r-1}}{(r-1)!}$ 且 $|R_{r-1}(x)|\le\frac{|x|^{r}}{r!}$。
        利用积分关系 $R_{r}(x)=i\int_{0}^{x}R_{r-1}(u)du$，取绝对值：

        \[
        |R_{r}(x)|\le\int_{0}^{|x|}|R_{r-1}(u)|du
        \]

        1. 积分归纳假设的第一部分：$|R_{r}(x)|\le\int_{0}^{|x|}\frac{2u^{r-1}}{(r-1)!}du = \frac{2|x|^{r}}{r!}$。
        2. 积分归纳假设的第二部分：$|R_{r}(x)|\le\int_{0}^{|x|}\frac{u^{r}}{r!}du = \frac{|x|^{r+1}}{(r+1)!}$。
        综上，对所有 $r\ge0$ 引理成立。

    <br>

    **定理证明：**
    由泰勒定理，对任意实数 $x$：

    \[
    e^{ix}=\sum_{j=0}^{r}\frac{(ix)^{j}}{j!}+R_{r}(x)
    \]

    代入 $x=tX$ 并取期望，得到特征函数：

    \[
    \phi_{X}(t)=E[e^{itX}]=\sum_{j=0}^{r}\frac{(it)^{j}}{j!}E[X^{j}]+E[R_{r}(tX)]
    \]

    需证明 $E[R_{r}(tX)]=o(|t|^{r})$，即当 $t\rightarrow0$ 时，$\frac{1}{|t|^{r}}E[R_{r}(tX)]\rightarrow0$。
    定义随机变量 $Z_{t}=\frac{R_{r}(tX)}{|t|^{r}}$。利用引理界限：

    \[
    |Z_{t}|=\frac{|R_{r}(tX)|}{|t|^{r}}\le \min \left( \frac{2|X|^{r}}{r!}, \frac{|t|\cdot|X|^{r+1}}{(r+1)!} \right)
    \]

    * **几乎处处收敛**：对于固定 $X$，当 $t\rightarrow0$ 时，$\frac{|t|\cdot|X|^{r+1}}{(r+1)!}\rightarrow0$，故 $Z_{t} \xrightarrow{a.s.} 0$。
    * **受控性**：对于所有 $t$，恒有 $|Z_{t}|\le\frac{2}{r!}|X|^{r}$。由于 $E|X|^{r}<\infty$，该上界是可积的。
    * **应用控制收敛定理 (DCT)**：

    \[
    \lim_{t\rightarrow0}E[Z_{t}]=E[\lim_{t\rightarrow0}Z_{t}]=E[0]=0
    \]

    因此 $E[R_{r}(tX)]=o(|t|^{r})$。结论证毕。

---

## 习题 2: Corollary of Liapounov CLT

**题目：**
设 $\{X_{nj}:1\le j\le k_n\}_{n\ge1}$ 为独立随机变量构成的三角阵列，记 $\sigma_{n}^{2}=Var(\sum_{j=1}^{k_n}X_{nj})$。若存在常数使得 $|X_{nj}/\sigma_{n}|\le M_{nj}$ a.e.，且满足 $\lim_{n\rightarrow\infty} \max_{1\le j \le k_n} M_{nj}=0$，证明 $S_{n}=\sum_{j=1}^{k_{n}}X_{nj}$ 满足：

\[
\frac{S_{n}-E(S_{n})}{\sigma_{n}} \xrightarrow{d} \mathcal{N}(0,1)
\]

??? success "解答 (点击展开)"

    **证明：**
    令 $M_{n}=\max_{1\le j\le k_{n}}M_{nj}$。已知当 $n\rightarrow\infty$ 时 $M_{n}\rightarrow0$。由条件得：

    \[
    |X_{nj}|\le M_{nj}\sigma_{n}\le M_{n}\sigma_{n} \quad \text{a.e.}
    \]

    为应用 Liapounov 中心极限定理，我们验证 $\delta=1$ 时的 Liapounov 条件，即证明三阶绝对中心矩之和趋于 0：

    \[
    \lim_{n\rightarrow\infty} \frac{1}{\sigma_{n}^{3}}\sum_{j=1}^{k_{n}}E|X_{nj}-E[X_{nj}]|^{3} = 0
    \]

    1.  **偏差估计**：由三角不等式，
        $|X_{nj}-E[X_{nj}]|\le|X_{nj}|+E|X_{nj}|\le M_{n}\sigma_{n}+M_{n}\sigma_{n}=2M_{n}\sigma_{n}$。

    2.  **三阶矩界限**：
        $E|X_{nj}-E[X_{nj}]|^{3}=E[|X_{nj}-E[X_{nj}]|\cdot|X_{nj}-E[X_{nj}]|^{2}] \le E[(2M_{n}\sigma_{n})\cdot|X_{nj}-E[X_{nj}]|^{2}] = 2M_{n}\sigma_{n}Var(X_{nj})$。

    3.  **求和化简**：
        由于三角阵列每行独立，$\sum_{j=1}^{k_{n}}Var(X_{nj})=\sigma_{n}^{2}$。则：

        \[
        \frac{1}{\sigma_{n}^{3}}\sum_{j=1}^{k_{n}}E|X_{nj}-E[X_{nj}]|^{3} \le \frac{2M_{n}\sigma_{n}}{\sigma_{n}^{3}} \sum_{j=1}^{k_{n}}Var(X_{nj}) = \frac{2M_{n}\sigma_{n}}{\sigma_{n}^{3}} \cdot \sigma_{n}^{2} = 2M_{n}
        \]

    由于 $M_{n}\rightarrow0$，Liapounov 条件满足。根据 Liapounov CLT，结论成立。

---

## 习题 3: Null array

**题目：**
设 $\{X_{nj}:1\le j\le k_{n}\}_{n\ge1}$ 为三角阵列，$E(X_{nj})=\alpha_{nj}$。证明如下条件的蕴含关系：$(d)\Rightarrow(c)\Rightarrow(b)\Rightarrow(a)$：

(a) $\lim_{n\rightarrow\infty}P(|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})=0$ 对各 $j$ 成立。

(b) $\lim_{n\rightarrow\infty} \max_{j} P(|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})=0$。

(c) $\lim_{n\rightarrow\infty}P(\max_{j} |X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})=0$。

(d) $\lim_{n\rightarrow\infty} \sum_{j=1}^{k_{n}}P(|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})=0$。

??? success "解答 (点击展开)"

    **证明过程：**

    * **$(d)\Rightarrow(c)$**：
        最大值超过阈值的事件是各分量超过阈值事件的并集。由 **Boole 不等式 (概率次可加性)**：

        \[
        P(\max_{1\le j\le k_{n}}|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n}) = P \left( \bigcup_{j=1}^{k_{n}} \{|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n}\} \right) \le \sum_{j=1}^{k_{n}} P(|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})
        \]

        若 (d) 成立，右端趋于 0，则 (c) 必然成立。

    * **$(c)\Rightarrow(b)$**：
        对于任意固定的 $j_{0}$，有：
        $\{|X_{nj_{0}}-\alpha_{nj_{0}}|>\epsilon\sigma_{n}\} \subseteq \{\max_{1\le j\le k_{n}}|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n}\}$。
        由概率的单调性：

        \[
        P(|X_{nj_{0}}-\alpha_{nj_{0}}|>\epsilon\sigma_{n}) \le P(\max_{1\le j\le k_{n}}|X_{nj}-\alpha_{nj}|>\epsilon\sigma_{n})
        \]

        对左侧取最大值，不等式依然保持。故 (c) 成立蕴含 (b) 成立。

    * **$(b)\Rightarrow(a)$**：
        此结论显而易见。各分量的概率被其最大值所控制，最大值趋于 0 则每个分量必趋于 0。

---

## 习题 4: Sufficiency of Lindeberg's condition

**题目：**
设 $X_{n}\sim\mathcal{N}(0,2^{-n})$ 为独立随机变量序列。证明该序列满足中心极限定理，但 **Lindeberg 条件失效**。

??? success "解答 (点击展开)"

    **证明：**

    **第一步：计算方差与极限分布**
    
    独立随机变量和的方差为方差之和：

    \[
    \sigma_{n}^{2} = Var \left( \sum_{j=1}^{n}X_{j} \right) = \sum_{j=1}^{n}2^{-j} = \frac{1/2(1-(1/2)^{n})}{1-1/2} = 1-(1/2)^{n}
    \]

    当 $n\rightarrow\infty$ 时，$\sigma_{n}^{2}\rightarrow1$。计算标准化和 $T_{n} = S_{n}/\sigma_{n}$ 的特征函数：

    \[
    \phi_{T_{n}}(t) = \prod_{j=1}^{n}\phi_{X_{j}}(t/\sigma_{n}) = \prod_{j=1}^{n} \exp \left( -\frac{1}{2}\cdot2^{-j}\cdot \left( \frac{t}{\sigma_{n}} \right) ^{2} \right) = \exp \left( -\frac{t^{2}}{2\sigma_{n}^{2}} \sum_{j=1}^{n}2^{-j} \right) = \exp \left( -\frac{t^{2}}{2} \right)
    \]

    这正是 $\mathcal{N}(0,1)$ 的特征函数。根据 Lévy 连续性定理，$S_{n}/\sigma_{n} \xrightarrow{d} \mathcal{N}(0,1)$，即满足中心极限定理。

    **第二步：验证 Lindeberg 条件的失效**
    
    Lindeberg 量定义为 $L_{n}(\epsilon) = \frac{1}{\sigma_{n}^{2}}\sum_{j=1}^{n}E[X_{j}^{2}I(|X_{j}|>\epsilon\sigma_{n})]$。对于任意 $n\ge1$：

    \[
    L_{n}(\epsilon) \ge \frac{1}{\sigma_{n}^{2}}E[X_{1}^{2}I(|X_{1}|>\epsilon\sigma_{n})]
    \]

    当 $n\rightarrow\infty$ 时，$\sigma_{n}\rightarrow1$。因此：

    \[
    \liminf_{n\rightarrow\infty} L_{n}(\epsilon) \ge E[X_{1}^{2}I(|X_{1}|>\epsilon)]
    \]

    由于 $X_{1}\sim\mathcal{N}(0,1/2)$，对于任何固定的 $\epsilon>0$，上述期望是一个正的常数 $C>0$。
    因此 $L_{n}(\epsilon)$ 不会收敛到 0，Lindeberg 条件不成立。这说明 Lindeberg 条件并非中心极限定理成立的必要条件。

## 习题 5: M-dependence

**题目：**
设 $X_{1}, X_{2}, \dots$ 是独立同分布的随机变量，均值为 $\mu$，方差为 $\sigma^{2}$，且具有有限的四阶矩。

(a) 求 $\overline{X}_{n}$ 与 $\overline{Z}_{n} = (1/n)\sum_{i=1}^{n} X_{i}X_{i+1}$ 的联合渐近分布。

(b) 获取样本自相关系数 $r_{n} = \frac{\overline{Z}_{n} - \overline{X}_{n}^{2}}{(1/n)\sum X_{i}^{2} - \overline{X}_{n}^{2}}$ 的渐近分布。

??? success "解答 (点击展开)"

    **(a) 联合渐近分布的证明**

    * **构造线性组合**：
        定义 $Y_{i} = aX_{i} + bX_{i}X_{i+1}$。由于 $X_{i}$ 是 i.i.d. 序列，因此 $\{Y_{i}\}$ 是一个 $m=1$ 的相依序列。
        令 $W_{n} = a\overline{X}_{n} + b\overline{Z}_{n} = \frac{1}{n}\sum_{i=1}^{n} Y_{i}$。根据 m-相依序列的中心极限定理：

    \[
    \sqrt{n}(W_{n} - E[Y_{1}]) \xrightarrow{d} \mathcal{N}(0, V)
    \]

    \[
    V = Var(Y_{1}) + 2Cov(Y_{1}, Y_{2})
    \]

    * **计算方差 $Var(Y_{1})$**：
        引入中心化变量 $q_{i} = X_{i} - \mu$。则 $E[q_{i}] = 0$，$E[q_{i}^{2}] = \sigma^{2}$。
        计算 $Y_{1}$ 的偏差：

    \[
    Y_{1} - E[Y_{1}] = a(X_{1} - \mu) + b(X_{1}X_{2} - \mu^{2})
    \]

    \[
    = aq_{1} + b[(q_{1} + \mu)(q_{2} + \mu) - \mu^{2}]
    \]

    \[
    = (a + b\mu)q_{1} + b\mu q_{2} + bq_{1}q_{2}
    \]

    利用 $q_{i}$ 的独立性及各阶矩：

    \[
    Var(Y_{1}) = E[(a + b\mu)^{2}q_{1}^{2} + (b\mu)^{2}q_{2}^{2} + b^{2}q_{1}^{2}q_{2}^{2}]
    \]

    \[
    = (a + b\mu)^{2}\sigma^{2} + b^{2}\mu^{2}\sigma^{2} + b^{2}\sigma^{4}
    \]

    * **计算协方差 $Cov(Y_{1}, Y_{2})$**：
        同理可得 $Y_{2} - E[Y_{2}] = (a + b\mu)q_{2} + b\mu q_{3} + bq_{2}q_{3}$。
        计算乘积期望，仅包含 $q_{2}^{2}$ 的项非零：

    \[
    Cov(Y_{1}, Y_{2}) = E[(b\mu q_{2}) \cdot ((a + b\mu)q_{2})] = ab\mu\sigma^{2} + b^{2}\mu^{2}\sigma^{2}
    \]

    * **合并计算 $V$**：

    \[
    V = (a^{2} + 2ab\mu + b^{2}\mu^{2})\sigma^{2} + b^{2}\mu^{2}\sigma^{2} + b^{2}\sigma^{4} + 2ab\mu\sigma^{2} + 2b^{2}\mu^{2}\sigma^{2}
    \]

    \[
    = a^{2}\sigma^{2} + 4ab\mu\sigma^{2} + b^{2}(4\mu^2\sigma^2 + \sigma^4)
    \]

    * **结论**：
        根据 Cramer-Wold 设备，联合分布对应的协方差矩阵 $\Sigma$ 为：

    \[
    \sqrt{n} \begin{pmatrix} \overline{X}_{n} - \mu \\ \overline{Z}_{n} - \mu^{2} \end{pmatrix} \xrightarrow{d} \mathcal{N} \left( 0, \begin{pmatrix} \sigma^{2} & 2\mu\sigma^{2} \\ 2\mu\sigma^{2} & 4\mu^{2}\sigma^{2} + \sigma^{4} \end{pmatrix} \right)
    \]

    **(b) $r_{n}$ 的渐近分布证明**

    * **分母分析**：
        由大数定律，$\frac{1}{n}\sum X_{i}^{2} \xrightarrow{p} E[X^{2}] = \mu^{2} + \sigma^{2}$ 且 $\overline{X}_{n} \xrightarrow{p} \mu$。
        故分母 $\frac{1}{n}\sum X_{i}^{2} - \overline{X}_{n}^{2} \xrightarrow{p} \sigma^{2}$。

    * **分子分析 (一阶 Delta 方法)**：
        令 $N_{n} = \overline{Z}_{n} - \overline{X}_{n}^{2}$。定义函数 $f(A, C) = C - A^{2}$。
        在 $(\mu, \mu^{2})$ 处的梯度为：

    \[
    \nabla f(\mu, \mu^{2}) = \begin{pmatrix} -2\mu \\ 1 \end{pmatrix}
    \]

    分子的渐近方差为 $\nabla f^{T} \Sigma \nabla f$：

    \[
    \begin{pmatrix} -2\mu & 1 \end{pmatrix} \begin{pmatrix} \sigma^{2} & 2\mu\sigma^{2} \\ 2\mu\sigma^{2} & 4\mu^{2}\sigma^{2} + \sigma^{4} \end{pmatrix} \begin{pmatrix} -2\mu \\ 1 \end{pmatrix} = \sigma^{4}
    \]

    * **最终结果**：
        结合 Slutsky 定理：

    \[
    \sqrt{n}r_{n} = \frac{1}{\sigma^{2} + o_{p}(1)} \sqrt{n}(\overline{Z}_{n} - \overline{X}_{n}^{2}) \xrightarrow{d} \frac{1}{\sigma^{2}} \mathcal{N}(0, \sigma^{4}) \sim \mathcal{N}(0, 1)
    \]

---

## 习题 6: First-Order Delta Method

**题目：**
(a) 寻找 $(\sqrt{n}(\overline{X} - \mu), \sqrt{n}(S^{2} - \sigma^{2}))$ 的联合极限分布，并讨论渐近独立的条件。

(b) 设 $X_{1}, \dots, X_{n}$ 是均值为 $\theta$ 的 Poisson 样本，寻找其样本均值的稳方差变换并构造置信区间。

??? success "解答 (点击展开)"

    **(a) 联合分布与独立性**

    * **多元中心极限定理**：
        令 $M_{n}^{(2)} = \frac{1}{n}\sum_{i=1}^{n}(X_{i} - \mu)^{2}$。其渐近协方差矩阵 $\Sigma^{*}$ 元素为：

    \[
    \Sigma_{11}^{*} = Var(X_{1}) = \sigma^{2}
    \]

    \[
    \Sigma_{22}^{*} = Var((X_{1} - \mu)^{2}) = E[(X_{1} - \mu)^{4}] - \sigma^{4}
    \]

    \[
    \Sigma_{12}^{*} = Cov(X_{1} - \mu, (X_{1} - \mu)^{2}) = E[(X_{1} - \mu)^{3}]
    \]

    * **渐近等价性**：
        由 $S_{n}^{2} = \frac{n}{n-1} [M_{n}^{(2)} - (\overline{X}_{n} - \mu)^{2}]$，项 $\sqrt{n}(\overline{X}_{n} - \mu)^{2} = o_{p}(1)$。
        因此 $\sqrt{n}(S_{n}^{2} - \sigma^{2})$ 与 $\sqrt{n}(M_{n}^{(2)} - \sigma^{2})$ 渐近等价。

    * **结论**：
        联合分布为：

    \[
    \sqrt{n} \begin{pmatrix} \overline{X}_{n} - \mu \\ S_{n}^{2} - \sigma^{2} \end{pmatrix} \xrightarrow{d} \mathcal{N} \left( 0, \begin{pmatrix} \sigma^{2} & E[(X_{1} - \mu)^{3}] \\ E[(X_{1} - \mu)^{3}] & E[(X_{1} - \mu)^{4}] - \sigma^{4} \end{pmatrix} \right)
    \]

    *渐近独立条件*：协方差项 $E[(X_{1} - \mu)^{3}] = 0$，即分布的三阶中心矩（偏度）为零（如对称分布）。

    **(b) 稳方差变换与置信区间**

    * **稳方差变换推导**：
        由 CLT，$\sqrt{n}(\overline{X}_{n} - \theta) \xrightarrow{d} \mathcal{N}(0, \theta)$。
        设变换为 $g(x)$，由 Delta 方法，渐近方差为 $[g'(\theta)]^{2}\theta$。
        令 $[g'(\theta)]^{2}\theta = 1$，则 $g'(\theta) = \theta^{-1/2}$。
        积分得 $g(\theta) = 2\sqrt{\theta}$。

    * **变换后的性质**：

    \[
    \sqrt{n}(2\sqrt{\overline{X}_{n}} - 2\sqrt{\theta}) \xrightarrow{d} \mathcal{N}(0, 1)
    \]

    * **构造置信区间**：
        对于 $(1-\alpha)$ 置信水平：

    \[
    P \left( -z_{\alpha/2} \le \sqrt{n}(2\sqrt{\overline{X}_{n}} - 2\sqrt{\theta}) \le z_{\alpha/2} \right) \approx 1 - \alpha
    \]

    整理得：

    \[
    \sqrt{\overline{X}_{n}} - \frac{z_{\alpha/2}}{2\sqrt{n}} \le \sqrt{\theta} \le \sqrt{\overline{X}_{n}} + \frac{z_{\alpha/2}}{2\sqrt{n}}
    \]

    两边平方得到 $\theta$ 的置信区间：

    \[
    \theta \in \left[ \left( \sqrt{\overline{X}_{n}} - \frac{z_{\alpha/2}}{2\sqrt{n}} \right)^{2}, \left( \sqrt{\overline{X}_{n}} + \frac{z_{\alpha/2}}{2\sqrt{n}} \right)^{2} \right]
    \]