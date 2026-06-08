# 📈 课后习题精解：大样本理论 Assignment 4

## 1. 单步估计量 (One-step estimator)

**题目：**
设 $\{P_{\theta}\}_{\theta\in\Theta}$ 是一个模型族，其中 $\Theta\subset\mathbb{R}^{d}$ 是开集，且 $X_{1},...,X_{n} \stackrel{i.i.d}{\sim} P_{\theta_{0}}$。对数似然函数 $l_{\theta}(x)=\log p_{\theta}(x)$ 满足给定的正则性条件。令 $\hat{\theta}_{n}$ 为一个 $n$-相合 (n-consistent) 估计量。单步估计量 $\delta_{n}$ 满足以下等式：

\[
\nabla L_{n}(\hat{\theta}_{n})+\nabla^{2}L_{n}(\hat{\theta}_{n})(\delta_{n}-\hat{\theta}_{n})=0
\]

(a) 假设模型族 $\{P_{\theta}\}_{\theta\in\mathbb{R}}$ 是柯西 (Cauchy) 族，其概率密度函数为 $p_{\theta}(x)= \frac{1}{\pi(1+(x-\theta)^{2})}$。令 $X_{1},...,X_{n} \stackrel{i.i.d}{\sim} P_{\theta}$ 并定义 $\hat{\theta}_{n}=\text{Median}(X_{1},...,X_{n})$。求 $\sqrt{n}(\hat{\theta}_{n}-\theta)$ 的极限分布。

(b) 令 $\delta_{n}$ 为该柯西族的单步估计量。它的渐近分布是什么？在 (a) 中，$\delta_{n}$ 是否比 $\hat{\theta}_{n}$ 更有效 (more efficient)？

??? success "解答 (点击展开)"

    **(a) 求样本中位数的渐近分布**

    样本中位数 $\hat{\theta}_{n}$ 对应于 $p=1/2$ 处的样本分位数。根据已建立的样本分位数渐近理论（例如，van der Vaart 中的定理 5.23），如果累积分布函数 $F$ 在真实总体分位数 $\theta$ 处具有连续且严格为正的导数（密度）$f$，则样本中位数是渐近正态的：

    \[
    \sqrt{n}(\hat{\theta}_{n}-\theta)\xrightarrow{d}\mathcal{N}\left(0,\frac{p(1-p)}{f(\theta)^{2}}\right)
    \]

    对于 $p=1/2,$ 方差的分子为 $1/4.$ 我们现在计算柯西密度函数在真实中位数 $\theta$ 处的值：

    \[
    f(\theta)=p_{\theta}(\theta)=\frac{1}{\pi(1+(\theta-\theta)^{2})}=\frac{1}{\pi}
    \]

    将 $f(\theta)=1/\pi$ 代入渐近方差公式：

    \[
    V=\frac{1/4}{(1/\pi)^{2}}=\frac{\pi^{2}}{4}
    \]

    因此，样本中位数的极限分布为：

    \[
    \sqrt{n}(\hat{\theta}_{n}-\theta)\xrightarrow{d}\mathcal{N}\left(0,\frac{\pi^{2}}{4}\right)
    \]

    **(b) 求单步估计量的渐近分布并比较有效性**

    根据单步估计量的定义，我们可以将其重写为：

    \[
    \sqrt{n}(\delta_{n}-\theta)=\sqrt{n}(\hat{\theta}_{n}-\theta)-\left[\frac{1}{n}\nabla^{2}L_{n}(\hat{\theta}_{n})\right]^{-1}\frac{1}{\sqrt{n}}\nabla L_{n}(\hat{\theta}_{n})
    \]

    我们利用泰勒定理 (Taylor's theorem) 将得分函数 $\frac{1}{\sqrt{n}}\nabla L_{n}(\hat{\theta}_{n})$ 在真实参数 $\theta$ 附近展开：

    \[
    \frac{1}{\sqrt{n}}\nabla L_{n}(\hat{\theta}_{n})=\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)+\frac{1}{n}\nabla^{2}L_{n}(\tilde{\theta})\sqrt{n}(\hat{\theta}_{n}-\theta)
    \]

    其中 $\tilde{\theta}$ 位于 $\hat{\theta}_{n}$ 和 $\theta$ 之间。因为 $\hat{\theta}_{n}$ 是 $\sqrt{n}$-相合的，$\tilde{\theta}\xrightarrow{P}\theta$。由弱大数定律 (Weak Law of Large Numbers) 和二阶导数的 Lipschitz 连续性，我们有：

    \[
    -\frac{1}{n}\nabla^{2}L_{n}(\hat{\theta}_{n})\xrightarrow{p}I(\theta) \quad \text{和} \quad -\frac{1}{n}\nabla^{2}L_{n}(\tilde{\theta})\xrightarrow{p}I(\theta)
    \]

    其中 $I(\theta)=-E[\nabla^{2}l_{\theta}(X)]$ 是 Fisher 信息量。将泰勒展开代回我们的估计量，得到：

    \[
    \sqrt{n}(\delta_{n}-\theta)=\sqrt{n}(\hat{\theta}_{n}-\theta)+I(\theta)^{-1}\left(\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)-I(\theta)\sqrt{n}(\hat{\theta}_{n}-\theta)\right)+o_{p}(1)
    \]

    \[
    =\sqrt{n}(\hat{\theta}_{n}-\theta)+I(\theta)^{-1}\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)-\sqrt{n}(\hat{\theta}_{n}-\theta)+o_{p}(1)
    \]

    \[
    =I(\theta)^{-1}\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)+o_{p}(1)
    \]

    由中心极限定理，标准化后的得分函数 $\frac{1}{\sqrt{n}}\nabla L_{n}(\theta)\xrightarrow{d}\mathcal{N}(0,I(\theta))$。由 Slutsky 定理，单步估计量达到了极大似然估计量的有效性：

    \[
    \sqrt{n}(\delta_{n}-\theta)\xrightarrow{d}\mathcal{N}(0,I(\theta)^{-1})
    \]

    现在，我们计算柯西分布的 Fisher 信息量 $I(\theta)$。对数密度为：

    \[
    l_{\theta}(x)=-\log\pi-\log(1+(x-\theta)^{2})
    \]

    一阶导数（得分函数）为：

    \[
    \nabla l_{\theta}(x)=\frac{2(x-\theta)}{1+(x-\theta)^{2}}
    \]

    Fisher 信息量是得分函数平方的期望：

    \[
    I(\theta)=E[(\nabla l_{\theta}(X))^{2}]=\int_{-\infty}^{\infty}\left(\frac{2(x-\theta)}{1+(x-\theta)^{2}}\right)^{2}\frac{1}{\pi(1+(x-\theta)^{2})}dx
    \]

    令 $y=x-\theta$。由于被积函数是偶函数：

    \[
    I(\theta)=\frac{4}{\pi}\int_{-\infty}^{\infty}\frac{y^{2}}{(1+y^{2})^{3}}dy=\frac{8}{\pi}\int_{0}^{\infty}\frac{y^{2}}{(1+y^{2})^{3}}dy
    \]

    我们可以通过三角代换 $y=\tan z$ 来求解，这给出 $dy=\sec^{2}z dz$ 和 $1+y^{2}=\sec^{2}z$。积分区间从 $[0,\infty)$ 变为 $[0,\pi/2)$：

    \[
    I(\theta)=\frac{8}{\pi}\int_{0}^{\pi/2}\frac{\tan^{2}z}{\sec^{6}z}\sec^{2}z dz
    \]

    \[
    =\frac{8}{\pi}\int_{0}^{\pi/2}\sin^{2}z \cos^{2}z dz
    \]

    \[
    =\frac{8}{\pi}\int_{0}^{\pi/2}\frac{1}{4}\sin^{2}(2z)dz
    \]

    \[
    =\frac{2}{\pi}\int_{0}^{\pi/2}\sin^{2}(2z)dz
    \]

    令 $u=2z$，所以 $du=2dz$，积分区间变为 $[0, \pi]$：

    \[
    I(\theta)=\frac{1}{\pi}\int_{0}^{\pi}\sin^{2}u du=\frac{1}{\pi}\cdot\frac{\pi}{2}=\frac{1}{2}
    \]

    因此，单步估计量的渐近方差为 $I(\theta)^{-1}=2$。

    \[
    \sqrt{n}(\delta_{n}-\theta)\xrightarrow{d}\mathcal{N}(0,2)
    \]

    **有效性比较**：由 (a) 部分，样本中位数 $\hat{\theta}_{n}$ 的渐近方差为 $\frac{\pi^{2}}{4}\approx2.467$。因为 $2<\frac{\pi^{2}}{4}$，单步估计量 $\delta_{n}$ 的渐近方差严格小于初始中位数估计量 $\hat{\theta}_{n}$ 的渐近方差。因此，$\delta_{n}$ 更有效。

---

## 2. 无偏估计量与 U-统计量 (Unbiased estimator and U-statistic)

**题目：**
令 $X_{1},...,X_{n}$ 为 i.i.d. 随机变量，具有有限的 $\mu=\mathbb{E}(X_{1})$ 和 $\overline{\mu}=\mathbb{E}(X_{1}^{-1})$。
寻找一个 U-统计量，它是 $\overline{\mu}\mu$ 的无偏估计量，并推导其方差和渐近分布。

??? success "解答 (点击展开)"

    **步骤 1：构造对称核和 U-统计量**

    我们想估计 $\theta=\mu\overline{\mu}$。考虑一对独立的观测值 $(X_{1},X_{2})$。由于它们是独立的，我们有：

    \[
    \mathbb{E}[X_{1}X_{2}^{-1}]=\mathbb{E}[X_{1}]\mathbb{E}[X_{2}^{-1}]=\mu\overline{\mu}
    \]

    因此，$g(x_{1},x_{2})=x_{1}x_{2}^{-1}$ 是 $\theta$ 的一个无偏估计量，但它是不对称的。我们通过排列参数构造一个次数 $m=2$ 的对称核：

    \[
    h(X_{1},X_{2})=\frac{1}{2}(g(X_{1},X_{2})+g(X_{2},X_{1}))=\frac{1}{2}\left(\frac{X_{1}}{X_{2}}+\frac{X_{2}}{X_{1}}\right)
    \]

    显然，$\mathbb{E}[h(X_{1},X_{2})]=\frac{1}{2}(\mu\overline{\mu}+\mu\overline{\mu})=\mu\overline{\mu}$。相应的 U-统计量为：

    \[
    U_{n}=\binom{n}{2}^{-1}\sum_{1\le i<j\le n}h(X_{i},X_{j})=\binom{n}{2}^{-1}\sum_{1\le i<j\le n}\frac{1}{2}\left(\frac{X_{i}}{X_{j}}+\frac{X_{j}}{X_{i}}\right)
    \]

    **步骤 2：推导 $U_{n}$ 的方差**

    根据 Hoeffding 定理，次数 $m=2$ 的 U-统计量的方差由下式给出：

    \[
    Var(U_{n})=\binom{n}{2}^{-1}\sum_{k=1}^{2}\binom{2}{k}\binom{n-2}{2-k}\zeta_{k}=\frac{2}{n(n-1)}(2(n-2)\zeta_{1}+\zeta_{2})
    \]

    首先，我们计算 $h_{1}(x_{1})=\mathbb{E}[h(X_{1},X_{2})|X_{1}=x_{1}]$：

    \[
    h_{1}(x_{1})=\mathbb{E}\left[\frac{1}{2}\left(\frac{x_{1}}{X_{2}}+\frac{X_{2}}{x_{1}}\right)\right]=\frac{1}{2}\overline{\mu}x_{1}+\frac{1}{2}\frac{\mu}{x_{1}}
    \]

    现在我们计算 $\zeta_{1}=Var(h_{1}(X_{1}))$。令 $\sigma^{2}=Var(X_{1})$ 和 $\overline{\sigma}^{2}=Var(X_{1}^{-1})$。利用性质 $Var(aX+bY)=a^{2}Var(X)+b^{2}Var(Y)+2abCov(X,Y)$：

    \[
    \zeta_{1}=Var\left(\frac{1}{2}\overline{\mu}X_{1}+\frac{1}{2}\mu X_{1}^{-1}\right)
    \]

    \[
    =\frac{1}{4}\overline{\mu}^{2}Var(X_{1})+\frac{1}{4}\mu^{2}Var(X_{1}^{-1})+2\left(\frac{1}{2}\overline{\mu}\right)\left(\frac{1}{2}\mu\right)Cov(X_{1},X_{1}^{-1})
    \]

    \[
    =\frac{1}{4}(\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(\mathbb{E}[X_{1}X_{1}^{-1}]-\mathbb{E}[X_{1}]\mathbb{E}[X_{1}^{-1}]))
    \]

    \[
    =\frac{1}{4}(\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(1-\mu\overline{\mu}))
    \]

    接下来，我们计算 $\zeta_{2}=Var(h(X_{1},X_{2}))=\mathbb{E}[h^{2}(X_{1},X_{2})]-(\mu\overline{\mu})^{2}$：

    \[
    \mathbb{E}[h^{2}(X_{1},X_{2})]=\mathbb{E}\left[\frac{1}{4}\left(\frac{X_{1}^{2}}{X_{2}^{2}}+\frac{X_{2}^{2}}{X_{1}^{2}}+2\frac{X_{1}X_{2}}{X_{2}X_{1}}\right)\right]
    \]

    \[
    =\frac{1}{4}(\mathbb{E}[X_{1}^{2}]\mathbb{E}[X_{2}^{-2}]+\mathbb{E}[X_{2}^{2}]\mathbb{E}[X_{1}^{-2}]+2)
    \]

    \[
    =\frac{1}{2}\mathbb{E}[X_{1}^{2}]\mathbb{E}[X_{1}^{-2}]+\frac{1}{2}
    \]

    回顾二阶矩和方差之间的关系：$\mathbb{E}[X_{1}^{2}]=\sigma^{2}+\mu^{2}$ 和 $\mathbb{E}[X_{1}^{-2}]=\overline{\sigma}^{2}+\overline{\mu}^{2}$，将它们代入方程：

    \[
    \zeta_{2}=\frac{1}{2}(\sigma^{2}+\mu^{2})(\overline{\sigma}^{2}+\overline{\mu}^{2})+\frac{1}{2}-\mu^{2}\overline{\mu}^{2}
    \]

    \[
    =\frac{1}{2}(\sigma^{2}\overline{\sigma}^{2}+\sigma^{2}\overline{\mu}^{2}+\mu^{2}\overline{\sigma}^{2}+\mu^{2}\overline{\mu}^{2})+\frac{1}{2}-\mu^{2}\overline{\mu}^{2}
    \]

    \[
    =\frac{1}{2}(\sigma^{2}\overline{\sigma}^{2}+\sigma^{2}\overline{\mu}^{2}+\mu^{2}\overline{\sigma}^{2}-\mu^{2}\overline{\mu}^{2}+1)
    \]

    将 $\zeta_{1}$ 和 $\zeta_{2}$ 代回公式即可得到精确方差：

    \[
    Var(U_{n})=\frac{4(n-2)}{n(n-1)}\zeta_{1}+\frac{2}{n(n-1)}\zeta_{2}
    \]

    \[
    =\frac{n-2}{n(n-1)}(\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(1-\mu\overline{\mu}))+\frac{1}{n(n-1)}(\sigma^{2}\overline{\sigma}^{2}+\sigma^{2}\overline{\mu}^{2}+\mu^{2}\overline{\sigma}^{2}-\mu^{2}\overline{\mu}^{2}+1)
    \]

    **步骤 3：推导渐近分布**

    根据 U-统计量的投影定理（Hoeffding's CLT），由于 $\mathbb{E}[h^{2}]<\infty$ 且 $\zeta_{1}>0$，U-统计量是渐近正态的，由一阶投影驱动：

    \[
    \sqrt{n}(U_{n}-\mu\overline{\mu})\rightarrow\mathcal{N}(0,m^{2}\zeta_{1})
    \]

    其中核的次数为 $m=2$。因此，渐近方差为 $4\zeta_{1}$：

    \[
    4\zeta_{1}=\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(1-\mu\overline{\mu})
    \]

    因此，渐近分布为：

    \[
    \sqrt{n}(U_{n}-\mu\overline{\mu})\rightarrow\mathcal{N}(0,\overline{\mu}^{2}\sigma^{2}+\mu^{2}\overline{\sigma}^{2}+2\mu\overline{\mu}(1-\mu\overline{\mu}))
    \]

---

## 3. Hoeffding 定理 (Hoeffding's theorem)

**题目：**
令 $X_{1},...,X_{n} \stackrel{i.i.d}{\sim} P\in\mathcal{P}$ 且 $\mathbb{E}_{P}h^{2}(X_{1},...,X_{m})<\infty$，则：

\[
Var_{P}(U_{n})=\binom{n}{m}^{-1}\sum_{k=1}^{m}\binom{m}{k}\binom{n-m}{m-k}\xi_{k}
\]

其中 $\xi_{k}=Var_{P}(h_{k}(X_{1},...,X_{k}))$。

(a) 证明 $\frac{m^{2}}{n}\xi_{1}\le Var_{P}(U_{n})\le\frac{m}{n}\xi_{m}$

(b) 证明 $(n+1)Var_{P}(U_{n+1})\le nVar_{P}(U_{n})$

??? success "解答 (点击展开)"

    **(a) 证明上下界**

    根据 U-统计量的 Hoeffding H-分解 (ANOVA decomposition)，我们可以将 $U_{n}$ 的方差重写为投影方差 $\eta_{j}$（第 $j$ 阶正交投影核的方差）的形式：

    \[
    Var_{P}(U_{n})=\sum_{j=1}^{m}\binom{m}{j}^{2}\binom{n}{j}^{-1}\eta_{j}
    \]

    根据定义，对于所有 $j=1,...,m$，$\eta_{j}\ge0$。$\xi_{k}$ 和 $\eta_{j}$ 之间的关系由下式给出：

    \[
    \xi_{k}=\sum_{j=1}^{k}\binom{k}{j}\eta_{j}
    \]

    特别地，对于 $k=1$，我们有 $\xi_{1}=\binom{1}{1}\eta_{1}=\eta_{1}$。对于 $k=m$，我们有 $\xi_{m}=\sum_{j=1}^{m}\binom{m}{j}\eta_{j}$。

    **下界：** 由于方差和中的所有项都是非负的 $(\eta_{j}\ge0)$，总方差有以下界，由它的第一项 $(j=1)$ 给出：

    \[
    Var_{P}(U_{n})\ge\binom{m}{1}^{2}\binom{n}{1}^{-1}\eta_{1}=\frac{m^{2}}{n}\eta_{1}=\frac{m^{2}}{n}\xi_{1}
    \]

    **上界：** 我们分析方差和中的组合系数。我们可以将其重写为：

    \[
    \binom{m}{j}^{2}\binom{n}{j}^{-1}=\binom{m}{j}\left[\binom{m}{j}\binom{n}{j}^{-1}\right]
    \]

    展开方括号中的项：

    \[
    \binom{m}{j}\binom{n}{j}^{-1}=\frac{m!}{j!(m-j)!}\frac{j!(n-j)!}{n!}=\frac{m(m-1)...(m-j+1)}{n(n-1)...(n-j+1)}=\prod_{i=0}^{j-1}\frac{m-i}{n-i}
    \]

    由于 $m\le n$，对于任意 $i\ge0$，我们有 $nn-in\le mn-im$，这蕴含着 $\frac{m-i}{n-i}\le\frac{m}{n}$。因此，该乘积的边界为：

    \[
    \prod_{i=0}^{j-1}\frac{m-i}{n-i}\le\frac{m}{n}\cdot1...1=\frac{m}{n}
    \]

    将此不等式代回方差公式，我们得到：

    \[
    Var_{P}(U_{n})=\sum_{j=1}^{m}\binom{m}{j}\left[\binom{m}{j}\binom{n}{j}^{-1}\right]\eta_{j}\le\sum_{j=1}^{m}\binom{m}{j}\left(\frac{m}{n}\right)\eta_{j}=\frac{m}{n}\sum_{j=1}^{m}\binom{m}{j}\eta_{j}
    \]

    由于 $\xi_{m}=\sum_{j=1}^{m}\binom{m}{j}\eta_{j}$，我们得出结论：

    \[
    Var_{P}(U_{n})\le\frac{m}{n}\xi_{m}
    \]

    **(b) 证明样本量递推不等式**

    使用 H-分解方差公式，我们写出样本大小 $n$ 和 $n+1$ 的按比例缩放的方差：

    \[
    nVar_{P}(U_{n})=\sum_{j=1}^{m}\binom{m}{j}^{2}\left[n\binom{n}{j}^{-1}\right]\eta_{j}
    \]

    \[
    (n+1)Var_{P}(U_{n+1})=\sum_{j=1}^{m}\binom{m}{j}^{2}\left[(n+1)\binom{n+1}{j}^{-1}\right]\eta_{j}
    \]

    为了证明该不等式，只需证明对于每个 $\eta_{j}$ 的系数，当样本大小从 $n$ 增加到 $n+1$ 时，它是不增的。我们比较方括号中的项。对于大小 $n+1$：

    \[
    (n+1)\binom{n+1}{j}^{-1}=(n+1)\frac{j!(n+1-j)!}{(n+1)!}=\frac{j!(n-j+1)!}{n!}
    \]

    我们可以从分子中提取出 $(n-j+1)$：

    \[
    \frac{j!(n-j+1)!}{n!}=(n-j+1)\frac{j!(n-j)!}{n!}=(n-j+1)\binom{n}{j}^{-1}
    \]

    由于对于和中的所有项 $j\ge1$，我们显然有 $n-j+1\le n$。因此：

    \[
    (n+1)\binom{n+1}{j}^{-1}=(n-j+1)\binom{n}{j}^{-1}\le n\binom{n}{j}^{-1}
    \]

    因为对所有 $j$，$\binom{m}{j}^{2}\eta_{j}\ge0$，将这个非负常数乘到不等式两边保持方向不变。对所有 $j=1,...,m$ 求和得出：

    \[
    \sum_{j=1}^{m}\binom{m}{j}^{2}\left[(n+1)\binom{n+1}{j}^{-1}\right]\eta_{j}\le\sum_{j=1}^{m}\binom{m}{j}^{2}\left[n\binom{n}{j}^{-1}\right]\eta_{j}
    \]

    这正是：

    \[
    (n+1)Var_{P}(U_{n+1})\le nVar_{P}(U_{n})
    \]

## 习题 4: U-统计量与 V-统计量 (U-statistics and V-statistics)

**题目：**
设 $X_1, X_2, \dots$ 为独立同分布的随机变量序列，且 $h$ 为对称核函数，满足 $0 < E[h^2(X_1, X_2)] < \infty$ 及 $E[h^2(X_1, X_1)] < \infty$。定义 $V_n$ 为：

\[
V_n := \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n h(X_i, X_j)
\]

(a) 利用 U-统计量的中心极限定理证明 $V_n$ 是渐近正态的。

(b) $V_n$ 能否写成一个 U-统计量？为什么？

(c) 现假设 $E X_1 = 0$ 且 $E X_1^4 < \infty$。取核函数 $h(x, y) = xy$，求 $V_n$ 相对于 U-统计量 $U_n := \frac{2}{n(n-1)} \sum_{1 \le i < j \le n} X_i X_j$ 的渐近相对效率 (ARE)，即 $\frac{Avar(U_n)}{Avar(V_n)}$。

??? success "解答 (点击展开)"

    **(a) 证明 $V_n$ 的渐近正态性**

    令 $\theta = E[h(X_1, X_2)]$ 为核函数对于独立观测值的期望。对应的 U-统计量为：

    \[
    U_n = \frac{2}{n(n-1)} \sum_{1 \le i < j \le n} h(X_i, X_j)
    \]

    我们可以将 V-统计量分解为 U-统计量部分（非对角线项）和对角线部分：

    \[
    V_n = \frac{1}{n^2} \sum_{i \ne j} h(X_i, X_j) + \frac{1}{n^2} \sum_{i=1}^n h(X_i, X_i)
    \]

    \[
    = \frac{n(n-1)}{n^2} U_n + \frac{1}{n^2} \sum_{i=1}^n h(X_i, X_i)
    \]

    \[
    = \left( 1 - \frac{1}{n} \right) U_n + \frac{1}{n^2} \sum_{i=1}^n h(X_i, X_i)
    \]

    为了分析其渐近分布，我们将中心化的 V-统计量乘以 $\sqrt{n}$ 进行缩放：

    \[
    \sqrt{n}(V_n - \theta) = \sqrt{n} \left( \left( 1 - \frac{1}{n} \right) U_n + \frac{1}{n^2} \sum_{i=1}^n h(X_i, X_i) - \theta \right)
    \]

    \[
    = \left( 1 - \frac{1}{n} \right) \sqrt{n}(U_n - \theta) - \frac{1}{\sqrt{n}} \theta + \frac{1}{n\sqrt{n}} \sum_{i=1}^n h(X_i, X_i)
    \]

    令 $\mu_D = E[h(X_1, X_1)]$。由弱大数定律 (WLLN) 可知，$\frac{1}{n} \sum_{i=1}^n h(X_i, X_i) \xrightarrow{p} \mu_D$。
    因此，上式中的最后一项为 $\frac{1}{n} (\frac{1}{n} \sum h(X_i, X_i)) = O_p(n^{-1})$，它是 $o_p(1)$ 的。
    同时，缩放因子 $(1 - 1/n) \rightarrow 1$。
    
    由 Slutsky 定理，V-统计量的缩放差趋于 U-统计量的缩放差：

    \[
    \sqrt{n}(V_n - \theta) = \sqrt{n}(U_n - \theta) + o_p(1)
    \]

    根据 U-统计量的中心极限定理，$\sqrt{n}(U_n - \theta) \xrightarrow{d} \mathcal{N}(0, 4\zeta_1)$。因此 $V_n$ 具有相同的渐近正态分布：

    \[
    \sqrt{n}(V_n - \theta) \xrightarrow{d} \mathcal{N}(0, 4\zeta_1)
    \]

    **(b) $V_n$ 能否写成 U-统计量**

    不能。根据定义，大小为 $n$ 的 U-统计量 $U_n$ 必须是其总体参数 $\theta$ 的无偏估计量，即对所有的 $n \ge m$ 满足 $E[U_n] = \theta$。
    计算 $V_n$ 的期望：

    \[
    E[V_n] = \left( 1 - \frac{1}{n} \right) \theta + \frac{1}{n} \mu_D = \theta + \frac{\mu_D - \theta}{n}
    \]

    由于 $E[V_n]$ 显式依赖于样本量 $n$（除非 $\theta = \mu_D$，这在一般情况下并不成立），因此 $V_n$ 是一个有偏估计量，不可能成为 U-统计量。

    **(c) 计算渐近相对效率 (ARE)**

    由于 $E X_1 = 0$ 且 $h(x, y) = xy$，我们计算 $\theta$ 和投影方差 $\zeta_1$：

    \[
    \theta = E[X_1 X_2] = E[X_1] E[X_2] = 0 \cdot 0 = 0
    \]

    投影函数为：

    \[
    h_1(x_1) = E[h(x_1, X_2)] = x_1 E[X_2] = 0
    \]

    因为 $h_1(x_1) \equiv 0$，所以 $\zeta_1 = Var(h_1(X_1)) = 0$。这意味着 $U_n$ 是一个退化的 (degenerate) U-统计量，标准 $\sqrt{n}$ 缩放后的极限为 0。
    我们需要使用 $n$ 的缩放因子。

    **分析 $n V_n$**：
    
    \[
    V_n = \frac{1}{n^2} \sum_{i=1}^n \sum_{j=1}^n X_i X_j = \left( \frac{1}{n} \sum_{i=1}^n X_i \right)^2 = \bar{X}_n^2
    \]

    由 CLT，$\sqrt{n}\bar{X}_n \xrightarrow{d} \mathcal{N}(0, \sigma^2)$。由连续映射定理 (CMT)：

    \[
    n V_n = (\sqrt{n} \bar{X}_n)^2 \xrightarrow{d} \sigma^2 \chi_1^2
    \]

    其渐近方差为 $Var(\sigma^2 \chi_1^2) = \sigma^4 Var(\chi_1^2) = \sigma^4 \cdot 2 = 2\sigma^4$。

    **分析 $n U_n$**：

    \[
    n U_n = \frac{n}{n(n-1)} \sum_{i \ne j} X_i X_j = \frac{n}{n-1} \left[ \left( \sum X_i \right)^2 - \sum X_i^2 \right] \cdot \frac{1}{n}
    \]

    \[
    = \frac{1}{n-1} \left[ n (\sqrt{n} \bar{X}_n)^2 - \sum X_i^2 \right] \xrightarrow{d} \sigma^2 \chi_1^2 - \sigma^2 = \sigma^2(\chi_1^2 - 1)
    \]

    其渐近方差为 $Var(\sigma^2(\chi_1^2 - 1)) = \sigma^4 Var(\chi_1^2) = 2\sigma^4$。

    因此，ARE 为：

    \[
    ARE(U_n, V_n) = \frac{Avar(n U_n)}{Avar(n V_n)} = \frac{2\sigma^4}{2\sigma^4} = 1
    \]

---

## 习题 5: 大数据下的分治 U-统计量 (Divide and Conquered U-statistics for big data)

**题目：**
设 $U_N$ 为核函数为 $h(x_1, \dots, x_m)$ 的 U-统计量。假设 $\theta = E[h(X_1, \dots, X_m)]$ 且 $E[|h|] < \infty$。定义 $h_1(x_1) = E[h(X_1, X_2, \dots, X_m) | X_1 = x_1]$ 并假设 $\zeta_1 = Var(h_1(X_1)) > 0$。假设总样本量 $N = nK$，数据被划分为 $K$ 个非重叠子集 $D_k$（每个大小为 $n$）。令 $U_{kn}$ 为基于第 $k$ 个子集 $D_k$ 的 U-统计量。定义 $AU_N = K^{-1} \sum_{k=1}^K U_{kn}$。证明：若 $K = o(N)$，则 $\sqrt{N}(AU_N - \theta) \xrightarrow{d} \mathcal{N}(0, m^2\zeta_1)$。

??? success "解答 (点击展开)"

    **证明过程：**

    我们要证明 $\sqrt{N}(AU_N - \theta) \xrightarrow{d} \mathcal{N}(0, m^2\zeta_1)$。我们将缩放后的平均 U-统计量重写为独立随机变量之和。
    因为 $N = nK$，我们有：

    \[
    \sqrt{N}(AU_N - \theta) = \sqrt{nK} \left( \frac{1}{K} \sum_{k=1}^K (U_{kn} - \theta) \right) = \sum_{k=1}^K \sqrt{\frac{n}{K}} (U_{kn} - \theta)
    \]

    对于给定的总样本量 $N$（从而 $n$ 和 $K$ 固定），定义 $Z_{N,k} = \sqrt{\frac{n}{K}} (U_{kn} - \theta)$。因为 $K$ 个子集是非重叠的，所以 $\{Z_{N,k}\}_{k=1}^K$ 构成了一个行间独立的随机变量三角阵列。我们使用 Lindeberg-Feller 中心极限定理。

    **步骤 1：行和方差的收敛性**

    令 $S_N = \sum_{k=1}^K Z_{N,k}$。由于 $Z_{N,k}$ 对固定的 $N$ 是独立同分布的：

    \[
    Var(S_N) = \sum_{k=1}^K Var(Z_{N,k}) = K \cdot Var\left( \sqrt{\frac{n}{K}} U_{kn} \right) = n Var(U_{kn})
    \]

    从标准 U-统计量理论可知，对于大小为 $n$ 的样本，$n Var(U_{kn}) \rightarrow m^2 \zeta_1$ (当 $n \rightarrow \infty$ 时)。因为 $K = o(N)$，即 $K/N \rightarrow 0$，这意味着 $n = N/K \rightarrow \infty$。因此：

    \[
    Var(S_N) \rightarrow m^2 \zeta_1
    \]

    **步骤 2：验证 Lindeberg 条件**

    我们需要证明对于任意 $\epsilon > 0$：

    \[
    \lim_{N \rightarrow \infty} \sum_{k=1}^K E[Z_{N,k}^2 I(|Z_{N,k}| > \epsilon)] = 0
    \]

    因为 $Z_{N,k}$ 是同分布的，该和简化为：

    \[
    K \cdot E[Z_{N,1}^2 I(|Z_{N,1}| > \epsilon)] = K \cdot E\left[ \frac{n}{K} (U_{1n} - \theta)^2 I\left(\sqrt{\frac{n}{K}}|U_{1n} - \theta| > \epsilon\right) \right]
    \]

    \[
    = E[n(U_{1n} - \theta)^2 I(| \sqrt{n}(U_{1n} - \theta) | > \epsilon \sqrt{K})]
    \]

    令 $W_n = \sqrt{n}(U_{1n} - \theta)$。上述期望变为 $E[W_n^2 I(|W_n| > \epsilon \sqrt{K})]$。
    由正则 U-统计量的 CLT，$W_n \xrightarrow{d} \mathcal{N}(0, m^2\zeta_1)$。
    此外，我们已知 $E[W_n^2] = n Var(U_{1n}) \rightarrow m^2 \zeta_1$。
    收敛分布加上二阶矩的收敛，意味着 $\{W_n^2\}$ 是**一致可积 (uniformly integrable, UI)** 的。
    对于任何一致可积序列，如果截断阈值趋于无穷大，则截断期望趋于 0。因为 $K \rightarrow \infty$，所以阈值 $\epsilon \sqrt{K} \rightarrow \infty$。因此：

    \[
    \lim_{N \rightarrow \infty} E[W_n^2 I(|W_n| > \epsilon \sqrt{K})] = 0
    \]

    这满足了 Lindeberg 条件。

    **步骤 3：结论**

    由于 Lindeberg 条件满足，且行和的方差收敛于 $m^2\zeta_1$，根据 Lindeberg-Feller 中心极限定理：

    \[
    \sum_{k=1}^K Z_{N,k} = \sqrt{N}(AU_N - \theta) \xrightarrow{d} \mathcal{N}(0, m^2 \zeta_1)
    \]

    证明完毕。 $\square$