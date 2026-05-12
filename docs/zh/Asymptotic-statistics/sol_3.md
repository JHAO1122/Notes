---
tags:
  - 大样本理论
  - 统计学
  - 课后习题
---

# 📈 课后习题精解：大样本理论 Assignment 3

## 习题 1: 二阶 Delta 方法 (Second-Order delta method)

**题目：**
设 $X_1, ..., X_n \stackrel{i.i.d}{\sim} Bernoulli(\theta)$。

(a) 当 $\theta = 1/2$ 时，求 $\hat{\theta}_n = \overline{X}_n$ 的渐近分布。

(b) 从 $P_\theta$ 到 $P_\eta$ 的分布的 Kullback-Leibler (KL) 散度定义为：

\[
D_{KL}(P_\eta || P_\theta) := -\mathbb{E}_\eta \log \frac{p_\theta}{p_\eta}(X) \quad , X \sim P_\eta
\]

其中 $p_\theta, p_\eta$ 分别为 $P_\theta$ 和 $P_\eta$ 的概率质量（密度）函数。请证明：

\[
n D_{KL}(P_\theta || P_{\hat{\theta}_n}) \xrightarrow{d} \frac{1}{2}\chi_1^2
\]

??? success "解答 (点击展开)"

    **(a) 求渐近分布**

    因为 $X_1, ..., X_n$ 是独立同分布的 Bernoulli 随机变量，所以：

    \[
    E[X_1] = \theta \quad \text{且} \quad Var(X_1) = \theta(1-\theta)
    \]

    估计量为样本均值 $\hat{\theta}_n = \overline{X}_n = \frac{1}{n}\sum_{i=1}^n X_i$。由经典的中心极限定理 (CLT)，有：

    \[
    \sqrt{n}(\overline{X}_n - \theta) \xrightarrow{d} \mathcal{N}(0, \theta(1-\theta))
    \]

    将已知条件 $\theta = 1/2$ 代入，均值为 $1/2$，方差为 $(1/2)(1 - 1/2) = 1/4$。因此得到：

    \[
    \sqrt{n} \left(\overline{X}_n - \frac{1}{2}\right) \xrightarrow{d} \mathcal{N} \left(0, \frac{1}{4}\right)
    \]

    或者可以等价地写成标准正态分布的形式：

    \[
    \sqrt{n}(2\overline{X}_n - 1) \xrightarrow{d} \mathcal{N}(0, 1)
    \]

    **(b) 证明 KL 散度的渐近分布**

    从 $P_\theta$ 到 $P_{\hat{\theta}_n}$ 的 KL 散度定义为：

    \[
    D_{KL}(P_\theta || P_{\hat{\theta}_n}) = \mathbb{E}_\theta \left[ \log \frac{p_\theta(X)}{p_{\hat{\theta}_n}(X)} \right]
    \]

    对于 Bernoulli 分布，其概率质量函数为 $p_t(x) = t^x (1-t)^{1-x}$。为了便于使用 Delta 方法，我们将 KL 散度视为变量 $t = \hat{\theta}_n$ 处的函数 $g(t)$：

    \[
    g(t) = D_{KL}(P_\theta || P_t) = \mathbb{E}_\theta \left[ X \log \frac{\theta}{t} + (1-X) \log \frac{1-\theta}{1-t} \right]
    \]

    由于是对真实参数 $\theta$ 求期望（即 $\mathbb{E}_\theta[X] = \theta$），上式可以写为：

    \[
    g(t) = \theta \log \left(\frac{\theta}{t}\right) + (1-\theta) \log \left(\frac{1-\theta}{1-t}\right)
    \]

    我们的目标是求 $n g(\hat{\theta}_n) = n(g(\hat{\theta}_n) - g(\theta))$ 的渐近分布（因为 $g(\theta) = D_{KL}(P_\theta || P_\theta) = 0$）。
    
    首先，计算 $g(t)$ 的一阶和二阶导数：

    \[
    g'(t) = -\frac{\theta}{t} + \frac{1-\theta}{1-t}
    \]

    在真实参数 $t=\theta$ 处取值：

    \[
    g'(\theta) = -\frac{\theta}{\theta} + \frac{1-\theta}{1-\theta} = -1 + 1 = 0
    \]

    由于一阶导数为零，我们必须使用**二阶 Delta 方法**。继续求二阶导数：

    \[
    g''(t) = \frac{\theta}{t^2} + \frac{1-\theta}{(1-t)^2}
    \]

    在 $t=\theta$ 处取值：

    \[
    g''(\theta) = \frac{\theta}{\theta^2} + \frac{1-\theta}{(1-\theta)^2} = \frac{1}{\theta} + \frac{1}{1-\theta} = \frac{1}{\theta(1-\theta)}
    \]

    由中心极限定理可知 $\sqrt{n}(\hat{\theta}_n - \theta) \xrightarrow{d} \mathcal{N}(0, \sigma^2)$，其中渐近方差 $\sigma^2 = \theta(1-\theta)$。
    根据二阶 Delta 方法的定理，因为 $g'(\theta) = 0$ 且 $g''(\theta) \neq 0$，我们有：

    \[
    n(g(\hat{\theta}_n) - g(\theta)) \xrightarrow{d} \frac{1}{2} g''(\theta) \sigma^2 \chi_1^2
    \]

    将 $\sigma^2$ 和 $g''(\theta)$ 的结果代入公式：

    \[
    g''(\theta)\sigma^2 = \left( \frac{1}{\theta(1-\theta)} \right) \cdot \theta(1-\theta) = 1
    \]

    这完美化简为 1。因此：

    \[
    n D_{KL}(P_\theta || P_{\hat{\theta}_n}) = n(g(\hat{\theta}_n) - 0) \xrightarrow{d} \frac{1}{2} \chi_1^2
    \]

---

## 习题 2: Hodge (1951) 超有效估计 (Superefficiency estimator)

**题目：**
假设 $X_1, X_2, \dots, X_n$ 是从 $\mathcal{N}(\theta, 1)$ 中抽取的独立同分布观测样本，其中 $\theta \in \mathbb{R}$。此时 $\overline{X}_n = n^{-1}\sum_{i=1}^n X_i$ 是 $\theta$ 的极大似然估计。请证明：当 $\theta \neq 0$ 时，

\[
P(|\overline{X}_n| \le n^{-1/4}) \rightarrow 0
\]

??? success "解答 (点击展开)"

    **证明过程：**

    由于 $X_1, X_2, \dots, X_n \stackrel{i.i.d}{\sim} \mathcal{N}(\theta, 1)$，样本均值 $\overline{X}_n$ 的精确分布为：

    \[
    \overline{X}_n \sim \mathcal{N} \left( \theta, \frac{1}{n} \right)
    \]

    我们将 $\overline{X}_n$ 标准化为一个标准正态分布随机变量 $Z$：

    \[
    Z = \sqrt{n}(\overline{X}_n - \theta) \sim \mathcal{N}(0, 1)
    \]

    现在，我们使用 $Z$ 重写目标概率：

    \[
    P(|\overline{X}_n| \le n^{-1/4}) = P(-n^{-1/4} \le \overline{X}_n \le n^{-1/4})
    \]

    \[
    = P(-n^{-1/4} - \theta \le \overline{X}_n - \theta \le n^{-1/4} - \theta)
    \]

    \[
    = P(-n^{1/4} - \sqrt{n}\theta \le \sqrt{n}(\overline{X}_n - \theta) \le n^{1/4} - \sqrt{n}\theta)
    \]

    \[
    = P(-n^{1/4} - \sqrt{n}\theta \le Z \le n^{1/4} - \sqrt{n}\theta)
    \]

    为了证明该概率在 $\theta \neq 0$ 时收敛于 0，我们分两种情况进行讨论：

    **情况 1：$\theta > 0$**
    
    考察不等式中 $Z$ 的上限：

    \[
    n^{1/4} - \sqrt{n}\theta = \sqrt{n}(n^{-1/4} - \theta)
    \]

    当 $n \rightarrow \infty$ 时，$n^{-1/4} \rightarrow 0$。因为 $\theta > 0$，所以项 $(n^{-1/4} - \theta) \rightarrow -\theta < 0$。因此：

    \[
    \lim_{n \rightarrow \infty} (n^{1/4} - \sqrt{n}\theta) = -\infty
    \]

    该概率受限于标准正态分布在上限处的累积分布函数 (CDF)：

    \[
    P(-n^{1/4} - \sqrt{n}\theta \le Z \le n^{1/4} - \sqrt{n}\theta) \le P(Z \le n^{1/4} - \sqrt{n}\theta)
    \]

    当 $n \rightarrow \infty$ 时，$P(Z \le n^{1/4} - \sqrt{n}\theta) \rightarrow \Phi(-\infty) = 0$。

    **情况 2：$\theta < 0$**
    
    考察不等式中 $Z$ 的下限：

    \[
    -n^{1/4} - \sqrt{n}\theta = -\sqrt{n}(n^{-1/4} + \theta)
    \]

    当 $n \rightarrow \infty$ 时，$n^{-1/4} \rightarrow 0$。由于 $\theta < 0$，项 $(n^{-1/4} + \theta) \rightarrow \theta < 0$。因此：

    \[
    \lim_{n \rightarrow \infty} -\sqrt{n}(n^{-1/4} + \theta) = +\infty
    \]

    该概率受限于标准正态分布的右尾概率：

    \[
    P(-n^{1/4} - \sqrt{n}\theta \le Z \le n^{1/4} - \sqrt{n}\theta) \le P(Z \ge -n^{1/4} - \sqrt{n}\theta)
    \]

    当 $n \rightarrow \infty$ 时，$P(Z \ge -n^{1/4} - \sqrt{n}\theta) \rightarrow 1 - \Phi(+\infty) = 0$。

    在这两种 $\theta \neq 0$ 的情况下，概率均随之消失。因此，当 $n \rightarrow \infty$ 时：

    \[
    P(|\overline{X}_n| \le n^{-1/4}) \rightarrow 0
    \]

---

## 习题 3: 非正则极大似然估计 (Irregular Maximum likelihood estimator)

**题目：**
设 $X_1, X_2, \dots, X_n$ 为来自均匀分布 $U(0, \theta)$（$\theta > 0$）的随机样本。

(a) 获得 $\theta$ 的极大似然估计量 $\hat{\theta}_n$，并证明 $\sqrt{n}(\hat{\theta}_n - \theta) \xrightarrow{p} 0$。

(b) 给出一个适当的缩放因子，使得其渐近分布是非退化的，并求出该渐近分布。

??? success "解答 (点击展开)"

    **(a) 求 MLE 并证明收敛性**

    样本的似然函数为：

    \[
    L(\theta) = \prod_{i=1}^n \frac{1}{\theta} I(0 \le X_i \le \theta) = \frac{1}{\theta^n} I(\max_{1 \le i \le n} X_i \le \theta)
    \]

    其中 $I(\cdot)$ 为指示函数。为了最大化 $L(\theta)$，我们需要在 $\theta \ge \max_{1 \le i \le n} X_i$ 的约束下，最小化分母 $\theta^n$。
    
    由于 $\theta^n$ 对 $\theta > 0$ 是严格递增函数，使似然函数最大的 $\theta$ 值即为满足约束的最小可能值，也就是最大次序统计量：

    \[
    \hat{\theta}_n = X_{(n)} = \max_{1 \le i \le n} X_i
    \]

    为了证明 $\sqrt{n}(\hat{\theta}_n - \theta) \xrightarrow{p} 0$，我们分析偏差 $\theta - \hat{\theta}_n$。显然 $\hat{\theta}_n \le \theta$ 几乎必然成立，故 $\theta - \hat{\theta}_n \ge 0$。我们可以使用马尔可夫 (Markov) 不等式。
    
    首先计算 $\hat{\theta}_n$ 的期望。$X_{(n)}$ 的概率密度函数 (PDF) 在 $0 < x < \theta$ 时为：

    \[
    f_{(n)}(x) = \frac{n x^{n-1}}{\theta^n}
    \]

    因此期望为：

    \[
    E[\hat{\theta}_n] = \int_0^\theta x \cdot \frac{n x^{n-1}}{\theta^n} dx = \frac{n}{\theta^n} \left[ \frac{x^{n+1}}{n+1} \right]_0^\theta = \frac{n}{n+1}\theta
    \]

    期望偏差为：

    \[
    E[\theta - \hat{\theta}_n] = \theta - \frac{n}{n+1}\theta = \frac{\theta}{n+1}
    \]

    对于任意 $\epsilon > 0$，由马尔可夫不等式：

    \[
    P(|\sqrt{n}(\hat{\theta}_n - \theta)| \ge \epsilon) = P(\sqrt{n}(\theta - \hat{\theta}_n) \ge \epsilon) \le \frac{E[\sqrt{n}(\theta - \hat{\theta}_n)]}{\epsilon} = \frac{\sqrt{n}\theta}{\epsilon(n+1)}
    \]

    当 $n \rightarrow \infty$ 时，$\frac{\sqrt{n}}{n+1} \rightarrow 0$。概率趋于 0，这意味着：

    \[
    \sqrt{n}(\hat{\theta}_n - \theta) \xrightarrow{p} 0
    \]

    **(b) 寻找适当的缩放因子及非退化渐近分布**

    从 (a) 问的结果 $\sqrt{n}(\theta - \hat{\theta}_n) \xrightarrow{p} 0$ 可以看出，$\hat{\theta}_n$ 收敛到 $\theta$ 的速度比常规的 $\sqrt{n}$ 还要快。此处适当的缩放因子是 $n$。
    
    我们要寻找非负随机变量 $Y_n = n(\theta - \hat{\theta}_n)$ 的渐近分布。我们通过推导其累积分布函数 (CDF) 的极限来求解。对于任意 $x > 0$：

    \[
    P(n(\theta - \hat{\theta}_n) \le x) = P \left( \theta - \hat{\theta}_n \le \frac{x}{n} \right)
    \]

    \[
    = P \left( \hat{\theta}_n \ge \theta - \frac{x}{n} \right) = 1 - P \left( \hat{\theta}_n < \theta - \frac{x}{n} \right)
    \]

    由于 $\hat{\theta}_n = \max_{1 \le i \le n} X_i$，样本最大值严格小于某一个数值，当且仅当每一个独立的观测值都严格小于该数值。因为 $X_i$ 是独立同分布于 $U(0, \theta)$，所以：

    \[
    P \left( \hat{\theta}_n < \theta - \frac{x}{n} \right) = \prod_{i=1}^n P \left( X_i < \theta - \frac{x}{n} \right) = \left( \frac{\theta - x/n}{\theta} \right)^n = \left( 1 - \frac{x/\theta}{n} \right)^n
    \]

    当 $n \rightarrow \infty$ 取极限时，利用指数函数的标准极限定义：

    \[
    \lim_{n \rightarrow \infty} \left( 1 - \frac{x/\theta}{n} \right)^n = e^{-x/\theta}
    \]

    因此，$n(\theta - \hat{\theta}_n)$ 的渐近 CDF 为：

    \[
    \lim_{n \rightarrow \infty} P(n(\theta - \hat{\theta}_n) \le x) = 1 - e^{-x/\theta} \quad (\text{对于 } x > 0)
    \]

    这正是均值为 $\theta$（或者说率参数为 $1/\theta$）的**指数分布 (Exponential Distribution)** 的 CDF。
    
    结论，非退化的渐近分布为：

    \[
    n(\theta - \hat{\theta}_n) \xrightarrow{d} Exp(mean = \theta)
    \]

## 习题 4: 极大似然估计的不同收敛速度 (Different convergence rate of MLEs)

**题目：**
假设 $X_1, ..., X_n \stackrel{i.i.d}{\sim} Exp(\mu, \sigma)$，其概率密度函数为 $f(x) = \frac{1}{\sigma}e^{-\frac{x-\mu}{\sigma}}$，其中 $x \ge \mu$。记参数 $(\mu, \sigma)$ 的极大似然估计 (MLE) 为 $(\hat{\mu}, \hat{\sigma})$。
请分别求出 $\hat{\mu}$ 和 $\hat{\sigma}$ 的渐近分布。（求边缘分布，非联合分布）

??? success "解答 (点击展开)"

    **第一步：求 $\mu$ 和 $\sigma$ 的极大似然估计 (MLE)**

    样本的似然函数为：

    \[
    L(\mu, \sigma) = \prod_{i=1}^n \frac{1}{\sigma}e^{-\frac{X_i-\mu}{\sigma}}I(X_i \ge \mu) = \frac{1}{\sigma^n} \exp \left\{ -\frac{\sum_{i=1}^n X_i - n\mu}{\sigma} \right\} I(\mu \le \min_{1 \le i \le n} X_i)
    \]

    在满足 $\mu \le X_{(1)}$ 的条件下，对数似然函数为：

    \[
    \ln L(\mu, \sigma) = -n \ln \sigma - \frac{\sum_{i=1}^n X_i - n\mu}{\sigma}
    \]

    首先关于 $\mu$ 求最大值。注意到：

    \[
    \frac{\partial \ln L}{\partial \mu} = \frac{n}{\sigma} > 0
    \]

    因此，似然函数关于 $\mu$ 是严格递增的。为了使其最大化，我们必须在约束条件 $\mu \le X_{(1)}$ 下使 $\mu$ 尽可能大。
    所以，$\mu$ 的极大似然估计为样本最小值：

    \[
    \hat{\mu} = X_{(1)} = \min_{1 \le i \le n} X_i
    \]

    接着，对 $\sigma$ 求偏导并令其为 0：

    \[
    \frac{\partial \ln L}{\partial \sigma} = -\frac{n}{\sigma} + \frac{\sum_{i=1}^n X_i - n\mu}{\sigma^2} = 0
    \]

    解得 $\hat{\sigma} = \overline{X}_n - \mu$。将 $\hat{\mu} = X_{(1)}$ 代入该等式，得到 $\sigma$ 的 MLE：

    \[
    \hat{\sigma} = \overline{X}_n - X_{(1)}
    \]

    **第二步：求 $\hat{\mu}$ 的渐近分布**

    我们利用累积分布函数 (CDF) 方法来求 $\hat{\mu} = X_{(1)}$ 的渐近分布。此处适当的缩放因子是 $n$。
    对于任意的 $x > 0$：

    \[
    P(n(\hat{\mu} - \mu) > x) = P \left( \hat{\mu} > \mu + \frac{x}{n} \right)
    \]

    \[
    = P \left( \min_{1 \le i \le n} X_i > \mu + \frac{x}{n} \right) = \prod_{i=1}^n P \left( X_i > \mu + \frac{x}{n} \right)
    \]

    因为 $X_i \sim Exp(\mu, \sigma)$，对于 $t \ge \mu$，有 $P(X_i > t) = e^{-(t-\mu)/\sigma}$。因此：

    \[
    P \left( X_i > \mu + \frac{x}{n} \right) = \exp \left\{ -\frac{\mu + x/n - \mu}{\sigma} \right\} = e^{-\frac{x}{n\sigma}}
    \]

    将其代回乘积中：

    \[
    P(n(\hat{\mu} - \mu) > x) = \left( e^{-\frac{x}{n\sigma}} \right)^n = e^{-x/\sigma}
    \]

    这意味着极限分布的 CDF 为 $1 - e^{-x/\sigma}$，这正好是位置参数为 $0$、尺度参数为 $\sigma$（均值为 $\sigma$）的**指数分布**的 CDF。因此：

    \[
    n(\hat{\mu} - \mu) \xrightarrow{d} Exp(mean = \sigma)
    \]

    **第三步：求 $\hat{\sigma}$ 的渐近分布**

    我们需要求 $\hat{\sigma} = \overline{X}_n - X_{(1)}$ 的渐近分布。
    首先，注意到 $X_1$ 的矩为：$E[X_1] = \mu + \sigma$ 且 $Var(X_1) = \sigma^2$。根据关于 $\overline{X}_n$ 的标准中心极限定理 (CLT)：

    \[
    \sqrt{n}(\overline{X}_n - (\mu + \sigma)) \xrightarrow{d} \mathcal{N}(0, \sigma^2)
    \]

    我们将标准化的估计量 $\sqrt{n}(\hat{\sigma} - \sigma)$ 重写为：

    \[
    \sqrt{n}(\hat{\sigma} - \sigma) = \sqrt{n}(\overline{X}_n - X_{(1)} - \sigma)
    \]

    \[
    = \sqrt{n}(\overline{X}_n - \mu - \sigma) - \sqrt{n}(X_{(1)} - \mu)
    \]

    由第二步可知，$n(X_{(1)} - \mu) \xrightarrow{d} Exp(\sigma)$，这意味着 $n(X_{(1)} - \mu) = O_p(1)$。
    因此：

    \[
    \sqrt{n}(X_{(1)} - \mu) = \frac{1}{\sqrt{n}} \cdot n(X_{(1)} - \mu) = o(1) \cdot O_p(1) = o_p(1)
    \]

    这说明 $\sqrt{n}(X_{(1)} - \mu) \xrightarrow{p} 0$。
    
    根据 Slutsky 定理，由于第二项依概率收敛于 0，渐近分布完全由第一项决定。因此：

    \[
    \sqrt{n}(\hat{\sigma} - \sigma) \xrightarrow{d} \mathcal{N}(0, \sigma^2)
    \]

---

## 习题 5: 不一致的极大似然估计与 Neyman-Scott 问题 (Inconsistent MLE and Neyman-Scott problem)

**题目：**
设 $X_{ij}$（$i=1,...,n$ 且 $j=1,...,k$）相互独立，且 $X_{ij} \sim \mathcal{N}(\mu_i, \sigma^2)$。
注意，这基本上是一个平衡的单因素 ANOVA 设计，我们假设 $k$ 是固定的而 $n \rightarrow \infty$，因此随着总样本量的增加，参数的维度也趋于无穷大。
证明公共方差 $\sigma^2$ 的极大似然估计 $\hat{\sigma}^2$ 满足：

\[
\hat{\sigma}^2 \xrightarrow{p} \left(1 - \frac{1}{k}\right)\sigma^2
\]

??? success "解答 (点击展开)"

    **第一步：导出 $\mu_i$ 和 $\sigma^2$ 的极大似然估计量**

    给定样本的似然函数为：

    \[
    L(\mu_1, ..., \mu_n, \sigma^2) = \prod_{i=1}^n \prod_{j=1}^k \frac{1}{\sqrt{2\pi\sigma^2}} \exp \left( -\frac{(X_{ij} - \mu_i)^2}{2\sigma^2} \right)
    \]

    取自然对数，得到对数似然函数：

    \[
    \ln L = -\frac{nk}{2} \ln(2\pi) - \frac{nk}{2} \ln(\sigma^2) - \frac{1}{2\sigma^2} \sum_{i=1}^n \sum_{j=1}^k (X_{ij} - \mu_i)^2
    \]

    为了寻找 $\mu_i$ 的 MLE，我们对 $\mu_i$ 求偏导并令其为 0：

    \[
    \frac{\partial \ln L}{\partial \mu_i} = \frac{1}{\sigma^2} \sum_{j=1}^k (X_{ij} - \mu_i) = 0 \implies \hat{\mu}_i = \frac{1}{k} \sum_{j=1}^k X_{ij} = \overline{X}_{i.}
    \]

    接着，我们对 $\sigma^2$ 求偏导并令其为 0 以寻找 $\sigma^2$ 的 MLE：

    \[
    \frac{\partial \ln L}{\partial (\sigma^2)} = -\frac{nk}{2\sigma^2} + \frac{1}{2(\sigma^2)^2} \sum_{i=1}^n \sum_{j=1}^k (X_{ij} - \mu_i)^2 = 0
    \]

    将 $\hat{\mu}_i = \overline{X}_{i.}$ 代入上式，我们得到公共方差的 MLE：

    \[
    \hat{\sigma}^2 = \frac{1}{nk} \sum_{i=1}^n \sum_{j=1}^k (X_{ij} - \overline{X}_{i.})^2
    \]

    **第二步：证明依概率收敛**

    令 $Y_i = \sum_{j=1}^k (X_{ij} - \overline{X}_{i.})^2$。对于每一个组 $i$，因为对于 $j=1,...,k$ 都有 $X_{ij} \sim \mathcal{N}(\mu_i, \sigma^2)$，根据标准统计性质，缩放后的离差平方和服从卡方分布：

    \[
    \frac{Y_i}{\sigma^2} = \frac{1}{\sigma^2} \sum_{j=1}^k (X_{ij} - \overline{X}_{i.})^2 \sim \chi_{k-1}^2
    \]

    因此，$Y_i \sim \sigma^2 \chi_{k-1}^2$。由于不同组（$i=1,...,n$）之间的样本是相互独立的，序列 $Y_1, Y_2, ..., Y_n$ 由独立同分布 (i.i.d.) 的随机变量组成。

    $Y_i$ 的期望值为 $E[Y_i] = \sigma^2(k-1)$。

    当 $n \rightarrow \infty$（而 $k$ 保持固定）时，我们对序列 $\{Y_i\}_{i=1}^n$ 应用弱大数定律 (WLLN)：

    \[
    \frac{1}{n} \sum_{i=1}^n Y_i \xrightarrow{p} E[Y_1] = (k-1)\sigma^2
    \]

    现在，我们将估计量 $\hat{\sigma}^2$ 用该序列重写：

    \[
    \hat{\sigma}^2 = \frac{1}{nk} \sum_{i=1}^n Y_i = \frac{1}{k} \left( \frac{1}{n} \sum_{i=1}^n Y_i \right)
    \]

    由连续映射定理（或依概率收敛的基本性质），乘以常数 $1/k$ 可得：

    \[
    \hat{\sigma}^2 \xrightarrow{p} \frac{1}{k} \cdot (k-1)\sigma^2 = \left( 1 - \frac{1}{k} \right)\sigma^2
    \]

    证明完毕。这个结果展示了著名的 **Neyman-Scott 问题**：公共方差的极大似然估计是不一致的，因为讨厌参数（nuisance parameters, $\mu_1, ..., \mu_n$）的数量以与样本量 $n$ 相同的速率增长。


## 习题 6: 样本分位数的相合性与渐近正态性 (Consistency and Asymptotic Normality of Sample Quantiles)

**题目：**
(a) 定义检验函数 (check function) $\rho_\tau(y) = y(\tau - I(y < 0))$ 作为损失函数，并设分布函数 $F$ 的分位数函数为 $F^{-1}(p) := \inf\{x \in \mathbb{R} : F(x) \ge p\}$。请证明：

\[
\theta_0 = \arg \min_\theta E\rho_\tau(X-\theta) = F^{-1}(\tau)
\]

(b) 验证由以下 Z-估计量定义的样本分位数的相合性：

\[
\Psi_n(\theta) := \frac{1}{n} \sum_{i=1}^n \left( (1-\tau)I(X_i < \theta) - \tau I(X_i > \theta) \right) = 0
\]

提示：可使用 van der Vaart [1] 中的定理 5.7 或引理 5.10（如有必要，可自行添加正则性条件）。

(c) (van der Vaart [1] Ex. 5.11) 在已知相合性的假设下，利用定理 5.23 证明样本分位数的渐近正态性。

??? success "解答 (点击展开)"

    **(a) 证明期望损失的极小值点为总体分位数**

    我们假设存在一个正则性条件：随机变量 $X$ 具有连续的概率密度函数 $f(x)$ 和累积分布函数 $F(x)$。
    定义目标函数为 $M(\theta) = E[\rho_\tau(X-\theta)]$。利用检验函数 $\rho_\tau(y) = y(\tau - I(y < 0))$ 的定义，我们可以将期望写为积分形式：

    \[
    M(\theta) = \int_{-\infty}^{\infty} (x-\theta)(\tau - I(x-\theta < 0)) f(x) dx
    \]

    将其拆分为两部分积分：

    \[
    M(\theta) = \int_{-\infty}^{\theta} (x-\theta)(\tau - 1) f(x) dx + \int_{\theta}^{\infty} (x-\theta)\tau f(x) dx
    \]

    为了求极小值，我们利用莱布尼茨积分法则 (Leibniz's rule) 对 $\theta$ 求导。注意到当 $x=\theta$ 时，边界项 $(x-\theta)$ 为 $0$，因此边界项求导后消失：

    \[
    M'(\theta) = \int_{-\infty}^{\theta} \frac{\partial}{\partial \theta} [(x-\theta)(\tau-1)] f(x) dx + \int_{\theta}^{\infty} \frac{\partial}{\partial \theta} [(x-\theta)\tau] f(x) dx
    \]

    \[
    = \int_{-\infty}^{\theta} (1-\tau) f(x) dx + \int_{\theta}^{\infty} (-\tau) f(x) dx
    \]

    计算积分得到：

    \[
    = (1-\tau)F(\theta) - \tau(1-F(\theta)) = F(\theta) - \tau F(\theta) - \tau + \tau F(\theta) = F(\theta) - \tau
    \]

    令导数为 $0$ 以求极小值：

    \[
    M'(\theta_0) = 0 \implies F(\theta_0) = \tau
    \]

    根据广义逆分布函数（分位数函数）的定义 $F^{-1}(\tau) = \inf\{x \in \mathbb{R} : F(x) \ge \tau\}$，我们得出结论：

    \[
    \theta_0 = F^{-1}(\tau)
    \]

    **(b) 验证样本分位数的相合性 (Consistency)**

    我们将利用 **引理 5.10** 来证明 $\hat{\theta}_n$ 的相合性。

    !!! info "引理 5.10 (van der Vaart)"

        设 $\Psi_n(\theta)$ 为随机函数，$\Psi(\theta)$ 为固定函数，且对每个 $\theta$，$\Psi_n(\theta)$ 依概率收敛于 $\Psi(\theta)$。
        假设映射 $\theta \mapsto \Psi_n(\theta)$ 是**非递减的 (nondecreasing)**，且 $\Psi_n(\hat{\theta}_n) = o_p(1)$。
        若 $\theta_0$ 满足对于任意 $\epsilon > 0$ 都有 $\Psi(\theta_0 - \epsilon) < 0 < \Psi(\theta_0 + \epsilon)$，则 $\hat{\theta}_n \xrightarrow{p} \theta_0$。

    我们需要逐一验证该引理的条件：

    * **条件 1：$\Psi_n(\theta)$ 的单调性**
    
    随着 $\theta$ 的增加，指示函数 $I(X_i < \theta)$ 是非递减的，而 $I(X_i > \theta)$ 是非递增的（从而 $-\tau I(X_i > \theta)$ 也是非递减的）。由于 $(1-\tau) > 0$ 且 $\tau > 0$，映射 $\theta \mapsto \Psi_n(\theta)$ 是多个非递减函数的正线性组合，因此整体上是**非递减的**。

    * **条件 2：依概率收敛到固定函数 $\Psi(\theta)$**
    
    由弱大数定律 (WLLN)，对于任何固定的 $\theta$，$\Psi_n(\theta)$ 依概率收敛于它的期望值。
    假设 $X$ 是连续随机变量（即 $P(X=\theta) = 0$）：

    \[
    \Psi(\theta) = E[\Psi_n(\theta)] = (1-\tau)P(X < \theta) - \tau P(X > \theta) = (1-\tau)F(\theta) - \tau(1-F(\theta)) = F(\theta) - \tau
    \]

    因此，$\Psi_n(\theta) \xrightarrow{p} \Psi(\theta)$ 成立。

    * **条件 3：唯一的零点穿越 (添加正则性条件)**
    
    引理要求对任意 $\epsilon > 0$，有 $\Psi(\theta_0 - \epsilon) < 0 < \Psi(\theta_0 + \epsilon)$。
    为此，我们必须添加一个正则性条件：**假设分布函数 $F(x)$ 在真实参数 $\theta_0 = F^{-1}(\tau)$ 的邻域内是严格递增的**。
    在该条件下，对于任意 $\epsilon > 0$：

    \[
    F(\theta_0 - \epsilon) < F(\theta_0) = \tau < F(\theta_0 + \epsilon)
    \]

    两边同时减去 $\tau$，正好得到：

    \[
    \Psi(\theta_0 - \epsilon) < 0 < \Psi(\theta_0 + \epsilon)
    \]

    由于 $\theta \mapsto \Psi_n(\theta)$ 是非递减的，且样本分位数使 $\Psi_n(\hat{\theta}_n) = 0$（或 $o_p(1)$），且真实参数 $\theta_0$ 唯一地穿过极限函数 $\Psi(\theta)$ 的零点，引理 5.10 的所有条件均已满足。由此得出样本分位数是相合的：

    \[
    \hat{\theta}_n \xrightarrow{p} \theta_0
    \]

    **(c) 证明样本分位数的渐近正态性**

    为了证明渐近正态性，我们将使用 **定理 5.23**。为了与教材的符号习惯保持一致，这里我们将分位数记为 $p$（即上一问中的 $\tau$）。

    !!! info "定理 5.23 (M-估计量的渐近正态性)"

        假设映射 $\theta \mapsto m_\theta(x)$ 在 $\theta_0$ 处对 $P$-几乎所有的 $x$ 可微，导数为 $\dot{m}_{\theta_0}(x)$。且满足 Lipschitz 条件：
        $|m_{\theta_1}(x) - m_{\theta_2}(x)| \le \dot{m}(x)|\theta_1 - \theta_2|$，其中包络函数满足 $P\dot{m}^2 < \infty$。
        进一步假设 $\theta \mapsto P m_\theta$ 在极大值点 $\theta_0$ 处存在二阶泰勒展开，其二阶导数矩阵 $V_{\theta_0}$ 是非奇异对称的。
        如果 $\hat{\theta}_n \xrightarrow{p} \theta_0$ 且近乎最大化经验目标函数，则有：
        $\sqrt{n}(\hat{\theta}_n - \theta_0) = -V_{\theta_0}^{-1} \frac{1}{\sqrt{n}} \sum_{i=1}^n \dot{m}_{\theta_0}(X_i) + o_p(1)$。

    样本 $p$-分位数 $\hat{\theta}_n = \mathbb{F}_n^{-1}(p)$ 最小化了经验目标函数 $\frac{1}{n}\sum \rho_p(X_i - \theta)$。为了将其转化为定理 5.23 中的“最大化”框架，我们定义目标函数为：

    \[
    m_\theta(x) = -\rho_p(x-\theta) = -(x-\theta)(p - I(x < \theta))
    \]

    我们逐步验证定理 5.23 在真实参数 $\theta_0 = F^{-1}(p)$ 处的各项条件：

    * **条件 1：$P$-几乎处处可微**
    
    映射 $\theta \mapsto m_\theta(x)$ 除了在 $\theta = x$ 处不可微之外处处可微。由于分布 $F$ 在 $\theta_0$ 处可微且存在密度 $f(\theta_0) > 0$，观测值恰好等于 $\theta_0$ 的概率为零，即 $P(X=\theta_0)=0$。因此 $m_\theta(x)$ 对 $P$-几乎所有的 $x$ 在 $\theta_0$ 处可微，其导数为：

    \[
    \dot{m}_{\theta_0}(x) = \frac{\partial}{\partial \theta} [-(x-\theta)p + (x-\theta)I(x < \theta)] \Big|_{\theta=\theta_0} = p - I(x \le \theta_0)
    \]

    * **条件 2：Lipschitz 连续性与平方可积包络**
    
    检验函数 $\rho_p(u)$ 具有 Lipschitz 连续性，其常数受限于 $\max(p, 1-p) \le 1$。
    因此，对于任意的 $\theta_1, \theta_2$：

    \[
    |m_{\theta_1}(x) - m_{\theta_2}(x)| = |\rho_p(x-\theta_1) - \rho_p(x-\theta_2)| \le 1 \cdot |\theta_1 - \theta_2|
    \]

    我们可以选择包络函数 $\dot{m}(x) = 1$。条件 $P\dot{m}^2 = E[1^2] = 1 < \infty$ 显然满足。

    * **条件 3：二阶泰勒展开与非奇异的 $V_{\theta_0}$**
    
    期望目标函数为 $Pm_\theta = -E[\rho_p(X-\theta)]$。由 (a) 问的结论，第一阶导数为：

    \[
    \frac{d}{d\theta} (Pm_\theta) = p - F(\theta)
    \]

    再次求导，我们得到 $\theta_0$ 处的二阶导数：

    \[
    V_{\theta_0} = \frac{d}{d\theta} (p - F(\theta)) \Big|_{\theta=\theta_0} = -f(\theta_0)
    \]

    题目假设在总体 $p$-分位数 $F^{-1}(p)$ 处，$F$ 具有正导数，这意味着 $f(\theta_0) > 0$。因此，标量 $V_{\theta_0} = -f(\theta_0)$ 严格为负，满足“非奇异”的条件。

    * **根据定理 5.23 得出结论**
    
    由于在 (b) 问中已经确立了相合性 $\hat{\theta}_n \xrightarrow{p} \theta_0$，所有条件均已满足。应用定理 5.23：

    \[
    \sqrt{n}(\hat{\theta}_n - \theta_0) = -V_{\theta_0}^{-1} \frac{1}{\sqrt{n}} \sum_{i=1}^n \dot{m}_{\theta_0}(X_i) + o_p(1)
    \]

    将 $V_{\theta_0}$ 和 $\dot{m}_{\theta_0}(X_i)$ 代入方程：

    \[
    \sqrt{n}(\hat{\theta}_n - \theta_0) = \frac{1}{f(\theta_0)} \frac{1}{\sqrt{n}} \sum_{i=1}^n (p - I(X_i \le \theta_0)) + o_p(1)
    \]

    令 $Y_i = p - I(X_i \le \theta_0)$。序列 $\{Y_i\}_{i=1}^n$ 由独立同分布 (i.i.d.) 随机变量组成，且满足：
    $E[Y_i] = p - P(X_i \le \theta_0) = p - F(\theta_0) = p - p = 0$
    $Var(Y_i) = Var(I(X_i \le \theta_0)) = p(1-p)$

    根据经典的中心极限定理 (CLT)：

    \[
    \frac{1}{\sqrt{n}} \sum_{i=1}^n Y_i \xrightarrow{d} \mathcal{N}(0, p(1-p))
    \]

    最后，根据 Slutsky 定理，乘以常数缩放因子 $1/f(\theta_0)$ 后，便得到了其渐近正态分布：

    \[
    \sqrt{n}(\hat{\theta}_n - \theta_0) \xrightarrow{d} \mathcal{N}\left(0, \frac{p(1-p)}{f(\theta_0)^2}\right)
    \]