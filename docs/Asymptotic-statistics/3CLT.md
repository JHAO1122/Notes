# 第三章：中心极限定理(一) 

与经典的经典中心极限定理（要求独立同分布 i.i.d.）不同，本章我们将探索更具一般性的中心极限定理，即随机变量是**独立但不同分布 (independent but not identically distributed, i.n.i.d.)** 的情况。

---

## 1. 双重数组与复数极限定理 (Double Arrays & A Lemma)

在处理不同分布的随机变量之和时，我们通常将其表示为“双重数组”（或三角阵列）的形式。

!!! abstract "定义 3.1：双重独立随机向量数组 (Double Array)"

    对于每一个 $n \ge 1$，设 $\{X_{n1}, X_{n2}, \dots, X_{nk_n}\}$ 是定义在概率空间 $(\Omega_n, \mathcal{F}_n, P_n)$ 上的一组随机向量，满足对于给定的 $n$，$X_{n1}, \dots, X_{nk_n}$ 相互独立，且当 $n \to \infty$ 时，$k_n \to \infty$。
    那么称 $\{X_{nj} : 1 \le j \le k_n\}_{n \ge 1}$ 为一个**双重独立随机向量数组**。

    **本章常用记号**：
    
    * 期望：$\alpha_{nj} = E(X_{nj})$，行的总期望 $\alpha_n = \sum_{j=1}^{k_n} E(X_{nj}) = \sum_{j=1}^{k_n} \alpha_{nj}$
    * 部分和：$S_n = \sum_{j=1}^{k_n} X_{nj}$
    * 方差：$\sigma_{nj}^2 = Var(X_{nj})$，行的总方差 $\sigma_n^2 = \sum_{j=1}^{k_n} \sigma_{nj}^2$

为了处理特征函数的乘积，我们需要引入一个复数序列乘积的极限引理。

!!! info "引理 3.2：复数序列的乘积极限"

    设 $\{\theta_{nj} : 1 \le j \le k_n\}_{n \ge 1}$ 为一个复数双重数组，满足当 $n \to \infty$ 时：
    
    **(i)** 均匀趋于 0：$\max_{1 \le j \le k_n} |\theta_{nj}| \to 0$
    **(ii)** 绝对和一致有界：$\sum_{j=1}^{k_n} |\theta_{nj}| \le M < \infty$（其中 $M$ 与 $n$ 无关）
    **(iii)** 和收敛：$\sum_{j=1}^{k_n} \theta_{nj} \to \theta$（其中 $\theta$ 为有限复数）
    
    那么其连乘积收敛于指数函数：
    
    $$
    \prod_{j=1}^{k_n} (1 + \theta_{nj}) \rightarrow e^\theta
    $$
    
    > *(注：这推广了微积分中经典的 $\lim_{n \to \infty} (1 + \theta/n)^n = e^\theta$ 结论，即将 $\theta_{nj} \equiv \theta/n$ 代入该引理。)*

    ??? proof "引理 3.2 的详细证明（点击展开）"

        对于非零复数 $z$，复对数主值定义为 $Log~z = \log|z| + i Arg~z$，其中 $Arg~z \in [-\pi, \pi]$。
        当 $|z| < 1$ 时，复对数存在如下泰勒级数展开：
        
        $$
        \log(1+z) = \sum_{m=1}^\infty (-1)^{m-1} \frac{z^m}{m}
        $$
        
        由条件 (i)，存在 $n_0$，使得对所有的 $n > n_0$ 有 $\max_{1 \le j \le k_n} |\theta_{nj}| \le 1/2$。
        此时 $|\theta_{nj}| < 1$ 且 $1+\theta_{nj} \neq 0$。考察其对数展开与一次项的截断误差：
        
        $$
        |\log(1+\theta_{nj}) - \theta_{nj}| = \left| \sum_{m \ge 2} (-1)^{m-1} \frac{\theta_{nj}^m}{m} \right| \le \sum_{m \ge 2} \frac{|\theta_{nj}|^m}{m}
        $$
        
        提取出二次项，并利用等比数列求和界定剩余项：
        
        $$
        \le \frac{|\theta_{nj}|^2}{2} \sum_{m=2}^\infty \left(\frac{1}{2}\right)^{m-2} = |\theta_{nj}|^2 < 1
        $$
        
        由于误差绝对值被界定在 $|\theta_{nj}|^2$ 内，我们可以将其写为：
        
        $$
        \log(1+\theta_{nj}) = \theta_{nj} + \Lambda_{nj}|\theta_{nj}|^2, \quad \text{其中 } |\Lambda_{nj}| < 1
        $$
        
        对一行求和：
        
        $$
        \sum_{j=1}^{k_n} \log(1+\theta_{nj}) = \sum_{j=1}^{k_n} \theta_{nj} + \sum_{j=1}^{k_n} \Lambda_{nj}|\theta_{nj}|^2
        $$
        
        利用条件 (i) 和 (ii)，估计误差项的总和：
        
        $$
        \left| \sum_{j=1}^{k_n} \Lambda_{nj}|\theta_{nj}|^2 \right| \le \max_{1 \le j \le k_n} |\theta_{nj}| \sum_{j=1}^{k_n} |\theta_{nj}| \le \left( \max_{1 \le j \le k_n} |\theta_{nj}| \right) \cdot M \rightarrow 0
        $$
        
        结合条件 (iii) $\sum_{j=1}^{k_n} \theta_{nj} \to \theta$，可得：
        
        $$
        \sum_{j=1}^{k_n} \log(1+\theta_{nj}) \rightarrow \theta
        $$
        
        两边同时取复指数映射，引理得证。$\square$

---

## 2. 李雅普诺夫中心极限定理 (Liapounov's CLT)

如果随机变量序列拥有高于二阶的矩，我们可以给出一个非常容易验证的充分条件。

!!! success "定理 3.3：李雅普诺夫 CLT (Liapounov's Theorem)"

    对于双重数组 $\{X_{nj} : 1 \le j \le k_n\}_{n \ge 1}$，定义其三阶中心绝对矩之和为 $\Gamma_n = \sum_{j=1}^{k_n} E|X_{nj} - \alpha_{nj}|^3$，假设对于每个 $n$ 该值均有限。
    如果满足**李雅普诺夫条件 (Liapounov's Condition)**：
    
    $$
    \frac{\Gamma_n}{\sigma_n^3} = \frac{1}{\sigma_n^3} \sum_{j=1}^{k_n} E|X_{nj} - \alpha_{nj}|^3 \rightarrow 0 \quad \text{as } n \to \infty
    $$
    
    那么标准化和将依分布收敛于标准正态：
    
    $$
    \frac{S_n - \alpha_n}{\sigma_n} \rightarrow N(0,1)
    $$
    
    > *(注：三阶矩可以被放宽到 $2+\delta$ 阶矩，其中 $\delta > 0$。)*

    ??? proof "定理 3.3 的严格证明（点击展开）"

        设 $\gamma_{nj} = E|X_{nj} - \alpha_{nj}|^3$。由 Liapounov 矩不等式可知：
        
        $$
        \sigma_{nj} = (E|X_{nj} - \alpha_{nj}|^2)^{1/2} \le (E|X_{nj} - \alpha_{nj}|^3)^{1/3} = \gamma_{nj}^{1/3}
        $$
        
        所以 $\sigma_{nj}^3 \le \gamma_{nj}$。进而有：
        
        $$
        \max_{1 \le j \le k_n} \sigma_{nj}^3 \le \max_{1 \le j \le k_n} \gamma_{nj} \le \sum_{j=1}^{k_n} \gamma_{nj} = \Gamma_n
        $$
        
        设 $\phi_{nj}(t)$ 是标准化变量 $(X_{nj} - \alpha_{nj})/\sigma_n$ 的特征函数。因为 $\gamma_{nj}$ 有限，特征函数可进行三阶泰勒展开：
        
        $$
        \phi_{nj}(t) = 1 - \frac{\sigma_{nj}^2 t^2}{2\sigma_n^2} + \frac{\Lambda_{nj}}{6} \frac{\gamma_{nj} t^3}{\sigma_n^3}, \quad \text{其中 } |\Lambda_{nj}| < 1
        $$
        
        为了应用引理 3.2，我们令 $\theta_{nj} = \phi_{nj}(t) - 1$，并验证它的三个条件：
        
        **验证 (i)**：
        
        $$
        \max_j |\phi_{nj}(t) - 1| \le \frac{t^2}{2\sigma_n^2} \max_j \sigma_{nj}^2 + \frac{|t|^3}{6\sigma_n^3} \max_j \gamma_{nj}
        $$
        
        由于 $\sigma_{nj}^2 = (\sigma_{nj}^3)^{2/3} \le (\max_j \sigma_{nj}^3)^{2/3}$，我们有：
        
        $$
        \frac{\max \sigma_{nj}^2}{\sigma_n^2} \le \left( \frac{\max \sigma_{nj}^3}{\sigma_n^3} \right)^{2/3} \le \left( \frac{\Gamma_n}{\sigma_n^3} \right)^{2/3} \rightarrow 0
        $$
        
        同时 $\max_j \gamma_{nj} / \sigma_n^3 \le \Gamma_n / \sigma_n^3 \rightarrow 0$。因此 $\max_j |\theta_{nj}| \to 0$，条件 (i) 成立。
        
        **验证 (ii)**：
        
        $$
        \sum_{j=1}^{k_n} |\phi_{nj}(t) - 1| \le \frac{t^2}{2\sigma_n^2} \sum_{j=1}^{k_n} \sigma_{nj}^2 + \frac{|t|^3}{6} \frac{\Gamma_n}{\sigma_n^3} = \frac{t^2}{2} + \frac{|t|^3}{6} \frac{\Gamma_n}{\sigma_n^3} \le M(t)
        $$
        
        此为有界量，条件 (ii) 成立。
        
        **验证 (iii)**：
        
        由于误差项之和满足：
        
        $$
        \left| \sum_{j=1}^{k_n} \frac{\Lambda_{nj} \gamma_{nj}}{\sigma_n^3} \right| \le \frac{\Gamma_n}{\sigma_n^3} \rightarrow 0
        $$
        
        因此特征函数偏移量之和的极限为：
        
        $$
        \sum_{j=1}^{k_n} (\phi_{nj}(t) - 1) = -\frac{t^2}{2} + t^3 \sum_{j=1}^{k_n} \frac{\Lambda_{nj} \gamma_{nj}}{\sigma_n^3} \rightarrow -\frac{t^2}{2}
        $$
        
        条件 (iii) 成立。
        
        综上所述，根据引理 3.2，独立标准化变量之和的特征函数满足：
        
        $$
        \prod_{j=1}^{k_n} \phi_{nj}(t) = \prod_{j=1}^{k_n} (1 + (\phi_{nj}(t) - 1)) \rightarrow e^{-t^2/2}
        $$
        
        根据 Lévy-Cramér 连续性定理，结果得证。$\square$

对于单一下标的一般序列，该推论同样适用：

!!! info "推论 3.4 (单序列)"

    设 $\{X_n\}_{n \ge 1}$ 是一列独立的随机向量序列。令 $\alpha_j = E(X_j)$，$\sigma_j^2 = Var(X_j)$ 且 $\gamma_j = E|X_j - \alpha_j|^3 < \infty$。
    令 $P_n = \sum_{j=1}^n \gamma_j$。如果 $P_n / \sigma_n^3 \rightarrow 0$，则：
    
    $$
    \frac{S_n - \sum_{j=1}^n \alpha_j}{\sigma_n} \rightarrow N(0,1)
    $$

---

## 3. 林德伯格嵌套法 (Lindeberg's Telescoping Method)

在不使用特征函数的情况下，直接运用解析技巧（Stein方法的前身）也可以证明 CLT。

假设 $\alpha_j = 0$。引入一列服从标准正态分布的辅助随机变量 $Y_1, \dots, Y_n$，使得它们与 $X_j$ 相互独立，且 $Y_j \sim N(0, \sigma_j^2)$ 匹配前两阶矩。
令 $Y_0 = \sum_{i=1}^n Y_i / \sigma_n \sim N(0,1)$。

我们目标是证明对所有有界且具备各阶有界导数的测试函数 $f \in C_B^\infty$：

$$
E\left[f\left(\frac{\sum_{i=1}^n X_i}{\sigma_n}\right)\right] - E[f(Y_0)] \rightarrow 0
$$

需要借助如下定理：

!!! success "定理 3.5 (Chung 6.1.6)"

    设 $\{\mu_n\}$ 为一列概率测度。如果对于所有无穷次可微且各阶导数均有界的测试函数 $f \in C_B^\infty$，满足：

    $$
    \int f(x) \mu_n(dx) \rightarrow \int f(x) \mu(dx)
    $$

    其中 $\mu$ 是一个概率测度，那么 $\mu_n$ 弱收敛于 $\mu$ 。

根据定理 3.5 (Chung 6.1.6)，如果上述期望对所有的测试函数收敛，则概率测度弱收敛。

!!! success "基于望远镜展开法的证明核心"

    构造混合部分和序列 $Z_j$：
    $Z_j = Y_1 + \dots + Y_{j-1} + X_{j+1} + \dots + X_n$ （对于 $2 \le j \le n-1$）。
    边界为 $Z_1 = X_2 + \dots + X_n$ 且 $Z_n = Y_1 + \dots + Y_{n-1}$。

    我们将全量差异写成了逐项替换的嵌套和 (Telescoping sum)：
    
    $$
    E\left[f\left(\frac{\sum X_i}{\sigma_n}\right)\right] - E\left[f\left(\frac{\sum Y_i}{\sigma_n}\right)\right] = \sum_{i=1}^n \left[ E\left[f\left(\frac{Z_i + X_i}{\sigma_n}\right)\right] - E\left[f\left(\frac{Z_i + Y_i}{\sigma_n}\right)\right] \right]
    $$
    
    对 $f$ 在 $Z_i/\sigma_n$ 处进行三阶泰勒展开：
    
    $$
    f\left(\frac{Z_i + X_i}{\sigma_n}\right) = f\left(\frac{Z_i}{\sigma_n}\right) + f'\left(\frac{Z_i}{\sigma_n}\right)\frac{X_i}{\sigma_n} + \frac{1}{2}f''\left(\frac{Z_i}{\sigma_n}\right)\frac{X_i^2}{\sigma_n^2} + \theta_i^{(1)}\frac{X_i^3}{3!\sigma_n^3}
    $$
    
    由于 $X_i$ 和 $Y_i$ 前两阶矩完全相同，取期望后一阶项与二阶项完美抵消。
    只剩下三阶误差项，受三阶导数上界 $M$ 限制：
    
    $$
    \le \sum_{i=1}^n \frac{M}{3!\sigma_n^3} (\gamma_i + \sqrt{\frac{8}{\pi}}\sigma_i^3) \le \frac{M_1}{6}\frac{\Gamma_n}{\sigma_n^3} \rightarrow 0
    $$
    
    此方法为 Stein 方法和高维随机向量的高斯逼近提供了基石。

推论 3.6 截断方法预备定理：
如果存在双重数组 $|X_{nj}/\sigma_n| \le M_{nj}$ a.s. 且 $\lim_{n \to \infty} \max_j M_{nj} = 0$，则标准化和收敛于正态分布。

---

## 4. 零阵列 (Null Arrays)

为了探讨 CLT 成立的充分必要条件，我们需要排除某个单一变量在总体方差中占据统治地位的病态情况。

!!! abstract "定义 3.7：零阵列"

    如果一个双重数组满足：对于任意 $\epsilon > 0$，
    
    $$
    \lim_{n \to \infty} \max_{1 \le j \le k_n} P(|X_{nj} - \alpha_{nj}| > \epsilon \sigma_n) = 0
    $$
    
    那么称该双重数组为**零阵列 (Null array)**。
    这等价于表示每个分量 $(X_{nj} - \alpha_j)/\sigma_n$ 在 $n \to \infty$ 时对于 $j$ 是一致依概率退化到 0 的。

利用特征函数，我们可以给出一个非常便于计算的等价刻画：

!!! info "命题 3.8：零阵列的特征函数等价形式"

    双重数组 $\{X_{nj}\}$ 是零阵列，当且仅当对于所有的 $t \in \mathbb{R}$：
    
    $$
    \lim_{n \to \infty} \max_{1 \le j \le k_n} |\phi_{nj}(t) - 1| = 0
    $$
    
    并且此收敛在任意有限区间 $[-K, K]$ 上是一致的。

    ??? proof "等价性证明（点击展开）"

        **($\Rightarrow$ 方向)**：不失一般性设 $\alpha_j = 0$。
        将期望分解为在阈值 $\epsilon\sigma_n$ 内和外的两部分：
        
        $$
        |\phi_{nj}(t) - 1| \le E\left[ |e^{itX_{nj}/\sigma_n} - 1| \mathbb{I}(|X_{nj}| > \epsilon\sigma_n) \right] + E\left[ |e^{itX_{nj}/\sigma_n} - 1| \mathbb{I}(|X_{nj}| \le \epsilon\sigma_n) \right]
        $$
        
        利用不等式 $|e^{itu} - 1| \le |tu|$ 和复指数模的界限 2：
        
        $$
        \le 2 P(|X_{nj}| > \epsilon\sigma_n) + |t| \epsilon
        $$
        
        因为是在 $[-K, K]$ 的有界闭集上：
        
        $$
        \sup_{|t|\le K} \max_j |\phi_{nj}(t) - 1| \le 2 \max_j P(|X_{nj}| > \epsilon\sigma_n) + K\epsilon
        $$
        
        当 $n \to \infty$ 第一项由零阵列定义趋于 0，第二项由于 $\epsilon$ 任意也可任意小，故一致收敛得证。
        
        **($\Leftarrow$ 方向)**：
        利用经典的特征函数不等式：
        
        $$
        P\left(\left|\frac{X_{nj}}{\sigma_n}\right| > \frac{2}{\delta}\right) \le \frac{1}{\delta} \int_{|t|\le\delta} (1 - \phi_{nj}(t)) dt \le \frac{1}{\delta} \int_{|t|\le\delta} |1 - \phi_{nj}(t)| dt
        $$
        
        两边取 max 后：
        
        $$
        \max_j P\left(\left|\frac{X_{nj}}{\sigma_n}\right| > \frac{2}{\delta}\right) \le \max_j \frac{1}{\delta} \int_{|t|\le\delta} |1 - \phi_{nj}(t)| dt
        $$
        
        由积分的有界收敛定理 (BCT) 及已知条件，上式极限趋于 0。命题得证 $\square$。

---

## 5. 林德伯格-费勒中心极限定理 (Lindeberg-Feller CLT)

当我们只知道二阶矩存在，而三阶矩不一定存在时，李雅普诺夫条件失效。此时**林德伯格条件 (Lindeberg Condition)** 成为了独立变量 CLT 成立的最精确的充要条件。

!!! abstract "定义 3.9：林德伯格条件 (LC)"

    对于双重数组 $\{X_{nj}\}$，如果对于任意 $\epsilon > 0$，其截尾方差比率满足：
    
    $$
    \lim_{n \to \infty} \frac{1}{\sigma_n^2} \sum_{j=1}^{k_n} E \left[ (X_{nj} - \alpha_{nj})^2 \mathbb{I}(|X_{nj} - \alpha_{nj}| > \epsilon \sigma_n) \right] = 0
    $$
    
    则称该数组满足林德伯格条件。

这结合下面的引理引出了独立变量渐近分布理论的著名定理：

!!! info "引理 3.10 (对角线构造法)"

    设 $u(m, n)$ 是定义在正整数 $m$ 和 $n$ 上的函数，满足对于每一个固定的 $m$，有：

    $$
    \lim_{n \to \infty} u(m, n) = 0
    $$

    那么存在一个单调递增并趋于无穷的序列 $\{m_n\}$ ($m_n \to \infty$)，使得：

    $$
    \lim_{n \to \infty} u(m_n, n) = 0
    $$

    ??? proof "引理 3.10 的证明（点击展开）"

        由于对于每个固定的 $m$，$\lim_{n \to \infty} u(m, n) = 0$，因此根据极限定义，必然存在一个下标 $n_m$，使得对于所有的 $n \ge n_m$，都有：

        $$
        u(m, n_m) \le \frac{1}{m}
        $$

        通过这种方式，我们可以构造出一个严格单调递增并趋于无穷的序列 $\{n_m\}_{m \ge 1}$。

        对于任何满足 $n_m \le n < n_{m+1}$ 的 $n$，我们令 $m_n \equiv m$。

        由此，当 $n \ge n_m$ 时，根据构造有：

        $$
        u(m_n, n) = u(m, n) \le \frac{1}{m}
        $$

        由于当 $n \to \infty$ 时，$n_m$ 所对应的索引 $m \to \infty$，因此 $m_n$ 也单调递增趋于无穷。
        根据上式夹逼可知 ：

        $$
        \lim_{n \to \infty} u(m_n, n) = 0
        $$

        证明完毕。$\square$ 

!!! success "定理 3.11：林德伯格-费勒 CLT (Lindeberg-Feller)"

    假设 $Var(X_{nj}) = \sigma_{nj}^2 < \infty$，$S_n = \sum_{j=1}^{k_n} X_{nj}$。那么以下两组命题等价：
    
    * **(i)** $\frac{S_n - E S_n}{\sigma_n} \rightarrow N(0,1)$ 且 **(ii)** 该双重数组是零阵列。
    * **$\Longleftrightarrow$** 该双重数组满足**林德伯格条件 (LC)**。

    ??? proof "定理核心证明推导（点击展开）"

        **(1) 充分性证明：LC $\Rightarrow$ CLT 且零阵列**
        
        假设 $E(X_{nj})=0$ 且 $\sigma_n^2 = 1$。定义基于截断点 $\eta \in (0, 1)$ 的截断随机变量：
        $X_{nj}' = X_{nj}$ (若 $|X_{nj}| < \eta$)，否则为 0。
        
        计算截断后的期望与方差：
        
        $$
        |E[X_{nj}']| = \left| \int_{|x|<\eta} x dF_{nj}(x) \right| = \left| -\int_{|x|\ge\eta} x dF_{nj}(x) \right| \le \frac{1}{\eta} \int_{|x|\ge\eta} x^2 dF_{nj}(x)
        $$
        
        求和后，由于 LC 条件，总期望趋于 0。同时，截断后的方差 $\sigma_n'^2 \to 1 = \sigma_n^2$。
        根据引理 3.10（对角线构造法则），我们可以选取一个单调递增趋于无穷的序列 $m_n \to \infty$，并令 $\eta_n = m_n^{-1} \to 0$。
        利用 $\eta_n$ 作为截断阈值，这保证了 $|X_{nj}'| \le \eta_n := M_{nj}$。由于 $\max M_{nj} = \eta_n \to 0$，根据前面推论 3.6 (有界变量 CLT)，我们有 $(S_n' - ES_n')/\sigma_n' \to N(0,1)$，故 $S_n' \to N(0,1)$。
        
        最后评估未截断和与截断和的差异：
        
        $$
        P(S_n \neq S_n') \le \sum_{j=1}^{k_n} P(|X_{nj}| \ge \eta_n) \le \sum_{j=1}^{k_n} \frac{1}{\eta_n^2} \int_{|x| > \eta_n} x^2 dF_{nj}(x) \rightarrow 0 \quad \text{(由 LC 保证)}
        $$
        
        由 Slutsky 定理，由于误差为 $o_p(1)$，最终可得 $S_n \to N(0,1)$。
        
        **(2) 必要性证明：CLT + 零阵列 $\Rightarrow$ LC**
        
        由 $S_n \xrightarrow{d} N(0,1)$ 可知连乘特征函数的对数：
        
        $$
        \lim_{n \to \infty} \sum_{j=1}^{k_n} \log \phi_{nj}(t) = -\frac{t^2}{2}
        $$
        
        零阵列保证了 $\max_j |\phi_{nj}(t) - 1| \to 0$。利用 $\log \phi_{nj}(t)$ 与 $\phi_{nj}(t) - 1$ 的等价关系：
        
        $$
        \lim_{n \to \infty} \sum_{j=1}^{k_n} (\phi_{nj}(t) - 1) = -\frac{t^2}{2}
        $$
        
        提取其实部：
        
        $$
        \lim_{n \to \infty} \sum_j \int_{-\infty}^{\infty} (1 - \cos tx) dF_{nj}(x) = \frac{t^2}{2}
        $$
        
        将积分拆分为 $|x| \le \eta$ 和 $|x| > \eta$ 的两部分，并利用不等式 $0 \le 1 - \cos tx \le (tx)^2/2$：
        
        $$
        \frac{t^2}{2} - \sum_j \int_{|x|\le\eta} (1 - \cos tx) dF_{nj}(x) \ge \frac{t^2}{2} \sum_j \int_{|x| \ge \eta} x^2 dF_{nj}(x)
        $$
        
        通过界定积分外部部分受切比雪夫不等式的控制（$\le 2/\eta^2 + \epsilon$），令 $t \to \infty$ 取极限，最后逼出 $\sum_{j=1}^{k_n} E[X_{nj}^2 \mathbb{I}(|X_{nj}| > \eta)] \to 0$，这恰好是林德伯格条件 LC。$\square$

---

## 6. 应用与判定扩展 (Applications & Further Conditions)

!!! tip "应用示例：普通最小二乘回归 (OLS Regression)"

    考虑经典线性回归模型 $y_j = x_j \beta + \epsilon_j$，其中误差项 $\epsilon_j \sim (0, \sigma_\epsilon^2)$ 且 i.i.d.。
    设计矩阵满足 $\max_{1 \le j \le n} \frac{|x_j|}{a_n} \to 0$，其中 $a_n^2 = \sum_{j=1}^n x_j^2$。
    OLS 估计量为 $\hat{\beta}_{LS} = \sum x_j y_j / a_n^2$。
    
    构造标准化双重数组 $X_{nj} = \frac{x_j \epsilon_j}{\sqrt{\sum x_j^2}}$。
    我们将林德伯格条件中的截断积分进行放缩：令 $m_n = \max_j |x_j/a_n|$，
    
    $$
    \frac{1}{\sigma_n^2 a_n^2} \sum_{j=1}^n x_j^2 E\left[ \epsilon_j^2 \mathbb{I}(|x_j \epsilon_j / a_n| > \delta \sigma_n) \right] \le \frac{1}{\sigma_\epsilon^2} E\left[ \epsilon_j^2 \mathbb{I}(|\epsilon_j| > \delta \sigma_\epsilon m_n^{-1}) \right]
    $$
    
    由于 $\epsilon_j$ 同分布且二阶矩有限，当 $m_n \to 0$ 时，该期望趋于 0。由 Theorem 3.11 立刻得到：
    
    $$
    a_n(\hat{\beta}_{LS} - \beta) \xrightarrow{d} N(0, \sigma_\epsilon^2)
    $$

验证 LC 条件时，除了利用有界截断，我们还可以利用高阶矩的存在性，这就是**李雅普诺夫类型条件的扩展**。

!!! info "命题 3.12：林德伯格条件的充分判定"

    对于双重数组 $\{X_{nj}\}$，如果存在某个实数 $\nu > 2$，满足：
    
    $$
    \sum_{j=1}^{k_n} E|X_{nj} - \mu_{nj}|^\nu = o(\sigma_n^\nu)
    $$
    
    那么，该数组必然满足林德伯格条件。

    ??? proof "证明推导（点击展开）"

        在截尾区域 $|t - \mu_{nj}| > \epsilon \sigma_n$ 内，我们对其二阶矩积分进行放缩：
        
        $$
        E\left[ (X_{nj} - \mu_{nj})^2 \mathbb{I}(|X_{nj} - \mu_{nj}| > \epsilon \sigma_n) \right] = \int_{|t - \mu_{nj}| > \epsilon \sigma_n} (t - \mu_{nj})^2 dF_{nj}(t)
        $$
        
        通过强行凑出 $\nu$ 次方，提出常数项：
        
        $$
        \le (\epsilon \sigma_n)^{2-\nu} \int_{|t - \mu_{nj}| > \epsilon \sigma_n} (t - \mu_{nj})^\nu dF_{nj}(t) \le (\epsilon \sigma_n)^{2-\nu} E|X_{nj} - \mu_{nj}|^\nu
        $$
        
        求和并除以 $\sigma_n^2$：
        
        $$
        \frac{1}{\sigma_n^2} \sum_{j=1}^{k_n} E\left[ \dots \right] \le \epsilon^{2-\nu} \frac{\sum_{j=1}^{k_n} E|X_{nj} - \mu_{nj}|^\nu}{\sigma_n^\nu} \rightarrow 0
        $$
        
        林德伯格条件得证。$\square$