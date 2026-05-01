# 第八章：正交分解与 Hilbert 空间理论

本章将深入研究 Hilbert 空间的核心拓扑性质，重点讨论正交分解定理、规范正交系、Fourier 展开以及可分 Hilbert 空间与 $l^2$ 空间的等距同构关系。

---

## 1. 正交分解定理

!!! success "定理 8.1 (正交分解定理)"

    设 $M$ 是 Hilbert 空间 $H$ 中的闭子空间，则 $H$ 中的任一元素 $x$，有下列唯一的正交分解：

    $$
    x = y + z, \quad y \in M, \quad z \in M^\perp
    $$

    $y$ 称为 $x$ 在 $M$ 中的正交投影。

??? proof "证明：正交分解的存在性与唯一性"

    **1. 存在性：**
    由假设 $M$ 是 $H$ 的闭子空间，故为凸闭集。 
    由上一章定理，$x$ 在 $M$ 中存在唯一的最佳逼近元 $y$。记 $\|x - y\| = \alpha$。

    * 对于实的 Hilbert 空间：由于 $y \in M$，于是对任一实数 $\lambda$ 以及任一元素 $u \in M$，有 $y + \lambda u \in M$。
    
    * 因此，函数 $\|-x + y + \lambda u\|^2$ 在 $\lambda = 0$ 时取到极小值。
    
    * 对 $\lambda$ 求导：

    $$
    0 = \frac{d}{d\lambda} \| -x + y + \lambda u \|^2 \Big|_{\lambda=0}
    $$

    $$
    = \frac{d}{d\lambda} \left[ \|y - x\|^2 + \lambda(y - x, u) + \lambda(u, y - x) + \lambda^2\|u\|^2 \right] \Big|_{\lambda=0}
    $$

    $$
    = (y - x, u) + (u, y - x) = 2(y - x, u)
    $$
    
    * 由于 $u$ 是 $M$ 中任一元素，故 $z = x - y \perp M$。便有 $x = y + z$，$y \in M$，$z \in M^\perp$。
    
    * 对于复的 Hilbert 空间：取 $\lambda = -\frac{(y - x, u)}{\|u\|^2}$，并注意到 $\|x - y\| = \alpha$ 得到：

    $$
    \alpha^2 \le \|-x + y + \lambda u\|^2
    $$

    $$
    = \|y - x\|^2 + \overline{\lambda}(y - x, u) + \lambda(u, y - x) + |\lambda|^2\|u\|^2
    $$

    $$
    = \|y - x\|^2 + \lambda(u, y - x) + \overline{\lambda}[(y - x, u) + \lambda\|u\|^2]
    $$

    * 将 $\lambda$ 代入化简得到：

    $$
    \alpha^2 \le \alpha^2 - \frac{|(y - x, u)|^2}{\|u\|^2}
    $$

    * 于是必然有 $(y - x, u) = 0$。 从而 $z = x - y \in M^\perp$。

    **2. 唯一性：**
    
    * 设另有分解 $x = y' + z'$，其中 $y' \in M$，$z' \in M^\perp$。
    
    * 由 $y + z = y' + z'$，可得 $y - y' = z' - z$。
    
    * 由于 $y - y' \in M$，$z' - z \in M^\perp$，且 $M \cap M^\perp = \{\theta\}$（若 $x \in M \cap M^\perp$，则 $(x, x) = 0$ 从而 $x = \theta$），
    
    * 故 $y - y' = z' - z = \theta$。因此 $y = y'$，$z = z'$。正交分解唯一。 $\square$

---

## 2. 规范正交系与傅里叶系数

!!! info "定义 8.1 (规范正交系)"

    设 $\{e_n\}$ 为内积空间中的元素系，满足：

    $$
    (e_n, e_m) = \begin{cases} 0, & n \neq m \\ 1, & n = m \end{cases}
    $$

    则称 $\{e_n\}$ 是 $H$ 中的一个**规范正交系**。 
    
    对任一元素 $x \in H$，称 $c_n = (x, e_n)$ 为 $x$ 关于 $\{e_n\}$ 的第 $n$ 个**傅里叶 (Fourier) 系数**，简称为 Fourier 系数，而 $\{(x, e_n)\}$ 为 $x$ 的 Fourier 系数集。

### 2.1 在有限维子空间上的正交投影

!!! success "推论"

    设 $n$ 是给定的自然数，$\{e_1, \dots, e_n\}$ 是内积空间中的一个规范正交系，$M$ 是 $\{e_1, \dots, e_n\}$ 张成的子空间，则对任给的 $x \in H$，$x$ 在 $M$ 上的正交投影为：

    $$
    y = \sum_{k=1}^n (x, e_k)e_k
    $$

??? proof "证明"

    * 由正交分解定理，$H$ 中的任一元素 $x$ 有唯一的正交分解 $x = y + z$，$y \in M$，$z \in M^\perp$。
    
    * 则有：

    $$
    y = \sum_{k=1}^n (y, e_k)e_k = \sum_{k=1}^n (x, e_k)e_k - \sum_{k=1}^n (z, e_k)e_k = \sum_{k=1}^n (x, e_k)e_k
    $$

    * （因为 $z \in M^\perp$，所以 $(z, e_k) = 0$） $\square$

### 2.2 典型规范正交系实例

* **空间 $L^2[0, 2\pi]$**：函数系 $\{ \frac{1}{\sqrt{2\pi}} e^{int} \}$ ($n = 0, \pm 1, \pm 2, \dots$) 是一个规范正交系。
    * 证明：当 $n \neq m$ 时：

    $$
    \left( \frac{1}{\sqrt{2\pi}} e^{int}, \frac{1}{\sqrt{2\pi}} e^{imt} \right) = \frac{1}{2\pi} \int_0^{2\pi} e^{i(n-m)t} dt = \frac{1}{2\pi i(n-m)} e^{i(n-m)t} \Big|_0^{2\pi} = 0
    $$

* **空间 $l^2$**：元素系 $e_n = (0, \dots, 0, 1, 0, \dots)$ （第 $n$ 个位置为1）是 $l^2$ 的一个规范正交系。

* **空间 $L^2([-1,1]; \frac{1}{\sqrt{1-t^2}})$ 与切比雪夫 (Chebyshev) 多项式**：考察函数系 $T_n(t) = \cos(n \arccos t)$，其中 $t \in [-1,1]$。
    * 在恒等式 $\cos n\theta + \cos(n-2)\theta = 2 \cos \theta \cos(n-1)\theta$ 中令 $\theta = \arccos t$，得到递推公式 $T_n(t) = 2t T_{n-1}(t) - T_{n-2}(t)$。
    
    * 由 $T_0(t)=1, T_1(t)=t$ 以及递推公式有 $T_2(t)=2t^2-1, T_3(t)=4t^3-3t$。用数学归纳法可证 $T_n(t)$ 确为 $n$ 次多项式（第一类 Chebyshev 多项式）。
    
    * 注意到 $t = \cos \theta$。当 $n \neq m$ 时：

    $$
    \int_{-1}^1 \frac{T_n(t)T_m(t)}{\sqrt{1-t^2}} dt = \int_0^\pi \cos n\theta \cdot \cos m\theta \, d\theta = 0
    $$

    * 当 $n = m$ 时：

    $$
    \int_{-1}^1 \frac{T_n(t)^2}{\sqrt{1-t^2}} dt = \int_0^\pi \cos^2 n\theta \, d\theta = \begin{cases} \frac{\pi}{2}, & n \neq 0 \\ \pi, & n = 0 \end{cases}
    $$

    * 从而可以归一化构造出规范正交系 $\tilde{T}_n(t)$，当 $n \ge 1$ 时系数为 $\sqrt{\frac{2}{\pi}}$，当 $n=0$ 时系数为 $\sqrt{\frac{1}{\pi}}$。

---

## 3. Bessel 不等式与 Riesz-Fischer 定理

!!! success "定理 8.2 (Bessel 不等式)"

    设 $\{e_n\}$ 是内积空间中的一个规范正交系，则对任意的 $x \in H$，下式成立：

    $$
    \sum_{k=1}^\infty |(x, e_k)|^2 \le \|x\|^2
    $$

??? proof "证明"

    * 任取 $x \in H$，由前面的推论知对任给的自然数 $n$，$y = \sum_{k=1}^n (x, e_k)e_k$ 是 $x$ 在 $\{e_1, \dots, e_n\}$ 张成的子空间上的正交投影。
    
    * 于是 $x = y + z$，其中 $y \perp z$。由此立即可知 $\|y\| \le \|x\|$，即：

    $$
    \sum_{k=1}^n |(x, e_k)|^2 \le \|x\|^2
    $$

    * 令 $n \to \infty$，得到 Bessel 不等式。 $\square$

由 Bessel 不等式可知，任一元素关于规范正交系的 Fourier 系数组成的序列必属于 $l^2$。 反之有如下定理：

!!! success "定理 8.3 (Riesz-Fischer 定理)"

    设 $\{e_n\}$ 是 Hilbert 空间中的一个规范正交系，数列 $\{c_n\} \in l^2$。那么存在 $H$ 中的元素 $x$，使得 $\{c_n\}$ 是这个元素关于 $\{e_n\}$ 的 Fourier 系数集，并且下述 **帕塞瓦尔 (Parseval) 公式** 成立：

    $$
    \sum_{k=1}^\infty |(x, e_k)|^2 = \|x\|^2
    $$

??? proof "证明"

    * 令 $x_n = \sum_{k=1}^n c_k e_k$。设自然数 $m, n$ 满足 $m > n$。

    $$
    \|x_m - x_n\|^2 = \left\| \sum_{k=n+1}^m c_k e_k \right\|^2 = \sum_{k=n+1}^m |c_k|^2
    $$

    * 由于 $\{c_n\} \in l^2$，当 $m, n \to \infty$ 时，尾部级数 $\sum_{k=n+1}^m |c_k|^2 \to 0$。
    
    * 故 $\|x_m - x_n\| \to 0$。因此 $\{x_n\}$ 是 $H$ 中的基本点列。
    
    * 由于 $H$ 完备，存在元素 $x \in H$ 使得 $x_n \to x$。由内积的连续性，对任给的自然数 $k_0$，有：

    $$
    (x_n, e_{k_0}) \to (x, e_{k_0})
    $$

    * 另一方面，当 $n \ge k_0$ 时，$(x_n, e_{k_0}) = \left(\sum_{k=1}^n c_k e_k, e_{k_0}\right) = c_{k_0}$。
    
    * 于是 $(x, e_{k_0}) = c_{k_0}$。由于 $k_0$ 是任意给定的，因此 $\{c_n\}$ 就是元素 $x$ 的 Fourier 系数集。
    
    * 再令：

    $$
    \sum_{k=1}^n |(x, e_k)|^2 = \|x_n\|^2
    $$

    * 其中 $n \to \infty$ 时，Parseval 等式成立。 $\square$

---

## 4. 规范正交系的完备性与完全性

Parseval 等式未必总是对所有元素成立，因为 $\{e_n\}$ 中的元素可能“不足够多”。如果元素“足够多”，Parseval 等式就有可能对一切元素成立。

!!! info "定义 8.2 (完备的与完全的)"

    * **完备的 (Complete)**：如果对每个 $x \in H$，Parseval 等式 $\sum_{k=1}^\infty |(x, e_k)|^2 = \|x\|^2$ 恒成立，则称 $\{e_n\}$ 是完备的。
    
    * **完全的 (Total)**：如果对每个 $x \in H$，由 $(x, e_k) = 0$ $(k=1,2,\dots)$ 可以导出 $x = \theta$，则称 $\{e_n\}$ 是完全的。

!!! success "定理 8.4"

    设 $\{e_n\}$ 是 Hilbert 空间中的一个规范正交系，则下列性质等价：
    * (i) $\{e_n\}$ 是完备的；

    * (ii) 对 $H$ 中任一元素 $x$，级数依范数收敛于 $x$，即 $x = \sum_{k=1}^\infty (x, e_k)e_k$；
    
    * (iii) 对 $H$ 中任意两个元素 $x, y$，有广义 Parseval 等式：$(x, y) = \sum_{k=1}^\infty (x, e_k)\overline{(y, e_k)}$，且右端绝对收敛；
    
    * (iv) $\{e_n\}$ 是完全的。

??? proof "证明：四大性质的等价性"

    * **(i) $\implies$ (ii)**：

    $$
    \left\|x - \sum_{k=1}^n (x, e_k)e_k\right\|^2 = \|x\|^2 - \sum_{k=1}^n |(x, e_k)|^2
    $$

    * 由完备性假设即 Parseval 等式，取极限得：

    $$
    \lim_{n \to \infty} \left\|x - \sum_{k=1}^n (x, e_k)e_k\right\|^2 = 0
    $$

    * **(ii) $\implies$ (iii)**：对任意两个元素 $x, y$，令 $x_n = \sum_{k=1}^n (x, e_k)e_k$ 和 $y_n = \sum_{k=1}^n (y, e_k)e_k$。
    
    * 注意到 $\{e_n\}$ 的正交性，我们有：

    $$
    (x_n, y_n) = \sum_{k=1}^n (x, e_k)\overline{(y, e_k)}
    $$

    * 由 (ii) 知 $x_n \to x, y_n \to y$。再由内积的连续性：

    $$
    (x, y) = \lim_{n \to \infty} (x_n, y_n) = \sum_{k=1}^\infty (x, e_k)\overline{(y, e_k)}
    $$

    * 由 Bessel 不等式可知，Fourier 系数构成的序列属于 $l^2$，故上式右端的级数绝对收敛。
    
    * **(iii) $\implies$ (iv)**：设 $x \in H$ 满足 $(x, e_k) = 0$ $(k=1,2,\dots)$。由 (iii) 可得对任意的 $y \in H$，有 $(x, y) = 0$。
    
    * 取 $y = x$，得 $(x, x) = 0$。所以 $x = \theta$。
    
    * **(iv) $\implies$ (i)**：由 Bessel 不等式可知 $\{(x, e_k)\} \in l^2$。由 Riesz-Fischer 定理知存在 $y$，使得 $\{(x, e_k)\}$ 为 $y$ 的 Fourier 系数集，且 $\sum |(x, e_k)|^2 = \|y\|^2$。
    
    * 注意到 $\{(x, e_k)\}$ 也是 $x$ 的 Fourier 系数，故 $(x-y, e_k) = 0$。由 (iv) 的完全性假设，导出 $x - y = \theta$，即 $x=y$。
    
    * 因此 (i) 的 Parseval 等式成立：$\sum_{k=1}^\infty |(x, e_k)|^2 = \|x\|^2$。 $\square$

* **推论 1**：定理中使 Parseval 等式成立的元素是唯一的。
    * （证明略述：若有 $x, x'$ 均满足，通过 (ii) 的级数展开可知 $x = x'$）。

* **推论 2**：设内积空间中存在完备的规范正交系，则该空间是可分的。
    * （证明略述：由于完备，正交系张成的子空间稠密。以有理数为系数的线性组合构成可列集且稠密，故可分）。

* **推论 3**：证明规范正交系完备的等价条件是：$\{e_k\}$ 所产生的闭子空间是整个空间。
    * （证明略述：利用正交投影是有限维空间最佳逼近元的性质放缩证明距离趋于 0）。

---

## 5. 施密特正交化定理与空间同构

!!! success "定理 8.5 (Gram-Schmidt 正交化定理)"

    设 $\{x_n\}$ 是内积空间中的一个可列子集，则由 $\{x_n\}$ 必可作一个规范正交系 $\{e_n\}$，使其张成的子空间与 $\{x_n\}$ 张成的子空间相同。

??? proof "证明：正交化构造过程"

    * 不妨设 $x_1$ 是第一个不等于零的元素，令 $e_1 = \frac{x_1}{\|x_1\|}$，则 $\|e_1\| = 1$。
    
    * 设 $x_2$ 是第一个与 $e_1$ 线性无关的元素，令 $h_2 = x_2 - (x_2, e_1)e_1$。
    
    * 则 $h_2 \neq 0$（否则线性相关），且易证 $h_2 \perp e_1$。令 $e_2 = \frac{h_2}{\|h_2\|}$。
    
    * 假设已作好了相互正交且范数均为 1 的元素 $e_1, \dots, e_{k-1}$。设 $x_k$ 是第一个与它们线性无关的元素，令：

    $$
    h_k = x_k - \sum_{j=1}^{k-1} (x_k, e_j)e_j
    $$

    * 则 $h_k \neq 0$ 且 $h_k \perp e_j$。再令 $e_k = \frac{h_k}{\|h_k\|}$。
    
    * 由归纳法，得到最多可列个元素的规范正交系。
    
    * 由构造过程，每个 $e_k$ 可由 $x_1, \dots, x_k$ 线性表示，反之每个 $x_k$ 亦可由 $e_1, \dots, e_k$ 线性表示。故两者张成的子空间完全相同。 $\square$

由上述正交化定理，可以直接得到一个重要结论：**任何可分内积空间必存在完备的规范正交系。** ### 5.1 $L^2(\mathbb{R}; e^{-t^2})$ 与埃尔米特 (Hermite) 多项式

* 函数族 $\{1, t, t^2, \dots\}$ 属于加权空间 $L^2(\mathbb{R}; e^{-t^2})$。记权函数 $\omega(t) = e^{-t^2}$。
    * 对 $\omega(t)$ 求高阶导数，由归纳法可证 $\omega^{(n)}(t) = y_n(t)e^{-t^2}$，其中 $y_n(t)$ 是最高次项系数为 $(-2)^n$ 的 $n$ 次多项式。
    
    * 对任一次数小于 $n$ 的多项式 $u(t)$，多次分部积分可得：

    $$
    \int_{-\infty}^{+\infty} y_n(t)e^{-t^2} u(t) dt = \int_{-\infty}^{+\infty} \omega^{(n)}(t)u(t) dt = \dots = (-1)^n \int_{-\infty}^{+\infty} \omega(t)u^{(n)}(t) dt = 0
    $$

    * 这表明 $y_n(t)$ 与所有低次多项式正交，即它就是 $\{1, t, t^2, \dots\}$ 正交化后的产物。
    
    * 归一化：计算范数得 $\|y_n\|^2 = 2^n n! \sqrt{\pi}$。
    
    * 令 $H_n(t) = \frac{y_n(t)}{\|y_n\|}$，即得到著名的 **Hermite 多项式** 规范正交系。 它可以表示为：

    $$
    H_n(t) = \frac{1}{(2^n n! \sqrt{\pi})^{\frac{1}{2}}} e^{t^2} \frac{d^n}{dt^n}(e^{-t^2})
    $$

    * 由于多项式系稠密，这构成一个完备的规范正交系。

### 5.2 空间的等距同构

!!! success "定理 8.6"

    **每一个实（或复）可分 Hilbert 空间必与实（或复）$l^2$ 空间等距同构。** 从而所有可分 Hilbert 空间彼此必然等距同构。

??? proof "证明：同构映射的构造"

    * 设 $H$ 为可分 Hilbert 空间，$\{e_n\}$ 是完备规范正交系。定义映射 $T: H \to l^2$ 为 $Tx = \{c_n\}$，其中 $\{c_n\}$ 为 $x$ 的 Fourier 系数集。
    
    * **线性性**：$T(\alpha x + \beta y) = \{\alpha c_n + \beta c_n'\} = \alpha Tx + \beta Ty$。
    
    * **等距性**：由 Parseval 等式，$\|x\|^2 = \sum |c_n|^2 = \|Tx\|^2$，所以 $T$ 是等距同态映射。 
    
    * **单射性**：若 $Tx = \{0,0,\dots\}$，由完全性推导 $x = \theta$。
    
    * **满射性**：对任意 $\{c_n\} \in l^2$，由 Riesz-Fischer 定理存在相应的 $x \in H$，故 $T$ 是满射。 $\square$

由于无穷维空间的性质与有限维截然不同，例如在 $l^2$ 中，有界闭球不能被有限个小开球覆盖（无穷维空间中的单位闭球不是紧致的）。 

---

## 6. 本章总结

* 内积可以导出范数，从而保留欧氏空间及赋范空间的几何与拓扑性质。内积定义了正交性，这是泛函分析中极其重要的结构特征。

* 正交分解定理保证了元素可以在任何闭子空间中找到唯一的投影，这是获得规范正交系存在性及相关重要结论的基石。

* “完备性”与“完全性”在可分 Hilbert 空间中是等价的；Bessel 不等式无条件成立，而 Parseval 等式仅当规范正交系完备时才对所有元素成立。

* 核心结论：可分 Hilbert 空间中存在完备规范正交系，所有可分 Hilbert 空间彼此等距同构，这是 Hilbert 空间理论中最美的基本事实之一。