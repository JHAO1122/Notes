# 2.1 可测函数的概念与性质

在测度空间 $(X, \mathcal{M})$ 中，可测函数是研究积分理论的核心对象。正如连续函数保持了拓扑空间中的开集结构，可测函数则保持了可测空间之间的 $\sigma$-代数结构。

---

## 1. 映射与逆像的性质

在定义可测函数之前，我们需要回顾集合逆像的基本性质。逆像映射 $f^{-1}: 2^Y \to 2^X$ 与集合的各种运算是完全相容的。

!!! info "引理 2.1.1 (逆像的性质)"

    设 $f: X \to Y$ 是一个映射，$E, E_{\lambda} \subset Y$，则：

    * **并集保持**：

    \[
    f^{-1}\left(\bigcup_{\lambda} E_{\lambda}\right) = \bigcup_{\lambda} f^{-1}(E_{\lambda})
    \]

    * **交集保持**：

    \[
    f^{-1}\left(\bigcap_{\lambda} E_{\lambda}\right) = \bigcap_{\lambda} f^{-1}(E_{\lambda})
    \]

    * **补集保持**：

    \[
    f^{-1}(E^c) = (f^{-1}(E))^c
    \]

由此性质可知，若 $\mathcal{N}$ 是 $Y$ 上的 $\sigma$-代数，则集族 $\mathcal{F} = \{E \subset Y : f^{-1}(E) \in \mathcal{M}\}$ 构成 $Y$ 上的一个 $\sigma$-代数。

---

## 2. 可测函数的定义与基本判定

!!! info "定义 2.1.1 ($(M, N)$-可测性)"

    设 $(X, \mathcal{M})$ 和 $(Y, \mathcal{N})$ 是两个可测空间。称映射 $f: X \to Y$ 是 **$(M, N)$-可测的**（或简称为可测的），如果对于任何 $E \in \mathcal{N}$，其逆像都满足：

    \[
    f^{-1}(E) \in \mathcal{M}
    \]

当目标空间是实数集 $\mathbb{R}$ 且配备 Borel $\sigma$-代数 $\mathcal{B}_{\mathbb{R}}$ 时，我们称 $f$ 为 **$\mathcal{M}$-可测函数**。

* **Borel 可测函数**：若 $f: (\mathbb{R}, \mathcal{B}_{\mathbb{R}}) \to (\mathbb{R}, \mathcal{B}_{\mathbb{R}})$ 是可测的。

* **Lebesgue 可测函数**：若 $f: (\mathbb{R}, \mathcal{L}) \to (\mathbb{R}, \mathcal{B}_{\mathbb{R}})$ 是可测的。

!!! success "命题 2.1.1 (生成元判定法)"

    若目标空间的 $\sigma$-代数 $\mathcal{N}$ 由集族 $\Sigma$ 生成（即 $\mathcal{N} = \mathcal{M}(\Sigma)$），则 $f$ 是 $(M, N)$-可测的，当且仅当：

    \[
    \forall E \in \Sigma, \quad f^{-1}(E) \in \mathcal{M}
    \]

??? proof "命题 2.1.1 的证明"

    **必要性**：如果 $f$ 是可测的，由定义，对所有 $E \in \mathcal{N}$ 都有 $f^{-1}(E) \in \mathcal{M}$。因为 $\Sigma \subset \mathcal{N}$，自然对所有 $E \in \Sigma$ 成立。

    **充分性**：构造集族 $\mathcal{F} = \{E \subset Y : f^{-1}(E) \in \mathcal{M}\}$。由引理 2.1.1 知，因为逆像运算与集合的并、交、补运算完全交换，所以 $\mathcal{F}$ 是 $Y$ 上的一个 $\sigma$-代数。
    根据条件，$\Sigma \subset \mathcal{F}$。由于 $\mathcal{N}$ 是包含 $\Sigma$ 的最小 $\sigma$-代数，必然有 $\mathcal{N} = \mathcal{M}(\Sigma) \subset \mathcal{F}$。这说明对任意 $E \in \mathcal{N}$，都有 $f^{-1}(E) \in \mathcal{M}$，即 $f$ 是可测的。

!!! success "推论 (连续映射的可测性)"

    若 $X, Y$ 为拓扑空间，$f: X \to Y$ 是连续映射，则 $f$ 是 $(\mathcal{B}_X, \mathcal{B}_Y)$-可测的（即 Borel 可测的）。

??? proof "推论的证明"

    拓扑空间 $Y$ 的 Borel $\sigma$-代数 $\mathcal{B}_Y$ 是由 $Y$ 中所有的开集族 $\mathcal{T}_Y$ 生成的。
    任取开集 $V \in \mathcal{T}_Y$，由于 $f$ 是连续映射，其逆像 $f^{-1}(V)$ 必为 $X$ 中的开集。
    而 $X$ 中的开集显然属于其 Borel $\sigma$-代数 $\mathcal{B}_X$。根据命题 2.1.1，连续函数必然是 Borel 可测的。

---

## 3. 实可测函数的等价条件

对于定义在实数轴上的函数，由于 Borel $\sigma$-代数可以由各种形式的射线生成，我们有以下常用的等价判定准则。

!!! success "命题 2.1.2 (实可测函数的等价条件)"

    对于映射 $f: X \to \mathbb{R}$，下列命题是等价的：

    * (a) $f$ 是可测的（即 $f^{-1}(B) \in \mathcal{M}, \forall B \in \mathcal{B}_{\mathbb{R}}$）。

    * (b) 对于任何 $a \in \mathbb{R}$，$f^{-1}((a, +\infty)) = \{x \in X : f(x) > a\} \in \mathcal{M}$。

    * (c) 对于任何 $a \in \mathbb{R}$，$f^{-1}([a, +\infty)) = \{x \in X : f(x) \ge a\} \in \mathcal{M}$。

    * (d) 对于任何 $a \in \mathbb{R}$，$f^{-1}((-\infty, a)) = \{x \in X : f(x) < a\} \in \mathcal{M}$。

    * (e) 对于任何 $a \in \mathbb{R}$，$f^{-1}((-\infty, a]) = \{x \in X : f(x) \le a\} \in \mathcal{M}$。

??? proof "命题 2.1.2 的证明"

    因为 $\mathbb{R}$ 上的 Borel $\sigma$-代数 $\mathcal{B}_{\mathbb{R}}$ 是由开射线族 $\Sigma = \{(a, +\infty) : a \in \mathbb{R}\}$ 生成的，根据命题 2.1.1，(a) 与 (b) 等价。
    接下来只需证明 (b), (c), (d), (e) 之间可以通过 $\sigma$-代数的封闭性相互推导：

    * **(b) $\Rightarrow$ (c)**：由于 $[a, +\infty) = \bigcap_{n=1}^\infty (a - \frac{1}{n}, +\infty)$，两边取逆像得：

    \[
    f^{-1}([a, +\infty)) = \bigcap_{n=1}^\infty f^{-1}\left(\left(a - \frac{1}{n}, +\infty\right)\right)
    \]

    由 (b) 知等号右边每个集合都在 $\mathcal{M}$ 中，根据 $\sigma$-代数对可数交封闭，得出左边集合也在 $\mathcal{M}$ 中。

    * **(c) $\Rightarrow$ (d)**：由于 $(-\infty, a) = ([a, +\infty))^c$，两边取逆像得：

    \[
    f^{-1}((-\infty, a)) = \left(f^{-1}([a, +\infty))\right)^c
    \]

    由 (c) 知 $f^{-1}([a, +\infty)) \in \mathcal{M}$，再由 $\sigma$-代数对补集封闭，得出原像属于 $\mathcal{M}$。

    * **(d) $\Rightarrow$ (e)**：类似于前面，利用 $(-\infty, a] = \bigcap_{n=1}^\infty (-\infty, a + \frac{1}{n})$ 取逆像即可。

    * **(e) $\Rightarrow$ (b)**：利用 $(a, +\infty) = ((-\infty, a])^c$ 取逆像和补集即可完成循环论证。

---

## 4. 局部可测性与复合、积空间

!!! info "定义 2.1.2 (局部可测性)"

    设 $E \in \mathcal{M}$ 为一个可测集。称函数 $f$ 在 $E$ 上是可测的，如果对于任意 $B \in \mathcal{B}_{\mathbb{R}}$，有：

    \[
    \{x \in E : f(x) \in B\} \in \mathcal{M}
    \]

    这等价于将全空间限制为 $E$，即 $f|_E$ 在子 $\sigma$-代数 $\mathcal{M}_E = \{A \cap E : A \in \mathcal{M}\}$ 意义下是可测的。

!!! success "定理 2.1.2 (复合性质)"

    若 $f: X \to Y$ 是 $(M, N)$-可测的，$g: Y \to Z$ 是 $(N, P)$-可测的，则复合映射 $g \circ f: X \to Z$ 是 $(M, P)$-可测的。

??? proof "定理 2.1.2 的证明"

    对于任意 $E \in P$，由于 $g$ 是 $(N, P)$-可测的，有 $g^{-1}(E) \in N$。
    又因为 $f$ 是 $(M, N)$-可测的，故 $f^{-1}(g^{-1}(E)) \in M$。
    显然 $(g \circ f)^{-1}(E) = f^{-1}(g^{-1}(E)) \in M$，得证。特别地，可测函数与连续函数的复合依然是可测的。

对于积空间 $Y = \prod_{\alpha \in A} Y_{\alpha}$ 及其上的乘积 $\sigma$-代数 $\mathcal{N} = \bigotimes_{\alpha \in A} \mathcal{N}_{\alpha}$，设 $\pi_\alpha : Y \to Y_\alpha$ 为投影映射。

!!! success "命题 2.1.3 (积空间的可测性)"

    映射 $f: X \to Y$ 是 $(M, N)$-可测的，当且仅当其每个坐标映射 $f_{\alpha} = \pi_{\alpha} \circ f$ 都是 $(M, N_{\alpha})$-可测的。

??? proof "命题 2.1.3 的证明"

    **必要性 ($\Rightarrow$)**：若 $f$ 可测，由于投影映射 $\pi_\alpha$ 总是 $(N, N_\alpha)$-可测的，根据复合性质，$f_\alpha = \pi_\alpha \circ f$ 必定是 $(M, N_\alpha)$-可测的。

    **充分性 ($\Leftarrow$)**：乘积 $\sigma$-代数 $\mathcal{N}$ 是由集族 $\Sigma = \{\pi_\alpha^{-1}(E_\alpha) : E_\alpha \in N_\alpha, \alpha \in A\}$ 生成的。
    对于 $\Sigma$ 中的任意元素，其逆像为：

    \[
    f^{-1}(\pi_\alpha^{-1}(E_\alpha)) = ( \pi_\alpha \circ f )^{-1}(E_\alpha) = f_\alpha^{-1}(E_\alpha)
    \]

    由假设 $f_\alpha$ 可测，所以上式属于 $\mathcal{M}$。再根据生成元判定法（命题 2.1.1），$f$ 是可测的。
