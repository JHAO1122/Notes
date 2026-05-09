# 2.2 可测函数序列的极限与分解

## 5. 复可测函数与广义实值函数

对于复值函数 $f: X \to \mathbb{C}$，我们利用复平面作为 $\mathbb{R}^2$ 的积空间拓扑结构，即 $\mathcal{B}_{\mathbb{C}} = \mathcal{B}_{\mathbb{R}} \otimes \mathcal{B}_{\mathbb{R}}$。

!!! success "命题 2.1.4 (复可测性)"

    复值函数 $f$ 是可测的，当且仅当其实部 $\text{Re } f$ 和虚部 $\text{Im } f$ 都是实可测函数。

在测度论中，我们经常允许函数取值 $\pm \infty$。定义广义实数集 $\bar{\mathbb{R}} = \mathbb{R} \cup \{-\infty, +\infty\}$，其 Borel $\sigma$-代数定义为：

\[
\mathcal{B}_{\bar{\mathbb{R}}} = \{E \subset \bar{\mathbb{R}} : E \cap \mathbb{R} \in \mathcal{B}_{\mathbb{R}}\}
\]

若 $f$ 为广义实可测函数，我们通常约定 $0 \cdot \infty = 0$ 且避免无意义的 $\infty - \infty$。

---

## 6. 代数运算与序列的极限

可测函数在各种运算和取极限操作下具有完美的封闭性。

!!! success "定理 2.1.3 (代数运算的封闭性)"

    设 $f, g: X \to \mathbb{C}$ 是可测函数，则 $f+g$ 和 $f \cdot g$ 都是可测函数。

??? proof "定理 2.1.3 的证明"

    我们定义映射 $F(x) = (f(x), g(x)): X \to \mathbb{C} \times \mathbb{C}$。
    由于 $f, g$ 均可测，根据积空间可测性（命题 2.1.3），联合映射 $F$ 是可测的。
    又因为复数加法映射 $S: \mathbb{C} \times \mathbb{C} \to \mathbb{C}$，即 $S(z, w) = z + w$ 是连续映射，因此是 Borel 可测的。
    根据复合性质，$f+g = S \circ F$ 必定是可测的。同理可证乘法 $f \cdot g$ 的可测性。

!!! success "定理 2.1.4 (函数列的上/下确界与极限)"

    设 $\{f_j\}_{j=1}^\infty$ 为一列定义在 $X$ 上的广义实值可测函数，则以下函数也是可测的：

    * $g_1(x) = \sup_{j} f_j(x)$
    * $g_2(x) = \inf_{j} f_j(x)$
    * $g_3(x) = \limsup_{j \to \infty} f_j(x)$
    * $g_4(x) = \liminf_{j \to \infty} f_j(x)$

??? proof "定理 2.1.4 的证明"

    我们仅证明上确界 $\sup_j f_j$ 的可测性，其余可由代数关系推导。
    对于任意 $a \in \mathbb{R}$，观察集合 $\{x : \sup_j f_j(x) > a\}$。
    上确界大于 $a$，当且仅当数列中至少存在某一项大于 $a$。因此：

    \[
    \{x : \sup_j f_j(x) > a\} = \bigcup_{j=1}^\infty \{x : f_j(x) > a\}
    \]

    由于每个 $f_j$ 都是可测的，等式右侧的每一个集合 $\{x : f_j(x) > a\}$ 都在 $\mathcal{M}$ 中。因为 $\mathcal{M}$ 对可数并封闭，所以上确界函数 $g_1$ 满足可测的等价条件（命题 2.1.2(b)），故可测。

    其余证明：
    $\inf_j f_j = -\sup_j (-f_j)$，由代数运算与 $\sup$ 的可测性即得。
    由于 $\limsup_{j \to \infty} f_j = \inf_{k \ge 1} (\sup_{j \ge k} f_j)$，经过一次可列 $\sup$ 和一次可列 $\inf$ 运算，自然保持可测性。$\liminf$ 同理。

!!! info "推论 (正负部与极分解)"

    对于任何实可测函数 $f$，我们可以将其分解为正部与负部：

    * **正部**：$f^+ = \max(f, 0)$
    * **负部**：$f^- = \max(-f, 0)$
    * 此时有 $f = f^+ - f^-$ 且 $|f| = f^+ + f^-$。这两个函数同样是可测的。

    对于复可测函数，我们可以进行极分解（符号分解）：
    
    \[
    f = (\text{sgn } f) |f|
    \]
    
    其中 $\text{sgn } z = z/|z|$ 当 $z \ne 0$，否则为 0。$\text{sgn } f$ 和 $|f|$ 同样是可测函数。