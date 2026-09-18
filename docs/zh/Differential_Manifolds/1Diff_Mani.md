# 第一章：微分流形

微分流形的基本思想是：一个空间在整体上可以具有复杂的拓扑结构，但在每个点附近都应当像某个 Euclidean 空间 $\mathbb{R}^n$。我们利用局部坐标把流形上的问题转化为 $\mathbb{R}^n$ 中的问题，再通过坐标变换保证不同局部描述彼此一致。

本章依次介绍拓扑流形、坐标卡、转移映射、$C^r$ 相容性、微分结构和极大图册，并以 $\mathbb{R}^n$ 与球面 $S^n$ 为基本例子。

---

## 1. 从坐标系到流形

在平面 $\mathbb{R}^2$ 中，一个点可以用直角坐标 $(x,y)$ 表示，也可以在适当区域内用极坐标 $(r,\theta)$ 表示。极坐标在原点处失效，角变量还需要选取分支，这表明一个坐标系通常只能定义在局部区域上。

类似地，我们会问：球面 $S^2$ 上能否存在一个覆盖整个球面的坐标系？答案是否定的。若存在整体坐标卡，则 $S^2$ 将与 $\mathbb{R}^2$ 的某个开子集同胚；但 $S^2$ 是紧致的，而 $\mathbb{R}^2$ 的非空开子集不可能紧致。因此，研究流形必须使用多个局部坐标卡。

!!! note "局部坐标的基本思想"

    对流形上的每个点 $p$，选取它的一个邻域 $U$，并用同胚映射

    \[
    \varphi:U\longrightarrow\varphi(U)\subseteq\mathbb{R}^n
    \]

    将 $U$ 中的点表示为 $n$ 个实数坐标。

---

## 2. 拓扑流形

!!! info "定义 2.1（拓扑流形）"

    设 $M$ 是一个 Hausdorff 拓扑空间。如果对任意 $p\in M$，都存在 $p$ 的一个开邻域 $U\subseteq M$ 以及一个同胚映射

    \[
    \varphi:U\longrightarrow\varphi(U),
    \]

    其中 $\varphi(U)$ 是 $\mathbb{R}^n$ 的开子集，则称 $M$ 是一个 **$n$ 维拓扑流形（topological manifold）**。

    本课程还约定流形具有可数拓扑基，即满足第二可数公理。

拓扑流形也可以概括为：**局部 Euclidean 的 Hausdorff 空间**。这里“局部 Euclidean”是指每个点都具有与 $\mathbb{R}^n$ 中某个开集同胚的邻域。

!!! note "流形的维数是良定义的"

    同一个非空开集不可能既与 $\mathbb{R}^m$ 的开集同胚，又与 $\mathbb{R}^n$ 的开集同胚而满足 $m\neq n$。这一事实可由拓扑中的区域不变性定理（invariance of domain）推出。因此，拓扑流形的维数 $n$ 不依赖于局部坐标的选择。

### 2.1 基本例子与反例

!!! example "例 1（Euclidean 空间）"

    $\mathbb{R}^n$ 本身是 $n$ 维拓扑流形。对每个点都可以直接取 $U=\mathbb{R}^n$，坐标映射为恒等映射

    \[
    \operatorname{id}_{\mathbb{R}^n}:\mathbb{R}^n\longrightarrow\mathbb{R}^n.
    \]

!!! example "例 2（球面与环面）"

    球面 $S^n$ 是 $n$ 维拓扑流形，但通常不能由一个坐标卡覆盖。环面

    \[
    T^2=S^1\times S^1
    \]

    也是二维拓扑流形。

!!! warning "反例（双锥面的锥点）"

    双锥面除去锥点后局部与 $\mathbb{R}^2$ 同胚，因此是二维拓扑流形；但锥点本身没有 Euclidean 邻域。

    直观上，删去锥点后，它的一个充分小邻域会分成两个连通分支；而从 $\mathbb{R}^2$ 的一个小圆盘中删去中心点后，所得集合仍然连通。因此，锥点附近不可能与平面开集同胚。

---

## 3. 坐标卡与局部坐标

!!! info "定义 3.1（坐标卡）"

    设 $M$ 是 $n$ 维拓扑流形。一个**坐标卡（coordinate chart）** 是一对 $(U,\varphi)$，其中 $U\subseteq M$ 是开集，而

    \[
    \varphi:U\longrightarrow\varphi(U)\subseteq\mathbb{R}^n
    \]

    是从 $U$ 到 $\mathbb{R}^n$ 中开集 $\varphi(U)$ 的同胚。

    集合 $U$ 称为坐标邻域，映射 $\varphi$ 称为坐标映射。

将坐标映射写成

\[
\varphi(p)=\bigl(x^1(p),\ldots,x^n(p)\bigr),
\]

则每个分量

\[
x^i:U\longrightarrow\mathbb{R},
\qquad
p\longmapsto x^i(p),
\qquad i=1,\ldots,n,
\]

称为坐标卡 $(U,\varphi)$ 的第 $i$ 个**局部坐标函数**。因此，$(x^1,\ldots,x^n)$ 也称为 $U$ 上的一个局部坐标系。

---

## 4. 转移映射

设 $(U,\varphi)$ 与 $(V,\psi)$ 是同一个 $n$ 维拓扑流形 $M$ 上的两个坐标卡，并假设

\[
U\cap V\neq\varnothing.
\]

在重叠区域中，同一个点既有 $x$-坐标，也有 $y$-坐标。由此得到**转移映射（transition map）**

\[
\psi\circ\varphi^{-1}:
\varphi(U\cap V)\longrightarrow\psi(U\cap V).
\]

若写成

\[
\varphi(p)=(x^1,\ldots,x^n),
\qquad
\psi(p)=(y^1,\ldots,y^n),
\]

则转移映射的分量可以表示为

\[
y^i=y^i(x^1,\ldots,x^n),
\qquad i=1,\ldots,n.
\]

反方向的坐标变换为

\[
\varphi\circ\psi^{-1}:
\psi(U\cap V)\longrightarrow\varphi(U\cap V),
\]

其分量形式为

\[
x^i=x^i(y^1,\ldots,y^n),
\qquad i=1,\ldots,n.
\]

---

## 5. 坐标卡的 $C^r$ 相容性

!!! info "定义 5.1（$C^r$ 相容）"

    设 $(U,\varphi)$ 与 $(V,\psi)$ 是 $M$ 上的两个坐标卡。

    * 若 $U\cap V=\varnothing$，约定这两个坐标卡 $C^r$ 相容；
    * 若 $U\cap V\neq\varnothing$，并且转移映射

    \[
    \psi\circ\varphi^{-1}:
    \varphi(U\cap V)\longrightarrow\psi(U\cap V)
    \]

    是一个 $C^r$ 微分同胚，则称 $(U,\varphi)$ 与 $(V,\psi)$ **$C^r$ 相容**。

因为

\[
(\psi\circ\varphi^{-1})^{-1}
=\varphi\circ\psi^{-1},
\]

所以要求转移映射为 $C^r$ 微分同胚，等价于要求两个方向的坐标变换都属于 $C^r$。

!!! warning "$C^r$ 相容一般不是等价关系"

    坐标卡的 $C^r$ 相容关系具有反身性和对称性，但在任意坐标卡的集合上不一定具有传递性。

    例如，若第二个坐标卡分别与第一个、第三个坐标卡不相交，则前两对坐标卡自动相容；然而第一个与第三个坐标卡可能存在重叠，并且它们之间的转移映射未必属于 $C^r$。因此，$C^r$ 相容关系一般不能直接视为等价关系。

---

## 6. 图册与微分结构

!!! info "定义 6.1（坐标图册）"

    $M$ 上的一族坐标卡

    \[
    \mathcal{A}=\{(U_\alpha,\varphi_\alpha):\alpha\in\Lambda\}
    \]

    称为一个**坐标图册（atlas）**，如果这些坐标邻域覆盖 $M$：

    \[
    \bigcup_{\alpha\in\Lambda}U_\alpha=M.
    \]

    若图册中的任意两个坐标卡都 $C^r$ 相容，则称 $\mathcal{A}$ 是一个 $C^r$ 图册。

!!! info "定义 6.2（$C^r$ 微分结构）"

    一个 $C^r$ 图册

    \[
    \mathcal{U}=\{(U_\alpha,\varphi_\alpha):\alpha\in\Lambda\}
    \]

    称为 $M$ 上的 **$C^r$ 微分结构**，如果它还满足极大性：对任意坐标卡 $(V,\psi)$，只要 $(V,\psi)$ 与 $\mathcal{U}$ 中的每个坐标卡都 $C^r$ 相容，就必有

    \[
    (V,\psi)\in\mathcal{U}.
    \]

    换言之，$C^r$ 微分结构就是一个极大 $C^r$ 图册。

!!! info "定义 6.3（$C^r$ 流形与光滑流形）"

    赋予了 $C^r$ 微分结构 $\mathcal{U}$ 的拓扑流形 $M$ 称为 **$C^r$ 微分流形**，记作 $(M,\mathcal{U})$。

    当 $r=\infty$ 时，$\mathcal{U}$ 称为 $M$ 上的**光滑结构（smooth structure）**，而 $(M,\mathcal{U})$ 称为**光滑流形（smooth manifold）**。

通常不需要把极大图册中的所有坐标卡逐一写出。只要给出一族彼此相容且覆盖 $M$ 的坐标卡，就可以唯一确定相应的极大图册。

!!! success "定理 6.4（相容图册的唯一极大扩张）"

    设 $M$ 是 $n$ 维拓扑流形，$\mathcal{A}$ 是一个覆盖 $M$ 的 $C^r$ 图册。则存在唯一的极大 $C^r$ 图册 $\mathcal{U}_{\max}$ 包含 $\mathcal{A}$。

??? proof "定理 6.4 的证明（点击展开）"

    定义

    \[
    \mathcal{U}_{\max}
    =\{(V,\psi):(V,\psi)\text{ 与 }\mathcal{A}
    \text{ 中每个坐标卡都 }C^r\text{ 相容}\}.
    \]

    显然 $\mathcal{A}\subseteq\mathcal{U}_{\max}$，因此 $\mathcal{U}_{\max}$ 覆盖 $M$。

    下面证明 $\mathcal{U}_{\max}$ 中任意两个坐标卡 $(V,\psi)$ 与 $(W,\chi)$ 都相容。固定 $p\in V\cap W$。由于 $\mathcal{A}$ 覆盖 $M$，存在 $(U_\alpha,\varphi_\alpha)\in\mathcal{A}$ 使 $p\in U_\alpha$。在 $p$ 的一个充分小邻域内，转移映射可以分解为

    \[
    \chi\circ\psi^{-1}
    =(\chi\circ\varphi_\alpha^{-1})
    \circ(\varphi_\alpha\circ\psi^{-1}).
    \]

    因为 $(V,\psi)$ 与 $(W,\chi)$ 分别同 $(U_\alpha,\varphi_\alpha)$ 相容，右侧两个映射都属于 $C^r$，所以 $\chi\circ\psi^{-1}$ 在 $p$ 附近属于 $C^r$。由于 $p$ 任意，而光滑性是局部性质，两个坐标卡 $C^r$ 相容。因此 $\mathcal{U}_{\max}$ 是一个 $C^r$ 图册。

    根据定义，任何与 $\mathcal{U}_{\max}$ 中所有坐标卡相容的坐标卡，尤其与 $\mathcal{A}$ 中所有坐标卡相容，因而已经属于 $\mathcal{U}_{\max}$。所以 $\mathcal{U}_{\max}$ 是极大的。

    最后，若另一个极大 $C^r$ 图册 $\mathcal{V}$ 也包含 $\mathcal{A}$，则 $\mathcal{V}$ 中每个坐标卡都与 $\mathcal{A}$ 相容，从而 $\mathcal{V}\subseteq\mathcal{U}_{\max}$。由 $\mathcal{V}$ 的极大性可得二者相等。因此极大扩张唯一。$\square$

---

## 7. 基本例子

### 7.1 $\mathbb{R}^n$ 的标准光滑结构

$\mathbb{R}^n$ 上的单个坐标卡

\[
(\mathbb{R}^n,\operatorname{id}_{\mathbb{R}^n})
\]

构成一个光滑图册。由定理 6.4，它唯一确定 $\mathbb{R}^n$ 的极大光滑图册，这称为 $\mathbb{R}^n$ 的**标准光滑结构**。

### 7.2 球面 $S^n$ 的光滑结构

将单位球面写成

\[
S^n
=\left\{(x^1,\ldots,x^{n+1})\in\mathbb{R}^{n+1}:
\sum_{j=1}^{n+1}(x^j)^2=1\right\}.
\]

球面带有从 $\mathbb{R}^{n+1}$ 继承的子空间拓扑。对 $i=1,\ldots,n+1$，定义开集

\[
U_i^+=\{x\in S^n:x^i>0\},
\qquad
U_i^-=\{x\in S^n:x^i<0\}.
\]

这些集合覆盖球面：

\[
S^n=\bigcup_{i=1}^{n+1}(U_i^+\cup U_i^-).
\]

令

\[
D^n=\left\{(u^1,\ldots,u^n)\in\mathbb{R}^n:
\sum_{j=1}^n(u^j)^2<1\right\}
\]

为 $\mathbb{R}^n$ 中的开单位球。定义坐标映射

\[
\varphi_i^\pm:U_i^\pm\longrightarrow D^n
\]

为删去第 $i$ 个坐标的投影：

\[
\varphi_i^\pm(x^1,\ldots,x^{n+1})
=(x^1,\ldots,\widehat{x^i},\ldots,x^{n+1}),
\]

其中符号 $\widehat{x^i}$ 表示省略第 $i$ 个分量。

其逆映射通过球面方程恢复第 $i$ 个坐标：

\[
(\varphi_i^\pm)^{-1}(u^1,\ldots,u^n)
=\left(u^1,\ldots,u^{i-1},
\pm\sqrt{1-\sum_{j=1}^n(u^j)^2},
u^i,\ldots,u^n\right).
\]

因此，每个 $(U_i^\pm,\varphi_i^\pm)$ 都是球面上的一个坐标卡。

!!! success "命题 7.1（球面坐标卡彼此光滑相容）"

    坐标卡族

    \[
    \mathcal{A}
    =\{(U_i^+,\varphi_i^+),(U_i^-,\varphi_i^-):
    i=1,\ldots,n+1\}
    \]

    是 $S^n$ 上的光滑图册，因而唯一确定 $S^n$ 的标准光滑结构。

??? proof "命题 7.1 的证明（点击展开）"

    已知这些坐标邻域覆盖 $S^n$，所以只需验证重叠区域上的转移映射光滑。一般情形完全类似，下面具体考察 $(U_1^+,\varphi_1^+)$ 与 $(U_2^-,\varphi_2^-)$。

    对

    \[
    u=(u^1,\ldots,u^n)\in
    \varphi_1^+(U_1^+\cap U_2^-),
    \]

    有

    \[
    (\varphi_1^+)^{-1}(u)
    =\left(\sqrt{1-\sum_{j=1}^n(u^j)^2},
    u^1,u^2,\ldots,u^n\right).
    \]

    由于交集内满足 $x^2<0$，转移映射的定义域为

    \[
    \varphi_1^+(U_1^+\cap U_2^-)
    =\{u\in D^n:u^1<0\}.
    \]

    再删去第二个坐标，得到

    \[
    \varphi_2^-\circ(\varphi_1^+)^{-1}(u^1,\ldots,u^n)
    =\left(\sqrt{1-\sum_{j=1}^n(u^j)^2},
    u^2,\ldots,u^n\right).
    \]

    该映射在其定义域上属于 $C^\infty$。

    反过来，对

    \[
    v=(v^1,\ldots,v^n)\in
    \varphi_2^-(U_1^+\cap U_2^-),
    \]

    有

    \[
    (\varphi_2^-)^{-1}(v)
    =\left(v^1,-\sqrt{1-\sum_{j=1}^n(v^j)^2},
    v^2,\ldots,v^n\right),
    \]

    其中定义域为

    \[
    \varphi_2^-(U_1^+\cap U_2^-)
    =\{v\in D^n:v^1>0\}.
    \]

    因此

    \[
    \varphi_1^+\circ(\varphi_2^-)^{-1}(v^1,\ldots,v^n)
    =\left(-\sqrt{1-\sum_{j=1}^n(v^j)^2},
    v^2,\ldots,v^n\right),
    \]

    也属于 $C^\infty$。其余坐标卡之间的转移映射具有完全相同的形式：重排若干坐标，并用

    \[
    \pm\sqrt{1-\sum_j(u^j)^2}
    \]

    恢复被省略的坐标。因此所有转移映射都是光滑的，$\mathcal{A}$ 是 $S^n$ 上的光滑图册。$\square$

---

## 8. 本章小结

1. 拓扑流形是局部与 $\mathbb{R}^n$ 的开集同胚的 Hausdorff 空间；本课程还要求第二可数性。
2. 坐标卡 $(U,\varphi)$ 为流形上的点提供局部 Euclidean 坐标。
3. 不同坐标卡通过转移映射 $\psi\circ\varphi^{-1}$ 联系。
4. 若所有转移映射均为 $C^r$ 微分同胚，则坐标图册给出一个 $C^r$ 微分结构。
5. 任意覆盖 $M$ 的相容 $C^r$ 图册都唯一确定一个极大 $C^r$ 图册。
6. $\mathbb{R}^n$ 与 $S^n$ 都具有自然的标准光滑结构。
