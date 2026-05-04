# 习题二 (Homework 03 & 04)

本部分包含点集拓扑课程第三、第四次作业相关练习。

---

## 问题 1

设 $S$ 为一个非空集合。

(i) 考虑离散度量 $d: S \times S \to \mathbb{R}_{\ge 0}$，定义为：当 $x=y$ 时 $d(x, y)=0$；当 $x \neq y$ 时 $d(x, y)=1$。证明该度量诱导的拓扑 $\mathcal{T}_d$ 正好是幂集 $\mathcal{P}(S)$，即 $S$ 的每一个子集在度量空间 $(S, d)$ 中都是开集。

(ii) 是否总能找到 $S$ 上的一个度量，使得其诱导的拓扑是平庸拓扑（即 $\{ \emptyset, S \}$）？

??? success "解答"

    **(i) 证明离散度量诱导离散拓扑**

    要证明 $\mathcal{T}_d = \mathcal{P}(S)$，只需证明 $S$ 中的每一个单点集 $\{x\}$ 都是开集即可。

    对于任意 $x \in S$，考虑以 $x$ 为中心，半径为 $\epsilon = 1/2$ 的开球：

    $$
    B(x, 1/2) = \{ y \in S \mid d(y, x) < 1/2 \}
    $$

    根据离散度量的定义，当 $y \neq x$ 时 $d(y, x) = 1 > 1/2$；只有当 $y = x$ 时 $d(y, x) = 0 < 1/2$。
    因此，$B(x, 1/2) = \{x\}$。

    由于度量空间中的开球必为开集，故 $\{x\}$ 是开集。
    因为 $S$ 的任意子集 $A$ 都可以表示为单点集的并集：

    $$
    A = \bigcup_{x \in A} \{x\}
    $$

    且拓扑对任意并运算封闭，所以 $A$ 必定是开集。故 $\mathcal{T}_d = \mathcal{P}(S)$。

    **(ii) 关于平庸拓扑的可度量化讨论**

    **结论**：当 $|S| \le 1$ 时可以；当 $|S| \ge 2$ 时不存在。

    **理由**：
    假设 $S$ 至少包含两个不同的点 $x$ 和 $y$。若存在度量 $d$ 诱导平庸拓扑，则 $d(x, y) = \epsilon > 0$。
    考虑开球 $B(x, \epsilon) = \{ z \in S \mid d(z, x) < \epsilon \}$。

    由于 $d(x, x) = 0 < \epsilon$，所以 $x \in B(x, \epsilon)$。
    由于 $d(y, x) = \epsilon \not< \epsilon$，所以 $y \notin B(x, \epsilon)$。

    这意味着 $B(x, \epsilon)$ 既不是空集（包含 $x$），也不是全空间 $S$（不包含 $y$）。
    但在平庸拓扑中，唯一的开集只有 $\emptyset$ 和 $S$。产生矛盾。
    因此，若 $|S| \ge 2$，不存在诱导平庸拓扑的度量。

---

## 问题 2

证明 $\mathcal{B} := \{ (a, b) \mid a < b \text{ 且 } a, b \in \mathbb{Q} \}$ 是 $\mathbb{R}$ 上标准拓扑（由度量 $d(x, y) = |x - y|$ 诱导）的一个基。并证明 $\text{Card}(\mathcal{B}) = \aleph_0$。

??? success "解答"

    **1. 证明 $\mathcal{B}$ 是标准拓扑的基**

    设 $\mathcal{T}$ 为 $\mathbb{R}$ 上的标准拓扑。要证明 $\mathcal{B}$ 是基，需证明对于任意开集 $U \in \mathcal{T}$ 及任意点 $x \in U$，都存在 $B \in \mathcal{B}$ 使得 $x \in B \subseteq U$。

    已知 $U$ 是开集，故存在 $\epsilon > 0$ 使得 $(x - \epsilon, x + \epsilon) \subseteq U$。
    根据有理数在实数中的**稠密性 (Density)**：
    在区间 $(x - \epsilon, x)$ 中必存在有理数 $a$；
    在区间 $(x, x + \epsilon)$ 中必存在有理数 $b$。

    由此我们得到 $a < x < b$，且 $a, b \in \mathbb{Q}$。
    显然 $(a, b) \in \mathcal{B}$。
    并且因为 $x - \epsilon < a$ 且 $b < x + \epsilon$，故：

    $$
    x \in (a, b) \subseteq (x - \epsilon, x + \epsilon) \subseteq U
    $$

    证毕。

    **2. 计算基的势**

    集合 $\mathcal{B}$ 的每一个元素由有序有理数对 $(a, b)$ 唯一确定，且满足 $a < b$。
    因此存在一个单射 $\Psi: \mathcal{B} \to \mathbb{Q} \times \mathbb{Q}$。
    由于 $\mathbb{Q}$ 是可数集，其笛卡尔积 $\mathbb{Q} \times \mathbb{Q}$ 也是可数集（即势为 $\aleph_0$）。
    因为 $\mathcal{B}$ 是无穷集（例如包含所有 $(0, n), n \in \mathbb{Q}^+$），且其势不大于 $\aleph_0$，由 Cantor-Bernstein 定理：

    $$
    \text{Card}(\mathcal{B}) = \aleph_0
    $$

---

## 问题 3

对于 $n \in \mathbb{N}_{\ge 1}$，考虑大小为 $n$ 的有限集 $X_n := \{1, 2, \dots, n\}$。
$X_n$ 上有多少个拓扑？此外，考虑 $X_n$ 上所有拓扑构成的集合 $\Omega$。定义 $\mathcal{T}_i < \mathcal{T}_j$ 为 $\mathcal{T}_i \subsetneq \mathcal{T}_j$。问 $(\Omega, <)$ 在何时构成线性序关系？

??? success "解答"

    **1. 拓扑的数量**

    有限集上的拓扑数量没有简单的通项公式，通常记为 $T(n)$。对于前几个 $n$，其数值为：

    * $n=1: T(1) = 1$（仅平庸拓扑）
    * $n=2: T(2) = 4$（平庸、离散、以及两个谢尔宾斯基拓扑）
    * $n=3: T(3) = 29$

    **2. 线性序关系的讨论**

    **结论**：$(\Omega, <)$ 构成线性序关系当且仅当 $n=1$。

    **理由**：
    线性序要求对于 $\Omega$ 中的任意两个元素 $\mathcal{T}_i, \mathcal{T}_j$，必有 $\mathcal{T}_i \subseteq \mathcal{T}_j$ 或 $\mathcal{T}_j \subseteq \mathcal{T}_i$。

    当 $n \ge 2$ 时，我们可以构造两个不可比的拓扑。设 $X_n$ 包含不同元素 $a$ 和 $b$。
    定义拓扑 $\mathcal{T}_1 = \{ \emptyset, X_n, \{a\} \}$；
    定义拓扑 $\mathcal{T}_2 = \{ \emptyset, X_n, \{b\} \}$。

    显然 $\{a\} \notin \mathcal{T}_2$ 且 $\{b\} \notin \mathcal{T}_1$。
    故 $\mathcal{T}_1 \not\subseteq \mathcal{T}_2$ 且 $\mathcal{T}_2 \not\subseteq \mathcal{T}_1$。
    因此，当 $n \ge 2$ 时，$(\Omega, <)$ 不是线性序。

---

## 问题 4

设 $S$ 为一个集合。

(i) 证明：若 $\{ \mathcal{T}_i \}_{i \in I}$ 是 $S$ 上的一族拓扑，则其交集 $\bigcap_{i \in I} \mathcal{T}_i$ 也是 $S$ 上的一个拓扑。
(ii) 给出两个拓扑 $\mathcal{T}_1$ 和 $\mathcal{T}_2$ 的例子，使得 $\mathcal{T}_1 \cup \mathcal{T}_2$ 不是拓扑。并找到包含 $\mathcal{T}_1 \cup \mathcal{T}_2$ 的唯一最小拓扑。

??? success "解答"

    **(i) 证明交集是拓扑**

    记 $\mathcal{T} = \bigcap_{i \in I} \mathcal{T}_i$。

    1. **包含空集与全空间**：因为对所有 $i \in I$ 都有 $\emptyset, S \in \mathcal{T}_i$，故 $\emptyset, S \in \mathcal{T}$。
    2. **对任意并封闭**：设 $\{ U_\alpha \}_{\alpha \in A} \subseteq \mathcal{T}$。则对每一个 $i \in I$，都有 $\{ U_\alpha \} \subseteq \mathcal{T}_i$。由于 $\mathcal{T}_i$ 是拓扑，故 $\bigcup U_\alpha \in \mathcal{T}_i$。因此 $\bigcup U_\alpha \in \bigcap \mathcal{T}_i = \mathcal{T}$。
    3. **对有限交封闭**：设 $U_1, \dots, U_k \in \mathcal{T}$。则对每一个 $i \in I$，都有 $U_1, \dots, U_k \in \mathcal{T}_i$，故 $\bigcap_{j=1}^k U_j \in \mathcal{T}_i$。因此其有限交属于 $\mathcal{T}$。

    **(ii) 并集不是拓扑的反例**

    设 $S = \{a, b, c\}$。
    取 $\mathcal{T}_1 = \{ \emptyset, S, \{a\} \}$，$\mathcal{T}_2 = \{ \emptyset, S, \{b\} \}$。
    则 $\mathcal{T}_1 \cup \mathcal{T}_2 = \{ \emptyset, S, \{a\}, \{b\} \}$。

    然而，$\{a\} \cup \{b\} = \{a, b\}$，但 $\{a, b\} \notin \mathcal{T}_1 \cup \mathcal{T}_2$。
    违反了对并运算封闭的公理，故并集不是拓扑。

    **最小拓扑**：
    包含 $\mathcal{T}_1 \cup \mathcal{T}_2$ 的唯一最小拓扑是由该并集作为**子基 (Subbasis)** 生成的拓扑。在这个例子中，该拓扑为：

    $$
    \mathcal{T}_{min} = \{ \emptyset, S, \{a\}, \{b\}, \{a, b\} \}
    $$

---

## 问题 5

$\mathbb{R}^2$ 的子集若满足：它在其每一个点处的每一个方向上都包含一条开线段，则称其为**径向开集** (Radially open)。记 $\mathcal{T}_{ro}$ 为所有径向开集构成的族。

证明 $\mathcal{T}_{ro}$ 是 $\mathbb{R}^2$ 上的一个拓扑。它与标准拓扑 $\mathcal{T}_{std}$ 的关系如何？

??? success "解答"

    **1. 证明是拓扑**

    重点在于证明对有限交封闭。
    设 $A, B \in \mathcal{T}_{ro}$。对于任何点 $x \in A \cap B$：

    * 因为 $x \in A$，在方向 $\theta$ 上存在开线段 $I_A(\theta) \subseteq A$；
    * 因为 $x \in B$，在方向 $\theta$ 上存在开线段 $I_B(\theta) \subseteq B$。

    取这两个线段的交集 $I_A \cap I_B$，由于它们共中心且都在同一射线上，其交集仍是一个包含 $x$ 的开线段，且 $I_A \cap I_B \subseteq A \cap B$。
    因此 $A \cap B$ 也是径向开的。结合 $\emptyset, S$ 的显见性和并运算的性质，$\mathcal{T}_{ro}$ 是拓扑。

    **2. 与标准拓扑的关系**

    **结论**：$\mathcal{T}_{std} \subsetneq \mathcal{T}_{ro}$。

    * **包含关系**：标准开集中的点包含一个开球 $B(x, \epsilon)$。该球在任何方向上显然都包含长度为 $2\epsilon$ 的开线段。故标准开集必为径向开集。
    * **严格性（反例）**：存在径向开集不是标准开集。
      构造集合：$U = \{ (r\cos\theta, r\sin\theta) \mid 0 \le r < f(\theta) \}$，其中 $f(\theta)$ 是一个在极角趋于某些值时迅速趋于 0 的正函数。
      例如，若使 $f(\theta)$ 在原点附近极其“尖锐”，使得原点虽在每个方向都有线段，但无法放入任何完整的圆盘。这种集合在原点处是径向开的，但在标准拓扑下原点不是内点。

---

## 问题 6

设 $X = \mathbb{R}^2$ 赋予标准拓扑。考虑其子集（拓扑学家正弦曲线）：

$$
A = \{ (x, \sin(1/x)) : x \in (0, 1] \} \subset X
$$

计算以下四个集合：$A^o$，$LP(A)$，$\overline{A}$ 以及 $\overline{A} \cap \overline{\mathbb{R}^2 - A}$。

??? success "解答"

    **1. 内部 $A^o$**

    由于 $A$ 是连续函数 $y = \sin(1/x)$ 在区间 $(0, 1]$ 上的图像，它在 $\mathbb{R}^2$ 中是一条曲线。
    对于 $A$ 中的任意一点 $P$，其任何开邻域（开圆盘）都会包含不属于该曲线的点（即纵坐标不满足 $y = \sin(1/x)$ 的点）。
    因此，$A$ 中不包含任何非空的开球。
    故：

    $$
    A^o = \emptyset
    $$

    **2. 极限点集 $LP(A)$ 与 闭包 $\overline{A}$**

    当 $x \to 0^+$ 时，函数 $\sin(1/x)$ 在 $-1$ 到 $1$ 之间进行无限次振荡。
    因此，纵轴上的线段 $\{0\} \times [-1, 1]$ 中的每一个点都是 $A$ 的极限点。
    集合 $A$ 自身的点显然也是其极限点。
    因此：

    $$
    LP(A) = A \cup ( \{0\} \times [-1, 1] )
    $$

    由于闭包 $\overline{A} = A \cup LP(A)$，我们可以得出：

    $$
    \overline{A} = A \cup ( \{0\} \times [-1, 1] )
    $$

    **3. 交集 $\overline{A} \cap \overline{\mathbb{R}^2 - A}$**

    因为 $A^o = \emptyset$，所以全空间中的每一个点要么属于 $\mathbb{R}^2 - A$ 的内部，要么属于 $A$ 的边界。
    由于 $\overline{\mathbb{R}^2 - A} = \mathbb{R}^2 - A^o = \mathbb{R}^2 - \emptyset = \mathbb{R}^2$。
    因此：

    $$
    \overline{A} \cap \overline{\mathbb{R}^2 - A} = \overline{A} \cap \mathbb{R}^2 = \overline{A}
    $$

---

## 问题 7

设 $X$ 为一个集合，$f: \mathcal{P}(X) \to \mathcal{P}(X)$ 是一个满足以下性质的函数（对任意 $A, B \subseteq X$）：
(i) $f(\emptyset) = \emptyset$；
(ii) $A \subseteq f(A)$；
(iii) $f(A) = f(f(A))$；
(iv) $f(A \cup B) = f(A) \cup f(B)$。

证明：
(1) 集合 $\mathcal{T}_f := \{ X - A : A \in \mathcal{P}(X) \text{ 且 } f(A) = A \}$ 是 $X$ 上的一个拓扑；
(2) 任意子集 $A \subseteq X$ 关于拓扑 $\mathcal{T}_f$ 的闭包 $\overline{A}$ 正好是 $f(A)$。

??? success "解答"

    **(1) 证明 $\mathcal{T}_f$ 是一个拓扑**

    * **空集与全空间**：
      由性质 (i) 知 $f(\emptyset) = \emptyset$，故 $X - \emptyset = X \in \mathcal{T}_f$。
      由性质 (ii) 知 $X \subseteq f(X)$，而 $f(X) \subseteq X$ 显然成立，故 $f(X) = X$。
      因此 $X - X = \emptyset \in \mathcal{T}_f$。

    * **任意并运算**：
      设 $\{ U_\alpha \}_{\alpha \in I}$ 是 $\mathcal{T}_f$ 中的一族开集，其中 $U_\alpha = X - A_\alpha$，且 $f(A_\alpha) = A_\alpha$。
      其并集 $\bigcup U_\alpha = X - \bigcap A_\alpha$。
      由性质 (iv) 可推导出 $f$ 是单调的：若 $A \subseteq B$，则 $f(A) \subseteq f(B)$。
      由于 $\bigcap A_\alpha \subseteq A_\alpha$，故 $f(\bigcap A_\alpha) \subseteq f(A_\alpha) = A_\alpha$。
      从而 $f(\bigcap A_\alpha) \subseteq \bigcap A_\alpha$。
      结合性质 (ii) 的反向包含，得 $f(\bigcap A_\alpha) = \bigcap A_\alpha$。
      故 $\bigcup U_\alpha \in \mathcal{T}_f$。

    * **有限交运算**：
      设 $U_1, U_2 \in \mathcal{T}_f$，则 $U_1 \cap U_2 = X - (A_1 \cup A_2)$。
      由性质 (iv) 知 $f(A_1 \cup A_2) = f(A_1) \cup f(A_2) = A_1 \cup A_2$。
      故 $U_1 \cap U_2 \in \mathcal{T}_f$。

    **(2) 证明闭包 $\overline{A} = f(A)$**

    根据闭包的定义，$\overline{A}$ 是包含 $A$ 的所有闭集的交集。在此拓扑下，闭集 $C$ 满足 $f(C) = C$。

    * 首先，由性质 (iii) 知 $f(f(A)) = f(A)$，故 $f(A)$ 本身是一个闭集。且由性质 (ii) 知 $A \subseteq f(A)$。
      因此 $f(A)$ 是一个包含 $A$ 的闭集，故 $\overline{A} \subseteq f(A)$。

    * 其次，设 $C$ 是任何包含 $A$ 的闭集（即 $A \subseteq C$ 且 $f(C) = C$）。
      由 $f$ 的单调性知 $f(A) \subseteq f(C) = C$。
      这意味着 $f(A)$ 包含在每一个包含 $A$ 的闭集中，因此 $f(A) \subseteq \bigcap C = \overline{A}$。

    综上所述，$\overline{A} = f(A)$。

---

## 问题 8

考虑函数 $f: [0, 1) \to S^1$，定义为 $f(t) = (\cos(2\pi t), \sin(2\pi t))$。
证明 $f$ 关于标准子空间拓扑不是同胚。是否可能找到 $[0, 1)$ 与 $S^1$ 之间的同胚映射？

??? success "解答"

    **1. 证明 $f$ 不是同胚**

    虽然 $f$ 是一个连续的双射，但它的逆映射 $f^{-1}$ 在点 $(1, 0)$ 处不连续。
    我们可以从拓扑不变量的角度给出更严谨的证明：

    * **方法一：紧致性**
      在标准欧氏拓扑下，$S^1$ 是 $\mathbb{R}^2$ 中的有界闭集，因此是紧致的。
      而 $[0, 1)$ 在 $\mathbb{R}$ 中不是闭集，因此不是紧致的。
      由于同胚映射必须保持紧致性，故 $f$ 不可能是同胚。

    * **方法二：割点（连通性）**
      若 $f$ 是同胚，则去掉一个点后空间的连通性应该保持不变。
      对于 $S^1$，去掉任意一点 $p$，剩余空间 $S^1 - \{p\}$ 仍然是连通的（同胚于一个开区间）。
      对于 $[0, 1)$，如果我们去掉点 $t = 1/2$，剩余空间 $[0, 1/2) \cup (1/2, 1)$ 是不连通的。
      因此，这两个空间不具有相同的局部连通结构，故 $f$ 不是同胚。

    **2. 结论**

    **不可能**。
    基于上述紧致性和割点的理由，$[0, 1)$ 与 $S^1$ 之间不存在任何同胚映射。

---

## 问题 9

子集 $A \subseteq X$ 若满足 $A = (\overline{A})^o$，则称其为**正则开集** (Regularly open)。

(i) 给出 $\mathbb{R}$ 中一个开集但不是正则开集的例子。你能刻画 $\mathbb{R}$ 中的正则开集吗？
(ii) 证明对于任何 $A \subseteq X$，集合 $(\overline{A})^o$ 都是正则开集。

??? success "解答"

    **(i) 例子与刻画**

    * **反例**：取 $A = (0, 1) \cup (1, 2)$。
      则 $\overline{A} = [0, 2]$，而 $(\overline{A})^o = (0, 2)$。
      显然 $A \neq (\overline{A})^o$，故 $A$ 不是正则开集。其原因在于 $A$ 中缺失了闭包内部的一个点 $\{1\}$。

    * **刻画**：$\mathbb{R}$ 中的正则开集是那些“没有多余缝隙”的开集。即它不能通过从一个更大的开集中挖去一些孤立点或闭集界点而得到。

    **(ii) 证明 $(\overline{A})^o$ 是正则开集**

    设 $U = (\overline{A})^o$。我们要证明 $U = (\overline{U})^o$。

    * 根据定义 $U \subseteq \overline{U}$，由于 $U$ 是开集，故 $U = U^o \subseteq (\overline{U})^o$。

    * 另一方面，由于 $U = (\overline{A})^o \subseteq \overline{A}$，且闭包运算具有单调性和幂等性：
      由此推得 $\overline{U} \subseteq \overline{(\overline{A})} = \overline{A}$。
      对两边取内部，得 $(\overline{U})^o \subseteq (\overline{A})^o = U$。

    综上所述，$U = (\overline{U})^o$，即 $(\overline{A})^o$ 是正则开集。

---

## 问题 10

设 $A$ 是拓扑空间 $X$ 的子集。证明等式 $\overline{A} = A \cup LP(A)$。

??? success "解答"

    **1. 证明 $A \cup LP(A) \subseteq \overline{A}$**

    * 显见 $A \subseteq \overline{A}$。
    
    * 设 $x \in LP(A)$。根据极限点的定义，对于 $x$ 的任何开邻域 $U$，都有 $U \cap (A - \{x\}) \neq \emptyset$。
      这意味着 $U \cap A \neq \emptyset$ 对所有开邻域成立。
      根据闭包的邻域刻画性质，这正是 $x \in \overline{A}$ 的充分必要条件。
      故 $LP(A) \subseteq \overline{A}$。
    
    因此 $A \cup LP(A) \subseteq \overline{A}$。

    **2. 证明 $\overline{A} \subseteq A \cup LP(A)$**

    设 $x \in \overline{A}$。

    * 若 $x \in A$，则结论显然成立。
    
    * 若 $x \notin A$，由于 $x \in \overline{A}$，对于 $x$ 的任何开邻域 $U$，都有 $U \cap A \neq \emptyset$。
      因为 $x \notin A$，集合 $U \cap A$ 中的点必然不是 $x$。
      也就是说，$U \cap (A - \{x\}) \neq \emptyset$。
      这符合极限点的定义，故 $x \in LP(A)$。
    
    因此 $x \in A \cup LP(A)$。

    综上所述，$\overline{A} = A \cup LP(A)$。