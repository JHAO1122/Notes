# 3.3 乘积测度与 Fubini 定理

在建立了一维（或单一抽象空间）的积分理论之后，我们自然需要将积分推广到多维或多个空间的乘积上（如从单变量积分过渡到重积分）。这就需要引入乘积测度的概念，并解决一个核心问题：**重积分是否可以化为累次积分计算？** 以及 **积分的顺序是否可以交换？**

这正是著名的 Tonelli 定理和 Fubini 定理所要回答的问题。在整个讨论中，我们默认涉及的测度空间都是 **$\sigma$-有限的 ( $\sigma$-finite )**。

---

## 1. 乘积可测空间与截面

为了在乘积空间 $X \times Y$ 上定义测度，我们首先需要建立其上的 $\sigma$-代数。最基本的构造单元是“可测长方形”。

!!! info "定义 3.3.1 (可测长方形与乘积 $\sigma$-代数)"

    设 $(X, \mathcal{M})$ 和 $(Y, \mathcal{N})$ 为两个可测空间。称集合 $A \times B$ 为一个**可测长方形**，如果 $A \in \mathcal{M}$ 且 $B \in \mathcal{N}$。
    
    由所有可测长方形构成的集族生成了 $X \times Y$ 上的**乘积 $\sigma$-代数**，记为：

    \[
    \mathcal{M} \otimes \mathcal{N} = \sigma(\{A \times B \mid A \in \mathcal{M}, B \in \mathcal{N}\})
    \]

在处理二元函数和二元集合时，固定其中一个变量所得到的“截面”是重要的研究对象。

!!! info "定义 3.3.2 (集合与函数的截面 Section)"

    对于 $E \subset X \times Y$，固定 $x \in X$ 或 $y \in Y$，定义其**集合截面**为：

    * $x$-截面：$E_x = \{y \in Y \mid (x, y) \in E\} \subset Y$

    * $y$-截面：$E^y = \{x \in X \mid (x, y) \in E\} \subset X$

    对于定义在 $X \times Y$ 上的函数 $f(x, y)$，定义其**函数截面**为：

    * $f_x(y) = f(x, y)$ （视为定义在 $Y$ 上的函数）

    * $f^y(x) = f(x, y)$ （视为定义在 $X$ 上的函数）

!!! success "定理 3.3.1 (截面的可测性)"

    如果 $E \in \mathcal{M} \otimes \mathcal{N}$，那么对于所有的 $x \in X$ 和 $y \in Y$，都有 $E_x \in \mathcal{N}$ 且 $E^y \in \mathcal{M}$。
    
    如果 $f$ 是 $X \times Y$ 上的 $\mathcal{M} \otimes \mathcal{N}$-可测函数，那么 $f_x$ 是 $\mathcal{N}$-可测的，$f^y$ 是 $\mathcal{M}$-可测的。

??? proof "集合截面可测性的证明"

    我们以 $E_x \in \mathcal{N}$ 的证明为例。构造集族：

    \[
    \mathcal{R} = \{E \subset X \times Y \mid E_x \in \mathcal{N}\}
    \]

    首先，如果 $E = A \times B$ 是一个可测长方形，那么当 $x \in A$ 时 $E_x = B \in \mathcal{N}$；当 $x \notin A$ 时 $E_x = \emptyset \in \mathcal{N}$。因此，所有可测长方形都在 $\mathcal{R}$ 中。

    其次，我们需要验证 $\mathcal{R}$ 是一个 $\sigma$-代数。
    
    * **对补集封闭**：$(E^c)_x = (E_x)^c$。因为 $\mathcal{N}$ 是 $\sigma$-代数，如果 $E_x \in \mathcal{N}$，则 $(E_x)^c \in \mathcal{N}$，故 $E^c \in \mathcal{R}$。
    
    * **对可数并封闭**：$(\bigcup_{j=1}^\infty E_j)_x = \bigcup_{j=1}^\infty (E_j)_x$。同样由 $\mathcal{N}$ 的封闭性可知并集也在 $\mathcal{R}$ 中。

    既然 $\mathcal{R}$ 是一个包含所有可测长方形的 $\sigma$-代数，而 $\mathcal{M} \otimes \mathcal{N}$ 是包含可测长方形的**最小** $\sigma$-代数，必然有 $\mathcal{M} \otimes \mathcal{N} \subset \mathcal{R}$。定理得证。

---

## 2. 单调类定理 (Monotone Class Theorem)

在证明乘积测度的性质时，由于 $\sigma$-代数非常复杂，我们往往先在代数上证明，然后通过“单调类定理”将其推广到整个 $\sigma$-代数。这是一个极其强大的证明工具。

!!! info "定义 3.3.3 (单调类 Monotone Class)"

    设 $\mathcal{C} \subset 2^X$ 为一集族。如果 $\mathcal{C}$ 满足以下条件，则称为**单调类**：

    * 对于递增集列 $E_1 \subset E_2 \subset \dots$，若 $E_j \in \mathcal{C}$，则 $\bigcup_{j=1}^\infty E_j \in \mathcal{C}$。

    * 对于递减集列 $E_1 \supset E_2 \supset \dots$，若 $E_j \in \mathcal{C}$，则 $\bigcap_{j=1}^\infty E_j \in \mathcal{C}$。

显然，任何 $\sigma$-代数都是单调类。反之，一个既是代数又是单调类的集族，必定是 $\sigma$-代数。

!!! success "定理 3.3.2 (单调类定理)"

    设 $\mathcal{A}$ 为一个代数。如果 $\mathcal{C}$ 是包含 $\mathcal{A}$ 的单调类，那么 $\mathcal{C}$ 必定包含由 $\mathcal{A}$ 生成的 $\sigma$-代数，即：

    \[
    \sigma(\mathcal{A}) \subset \mathcal{C}
    \]

---

## 3. 乘积测度的构造

现在我们可以定义两个测度的乘积了。

!!! success "定理 3.3.3 (截面积分与乘积测度)"

    设 $(X, \mathcal{M}, \mu)$ 和 $(Y, \mathcal{N}, \nu)$ 为 $\sigma$-有限测度空间。对于任意 $E \in \mathcal{M} \otimes \mathcal{N}$：

    * 映射 $x \mapsto \nu(E_x)$ 是 $\mathcal{M}$-可测的；

    * 映射 $y \mapsto \mu(E^y)$ 是 $\mathcal{N}$-可测的；

    * 且它们的积分相等，共同定义了 $X \times Y$ 上的**乘积测度 $\mu \times \nu$**：

    \[
    (\mu \times \nu)(E) = \int_X \nu(E_x) \, d\mu(x) = \int_Y \mu(E^y) \, d\nu(y)
    \]

??? proof "定理 3.3.3 的证明简述"

    证明的核心利用了**单调类定理**。
    
    首先，如果 $E = A \times B$ 是可测长方形，$\nu(E_x) = \nu(B) \cdot \chi_A(x)$。这是一个简单函数，显然可测，且其积分为 $\mu(A)\nu(B)$。对称地，$\int \mu(E^y)d\nu(y) = \mu(A)\nu(B)$ 也是成立的。
    
    我们将满足上述性质的所有集合 $E \in \mathcal{M} \otimes \mathcal{N}$ 收集起来，构成集族 $\mathcal{C}$。由上述讨论，有限不相交的长方形并集（即由长方形生成的代数 $\mathcal{A}$）都在 $\mathcal{C}$ 中。

    由于测度空间是 $\sigma$-有限的，我们可以利用测度的自下连续性（针对递增集列）和自上连续性（结合控制收敛定理 DCT 针对递减集列），证明 $\mathcal{C}$ 具有**单调类**的性质。
    
    根据单调类定理，包含代数 $\mathcal{A}$ 的单调类 $\mathcal{C}$ 必定包含 $\sigma(\mathcal{A}) = \mathcal{M} \otimes \mathcal{N}$。因此，对所有的可测集，该性质都成立。

---

## 4. Tonelli 定理与 Fubini 定理

这是重积分计算的两大核心定理。**Tonelli 定理**处理非负函数，而 **Fubini 定理**处理一般可积函数。它们通常配合使用。

!!! success "定理 3.3.4 (Tonelli 定理：针对非负函数)"

    设 $(X, \mathcal{M}, \mu)$ 和 $(Y, \mathcal{N}, \nu)$ 为 $\sigma$-有限测度空间。如果 $f \in L^+(X \times Y)$（即 $f$ 是非负可测函数），那么：

    * 函数 $g(x) = \int_Y f_x \, d\nu$ 和 $h(y) = \int_X f^y \, d\mu$ 分别是 $X$ 和 $Y$ 上的非负可测函数。

    * $f$ 的重积分可以通过任意顺序的累次积分计算（即使结果为 $+\infty$）：

    \[
    \int_{X \times Y} f \, d(\mu \times \nu) = \int_X \left( \int_Y f(x, y) \, d\nu(y) \right) d\mu(x) = \int_Y \left( \int_X f(x, y) \, d\mu(x) \right) d\nu(y)
    \]

??? proof "Tonelli 定理的推导过程"

    1. **特征函数**：若 $f = \chi_E$，由乘积测度的定义（定理 3.3.3），定理直接成立。
    
    2. **非负简单函数**：若 $f = \sum a_j \chi_{E_j}$，由积分的线性性质，定理成立。
    
    3. **非负可测函数**：由简单函数逼近定理，存在单调递增的非负简单函数列 $\phi_n \uparrow f$。
       对于每一个 $\phi_n$，累次积分等于重积分。对等式两边应用**单调收敛定理 (MCT)**，极限与积分号交换，即可得到任意非负函数 $f$ 满足 Tonelli 定理。

!!! success "定理 3.3.5 (Fubini 定理：针对绝对可积函数)"

    设 $(X, \mathcal{M}, \mu)$ 和 $(Y, \mathcal{N}, \nu)$ 为 $\sigma$-有限测度空间。如果 $f \in L^1(\mu \times \nu)$（即 $f$ 是乘积空间上的**绝对可积函数**），那么：

    * 对于几乎所有的 $x \in X$，$f_x \in L^1(\nu)$（即截面函数是可积的）；对于几乎所有的 $y \in Y$，$f^y \in L^1(\mu)$。

    * 对截面积分得到的函数 $g(x) = \int_Y f_x \, d\nu$ 和 $h(y) = \int_X f^y \, d\mu$ 分别在 $X$ 和 $Y$ 上是可积的。

    * 重积分等于累次积分：

    \[
    \int_{X \times Y} f \, d(\mu \times \nu) = \int_X \left( \int_Y f(x, y) \, d\nu(y) \right) d\mu(x) = \int_Y \left( \int_X f(x, y) \, d\mu(x) \right) d\nu(y)
    \]

??? proof "Fubini 定理的证明思路"

    因为 $f \in L^1$，可知 $\int |f| \, d(\mu \times \nu) < \infty$。
    
    由于 $|f|$ 是非负可测函数，我们可以对其应用 **Tonelli 定理**：

    \[
    \int_X \left( \int_Y |f(x, y)| \, d\nu(y) \right) d\mu(x) = \int_{X \times Y} |f| \, d(\mu \times \nu) < \infty
    \]

    这说明外部积分是有限的。由非负函数积分的“几乎处处”性质，被积函数 $\int_Y |f_x| d\nu$ 必须在 $X$ 上几乎处处有限。这就证明了 $f_x \in L^1(\nu)$ a.e.。
    
    最后，将实函数 $f$ 分解为正部和负部 $f = f^+ - f^-$（若为复函数则分解实虚部及正负部）。由于 $|f|$ 可积，可知 $f^+$ 和 $f^-$ 都在 $X \times Y$ 上可积。
    分别对 $f^+$ 和 $f^-$ 应用 Tonelli 定理，然后利用积分的线性性质相减，即证明了 Fubini 定理。

**总结使用策略**：在实际计算中，我们经常先对 $|f|$ 使用 Tonelli 定理来验证积分是否有限。如果 $\int \int |f| d\mu d\nu < \infty$，则 $f \in L^1(\mu \times \nu)$，此时就可以放心大胆地使用 Fubini 定理交换 $f$ 本身的积分顺序了。

---

## 5. $\mathbb{R}^n$ 上的 Lebesgue 积分简述

基于上述乘积测度的理论，高维欧氏空间 $\mathbb{R}^n$ 上的积分理论变得非常自然。

* **Lebesgue 测度的乘积**：$\mathbb{R}^n$ 上的 Borel 测度可视为 $n$ 个一维 Borel 测度的乘积。将其进行完备化处理后，就得到了 $\mathbb{R}^n$ 上的 Lebesgue 测度 $m_n$。

* **外测度逼近**：对于 $\mathbb{R}^n$ 中任何可测集 $E$ 及其 Lebesgue 测度 $m_n(E)$，我们既可以通过开集从外部逼近（$\inf \{m(U) \mid U \supset E\}$），也可以通过紧集从内部逼近（$\sup \{m(K) \mid K \subset E\}$）。

* **几何不变性**：在一维情况中，我们知道 Lebesgue 测度具有平移不变性。在 $\mathbb{R}^n$ 中，Lebesgue 测度 $m_n$ 不仅对于平移保持不变，对于空间中的**旋转变换 (Rotations)** 同样保持不变。这使得 Lebesgue 积分在处理几何与物理问题时具有天然的优势。