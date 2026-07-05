# 📏 Measure Theory Basics

This module covers the core content of set theory and measure theory, including Lebesgue measure, measurable functions, Lebesgue integral theory, and the three major convergence theorems.

---

## I. Set Algebras and Cardinality

!!! info "Equipotence and Cardinality"
    If there exists a bijection between sets \(A\) and \(B\), then \(A\) and \(B\) are said to be equipotent, denoted \(A \sim B\).

    * **Countable set**: A set equipotent to the positive integers \(\mathbb{N}\) is called a countable set, and its cardinality is denoted \(\aleph_0\).

    * **Continuum cardinality**: The cardinality of a set equipotent to the real numbers \(\mathbb{R}\) is denoted \(\mathfrak{c}\). It is known that \(\mathfrak{c} = 2^{\aleph_0}\).

??? note "Bernstein's Theorem"
    Let \(A\) and \(B\) be two sets. If there exists an injection from \(A\) to \(B\), and also an injection from \(B\) to \(A\), then \(A \sim B\).

!!! success "Open Sets, Closed Sets, and the Cantor Set"
    * **Structure of open sets**: Any non‑empty open set in \(\mathbb{R}\) can be uniquely represented as a union of at most countably many disjoint open intervals.

    * **Cantor ternary set (\(C\))**: Constructed by repeatedly removing the middle third open interval from the closed interval \([0, 1]\).
        * **Properties**: \(C\) is a compact, perfect, nowhere dense set; its cardinality is \(\mathfrak{c}\), but its Lebesgue measure is 0.

---

## II. Lebesgue Measure Theory

!!! info "Outer Measure"
    For any subset \(E\) of \(\mathbb{R}^n\), the Lebesgue outer measure \(m^*(E)\) is defined as:

    \[
    m^*(E) = \inf \left\{ \sum_{i=1}^\infty |I_i| \ \middle|\  E \subset \bigcup_{i=1}^\infty I_i \right\}
    \]

    where \(\{I_i\}\) is an at most countable sequence of open rectangles covering \(E\), and \(|I_i|\) denotes the volume of the rectangle.

!!! success "Carathéodory's Definition of Measurability"
    A set \(E \subset \mathbb{R}^n\) is called **Lebesgue measurable** if for every “test set” \(T \subset \mathbb{R}^n\),

    \[
    m^*(T) = m^*(T \cap E) + m^*(T \cap E^c)
    \]

    holds. In this case, its outer measure is called the Lebesgue measure, denoted \(m(E) = m^*(E)\). The collection of all measurable sets forms a \(\sigma\)-algebra.

??? note "Continuity of Measure"
    1. **Continuity of increasing sequence**: If \(E_1 \subset E_2 \subset \dots \subset E_n \subset \dots\) are all measurable, then

    \[
    m\left( \bigcup_{n=1}^\infty E_n \right) = \lim_{n \to \infty} m(E_n)
    \]

    2. **Continuity of decreasing sequence**: If \(E_1 \supset E_2 \supset \dots \supset E_n \supset \dots\) are all measurable **and \(m(E_1) < \infty\)**, then

    \[
    m\left( \bigcap_{n=1}^\infty E_n \right) = \lim_{n \to \infty} m(E_n)
    \]

---

## III. Measurable Functions and Their Convergence

!!! abstract "Definition of Measurable Function"
    Let \(E \subset \mathbb{R}^n\) be a measurable set, and let \(f\) be an extended real‑valued function defined on \(E\). If for every real number \(a\) the set

    \[
    \{x \in E \mid f(x) > a\}
    \]

    is Lebesgue measurable, then \(f\) is called a **measurable function** on \(E\).

!!! info "Almost Everywhere (a.e.)"
    If a property holds at all points of \(E\) except possibly on a subset of measure zero, then it is said to hold **almost everywhere** on \(E\).

??? tip "Relations among Four Types of Convergence"
    Let \(\{f_n\}\) and \(f\) be finite measurable functions on \(E\).

    1. **Almost everywhere convergence (\(f_n \xrightarrow{a.e.} f\))**: \(m(\{x \in E \mid \lim_{n\to\infty} f_n(x) \neq f(x)\}) = 0\).

    2. **Convergence in measure (\(f_n \xrightarrow{m} f\))**: For every \(\epsilon > 0\), \(\lim_{n \to \infty} m(\{x \in E \mid |f_n(x) - f(x)| > \epsilon\}) = 0\).

    * **Implication chain and conversion**:
        * If \(m(E) < \infty\), then \(f_n \xrightarrow{a.e.} f \implies f_n \xrightarrow{m} f\) (Egorov’s theorem).
    
        * If \(f_n \xrightarrow{m} f\), then there exists a **subsequence** \(\{f_{n_k}\}\) such that \(f_{n_k} \xrightarrow{a.e.} f\) (Riesz’s theorem).

!!! success "Two Core Theorems (Egorov & Lusin)"
    * **Egorov’s theorem**: Let \(m(E) < \infty\). If \(f_n \xrightarrow{a.e.} f\), then for any \(\delta > 0\) there exists a measurable subset \(E_0 \subset E\) with \(m(E \setminus E_0) < \delta\) such that \(\{f_n\}\) converges **uniformly** to \(f\) on \(E_0\).

    * **Lusin’s theorem**: Let \(f\) be a measurable function on \(E\) that is finite almost everywhere. Then for any \(\delta > 0\) there exists a closed set \(F \subset E\) with \(m(E \setminus F) < \delta\) such that the restriction of \(f\) to \(F\) is **continuous**.

---

## IV. Lebesgue Integral and Dominated Convergence Theorems

!!! info "Integral of Non‑negative Simple Functions"
    Let \(\phi(x) = \sum_{i=1}^k c_i \chi_{A_i}(x)\) be a non‑negative simple function on \(E\), where the \(A_i\) are pairwise disjoint measurable sets. Its Lebesgue integral is defined as

    \[
    \int_E \phi(x) dx = \sum_{i=1}^k c_i m(A_i)
    \]

!!! success "Lebesgue Integral of General Functions"
    * **Non‑negative measurable function**: Define \(\int_E f(x) dx = \sup \{ \int_E \phi(x) dx \mid 0 \le \phi \le f,\ \phi \text{ is simple} \}\).

    * **General measurable function**: Introduce the positive part \(f^+ = \max(f, 0)\) and the negative part \(f^- = \max(-f, 0)\). If at least one of \(\int_E f^+ dx\) and \(\int_E f^- dx\) is finite, then define

    \[
    \int_E f(x) dx = \int_E f^+(x) dx - \int_E f^-(x) dx
    \]

    If both parts are finite, \(f\) is said to be **Lebesgue integrable** on \(E\).

??? success "Three Core Limit Theorems for Integrals"
    * **1. Levi’s Monotone Convergence Theorem**: If \(0 \le f_1(x) \le f_2(x) \le \dots \le f_n(x) \le \dots\) and \(f_n(x) \xrightarrow{a.e.} f(x)\), then

    \[
    \lim_{n \to \infty} \int_E f_n(x) dx = \int_E f(x) dx
    \]

    * **2. Fatou’s Lemma**: If \(\{f_n\}\) is a sequence of non‑negative measurable functions on \(E\), then

    \[
    \int_E \liminf_{n \to \infty} f_n(x) dx \le \liminf_{n \to \infty} \int_E f_n(x) dx
    \]

    * **3. Lebesgue’s Dominated Convergence Theorem (LDCT)**: Let \(\{f_n\}\) be a sequence of measurable functions and \(f_n(x) \xrightarrow{a.e.} f(x)\). If there exists a non‑negative **integrable** function \(F(x)\) such that \(|f_n(x)| \le F(x)\) a.e. for all \(n\), then \(f\) is also integrable, and

    \[
    \lim_{n \to \infty} \int_E f_n(x) dx = \int_E f(x) dx
    \]

---

## V. Integration Spaces and Differentiation Theory

!!! info "Differentiation under the Integral Sign and Fubini’s Theorem"
    * **Fubini’s theorem**: If \(f(x, y)\) is non‑negative measurable (or integrable) on the measurable set \(E_X \times E_Y\) in \(\mathbb{R}^p \times \mathbb{R}^q\), then the double integral equals the iterated integrals:

    \[
    \int_{E_X \times E_Y} f(x, y) dxdy = \int_{E_X} \left( \int_{E_Y} f(x, y) dy \right) dx = \int_{E_Y} \left( \int_{E_X} f(x, y) dx \right) dy
    \]

!!! success "Functions of Bounded Variation and Absolutely Continuous Functions"
    * **Functions of bounded variation (\(BV\))**: If the total variation \(V_a^b(f)\) of \(f\) on \([a, b]\) is finite, then \(f\) can be expressed as the difference of two increasing functions (Jordan decomposition), and \(f\) is differentiable almost everywhere.

    * **Absolutely continuous functions (\(AC\))**: For any \(\epsilon > 0\), there exists \(\delta > 0\) such that for any finite collection of pairwise disjoint open intervals \(\{(a_i, b_i)\}_{i=1}^m\) in \([a, b]\) with \(\sum_{i=1}^m (b_i - a_i) < \delta\), we have

    \[
    \sum_{i=1}^m |f(b_i) - f(a_i)| < \epsilon
    \]

    * **Necessary and sufficient condition for the Fundamental Theorem of Calculus**: The formula \(f(x) - f(a) = \int_a^x f'(t) dt\) holds **if and only if** \(f(x)\) is absolutely continuous on \([a, b]\).