# Chapter 1: Differentiable Manifolds

The basic idea of a differentiable manifold is: a space may globally have a complicated topological structure, but near each point it should look like some Euclidean space $\mathbb{R}^n$. We use local coordinates to transform problems on the manifold into problems in $\mathbb{R}^n$, and then use coordinate transformations to ensure that different local descriptions are consistent with each other.

This chapter successively introduces topological manifolds, coordinate charts, transition maps, $C^r$ compatibility, differentiable structures, and maximal atlases, and takes $\mathbb{R}^n$ and the sphere $S^n$ as basic examples.

---

## 1. From Coordinate Systems to Manifolds

In the plane $\mathbb{R}^2$, a point can be represented by rectangular coordinates $(x,y)$, and in an appropriate region it can also be represented by polar coordinates $(r,\theta)$. Polar coordinates fail at the origin, and the angular variable also requires a choice of branch, which shows that a coordinate system can usually only be defined on a local region.

Similarly, we may ask: can there exist a coordinate system on the sphere $S^2$ that covers the entire sphere? The answer is no. If a global coordinate chart existed, then $S^2$ would be homeomorphic to some open subset of $\mathbb{R}^2$; but $S^2$ is compact, whereas no nonempty open subset of $\mathbb{R}^2$ can be compact. Therefore, the study of manifolds must use multiple local coordinate charts.

!!! note "Basic Idea of Local Coordinates"

    For each point $p$ on the manifold, choose a neighborhood $U$ of it, and use a homeomorphism

    \[
    \varphi:U\longrightarrow\varphi(U)\subseteq\mathbb{R}^n
    \]

    to represent points in $U$ by $n$ real coordinates.

---

## 2. Topological Manifolds

!!! info "Definition 2.1 (Topological Manifold)"

    Let $M$ be a Hausdorff topological space. If for every $p\in M$, there exist an open neighborhood $U\subseteq M$ of $p$ and a homeomorphism

    \[
    \varphi:U\longrightarrow\varphi(U),
    \]

    where $\varphi(U)$ is an open subset of $\mathbb{R}^n$, then $M$ is called an **$n$-dimensional topological manifold**.

    This course also stipulates that manifolds have a countable topological basis, that is, they satisfy the second countability axiom.

A topological manifold can also be summarized as: a **locally Euclidean Hausdorff space**. Here “locally Euclidean” means that every point has a neighborhood homeomorphic to some open subset of $\mathbb{R}^n$.

!!! note "The Dimension of a Manifold Is Well Defined"

    The same nonempty open set cannot be homeomorphic both to an open subset of $\mathbb{R}^m$ and to an open subset of $\mathbb{R}^n$ with $m\neq n$. This fact follows from the invariance of domain theorem in topology. Therefore, the dimension $n$ of a topological manifold does not depend on the choice of local coordinates.

### 2.1 Basic Examples and Counterexamples

!!! example "Example 1 (Euclidean Space)"

    $\mathbb{R}^n$ itself is an $n$-dimensional topological manifold. For every point one can directly take $U=\mathbb{R}^n$, with coordinate map the identity map

    \[
    \operatorname{id}_{\mathbb{R}^n}:\mathbb{R}^n\longrightarrow\mathbb{R}^n.
    \]

!!! example "Example 2 (Sphere and Torus)"

    The sphere $S^n$ is an $n$-dimensional topological manifold, but usually cannot be covered by a single coordinate chart. The torus

    \[
    T^2=S^1\times S^1
    \]

    is also a two-dimensional topological manifold.

!!! warning "Counterexample (The Cone Point of a Double Cone)"

    After removing the cone point, a double cone is locally homeomorphic to $\mathbb{R}^2$, and hence is a two-dimensional topological manifold; but the cone point itself has no Euclidean neighborhood.

    Intuitively, after deleting the cone point, a sufficiently small neighborhood of it splits into two connected components; whereas after deleting the center point from a small disk in $\mathbb{R}^2$, the resulting set is still connected. Therefore, near the cone point it cannot be homeomorphic to an open subset of the plane.

---

## 3. Coordinate Charts and Local Coordinates

!!! info "Definition 3.1 (Coordinate Chart)"

    Let $M$ be an $n$-dimensional topological manifold. A **coordinate chart** is a pair $(U,\varphi)$, where $U\subseteq M$ is an open set, and

    \[
    \varphi:U\longrightarrow\varphi(U)\subseteq\mathbb{R}^n
    \]

    is a homeomorphism from $U$ to the open set $\varphi(U)$ in $\mathbb{R}^n$.

    The set $U$ is called the coordinate neighborhood, and the map $\varphi$ is called the coordinate map.

Write the coordinate map as

\[
\varphi(p)=\bigl(x^1(p),\ldots,x^n(p)\bigr),
\]

then each component

\[
x^i:U\longrightarrow\mathbb{R},
\qquad
p\longmapsto x^i(p),
\qquad i=1,\ldots,n,
\]

is called the $i$-th **local coordinate function** of the coordinate chart $(U,\varphi)$. Therefore, $(x^1,\ldots,x^n)$ is also called a local coordinate system on $U$.

---

## 4. Transition Maps

Let $(U,\varphi)$ and $(V,\psi)$ be two coordinate charts on the same $n$-dimensional topological manifold $M$, and assume

\[
U\cap V\neq\varnothing.
\]

In the overlap region, the same point has both $x$-coordinates and $y$-coordinates. This yields the **transition map**

\[
\psi\circ\varphi^{-1}:
\varphi(U\cap V)\longrightarrow\psi(U\cap V).
\]

If we write

\[
\varphi(p)=(x^1,\ldots,x^n),
\qquad
\psi(p)=(y^1,\ldots,y^n),
\]

then the components of the transition map can be expressed as

\[
y^i=y^i(x^1,\ldots,x^n),
\qquad i=1,\ldots,n.
\]

The coordinate transformation in the reverse direction is

\[
\varphi\circ\psi^{-1}:
\psi(U\cap V)\longrightarrow\varphi(U\cap V),
\]

whose component form is

\[
x^i=x^i(y^1,\ldots,y^n),
\qquad i=1,\ldots,n.
\]

---

## 5. $C^r$ Compatibility of Coordinate Charts

!!! info "Definition 5.1 ($C^r$ Compatible)"

    Let $(U,\varphi)$ and $(V,\psi)$ be two coordinate charts on $M$.

    * If $U\cap V=\varnothing$, the two coordinate charts are by convention $C^r$ compatible;
    * If $U\cap V\neq\varnothing$, and the transition map

    \[
    \psi\circ\varphi^{-1}:
    \varphi(U\cap V)\longrightarrow\psi(U\cap V)
    \]

    is a $C^r$ diffeomorphism, then $(U,\varphi)$ and $(V,\psi)$ are called **$C^r$ compatible**.

Because

\[
(\psi\circ\varphi^{-1})^{-1}
=\varphi\circ\psi^{-1},
\]

requiring the transition map to be a $C^r$ diffeomorphism is equivalent to requiring the coordinate transformations in both directions to be $C^r$.

!!! warning "$C^r$ Compatibility Is Generally Not an Equivalence Relation"

    The $C^r$ compatibility relation for coordinate charts is reflexive and symmetric, but it is not necessarily transitive on an arbitrary set of coordinate charts.

    For example, if the second coordinate chart is disjoint from the first and third coordinate charts, respectively, then the first two pairs of coordinate charts are automatically compatible; however, the first and third coordinate charts may overlap, and the transition map between them need not be $C^r$. Therefore, the $C^r$ compatibility relation generally cannot be directly regarded as an equivalence relation.

---

## 6. Atlases and Differentiable Structures

!!! info "Definition 6.1 (Coordinate Atlas)"

    A family of coordinate charts on $M$

    \[
    \mathcal{A}=\{(U_\alpha,\varphi_\alpha):\alpha\in\Lambda\}
    \]

    is called a **coordinate atlas** if these coordinate neighborhoods cover $M$:

    \[
    \bigcup_{\alpha\in\Lambda}U_\alpha=M.
    \]

    If any two coordinate charts in the atlas are $C^r$ compatible, then $\mathcal{A}$ is called a $C^r$ atlas.

!!! info "Definition 6.2 ($C^r$ Differentiable Structure)"

    A $C^r$ atlas

    \[
    \mathcal{U}=\{(U_\alpha,\varphi_\alpha):\alpha\in\Lambda\}
    \]

    is called a **$C^r$ differentiable structure** on $M$ if it also satisfies maximality: for any coordinate chart $(V,\psi)$, as long as $(V,\psi)$ is $C^r$ compatible with every coordinate chart in $\mathcal{U}$, we must have

    \[
    (V,\psi)\in\mathcal{U}.
    \]

    In other words, a $C^r$ differentiable structure is a maximal $C^r$ atlas.

!!! info "Definition 6.3 ($C^r$ Manifolds and Smooth Manifolds)"

    A topological manifold $M$ endowed with a $C^r$ differentiable structure $\mathcal{U}$ is called a **$C^r$ differentiable manifold**, denoted by $(M,\mathcal{U})$.

    When $r=\infty$, $\mathcal{U}$ is called a **smooth structure** on $M$, and $(M,\mathcal{U})$ is called a **smooth manifold**.

Usually one does not need to list all coordinate charts in the maximal atlas one by one. As long as one gives a family of mutually compatible coordinate charts covering $M$, the corresponding maximal atlas is uniquely determined.

!!! success "Theorem 6.4 (Unique Maximal Extension of a Compatible Atlas)"

    Let $M$ be an $n$-dimensional topological manifold, and let $\mathcal{A}$ be a $C^r$ atlas covering $M$. Then there exists a unique maximal $C^r$ atlas $\mathcal{U}_{\max}$ containing $\mathcal{A}$.

??? proof "Proof of Theorem 6.4 (click to expand)"

    Define

    \[
    \mathcal{U}_{\max}
    =\{(V,\psi):(V,\psi)\text{ is }C^r
    \text{ compatible with every coordinate chart in }\mathcal{A}\}.
    \]

    Clearly $\mathcal{A}\subseteq\mathcal{U}_{\max}$, so $\mathcal{U}_{\max}$ covers $M$.

    Next we prove that any two coordinate charts $(V,\psi)$ and $(W,\chi)$ in $\mathcal{U}_{\max}$ are compatible. Fix $p\in V\cap W$. Since $\mathcal{A}$ covers $M$, there exists $(U_\alpha,\varphi_\alpha)\in\mathcal{A}$ such that $p\in U_\alpha$. In a sufficiently small neighborhood of $p$, the transition map factors as

    \[
    \chi\circ\psi^{-1}
    =(\chi\circ\varphi_\alpha^{-1})
    \circ(\varphi_\alpha\circ\psi^{-1}).
    \]

    Because $(V,\psi)$ and $(W,\chi)$ are respectively compatible with $(U_\alpha,\varphi_\alpha)$, the two maps on the right are $C^r$, so $\chi\circ\psi^{-1}$ is $C^r$ near $p$. Since $p$ is arbitrary and smoothness is a local property, the two coordinate charts are $C^r$ compatible. Therefore $\mathcal{U}_{\max}$ is a $C^r$ atlas.

    By definition, any coordinate chart compatible with all coordinate charts in $\mathcal{U}_{\max}$ is, in particular, compatible with all coordinate charts in $\mathcal{A}$, and hence already belongs to $\mathcal{U}_{\max}$. Therefore $\mathcal{U}_{\max}$ is maximal.

    Finally, if another maximal $C^r$ atlas $\mathcal{V}$ also contains $\mathcal{A}$, then every coordinate chart in $\mathcal{V}$ is compatible with $\mathcal{A}$, so $\mathcal{V}\subseteq\mathcal{U}_{\max}$. By the maximality of $\mathcal{V}$, the two are equal. Hence the maximal extension is unique. $\square$

---

## 7. Basic Examples

### 7.1 The Standard Smooth Structure on $\mathbb{R}^n$

The single coordinate chart on $\mathbb{R}^n$

\[
(\mathbb{R}^n,\operatorname{id}_{\mathbb{R}^n})
\]

constitutes a smooth atlas. By Theorem 6.4, it uniquely determines the maximal smooth atlas on $\mathbb{R}^n$, which is called the **standard smooth structure** on $\mathbb{R}^n$.

### 7.2 The Smooth Structure on the Sphere $S^n$

Write the unit sphere as

\[
S^n
=\left\{(x^1,\ldots,x^{n+1})\in\mathbb{R}^{n+1}:
\sum_{j=1}^{n+1}(x^j)^2=1\right\}.
\]

The sphere carries the subspace topology inherited from $\mathbb{R}^{n+1}$. For $i=1,\ldots,n+1$, define the open sets

\[
U_i^+=\{x\in S^n:x^i>0\},
\qquad
U_i^-=\{x\in S^n:x^i<0\}.
\]

These sets cover the sphere:

\[
S^n=\bigcup_{i=1}^{n+1}(U_i^+\cup U_i^-).
\]

Let

\[
D^n=\left\{(u^1,\ldots,u^n)\in\mathbb{R}^n:
\sum_{j=1}^n(u^j)^2<1\right\}
\]

be the open unit ball in $\mathbb{R}^n$. Define the coordinate maps

\[
\varphi_i^\pm:U_i^\pm\longrightarrow D^n
\]

to be the projection that deletes the $i$-th coordinate:

\[
\varphi_i^\pm(x^1,\ldots,x^{n+1})
=(x^1,\ldots,\widehat{x^i},\ldots,x^{n+1}),
\]

where the notation $\widehat{x^i}$ indicates that the $i$-th component is omitted.

Its inverse recovers the $i$-th coordinate through the sphere equation:

\[
(\varphi_i^\pm)^{-1}(u^1,\ldots,u^n)
=\left(u^1,\ldots,u^{i-1},
\pm\sqrt{1-\sum_{j=1}^n(u^j)^2},
u^i,\ldots,u^n\right).
\]

Therefore, each $(U_i^\pm,\varphi_i^\pm)$ is a coordinate chart on the sphere.

!!! success "Proposition 7.1 (The Sphere Coordinate Charts Are Mutually Smoothly Compatible)"

    The family of coordinate charts

    \[
    \mathcal{A}
    =\{(U_i^+,\varphi_i^+),(U_i^-,\varphi_i^-):
    i=1,\ldots,n+1\}
    \]

    is a smooth atlas on $S^n$, and hence uniquely determines the standard smooth structure on $S^n$.

??? proof "Proof of Proposition 7.1 (click to expand)"

    It is known that these coordinate neighborhoods cover $S^n$, so it suffices to verify that the transition maps on overlap regions are smooth. The general case is entirely analogous; below we specifically examine $(U_1^+,\varphi_1^+)$ and $(U_2^-,\varphi_2^-)$.

    For

    \[
    u=(u^1,\ldots,u^n)\in
    \varphi_1^+(U_1^+\cap U_2^-),
    \]

    we have

    \[
    (\varphi_1^+)^{-1}(u)
    =\left(\sqrt{1-\sum_{j=1}^n(u^j)^2},
    u^1,u^2,\ldots,u^n\right).
    \]

    Since $x^2<0$ in the intersection, the domain of the transition map is

    \[
    \varphi_1^+(U_1^+\cap U_2^-)
    =\{u\in D^n:u^1<0\}.
    \]

    Deleting the second coordinate, we obtain

    \[
    \varphi_2^-\circ(\varphi_1^+)^{-1}(u^1,\ldots,u^n)
    =\left(\sqrt{1-\sum_{j=1}^n(u^j)^2},
    u^2,\ldots,u^n\right).
    \]

    This map is $C^\infty$ on its domain.

    Conversely, for

    \[
    v=(v^1,\ldots,v^n)\in
    \varphi_2^-(U_1^+\cap U_2^-),
    \]

    we have

    \[
    (\varphi_2^-)^{-1}(v)
    =\left(v^1,-\sqrt{1-\sum_{j=1}^n(v^j)^2},
    v^2,\ldots,v^n\right),
    \]

    where the domain is

    \[
    \varphi_2^-(U_1^+\cap U_2^-)
    =\{v\in D^n:v^1>0\}.
    \]

    Therefore

    \[
    \varphi_1^+\circ(\varphi_2^-)^{-1}(v^1,\ldots,v^n)
    =\left(-\sqrt{1-\sum_{j=1}^n(v^j)^2},
    v^2,\ldots,v^n\right),
    \]

    is also $C^\infty$. The transition maps between the remaining coordinate charts have exactly the same form: rearrange some coordinates and use

    \[
    \pm\sqrt{1-\sum_j(u^j)^2}
    \]

    to recover the omitted coordinate. Therefore all transition maps are smooth, and $\mathcal{A}$ is a smooth atlas on $S^n$. $\square$

---

## 8. Summary of This Chapter

1. A topological manifold is a Hausdorff space that is locally homeomorphic to an open subset of $\mathbb{R}^n$; this course also requires second countability.
2. A coordinate chart $(U,\varphi)$ provides local Euclidean coordinates for points on the manifold.
3. Different coordinate charts are related through the transition map $\psi\circ\varphi^{-1}$.
4. If all transition maps are $C^r$ diffeomorphisms, then the coordinate atlas gives a $C^r$ differentiable structure.
5. Any compatible $C^r$ atlas covering $M$ uniquely determines a maximal $C^r$ atlas.
6. Both $\mathbb{R}^n$ and $S^n$ have natural standard smooth structures.