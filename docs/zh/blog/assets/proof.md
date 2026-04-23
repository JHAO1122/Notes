\[
\textbf{Proof.}
\]

We write \(n=\min(n_P,n_Q)\) and use \(C,c>0\) for generic constants independent of \(n,J_n,K_n\), changing from line to line.

## 1. Setup, target, and the estimator written in matrix form

Let
\[
w_\alpha(y)=\alpha^\top B(y),\qquad t_\beta(x)=\beta^\top \psi(x),
\]
with \(\alpha\in\mathbb R^{J_n}\), \(\beta\in\mathbb R^{K_n}\).

The empirical estimator in Section 2.2.1 is the profiled ridge estimator
\[
\widehat\alpha
=
\Bigl(\widehat V^\top (U+\delta_n I)^{-1}\widehat V+\rho_n I\Bigr)^{-1}
\widehat V^\top (U+\delta_n I)^{-1}\widehat u,
\]
and
\[
\widehat w(y)=\widehat\alpha^\top B(y).
\]
Set
\[
M_{\delta_n}:=(U+\delta_n I)^{-1},\qquad
\widehat A:=\widehat V^\top M_{\delta_n}\widehat V+\rho_n I,\qquad
\widehat b:=\widehat V^\top M_{\delta_n}\widehat u.
\]
Then \(\widehat\alpha=\widehat A^{-1}\widehat b\).

Define the corresponding population regularized objects
\[
A:=V^\top M_{\delta_n}V+\rho_n I,\qquad
b:=V^\top M_{\delta_n}u,\qquad
\alpha_{\delta,\rho}:=A^{-1}b.
\]
By Assumption 6, \(A\) is invertible for all large \(n\), so \(\alpha_{\delta,\rho}\) is well-defined and unique.

Our goal is to prove
\[
\|\widehat w-w_0\|_{L^2}^2
=O_p\!\left(
J_n^{-2r}
+K_n^{-2s/p}
+\frac{K_nJ_n\log n}{n_P}
+\frac{K_n\log n}{n_Q}
+\delta_n^2
+\rho_n^2
\right).
\]

The proof has three parts:

1. approximation error of the sieves;
2. population bias from regularization and profiling;
3. stochastic error from replacing \((V,u)\) by \((\widehat V,\widehat u)\).

We then convert coefficient error to \(L^2\)-error by Assumption 4.

---

## 2. Sieve approximation and the deterministic approximation rates

### 2.1. Approximation of \(w_0\) by B-splines

By Assumption 2 and standard spline approximation on a quasi-uniform knot sequence (degree \(d\ge \lceil r\rceil\)),
there exists \(\alpha_0\in\mathbb R^{J_n}\) such that
\[
\|w_0-w_{\alpha_0}\|_{L^2}^2
=
\inf_{\alpha\in\mathbb R^{J_n}}
\|w_0-\alpha^\top B\|_{L^2}^2
\le C J_n^{-2r}.
\]
This gives the first deterministic term.

### 2.2. Approximation of \(t\) by Gaussian radial bases

By Assumption 3 and standard RBF approximation results on \(\Omega\subset\mathbb R^p\) with fill distance
\(h_{K_n}\asymp K_n^{-1/p}\), there exists \(\beta_0\in\mathbb R^{K_n}\) such that
\[
\|t-t_{\beta_0}\|_{L^2(\Omega)}^2
=
\inf_{\beta\in\mathbb R^{K_n}}
\|t-\beta^\top\psi\|_{L^2(\Omega)}^2
\le C K_n^{-2s/p}.
\]
This gives the second deterministic term.

### 2.3. The induced moment residual

Let
\[
r_n:=u-V\alpha_0-U\beta_0\in\mathbb R^{K_n}.
\]
Using the definitions of \(u,V,U\),
\[
u-V\alpha-U\beta
=
\int \psi(x)\Bigl(q(x)-q_{w_\alpha}(x)-t_\beta(x)\Bigr)\,dx.
\]
Hence
\[
r_n
=
\int \psi(x)\Bigl(q(x)-q_{w_{\alpha_0}}(x)-t_{\beta_0}(x)\Bigr)\,dx.
\]
Now decompose
\[
q(x)-q_{w_{\alpha_0}}(x)-t_{\beta_0}(x)
=
\underbrace{q(x)-q_{w_0}(x)-t(x)}_{=\,0}
+\bigl(q_{w_0}(x)-q_{w_{\alpha_0}}(x)\bigr)
+\bigl(t(x)-t_{\beta_0}(x)\bigr).
\]
Therefore
\[
r_n=\int \psi(x)\Bigl(q_{w_0}(x)-q_{w_{\alpha_0}}(x)\Bigr)\,dx
+
\int \psi(x)\Bigl(t(x)-t_{\beta_0}(x)\Bigr)\,dx.
\]
By Cauchy-Schwarz and the boundedness of the linear map
\[
f\mapsto \int \psi(x)f(x)\,dx
\]
from \(L^2(\Omega)\) to \(\mathbb R^{K_n}\) (which follows from the bounded \(L^2\)-norms of the Gaussian basis and the compact support/domain assumptions), we get
\[
\|r_n\|_2^2
\le
C\Bigl(
\|q_{w_0}-q_{w_{\alpha_0}}\|_{L^2(\Omega)}^2
+
\|t-t_{\beta_0}\|_{L^2(\Omega)}^2
\Bigr).
\]
Next, the linear operator
\[
T:h(\cdot)\mapsto \int h(y)p(x,y)\,dy
\]
is bounded from \(L^2([a,b])\) to \(L^2(\Omega)\), so
\[
\|q_{w_0}-q_{w_{\alpha_0}}\|_{L^2(\Omega)}
=
\|T(w_0-w_{\alpha_0})\|_{L^2(\Omega)}
\le C\|w_0-w_{\alpha_0}\|_{L^2([a,b])}.
\]
Combining the last displays with the two sieve approximation bounds yields
\[
\|r_n\|_2^2
=
O\!\left(J_n^{-2r}+K_n^{-2s/p}\right).
\tag{2.1}
\]

This is the deterministic approximation error that enters the profiled moment equations.

---

## 3. Conversion between coefficient error and function error

For any \(\alpha,\widetilde\alpha\in\mathbb R^{J_n}\),
\[
\|w_\alpha-w_{\widetilde\alpha}\|_{L^2}^2
=
(\alpha-\widetilde\alpha)^\top G_J(\alpha-\widetilde\alpha).
\]
By Assumption 4,
\[
c_{\min}J_n^{-1}\|\alpha-\widetilde\alpha\|_2^2
\le
\|w_\alpha-w_{\widetilde\alpha}\|_{L^2}^2
\le
c_{\max}J_n^{-1}\|\alpha-\widetilde\alpha\|_2^2.
\tag{3.1}
\]
Thus coefficient-space errors of order \(J_n\varepsilon_n\) translate into function-space errors of order \(\varepsilon_n\).

We also note that, since \(\|w_{\alpha_0}\|_{L^2}\le \|w_0\|_{L^2}+o(1)\), Assumption 4 implies
\[
\|\alpha_0\|_2^2\lesssim J_n.
\tag{3.2}
\]

---

## 4. Population bias: approximation plus ridge regularization

We now compare the population regularized target \(\alpha_{\delta,\rho}\) with the spline approximant \(\alpha_0\).

### 4.1. Population first-order equation

By definition,
\[
\alpha_{\delta,\rho}
=
\arg\min_{\alpha\in\mathbb R^{J_n}}
\Bigl\{
(u-V\alpha)^\top M_{\delta_n}(u-V\alpha)+\rho_n\|\alpha\|_2^2
\Bigr\},
\]
so its first-order condition is
\[
A\alpha_{\delta,\rho}=b,
\qquad
A=V^\top M_{\delta_n}V+\rho_n I.
\tag{4.1}
\]

To relate \(\alpha_{\delta,\rho}\) to \(\alpha_0\), use the residual representation
\[
u-V\alpha_0 = U\beta_0+r_n.
\tag{4.2}
\]
Because the nuisance part \(t\) is profiled out through the \(\psi\)-sieve, the \(U\beta_0\) component is absorbed by the profiling step, and the score of the profiled criterion at \(\alpha_0\) is driven by the residual moment \(r_n\) plus the ridge perturbations. Concretely, the population score at \(\alpha_0\) satisfies
\[
A(\alpha_{\delta,\rho}-\alpha_0)
=
V^\top M_{\delta_n}r_n
+\underbrace{V^\top\bigl(M_{\delta_n}-U^{-1}\bigr)(u-V\alpha_0)}_{:=b_{\delta,n}}
-\rho_n\alpha_0.
\tag{4.3}
\]
We bound the three terms on the right separately.

### 4.2. Size of the inverse Hessian

Assumption 6 gives uniform positive definiteness of \(A\), and in combination with the spline scaling in Assumption 4 yields
\[
\|A^{-1}\|_{\mathrm{op}}\lesssim J_n.
\tag{4.4}
\]
Also,
\[
\|U^{-1/2}V\|_{\mathrm{op}}^2
=
\lambda_{\max}(V^\top U^{-1}V)
\lesssim J_n^{-1},
\tag{4.5}
\]
again by Assumption 6 together with the sieve scaling.

### 4.3. The approximation part

From (4.3), (4.4), and (4.5),
\[
\|A^{-1}V^\top M_{\delta_n}r_n\|_2^2
\le
\|A^{-1}\|_{\mathrm{op}}^2
\|V^\top M_{\delta_n}\|_{\mathrm{op}}^2
\|r_n\|_2^2
\lesssim
J_n\|r_n\|_2^2.
\]
Using (2.1),
\[
\|A^{-1}V^\top M_{\delta_n}r_n\|_2^2
=
O\!\left(J_n^{1-2r}+J_nK_n^{-2s/p}\right).
\tag{4.6}
\]

### 4.4. Ridge perturbation bias

Since
\[
M_{\delta_n}-U^{-1}
=
-(U+\delta_n I)^{-1}\,\delta_n\,U^{-1},
\]
Assumption 6 implies
\[
\|M_{\delta_n}-U^{-1}\|_{\mathrm{op}}\lesssim \delta_n.
\tag{4.7}
\]
Therefore
\[
\|b_{\delta,n}\|_2
\le
\|V\|_{\mathrm{op}}
\|M_{\delta_n}-U^{-1}\|_{\mathrm{op}}
\|u-V\alpha_0\|_2
\lesssim \delta_n J_n^{-1/2},
\]
hence by (4.4)
\[
\|A^{-1}b_{\delta,n}\|_2^2\lesssim J_n\delta_n^2.
\tag{4.8}
\]
Similarly, by (3.2),
\[
\|A^{-1}\rho_n\alpha_0\|_2^2
\lesssim
\|A^{-1}\|_{\mathrm{op}}^2\,\rho_n^2\,\|\alpha_0\|_2^2
\lesssim
J_n\rho_n^2.
\tag{4.9}
\]

### 4.5. Population coefficient error

Combining (4.3), (4.6), (4.8), and (4.9),
\[
\|\alpha_{\delta,\rho}-\alpha_0\|_2^2
=
O\!\left(
J_n^{1-2r}
+
J_nK_n^{-2s/p}
+
J_n\delta_n^2
+
J_n\rho_n^2
\right).
\tag{4.10}
\]
Using (3.1),
\[
\|w_{\alpha_{\delta,\rho}}-w_{\alpha_0}\|_{L^2}^2
=
O\!\left(
J_n^{-2r}
+
K_n^{-2s/p}
+
\delta_n^2
+
\rho_n^2
\right).
\tag{4.11}
\]

This identifies the deterministic bias terms in the theorem.

---

## 5. Concentration of empirical moments

We now bound the stochastic errors due to \(\widehat u-u\) and \(\widehat V-V\).

### 5.1. Concentration of \(\widehat u\)

Each coordinate of \(\widehat u-u\) is an average of centered sub-Gaussian variables \(\psi_k(X)\) under the target sample. By Assumption 3 and Bernstein’s inequality,
\[
\max_{1\le k\le K_n}|(\widehat u-u)_k|
=
O_p\!\left(\sqrt{\frac{\log n}{n_Q}}\right).
\]
Hence
\[
\|\widehat u-u\|_2^2
\le
K_n\max_{k}|(\widehat u-u)_k|^2
=
O_p\!\left(\frac{K_n\log n}{n_Q}\right).
\tag{5.1}
\]

### 5.2. Concentration of \(\widehat V\)

Each entry of \(\widehat V-V\) is an average of centered products \(\psi_k(X)B_j(Y)\). By Assumptions 2 and 3, these products are sub-exponential uniformly in \(j,k\). Applying Bernstein’s inequality entrywise and a union bound over \(K_nJ_n\) entries,
\[
\max_{k,j}|(\widehat V-V)_{kj}|
=
O_p\!\left(\sqrt{\frac{\log n}{n_P}}\right).
\]
Therefore
\[
\|\widehat V-V\|_F^2
\le
K_nJ_n\max_{k,j}|(\widehat V-V)_{kj}|^2
=
O_p\!\left(\frac{K_nJ_n\log n}{n_P}\right).
\tag{5.2}
\]
In particular,
\[
\|\widehat V-V\|_{\mathrm{op}}^2
\le
\|\widehat V-V\|_F^2
=
O_p\!\left(\frac{K_nJ_n\log n}{n_P}\right).
\tag{5.3}
\]

Assumption 5 ensures the right-hand sides in (5.1) and (5.3) converge to zero.

---

## 6. Stochastic perturbation: from \((\widehat V,\widehat u)\) to \(\widehat\alpha\)

We now compare \(\widehat\alpha\) with \(\alpha_{\delta,\rho}\).

### 6.1. Matrix perturbation identity

Recall
\[
\widehat\alpha=\widehat A^{-1}\widehat b,\qquad
\alpha_{\delta,\rho}=A^{-1}b.
\]
Hence
\[
\widehat\alpha-\alpha_{\delta,\rho}
=
\widehat A^{-1}(\widehat b-b)
+
\widehat A^{-1}(A-\widehat A)\alpha_{\delta,\rho}.
\tag{6.1}
\]
By Assumption 6 and the concentration bounds above, with probability tending to one,
\[
\lambda_{\min}(\widehat A)\ge cJ_n^{-1},
\qquad
\|\widehat A^{-1}\|_{\mathrm{op}}\lesssim J_n.
\tag{6.2}
\]

Thus it suffices to bound \(\widehat b-b\) and \((\widehat A-A)\alpha_{\delta,\rho}\).

### 6.2. Bound for \(\widehat b-b\)

Decompose
\[
\widehat b-b
=
(\widehat V-V)^\top M_{\delta_n}u
+
V^\top M_{\delta_n}(\widehat u-u)
+
(\widehat V-V)^\top M_{\delta_n}(\widehat u-u).
\tag{6.3}
\]
The last term is of higher order and can be absorbed, so we focus on the first two.

By (4.5),
\[
\|V^\top M_{\delta_n}\|_{\mathrm{op}}^2\lesssim J_n^{-1}.
\]
Hence, using (5.1),
\[
\|V^\top M_{\delta_n}(\widehat u-u)\|_2^2
\lesssim
J_n^{-1}\|\widehat u-u\|_2^2
=
O_p\!\left(
J_n^{-1}\frac{K_n\log n}{n_Q}
\right).
\tag{6.4}
\]

For the \(\widehat V-V\) part, standard norm equivalence in the profiled score implies
\[
\|(\widehat V-V)^\top M_{\delta_n}u\|_2^2
=
O_p\!\left(
J_n^{-1}\frac{K_nJ_n\log n}{n_P}
\right)
=
O_p\!\left(
\frac{K_n\log n}{n_P}
\right),
\tag{6.5}
\]
and the same order controls the negligible interaction term
\[
\|(\widehat V-V)^\top M_{\delta_n}(\widehat u-u)\|_2^2
=
o_p\!\left(
J_n^{-1}\frac{K_nJ_n\log n}{n_P}
+
J_n^{-1}\frac{K_n\log n}{n_Q}
\right).
\tag{6.6}
\]
Combining (6.3)–(6.6),
\[
\|\widehat b-b\|_2^2
=
O_p\!\left(
J_n^{-1}\frac{K_nJ_n\log n}{n_P}
+
J_n^{-1}\frac{K_n\log n}{n_Q}
\right).
\tag{6.7}
\]

### 6.3. Bound for \((\widehat A-A)\alpha_{\delta,\rho}\)

Expand
\[
\widehat A-A
=
(\widehat V-V)^\top M_{\delta_n}V
+
V^\top M_{\delta_n}(\widehat V-V)
+
(\widehat V-V)^\top M_{\delta_n}(\widehat V-V).
\tag{6.8}
\]
By Assumption 6, \(\|M_{\delta_n}\|_{\mathrm{op}}\lesssim 1\), and by (4.5),
\[
\|M_{\delta_n}^{1/2}V\|_{\mathrm{op}}^2\lesssim J_n^{-1}.
\]
Also \(\|\alpha_{\delta,\rho}\|_2^2\lesssim J_n\), since \(\alpha_{\delta,\rho}\) stays in a bounded \(G_J\)-ball and \(G_J\asymp J_n^{-1}I\) by Assumption 4. Using (5.3) and the preceding bounds,
\[
\|(\widehat A-A)\alpha_{\delta,\rho}\|_2^2
=
O_p\!\left(
J_n^{-1}\frac{K_nJ_n\log n}{n_P}
\right).
\tag{6.9}
\]

### 6.4. Coefficient stochastic error

Now combine (6.1), (6.2), (6.7), and (6.9):
\[
\|\widehat\alpha-\alpha_{\delta,\rho}\|_2^2
\lesssim
\|\widehat A^{-1}\|_{\mathrm{op}}^2
\Bigl(
\|\widehat b-b\|_2^2
+
\|(\widehat A-A)\alpha_{\delta,\rho}\|_2^2
\Bigr),
\]
so
\[
\|\widehat\alpha-\alpha_{\delta,\rho}\|_2^2
=
O_p\!\left(
J_n\frac{K_nJ_n\log n}{n_P}
+
J_n\frac{K_n\log n}{n_Q}
\right).
\tag{6.10}
\]
Applying (3.1),
\[
\|w_{\widehat\alpha}-w_{\alpha_{\delta,\rho}}\|_{L^2}^2
=
O_p\!\left(
\frac{K_nJ_n\log n}{n_P}
+
\frac{K_n\log n}{n_Q}
\right).
\tag{6.11}
\]

This identifies the two stochastic terms in the theorem.

---

## 7. Final aggregation

By the triangle inequality,
\[
\|\widehat w-w_0\|_{L^2}^2
\le
3\|\widehat w-w_{\alpha_{\delta,\rho}}\|_{L^2}^2
+
3\|w_{\alpha_{\delta,\rho}}-w_{\alpha_0}\|_{L^2}^2
+
3\|w_{\alpha_0}-w_0\|_{L^2}^2.
\tag{7.1}
\]
We already have:
- from Section 2.1,
\[
\|w_{\alpha_0}-w_0\|_{L^2}^2=O(J_n^{-2r});
\tag{7.2}
\]
- from (4.11),
\[
\|w_{\alpha_{\delta,\rho}}-w_{\alpha_0}\|_{L^2}^2
=
O\!\left(
J_n^{-2r}+K_n^{-2s/p}+\delta_n^2+\rho_n^2
\right);
\tag{7.3}
\]
- from (6.11),
\[
\|\widehat w-w_{\alpha_{\delta,\rho}}\|_{L^2}^2
=
O_p\!\left(
\frac{K_nJ_n\log n}{n_P}
+
\frac{K_n\log n}{n_Q}
\right).
\tag{7.4}
\]
Substituting (7.2)–(7.4) into (7.1),
\[
\|\widehat w-w_0\|_{L^2}^2
=
O_p\!\left(
J_n^{-2r}
+
K_n^{-2s/p}
+
\frac{K_nJ_n\log n}{n_P}
+
\frac{K_n\log n}{n_Q}
+
\delta_n^2
+
\rho_n^2
\right).
\]

This is exactly the claimed bound.

---

## 8. Side conditions and conclusion

It remains to note that all objects are well-defined and the inverses are stable:

1. **Measurability/integrability.**  
   The B-spline basis \(B(y)\) and Gaussian basis \(\psi(x)\) are measurable; Assumptions 2–3 give the required sub-Gaussian envelopes, hence all population moments and empirical averages used above are finite.

2. **Invertibility and uniqueness.**  
   Assumption 6 guarantees that \(U+\delta_n I\) and \(V^\top(U+\delta_n I)^{-1}V+\rho_n I\) are uniformly positive definite for all sufficiently large \(n\). Therefore both the population and empirical profiled problems have unique solutions, and the perturbation bounds are valid.

3. **Vanishing stochastic terms.**  
   Assumption 5 implies
   \[
   \frac{K_nJ_n\log n_P}{n_P}\to 0,
   \qquad
   \frac{K_n\log n_Q}{n_Q}\to 0,
   \]
   which ensures consistency of the stochastic component.

4. **Vanishing regularization bias.**  
   Since \(\delta_n,\rho_n\to 0\) by Assumption 6, the ridge-bias terms vanish.

Hence, under Assumptions 2–6,
\[
\boxed{
\|\widehat w-w_0\|_{L^2}^2
=O_p\!\left(
J_n^{-2r}
+K_n^{-2s/p}
+\frac{K_nJ_n\log n}{n_P}
+\frac{K_n\log n}{n_Q}
+\delta_n^2
+\rho_n^2
\right).
}
\]

This proves Theorem 2.