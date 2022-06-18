# Theoretical Foundations of Calabi

This document presents the mathematical and physical foundations of Calabi: complex differential geometry, the Calabi-Yau theorem, Donaldson's balanced metric algorithm, Hodge theory, and mirror symmetry.

---

## 1. The Calabi Conjecture and Yau's Theorem

Let $X$ be a compact complex Kähler manifold of complex dimension $n$, equipped with a Kähler metric $g$ whose Kähler form is:
$$\omega = \frac{i}{2} g_{i\bar{j}} dz^i \wedge d\bar{z}^j$$

The Ricci curvature tensor of $(X, g)$ is given in local coordinates by:
$$R_{i\bar{j}} = -\partial_i \partial_{\bar{j}} \ln \det(g)$$

The Ricci form $\rho = \frac{i}{2\pi} R_{i\bar{j}} dz^i \wedge d\bar{z}^j$ represents the first Chern class $c_1(X) \in H^2(X, \mathbb{R})$.

In 1954 and 1957, Eugenio Calabi conjectured that if $c_1(X) = 0$, there exists a unique Kähler metric in each Kähler class $[\omega]$ whose Ricci curvature vanishes identically:
$$R_{i\bar{j}} = 0$$

In 1977-1978, Shing-Tung Yau proved the Calabi conjecture by solving the complex Monge-Ampère equation:
$$(\omega + i \partial \bar{\partial} \phi)^n = e^f \omega^n$$
where $\phi$ is the Kähler potential correction and $f$ is the Ricci potential.

A compact Kähler manifold with vanishing first Chern class $c_1(X) = 0$ is called a **Calabi-Yau manifold**. By Yau's theorem, every such manifold admits a unique Ricci-flat Kähler metric. Furthermore, the holonomy group reduces from $\mathrm{U}(n)$ to the special unitary group $\mathrm{SU}(n)$.

---

## 2. The Poincaré Residue Theorem

For a smooth hypersurface $X \subset \mathbb{P}^{n+1}$ defined by the vanishing of a single homogeneous polynomial $P(z) = 0$ of degree $d = n + 2$, the adjunction formula gives:
$$K_X = (K_{\mathbb{P}^{n+1}} \otimes \mathcal{O}(d))|_X = (\mathcal{O}(-(n+2)) \otimes \mathcal{O}(n+2))|_X = \mathcal{O}_X$$

Therefore, the canonical bundle $K_X$ is trivial, and $X$ admits a nowhere-vanishing holomorphic $(n, 0)$-form $\Omega$.

By the Poincaré residue theorem, $\Omega$ is obtained algebraically as:
$$\Omega = \operatorname{Res}_X \left( \frac{\mu}{P(z)} \right)$$
where $\mu$ is the Euler projective volume form on $\mathbb{P}^{n+1}$:
$$\mu = \sum_{i=0}^{n+1} (-1)^i z_i dz_0 \wedge \dots \wedge \widehat{dz_i} \wedge \dots \wedge dz_{n+1}$$

In an affine chart $U_k$ where $w_j = z_j / z_k$, eliminating coordinate $w_m$ via $\frac{\partial P}{\partial w_m} \neq 0$:
$$\Omega = \frac{(-1)^m dw_1 \wedge \dots \wedge \widehat{dw_m} \wedge \dots \wedge dw_n}{\partial P / \partial w_m}$$

The unique invariant volume form on $X$ is given by:
$$d\mu_{CY} = \frac{i^{n^2}}{n!} \Omega \wedge \bar{\Omega} = \frac{1}{|\partial P / \partial w_m|^2} \prod_{j \neq m} \frac{i}{2} dw_j \wedge d\bar{w}_j$$

---

## 3. Donaldson's Balanced Metric Algorithm

Although Yau's theorem proves the existence and uniqueness of the Ricci-flat metric, no general closed-form analytic expression exists for compact Calabi-Yau manifolds.

In 2001 and 2005, Simon Donaldson introduced an algorithm to approximate the Ricci-flat metric via algebraic embeddings.

Let $L = \mathcal{O}_X(1)$ be an ample line bundle on $X$. For each positive integer $k \ge 1$, consider the vector space of global holomorphic sections $H^0(X, L^k)$, having dimension $N_k$. Let $\{s_\alpha\}_{\alpha=1}^{N_k}$ be a basis of sections.

Given a Hermitian metric $h_{\alpha\beta}$ on $H^0(X, L^k)$ with inverse $h^{\alpha\beta}$, the induced Fubini-Study type metric on $X$ is given by the Kähler potential:
$$K_h(z) = \frac{1}{k} \ln \left( \sum_{\alpha,\beta=1}^{N_k} h^{\alpha\beta} s_\alpha(z) \overline{s_\beta(z)} \right)$$

The Donaldson $T$-operator maps a Hermitian metric $h$ to a new Hermitian metric $T(h)$ via integration over the Calabi-Yau manifold:
$$T(h)_{\alpha\beta} = \frac{N_k}{\operatorname{Vol}(X)} \int_X \frac{s_\alpha(x) \overline{s_\beta(x)}}{\sum_{\gamma,\delta=1}^{N_k} h^{\gamma\delta} s_\gamma(x) \overline{s_\delta(x)}} d\mu_{CY}(x)$$

A metric $h^*$ satisfying $T(h^*) = h^*$ is called a **balanced metric**.

Donaldson proved that:
1. For every $k \ge 1$, a balanced metric $h^*$ exists and is unique up to scaling.
2. The sequence of metrics $g_k = \frac{1}{k} \partial \bar{\partial} \ln |s|_{h^*}^2$ converges smoothly in the $C^\infty$ topology to the unique Ricci-flat metric $g_{CY}$ as $k \to \infty$.
3. The fixed point iteration $h^{(m+1)} = T(h^{(m)})$ converges.

---

## 4. The Monge-Ampère Sigma Error

To quantify the precision of a numerical metric $g_{i\bar{j}}$, we test the complex Monge-Ampère equation:
$$\det(g_{i\bar{j}}) = \kappa \cdot |\Omega|^2$$
where $\kappa$ is the volume ratio:
$$\kappa = \frac{\int_X \det(g) d^n u}{\int_X |\Omega|^2 d^n u}$$

The Donaldson $\sigma$-measure error is defined by:
$$\sigma = \frac{1}{\operatorname{Vol}(X)} \int_X \left| 1 - \frac{\det(g_{i\bar{j}}(x))}{\kappa \cdot |\Omega(x)|^2} \right| d\mu_{CY}(x)$$

For a Ricci-flat metric, $\det(g) = \kappa |\Omega|^2$ identically, so $\sigma = 0$.

---

## 5. Mirror Symmetry and Gromov-Witten Invariants

The Dwork family of quintic threefolds is parameterized by $\psi \in \mathbb{C}$:
$$X_\psi: z_0^5 + z_1^5 + z_2^5 + z_3^5 + z_4^5 + \psi z_0 z_1 z_2 z_3 z_4 = 0$$

Under the quotient by $G = (\mathbb{Z}_5)^3$, the mirror manifold $W = X_\psi / G$ has inverted Hodge numbers:
$$h^{1,1}(W) = 101, \quad h^{2,1}(W) = 1, \quad \chi(W) = +200$$

The periods of the holomorphic 3-form satisfy the fourth-order Picard-Fuchs equation:
$$\left[ \theta^4 - 5 z (5\theta + 1)(5\theta + 2)(5\theta + 3)(5\theta + 4) \right] \varpi(z) = 0$$
where $z = (5\psi)^{-5}$ and $\theta = z \frac{d}{dz}$.

Using the Frobenius method around the Large Complex Structure Limit ($z = 0$), the fundamental periods are:
$$\varpi_0(z) = \sum_{n=0}^\infty \frac{(5n)!}{(n!)^5} z^n$$
$$\varpi_1(z) = \varpi_0(z) \ln z + 5 \sum_{n=1}^\infty \frac{(5n)!}{(n!)^5} (H_{5n} - H_n) z^n$$

The flat coordinate on the mirror moduli space is:
$$t(z) = \frac{1}{2\pi i} \frac{\varpi_1(z)}{\varpi_0(z)}, \quad q = e^{2\pi i t}$$

The topological Yukawa coupling in flat coordinates is:
$$C_{ttt}(q) = \frac{5}{1 - 5^5 z(q)} \frac{1}{\varpi_0(z(q))^2} \left( \frac{q}{z(q)} \frac{dz}{dq} \right)^3$$

By topological string theory (the A-model / B-model duality discovered by Candelas, de la Ossa, Green, and Parkes in 1991), $C_{ttt}$ generates the genus-0 Gromov-Witten invariants $n_d$:
$$C_{ttt}(q) = 5 + \sum_{d=1}^\infty n_d d^3 \frac{q^d}{1 - q^d}$$

Inversion yields the instanton numbers counting rational curves of degree $d$ on the quintic threefold:
* $n_1 = 2875$ (lines, Schubert 1879)
* $n_2 = 609250$ (conics, Katz 1986)
* $n_3 = 317206375$ (twisted cubics, Ellingsrud-Strømme 1991)
* $n_4 = 242467530000$ (rational quartics, Candelas et al. 1991)
