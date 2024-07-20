---
title: Circle STARK notes

---

# Circle STARK notes
In Starks, we typically want a small prime $p$ with $p + 1$ smooth. We choose the Mersenne prime M31.

# The Circle curve
Let the *Circle Curve* be defined as $C = C(\mathbb{F}_p)$a is a smooth algebraic variety over $\mathbb{F}_p$ defined by 

$$C: x^2 + y^2 = 1$$

**The arithmetization of the circle STARK takes place over the circle curve domain. Witnesses are encoded into low-degree polynomials over that domain. They are then (as usual) subject to a set of selected constraints.**

Another useful note is that the circle curve $C(\mathbb{F}_p) \cong P^1(\mathbb{F}_p)$ (the projective line). This also extends to all field extensions of $\mathbb{F}_p$ including $\overline{\mathbb{F}_p}$. This is constructed via a stereographic projection onto the y-axis.

## The circle curve as a group
The $p + 1$ points of $C(\mathbb{F}_p)$ form a group with the group operation
$$ (x_0, y_0) \cdot (x_1, y_1) := (x_0 \cdot x_1 - y_0 \cdot y_1, x_0 \cdot y_1 - y_0 \cdot x_1) $$

The identity element is $(1, 0)$.

Check: 
$(1, 0) \cdot (x_1, y_1) = (x_1 - 0 \cdot y_1, y_1, - 0 \cdot x_1) \stackrel{\checkmark}{=} (x_1, y_1)$

![Screenshot 2024-07-19 at 3.00.51 AM](https://hackmd.io/_uploads/SJRUtnPdA.png)

#### Couple things:
- Note that $J(J(P)) = P,  \forall P \in C(\mathbb{F}_p)$
- Also, $\pi(J(P)) = J(P) \cdot J(P) = (P_x, -P_y) \cdot (P_x, -P_y) = J(\pi(P))$
- $T_P$ (translations/rotations by P) extends over the entire projective variety of $\overline{\mathbb{F}_p}$ so that $\infty = (1: i: 0)$ and $\bar{\infty} = (1: -i: 0)$  points at infinity are **fixed under the action of the circle group.** 

Under the isomorphisim due to the stereographic projection, we can write the circle group law on a projective line as $(t_1 : s_1) \oplus (t_2: s_2) = (t_1 \cdot s_2 + t_2 \cdot s_1 : s_1 \cdot s_2 - t_1 \cdot t_2 )$ and similarly in affine coordinates. 
- So, **the circle group is isomorphic to a subgroup of linear automorphisms of the projective line**
- Another thing: $C(\mathbb{F}_p)$ is cyclic (order $p + 1$).
    - Thus, $\forall N | (p + 1)$, there exists a unique cyclic subgroup of size $N$. (Might be useful for doing a FFT kind of thing :eyes: )
    

#### Definitions:
- Any disjoint union (in **Set**) $D$ of a (cyclic) subgroup $G_{n - 1}$ of order $2^{n - 1}$ is called a **twin-coset** of size $2^n$.
    - The twin-cosets are sets with pairs that are related by the map $J$. (Note that fixed points are in the intersection, so they are **excluded form all twin-cosets** since they satisfy $Q \cdot G_{n - 1} = Q^{-1} \cdot G_{n - 1}$).
- If the twin-coset is a coset $D$ of the subgroup $G_n$, then $D$ is a **standard position coset** of size $N$.

 Twin-sets are interesting generally because they can be used to work around the non-smooth behaviour of the circle FFT under rotation. Also, they can be decomposed into twin-cosets of smaller size.
 
 ### Proposition 1. 
 The existence of standard position cosets (of size $2^n$) is *equivalent* to "$p$ is CFFT-friendly" supporting the order $n$. Note that the standard position cosets can be written as 
 $$D = Q \cdot G_n = Q \cdot G_{n - 1} \cup Q^{-1} \cdot G_{n - 1}$$ 
where $Q \in C(\mathbb{F_p})$. 
 
**Proof**: The quick reason is because the quotient group $C(\mathbb{F}_p)/G_{n - 1}$ is cyclic so there is at most one $J$-invariant pair $(Q \cdot G_{n - 1},\  Q^{-1} \cdot G_{n - 1})$ where $Q \cdot G_{n - 1} \neq  Q^{-1} \cdot G_{n - 1}$  (also forms a coset in the quotient group). Then, since the order of $Q$ is $2^{n + 1}$ we have that $Q^4$ has order $2^{n - 1}$ which is also the order of $G_{n - 1}$. So $Q \cdot G_{n -1}$ has order 4. Finally, the existence of order 4 elements in $C(\mathbb{F}_p)/G_{n - 1}$ is equivalent to the claim that $2^{n + 1} | (p + 1). \blacksquare$
 
 ### Lemma 2
 Consider CFFT-friendly prime and subgroup $G_k, k \leq m$. Then any subset $D \subseteq C(\mathbb{F}_p) - G_m$ invariant under $G_{m - 1}$ and $J$ **can be decomposed into twin-cosets** of size $2^n, n \leq m$.
 
![Screenshot 2024-07-19 at 4.01.10 AM](https://hackmd.io/_uploads/B1GYDaDO0.png)
**Proof:**  For the first assertion, D is also invariant under $_{n - 1}$ and disjoint from $G_n$ since $n \leq m$ so then $C(\mathbb{F}_p)/G_{n - 1}$ is teh union of $J$-invariant pairs. 

The long expansion for $D$ comes from $D = Q \cdot G_{m - 1} \cup Q^{-1} \cdot G_{m - 1}$.

### Lemma 3
The image $\pi(D)$ of a twin-coset $D$ of size $N = 2^n, n \geq 2$ under the squaring map $\pi$ **is also a twin-coset** of size $N/2$.
- If $D$ is a standard position coset, $\pi(D)$ also is one

Challenge for reader. Hint: follows from considering $D = Q \cdot G_{m - 1} \cup Q^{-1} \cdot G_{m - 1}$![Screenshot 2024-07-19 at 4.06.04 AM](https://hackmd.io/_uploads/Syood6wOA.png)


## The space of polynomials and circle codes
### Definition: 

Let $\mathcal{L}_N(F)$ be the **space of all bivariate polynomials with coefficients in F and of total degree at most $N / 2$ over the circle curve**
$$ \mathcal{L}_N (F) = \{ p(x, y) \in F[x, y]/(x^2 + y^2 - 1): \deg p \leq N / 2 \}$$

The important properties of $\mathcal{L}_N (F)$ are:
1) Rotation invariance (for "next-neighbor relation" and efficient encoding)
2) good separability (maximum distance separable codes)

![Screenshot 2024-07-20 at 2.28.23 AM](https://hackmd.io/_uploads/BkMrQZKd0.png)

## Circle Code
The circle code $C_N (F, D)$ (values in $F$ and evaluation in $D$) can be defined as 
$$ C_N (F, D) = \{ f (P)|_{P \in D} : f \in \mathcal{L}_N (F) \} $$

The evaluation domain is a set of $\mathbb{F}_p$-rational points. And by the isomorphism between projective lines and circles (via stereographic projection), this is actually **equivalent to a Reed-Solomon codeword**. The block length is the size of the eval domain $D$. Some interesting properties include:
- The code is invariant under rotation by $Q \in C(\mathbb{F}_p)$ and involution when $J(D) = D$
- The domain eval map $\mathcal{L}_N (F) \rightarrow F^D$ which sends a polynomial $f \rightarrow f(P)|_{P \in D}$ is **linear and injective**. 
- There is a minimum distance of a code defined by $k = \dim \mathcal{C}_N (F, D) = N + 1$ which creates a **minimum distance of a code** = $|D| - k + 1$. 

See the paper for the full proof for why Circle Reed-Solomon codes and RS codes are equivalent (Deals with using teh stereographic isomorphism to claim that the space of polynomials $\mathcal{L}_N (F)$ is equivalent to a univariate space).


For circle STARKS, we will consider **purely two-adic $N = 2^n$** and encoding will be done using FFT-based extrapolation. THis will extend a given set of values over a *standard position coset of size $N$* to be over an evaluation domain $D$. The **dimension of the code is $N + 1$**, which will need to be accounted for. 

### Vanishing polynomials and quotients
 The set of vanishing polynomials from $\mathcal{L}_N$ that vanish over $D$ form a one-dimensional subspace of $\mathcal{L}_N$ defined as 
 
 $$\nu (D) = \{ v \in \mathcal{L}_N : v |_D = 0 \}$$
 
 