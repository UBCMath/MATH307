# QR by Reflectors

```{div} definition
An **elementary reflector** (or **Householder transformation**) is matrix of the form

$$
H = I - \frac{2}{\| \boldsymbol{u} \|^2} \boldsymbol{u} \boldsymbol{u}^T
$$

for some nonzero vector $\boldsymbol{u} \in \mathbb{R}^n$. Note that a reflector is an orthogonal matrix since $H^T = H$ and $H^2 = I$. Note also that if $P$ is the orthogonal projection onto $\boldsymbol{u}$ then $H = I - 2P$ and $H$ is the reflection through the hyperplane orthogonal to $\boldsymbol{u}$.
```

```{div} definition
Let $\{ \boldsymbol{e}_1, \dots, \boldsymbol{e}_n \}$ be the standard orthonormal basis of $\mathbb{R}^n$

$$
\boldsymbol{e}_1 = \begin{bmatrix} 1 \\ 0 \\ \vdots \\ 0 \end{bmatrix}
\hspace{5mm}
\boldsymbol{e}_2 = \begin{bmatrix} 0 \\ 1 \\ \vdots \\ 0 \end{bmatrix}
\hspace{5mm}
\cdots
\hspace{5mm}
\boldsymbol{e}_n = \begin{bmatrix} 0 \\ 0 \\ \vdots \\ 1 \end{bmatrix}
$$
```

```{div} theorem
Let $\boldsymbol{a} \in \mathbb{R}^n$ and $\alpha = -\mathrm{sign}(a_1) \| \boldsymbol{a} \|$. Let  $\boldsymbol{u} = \boldsymbol{a} - \alpha \boldsymbol{e}_1$, let $P$ be the orthogonal projector onto $\boldsymbol{u}$ and let $H = I - 2P$ be the corresponding elementary reflector. Then

$$
H \boldsymbol{a} = \alpha \boldsymbol{e}_1 = \begin{bmatrix} \alpha \\ 0 \\ \vdots \\ 0 \end{bmatrix}
$$

More generally, let $\boldsymbol{a} \in \mathbb{R}^n$ and partition the vector

$$
\boldsymbol{a} = \begin{bmatrix} \boldsymbol{a}_1 \\ \boldsymbol{a}_2 \end{bmatrix}
\ \ \text{where} \ \
\boldsymbol{a}_1 = \begin{bmatrix} a_1 \\ \vdots \\ a_{k-1} \end{bmatrix} \in \mathbb{R}^{k-1}
 \ \ \text{and} \ \ \boldsymbol{a}_2 = \begin{bmatrix} a_k \\ \vdots \\ a_n \end{bmatrix} \in \mathbb{R}^{n-k+1}
$$

Let $\alpha = -\mathrm{sign}(a_k) \| \boldsymbol{a}_2 \|$ and let

$$
\boldsymbol{u} = \begin{bmatrix} \boldsymbol{0} \\ \boldsymbol{a}_2 \end{bmatrix} - \alpha \boldsymbol{e}_k = \begin{bmatrix} 0 \\ \vdots \\ 0 \\ a_k - \alpha \\ a_{k+1} \\ \vdots \\ a_n \end{bmatrix}
$$

and let $H$ be the corresponding elementary reflector. Then

$$
H \boldsymbol{a} = \begin{bmatrix} a_1 \\ \vdots \\ a_{k-1} \\ \alpha \\ 0 \\ \vdots \\ 0 \end{bmatrix}
$$
```

```{div} theorem
Let $A$ be an $n \times m$ matrix with $n > m$. There exists a sequence of elementary reflectors $H_1,\dots,H_m$ such that $H_m\cdots H_1A = R$ is upper triangular and therefore

$$
A = QR
$$

where $Q = H_1 \cdots H_m$.

---

*Proof*. For each column, construct an elementary reflector to annihilate the entries below the diagonal. For example, if $A$ has 3 columns and 4 rows then

$$
A = \begin{bmatrix} * & * & * \\ * & * & * \\ * & * & * \\ * & * & * \end{bmatrix}
\ \
H_1A = \begin{bmatrix} * & * & * \\ 0 & * & * \\ 0 & * & * \\ 0 & * & * \end{bmatrix}
\ \
H_2H_1A = \begin{bmatrix} * & * & * \\ 0 & * & * \\ 0 & 0 & * \\ 0 & 0 & * \end{bmatrix}
\ \
H_3H_2H_1A = \begin{bmatrix} * & * & * \\ 0 & * & * \\ 0 & 0 & * \\ 0 & 0 & 0 \end{bmatrix}
$$

Since each $H$ is an elementary reflector, we have $A=H_1^{-1}H_2^{-1}H_3^{-1}R=H_1H_2H_3R$.
```

```{div} example
Find the QR decomposition of $A = \begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 1 & 1 & 0 \end{bmatrix}$ by elementary reflectors.

Construct the vector

$$
\boldsymbol{u}_1 =  \boldsymbol{a}_1 - \alpha \boldsymbol{e}_1
= \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix} + \sqrt{2} \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix}
= \begin{bmatrix} 1 + \sqrt{2} \\ 0 \\ 1 \end{bmatrix} \\
$$

Compute the norm squared

$$
\| \boldsymbol{u}_1 \|^2 = (1 + \sqrt{2})^2 + 1 = 4 + 2\sqrt{2}
$$

and construct the elementary reflector

$$
\begin{align*}
H_1 &= I - \frac{2}{\| \boldsymbol{u}_1 \|^2} \boldsymbol{u}_1 \boldsymbol{u}_1^T
= I - \frac{2}{4 + 2\sqrt{2}} \begin{bmatrix} 1 + \sqrt{2} \\ 0 \\ 1 \end{bmatrix} \begin{bmatrix} 1 + \sqrt{2} & 0 & 1 \end{bmatrix} \\
&= I - \frac{1}{2 + \sqrt{2}} \begin{bmatrix} 3 + 2\sqrt{2} & 0 & 1 + \sqrt{2} \\ 0 & 0 & 0 \\ 1 + \sqrt{2} & 0 & 1 \end{bmatrix}
= \frac{1}{2 + \sqrt{2}} \begin{bmatrix} - 1 - \sqrt{2} & 0 & - 1 - \sqrt{2} \\ 0 & 2 + \sqrt{2} & 0 \\ - 1 - \sqrt{2} & 0 & 1 + \sqrt{2} \end{bmatrix} \\
&= \frac{1}{\sqrt{2}} \left[ \begin{array}{rrr} - 1 & 0 & - 1 \\ 0 & \sqrt{2} & 0 \\ - 1 & 0 & 1 \end{array} \right]
\end{align*}
$$

Compute

$$
H_1A
=
\frac{1}{\sqrt{2}} \left[ \begin{array}{rrr} - 1 & 0 & - 1 \\ 0 & \sqrt{2} & 0 \\ - 1 & 0 & 1 \end{array} \right]
\begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 1 & 1 & 0 \end{bmatrix} \\
= \frac{1}{\sqrt{2}} \left[ \begin{array}{rrr} - 2 & -2 & - 1 \\ 0 & \sqrt{2} & \sqrt{2} \\ 0 & 0 & -1 \end{array} \right]
$$

Since the result is already upper triangular we have $A = QR$ where

$$
Q = \frac{1}{\sqrt{2}} \left[ \begin{array}{rrr} - 1 & 0 & - 1 \\ 0 & \sqrt{2} & 0 \\ - 1 & 0 & 1 \end{array} \right]
\hspace{10mm}
R = \frac{1}{\sqrt{2}} \left[ \begin{array}{rrr} - 2 & -2 & - 1 \\ 0 & \sqrt{2} & \sqrt{2} \\ 0 & 0 & -1 \end{array} \right]
$$
```

## Exercises

**Exercise 2.** Find the elementary reflector $H$ corresponding to the vector

$$
\boldsymbol{u} = \begin{bmatrix} 2 \\ 1 \\ 0 \\ 1 \end{bmatrix}
$$

**Exercise 3.** Find an elementary reflector $H$ such that

$$
HA = \displaystyle \begin{bmatrix} * & * & * \\ 0 & * & * \\ 0 & * & * \\ 0 & * & * \end{bmatrix}
$$

where

$$
A = \begin{bmatrix} 1 & 2 & 1 \\ 1 & 0 & 1 \\ 1 & 2 & 1 \\ 1 & 1 & 0 \end{bmatrix}
$$