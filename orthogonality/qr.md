# QR Decomposition

```{div} bigidea
The QR decomposition of a matrix $A$ is $A = QR$ where $Q$ is an orthogonal matrix and $R$ is an upper triangular matrix. There are several ways to compute the QR decomposition including Gram-Schmidt orthogonalization and elementary reflectors.
```

## Orthogonal Matrices

```{div} definition
A matrix $A$ is **orthogonal** if $A^TA = AA^T = I$.
```

```{div} note
 If $A$ is an orthogonal matrix, then:

* $\| A \boldsymbol{x} \| = \| \boldsymbol{x} \|$ for all $\boldsymbol{x} \in \mathbb{R}^n$ since $\| A \boldsymbol{x} \|^2 = (A \boldsymbol{x})^T A \boldsymbol{x} = \boldsymbol{x}^T A^T A \boldsymbol{x} = \boldsymbol{x}^T \boldsymbol{x} = \| \boldsymbol{x} \|^2$
* the columns of $A$ are orthonormal
* the rows of $A$ are orthonormal
```

```{div} example
Rotations and reflections are examples of orthogonal matrices.
```

```{div} note
An orthogonal matrix and an orthogonal projector are *not* the same thing but they are related. If $P$ is an orthogonal projector then $Q = I - 2P$ is an orthogonal (and symmetric) matrix. In fact, if $P$ projects onto a subspace $U$ then $Q$ is the reflection through $U^{\perp}$.
```

## QR by Gram-Schmidt

```{div} definition
Let $A$ be an $n \times m$ matrix with $\mathrm{rank}(A) = m$ and let $\boldsymbol{a}_1,\dots,\boldsymbol{a}_m$ be the columns of $A$. Apply the Gram-Schmidt algorithm to the columns and construct an orthonormal basis $\{ \boldsymbol{w}_1,\dots,\boldsymbol{w}_m \}$ of the column space $R(A)$. Project the columns onto the basis

$$
\begin{align*}
\boldsymbol{a}_1 &= \langle \boldsymbol{w}_1 , \boldsymbol{a}_1 \rangle \boldsymbol{w}_1 \\
\boldsymbol{a}_2 &= \langle \boldsymbol{w}_1 , \boldsymbol{a}_2 \rangle \boldsymbol{w}_1 + \langle \boldsymbol{w}_2 , \boldsymbol{a}_2 \rangle \boldsymbol{w}_2 \\
& \ \ \vdots \\
\boldsymbol{a}_m &= \langle \boldsymbol{w}_1 , \boldsymbol{a}_m \rangle \boldsymbol{w}_1 + \langle \boldsymbol{w}_2 , \boldsymbol{a}_m \rangle \boldsymbol{w}_2 + \cdots + \langle \boldsymbol{w}_m , \boldsymbol{a}_m \rangle \boldsymbol{w}_m
\end{align*}
$$

where $\boldsymbol{a}_k \in \mathrm{span} \{ \boldsymbol{w}_1 , \dots , \boldsymbol{w}_k \}$ by construction. Write as matrix multiplication

$$
A = Q_1R_1
$$

where

$$
Q_1 = \begin{bmatrix} & & \\ \boldsymbol{w}_1 & \cdots & \boldsymbol{w}_m \\ & & \end{bmatrix}
\hspace{5mm}
R_1 = \begin{bmatrix}
\langle \boldsymbol{w}_1 , \boldsymbol{a}_1 \rangle & \langle \boldsymbol{w}_1 , \boldsymbol{a}_2 \rangle & \cdots & \langle \boldsymbol{w}_1 , \boldsymbol{a}_m \rangle \\
& \langle \boldsymbol{w}_2 , \boldsymbol{a}_2 \rangle & \cdots & \langle \boldsymbol{w}_2 , \boldsymbol{a}_m \rangle \\
& & \ddots & \vdots \\
& & & \langle \boldsymbol{w}_m , \boldsymbol{a}_m \rangle
\end{bmatrix}
$$

Extend the basis to an orthonormal basis $\{ \boldsymbol{w}_1 , \dots , \boldsymbol{w}_m , \boldsymbol{w}_{m+1} , \dots , \boldsymbol{w}_n \}$ of $\mathbb{R}^n$ where $\{ \boldsymbol{w}_{m+1} , \dots , \boldsymbol{w}_n \}$ is *any* orthonormal basis of the orthogonal complement $R(A)^{\perp}$ and let

$$
Q_2 = \begin{bmatrix} & & \\ \boldsymbol{w}_{m+1} & \cdots & \boldsymbol{w}_n \\ & & \end{bmatrix}
$$

Finally, the **QR decomposition** of $A$ is

$$
A = QR =
\begin{bmatrix} Q_1 & Q_2 \end{bmatrix}
\begin{bmatrix} R_1 \\ 0 \end{bmatrix}
$$

where $Q$ is a $n \times n$ orthogonal matrix and $R$ is a $n \times m$ upper triangular matrix. The decomposition $A = Q_1 R_1$ is called the **thin QR decomposition**. See [Wikipedia:QR decomposition](https://en.wikipedia.org/wiki/QR_decomposition).
```

```{div} example
Compute the QR decomposition for the matrix

$$
A = \begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\  1 & 1 & 0 \\  0 & 0 & 0 \\ \end{bmatrix}
$$

In a previous example, we found an orthonormal basis of the column space

$$
\boldsymbol{w}_1 = \frac{1}{\sqrt{2}} \begin{bmatrix} 1 \\ 0 \\ 1 \\ 0 \end{bmatrix}
\hspace{5mm}
\boldsymbol{w}_2 = \begin{bmatrix} 0 \\ 1 \\ 0 \\ 0 \end{bmatrix}
\hspace{5mm}
\boldsymbol{w}_3 = \frac{1}{\sqrt{2}} \left[ \begin{array}{r} 1 \\ 0 \\ -1 \\ 0 \end{array} \right]
$$

Extend to an orthonormal basis of $\mathbb{R}^4$ by

$$
\boldsymbol{w}_4 = \begin{bmatrix} 0 \\ 0 \\ 0 \\ 1 \end{bmatrix}
$$

Therefore we have

$$
Q =
\begin{bmatrix}
1/\sqrt{2} & 0 & 1/\sqrt{2} & 0 \\
0 & 1 & 0 & 0 \\
1/\sqrt{2} & 0 & -1/\sqrt{2} & 0 \\
0 & 0 & 0 & 1
\end{bmatrix}
$$

and

$$
R =
\begin{bmatrix}
\langle \boldsymbol{w}_1 , \boldsymbol{a}_1 \rangle & \langle \boldsymbol{w}_1 , \boldsymbol{a}_2 \rangle & \langle \boldsymbol{w}_1 , \boldsymbol{a}_3 \rangle \\
0 & \langle \boldsymbol{w}_2 , \boldsymbol{a}_2 \rangle & \langle \boldsymbol{w}_2 , \boldsymbol{a}_3 \rangle \\
0 & 0 & \langle \boldsymbol{w}_3 , \boldsymbol{a}_3 \rangle \\
0 & 0 & 0
\end{bmatrix}
=
\begin{bmatrix}
\sqrt{2} & \sqrt{2} & 1/\sqrt{2} \\
0 & 1 & 1 \\
0 & 0 & 1/\sqrt{2} \\
0 & 0 & 0
\end{bmatrix}
$$
```

```{div} note
The Gram-Schmidt algorithm shows that the QR decomposition exists but it is not the most efficient way to compute the QR decomposition. Software such as the MATLAB function `qr` (see [documentation](https://www.mathworks.com/help/matlab/ref/qr.html)) and the SciPy function `scipy.linalg.qr` (see [documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr.html)), and LAPACK (see [documentation](https://www.netlib.org/lapack/lug/node128.html#secorthog)) use elementary reflectors to construct the matrices $Q$ and $R$.
```

## QR by Reflectors

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

## Orthogonal Projection by QR

```{div} theorem
Let $A = Q_1 R_1$ be the thin QR decomposition of $A$. The orthogonal projection of $\boldsymbol{x}$ onto $R(A)$ is given by

$$
\mathrm{proj}_{R(A)}(\boldsymbol{x}) = Q_1 Q_1^T \boldsymbol{x}
$$
```

## Exercises

1. Determine whether the statement is **True** or **False**.

   * Let $\boldsymbol{u} \in \mathbb{R}^n$ be a nonzero vector and let $\displaystyle H = I - \frac{2}{\| \boldsymbol{u} \|^2}\boldsymbol{u} \boldsymbol{u}^T$ be the corresponding elementary reflector. There is a unique unit vector $\boldsymbol{v}$ such that $H \boldsymbol{v} = \boldsymbol{v}$.
   * Let $\boldsymbol{u} \in \mathbb{R}^n$ be a nonzero vector and let $\displaystyle H = I - \frac{2}{\| \boldsymbol{u} \|^2}\boldsymbol{u} \boldsymbol{u}^T$ be the corresponding elementary reflector. Then  $H \boldsymbol{v} = \boldsymbol{v}$ for all $\boldsymbol{v} \in \mathrm{span} \{ \boldsymbol{u} \}^{\perp}$.
   * Let $A=QR$ where $Q$ is an orthogonal matrix and $R$ is upper triangular. Then $\| A \|_p = \| R \|_p$ for any $p \geq 1$.

2. Find the elementary reflector $H$ corresponding to the vector

   $$
   \boldsymbol{u} = \begin{bmatrix} 2 \\ 1 \\ 0 \\ 1 \end{bmatrix}
   $$

3. Find an elementary reflector $H$ such that

   $$
   HA = \displaystyle \begin{bmatrix} * & * & * \\ 0 & * & * \\ 0 & * & * \\ 0 & * & * \end{bmatrix}
   $$

   where

   $$
   A = \begin{bmatrix} 1 & 2 & 1 \\ 1 & 0 & 1 \\ 1 & 2 & 1 \\ 1 & 1 & 0 \end{bmatrix}
   $$
