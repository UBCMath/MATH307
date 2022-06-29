# QR Decomposition

```{div} bigidea
If $A$ is a $m \times n$ matrix with $\mathrm{rank}(A) = n$, then the decomposition $A = QR$ provides orthonormal bases of both $R(A)$ and $R(A)^{\perp}$. In particular, write $Q = [Q_1 \ Q_2]$ where $Q_1$ is the first $n$ columns of $Q$, then the columns of $Q_1$ provide an orthonormal basis of $R(A)$ and the columns of $Q_2$ provide an orthonormal basis of $R(A)^{\perp}$. 
```

```{image} /img/02_04_01.png
:width: 100%
:align: center
```

## Orthogonal Matrices

```{div} definition
A matrix $A$ is **orthogonal** if $A^TA = AA^T = I$.
```

```{div} note
The condition $A^TA = I$ implies that the columns of $A$ are orthonormal, and $AA^T = I$ implies the rows of $A$ are orthonormal.
```

```{div} theorem
 If $A$ is an orthogonal matrix, then $\| A \boldsymbol{x} \| = \| \boldsymbol{x} \|$ for all $\boldsymbol{x} \in \mathbb{R}^n$.

---

*Proof.* Compute

$$
\| A \boldsymbol{x} \|^2 = (A \boldsymbol{x})^T A \boldsymbol{x} = \boldsymbol{x}^T A^T A \boldsymbol{x} = \boldsymbol{x}^T \boldsymbol{x} = \| \boldsymbol{x} \|^2
$$
```

```{div} example
Rotations and reflections are examples of orthogonal matrices.
```

```{div} note
An orthogonal matrix and an orthogonal projector are *not* the same thing but they are related. If $P$ is an orthogonal projector then $Q = I - 2P$ is an orthogonal (and symmetric) matrix. In fact, if $P$ projects onto a subspace $U$ then $Q$ is the reflection through $U^{\perp}$.
```

## QR by Gram-Schmidt

```{div} definition
Let $A$ be an $m \times n$ matrix with $\mathrm{rank}(A) = n$ and let $\boldsymbol{a}_1,\dots,\boldsymbol{a}_n$ be the columns of $A$. Apply the Gram-Schmidt algorithm to the columns and construct an orthonormal basis $\{ \boldsymbol{w}_1,\dots,\boldsymbol{w}_n \}$ of $R(A)$. Project the columns onto the basis

$$
\begin{align*}
\boldsymbol{a}_1 &= \langle \boldsymbol{w}_1 , \boldsymbol{a}_1 \rangle \boldsymbol{w}_1 \\
\boldsymbol{a}_2 &= \langle \boldsymbol{w}_1 , \boldsymbol{a}_2 \rangle \boldsymbol{w}_1 + \langle \boldsymbol{w}_2 , \boldsymbol{a}_2 \rangle \boldsymbol{w}_2 \\
& \ \ \vdots \\
\boldsymbol{a}_n &= \langle \boldsymbol{w}_1 , \boldsymbol{a}_n \rangle \boldsymbol{w}_1 + \langle \boldsymbol{w}_2 , \boldsymbol{a}_n \rangle \boldsymbol{w}_2 + \cdots + \langle \boldsymbol{w}_n , \boldsymbol{a}_n \rangle \boldsymbol{w}_n
\end{align*}
$$

where $\boldsymbol{a}_k \in \mathrm{span} \{ \boldsymbol{w}_1 , \dots , \boldsymbol{w}_k \}$ by construction. Write as matrix multiplication

$$
A = Q_1R_1
$$

where

$$
Q_1 = \begin{bmatrix} & & \\ \boldsymbol{w}_1 & \cdots & \boldsymbol{w}_n \\ & & \end{bmatrix}
\hspace{5mm}
R_1 = \begin{bmatrix}
\langle \boldsymbol{w}_1 , \boldsymbol{a}_1 \rangle & \langle \boldsymbol{w}_1 , \boldsymbol{a}_2 \rangle & \cdots & \langle \boldsymbol{w}_1 , \boldsymbol{a}_n \rangle \\
& \langle \boldsymbol{w}_2 , \boldsymbol{a}_2 \rangle & \cdots & \langle \boldsymbol{w}_2 , \boldsymbol{a}_n \rangle \\
& & \ddots & \vdots \\
& & & \langle \boldsymbol{w}_n , \boldsymbol{a}_n \rangle
\end{bmatrix}
$$

Extend the basis to an orthonormal basis $\{ \boldsymbol{w}_1 , \dots , \boldsymbol{w}_n , \boldsymbol{w}_{n+1} , \dots , \boldsymbol{w}_m \}$ of $\mathbb{R}^m$ where $\{ \boldsymbol{w}_{n+1} , \dots , \boldsymbol{w}_m \}$ is *any* orthonormal basis of the orthogonal complement $R(A)^{\perp}$ and let

$$
Q_2 = \begin{bmatrix} & & \\ \boldsymbol{w}_{n+1} & \cdots & \boldsymbol{w}_m \\ & & \end{bmatrix}
$$

Finally, the **QR decomposition** of $A$ is

$$
A = QR =
\begin{bmatrix} Q_1 & Q_2 \end{bmatrix}
\begin{bmatrix} R_1 \\ \mathbf{0} \end{bmatrix}
$$

where $Q$ is a $m \times m$ orthogonal matrix and $R$ is a $m \times n$ upper triangular matrix. The decomposition $A = Q_1 R_1$ is called the **thin QR decomposition**. See [Wikipedia:QR decomposition](https://en.wikipedia.org/wiki/QR_decomposition).
```

```{div} example
Compute the QR decomposition for the matrix

$$
A = \begin{bmatrix} 1 & 1 & 1 \\ 0 & 1 & 1 \\ 1 & 1 & 0 \end{bmatrix}
$$

Apply Gram-Schmidt to find an orthonormal basis of the column space

$$
\boldsymbol{w}_1 = \frac{1}{\sqrt{2}} \begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}
\hspace{5mm}
\boldsymbol{w}_2 = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix}
\hspace{5mm}
\boldsymbol{w}_3 = \frac{1}{\sqrt{2}} \left[ \begin{array}{r} 1 \\ 0 \\ -1 \end{array} \right]
$$

Therefore

$$
Q =
\begin{bmatrix}
1/\sqrt{2} & 0 & 1/\sqrt{2} \\
0 & 1 & 0 \\
1/\sqrt{2} & 0 & -1/\sqrt{2}
\end{bmatrix}
$$

and

$$
R =
\begin{bmatrix}
\langle \boldsymbol{w}_1 , \boldsymbol{a}_1 \rangle & \langle \boldsymbol{w}_1 , \boldsymbol{a}_2 \rangle & \langle \boldsymbol{w}_1 , \boldsymbol{a}_3 \rangle \\
0 & \langle \boldsymbol{w}_2 , \boldsymbol{a}_2 \rangle & \langle \boldsymbol{w}_2 , \boldsymbol{a}_3 \rangle \\
0 & 0 & \langle \boldsymbol{w}_3 , \boldsymbol{a}_3 \rangle
\end{bmatrix}
=
\begin{bmatrix}
\sqrt{2} & \sqrt{2} & 1/\sqrt{2} \\
0 & 1 & 1 \\
0 & 0 & 1/\sqrt{2}
\end{bmatrix}
$$
```

```{div} note
The Gram-Schmidt algorithm shows that the QR decomposition exists but it is not the most efficient way to compute the QR decomposition. Software such as the MATLAB function `qr` (see [documentation](https://www.mathworks.com/help/matlab/ref/qr.html)) and the SciPy function `scipy.linalg.qr` (see [documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.qr.html)), and LAPACK (see [documentation](https://www.netlib.org/lapack/lug/node128.html#secorthog)) use elementary reflectors to construct the matrices $Q$ and $R$.
```

## Orthogonal Projection by QR

```{div} theorem
Let $A$ be a $m \times n$ matrix such that $\mathrm{rank}(A) = n$, and let $A = Q R$ be the QR decomposition. The projection of $\boldsymbol{x} \in \mathbb{R}^m$ onto $R(A)$ is given by

$$
\mathrm{proj}_{R(A)}(\boldsymbol{x}) = Q_1 Q_1^T \boldsymbol{x}
$$

and the projection onto $R(A)^{\perp}$ is given by

$$
\mathrm{proj}_{R(A)^{\perp}}(\boldsymbol{x}) = Q_2 Q_2^T \boldsymbol{x}
$$
```

## Exercises

````{div} exercise
Compute the thin QR decomposition of the matrix

$$
A = \left[ \begin{array}{rr} 1 & 1 \\ 1 & -1 \\ 1 & 1 \\ 1 & 1 \end{array} \right]
$$

```{dropdown} Solution
Apply Gram-Schmidt to the columns of $A$

$$
\boldsymbol{v}_1 = \boldsymbol{a}_1
\hspace{10mm}
\boldsymbol{v}_2 = \boldsymbol{a}_2 - \frac{\langle \boldsymbol{v}_1 , \boldsymbol{a}_2 \rangle}{\langle \boldsymbol{v}_1 , \boldsymbol{v}_1 \rangle} \boldsymbol{v}_1 = \frac{1}{2} \left[ \begin{array}{r} 1 \\ -3 \\ 1 \\ 1 \end{array} \right]
$$

Therefore

$$
Q_1 = \begin{bmatrix} & \\ \boldsymbol{w}_1 & \boldsymbol{w}_2 \\ & \end{bmatrix} = \left[ \begin{array}{rr} 1/2 & 1/\sqrt{12} \\ 1/2 & -3/\sqrt{12} \\ 1/2 & 1/\sqrt{12} \\ 1/2 & 1/\sqrt{12} \end{array} \right]
$$

and

$$
R_1 =
\begin{bmatrix}
\langle \boldsymbol{w}_1 , \boldsymbol{a}_1 \rangle & \langle \boldsymbol{w}_1 , \boldsymbol{a}_2 \rangle \\
0 & \langle \boldsymbol{w}_2 , \boldsymbol{a}_2 \rangle
\end{bmatrix}
=
\begin{bmatrix}
2 & 1 \\
0 & \sqrt{3}
\end{bmatrix}
$$

```
````