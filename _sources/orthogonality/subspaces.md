# Subspaces

## Subspaces and Bases

```{div} bigidea
Subspaces of $\mathbb{R}^n$ include lines, planes and hyperplanes through the origin. A basis of a subspace is a linearly independent set of spanning vectors.
```

```{div} definition
A subset $U \subseteq \mathbb{R}^n$ is a **subspace** if:

1. $U$ contains the zero vector $\boldsymbol{0}$
2. $\boldsymbol{u}_1 + \boldsymbol{u}_2 \in U$ for all $\boldsymbol{u}_1,\boldsymbol{u}_2 \in U$
3. $c \boldsymbol{u} \in U$ for all $c \in \mathbb{R},\boldsymbol{u} \in U$
```

```{div} example
<p>

1. The zero subspace $\{ \boldsymbol{0} \}$ and the entire space $\mathbb{R}^n$ are both subspaces of $\mathbb{R}^n$.
2. Subspaces of $\mathbb{R}^2$ include any line through the origin.
3. Subspaces of $\mathbb{R}^3$ include any line or plane through the origin.
4. In general, subspaces of $\mathbb{R}^n$ are hyperplanes of any dimension through the origin.

</p>
```

```{div} definition
A **linear combination** of vectors $\boldsymbol{u}_1,\dots,\boldsymbol{u}_m \in \mathbb{R}^n$ is a vector

$$
c_1 \boldsymbol{u}_1 + \cdots + c_m \boldsymbol{u}_m
$$

where $c_1,\dots,c_m \in \mathbb{R}$. The **span** of vectors $\boldsymbol{u}_1,\dots,\boldsymbol{u}_m \in \mathbb{R}^n$ is the set of all linear combinations

$$
\mathrm{span} \{ \boldsymbol{u}_1 , \dots , \boldsymbol{u}_m \} = \{ c_1 \boldsymbol{u}_1 + \cdots + c_m \boldsymbol{u}_m \in \mathbb{R}^n : c_1,\dots,c_m \in \mathbb{R} \}
$$
```

```{div} theorem
Let $\boldsymbol{u}_1 , \dots , \boldsymbol{u}_m \in \mathbb{R}^n$. Then $\mathrm{span} \{ \boldsymbol{u}_1 , \dots , \boldsymbol{u}_m \}$ is a subspace of $\mathbb{R}^n$.
```

```{div} example
The span of a single nonzero vector $\boldsymbol{u}$ is a line with direction $\boldsymbol{u}$. The span of two nonzero vectors $\boldsymbol{u}$ and $\boldsymbol{v}$ is a plane as long as $\boldsymbol{u}$ and $\boldsymbol{v}$ are not colinear.
```

```{div} definition
A set of vectors $\{ \boldsymbol{u}_1,\dots,\boldsymbol{u}_m \} \subset \mathbb{R}^n$ forms a **linearly independent** set if the vectors satisfy the property:

$$
c_1 \boldsymbol{u}_1 + \cdots + c_m \boldsymbol{u}_m = \boldsymbol{0} \hspace{5mm} \text{if and only if} \hspace{5mm} c_1 = \cdots = c_m = 0
$$

In other words, $\{ \boldsymbol{u}_1,\dots,\boldsymbol{u}_m \}$ is a linearly independent set if no vector $\boldsymbol{u}_i$ in the set can be expressed as a linear combination of the others.
```

```{div} definition
Let $U \subseteq \mathbb{R}^n$ be a subspace. A set of vectors $\{ \boldsymbol{u}_1 , \dots , \boldsymbol{u}_m \}$ forms a **basis** of $U$ if:

1. $\{ \boldsymbol{u}_1 , \dots , \boldsymbol{u}_m \}$ is a linearly independent set
2. $\mathrm{span} \{ \boldsymbol{u}_1 , \dots , \boldsymbol{u}_m \} = U$

The **dimension** of $U$ is the number $m$ of vectors in a basis.
```

```{div} note
How do we know if a set of vectors $\{ \boldsymbol{u}_1,\dots,\boldsymbol{u}_m \}$ is linearly independent? Create a matrix where the columns are the given vectors

$$
A = \begin{bmatrix} & & \\ \boldsymbol{u_1} & \cdots & \boldsymbol{u_m} \\ & & \end{bmatrix}
$$

Then $\{ \boldsymbol{u}_1,\dots,\boldsymbol{u}_m \}$ is a linearly independent set if and only if the linear system $A \boldsymbol{x} = \boldsymbol{0}$ has only the trivial solution $\boldsymbol{x} = \boldsymbol{0}$.
```

## Orthogonal Vectors

```{div} bigidea
Vectors $\boldsymbol{x},\boldsymbol{y} \in \mathbb{R}^n$ are orthogonal if $\langle \boldsymbol{x} , \boldsymbol{y} \rangle = 0$. Subspaces $S_1 \subseteq \mathbb{R}^n$ and $S_2 \subseteq \mathbb{R}^n$ are orthogonal if $\langle \boldsymbol{x}_1 , \boldsymbol{x}_2 \rangle = 0$ for all $\boldsymbol{x}_1 \in S_1$ and $\boldsymbol{x}_2 \in S_2$.
```

```{div} definition
The **inner product** of vectors $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^n$ is

$$
\langle \boldsymbol{x} , \boldsymbol{y} \rangle = \sum_{k=1}^n x_k y_k = x_1y_1 + \cdots + x_ny_n
$$
```

```{div} note
<p>

* The inner product is symmetric: $\langle \boldsymbol{x} , \boldsymbol{y} \rangle = \langle \boldsymbol{y} , \boldsymbol{x} \rangle$ for all $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^n$
* The inner product of column vectors is the same as matrix multiplication:

  $$
  \langle \boldsymbol{x} , \boldsymbol{y} \rangle = \boldsymbol{x}^T \boldsymbol{y} =
  \begin{bmatrix} x_1 & \cdots & x_n \end{bmatrix} \begin{bmatrix} y_1 \\ \vdots \\ y_n \end{bmatrix}
  $$

* The inner product satisfies the usual distributive rules of multiplication:

  $$
  \langle \boldsymbol{x} ,  c \boldsymbol{y} + d \boldsymbol{z} \rangle = c \langle \boldsymbol{x} , \boldsymbol{y} \rangle + d \langle \boldsymbol{x} , \boldsymbol{z} \rangle
  $$

  for all $c,d \in \mathbb{R}$ and $\boldsymbol{x} , \boldsymbol{y} , \boldsymbol{z} \in \mathbb{R}^n$.

* The square root of the inner product of a vector $\boldsymbol{x}$ with itself is equal to the 2-norm

  $$
  \sqrt{ \langle \boldsymbol{x} , \boldsymbol{x} \rangle } = \| \boldsymbol{x} \|
  $$

* We can also write the inner product in terms of the angle between vectors

  $$
  \langle \boldsymbol{x} , \boldsymbol{y} \rangle = \| \boldsymbol{x} \| \| \boldsymbol{y} \| \cos \theta \hspace{10mm} 0 \leq \theta \leq \pi
  $$

* Let $A$ be a $m \times n$ matrix, let $\boldsymbol{u} \in \mathbb{R}^n$ and let $\boldsymbol{v} \in \mathbb{R}^m$. Then

  $$
  \langle A \boldsymbol{u} , \boldsymbol{v} \rangle = \langle \boldsymbol{u} , A^T \boldsymbol{v} \rangle
  $$

</p>
```

```{div} definition
Vectors $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^n$ are **orthogonal** if $\langle \boldsymbol{x} , \boldsymbol{y} \rangle = 0$. More generally, vectors $\boldsymbol{x}_1, \dots, \boldsymbol{x}_m \in \mathbb{R}^n$ are **orthogonal** if $\langle \boldsymbol{x}_i , \boldsymbol{x}_j \rangle = 0$ for all $i \not= j$. In other words, each $\boldsymbol{x}_i$ is orthogonal to every other vector $\boldsymbol{x}_j$ in the set. Furthermore, vectors $\boldsymbol{x}_1, \dots, \boldsymbol{x}_m \in \mathbb{R}^n$ are **orthonormal** if they are orthogonal and each is a unit vector, $\| \boldsymbol{x}_k \| = 1$, $k=1,\dots,m$.
```

```{div} note
Vectors $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{R}^n$ are orthogonal if and only if the acute angle between $\boldsymbol{x}$ and $\boldsymbol{y}$ is $\pi/2$ radians (or 90 degrees).
```

```{div} theorem
Let $\boldsymbol{x}_1, \dots, \boldsymbol{x}_m \in \mathbb{R}^n$ be orthogonal. Then

$$
\| \boldsymbol{x}_1 + \cdots + \boldsymbol{x}_m \|^2 = \| \boldsymbol{x}_1 \|^2 + \cdots + \| \boldsymbol{x}_m \|^2
$$

This is called the **Pythagorean theorem**.

---

*Proof*. Compute the left side of the equation using orthogonality $\langle \boldsymbol{x}_i , \boldsymbol{x}_i \rangle = 0$ if $i \not= j$

$$
\begin{align*}
\| \boldsymbol{x}_1 + \cdots + \boldsymbol{x}_m \|^2 &= \langle \boldsymbol{x}_1 + \cdots + \boldsymbol{x}_m ,\boldsymbol{x}_1 + \cdots + \boldsymbol{x}_m \rangle \\
&= \sum_{i=1}^m \sum_{j=1}^m \langle \boldsymbol{x}_i , \boldsymbol{x}_j \rangle \\
&= \sum_{i=1}^m \langle \boldsymbol{x}_i , \boldsymbol{x}_i \rangle \\
&= \| \boldsymbol{x}_1 \|^2 + \cdots + \| \boldsymbol{x}_m \|^2
\end{align*}
$$

```

## Orthogonal Subspaces

```{div} definition
Let $S_1 \subseteq \mathbb{R}^n$ and $S_2 \subseteq \mathbb{R}^n$ be subspaces. Then $S_1$ and $S_2$ are **orthogonal** if $\langle \boldsymbol{x}_1 , \boldsymbol{x}_2 \rangle = 0$ for all $\boldsymbol{x}_1 \in S_1$ and $\boldsymbol{x}_2 \in S_2$. If $S_1$ and $S_2$ are orthogonal subspaces, then we write $S_1 \perp S_2$.
```

```{div} theorem
Let $\{ \boldsymbol{u}_1,\dots,\boldsymbol{u}_m \}$ be a basis of a subspace $S_1 \subseteq \mathbb{R}^n$ and let $\{ \boldsymbol{v}_1,\dots,\boldsymbol{v}_{\ell} \}$ be a basis of a subspace $S_2 \subseteq \mathbb{R}^n$. Then $S_1 \perp S_2$ if and only if $\langle \boldsymbol{u}_i , \boldsymbol{v}_j \rangle = 0$ for all $i,j$. In other words, every $\boldsymbol{u}_i$ in the basis of $S_1$ is orthogonal to each $\boldsymbol{v}_j$ in the basis of $S_2$.
```

```{div} example
Let $S_1 \subset \mathbb{R}^3$ and $S_2 \subset \mathbb{R}^3$ be 2-dimensional subspaces (planes). Is it possible that $S_1 \perp S_2$? No!
```

```{div} definition
Let $S \subseteq \mathbb{R}^n$ be a subspace. The **orthogonal complement of $S$** is

$$
S^{\perp} = \{ \boldsymbol{x} \in \mathbb{R}^n : \langle \boldsymbol{x} , \boldsymbol{y} \rangle = 0 \text{ for all } \boldsymbol{y} \in S \}
$$
```

```{div} note
<p>

* If $S \subseteq \mathbb{R}^n$ is any subspace then $S = (S^{\perp})^{\perp}$ and also $S \cap S^{\perp} = \{ \boldsymbol{0} \}$.
* $\{ \boldsymbol{0} \}^{\perp} = \mathbb{R}^n$.

</p>

```

```{div} theorem
Let $S \subseteq \mathbb{R}^n$ be a subspace. Then $S^{\perp} \subseteq \mathbb{R}^n$ is a subspace.

---

*Proof*. Let us verify that $S^{\perp}$ satisfies the properties of a subspace.

Clearly $\langle \boldsymbol{0} , \boldsymbol{x} \rangle = 0$ for all $\boldsymbol{x} \in S$ therefore $\boldsymbol{0} \in S^{\perp}$.

Let $\boldsymbol{x}_1,\boldsymbol{x}_2 \in S^{\perp}$. Then

$$
\langle \boldsymbol{x}_1 + \boldsymbol{x}_2) , \boldsymbol{y} \rangle = \langle \boldsymbol{x}_1 , \boldsymbol{y} \rangle + \langle \boldsymbol{x}_2 , \boldsymbol{y} \rangle = 0 + 0 = 0
$$

for all $\boldsymbol{y} \in S$ therefore $\boldsymbol{x}_1 + \boldsymbol{x}_2 \in S^{\perp}$.

Let $c \in \mathbb{R}$ and $\boldsymbol{x} \in S^{\perp}$. Then

$$
\langle c\boldsymbol{x} , \boldsymbol{y} \rangle = c \langle \boldsymbol{x} , \boldsymbol{y} \rangle = c(0) = 0
$$

for all $\boldsymbol{y} \in S$ therefore $c \boldsymbol{x} \in S^{\perp}$.

Therefore $S^{\perp}$ is a subspace.
```

```{div} theorem
Let $S \subseteq \mathbb{R}^n$ be a subspace. Then

$$
\dim(S) + \dim(S^{\perp}) = n
$$
```

## Fundamental Subspaces of a Matrix

```{div} bigidea
Any $m \times n$ matrix $A$ naturally decomposes $\mathbb{R}^n$ into $N(A)$ and $R(A^T)$ and $\mathbb{R}^m$ into $N(A^T)$ and $R(A)$.
```

```{div} definition
Let $A$ be a $m \times n$ matrix. The four **fundamental subspaces** of $A$ are:

$$
\begin{align*}
N(A) &= \{ \boldsymbol{x} \in \mathbb{R}^n : A \boldsymbol{x} = \boldsymbol{0} \} \\
R(A) &= \{ \boldsymbol{y} \in \mathbb{R}^m : \text{there exists } \boldsymbol{x} \in \mathbb{R}^n \text{ such that } A \boldsymbol{x} = \boldsymbol{y} \} \\
N(A^T) &= \{ \boldsymbol{y} \in \mathbb{R}^m : A^T \boldsymbol{y} = \boldsymbol{0} \} \\
R(A^T) &= \{ \boldsymbol{x} \in \mathbb{R}^n : \text{there exists } \boldsymbol{y} \in \mathbb{R}^m \text{ such that } A^T \boldsymbol{y} = \boldsymbol{x} \}
\end{align*}
$$

The subspace $N(A)$ is the **nullspace** of $A$, $N(A^T)$ is the nullspace of $A^T$, $R(A)$ is the **range** of $A$ (or the **image** or **column space** of $A$), and $R(A^T)$ is the range of $A^T$ (or the **row space** of $A$).
```

```{div} theorem
Let $A$ be a $m \times n$ matrix. The nullspace $N(A)$ is a subspace of $\mathbb{R}^n$ and the range $R(A)$ is a subspace of
$\mathbb{R}^m$.
```

```{div} theorem
Let $A$ be a $m \times n$ matrix. Then $N(A) = R(A^T)^{\perp}$ and $R(A) = N(A^T)^{\perp}$.

---

*Proof*. The second equality follows from the first by replacing $A$ with $A^T$ therefore it is sufficient to prove $N(A) = R(A^T)^{\perp}$. A general strategy to prove equality of sets is to show that each set contains the other therefore let us prove  $N(A) \subseteq R(A^T)^{\perp}$ and then prove the reverse $R(A^T)^{\perp} \subseteq N(A)$.

Let $\boldsymbol{x} \in N(A)$. Then $A \boldsymbol{x} = \boldsymbol{0}$ therefore $\langle A \boldsymbol{x} , \boldsymbol{y} \rangle = 0$ for all $\boldsymbol{y} \in \mathbb{R}^m$. Using properties of the inner product we see that $\langle \boldsymbol{x} , A^T \boldsymbol{y} \rangle = 0$ for all $\boldsymbol{y} \in \mathbb{R}^m$ therefore $\boldsymbol{x} \in R(A^T)^{\perp}$.

Now let $\boldsymbol{x} \in R(A^T)^{\perp}$. Then $\langle \boldsymbol{x} , A^T \boldsymbol{y} \rangle = 0$ and so $\langle A \boldsymbol{x} , \boldsymbol{y} \rangle = 0$ for all $\boldsymbol{y} \in \mathbb{R}^m$. Choose $\boldsymbol{y} = A\boldsymbol{x} \in \mathbb{R}^m$ and then $\langle A \boldsymbol{x} , A \boldsymbol{x} \rangle = 0$. Therefore $\| A \boldsymbol{x} \| = 0$ and so $A \boldsymbol{x} = \boldsymbol{0}$ and finally $\boldsymbol{x} \in N(A)$.

Since $N(A) \subseteq R(A^T)^{\perp}$ and $R(A^T)^{\perp} \subseteq N(A)$ we have $N(A) = R(A^T)^{\perp}$.
```

```{div} example
Let $A$ be a matrix such that its LU decomposition is of the form

$$
A = LU =
\begin{bmatrix} 1 & 0 & 0 \\ * & 1 & 0 \\ * & * & 1 \end{bmatrix}
\begin{bmatrix} * & * & * & * \\ 0 & * & * & * \\ 0 & 0 & * & * \end{bmatrix}
$$

where $*$ denotes a nonzero number. Find the dimension of each subspace $N(A)$, $R(A)$, $R(A^T)$ and $N(A^T)$.
```

## Exercises
1. Let $A = LU$ be the LU decomposition of $A$. Determine whether the statement is **True** or **False**.
   * $N(A) = N(U)$
   * $N(A^T) = N(U^T)$
   * $R(A) = R(U)$
   * $R(A^T) = R(U^T)$
   * $\mathrm{dim} (R(A)) = \mathrm{dim} (R(U))$
2. Let $A$ be a $m \times n$ matrix and let $\{ \boldsymbol{u}_1,\boldsymbol{u}_2 \} \subset \mathbb{R}^n$ be a basis of the nullspace $N(A)$. Determine $\mathrm{dim}(R(A^T))$ and $\mathrm{dim}(N(A^T))$.
3. Let $A$ be a $4 \times 4$ matrix such that

   $$
   A = LU =
   \left[ \begin{array}{rrrrrr} 1 & 0 & 0 & 0 \\ 1 & 1 & 0 & 0 \\ 0 &  1 & 1 & 0 \\ 0 & 2 & 1 & 1 \end{array} \right]
   \left[ \begin{array}{rrrrrr} 1 & -1 & 2 & -1 \\ 0 & 1 & -3 & 4 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 \end{array} \right]
   $$

   Find a basis of $N(A^T)$ and find a basis of $R(A^T)$.
4. Let $A$ be a matrix such that its LU decomposition is of the form

   $$
   A = LU = \begin{bmatrix} 1 & 0 & 0 \\ * & 1 & 0 \\ * & * & 1 \end{bmatrix}
   \begin{bmatrix} * & * & * & * \\ 0 & * & * & * \\ 0 & 0 & 0 & * \end{bmatrix}
   $$

   where $*$ denotes a nonzero number. Determine the dimension of $R(A^T)$ and the dimension of $N(A^T)$.
5. Determine whether the statement is **True** or **False**.
   * Let $S \subseteq \mathbb{R}^n$ be a subspace. If $\boldsymbol{u} \in \mathbb{R}^n$ such that $\boldsymbol{u} \not= 0$ then either $\boldsymbol{u} \in S$ or $\boldsymbol{u} \in S^{\perp}$.
   * Let $L_1 \subset \mathbb{R}^2$ be a line through the origin. There is a unique line $L_2$ through the origin such that $L_1 \perp L_2$.
   * Let $L_1 \subset \mathbb{R}^3$ be a line through the origin. There is a unique line $L_2$ through the origin such that $L_1 \perp L_2$.
   * Let $U_1 \subset \mathbb{R}^4$ be a 2-dimensional subspace. There is a unique plane $U_2$ through the origin such that $U_1 \perp U_2$.
