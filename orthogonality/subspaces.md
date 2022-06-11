# Subspaces

```{div} bigidea
Subspaces of $\mathbb{R}^n$ include lines, planes and hyperplanes through the origin. A basis of a subspace is a linearly independent set of spanning vectors. The Rank-Nullity Theorem describes the dimensions of the nullspace and range of a matrix.
```

## Subspaces

```{div} definition
A subset $U \subseteq \mathbb{R}^n$ is a **subspace** if:

1. $U$ contains the zero vector $\boldsymbol{0}$
2. $\boldsymbol{u}_1 + \boldsymbol{u}_2 \in U$ for all $\boldsymbol{u}_1,\boldsymbol{u}_2 \in U$
3. $c \boldsymbol{u} \in U$ for all $c \in \mathbb{R},\boldsymbol{u} \in U$
```

```{div} note
<p>

* Condition 2 is called **closed under addition**
* Condition 3 is called **closed under scalar multiplication**
* Condition 3 with $c=0$ implies Condition 1
</p>
```

```{div} example
<p>

1. The zero subspace $\{ \boldsymbol{0} \}$ and the entire space $\mathbb{R}^n$ are both subspaces of $\mathbb{R}^n$.
2. Subspaces of $\mathbb{R}^2$ include any line through the origin.
3. Subspaces of $\mathbb{R}^3$ include any line or plane through the origin.
4. In general, subspaces of $\mathbb{R}^n$ are hyperplanes of any dimension through the origin.

</p>
```

```{div} example
Consider the set

$$
U = \left\{ \begin{bmatrix} x \\ y \end{bmatrix} : y \geq 0 \right\}
$$

Then $U$ contains the zero vector

$$
\boldsymbol{0} = \begin{bmatrix} 0 \\ 0 \end{bmatrix} \in U
$$

and $U$ is closed under addition since

$$
\boldsymbol{u}_1 + \boldsymbol{u}_2 = \begin{bmatrix} x_1 \\ y_1 \end{bmatrix} + \begin{bmatrix} x_2 \\ y_2 \end{bmatrix} = \begin{bmatrix} x_1 + x_2 \\ y_1 + y_2 \end{bmatrix} \in U
$$

because $y_1 + y_2 \geq 0$ since $y_1 \geq 0$ and $y_2 \geq 0$. However $U$ is not closed under scalar multiplication because

$$
(-1) \begin{bmatrix} 1 \\ 1 \end{bmatrix} = \begin{bmatrix} -1 \\ -1 \end{bmatrix} \not\in U
$$

Therefore $U$ is not a subspace of $\mathbb{R}^2$.
```

## Linear Independence and Span

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

In other words, $\{ \boldsymbol{u}_1,\dots,\boldsymbol{u}_m \}$ is a linearly independent set if no vector in the set can be expressed as a linear combination of the others.
```

```{div} note
How do we know if a set of vectors $\{ \boldsymbol{u}_1,\dots,\boldsymbol{u}_m \}$ is linearly independent? Create a matrix where the columns are the given vectors

$$
A = \begin{bmatrix} & & \\ \boldsymbol{u_1} & \cdots & \boldsymbol{u_m} \\ & & \end{bmatrix}
$$

Then $\{ \boldsymbol{u}_1,\dots,\boldsymbol{u}_m \}$ is a linearly independent set if and only if the linear system $A \boldsymbol{x} = \boldsymbol{0}$ has only the trivial solution $\boldsymbol{x} = \boldsymbol{0}$.
```

## Basis and Dimension

```{div} definition
Let $U \subseteq \mathbb{R}^n$ be a subspace. A set of vectors $\{ \boldsymbol{u}_1 , \dots , \boldsymbol{u}_m \}$ forms a **basis** of $U$ if:

1. $\{ \boldsymbol{u}_1 , \dots , \boldsymbol{u}_m \}$ is a linearly independent set
2. $\mathrm{span} \{ \boldsymbol{u}_1 , \dots , \boldsymbol{u}_m \} = U$

The **dimension** of $U$ is the number $m$ of vectors in a basis.
```

## Nullspace

```{div} definition
The **nullspace** of a $m \times n$ matrix $A$ is

$$
N(A) = \{ \boldsymbol{x} \in \mathbb{R}^n : A\boldsymbol{x} = \boldsymbol{0} \}
$$
```

```{div} theorem
Let $A$ be a $m \times n$ matrix. The nullspace $N(A)$ is a subspace of $\mathbb{R}^n$.
```

```{div} theorem
Let $A$ be a $m \times n$ matrix and let $A = LU$ be the LU decomposition (if it exists). Then $N(A) = N(U)$.
```

## Range

```{div} definition
The **range** of a $m \times n$ matrix $A$ is:

$$
R(A) = \{ A \boldsymbol{x} : \boldsymbol{x} \in \mathbb{R}^n \}
$$

The range of $A$ is also called **image** or **column space** of $A$.
```

```{div} note
Matrix multiplication can be written as

$$
A \boldsymbol{x} = \begin{bmatrix} & & \\ \boldsymbol{a}_1 & \cdots & \boldsymbol{a}_n \\ & & \end{bmatrix} \begin{bmatrix} x_1 \\ \vdots \\ x_n \end{bmatrix} = x_1 \boldsymbol{a}_1 + \cdots + x_n \boldsymbol{a}_n
$$

Therefore the range of $A$ is the equal to the span of the columns

$$
R(A) = \mathrm{span} \{ \boldsymbol{a}_1 , \dots, \boldsymbol{a}_n \}
$$

and that's why $R(A)$ is sometimes called the column space.
```

```{div} theorem
Let $A$ be a $m \times n$ matrix. The range $R(A)$ is a subspace of $\mathbb{R}^m$.
```

```{div} theorem
Let $A$ be an $m \times n$ matrix. Then

$$
\dim (R(A)) = \mathrm{rank}(A)
$$

---

*Proof*. The rank of $A$ is the number of nonzero rows in the row echelon form of $A$. The dimension of $A$ is the number of linearly independent columns in $A$ which is also equal to the number of nonzero rows in $A$.
```

```{div} theorem
Let $A = LU$ be the LU decomposition of $A$ (if it exists) and let $r = \mathrm{rank}(A)$. Then

$$
R(A) = \mathrm{span} \{ \boldsymbol{\ell}_1 , \dots , \boldsymbol{\ell}_r \}
$$

where $\boldsymbol{\ell}_1 , \dots , \boldsymbol{\ell}_r$ are the first $r$ columns of $L$. In particular, $\boldsymbol{\ell}_1 , \dots , \boldsymbol{\ell}_r$ is a basis of $R(A)$.

---

*Proof*. If $\mathrm{rank}(A) = r$ then only the first $r$ entries of the vector $U \boldsymbol{x}$ are nonzero

$$
U \boldsymbol{x} = \begin{bmatrix} * & * & \cdots & * \\ 0 & \ddots & \ddots & \vdots \\ \vdots & \ddots & * & * \\ 0 & \cdots & 0 & 0 \end{bmatrix} \boldsymbol{x} = \begin{bmatrix} * \\ \vdots \\ * \\ 0 \end{bmatrix}
$$

Therefore

$$
LU \boldsymbol{x} = \begin{bmatrix} & & \\ \boldsymbol{\ell}_1 & \cdots & \boldsymbol{\ell}_n \\ & & \end{bmatrix} \begin{bmatrix} * \\ \vdots \\ * \\ 0 \end{bmatrix} = (*) \boldsymbol{\ell}_1 + \cdots + (*) \boldsymbol{\ell}_r
$$
```

## Rank-Nullity Theorem

```{div} theorem
Let $A$ be an $m \times n$ matrix. Then

$$
\mathrm{rank}(A) + \dim(N(A)) = n
$$

---

*Proof*. The dimension of $N(A)$ is equal to the number of columns of the row echelon form of $A$ *without* a leading nonzero entry, and $\mathrm{rank}(A)$ is equal to the number of columns of the row echelon form of $A$ *with* a leading nonzero entry, and there are $n$ total columns.
```

## Exercises

**Exercise 1.** Let $A = LU$ be the LU decomposition of $A$. Determine whether the statement is **True** or **False**.

  * $N(A) = N(U)$
  * $\dim (N(A)) = \dim (N(U))$
  * $R(A) = R(U)$
  * $\dim (R(A)) = \dim (R(U))$