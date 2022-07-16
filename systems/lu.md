# LU Decomposition

```{div} bigidea
Record the row operations of the Gaussian elimination algorithm in the LU decomposition and use the decomposition $A = LU$ in backward/forward subsitution to efficiently solve a system of linear equations $A \boldsymbol{x} = \boldsymbol{b}$.
```

```{image} /img/01_01_01.png
:width: 100%
:align: center
```

## Gaussian Elimination

```{div} definition
A **linear system of equations** is a collection of equations of the form

$$
\begin{array}{ccccccccc}
a_{1,1} x_1 & + & a_{1,2} x_2 & + & \cdots & + & a_{1,n} x_n & = & b_1 \\
a_{2,1} x_1 & + & a_{2,2} x_2 & + & \cdots & + & a_{2,n} x_n & = & b_2 \\
& \vdots & & & & \vdots & & \vdots & \\
a_{m,1} x_1 & + & a_{m,2} x_2 & + & \cdots & + & a_{m,n} x_n & = & b_m
\end{array}
$$

where the coefficients $a_{i,j}$ and $b_i$ are known constants and $x_i$ are unknown variables. In matrix notation, a linear system is $A \boldsymbol{x} = \boldsymbol{b}$ where

$$
A =
\begin{bmatrix}
a_{1,1} & a_{1,2} & \cdots & a_{1,n} \\
a_{2,1} & a_{2,2} & \cdots & a_{2,n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{m,1} & a_{m,2} & \cdots & a_{m,n}
\end{bmatrix}
\hspace{10mm}
\boldsymbol{x} = \begin{bmatrix} x_1 \\ x_2 \\ \vdots \\ x_n \end{bmatrix}
\hspace{10mm}
\boldsymbol{b} = \begin{bmatrix} b_1 \\ b_2 \\ \vdots \\ b_m \end{bmatrix}
$$
```

```{div} definition
There are three types of **elementary row operations** are:

1. Interchange two rows.
2. Multiply a row by a nonzero number.
3. Add a multiple of a row to another row.
```

```{div} definition
A matrix is in **row echelon form** if:

1. All zero rows are at the bottom.
2. The first nonzero entry in a row (from the left) is located to the right of the first nonzero entry in every row above.
```

```{div} example
The following matrices are in row echelon form:

$$
\begin{bmatrix} \times & * & * \\ 0 & \times & * \\ 0 & 0 & \times \end{bmatrix} \hspace{10mm}
\begin{bmatrix} \times & * & * & * & * \\ 0 & \times & * & * & * \\ 0 & 0 & \times & * & * \\ 0 & 0 & 0 & 0 & 0 \end{bmatrix} \hspace{10mm}
\begin{bmatrix} \times & * & * & * & * \\ 0 & 0 & \times & * & * \\ 0 & 0 & 0 & 0 & \times \\ 0 & 0 & 0 & 0 & 0 \end{bmatrix} \hspace{10mm}
\begin{bmatrix} \times & * \\ 0 & 0 \end{bmatrix}
$$

Symbol $\times$ denotes a nonzero entry and $*$ is any value.
```

```{div} definition
The **Gaussian elimination algorithm** applies a sequence of elementary row operations to transform a matrix to row echelon form. See [Wikipedia:Gaussian Elimination](https://en.wikipedia.org/wiki/Gaussian_elimination) to review how the algorithm works.
```

```{div} definition
Let $A$ be a matrix and let $U$ be a matrix in row echelon form obtained from $A$ through a sequence of elementary row operations. The **rank** of $A$ is the number of nonzero rows in $U$. Denote the rank of a matrix $A$ by $\mathrm{rank}(A)$.
```

```{div} definition
Let $A$ be a $m \times n$ matrix with $\mathrm{rank}(A) = r$.

The system $A \boldsymbol{x} = \boldsymbol{b}$ is **inconsistent** (ie. no solution) if $\mathrm{rank}(A) < \mathrm{rank}([ \, A \ \boldsymbol{b} \, ])$. In other words, the row echelon form of the augmented matrix $[ \, A \ \boldsymbol{b} \, ]$ has a row of the form

$$
\begin{array}{rrrr|r} 0 & 0 & \cdots & 0 & 1 \end{array}
$$

which implies $0 = 1$.

The system has a **unique solution** when $\mathrm{rank}(A) = \mathrm{rank}([ \, A \ \boldsymbol{b} \, ])$ and $\mathrm{rank}(A) = n$. In other words, the system is consistent and the rank of $A$ is equal to the number of variables in the system.

The system has **infinitely many solutions** when $\mathrm{rank}(A) = \mathrm{rank}([ \, A \ \boldsymbol{b} \, ])$ and $\mathrm{rank}(A) < n$.
```

```{div} example
Find all solutions of $A \boldsymbol{x} = \boldsymbol{b}$ where

$$
A = \begin{bmatrix} 1 & 1 & 1 \\ 1 & 0 & 1 \\ 2 & 5 & 2 \end{bmatrix} \hspace{10mm} \boldsymbol{b} = \begin{bmatrix} 2 \\ 1 \\ 7 \end{bmatrix}
$$

Add -1 times row 1 to row 2, and add -2 times row 1 to row 3:

$$
\left[ \begin{array}{rrr|r} 1 & 1 & 1 & 2 \\ 1 & 0 & 1 & 1 \\ 2 & 5 & 2 & 7 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrr|r} 1 & 1 & 1 & 2 \\ 0 & -1 & \phantom{+}0 & -1 \\ 0 & 3 & 0 & 3 \end{array} \right]
$$

Add 3 times row 2 to row 3:

$$
\left[ \begin{array}{rrr|r} 1 & 1 & 1 & 2 \\ 0 & -1 & \phantom{+}0 & -1 \\ 0 & 3 & 0 & 3 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrr|r} 1 & 1 & 1 & 2 \\ 0 & -1 & \phantom{+}0 & -1 \\ 0 & 0 & 0 & 0 \end{array} \right]
$$

Since there is no leading one in column 3, assign the free variable $t$ to $x_3$. The second row implies $x_2 = 1$ and the first row gives us

$$
x_1 = 2 - x_2 - x_3 = 1 - t
$$

There are infinitely many solutions

$$
\begin{bmatrix} 1-t \\ 1 \\ t \end{bmatrix} \ , \ \ t \in \mathbb{R}
$$
```

```{div} example
Find all solutions of $A \boldsymbol{x} = \boldsymbol{b}$ where

$$
A = \left[ \begin{array}{rrrr} 1 & -1 & \phantom{+}1 & -2 \\ -1 & 1 & 1 & 1 \\ -1 & 2 & 3 & 1 \\ 1 & -1 & 2 & 1 \end{array} \right] \hspace{10mm} \boldsymbol{b} = \left[ \begin{array}{r} 1 \\ -1 \\ 2 \\ 1 \end{array} \right]
$$

Add row 1 to row 2, add row 2 to row 3 and add -1 times row 1 to row 4:

$$
\left[ \begin{array}{rrrr|r} 1 & -1 & \phantom{+}1 & -2 & 1 \\ -1 & 1 & 1 & 1 & -1 \\ -1 & 2 & 3 & 1 & 2 \\ 1 & -1 & 2 & 1 & 1 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrrr|r} 1 & -1 & \phantom{+}1 & -2 & \phantom{+}1 \\ 0 & 0 & 2 & -1 & 0 \\ 0 & 1 & 4 & -1 & 3 \\ 0 & 0 & 1 & 3 & 0 \end{array} \right]
$$

Move row 2 to the bottom:

$$
\left[ \begin{array}{rrrr|r} 1 & -1 & \phantom{+}1 & -2 & \phantom{+}1 \\ 0 & 0 & 2 & -1 & 0 \\ 0 & 1 & 4 & -1 & 3 \\ 0 & 0 & 1 & 3 & 0 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrrr|r} 1 & -1 & \phantom{+}1 & -2 & \phantom{+}1 \\ 0 & 1 & 4 & -1 & 3 \\ 0 & 0 & 1 & 3 & 0 \\ 0 & 0 & 2 & -1 & 0 \end{array} \right]
$$

Add -2 times row 3 to row 4:

$$
\left[ \begin{array}{rrrr|r} 1 & -1 & \phantom{+}1 & -2 & \phantom{+}1 \\ 0 & 1 & 4 & -1 & 3 \\ 0 & 0 & 1 & 3 & 0 \\ 0 & 0 & 2 & -1 & 0 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrrr|r} 1 & -1 & \phantom{+}1 & -2 & \phantom{+}1 \\ 0 & 1 & 4 & -1 & 3 \\ 0 & 0 & 1 & 3 & 0 \\ 0 & 0 & 0 & -7 & 0 \end{array} \right]
$$

There is a unique solution

$$
\boldsymbol{x} = \begin{bmatrix} 4 \\ 3 \\ 0 \\ 0 \end{bmatrix}
$$
```

```{div} example
Find all solutions of $A \boldsymbol{x} = \boldsymbol{b}$ where

$$
A = \begin{bmatrix} 1 & 1 & 2 & 1 \\ 1 & 0 & 1 & 1 \\ 0 & 1 & 1 & 0 \end{bmatrix} \hspace{10mm} \boldsymbol{b} = \begin{bmatrix} 1 \\ 1 \\ 2 \end{bmatrix}
$$

Add -1 times row 1 to row 2:

$$
\left[ \begin{array}{rrrr|r} 1 & 1 & 2 & 1 & 1 \\ 1 & 0 & 1 & 1 & 1 \\ 0 & 1 & 1 & 0 & 2 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrrr|r} 1 & 1 & 2 & 1 & 1 \\ 0 & -1 & -1 & \phantom{+}0 & 0 \\ 0 & 1 & 1 & 0 & 2 \end{array} \right]
$$

Add row 2 to row 3:

$$
\left[ \begin{array}{rrrr|r} 1 & 1 & 2 & 1 & 1 \\ 0 & -1 & -1 & \phantom{+}0 & 0 \\ 0 & 1 & 1 & 0 & 2 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrrr|r} 1 & 1 & 2 & 1 & 1 \\ 0 & -1 & -1 & \phantom{+}0 & 0 \\ 0 & 0 & 0 & 0 & 2 \end{array} \right]
$$

The systems is inconsistent. In other words, there is no solution.
```

## LU Decomposition

```{div} definition
A **lower triangular matrix** is a matrix with zeros above the diagonal. For example:

$$
\begin{bmatrix} * &  &  &  \\ * & * &  &  \\ * & * & * &  \\ * & * & * & * \end{bmatrix}
\hspace{10mm}
\begin{bmatrix} * &  &  \\ * & * &  \\ * & * & * \\ * & * & * \\ * & * & * \end{bmatrix}
\hspace{10mm}
\begin{bmatrix} * &  &  &  &  \\ * & * &  &  &  \\ * & * & * &  &  \end{bmatrix}
$$

A **unit lower triangular matrix** is a *square* matrix with ones on the diagonal and zeros above the diagonal. For example:

$$
\begin{bmatrix} 1 & & \\ * & 1 & & \\ * & * & 1 & \\ * & * & * & 1 \end{bmatrix}
$$

An **upper triangular matrix** is a matrix with zeros below the diagonal. For example:

$$
\begin{bmatrix} * & * & * & * \\  & * & * & * \\  &  & * & * \\  &  &  & * \end{bmatrix}
\hspace{10mm}
\begin{bmatrix} * & * & * \\  & * & * \\  &  & * \\  &  &  \\  &  &  \end{bmatrix}
\hspace{10mm}
\begin{bmatrix} * & * & * & * & * \\  & * & * & * & * \\  &  & * & * & * \end{bmatrix}
$$

A **unit upper triangular matrix** is a *square* matrix with ones on the diagonal and zeros below the diagonal. For example:

$$
\begin{bmatrix} 1 & * & * & * \\ & 1 & * & * \\ & & 1 & * \\ & & & 1 \end{bmatrix}
$$
```

```{div} theorem
Let $E$ be the $m \times m$ matrix with ones along the diagonal, $c$ in the entry at row $i$ and column $j$ with $i > j$ and all other entries are zeros

$$
E = \begin{bmatrix} 1 & & & & \\ & \ddots & & & \\ & & \ddots & & \\ & c & & \ddots & \\ & & & & 1 \end{bmatrix}
$$

Then, for any $m \times n$ matrix $A$, matrix multiplication $EA$ applies to $A$ the elementary row operation: add $c$ times row $j$ to row $i$.

Furthermore, the inverse of $E$ is given by

$$
E^{-1} = \left[ \begin{array}{rrrrr} 1 & & & & \\ & \ddots & & & \\ & & \ddots & & \\ & -c & & \ddots & \\ & & & & 1 \end{array} \right]
$$

where $-c$ is the entry at row $i$ and column $j$.
```

```{div} example
Consider the matrix

$$
A = \left[ \begin{array}{rrrr} 1 & -1 & \phantom{+}1 & -2 \\ -1 & 1 & 1 & 1 \\ -1 & 2 & 3 & 1 \\ 1 & -1 & 2 & 1 \end{array} \right]
$$

The elementary matrix which adds -1 times row 1 to row 4 is

$$
E = \left[ \begin{array}{rrrr|r} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ -1 & \phantom{+}0 & \phantom{+}0 & \phantom{+}1 \end{array} \right]
$$

Perform matrix multiplication to verify

$$
EA =
\left[ \begin{array}{rrrr|r} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ -1 & \phantom{+}0 & \phantom{+}0 & \phantom{+}1 \end{array} \right]
\left[ \begin{array}{rrrr} 1 & -1 & \phantom{+}1 & -2 \\ -1 & 1 & 1 & 1 \\ -1 & 2 & 3 & 1 \\ 1 & -1 & 2 & 1 \end{array} \right]
=
\left[ \begin{array}{rrrr} 1 & -1 & \phantom{+}1 & -2 \\ -1 & 1 & 1 & 1 \\ -1 & 2 & 3 & 1 \\ 0 & 0 & 1 & 3 \end{array} \right]
$$
```

````{div} theorem
If $A$ can be reduced by Gaussian elimination to row echelon form *only* with operations "add $c$ times row $j$ to row $i$" (in other words, without scaling rows and without interchanging rows), then $A$ has an **LU decomposition** of the form

$$A = LU$$

where $L$ is a unit lower triangular matrix and $U$ is an upper triangular matrix. In particular, after performing Gaussian elimination on $A$, the matrix $U$ is the corresponding row echelon form of $A$ and $L$ is given by

$$
L =
\begin{bmatrix}
1 & & & & \\
-c_{2,1} & 1 & & & \\
-c_{3,1} & -c_{3,2} & 1 & & \\
\vdots & \vdots & \ddots & \ddots & \\
-c_{m,1} & -c_{m,2} & \cdots & -c_{m,m-1} & 1
\end{bmatrix}
$$

where each entry corresponds to the elementary row operation "add $c_{i,j}$ times row $j$ to row $i$" performed during Gaussian elimination.

```{dropdown} Proof
Gaussian elimination (without scaling rows and without interchanging rows) yields a sequence of elementary row operations of the form "add $c_{i,j}$ times row $j$ to row $i$" for each pair $(i,j)$ with $i > j$. Let $E_{i,j}$ be the corresponding elementary matrix for each $(i,j)$. Then

$$
E_{m,m-1} E_{m,m-2} E_{m-1,m-2} \cdots E_{m,2} \cdots E_{3,2} E_{m,1} \cdots E_{2,1} A = U
$$

where $U$ is in row echelon form. Rearrange the equation to get

$$
A = E_{2,1}^{-1} \cdots E_{m,1}^{-1} E_{3,2}^{-1} \cdots E_{m,2}^{-1} \cdots E_{m-1,m-2}^{-1} E_{m,m-2}^{-1} E_{m,m-1}^{-1} U
$$

Each matrix $E_{i,j}^{-1}$ is the unit lower triangular matrix with entry $-c_{i,j}$ at index $(i,j)$. In particular, we have

$$
E_{m,m-1}^{-1} =
\begin{bmatrix}
1 & & & & \\
  & 1 & & & \\
  &   & \ddots & & \\
  &   &   & 1 & \\
 & & & -c_{m,m-1} & 1
\end{bmatrix}
$$

Then $E_{m,m-2}^{-1}E_{m,m-1}^{-1}$ results in add $-c_{m,m-2}$ times row $m-2$ to row $m$ and so

$$
E_{m,m-2}^{-1} E_{m,m-1}^{-1} =
\begin{bmatrix}
1 & & & & \\
  & \ddots & & & \\
  &   & 1 & & \\
  &   &   & 1 & \\
 & & -c_{m,m-2} & -c_{m,m-1} & 1
\end{bmatrix}
$$

Applying all the corresponding rows operations in order yields the result

$$
E_{2,1}^{-1} \cdots E_{m,1}^{-1} E_{3,2}^{-1} \cdots E_{m,2}^{-1} \cdots E_{m-1,m-2}^{-1} E_{m,m-2}^{-1} E_{m,m-1}^{-1}
$$

$$
=
\begin{bmatrix}
1 & & & & \\
-c_{2,1} & 1 & & & \\
-c_{3,1} & -c_{3,2} & 1 & & \\
\vdots & \vdots & \ddots & \ddots & \\
-c_{m,1} & -c_{m,2} & \cdots & -c_{m,m-1} & 1
\end{bmatrix}
$$
```
````

```{div} example
Compute the LU decomposition of

$$
A = \begin{bmatrix} 2 & 1 & 1 \\ 2 & 0 & 2 \\ 4 & 3 & 4 \end{bmatrix}
$$

Add -1 times row 1 to row 2 and add -2 times row 1 to row 3

$$
\begin{bmatrix} 2 & 1 & 1 \\ 2 & 0 & 2 \\ 4 & 3 & 4 \end{bmatrix}
\longrightarrow
\left[ \begin{array}{rrr} 2 & 1 & 1 \\ 0 & -1 & \phantom{+}1 \\ 0 & 1 & 2 \end{array} \right]
$$

Add row 2 to row 3

$$
\left[ \begin{array}{rrr} 2 & 1 & 1 \\ 0 & -1 & \phantom{+}1 \\ 0 & 1 & 2 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrr} 2 & 1 & 1 \\ 0 & -1 & \phantom{+}1 \\ 0 & 0 & 3 \end{array} \right]
$$

In terms of elementary matrices, we have just shown that

$$
\begin{align*}
E_{3,2} E_{3,1} E_{2,1} A &= U \\
& \\
\left[ \begin{array}{rrr} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & 1 & 1 \end{array} \right]
\left[ \begin{array}{rrr} 1 & 0 & 0 \\ 0 & 1 & 0 \\ -2 & 0 & 1 \end{array} \right]
\left[ \begin{array}{rrr} 1 & 0 & 0 \\ -1 & 1 & 0 \\ 0 & 0 & 1 \end{array} \right]
\left[ \begin{array}{rrr} 2 & 1 & 1 \\ 2 & 0 & 2 \\ 4 & 3 & 4 \end{array} \right]
&=
\left[ \begin{array}{rrr} 2 & 1 & 1 \\ 0 & -1 & \phantom{+}1 \\ 0 & 0 & 3 \end{array} \right]
\end{align*}
$$

Rearrange the equation to get

$$
\begin{align*}
A &= E_{2,1}^{-1} E_{3,1}^{-1} E_{3,2}^{-1} U \\
& \\
\begin{bmatrix} 2 & 1 & 1 \\ 2 & 0 & 2 \\ 4 & 3 & 4 \end{bmatrix}
&=
\begin{bmatrix} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 0 & 0 & 1 \end{bmatrix}
\begin{bmatrix} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 2 & 0 & 1 \end{bmatrix}
\left[ \begin{array}{rrr} 1 & 0 & 0 \\ 0 & 1 & 0 \\ 0 & -1 & \phantom{+}1 \end{array} \right]
\left[ \begin{array}{rrr} 2 & 1 & 1 \\ 0 & -1 & \phantom{+}1 \\ 0 & 0 & 3 \end{array} \right]
\end{align*}
$$

Combine the matrices as in the proof of the LU decomposition to find

$$
\begin{bmatrix} 2 & 1 & 1 \\ 2 & 0 & 2 \\ 4 & 3 & 4 \end{bmatrix}
=
\left[ \begin{array}{rrr} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 2 & -1 & \phantom{+}1 \end{array} \right]
\left[ \begin{array}{rrr} 2 & 1 & 1 \\ 0 & -1 & \phantom{+}1 \\ 0 & 0 & 3 \end{array} \right]
$$

Therefore $A = LU$ where

$$
L =
\left[ \begin{array}{rrr} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 2 & -1 & \phantom{+}1 \end{array} \right]
\hspace{10mm}
U = \left[ \begin{array}{rrr} 2 & 1 & 1 \\ 0 & -1 & \phantom{+}1 \\ 0 & 0 & 3 \end{array} \right] \\
$$

Note that we can construct the matrix $L$ directly from the list of operations:

1. Add -1 times row 1 to row 2.
2. Add -2 times row 1 to row 3.
3. Add 1 times row 2 to row 3.
```

```{div} example
Compute the LU decomposition of

$$
A = \left[ \begin{array}{rrrr} 1 & 0 & 2 & 1 \\ 2 & 1 & 5 & 3 \\ 1 & 0 & 0 & 2 \\ 0 & -1 & \phantom{+}1 & \phantom{+}1 \end{array} \right]
$$

Add -2 times row 1 to row 2 and add -1 times row 1 to row 3:

$$
\left[ \begin{array}{rrrr} 1 & 0 & 2 & 1 \\ 2 & 1 & 5 & 3 \\ 1 & 0 & 0 & 2 \\ 0 & -1 & \phantom{+}1 & \phantom{+}1 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrrr} 1 & 0 & 2 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & -2 & 1 \\ 0 & -1 & 1 & \phantom{+}1 \end{array} \right]
$$

Add row 2 to row 4:

$$
\left[ \begin{array}{rrrr} 1 & 0 & 2 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & -2 & 1 \\ 0 & -1 & 1 & \phantom{+}1 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrrr} 1 & 0 & 2 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & -2 & 1 \\ 0 & \phantom{+}0 & 2 & \phantom{+}2 \end{array} \right]
$$

Add row 3 to row 4:

$$
\left[ \begin{array}{rrrr} 1 & 0 & 2 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & -2 & 1 \\ 0 & \phantom{+}0 & 2 & \phantom{+}2 \end{array} \right]
\longrightarrow
\left[ \begin{array}{rrrr} 1 & 0 & 2 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & -2 & 1 \\ 0 & \phantom{+}0 & 0 & \phantom{+}3 \end{array} \right]
$$

Therefore $A = LU$ where

$$
L = \left[ \begin{array}{rrrr} 1 & 0 & 0 & 0 \\ 2 & 1 & 0 & 0 \\ 1 & 0 & 1 & 0 \\ 0 & -1 & -1 & \phantom{+}1 \end{array} \right]
\hspace{10mm}
U = \left[ \begin{array}{rrrr} 1 & 0 & 2 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & -2 & 1 \\ 0 & \phantom{+}0 & 0 & \phantom{+}3 \end{array} \right]
$$
```

```{div} theorem
Suppose $A$ has an LU decomposition $A = LU$.

1. $\mathrm{rank}(A) = \mathrm{rank}(U)$
2. $\det(A) = \det(U) = u_{1,1} \cdots u_{m,m}$ where $u_{1,1} , \dots , u_{m,m}$ are the diagonal entries of $U$.
```

```{div} note
Not all matrices have an LU decomposition. For example,

$$
A = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}
$$

does not have an LU decomposition (*why not?*). However, if we allow partial pivoting (ie. interchanging rows during Gaussian elimination), then Gaussian elimination with partial pivoting computes for *any* matrix $A$ a decomposition $A = PLU$ where $P$ is a permutation matrix, $L$ is unit lower triangular and $U$ is upper triangular. This is called the **LU decomposition with partial pivoting** and has similar computational advantages as the LU decomposition. See [Wikipedia:LU decomposition](https://en.wikipedia.org/wiki/LU_decomposition).
```

## Forward and Backward Substitution

```{div} definition
Let $A = LU$ be the LU decomposition of $A$, let $\ell_{i,j}$ denote the entries of $L$ and let $u_{i,j}$ denote the entries of $U$. Consider the system $A \boldsymbol{x} = \boldsymbol{b}$ and let $\boldsymbol{y} = U \boldsymbol{x}$.

**Forward substitution** is the process of solving the lower triangular system $L \boldsymbol{y} = \boldsymbol{b}$ from top to bottom:

$$
\begin{align*}
y_1 &= b_1 \\
y_2 &= b_2 - \ell_{2,1} y_1 \\
 & \vdots \\
y_n &= b_n - \ell_{n,1} y_1 - \cdots - \ell_{n,n-1} y_{n-1}
\end{align*}
$$

**Backward substitution** is the process of solving the upper triangular system $U \boldsymbol{x} = \boldsymbol{y}$ from bottom to top:

$$
\begin{align*}
x_n &= y_n / u_{n,n} \\
x_{n-1} &= (y_{n-1} - u_{n-1,n} x_n) / u_{n-1,n-1} \\
 & \vdots \\
x_1 &= (y_1 - u_{1,2} x_2 - \cdots - u_{1,n} x_n) / u_{1,1}
\end{align*}
$$
```

```{div} example
Solve the system $A \boldsymbol{x} = \boldsymbol{b}$ where

$$
A = LU =
\begin{bmatrix} 1 & 0 & 0 \\ 2 & 1 & 0 \\ 1 & 1 & 1 \end{bmatrix}
\begin{bmatrix} 2 & 4 & 2 \\ 0 & 1 & 1 \\ 0 & 0 & 1 \end{bmatrix}
\hspace{10mm}
\boldsymbol{b} = \left[ \begin{array}{r} -1 \\ 1 \\ 2 \end{array} \right]
$$

Solve $L \boldsymbol{y} = \boldsymbol{b}$

$$
\begin{align*}
y_1 &= -1 \\
y_2 &= 1 - 2(-1) = 3 \\
y_3 &= 2 - (-1) - 3 = 0
\end{align*}
$$

and then solve $U \boldsymbol{x} = \boldsymbol{y}$

$$
\begin{align*}
x_3 &= 0 \\
x_2 &= 3 \\
x_1 &= (-1 - 4(3) - 0) / 2 = -13/2
\end{align*}
$$

Therefore

$$
\boldsymbol{x} = \begin{bmatrix} -13/2 \\ 3 \\ 0 \end{bmatrix}
$$
```

```{div} note
The LU decomposition is especially useful when solving many different systems with the *same* coefficient matrix $A$. For example, to compute the inverse $A^{-1}$ of a square matrix of size $n$ we need to solve $n$ different systems $A \boldsymbol{x}_k = \boldsymbol{e}_k$ for $k=1,\dots,n$ where $\boldsymbol{e}_k$ is the $k$th column of the identity matrix $I$. The result is $A^{-1} = [\boldsymbol{x}_1 \cdots \boldsymbol{x}_n]$. In other words, the columns of $A^{-1}$ are given by $\boldsymbol{x}_1, \dots, \boldsymbol{x}_n$.
```

## Exercises

````{div} exercise
Let $A$ be a $m$ by $n$ matrix. Determine whether the statement is **True** or **False**.

  * If $m > n$ and $\mathrm{rank}(A) = n$, then here is a unique solution of $A \boldsymbol{x} = \boldsymbol{b}$ for any $\boldsymbol{b}$.
  * If $m < n$ and $\mathrm{rank}(A) = m$, then there are infinitely many solutions of $A \boldsymbol{x} = \boldsymbol{b}$ for any $\boldsymbol{b}$.
  * If $m > n$ and $\mathrm{rank}(A) = n$, then if the system $A \boldsymbol{x} = \boldsymbol{b}$ has one solution then there is only one solution.
  * If $m > n$ and $\mathrm{rank}(A) < n$, then if the system $A \boldsymbol{x} = \boldsymbol{b}$ has one solution then there are infinitely many solutions.
  * If $A = LU$ is the LU decomposition of $A$ then $\det(L) \not= 0$.

```{dropdown} Solution
* False
* True
* True
* True
* True
```

````

````{div} exercise
Determine whether the statement is **True** or **False**. If $A$ is of the form

$$
A = \begin{bmatrix} * & * & 0 & 0 \\ * & * & * & 0 \\ 0 & * & * & * \\ 0 & 0 & * & * \end{bmatrix}
$$
and the $LU$ decomposition $A = LU$ exists, then $L$ and $U$ are of the form

$$
L = \begin{bmatrix} 1 & 0 & 0 & 0 \\ * & 1 & 0 & 0 \\ 0 & * & 1 & 0 \\ 0 & 0 & * & 1 \end{bmatrix}
\hspace{5mm}
U = \begin{bmatrix} * & * & 0 & 0 \\ 0 & * & * & 0 \\ 0 & 0 & * & * \\ 0 & 0 & 0 & * \end{bmatrix}
$$

```{dropdown} Solution
True
```

````

````{div} exercise
Let $I$ be the identity matrix of size $n$ and let $R$ be the $n$ by $n$ matrix with all zeros except for the nonzero scalar $c$ at row $i$ and column $j$ where $i \not= j$. Let $E = I + R$ and let $A$ be any $n$ by $n$ matrix.
  * Matrix multiplication $EA$ is equivalent to which elementary row/column operation on $A$?
  * Matrix multiplication $AE$ is equivalent to which elementary row/column operation on $A$?

```{dropdown} Solution
* Add $c$ times row $j$ to row $i$.
* Add $c$ times column $i$ to column $j$.
```

````

````{div} exercise
Find a value $c$ such that the system $A \boldsymbol{x} = \boldsymbol{b}$ has infinitely many solutions where

$$
A = \left[ \begin{array}{rrr} 3 & -1 & 2 \\ 1 & 1 & -1 \\ 2 & -2 & 3 \end{array} \right]
\hspace{5mm}
\boldsymbol{b} = \begin{bmatrix} 3 \\ 2 \\ c \end{bmatrix}
$$


```{dropdown} Solution
$c = 1$
```

````

````{div} exercise
Compute the LU decomposition of the matrix

$$
A = \left[ \begin{array}{rrrr} 2 & 0 & 1 & 1 \\ -2 & -1 & 2 & -1 \\ 0 & 2 & -5 & -1 \\ 4 & 0 & 6 & 0 \end{array} \right]
$$

```{dropdown} Solution
$$
L = \left[ \begin{array}{rrrr} 1 & 0 & \phantom{+}0 & \phantom{+}0 \\ -1 & 1 & 0 & 0 \\ 0 & -2 & 1 & 0 \\ 2 & 0 & 4 & 1 \end{array} \right]
\hspace{20mm}
U = \left[ \begin{array}{rrrr} 2 & 0 & \phantom{+}1 & 1 \\ 0 & -1 & 3 & 0 \\ 0 & 0 & 1 & -1 \\ 0 & 0 & 0 & 2 \end{array} \right]
$$
```

````

````{div} exercise
Consider the matrix

$$
A = \left[ \begin{array}{rrrr}
-3 & \phantom{+}1 & 2 & \phantom{+}0 \\
3 & 1 & -2 & 1 \\
-6 & 2 & 5 & 1 \\
-9 & 3 & 4 & 2
\end{array} \right]
$$

Find the $LU$ decomposition of $A$ and compute $\mathrm{det}(A)$.

```{dropdown} Solution
$$
L = \left[ \begin{array}{rrrr} 1 & \phantom{+}0 & 0 & \phantom{+}0 \\ -1 & 1 & 0 & 0 \\ 2 & 0 & 1 & 0 \\ 3 & 0 & -2 & 1 \end{array} \right]
\hspace{10mm}
U = \left[ \begin{array}{rrrr} -3 & 1 & 2 & 0 \\ 0 & 2 & 0 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 4 \end{array} \right]
$$

$$
\mathrm{det}(A) = -24
$$

```

````

````{div} exercise
Find the solution of the system $A \boldsymbol{x} = \boldsymbol{b}$ for

$$
\boldsymbol{b} = \left[ \begin{array}{r} 2 \\ -1 \\ 2 \end{array} \right]
$$

given the LU decomposition

$$
A = LU =
\left[ \begin{array}{rrr} 1 & \phantom{+}0 & \phantom{+}0 \\ -1 & 1 & 0 \\ 1 & 1 & 1 \end{array} \right]
\left[ \begin{array}{rrr} 3 & 0 & \phantom{+}1 \\ 0 & -1 & 1 \\ 0 & 0 & 1 \end{array} \right]
$$

```{dropdown} Solution
$$
\boldsymbol{x} = \left[ \begin{array}{r} 1 \\ -2 \\ -1 \end{array} \right]
$$
```

````

````{div} exercise
Suppose we compute a decomposition $A = L_0U_0$ such that $U_0$ is *unit* upper triangular and $L_0$ is lower triangular. Describe a method to derive a decomposition $A = LU$ such that $L$ is *unit* lower triangular and $U$ is upper triangular.

```{dropdown} Solution
Factor of the diagonal entries of $L_0$ into a diagonal matrix $D$ such that $L_0 = LD$ where $L$ is unit lower triangular then let $U = DU_0$ to find $A = LU$. For example:

$$
A = \begin{bmatrix} 2 & 0 \\ 4 & 3 \end{bmatrix} \left[ \begin{array}{rr} 1 & -1 \\ 0 & 1 \end{array} \right]
= \begin{bmatrix} 1 & 0 \\ 2 & 1 \end{bmatrix} \begin{bmatrix} 2 & 0 \\ 0 & 3 \end{bmatrix} \left[ \begin{array}{rr} 1 & -1 \\ 0 & 1 \end{array} \right]
= \begin{bmatrix} 1 & 0 \\ 2 & 1 \end{bmatrix} \left[ \begin{array}{rr} 2 & -2 \\ 0 & 3 \end{array} \right]
$$
```

````