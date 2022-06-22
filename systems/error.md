# Error Analysis

```{div} bigidea
The condition number of a nonsingular matrix $A$ is $\mathrm{cond}(A) = \| A \| \| A^{-1} \|$. Given a linear system $A \boldsymbol{x} = \boldsymbol{b}$, the condition number of $A$ quantifies how sensitive the solution $\boldsymbol{x}$ is relative to changes in $\boldsymbol{b}$.
```

## Matrix Norm

```{div} definition
The **norm** of a vector $\boldsymbol{x} \in \mathbb{R}^n$ is given by

$$
\| \boldsymbol{x} \| = \sqrt{ | x_1|^2 + \cdots + | x_n |^2 } = \sqrt{ \sum_{k=1}^n | x_k |^2 }
$$
```

```{div} definition
The vector norm $\| \cdot \|$ satisfies the following properties.

1. $\| \boldsymbol{x} \| \geq 0$ for all $\boldsymbol{x} \in \mathbb{R}^n$
2. $\| \boldsymbol{x} \| = 0$ if and only if $\boldsymbol{x} = \boldsymbol{0}$
3. $\| c \boldsymbol{x} \| = |c| \| \boldsymbol{x} \|$ for any $c \in \mathbb{R}$ and $\boldsymbol{x} \in \mathbb{R}^n$
4. $\| \boldsymbol{x} + \boldsymbol{y} \| \leq \| \boldsymbol{x} \| + \| \boldsymbol{y} \|$ for all $\boldsymbol{x} , \boldsymbol{y} \in \mathbb{R}^n$

Condition 4 is called the **triangle inequality**. See [Wikipedia:Norm](https://en.wikipedia.org/wiki/Norm_(mathematics)).
```

```{div} definition
The **norm** of a matrix $A$ is

$$
\| A \| = \max_{\boldsymbol{x} \not= \boldsymbol{0} } \frac{ \| A \boldsymbol{x} \| }{ \| \boldsymbol{x}  \|}
$$

The matrix norm is also called the **operator norm**.
```

```{div} theorem
Let $A$ be a nonsingular matrix. Then

$$
\| A \| = \max_{ \| \boldsymbol{x} \| = 1 } \| A \boldsymbol{x} \|
\hspace{10mm} \text{and} \hspace{10mm}
\| A^{-1} \| = \frac{1}{ \displaystyle \min_{ \| \boldsymbol{x} \| = 1 } \| A \boldsymbol{x} \|}
$$

In other words, $\| A \|$ is the *maximum stretch of a unit vector* by the linear transformation $A$, and $\| A^{-1} \|$ is the *reciprocal* of the *minimum stretch of a unit vector* by the linear transformation $A$.

---

*Proof*. Note that $\| A \boldsymbol{x} \| / \| \boldsymbol{x} \|= \| A ( \boldsymbol{x} / \| \boldsymbol{x} \| ) \|$ therefore

$$
\| A \| = \max_{ \| \boldsymbol{x} \| = 1 } \| A \boldsymbol{x} \|
$$

Similarly, we can rearrange the definition of $\| A^{-1} \|$ to find:

$$
\begin{align*}
\| A^{-1} \| &= \max_{ \boldsymbol{x} \not= 0} \frac{\| A^{-1} \boldsymbol{x} \|}{\| \boldsymbol{x} \|} \\
&= \max_{ \boldsymbol{x} \not= 0} \frac{\| A^{-1} A \boldsymbol{x} \|}{\| A \boldsymbol{x} \|} \\
&= \max_{ \boldsymbol{x} \not= 0} \frac{\| \boldsymbol{x} \|}{\| A \boldsymbol{x} \|} \\
&= \max_{ \| \boldsymbol{x} \| = 1} \frac{1}{\| A \boldsymbol{x} \|} \\
&= \frac{1}{\displaystyle \min_{ \| \boldsymbol{x} \| = 1} \| A \boldsymbol{x} \|}
\end{align*}
$$
```

```{div} proposition
Let $D$ be a diagonal matrix and let $\boldsymbol{d}$ be the vector of diagonal entries of $D$:

$$
D = \begin{bmatrix} d_1 & & & \\ & d_2 & & \\ & & \ddots & \\ & & & d_n \end{bmatrix}
\hspace{10mm}
\boldsymbol{d} = \begin{bmatrix} d_1 \\ d_2 \\ \vdots \\ d_n \end{bmatrix}
$$

Then $\| D \| = \max \{ |d_1| , \dots, |d_n| \}$.

---

*Proof*. Compute

$$
\begin{align*}
\| D \| &= \max_{\| \boldsymbol{x} \| = 1} \| D \boldsymbol{x} \| \\
&= \max_{\| \boldsymbol{x} \| = 1} \sqrt{d_1^2 x_1^2 + \cdots + d_n^2 x_n^2} \\
&\leq \max_{\| \boldsymbol{x} \| = 1} \left( \max_i \{ | d_i | \} \right) \sqrt{ x_1^2 + \cdots + x_n^2} \\
&= \max_i \{ | d_i | \}
\end{align*}
$$
```

```{div} note
How do we compute the matrix norm $\| A \|$ for a general matrix? This is a nontrivial problem and we will see later how to use the singular values of $A$ to determine the matrix norm.
```

```{div} theorem
The matrix norm satisfies the properties:

1. $\| A \| > 0$ for all $A \not= 0$
2. $\| A \| = 0$ if and only $A = 0$
3. $\| c A \| = |c| \| A \|$ for any $c \in \mathbb{R}$
4. $\| A + B \| \leq \| A \| + \| B \|$
5. $\| A B \| \leq \| A \| \| B \|$ (if $A$ and $B$ are square matrices)
6. $\| A \boldsymbol{x} \| \leq \| A \| \| \boldsymbol{x} \|$ for any $\boldsymbol{x} \in \mathbb{R}^n$

See [Wikipedia: Matrix norm](https://en.wikipedia.org/wiki/Matrix_norm).
```

## Condition Number

```{div} definition
The **condition number** of a nonsingular square matrix $A$ is

$$
\mathrm{cond}(A) = \| A \| \| A^{-1} \|
$$

By convention, we define $\mathrm{cond}(A) = \infty$ if $\det(A) = 0$.
```

```{div} note
If $A$ is nonsingular, we have

$$
\mathrm{cond}(A) = \| A \| \| A^{-1} \| = \frac{\text{maximum stretch of a unit vector}}{\text{minimum stretch of a unit vector}}
$$
```

```{div} example
The image below shows the unit circle and its image under the linear transformation defined by a $2 \times 2$ matrix $A$. Determine $\| A \|$, $\| A^{-1} \|$ and $\mathrm{cond}(A)$.

![img/01_02_01.png](/img/01_02_01.png)

Observe the maximum stretch of a unit vector is $\| A \| =  3 / \sqrt{2}$, the minimum stretch is $1/\sqrt{2}$ therefore $\| A^{-1} \| = \sqrt{2}$ and the condition number is $\mathrm{cond}(A) = 3$.
```

## Relative Errors

```{div} theorem
Let $A$ be a nonsingular matrix and consider the linear system $A \boldsymbol{x} = \boldsymbol{b}$. If a small change $\Delta \boldsymbol{b}$ corresponds to a change $\Delta \boldsymbol{x}$ in the sense that $A(\boldsymbol{x} + \Delta \boldsymbol{x}) = \boldsymbol{b} + \Delta \boldsymbol{b}$, then

$$
\frac{\| \Delta \boldsymbol{x} \|}{\| \boldsymbol{x} \|} \leq \mathrm{cond}(A) \frac{\| \Delta \boldsymbol{b} \|}{\| \boldsymbol{b} \|}
$$

---

*Proof*. Since $A \boldsymbol{x} = \boldsymbol{b}$, we have $\Delta x = A^{-1} \Delta \boldsymbol{b}$. Computing norms we find

$$
\begin{align*}
\| \boldsymbol{b} \| &= \| A \boldsymbol{x} \| \\ & \\
\| \Delta \boldsymbol{x} \| \| \boldsymbol{b} \| &= \| A^{-1} \Delta \boldsymbol{b} \| \| A \boldsymbol{x} \| \\ & \\
\| \Delta \boldsymbol{x} \| \| \boldsymbol{b} \| &\leq \| A^{-1} \| \| \Delta\boldsymbol{b} \| \| A \| \| \boldsymbol{x} \| \\ & \\
\frac{\| \Delta \boldsymbol{x} \|}{ \| \boldsymbol{x} \|}  &\leq  \| A \| \| A^{-1} \| \frac{\| \Delta \boldsymbol{b} \|}{\| \boldsymbol{b} \|}
\end{align*}
$$
```

```{div} definition
Given a vector $\boldsymbol{b}$ and small change $\Delta \boldsymbol{b}$, the **relative change** (or **relative error**) is

$$
\frac{\| \Delta \boldsymbol{b} \|}{\| \boldsymbol{b} \|}
$$
```

```{div} note
The error bound

$$
\frac{\| \Delta \boldsymbol{x} \|}{\| \boldsymbol{x} \|} \leq \mathrm{cond}(A) \frac{\| \Delta \boldsymbol{b} \|}{\| \boldsymbol{b} \|}
$$

implies that if $A$ has a large condition number then small changes in $\boldsymbol{b}$ may result in *very* large changes in the solution $\boldsymbol{x}$. In other words, the solution $\boldsymbol{x}$ is sensitive to errors in $\Delta \boldsymbol{b}$.
```

## Exercises

**Exercise 1.** Determine whether the statement is **True** or **False**: If $\| A \| = 1$ then $A = I$.

**Exercise 2.** Suppose $A$ is a 2 by 2 matrix such that the image of the unit circle under the linear transformation $A$ is:

![/img/01_02_02.png](/img/01_02_02.png)

Determine $\mathrm{cond}(A)$.
