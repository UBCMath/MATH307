# Vector and Matrix Norms

```{div} bigidea
Norms provide different ways of describing the magnitude of a vector or matrix. The length of a vector $\| \boldsymbol{x} \|_2 = \sqrt{x_1^2 + \cdots + x_n^2}$ is called the 2-norm but there are infinitely many other norms like the 1-norm and $\infty$-norm. The norm $\| A \|$ of a matrix $A$ quantifies how much $A$ stretches vectors.
```

```{div} definition
A **norm** on the vector space $\mathbb{R}^n$ is a function $\| \cdot \|$ such that:

1. $\| \boldsymbol{x} \| \geq 0$ for all $\boldsymbol{x} \in \mathbb{R}^n$
2. $\| \boldsymbol{x} \| = 0$ if and only if $\boldsymbol{x} = \boldsymbol{0}$
3. $\| c \boldsymbol{x} \| = |c| \| \boldsymbol{x} \|$ for any $c \in \mathbb{R}$ and $\boldsymbol{x} \in \mathbb{R}^n$
4. $\| \boldsymbol{x} + \boldsymbol{y} \| \leq \| \boldsymbol{x} \| + \| \boldsymbol{y} \|$ for all $\boldsymbol{x} , \boldsymbol{y} \in \mathbb{R}^n$

Condition 4 is called the **triangle inequality**. See [Wikipedia:Norm](https://en.wikipedia.org/wiki/Norm_(mathematics)).
```

```{div} example
Let $\boldsymbol{x} \in \mathbb{R}^n$.

1. The 2-norm (also called the $\ell_2$-norm) is given by the familiar formula

   $$
   \| \boldsymbol{x} \|_2 = \sqrt{ | x_1|^2 + \cdots + | x_n |^2 } = \sqrt{ \sum_{k=1}^n | x_k |^2 }
   $$

2. More generally, the $p$-norm (also called the $\ell_p$-norm) is given by

   $$
   \| \boldsymbol{x} \|_p = \left( \sum_{k=1}^n | x_k |^p \right)^{1/p}
   $$

   For example, a commonly used norm is the 1-norm

   $$
   \| \boldsymbol{x} \|_1 =  | x_1| + \cdots + | x_n | = \sum_{k=1}^n | x_k |
   $$

3. The $\infty$-norm (also called the $\ell_\infty$-norm) is given by

   $$
   \| \boldsymbol{x} \|_{\infty} = \max_k | x_k |
   $$
```

```{div} definition
The **Frobenius norm** of a matrix $A$ (also called the **Hilbert-Schmidt norm**) is given by

$$
\|A\|_{\rm FR} = \sqrt{ \sum_{i=1}^m \sum_{j=1}^n |a_{i,j}|^2 }
$$
```

```{div} example
Let $D$ be a diagonal matrix and let $\boldsymbol{d}$ be the vector of diagonal entries of $D$:

$$
D = \begin{bmatrix} d_1 & & & \\ & d_2 & & \\ & & \ddots & \\ & & & d_n \end{bmatrix}
\hspace{10mm}
\boldsymbol{d} = \begin{bmatrix} d_1 \\ d_2 \\ \vdots \\ d_n \end{bmatrix}
$$

Then $\|D\|_{\rm FR} = \| \boldsymbol{d} \|_2$.
```

```{div} definition
Choose a vector norm $\| \cdot \|$. The **matrix norm** (with respect to the vector norm $\| \cdot \|$) is

$$
\| A \| = \max_{\boldsymbol{x} \not= \boldsymbol{0} } \frac{\| A \boldsymbol{x} \|}{ \| \boldsymbol{x}  \|}
$$

Note that $\| A \boldsymbol{x} \| / \| \boldsymbol{x} \|= \| A ( \boldsymbol{x} / \| \boldsymbol{x} \| ) \|$ therefore

$$
\| A \| = \max_{ \| \boldsymbol{x} \| = 1 } \| A \boldsymbol{x} \|
$$

In other words, the matrix norm is the maximum stretch of a unit vector by the linear transformation $A$. The matrix norm is also called the **operator norm**.
```

```{div} note
We *almost always* use the 2-norm when defining the matrix norm. Therefore

$$
\| A \| = \| A \|_2 = \max_{\boldsymbol{x} \not= \boldsymbol{0} } \frac{\| A \boldsymbol{x} \|_2}{ \| \boldsymbol{x}  \|_2}
$$

is the matrix norm with respect to the 2-norm unless a different vector norm is explicitly stated.
```

```{div} example
Let $D$ be a diagonal matrix and let $\boldsymbol{d}$ be the vector of diagonal entries of $D$:

$$
D = \begin{bmatrix} d_1 & & & \\ & d_2 & & \\ & & \ddots & \\ & & & d_n \end{bmatrix}
\hspace{10mm}
\boldsymbol{d} = \begin{bmatrix} d_1 \\ d_2 \\ \vdots \\ d_n \end{bmatrix}
$$

Then $\|D\|_2 = \| \boldsymbol{d} \|_{\infty}$.
```

```{div} theorem
A matrix norm satisfies the properties:

1. $\| A \| > 0$ for all $A \not= 0$
2. $\| A \| = 0$ if and only $A = 0$
3. $\| c A \| = |c| \| A \|$ for any $c \in \mathbb{R}$
4. $\| A + B \| \leq \| A \| + \| B \|$
5. $\| A B \| \leq \| A \| \| B \|$
6. $\| A \boldsymbol{x} \| \leq \| A \| \| \boldsymbol{x} \|$ for any $\boldsymbol{x} \in \mathbb{R}^n$
```

## Exercises

1. Prove that the $\infty$-norm satisfies the required properties of a norm.
2. Sketch the "unit ball" in $\mathbb{R}^2$ for each norm:

   $$
   \begin{align*}
   B_1 &= \{ \mathbf{x} \in \mathbb{R}^2 : \| \boldsymbol{x} \|_1 = 1 \} \\
   B_2 &= \{ \mathbf{x} \in \mathbb{R}^2 : \| \boldsymbol{x} \|_2 = 1 \} \\
   B_{\infty} &= \{ \mathbf{x} \in \mathbb{R}^2 : \| \boldsymbol{x} \|_{\infty} = 1 \}
   \end{align*}
   $$

3. Which set of inequalities is always true? Explain.

   $$
   \| \boldsymbol{x} \|_1 \leq \| \boldsymbol{x} \|_2 \leq \| \boldsymbol{x} \|_{\infty}
   \hspace{5mm} \text{or} \hspace{5mm}
   \| \boldsymbol{x} \|_1 \geq \| \boldsymbol{x} \|_2 \geq \| \boldsymbol{x} \|_{\infty}
   $$

4. Determine whether the statement is **True** or **False**.
   * If $\boldsymbol{x} \in \mathbb{R}^n$ such that $\| \boldsymbol{x} \|_1 = \| \boldsymbol{x} \|_{\infty} = \lambda$ then $\| \boldsymbol{x} \|_p = \lambda$ for any $p > 1$.
   * Define $\| \boldsymbol{x} \|_0 = \sum_{k=1}^n x_k^2$ for any $\boldsymbol{x} = \begin{bmatrix} x_1 & \cdots & x_n \end{bmatrix}^T \in \mathbb{R}^n$. Then $\| \boldsymbol{x} \|_{0}$ is a norm.
   * Define $\| \boldsymbol{x} \|_{\min} = \min_k | x_k |$ for any $\boldsymbol{x} = \begin{bmatrix} x_1 & \cdots & x_n \end{bmatrix}^T \in \mathbb{R}^n$. Then $\| \boldsymbol{x} \|_{\min}$ is a norm.

5. Suppose $\boldsymbol{x} = \begin{bmatrix} 1 & 2 & 2 & c & 1 & 1 \end{bmatrix}^T \in \mathbb{R}^6$ such that $\| \boldsymbol{x} \|_{3} = 3$. Find all possible values $c$.

6. Suppose $\boldsymbol{x} = \begin{bmatrix} 1 & 3 & 2 & 2 & c & 2 & 1 & 2 \end{bmatrix}^T \in \mathbb{R}^8$ such that $\| \boldsymbol{x} \|_{3} = 5$. Find all possible values $c$.

7. Suppose $\boldsymbol{x} = \begin{bmatrix} 1 & 0 & 2 & c & -1 \end{bmatrix}^T \in \mathbb{R}^5$ such that $\| \boldsymbol{x} \|_1 = 5$. Find all possible values $c$.

8. Suppose $\boldsymbol{x} = \begin{bmatrix} 1 & 0 & 2 & c & -1 \end{bmatrix}^T \in \mathbb{R}^5$ such that $\| \boldsymbol{x} \|_{\infty} = 5$. Find all possible values $c$.

9. Suppose $A$ is a 2 by 2 matrix such that the image of the unit square under the linear transformation $A$ is:

   ![/img/img02.png](/img/img02.png)

   Determine $\mathrm{cond}_{\infty}(A)$ (the condition number with respect to the $\infty$-norm).

10. Suppose $A$ is a 2 by 2 matrix such that the image of the unit "diamond" under the linear transformation $A$ is:

    ![/img/img03.png](/img/img03.png)

    Determine $\mathrm{cond}_1(A)$ (the condition number with respect to the 1-norm).
