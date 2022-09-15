# Matrix Multiplication

```{div} theorem
Let $A$ and $B$ be $m \times n$ matrices. If $A \boldsymbol{x} = B \boldsymbol{x}$ for all $\boldsymbol{x} \in \mathbb{R}^n$, then $A = B$.

---

*Proof.* The standard basis of $\mathbb{R}^n$ is

$$
\boldsymbol{e}_1 = \begin{bmatrix} 1 \\ 0 \\ \vdots \\ 0 \end{bmatrix}
\ \ , \ \
\boldsymbol{e}_2 = \begin{bmatrix} 0 \\ 1 \\ \vdots \\ 0 \end{bmatrix}
\ \ , \ \ \dots \ \ , \ \
\boldsymbol{e}_n = \begin{bmatrix} 0 \\ 0 \\ \vdots \\ 1 \end{bmatrix}
$$

In other words, $\boldsymbol{e}_k$ is the vector with 1 at index $k$ and 0 everywhere else. Then $A \boldsymbol{e}_k$ is equal to the $k$th column of $A$. Since $A \boldsymbol{e}_k = B \boldsymbol{e}_k$ for each $k=1,\dots,n$ we see that the columns of $A$ and $B$ are equal therefore $A = B$.
```

```{div} theorem
Let $A$ be a $p \times m$ matrix, let $B$ be a $m \times n$ matrices and let $\boldsymbol{b}_1,\dots,\boldsymbol{b}_n$ be the columns of $B$

$$
B = \begin{bmatrix} & & \\ \boldsymbol{b}_1 & \cdots & \boldsymbol{b}_n \\ & & \end{bmatrix}
$$

Then the $k$th column of $AB$ is $A \boldsymbol{b}_k$. In other words, matrix multiplication $AB$ can be written as

$$
AB = \begin{bmatrix} & & \\ A \boldsymbol{b}_1 & \cdots & A \boldsymbol{b}_n \\ & & \end{bmatrix}
$$

```