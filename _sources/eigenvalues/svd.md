# Singular Value Decomposition

```{div} bigidea
<p>

* Any $m \times n$ matrix $A$ has a singular value decomposition $A = P \Sigma Q^T$ where $P$ and $Q$ are orthogonal matrices and $\Sigma$ is a diagonal $m \times n$ matrix.
* An $n \times p$ data matrix $X$ represents a set of $n$ samples in $\mathbb{R}^p$ and projecting the data onto the first $k$ principal components allows us to view the data in $\mathbb{R}^k$. Usually $k=2$ such that we can visualize the data in 2D.

</p>
```

```{div} note
All matrices in this chapter are *real* unless explicitly stated otherwise.
```

## SVD Construction

```{div} note
If $A$ is any $m \times n$ matrix, then $AA^T$ and $A^TA$ are both symmetric therefore both are orthogonally diagonalizable

$$
AA^T = PD_1P^T \hspace{5mm} A^TA = QD_2Q^T
$$
```

```{div} theorem
Let $A$ be an $m \times n$ matrix.

1. If $\lambda$ is a nonzero eigenvalue of $AA^T$ then $\lambda$ is an eigenvalue of $A^TA$, and vice versa.
2. All eigenvalues of $AA^T$ (and $A^TA$) are non-negative (that is, $\lambda \geq 0$).
```

```{div} definition
The matrices $AA^T$ and $A^TA$ have the same set of positive eigenvalues. Label the eigenvalues in decreasing order $\lambda_1 \geq \lambda_2 \geq \cdots \geq \lambda_r > 0$. The **singular values** of $A$ are

$$
\sigma_i = \sqrt{\lambda_i} \ \  , \ \ i=1,\dots,r
$$
```

```{div} theorem
Let $A$ be an $m \times n$ matrix and let $\sigma_1 \geq \sigma_2 \geq \cdots \geq \sigma_r > 0$ be the singular values of $A$. Then there are orthogonal matrices $P$ and $Q$ such that

$$
A = P \Sigma Q^T
\hspace{5mm}
\text{where}
\hspace{5mm}
\Sigma =
\left[
\begin{array}{ccc|c}
\sigma_1 & & & \\
& \ddots & & \boldsymbol{0} \\
& & \sigma_r & \\ \hline
& \boldsymbol{0} & & \boldsymbol{0}
\end{array} \right]_{m \times n}
$$

This is called the **singular value decomposition** of $A$.

---

*Proof*. Let $\boldsymbol{q}_1,\dots,\boldsymbol{q}_n$ be orthonormal eigenvectors of $A^TA$ chosen in order such that

$$
A^TA \boldsymbol{q}_i = \sigma_i^2 \boldsymbol{q}_i \ \ , \ \ i =1,\dots,r
\hspace{15mm}
A^TA \boldsymbol{q}_i = \boldsymbol{0} \ \ , \ \ i =r+1,\dots,n
$$

Note that in fact $A\boldsymbol{q}_i = \boldsymbol{0}$ for $i=r+1,\dots,n$ since

$$
\| A\boldsymbol{q}_i \|^2 = \boldsymbol{q}_i^T A^T A\boldsymbol{q}_i = 0 \ \ , \ \ i=r+1,\dots,n
$$

Let $Q$ be the orthogonal matrix

$$
Q = \begin{bmatrix} & & \\ \boldsymbol{q}_1 & \cdots & \boldsymbol{q}_n \\ & & \end{bmatrix}
$$

Now construct the matrix $P$. Let

$$
\boldsymbol{p}_i = \frac{1}{\sigma_i} A \boldsymbol{q}_i \ \ , \ \ i = 1,\dots,r
$$

Note that

$$
AA^T \boldsymbol{p}_i =  AA^T \left( \frac{1}{\sigma_i} A \boldsymbol{q}_i \right)
= \frac{1}{\sigma_i} A \left( A^TA \boldsymbol{q}_i \right)
= \sigma_i A \boldsymbol{q}_i
= \sigma_i^2 \boldsymbol{p}_i
$$

therefore each $\boldsymbol{p}_i$ is an eigenvector for $AA^T$ with eigenvalue $\sigma_i^2$. Note also that

$$
\| \boldsymbol{p}_i \|^2 = \boldsymbol{p}_i^T \boldsymbol{p}_i
= \frac{1}{\sigma_i^2}  \boldsymbol{q}_i^T A^T A \boldsymbol{q}_i
= \boldsymbol{q}_i^T \boldsymbol{q}_i = 1
$$

therefore each $\boldsymbol{p}_i$ is a unit vector. Extend (by Gram-Schmidt algorithm) to an orthonormal basis $
\boldsymbol{p}_1,\dots,\boldsymbol{p}_r,\boldsymbol{p}_{r+1},\dots,\boldsymbol{p}_m$ of $\mathbb{R}^m$. Define the orthogonal matrix

$$
P = \begin{bmatrix} & & \\ \boldsymbol{p}_1 & \cdots & \boldsymbol{p}_m \\ & & \end{bmatrix}
$$

Compute

$$
AQ = \begin{bmatrix} & & & & & \\ A \boldsymbol{q}_1 & \cdots & A\boldsymbol{q}_r & A\boldsymbol{q}_{r+1} & \cdots & A\boldsymbol{q}_n \\ & & & & & \end{bmatrix} = \begin{bmatrix} & & & & & \\ \sigma_1 \boldsymbol{p}_1 & \cdots & \sigma_r \boldsymbol{p}_r & \boldsymbol{0} & \cdots & \boldsymbol{0} \\ & & & & & \end{bmatrix}
$$

Finally, compute

$$
P \Sigma = \begin{bmatrix} & & & & & \\ \boldsymbol{p}_1 & \cdots & \boldsymbol{p}_r & \boldsymbol{p}_{r+1} & \cdots & \boldsymbol{p}_n \\ & & & & & \end{bmatrix} \left[
\begin{array}{ccc|c}
\sigma_1 & & & \\
& \ddots & & \boldsymbol{0} \\
& & \sigma_r & \\ \hline
& \boldsymbol{0} & & \boldsymbol{0}
\end{array} \right]
=
\begin{bmatrix} & & & & & \\ \sigma_1 \boldsymbol{p}_1 & \cdots & \sigma_r \boldsymbol{p}_r & \boldsymbol{0} & \cdots & \boldsymbol{0} \\ & & & & & \end{bmatrix}
$$

Therefore $AQ = P \Sigma$ and so $A = P \Sigma Q^T$.
```

```{div} note
In the construction of the SVD, we may chose to first construct either $P$ or $Q$. The connection between the columns for $i = 1,\dots,r$ are given by the equations:

$$
\begin{align*}
\boldsymbol{q}_i &= \frac{1}{\sigma_i} A^T \boldsymbol{p}_i & A^TA \boldsymbol{q}_i &= \sigma_i^2 \boldsymbol{q}_i & \| \boldsymbol{q}_i \| &= 1 \\
\boldsymbol{p}_i &= \frac{1}{\sigma_i} A \boldsymbol{q}_i & AA^T \boldsymbol{p}_i &= \sigma_i^2 \boldsymbol{p}_i & \| \boldsymbol{p}_i \| &= 1
\end{align*}
$$
```

```{div} example
Construct the SVD for

$$
A = \left[ \begin{array}{rr} 1 & 1 \\ 1 & -1 \\ 0 & 1 \end{array} \right]
$$

Since $A^TA$ is a smaller matrix, let us first construct $Q$. Compute

$$
A^TA = \begin{bmatrix} 2 & 0 \\ 0 & 3 \end{bmatrix}
$$

Therefore $\sigma_1 = \sqrt{3}$ and $\sigma_2 = \sqrt{2}$. By inspection, we find

$$
\boldsymbol{q}_1 = \begin{bmatrix} 0 \\ 1 \end{bmatrix}
\hspace{5mm}
\boldsymbol{q}_2 = \begin{bmatrix} 1 \\ 0 \end{bmatrix}
\hspace{5mm}
\Rightarrow
\hspace{5mm}
Q = \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}
$$

Construct the matrix $P$

$$
\boldsymbol{p}_1 = \frac{1}{\sigma_1} A\boldsymbol{q}_1 = \frac{1}{\sqrt{3}} \left[ \begin{array}{r} 1 \\ -1 \\ 1 \end{array} \right]
\hspace{5mm}
\boldsymbol{p}_2 = \frac{1}{\sigma_2} A\boldsymbol{q}_2 = \frac{1}{\sqrt{2}} \left[ \begin{array}{r} 1 \\ 1 \\ 0 \end{array} \right]
$$

Extend to an orthonormal basis of $\mathbb{R}^3$ by finding $\boldsymbol{p}_3$ orthogonal to $\boldsymbol{p}_1$ and $\boldsymbol{p}_2$. There are different ways of doing this. Setup equations $\langle \boldsymbol{p}_1 , \boldsymbol{p}_3 \rangle = 0$ and $\langle \boldsymbol{p}_2 , \boldsymbol{p}_3 \rangle = 0$ in a linear system and solve

$$
\left[ \begin{array}{rrr|r} 1 & -1 & \phantom{+}1 & 0 \\ 1 & 1 & 0 & 0 \end{array} \right]
\hspace{5mm}
\Rightarrow
\hspace{5mm}
\boldsymbol{p}_3 = \frac{1}{\sqrt{6}} \left[ \begin{array}{r} -1 \\ 1 \\ 2 \end{array} \right]
$$

Therefore the SVD is given by

$$
A = P \Sigma Q^T
=
\left[ \begin{array}{rcr} 1/\sqrt{3} & 1/\sqrt{2} & -1/\sqrt{6} \\ -1/\sqrt{3} & 1/\sqrt{2} & 1/\sqrt{6} \\ 1/\sqrt{3} & 0 & 2/\sqrt{6} \end{array} \right]
\left[ \begin{array}{rr} \sqrt{3} & 0 \\ 0 & \sqrt{2} \\ 0 & 0 \end{array} \right]  \begin{bmatrix} 0 & 1 \\ 1 & 0 \end{bmatrix}^T
$$
```

```{div} theorem
Let $A = P \Sigma Q^T$.

* $\mathrm{rank}(A) = r$
* $\| A \| = \sigma_1$
* $\| A^{-1} \| = 1/\sigma_r$
* $\mathrm{cond}(A) = \sigma_1 / \sigma_r$
* $\mathrm{null}(A) = \mathrm{span} \{ \boldsymbol{q}_{r+1},\dots,\boldsymbol{q}_n \}$
* $\mathrm{range}(A) = \mathrm{span} \{ \boldsymbol{p}_1,\dots,\boldsymbol{p}_r \}$
```

## Principal Component Analysis

```{div} definition
Let $\boldsymbol{x}_1 , \dots , \boldsymbol{x}_n \in \mathbb{R}^p$ (viewed as row vectors) and let $X$ be the $n \times p$ data matrix where row $k$ is given by $\boldsymbol{x}_k$. Assume the data is **normalized** such that the mean value of each column of $X$ is 0. The unit vector $\boldsymbol{w}_1$ which maximizes the sum

$$
\sum_{i=1}^n \langle \boldsymbol{x}_i , \boldsymbol{w}_1 \rangle^2 = \| X \boldsymbol{w}_1 \|^2
$$

is called the **first weight vector** of $X$ (see [Wikipedia:Principal component analysis](https://en.wikipedia.org/wiki/Principal_component_analysis)). More generally, given weight vectors $\boldsymbol{w}_1 , \dots, \boldsymbol{w}_{k-1}$, the **$k$th weight vector** of $X$ is the unit vector $\boldsymbol{w}_k$ which maximizes

$$
\| X_k \boldsymbol{w}_k \|^2
$$

where $X_k$ is the projection of the data matrix $X$ onto $\mathrm{span} \{ \boldsymbol{w}_1 , \dots, \boldsymbol{w}_{k-1} \}^{\perp}$

$$
X_k = X - \sum_{j=1}^{k-1} X \boldsymbol{w}_j \boldsymbol{w}_j^T
$$

The projection coefficient $\langle \boldsymbol{x}_i , \boldsymbol{w}_k \rangle$ is called the **$k$th principal component** of a data vector $\boldsymbol{x}_i$.
```

```{div} note
Each $\langle \boldsymbol{x}_k , \boldsymbol{w}_1 \rangle^2$ is the length squared of the orthogonal projection of $\boldsymbol{x}_k$ onto $\boldsymbol{w}_1$. Therefore the first weight vector $\boldsymbol{w}_1$ points in the direction which captures the most information (ie. the maximum variance) of the data, and the second weight vector $\boldsymbol{w}_2$ is orthogonal to $\boldsymbol{w}_1$.

```{image} /img/03_02_01.png
:width: 500px
:align: center
```

```{div} theorem
The weight vectors $\boldsymbol{w}_1 , \dots , \boldsymbol{w}_p$ are the right singular vectors of the matrix $X$. In other words, let $X = P \Sigma Q^T$ be a singular value decomposition of $X$ and let $\boldsymbol{q}_1, \dots, \boldsymbol{q}_p$ be the columns of $Q$ corresponding to singular values $\sigma_1 > \cdots > \sigma_p > 0$. Then $\boldsymbol{w}_1 = \boldsymbol{q}_1 , \dots , \boldsymbol{w}_p = \boldsymbol{q}_p$.

---

*Proof*. Let $X = P \Sigma Q^T$ be a singular value decomposition of $X$ and note

$$
\| X \boldsymbol{w} \|^2 = \| P \Sigma Q^T \boldsymbol{w} \|^2 = \| \Sigma Q^T \boldsymbol{w} \|^2
$$

since $P$ is orthogonal. Since $\Sigma$ is diagonal with diagonal entires $\sigma_1 > \cdots > \sigma_p$, the maximum value of $\| X \boldsymbol{w} \|^2$ occurs when $Q^T \boldsymbol{w} = \begin{bmatrix} 1 & 0 & \cdots & 0 \end{bmatrix}^T$ therefore $\boldsymbol{w}_1 = \boldsymbol{q}_1$. For general $k$, note that the singular value decomposition $X_k = P_k \Sigma_k Q^T_k$ is obtained from $X$ by removing the singular values $\sigma_1,\dots,\sigma_{k-1}$. Therefore the largest singular value of $X_k$ is $\sigma_k$ with corresponding right singular vector $\boldsymbol{q}_k$ and therefore $\boldsymbol{w}_k = \boldsymbol{q}_k$.
```

````{div} example
Find the first weight vector for the data given in the image below.

```{image} /img/03_02_02.png
:width: 500px
:align: center
```

We expect $\boldsymbol{w}_1 = \begin{bmatrix} 1/\sqrt{2} & 1/\sqrt{2} \end{bmatrix}^T$ since that direction clearly captures the most information. Form the data matrix

$$
X^T = \left[ \begin{array}{rrrrrrr} -2 & -1 & -1 & \phantom{+}0 & 1 & \phantom{+}1 & 2 \\ -2 & -1 & 1 & 0 & -1 & 1 & 2 \end{array} \right]
$$

We don't need to compute the full SVD of $X$ but just the first right singular vector. Compute

$$
X^T X = \begin{bmatrix} 12 & 8 \\ 8 & 12 \end{bmatrix}
$$

The characteristic polynomial of $X^TX$ is

$$
\det (xI - X^TX) = (x-12)^2 - 8^2 = x^2 - 24x + 80 = (x - 4)(x - 20)
$$

The right singular vector $\boldsymbol{q}_1$ for $X$ is a unit eigenvector for $X^TX$ for the eigenvalue $\lambda_1 = 20$. Compute

$$
( X^T X - 20 I)\boldsymbol{w}_1 = \boldsymbol{0} \ \ \Rightarrow \ \
\left[ \begin{array}{rr|r} -8 & 8 & 0 \\ 8 & -8 & \phantom{+}0 \end{array} \right]
\ \ \Rightarrow \ \
\boldsymbol{w}_1 = \frac{1}{\sqrt{2}} \begin{bmatrix} 1 \\ 1 \end{bmatrix}
$$
````

````{div} example
The [digits dataset from sklearn](https://scikit-learn.org/stable/modules/generated/sklearn.datasets.load_digits.html#sklearn.datasets.load_digits) is a $1797 \times 64$ data matrix $X$ such that each row represents an $8 \times 8$ pixel image of a handwritten number. The first 10 rows of X (reshaped from vectors of length 64 to $8 \times 8$ matrices to visualize) are:

```{image} /img/03_02_03.png
:width: 500px
:align: center
```

Compute the first 2 weight vectors and find (again $\boldsymbol{w}_1,\boldsymbol{w}_2$ reshaped from vectors of length 64 to $8 \times 8$ matrices to visualize)

```{image} /img/03_02_04.png
:width: 500px
:align: center
```

We can see $\boldsymbol{w}_1$ looks like a 3 and $\boldsymbol{w}_2$ looks like 0. Project the entire dataset onto these weight vectors and label each data point by a color according to the digit:

```{image} /img/03_02_05.png
:width: 500px
:align: center
```

We can see that the 3s are to the right in the horizontal direction since these points most similar to $\boldsymbol{w}_1$, and the 0s are at the top in the vertical direction since these points most similar to $\boldsymbol{w}_2$. We can make other interesting observations such as the 4s are opposite to the 3s and orthogonal to 0s, and 7s and 1s are opposite to 0s and orthogonal to 3s.
````

## Pseudoinverse, Least Squares and the SVD Expansion

```{div} bigidea
The pseudoinverse $A^+$ solves the least squares problem $A \boldsymbol{x} \cong \boldsymbol{b}$ by $\boldsymbol{x} = A^+ \boldsymbol{b}$.
```

```{div} definition
Let $A$ be a $m \times n$ matrix of rank $r$ and let $A = P \Sigma Q^T$ be a SVD of $A$. The **pseudoinverse** of $A$ is

$$
A^+ = Q \Sigma^+ P^T
$$

where

$$
\Sigma^+ =
\renewcommand{\arraystretch}{1.25}
\left[ \begin{array}{ccc|c}
\sigma_1^{-1} & & & \\
& \ddots & & \boldsymbol{0} \\
& & \sigma_r^{-1} & \\ \hline
& \boldsymbol{0} & & \boldsymbol{0}
\end{array} \right]_{n \times m}
\renewcommand{\arraystretch}{1}
$$
```

```{div} note
If $A$ is invertible, then $A^+ = A^{-1}$.
```

```{div} theorem
Let $A$ be an $m \times n$ matrix with $\mathrm{rank}(A) = n$ and let $\boldsymbol{b} \in \mathbb{R}^m$. The least squares approximation of the system $A \boldsymbol{x} \cong \boldsymbol{b}$ is given by $\boldsymbol{x} = A^+ \boldsymbol{b}$.

---

*Proof*. Let $A = P \Sigma Q^T$ be a SVD of $A$. Let $P^T \boldsymbol{b} = \boldsymbol{c}$ and write

$$
\boldsymbol{c} = \begin{bmatrix}
c_1 \\ \vdots \\ c_m \end{bmatrix}
\in \mathbb{R}^m
$$

Since $P$ is an orthogonal matrix we have

$$
\| A \boldsymbol{x} - \boldsymbol{b} \| = \| P \Sigma Q^T \boldsymbol{x} - P \boldsymbol{c} \| = \| \Sigma Q^T \boldsymbol{x} - \boldsymbol{c} \|
$$

The matrix $\Sigma$ is of the form

$$
\Sigma =
\renewcommand{\arraystretch}{1.25}
\left[ \begin{array}{ccc}
\sigma_1 & & \\
& \ddots & \\
& & \sigma_n \\ \hline
\boldsymbol{0} & \cdots & \boldsymbol{0}
\end{array} \right]_{m \times n}
\renewcommand{\arraystretch}{1}
$$

and so only the first $n$ entries of $\Sigma \boldsymbol{v}$ are nonzero for any vector $\boldsymbol{v} \in \mathbb{R}^n$. Therefore the minimum value $\| A \boldsymbol{x} - \boldsymbol{b} \| = \| \Sigma Q^T \boldsymbol{x} - \boldsymbol{c} \|$ occurs when

$$
\Sigma Q^T \boldsymbol{x}
= \begin{bmatrix}
c_1 \\ \vdots \\ c_n \\ \boldsymbol{0} \end{bmatrix}
$$

and so $\boldsymbol{x} = Q \Sigma^+ \boldsymbol{c}$. Altogether, we have

$$
\boldsymbol{x} = Q \Sigma^+ P^T \boldsymbol{b} = A^+ \boldsymbol{b}
$$
```

```{div} definition
Let $A = P \Sigma Q^T$. The **SVD expansion** of $A$ is

$$
A = \sum_{i=1}^r \sigma_i \boldsymbol{p}_i \boldsymbol{q}_i^T = \sigma_1 \boldsymbol{p}_1 \boldsymbol{q}_1^T + \cdots + \sigma_r \boldsymbol{p}_r \boldsymbol{q}_r^T
$$

Note that each product $\boldsymbol{p}_i \boldsymbol{q}_i^T$ is a $m \times n$ matrix of rank 1.
```

```{div} theorem
Let $A$ be a $m \times n$ matrix and let $B$ be a $n \times \ell$ matrix with SVD expansions

$$
A = \sum_{i=1}^r \sigma_i \boldsymbol{p}_i \boldsymbol{q}_i^T
\hspace{5mm} \text{and} \hspace{5mm}
B = \sum_{j=1}^s \mu_j \boldsymbol{u}_j \boldsymbol{v}_j^T
$$

If $\{ \boldsymbol{q}_1,\dots,\boldsymbol{q}_r \}$ and $\{ \boldsymbol{u}_1,\dots,\boldsymbol{u}_s \}$ are orthogonal sets of vectors, then $AB = 0$.
```

```{div} definition
Let $A = P \Sigma Q^T$. The **truncated SVD expansion of rank $k$** of $A$ is

$$
A_k = \sum_{i=1}^k \sigma_i \boldsymbol{p}_i \boldsymbol{q}_i^T = \sigma_1 \boldsymbol{p}_1 \boldsymbol{q}_1^T + \cdots + \sigma_k \boldsymbol{p}_k \boldsymbol{q}_k^T
$$
```

```{div} note
Suppose we want to solve the system $A \boldsymbol{x} = \boldsymbol{b}$ however the right side is corrupted by noise $\boldsymbol{e}$ and we must work with the system

$$
A \hat{\boldsymbol{x}} = \boldsymbol{b} + \boldsymbol{e}
$$

Solving directly we get

$$
\hat{\boldsymbol{x}} =  A^{-1} \left( \boldsymbol{b} + \boldsymbol{e} \right)= A^{-1}\boldsymbol{b} + A^{-1}\boldsymbol{e}
$$

and the term $A^{-1}\boldsymbol{e}$ is called the **inverted noise** which may dominate the true solution $\boldsymbol{x} = A^{-1} \boldsymbol{b}$. From the SVD expansion, we see that most of $A$ is composed of the terms $\sigma_i  \boldsymbol{p}_i \boldsymbol{q}_i^T$ for *large* singular values $\sigma_i$. If we know that the error $\boldsymbol{e}$ is unrelated to $A$ in the sense that $\boldsymbol{e}$ is (mostly) orthogonal to the singular vectors $\boldsymbol{p}_i$ of $A$ corresponding to large singular values, then the truncated SVD expansion of the pseudoinverse

$$
A_k^+ = \sum_{i=1}^k \frac{1}{\sigma_i} \boldsymbol{q}_i \boldsymbol{p}_i^T
$$

gives a better solution

$$
\hat{\boldsymbol{x}} =  A^+_k \left( \boldsymbol{b} + \boldsymbol{e} \right)= A^+_k \boldsymbol{b} + A^+_k\boldsymbol{e}
$$

since the term $A^+_k\boldsymbol{e}$ will be smaller. In other words, we avoid terms $\sigma_i^{-1} \boldsymbol{p}_i \boldsymbol{q}_i^T$ in the SVD expansion of $A^{-1}$ for small singular values $\sigma_i$ which produce large values $\sigma_i^{-1}$ which may amplify the error. This is the strategy for image deblurring and computed tomography in the next sections.
```

## Exercises

**Exercise 1.** Find the singular value decomposition of the matrix

$$
A = \left[ \begin{array}{rrr} 1 & \ \, 2 & -1 \\ 2 & 1 & 4 \end{array} \right]
$$

**Exercise 2.** Find the singular value decomposition of the matrix

$$
A = \left[ \begin{array}{rrr} 1 & 1 & 1 \\ -1 & 2 & -1 \\ 1 & 0 & -1 \end{array} \right]
$$

**Exercise 3.** Determine whether the statement is **True** or **False**. (Assume all data matrices are normalized.)

* Let $X$ be a $n \times p$ data matrix and let $\boldsymbol{x}_i , \boldsymbol{x}_j \in \mathbb{R}^p$ be two different rows of $X$ such that $\| \boldsymbol{x}_i \| < \| \boldsymbol{x}_j \|$. If $\boldsymbol{w}_1$ is the first weight vector of $X$, then $| \langle \boldsymbol{x}_i , \boldsymbol{w}_1 \rangle | < | \langle \boldsymbol{x}_j , \boldsymbol{w}_1 \rangle |$.
* Let $X$ be a $n \times p$ data matrix and let $\boldsymbol{x}_i , \boldsymbol{x}_j \in \mathbb{R}^p$ be two different rows of $X$ such that $\langle \boldsymbol{x}_i , \boldsymbol{x}_j \rangle = 0$. If $\boldsymbol{w}_1$ is the first weight vector of $X$ and $\langle \boldsymbol{x}_i , \boldsymbol{w}_1 \rangle = 0$ then $\langle \boldsymbol{x}_j , \boldsymbol{w}_1 \rangle = 0$.
* Let $X$ be a $n \times 2$ data matrix and let $Y$ be the matrix with the same columns as $X$ but switched. (In other words, the first column of $Y$ is the same as the second column of $X$, and the second column of $Y$ is the first column of $X$.) If $X$ and $Y$ represent the same set of data points, then all the singular values of $X$ equal.
* Let $X$ be a $n \times 2$ data matrix and let $Y$ be the matrix with the same columns as $X$ but switched. (In other words, the first column of $Y$ is the same as the second column of $X$, and the second column of $Y$ is the first column of $X$.) If $X$ and $Y$ represent the same set of data points, then $\boldsymbol{w}_1 = \begin{bmatrix} 1/\sqrt{2} & 1/\sqrt{2} \end{bmatrix}^T$.

**Exercise 4.** Find the weight vectors for the data matrix $X$ representing the points:

```{image} /img/03_02_06.png
:width: 500px
:align: center
```

**Exercise 5.** Suppose $X$ is a $100 \times 4$ data matrix such that

$$
X^T X = \begin{bmatrix} 2 & 0 & 0 & 0 \\ 0 & 1.5 & 0 & 0 \\ 0 & 0 & 2 & 1 \\ 0 & 0 & 1 & 2 \end{bmatrix}
$$

Find all the weight vectors of $X$.

**Exercise 6.** Suppose we want to solve a system $A \boldsymbol{x} = \boldsymbol{b}$. A small change $\Delta \boldsymbol{b}$ produces a change in the solution

$$
A(\boldsymbol{x} + \Delta \boldsymbol{x}) = \boldsymbol{b} + \Delta \boldsymbol{b}
$$

Describe the unit vector $\Delta \boldsymbol{b}$ that will produce the largest change $\| \Delta \boldsymbol{x} \|$.

**Exercise 7.** Find the rank 2 pseudo inverse

$$
A_2^+ = \frac{1}{\sigma_1} \boldsymbol{q}_1 \boldsymbol{p}_1^T + \frac{1}{\sigma_2} \boldsymbol{q}_2 \boldsymbol{p}_2^T
$$

of the matrix

$$
A = \left[ \begin{array}{rrr} 1 & 1 & 1 \\ 1 & 0 & -2 \\ 1 & -1 & 1 \end{array} \right]
$$

(Note: the columns of $A$ are orthogonal.)

**Exercise 8.** Let $A$ be a $m \times n$ matrix with singular value decomposition $A = P \Sigma Q^T$. Let $k < \min\{m,n\}$ and let

$$
A_k = \sum_{i=1}^k \sigma_i \boldsymbol{p}_i \boldsymbol{q}_i^T
$$

Describe the singular value decomposition of $A - A_k$.
