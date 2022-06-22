# Interpolation

```{div} bigidea
An interpolating function for a set of data provides information about values between data points and beyond the range of the data. There are infinitely many different ways to interpolate a set of data. Polynomial interpolation and cubic spline interpolation are simple methods of interpolation.
```

## Interpolating Functions

```{div} definition
An **interpolant** (or **interpolating function**) for data points $(t_0,y_0),\dots,(t_d,y_d)$ is a function $f(t)$ such that $f(t_k) = y_k$ for $k=0,\dots,d$.
```

```{div} note
There are *infinitely* many different interpolating functions for any set of data. We impose different kinds of constraints on the form of the interpolant to ensure there is a unique solution which is meaningful for the data and can be computed efficiently. In the following sections we describe two kinds of interpolation: *polynomial interpolation* and *cubic spline interpolation*.  
```

## Polynomial Interpolation

```{div} bigidea
Given points $(t_0,y_0), \dots , (t_d,y_d)$, there exists a unique polynomial $p(t)$ of degree (at most) $d$ such that $p(t_k) = y_k$ for each $k=0,\dots,d$. Monomial, Lagrange and Newton interpolation provide different ways to compute the exact same interpolating polynomial.
```

```{div} definition
A **polynomial** of degree (at most) $d$ is a function of the form

$$
p(t)= c_0 + c_1 t + \cdots + c_d t^d
$$

where $c_0,c_1,\dots,c_d \in \mathbb{R}$. The collection of all polynomials of degree (at most) $d$ is denoted by

$$
\mathbb{P}_d = \{ c_0 + c_1 t + \cdots + c_d t^d : c_0,c_1,\dots,c_d \in \mathbb{R} \}
$$

Note that $\mathbb{P}_d$ is a vector space of dimension $d+1$.
```

### Monomial Basis

```{div} definition
Fix a positive integer $d$ and let $\phi_k(t) = t^k$ for $k=0,\dots,d$. In other words, consider the polynomials

$$
\phi_0(t) = 1 \ , \ \phi_1(t) = t \ , \ \dots \ , \ \phi_d(t) = t^d
$$

Since a polynomial $p(t)$ of degree (at most) $d$ is of the form

$$
p(t) = c_0 \phi_0(t) + c_1 \phi_1(t) + \cdots + c_d \phi_d(t)
$$

we see that the polynomials $\phi_k(t)$ form a basis of $\mathbb{P}_d$. The set $\{ \phi_0(t),\dots,\phi_d(t) \}$ is called the **monomial basis** of $\mathbb{P}_d$.
```

```{div} definition
Given $d+1$ points $(t_0,y_0), \dots,(t_d,y_d)$, **polynomial interpolation with respect to the monomial basis** seeks a polynomial of the form

$$
p(t) = c_0 + c_1 t + \cdots + c_d t^d
$$

such that $p(t_k) = y_k$ for each $k=0,\dots,d$. Each point defines an equation

$$
\begin{align*}
c_0 + c_1t_0 + \cdots + c_d t_0^d &= y_0 \\
c_0 + c_1t_1 + \cdots + c_d t_1^d &= y_1 \\
& \ \vdots \\
c_0 + c_1t_d + \cdots + c_d t_d^d &= y_d
\end{align*}
$$

This yields a linear system of equations $A \boldsymbol{c} = \boldsymbol{y}$

$$
\begin{bmatrix}
1 & t_0 & \cdots & t_0^d \\
1 & t_1 & \cdots & t_1^d \\
\vdots & \vdots & \ddots & \vdots \\
 1 & t_d & \cdots & t_d^d
\end{bmatrix}
\begin{bmatrix} c_0 \\ c_1 \\ \vdots \\ c_d \end{bmatrix}
=
\begin{bmatrix} y_0 \\ y_1 \\ \vdots \\ y_d \end{bmatrix}
$$

The matrix

$$
A = \begin{bmatrix}
1 & t_0 & \cdots & t_0^d \\
1 & t_1 & \cdots & t_1^d \\
\vdots & \vdots & \ddots & \vdots \\
 1 & t_d & \cdots & t_d^d
\end{bmatrix}
$$

is called a **Vandermonde matrix**. Solve the system $A \boldsymbol{c} = \boldsymbol{y}$ to find the coefficients

$$
\boldsymbol{c} = \begin{bmatrix} c_0 \\ c_1 \\ \vdots \\ c_d \end{bmatrix}
$$
```

```{div} note
The condition number of a Vandermonde matrix gets *very* large as the size of the matrix increases. This means that interpolation by the monomial basis is very sensitive to changes in the data for polynomials of large degree. For example, for $11$ equally spaced points $t_0=0,\dots,t_{10}=10$, the Vandermonde matrix $A$ is 11 by 11 and has condition number larger than $10^{12}$. Yikes!

```{image} /img/01_03_01.png
:width: 500px
:align: center
```

```{div} theorem
Let $A$ be the Vandermonde matrix for points $(t_0,y_0), \dots,(t_d,y_d)$. Then

$$
\mathrm{det}(A) = \prod_{0 \leq i < j \leq d} (t_j - t_i)
$$
```

```{div} example
Let's compute the formula for the determinant of the Vandermonde matrix for the cases $d = 1$ and $d = 2$. Compute for $d = 1$

$$
\mathrm{det} \left( \begin{bmatrix} 1 & t_0 \\ 1 & t_1 \end{bmatrix} \right) = t_1 - t_0
$$

For $d=2$, use $\mathrm{det}(A) = \mathrm{det}(U)$ from the LU decomposition. Compute

$$
\begin{bmatrix} 1 & t_0 & t_0^2 \\ 1 & t_1 & t_1^2 \\ 1 & t_2 & t_2^2 \end{bmatrix}
\longrightarrow
\begin{bmatrix} 1 & t_0 & t_0^2 \\ 0 & t_1 - t_0 & t_1^2 - t_0^2 \\ 0 & t_2 - t_0 & t_2^2 - t_0^2 \end{bmatrix}
\longrightarrow
\begin{bmatrix} 1 & t_0 & t_0^2 \\ 0 & t_1 - t_0 & t_1^2 - t_0^2 \\ 0 & 0 & * \end{bmatrix}
$$

where

$$
* = t_2^2 - t_0^2 - \frac{t_2 - t_0}{t_1 - t_0} (t_1^2 - t_0^2) = (t_2 - t_0)(t_2 - t_1)
$$

The determinant is the product of the diagonal entries of $U$

$$
\mathrm{det} \left( \begin{bmatrix} 1 & t_0 & t_0^2 \\ 1 & t_1 & t_1^2 \\ 1 & t_2 & t_2^2 \end{bmatrix}
 \right)
= (t_1 - t_0)(t_2 - t_0)(t_2 - t_1)
$$
```

```{div} theorem
Consider $d+1$ data points $(t_0,y_0), \dots , (t_d,y_d)$ (such that $t_i \not= t_j$ for $i \not= j$). There exists a unique polynomial $p(t)$ of degree (at most) $d$ such that $p(t_k) = y_k$ for each $k=0,\dots,d$.

---

*Proof*. The determinant of the Vandermonde matrix is nonzero when the values $t_k$ are distinct therefore it is invertible and so there is a unique solution of the system $A \boldsymbol{c} = \boldsymbol{y}$.
```

### Lagrange Basis

```{div} definition
Given $d+1$ points $(t_0,y_0),\dots,(t_d,y_d)$, the **Lagrange basis** of $\mathbb{P}_d$ (with respect to $t_0,\dots,t_d$) is given by $\{ \ell_0(t),\dots,\ell_d(t) \}$ where

$$
\ell_k(t) = \frac{\prod_{j=0, j\not=k}^d (t - t_j)}{\prod_{j=0, j\not=k}^d (t_k - t_j)}
$$

The essential property of these polynomials is

$$
\ell_k(t_j) = \left\{ \begin{array}{cc} 1 & \text{ if } k = j \\ 0 & \text{ if } k \not= j \end{array} \right. \ , \ \ k,j=0,\dots,d
$$
```

```{div} definition
Given $d+1$ points $(t_0,y_0), \dots,(t_d,y_d)$, **polynomial interpolation with respect to the Lagrange basis** seeks a polynomial of the form

$$
p(t) = c_0 \ell_0(t) + c_1 \ell_1(t) + \cdots + c_d \ell_d(t)
$$

such that $p(t_k) = y_k$ for each $k=0,\dots,d$. Each point defines an equation and we get a linear system $A \boldsymbol{c} = \boldsymbol{y}$

$$
\begin{bmatrix}
\ell_0(t_0) & \ell_1(t_0) & \cdots & \ell_d(t_0) \\
\ell_0(t_1) & \ell_1(t_1) & \cdots & \ell_d(t_1) \\
\vdots & \vdots & \ddots & \vdots \\
\ell_0(t_d) & \ell_1(t_d) & \cdots & \ell_d(t_d)
\end{bmatrix}
\begin{bmatrix} c_0 \\ c_1 \\ \vdots \\ c_d \end{bmatrix}
=
\begin{bmatrix} y_0 \\ y_1 \\ \vdots \\ y_d \end{bmatrix}
$$

The essential property of Lagrange polynomials implies that $A$ is the identity matrix and so $\boldsymbol{c} = \boldsymbol{y}$ therefore

$$
p(t) = y_0 \ell_0(t) + y_1 \ell_1(t) + \cdots + y_d \ell_d(t)
$$
```

### Newton Basis

```{div} definition
Given $d+1$ points $(t_0,y_0),\dots,(t_d,y_d)$, the **Newton basis** of $\mathbb{P}_d$ (with respect to $t_0,\dots,t_d$) is $\{ \pi_0(t),\dots,\pi_d(t) \}$ where

$$
\begin{align*}
\pi_0(t) &= 1 \\
\pi_1(t) &= t - t_0 \\
\pi_2(t) &= (t - t_0)(t - t_1) \\
& \ \ \vdots \\
\pi_d(t) &= (t - t_0)(t - t_1)(t - t_2) \cdots (t - t_{d-1}) \\
\end{align*}
$$

An important property of these polynomials is

$$
\pi_k(t_j) = 0 \ , \ \ j=0,\dots,k-1
$$
```

```{div} definition
Given $d+1$ points $(t_0,y_0), \dots,(t_d,y_d)$, **polynomial interpolation with respect to the Newton basis** seeks a polynomial of the form

$$
p(t) = c_0 \pi_0(t) + c_1 \pi_1(t) + \cdots + c_d \pi_d(t)
$$

such that $p(t_k) = y_k$ for each $k=0,\dots,d$. Each point defines an equation and we get a linear system $A \boldsymbol{c} = \boldsymbol{y}$

$$
\begin{bmatrix}
\pi_0(t_0) & \pi_1(t_0) & \cdots & \pi_d(t_0) \\
\pi_0(t_1) & \pi_1(t_1) & \cdots & \pi_d(t_1) \\
\vdots & \vdots & \ddots & \vdots \\
\pi_0(t_d) & \pi_1(t_d) & \cdots & \pi_d(t_d)
\end{bmatrix}
\begin{bmatrix} c_0 \\ c_1 \\ \vdots \\ c_d \end{bmatrix}
=
\begin{bmatrix} y_0 \\ y_1 \\ \vdots \\ y_d \end{bmatrix}
$$

The property of Newton polynomials implies that $A$ is a lower triangular matrix

$$
A =
\begin{bmatrix}
\pi_0(t_0) & 0 & \cdots & 0 \\
\pi_0(t_1) & \pi_1(t_1) & \ddots & \vdots \\
\vdots & \vdots & \ddots & 0 \\
\pi_0(t_d) & \pi_1(t_d) & \cdots & \pi_d(t_d)
\end{bmatrix}
$$

Consequently, the system $A \boldsymbol{c} = \boldsymbol{y}$ is solved by forward substitution.
```

```{div} note
An advantage of the Newton basis is that updating the interpolating function with additional data is simple. In particular, suppose we have data points $(t_0,y_0),\dots,(t_d,y_d)$ and compute the coefficients of the interpolating polynomial

$$
p_d(t)=c_0 + c_1 \pi_1(t)+\dots +c_d \pi_d(t)
$$

Now suppose we get another data point $(t_{d+1},y_{d+1})$. We can find the interpolating polynomial $p_{d+1}(t)$ for the full data set by

$$
p_{d+1}(t)=p_d(t)+c_{d+1} \pi_{d+1}(t)
$$

where the only new coefficient is

$$
c_{d+1} = \frac{y_{d+1}-p_d(t_{d+1})}{\pi_{d+1}(t_{d+1})}.
$$

This observation also leads to an iterative algorithm to compute Newton interpolation in general: given  $(t_0,y_0),\dots,(t_d,y_d)$, find $p_0(t)$ using only $(t_0,y_0)$ and then find $p_{k+1}(t)$ from $p_k(t)$ for each $k$ using the formula above.
```

### Example

```{div} note
Recall that there is a unique polynomial $p(t)$ of degree (at most) $d$ which interpolates $d+1$ points $(t_0,y_0), \dots , (t_d,y_d)$ if the values $t_0,\dots,t_d$ are different. Therefore the monomial, Lagrange and Newton bases all produce the *same* result but computed and represented differently.
```

```{div} example
Find the interpolating polynomial for $(-1,1),(0,0),(1,1)$ using each of the monomial, Lagrange and Newton bases. We know the result is $p(t) = t^2$. Begin with monomial interpolation and setup the Vandermonde matrix and solve the linear system $A \boldsymbol{c} = \boldsymbol{y}$

$$
\left[ \begin{array}{rrr} 1 & -1 & \phantom{+}1 \\ 1 & 0 & 0 \\ 1 & 1 & 1 \end{array} \right]
\begin{bmatrix} c_0 \\ c_1 \\ c_2 \end{bmatrix}
=
\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}
\ \
\Rightarrow
\ \
\boldsymbol{c} = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}
$$

and therefore $c_0 = c_1 = 0$ and $c_2 = 1$ and so $p(t) = t^2$. Now construct the Lagrange basis

$$
\begin{align*}
\ell_0(t) &= \frac{(t - 0)(t - 1)}{(-1 - 0)(-1 - 1)} = \frac{t(t-1)}{2} \\
\ell_1(t) &= \frac{(t - (-1))(t - 1)}{(0 - (-1))(0 - 1)} = 1-t^2 \\
\ell_2(t) &= \frac{(t - (-1))(t - 0)}{(1 - (-1))(1 - 0)} = \frac{t(t+1)}{2}
\end{align*}
$$

and the interpolating polynomial

$$
p(t) = y_0 \ell_0(t) + y_1 \ell_1(t) + y_2 \ell_2(t) = (1)\frac{t(t-1)}{2} + (0) (1-t^2) + (1) \frac{t(t+1)}{2} = t^2
$$

Now construct the Newton basis

$$
\begin{align*}
p_0(t) &= 1 \\
p_1(t) &= t - (-1) = t + 1 \\
p_2(t) &= (t - (-1))(t - 0) = t^2 + t
\end{align*}
$$

Each point yields an equation $p(t_k) = y_k$ for $k=0,1,2$ and so we solve the linear system

$$
\begin{align*}
\left[ \begin{array}{ccc} 1 & 0 & 0 \\ 1 & t_1 - t_0 & 0 \\ 1 & t_2 - t_0 & (t_2-t_0)(t_2-t_1) \end{array} \right]
\begin{bmatrix} c_0 \\ c_1 \\ c_2 \end{bmatrix}
&=
\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}
\\
\left[ \begin{array}{rrr} 1 & 0 & 0 \\ 1 & 1 & 0 \\ 1 & 2 & 2 \end{array} \right]
\begin{bmatrix} c_0 \\ c_1 \\ c_2 \end{bmatrix}
&=
\begin{bmatrix} 1 \\ 0 \\ 1 \end{bmatrix}
\ \
\Rightarrow
\ \
\begin{bmatrix} c_0 \\ c_1 \\ c_2 \end{bmatrix} = \left[ \begin{array}{r} 1 \\ -1 \\ 1 \end{array} \right]
\end{align*}
$$

The interpolating polynomial is again

$$
p(t) = c_0 p_0(t) + c_1 p_1(t) + c_2 p_2(t) = 1 - (t+1) + t^2 + t = t^2
$$
```

## Cubic Spline Interpolation

```{div} bigidea
Given $N+1$ points $(t_0,y_0),\dots,(t_N,y_N)$, a cubic spline is a piecewise cubic polynomial defined by a different polynomial $p_k(t)$ for each subinterval $[t_{k-1},t_k]$, $k=1,\dots,N$.
```

```{div} definition
Consider $N+1$ points $(t_0,y_0),\dots,(t_N,y_N)$. A **cubic spline** is a function $p(t)$ defined piecewise by $N$ cubic polynomials $p_1(t),\dots,p_N(t)$ where

$$
p_k(t) = a_k(t - t_{k-1})^3 + b_k(t - t_{k-1})^2 + c_k(t - t_{k-1}) + d_k \ \ , \ \ t \in [t_{k-1},t_k]
$$

such that $p(t)$, $p'(t)$ and $p''(t)$ are continuous functions.
```

```{div} note
Each polynomial $p_k(t)$ has 4 coefficients $a_k,b_k,c_k,d_k$ therefore we require $4N$ equations to specify the $4N$ unknowns:

1. Interpolation at left endpoints: $p_k(t_{k-1}) = y_{k-1}$ for $k=1,\dots,N$ yields $N$ equations.
2. Interpolation at right endpoints: $p_k(t_k) = y_k$ for $k=1,\dots,N$ yields $N$ equations.
3. Continuity of $p'(t)$: $p_k'(t_k) = p_{k+1}'(t_k)$ for $k=1,\dots,N-1$ yields $N-1$ equations.
4. Continuity of $p''(t)$: $p_k''(t_k) = p_{k+1}''(t_k)$ for $k=1,\dots,N-1$ yields $N-1$ equations.

The conditions impose only $4N-2$ equations therefore we need 2 more to determine the cubic spline uniquely. There are different choices such as the natural spline and the ``not-a-knot" condition.
```

```{div} definition
The **natural cubic spline** requires $p''_1(t_0) = p''_N(t_N) = 0$.
```

```{div} definition
Represent a cubic spline $p(t)$ by the **coefficient matrix**

$$
C=
\begin{bmatrix}
a_1 & a_2 & \cdots & a_N \\
b_1 & b_2 & \cdots & b_N \\
c_1 & c_2 & \cdots & c_N \\
d_1 & d_2 & \cdots & d_N
\end{bmatrix}
$$

where the $k$th column of $C$ consists of the coefficients for the $k$th cubic polynomial in the spline

$$
p_k(t) = a_k(t - t_{k-1})^3 + b_k(t - t_{k-1})^2 + c_k(t - t_{k-1}) + d_k \ \ , \ \ t \in [t_{k-1},t_k]
$$
```

```{div} theorem
Consider $N+1$ points $(t_0,y_0),\dots,(t_N,y_N)$ (with $t_i \not= t_j$ for $i \not= j$). The unique natural cubic spline $p(t)$ which interpolates the points is given by the coefficient matrix

$$
C=
\begin{bmatrix}
a_1 & a_2 & \cdots & a_N \\
b_1 & b_2 & \cdots & b_N \\
c_1 & c_2 & \cdots & c_N \\
d_1 & d_2 & \cdots & d_N
\end{bmatrix}
$$

where $d_k = y_{k-1}$ for $k=1,\dots,N$ and the coefficients $a_1,b_1,c_1,\dots,a_N,b_N,c_N$ are the solution of the linear system

$$
\renewcommand{\arraystretch}{1.5}
\left[ \begin{array}{c|c|c|c|c}
A(L_1) & B & & & \\ \hline
& A(L_2) & B & & \\ \hline
& & \ddots & \ddots & \phantom{A(L_{N-1})} \\ \hline
\phantom{A(L_{N-1})} & \phantom{A(L_{N-1})} & \phantom{A(L_{N-1})} & A(L_{N-1}) & B \\ \hline
T & & & & V
\end{array} \right]
\renewcommand{\arraystretch}{1}
\begin{bmatrix} a_1 \\ b_1 \\ c_1 \\ \vdots \\ a_N \\ b_N \\ c_N \end{bmatrix}
=
\begin{bmatrix} y_1 - y_0 \\ 0 \\ 0 \\ \vdots \\ y_N - y_{N-1} \\ 0 \\ 0 \end{bmatrix}
$$

where $L_k = t_k - t_{k-1}$ is the length of the subinterval $[t_{k-1},t_k]$ and

$$
A(L) = \begin{bmatrix} L^3 & L^2 & L \\ 3L^2 & 2L & 1 \\ 6L & 2 & 0 \end{bmatrix}
\ \
B = \left[ \begin{array}{rrr} 0 & 0 & 0 \\ 0 & 0 & -1 \\ 0 & -2 & 0 \end{array} \right]
\ \
T = \begin{bmatrix} 0 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 0 \end{bmatrix}
\ \
V = \begin{bmatrix} L_N^3 & L_N^2 & L_N \\ 0 & 0 & 0 \\ 6L_N & 2 & 0 \end{bmatrix}
$$
```

```{div} example
Construct the natural cubic spline interpolating points

$$
(0,0),(1,2),(2,1),(3,3),(4,2),(5,4),(6,2),(7,3),(8,1),(9,2),(10,0)
$$

```{image} /img/01_03_02.png
:width: 500px
:align: center
```

```{div} note
The condition number of the matrix for constructing the natural cubic spline does not increase as drastically with the number of points $N+1$ as compared with the Vandermonde matrix. For example, for $11$ equally spaced points $t_0=0,\dots,t_{10}=10$, the Vandermonde matrix is 11 by 11 and has $\mathrm{cond}(A) \approx 10^{12}$ whereas the cubic spline matrix is 30 by 30 and the condition number is only around $33$.

```{image} /img/01_03_03.png
:width: 500px
:align: center
```

## Exercises

**Exercise 1.** Determine whether the statement is **True** or **False**.

  * Consider $N+1$ data points $(t_0,y_0),\dots,(t_N,y_N)$. Let $A_1 \boldsymbol{c}_1 = \boldsymbol{b}_1$ be the linear system such that the solution $\boldsymbol{c}_1$ consists of the coefficients of the interpolating polynomial with respect to the monomial basis. Let $A_2 \boldsymbol{c}_2 = \boldsymbol{b}_2$ be the linear system such that the solution $\boldsymbol{c}_2$ consists of the coefficients of the interpolating natural cubic spline. Then we expect $\mathrm{cond}(A_1) < \mathrm{cond}(A_2)$ for large values of $N$.
  * Consider $d$ data points $(t_1,y_1),\dots,(t_d,y_d)$ (such that $t_i \not= t_j$). There exists a unique polynomial of degree (at most) $d-1$ which interpolates the data.
  * Consider $d$ data points $(t_1,y_1),\dots,(t_d,y_d)$ (such that $t_i \not= t_j$). There exists a unique polynomial $p(t)$ of degree (at most) $d$ which interpolates the data and also satisfies $p'(t_1)=0$ and $p''(t_1)=0$.
  * Consider $d$ data points $(t_1,y_1),\dots,(t_d,y_d)$ (such that $t_i \not= t_j$). There exists a unique polynomial $p(t)$ of degree (at most) $d$ which interpolates the data and $p'(t_1)=0$.
  * Consider $d+1$ data points $(t_0,y_0),(t_1,y_1),\dots,(t_d,y_d)$. Let $p(t)$ be the interpolating polynomial constructed using the monomial basis. Let $q(t)$ be the interpolating polynomial constructed using the Lagrange basis. Let $r(t)$ be the interpolating polynomial constructed using the Newton basis. Then $p(t) = q(t) = r(t)$ for all $t$.

**Exercise 2.** Suppose we have 4 points $(0,y_0),(1,y_1),(2,y_2),(3,y_3)$ and we want to interpolate the data using a spline $p(t)$ constructed from $3$ degree 2 polynomials $p_1,p_2,p_3$ where

$$
p_k(t) = a_k(t - t_{k-1})^2 + b_k(t - t_{k-1}) + c_k \ \ , \ \ t \in [t_{k-1},t_k]
$$

We require that $p(t)$ and $p'(t)$ are continuous and $p''(t_0)=0$. Setup a linear system $A \boldsymbol{x} = \boldsymbol{b}$ where the solution is

$$
\boldsymbol{x} = \begin{bmatrix} a_1 & b_1 & a_2 & b_2 & a_3 & b_3 \end{bmatrix}^T
$$

Note: the system depends on the unspecified values $y_0,y_1,y_2,y_3$.

**Exercise 3.** Suppose a cubic spline $p(t)$ interpolates the data

$$
(0.0,0.0),(1.0,2.0),(2.0,1.0),(3.0,3.0),(4.0,-1.0),(5.0,1.0)
$$

and $p(t)$ has coefficient matrix

$$
\left[ \begin{array}{rrrrr}
1.9 & 1.9 & a_3 & 3.1 & 3.1 \\
-7.2 & -1.5 & b_3 & -6.3 & 3.0 \\
7.3 & -1.4 & c_3 & -0.8 & -4.1 \\
0.0 & 2.0 & 1.0 & 3.0 & -1.0
\end{array}\right]
$$

Determine the coefficients $a_3,b_3,c_3$.

**Exercise 4.** Suppose a cubic spline $p(t)$ interpolates the data

$$
(0.0,1.0),(1.0,3.0),(2.0,1.0),(3.0,1.0),(4.0,2.0),(5.0,1.0)
$$

and $p(t)$ has coefficient matrix (rounded to 2 decimal places)

$$
\left[\begin{array}{rrrrr}
-1.19 & 1.93 & a_3 & -0.75 & 0.55 \\
0.00 & -3.56 & b_3 & 0.60 & -1.65 \\
3.19 & -0.37 & c_3 & 1.15 & 0.10\\
1.00 & 3.00 & 1.00 & 1.00 & 2.00
\end{array} \right]
$$

Determine the coefficients $a_3,b_3,c_3$.
