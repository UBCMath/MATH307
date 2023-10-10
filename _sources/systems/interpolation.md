# Interpolation

```{div} bigidea
An interpolating function provides information about values between points and beyond the range of the data. There are infinitely many different ways to interpolate a set of data! Polynomial interpolation is the simplest method whereas cubic spline interpolation provides much more flexibility.
```

```{image} /img/01_03_01.png
:width: 100%
:align: center
```

## Interpolating Functions

```{div} definition
An **interpolating function** (or **interpolant**) for points $(t_0,y_0),\dots,(t_d,y_d)$ is a function $f(t)$ such that $f(t_k) = y_k$ for $k=0,\dots,d$.
```

```{div} note
There are *infinitely* many different interpolating functions for any set of data. We impose different kinds of constraints on the form of the interpolant to ensure there is a unique solution which is meaningful for the data and can be computed efficiently. In the following sections we describe two kinds of interpolation: *polynomial interpolation* and *cubic spline interpolation*.  
```

## Polynomial Interpolation

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
Given $d+1$ points $(t_0,y_0), \dots,(t_d,y_d)$, **polynomial interpolation** (with respect to the monomial basis) is a polynomial of the form

$$
p(t) = c_0 + c_1 t + \cdots + c_d t^d
$$

such that $p(t_k) = y_k$ for each $k=0,\dots,d$. Each point defines an equation

$$
p(t_k) = y_k \ \ , \ \ \ k=0,\dots,d
$$

therefore the coefficients $c_0,c_1,\dots,c_d$ satisfy the system of linear equations

$$
\begin{align*}
c_0 + c_1t_0 + \cdots + c_d t_0^d &= y_0 \\
c_0 + c_1t_1 + \cdots + c_d t_1^d &= y_1 \\
& \ \vdots \\
c_0 + c_1t_d + \cdots + c_d t_d^d &= y_d
\end{align*}
$$

In matrix notation, the system is $A \boldsymbol{c} = \boldsymbol{y}$ where

$$
A = 
\begin{bmatrix}
1 & t_0 & \cdots & t_0^d \\
1 & t_1 & \cdots & t_1^d \\
\vdots & \vdots & \ddots & \vdots \\
 1 & t_d & \cdots & t_d^d
\end{bmatrix}
\hspace{10mm}
\boldsymbol{c} = \begin{bmatrix} c_0 \\ c_1 \\ \vdots \\ c_d \end{bmatrix}
\hspace{10mm}
\boldsymbol{y} = \begin{bmatrix} y_0 \\ y_1 \\ \vdots \\ y_d \end{bmatrix}
$$

The matrix $A$ is called the **Vandermonde matrix** for $t_0,\dots,t_d$.
```

```{div} theorem
Let $A$ be the Vandermonde matrix for points $(t_0,y_0), \dots,(t_d,y_d)$. Then

$$
\mathrm{det}(A) = \prod_{0 \leq i < j \leq d} (t_j - t_i)
$$

---

*Proof*. See [Wikipedia:Vandermonde matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix) for a proof. Let's compute the formula for the cases $d = 1$ and $d = 2$. Compute for $d = 1$

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
Consider $d+1$ data points $(t_0,y_0), \dots , (t_d,y_d)$ such that $t_i \not= t_j$ for $i \not= j$. There exists a unique polynomial $p(t)$ of degree (at most) $d$ such that $p(t_k) = y_k$ for each $k=0,\dots,d$.

---

*Proof*. The formula for the determinant of the Vandermonde matrix shows that $\mathrm{det}(A) \ne 0$ when the values $t_k$ are distinct therefore it is invertible and so there is a unique solution of the system $A \boldsymbol{c} = \boldsymbol{y}$.
```

````{div} example
Find the unique polynomial of degree 3 which interpolates the points

$$
(0,-1),(1,-1),(2,1),(3,-1)
$$

Solve the system $A \boldsymbol{c} = \boldsymbol{y}$

$$
A =  
\begin{bmatrix}
1 & 0 & 0 &  0 \\
1 & 1 & 1 &  1 \\
1 & 2 & 4 &  8 \\
1 & 3 & 9 & 27
\end{bmatrix}
\hspace{20mm}
\boldsymbol{y} = \left[ \begin{array}{r} -1 \\ -1 \\ 1 \\ -1 \end{array} \right]
$$

Let's compute the LU decomposition of the Vandermonde matrix $A$ so that we can reuse the result in other examples

$$
L =  
\begin{bmatrix}
1 & 0 & 0 & 0 \\
1 & 1 & 0 & 0 \\
1 & 2 & 1 & 0 \\
1 & 3 & 3 & 1
\end{bmatrix}
\hspace{10mm}
U =
\begin{bmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 1 & 1 \\
0 & 0 & 2 & 6 \\
0 & 0 & 0 & 6
\end{bmatrix}
$$

Solve the systems $L \boldsymbol{x} = \boldsymbol{y}$ and $U \boldsymbol{c} = \boldsymbol{x}$ to find

$$
\boldsymbol{c} = \left[ \begin{array}{r} -1 \\ -3 \\ 4 \\ -1 \end{array} \right]
$$

Therefore the interpolating polynomial is

$$
p(t) = -1 - 3t + 4t^2 - t^3
$$

```{image} /img/01_03_02.png
:width: 500px
:align: center
```

````

````{div} note
The condition number of a Vandermonde matrix gets *very* large as the size of the matrix increases. This means that interpolation by the monomial basis is very sensitive to changes in the data for polynomials of large degree. For example, for $11$ equally spaced points $t_0=0,\dots,t_{10}=10$, the Vandermonde matrix $A$ is 11 by 11 and has condition number larger than $10^{12}$. Yikes!

```{image} /img/01_03_03.png
:width: 500px
:align: center
```

````

## Cubic Spline Interpolation

```{div} definition
Consider $N+1$ points $(t_0,y_0),\dots,(t_N,y_N)$. A **cubic spline** is a function $p(t)$ defined piecewise by $N$ cubic polynomials $p_1(t),\dots,p_N(t)$ where

$$
p_k(t) = a_k(t - t_{k-1})^3 + b_k(t - t_{k-1})^2 + c_k(t - t_{k-1}) + d_k \ \ , \ \ t \in [t_{k-1},t_k]
$$

such that $p(t)$, $p'(t)$ and $p''(t)$ are continuous functions.
```

```{div} note
Each polynomial $p_k(t)$ is defined by four coefficients $a_k,b_k,c_k,d_k$ therefore we require $4N$ equations to specify the $4N$ unknowns.

Interpolation at left endpoints yields $N$ equations:

$$
p_k(t_{k-1}) = y_{k-1} \ \ , \ \ \ k=1,\dots,N
$$

Interpolation at right endpoints yields $N$ equations:

$$
p_k(t_k) = y_k \ \ , \ \ \ k=1,\dots,N
$$

Continuity of $p'(t)$ yields $N-1$ equations:

$$
p_k'(t_k) = p_{k+1}'(t_k) \ \ , \ \ \ k=1,\dots,N-1
$$

Continuity of $p''(t)$ yields $N-1$ equations:

$$
p_k''(t_k) = p_{k+1}''(t_k) \ \ , \ \ \ k=1,\dots,N-1
$$

The conditions impose only $4N-2$ equations therefore we need 2 more to determine the cubic spline uniquely. There are different choices such as the natural spline and the "not-a-knot" condition.
```

```{div} definition
A **natural cubic spline** satisfies $p''_1(t_0) = p''_N(t_N) = 0$.
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
$$
$$
T = \begin{bmatrix} 0 & 0 & 0 \\ 0 & 2 & 0 \\ 0 & 0 & 0 \end{bmatrix}
\ \
V = \begin{bmatrix} L_N^3 & L_N^2 & L_N \\ 0 & 0 & 0 \\ 6L_N & 2 & 0 \end{bmatrix}
$$
```

````{div} example
Construct the natural cubic spline interpolating points

$$
(0,0),(1,2),(2,1),(3,3),(4,2),(5,4),(6,2),(7,3),(8,1),(9,2),(10,0)
$$

```{image} /img/01_03_04.png
:width: 500px
:align: center
```

````

````{div} example
Let $p(t)$ be the natural cubic spline which interpolates the data

$$
(0, 1) \ , \ (1, 3) \ , \ (2, 8) \ , \ (3, 10) \ , \ (4, 9) \ , \ (5,-1) \ , \ (6,-17)
$$

Suppose the coefficient matrix of $p(t)$ is 

$$
C = \left[ \begin{array}{rrrrrr}
1 & -2 & 1 & \phantom{+}a_4 & 1 & 1 \\
0 & 3 & -3 & b_4 & -6 & -3 \\
1 & 4 & 4 & c_4 & -5 & -14 \\
1 & 3 & 8 & 10 & 9 & -1
\end{array} \right]
$$
% a4 = -2 , b4 = 0 , c4 = 1
Determine the coefficients $a_4, b_4, c_4$, then compute the value $p''(2.5)$.

First note that the $t$ values are equally spaced thefore $L_k = t_k - t_{k-1} = 1$ for each $k$. Continuity of $p(t)$, $p'(t)$ and $p''(t)$ yields equations

$$
p_k(t_k) = p_{k+1}(t_k) \ \ \Rightarrow \ \ a_k + b_k + c_k + d_k = d_{k+1}
$$

$$
p'_k(t_k) = p'_{k+1}(t_k) \ \ \Rightarrow \ \ 3a_k + 2b_k + c_k = c_{k+1}
$$

$$
p''_k(t_k) = p''_{k+1}(t_k) \ \ \Rightarrow \ \ 6a_k + 2b_k = 2b_{k+1}
$$

Therefore $6a_3 + 2b_3 = 2b_4$ implies $b_4 = 0$, the equation $3a_3 + 2b_3 + c_3 = c_4$ implies $c_4 = 1$ and finally $6a_4 + 2b_4 = 2b_5$ implies $a_4 = -2$.

The value $t=2.5$ lies in the interval $[t_2,t_3]$ therefore we compute

$$
p''(2.5) = p_3''(2.5) = 6a_3(2.5 - 2) + 2b_3 = -3
$$

```{image} /img/01_03_05.png
:width: 500px
:align: center
```

````

````{div} note
The condition number of the matrix for constructing the natural cubic spline does not increase as drastically with the number of points $N+1$ as compared with the Vandermonde matrix. For example, for $11$ equally spaced points $t_0=0,\dots,t_{10}=10$, the Vandermonde matrix is 11 by 11 and has $\mathrm{cond}(A) \approx 10^{12}$ whereas the cubic spline matrix is 30 by 30 and the condition number is only around $33$.

```{image} /img/01_03_06.png
:width: 500px
:align: center
```

````

## Exercises

````{div} exercise
Consider $d$ data points $(t_1,y_1),\dots,(t_d,y_d)$ such that $t_i \not= t_j$ for $i \ne j$. Determine whether the statement is **True** or **False**.

* There exists a unique polynomial of degree (at most) $d-1$ which interpolates the data.
* There exists a unique polynomial $p(t)$ of degree (at most) $d$ which interpolates the data and also satisfies $p'(t_1)=0$ and $p''(t_1)=0$.
* There exists a unique polynomial $p(t)$ of degree (at most) $d$ which interpolates the data and $p'(t_1)=0$.

```{dropdown} Solution
* True
* False
* True
```

````

````{div} exercise
Consider $N+1$ data points $(t_0,y_0),\dots,(t_N,y_N)$. Let $A_1 \boldsymbol{c}_1 = \boldsymbol{b}_1$ such that the solution $\boldsymbol{c}_1$ consists of the coefficients of the interpolating polynomial with respect to the monomial basis. Let $A_2 \boldsymbol{c}_2 = \boldsymbol{b}_2$ such that the solution $\boldsymbol{c}_2$ consists of the coefficients of the interpolating natural cubic spline. Do we expect $\mathrm{cond}(A_1) < \mathrm{cond}(A_2)$ or $\mathrm{cond}(A_1) > \mathrm{cond}(A_2)$ for large values of $N$? Explain.

```{dropdown} Solution
We expect $\mathrm{cond}(A_1) > \mathrm{cond}(A_2)$.
```

````

````{div} exercise
Suppose we have 4 points $(0,y_0),(1,y_1),(2,y_2),(3,y_3)$ and we want to interpolate the data using a spline $p(t)$ constructed from $3$ degree 2 polynomials $p_1,p_2,p_3$ where

$$
p_k(t) = a_k(t - t_{k-1})^2 + b_k(t - t_{k-1}) + c_k \ \ , \ \ t \in [t_{k-1},t_k]
$$

We require that $p(t)$ and $p'(t)$ are continuous and $p''(t_0)=0$. Setup a linear system $A \boldsymbol{x} = \boldsymbol{b}$ where the solution is

$$
\boldsymbol{x} = \begin{bmatrix} a_1 & b_1 & a_2 & b_2 & a_3 & b_3 \end{bmatrix}^T
$$

Note: the system depends on the unspecified values $y_0,y_1,y_2,y_3$.

```{dropdown} Solution
$$
A = \left[ \begin{array}{rrrrrr}
1 & \phantom{+}1 & \phantom{+}0 &  0 & \phantom{+}0 &  0 \\
0 & 0 & 1 &  1 & 0 &  0 \\
0 & 0 & 0 &  0 & 1 &  1 \\
2 & 1 & 0 & -1 & 0 &  0 \\
0 & 0 & 2 &  1 & 0 & -1 \\
1 & 0 & 0 &  0 & 0 &  0
\end{array} \right]
\hspace{10mm}
\boldsymbol{b} = \begin{bmatrix} y_1 - y_0 \\ y_2 - y_1 \\ y_3 - y_2 \\ 0 \\ 0 \\ 0 \end{bmatrix}
$$
```

````

````{div} exercise
Suppose a cubic spline $p(t)$ interpolates the data

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

```{dropdown} Solution
$a_3 = -0.55$, $b_3 = 2.24$, $c_3 = -1.69$
```

````

````{div} exercise
Consider the natural cubic spline $p(t)$ represented by the coefficient matrix

$$
C = \left[ \begin{array}{rrrrrrr}
3 & a_2 & 2 & -4 & -2 & -2 & 8 \\
0 & 9 & b_3 & b_4 & b_5 & -18 & -24 \\
0 & 9 & 12 & 6 & c_5 & -36 & -78 \\
0 & 3 & 16 & 24 & d_5 & 6 & -50
\end{array} \right]
$$     

where $t_0 = 0,t_1 = 1,t_2 = 2,t_3 = 3,t_4 = 4,t_5 = 5,t_6 = 6,t_7 = 7$. Find the value $b_4$ and compute $p''(5.5)$.

```{dropdown} Solution
$b_4 = 0$, $p''(5.5) = -42$
```

````