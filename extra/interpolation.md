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