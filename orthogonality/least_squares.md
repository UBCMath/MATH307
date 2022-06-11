# Least Squares Approximation

```{div} bigidea
Find the best approximation of the system $A \boldsymbol{x} \cong \boldsymbol{b}$ by minimizing the distance $\| A \boldsymbol{x} - \boldsymbol{b}\|$. There are several methods to find the approximation including the normal equations and the QR decomposition.
```

## Normal Equations

```{div} definition
Let $A$ be an $m \times n$ matrix with $m > n$ and $\mathrm{rank}(A) = n$. The best approximation of the system $A \boldsymbol{x} \cong \boldsymbol{b}$ is the vector $\boldsymbol{x}$ which minimizes the distance $\| A\boldsymbol{x} - \boldsymbol{b} \|$. Since $\| \cdot \|$ is the 2-norm, the best approximation is called the **least squares approximation**.
```

```{div} theorem
Let $A$ be an $m \times n$ matrix with $m > n$ and $\mathrm{rank}(A) = n$. The least squares approximation of the system $A \boldsymbol{x} \cong \boldsymbol{b}$ is the solution of the **normal equations**

$$
A^TA\boldsymbol{x} = A^T\boldsymbol{b}
$$

---

*Proof*. If $\boldsymbol{x} \in \mathbb{R}^n$, then $A \boldsymbol{x} \in R(A)$. The projection theorem states that the point in $R(A)$ nearest to $\boldsymbol{b} \in \mathbb{R}^m$ is the orthogonal projection of $\boldsymbol{b}$ onto $R(A)$. If $\boldsymbol{x}$ is the vector such that $A\boldsymbol{x} = \mathrm{proj}_{R(A)}(\boldsymbol{b})$, then $A\boldsymbol{x} - \boldsymbol{b}$ is in $R(A)^{\perp}$ and therefore

$$
A^T(A\boldsymbol{x} - \boldsymbol{b}) = 0 \ \ \Rightarrow \ \ A^TA\boldsymbol{x} = A^T\boldsymbol{b}
$$

We assume $\mathrm{rank}(A) = n$, therefore $A^TA$ is nonsingular and the solution exists and is unique.
```

## QR Equations

```{div} theorem
Let $A$ be an $m \times n$ matrix with $m > n$ and $\mathrm{rank}(A) = n$. The least squares approximation of the system $A \boldsymbol{x} \cong \boldsymbol{b}$ is the solution of the system of equations

$$
R_1\boldsymbol{x} = \boldsymbol{c}_1 \hspace{5mm} \text{where} \hspace{2mm} A = QR = [ Q_1 \ \ Q_2 ] \begin{bmatrix} R_1 \\ 0 \end{bmatrix} \ \ \text{and} \ \ Q^T\boldsymbol{b} = \begin{bmatrix} \boldsymbol{c}_1 \\ \boldsymbol{c}_2 \end{bmatrix}
$$

---

*Proof*. The matrix $Q$ is orthogonal therefore

$$
\| A \boldsymbol{x} - \boldsymbol{b} \|^2 = \| Q(R \boldsymbol{x} - Q^T\boldsymbol{b}) \|^2 = \| R \boldsymbol{x} - Q^T\boldsymbol{b} \|^2
= \left\| \begin{bmatrix} R_1\boldsymbol{x} \\ 0 \end{bmatrix} - \begin{bmatrix} \boldsymbol{c}_1 \\ \boldsymbol{c}_2 \end{bmatrix} \right\|^2
= \| R_1 \boldsymbol{x} - \boldsymbol{c}_1 \|^2 + \| \boldsymbol{c}_2 \|^2
$$

where we use the Pythagoras theorem in the last equality. The vector $\boldsymbol{c}_2$ does not depend on $\boldsymbol{x}$ therefore the minimum value of $\| A \boldsymbol{x} - \boldsymbol{b} \|$ occurs when $R_1 \boldsymbol{x} = \boldsymbol{c}_1$.
```

## Fitting Models to Data

```{div} bigidea
Least squares data fitting computes coefficients $c_1,\dots,c_n$ such that the model function

$$
f(t,\boldsymbol{c}) = c_1 f_1(t) + \cdots + c_n f_n(t)
$$

best fits the data $(t_1,y_1),\dots,(t_m,y_m)$.
```

```{div} definition
Suppose we have $m$ points

$$
(t_1,y_1) , \dots , (t_m,y_m)
$$

and we want to find a line

$$
y=c_1 + c_2t
$$

that "best fits" the data. There are different ways to quantify what "best fit" means but the most common method is called **least squares linear regression**. In least squares linear regression, we want to minimize the sum of squared errors

$$
SSE = \sum_i (y_i - (c_1 + c_2 t_i))^2
$$

In matrix notation, the sum of squared errors is

$$
SSE = \Vert \boldsymbol{y} - A \boldsymbol{c} \Vert^2
$$

where

$$
\boldsymbol{y} = \begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_m \end{bmatrix}
\ \ \
A = \begin{bmatrix} 1 & t_1 \\ 1 & t_2 \\ \vdots & \vdots \\ 1 & t_m \end{bmatrix}
\ \ \
\boldsymbol{c} = \begin{bmatrix} c_1 \\ c_2 \end{bmatrix}
$$

We assume that $m \geq 2$ and $t_i \not= t_j$ for all $i \not= j$ (which implies $\mathrm{rank}(A) = 2$). Therefore the vector of coefficients

$$
\boldsymbol{c} = \begin{bmatrix} c_1 \\ c_2 \end{bmatrix}
$$

is the least squares approximation of the system $A \boldsymbol{c} \cong \boldsymbol{y}$. See [Wikipedia:Simple linear regression](https://en.wikipedia.org/wiki/Simple_linear_regression).
```

```{div} definition
More generally, given $m$ data points

$$
(t_1,y_1) , \dots , (t_m,y_m)
$$

and a **model function** $f(t,\boldsymbol{c})$ which depends on parameters $c_1,\dots,c_n$, the **least squares data fitting problem** consists of computing parameters $c_1,\dots,c_n$ which minimize the sum of squared errors

$$
SSE = \sum_i (y_i - f(t_i,\boldsymbol{c}))^2
$$

If the model function is of the form

$$
f(t,\boldsymbol{c}) = c_1 f_1(t) + \cdots + c_n f_n(t)
$$

for some functions $f_1(t),\dots,f_n(t)$ then we say the data fitting problem is **linear** (but note the function $f_1,\dots,f_n$ are not necessarily linear). In the linear case, use matrix notation to write the sum of squared errors as

$$
SSE = \Vert \boldsymbol{y} - A \boldsymbol{c} \Vert^2
$$

where

$$
\boldsymbol{y} = \begin{bmatrix} y_1 \\ y_2 \\ \vdots \\ y_m \end{bmatrix}
\hspace{10mm}
A = \begin{bmatrix}
f_1(t_1) & f_2(t_1) & \cdots & f_n(t_1) \\
f_1(t_2) & f_2(t_2) & \cdots & f_n(t_2) \\
\vdots & & & \vdots \\
f_1(t_m) & f_2(t_m) & \cdots & f_n(t_m)
\end{bmatrix}
\hspace{10mm}
\boldsymbol{c} = \begin{bmatrix} c_1 \\ c_2 \\ \vdots \\ c_n \end{bmatrix}
$$

We assume that $m \geq n$ and $f_1,\dots,f_n$ are linearly independenty (which implies $\mathrm{rank}(A) = n$). Therefore the vector of coefficients $\boldsymbol{c}$ is the least squares approximation of the system $A \boldsymbol{c} \cong \boldsymbol{y}$.
```

## Exercises

**Exercise 1.** **True** or **False**. Let $A$ be a $m \times n$ matrix with $m \geq n$ and let $\boldsymbol{b} \in \mathbb{R}^m$. There is a unique vector $\boldsymbol{x} \in \mathbb{R}^n$ which minimizes the norm of the residual $\| A \boldsymbol{x} - \boldsymbol{b} \|$.

**Exercise 2.** Let $A = QR$ where

$$
Q = \left[ \begin{array}{rrrrr} 0 & 0 & 0 & 1 & 0 \\ 0 & 1 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 & 0 \\ 1 & 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 & 1 \end{array} \right]
\ \
R = \left[ \begin{array}{rrrr} 1 & 1 & 1 & 1 \\ 0 & 1 & 1 & 1 \\ 0 & 0 & 1 & 1 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 0 & 0 \end{array} \right]
$$

Find the least squares approximation $A\boldsymbol{x} \approx \boldsymbol{b}$ for:

$$
\boldsymbol{b} = \begin{bmatrix} -2 \\ -1 \\ 0 \\ 1 \\ 2 \end{bmatrix}
$$

**Exercise 3.** Setup (but do not solve) a linear system $B \boldsymbol{c} = \boldsymbol{y}$ where the solution is the coefficient vector

$$
\boldsymbol{c} = \begin{bmatrix} c_0 \\ c_1 \\ c_2 \end{bmatrix}
$$

such that the function

$$
f(t) = c_0  + c_1\cos(2 \pi t) + c_2 \sin(2 \pi t)
$$

bests fits the data $(0,1),(1/4,3),(1/2,2),(3/4,-1),(1,0)$.