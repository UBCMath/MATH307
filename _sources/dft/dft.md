# Discrete Fourier Transform

```{div} bigidea
The discrete Fourier transform (DFT) is the orthogonal projection onto the Fourier basis vectors $\boldsymbol{f}_0 , \dots, \boldsymbol{f}_{N-1}$.
```

## Roots of Unity

```{div} definition
An **$N$th root of unity** is a complex number $\omega$ such that $\omega^N = 1$.
```

```{div} proposition
Let $\omega_N = e^{2 \pi i / N}$. Then $\omega_N$ is an $N$th root of unity and $1,\omega_N,\omega_N^2,\dots,\omega_N^{N-1}$ are all the $N$th roots of unity.
```

```{div} proposition
Let $\omega_N = e^{2 \pi i / N}$.

1. $\overline{\omega}_N = \omega_N^{-1} = \omega_N^{N-1}$
2. $\displaystyle \omega_N^k = \cos\left( \frac{2 \pi k}{N} \right) + i \sin \left( \frac{2 \pi k}{N} \right)$
```

```{div} proposition
Let $\omega_N = e^{2 \pi i / N}$ and let $k$ be an integer such that $0<k<N$. Then

$$
\sum_{n=0}^{N-1} \omega_N^{nk} = 0
$$

---

*Proof*. The sum is a geometric series

$$
\sum_{n=0}^{N-1} r^n = \frac{1 - r^N}{1 - r}
$$

and so with $r = \omega_N^k$ we have

$$
\sum_{n=0}^{N-1} \omega_N^{nk} =  \frac{1 - \omega_N^{kN}}{1 - \omega_N^k} = 0
$$

since $\omega_N^{kN} = 1$ and $\omega_N^k \not= 1$.
```

## Fourier Basis

```{div} definition
The **standard basis** of $\mathbb{C}^N$ is $\boldsymbol{e}_0 , \dots, \boldsymbol{e}_{N-1}$ where $\boldsymbol{e}_k$ is the vector with all 0s except 1 in index $k$

$$
\boldsymbol{e}_k = \begin{bmatrix} \boldsymbol{0} \\ 1 \\ \boldsymbol{0} \end{bmatrix} \leftarrow \text{index } k
$$

Use 0-indexing (as in Python) such that the first entry is at index 0. For example, for $N=3$,

$$
\boldsymbol{e}_0 = \begin{bmatrix} 1 \\ 0 \\ 0 \end{bmatrix} \hspace{10mm}
\boldsymbol{e}_1 = \begin{bmatrix} 0 \\ 1 \\ 0 \end{bmatrix} \hspace{10mm}
\boldsymbol{e}_2 = \begin{bmatrix} 0 \\ 0 \\ 1 \end{bmatrix}
$$
```

```{div} definition
Let $N$ be a positive integer and let $\omega_N = e^{2 \pi i / N}$. The **Fourier basis** of $\mathbb{C}^N$ is $\boldsymbol{f}_0 , \dots, \boldsymbol{f}_{N-1}$ where

$$
\renewcommand{\arraystretch}{1.5}
\boldsymbol{f}_k = \begin{bmatrix} 1 \\ \omega^k_N \\ \omega^{2k}_N \\ \vdots \\ \omega^{(N-1)k}_N \end{bmatrix}
\renewcommand{\arraystretch}{1}
$$
```

```{div} example
For $N=2$, $\omega_2 = -1$ and the Fourier basis of $\mathbb{C}^2$ is

$$
\boldsymbol{f}_0 = \left[ \begin{array}{r} 1 \\ 1 \end{array} \right]
\hspace{10mm}
\boldsymbol{f}_1 = \left[ \begin{array}{r} 1 \\ -1 \end{array} \right]
$$

For $N=3$, $\omega_3 = e^{2 \pi i/3} = (-1 + \sqrt{3}i)/2$ and the Fourier basis of $\mathbb{C}^3$ is

$$
\boldsymbol{f}_0 = \left[ \begin{array}{r} 1 \\ 1 \\ 1 \end{array} \right]
\hspace{10mm}
\boldsymbol{f}_1 = \left[ \begin{array}{c} 1 \\ (-1 + \sqrt{3}i)/2 \\ (-1 - \sqrt{3}i)/2 \end{array} \right]
\hspace{10mm}
\boldsymbol{f}_2 = \left[ \begin{array}{c} 1 \\ (-1 - \sqrt{3}i)/2 \\ (-1 + \sqrt{3}i)/2 \end{array} \right]
$$

For $N=4$, $\omega_4 = i$ and the Fourier basis of $\mathbb{C}^4$ is

$$
\boldsymbol{f}_0 = \left[ \begin{array}{r} 1 \\ 1 \\ 1 \\ 1 \end{array} \right]
\hspace{10mm}
\boldsymbol{f}_1 = \left[ \begin{array}{r} 1 \\ i \\ -1 \\ -i \end{array} \right]
\hspace{10mm}
\boldsymbol{f}_2 = \left[ \begin{array}{r} 1 \\ -1 \\ 1 \\ -1 \end{array} \right]
\hspace{10mm}
\boldsymbol{f}_3 = \left[ \begin{array}{r} 1 \\ -i \\ -1 \\ i \end{array} \right]
$$
```

```{div} proposition
The Fourier basis $\boldsymbol{f}_0 , \dots, \boldsymbol{f}_{N-1}$ satisfies

$$
\langle \boldsymbol{f}_k , \boldsymbol{f}_{\ell} \rangle = \left\{ \begin{array}{cl} N & \text{if } k = \ell \\ 0 & \text{otherwise} \end{array} \right.
$$

Therefore the Fourier basis is an orthogonal basis of $\mathbb{C}^N$.

---

*Proof*. Compute

$$
\langle \boldsymbol{f}_k , \boldsymbol{f}_{\ell} \rangle
= \sum_{n=0}^{N-1} \omega_N^{nk} \omega_N^{-n\ell}
= \sum_{n=0}^{N-1} \omega_N^{n(k -\ell)}
$$

We showed in a previous proposition that the sum is equal to 0 if $k \not= \ell$. If $k = \ell$ then clearly the sum is equal to $N$.
```

```{div} proposition
Let $0<k<N$. Then

$$
\overline{\boldsymbol{f}}_k = \boldsymbol{f}_{N-k}
$$

---

*Proof*. By definition and using $\omega_N^N = 1$ we have

$$
\renewcommand{\arraystretch}{1.5}
\overline{\boldsymbol{f}}_k = \begin{bmatrix} 1 \\ \overline{\omega}^k_N \\ \overline{\omega}^{2k}_N \\ \vdots \\ \overline{\omega}^{(N-1)k}_N \end{bmatrix}
= \begin{bmatrix} 1 \\ \omega^{-k}_N \\ \omega^{-2k}_N \\ \vdots \\ \omega^{-(N-1)k}_N \end{bmatrix}
= \begin{bmatrix} 1 \\ \omega^{N-k}_N \\ \omega^{2(N-k)}_N \\ \vdots \\ \omega^{(N-1)(N-k)}_N \end{bmatrix}
= \boldsymbol{f}_{N-k}
\renewcommand{\arraystretch}{1}
$$
```

## Discrete Fourier Transform

```{div} definition
Let $\boldsymbol{x} \in \mathbb{C}^N$. The **discrete Fourier transform** of $\boldsymbol{x}$ is

$$
\mathrm{DFT}(\boldsymbol{x}) = F_N \boldsymbol{x}
$$

where $F_N$ is the **Fourier matrix**

$$
\renewcommand{\arraystretch}{1.5}
F_N =
\begin{bmatrix} & & \overline{\boldsymbol{f}}^T_0 & & \\ & & \overline{\boldsymbol{f}}^T_1 & & \\ & & \vdots & & \\ & & \overline{\boldsymbol{f}}^T_{N-1} & & \end{bmatrix}
=
\renewcommand{\arraystretch}{1.25}
\begin{bmatrix}
1 & 1 & 1 & \cdots & 1 \\
1 & \overline{\omega}_N & \overline{\omega}_N^2 & \cdots & \overline{\omega}_N^{N-1} \\
1 & \overline{\omega}_N^2 & \overline{\omega}_N^4 & \cdots & \overline{\omega}_N^{2(N-1)} \\
1 & \vdots & \vdots & \ddots & \vdots \\
1 & \overline{\omega}_N^{N-1} & \overline{\omega}_N^{2(N-1)} & \cdots & \overline{\omega}_N^{(N-1)^2}
\end{bmatrix}
\renewcommand{\arraystretch}{1}
$$
```

```{div} note
Expand $\boldsymbol{x}$ in terms of the Fourier basis

$$
\boldsymbol{x} = \frac{\langle \boldsymbol{x} , \boldsymbol{f}_0 \rangle}{\langle \boldsymbol{f}_0 , \boldsymbol{f}_0 \rangle} \boldsymbol{f}_0 + \cdots + \frac{\langle \boldsymbol{x} , \boldsymbol{f}_{N-1} \rangle}{\langle \boldsymbol{f}_{N-1} , \boldsymbol{f}_{N-1} \rangle} \boldsymbol{f}_{N-1}
$$

Note that $\langle \boldsymbol{f}_k , \boldsymbol{f}_k \rangle = N$ for each $k=0,\dots,N-1$ and write as matrix multiplication

$$
\boldsymbol{x} = \frac{1}{N}
\begin{bmatrix} & & \\ \boldsymbol{f}_0 & \cdots & \boldsymbol{f}_{N-1} \\ & & \end{bmatrix}
\begin{bmatrix} \langle \boldsymbol{x} , \boldsymbol{f}_0 \rangle \\ \langle \boldsymbol{x} , \boldsymbol{f}_1 \rangle \\ \vdots \\ \langle \boldsymbol{x} , \boldsymbol{f}_{N-1} \rangle \end{bmatrix}
= \frac{1}{N}
\begin{bmatrix} & & \\ \boldsymbol{f}_0 & \cdots & \boldsymbol{f}_{N-1} \\ & & \end{bmatrix}
\renewcommand{\arraystretch}{1.5}
\begin{bmatrix} & & \overline{\boldsymbol{f}}^T_0 & & \\ & & \overline{\boldsymbol{f}}^T_1 & & \\ & & \vdots & & \\ & & \overline{\boldsymbol{f}}^T_{N-1} & & \end{bmatrix}
\renewcommand{\arraystretch}{1}
\boldsymbol{x}
$$

Therefore $DFT(\boldsymbol{x})$ is the vector of coefficients of $\boldsymbol{x}$ with respect to the Fourier basis (up to multiplication by $N$)

$$
DFT(\boldsymbol{x})
=
\begin{bmatrix} \langle \boldsymbol{x} , \boldsymbol{f}_0 \rangle \\ \langle \boldsymbol{x} , \boldsymbol{f}_1 \rangle \\ \vdots \\ \langle \boldsymbol{x} , \boldsymbol{f}_{N-1} \rangle \end{bmatrix}
$$
```

```{div} definition
The DFT is used to study sound, images and any kind of information that can be represented by a vector $\boldsymbol{x} \in \mathbb{C}^N$. Therefore, in the context of the DFT, we use the term **signal** to refer to a (column) vector $\boldsymbol{x} \in \mathbb{C}^N$ and we use the notation

$$
\boldsymbol{x} = \begin{bmatrix} x_0 \\ x_1 \\ x_2 \\ \vdots \\ x_{N-1} \end{bmatrix}
\hspace{10mm}
\boldsymbol{x}[n] = x_n
$$
```

```{div} proposition
Let $\boldsymbol{x}$ be a real signal (that is, $\boldsymbol{x}[k] \in \mathbb{R}$ for each $k=0,\dots,N-1$) and let $\boldsymbol{y} = \mathrm{DFT}(\boldsymbol{x})$. Then

$$
\overline{\boldsymbol{y}[k]} = \boldsymbol{y}[N-k]
$$

---

*Proof*. Compute from the definition

$$
\begin{align*}
\overline{\boldsymbol{y}[k]} &= \overline{\langle \boldsymbol{x} , \boldsymbol{f}_k \rangle} =  \langle \boldsymbol{f}_k , \boldsymbol{x} \rangle = \sum_{n=0}^{N-1} \omega_N^{nk} \overline{x}_n \\
&= \sum_{n=0}^{N-1} \omega_N^{nk-nN} x_n \\
&= \sum_{n=0}^{N-1} \omega_N^{-(N-k)n} x_n = \boldsymbol{y}[N-k]
\end{align*}
$$
```

```{div} proposition
For each $k=0,\dots,N-1$, we have

$$
\mathrm{DFT}(\boldsymbol{f}_k) = N \boldsymbol{e}_k
$$

where $\boldsymbol{e}_k$ is the $k$th standard basis vector.

---

*Proof*. By definition of DFT and the fact $\langle \boldsymbol{f}_k , \boldsymbol{f}_{\ell} \rangle = 0$ when $k \not= \ell$ and $\langle \boldsymbol{f}_k , \boldsymbol{f}_k \rangle = N$, compute

$$
\mathrm{DFT}(\boldsymbol{f}_k) =
\renewcommand{\arraystretch}{1.5}
\begin{bmatrix} & & \overline{\boldsymbol{f}}^T_0 & & \\ & & \vdots & & \\ & & \overline{\boldsymbol{f}}^T_{N-1} & & \end{bmatrix} \boldsymbol{f}_k
=
\begin{bmatrix} \langle \boldsymbol{f}_k , \boldsymbol{f}_0 \rangle \\ \vdots \\ \langle \boldsymbol{f}_k , \boldsymbol{f}_{N-1} \rangle \end{bmatrix}
=
\begin{bmatrix} \boldsymbol{0} \\ N \\ \boldsymbol{0} \end{bmatrix} \leftarrow \text{index } k
\renewcommand{\arraystretch}{1}
$$

and so $\mathrm{DFT}(\boldsymbol{f}_k) = N \boldsymbol{e}_k$.
```

```{div} definition
Let $\boldsymbol{y} \in \mathbb{C}^N$. The **inverse discrete Fourier transform** of $\boldsymbol{y}$ is

$$
\mathrm{IDFT}(\boldsymbol{y}) = \frac{1}{N} \overline{F}^T_N \boldsymbol{y}
$$
```

```{div} note
The Fourier matrix $F_N$ is *not* unitary however the matrix $\frac{1}{\sqrt{N}} F_N$ is unitary.
```

## Exercises

1. Determine whether the statement is **True** or **False**.

   * If $N$ is an even integer then the vector $\boldsymbol{f}_{N/2}$ in the Fourier basis of $\mathbb{C}^N$ has real entries.
   * Let $\boldsymbol{x} \in \mathbb{R}^N$ and let $\boldsymbol{y} = \mathrm{DFT}(\boldsymbol{x})$. Then $\overline{\boldsymbol{x}[k]} = \boldsymbol{x}[N-k]$ for all $0<k<N$.

2. Suppose a signal $\boldsymbol{x} \in \mathbb{R}^9$ of length 9 has real values and let $\boldsymbol{y} = \mathrm{DFT}(\boldsymbol{x})$. Determine all the values of $\boldsymbol{y}$ given the values at even indices

   $$
   \boldsymbol{y}[0] = 1 \hspace{5mm}
   \boldsymbol{y}[2] = 2+i \hspace{5mm}
   \boldsymbol{y}[4] = 1+2i \hspace{5mm}
   \boldsymbol{y}[6] = 1-3i \hspace{5mm}
   \boldsymbol{y}[8] = 1-i
   $$

3. Let $N$ be an even integer and let $\boldsymbol{x} \in \mathbb{R}^N$ such that $\boldsymbol{x}[n] = 1$ if $n$ is even and $\boldsymbol{x}[n] = 0$ if $n$ is odd. Find $\mathrm{DFT}(\boldsymbol{x})$.

4. Let $N$ be an even integer and let $\boldsymbol{x} \in \mathbb{R}^N$ such that $\boldsymbol{x}[n] = 1$ if $n$ is even and $\boldsymbol{x}[n] = -1$ if $n$ is odd. Find $\mathrm{DFT}(\boldsymbol{x})$.

5. Let $N$ be an integer and let $\boldsymbol{x} \in \mathbb{R}^N$ such that $\boldsymbol{x}[0] = 0$ and $\boldsymbol{x}[n] = 1$ for $0<n<N$. Find $\mathrm{DFT}(\boldsymbol{x})$.
