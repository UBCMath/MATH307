# Convolution and Filtering

```{div} bigidea
The discrete Fourier transform of the convolution of two signals is equal to the elementwise product of the discrete Fourier transforms of those signals. In other words, convolution in the time domain corresponds via DFT to elementwise multiplication in the frequency domain.
```

```{div} definition
Let $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{C}^N$. The **convolution** of $\boldsymbol{x}$ and $\boldsymbol{y}$ is the vector $\boldsymbol{x} * \boldsymbol{y} \in \mathbb{C}^N$ given by

$$
( \boldsymbol{x} * \boldsymbol{y} )[n] = \sum_{m=0}^{N-1} \boldsymbol{x}[m] \boldsymbol{y}[n - m]
$$

We interpret $\boldsymbol{y}[n - m]$ using modular arithmetic: if $n - m$ is outside the interval $[0,N-1]$, then we add/subtract multiples of $N$ until the index is inside the interval $[0,N-1]$. See [Wikipedia:Convolution](https://en.wikipedia.org/wiki/Convolution).
```

```{div} example
Let $\boldsymbol{x} = \begin{bmatrix} 1 & 2 & 1 & -1 \end{bmatrix}$ and $\boldsymbol{y} = \begin{bmatrix} 0 & 1/2 & 1/2 & 0 \end{bmatrix}$. Compute

$$
\begin{align*}
( \boldsymbol{x} * \boldsymbol{y} )[0] = \sum_{m=0}^3 \boldsymbol{x}[m] \boldsymbol{y}[- m] &= \boldsymbol{x}[0] \boldsymbol{y}[0] + \boldsymbol{x}[1] \boldsymbol{y}[-1] + \boldsymbol{x}[2] \boldsymbol{y}[-2] + \boldsymbol{x}[3] \boldsymbol{y}[-3] \\
&= \boldsymbol{x}[0] \boldsymbol{y}[0] + \boldsymbol{x}[1] \boldsymbol{y}[3] + \boldsymbol{x}[2] \boldsymbol{y}[2] + \boldsymbol{x}[3] \boldsymbol{y}[1] \\
&= (1)(0) + (2)(0) + (1)(1/2) + (-1)(1/2) \\
&= 0
\end{align*}
$$

$$
\begin{align*}
( \boldsymbol{x} * \boldsymbol{y} )[1] = \sum_{m=0}^3 \boldsymbol{x}[m] \boldsymbol{y}[1 - m] &= \boldsymbol{x}[0] \boldsymbol{y}[1] + \boldsymbol{x}[1] \boldsymbol{y}[0] + \boldsymbol{x}[2] \boldsymbol{y}[-1] + \boldsymbol{x}[3] \boldsymbol{y}[-2] \\
&= \boldsymbol{x}[0] \boldsymbol{y}[1] + \boldsymbol{x}[1] \boldsymbol{y}[0] + \boldsymbol{x}[2] \boldsymbol{y}[3] + \boldsymbol{x}[3] \boldsymbol{y}[2] \\
&= (1)(1/2) + (2)(0) + (1)(0) + (-1)(1/2) \\
&= 0
\end{align*}
$$

$$
\begin{align*}
( \boldsymbol{x} * \boldsymbol{y} )[2] = \sum_{m=0}^3 \boldsymbol{x}[m] \boldsymbol{y}[2 - m] &= \boldsymbol{x}[0] \boldsymbol{y}[2] + \boldsymbol{x}[1] \boldsymbol{y}[1] + \boldsymbol{x}[2] \boldsymbol{y}[0] + \boldsymbol{x}[3] \boldsymbol{y}[-1] \\
&= \boldsymbol{x}[0] \boldsymbol{y}[2] + \boldsymbol{x}[1] \boldsymbol{y}[1] + \boldsymbol{x}[2] \boldsymbol{y}[0] + \boldsymbol{x}[3] \boldsymbol{y}[3] \\
&= (1)(1/2) + (2)(1/2) + (1)(0) + (-1)(0) \\
&= 3/2
\end{align*}
$$

$$
\begin{align*}
( \boldsymbol{x} * \boldsymbol{y} )[3] = \sum_{m=0}^3 \boldsymbol{x}[m] \boldsymbol{y}[3 - m] &= \boldsymbol{x}[0] \boldsymbol{y}[3] + \boldsymbol{x}[1] \boldsymbol{y}[2] + \boldsymbol{x}[2] \boldsymbol{y}[1] + \boldsymbol{x}[3] \boldsymbol{y}[0] \\
&= (1)(0) + (2)(1/2) + (1)(1/2) + (-1)(0) \\
&= 3/2
\end{align*}
$$
```

```{div} definition
Let $\boldsymbol{u}, \boldsymbol{v} \in \mathbb{C}^N$. The **elementwise product** (or **Hadamard product**) of $\boldsymbol{u}$ and $\boldsymbol{v}$ is the vector $\boldsymbol{u} \circ \boldsymbol{v} \in \mathbb{C}^N$ given by

$$
\boldsymbol{u} \circ \boldsymbol{v} =
\begin{bmatrix} u_0 \\ u_1 \\ \vdots \\ u_{N-1} \end{bmatrix}
\circ
\begin{bmatrix} v_0 \\ v_1 \\ \vdots \\ v_{N-1} \end{bmatrix}
=
\begin{bmatrix} u_0v_0 \\ u_1v_1 \\ \vdots \\ u_{N-1}v_{N-1} \end{bmatrix}
$$

See [Wikipedia:Hadamard product](https://en.wikipedia.org/wiki/Hadamard_product_(matrices)).
```

```{div} theorem
Let $\boldsymbol{x}, \boldsymbol{y} \in \mathbb{C}^N$. Then

$$
\mathrm{DFT}( \boldsymbol{x} * \boldsymbol{y} ) = \mathrm{DFT}(\boldsymbol{x}) \circ \mathrm{DFT}(\boldsymbol{y})
$$

where $\circ$ denotes elementwise multiplication. See [Wikipedia:Convolution](https://en.wikipedia.org/wiki/Convolution).

---

*Proof*. Compute from the definitions

$$
\begin{align*}
(\mathrm{DFT}(\boldsymbol{x} * \boldsymbol{y}))[k] &= \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} \boldsymbol{x}[m] \boldsymbol{y}[n - m] \omega_N^{-n k} \\
&= \sum_{n=0}^{N-1} \sum_{m=0}^{N-1} \boldsymbol{x}[m] \boldsymbol{y}[n - m] \omega_N^{-n k} \omega_N^{m k} \omega_N^{-m k} \\
&= \sum_{m=0}^{N-1} \boldsymbol{x}[m] \omega_N^{-m k} \sum_{n=0}^{N-1}  \boldsymbol{y}[n - m] \omega_N^{-(n-m) k} \\
&= \sum_{m=0}^{N-1} \boldsymbol{x}[m] \omega_N^{-m k} \sum_{n=0}^{N-1}  \boldsymbol{y}[n] \omega_N^{-n k} \\
&= \mathrm{DFT}(\boldsymbol{x})[k] \ \mathrm{DFT}(\boldsymbol{y})[k]
\end{align*}
$$

Therefore $\mathrm{DFT}( \boldsymbol{x} * \boldsymbol{y} ) = \mathrm{DFT}(\boldsymbol{x}) \cdot \mathrm{DFT}(\boldsymbol{y})$.
```
