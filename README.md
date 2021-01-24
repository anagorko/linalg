# Interactive Linear Algebra with Sage

This module, meant for **educational purposes only**, supports learning of basics of linear algebra.

It was created to supplement a "Linear algebra and geometry I" course taught during winter semester
at the University of Warsaw Mathematics Department.

## Installation

Just paste and run
```
!wget -N https://raw.githubusercontent.com/anagorko/linalg/main/linalg.py
from linalg import IMatrix
```
in your Sage notebook. Thanks to @samorajp for the tip.

## Usage

### Constructors

A basic matrix notation.

```python
A = IMatrix([[1, 2, 3], [4, 5, 6]], separate=1, names=['x', 'y'])
show(A)
```

<img src="https://render.githubusercontent.com/render/math?math=%5Cleft%5B%5Cbegin%7Barray%7D%7Brr%7Cr%7D1%20%26%202%20%26%203%5C%5C4%20%26%205%20%26%206%5C%5C%5Cend%7Barray%7D%5Cright%5D">

We can render it as a system of linear equations.

```python
show(A.as_equations())
```

<img src="https://render.githubusercontent.com/render/math?math=%5Cleft%5C%7B%5Cbegin%7Barray%7D%7Bccccl%7D%0Ax%26%2B%262%20y%26%3D%263%5C%5C%0A4%20x%26%2B%265%20y%26%3D%266%5C%5C%0A%5Cend%7Barray%7D%5Cright.">

We can render it as a linear combination of column vectors.

```python
show(A.as_combination())
```

<img src="https://render.githubusercontent.com/render/math?math=x%5Cleft%5B%5Cbegin%7Barray%7D%7Bc%7D%0A1%20%5C%5C%0A4%20%5C%5C%0A%5Cend%7Barray%7D%5Cright%5D%2By%5Cleft%5B%5Cbegin%7Barray%7D%7Bc%7D%0A2%20%5C%5C%0A5%20%5C%5C%0A%5Cend%7Barray%7D%5Cright%5D%20%3D%20%5Cleft%5B%5Cbegin%7Barray%7D%7Bc%7D%0A3%20%5C%5C%0A6%20%5C%5C%0A%5Cend%7Barray%7D%5Cright%5D">

A square matrix can be interpreted as a determinant.

We can do symbolic expressions as well. `FractionField` is preferred over `SymbolicRing` because `SR` doesn't work over finite fields (so for example we can't mix parameters and Z_5 in it.)

```python
a, b, c = var('a b c')
F = FractionField(QQ[a, b, c])

B = IMatrix(matrix(F, [[1,a,a^2], [1,b,b^2], [1,c,c^2]]))
show(B)
```
<img src="https://render.githubusercontent.com/render/math?math=%5Cleft%5B%5Cbegin%7Barray%7D%7Brrr%7D%0A1%20%26%20a%20%26%20a%5E%7B2%7D%5C%5C%0A1%20%26%20b%20%26%20b%5E%7B2%7D%5C%5C%0A1%20%26%20c%20%26%20c%5E%7B2%7D%5C%5C%0A%5Cend%7Barray%7D%5Cright%5D">


### Basics

```python
A = IMatrix([[2, 3, 1], [3,1,0]], separate=1, names=['x', 'y'])
show(A.as_equations())
```

```python
A.as_equations().rescale_row(0, -3)
```

```python
A.as_equations().rescale_row(1, 2)
```

```python
A.as_equations().add_multiple_of_row(0, 1, 1)
```

```python
A.as_equations().rescale_row(0, -1/7)
```

```python
A.as_equations().add_multiple_of_row(1, 0, -2)
```

```python
A.as_equations().rescale_row(1, 1/6)
```

```python
A.as_equations().swap_rows(0, 1)
```

### Gaussian elimination

```python
A.to_echelon_form()
```

```python
A.to_reduced_form()
```

### Parameters

```python
t, x1, x2, x3, x4 = var('t x1 x2 x3 x4')
F = PolynomialRing(QQ, [t, x1, x2, x3, x4])

A = IMatrix(matrix(F, [[1, 3, x1], [1, 2, x2], [t, 1, x3], [3, 2, x4]]), separate=1, names=['a_1', 'a_2'])
A.to_echelon_form()
```

## Note

Github doesn't support LaTeX directly so we have to use [a hack](https://gist.github.com/a-rodin/fef3f543412d6e1ec5b6cf55bf197d7b), which negatively affects quality.
