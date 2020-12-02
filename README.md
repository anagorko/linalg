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

## Examples
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
