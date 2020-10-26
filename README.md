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
into your Sage notebook. Thanks to @samorajp for the tip.

## Examples


```python
A = IMatrix([[2, 3, 1], [3,1,0]], separate=1, var=['x', 'y'])
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
