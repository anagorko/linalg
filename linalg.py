"""
This module, meant for **educational purposes only**, supports learning of basics of linear algebra.

It was created to supplement a "Linear algebra and geometry I" course taught during winter semester 2020
at the University of Warsaw Mathematics Department.

To update, run:

!wget -N https://raw.githubusercontent.com/anagorko/linalg/main/linalg.py

AUTHORS:
    Andrzej Nagórko, Jarosław Wiśniewski

VERSION:
    26/10/2020
"""

from typing import Dict, List, Optional, Callable
from itertools import combinations

import sage.all
import sage.structure.sage_object

# noinspection PyUnresolvedReferences
from sage.symbolic.expression import Expression
from sage.misc.html import HtmlFragment
from sage.repl.ipython_kernel.interact import sage_interactive

from ipywidgets import SelectionSlider

import matplotlib.pyplot as plt
import numpy as np

get_ipython().run_line_magic('matplotlib', 'inline')


class IMatrix(sage.structure.sage_object.SageObject):
    """Interactive matrix class."""

    def __init__(self, M, separate=0, copy=True, var: Optional[List[str]] = None):
        """Instantiate IMatrix from a matrix-like object M. Entries of M are copied."""

        if copy:
            self.M = sage.all.matrix(sage.all.SR, M)
        else:
            self.M = M
            """Matrix entries."""
        self.column_separator = separate
        """Column separator placement, counting from the right side."""

        if var is None:
            self.var = [f'x_{i+1}' for i in range(self.M.ncols() - self.column_separator)]
        elif isinstance(var, str):
            self.var = [f'{var}_{i+1}' for i in range(self.M.ncols() - self.column_separator)]
        else:
            assert len(var) == self.M.ncols() - self.column_separator
            self.var = var
        """Variable names"""

        self.var_sr = list()
        for v in self.var:
            self.var_sr.append(sage.all.SR.var(v))
        """Variables."""

    def _repr_(self) -> str:
        return f'IMatrix({repr(list(self.M))}, separate={self.column_separator})'

    def _latex_(self) -> str:
        output = list()

        if self.column_separator > 0:
            column_format = 'r' * (self.M.ncols() - self.column_separator) + '|' + \
                            'r' * self.column_separator
        else:
            column_format = 'r' * self.M.ncols()

        output.append(r'\left[\begin{array}{'f'{column_format}''}')
        for row in self.M:
            output.append(' & '.join([sage.all.latex(el) for el in row]) + r'\\')
        output.append(r'\end{array}\right]')

        return '\n'.join(output)

    @staticmethod
    def _format_coefficient(s: sage.all.SR, for_addition=True) -> str:
        if isinstance(s, Expression):
            if s.args():
                return f'+({sage.all.latex(s)})'

        if s > 0 and for_addition:
            return f'+{sage.all.latex(s)}'

        return f'{sage.all.latex(s)}'

    @staticmethod
    def _format_coefficient_separated(s: sage.all.SR) -> str:
        if isinstance(s, Expression):
            if s.args():
                return f'+ & ({sage.all.latex(s)})'

        if s == 1:
            return '+ & '
        if s >= 0:
            return f'+ & {sage.all.latex(s)}'

        return f'- & {sage.all.latex(-s)}'

    def _format_row_operations(self, op: Dict[int, str]) -> str:
        output = list()

        operations = [r'\ '] * self.M.nrows()
        for i, operation in op.items():
            operations[i] = operation

        output.append(r'\begin{array}{c}')
        output.append(r'\\'.join(operations))
        output.append(r'\end{array}')

        return '\n'.join(output)

    def _add_multiple_of_row(self, i: int, j: int, s: sage.all.SR) -> Dict[int, str]:
        """Add s times row j to row i. The operation is done in place."""

        self.M.add_multiple_of_row(i, j, s)

        return {i: fr'{IMatrix._format_coefficient(s)} \cdot w_{j+1}'}

    def add_multiple_of_row(self, i: int, j: int, s: sage.all.SR) -> HtmlFragment:
        """Add s times row j to row i. The operation is done in place."""

        return HtmlFragment(r'\[' + self._latex_() + self._format_row_operations(self._add_multiple_of_row(i, j, s))
                            + r'\rightarrow' + self._latex_() + r'\]')

    def _swap_rows(self, r1: int, r2: int) -> Dict[int, str]:
        """Swap rows r1 and r2 of self. The operation is done in place."""

        self.M.swap_rows(r1, r2)

        return {r1: fr'\leftarrow w_{r2+1}', r2: fr'\leftarrow w_{r1+1}'}

    def swap_rows(self, r1: int, r2: int) -> HtmlFragment:
        """Swap rows r1 and r2 of self. The operation is done in place."""

        return HtmlFragment(r'\[' + self._latex_() + self._format_row_operations(self._swap_rows(r1, r2))
                            + r'\rightarrow' + self._latex_() + r'\]')

    def _rescale_row(self, i: int, s: sage.all.SR) -> Dict[int, str]:
        """Replace i-th row of self by s times i-th row of self. The operation is done in place."""

        self.M.rescale_row(i, s)

        return {i: fr'\cdot {IMatrix._format_coefficient(s, False)}'}

    def rescale_row(self, i: int, s: sage.all.SR) -> HtmlFragment:
        """Replace i-th row of self by s times i-th row of self. The operation is done in place."""

        return HtmlFragment(r'\[' + self._latex_() + self._format_row_operations(self._rescale_row(i, s))
                            + r'\rightarrow' + self._latex_() + r'\]')

    def to_echelon_form(self) -> HtmlFragment:
        """Transform self to the echelon form of self."""

        output = list()

        # Gaussian elimination algorithm is a modified version of
        # https://ask.sagemath.org/question/8840/how-to-show-the-steps-of-gauss-method/

        col = 0  # all cols before this are already done
        for row in range(0, self.M.nrows()):
            # Need to swap in a nonzero entry from below
            while col < self.M.ncols() and self.M[row][col] == 0:
                for i in self.M.nonzero_positions_in_column(col):
                    if i > row:
                        output.append(r'\[')
                        output.append(self._latex_())
                        output.append(self._format_row_operations(self._swap_rows(row, i)))
                        output.append(r'\rightarrow')
                        output.append(self._latex_())
                        output.append(r'\]')
                        break
                else:
                    col += 1

            if col >= self.M.ncols() - self.column_separator:
                break

            # Now guaranteed M[row][col] != 0
            if self.M[row][col] != 1:
                if self.M[row][col].args():
                    output.append(f'<br>Przerywam eliminację bo nie wiem, czy wyrażenie'
                                  f'${sage.all.latex(self.M[row][col])}$ jest niezerowe.')
                    break
                else:
                    output.append(r'\[')
                    output.append(self._latex_())
                    output.append(self._format_row_operations(self._rescale_row(row, 1 / self.M[row][col])))
                    output.append(r'\rightarrow')
                    output.append(self._latex_())
                    output.append(r'\]')

            change_flag = False
            unchanged = self._latex_()
            operations = dict()
            for changed_row in range(row + 1, self.M.nrows()):
                if self.M[changed_row][col] != 0:
                    change_flag = True
                    factor = -1 * self.M[changed_row][col] / self.M[row][col]
                    operations.update(self._add_multiple_of_row(changed_row, row, factor))

            if change_flag:
                output.append(r'\[')
                output.append(unchanged)
                output.append(self._format_row_operations(operations))
                output.append(r'\rightarrow')
                output.append(self._latex_())
                output.append(r'\]')

            col += 1

        return HtmlFragment('\n'.join(output))

    def to_reduced_form(self):
        """Transform self to the reduced echelon form of self."""
        output = list()

        for changed_row in range(1, self.M.nrows()):
            operations = {i: r'\ ' for i in range(self.M.nrows())}
            unchanged = self._latex_()

            for row in range(0, changed_row):
                factor = -self.M[row][changed_row]
                operations.update(self._add_multiple_of_row(row, changed_row, factor))

            output.append(r'\[')
            output.append(unchanged)
            output.append(self._format_row_operations(operations))
            output.append(r'\rightarrow')
            output.append(self._latex_())
            output.append(r'\]')

        return HtmlFragment('\n'.join(output))

    def as_equations(self) -> 'SoLE':
        return SoLE(self)

    def as_combination(self) -> 'LinearCombination':
        return LinearCombination(self)

    def plot(self):
        return self.as_equations().plot()


class SoLE(IMatrix):
    """System of linear equations."""

    def __init__(self, M: IMatrix):
        super().__init__(M.M, separate=M.column_separator, copy=False, var=M.var)

        if self.column_separator != 1:
            print('Uwaga: macierz nie wygląda na układ równań.')

        self.x_min = None
        self.x_max = None
        self.y_min = None
        self.y_max = None
        """Plot range"""

    def _format_row_operations(self, op: Dict[int, str]) -> str:
        output = list()

        operations = [r'\ '] * self.M.nrows()
        for i, operation in op.items():
            operations[i] = '/' + operation

        output.append(r'\begin{array}{c}')
        output.append(r'\\'.join(operations))
        output.append(r'\end{array}')

        return '\n'.join(output)

    def _latex_(self) -> str:
        output = list()

        column_format = 'c' * (self.M.ncols() - self.column_separator) * 2 + 'l' * self.column_separator

        output.append(r'\left\{\begin{array}{'f'{column_format}''}')
        output += ['&'.join(self._format_row(row)) + r'\\' for row in self.M]
        output.append(r'\end{array}\right.')

        return '\n'.join(output)

    def _format_row(self, row: List) -> List[str]:
        """A latex representation of a row."""

        empty_so_far = True
        row_output = list()
        for i, coefficient in enumerate(row):
            sign = ''
            if i > 0 and not empty_so_far:
                if i == self.M.ncols() - self.column_separator:
                    sign = '='
                elif i < self.M.ncols() - self.column_separator:
                    if coefficient == 0.0:
                        sign = ''
                    elif coefficient < 0.0:
                        sign = '-'
                        coefficient = -coefficient
                    else:
                        sign = '+'

            variable = 1
            if coefficient != 0.0 and i < self.M.ncols() - self.column_separator:
                variable = self.var_sr[i]
            term = variable * coefficient

            if term == 0.0 and i < self.M.ncols() - self.column_separator:
                term = ''

            if i == 0:
                row_output.append(sign + sage.all.latex(term))
            else:
                row_output.append(sign)
                row_output.append(sage.all.latex(term))

            if coefficient != 0.0:
                empty_so_far = False

        return row_output

    """
    Methods related to plot() function.
    """

    MAX_ASPECT_RATIO = 2.0
    SLIDER_RANGE = [sage.all.QQ(i - 125) / 50 for i in range(251)]

    @staticmethod
    def _row_to_function(coefficients: List) -> Callable:
        """Converts a list of coefficients into a linear function."""

        def f(*args):
            return (1 / coefficients[-2]) * (coefficients[-1] - sum(v * coefficients[i] for i, v in enumerate(args)))

        return f

    def _equation(self, coefficients: List) -> Expression:
        return sum(self.var_sr[i] * coefficients[i] for i in range(len(coefficients) - 1)) == coefficients[-1]
#        return self.var_sr[-1] == coefficients[-1] / coefficients[-2] - \
#               sum(self.var_sr[i] * coefficients[i] / coefficients[-2] for i in range(len(coefficients) - 2))

    def _format_equation(self, row: List) -> str:
        return f'${"".join(self._format_row(row))}$'

    def plot(self):
        """Interactive plot of self. Supported in two and three dimensional cases."""

        free_variables = list(sum(sum(self.M)).free_variables())

        var_sliders = {str(var): SelectionSlider(options=SoLE.SLIDER_RANGE, continuous_update=False, value=0)
                       for var in free_variables}

        def f(**_var_sliders):
            M = self.M.subs({sage.all.SR.var(v): _var_sliders[v] for v in _var_sliders})

            avg_y = 0
            """Average y-value of subsystems consisting of two equations."""

            """Determine x- and y- range of the plot."""
            for i, j in combinations(range(M.nrows()), 2):
                solution = sage.all.solve([self._equation(M[i]), self._equation(M[j])],
                                          *self.var_sr, solution_dict=True)
                if solution:
                    solution = solution[0]
                    if solution[self.var_sr[0]].free_variables() or solution[self.var_sr[1]].free_variables():
                        continue

                    avg_y += solution[self.var_sr[1]]

                    if self.x_min is None or solution[self.var_sr[0]] - 1 < self.x_min:
                        self.x_min = solution[self.var_sr[0]] - 1

                    if self.x_max is None or solution[self.var_sr[0]] + 1 > self.x_max:
                        self.x_max = solution[self.var_sr[0]] + 1

            if self.x_min is None:
                self.x_min = -10.0

            if self.x_max is None:
                self.x_max = 10.0

            avg_y = avg_y / M.nrows()
            self.y_max = avg_y + SoLE.MAX_ASPECT_RATIO * (self.x_max - self.x_min) / 2
            self.y_min = avg_y - SoLE.MAX_ASPECT_RATIO * (self.x_max - self.x_min) / 2

            x = np.arange(self.x_min, self.x_max, sage.all.QQ((self.x_max - self.x_min)/100))

            fig, ax = plt.subplots()

            for i in range(M.nrows()):
                if M[i][-2] == 0.0:
                    # vertical line
                    if M[i][0] != 0.0:
                        xv = M[i][-1] / M[i][0]
                        ax.plot([xv, xv], [self.y_min, self.y_max], label=self._format_equation(M[i]))
                else:
                    row_f = SoLE._row_to_function(M[i])
                    y = [row_f(_x) for _x in x]

                    x_clip, y_clip = list(), list()

                    for _x, _y in zip(x, y):
                        if self.y_min <= _y <= self.y_max:
                            x_clip.append(_x)
                            y_clip.append(_y)

                    ax.plot(x_clip, y_clip, label=self._format_equation(M[i]))

            ax.set(xlabel=self.var[0], ylabel=self.var[1])
            ax.grid()

            plt.legend()
            plt.show()

        w = sage_interactive(f, **var_sliders)
        # output = w.children[-1]
        # output.layout.height = '600px'

        default_dpi = plt.rcParamsDefault['figure.dpi']
        plt.rcParams['figure.dpi'] = default_dpi * 1.5

        display(w)

    def _repr_(self) -> str:
        return f'IMatrix({repr(list(self.M))}, separate={self.column_separator}).as_equations()'


class LinearCombination(IMatrix):
    """System of linear equations interpreted as linear combination."""

    def __init__(self, M: IMatrix):
        super().__init__(M.M, separate=M.column_separator, copy=False, var=M.var)

        if self.column_separator != 1:
            print('Uwaga: macierz nie wygląda na układ równań.')

    def _format_column(self, col_n: int) -> str:
        """Format column as a column vector."""

        output = list()

        output.append(r'\left[\begin{array}{c}')
        output += [sage.all.latex(self.M[i][col_n]) + r'\\' for i in range(self.M.ncols())]
        output.append(r'\end{array}\right]')

        return '\n'.join(output)

    def _latex_(self) -> str:
        lhs = list()
        for i in range(self.M.ncols() - self.column_separator):
            lhs.append(self.var[i] + self._format_column(i))

        output = ['+'.join(lhs), '=', self._format_column(-1)]

        return ' '.join(output)


def main():
    """Basic functionality check, ran when script is invoked from the command line."""

    M = IMatrix([[1, 2, 3], [4, 5, 6]])
    print(M.add_multiple_of_row(0, 1, 3))


if __name__ == '__main__':
    main()
else:
    print(__doc__)
    from IPython.core.display import HTML
    display(HTML("<style>.container { width:99% !important; }</style>"))
