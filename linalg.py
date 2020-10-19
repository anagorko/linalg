"""
This module, meant for **educational purposes only**, supports learning of basics of linear algebra.

It was created to supplement a "Linear algebra and geometry I" course taught during winter semester 2020
at the University of Warsaw Mathematics Department.

AUTHORS:
    Andrzej Nagórko, Jarosław Wiśniewski
"""

from typing import Dict, List, Optional, Union

import sage.all
import sage.structure.sage_object

# noinspection PyUnresolvedReferences
from sage.symbolic.expression import Expression
from sage.misc.html import HtmlFragment


class IMatrix(sage.structure.sage_object.SageObject):
    """Interactive matrix class."""

    def __init__(self, M, separate=0, copy=True, var: Optional[Union[str, List[str]]] = None):
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

    def _repr_(self) -> str:
        return repr(self.M)

    def _latex_(self) -> str:
        output = list()

        if self.column_separator > 0:
            column_format = 'r' * (self.M.ncols() - self.column_separator) + '|' + \
                            'r' * self.column_separator
        else:
            column_format = 'r' * self.M.ncols()

        output.append('\\left[\\begin{array}{'f'{column_format}''}')
        for row in self.M:
            output.append(' & '.join([sage.all.latex(el) for el in row]) + '\\\\')
        output.append('\\end{array}\\right]')

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

        operations = ['\\ '] * self.M.nrows()
        for i, operation in op.items():
            operations[i] = operation

        output.append('\\begin{array}{c}')
        output.append('\\\\'.join(operations))
        output.append('\\end{array}')

        return '\n'.join(output)

    def _add_multiple_of_row(self, i: int, j: int, s: sage.all.SR) -> Dict[int, str]:
        """Add s times row j to row i. The operation is done in place."""

        self.M.add_multiple_of_row(i, j, s)

        return {i: f'{IMatrix._format_coefficient(s)} \\cdot w_{j+1}'}

    def add_multiple_of_row(self, i: int, j: int, s: sage.all.SR) -> HtmlFragment:
        """Add s times row j to row i. The operation is done in place."""

        return HtmlFragment('\\[' + self._latex_() + self._format_row_operations(self._add_multiple_of_row(i, j, s))
                            + '\\rightarrow' + self._latex_() + '\\]')

    def _swap_rows(self, r1: int, r2: int) -> Dict[int, str]:
        """Swap rows r1 and r2 of self. The operation is done in place."""

        self.M.swap_rows(r1, r2)

        return {r1: f'\\leftarrow w_{r2+1}', r2: f'\\leftarrow w_{r1+1}'}

    def swap_rows(self, r1: int, r2: int) -> HtmlFragment:
        """Swap rows r1 and r2 of self. The operation is done in place."""

        return HtmlFragment('\\[' + self._latex_() + self._format_row_operations(self._swap_rows(r1, r2))
                            + '\\rightarrow' + self._latex_() + '\\]')

    def _rescale_row(self, i: int, s: sage.all.SR) -> Dict[int, str]:
        """Replace i-th row of self by s times i-th row of self. The operation is done in place."""

        self.M.rescale_row(i, s)

        return {i: f'\\cdot {IMatrix._format_coefficient(s, False)}'}

    def rescale_row(self, i: int, s: sage.all.SR) -> HtmlFragment:
        """Replace i-th row of self by s times i-th row of self. The operation is done in place."""

        return HtmlFragment('\\[' + self._latex_() + self._format_row_operations(self._rescale_row(i, s))
                            + '\\rightarrow' + self._latex_() + '\\]')

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
                        output.append('\\[')
                        output.append(self._latex_())
                        output.append(self._format_row_operations(self._swap_rows(row, i)))
                        output.append('\\rightarrow')
                        output.append(self._latex_())
                        output.append('\\]')
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
                    output.append('\\[')
                    output.append(self._latex_())
                    output.append(self._format_row_operations(self._rescale_row(row, 1 / self.M[row][col])))
                    output.append('\\rightarrow')
                    output.append(self._latex_())
                    output.append('\\]')

            change_flag = False
            unchanged = self._latex_()
            operations = dict()
            for changed_row in range(row + 1, self.M.nrows()):
                if self.M[changed_row][col] != 0:
                    change_flag = True
                    factor = -1 * self.M[changed_row][col] / self.M[row][col]
                    operations.update(self._add_multiple_of_row(changed_row, row, factor))

            if change_flag:
                output.append('\\[')
                output.append(unchanged)
                output.append(self._format_row_operations(operations))
                output.append('\\rightarrow')
                output.append(self._latex_())
                output.append('\\]')

            col += 1

        return HtmlFragment('\n'.join(output))

    def to_reduced_form(self):
        """Transform self to the reduced echelon form of self."""
        output = list()

        for changed_row in range(1, self.M.nrows()):
            operations = {i: '\\ ' for i in range(self.M.nrows())}
            unchanged = self._latex_()

            for row in range(0, changed_row):
                factor = -self.M[row][changed_row]
                operations.update(self._add_multiple_of_row(row, changed_row, factor))

            output.append('\\[')
            output.append(unchanged)
            output.append(self._format_row_operations(operations))
            output.append('\\rightarrow')
            output.append(self._latex_())
            output.append('\\]')

        return HtmlFragment('\n'.join(output))

    def as_equations(self) -> 'SoLE':
        return SoLE(self)


class SoLE(IMatrix):
    """System of linear equations."""

    def __init__(self, M: IMatrix):
        super().__init__(M.M, separate=M.column_separator, copy=False)

        if self.column_separator != 1:
            print('Uwaga: macierz nie wygląda na układ równań.')

        self.var = M.var
        """Variable names"""

    def _format_row_operations(self, op: Dict[int, str]) -> str:
        output = list()

        operations = ['\\ '] * self.M.nrows()
        for i, operation in op.items():
            operations[i] = '/' + operation

        output.append('\\begin{array}{c}')
        output.append('\\\\'.join(operations))
        output.append('\\end{array}')

        return '\n'.join(output)

    def _latex_(self) -> str:
        output = list()

        column_format = 'c' * (self.M.ncols() - self.column_separator) * 2 + 'l' * self.column_separator

        output.append('\\left\\{\\begin{array}{'f'{column_format}''}')
        for row in self.M:
            empty_so_far = True
            row_output = ''
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
                        else:
                            sign = '+'

                if isinstance(coefficient, Expression) and coefficient.args():
                    coeff = f'({sage.all.latex(coefficient)})'
                else:
                    if coefficient < 0:
                        if empty_so_far:
                            if coefficient == -1.0:
                                coeff = '-'
                            else:
                                coeff = sage.all.latex(coefficient)
                        elif i >= self.M.ncols() - self.column_separator:
                            coeff = sage.all.latex(coefficient)
                        else:
                            if coefficient == -1.0:
                                coeff = ''
                            else:
                                coeff = sage.all.latex(-coefficient)
                    elif coefficient in {0.0, 1.0} and i < self.M.ncols() - self.column_separator:
                        coeff = ''
                    else:
                        coeff = sage.all.latex(coefficient)

                variable = ''
                if coefficient != 0.0 and i < self.M.ncols() - self.column_separator:
                    variable = self.var[i]

                if i == 0:
                    row_output += sign + coeff + ' ' + variable
                else:
                    row_output += '&' + sign + '&' + coeff + ' ' + variable

                if coefficient != 0.0:
                    empty_so_far = False

            output.append(row_output + '\\\\')
        output.append('\\end{array}\\right.')

        return '\n'.join(output)


def main():
    """Basic functionality check, ran when script is invoked from the command line."""

    M = IMatrix([[1, 2, 3], [4, 5, 6]])
    print(M.add_multiple_of_row(0, 1, 3))


if __name__ == '__main__':
    main()
else:
    print(__doc__)
