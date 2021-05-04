/*************************************************************************************
 *  Copyright (c) 2021 Stefan Watt
 *  Robert Gordon University, School of Computing
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 3 of the License, or (at
 *  your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 ************************************************************************/

using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace GalaxySim.math
{

    //templates???? not generics _T

    public class matrix<_T>
    {
        protected int m_size;   //  matrix m_size; can be computed at any time, but this should be faster
        protected int m_rows, m_columns;//  number of m_rows & m_columns in matrix
        protected bool m_bClassStorage;//  is (array of m_values) created on heap?
        protected _T[] m_values;

        //  default constructor
        public matrix(int n_mathRows = 16, int n_mathColumns = 16)
        {
            //  create "internal" storage for m_values on heap
            m_values = < _T > (n_mathRows * n_mathColumns);
            m_bClassStorage = true;

            this.m_rows = n_mathRows;
            this.m_columns = n_mathColumns;
            this.m_size = n_mathRows * n_mathColumns;
        }

        private matrix(int n_mathRows, int n_columns, object n_values, bool n_bClassStorage = false)
        {
            //  "external" m_values are used (and modified)
            this.m_values = (_T)n_values;
            this.m_bClassStorage = n_bClassStorage;

            this.m_rows = n_mathRows;
            this.m_columns = n_columns;
            m_size = n_mathRows * n_columns;
        }

        private matrix(matrix<_T> original)
        {
            //  create "internal" storage for m_values on heap
            m_values = < _T > (original.m_rows * original.m_columns);
            m_bClassStorage = true;

            this.m_rows = original.m_rows;
            this.m_columns = original.m_columns;
            this.m_size = m_rows * m_columns;

			// calls Matrix assigment operator (it's the same code)
			operator = (original);
        }

        //  "virtual" to call appropriate destructor when using "late-binding"
        //  or polymorphism
        public virtual void Dispose()
        {
            if (m_bClassStorage)
            {
                m_values;
            }
        }

        //  Assignment operators are not inherited by derived classes.
        //  Performs deep copy. Sometimes called by Matrix constructor.
        //  eg.
        //	Math::Vector<> pos1, pos2;
        //  pos1 = pos2;
        private matrix<_T> CopyFrom(matrix<_T> right_side)
        {
            if (this == right_side)
            {
                return this;
            }
            Assign(this, right_side);
            return this;
        }

        //  eg.
        //	Math::Vector<> pos1, pos2;
        //  pos1 = pos2;		//  operator= internally calls Assign

        private void Assign(matrix<_T> left, matrix<_T> right)
        {
            for (int i = 0; i < right.m_size; i++)
            {
                if ((left).m_size != (right).m_size)
                {
                    throw new System.Exception("These two matrices aren't equal.");
                }
                else
                {
                    left[i] = right[i];
                }
            }

        }

        private virtual _T GetMatrix()
        {
            return m_values;
        }

        private _T GetValue(int column, int row)
        {
            if (((this).m_columns <= column) || ((this).m_rows <= row))
            {
                throw new System.Exception("Index out of matrix bounds.");
            }
            return m_values[row * m_columns + column];
        }

        private void SetValue(int column, int row, _T val)
        {
            if (((this).m_columns <= column) || ((this).m_rows <= row))
            {
                throw new System.Exception("Index out of matrix bounds.");
            }
            m_values[row * m_columns + column] = val;
        }

        private uint GetRows()
        {
            return m_rows;
        }
        private uint GetColumns()
        {
            return m_columns;
        }

        private void LoadIdentity()
        {
            memset(m_values, 0x0, m_size * sizeof(_T));
            for (int i = 0; i < m_size; i += (m_columns + 1))
            {
                m_values[i] = 1;
            }
        }

        //  transpose this matrix
        private void Transpose()
        {
            TransposeInto(this);
        }

        //  transpose this matrix into targetMatrix; this stays the same
        private Matrix<_T> TransposeInto(Matrix<_T> targetMatrix)
        {
            CHECK_SIZES(targetMatrix, this) Matrix<_T>* resultMatrix = targetMatrix;
            if (targetMatrix == this)
            {
                resultMatrix = new Matrix<_T>(m_rows, m_columns);
            }

            //  Main transpose routine
            for (int r = 0; r < m_rows; r++)
            {
                for (int c = 0; c < m_columns; c++)
                {
                    resultMatrix.m_values[c * m_rows + r] = m_values[r * m_columns + c];
                }
            }

            //  If current & target matrices were the same, copy values from
            //  temporary into current matrix.

            if (targetMatrix == this)
            {
                Assign(targetMatrix, *resultMatrix);
                resultMatrix = null;
            }

            return targetMatrix;
        }

        //  compute inversion matrix; inverse this matrix
        private void Inverse()
        {
            InverseInto(this);
        }

        //  inverse this matrix into targetMarix; this matrix stays the same
        private Matrix<_T> InverseInto(Matrix<_T> targetMatrix)
        {
            CHECK_SQUARE(targetMatrix) CHECK_SIZES(this, targetMatrix) Matrix<_T> matrixWithIdentity(targetMatrix.m_rows, 2 * targetMatrix.m_columns);

            //  fill the new matrix with values from current matrix and with identity matrix
            for (int i = 0; i < targetMatrix.m_columns; i++)
            {
                //  i is a row
                for (int j = 0; j < targetMatrix.m_rows; j++)
                {
                    //  j is a column

                    //  fill with current matrix (left side)
                    matrixWithIdentity.m_values[i * matrixWithIdentity.m_columns + j] = m_values[i * targetMatrix.m_columns + j];

                    //  fill with identity matrix (right side)
                    matrixWithIdentity.m_values[i * matrixWithIdentity.m_columns + targetMatrix.m_columns + j] = (i == j ? 1 : 0);
                }
            }

            //  perform GEM
            GEM_PC(matrixWithIdentity);

            //  extract right side of matrix
            for (int i = 0; i < targetMatrix.m_columns; i++)
            {
                //  i is a row
                for (int j = 0; j < targetMatrix.m_rows; j++)
                {
                    //  j is a column
                    targetMatrix.m_values[i * targetMatrix.m_columns + j] = matrixWithIdentity.m_values[i * matrixWithIdentity.m_columns + targetMatrix.m_columns + j];
                }
            }

            return targetMatrix;

        }

        //  Gauss eliminating method (GEM)
        private void GEM(Matrix<_T> targetMatrix)
        {
            for (int i = 0; i < targetMatrix.m_rows; i++)
            {
                //  i is a current phase of GEM otherwise a row, that is taken as
                //  a reference row and isn't modified. This row is called a pivot row.
                for (int j = 0; j < targetMatrix.m_rows; j++)
                {
                    if (i == j)
                    {
                        continue;
                    }
                    //  j is a currently modified row
                    _T coeff = (-1.0) * targetMatrix.m_values[j * targetMatrix.m_columns + i] / targetMatrix.m_values[i * targetMatrix.m_columns + i];

                    for (int k = 0; k < targetMatrix.m_columns; k++)
                    {
                        //  k is a currently modified column in j-row
                        targetMatrix.m_values[j * targetMatrix.m_columns + k] = coeff * targetMatrix.m_values[i * targetMatrix.m_columns + k] + targetMatrix.m_values[j * targetMatrix.m_columns + k];
                    }
                }
            }
        }

        //  Gauss eliminating method with partial search in Column (GEM-PC)
        private void GEM_PC(matrix<_T> targetMatrix)
        {
            for (int i = 0; i < targetMatrix.m_rows; i++)
            {
                //  i is a current phase of GEM otherwise a row, that is taken as
                //  a reference row and isn't modified. This row is called a pivot row.

                //  Look for largest element in column on rows, that does not act
                //  as pivot row yet.
                int swap = i;
                _T val = targetMatrix.m_values[i * targetMatrix.m_columns + i];
                for (int j = i + 1; j < targetMatrix.m_rows; j++)
                {
                    if (Math.Abs(targetMatrix.m_values[j * targetMatrix.m_columns + i]) > Math.Abs(val))
                    {
                        swap = j;
                        val = targetMatrix.m_values[j * targetMatrix.m_columns + i];
                    }
                }

                if (swap != i)
                {
                    //  Swap rows if needed
                    for (int j = 0; j < targetMatrix.m_columns; j++)
                    {
                        /*								std::swap<_T>( targetMatrix.m_values[swap*targetMatrix.m_columns+j],
													targetMatrix.m_values[i*targetMatrix.m_columns+j] );
													*/
                        _T tmp = targetMatrix.m_values[swap * targetMatrix.m_columns + j];
                        targetMatrix.m_values[swap * targetMatrix.m_columns + j] = targetMatrix.m_values[i * targetMatrix.m_columns + j];
                        targetMatrix.m_values[i * targetMatrix.m_columns + j] = tmp;
                    }
                }

                //  Or "== 0" can be sustitute with "<= MACHINE_EPS"
                if (targetMatrix.m_values[i * targetMatrix.m_columns + i] == 0)
                {
                    //  No non-zero pivot. The matrix is singular.
                    throw new System.Exception("Passed singular matrix.");
				}

                //  This assure better accurracy of GEM
                _T temp = targetMatrix.m_values[i * targetMatrix.m_columns + i];
                for (int j = 0; j < targetMatrix.m_columns; j++)
                {
                    targetMatrix.m_values[i * targetMatrix.m_columns + j] /= temp;
                }

                for (int j = 0; j < targetMatrix.m_rows; j++)
                {
                    if (i == j)
                    {
                        continue;
                    }
                    //  j is a currently modified row

                    temp = targetMatrix.m_values[j * targetMatrix.m_columns + i];

                    for (int k = 0; k < targetMatrix.m_columns; k++)
                    {
                        //  k is a currently modified column in j-row
                        targetMatrix.m_values[j * targetMatrix.m_columns + k] -= temp * targetMatrix.m_values[i * targetMatrix.m_columns + k];
                    }
                }
            }
        }
    }


}
