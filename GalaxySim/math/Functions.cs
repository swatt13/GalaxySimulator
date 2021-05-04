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

namespace GalaxySim.math
{
    public delegate double Function(double func);
    public delegate double Integral(double integral);
    class Functions
    {
        public static double Sin_angle(double angle)
        {
            double res;
            if ((angle == 360) || (angle == 180))
            {
                res = 0;
            }
            else
            {
                res = Math.Sin(angle * Math.PI / 180.0);
            }
            return res;
        }

        public static double Cos_angle(double angle)
        {
            double res;
            if ((angle == 90) || (angle == 270) || (angle == -90) || (angle == -270))
            {
                res = 0;
            }
            else
            {
                res = Math.Cos(angle * Math.PI / 180.0);
            }
            return res;
        }

        public static double Tan_angle(double angle)
        {
            return Math.Tan(angle + Math.PI / 180.0);
        }

        public static double Asin_angle(double x)
        {
            double res;
            if (Math.Ceiling(x) == -1.0) res = -90.0;
            else
            {
                if (Math.Floor(x) == 1.0) res = 90.0;
                else
                    res = Math.Asin(x) * 180.0 / Math.PI;
            }
            return res;
        }

        public static double Acos_angle(double x)
        {
            double res;
            if (Math.Floor(x) == 1.0) res = 0;
            else
            {
                if (Math.Ceiling(x) == -1.0) res = 100.0;
                else
                    res = Math.Acos(x) * 180.0 / Math.PI;
            }
            return res;
        }

        public static double Atan_angle(double x)
        {
            return Math.Atan(x) * 180.0 / Math.PI;
        }

        public static double Atan_angle2(double y, double x)
        {
            return Math.Atan2(y, x) * 180.0 / Math.PI;
        }

        public static double Random(double idum)
        {
            Random random = new Random(1047239);
            ///////////////////////////
            return random.NextDouble();
        }

        public static double TrapezoidalIntegration(Integral func, double lower, double upper, int steps)
        {
            double integral = 0;
            double h = (upper - lower) / steps;
            double xi = lower;

            for (int k = 0; k < steps; k++)
            {
                integral = integral + h * (func(xi) + func(xi + h)) / 2;
                if ((integral == integral) == false)        //  integral is NaN
                {
                    throw new Exception();
                }

                xi = xi + h;
            }

            return integral;
        }

        /* Returns the integral of the function func from a to b. The parameters EPS can be set to the
        desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed
        number of steps. Integration is performed by Simpson's rule.
        */
        public static double SimpsonIntegration(Integral f, double lower,
            double upper, int steps)
        {
            double EPS = 1.0e-6;
            double JMAX = 20;

            int j;
            double s, st, ost, os;
            ost = os = -1.0e30;
            for (j = 1; j <= JMAX; j++)
            {
                st = TrapezoidalIntegration(f, lower, upper, j);
                s = (4.0 * st - ost) / 3.0; //  Compare equation (4.2.4), above.
                if (j > 5) //  Avoid spurious early convergence.
                    if (Math.Abs(s - os) < EPS * Math.Abs(os) || (s == 0.0 && os == 0.0)) return s;
                os = s;
                ost = st;
            }

            throw new Exception();
            return 0.0; ///  Never get here.
        }
        public delegate double Point(double f);
        public static double FindMaximum(Point f, double L, double U, double tol)
        {

            double t_L, t_R;

            t_L = L + 0.382 * (U - L);
            t_R = L + 0.618 * (U - L);

            while ((U - L) > tol)
            {

                if (f(t_L) > f(t_R))
                {

                    U = t_R;
                    t_R = t_L;
                    t_L = L + 0.382 * (U - L);

                }
                else
                {

                    L = t_L;
                    t_L = t_R;
                    t_R = L + 0.618 * (U - L);
                }
            }

            return (U + L) / 2.0;
        }

        public static double FindMinimum(Point f, double L, double U, double tol)
        {

            double t_L, t_R;

            t_L = L + 0.382 * (U - L);
            t_R = L + 0.618 * (U - L);

            while ((U - L) > tol)
            {

                if (f(t_L) < f(t_R))
                {

                    U = t_R;
                    t_R = t_L;
                    t_L = L + 0.382 * (U - L);

                }
                else
                {

                    L = t_L;
                    t_L = t_R;
                    t_R = L + 0.618 * (U - L);
                }
            }

            return (U + L) / 2.0;
        }

        public static double semigausdistr(double x)
        {
            // y = exp(- x^2)
            return Math.Exp(-Math.Pow(x, 2));
        }

        public static double erf(double q)
        {
            return ((2 / Math.Sqrt(Math.PI)) * TrapezoidalIntegration(semigausdistr, 0, q, 5000));
        }

        //  Lagrange interpolation
        //  Appeared in astronomical computing, Sky & Telescope, April, 1984
        public static double? LagrangeInterpolation(int n, double[] x, double[] f, double val)
        {
            double[] L = new double[n];
            for (int i = 0; i < n; i++)
            {
                L[i] = 1;

                for (int j = 0; j < n; j++)
                {
                    if (j == i) continue;
                    L[i] = L[i] * (x[i] - x[j]);

                }

                L[i] = f[i] / L[i];
            }

            double F1 = 0;
            double? F = null;
            for (int i = 0; i < n; i++)
            {
                if (val != x[i]) continue;
                F = f[i];
                F1 = 1;
            }

            if (F1 != 1)
            {
                double T = 1;
                F = 0;

                for (int i = 0; i < n; i++)
                {
                    T = T * (val - x[i]);
                }

                for (int i = 0; i < n; i++)
                {
                    F = F + L[i] * T / (val - x[i]);
                }
            }
            if (F.Equals(null))
                throw new Exception();
            else { return F; }
        }

        //  Interpolate a set of N points by fitting a polynomial of degree N-1
        //  Adapted from algorithm in Numerical Recipes, Press et al. (1992), 
        //  Section 3.1. 

        public static void PolynomialInterpolation(double[] xa, double[] ya, int n, double x, double y/*out*/, double dy/*out*/)
        {
            int i, m, ns = 1;
            double den, dift, ho, hp, w;

            double[] c = new double[n + 1];  //  bacause C recipes arrays are from 1 not 0
            double[] d = new double[n + 1];

            double dif = Math.Abs(x - xa[1]);

            for (i = 1; i <= n; i++)
            {
                //  find the index of xa which is closest to x
                if ((dift = Math.Abs(x - xa[i])) < dif)
                {
                    ns = i;
                    dif = dift;
                }

                //  and initialize the tableau of c's and d's
                c[i] = ya[i];
                d[i] = ya[i];
            }

            y = ya[ns--]; //  this is the initial approximation to y

            for (m = 1; m < n; m++)
            {
                //  for each column of the tableau
                for (i = 1; i <= n - m; i++)
                {
                    //  we loop over the current c's and d's and update them
                    ho = xa[i] - x;
                    hp = xa[i + m] - x;
                    w = c[i + 1] - d[i];
                    //  This error can occur only if two input xa's are (to within roundof) identical
                    if ((den = ho - hp) == 0.0) throw new Exception();  //std::runtime_error("Error in routine polint");
                    den = w / den;
                    //  Here the c's and d's are updated.
                    d[i] = hp * den;
                    c[i] = ho * den;
                }

                y += (dy = (2 * ns < (n - m) ? c[ns + 1] : d[ns--]));
            }
        }
    }
    class AdditiveIntegral
    {
        protected double m_integral;
        protected double m_xi;
        public AdditiveIntegral()
        {
            m_integral = 0;
            m_xi = 0;
        }
        public double Integrate(Integral func, double upper, int steps)
        {
            double h = (upper - m_xi) / (double)steps;

            for (int k = 0; k < steps; k++)
            {
                m_integral += h * (func(m_xi) + func(m_xi + h) )/ 2;
                if (double.IsNaN(m_integral))
                {
                    throw new Exception();
                }
                m_xi = m_xi + h;
            }

            m_xi = m_xi - h;
            return m_integral;
        }
    }
}
