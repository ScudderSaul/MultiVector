using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace MvLib
{
    public class Mv
    {
        List<double> _mv = new();

        const double MLARGE = 9999.9999;
        public const int Mvs = 0;
        public  const int Mvi = 1;
        public const int Mvj = 2;
        public const int Mvk = 3;
        public const int Mvij = 4;
        public const int Mvjk = 5;
        public const int Mvki = 6;
        public const int Mvp = 7;

        public Mv()
        {
            _mv = new List<double>();
            for (int ii = 0; ii < 8; ii++)
            {
                _mv.Add(0.0);
            }
        }

        public Mv(double[] mn)
        {
            Init();
            
            for (int ii = 0; ii < 8; ii++)
            {
                try
                {
                    _mv[ii] = mn[ii];
                }
                catch (Exception)
                {
                    _mv[ii] = 0.0;
                }
            }
        }

        public Mv(double s, double i, double j, double k, double ij, double jk, double ki, double p)
        {
            Init();

            _mv[Mvs] = s;
            _mv[Mvi] = i;
            _mv[Mvj] = j;
            _mv[Mvk] = k;
            _mv[Mvij] = ij;
            _mv[Mvjk] = jk;
            _mv[Mvki] = ki;
            _mv[Mvp] = p;
            
        }

        public Mv(double s)
        {
            Init();

            _mv[Mvs] = s;
        }

        public Mv Init()
        {
            _mv = new List<double>();
            for (int ii = 0; ii < 8; ii++)
            {
                _mv.Add(0.0);
            }
            return this;
        }

        public Mv(Mv mm)
        {

            Init();

            for (int qq = 0; qq < 8; qq++)
            {
                _mv.Add(mm._mv[qq]);
            }
        }

        public Mv(double ang, double xx, double yy, double zz)
        {
            double ang2;
            double sang;
            double den;
            Mv q = new();
           

            Init();

            ang2 = ang / 2.0;
            sang = Math.Sin(ang2);
            den = Math.Sqrt((xx * xx) + (yy * yy) + (zz * zz));

            q._mv[Mvij] = (xx * sang) / den;
            q._mv[Mvjk] = (yy * sang) / den;
            q._mv[Mvki] = (zz * sang) / den;
            q._mv[0] = Math.Cos(ang2);


            _mv[Mvs] = q._mv[Mvs];
            _mv[Mvi] = 0.0;
            _mv[Mvj] = 0.0;
            _mv[Mvk] = 0.0;
            _mv[Mvij] = q._mv[Mvij];
            _mv[Mvjk] = q._mv[Mvjk];
            _mv[Mvki] = q._mv[Mvki];
            _mv[Mvp] = 0.0;
        }

      public  Mv(double ang, Mv vv)
        {
            double ang2;
            double sang;
            double den;
            Mv q = new();

            ang2 = ang / 2.0;
            sang = Math.Sin(ang2);
            den = vv.Mag();

            q._mv[Mvij] = (vv[Mvi] * sang) / den;
            q._mv[Mvjk] = (vv[Mvj] * sang) / den;
            q._mv[Mvki] = (vv[Mvk] * sang) / den;

            q._mv[0] = Math.Cos(ang2);

            _mv[Mvs] = q._mv[Mvs];
            _mv[Mvi] = 0.0;
            _mv[Mvj] = 0.0;
            _mv[Mvk] = 0.0;
            _mv[Mvij] = q._mv[Mvij];
            _mv[Mvjk] = q._mv[Mvjk];
            _mv[Mvki] = q._mv[Mvki];
            _mv[Mvp] = 0.0;
        }


        public Mv Inverse()
        {
            Mv r = new();
            int ii;

            for (ii = 0; ii < 8; ii++)
            {
                r._mv[ii] = _mv[ii];
            }

            Mv t = Conj(r);             // t = conj r
            Mv q = (r * t);          // q = r * conj r

            if (q._mv[0] != 0.0)
            {
                for (ii = 0; ii < 8; ii++)
                {
                    r._mv[ii] = t._mv[ii] / q._mv[0];
                }
            }
            else
            {
                for (ii = 0; ii < 8; ii++)
                {
                    r._mv[ii] = MLARGE;
                }
            }
            return (r);
        }


        public static Mv Conj(Mv m)            // conj
        {
            double a;
            Mv vv = new();

            for (int ii = 0; ii < 8; ii++)
            {
                if (ii < 4)
                {
                    vv._mv[ii] = m._mv[ii];
                }
                else
                {
                    a = m._mv[ii];
                    vv._mv[ii] = (-a);
                }
            }
            return (vv);
        }

        public double Mag()
        {
            Mv r = new();
            double rt;
            int ii;

            for (ii = 0; ii < 8; ii++)
            {
                r._mv[ii] = _mv[ii];
            }

            Mv t = Conj(r);
            Mv q = (r * t);

            if (q._mv[0] != 0.0)
            {
                rt = q._mv[0];
                rt = Math.Sqrt(Math.Abs(rt));
            }
            else
            {
                rt = 0.0;
            }
            return (rt);
        }

        public double MagSq()
        {
            Mv r = new();
            int ii;

            for (ii = 0; ii < 8; ii++)
            {
                r._mv[ii] = _mv[ii];
            }

            Mv t = Conj(r);
            Mv q = (r * t);

            return (q._mv[0]);
        }

        public static Mv operator *(double a, Mv b)
        {
            Mv vv = new();

            vv._mv[Mvs] = b._mv[Mvs] * a;
            vv._mv[Mvs] = b._mv[Mvi] * a;
            vv._mv[Mvs] = b._mv[Mvj] * a;
            vv._mv[Mvs] = b._mv[Mvk] * a;
            vv._mv[Mvs] = b._mv[Mvij] * a;
            vv._mv[Mvs] = b._mv[Mvjk] * a;
            vv._mv[Mvs] = b._mv[Mvki] * a;
            vv._mv[Mvs] = b._mv[Mvp] * a;
            return (vv);
        }

        public static Mv operator *(Mv b, double a)
        {
            Mv vv = new();

            vv._mv[Mvs] = b._mv[Mvs] * a;
            vv._mv[Mvs] = b._mv[Mvi] * a;
            vv._mv[Mvs] = b._mv[Mvj] * a;
            vv._mv[Mvs] = b._mv[Mvk] * a;
            vv._mv[Mvs] = b._mv[Mvij] * a;
            vv._mv[Mvs] = b._mv[Mvjk] * a;
            vv._mv[Mvs] = b._mv[Mvki] * a;
            vv._mv[Mvs] = b._mv[Mvp] * a;
            return (vv);
        }

        public static Mv operator /(Mv b, double a)
        {
            if (a == 0.0)
            {
                throw new DivideByZeroException();
            }
            
            Mv vv = new();

            vv._mv[Mvs] = b._mv[Mvs] / a;
            vv._mv[Mvs] = b._mv[Mvi] / a;
            vv._mv[Mvs] = b._mv[Mvj] / a;
            vv._mv[Mvs] = b._mv[Mvk] / a;
            vv._mv[Mvs] = b._mv[Mvij] / a;
            vv._mv[Mvs] = b._mv[Mvjk] / a;
            vv._mv[Mvs] = b._mv[Mvki] / a;
            vv._mv[Mvs] = b._mv[Mvp] / a;
            return (vv);
        }

        public Mv UnitMv()
        {
            Mv nv = new();
            Mv qv = new();
            double db;
            double a;
            int ii;

            for (ii = 0; ii < 8; ii++)
            {
                qv._mv[ii] = _mv[ii];
            }

            db = qv.Mag();

            if (db != 0.0)
            {
                for (ii = 0; ii < 8; ii++)
                {
                    a = _mv[ii];
                    nv._mv[ii] = a / db;
                }
            }
            else
            {
                for (ii = 0; ii < 8; ii++)
                {
                    nv._mv[ii] = MLARGE;
                }
            }
            return (nv);
        }

        public double this[int index]
        {
            get
            {
                if (index < 0 || index > 7)
                    throw new IndexOutOfRangeException("index must be between 0 and 7");
                return _mv[index];
            }

            set
            {
                if (index < 0 || index > 7)
                    throw new IndexOutOfRangeException("index must be between 0 and 7");
                _mv[index] = value;
            }
        }

        public static Mv operator *(Mv a, Mv m)
        {
            Mv nv = new();

            nv._mv[Mvs] = a._mv[Mvs] * m._mv[Mvs];     //MNA
            nv._mv[Mvi] = a._mv[Mvi] * m._mv[Mvs];    //MNA
            nv._mv[Mvj] = a._mv[Mvj] * m._mv[Mvs];     //MNA
            nv._mv[Mvk] = a._mv[Mvk] * m._mv[Mvs];     //MNA
            nv._mv[Mvij] = a._mv[Mvij] * m._mv[Mvs];   //MNA
            nv._mv[Mvjk] = a._mv[Mvjk] * m._mv[Mvs];   //MNA
            nv._mv[Mvki] = a._mv[Mvki] * m._mv[Mvs];   //MNA
            nv._mv[Mvp] = a._mv[Mvp] * m._mv[Mvs];     //MNA

             nv._mv[Mvi] += ( a._mv[Mvs] * m._mv[Mvi]);   //MA
             nv._mv[Mvs] += ( a._mv[Mvi] *m._mv[Mvi]);   //MA
             nv._mv[Mvij] -= ( a._mv[Mvj] * m._mv[Mvi]);  //MS
             nv._mv[Mvki] += ( a._mv[Mvk] * m._mv[Mvi]);  //MA
             nv._mv[Mvj] -= ( a._mv[Mvij] * m._mv[Mvi]);  //MS
             nv._mv[Mvp] += ( a._mv[Mvjk] * m._mv[Mvi]);  //MA
             nv._mv[Mvk] += ( a._mv[Mvki] * m._mv[Mvi]);  //MA
             nv._mv[Mvjk] += ( a._mv[Mvp] * m._mv[Mvi]);  //MA

             nv._mv[Mvj] += ( a._mv[Mvs] * m._mv[Mvj]);  //MA
             nv._mv[Mvij] += ( a._mv[Mvi] * m._mv[Mvj]);  //MA
             nv._mv[Mvs] += ( a._mv[Mvj] * m._mv[Mvj]);  //MA
             nv._mv[Mvjk] -= ( a._mv[Mvk] * m._mv[Mvj]);  //MS
             nv._mv[Mvi] += ( a._mv[Mvij] * m._mv[Mvj]);  //MA
             nv._mv[Mvk] -= ( a._mv[Mvjk] * m._mv[Mvj]);   //MS
             nv._mv[Mvp] += (a._mv[Mvki] * m._mv[Mvj]);  //MA
             nv._mv[Mvki] += ( a._mv[Mvp] * m._mv[Mvj]);  //MA

             nv._mv[Mvk] += ( a._mv[Mvs] * m._mv[Mvk]);  //MA
             nv._mv[Mvki] -= ( a._mv[Mvi] * m._mv[Mvk]);  //MS
             nv._mv[Mvjk] += ( a._mv[Mvj] * m._mv[Mvk]);  //MA
             nv._mv[Mvs] += ( a._mv[Mvk] * m._mv[Mvk]);  //MA
             nv._mv[Mvp] += ( a._mv[Mvij] * m._mv[Mvk]);  //MA
             nv._mv[Mvj] += (a._mv[Mvjk] * m._mv[Mvk]);  //MA
             nv._mv[Mvi] -= ( a._mv[Mvki] * m._mv[Mvk]);  //MS
             nv._mv[Mvij] += ( a._mv[Mvp] * m._mv[Mvk]);  //MA

             nv._mv[Mvij] += ( a._mv[Mvs] * m._mv[Mvij]);  //MA
             nv._mv[Mvj] += ( a._mv[Mvi] * m._mv[Mvij]);  //MA
             nv._mv[Mvi] -= ( a._mv[Mvj] * m._mv[Mvij]);  //MS
             nv._mv[Mvp] += ( a._mv[Mvk] * m._mv[Mvij]);  //MA
             nv._mv[Mvs] -= ( a._mv[Mvij] * m._mv[Mvij]);  //MS
             nv._mv[Mvki] += ( a._mv[Mvjk] * m._mv[Mvij]);  //MA
             nv._mv[Mvjk] -= (a._mv[Mvki] * m._mv[Mvij]);  //MS
             nv._mv[Mvk] -= (a._mv[Mvp] * m._mv[Mvij]);  //MS

             nv._mv[Mvjk] += ( a._mv[Mvs] *  m._mv[Mvjk]);  //MA
             nv._mv[Mvp] += ( a._mv[Mvi] *  m._mv[Mvjk]);  //MA
             nv._mv[Mvk] += (a._mv[Mvj] *  m._mv[Mvjk]);  //MA
             nv._mv[Mvj] -= ( a._mv[Mvk] *  m._mv[Mvjk]);  //MS
             nv._mv[Mvki] -= ( a._mv[Mvij] *  m._mv[Mvjk]);  //MS
             nv._mv[Mvs] -= ( a._mv[Mvjk] *  m._mv[Mvjk]);  //MS
             nv._mv[Mvij] += ( a._mv[Mvki] *  m._mv[Mvjk]);  //MA
             nv._mv[Mvi] -= ( a._mv[Mvp] *  m._mv[Mvjk]);  //MS

             nv._mv[Mvki] += ( a._mv[Mvs] *  m._mv[Mvki]);  //MA
             nv._mv[Mvk] -= ( a._mv[Mvi] *  m._mv[Mvki]);  //MS
             nv._mv[Mvp] += ( a._mv[Mvj] *  m._mv[Mvki]);  //MA
             nv._mv[Mvi] += ( a._mv[Mvk] *  m._mv[Mvki]);  //MA
             nv._mv[Mvjk] += ( a._mv[Mvij] *  m._mv[Mvki]);  //MA
             nv._mv[Mvij] -= ( a._mv[Mvjk] *  m._mv[Mvki]);  //MS
             nv._mv[Mvs] -= ( a._mv[Mvki] *  m._mv[Mvki]);  //MS
             nv._mv[Mvj] -= ( a._mv[Mvp] *  m._mv[Mvki]);  //MS

             nv._mv[Mvp] += ( a._mv[Mvs] *  m._mv[Mvp]);  //MA
             nv._mv[Mvjk] += ( a._mv[Mvi] *  m._mv[Mvp]);  //MA
             nv._mv[Mvki] += ( a._mv[Mvj] *  m._mv[Mvp]);  //MA
             nv._mv[Mvij] += ( a._mv[Mvk] *  m._mv[Mvp]);  //MA
             nv._mv[Mvk] -= ( a._mv[Mvij] *  m._mv[Mvp]);  //MS
             nv._mv[Mvi] -= ( a._mv[Mvjk] *  m._mv[Mvp]);  //MS
             nv._mv[Mvj] -= ( a._mv[Mvki] *  m._mv[Mvp]);  //MS
             nv._mv[Mvs] -= ( a._mv[Mvp] *  m._mv[Mvp]);  //MS
            return (nv);
        }

        //public static Mv operator *(Mv a, Mv m)
        //{
        //    Mv nv = new Mv();

        //    MNA(nv._mv[Mvs], a._mv[Mvs], m._mv[Mvs]);
        //    MNA(nv._mv[Mvi], a._mv[Mvi], m._mv[Mvs]);
        //    MNA(nv._mv[Mvj], a._mv[Mvj], m._mv[Mvs]);
        //    MNA(nv._mv[Mvk], a._mv[Mvk], m._mv[Mvs]);
        //    MNA(nv._mv[Mvij], a._mv[Mvij], m._mv[Mvs]);
        //    MNA(nv._mv[Mvjk], a._mv[Mvjk], m._mv[Mvs]);
        //    MNA(nv._mv[Mvki], a._mv[Mvki], m._mv[Mvs]);
        //    MNA(nv._mv[Mvp], a._mv[Mvp], m._mv[Mvs]);

        //    MA(nv._mv[Mvi], a._mv[Mvs], m._mv[Mvi]);
        //    MA(nv._mv[Mvs], a._mv[Mvi], m._mv[Mvi]);
        //    MS(nv._mv[Mvij], a._mv[Mvj], m._mv[Mvi]);
        //    MA(nv._mv[Mvki], a._mv[Mvk], m._mv[Mvi]);
        //    MS(nv._mv[Mvj], a._mv[Mvij], m._mv[Mvi]);
        //    MA(nv._mv[Mvp], a._mv[Mvjk], m._mv[Mvi]);
        //    MA(nv._mv[Mvk], a._mv[Mvki], m._mv[Mvi]);
        //    MA(nv._mv[Mvjk], a._mv[Mvp], m._mv[Mvi]);

        //    MA(nv._mv[Mvj], a._mv[Mvs], m._mv[Mvj]);
        //    MA(nv._mv[Mvij], a._mv[Mvi], m._mv[Mvj]);
        //    MA(nv._mv[Mvs], a._mv[Mvj], m._mv[Mvj]);
        //    MS(nv._mv[Mvjk], a._mv[Mvk], m._mv[Mvj]);
        //    MA(nv._mv[Mvi], a._mv[Mvij], m._mv[Mvj]);
        //    MS(nv._mv[Mvk], a._mv[Mvjk], m._mv[Mvj]);
        //    MA(nv._mv[Mvp], a._mv[Mvki], m._mv[Mvj]);
        //    MA(nv._mv[Mvki], a._mv[Mvp], m._mv[Mvj]);

        //    MA(nv._mv[Mvk], a._mv[Mvs], m._mv[Mvk]);
        //    MS(nv._mv[Mvki], a._mv[Mvi], m._mv[Mvk]);
        //    MA(nv._mv[Mvjk], a._mv[Mvj], m._mv[Mvk]);
        //    MA(nv._mv[Mvs], a._mv[Mvk], m._mv[Mvk]);
        //    MA(nv._mv[Mvp], a._mv[Mvij], m._mv[Mvk]);
        //    MA(nv._mv[Mvj], a._mv[Mvjk], m._mv[Mvk]);
        //    MS(nv._mv[Mvi], a._mv[Mvki], m._mv[Mvk]);
        //    MA(nv._mv[Mvij], a._mv[Mvp], m._mv[Mvk]);

        //    MA(nv._mv[Mvij], a._mv[Mvs], m._mv[Mvij]);
        //    MA(nv._mv[Mvj], a._mv[Mvi], m._mv[Mvij]);
        //    MS(nv._mv[Mvi], a._mv[Mvj], m._mv[Mvij]);
        //    MA(nv._mv[Mvp], a._mv[Mvk], m._mv[Mvij]);
        //    MS(nv._mv[Mvs], a._mv[Mvij], m._mv[Mvij]);
        //    MA(nv._mv[Mvki], a._mv[Mvjk], m._mv[Mvij]);
        //    MS(nv._mv[Mvjk], a._mv[Mvki], m._mv[Mvij]);
        //    MS(nv._mv[Mvk], a._mv[Mvp], m._mv[Mvij]);

        //    MA(nv._mv[Mvjk], a._mv[Mvs], m._mv[Mvjk]);
        //    MA(nv._mv[Mvp], a._mv[Mvi], m._mv[Mvjk]);
        //    MA(nv._mv[Mvk], a._mv[Mvj], m._mv[Mvjk]);
        //    MS(nv._mv[Mvj], a._mv[Mvk], m._mv[Mvjk]);
        //    MS(nv._mv[Mvki], a._mv[Mvij], m._mv[Mvjk]);
        //    MS(nv._mv[Mvs], a._mv[Mvjk], m._mv[Mvjk]);
        //    MA(nv._mv[Mvij], a._mv[Mvki], m._mv[Mvjk]);
        //    MS(nv._mv[Mvi], a._mv[Mvp], m._mv[Mvjk]);

        //    MA(nv._mv[Mvki], a._mv[Mvs], m._mv[Mvki]);
        //    MS(nv._mv[Mvk], a._mv[Mvi], m._mv[Mvki]);
        //    MA(nv._mv[Mvp], a._mv[Mvj], m._mv[Mvki]);
        //    MA(nv._mv[Mvi], a._mv[Mvk], m._mv[Mvki]);
        //    MA(nv._mv[Mvjk], a._mv[Mvij], m._mv[Mvki]);
        //    MS(nv._mv[Mvij], a._mv[Mvjk], m._mv[Mvki]);
        //    MS(nv._mv[Mvs], a._mv[Mvki], m._mv[Mvki]);
        //    MS(nv._mv[Mvj], a._mv[Mvp], m._mv[Mvki]);

        //    MA(nv._mv[Mvp], a._mv[Mvs], m._mv[Mvp]);
        //    MA(nv._mv[Mvjk], a._mv[Mvi], m._mv[Mvp]);
        //    MA(nv._mv[Mvki], a._mv[Mvj], m._mv[Mvp]);
        //    MA(nv._mv[Mvij], a._mv[Mvk], m._mv[Mvp]);
        //    MS(nv._mv[Mvk], a._mv[Mvij], m._mv[Mvp]);
        //    MS(nv._mv[Mvi], a._mv[Mvjk], m._mv[Mvp]);
        //    MS(nv._mv[Mvj], a._mv[Mvki], m._mv[Mvp]);
        //    MS(nv._mv[Mvs], a._mv[Mvp], m._mv[Mvp]);
        //    return (nv);

        //}

        public static Mv operator /(Mv a, Mv m)
        {
            Mv av = new();
            Mv qv = new();
            int ii;

            for (ii = 0; ii < 8; ii++)
            {
                qv._mv[ii] = a._mv[ii];
                av._mv[ii] = m._mv[ii];
            }

            Mv tv = av.Inverse();
            Mv nv = (qv * tv);
            return (nv);
        }

        
        public static Mv operator +(Mv a, Mv m)
        {
            Mv nv = new();

            for (int ii = 0; ii < 8; ii++)
            {
                nv._mv[ii] = a._mv[ii] + m._mv[ii];
            }
            return nv;
        }

        public static Mv operator -(Mv a, Mv m)
        {
            Mv nv = new();

            for (int ii = 0; ii < 8; ii++)
            {
                nv._mv[ii] = a._mv[ii] - m._mv[ii];
            }
            return nv;
        }

        public static bool operator ==(Mv a, Mv m)
        {
            for (int ii = 0; ii < 8; ii++)
            {
                if (a._mv[ii] != m._mv[ii])
                {
                    return (false);
                }
            }
            return (true);
        }

        public static bool operator !=(Mv a, Mv m)
        {
            for (int ii = 0; ii < 8; ii++)
            {
                if (a._mv[ii] != m._mv[ii])
                {
                    return (true);
                }
            }
            return (false);
        }

        public override bool Equals(object obj)
        {
            return (this == (Mv)obj);
        }

        public override int GetHashCode()
        {
            int ii = _mv[0].GetHashCode();
            for (int kk = 1; kk < 8; kk++)
            {
                ii ^= (int)_mv[kk].GetHashCode();
            }

            return ii;
        }

        public Mv VectorOf()
        {
            Mv vv = new();

            vv[Mvj] = _mv[Mvi];
            vv[Mvj] = _mv[Mvj];
            vv[Mvj] = _mv[Mvk];
            return (vv);
        }

        public static double AngVect(Mv va, Mv vb)
        {
            double dotvect = (va[Mvi] * vb[Mvi]) + (va[Mvj] * vb[Mvj]) + (va[Mvk] * vb[Mvk]);
            double vamag = Math.Sqrt(va[Mvi] * va[Mvi]) + (va[Mvj] * va[Mvj]) + (va[Mvk] * va[Mvk]);
            double vbmag = Math.Sqrt(vb[Mvi] * vb[Mvi]) + (vb[Mvj] * vb[Mvj]) + (vb[Mvk] * vb[Mvk]);
            double ang = Math.Acos(dotvect / (vamag * vbmag));
            return (ang);
        }

        public static Mv VectCross(Mv a, Mv b)
        {
            Mv vv = new();

            vv._mv[Mvi] = (a._mv[Mvj] * b._mv[Mvk]) - (b._mv[Mvj] * a._mv[Mvk]);
            vv._mv[Mvj] = -(a._mv[Mvi] * b._mv[Mvk]) + (b._mv[Mvi] * a._mv[Mvk]);
            vv._mv[Mvk] = (a._mv[Mvi] * b._mv[Mvj]) - (b._mv[Mvi] * a._mv[Mvj]);

            return (vv);
        }

        public static Mv GetSpinor(Mv va, Mv vb)
        {
            double ang;
            double ang2;
            double sang;
            double den;
            Mv q = new();

            ang = AngVect(va, vb);
            ang2 = ang / 2.0;
            sang = Math.Sin(ang2);
            Mv vv = VectCross(va, vb);
            den = vv.Mag();

            q._mv[Mvij] = (vv[Mvi] * sang) / den;
            q._mv[Mvjk] = (vv[Mvj] * sang) / den;
            q._mv[Mvki] = (vv[Mvk] * sang) / den;
            q._mv[Mvs] = Math.Cos(ang2);

            return (q);
        }

        public Mv Spin(Mv vv)
        {
            Mv conjth = Conj(this);
            Mv rmv = this * (vv * conjth);
            return (rmv.VectorOf());
        }

        public override string ToString()
        {
            string rt = "";
            for (int ii = 0; ii < 8; ii++)
            {
                rt += (_mv[ii].ToString() + " ");
            }
            return (rt);
        }

     public static Mv Exp(Mv tt)
        {
            double fact = 1.0;
            Mv sumb = new();
            Mv cv = new();
            cv[Mvs] = 1.0;
            sumb[Mvs] = 1.0;

            for (int ii = 1; ii <= 10; ii++)
            {
                fact *= ((double)ii);
                cv = cv * tt;
                sumb = sumb + (cv * (1.0/fact));
            }
            return (sumb);
        }

        public static Mv Sinh(Mv tt)
        {
            double fact = 1.0;
            Mv sumb = new();
            Mv cv = new();
            cv = tt;
            sumb = tt;

            for (int ii = 2; ii <= 11; ii++)
            {
                fact *= ((double)ii);
                cv = cv * tt;
                if ((ii % 2) == 1)
                {
                    sumb = sumb + (cv * (1.0/fact));
                }
            }
            return (sumb);
        }


        public static Mv Cosh(Mv tt)
        {
            double fact = 1.0;
            Mv sumb = new();
            Mv cv = new();
            sumb[Mvs] = 1.0;

            for (int ii = 2; ii <= 10; ii++)
            {
                fact *= ((double)ii);
                cv = cv * tt;
                if ((ii % 2) == 0)
                {
                    sumb = sumb + (cv * (1.0/fact));
                }
            }
            return (sumb);
        }

        public static  Mv Sin(Mv tt)
        {
            double fact = 1.0;
            Mv sumb = new();
            Mv cv = new();
            cv = tt;
            sumb = tt;

            for (int ii = 2; ii <= 11; ii++)
            {
                fact *= ((double)ii);
                cv = cv * tt;
                if ((ii % 2) == 1)
                {
                    if ((ii % 4) == 1)
                    {
                        sumb = sumb + (cv * (1.0/fact));
                    }
                    else
                    {
                        sumb = sumb - (cv * (1.0/fact));
                    }
                }
            }
            return (sumb);
        }

        public static Mv Cos(Mv tt)
        {
            double fact = 1.0;
            Mv sumb = new();
            Mv cv = new();
            sumb[Mvs] = 1.0;

            for (int ii = 2; ii <= 10; ii++)
            {
                fact *= ((double)ii);
                cv = cv * tt;
                if ((ii % 2) == 0)
                {
                    if ((ii % 4) == 0)
                    {
                        sumb = sumb + (cv * (1.0/fact));
                    }
                    else
                    {
                        sumb = sumb - (cv * (1.0/fact));
                    }
                }
            }
            return (sumb);
        }


    }
}
