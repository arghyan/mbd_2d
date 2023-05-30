# mbd_2d_1
Multibody Dynamics - 2d example problems

Example 1 - Double Pendulum

            1a. Lagrange's equations, q1, q2 measured from inertial x axis - Lagrange_double_pend1
            
            1b.  Lagrange's equations, q1 from inertial x axis, q2 from pendulum 1 axis - Lagrange_double_pend2
            
            1c.  sympy, q1, q2 measured from inertial x axis - Kane_double_pend_sympy1
            
            1d.  sympy, q1 from inertial x axis, q2 from pendulum 1 axis - Kane_double_pend_sympy2
            
            1e.  Amirouche's approach, q1, q2 measured from inertial x axis - dbp_mbd1
            
            1f.  Amirouche's approach, q1 from inertial x axis, q2 from pendulum 1 axis - dbp_mbd2
            
            1g. Lagrange's equations, use of lambify to get numerical equations, use of odeint of scipy and plotting 
                using matplotlib - Lagrange_double_pend3
                
             1h & 1i. The ode of 1g integrated using ode45 of octave - dydt_dbp.m and db_pend_rkutta.m
            
            Youtube link for Amirouche approach - https://youtu.be/bxrilrDBvZw?list=PLboiG5KX48xiiuED1r5YcBWlwjZDw6fOL
            
 Example 2  -  Amirouche Book Page 338 Ex 7.7.2. For checking application of constraints - amirouche_ex2.jpg
            
            2a & b amirouche_ex2.ipynb, amirouche_ex_kane.ipynb - Lagrange's method with numerical integration and Kane's
            method respectively. Lagrange's method of sympy uses augmented constraint equations. Lagrange multipliers along with
            generalised accelerations are computed. In Kane's method only equations are derived. Kane's method in sympy uses
            substitution of lambdas.
            
            2c Amirouche_ex_kane3.py integrates sympy equations. Results compared with amirouche_ex2.ipynb, which uses 
            Lagrange's method with augmented constraints.Screenshots amirouche_ex2_ipynb_results.png &
            amirouche_ex_kane3_py_results.png
            
            2d & e Amirouche_ex_kane1.py - same as Amirouche_ex_kane3.py but without integration
                   Amirouche_ex_kane2.py - Unconstrained Kane's equation from sympy and incorporaating constraints
                   using own code segment. 2d & e are compared after numerical evaluation of matrices. For
                   UNDERSTANDING purpose and COMPARISON of mass matrix and forcing vector.
