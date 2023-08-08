This is a basic and general program for element differential method (EDM), a numerical algorithm proposed by Xiao-Wei Gao et al.
The program is developed by Xiao-Wei Gao, Yong-Tong Zheng et al.

The program can solve linear steady heat transfer problems (when ndf=1, istrans=0, isnonli=0),
                      linear elasticity problems (when ndf=2 or 3, istrans=0, isnonli=0),
                      linear transient heat transfer problems (when ndf=1, istrans=1, isnonli=0),
                      linear elastodynamic response problems (when ndf=2 or 3, istrans=1, isnonli=0),
                      non-linear steady heat transfer problems (when ndf=1, istrans=0, isnonli=1),
                  and non-linear transient heat transfer problems (when ndf=1, istrans=1, isnonli=1).

Before using the program EDM.exe:
  1. Some basic knowledge (governing equation and boundary condition etc.) about heat conduction and mechanics should be known.

  2. A FORTRAN compiler should be installed in the computer and it had better to be Intel Visual Fortran or Intel One API, otherwise the linear algebraic solver PARDISO, DSS,
      FGMRES_ILU0, FGMRES_ILUT and some other subroutines about matrix adding or multiplying may not work. It is highly recommended to use:
      Microsoft Visual Studio Community (MVS) (https://visualstudio.microsoft.com/downloads),
      Intel OneAPI Base Toolkit and HPC Toolkit (https://www.intel.com/content/www/us/en/developer/tools/oneapi/toolkits.html#base-kit) or Intel visual FORTRAN (IVF).
      If one want to use some other compilers, It is recommended to find some other sparse matrix solvers and replace the origin solvers with them by modifying the codes. 

  3. It is highly recommended to read the EDM_User_Manual.pdf or EDM_User_MANUAL_English.pdf to learn how to make a input file (EDM.inp)
      and what files will be generated after running the program.

Before reading the source codes, it is highly recommended to refer to the following five papers:
    International Journal of Heat and Mass Transfer 115 (2017) 882–894 ----- linear steady heat transfer problems
    International Journal for Numerical Methods in Engineering vol.: 113, issue: 1,(2018) 82–108 ----- linear elasticity problems
    International Journal of Heat and Mass Transfer 127 (2018) 1189–1197 ----- linear transient heat transfer problems
    International Journal of Mechanical Sciences 151 (2019) 828–841 ----- linear elastodynamic response problems
    International Journal of Heat and Mass Transfer 126 (2018) 1111–1119 ----- non-linear steady and transient heat transfer problems

How to generate EDM.exe by the source codes?
    If using Intel OneAPI or IVF, to compile the codes into a executable program named EDM.exe, one can use the following compile command:
    IFORT /Qmkl:parallel /exe:EDM.exe *.F *.f90
    If compiling them in IDE of MVS, one still need to choose 'Use Intel Math Kernel Library' in 'Property' (not to be 'No').

If one want to get more details or give some suggestions, please contact: prof. Gao via email: xwgao@dlut.edu.cn