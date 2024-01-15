c=======================================================================
c
c                      INCLUDE FILE common.inc
c
c=======================================================================

        implicit none

c---    pi 
        REAL  pi
        PARAMETER (pi=3.141592654)

c---  I/O related parameters ...

        integer :: Ngaldim, Ngalmax  
        COMMON/inputnumbers/Ngaldim,Ngalmax

c=======================================================================
c                             END
c=======================================================================
