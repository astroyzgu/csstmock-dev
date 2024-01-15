      subroutine rangeLM(rg,Nin)
c================================================================
      include 'common.inc'
c================================================================
      real rg(Nin)
      real xMmax,xMmin,xMst
      INTEGER*8 Nin,i1

      xMmax=0.
      xMmin=100.
      do i1=1,Nin
          xMst=rg(i1) 
          xMmax=max(xMst,xMmax)
          xMmin=min(xMst,xMmin)
          IF(abs(xMst).gt.100.) print*,i1,rg(i1)
      enddo
      xlum_1=xMmin         !! gax lum min

      xlum_1=max(xlum_cut,xMmin)

      xlum_2=xMmax        !! gax lum max
      dxlum=(xlum_2-xlum_1)/float(N_GLF-1)
      WRITE(*,*) 'Min and Max log Mst or L:', xlum_1,xlum_2

      RETURN
      END

