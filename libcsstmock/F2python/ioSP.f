      SUBROUTINE inputSP(rg,Nd1,Nd2max,idrg,Nin)
c================================================================
c Input spectroscopic redshift data 
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2max,Nin,idb1

      REAL*8 rg(Nd1,Nd2max)
      INTEGER*8 idrg(Nd2max)
      REAL*8 xxx(20)

      INTEGER, external:: lblnk

      OPEN(66,FILE=indir(1:lblnk(indir))//
     &  infile(1:lblnk(infile))//file_sp(1:lblnk(file_sp)),
     & STATUS='OLD')
        do idb1=1,Nd2max
         Nin=Nin+1
         READ(66,*,end=121) idrg(Nin),xxx(1:13)
         rg(1,Nin)=xxx(3)
         rg(2,Nin)=xxx(8)
         rg(3,Nin)=xxx(10)
        enddo
121   CLOSE(66)
      
      Print*,'Galaxies with spectro redshift:', Nin

      RETURN
      END


