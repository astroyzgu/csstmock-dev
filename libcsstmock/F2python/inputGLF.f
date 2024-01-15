c*************************************************************
      subroutine getGLFfromLFdata
c================================================================
c  Get the CUMULATIVE galaxy luminosity functions from LF data,
c  which are used to get the stellar-luminosity relations.
c  Writen by: Xiaohu Yang
c================================================================ 
      include 'common.inc'
c================================================================
    
      real xl_LF(N_GLF),all_LF(N_GLF),cent_LF(N_GLF),sate_LF(N_GLF)
      integer, external:: lblnk
      integer NLFobs

      IF(CLFfile(1:lblnk(CLFfile)).eq.'LF_sdss_Yang2009') THEN
         call inputLF_Yang2009(xl_LF,all_LF,cent_LF,sate_LF,NLFobs)
         call cumulative(xl_LF,all_LF,NLFobs,'a')
         call cumulative(xl_LF,cent_LF,NLFobs,'c')
         call cumulative(xl_LF,sate_LF,NLFobs,'s')
      ENDIF 


      open(3,FILE=outLF(1:lblnk(outLF)) ,status='replace')
      do i=1,N_GLF
         write(3,142) GLFa(1:2,i),GLFc(2,i),GLFs(2,i)
      enddo
      close(3)
142   format(1x,f9.5,6(1x,e11.4))

      print*,'Input galaxy luminosity function OK.'

      return
      end

c*************************************************************

      subroutine inputLF_Yang2009(xl_LF,all_LF,cent_LF,sate_LF,NLFobs)
c================================================================ 
      include 'common.inc'
c================================================================

      real xl_LF(N_GLF),all_LF(N_GLF),cent_LF(N_GLF),sate_LF(N_GLF)
      integer, external:: lblnk
      integer NLFobs
      real xxx(20),xxxmin

      xxxmin=1.0e-6

      NLFobs=43
      open(63,FILE='./MF/obsLFs/'//CLFfile(1:lblnk(CLFfile))
     x   ,status='old')
      do i=1,NLFobs
          read(63,*) xl_LF(i),xxx(1:18)
          all_LF(i)=log10(xxx(1)+xxxmin)-2.0
          cent_LF(i)=log10(xxx(7)+xxxmin)-2.0
          sate_LF(i)=log10(xxx(13)+xxxmin)-2.0
c          print*,xl_LF(i),all_LF(i),cent_LF(i),sate_LF(i)
      enddo
      close(63)

      return
      end
     

c*************************************************************

      subroutine cumulative(xl_LF,data_LF,NLFobs,type_LF)
c================================================================ 
      include 'common.inc'
c================================================================

      real xl_LF(N_GLF),data_LF(N_GLF),xlmbin(N_GLF),temp_LF(N_GLF)
      integer NLFobs
      character type_LF
      real xlum_a,dxluma,xlum_i,xLF_i
      REAL     r1,r2,z1,z2

      xlum_a=6.0
      dxluma=0.08
      do i=1,N_GLF
         xlmbin(i)=xlum_a+dxluma*(i-1)
         temp_LF(i)=0.
      enddo

      do i=N_GLF,1,-1
         xlum_i=xlmbin(i)+0.5*dxluma
         CALL locate(xl_LF(1:NLFobs),NLFobs,xlum_i,j)

        IF (j.EQ.0) THEN
           xLF_i=data_LF(1)
           xLF_i=10.**xLF_i
           GOTO 111
        END IF

        IF (j.EQ.NLFobs) THEN
           xLF_i=0.
           GOTO 111
        END IF
      
        r1 = xl_LF(j)
        r2 = xl_LF(j+1)

        z1 = data_LF(j)
        z2 = data_LF(j+1)

        xLF_i = z1 + ((xlum_i-r1)/(r2-r1)) * (z2-z1)
        xLF_i=10.**xLF_i

111     do k=1,i
        temp_LF(k)=temp_LF(k)+xLF_i*dxluma
        enddo
            
        IF(type_LF.eq.'a') THEN
           GLFa(1,i)=xlmbin(i)
           GLFa(2:Nzbin+1,i)=temp_LF(i)**(-1./3.)  !!! average distance
        ENDIF

        IF(type_LF.eq.'c') THEN
           GLFc(1,i)=xlmbin(i)
           GLFc(2:Nzbin+1,i)=temp_LF(i)**(-1./3.)  !!! average distance
        ENDIF

        IF(type_LF.eq.'s') THEN
           GLFs(1,i)=xlmbin(i)
           GLFs(2:Nzbin+1,i)=temp_LF(i)**(-1./3.)  !!! average distance
        ENDIF

      enddo


      return
      end


