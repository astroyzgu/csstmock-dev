      SUBROUTINE initzbins
c================================================================

      include 'common.inc'

c================================================================
      REAL get_r,get_z,reff
      EXTERNAL get_r,get_z

      fact=skycover/(180./pi)**2/3
      dzbin=(z_cut2-z_cut1)/Nzbin
      do i=1,Nzbin
         z_min(i)=z_cut1+dzbin*(i-1.0)
         z_max(i)=z_min(i)+dzbin
         z_min(1)=z_cut0   !!! set differnt cut redshift
         v_min(i)=fact*(get_r(z_min(i)))**3
         v_max(i)=fact*(get_r(z_max(i)))**3
         v_eff(i)=v_max(i)-v_min(i)
         reff=((v_min(i)+v_eff(i)/2.)/fact)**(1./3.)
         z_eff(i)=get_z(reff)  
      enddo

      print*,'z_eff in Nzbin:',z_cut1,z_cut2, z_eff(1:Nzbin)

      print*,' '
      
      RETURN
      END


