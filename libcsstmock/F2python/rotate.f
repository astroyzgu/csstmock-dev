c*****************************************************************
      subroutine randomize(gaxsim,Nd1,Nd2,Ntreat)
c================================================================
c  Randomize the cordinates of the galaxies.
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2,Ntreat,i8

      REAL    gaxsim(Nd1,Nd2)
      REAL ran1
      External ran1

      do i8=1,Ntreat
        do j=1,3
          gaxsim(j,i8)=ran1(iseed)*rLbox
        enddo
      enddo

      return
      end


c*****************************************************************
      subroutine rotate(r,Nd1,Nd2,Ntreat)
c================================================================
c  Randomize the cordinates of the galaxies.
c================================================================
      include 'common.inc'
c================================================================
      INTEGER Nd1
      INTEGER*8 Nd2,Ntreat,i8

      REAL    r(Nd1,Nd2)
      REAL x1,x2,x3,v1,v2,v3

      do idble=1,Ntreat
        x1 = r(1,idble)
        x2 = r(2,idble)
        x3 = r(3,idble)
        v1 = r(4,idble)
        v2 = r(5,idble)
        v3 = r(6,idble)
        r(1,idble) = mod(rLbox+x2,rLbox)
        r(2,idble) = mod(rLbox+x3,rLbox)
        r(3,idble) = mod(rLbox*0.4+x1,rLbox)
        r(4,idble) = v2
        r(5,idble) = v3
        r(6,idble) = v1
      enddo

      return
      end

     


