        PROGRAM test 
        use, intrinsic :: iso_c_binding ! 和C语言类型互通的模块iso_c_binding
        implicit none
        integer(4)  :: snapnum, ifile, filenum  
        integer(8)  :: nmax, nsub
        integer(c_long), allocatable :: larr(:,:)
        real(c_float),   allocatable :: darr(:,:)
        nmax    = 500000000
        allocate( larr( 5, nmax) )
        allocate( darr(20, nmax) )
        do snapnum = 127, 120, -1 

            !>>> to know how many files at the given snapshot number. 
            call fullsky_z2_filenum_jiutian(snapnum, filenum) 
	    print*, snapnum, filenum 

            do ifile = 1,filenum
               nsub    = nmax 
               !>>> to read the brach file one by one. 
               call fullsky_z2_jiutian(snapnum, ifile-1, 
     &                                 larr, darr, nsub)
               !>>> operation herein <<<!
               ! print*, snapnum, ifile, filenum, nsub  
               ! write(*,*)(larr(:,1) )  
            end do 
       end do 
       END 
