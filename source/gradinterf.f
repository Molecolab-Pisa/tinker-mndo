c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program gradintergace  --  Output QM/MM gradients to      ##
c     ##  an interface file                                         ##
c     ##                                                            ##
c     ################################################################
c
c
      program gradinterface
      use sizes
      use atoms
      use deriv
      use energi
      use files
      use inform
      use inter
      use iounit
      use mndo 
      use solute
      use usage
      implicit none
      integer i,j,ixyz, iintf
      integer freeunit
      real*8 etot, e0
      real*8, allocatable :: detot(:,:)

      character*240 xyzfile, intffile
c
c     set up the structure and mechanics calculation
c
      call initial
      call getxyz
      call mechanic
c
c     reopen the coordinates file and read the first structure
c
      ixyz = freeunit ()
      xyzfile = filename
      call suffix (xyzfile,'xyz','old')
      open (unit=ixyz,file=xyzfile,status ='old')
      rewind (unit=ixyz)
      call readxyz (ixyz)
c
c     perform dynamic allocation of some local arrays
c
      allocate (detot(3,n))
c
c     compute the analytical gradient components
c
      call gradient (etot,detot)
c
c     Output an interface file.
c
      iintf = freeunit()
      intffile = 'interface.dat'

      open(unit=iintf, file=intffile)
      
      if(mndo_nstates .gt. 0) then
        write(iintf, '(a, i10)') "POTENTIAL", mndo_nstates
        do i=1, mndo_nstates
          write(iintf, '(f20.8)') emndo_tmp(i)
        end do
      else
        write(iintf, '(a, i10)') "POTENTIAL", 1
        write(iintf, '(f20.8)') etot
      end if 

      write(iintf, '(a, 2i10)') "GRADIENT", mndo_currentstate
      do i=1, n
        write(iintf, '(3f20.8)') detot(:,i)
      end do

      close(unit=iintf)
c
c     perform deallocation of some local arrays
c
      deallocate (detot)
c
c     perform any final tasks before program exit
c
      call final
      end
