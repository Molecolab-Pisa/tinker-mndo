c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ################################################################
c     ##                                                            ##
c     ##  program gradintergace  --  Output QM/MM potential to      ##
c     ##  an interface file                                         ##
c     ##                                                            ##
c     ################################################################
c
c
      program potinterface
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
      real*8 energy

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
c     compute the analytical gradient components
c
      etot = energy()
c
c     Output an interface file.
c

      iintf = freeunit()
      intffile = 'interface.dat'

      open(unit=iintf, file=intffile)
      
      if(mndo_nstates .gt. 0) then
        write(iintf, '(a, i10)') "POTENTIAL", mndo_nstates
        do i=1, mndo_nstates
          write(iintf, '(f20.8)') etot - emndo_tmp(mndo_currentstate) 
     &                            + emndo_tmp(i)
        end do
      else
        write(iintf, '(a, i10)') "POTENTIAL", 1
        write(iintf, '(f20.8)') etot
      end if 

      close(unit=iintf)
c
c     perform any final tasks before program exit
c
      call final
      end
