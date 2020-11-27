      module mndo
      use sizes

      implicit none
      integer nqmatoms
      integer qmlist(maxatm)
      
      integer mndot_nline, mndot_oline, mndot_eline
      character(len=1024), allocatable :: mndo_template(:) 
      save
      end
