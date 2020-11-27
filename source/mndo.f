      module mndo
      use sizes

      implicit none
      character(len=1024), parameter ::
     & mndo_exe = "/home/mattia/work_tmp/MNDO2020/mndo2020/mndo2020",
     & mndo_in = "tmp_mndo.inp",
     & mndo_out = "tmp_mndo.out" 
      integer, parameter :: mndo_in_unit = 998
      
      integer :: nqmatoms
      integer :: qmlist(maxatm)
      
      integer :: mndot_nline, mndot_oline, mndot_eline
      character(len=1024), allocatable :: mndo_template(:) 
      save
      end
