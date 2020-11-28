      module mndo
      use sizes

      implicit none
      character(len=1024), parameter ::
     & mndo_exe = "/home/mattia/work_tmp/MNDO2020/mndo2020/mndo2020",
     & mndo_in = "tmp_mndo.inp",
     & mndo_out = "tmp_mndo.out", 
     & mndo_intfi = "fort.15",
     & template_fname = 'template.inp'

      integer, parameter :: mndo_in_unit = 998, mndo_out_unit = 997,
     & temp_unit=999

      integer :: nqmatoms
      integer :: qmlist(maxatm), mmlist(maxatm)
      logical :: isqm(maxatm), ismndoinit
      
      integer :: mndo_nwk, mndo_neline
      character(len=128), allocatable :: mndo_keyword(:), mndo_eline(:)
     
      contains

      subroutine mndo_out_keys(unitio)
        implicit none

        integer :: unitio

        integer :: nchl, i, l

        integer trimtext

        nchl = 1
        do i=1, mndo_nwk
          l = trimtext(mndo_keyword(i))

c         Create a line every 80 char
          if(nchl + l .gt. 50 .and. nchl .gt. 1) then
            write(unitio, *) '+'
            nchl = 1
          end if
          
          write(unitio, "(A)", advance="no") mndo_keyword(i)(:l)
          write(unitio, "(A)", advance="no") " " 
          nchl = nchl + l +  1
        end do

        write(unitio, *) ""
      end

      subroutine mndo_parse_key(optline, noptline)
        implicit none

        character(len=1024) :: optline(:)
        integer :: noptline
        
        character(len=128) :: original_keyword(1024)
        integer :: ikey, i, j, l, kb, ke

        integer trimtext

        ikey = 1
        do i=1, noptline
          l = trimtext(optline(i))
          kb = 1
          ke = 1
          do j=1, l
            if(optline(i)(j:j) .eq. ' ' .or. j .eq. l) then
              ke = j - 1
              if(j .eq. l) ke = ke + 1
              if(optline(i)(ke:ke) .eq. '+') ke = ke - 1
              if(ke-kb .gt. 0) then
                original_keyword(ikey) = optline(i)(kb:ke)
                ikey = ikey + 1
              end if
              kb = ke + 2
            end if
          end do
        end do

c       A smart check should be done TODO
        mndo_nwk = ikey
        allocate(mndo_keyword(mndo_nwk))
        do i=1, ikey
          l = trimtext(original_keyword(i))
          mndo_keyword(i) = original_keyword(i)(:l)
        end do
      end
      
      end module mndo
