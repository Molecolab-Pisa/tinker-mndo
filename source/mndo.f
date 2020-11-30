      module mndo
      use sizes

      implicit none
      character(len=1024), parameter ::
     & mndo_in = "tmp_mndo.inp",
     & mndo_out = "tmp_mndo.out", 
     & mndo_intfi = "fort.15",
     & mndo_default_exe = "mndo2020",
     & mndo_default_template = 'template.inp' 
      integer, parameter :: mndo_in_unit = 998, mndo_out_unit = 997,
     & temp_unit=999
      logical, parameter :: iter_guess = .true., mndo_debug = .true.
     
      character(len=1024) :: mndo_exe, template_fname, mndo_postexe


      integer :: nqmatoms
      integer :: qmlist(maxatm), mmlist(maxatm)
      logical :: isqm(maxatm), ismndoinit, mndo_dope
      
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

c         Create a line every 50 char
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

c       A smart check should be done
        call mndo_check_keys(original_keyword, ikey)
      end

      subroutine mndo_check_keys(orig_kw, nkw)
        use atoms

        implicit none

        character(len=128) :: orig_kw(:)
        integer :: nkw

        character(len=128) :: key_buffer(1024)

        integer, parameter :: nauto = 8
        character(len=128) :: automatic_kwd(nauto) = (/
c       1         2         3         4         5         6
     &  "iform ", "jop   ", "mminp ", "numatm", "mmcoup", "mmskip",
c       7         8
     &  "nsav15", "igeom "
     &  /)
        integer :: automatic_prm(nauto) = (/ 1, -2, 2, -1, 2, 0, 3, 1/)
        
        integer :: i, j, k, l, prm
        character(len=128) :: rch, kw
        logical :: isauto

        integer trimtext

        automatic_prm(4) = n - nqmatoms
        if (n - nqmatoms .eq. 0) then
          automatic_prm(3) = 0
          automatic_prm(5) = 0
        end if

c       Check and remove unneeded automatic keyword

        mndo_nwk = 1
        do i=1, nkw
          l = trimtext(orig_kw(i))
          if(l.eq.0) cycle

          do j=1, l
            if(orig_kw(i)(j:j) .eq. '=') exit
          end do

          if(j.eq.l+1) then
            write(6, *) "A keyword provided in template file (",
     &      orig_kw(i)(:l), ") does not contain '='."
            write(6, *) "This is not allowed, please only use MNDO ",
     &      "standard input format."
            call fatal
          end if

          kw = orig_kw(i)(:j-1)
          kw(j:) = ' '
          read(orig_kw(i)(j+1:l),'(I10)') prm
          l = trimtext(kw)
          
c         Check if it is not in an automatic keyword
          isauto = .false.
          do j=1, nauto
            if(kw(:8) .eq. automatic_kwd(j)) then
              write(6, *) "Keyword ", kw(:l), " is handled by Tinker-",
     &        "MNDO interface. The value found in template will be ",
     &        "ignored."
              isauto = .true.
              exit
            endif
          end do

          if(.not. isauto) then
            write(rch, *) prm
            rch = adjustl(rch)
            write(key_buffer(mndo_nwk), "(A,'=',A)") 
     &      kw(:l),
     &      rch(:trimtext(rch))     
            mndo_nwk = mndo_nwk + 1
          end if
        end do
        
        do i=1, nauto
          write(rch, *) automatic_prm(i)
          rch = adjustl(rch)
          write(key_buffer(mndo_nwk), "(A,'=',A)") 
     &    automatic_kwd(i)(:trimtext(automatic_kwd(i))),
     &    rch(:trimtext(rch))     
          mndo_nwk = mndo_nwk + 1
        end do

        allocate(mndo_keyword(mndo_nwk))
        do i=1, mndo_nwk
          mndo_keyword(i) = key_buffer(i)
        end do
      end

      end module mndo
