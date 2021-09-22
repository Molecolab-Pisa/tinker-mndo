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
     & temp_unit=999, mndo_maxs=5
      logical, parameter :: mndo_debug = .false., mndo_la13 = .true.
      
      real*8, parameter :: distqmla = 1.10, mndo_l3t=1.0e-3

      character(len=1024) :: mndo_exe, template_fname, mndo_postexe
      logical :: mndo_iterguess


      integer :: nqmatoms, mndo_nla, mndo_nconjat
      integer :: qmlist(maxatm), mmlist(maxatm), mndo_laqm(maxatm),
     &           mndo_lamm(maxatm), mndo_conjlist(maxatm), 
     &           mndo_qmconjl(maxatm)
      logical :: isqm(maxatm), ismndoinit, mndo_dope, mndo_usela, 
     &           mndo_delchg(maxatm) 

      integer :: mndo_nstates, mndo_currentstate, 
     &           mndo_kci, mndo_icross
      logical :: mndo_multistate, mndo_ms_all
      
      integer :: mndo_nwk, mndo_neline
      character(len=128), allocatable :: mndo_keyword(:), mndo_eline(:)

      real*8, allocatable :: demndo_tmp(:,:,:), emndo_tmp(:)
     
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
          
          if( mndo_keyword(i)(:5) .eq. 'imult' ) nchl = 51
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

        integer, parameter :: nauto = 8, nskip = 1, nitgu = 3, nla = 2
        character(len=128) :: automatic_kwd(nauto) = (/
c       1         2         3         4         5         6
     &  "iform ", "mminp ", "numatm", "mmcoup", "mmskip", "nsav15",
c       7         8         9         10        11        12
     &  "igeom ", "nconj "
     &  /)
        integer :: automatic_prm(nauto) = (/ 1, 2, -1, 2, 1, 3, 1, 0/)

        character(len=128) :: skip_kwd(nskip) = (/
c       1         2         3         4         5         6
     &  "jop   "
     &  /)

        character(len=128) :: itgu_kwd(nitgu) = (/
c       1         2         3         4         5         6
     &  "ipubo ", "ktrial", "imomap"
     &  /)
        integer :: itgu_prm(nitgu) = (/ 1, 11, 3 /)
        
        character(len=128) :: la_kwd(nla) = (/
c       1         2         3         4         5         6
     &  "mmlink", "nlink "
     &  /)
        integer :: la_prm(nla) = (/ 2, -1 /)
        
        integer :: i, j, k, l, prm
        character(len=128) :: rch, kw
        logical :: toadd

        integer trimtext

        automatic_prm(3) = n - nqmatoms
        if (n - nqmatoms .eq. 0) then
          automatic_prm(2) = 0
          automatic_prm(4) = 0
        end if
        if(mndo_usela) la_prm(2) = mndo_nla
        
        if(mndo_nconjat .gt. 0) then
          automatic_prm(8) = mndo_nconjat
        end if
        
        mndo_kci = 0
        mndo_icross = 0

c       Check and remove unneeded automatic and skip keyword

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
         
          if(kw(:8) .eq. 'kci   ') mndo_kci = prm
          if(kw(:8) .eq. 'icross') mndo_icross = prm

c         Check if it is not in an automatic keyword
          toadd = .true.
          do j=1, nauto
            if(kw(:8) .eq. automatic_kwd(j)) then
              write(6, *) "Keyword ", kw(:l), " is handled by Tinker-",
     &        "MNDO interface. The value found in template will be ",
     &        "ignored."
              toadd = .false.
              exit
            endif
          end do
          
          do j=1, nitgu
            if(kw(:8) .eq. itgu_kwd(j)) then
              write(6, *) "Keyword ", kw(:l), " is handled by Tinker-",
     &        "MNDO interface. The value found in template will be ",
     &        "ignored."
              toadd = .false.
              exit
            endif
          end do

          if(mndo_usela) then
            do j=1, nla
              if(kw(:8) .eq. la_kwd(j)) then
                write(6, *) "Keyword ",kw(:l), " is handled by Tinker-",
     &          "MNDO interface. The value found in template will be ",
     &          "ignored."
                toadd = .false.
                exit
              endif
            end do
          end if

          do j=1, nskip
            if(kw(:8) .eq. skip_kwd(j)) then
              write(6, *) "Keyword ", kw(:l), " is handled by Tinker-",
     &        "MNDO interface. The value found in template will be ",
     &        "ignored."
              toadd = .false.
              exit
            endif
          end do

          if(toadd) then
            write(rch, *) prm
            rch = adjustl(rch)
            write(key_buffer(mndo_nwk), "(A,'=',A)") 
     &      kw(:l),
     &      rch(:trimtext(rch))     
            mndo_nwk = mndo_nwk + 1
          end if
        end do

c       Now add each automatic keyword
        do i=1, nauto
          write(rch, *) automatic_prm(i)
          rch = adjustl(rch)
          write(key_buffer(mndo_nwk), "(A,'=',A)") 
     &    automatic_kwd(i)(:trimtext(automatic_kwd(i))),
     &    rch(:trimtext(rch))     
          mndo_nwk = mndo_nwk + 1
        end do
        
        if(mndo_iterguess) then
          do i=1, nitgu
            write(rch, *) itgu_prm(i)
            rch = adjustl(rch)
            write(key_buffer(mndo_nwk), "(A,'=',A)") 
     &      itgu_kwd(i)(:trimtext(itgu_kwd(i))),
     &      rch(:trimtext(rch))     
            mndo_nwk = mndo_nwk + 1
          end do
        end if

        if(mndo_usela) then
          do i=1, nla
            write(rch, *) la_prm(i)
            rch = adjustl(rch)
            write(key_buffer(mndo_nwk), "(A,'=',A)") 
     &      la_kwd(i)(:trimtext(la_kwd(i))),
     &      rch(:trimtext(rch))     
            mndo_nwk = mndo_nwk + 1
          end do
        end if

        mndo_nwk = mndo_nwk - 1
        allocate(mndo_keyword(mndo_nwk + nskip))
        do i=1, mndo_nwk
          mndo_keyword(i) = key_buffer(i)
          if(mndo_debug) 
     &      write(6, '("(",I3,")  ", A)')
     &      i, mndo_keyword(i)(:trimtext(mndo_keyword(i)))      
        end do
      end

      subroutine mndomskey(new_icross) 
        implicit none

        integer :: new_icross

        character(len=128) :: key_buffer(1024)

        integer, parameter :: nauto = 4, nskip = 1
        character(len=128) :: automatic_kwd(nauto) = (/
c       1         2         3         4         5         6
     &  "icross", "ncigrd", "iroot ", "lroot "
     &  /)
        integer :: automatic_prm(nauto) = (/ 0, 0, 0, 1 /)

        integer :: i, j, k, l, nwk, prm
        character(len=128) :: rch, kw
        logical :: toadd

        integer trimtext

        automatic_prm(1) = new_icross
        if(mndo_ms_all) then
          automatic_prm(2) = mndo_nstates
        else
          automatic_prm(2) = 1
        end if
        automatic_prm(3) = mndo_nstates
        
        if(mndo_currentstate .gt. 0) then
          automatic_prm(4) = mndo_currentstate
        end if

        nwk = 0
        do i=1, mndo_nwk
          l = trimtext(mndo_keyword(i))

          do j=1, l
            if(mndo_keyword(i)(j:j) .eq. '=') exit
          end do

          kw = mndo_keyword(i)(:j-1)
          kw(j:) = ' '
          read(mndo_keyword(i)(j+1:l),'(I10)') prm
          l = trimtext(kw)
         
c         Check if it is not in an automatic keyword
          toadd = .true.
          do j=1, nauto
            if(kw(:8) .eq. automatic_kwd(j)) then
              write(6, *) "Keyword ", kw(:l), " is handled by Tinker-",
     &        "MNDO interface. The value found in template will be ",
     &        "ignored."
              toadd = .false.
              exit
            endif
          end do
          
          if(toadd) then
            nwk = nwk + 1
            write(rch, *) prm
            rch = adjustl(rch)
            write(key_buffer(nwk), "(A,'=',A)") 
     &      kw(:l),
     &      rch(:trimtext(rch))     
          end if
        end do

c       Now add each automatic keyword
        do i=1, nauto
          nwk = nwk + 1
          write(rch, *) automatic_prm(i)
          rch = adjustl(rch)
          write(key_buffer(nwk), "(A,'=',A)") 
     &    automatic_kwd(i)(:trimtext(automatic_kwd(i))),
     &    rch(:trimtext(rch))     
        end do
        
        mndo_nwk = nwk
        deallocate(mndo_keyword)
        allocate(mndo_keyword(mndo_nwk + nskip))
        do i=1, mndo_nwk
          mndo_keyword(i) = key_buffer(i)
          if(mndo_debug) 
     &      write(6, '("(",I3,")  ", A)')
     &      i, mndo_keyword(i)(:trimtext(mndo_keyword(i)))      
        end do
      end
  
      subroutine mndo_lapos(ila, pos)
        use atoms

        implicit none

        integer :: ila
        real*8 :: pos(3)

        integer :: imm, iqm 
        real*8 :: mmp(3), qmp(3)

        imm = mndo_lamm(ila)
        iqm = mndo_laqm(ila)

        qmp(1) = x(iqm)
        qmp(2) = y(iqm)
        qmp(3) = z(iqm)

        mmp(1) = x(imm)
        mmp(2) = y(imm)
        mmp(3) = z(imm)
        
        call mndo_lapos_def(mmp, qmp, pos)
      end
      
      subroutine mndo_lapos_def(mmp, qmp, pos)
        implicit none
        real*8 :: pos(3), mmp(3), qmp(3), dp(3)
        
        dp = mmp - qmp
        pos = qmp + distqmla/norm2(dp)*dp
      end
      
      subroutine mndo_laproj_tensor(qmp, mmp, tensor, qm)
        implicit none

        real*8 :: qmp(3), mmp(3), tensor(3,3), delta, p0(3), p1(3)
        integer :: i, j
        logical :: qm

        delta = 1e-5

        do i=1, 3
          if(qm) then
            qmp(i) = qmp(i) + delta
          else
            mmp(i) = mmp(i) + delta
          end if

          call mndo_lapos_def(mmp, qmp, p0)

          if(qm) then
            qmp(i) = qmp(i) - 2*delta
          else
            mmp(i) = mmp(i) - 2*delta
          end if

          call mndo_lapos_def(mmp, qmp, p1)

          tensor(:,i) = (p0 - p1)/(2*delta)
          
          if(qm) then
            qmp(i) = qmp(i) + delta
          else
            mmp(i) = mmp(i) + delta
          end if
        end do
      end 
      
      subroutine num_mndo_laproj(laforces, forces)
        use atoms

        implicit none

        real*8 :: laforces(3,mndo_nla), forces(3,n)

        integer :: ila, imm, iqm, i, j
        real*8 :: mmp(3), qmp(3), tensor(3,3)

        if(mndo_debug) write(6, *) "Entering mndo_laproj_numerical"

        do ila=1, mndo_nla
          imm = mndo_lamm(ila)
          iqm = mndo_laqm(ila)

          qmp(1) = x(iqm)
          qmp(2) = y(iqm)
          qmp(3) = z(iqm)

          mmp(1) = x(imm)
          mmp(2) = y(imm)
          mmp(3) = z(imm)
          
          call mndo_laproj_tensor(qmp, mmp, tensor, .true.)
          
          do i=1, 3
            do j=1, 3
              forces(i,iqm) = forces(i,iqm) +
     &                        laforces(j,ila)*tensor(j,i)
            end do
          end do
          
          call mndo_laproj_tensor(qmp, mmp, tensor, .false.)
          
          do i=1, 3
            do j=1, 3
              forces(i,imm) = forces(i,imm) +
     &                        laforces(j,ila)*tensor(j,i)
            end do
          end do

        end do
        
        if(mndo_debug) write(6, *) "Exiting mndo_laproj_numerical"
      end

      subroutine mndo_laproj(laforces, forces)
        use atoms

        implicit none

        real*8 :: laforces(3,mndo_nla), forces(3,n)

        integer :: ila, imm, iqm, i, j
        real*8 :: mmp(3), qmp(3), g, np(3), proj

        if(mndo_debug) write(6, *) "Entering mndo_laproj"

        do ila=1, mndo_nla
          imm = mndo_lamm(ila)
          iqm = mndo_laqm(ila)

          qmp(1) = x(iqm)
          qmp(2) = y(iqm)
          qmp(3) = z(iqm)

          mmp(1) = x(imm)
          mmp(2) = y(imm)
          mmp(3) = z(imm)
          
          np = (mmp - qmp) / norm2(mmp - qmp)
          proj = sum(laforces(:,ila) * np)
          g = distqmla/norm2(mmp - qmp)

          forces(:,iqm) = forces(:,iqm) + g*proj*np +
     &                    (1.0-g)*laforces(:,ila)
          forces(:,imm) = forces(:,imm) - g*proj*np + g*laforces(:,ila)
        end do
        
        if(mndo_debug) write(6, *) "Exiting mndo_laproj"
      end

      end module mndo
