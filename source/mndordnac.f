      subroutine mndordnac()
        use atoms
        use mndo

        implicit none

        integer :: i, j, isec, ios, if_natqm, if_natmm, xx, isla, ila,
     &             mndo_nnac
        
        character(len=1024) :: line
        character(len=1024), parameter :: 
     & QM_NAD_H = 'CARTESIAN INTERSTATE COUPLING GRADIENT FOR STATES ', 
     & MM_NAD_H = 
     & 'CARTESIAN INTERSTATE COUPLING GRADIENT OF MM ATOMS FOR STATES '
        
        real*8 :: law3(3), yy
        real*8, allocatable :: nac(:,:,:), nac_la(:,:,:)

        integer trimtext

  10    format(I5)
  30    format(I5,5x,3F20.10)
  40    format(2I5,3F20.10,I5)
        
        if(mndo_debug) write(6, *) "Entering subroutine mndordnac"

        inquire(file=mndo_intfi, exist=intfexists)
        if (.not. intfexists) then
          write(6, *) "MNDO interface file was not found."
          write(6, *) "Interface file path ",
     &    mndo_intfi(:trimtext(mndo_intfi))
          call fatal
        end if
        
        mndo_nnac = mndo_nstates * (mndo_nstates-1) / 2
       
        allocate(nac(3,n,mndo_nnac)
        if(mndo_usela) then
          allocate(demndo_la(3,mndo_nla,mndo_nstates))
        end if

        open(unit=mndo_out_unit, file=mndo_intfi)
       
        isec = 0
        do 
          read(mndo_out_unit, '(A)', iostat=ios) line
          if (ios .ne. 0) exit

          else if(isec .eq. 8) then
c         QM gradients values
            if(trimtext(line) .eq. 0) then
              if(n-nqmatoms .eq. 0) then
c               Case full QM
                if(istate .eq. if_nstate) then
                  isec = 11
                else
                  isec = 7
                end if
              else
                isec = 9
              end if
            else
              if(j .le. nqmatoms) then
                read(line, 40) xx, xx, 
     $          demndo_tmp(1, qmlist(j),if_istate), 
     $          demndo_tmp(2, qmlist(j),if_istate), 
     $          demndo_tmp(3, qmlist(j),if_istate), isla
                if(isla .ne. 0) then
                  write(6, *) "Reading a QM atom gradient as a",
     &            " link atom's one. This is a bug."
                  call fatal
                end if
              else if(j .le. nqmatoms + mndo_nla .and. mndo_usela) then
                read(line, 40) xx, xx, demndo_la(1,ila,if_istate), 
     $          demndo_la(2,ila,if_istate), 
     $          demndo_la(3,ila,if_istate),isla
                ila = ila + 1
                if(isla .ne. 1) then
                  write(6, *) "Reading a link atom gradient as a",
     &            " regular QM atom's one. This is a bug."
                  call fatal
                end if
              else
                write(6, *) "Wrong number of QM atoms' gradient in ",
     &          mndo_intfi(:trimtext(mndo_intfi))
                call fatal
              end if
              j = j + 1
            end if
          else if(isec .eq. 9) then
c           MM gradients header
            isec = 10
            read(line(42:), *) if_istate 
            if(if_states(istate) .ne. if_istate) then
              write(6, *) "Error while reading MM gradient section.",
     $        "Expected gradient for a different state."
              call fatal
            end if
            j = 1
          else if(isec .eq. 10) then
            if(trimtext(line) .eq. 0) then
              if(istate.eq.if_nstate) then
                isec = 11
              else
                isec = 7
              end if
            else
              read(line, 40) xx, xx, demndo_tmp(1, mmlist(j),if_istate), 
     $        demndo_tmp(2, mmlist(j),if_istate), 
     $        demndo_tmp(3, mmlist(j),if_istate), xx 
              j = j + 1
            end if
          else
            continue
          end if
        end do
        close(unit=mndo_out_unit)

c     Insert some sanity check here
        if(dosck) then
          sck_passed = .true.
c         Check if you found the correct number of atoms
          if(if_natqm .ne. nqmatoms + mndo_nla) then
            sck_passed = .false.
            write(6, *) "Wrong number of QM atoms in ", 
     &      mndo_intfi(:trimtext(mndo_intfi))
          end if
          if(if_natmm + nqmatoms .ne. n .and. n .ne. nqmatoms) then
            sck_passed = .false. 
            write(6, *) "Wrong number of MM atoms in ", 
     &      mndo_intfi(:trimtext(mndo_intfi))
          end if

c         Check if the norm of QM atoms' gradiend is equal to the one in
c         output file
          do istate=1, if_nstate
            mygn = 0.0
            do i=1, nqmatoms 
              do j=1, 3
              mygn = mygn + demndo_tmp(j,qmlist(i),if_states(istate))**2
              end do
            end do

            if(mndo_usela) then
              do i=1, mndo_nla
                do j=1, 3
                mygn = mygn + demndo_la(j, i,if_states(istate))**2
                end do
              end do
            end if
            mygn = sqrt(mygn)
            
            if(abs(mygn - cgn(if_states(istate))) .gt. 1.0e-7) then
              sck_passed = .false.
              write(6, *)"Difference between compute and expected norm",
     $        " of QM atoms' gradient > 1.0e-7; State = ",
     $        if_states(istate)
              write(6, "('Found= ', F12.6, '  Computed= ', F12.6)") 
     $        cgn(if_states(istate)), mygn
            end if
          end do
        end if

        demndo = demndo_tmp(:,:,mndo_currentstate)
        emndo = emndo_tmp(mndo_currentstate)

c       Project LA gradients on QM and MM atoms
        call mndo_laproj(demndo_la(:,:,mndo_currentstate))

        if(dosck .and. count(use) .eq. n) then
c         Check if the computed gradients are OK with Newton 3rd law
          law3 = 0.0
          do i=1, n
            do j=1, 3
              law3(j) = law3(j) + demndo(j, i)
            end do
          end do
          
          if(norm2(law3) .gt. mndo_l3t) then
            sck_passed = .false.
            write(6, *) "Computed forces does not respect Newton 3rd", 
     &      " law."
            write(6, '("Ftot(X,Y,Z) = ", 3F12.6)') law3
          end if
        end if
        
        if(dosck) then
          if(.not. sck_passed) then
            write(6, *) "MNDO sanity check not passed!"
            write(6, *) "This is probably due to wrong input or bugs."
            call fatal
          end if
        end if
        
        if(mndo_debug) write(6, *) "Exiting subroutine mndordmsout"

      end
