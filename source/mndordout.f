      subroutine mndordout(doforce)
        use atoms
        use energi
        use deriv
        use mndo

        implicit none

        logical :: doforce

        logical, parameter :: dosck = .true.

        integer :: i, j, isec, ios, if_natqm, if_natmm, xx, isla, ila
        logical :: intfexists, sck_passed
        character(len=1024) :: line
        real*8 :: cgn, ign, mygn, law3(3)
        real*8, allocatable :: demndo_la(:,:)

        integer trimtext

  10    format(I5)
  30    format(10x,3F20.10)
  40    format(2I5,3F20.10,I5)

        inquire(file=mndo_intfi, exist=intfexists)
        if (.not. intfexists) then
          write(6, *) "MNDO interface file was not found."
          write(6, *) "Interface file path ",
     &    mndo_intfi(:trimtext(mndo_intfi))
          call fatal
        end if

        if(mndo_usela) then
          allocate(demndo_la(3,mndo_nla))
        end if

        open(unit=mndo_out_unit, file=mndo_intfi)
       
        isec = 1
        do 
          read(mndo_out_unit, '(A)', iostat=ios) line
          if (ios .ne. 0) exit

          if (isec .eq. 1) then
            read(line(32:), 10) if_natqm
            isec = 2
          else if(isec .eq. 2) then
            if( trimtext(line) .eq. 0) then
              isec = 3
              if(n-nqmatoms .eq. 0) isec = 5
            end if
          else if(isec .eq. 3) then
            read(line(58:), 10) if_natmm
            isec = 4
          else if(isec .eq. 4) then
            if( trimtext(line) .eq. 0) isec = 5
          else if(isec .eq. 5) then
            isec = 6
          else if(isec .eq. 6) then
            if(trimtext(line) .eq. 0) then
              if(doforce) then
                isec = 7
              else
                isec = 11
              end if
            else
              read(line, 30) emndo, cgn, ign
            end if
          else if(isec .eq. 7) then
c         Header of QM gradients
            isec = 8
            j = 1
            ila = 1
          else if(isec .eq. 8) then
c         QM gradients values
            if(trimtext(line) .eq. 0) then
              isec = 9
              if(n-nqmatoms .eq. 0) isec = 11
            else
              if(j .le. nqmatoms) then
                read(line, 40) xx, xx, demndo(1, qmlist(j)), 
     $          demndo(2, qmlist(j)), demndo(3, qmlist(j)), isla
                if(isla .ne. 0) then
                  write(6, *) "Reading a QM atom gradient as a",
     &            " link atom's one. This is a bug."
                  call fatal
                end if
              else if(j .le. nqmatoms + mndo_nla .and. mndo_usela) then
                read(line, 40) xx, xx, demndo_la(1, ila), 
     $          demndo_la(2, ila), demndo_la(3, ila), isla
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
            isec = 10
            j = 1
          else if(isec .eq. 10) then
            if(trimtext(line) .eq. 0) then
              isec = 11
            else
              read(line, 40) xx, xx, demndo(1, mmlist(j)), 
     $        demndo(2, mmlist(j)), demndo(3, mmlist(j)), xx 
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

          if(doforce) then
c           Check if the norm of QM atoms' gradiend is equal to the one in
c           output file
            mygn = 0.0
            do i=1, nqmatoms 
              do j=1, 3
                mygn = mygn + demndo(j,qmlist(i))**2
              end do
            end do

            if(mndo_usela) then
              do i=1, mndo_nla
                do j=1, 3
                  mygn = mygn + demndo_la(j, i)**2
                end do
              end do
            end if
            mygn = sqrt(mygn)

            if(abs(mygn - cgn) .gt. 1.0e-7) then
              sck_passed = .false.
              write(6, *)"Difference between compute and expected norm",
     $        " of QM atoms' gradient > 1.0e-7"
              write(6, "('Found= ', F12.6, '  Computed= ', F12.6)") cgn,
     $        mygn
            end if
          end if
        end if

c       Project LA gradients on QM and MM atoms
        call mndo_laproj(demndo_la)

        if(dosck .and. doforce) then
c         Check if the computed gradients are OK with Newton 3th law
          law3 = 0.0
          do i=1, n
            do j=1, 3
              law3(j) = law3(j) + demndo(j, i)
            end do
          end do
          
          if(norm2(law3) .gt. mndo_l3t) then
            sck_passed = .false.
            write(6, *) "Computed forces does not respect Newton 3th", 
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

      end
