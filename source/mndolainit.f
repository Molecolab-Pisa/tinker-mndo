      subroutine mndolainit()
        use atoms
        use couple
        use mndo

        implicit none

        integer :: i, j, k, iatqm

c       Find out if a link atom is needed
        mndo_usela = .false.
        mndo_nla = 0
        do i=1, nqmatoms
          iatqm = qmlist(i)
          do j=1, n12(iatqm)
            do k=1, n-nqmatoms
              if(i12(j,iatqm) .eq. mmlist(k)) then
                mndo_usela = .true.
                mndo_nla = mndo_nla + 1
                mndo_laqm(mndo_nla) = iatqm
                mndo_lamm(mndo_nla) = i12(j,iatqm)
              end if
            end do
          end do
        end do

c      Find out which atoms should have charge removed
       mndo_delchg(:) = .false.
       do i=1, mndo_nla
        mndo_delchg(mndo_lamm(i)) = .true.
        if(mndo_la13) then
          do j=1, n12(mndo_lamm(i))
            mndo_delchg(i12(j,mndo_lamm(i))) = .true.
          end do
        end if
       end do
       
        if(mndo_debug .and. mndo_nla .gt. 0) then
          write(6, *) ""
          write(6, *) "Link atom table"
          write(6, *) "|  QM  |  MM  |"
          write(6, *) "---------------"
          do i=1, mndo_nla
            write(6, "(' | ', I4,' | ', I4,' |')") mndo_laqm(i),
     &       mndo_lamm(i)
          end do
          write(6, *) "---------------"
          write(6, *) ""
          write(6, *) "Charge will be removed on the following centers:"
          do i=1, n
            if(mndo_delchg(i)) write(6, "('ATOM : ', I5)") i
          end do
        end if

      end  
