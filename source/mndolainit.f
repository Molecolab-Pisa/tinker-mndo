      subroutine mndolainit()
        use atoms
        use couple
        use mndo

        implicit none

        integer :: i, j, k, iatqm

c       Find out if a link atom is needed
        write(6,*) 'SCHIANTA FORSE QUA?'
        flush(6)
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
        end if

        if(mndo_usela) then
          write(6, *) "Link atom are still not implemented."
c          call fatal
        end if
      end  
