      subroutine mndomkin()
        use atomid
        use atoms
        use mndo
        use charge
        use usage

        implicit none

        integer :: i, j, l, elb, iselct
        logical :: delchg
        real*8 :: lapos(3)

        integer trimtext
  
  10   format(i2, 8x, 3(F10.5,2x,i2,6x))
  20   format(3F12.4,F8.4,I3)

        if(mndo_debug) write(6, *) "Entering subroutine mndomkin"

        open(unit=mndo_in_unit, file=mndo_in)

c       Write input options
        call mndo_out_keys(mndo_in_unit)
        if(mndo_debug) then
          write(6, *) "INPUT SECTION"
          call mndo_out_keys(6)
          write(6, *)
        end if
        
c       Write header comments
        write(mndo_in_unit, *) "Generated by Tinker-MNDO interface."
        write(mndo_in_unit, *) "Finger crossed."

        do i=1, nqmatoms
          write(mndo_in_unit, 10) atomic(qmlist(i)), x(qmlist(i)), 1,
     &    y(qmlist(i)), 1, z(qmlist(i)), 1
        end do

        if(mndo_usela) then
          do i=1, mndo_nla
c           Compute coordinates of LAs (fixed distance)
            call mndo_lapos(i, lapos)

c           write down coordinates of LAs
            write(mndo_in_unit, 10) 1, lapos(1), 1, lapos(2), 1, 
     &      lapos(3), 1
          end do
        end if

        write(mndo_in_unit, *) ""
     
        do i=1, mndo_neline
          l = trimtext(mndo_eline(i))
          write(mndo_in_unit, *) mndo_eline(i)(:l)
        end do

        if(mndo_nconjat > 0) then
c         write down conjugate atoms list format 20i4
          do i=1, mndo_nconjat
            if(mod(i,20) .eq. 0) write(mndo_in_unit, *) ""
            write(mndo_in_unit, "(I4)", advance="no") mndo_qmconjl(i)
          end do
          write(mndo_in_unit, *) ""
        end if
        
        if(mndo_multistate) then
          if(mndo_ms_all) then
            do i=1, mndo_nstates
              write(mndo_in_unit, "(I4)", advance="no") i
            end do
            write(mndo_in_unit, *) ""
          else
              write(mndo_in_unit, "(I4)") mndo_currentstate
          end if
        end if

        if(mndo_usela) then
c         write down link atoms list
          do i=1, mndo_nla
            if(mod(i,16) .eq. 0) write(mndo_in_unit, *) ""
            write(mndo_in_unit, "(I5)", advance="no") nqmatoms + i
          end do
          write(mndo_in_unit, *) ""
        end if

        do i=1, n-nqmatoms
          if( use(mmlist(i)) ) then
            iselct = 0
          else
            iselct = 1
          end if
          
          delchg = .false.
          if( mndo_usela ) then
            do j=1, mndo_nla
              if(mndo_lamm(j) .eq. mmlist(i)) delchg = .true.
            end do
          end if
          
          if(delchg) then
c           This charge is just here as place holder!
            if(mndo_debug) write(*, *) "Charge on atom ", mmlist(i), 
     $         " not included in MNDO."
            write(mndo_in_unit, 20) x(mmlist(i)), y(mmlist(i)), 
     $      z(mmlist(i)), 0.00, iselct
          else
            write(mndo_in_unit, 20) x(mmlist(i)), y(mmlist(i)), 
     $      z(mmlist(i)), pchg(mmlist(i)), iselct
          end if
        end do

        close(mndo_in_unit)

        if(mndo_debug) write(6, *) "Exiting subroutine mndomkin"
      end
