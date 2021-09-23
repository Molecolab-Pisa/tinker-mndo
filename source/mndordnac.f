      subroutine mndordnac()
        use atoms
        use mndo

        implicit none

        integer :: i, j, isec, ios, mndo_nnac, nac_st_a, nac_st_b, nc,
     &            ila, isla, xx, if_natqm, if_natmm, iintf

        logical :: intfexists, dosck, sck_passed
        
        character(len=1024) :: line, nac_intffile
        character(len=1024), parameter :: 
     & QM_NAD_H = 'CARTESIAN INTERSTATE COUPLING GRADIENT FOR STATES ', 
     & MM_NAD_H = 
     & 'CARTESIAN INTERSTATE COUPLING GRADIENT OF MM ATOMS FOR STATES ',
     & N_NAD_H = 'CARTESIAN AND INTERNAL INTERSTATE COUPLING NORMS'
        
        real*8 :: law3(3), yy, tmp_cart, tmp_int, mygn
        real*8, allocatable :: nac(:,:,:), nac_la(:,:,:), nac_norm(:)

        integer trimtext, freeunit

  40    format(2I5,3F20.10,I5)
  50    format(2I5,20X,2F20.10) 
  60    format("CARTESIAN INTERSTATE COUPLING GRADIENT FOR STATES",2I5)
  70    format(I5, 3F20.10)

        if(mndo_debug) write(6, *) "Entering subroutine mndordnac"
        
        dosck = .true.
        
        inquire(file=mndo_intfi, exist=intfexists)
        if (.not. intfexists) then
          write(6, *) "MNDO interface file was not found."
          write(6, *) "Interface file path ",
     &    mndo_intfi(:trimtext(mndo_intfi))
          call fatal
        end if
        
        mndo_nnac = mndo_nstates * (mndo_nstates-1) / 2
       
        allocate(nac(3,n,mndo_nnac))
        allocate(nac_norm(mndo_nnac))
        if(mndo_usela) then
          allocate(nac_la(3,mndo_nla,mndo_nnac))
        end if

c        if(mndo_debug) write(*, *) "Reading fort.15"
        open(unit=mndo_out_unit, file=mndo_intfi)
       
        isec = 0
        do 
          read(mndo_out_unit, '(A)', iostat=ios) line
          if (ios .ne. 0) exit
          
          if (isec .eq. 0 .and. 
     &        line(2:trimtext(QM_NAD_H)+2) .eq. QM_NAD_H) then
            read(line(trimtext(QM_NAD_H)+2:), *) nac_st_a, nac_st_b
            nc = (nac_st_a - 2)*(nac_st_a - 1) / 2 + nac_st_b 
c            if(mndo_debug) write(*, *) "Reading QM", nac_st_a, nac_st_b,
c     &        nc
            isec = 1
            j = 1
            ila = 1
          else if(isec .eq. 0 .and.
     &            line(2:trimtext(MM_NAD_H)+2) .eq. MM_NAD_H) then
            read(line(trimtext(MM_NAD_H)+2:), *) nac_st_a, nac_st_b
            nc = (nac_st_a - 2)*(nac_st_a - 1) / 2 + nac_st_b 
c            if(mndo_debug) write(*, *) "Reading MM", nac_st_a, nac_st_b,
c     &        nc
            isec = 2
            j = 1
          else if(isec .eq. 0 .and.
     &            line(2:trimtext(N_NAD_H)+2) .eq. N_NAD_H) then 
            isec = 3

          else if(isec .eq. 1) then
c         QM NAC values
            if(trimtext(line) .eq. 0) then
c             End of QM atoms
c              if(mndo_debug) write(*, *) "Done"
              if_natqm = j - 1
              isec = 0
            else
              if(j .le. nqmatoms) then
                read(line, 40) xx, xx, 
     $          nac(1, qmlist(j), nc), 
     $          nac(2, qmlist(j), nc), 
     $          nac(3, qmlist(j), nc), isla
                if(isla .ne. 0) then
                  write(6, *) "Reading a QM atom gradient as a",
     &            " link atom's one. This is a bug."
                  call fatal
                end if
              else if(j .le. nqmatoms + mndo_nla .and. mndo_usela) then
                read(line, 40) xx, xx, nac_la(1,ila,nc), 
     $          nac_la(2,ila,nc), 
     $          nac_la(3,ila,nc),isla
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
          else if(isec .eq. 2) then
            if(trimtext(line) .eq. 0) then
c           end of MM atoms
c              if(mndo_debug) write(*, *) "Done"
              if_natmm = j - 1
              isec = 0
            else
              read(line, 40) xx, xx, nac(1, mmlist(j),nc), 
     $        nac(2, mmlist(j),nc), 
     $        nac(3, mmlist(j),nc), xx 
              j = j + 1
            end if
          else if(isec .eq. 3) then
c         NAC norms values
            if(trimtext(line) .eq. 0) then
c             End of QM atoms
              isec = 0
            else
              read(line, 50) nac_st_a, nac_st_b, tmp_cart, tmp_int
              nc = (nac_st_a - 2)*(nac_st_a - 1) / 2 + nac_st_b 
              nac_norm(nc) = tmp_cart
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
          do nc=1, mndo_nnac
            mygn = 0.0
            do i=1, nqmatoms 
              do j=1, 3
              mygn = mygn + nac(j, qmlist(i), nc)**2
              end do
            end do

            if(mndo_usela) then
              do i=1, mndo_nla
                do j=1, 3
                mygn = mygn + nac_la(j, i, nc)**2
                end do
              end do
            end if
            mygn = sqrt(mygn)
            
            if(abs(mygn - nac_norm(nc)) .gt. 1.0e-7) then
              sck_passed = .false.
              write(6, *)"Difference between compute and expected norm",
     $        " of QM atoms' NAC > 1.0e-7; NC = ",
     $        nc
              write(6, "('Found= ', F12.6, '  Computed= ', F12.6)") 
     $        nac_norm(nc), mygn
            end if
          end do
        end if

c       Project LA gradients on QM and MM atoms
        if(mndo_usela) then
          do i=1, mndo_nnac
            call mndo_laproj(nac_la(:,:,i),
     &                       nac(:,:,i))
          end do
        end if

c       Write down the nac in the xyz order and with LA projected out!
        nac_intffile = 'interface_nac.dat'
        iintf = freeunit()

        open(unit=iintf, file=nac_intffile)
        
        do nac_st_b=1, mndo_nstates-1
          do nac_st_a=nac_st_b+1, mndo_nstates
            nc = (nac_st_a - 2)*(nac_st_a - 1) / 2 + nac_st_b
c           write header
            write(iintf, 60) nac_st_a, nac_st_b
            do i=1, n
              write(iintf, 70) i, nac(:,i,nc)
            end do
            write(iintf, *) ''
          end do
        end do

        close(unit=iintf)
        
        deallocate(nac)
        if(mndo_usela) deallocate(nac_la)

        if(mndo_debug) write(6, *) "Exiting subroutine mndordnac"

      end
