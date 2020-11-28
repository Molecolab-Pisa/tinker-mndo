      subroutine mndoinit()
        use atoms
        use mndo
        
        implicit none
      
        integer :: i, j, imm, l, ios, isec, temp_nat, mndot_nline,
     $  mndot_oline, mndot_eline
        logical :: file_exists
        character(len=1024), allocatable :: mndo_template(:)

        integer trimtext
        
        imm = 1
        do i=1, n
          isqm(i) = .false.
          do j=1, nqmatoms
            if(qmlist(j) .eq. i) isqm(i) = .true.
          end do
          if(.not. isqm(i)) then
            mmlist(imm) = i
            imm = imm + 1
          end if
        end do

        inquire(file=template_fname, exist=file_exists)
        if(.not. file_exists) then
          write(6, *) "Template file for MNDO does not exist."
          call fatal
        end if

        open(unit=temp_unit, file=template_fname, iostat=ios)
        
        mndot_nline = 0
        allocate(mndo_template(1))
        do
          read(temp_unit, '(A)', iostat=ios) mndo_template(1)
          if (ios .ne. 0) exit
          mndot_nline = mndot_nline + 1
        end do

        deallocate(mndo_template)
        allocate(mndo_template(mndot_nline))

        rewind(temp_unit)
        do i=1, mndot_nline
          read(temp_unit, '(A)') mndo_template(i)
        end do
          
        close(temp_unit)

        mndot_eline = 0
        mndot_oline = 0
        temp_nat = 0
        isec = 0

        do i=1, mndot_nline
          l = trimtext(mndo_template(i))
          if (isec .eq. 0) then
            write(6, *) '(O)  '
            mndot_oline = mndot_oline + 1
            if (mndo_template(i)(l:l) .ne. '+') isec = 1
          else if (isec .eq. 1) then
            write(6, *) '(H)  '
            if (i .eq. mndot_oline + 2) isec = 2
          else if (isec .eq. 2) then
            write(6, *) '(C)  '
            if (l .eq. 0) then
              isec = 3
            else 
              temp_nat = temp_nat + 1
            end if
          else if (isec .eq. 3) then
            write(6, *) '(E) '
            mndot_eline = mndot_eline + 1
          else 
            write(6, *) '(?)  '
          end if
          write(6,*) mndo_template(i)(:l)
        end do

c       Check keyword
        call mndo_parse_key(mndo_template, mndot_oline)

c       Copy "extra" lines
        do i=1, mndot_eline
          mndo_eline(i) = mndo_template(mndot_oline+nqmatoms+3+i)
        end do
        mndo_neline = mndot_eline

c       Sanity check
        if (temp_nat .ne. nqmatoms) then
          write(6, *) "Wrong number of QM atoms in template"
          write(6, *) "TEMPLATE ", temp_nat, "KEYFILE ", nqmatoms
          call fatal
        end if

        ismndoinit = .true.
      end
