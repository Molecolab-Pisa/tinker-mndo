      subroutine mndoinit()
        use atoms
        use mndo
        
        implicit none
      
        integer :: i, j, imm, l, ios, isec, temp_nat, mndot_nline,
     $  mndot_oline
        logical :: file_exists
        character(len=1024), allocatable :: mndo_template(:)

        integer trimtext
       
        write(6, *) "+++ MNDO INITIALIZATION +++"

c       Build MM/QM atoms list
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

        call mndolainit

c       Read template file and save keyword
        inquire(file=template_fname, exist=file_exists)
        if(.not. file_exists) then
          write(6, *) "Template file for MNDO (",
     &    template_fname(:trimtext(template_fname)), ") does not exist."
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

        mndo_neline = 0
        mndot_oline = 0
        temp_nat = 0
        isec = 0

        do i=1, mndot_nline
          l = trimtext(mndo_template(i))
          if (isec .eq. 0) then
            if(mndo_debug) write(6, "(A)", advance="no" ) '(O)  '
            mndot_oline = mndot_oline + 1
            if (mndo_template(i)(l:l) .ne. '+') isec = 1
          else if (isec .eq. 1) then
            if(mndo_debug) write(6, "(A)", advance="no" ) '(H)  '
            if (i .eq. mndot_oline + 2) isec = 2
          else if (isec .eq. 2) then
            if(mndo_debug) write(6, "(A)", advance="no" ) '(C)  '
            if (l .eq. 0) then
              isec = 3
            else 
              temp_nat = temp_nat + 1
            end if
          else if (isec .eq. 3) then
            if(mndo_debug) write(6, "(A)", advance="no" ) '(E)  '
            mndo_neline = mndo_neline + 1
            mndo_eline(mndo_neline) = mndo_template(i)
          else 
            if(mndo_debug) write(6, "(A)", advance="no" ) '(?)  '
          end if
            if(mndo_debug) write(6, "(A)") mndo_template(i)(:l) 
        end do

c       Check keyword
        call mndo_parse_key(mndo_template, mndot_oline)

c       Sanity check of template
c        if (temp_nat .ne. nqmatoms+mndo_nla) then
c          write(6, *) "Wrong number of QM atoms in template"
c          write(6, *) "TEMPLATE ", temp_nat, "KEYFILE ", nqmatoms
c          call fatal
c        end if

c       Debug information        
        if(mndo_debug) then
          write(6, *) "=== MNDO OPTIONS ==="
          write(6, *) "  MNDO TEMPLATE: ",
     &      template_fname(:trimtext(template_fname))
          write(6, *) "  MNDO EXECUTABLE: ",
     &      mndo_exe(:trimtext(mndo_exe))
          write(6, "('   QM ATOMS: ', I5)") nqmatoms
          write(6, "('   MM ATOMS: ', I5)") n - nqmatoms
          if(mndo_dope) 
     &      write(6, *) "  MNDO POST EXECUTION SCRIPT: ",
     &      mndo_postexe(:trimtext(mndo_postexe))
          if(mndo_iterguess) write(6, *) "  MNDO ITER GUESS IS ON"
          write(6, *) "===================="
        end if
        write(6, *) "+++ MNDO INITIALIZATION FINISHED +++"

        ismndoinit = .true.
      end
