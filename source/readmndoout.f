      subroutine readmndoout(filename,energy,deqm,demm,nqm,nmm)
      implicit none
      character*(*) filename 
      character*(120) line
      integer nqm, nmm
      real*8 energy, deqm(3,nqm), demm(3,nmm)
      real*8 fx, fy, fz
      integer read_unit, n, istatus, i
      logical deqm_section, demm_section, demm_header
      read_unit = 8
      open(unit=read_unit, file=filename, iostat=istatus)
      if (istatus.ne.0) then
        write(6,*) 'Failed to open MNDO log'
        call fatal
      end if
c
      deqm_section = .false.
      demm_section = .false.
      demm_header = .false.
c
      do
        read(read_unit, '(A)', iostat=istatus) line
        if (istatus.ne.0) exit
c
c       parse energy
c
        if (line(2:21).eq.'CI HEAT OF FORMATION') then
          read(line(22:37),'(F10.5)') energy
        end if
c
c       parse QM gradient
c
        if (deqm_section .and. (n.ge.5)) then
          if (len(trim(line)).eq.0) then
            if (i.le.nqm) then
              write(6,*) 'Not enough QM atoms found in log file'
              call fatal
            end if
            deqm_section = .false.
          else
            read(line(59:94),'(3F12.5)') fx, fy, fz
            if (i.gt.nqm) then
              write(6,*) 'Too many QM atoms found in log file'
              call fatal
            end if
            deqm(1,i) = fx
            deqm(2,i) = fy
            deqm(3,i) = fz 
            i = i + 1
          end if 
        end if
c
c       parse MM gradient
c
        if (demm_section .and. (n.ge.7)) then
          if (len(trim(line)).eq.0) then
            if (i.le.nmm) then
              write(6,*) "Not enough MM atoms found in log file"
              call fatal
            end if
            demm_section = .false.
          else
            read(line(59:94),'(3F12.5)') fx, fy, fz
            if (i.gt.nmm) then
              write(6,*) "Too many MM atoms found in log file"
              call fatal
            end if
            demm(1,i) = fx
            demm(2,i) = fy
            demm(3,i) = fz 
            i = i + 1
          end if
        end if
c
c       open sections
c
        if (line(27:48).eq.'COORDINATES (ANGSTROM)' .and.
     $      line(64:94).eq.'GRADIENTS (KCAL/(MOL*ANGSTROM))' .and.
     $      .not.demm_section) then
          deqm_section = .true.
          n = 1
          i = 1
        end if
c
        if (line(6:53).eq.'GRADIENT CONTRIBUTIONS TO EXTERNAL POINT CHAR
     $GES') then
          demm_section = .true.
          n = 1
          i = 1
        end if
        n = n + 1
      end do
      close(unit=8)
      end subroutine
