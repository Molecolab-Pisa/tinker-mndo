      subroutine emndoqm1()
        use atoms
        use energi
        use mndo
       
        character(len=2048) :: command
        integer :: status
        real*8, allocatable :: deqm(:,:), demm(:,:)
        
        integer trimtext

        write(6, *) "DEBUG: EMNDOQM1"
        call mndomkin
        write(command, *) mndo_exe(:trimtext(mndo_exe)), 
     &  "<", mndo_in(:trimtext(mndo_in)), ">",
     &  mndo_out(:trimtext(mndo_out))

        status = system(command)
        write(6, *) status

        allocate(deqm(3,nqmatoms))
        allocate(demm(3,n-nqmatoms))
        call readmndoout(mndo_out(:trimtext(mndo_out)), 
     &  emndo, deqm, demm, nqmatoms, n-nqmatoms)
        write(6, *) "EMNDO: ", emndo
        write(6, *) "FMM(1,1): ", demm(1,1)
        write(6, *) "FQM(1,1): ", deqm(1,1)

      end
