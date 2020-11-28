      subroutine emndoqm1()
        use atoms
        use energi
        use mndo
       
        character(len=2048) :: command
        integer :: status
        
        integer trimtext

        if(.not. ismndoinit) call mndoinit
        
        write(6, *) "DEBUG: EMNDOQM1"

c       Before each call remove old temp files.
        write(command, *) "rm -f ", mndo_in(:trimtext(mndo_in)), " ",
     &  mndo_out(:trimtext(mndo_out)), " ", 
     &  mndo_intfi(:trimtext(mndo_intfi))
        status = system(command)

c       Create new input file
        call mndomkin

c       Run MNDO
        write(command, *) mndo_exe(:trimtext(mndo_exe)), 
     &  "<", mndo_in(:trimtext(mndo_in)), ">&",
     &  mndo_out(:trimtext(mndo_out))
        status = system(command)
        write(6, *) "STATUS: ", status

c       Read the output to populate emndo demndo
        call mndordout
        write(6, *) "EMNDO: ", emndo

c       TODO Post-execution script
      end
