      subroutine emndoqm1()
        use atoms
        use energi
        use mndo
       
        character(len=2048) :: command
        integer :: status
        
        integer trimtext

        if(.not. ismndoinit) call mndoinit
        
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

c       Read the output to populate emndo demndo
        call mndordout

c       TODO Post-execution script
      end
