      subroutine emndoqm1()
        use atoms
        use energi
        use mndo
       
        character(len=2048) :: command
        integer :: status
        
        integer trimtext

        if(.not. ismndoinit) call mndoinit
        
        write(6, *) "DEBUG: EMNDOQM1"
        
        call mndomkin
        
        write(command, *) mndo_exe(:trimtext(mndo_exe)), 
     &  "<", mndo_in(:trimtext(mndo_in)), ">",
     &  mndo_out(:trimtext(mndo_out))

        status = system(command)
        write(6, *) "STATUS: ", status

        call mndordout
        write(6, *) "EMNDO: ", emndo

      end
