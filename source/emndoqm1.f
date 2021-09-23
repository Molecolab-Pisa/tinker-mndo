      subroutine emndoqm1()
        use atoms
        use energi
        use mndo
       
        character(len=2048) :: command
        integer :: status
        
        integer trimtext

        if(mndo_debug) write(6, *) "Entering EMNDOQM1"

c       Before each call remove old temp files.
        write(command, *) "rm -f ", mndo_in(:trimtext(mndo_in)), " ",
     &  mndo_out(:trimtext(mndo_out)), " ", 
     &  mndo_intfi(:trimtext(mndo_intfi))
        if(mndo_debug) write(6, *) command(:trimtext(command))
        status = system(command)

c       Create new input file with jop=-2 for gradients calculation
        write(mndo_keyword(mndo_nwk+1), "(A,'=',A)") "jop", "-2"
        mndo_nwk = mndo_nwk + 1
        
        call mndomkin
        
        mndo_nwk = mndo_nwk - 1

c       Run MNDO
        write(command, *) mndo_exe(:trimtext(mndo_exe)), 
     &  " < ", mndo_in(:trimtext(mndo_in)), " > ",
     &  mndo_out(:trimtext(mndo_out)), " 2>&1 "
        if(mndo_debug) write(6, *) command(:trimtext(command))
        status = system(command)

c       Read the output to populate emndo demndo
        if(mndo_multistate) then
          call mndordmsout()
        else
          call mndordout(.true.)
        end if

c       Post-execution script
        if(mndo_dope) then
          write(command, *) mndo_postexe
          if(mndo_debug) write(6, *) command(:trimtext(command))
          status = system(command)
        end if
      end
