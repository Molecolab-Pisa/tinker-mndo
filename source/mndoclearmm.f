c
c     strip bonded lists from terms not required in a QM/MM calculation
c
      subroutine mndoclearmm
      use mndo
      use sizes
      use atoms
      use bndstr
      use atmlst
      implicit none
      integer i, iat, jat, istatus
      integer nqm, maxbnd
      integer new_nbond
      integer, allocatable :: iscratch(:,:), jscratch(:,:)
c     strip stretchings: nbond, ibnd, bndlist
      maxbnd = 4*n
c     for the moment it's better to make a scratch copy of ibnd
      allocate(iscratch(2,maxbnd)
      new_nbond = 0
      jscratch = 0
      do i = 1, nbond
        iat = ibnd(1,i)
        jat = ibnd(2,i)
        nqm = 0
        if (isqm(iat)) nqm = nqm + 1
        if (isqm(jat)) nqm = nqm + 1
        if (nqm.eq.2 .and. mndo_debug) write(6,*) 'Found a bond between 
     $two QM atoms: removing it'
        if (nqm.lt.2) then
          new_nbond = new_nbond + 1
          iscratch(1,new_nbond) = iat
          iscratch(2,new_nbond) = jat
        else
          bndlist(:,iat) = 0
        end if
      end do
      nbond = new_nbond
      ibnd = iscratch
      deallocate(iscratch)
c
c     strip angles: nangle, iang, anglist, balist
c
      maxang = 6*n
      allocate(iscratch(4,maxang))
      do i = 1, nangle
      end do
      end
