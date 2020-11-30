c
c     strip bonded lists from terms not required in a QM/MM calculation
c
      subroutine mndoclearmm
      use mndo
      use sizes
      use atoms
      use bndstr
      use atmlst
      use angbnd
      implicit none
      integer i, iat, jat, kat, lat, istatus
      integer nqm, maxbnd, maxang
      integer new_index
      integer, allocatable :: iscratch(:,:), jscratch(:,:)
c
c     strip stretchings: nbond, ibnd, bndlist
c
      maxbnd = 4*n
c     for the moment it's better to make a scratch copy of ibnd
      allocate(iscratch(2,maxbnd))
      new_index = 0
      do i = 1, nbond
        iat = ibnd(1,i)
        jat = ibnd(2,i)
        nqm = 0
        if (isqm(iat)) nqm = nqm + 1
        if (isqm(jat)) nqm = nqm + 1
        if (nqm.ge.2 .and. mndo_debug) write(6,*) 'Found a bond between 
     $two QM atoms: removing it'
        if (nqm.lt.2) then
          new_index = new_index + 1
          iscratch(1,new_index) = iat
          iscratch(2,new_index) = jat
        else
          bndlist(:,iat) = 0
        end if
      end do
      nbond = new_index
      ibnd = iscratch
      deallocate(iscratch)
c
c     strip angles: nangle, iang, anglist, balist
c
      new_index = 0
      maxang = 6*n
      allocate(iscratch(4,maxang))
      do i = 1, nangle
        iat = iang(1,i)
        jat = iang(2,i)
        kat = iang(3,i)
c       lat can either be 0 for plain angles or an index for 
c       pyramidalizations, in either case it's okay don't checking it
        lat = iang(4,i)

        nqm = 0
        if (isqm(iat)) nqm = nqm + 1
        if (isqm(jat)) nqm = nqm + 1
        if (isqm(kat)) nqm = nqm + 1
        if (nqm.ge.2 .and. mndo_debug) write(6,*) 'Found an angle involv
     $ing more than one QM atom: removing it'
        if (nqm.lt.2) then
          new_index = new_index + 1
          iscratch(1,new_index) = iat
          iscratch(2,new_index) = jat
          iscratch(3,new_index) = kat
          iscratch(4,new_index) = lat
          anglist(:,iat) = 0
        end if
      end do
      nangle = new_index
      iang = iscratch
      deallocate(iscratch)
      end
