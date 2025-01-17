c
c
c     ###################################################
c     ##  COPYRIGHT (C)  1995  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     ###########################################################
c     ##                                                       ##
c     ##  subroutine eopbend  --  out-of-plane bending energy  ##
c     ##                                                       ##
c     ###########################################################
c
c
c     "eopbend" computes the out-of-plane bend potential energy at
c     trigonal centers via a Wilson-Decius-Cross or Allinger angle
c
c
      subroutine eopbend
      use angbnd
      use angpot
      use atoms
      use bound
      use energi
      use fields
      use group
      use math
      use opbend
      use usage
      use mndo
      implicit none
      integer i,iopbend,nqm
      integer ia,ib,ic,id
      real*8 e,force,angle
      real*8 sine,fgrp
      real*8 dt,dt2,dt3,dt4
      real*8 xia,yia,zia
      real*8 xib,yib,zib
      real*8 xic,yic,zic
      real*8 xid,yid,zid
      real*8 xab,yab,zab
      real*8 xcb,ycb,zcb
      real*8 xdb,ydb,zdb
      real*8 xad,yad,zad
      real*8 xcd,ycd,zcd
      real*8 rdb2,rad2,rcd2
      real*8 rab2,rcb2
      real*8 cc,ee
      logical proceed
c
c
c     zero out the out-of-plane bending energy component
c
      eopb = 0.0d0
      if (nopbend .eq. 0)  return
c
c     OpenMP directives for the major loop structure
c
!$OMP PARALLEL default(private) shared(nopbend,iopb,iang,opbk,use,
!$OMP& x,y,z,opbtyp,copb,qopb,popb,sopb,opbunit,use_group,use_polymer)
!$OMP& shared(eopb,isqm)
!$OMP DO reduction(+:eopb) schedule(guided)
c
c     calculate the out-of-plane bending energy term
c
      do iopbend = 1, nopbend
         i = iopb(iopbend)
         ia = iang(1,i)
         ib = iang(2,i)
         ic = iang(3,i)
         id = iang(4,i)
c
c     skip interaction if mor than one of the atoms is qm
c
         nqm = 0
         if (isqm(ia)) nqm = nqm + 1
         if (isqm(ib)) nqm = nqm + 1
         if (isqm(ic)) nqm = nqm + 1
         if (isqm(id)) nqm = nqm + 1
         if (nqm.gt.1) cycle
         force = opbk(iopbend)
c
c     decide whether to compute the current interaction
c
         proceed = .true.
         if (use_group)  call groups (proceed,fgrp,ia,ib,ic,id,0,0)
         if (proceed)  proceed = (use(ia) .or. use(ib) .or.
     &                              use(ic) .or. use(id))
c
c     get the coordinates of the atoms at trigonal center
c
         if (proceed) then
            xia = x(ia)
            yia = y(ia)
            zia = z(ia)
            xib = x(ib)
            yib = y(ib)
            zib = z(ib)
            xic = x(ic)
            yic = y(ic)
            zic = z(ic)
            xid = x(id)
            yid = y(id)
            zid = z(id)
c
c     compute the out-of-plane bending angle
c
            xab = xia - xib
            yab = yia - yib
            zab = zia - zib
            xcb = xic - xib
            ycb = yic - yib
            zcb = zic - zib
            xdb = xid - xib
            ydb = yid - yib
            zdb = zid - zib
            xad = xia - xid
            yad = yia - yid
            zad = zia - zid
            xcd = xic - xid
            ycd = yic - yid
            zcd = zic - zid
            if (use_polymer) then
               call image (xab,yab,zab)
               call image (xcb,ycb,zcb)
               call image (xdb,ydb,zdb)
               call image (xad,yad,zad)
               call image (xcd,ycd,zcd)
            end if
c
c     W-D-C angle between A-B-C plane and B-D vector for D-B<AC
c
            if (opbtyp .eq. 'W-D-C') then
               rab2 = xab*xab + yab*yab + zab*zab
               rcb2 = xcb*xcb + ycb*ycb + zcb*zcb
               cc = rab2*rcb2 - (xab*xcb+yab*ycb+zab*zcb)**2
c
c     Allinger angle between A-C-D plane and D-B vector for D-B<AC
c
            else if (opbtyp .eq. 'ALLINGER') then
               rad2 = xad*xad + yad*yad + zad*zad
               rcd2 = xcd*xcd + ycd*ycd + zcd*zcd
               cc = rad2*rcd2 - (xad*xcd+yad*ycd+zad*zcd)**2
            end if
c
c     find the out-of-plane angle bending energy
c
            ee = xdb*(yab*zcb-zab*ycb) + ydb*(zab*xcb-xab*zcb)
     &              + zdb*(xab*ycb-yab*xcb)
            rdb2 = xdb*xdb + ydb*ydb + zdb*zdb
            if (rdb2.ne.0.0d0 .and. cc.ne.0.0d0) then
               sine = abs(ee) / sqrt(cc*rdb2)
               sine = min(1.0d0,sine)
               angle = radian * asin(sine)
               dt = angle
               dt2 = dt * dt
               dt3 = dt2 * dt
               dt4 = dt2 * dt2
               e = opbunit * force * dt2
     &                * (1.0d0+copb*dt+qopb*dt2+popb*dt3+sopb*dt4)
c
c     scale the interaction based on its group membership
c
               if (use_group)  e = e * fgrp
c
c     increment the total out-of-plane bending energy
c
               eopb = eopb + e
            end if
         end if
      end do
c
c     OpenMP directives for the major loop structure
c
!$OMP END DO
!$OMP END PARALLEL
      return
      end
