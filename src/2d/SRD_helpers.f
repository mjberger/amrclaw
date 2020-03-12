c
c ----------------------------------------------------------------------------
c
      double precision function bigconck(q,irr,mitot,mjtot,lwidth,nvar)
c
      use amr_module
      implicit double precision (a-h,o-z)

      dimension q(nvar,mitot,mjtot), irr(mitot,mjtot)

      totmass = 0.d0

      do 10 j = lwidth+1, mjtot-lwidth
      do 10 i = lwidth+1, mitot-lwidth

         k = irr(i,j)
         if (k .eq. -1) cycle 
         totmass = totmass + q(1,i,j)*ar(k)

 10   continue

      bigconck = totmass

      return
      end
c
c ------------------------------------------------------------------------
c 
      subroutine getAdjCell(i,j,k,ioff,joff,ffluxlen,gfluxlen,
     .                      mitot,mjtot,dx,dy)

      use amr_module
      implicit double precision (a-h, o-z)
      dimension ffluxlen(mitot+1,mjtot+1),gfluxlen(mitot+1,mjtot+1)
      dimension rl(4),iset(4),jset(4)
      data iset/-1, 1, 0, 0/
      data jset/0, 0, -1, 1/


c  for cell i,j return offsets for most normal neighbor
c  either look at normal direction, or full edge lengths
      
      rl(1) = ffluxlen(i,j)/dy
      rl(2) = ffluxlen(i+1,j)/dy
      rl(3) = gfluxlen(i,j)/dx
      rl(4) = gfluxlen(i,j+1)/dx

c take max over all and thats the neighbor
      ispot = 1
      do index = 2, 4
         if (rl(index) .gt. rl(ispot)) ispot = index  ! always positive
      end do
      ioff = iset(ispot)
      joff = jset(ispot)


!!    will need to check that i+ioff,j+joff is inside
!!    domain and valid, but for cylinder case for now
c     skip it

!     test if above algorithm is same as max component of normal
      do 20 kside=1,6
         if (poly(kside+2,1,k).eq.-11.) then
            x1 = poly(kside,1,k)
            y1 = poly(kside,2,k)
            x2 = poly(kside+1,1,k)
            y2 = poly(kside+1,2,k)
            go to 25
         endif
 20   continue
 25   continue
      rlen = dsqrt((y2-y1)**2+(x2-x1)**2)
      alf = (y2-y1)/rlen
      beta = (x1-x2)/rlen
      if (abs(alf) .gt. abs(beta)) then
        if (alf .lt. 0.d0) then 
           iin = 1
        else
           iin = 2
        endif
      else
        if (beta .lt. 0.d0) then
           iin = 3
        else
           iin = 4
        endif
      endif
      if (ispot .ne. iin) then
         write(*,*)" tests for adjacent cell differ for i,j = ",i,j
         ! trust this latter one
         ioff = iset(iin)
         joff = jset(iin)
      endif

      return
      end

c
c -------------------------------------------------------------------------
c
      subroutine limitTileCellBJ(qmerge,gradmx,gradmy,xcentMerge,
     &                           ycentMerge,irr,nvar,mitot,mjtot,
     &                           nborList,nborCount,i,j)

       use amr_module
       implicit double precision (a-h,o-z)

       dimension qmerge(nvar,mitot,mjtot), irr(mitot,mjtot)
       dimension gradmx(nvar,mitot,mjtot), gradmy(nvar,mitot,mjtot)
       dimension xcentMerge(mitot,mjtot),ycentMerge(mitot,mjtot)

       dimension dumax(nvar),dumin(nvar),phimin(nvar),graddot(nvar)
       dimension alpha(nvar),recon(nvar)
       dimension nborList(25,2)

          
            ! find max and min needed for BJ limiting. Use all nbors used in gradient comp
       eps = 1d-4
       dumax = 0.d0 
       dumin = 0.d0
       do 31 ico = 1, nborCount
          inbor = nborList(ico,1)
          jnbor = nborList(ico,2)
          dumax = max(dumax,qmerge(:,inbor,jnbor)-qmerge(:,i,j))
          dumin = min(dumin,qmerge(:,inbor,jnbor)-qmerge(:,i,j))
 31    continue

       phimin = 1.d0
       do 32 ico = 1, nborCount
          inbor = nborList(ico,1)
          jnbor = nborList(ico,2)
          diffx = xcentMerge(inbor,jnbor)-xcentMerge(i,j)
          diffy = ycentMerge(inbor,jnbor)-ycentMerge(i,j)
          graddot  = gradmx(:,i,j)*diffx + gradmy(:,i,j)*diffy
          do m = 1,nvar
             if (graddot(m) > 0.d0) then
                alpha(m) = min(1.d0, dumax(m)/graddot(m))
             else if (graddot(m) < 0.d0) then
                alpha(m) = min(1.d0, dumin(m)/graddot(m))
             else
                alpha(m) = 1.d0
             endif
          end do

! one last check for positivity
          recon = qMerge(:,i,j) + graddot  
          xymomsq = recon(2)**2+recon(3)**2
          press = .4d0*(recon(4)-0.5d0*xymomsq/recon(1))
          if (recon(1).le.eps .or. press.le.eps) alpha = 0.d0
          phimin = min(phimin, alpha)
 32    continue
       gradmx(:,i,j) = gradmx(:,i,j)*phimin(:)
       gradmy(:,i,j) = gradmy(:,i,j)*phimin(:)

       return 
       end
c
c -------------------------------------------------------------------------
c
      subroutine limitTileGradientBJ(qmerge,gradmx,gradmy,xcentMerge,
     &                               ycentMerge,xlow,ylow,dx,dy,irr,
     &                               lwidth,nvar,mitot,mjtot,lstgrd,
     &                               areaMin,nborList,nborCount)

       use amr_module
       implicit double precision (a-h,o-z)

       dimension qmerge(nvar,mitot,mjtot), irr(mitot,mjtot)
       dimension gradmx(nvar,mitot,mjtot), gradmy(nvar,mitot,mjtot)
       dimension xcentMerge(mitot,mjtot),ycentMerge(mitot,mjtot)

       dimension dumax(nvar),dumin(nvar),phimin(nvar),graddot(nvar)
       dimension alpha(nvar),recon(nvar)
       dimension nborList(25,2)
       logical  IS_OUTSIDE

       IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                    y .lt. ylower .or. y .gt. yupper)

      
       eps = 1d-4
c       nco = 1   ! size of nhood for limiting tiles
        nco = 2 ! but throw away l1 distance 3 points

c      do 30 j = lwidth-1, mjtot-lwidth+1
c      do 30 i = lwidth-1, mitot-lwidth+1
       do 30 j = lwidth+1, mjtot-lwidth
       do 30 i = lwidth+1, mitot-lwidth

          k = irr(i,j)
          if (k .eq. -1) go to 30  ! solid cells have no gradient
          if (k .eq. lstgrd) go to 30 ! full cells dont use a gradient
          if (ar(k) .gt. areaMin) go to 30 ! no gradient needed

            ! find max and min needed for BJ limiting
            dumax = 0.d0 
            dumin = 0.d0
            do 31 joff = -nco, nco
            do 31 ioff = -nco, nco
              if (ioff .eq. 0 .and. joff .eq. 0) go to 31
              if (iabs(ioff)+iabs(joff) .gt. 2)  go to 31 ! too far
              koff = irr(i+ioff,j+joff)
              if (koff .eq. -1) go to 31
              call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                             xlow,ylow,dx,dy,koff)
              if (IS_OUTSIDE(xcn,ycn)) go to 31
              dumax = max(dumax,qmerge(:,i+ioff,j+joff)-qmerge(:,i,j))
              dumin = min(dumin,qmerge(:,i+ioff,j+joff)-qmerge(:,i,j))
 31         continue

           phimin = 1.d0
            do 32 joff = -nco, nco
            do 32 ioff = -nco, nco
                !! there may be an equation here since have to reconstruct
                !! to cell center, not merge center
                !! but presumably one would only check positivity since nothing
                !! to limit using same qMerge nhood
                if (ioff .eq. 0 .and. joff .eq. 0) go to 32 
                koff = irr(i+ioff,j+joff)
                if (koff .eq. -1) go to 32
                call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                               xlow,ylow,dx,dy,koff)
                if (IS_OUTSIDE(xcn,ycn)) go to 32
                ! use merged val at ioff,joff to limit
                !diffx = xcentMerge(i+ioff,j+joff)-xcentMerge(i,j)
                !diffy = ycentMerge(i+ioff,j+joff)-ycentMerge(i,j)
                ! these next lines limit at cell where will evaluate
                ! previous lines limit at cells used in making gradient
c                diffx = xcn-xcentMerge(i,j)
c                diffy = ycn-ycentMerge(i,j)
                diffx = xcentMerge(i+ioff,j+joff)-xcentMerge(i,j)
                diffy = ycentMerge(i+ioff,j+joff)-ycentMerge(i,j)
                graddot  = gradmx(:,i,j)*diffx + gradmy(:,i,j)*diffy
                  do m = 1,4
                     if (graddot(m) > 0.d0) then
                        alpha(m) = min(1.d0, dumax(m)/graddot(m))
                     else if (graddot(m) < 0.d0) then
                        alpha(m) = min(1.d0, dumin(m)/graddot(m))
                     else
                        alpha(m) = 1.d0
                     endif
                  end do
                  ! one last check for positivity
                  recon = qMerge(:,i,j) + graddot  
                  xymomsq = recon(2)**2+recon(3)**2
                  press = .4d0*(recon(4)-0.5d0*xymomsq/recon(1))
                  if (recon(1).le.eps .or. press.le.eps) alpha = 0.d0
                  phimin = min(phimin, alpha)
 32         continue
            gradmx(:,i,j) = gradmx(:,i,j)*phimin(:)
            gradmy(:,i,j) = gradmy(:,i,j)*phimin(:)
 30     continue

        return 
        end
c
c -------------------------------------------------------------------------
c
      subroutine limitTileGradientLP(qmerge,gradmx,gradmy,xcentMerge,
     &                               ycentMerge,xlow,ylow,dx,dy,irr,
     &                               lwidth,nvar,mitot,mjtot,
     &                               lstgrd,areaMin,mptr,lpChoice,
     &                               nborList,nborCount)

       use amr_module
       implicit double precision (a-h,o-z)

       dimension qmerge(nvar,mitot,mjtot), irr(mitot,mjtot)
       dimension gradmx(nvar,mitot,mjtot), gradmy(nvar,mitot,mjtot)
       dimension xcentMerge(mitot,mjtot),ycentMerge(mitot,mjtot)

       dimension dumax(nvar),dumin(nvar),phimin(nvar),graddot(nvar)
       dimension alpha(nvar),recon(nvar), nborList(25,2)
       logical  IS_OUTSIDE, INDEX_OUTSIDE

c      sandras variables
       parameter (max_num_nbor=22)
       logical  LP_success, output_debug
       double precision neighb_x(max_num_nbor)
       double precision neighb_y(max_num_nbor)
       double precision neighb_u(max_num_nbor), allu(nvar,max_num_nbor)
       double precision Phi_x, Phi_y, LP_Dx,LP_Dy, LS_Dx, LS_Dy
       integer row_A
       logical do_pos_constr ! are we really a cutcell?
       integer my_iters, tot_iters, max_iters   ! maximum no of iterations needed in Simplex
       double precision avg_iters    ! avg number of iterations needed in Simplex
       integer count_cells
       logical added_neighbor


       IS_OUTSIDE(x,y) = (x .lt. xlower .or. x .gt. xupper .or.
     .                    y .lt. ylower .or. y .gt. yupper)
       INDEX_OUTSIDE(i,j) = (i .lt. 1 .or. i .gt. mitot .or.
     .                       j .lt. 1 .or. j .gt. mjtot)

       eps = 1d-4
       nco = 1   ! size of nhood for limiting tiles

c       do 30 j = lwidth-1, mjtot-lwidth+1
c       do 30 i = lwidth-1, mitot-lwidth+1
       do 30 j = lwidth+1, mjtot-lwidth
       do 30 i = lwidth+1, mitot-lwidth

          k = irr(i,j)
          if (k .eq. -1) go to 30  ! solid cells have no gradient
          if (k .eq. lstgrd) go to 30 ! full cells dont use a gradient
          if (ar(k) .gt. areaMin) go to 30 ! no gradient needed

            ! convert to sandra lp routine
            call getBndryPt(bxpt,bypt,k)
            center_x = xcentMerge(i,j)
            center_y = ycentMerge(i,j)
            num_neighb = 0

            do 31 joff = -nco, nco
            do 31 ioff = -nco, nco
              if (ioff .eq. 0 .and. joff .eq. 0) go to 31
              if (INDEX_OUTSIDE(i+ioff,j+joff)) go to 31
              koff = irr(i+ioff,j+joff)
              if (koff .eq. -1) go to 31
              call getCellCentroid(lstgrd,i+ioff,j+joff,xcn,ycn,
     &                             xlow,ylow,dx,dy,koff)
             xcn = xcentMerge(i+ioff,j+joff)
             ycn = ycentMerge(i+ioff,j+joff)
              if (IS_OUTSIDE(xcn,ycn)) go to 31
              num_neighb = num_neighb + 1
              output_debug = .false.
              neighb_x(num_neighb) = xcn
              neighb_y(num_neighb) = ycn
              allu(:,num_neighb) = qmerge(:,i+ioff,j+joff)
 31         continue
            row_A     = 2*num_neighb + 4
            do_pos_constr = .false.
              do m = 1, nvar
                 neighb_u(1:num_neighb) = allu(m,1:num_neighb)
                 center_u = qmerge(m,i,j)
                 LS_Dx = gradmx(m,i,j)
                 LS_Dy = gradmy(m,i,j)
                 LP_Success = .false.
                 call vector_lim_allineq(num_neighb,neighb_x,neighb_y,
     &               neighb_u,center_x,center_y,center_u,output_debug,
     &               row_A,Phi_x,Phi_y, LP_Dx,LP_Dy,
     &               LP_success,LS_Dx,LS_Dy,bxpt,bypt,do_pos_constr,
     &               my_iters,lpChoice)
                if (LP_Success) then
                   gradmx(m,i,j) = LP_Dx
                   gradmy(m,i,j) = LP_Dy
                else
                   write(*,*)"LP failure for grid tile var ",
     &                       mptr,i,j,m
                   gradmx(m,i,j) = 0.d0
                   gradmy(m,i,j) = 0.d0
                endif
              end do
 30     continue

        return 
        end
c
c ------------------------------------------------------------------
c
       subroutine getBndryPt(bxpt,bypt,k)

       use amr_module
       implicit double precision (a-h,o-z)

      do 20 kside=1,6
        if (poly(kside+2,1,k).eq.-11.) then
           x1 = poly(kside,1,k)
           y1 = poly(kside,2,k)
           x2 = poly(kside+1,1,k)
           y2 = poly(kside+1,2,k)
           go to 25
        endif
 20   continue

 25   continue

      bxpt = .5d0*(x2+x1)
      bypt = .5d0*(y2+y1)

      return
      end
c
c ------------------------------------------------------------------
c
       subroutine getBndryInfo(alf,beta,k,bxpt,bypt)

       use amr_module
       implicit double precision (a-h,o-z)


      do 20 kside=1,6
        if (poly(kside+2,1,k).eq.-11.) then
           x1 = poly(kside,1,k)
           y1 = poly(kside,2,k)
           x2 = poly(kside+1,1,k)
           y2 = poly(kside+1,2,k)
           go to 25
        endif
 20   continue

 25   continue

      bxpt = .5d0*(x2+x1)
      bypt = .5d0*(y2+y1)

      rlenb = dsqrt((y1-y2)**2 + (x1-x2)**2)
      alf  = (y2-y1)/rlenb
      beta = (x1-x2)/rlenb

      return
      end
c
c ------------------------------------------------------------------
c