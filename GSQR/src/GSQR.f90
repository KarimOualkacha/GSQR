

      SUBROUTINE sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        IMPLICIT NONE
        INTEGER:: nobs
        INTEGER:: nvars
        INTEGER:: i
        DOUBLE PRECISION:: taux
        DOUBLE PRECISION:: k
        DOUBLE PRECISION:: dl(nobs)
        DOUBLE PRECISION:: r(nobs)
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION:: vl(nvars)
        vl = 0.0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        ENDDO
        vl = matmul(dl, x) / nobs
        END SUBROUTINE sqrdrv2

        SUBROUTINE sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        IMPLICIT NONE
        INTEGER:: nobs
        INTEGER:: nvars
        INTEGER:: i
        DOUBLE PRECISION:: taux
        DOUBLE PRECISION:: c
        DOUBLE PRECISION:: dl(nobs)
        DOUBLE PRECISION:: r(nobs)
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION:: vl(nvars)
        vl = 0.0
        DO i = 1, nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        ENDDO
        vl = matmul(dl, x) / nobs
        END SUBROUTINE sqrdrv1

! --------------------------------------------------! -------------------------
! -------------------------------------------------! --------------------------
! --------quantile regression  with group lasso penalite--------------------
! ---------------------------------------------f1 and f2-----------------------
! --------------------------------------------------! -------------------------

SUBROUTINE sqr1lasso (taux,c,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::c
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::delta
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
          delta=c/max(taux,1-taux)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)

        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
!----------------------------------------------------------------------------------------

        DO l=1,nlam
        
        write(*,*) l

        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0

        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO

        al = al0 * alf
        ENDIF
        ENDIF

        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))

        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
		dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.50 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
  
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE sqr1lasso
! ----------------------------------------------------------------------------
! -----------------------------------------------------------------------------
SUBROUTINE sqr2lasso (taux,k,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                                eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::k
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::delta
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION:: dd1
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
          delta=k
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)

        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
!----------------------------------------------------------------------------------------

        DO l=1,nlam

        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0

        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO

        al = al0 * alf
        ENDIF
        ENDIF

        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb

        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))

        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
		dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.50 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
        jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
  
        IF (jx /= 0) CYCLE
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE sqr2lasso
! --------------------------------------------------! -------------------------
! -------------------------------------------------! --------------------------
! -------------quantile regression  with group mcp penalite--------------------
! ---------------------------------------------f1 and f2-----------------------
! --------------------------------------------------! -------------------------
! --------------------------------------------------! -------------------------
SUBROUTINE sqr1mcp (gamm, taux,c,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
    eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
        ! --------------------------------------------------
          IMPLICIT NONE
        ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::c
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION::delta
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
          DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
          DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
          delta=c/max(taux,1-taux)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
          IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
        DO l=1,nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0
        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO
        al = al0 * alf
        ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
          DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))

        t=unorm-pf(g)*al

        IF (t>0.0D0) THEN
				tt=unorm-gam(g)*gamm*pf(g)*al
				IF (tt>0.0D0) THEN
						b(start:end) = u/gam(g)
				ELSE
						b(start:end)=(u*t/((gam(g)-1/gamm)*unorm))
				END IF
        ELSE
				b(start:end)=0.0D0

	END IF

        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))


        t=unorm-pf(g)*al

        IF (t>0.0D0) THEN
				tt=unorm-gam(g)*gamm*pf(g)*al
				IF (tt>0.0D0) THEN
						b(start:end) = u/gam(g)
				ELSE
						b(start:end)=(u*t/((gam(g)-1/gamm)*unorm))
				END IF
        ELSE
				b(start:end)=0.0D0

	END IF


        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
         jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE sqr1mcp
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
SUBROUTINE sqr2mcp (gamm, taux,k,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
    eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
        ! --------------------------------------------------
          IMPLICIT NONE
        ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::k
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION::delta
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
          DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
          DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
          delta=k
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
          IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
        DO l=1,nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0
        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO
        al = al0 * alf
        ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
          DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))

        t=unorm-pf(g)*al

        IF (t>0.0D0) THEN
				tt=unorm-gam(g)*gamm*pf(g)*al
				IF (tt>0.0D0) THEN
						b(start:end) = u/gam(g)
				ELSE
						b(start:end)=(u*t/((gam(g)-1/gamm)*unorm))
				END IF
        ELSE
				b(start:end)=0.0D0

	END IF

        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))


        t=unorm-pf(g)*al

        IF (t>0.0D0) THEN
				tt=unorm-gam(g)*gamm*pf(g)*al
				IF (tt>0.0D0) THEN
						b(start:end) = u/gam(g)
				ELSE
						b(start:end)=(u*t/((gam(g)-1/gamm)*unorm))
				END IF
        ELSE
				b(start:end)=0.0D0

	END IF


        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
         jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE sqr2mcp

! --------------------------------------------------! -------------------------
! -------------------------------------------------! --------------------------
! -----------------------quantile regression  with group scad penalite----
! ---------------------------------------------f1 and f2-----------------------
! --------------------------------------------------! -------------------------
! --------------------------------------------------! -------------------------

SUBROUTINE sqr1scad (gamm, taux,c,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
      eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)

        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::c
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION::delta
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::ts
        DOUBLE PRECISION::tg
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
          DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
          delta=c/max(taux,1-taux)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
          IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
        DO l=1,nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0
        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO
        al = al0 * alf
        ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
          DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (unorm <= ((pf (g)+gam(g))*al)) THEN
        	t=unorm-pf(g)*al
        	IF(t>0.0D0) THEN
        		b(start:end)=u*t/(gam(g)*unorm)
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF
        ELSE IF (unorm <= (gam(g)*al*gamm)) THEN

        	tt=unorm-pf(g)*al*(gamm/(gamm-1))

        	IF(tt>0.0D0) THEN
        		b(start:end)=u*tt/((unorm)*(gam(g)-pf(g)/(gamm-1.0D0)))
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF                     

        ELSE
        	b(start:end) = u/gam(g)
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (unorm <= ((pf (g)+gam(g))*al)) THEN
        	t=unorm-pf(g)*al
        	IF(t>0.0D0) THEN
        		b(start:end)=u*t/(gam(g)*unorm)
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF
        ELSE IF (unorm <= (gam(g)*al*gamm)) THEN

        	tt=unorm-pf(g)*al*(gamm/(gamm-1))

        	IF(tt>0.0D0) THEN
        		b(start:end)=u*tt/((unorm)*(gam(g)-pf(g)/(gamm-1.0D0)))
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF                     

        ELSE
        	b(start:end) = u/gam(g)
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
          jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE sqr1scad
! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
SUBROUTINE sqr2scad (gamm, taux,k,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
      eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)

        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::k
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION::delta
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::tt
        DOUBLE PRECISION::ts
        DOUBLE PRECISION::tg
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
          DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
          delta=k
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
          IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO
        DO l=1,nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0
        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO
        al = al0 * alf
        ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
          DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (unorm <= ((pf (g)+gam(g))*al)) THEN
        	t=unorm-pf(g)*al
        	IF(t>0.0D0) THEN
        		b(start:end)=u*t/(gam(g)*unorm)
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF
        ELSE IF (unorm <= (gam(g)*al*gamm)) THEN

        	tt=unorm-pf(g)*al*(gamm/(gamm-1))

        	IF(tt>0.0D0) THEN
        		b(start:end)=u*tt/((unorm)*(gam(g)-pf(g)/(gamm-1.0D0)))
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF                     

        ELSE
        	b(start:end) = u/gam(g)
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF (unorm <= ((pf (g)+gam(g))*al)) THEN
        	t=unorm-pf(g)*al
        	IF(t>0.0D0) THEN
        		b(start:end)=u*t/(gam(g)*unorm)
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF
        ELSE IF (unorm <= (gam(g)*al*gamm)) THEN

        	tt=unorm-pf(g)*al*(gamm/(gamm-1))

        	IF(tt>0.0D0) THEN
        		b(start:end)=u*tt/((unorm)*(gam(g)-pf(g)/(gamm-1.0D0)))
        	ELSE
        		b(start:end)=0.0D0
        	ENDIF                     

        ELSE
        	b(start:end) = u/gam(g)
        END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
          jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE sqr2scad

! --------------------------------------------------! -------------------------
! -------------------------------------------------! --------------------------
! ------------quantile regression  with group mcp penalite aproximation-------
! ---------------------------------------------f1 and f2-----------------------
! --------------------------------------------------! -------------------------
! --------------------------------------------------! -------------------------
SUBROUTINE asqr1mcp (gamm, taux,c,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                        eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::c
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::delta
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
          DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
		DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
		DOUBLE PRECISION::pfHis(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
		pfHis = pf
        ! - - - some initial setup - - -
          delta=c/max(taux,1-taux)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
          IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO

        DO l=1,nlam
        pf = pfHis
        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0
        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO
        al = al0 * alf
        ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --------- ----------------------------
        !pf = pfHis
        IF(l>1) THEN
				DO j=1,bn
					g=idx(j)
					bnorm = sqrt(dot_product(b(ix(g):iy(g)),b(ix(g):iy(g))))
					IF  (pf(g)*al*gamm>bnorm)   THEN
							pf(g)=pf(g)-bnorm/(al*gamm)
					ELSE
							pf(g)=0.0D0
					END IF
				ENDDO
        END IF
        ! -------------------------------------
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
          jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE asqr1mcp

! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
SUBROUTINE asqr2mcp (gamm, taux,k,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                        eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::k
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::delta
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
          DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
		DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
		DOUBLE PRECISION::pfHis(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
		pfHis = pf
        ! - - - some initial setup - - -
          delta=k
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
          IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO

        DO l=1,nlam
        
        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0
        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO
        al = al0 * alf
        ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO

        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --------- ----------------------------
        !pf = pfHis
        IF(l>1) THEN
				DO j=1,bn
					g=idx(j)
					bnorm = sqrt(dot_product(b(ix(g):iy(g)),b(ix(g):iy(g))))
					IF  (pf(g)*al*gamm>bnorm)   THEN
							pf(g)=pf(g)-bnorm/(al*gamm)
					ELSE
							pf(g)=0.0D0
					END IF
				ENDDO
		END IF
        ! -------------------------------------
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
          jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE asqr2mcp

! --------------------------------------------------! -------------------------
! -------------------------------------------------! --------------------------
! --------quantile regression  with group scad penalite aproximation-----------
! ---------------------------------------------f1 and f2-----------------------
! --------------------------------------------------! -------------------------
! --------------------------------------------------! -------------------------
SUBROUTINE asqr1scad (gamm, taux,c,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                        eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::c
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::delta
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
          DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
		DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
		DOUBLE PRECISION::pfHis(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
		pfHis = pf
        ! - - - some initial setup - - -
          delta=c/max(taux,1-taux)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
          IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO

        DO l=1,nlam
        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0
        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO
        al = al0 * alf
        ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO

        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --------- ----------------------------
        !pf = pfHis
        IF(l>1) THEN
				DO j=1,bn
					g=idx(j)
					bnorm = sqrt(dot_product(b(ix(g):iy(g)),b(ix(g):iy(g))))
					IF  (pf(g)*al>bnorm)   THEN
							pf(g)=pf(g)
					ELSE IF ((pf(g)*al*gamm>bnorm)) THEN
							pf(g)=(gamm/(gamm-1))*pf(g)-bnorm/(al*(gamm-1))
					ELSE
							pf(g)=0.0D0
					END IF
				ENDDO
		END IF
        ! -------------------------------------
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
         jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE asqr1scad

! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
SUBROUTINE asqr2scad (gamm, taux,k,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,flmin,ulam,&
                        eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)
      ! --------------------------------------------------
        IMPLICIT NONE
      ! - - - arg types - - -
        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION:: flmin
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::k
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::delta
        DOUBLE PRECISION::gamm
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
          DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
		DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
        DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
		DOUBLE PRECISION::pfHis(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
		pfHis = pf
        ! - - - some initial setup - - -
          delta=k
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        IF(flmin < 1.0D0) THEN
        flmin = Max (mfl, flmin)
        alf = flmin ** (1.0D0 / (nlam - 1.0D0))
        ENDIF
        vl = 0.0
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO

        DO l=1,nlam
        
        al0 = al
        IF(flmin>=1.0D0) THEN
        al=ulam(l)
        ELSE
        IF(l > 2) THEN
        al=al*alf
        ELSE IF(l==1) THEN
        al=big
        ELSE IF(l==2) THEN
        al0 = 0.0D0
        DO g = 1,bn
        IF(pf(g)>0.0D0) THEN
        al0 = max(al0, ga(g) / pf(g))
        ENDIF
        END DO
        al = al0 * alf
        ENDIF
        ENDIF
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO

        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --------- ----------------------------
        !pf = pfHis
        IF(l>1) THEN
				DO j=1,bn
					g=idx(j)
					bnorm = sqrt(dot_product(b(ix(g):iy(g)),b(ix(g):iy(g))))
					IF  (pf(g)*al>bnorm)   THEN
							pf(g)=pf(g)
					ELSE IF ((pf(g)*al*gamm>bnorm)) THEN
							pf(g)=(gamm/(gamm-1))*pf(g)-bnorm/(al*(gamm-1))
					ELSE
							pf(g)=0.0D0
					END IF
				ENDDO
		END IF
        ! -------------------------------------
        ! --middle loop-------------------------------------
        DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO

        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u
        unorm=sqrt(dot_product(u,u))
        t=unorm-pf(g)*al
        IF(t>0.0D0) THEN
        b(start:end)=u*t/(gam(g)*unorm)
        ELSE
        b(start:end)=0.0D0
        ENDIF
        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
         jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
        END SUBROUTINE asqr2scad
! --------------------------------------------------! -------------------------
! -------------------------------------------------! --------------------------
! --quantile regression  with group  sparse lasso penalite---------------
! ---------------------------------------------f1 and f2-----------------------
! --------------------------------------------------! -------------------------
! --------------------------------------------------! -------------------------
SUBROUTINE sqr1Slasso (alpha,taux,c,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,ulam,&
         eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)

        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::c
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::alpha
        DOUBLE PRECISION::delta
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER:: n
        INTEGER:: m
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
          DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        DOUBLE PRECISION:: v
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
          delta=c/max(taux,1-taux)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        vl = 0.0
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO

        al=big
        DO l=1,nlam
        al0 = al
        al=ulam(l)
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO

        ! --------- outer loop ----------------------------
        DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u

        unorm=sqrt(dot_product(u,u))
        
        bnorm=sqrt(dot_product(b(start:end),b(start:end)))
        if (bnorm == 0) THEN
              bnorm = bnorm + 1.0E-10
        ENDIF
        
        t=unorm-pf(g)*al*(1-alpha) 
        IF(t>0.0D0) THEN
            n=0
            DO m = start,end
                n=n+1
                v= alpha*al
                v = Abs (u(n)) - v
                IF(v>0.0D0) THEN
                      b(m)=  sign (v, u(n))/(gam(g)+pf(g)*al*(1-alpha)/bnorm)
                ELSE
                      b(m)= 0.0D0
                ENDIF
            ENDDO
        ELSE
                b(start:end)=0.0D0
        ENDIF

        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO

        u=gam(g)*b(start:end) + u

        unorm=sqrt(dot_product(u,u))
        
        bnorm=sqrt(dot_product(b(start:end),b(start:end)))
        if (bnorm == 0) THEN
              bnorm = bnorm + 1.0E-10
        ENDIF
        
        t=unorm-pf(g)*al*(1-alpha) 
        IF(t>0.0D0) THEN
            n=0
            DO m = start,end
                n=n+1
                v= alpha*al
                v = Abs (u(n)) - v
                IF(v>0.0D0) THEN
                      b(m)=  sign (v, u(n))/(gam(g)+pf(g)*al*(1-alpha)/bnorm)
                ELSE
                      b(m)= 0.0D0
                ENDIF
            ENDDO
        ELSE
                b(start:end)=0.0D0
        ENDIF


        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-c)) THEN
        dl (i) = taux - 1.0D0
        ELSE IF (r(i) < 0.0D0) THEN
        dl (i) =  (1.0D0-taux) * r(i) / c
        ELSE IF (r(i) < c) THEN
        dl (i) = taux * r(i) / c
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
          jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv1(taux,c,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
END SUBROUTINE sqr1Slasso

! ----------------------------------------------------------------------------
! ----------------------------------------------------------------------------
SUBROUTINE sqr2Slasso (alpha,taux,k,bn,bs,ix,iy,gam,nobs,nvars,x,y,pf,dfmax,pmax,nlam,ulam,&
         eps,maxit,intr,nalam,b0,beta,idx,nbeta,alam,npass,jerr)

        DOUBLE PRECISION, PARAMETER :: big=9.9E30
        DOUBLE PRECISION, PARAMETER :: mfl = 1.0E-6
        INTEGER, PARAMETER :: mnlam = 6
        INTEGER:: mnl
        INTEGER:: bn
        INTEGER::bs(bn)
        INTEGER::ix(bn)
        INTEGER::iy(bn)
        INTEGER:: nobs
        INTEGER::nvars
        INTEGER::dfmax
        INTEGER::pmax
        INTEGER::nlam
        INTEGER::nalam
        INTEGER::npass
        INTEGER::jerr
        INTEGER::maxit
        INTEGER::intr
        INTEGER::idx(pmax)
        INTEGER::nbeta(nlam)
        DOUBLE PRECISION::eps
        DOUBLE PRECISION::k
        DOUBLE PRECISION::taux
        DOUBLE PRECISION::alpha
        DOUBLE PRECISION::delta
        DOUBLE PRECISION:: x(nobs,nvars)
        DOUBLE PRECISION::y(nobs)
        DOUBLE PRECISION::pf(bn)
        DOUBLE PRECISION::ulam(nlam)
        DOUBLE PRECISION::gam(bn)
        DOUBLE PRECISION:: b0(nlam)
        DOUBLE PRECISION::beta(nvars,nlam)
        DOUBLE PRECISION::alam(nlam)
        ! - - - local declarations - - -
        DOUBLE PRECISION:: max_gam
        DOUBLE PRECISION::d
        DOUBLE PRECISION::t
        DOUBLE PRECISION::dif
        DOUBLE PRECISION::unorm
        DOUBLE PRECISION::bnorm
        DOUBLE PRECISION::al
        DOUBLE PRECISION::alf
        DOUBLE PRECISION::dl(nobs)
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE :: b
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldbeta
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::r
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::oldb
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::u
        DOUBLE PRECISION, DIMENSION (:), ALLOCATABLE ::dd
        INTEGER, DIMENSION (:), ALLOCATABLE :: oidx
        INTEGER:: i
        INTEGER:: n
        INTEGER:: m
        INTEGER::g
        INTEGER::j
        INTEGER::l
        INTEGER::ierr
        INTEGER::ni
        INTEGER::me
        INTEGER::start
        INTEGER::end
        ! - - - local declarations - - -
          DOUBLE PRECISION:: tlam
        INTEGER:: jx
        INTEGER:: jxx(bn)
        DOUBLE PRECISION:: ga(bn)
        DOUBLE PRECISION:: vl(nvars)
        DOUBLE PRECISION:: al0
        DOUBLE PRECISION:: v
        ! - - - begin - - -
          ! - - - allocate variables - - -
          ALLOCATE(b(0:nvars), STAT = jerr)
        ALLOCATE(oldbeta(0:nvars), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(r(1:nobs), STAT = ierr)
        jerr=jerr + ierr
        ALLOCATE(oidx(1:bn), STAT = ierr)
        jerr=jerr + ierr
        IF(jerr /= 0) RETURN
        ! - - - checking pf - - -
          IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
        ENDIF
        pf = max(0.0D0, pf)
        ! - - - some initial setup - - -
          delta=c/max(taux,1-taux)
        jxx = 0
        al = 0.0D0
        mnl = Min (mnlam, nlam)
        r = y
        b = 0.0D0
        oldbeta = 0.0D0
        idx = 0
        oidx = 0
        npass = 0
        ni = npass
        ! --------- lambda loop ----------------------------
        vl = 0.0
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1,bn
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        DEALLOCATE(u)
        END DO

        al=big
        DO l=1,nlam
        al0 = al
        al=ulam(l)
        tlam = (2.0*al-al0)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        IF(ga(g) > pf(g) * tlam) jxx(g) = 1
        ENDDO
        ! --------- outer loop ----------------------------
          DO
        oldbeta(0)=b(0)
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        oldbeta(ix(g):iy(g))=b(ix(g):iy(g))
        ENDDO
        ENDIF
        ! --middle loop-------------------------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO g=1,bn
        IF(jxx(g) == 0) CYCLE
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO
        u=gam(g)*b(start:end) + u

        unorm=sqrt(dot_product(u,u))
        bnorm=sqrt(dot_product(b(start:end),b(start:end)))
        if (bnorm == 0) THEN
              bnorm = bnorm + 1.0E-10
        ENDIF
        t=unorm-pf(g)*al*(1-alpha) 
        IF(t>0.0D0) THEN
            n=0
            DO m = start,end
                n=n+1
                v= alpha*al
                v = Abs (u(n)) - v
                IF(v>0.0D0) THEN
                      b(m)=  sign (v, u(n))/(gam(g)+pf(g)*al*(1-alpha)/bnorm)
                ELSE
                      b(m)= 0.0D0
                ENDIF
            ENDDO
        ELSE
                b(start:end)=0.0D0
        ENDIF

        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        IF(oidx(g)==0) THEN
        ni=ni+1
        IF(ni>pmax) EXIT
        oidx(g)=ni
        idx(ni)=g
        ENDIF
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF (ni > pmax) EXIT
        IF (dif < eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ! --inner loop----------------------
          DO
        npass=npass+1
        dif=0.0D0
        DO j=1,ni
        g=idx(j)
        start=ix(g)
        end=iy(g)
        ALLOCATE(u(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(dd(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        ALLOCATE(oldb(bs(g)),STAT=ierr)
        jerr=jerr+ierr
        IF(jerr/=0) RETURN
        oldb=b(start:end)
        u = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        u = u + dl(i)*x(i,start:end)/nobs
        ENDDO

        u=gam(g)*b(start:end) + u

        unorm=sqrt(dot_product(u,u))
        bnorm=sqrt(dot_product(b(start:end),b(start:end)))
        if (bnorm == 0) THEN
              bnorm = bnorm + 1.0E-10
        ENDIF
        t=unorm-pf(g)*al*(1-alpha) 
        IF(t>0.0D0) THEN
            n=0
            DO m = start,end
                n=n+1
                v= alpha*al
                v = Abs (u(n)) - v
                IF(v>0.0D0) THEN
                      b(m)=  sign (v, u(n))/(gam(g)+pf(g)*al*(1-alpha)/bnorm)
                ELSE
                      b(m)= 0.0D0
                ENDIF
            ENDDO
        ELSE
                b(start:end)=0.0D0
        ENDIF


        dd=b(start:end)-oldb
        IF(any(dd/=0.0D0)) THEN
        dif=max(dif,gam(g)**2*dot_product(dd,dd))
        r=r-matmul(x(:,start:end),dd)
        ENDIF
        DEALLOCATE(u,dd,oldb)
        ENDDO
        IF(intr /= 0) THEN
        d = 0.0D0
        DO i = 1,nobs
        IF (r(i) < (-(1-taux)*k)) THEN
        dl (i) = taux - 1
        ELSE IF (r(i) < (taux)*k) THEN
        dl (i) =  r(i) / k
        ELSE
        dl (i) = taux
        END IF
        d = d + dl(i)
        ENDDO
        d = 0.5 * delta * d / nobs
        IF(d/=0.0D0) THEN
        b(0)=b(0)+d
        r=r-d
        dif=max(dif, d**2)
        ENDIF
        ENDIF
        IF(dif<eps) EXIT
        IF(npass > maxit) THEN
        jerr=-l
        RETURN
        ENDIF
        ENDDO
        ENDDO
        IF(ni>pmax) EXIT
        !--- final check ------------------------
          jx = 0
        max_gam = maxval(gam)
	IF(any((max_gam*(b-oldbeta)/(1+abs(b)))**2 >= eps)) jx = 1
        IF (jx /= 0) CYCLE
        CALL sqrdrv2(taux,k,nobs,nvars,x,r,vl)
        DO g = 1, bn
        IF(jxx(g) == 1) CYCLE
        ALLOCATE(u(bs(g)),STAT=ierr)
        IF(ierr/=0) RETURN
        u = vl(ix(g):iy(g))
        ga(g) = sqrt(dot_product(u,u))
        IF(ga(g) > al*pf(g))THEN
        jxx(g) = 1
        jx = 1
        ENDIF
        DEALLOCATE(u)
        ENDDO
        IF(jx == 1) CYCLE
        EXIT
        ENDDO
        !---------- final update variable and save results------------
          IF(ni>pmax) THEN
        jerr=-10000-l
        EXIT
        ENDIF
        IF(ni>0) THEN
        DO j=1,ni
        g=idx(j)
        beta(ix(g):iy(g),l)=b(ix(g):iy(g))
        ENDDO
        ENDIF
        nbeta(l)=ni
        b0(l)=b(0)
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        me=0
        DO j=1,ni
        g=idx(j)
        IF(any(beta(ix(g):iy(g),l)/=0.0D0)) me=me+1
        ENDDO
        IF(me>dfmax) EXIT
        ENDDO
        DEALLOCATE(b,oldbeta,r,oidx)
        RETURN
END SUBROUTINE sqr2Slasso


