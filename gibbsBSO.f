       subroutine gibbs( iseed )  
       implicit real*8 (a-h,o-z)
       parameter(n=66357,ngmax=110) !capK=5
       real*8 ts(n),dold(n),d(n),x(8,n),xs(8,n),xcure(8,n)
       real*8 beta(5,8),betaC(5,9),alam(5,ngmax)
       real*8 cureD(n,5)
       real*8 bo(5,2)
       real*8 s(5,ngmax)
       integer icount1(5),ntt(5)
       integer iseed,iaccept(n)
       integer ng(5)
c      data
       common  /vects/ts
       common  /vecd/dold ! taking a value of 0,...,6 with 6 denoting unknown cause
       common  /vecx/x
       common  /vecxs/xs
       common  /vecxcure/xcure
       common  /vecs/s
       common  /vecng/ng
c      quantities to be updated
       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
       common  /vbeta/beta ! coefficients of fixed covariates in surv model
       common  /valam/alam ! constant piecewise hazard parameter
       common  /vecbo/bo    ! hyper-parameter for lam
       common  /dummy3/idum3
       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
       common  /vcureD/cureD !! coefficients of fixed covariates in surv model

c......Generate "d" by imputing for unknown causes using Gibbs sampling
       call Gend( iseed )
      write(*,*) 'after Gend'

c......generate thetaR 
c        call GthetaR( iseed )
c      write(*,*) 'after GthetaR'

c......Generate "cureD" by imputing for unknown cure rate status (1 = cured)
       call GencureD( iseed )
      write(*,*) 'after GencureD'

c......generate theta and gamma
c       call Gthetagamma( iseed )
c      write(*,*) 'after Gthetagam'

c......generate Omega
c       call GOmega( iseed )
c      write(*,*) 'after GOmega'

c......generate Sig2 
c       call Gsig2( iseed )
c      write(*,*) 'after Gsig2'

c......generate alpha1
c        call Galpha( iseed )
c      write(*,*) 'after Galpha'

c......generate alpha2
c        call Galpha2( iseed )
c      write(*,*) 'after Galpha2'

c......generate alam
        call Galam( iseed )
      write(*,*) 'after Galam'

c......generate beta (this is betaNC)
        call Gbeta( iseed )
      write(*,*) 'after Gbeta'

c......generate betaC1   !!!
        call GbetaC( iseed )
      write(*,*) 'after GbetaC'

 
       return
       end


       subroutine Gend( iseed )
c      generate censoring indicator for missing causes of failure 
       implicit real*8 (a-h,o-z)
       parameter(n=66357,ngmax=110) !capK=5
       real*8 ts(n),dold(n),d(n)
       real*8 x(8,n),xs(8,n),xcure(8,n)
       real*8 beta(5,8),betaC(5,9),alam(5,ngmax)
       real*8 cureD(n,5)
       real*8 bo(5,2)
       real*8 s(5,ngmax)
       real*8 sum1(5),pf(5),ps(5)
       integer iseed,iaccept(n)
       integer ng(5)
c      data
       common  /vects/ts
       common  /vecd/dold
       common  /vecx/x
       common  /vecxs/xs
       common  /vecxcure/xcure 
       common  /vecs/s
       common  /vecng/ng
c      quantities to be updated
       common  /vecd/d
       common  /vbeta/beta
       common  /valam/alam
       common  /vecbo/bo
       common  /dummy3/idum3
       common  /vbetaC/betaC
       common  /vcureD/cureD
       external  DRNUNF

       do 100 i=1,n
        if (dold(i) .eq. 6.0d0) then 
c........ein=\alpha_k'\theta_i^R+x'\beta_k      
         do 101 k=1,5
          ein=0.0d0
          do j=1,8
           ein=ein+xs(j,i)*beta(k,j)
          enddo
          if (ts(i) .le. s(k,1)) jg=1
          do j1=1,ng(k)-1
           if ( (ts(i) .gt. s(k,j1)) .and.
     1        (ts(i) .le. s(k,j1+1)) ) then
            jg=j1+1
           endif
          enddo
          if (jg .eq. 1) then
           ein1=alam(k,1)*ts(i)
          else
           ein1=alam(k,1)*s(k,1)
           do j2=2,jg-1
            ein1=ein1+alam(k,j2)*(s(k,j2)-s(k,j2-1))
           enddo
           ein1=ein1+alam(k,jg)*(ts(i)-s(k,jg-1))
          endif
          ss=-ein1*dexp(ein)
c.........eincure=x_{cure}'\beta_{k}^c
          eincure=betaC(k,1)
          do j=1,8
           eincure=eincure+xcure(j,i)*betaC(k,j+1)
          enddo
          if(eincure .ge. 0.0d0) then
           einCR=1.0d0/(1.0d0+dexp(-eincure))
          else
           einCR=dexp(eincure)/(1.0d0+dexp(eincure))
          endif
          ps(k)=einCR+(1.0d0-einCR)*dexp(ss) 
          ein2=ss+ein+dlog(alam(k,jg))
          pf(k)=(1.0d0-einCR)*dexp(ein2)
101      continue
c........Start of previous logic for k=2
c         pp=pf(1)*ps(2)/(pf(1)*ps(2)+pf(2)*ps(1)) !!!!!!!!!!!!!!!!up to here
c         call rnset( iseed )
c          u=DRNUNF()
c         call rnget( iseed )
c         if (u .le. pp) then 
c          d(i)=1.0d0 
c         else
c          d(i)=2.0d0
c         endif
c........End of previous logic for k=2
         ppdenom=pf(1)*ps(2)*ps(3)*ps(4)*ps(5)
     1     + pf(2)*ps(1)*ps(3)*ps(4)*ps(5)
     1     + pf(3)*ps(1)*ps(2)*ps(4)*ps(5)
     1     + pf(4)*ps(1)*ps(2)*ps(3)*ps(5)
     1     + pf(5)*ps(1)*ps(2)*ps(3)*ps(4)
         ppnum1=pf(1)*ps(2)*ps(3)*ps(4)*ps(5) 
         ppnum2=pf(2)*ps(1)*ps(3)*ps(4)*ps(5) 
         ppnum3=pf(3)*ps(1)*ps(2)*ps(4)*ps(5) 
         ppnum4=pf(4)*ps(1)*ps(2)*ps(3)*ps(5) 
         !don't need ppnum5
         pp1=ppnum1/ppdenom
         pp2=ppnum2/ppdenom
         pp3=ppnum3/ppdenom
         pp4=ppnum4/ppdenom
         !don't need pp5
         call rnset( iseed )
          u1=DRNUNF()
         call rnget( iseed )
         if (u1 .le. pp1) then 
          d(i)=1.0d0 
         else
          u2=u1-pp1
          if (u2 .le. pp2) then
           d(i)=2.0d0
          else 
           u3=u2-pp2
           if (u3 .le. pp3) then
            d(i)=3.0d0
           else 
            u4=u3-pp3
            if (u4 .le. pp4) then
             d(i)=4.0d0
            else
             d(i)=5.0d0
            endif
           endif
          endif 
         endif
        endif
100    continue  
       return
       end 

       subroutine GencureD( iseed )
c      generate cure indicator (1 = cured)
       implicit real*8 (a-h,o-z)
       parameter(n=66357,ngmax=110) !capK=5
       real*8 ts(n),dold(n),d(n),x(8,n),xs(8,n)
       real*8 xcure(8,n)
       real*8 beta(5,8)
       real*8 betaC(5,9) 
       real*8 cureD(n,5)
       real*8 alam(5,ngmax)
       real*8 bo(5,2)
       real*8 s(5,ngmax)
       real*8 sum1(5),ein
c       real*8 aloghr(2,500)
       real*8 csum(5)
       integer icount1(5),ntt(5)
       integer iseed,iaccept(n)
       integer ng(5)
c      data
       common  /vects/ts
       common  /vecd/dold ! taking a value of 0,...,6 with 6 denoting unknown cause
       common  /vecx/x
       common  /vecxs/xs
       common  /vecxcure/xcure 
       common  /vecs/s
       common  /vecng/ng
c      quantities to be updated
       common  /vecd/d    ! indicator for censoring (0) or cause of death (1,...,5)  
       common  /vbeta/beta ! coefficients of fixed covariates in surv model
       common  /valam/alam ! constant piecewise hazard parameter
       common  /vecbo/bo    ! hyper-parameter for lam
       common  /dummy3/idum3
c       common /valoghr/aloghr
       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
       external  DRNUNF
!!!!!!!!!!!!! Logic     
       ! Only sample when k neq d(i)
       ! General logic of the imputation:
       ! - Step1: first calculate pik = pi / (pi + (1- pi)Sk^*)
       ! - Step2: Generate U ~ U(0,1)
       ! If u ≤ pik then dik = 1, else dik = 0 endif
       
       !Complexity from floating control on calculation of pi
       ! - first calculate linear predictor
       ! - calculate numerator of pi
       ! - calculate Lambda
       ! - calculate pi

!!!!!!!!!!!!!!!!!!!!!!!
       icount=0
       csum(1)=0.0d0
       csum(2)=0.0d0
       csum(3)=0.0d0
       csum(4)=0.0d0
       csum(5)=0.0d0
       do 100 i=1,n !this calculates the linear predictor
        icount=icount+1
        do k=1,5 !This creates the linear predictor for non-cure part
         if (k .eq. int(d(i))) then 
          cureD(i,k) = 0.0d0 ! d(i) = k implies d_ik = 0 w.p.1. since k > 0 
c             (k=int(3.0d0 - d(i))) then ! d(i) ≠ k
         else
          ein=0.0d0 !initializes the value
          eincure=betaC(k,1)
          do j=1,8
           ein=ein+xs(j,i)*beta(k,j)
           eincure = eincure+xcure(j,i)*betaC(k,j+1)
          enddo
c.........calculate lambda by first creating the "s" intervals
          if (ts(i) .le. s(k,1)) jg=1
          do j1=1,ng(k)-1
           if ( (ts(i) .gt. s(k,j1)) .and.
     1          (ts(i) .le. s(k,j1+1)) ) then
            jg=j1+1
           endif
           enddo
           if (jg .eq. 1) then
            ein1=alam(k,1)*ts(i)
           else
            ein1=alam(k,1)*s(k,1)
            do j2=2,jg-1
             ein1=ein1+alam(k,j2)*(s(k,j2)-s(k,j2-1))
            enddo
            ein1=ein1+alam(k,jg)*(ts(i)-s(k,jg-1))
           endif !the final outcome of this "do" is Big Lambda as "ein1"
           ss=-ein1*dexp(ein) !this is log(S_k^*)
c..........calculate denominator like this: 
           if (eincure .le. 0.0d0) then
            ein3 = eincure - dlog(1.0d0 + dexp(eincure)) !this is log(pi)
            ein4 = -dlog(1.0d0 + dexp(eincure)) + ss !this is log((1-pi)Sk^*)
           else
            ein3 = - dlog(1.0d0 + dexp(-eincure)) !this is log(pi)
            ein4 = -eincure - dlog(1.0d0 + dexp(-eincure)) + ss
           endif
c       
           emax = ein3
           if (emax .lt. ein4) emax = ein4
c       
           prob = dexp(ein3 - emax) / (dexp(ein3 - emax) 
     1              + dexp(ein4 - emax)) !this prob is for a particular k
           call rnset( iseed )
            u = DRNUNF() !generate random number from U(0,1)
           call rnget( iseed )
           if (u .le. prob) then
            cureD(i,k) = 1.0d0
           else
            cureD(i,k) = 0.0d0
           endif
          endif
          csum(k)=csum(k)+cureD(i,k)
         enddo
100     continue  
        do k=1,5
         cr=csum(k)/real(icount)
c        write(*,*) 'k=',k,' cr=',cr
        enddo
       return
       end 


c       subroutine Gthetagamma( iseed )
cc      generate censoring indicator for missing causes of failure 
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n), 
c     1        xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7), 
c     1        betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax),Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       real*8 xx(9)
c       real*8 ghat(9),ghat1(9)
c       real*8 Bm(9,9),Bminv(9,9)
c       real*8 tol,Rsig(9,9),RV(1,9)  
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy3/idum3
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c       external DMACH,DCHFAC,DRNMVN,DLINDS 
c       do j1=1,9
c        ghat(j1)=0.0d0
c        ghat1(j1)=0.0d0
c        do j2=1,9
c          Bm(j1,j2)=0.0d0
c        enddo
c        Bm(j1,j1)=0.001d0
c       enddo
c
c       do 100 i=1,n
c         if (m(i) .eq. 1) goto 100
c         do j2=1,m(i)
c          ein3=Omega(1,1)+Omega(2,1)*ty(i,j2) 
c          ein4=Omega(2,2)*ty(i,j2)
c          ein=y(i,j2)-(ein3*thetaR(i,1)+ein4*thetaR(i,2))
c          xx(1)=1.0d0
c          xx(2)=ty(i,j2)
c          do k=1,7
c           k1=k+2
c           xx(k1)=x(k,i)
c          enddo
c          do j1=1,9
c           ghat1(j1)=ghat1(j1)+xx(j1)*ein/sig2
c           do j3=1,9
c             Bm(j1,j3)=Bm(j1,j3)+xx(j1)*xx(j3)/sig2
c           enddo
c          enddo
c         enddo
c100    continue 

cc......compute inverse of ( V1{-1} 
cc.....+ 1/(sigma^2) sum{i=1}^n Z'i Z_i )
c       idim=9
c       call DLINDS(idim,Bm,idim,Bminv,idim)
c
cc......compute conditional mean ghat(3)
c       do j1=1,9
c         ein=0.0d0
c         do j2=1,9
c           ein=ein+Bminv(j1,j2)*ghat1(j2)
c         enddo
c         ghat(j1)=ein
c       enddo
cc......generate multivariate nomal
c       lda=9
c       ldr=9
c       tol=100.0d0*DMACH(4)
c       call DCHFAC(9,Bminv,lda,tol,irank,Rsig,ldr)
cc......generate multivariate normal
c       nr=1
c       call rnset( iseed )
c       call DRNMVN(nr,9,Rsig,9,RV,1)
c       call rnget( iseed )
c       do j=1,9
c        if (j .le. 2) then 
c         theta(j)=RV(1,j)+ghat(j)
c        else
c         j1=j-2
c         gam(j1)=RV(1,j)+ghat(j)
c        endif 
c       enddo
c   
c        return
c       end 
   
c       subroutine GOmega( iseed )
cc      new code
cc      \Gamma=\bgin{bmatrix}
cc             Omega(1,1)  &  0 \\
cc             Omega(2,1)  & Omega(2,2)
cc             \end{bmatrix}
cc      generate Omega(1,1) and Omega(2,2) from truncated normals 
cc      generate Omega(2,1) from a normal distribution
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n), 
c     1        xcure(7,n)
c       real*8 gam(7), theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7), 
c     1        betaC(2,8) !! 
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax),Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       logical la,lb
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecnstar/nstar
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c!	   external  DRNUNF
c       external DRNNOF,YTUVN
cc      generate Omega(2,1)
c       sum1=0.0d0 
c       sum2=0.0005d0 ! prior precision for Omega(2,1) 
c       do 400 i=1,n
c        if (m(i) .eq. 1) goto 400
c        ein1=0.0d0
c        do j1=1,7
c         ein1=ein1+x(j1,i)*gam(j1)
c        enddo
c        do j2=1,m(i)
c         ein3=theta(1)+Omega(1,1)*thetaR(i,1) 
c         ein4=(theta(2)+Omega(2,2)*thetaR(i,2))*ty(i,j2)
c         einy=y(i,j2)-(ein3+ein4+ein1)
c         sum1=sum1+einy*ty(i,j2)*thetaR(i,1)/sig2
c         sum2=sum2+(ty(i,j2)*thetaR(i,1))**2/sig2
c        enddo
c 400   continue
c       amu=sum1/sum2
c       std=1.0d0/dsqrt(sum2)
c       call rnset( iseed )
c         rn=DRNNOF()
c       call rnget( iseed )
c       Omega(2,1)=amu+rn*std
cc
cc      generate Omega(1,1)
c       sum1=0.0d0 
c       sum2=0.0005d0 ! prior precision for Omega(1,1)  
c                    ! using a truncated normal prior 
c       do 401 i=1,n
c        if (m(i) .eq. 1) goto 401
c        ein1=0.0d0
c        do j1=1,7
c         ein1=ein1+x(j1,i)*gam(j1)
c        enddo
c        do j2=1,m(i)
c         ein3=theta(1)+Omega(2,1)*ty(i,j2)*thetaR(i,1) 
c         ein4=(theta(2)+Omega(2,2)*thetaR(i,2))*ty(i,j2)
c         einy=y(i,j2)-(ein3+ein4+ein1)
c         sum1=sum1+einy*thetaR(i,1)/sig2
c         sum2=sum2+(thetaR(i,1))**2/sig2
c        enddo
c 401   continue
c       amu=sum1/sum2
c       std=1.0d0/dsqrt(sum2)
c       lb=.true.
c       aleft=-amu/std
c       call rnset( iseed )
c       tn=YTUVN(aleft,aright,la,lb,iseed)
c       call rnget( iseed )
c       lb=.false.
c       Omega(1,1)=amu+std*tn
cc
cc      generate Omega(2,2)
c       sum1=0.0d0 
c       sum2=0.0005d0 ! prior precision for Omega(2,2)  
c                    ! using a truncated normal prior 
c       do 402 i=1,n
c        if (m(i) .eq. 1) goto 402
c        ein1=0.0d0
c        do j1=1,7
c         ein1=ein1+x(j1,i)*gam(j1)
c        enddo
c        do j2=1,m(i)
c         ein3=theta(1)+(Omega(1,1)+Omega(2,1)*ty(i,j2))*thetaR(i,1) 
c         ein4=theta(2)*ty(i,j2)
c         einy=y(i,j2)-(ein3+ein4+ein1)
c         sum1=sum1+einy*thetaR(i,2)*ty(i,j2)/sig2
c         sum2=sum2+(thetaR(i,2)*ty(i,j2))**2/sig2
c        enddo
c 402   continue
c       amu=sum1/sum2
c       std=1.0d0/dsqrt(sum2)
c       lb=.true.
c       aleft=-amu/std
c       call rnset( iseed )
c       tn=YTUVN(aleft,aright,la,lb,iseed)
c       call rnget( iseed )
c       lb=.false.
c       Omega(2,2)=amu+std*tn
c       ! testing with \Omega=\Gamma \Gamma'
c       ! \Omega2 calculation is not correct with this setting
c!       a11=Omega(1,1)**2
c!       ein5=Omega(2,1)**2+Omega(2,2)**2
c!       ein6=Omega(2,1)*Omega(1,1)
c!       Omega(1,1)=a11
c!       Omega(2,1)=ein5
c!       Omega(1,2)=ein5
c!       Omega(2,2)=ein6
c
c       return
c       end

   
c       subroutine Gsig2( iseed )
cc      generate censoring indicator for missing causes of failure 
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261, mmax=10, ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n)
c       real*8 xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7)
c       real*8 betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2), r(1)
c       real*8 s(2,ngmax)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy3/idum3
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c       external DRNGAM
c       shape1=0.001d0
c       scale1=0.001d0
c       do 100 i=1,n
c        if (m(i) .eq. 1) goto 100
c         ein=real(m(i))
c         shape1=shape1+ein/2.0d0 
c         ein=0.0d0
c         do j=1,7
c          ein=ein+x(j,i)*gam(j)
c         enddo
c         do j=1,m(i)
c          ein3=Omega(1,1)+Omega(2,1)*ty(i,j) 
c          ein4=Omega(2,2)*ty(i,j)
c          ein1=y(i,j)-(theta(1)+ein3*thetaR(i,1)
c     1    +theta(2)*ty(i,j)+ein4*thetaR(i,2)+ein)
c          scale1=scale1+ein1**2/2.0d0
c         enddo
c100    continue
c       call rnset( iseed )
c        call DRNGAM(1,shape1,r)
c       call rnget( iseed )
c       sig2=scale1/r(1)
c        return
c       end 
   
c       subroutine Galpha( iseed )
cc      generate censoring indicator for missing causes of failure 
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       PARAMETER(maxtries=25,im=2,NS=10,np=7,in=1000)
c       real*8  emax, step
cC ------------------------ Working variables ------------------------------
c       real*8  val,a(im+in),ha(im+in),hpa(im+in),
c     1         xlb,xub,rwv(6*NS+15+in)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n)
c       real*8 xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7)
c       real*8 betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       real*8 start(2),xmin(2),ynewlo,reqmin,step1(2)
c       real*8 xiold(2),xinew(2),ximean(2)
c       real*8 tol,Rsig(2,2),RV(1,2)
c       real*8 sig2inv(2,2),sig2inv2(2,2),falphadev2(2,2)
c       integer ifault1,icount
c!       integer iseed
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy3/idum3
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c       external falpha

c        nopt = 2
c        reqmin = 1.0d-10
c        konvge = 5
c        kcount = 1000
c        step1(1)=0.2d0
c        step1(2)=0.2d0
c        start(1)=alpha(1,1)
c        start(2)=alpha(1,2)
c        xiold(1)=alpha(1,1)
c        xiold(2)=alpha(1,2)
c        call nelmin(falpha, nopt, start, xmin, ynewlo, reqmin, step1,
c     1             konvge, kcount, icount, numres, ifault1)
c!       write(*,*) 'xmin=',xmin(1)
c        ximean(1)=xmin(1)
c        ximean(2)=xmin(2)
cc.......calculate Hessian matrix
cc        write(*,*) 'ximean',ximean(j)
c         call alphadev2(ximean(1),ximean(2),falphadev2)
cc        write(*,*) 'alphadev2',falphadev2
c         sig2inv=-falphadev2
c         call DLINDS(2,sig2inv,2,sig2inv2,2)
cc........generate multivariate normal
c         lda=2
c         ldr=2
c         tol=100.0d0*DMACH(4)
c         call DCHFAC(2,sig2inv2,lda,tol,irank,Rsig,ldr)
c
c
c         ein=0.0d0
c         do j1=1,2
c          do j2=1,2
c           ein=ein+(xiold(j2)-ximean(j2))*sig2inv(j2,j1)
c     1            *(xiold(j1)-ximean(j1))/2.0d0
c          enddo
c         enddo
c
c         sumold=-falpha(xiold)+ein
c
c
cc.......perform metropolis
c         do imetr=1,10
c          nr=1
c          call rnset( iseed )
c          call DRNMVN(nr,2,Rsig,2,RV,1)
c          call rnget( iseed )
c          do j=1,2
c           xinew(j)=ximean(j)+RV(1,j)
c          enddo
c
c
c          ein=0.0d0
c          do j1=1,2
c           do j2=1,2
c            ein=ein+(xinew(j2)-ximean(j2))*sig2inv(j2,j1)
c     1             *(xinew(j1)-ximean(j1))/2.0d0
c           enddo
c          enddo
c
c          sumnew=-falpha(xinew)+ein
c          ratio=sumnew-sumold
cc         write(*,*) 'ratio',ratio
c          if (ratio .ge. 0.0d0) then
c           xiold=xinew
c           alpha(1,1)=xinew(1)
c           alpha(1,2)=xinew(2)
c           sumold=sumnew
cc          iaccept=iaccept+1
c          else
c           call rnset( iseed )
c            u=DRNUNF()
c           call rnget( iseed )
c           if (dlog(u) .le. ratio) then
c            xiold=xinew
c            sumold=sumnew
cc           iaccept=iaccept+1
c            alpha(1,1)=xinew(1)
c            alpha(1,2)=xinew(2)
c           endif
c          endif
c         enddo  
c
c       return
c       end 
    
c       real*8 function falpha(z1)
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n)
c       real*8 xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7)
c       real*8 betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       real*8 start(2,2),xmin(2,2),ynewlo,reqmin,step1(2,2)
c       real*8 xiold(2,2),xinew(2,2),ximean(2,2)
c       real*8 tol,Rsig(2,2),RV(1,2)
c       real*8 sig2inv(2,2),sig2inv2(2,2),falphadev2(2,2)
c       real*8 z1(2)
c       integer ifault1,icount
c!       integer iseed
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy3/idum3
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model

c       alpha(1,1)=z1(1)
c       alpha(1,2)=z1(2)
cc......calculate log likelihood
c       sum1=0.0d0
c       do 300 i=1,n
c        if (m(i) .eq. 1) goto 300
c        if (cureD(i,1) .eq. 1.0d0) goto 300
c
c        ein=alpha(1,1)*thetaR(i,1)
c     1       +alpha(1,2)*thetaR(i,2)
c        do j=1,7
c         ein=ein+xs(j,i)*beta(1,j)
c        enddo
c        if (ts(i) .le. s(1,1)) jg=1
c        do j1=1,ng(1)-1
c          if ( (ts(i) .gt. s(1,j1)) .and.
c     1         (ts(i) .le. s(1,j1+1)) ) then
c           jg=j1+1
c          endif
c        enddo
c        if (jg .eq. 1) then
c          ein1=alam(1,1)*ts(i)
c        else
c          ein1=alam(1,1)*s(1,1)
c          do j2=2,jg-1
c           ein1=ein1+alam(1,j2)*(s(1,j2)-s(1,j2-1))
c          enddo
c          ein1=ein1+alam(1,jg)*(ts(i)-s(1,jg-1))
c        endif
c        ss=-ein1*dexp(ein)
c        sum1=sum1+ss
c        if (d(i) .eq. 1.0d0) then
c         sum1=sum1+ein
c        endif
c 300   continue
c       sum1=sum1-0.0001d0*z1(1)**2/2.0d0-0.0001d0*z1(2)**2/2.0d0
c       falpha=-sum1

c       return
c       end

c       subroutine alphadev2(z1,z2,falphadev2)
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n)
c       real*8 xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7)
c       real*8 betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       real*8 start(2,2),xmin(2,2),ynewlo,reqmin,step1(2,2)
c       real*8 xiold(2,2),xinew(2,2),ximean(2,2)
c       real*8 tol,Rsig(2,2),RV(1,2)
c       real*8 sig2inv(2,2),sig2inv2(2,2),falphadev2(2,2)
c       real*8 sum2(2,2),z1,z2
c       integer ifault1,icount
c!       integer iseed
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy3/idum3
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model

c       alpha(1,1)=z1
c       alpha(1,2)=z2
cc......calculate the second derivative
c        do j1=1,2
c         do j2=1,2
c          sum2(j1,j2)=0.0d0
c         enddo
c        enddo
c
c       do 300 i=1,n
c        if (m(i) .eq. 1) goto 300
c        if (cureD(i,1) .eq. 1.0d0) goto 300
c
c        ein=alpha(1,1)*thetaR(i,1)
c     1       +alpha(1,2)*thetaR(i,2)
c        do j=1,7
c         ein=ein+xs(j,i)*beta(1,j)
c        enddo
c        if (ts(i) .le. s(1,1)) jg=1
c        do j1=1,ng(1)-1
c          if ( (ts(i) .gt. s(1,j1)) .and.
c     1         (ts(i) .le. s(1,j1+1)) ) then
c           jg=j1+1
c          endif
c        enddo
c        if (jg .eq. 1) then
c          ein1=alam(1,1)*ts(i)
c        else
c          ein1=alam(1,1)*s(1,1)
c          do j2=2,jg-1
c           ein1=ein1+alam(1,j2)*(s(1,j2)-s(1,j2-1))
c          enddo
c          ein1=ein1+alam(1,jg)*(ts(i)-s(1,jg-1))
c        endif
c        ss=-ein1*dexp(ein)
c        sum2(1,1)=sum2(1,1)+thetaR(i,1)**2*ss
c        sum2(1,2)=sum2(1,2)+thetaR(i,1)*thetaR(i,2)*ss
c        sum2(2,1)=sum2(2,1)+thetaR(i,2)*thetaR(i,1)*ss
c        sum2(2,2)=sum2(2,2)+thetaR(i,2)**2*ss
c 300   continue
c        sum2(1,1)=sum2(1,1)-0.0001d0
c        sum2(2,2)=sum2(2,2)-0.0001d0
c
c       falphadev2=sum2
c
c       return
c       end


c       subroutine Galpha2( iseed )
cc      generate censoring indicator for missing causes of failure 
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       PARAMETER(maxtries=25,im=2,NS=10,np=7,in=1000)
c       real*8  emax, step
cC ------------------------ Working variables ------------------------------
c       real*8  val,a(im+in),ha(im+in),hpa(im+in),
c     1         xlb,xub,rwv(6*NS+15+in)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n)
c       real*8 xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7)
c       real*8 betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       real*8 start(2),xmin(2),ynewlo,reqmin,step1(2)
c       real*8 xiold(2),xinew(2),ximean(2)
c       real*8 tol,Rsig(2,2),RV(1,2)
c       real*8 sig2inv(2,2),sig2inv2(2,2),falphadev2(2,2)
c       integer ifault1,icount
c!       integer iseed
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy3/idum3
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c       external falpha2

c        nopt = 2
c        reqmin = 1.0d-10
c        konvge = 5
c        kcount = 1000
c        step1(1)=0.2d0
c        step1(2)=0.2d0
c        start(1)=alpha(2,1)
c        start(2)=alpha(2,2)
c        xiold(1)=alpha(2,1)
c        xiold(2)=alpha(2,2)
c
c        call nelmin(falpha2, nopt, start, xmin, ynewlo, reqmin, step1,
c     1             konvge, kcount, icount, numres, ifault1)
c!       write(*,*) 'xmin=',xmin(1)
c        ximean(1)=xmin(1)
c        ximean(2)=xmin(2)
cc.......calculate Hessian matrix
cc        write(*,*) 'ximean',ximean(j)
c         call alpha2dev2(ximean(1),ximean(2),falphadev2)
cc        write(*,*) 'alphadev2',falpha2dev2
c         sig2inv=-falphadev2
c
c         call DLINDS(2,sig2inv,2,sig2inv2,2)
cc........generate multivariate normal
c         lda=2
c         ldr=2
c         tol=100.0d0*DMACH(4)
c         call DCHFAC(2,sig2inv2,lda,tol,irank,Rsig,ldr)
c
c
c         ein=0.0d0
c         do j1=1,2
c          do j2=1,2
c           ein=ein+(xiold(j2)-ximean(j2))*sig2inv(j2,j1)
c     1            *(xiold(j1)-ximean(j1))/2.0d0
c          enddo
c         enddo
c
c         sumold=-falpha2(xiold)+ein
c
c
cc.......perform metropolis
c         do imetr=1,10
c          nr=1
c          call rnset( iseed )
c          call DRNMVN(nr,2,Rsig,2,RV,1)
c          call rnget( iseed )
c          do j=1,2
c           xinew(j)=ximean(j)+RV(1,j)
c          enddo
c
c
c          ein=0.0d0
c          do j1=1,2
c           do j2=1,2
c            ein=ein+(xinew(j2)-ximean(j2))*sig2inv(j2,j1)
c     1             *(xinew(j1)-ximean(j1))/2.0d0
c           enddo
c          enddo
c
c          sumnew=-falpha2(xinew)+ein
c          ratio=sumnew-sumold
cc         write(*,*) 'ratio',ratio
c          if (ratio .ge. 0.0d0) then
c           xiold=xinew
c           alpha(2,1)=xinew(1)
c           alpha(2,2)=xinew(2)
c           sumold=sumnew
cc          iaccept=iaccept+1
c          else
c           call rnset( iseed )
c            u=DRNUNF()
c           call rnget( iseed )
c           if (dlog(u) .le. ratio) then
c            xiold=xinew
c            sumold=sumnew
cc           iaccept=iaccept+1
c            alpha(2,1)=xinew(1)
c            alpha(2,2)=xinew(2)
c           endif
c          endif
c         enddo  
c
c       return
c       end 
    
c       real*8 function falpha2(z1)
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n)
c       real*8 xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7)
c       real*8  betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       real*8 start(2,2),xmin(2,2),ynewlo,reqmin,step1(2,2)
c       real*8 xiold(2,2),xinew(2,2),ximean(2,2)
c       real*8 tol,Rsig(2,2),RV(1,2)
c       real*8 z1(2)
c       real*8 sig2inv(2,2),sig2inv2(2,2),falpha2dev2(2,2)
c       integer ifault1,icount
c!       integer iseed
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy3/idum3
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model

c       alpha(2,1)=z1(1)
c       alpha(2,2)=z1(2)
cc......calculate log likelihood
c       sum1=0.0d0
c       do 300 i=1,n
c        if (m(i) .eq. 1) goto 300
c        if (cureD(i,2) .eq. 1.0d0) goto 300
c        ein=alpha(2,1)*thetaR(i,1)
c     1       +alpha(2,2)*thetaR(i,2)
c        do j=1,7
c         ein=ein+xs(j,i)*beta(2,j)
c        enddo
c        if (ts(i) .le. s(2,1)) jg=1
c        do j1=1,ng(2)-1
c          if ( (ts(i) .gt. s(2,j1)) .and.
c     1         (ts(i) .le. s(2,j1+1)) ) then
c           jg=j1+1
c          endif
c        enddo
c        if (jg .eq. 1) then
c          ein1=alam(2,1)*ts(i)
c        else
c          ein1=alam(2,1)*s(2,1)
c          do j2=2,jg-1
c           ein1=ein1+alam(2,j2)*(s(2,j2)-s(2,j2-1))
c          enddo
c          ein1=ein1+alam(2,jg)*(ts(i)-s(2,jg-1))
c        endif
c        ss=-ein1*dexp(ein)
c        sum1=sum1+ss
c        if (d(i) .eq. 2.0d0) then
c         sum1=sum1+ein
c        endif
c 300   continue
c       sum1=sum1-0.0001d0*z1(1)**2/2.0d0-0.0001d0*z1(2)**2/2.0d0
c       falpha2=-sum1
c
c       return
c       end

c       subroutine alpha2dev2(z1,z2,falpha2dev2)
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n)
c       real*8 xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7)
c       real*8 betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       real*8 start(2,2),xmin(2,2),ynewlo,reqmin,step1(2,2)
c       real*8 xiold(2,2),xinew(2,2),ximean(2,2)
c       real*8 tol,Rsig(2,2),RV(1,2)
c       real*8 sig2inv(2,2),sig2inv2(2,2),falpha2dev2(2,2)
c       real*8 sum2(2,2),z1,z2
c       integer ifault1,icount
c!       integer iseed
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy3/idum3
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model

c       alpha(2,1)=z1
c       alpha(2,2)=z2
cc......calculate the second derivative
c        do j1=1,2
c         do j2=1,2
c          sum2(j1,j2)=0.0d0
c         enddo
c        enddo
c
c       do 300 i=1,n
c        if (m(i) .eq. 1) goto 300
c        if (cureD(i,2) .eq. 1.0d0) goto 300
cc
c        ein=alpha(2,1)*thetaR(i,1)
c     1       +alpha(2,2)*thetaR(i,2)
c        do j=1,7
c         ein=ein+xs(j,i)*beta(2,j)
c        enddo
c        if (ts(i) .le. s(2,1)) jg=1
c        do j1=1,ng(2)-1
c          if ( (ts(i) .gt. s(2,j1)) .and.
c     1         (ts(i) .le. s(2,j1+1)) ) then
c           jg=j1+1
c          endif
c        enddo
c        if (jg .eq. 1) then
c          ein1=alam(2,1)*ts(i)
c        else
c          ein1=alam(2,1)*s(2,1)
c          do j2=2,jg-1
c           ein1=ein1+alam(2,j2)*(s(2,j2)-s(2,j2-1))
c          enddo
c          ein1=ein1+alam(2,jg)*(ts(i)-s(2,jg-1))
c        endif
c        ss=-ein1*dexp(ein)
c        sum2(1,1)=sum2(1,1)+thetaR(i,1)**2*ss
c        sum2(1,2)=sum2(1,2)+thetaR(i,1)*thetaR(i,2)*ss
c        sum2(2,1)=sum2(2,1)+thetaR(i,2)*thetaR(i,1)*ss
c        sum2(2,2)=sum2(2,2)+thetaR(i,2)**2*ss
c 300   continue
c        sum2(1,1)=sum2(1,1)-0.0001d0
c        sum2(2,2)=sum2(2,2)-0.0001d0
c
c       falpha2dev2=sum2
c
c       return
c       end

   
       subroutine Gbeta( iseed )
c      generate censoring indicator for missing causes of failure 
       implicit real*8 (a-h,o-z)
       parameter(n=66357, ngmax=110) !capK=5
       PARAMETER(maxtries=25,im=2,NS=10,in=1000)
       real*8 emax,step
C ------------------------ Working variables ------------------------------
       real*8 val,a(im+in),ha(im+in),hpa(im+in),
     1         xlb,xub,rwv(6*NS+15+in)
       real*8 ts(n),dold(n),d(n),x(8,n),xs(8,n),xcure(8,n)
       real*8 beta(5,8),betaC(5,9) 
       real*8 cureD(n,5)
       real*8 alam(5,ngmax)
       real*8 bo(5,2)
       real*8 s(5,ngmax)
       integer iaccept(n)
       integer ng(5)
       integer ifault,iwv(NS+7+in),ifcount
       integer iseed,idum2
c      data
       common  /vects/ts
       common  /vecd/dold ! taking a value of 0,...,6 with 6 denoting unknown cause
       common  /vecx/x
       common  /vecxs/xs
       common  /vecxcure/xcure
       common  /vecs/s
       common  /vecng/ng
c      quantities to be updated
       common  /vecd/d    ! indicator for censoring (0) or cause of death (1,...,5)  
       common  /vbeta/beta ! coefficients of fixed covariates in surv model
       common  /valam/alam ! constant piecewise hazard parameter
       common  /vecbo/bo    ! hyper-parameter for lam
       common  /dummy1/idum1
       common  /dummy2/idum2
       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c      external  DRNUNF
       external eval
       emax=60.0d0
       do 79 k=1,5
              idum1=k
       DO 80 j=1,8
        idum2=j
        step=1.0d0
c       step=0.5d0
        a(1)=beta(k,j)-step
        a(2)=beta(k,j)+step
20      CALL eval(a(1), ha(1), hpa(1))
        CALL eval(a(2), ha(2), hpa(2))
c        write(*,*) 'a(1)=',a(1),' ha(1)=',ha(1)
c        write(*,*) ' hpa(1)=',hpa(1)
c ------------------ Get starting points for Gilks -----------------------
         IF (hpa(2) .gt. 0.0) THEN
c        Both points to the left of the mode; push a(2) to the right
30        a(2)=a(2)+step
          CALL EVAL(a(2), ha(2), hpa(2))
c         If a(2) still to the left of the mode repeat, else continue
          IF (hpa(2).gt.0.0) THEN
           GOTO 30
          ELSE
           a(1)=a(2)-step
          ENDIF
         ELSE
          IF (hpa(1) .le. 0.0) THEN
c        Both points to the right of the mode
c        Insert new point to the left of a(1)
           a(2)=a(1)
           ha(2)=ha(1)
           hpa(2)=hpa(1)
40         a(1)=a(1)-step
           CALL EVAL(a(1), ha(1), hpa(1))
           IF (hpa(1) .le. 0.0) THEN
            GOTO 40
c          ELSE
c           a(2)=a(1)+step
           ENDIF
          ENDIF
         ENDIF
c -----------------------------------------------------------------------
         ifcount=0
50       CALL INITIAL(ns, im, emax, a, ha, hpa, lb, xlb, ub, xub,
     +                ifault, iwv, rwv, eval)
         IF((ifault.eq.3) .or. (ifault.eq.4)) GOTO 20
C
         CALL SAMPLE(iwv,rwv,eval,val,ifault,iseed)
         IF(ifault.eq.5) THEN
c          PRINT*, 'ifault =', ifault, ' trying again'
           a(1)=a(1)-step
           a(2)=a(2)+step
           ifcount=ifcount+1
           IF(ifcount.le.maxtries) THEN
             GOTO 50
           ELSE
c             PRINT*, 'Too many tries resetting starting values'
              a(1)=0.0
              a(2)=0.0
              step=1.2d0*step
              GOTO 20
           ENDIF
         ENDIF
         IF(ifault.eq.7) THEN
           a(1)=a(1)-step
           a(2)=a(2)+step
           ifcount=ifcount+1
           IF(ifcount.le.maxtries) THEN
             GOTO 50
           ELSE
c            PRINT*, 'Too many tries resetting starting values'
             a(1)=0.0
             a(2)=0.0
             step=1.2d0*step
             GOTO 20
           ENDIF
         ENDIF
         IF (IFAULT .NE. 0) THEN
c           PRINT*, ' The sampling status is ', IFAULT
         ENDIF
C
            beta(k,j)=val
C
80     continue
79     continue
   
       return
       end 
              

       subroutine eval(z1,ahz,ahpz)
       implicit real*8 (a-h,o-z)
       parameter(n=66357, ngmax=110) !capK=5
       real*8 ts(n),dold(n),d(n),x(8,n),xs(8,n),xcure(8,n)
       real*8 beta(5,8),betaC(5,9)  
       real*8 cureD(n,5) 
       real*8 alam(5,ngmax)
       real*8 bo(5,2)
       real*8 s(5,ngmax)
       real*8 z1,ahz,ahpz
       integer iseed,iaccept(n)
       integer ng(5)
       common  /vects/ts
       common  /vecd/dold ! taking a value of 0,...,6 with 6 denoting unknown cause
       common  /vecx/x
       common  /vecxs/xs
       common  /vecxcure/xcure
       common  /vecs/s
       common  /vecng/ng
c      quantities to be updated
       common  /vecd/d    ! indicator for censoring (0) or cause of death (1,...,5)  
       common  /vbeta/beta ! coefficients of fixed covariates in surv model
       common  /valam/alam ! constant piecewise hazard parameter
       common  /vecbo/bo    ! hyper-parameter for lam
       common  /dummy1/idum1
       common  /dummy2/idum2
       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
       external  DRNUNF
       k=idum1
       beta(k,idum2)=z1
c......calculate log likelihood and its derivative
       sum1=0.0d0
       sum2=0.0d0
       do 300 i=1,n
        if (cureD(i,k) .eq. 1.0d0) goto 300
        ein=0.0d0
        do j=1,8
         ein=ein+xs(j,i)*beta(k,j)
        enddo
        if (ts(i) .le. s(k,1)) jg=1
        do j1=1,ng(k)-1
          if ( (ts(i) .gt. s(k,j1)) .and.
     1         (ts(i) .le. s(k,j1+1)) ) then
           jg=j1+1
          endif
        enddo
        if (jg .eq. 1) then
          ein1=alam(k,1)*ts(i)
        else
          ein1=alam(k,1)*s(k,1)
          do j2=2,jg-1
           ein1=ein1+alam(k,j2)*(s(k,j2)-s(k,j2-1))
          enddo
          ein1=ein1+alam(k,jg)*(ts(i)-s(k,jg-1))
        endif
        ss=-ein1*dexp(ein)
        sum1=sum1+ss
        sum2=sum2+xs(idum2,i)*ss
        if (d(i) .eq. k) then
         sum1=sum1+ein
         sum2=sum2+xs(idum2,i)
        endif
 300   continue
       sum1=sum1-0.001d0*z1**2/2.0d0
       sum2=sum2-0.001d0*z1 !derivative
       ahz=sum1
       ahpz=sum2   
       return
       end

c       subroutine Gbeta2( iseed )
c      generate censoring indicator for missing causes of failure 
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261, mmax=10, ngmax=110)
c       PARAMETER(maxtries=25,im=2,NS=10,in=1000)
c       real*8  emax, step
C ------------------------ Working variables ------------------------------
c       real*8  val,a(im+in),ha(im+in),hpa(im+in),
c     1         xlb,xub,rwv(6*NS+15+in)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n), 
c     1        xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7), 
c     1        betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax),Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       real*8 start(1),xmin(1),ynewlo,reqmin,step1
c       integer iaccept(n)
c       integer ifault,iwv(NS+7+in),ifcount
c       integer iseed,idum2
c       integer m(n),ng(2)
c      data
c       common  /vecy/y
c       common  /vecty/ty
c      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
c      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy2/idum2
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c      external  fnbeta2

c       external eval2
c       emax=60.0d0
c       DO 80 j=1,7
c        idum2=j
c		find the mode of the conditional posterior distribution 
c       nopt = 1
c       reqmin = 1.0d-10
c       konvge = 5
c       kcount = 1000
c       step1=0.2d0
c       start(1)=beta(2,j)
c       call nelmin(fnbeta2, nopt, start, xmin, ynewlo, reqmin, step1,
c    1             konvge, kcount, icount, numres, ifault1)
c       write(*,*) 'j=',j, ' xmin=',xmin(1),' beta2=',beta(2,j)
c        ximean=xmin(1)				
c        step=1.0d0
c       step=0.5d0
c        xmin(1)=beta(2,j)
c        a(1)=xmin(1)-step
c        a(2)=xmin(1)+step
c20      CALL eval2(a(1), ha(1), hpa(1))
c        CALL eval2(a(2), ha(2), hpa(2))
c        write(*,*) 'a(1)=',a(1),' ha(1)=',ha(1)
c        write(*,*) ' hpa(1)=',hpa(1)
c ------------------ Get starting points for Gilks -----------------------
c         IF (hpa(2) .gt. 0.0) THEN
c        Both points to the left of the mode; push a(2) to the right
c30          a(2)=a(2)+step
c            CALL EVAL2(a(2), ha(2), hpa(2))
c        If a(2) still to the left of the mode repeat, else continue
c            IF (hpa(2).gt.0.0) THEN
c               GOTO 30
c            ELSE
c               a(1)=a(2)-step
c            ENDIF
c         ELSE
c            IF (hpa(1) .le. 0.0) THEN
c        Both points to the right of the mode
c        Insert new point to the left of a(1)
c               a(2)=a(1)
c               ha(2)=ha(1)
c               hpa(2)=hpa(1)
c40             a(1)=a(1)-step
c               CALL EVAL2(a(1), ha(1), hpa(1))
c               IF (hpa(1) .le. 0.0) THEN
c                  GOTO 40
c              ELSE
c                 a(2)=a(1)+step
c               ENDIF
c             ENDIF
c         ENDIF
c -----------------------------------------------------------------------
c         ifcount=0
c50       CALL INITIAL(ns, im, emax, a, ha, hpa, lb, xlb, ub, xub,
c     +                ifault, iwv, rwv, eval2)
c         IF((ifault.eq.3) .or. (ifault.eq.4)) GOTO 20
C
c         CALL SAMPLE(iwv,rwv,eval2,val,ifault,iseed)
c         IF(ifault.eq.5) THEN
c          PRINT*, 'ifault =', ifault, ' trying again'
c           a(1)=a(1)-step
c           a(2)=a(2)+step
c           ifcount=ifcount+1
c           IF(ifcount.le.maxtries) THEN
c             GOTO 50
c           ELSE
c             PRINT*, 'Too many tries resetting starting values'
c              a(1)=0.0
c              a(2)=0.0
c              step=1.2d0*step
c              GOTO 20
c           ENDIF
c         ENDIF
c         IF(ifault.eq.7) THEN
c           a(1)=a(1)-step
c           a(2)=a(2)+step
c           ifcount=ifcount+1
c           IF(ifcount.le.maxtries) THEN
c             GOTO 50
c           ELSE
c            PRINT*, 'Too many tries resetting starting values'
c             a(1)=0.0
c             a(2)=0.0
c             step=1.2d0*step
c             GOTO 20
c           ENDIF
c         ENDIF
c         IF (IFAULT .NE. 0) THEN
c           PRINT*, ' The sampling status is ', IFAULT
c         ENDIF
C
c            beta(2,j)=val
C
c 80    continue
   
c       return
c       end 

c       subroutine eval2(z1,ahz,ahpz)
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261, mmax=10, ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n), 
c     1        xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7), 
c     1       betaC(2,8) !! 
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax),Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       real*8 z1,ahz,ahpz
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       common  /vecy/y
c       common  /vecty/ty
c      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /xcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
c      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy2/idum2
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c       external  DRNUNF
c       beta(2,idum2)=z1
c......calculate log likelihood and its derivative
c       sum1=0.0d0
c       sum2=0.0d0
c       do 300 i=1,n
c        if (m(i) .eq. 1) goto 300
c        if (cureD(i,2) .eq. 1.0d0) goto 300
c        ein=alpha(2,1)*thetaR(i,1)
c     1       +alpha(2,2)*thetaR(i,2)
c        do j=1,7
c         ein=ein+xs(j,i)*beta(2,j)
c        enddo
c        if (ts(i) .le. s(2,1)) jg=1
c        do j1=1,ng(2)-1
c          if ( (ts(i) .gt. s(2,j1)) .and.
c     1         (ts(i) .le. s(2,j1+1)) ) then
c           jg=j1+1
c          endif
c        enddo
c        if (jg .eq. 1) then
c          ein1=alam(2,1)*ts(i)
c        else
c          ein1=alam(2,1)*s(2,1)
c          do j2=2,jg-1
c           ein1=ein1+alam(2,j2)*(s(2,j2)-s(2,j2-1))
c          enddo
c          ein1=ein1+alam(2,jg)*(ts(i)-s(2,jg-1))
c        endif
c        ss=-ein1*dexp(ein)
c        sum1=sum1+ss
c        sum2=sum2+xs(idum2,i)*ss
c        if (d(i) .eq. 2.0d0) then
c         sum1=sum1+ein
c         sum2=sum2+xs(idum2,i)
c        endif
c 300   continue
c       sum1=sum1-0.001d0*z1**2/2.0d0
c       sum2=sum2-0.001d0*z1
c       ahz=sum1
c       ahpz=sum2   
c       return
c       end

       subroutine GbetaC( iseed )
c      generate censoring indicator for missing causes of failure 
       implicit real*8 (a-h,o-z)
       parameter(n=66357, ngmax=110) !capK=5
       real*8 step
C ------------------------ Working variables ------------------------------
       real*8 ts(n),dold(n),d(n),x(8,n),xs(8,n),xcure(8,n)
       real*8 beta(5,8),betaC(5,9) 
       real*8 cureD(n,5) 
       real*8 alam(5,ngmax)
       real*8 bo(5,2)
       real*8 s(5,ngmax)
       real*8 start(1),xmin(1),step1(1)
       real*8 xa(5), ya(5)
       integer iaccept(n)
       integer ng(5)
       integer ifault,ifcount
       integer iseed,idum2
c      data
       common  /vects/ts
       common  /vecd/dold ! taking a value of 0,...,6 with 6 denoting unknown cause
       common  /vecx/x
       common  /vecxs/xs
       common  /vecxcure/xcure
       common  /vecs/s
       common  /vecng/ng
c      quantities to be updated
       common  /vecd/d    ! indicator for censoring (0) or cause of death (1,...,5)  
       common  /vbeta/beta ! coefficients of fixed covariates in surv model
       common  /valam/alam ! constant piecewise hazard parameter
       common  /vecbo/bo    ! hyper-parameter for lam
       common  /dummy1C/idumC1
       common  /dummy2C/idumC2
       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c      external  DRNUNF
        external fbetaC!,evalC
!Step 1: Let Bj be the current value
c       betacurrent=0.0d0 
!Step 2A: Find the mode of the log of the conditional posterior
       do 79 k=1,5
            idumC1=k
       DO 80 j=1,9
              idumC2=j
              step=0.5d0
         nopt = 1
         reqmin = 1.0d-10
         konvge = 5
         kcount = 1000
         step1(1)=0.2d0
         start(1)=betaC(k,j)
		 betacurrent=betaC(k,j)
         call nelmin(fbetaC, nopt, start, xmin, ynewlo, reqmin, step1,
     1         konvge, kcount, icount, numres, ifault1) ! need to know the function 
C
         betaCmode=xmin(1) !xmin(1) is mode
         amode=betaCmode
         if (amode .lt. 0.0d0) amode=-betaCmode !abs value
         stepa=0.30d0*amode 
!Step 2B: Use quadratic regression to compute "a"
         xa(1)=betaCmode
         xa(2)=betaCmode-stepa
         xa(3)=betaCmode+stepa
         xa(4)=betaCmode-2.0d0*stepa
         xa(5)=betaCmode+2.0d0*stepa
         ya(1)= -fbetaC(xa(1))
         ya(2)= -fbetaC(xa(2))
         ya(3)= -fbetaC(xa(3))
         ya(4)= -fbetaC(xa(4))
         ya(5)= -fbetaC(xa(5))
c         Ma=0.0d0
c         Mb=0.0d0
c         Mc=0.0d0
c         Mf=0.0d0
c         Mi=0.0d0
c         respy=0.0d0
c         respxy=0.0d0
c         respxxy=0.0d0
c         do ibeta=1,5
c          Ma=Ma+1.0d0  !Ma=5
c          Mb=Mb+xa(ibeta) !sum of x (Md=Mb)
c          Mc=Mc+xa(ibeta)**2 !sum of x^2 (Me=Mc=Mg)
c          Mf=Mf+xa(ibeta)**3 !sum of x^3 (Mh=Mf)
c          Mi=Mi+xa(ibeta)**4 !sum of x^4
c          respy=respy+ya(ibeta) !sum of y
c          respxy=respxy+xa(ibeta)*ya(ibeta) !sum of xy
c          respxxy=respxxy+(xa(ibeta)**2)*ya(ibeta) !sum of (x^2)y
c         enddo
c         aestnum=(Mb*Mf-Mc**2)*respy+(Mb*Mc-Ma*Mf)*respxy
c     1                +(Ma*Mc-Mb**2)*respxxy
c         aestdenom=Ma*(Mc*Mi-Mf**2)-Mb*(Mb*Mi-Mf*Mc)+Mc*(Mb*Mf-Mc**2)
c         aest=aestnum/aestdenom

         xsum=0.0d0
         xxsum=0.0d0
         ysum=0.0d0
         do imean=1,5
          xsum=xsum+xa(imean) !sum of x
          xxsum=xxsum+xa(imean)**2  !sum of x^2
          ysum=ysum+ya(imean) !sum of y
         enddo

         xmean=xsum/5.0d0
         xxmean=xxsum/5.0d0
         ymean=ysum/5.0d0

         Sxx=0.0d0 !S_{X,X}
         Sxsqy=0.0d0 !S_{X^2,Y}
         Sxy=0.0d0 !S_{X,Y}
         Sxxsq=0.0d0 !S_{X,X^2} 
         Sxsqxsq=0.0d0 !S_{X^2,X^2}
         do ibeta=1,5
          Sxx=Sxx+(xa(ibeta)-xmean)**2
          Sxsqy=Sxsqy+(xa(ibeta)**2-xxmean)*(ya(ibeta)-ymean)
          Sxy=Sxy+(xa(ibeta)-xmean)*(ya(ibeta)-ymean)
          Sxxsq=Sxxsq+(xa(ibeta)-xmean)*(xa(ibeta)**2-xxmean)
          Sxsqxsq=Sxsqxsq+(xa(ibeta)**2-xxmean)**2
         enddo
         aestnum=Sxx*Sxsqy-Sxy*Sxxsq
         aestdenom=Sxx*Sxsqxsq-Sxxsq**2
         aest=aestnum/aestdenom

!Step 3: Draw a new beta (beta*) value from N(mode, -1/2a)

         sd=dsqrt(-1.0d0/(2.0d0*aest))
         rold=-fbetaC(betacurrent)-aest*(betacurrent-betaCmode)**2 
         do imetro=1,10		  
         call rnset( iseed ) !call RNSET to initialize the seed
           betanew=DRNNOF() !generates standard normal
         call rnget( iseed )
         betanew=betanew*sd+betaCmode !transform(betanew=beta*)
!Step 4: Move from current to betanew with probability alpha
          rnew=-fbetaC(betanew)-aest*(betanew-betaCmode)**2
         ratioC=rnew-rold
c		 
c		 -fbetaC(betanew) !fbetaC=-log(cond post)
c     1     +aest*(betacurrent-betaCmode)**2
c     1     +fbetaC(betacurrent)
c     1     -aest*(betanew-betaCmode)**2
         if (ratioC .ge. 0.0d0) then ! case where alpha=1 (auto move to new)
              betacurrent=betanew
			  rold=rnew
           betaC(k,j)=betanew
         else
          call rnset( iseed )
           u=DRNUNF() !generate random number from U(0,1)
          call rnget( iseed )
           if (dlog(u) .le. ratioC) then !if < log(alpha) then use new val
              betacurrent=betanew
			  rold=rnew
            betaC(k,j)=betanew
           endif
         endif
		 enddo
         write(*,*) 'betaC(',k,',',j,')=',betaC(k,j)
         
c            !Step 3: draw beta star from normal based on mode
c             !Box-Muller transformation used
c              call rnset( iseed )
c               u1 = DRNUNF() !generate random number from U(0,1)
c              call rnget( iseed )
c              call rnset( iseed )
c               u2 = DRNUNF() !generate random number from U(0,1)
c              call rnget( iseed )
c              PI=4.D0*DATAN(1.D0) !this is pi with maximum precision
c              bstar=(sqrt(-2.0d0*dlog(u1))*cos(2.0d0*PI*u2))
c     1         *sqrt(-1.0d0/(2.0d0*aest))+betaC(k,j)
              

C
80        continue
79        continue
                 
       return
       end 
                            
       real*8 function fbetaC(z1)
       implicit real*8 (a-h,o-z)
       parameter(n=66357, ngmax=110) !capK=5
       real*8 ts(n),dold(n),d(n),x(8,n),xs(8,n),xcure(8,n)
       real*8 beta(5,8),betaC(5,9) 
       real*8 cureD(n,5) 
       real*8 alam(5,ngmax)
       real*8 bo(5,2)
       real*8 s(5,ngmax)
       real*8 z1,ahz,ahpz
       real*8 xx(9)
       real*8 A(9,9), detm1, detm2, fact(9,9)

       integer iseed,iaccept(n)
       integer ng(5)
       common  /vects/ts
       common  /vecd/dold ! taking a value of 0,...,6 with 6 denoting unknown cause
       common  /vecx/x
       common  /vecxs/xs
       common  /vecxcure/xcure
       common  /vecs/s
       common  /vecng/ng
c      quantities to be updated
       common  /vecd/d    ! indicator for censoring (0) or cause of death (1,...,5)  
       common  /vbeta/beta ! coefficients of fixed covariates in surv model
       common  /valam/alam ! constant piecewise hazard parameter
       common  /vecbo/bo    ! hyper-parameter for lam
       common  /dummy1C/idumC1
       common  /dummy2C/idumC2
       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
       external  DRNUNF
        k=idumC1
        betaC(k,idumC2)=z1
c......calculate log likelihood and its derivative
       !initialize the matrix for prior
       sum1=0.0d0
       do k1a = 1,9
        do k2a=1,9
         A(k1a,k2a)=0.0d0
        enddo
       enddo
       !calculate log likelihood
       do 300 i=1,n
        xx(1)=1.0d0
        do j=2,9
         xx(j)=xcure(j-1,i)
        enddo
        ein=0.0d0
        do j=1,9
         ein=ein+betaC(k,j)*xx(j)
        enddo
		
		if (ein .le. 0.0d0) then
         wi=dexp(ein)/((1.0d0+dexp(ein))**2) !part of prior
        else
         wi=dexp(-ein)/((1.0d0+dexp(-ein))**2) !part of prior
        endif
		
		do k1a=1,9
          do k2a=1,9
           A(k1a,k2a)=A(k1a,k2a)+xx(k1a)*wi*xx(k2a)
          enddo
        enddo
	     
		ak=real(k)
c		if (d(i) .eq. ak) then 
        sum1=sum1+ein*cureD(i,k) 
c        j=idumC2
         if (ein .le. 0.0d0) then
          sum1=sum1-dlog(1.0d0+dexp(ein)) !log likelihood
         else
          sum1=sum1-ein-dlog(1.0d0+dexp(-ein)) !log likelihood
         endif
c	    else
c         eink=alpha(k,1)*thetaR(i,1)
c    1       +alpha(k,2)*thetaR(i,2)
c         do j=1,7
c          eink=eink+xs(j,i)*beta(k,j)
c         enddo
c         if (ts(i) .le. s(k,1)) jg=1
c         do j1=1,ng(k)-1
c          if ( (ts(i) .gt. s(k,j1)) .and.
c     1         (ts(i) .le. s(k,j1+1)) ) then
c           jg=j1+1
c          endif
c         enddo
c         if (jg .eq. 1) then
c          ein1=alam(k,1)*ts(i)
c         else
c          ein1=alam(k,1)*s(k,1)
c          do j2=2,jg-1
c           ein1=ein1+alam(k,j2)*(s(k,j2)-s(k,j2-1))
c          enddo
c          ein1=ein1+alam(k,jg)*(ts(i)-s(k,jg-1))
c         endif
c         ss=-ein1*dexp(eink)
c		 if (ein .le. 0.0d0) then
c           eini=dexp(ein)/((1.0d0+dexp(ein))**2) !part of prior
c         else
c          eini=dexp(-ein)/((1.0d0+dexp(-ein))**2) !part of prior
c         endif
c 		 sum1=sum1+dlog(eini+(1.0d0-eini)*dexp(ss))
c		endif
300    continue
c       if (idumC2 .eq. 1) then
c        sum1=sum1-z1**2/(2.0d0*1.0d0) !Jeffreys prior
c       else
c        sum1=sum1-z1**2/(2.0d0*3.0d0) !change
c       endif
       call DLFTDS (9,A,9,fact,9) !nrow, matrix, ncol, output of factored matrix, nrow of fact
       call DLFDDS(9,fact,9,detm1,detm2) !outputs determinant
       sum1=sum1+0.5d0*dlog(detm1*10.0d0**detm2) !log of Jeffreys-type prior (this time it's plus!)
       fbetaC=-sum1
       return
       end

              
       !subroutine evalC(z1,ahz,ahpz)
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261, mmax=10, ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n), 
c     1        xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7), 
c     1        betaC(2,8) !! 
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax),Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       real*8 z1,ahz,ahpz
c       real*8 xx(8)
c       real*8 A(8,8), detm1, detm2, fact(8,8), Ainv(8,8), Ader(8,8)
c       real*8 interior(8,8)
c       integer iseed,iaccept(n)
c       integer m(n),ng(2)
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy1C/idumC1
c       common  /dummy2C/idumC2
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c       external  DRNUNF
c        k=idumC1
c        betaC(k,idumC2)=z1
cc......calculate log likelihood and its derivative
c       sum1=0.0d0
c       sum2=0.0d0
c       do k1a = 1,8
c        do k2a=1,8
c         A(k1a,k2a)=0.0d0
c         Ader(k1a,k2a)=0.0d0
c        enddo
c       enddo
c       do 300 i=1,n
c        if (m(i) .eq. 1) goto 300
c        xx(1)=1.0d0
c        do j=2,8
c         xx(j)=xcure(j-1,i)
c        enddo
c        ein=0.0d0
c        do j=1,8
c         ein=ein+betaC(k,j)*xx(j)
c        enddo
c        sum1=sum1+ein*cureD(i,k) 
c        j=idumC2
c        sum2=sum2+xx(j)*cureD(i,k) 
c        if (ein .le. 0.0d0) then
c         sum1=sum1-dlog(1.0d0+dexp(ein))
c         sum2=sum2-xx(j)*dexp(ein)/(1.0d0+dexp(ein))
c         wi=dexp(ein)/((1.0d0+dexp(ein))**2) !part of prior
c        else
c         sum1=sum1-ein-dlog(1.0d0+dexp(-ein))
c         sum2=sum2-xx(j)/(1.0d0+dexp(-ein))
c         wi=dexp(-ein)/((1.0d0+dexp(-ein))**2) !part of prior
c        endif
c        dwi=xx(idumC2)*dexp(ein)*(1.0d0-dexp(ein)) !issue - depends on beta
c     1   /((1.0d0+dexp(ein))**3)!part of prior derivative(no floating control)
c        do k1a=1,8
c         do k2a=1,8
c          A(k1a,k2a)=A(k1a,k2a)+xx(k1a)*wi*xx(k2a)
c          Ader(k1a,k2a)=Ader(k1a,k2a)+xx(k1a)*dwi*xx(k2a)
c         enddo
c        enddo
c 300    continue
cc       if (idumC2 .eq. 1) then
cc        sum1=sum1-z1**2/(2.0d0*1.0d0) !kernel of N(0,1) on log scale
cc        sum2=sum2-z1/1.0d0 !derivative on log scale
cc       else
cc        sum1=sum1-z1**2/(2.0d0*3.0d0) !kernel of N(0,3) on log scale
cc        sum2=sum2-z1/3.0d0 !derivative on log scale
cc       endif
c        call DLFTDS (8,A,8,fact,8) !nrow, matrix, ncol, output of factored matrix, nrow of fact (factors the matrix)
c        call DLFDDS(8,fact,8,detm1,detm2) !outputs determinant
c        sum1=sum1+dlog(detm1*10.0d0**detm2) !log of Jeffreys-type prior (this time it's plus!)
c        !need pieces of Jacobi formula for derivative of Jeffreys-type prior
c        call DLINDS(8,A,8,Ainv,8) !outputs inverse of "A" matrix as "Ainv"
c        !multiplying Ainv by Ader (call it "interior")
c        do jder=1,8
c         do ider=1,8
c          tmp=0.0d0
c          do kder=1,8
c              tmp=tmp+Ainv(ider,kder)*Ader(kder,jder)
c          enddo
c          interior(ider,jder)=tmp
c         enddo
c        enddo
c        !take trace of Ainv*Ader
c        inttrace=0.0d0
c        do itr=1,8
c         inttrace=inttrace+interior(itr,itr)
c        enddo
cc        jaconum=(detm1*10.0d0**detm2)* !multiply |A|*tr(Ainv*Ader)
c        sum2=sum2+inttrace !derivative on log scale of Jeffreys-type prior (this time it's plus!)
c       ahz=sum1
c       ahpz=sum2   
c       return
c       end
           
       


       subroutine Galam( iseed )
c      generate censoring indicator for missing causes of failure 
       implicit real*8 (a-h,o-z)
       parameter(n=66357, ngmax=110) !capK=5
       real*8 ts(n),dold(n),d(n),x(8,n),xs(8,n), 
     1        xcure(8,n)
       real*8 beta(5,8), betaC(5,9) !!
       real*8 cureD(n,5) !!
       real*8 alam(5,ngmax)
       real*8 bo(5,2)
       real*8 s(5,ngmax)
       real*8 shape(5,ngmax),scale(5,ngmax),rg(1)
       integer iseed,iaccept(n)
       integer ng(5)
c      data
       common  /vects/ts
       common  /vecd/dold ! taking a value of 0,...,6 with 6 denoting unknown cause
       common  /vecx/x
       common  /vecxs/xs
       common  /vecxcure/xcure
       common  /vecs/s
       common  /vecng/ng
c      quantities to be updated
       common  /vecd/d    ! indicator for censoring (0) or cause of death (1,...,5)  
       common  /vbeta/beta ! coefficients of fixed covariates in surv model
       common  /valam/alam ! constant piecewise hazard parameter
       common  /vecbo/bo    ! hyper-parameter for lam
       common  /dummy3/idum3
       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
       external  DRNGAM
       do k=1,5
        do j=1,ng(k)
         shape(k,j)=bo(k,1)
         scale(k,j)=bo(k,2)
        enddo
       enddo 
       
       do k=1,5
        ak=real(k)
        do 100 i=1,n
         if (cureD(i,k) .eq. 1.0d0) goto 100
         ein=0.0d0
         do j=1,8
          ein=ein+xs(j,i)*beta(k,j)
         enddo
         if (ts(i) .le. s(k,1)) jg=1
         do j1=1,ng(k)-1
          if ( (ts(i) .gt. s(k,j1)) .and.
     1       (ts(i) .le. s(k,j1+1)) ) then
             jg=j1+1
          endif
         enddo
         if (jg .eq. 1) then
          scale(k,1)=scale(k,1)+ts(i)*dexp(ein)
          if (d(i) .eq. ak) then
            shape(k,1)=shape(k,1)+1.0d0
          endif 
         else
          scale(k,1)=scale(k,1)+s(k,1)*dexp(ein)
         if (d(i) .eq. ak) then
           shape(k,jg)=shape(k,jg)+1.0d0
         endif 
          do j=2,jg-1
           scale(k,j)=scale(k,j)+(s(k,j)-s(k,j-1))*dexp(ein)
          enddo
          scale(k,jg)=scale(k,jg)+(ts(i)-s(k,jg-1))*dexp(ein)
         endif
 100    continue
c
        do j=1,ng(k)
         shape1=shape(k,j)
         call rnset( iseed )
         call DRNGAM(1,shape1,rg)
         call rnget( iseed )
         alam(k,j)=rg(1)/scale(k,j)
        enddo
       enddo  
   
       return
       end 
 
c       subroutine GthetaR( iseed )
cc      generate censoring indicator for missing causes of failure 
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n), 
c     1        xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7), 
c     1        betaC(2,8) !! 
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 ynewlo,ysec,reqmin
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       real*8 start(2),xmin(2),step1(2)
c       real*8 xiold(2),xinew(2),ximean(2)
c       real*8 sig2inv(2,2),sig2inv2(2,2),fthetaRdev2(2,2)
c       real*8 tol,Rsig(2,2),RV(1,2)
c       integer iseed,iaccept(n),idum1
c       integer m(n),ng(2)
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy1/idum1
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c       external fthetaR,DRNUNF
       
c      write(*,*) 'before'
       
c        do 81 i=1,n
c         if (m(i) .eq. 1) goto 81
c         idum1=i
c         nopt = 2
c         reqmin = 1.0d-10
c         konvge = 5
c         kcount = 1000
c         step1(1)=0.2d0
c         step1(2)=0.2d0
c         start(1)=thetaR(i,1)
c         start(2)=thetaR(i,2)
c         xiold(1)=thetaR(i,1)
c         xiold(2)=thetaR(i,2)
c        
c         call nelmin(fthetaR, nopt, start, xmin, ynewlo, reqmin, step1,
c     1         konvge, kcount, icount, numres, ifault1) ! need to know the function 
c         ximean(1)=xmin(1)
c         ximean(2)=xmin(2)
cc........calculate Hessian matrix
cc        write(*,*) 'ximean',ximean(1),ximean(2)
c         call dev2(ximean(1),ximean(2),fthetaRdev2)
cc        write(*,*) 'dev2',fthetaRdev2
c         sig2inv=-fthetaRdev2
cc         write(*,*) 'i=',i,'sig2inv=',sig2inv 
c         call DLINDS(2,sig2inv,2,sig2inv2,2)
cc........generate multivariate normal
c         lda=2
c         ldr=2
c         tol=100.0d0*DMACH(4)
c         call DCHFAC(2,sig2inv2,lda,tol,irank,Rsig,ldr)
c         
c         ein=0.0d0
c         do j1=1,2
c          do j2=1,2
c           ein=ein+(xiold(j2)-ximean(j2))*sig2inv(j2,j1)
c     1            *(xiold(j1)-ximean(j1))/2.0d0
c          enddo
c         enddo
c         
c         sumold=-fthetaR(xiold)+ein
c!		 write(*,*) 'before'
cc.........perform metropolis
c         do imetr=1,10
c          nr=1
c          call rnset( iseed ) !call RNSET to initialize the seed
c          call DRNMVN(nr,2,Rsig,2,RV,1)
c          call rnget( iseed )
c          do j=1,2
c           xinew(j)=ximean(j)+RV(1,j)
c          enddo
c
c          ein=0.0d0
c          do j1=1,2
c           do j2=1,2
c            ein=ein+(xinew(j2)-ximean(j2))*sig2inv(j2,j1)
c     1             *(xinew(j1)-ximean(j1))/2.0d0
c           enddo
c          enddo
c
c          sumnew=-fthetaR(xinew)+ein
c          ratio=sumnew-sumold
cc         write(*,*) 'ratio',ratio
c          if (ratio .ge. 0.0d0) then
c           xiold=xinew
c           thetaR(i,1)=xinew(1)
c           thetaR(i,2)=xinew(2)
c           sumold=sumnew
c           iaccept(i)=iaccept(i)+1
c          else
c           call rnset( iseed )
c            u=DRNUNF() !generate random number from U(0,1)
c           call rnget( iseed )
c           if (dlog(u) .le. ratio) then
c            xiold=xinew
c            sumold=sumnew
c            iaccept(i)=iaccept(i)+1
c            thetaR(i,1)=xinew(1)
c            thetaR(i,2)=xinew(2)
c           endif
c          endif
c         enddo
cc       write(*,*) 'Running', i
c 81    continue
c   
c       return
c       end 
   
c       real*8 function fthetaR(z1)
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n), 
c     1        xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7), 
c     1        betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       real*8 start(2),xmin(2),step1(2)
c       real*8 xiold(2),xinew(2),ximean(2)
c       real*8 sig2inv(2,2),sig2inv2(2,2),fthetaRdev2(2,2)
c       real*8 tol,Rsig(2,2),RV(1,2)   
c       real*8 z1(2)
c       integer iseed,iaccept(n),idum1
c       integer m(n),ng(2)
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy1/idum1
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c       external DRNUNF

c       thetaR(idum1,1)=z1(1)
c       thetaR(idum1,2)=z1(2)
c       i=idum1
c       sum1=0.0d0
c   
c       do k=1,2
c        if (cureD(i,k) .eq. 0.0d0) then
c         ein=alpha(k,1)*thetaR(i,1)
c     1       +alpha(k,2)*thetaR(i,2)
c         ak=real(k)
c        do j=1,7
c         ein=ein+xs(j,i)*beta(k,j)
c        enddo
c        if (ts(i) .le. s(k,1)) jg=1
c        do j1=1,(ng(k)-1)
c         if ( (ts(i) .gt. s(k,j1)) .and.
c     1        (ts(i) .le. s(k,j1+1)) ) then
c           jg=j1+1
c         endif
c        enddo
c        if (jg .eq. 1) then
c          ein1=alam(k,1)*ts(i)
c        else
c         ein1=alam(k,1)*s(k,1)
c         do j2=2,jg-1
c          ein1=ein1+alam(k,j2)*(s(k,j2)-s(k,j2-1))
c         enddo
c         ein1=ein1+alam(k,jg)*(ts(i)-s(k,jg-1))
c        endif
c        ss=-ein1*dexp(ein)
c        sum1=sum1+ss
c        if (d(i) .eq. ak) then
c         sum1=sum1+ein
c        endif
c       endif
c      enddo 
cc......longitudinal part
c       ss2=0.0d0
c       ein1=0.0d0
c        do j1=1,7
c         ein1=ein1+x(j1,i)*gam(j1)
c        enddo
c        do j=1,m(i)
c         ein3=Omega(1,1)+Omega(2,1)*ty(i,j) 
c         ein2=y(i,j)-(theta(1)+ein3*thetaR(i,1)
c     1    +(theta(2)+Omega(2,2)*thetaR(i,2))*ty(i,j)+ein1)
c         ss2=ss2-ein2**2/(2.0d0*sig2)
c        enddo
c        sum1=sum1+ss2
c        ss3=0.0d0
c        do j1=1,2
c         ss3=ss3-thetaR(i,j1)**2/2.0d0
cc        do j2=1,2
cc         ss3=ss3-thetaR(i,j2)*Omegainv(j2,j1)
cc    1       *thetaR(i,j1)/2.0d0
cc        enddo
c        enddo
c        sum1=sum1+ss3
c        fthetaR=-sum1
c
c       return
c       end
   
   
c       subroutine dev2(z1,z2,fthetaRdev2)
c       implicit real*8 (a-h,o-z)
c       parameter(n=35261,mmax=10,ngmax=110)
c       real*8 y(n,mmax),ty(n,mmax),ts(n),dold(n),d(n),x(7,n),xs(7,n), 
c     1         xcure(7,n)
c       real*8 gam(7),theta(2),sig2,Omega(2,2),alpha(2,2),beta(2,7), 
c     1         betaC(2,8) !!
c       real*8 thetaR(n,2), cureD(n,2) !!
c       real*8 alam(2,ngmax), Omegainv(2,2)
c       real*8 bo(2,2)
c       real*8 s(2,ngmax)
c       real*8 start(2),xmin(2),step1(2)
c       real*8 xiold(2),xinew(2),ximean(2)
c       real*8 sig2inv(2,2),sig2inv2(2,2),fthetaRdev2(2,2)
c       real*8 tol,Rsig(2,2),RV(1,2)
c       real*8 sum2(2,2),z1,z2  
c       integer iseed,iaccept(n),idum1
c       integer m(n),ng(2)
cc      data
c       common  /vecy/y
c       common  /vecty/ty
cc      common  /vectymean/tymean
c       common  /vects/ts
c       common  /vecd/dold ! taking a value of 0,1,2,3 with 3 denoting unknown cause
c       common  /vecx/x
c       common  /vecxs/xs
c       common  /vecxcure/xcure
c       common  /vecm/m
c       common  /vecs/s
c       common  /vecng/ng
cc      quantities to be updated
c       common  /vecd/d    ! indicator for censoring (0) or cause of death (1 or 2)  
c       common  /vgam/gam   ! coefficients of fixed covariates in long model
c       common  /vtheta/theta ! the mean of random effects
c       common  /vsig2/sig2    ! variance of long model
c       common  /vthetaR/thetaR ! random effects
c       common  /vOmega/Omega  ! variance of random effects
c       common  /vOmegainv/Omegainv  ! inverse of the variance of random effects
c       common  /valpha/alpha ! coefficients for random effects in surv model
c       common  /vbeta/beta ! coefficients of fixed covariates in surv model
c       common  /valam/alam ! constant piecewise hazard parameter
c       common  /vecbo/bo    ! hyper-parameter for lam
c       common  /dummy1/idum1
c       common  /vbetaC/betaC !! coefficients of fixed covariates in cure rate model
c       common  /vcureD/cureD !! coefficients of fixed covariates in surv model
c       external DRNUNF

c       thetaR(idum1,1)=z1
c       thetaR(idum1,2)=z2
c       i=idum1     
cc......calculate the second derivative
c        do j1=1,2
c         do j2=1,2
c           sum2(j1,j2)=0.0d0
c         enddo
c        enddo
c
cc......survival part
c        do k=1,2
c        if (cureD(i,k) .eq. 0.0d0) then
c          ein=alpha(k,1)*thetaR(i,1)
c     1      +alpha(k,2)*thetaR(i,2)
c         do j1=1,7
c          ein=ein+beta(k,j1)*xs(j1,i)
c         enddo
c         if (ts(i) .le. s(k,1)) jg=1
c         do j1=1,ng(k)-1
c          if ( (ts(i) .gt. s(k,j1)) .and.
c     1         (ts(i) .le. s(k,j1+1)) ) then
c           jg=j1+1
c          endif
c         enddo
c         if (jg .eq. 1) then
c          ein1=alam(k,1)*ts(i)
c         else
c          ein1=alam(k,1)*s(k,1)
c          do j2=2,jg-1
c           ein1=ein1+alam(k,j2)*(s(k,j2)-s(k,j2-1))
c          enddo
c          ein1=ein1+alam(k,jg)*(ts(i)-s(k,jg-1))
c         endif
c         ss1=-ein1*dexp(ein)
c         sum2(1,1)=sum2(1,1)+alpha(k,1)**2*ss1
c         sum2(1,2)=sum2(1,2)+alpha(k,1)*alpha(k,2)*ss1
c         sum2(2,1)=sum2(2,1)+alpha(k,2)*alpha(k,1)*ss1
c         sum2(2,2)=sum2(2,2)+alpha(k,2)**2*ss1
c        endif
c        enddo 
cc......longitudinal part
c        do j=1,m(i)
c         ein3=Omega(1,1)+Omega(2,1)*ty(i,j) 
c         ein4=Omega(2,2)*ty(i,j)
c         sum2(1,1)=sum2(1,1)-ein3**2/sig2
c         sum2(1,2)=sum2(1,2)-ein3*ein4/sig2
c         sum2(2,1)=sum2(2,1)-ein3*ein4/sig2
c         sum2(2,2)=sum2(2,2)-ein4**2/sig2
cc        sum2(1,1)=sum2(1,1)-1.0d0/sig2
cc        sum2(1,2)=sum2(1,2)-ty(i,j)/sig2
cc        sum2(2,1)=sum2(2,1)-ty(i,j)/sig2
cc        sum2(2,2)=sum2(2,2)-ty(i,j)**2/sig2
c        enddo
c        sum2(1,1)=sum2(1,1)-1.0d0
c        sum2(2,2)=sum2(2,2)-1.0d0
cc       sum2(1,1)=sum2(1,1)-Omegainv(1,1)
cc       sum2(1,2)=sum2(1,2)-(Omegainv(1,2)+Omegainv(2,1))/2.0d0
cc       sum2(2,1)=sum2(2,1)-(Omegainv(1,2)+Omegainv(2,1))/2.0d0
cc       sum2(2,2)=sum2(2,2)-Omegainv(2,2)
c        fthetaRdev2=sum2 
c
c       return
c       end   