       program Bayessurvonly
c......Bayesian Modeling of Competing Risk Survival Data 
c......in the Presence of Cure Fractions (No Cure)
c......Austin Menger, Md Tuhin Sheikh, and Ming-Hui Chen, UConn
       implicit real*8 (a-h,o-z) 
       parameter(n=66357, ngmax=110,np=8) !capK=5
       parameter(iprint=3,maxlag=20)
       real*8 ts(n),dold(n),d(n),x(8,n),xs(8,n)!,xcure(8,n),xx(9) 
       real*8 xold(8,n),xbar(8),xbar2(8),xsd(8)
   
       real*8 beta(5,8)!,betaC(5,9)!beta=betaNC
       real*8 alam(5,ngmax)!,cureD(n,5)
    
       real*8 ealam(5,ngmax),sealam(5,ngmax)
       real*8 ebeta(5,8),sebeta(5,8)
!       real*8 ebetaC(5,9)
c       real*8 ecr(n,5)
       real*8 ecindex(2),secindex(2)
c       
       real*8 salam(5,ngmax),s2alam(5,ngmax)
       real*8 sbeta(5,8),s2beta(5,8)
!       real*8 sbetaC(5,9),s2betaC(5,9),sebetaC(5,9)
!       real*8 scureD(n,5)
       real*8 scindex(5),s2cindex(5)
       real*8 swcindex(5),s2wcindex(5)

       real*8 bo(5,2)
   
       real*8 alamupp(5,ngmax),alamlow(5,ngmax)
       real*8 cindexupp(5),cindexlow(5)
       real*8 betaupp(5,8),betalow(5,8)
!       real*8 betaClow(5,9),betaCupp(5,9)
       real*8 alow(2),aupp(2)  !this is 2 for two types of intervals
c
!       real*8 ecureD(n,5)
       real*8 seqlam(5,ngmax,50001)
       real*8 seqbeta(5,8,50001)!,seqbetaC(5,9,50001)
       real*8 seqcindex(5,50001)

       real*8 ac(maxlag+1),acv(maxlag+1),seac(maxlag+1),
     1         xacmean
       real*8 ahpd(50001)
c
       integer ide(n),ide2(5,n),ne(5) ! ids with observed/unobserved events
       real*8 comp(5),tie(5),conc(5),cindex(5)       

       integer now(3)
       integer nowb(3),nowe(3)
       integer ntoday(3),ntodayb(3)
       real etime
       real elapsed(2)
c
       real*8 s(5,ngmax),tt(5,n),p(ngmax),sall(ngmax)
       real*8 tt1(5,n)
       real*8  s6(ngmax),tt3(n)
c       real*8  sumC(5) !added for cure DEV
       real*8  sumNC(5) ! added for noncure DEV
       real*8  ss(5), h(5), cpr(5), einL(5)!,fbetaC(5)
       real*8  prodint(5),cprnum2(5)
       integer imean,iseopt
       integer icount1(5),ntt(5)
       integer iseed,iaccept(n)
       
       integer idc(n),ng(5)
c      data
       common  /vects/ts
       common  /vecd/dold ! taking a value of 0,...,6 with 6 denoting unknown cause
       common  /vecx/x
       common  /vecxs/xs
!       common  /vecxcure/xcure ! cure rate model covariates
c       common  /vecnstar/nstar
       common  /vecs/s
       common  /vecng/ng
c      quantities to be updated
       common  /vecd/d    ! indicator for censoring (0) or cause of death (1,...,5)  
       common  /vbeta/beta ! coefficients of fixed covariates in surv model
       common  /valam/alam ! constant piecewise hazard parameter
       common  /vecbo/bo    ! hyper-parameter for lam
       common  /dummy3/idum3
       common  /veccindex/cindex
!       common  /vbetaC/betaC ! coefficients of fixed covariates in cure rate model
!       common  /vcureD/cureD ! cure indicators
       external dacf,DRNUNF,DRNNOF,DLINDS 
       call idate(ntoday)
       ntodayb(1)=ntoday(1)
       ntodayb(2)=ntoday(2)
       ntodayb(3)=ntoday(3)
       call itime(now)
       nowb(1)=now(1)
       nowb(2)=now(2)
       nowb(3)=now(3)
!      number of warm-up iterations
       nbin=2500 !2500 to fully run, was 100 for tuning and for first paper submission
c      nbin=25
!      number of thinning
       nthin=5 !5 to fully run, was 1 for tuning and for first paper submission
c......set pieces for baseline hazard function, ng 
       ng(1)=80 !was 50
       ng(2)=49 !was 50
       ng(3)=12 !was 50
       ng(4)=11 !was 50
       ng(5)=100 !was 50
       ng6=15
c......prior parameter for lambda
       bo(1,1)=0.001d0  ! shape parameter in gamma prior
       bo(1,2)=0.001d0  ! scale parameter in gamma prior
       bo(2,1)=0.001d0  ! shape parameter in gamma prior
       bo(2,2)=0.001d0  ! scale parameter in gamma prior
       bo(3,1)=0.001d0  ! shape parameter in gamma prior
       bo(3,2)=0.001d0  ! scale parameter in gamma prior
       bo(4,1)=0.001d0  ! shape parameter in gamma prior
       bo(4,2)=0.001d0  ! scale parameter in gamma prior
       bo(5,1)=0.001d0  ! shape parameter in gamma prior
       bo(5,2)=0.001d0  ! scale parameter in gamma prior
c......read surv data
c    "id" "t" "delta" "age" "race_black" "race_other" "diag_year"
c     "hisp_ind" "stage" "grade" "stage_grade_INT"
       open(unit=16,file='prostate_trim_sorted.txt',status='old')
c       open(unit=16,file='prostate_trim.txt',status='old')
       do i=1,n
        read(16,*) idc(i),ts(i),dold(i),(xold(j,i),j=1,8)
        if (i .le. 10) then
           write(*,65) idc(i),ts(i),nint(dold(i)),
     1         (xold(j,i),j=1,8)
 65        format(I6,f16.6,I6,8f16.6)
        endif
       enddo
       close(16)

c......convert month to year
       do isurv=1,n
        ts(isurv)=ts(isurv)/12.0d0
       enddo
       
!      standardizing the covariates
       do j=1,8
        xbar(j)=0.0d0
        xbar2(j)=0.0d0 
       enddo
       n1=0.0d0
       do 175 i=1,n
        n1=n1+1 !no. of ids 
        do j=1,8
         xbar(j)=xbar(j)+xold(j,i)
         xbar2(j)=xbar2(j)+xold(j,i)**2
        enddo
175    continue 

       do j=1,8
        xbar(j)=xbar(j)/real(n1)
        xsd(j)=dsqrt((xbar2(j)-real(n1)*xbar(j)**2)/real(n1-1))
       enddo

       do 176 i=1,n
        do j=1,8
         xs(j,i)=(xold(j,i)-xbar(j))/xsd(j)
         x(j,i)=(xold(j,i)-xbar(j))/xsd(j)
!         xcure(j,i)=x(j,i)
        enddo
176    continue 

c......set initial values of s
        do k=1,5
         do j2=1,ngmax
          s(k,j2)=0.0d0
         enddo        
        enddo

      do k=1,5
       icount1(k) = 0
      enddo
      icount6=0

      do 335 i =1,n
        if ((dold(i) .eq. 1.0d0) .or. (dold(i) .eq. 6.0d0)) then 
          icount1(1) = icount1(1) + 1
          tt1(1,icount1(1)) = ts(i)
        endif
        if ((dold(i) .eq. 2.0d0) .or. (dold(i) .eq. 6.0d0)) then
          icount1(2) = icount1(2) + 1 
          tt1(2,icount1(2)) = ts(i)
        endif
        if ((dold(i) .eq. 3.0d0) .or. (dold(i) .eq. 6.0d0)) then
          icount1(3) = icount1(3) + 1 
          tt1(3,icount1(3)) = ts(i)
        endif
        if ((dold(i) .eq. 4.0d0) .or. (dold(i) .eq. 6.0d0)) then
          icount1(4) = icount1(4) + 1 
          tt1(4,icount1(4)) = ts(i)
        endif
        if ((dold(i) .eq. 5.0d0) .or. (dold(i) .eq. 6.0d0)) then
          icount1(5) = icount1(5) + 1 
          tt1(5,icount1(5)) = ts(i)
        endif
        if (dold(i) .eq. 6.0d0) then !masked cause
          icount6 = icount6 + 1 
          tt3(icount6) = ts(i)
        endif
335   continue

      do k=1,5
       ntt(k) = icount1(k) 
      enddo
      ntt6=icount6


c.....Sorting the survival time
      do k=1,5
       do i1=1,ntt(k)-1
        do i2=i1+1,ntt(k)
         if (tt1(k,i1) .gt. tt1(k,i2)) then
          einn=tt1(k,i1)
          tt1(k,i1)=tt1(k,i2)
          tt1(k,i2)=einn
         endif
        enddo
       enddo
      enddo

c.....Getting distinct failure times
      do k=1,5
       tt(k,1)=tt1(k,1)
      enddo

      do k=1,5
       icount=1.0d0  
       do i1=2,ntt(k)
        if (tt1(k,i1) .gt. tt(k,icount)) then
          tt(k,icount+1)=tt1(k,i1)
          icount=icount+1
        endif
       enddo
       ntt(k)=icount
      enddo 

      do i1=1,ntt6-1
       do i2=i1+1,ntt6
        if (tt3(i1) .gt. tt3(i2)) then
         einn=tt3(i1)
         tt3(i1)=tt3(i2)
         tt3(i2)=einn
        endif
       enddo
      enddo

      do k1=1,5
       ang=real(ng(k1))
       k=int(dlog(ang)/dlog(2.0d0))
       mm=ng(k1)-2**k
       kk=2**k-1
       do i=1,kk
        p(i)=real(i)/real(kk+1)
       enddo
       do j=1,mm
        p(kk+j)=real(2*j-1)/real(2**(k+1))
       enddo
       do i=1,kk+mm-1
        do j=i+1,kk+mm
         if (p(i) .gt. p(j)) then
          temp=p(i)
          p(i)=p(j)
          p(j)=temp
         endif
        enddo
       enddo
       do i=1,kk+mm
        write(*,*) 'j=',i,' p=',p(i)
       enddo
c.....compute cut points for td
       do j=1,ng(k1)-1
        antt=real(ntt(k1))
        k=int(antt*p(j))
        einn=antt*p(j)!ein=antt*p(j)
        pp=einn-real(k)!ein-real(k)
        if (pp .eq. 0.0d0) then
         s(k1,j)=(tt(k1,k)+tt(k1,k+1))/2.0d0
        else
         s(k1,j)=tt(k1,k+1)
        endif
        write(*,*) 'k1=',k1,' ntt=',ntt,' j=',j,' s=',s(k1,j)
       enddo
       s(k1,ng(k1))=10000.0d0
      enddo

       ang=real(ng6)
       k=int(dlog(ang)/dlog(2.0d0))
       mm=ng6-2**k
       kk=2**k-1
       do i=1,kk
         p(i)=real(i)/real(kk+1)
       enddo
       do j=1,mm
        p(kk+j)=real(2*j-1)/real(2**(k+1))
       enddo
       do i=1,kk+mm-1
        do j=i+1,kk+mm
         if (p(i) .gt. p(j)) then
          temp=p(i)
          p(i)=p(j)
          p(j)=temp
         endif
        enddo
       enddo
       do i=1,kk+mm
        write(*,*) 'j=',i,' p=',p(i)
       enddo
c.....compute cut points for td
       do j=1,ng6-1
        antt=real(ntt6)
        k=int(antt*p(j))
        einn=antt*p(j)!ein=antt*p(j)
        pp=einn-real(k)!ein-real(k)
        if (pp .eq. 0.0d0) then
         s6(j)=(tt3(k)+tt3(k+1))/2.0d0
        else
         s6(j)=tt3(k+1)
        endif
        write(*,*) 'ntt6=',ntt6,' j=',j,' s=',s6(j)
       enddo
       s6(ng6)=10000.0d0

c......ide: contains the ids for whom event occurred
c...........with at least 2 long data
       icount=0
       do 333 i=1,n
c        following condition is needed for long data              
c         if(m(i) .eq. 1) goto 333
         if(dold(i) .ne. 0.0d0) then
          icount=icount+1         
          ide(icount)=i 
         endif  
333    continue      
       nide=icount 

       
c......initialize random seed
       iseed=9999999
 
       do i=1,n
        iaccept(i)=0
       enddo

c......initialize surv params

c       beta(2,1)=-0.26652d0
c       beta(2,2)=0.12086d0
c       beta(2,3)=-0.80538d0
c       beta(2,4)=0.32776d0
c       beta(2,5)=0.53856d0
c       beta(2,6)=-0.37988d0
c       beta(2,7)=-0.17423d0
c       beta(1,1)=0.02422d0
c       beta(1,2)=-0.07167d0
c       beta(1,3)=-0.72296d0
c       beta(1,4)=0.10791d0
c       beta(1,5)=0.71234d0
c       beta(1,6)=-0.34808d0
c       beta(1,7)=-0.14978d0
       
       !initializing beta with coxph in R
       do k=1,5
        beta(k,1)=0.78268d0 !age
        beta(k,2)=0.09915d0 !race_black
        beta(k,3)=-0.02784d0 !race_other
        beta(k,4)=-0.13430d0 !diag_year
        beta(k,5)=0.01350d0 !hisp_ind
        beta(k,6)=0.10655d0 !stage
        beta(k,7)=0.12668d0 !grade
        beta(k,8)=0.13408d0 !stage_grade_int
       enddo

       do k=1,5
        do j=1,ngmax
         alam(k,j)=0.001d0
        enddo
       enddo

!       do k=1,5
!        do j=1,9
!         betaC(k,j)=0.0d0
!        enddo 
!       enddo

c.......initialize d(i) (I adjusted the logic here a bit)
       do 224 i=1,n
        if (dold(i) .lt. 6.0d0) d(i)=dold(i)
         if (dold(i) .eq. 6.0d0) then 
          call rnset( iseed )
           u=DRNUNF() !generate random number from U(0,1)
          call rnget( iseed )
          if (u .le. 0.2d0) d(i)=1.0d0
          if ((u .gt. 0.2d0) .and. (u .le. 0.4d0)) d(i)=2.0d0
          if ((u .gt. 0.4d0) .and. (u .le. 0.6d0)) d(i)=3.0d0
          if ((u .gt. 0.6d0) .and. (u .le. 0.8d0)) d(i)=4.0d0
          if ((u .gt. 0.8d0) .and. (u .le. 1.0d0)) d(i)=5.0d0
         endif
 224   continue

c......initialize cureD(i) for when we don't know an event occurred
!       do 225 i=1,n
!        if (d(i) .eq. 0.0d0) then 
!         do k=1,5 
!          call rnset( iseed )
!           u1=DRNUNF() !generate random number from U(0,1)
!          call rnget( iseed )
!          if (u1 .le. 0.5d0) then 
!           cureD(i,k)=0.0d0 
!          else
!           cureD(i,k)=1.0d0
!          endif 
!         enddo
!        else 
!         k1 = int(d(i))  !known, actually die from this cause
!         cureD(i,k1) = 0.0d0
!         ! have to adjust the logic here a bit for d(i) ≠ k
!         do k=1,5
!          if (int(d(i)) .ne. k) then
!c           k=int(3.0d0 - d(i)) ! d(i) ≠ k
!           call rnset( iseed ) !call again for k=2
!            u2=DRNUNF() !generate random number from U(0,1)
!           call rnget( iseed )
!           if (u2 .le. 0.5d0) then 
!            cureD(i,k)=0.0d0 
!           else
!            cureD(i,k)=1.0d0
!           endif
!          endif
!         enddo 
!        endif      
! 225   continue 
c......read the number of Gibbs iterates
       nrep=5000 !5000 to run, was 250 for tuning and 500 for initial paper

c......warm up gibbs with longitudinal
       call rnset( iseed )

       icount=1
       do 51 i=1,nbin
        call gibbs( iseed )
        if ( (i .eq. 1) .or. (i .eq. icount) ) then
         do k = 1,5
          write(*,*) 'gibbs warm up i=',i,' k=',k 
          write(*,980) (beta(k,j),j=1,8)
980       format('beta=',8f12.6)
!          write(*,972) (betaC(k,j),j=1,9)
!972       format('betaC=',9f12.6)
          write(*,984) (alam(k,j),j=1,5)
984       format('alam1--alam5=',5f12.6)
         enddo
          icount2=0
          do 200 i2=1,1000
           if (icount2 .eq. 5) goto 200
c           if ((m(i2) .gt. 1) .and. (dold(i2) .eq. 3.0d0)) then
           if (dold(i2) .eq. 6.0d0) then
            icount2=icount2+1
            write(*,*) 'i=',i2,' d(i)=',d(i2)
           endif
200       continue
          icount=icount+1
        endif
51     continue
!......initialize cumulative variables

       do k=1,5
        do j=1,8
         sbeta(k,j)=0.0d0
         s2beta(k,j)=0.0d0
        enddo
       enddo

!       do k=1,5
!        do j=1,9
!         sbetaC(k,j)=0.0d0
!         s2betaC(k,j)=0.0d0
!        enddo
!       enddo

!       do i=1,n
!        do k=1,5
!         scureD(i,k)=0.0d0
!        enddo
!       enddo

       do k=1,5
        do j=1,ng(k)
         salam(k,j)=0.0d0
         s2alam(k,j)=0.0d0
        enddo
       enddo

       !for C-index
       do k=1,5
        scindex(k)=0.0d0
        s2cindex(k)=0.0d0
        !swcindex(k)=0.0d0
        !s2wcindex(k)=0.0d0
       enddo
c       do 123 i=1,n
c        do k=1,5
c         ecr(i,k)=0.0d0
c        enddo
c123    continue  


c......gibbs sampling
       icount=0
       sDEVsurv1=0.0d0
       sDEVsurv2=0.0d0
       sDEVsurv3=0.0d0
       sDEVsurv4=0.0d0
       sDEVsurv5=0.0d0
       sDEVsurvTot=0.0d0
  
       do 111 i=1,nrep
        do ithin=1,nthin
         call gibbs( iseed )
        enddo

       if ( (i .eq. 1) .or. (i .eq. icount) ) then
        do k=1,5
         write(*,*) 'gibbs  i=',i,' k=',k 
         write(*,980) (beta(k,j),j=1,8)
!         write(*,972) (betaC(k,j),j=1,9)
         write(*,984) (alam(k,j),j=1,5)
        enddo
        icount2=0
        do 201 i2=1,1000
         if (icount2 .eq. 5) goto 201
c          if ((m(i2) .gt. 1) .and. (dold(i2) .eq. 3.0d0)) then
          if (dold(i2) .eq. 6.0d0) then 
           icount2=icount2+1
           write(*,*) 'i=',i2,' d(i)=',d(i2)
          endif
201      continue
        icount=icount+10
       endif
       

c......codes for c-index starts----
c......ide2(k,): event ids for kth cause 
c......ne(k): total number of observed cases for kth cause 
       do k=1,5
        ne(k)=0.0d0 
        do j=1,nide 
         if(dold(ide(j)) .eq. k) then 
          ne(k)=ne(k)+1
          ide2(k,ne(k))=ide(j)
         endif
        enddo
       enddo

c......indicator whether to compute cause-specific C-index
       compcind=1.0d0
c......if compcind==1, it will compute
       if(compcind .eq. 1.0d0) then
        do 108 k=1,5
         comp(k)=0.0d0
         conc(k)=0.0d0 
         tie(k)=0.0d0 
         do 109 i1=1,ne(k)
          id1=ide2(k,i1)
          ss1=0.0d0
c         following two lines need to uncomment when long data is used          
c          ss1=alpha(k,1)*thetaR(id1,1)
c     1      +alpha(k,2)*thetaR(id1,2)
           do j=1,np !np=8 defined in the parameter (top)
            ss1=ss1+xs(j,id1)*beta(k,j)
           enddo  
          do 110 i2=ide2(k,i1)+1,n !the data needs to be ordered by event time
c         following condition needs to uncomment when long data is used
c           if(m(i2) .eq. 1) goto 110 
c!         do 110 i2=1,n
           if(ts(i2) .gt. ts(id1)) then 
            comp(k)=comp(k)+1.0d0 
            sum1=0.0d0 
            ss2=0.0d0
c           following line needs to uncomment when long data is used
c            ss2=alpha(k,1)*thetaR(i2,1)+alpha(k,2)*thetaR(i2,2) 
            do j=1,np
             ss2=ss2+xs(j,i2)*beta(k,j)
            enddo
            if(ss2 .lt. ss1) then 
             conc(k)=conc(k)+1.0d0
            else if(ss2 .eq. ss1) then 
             tie(k)=tie(k)+1.0d0
            endif 
           endif
110       continue 
109      continue
         cindex(k)=(conc(k)+(0.5*tie(k)))/comp(k)
         write(*,*) 'k=',k,'c-index(k)',cindex(k)
108     continue
       endif       
      

       do k=1,5
        do j=1,8
         sbeta(k,j)=sbeta(k,j)+beta(k,j)
         s2beta(k,j)=s2beta(k,j)+beta(k,j)**2
         seqbeta(k,j,i)=beta(k,j)
        enddo
       enddo

!       do k=1,5
!        do j=1,9
!         sbetaC(k,j)=sbetaC(k,j)+betaC(k,j)
!         s2betaC(k,j)=s2betaC(k,j)+betaC(k,j)**2
!         seqbetaC(k,j,i)=betaC(k,j)
!        enddo
!       enddo

!       do 337 i1=1,n
!        do k=1,5
!         scureD(i1,k)=scureD(i1,k)+cureD(i1,k)
!        enddo
!337    continue 

       do k=1,5
        do j=1,ng(k)
         salam(k,j)=salam(k,j)+alam(k,j)
         s2alam(k,j)=s2alam(k,j)+alam(k,j)**2
         seqlam(k,j,i)=alam(k,j)
        enddo
       enddo

! for c-index
       do k=1,5
        scindex(k)=scindex(k)+cindex(k)
        s2cindex(k)=s2cindex(k)+cindex(k)**2
        seqcindex(k,i)=cindex(k)
c.......weighted
c              swcindex(k)=swcindex(k)+wcindex(k)
c              s2wcindex(k)=s2wcindex(k)+wcindex(k)**2
       enddo

c......storing ecr
c       do 211 i1=1,n
c        do k=1,2
c         eincure=betaC(k,1)
c         do j=1,7 !may not be right in original code
c          eincure=eincure+xcure(j,i1)*betaC(k,j+1)
c         enddo
c         if(eincure .ge. 0.0d0) then
c          ein=1.0d0/(1.0d0+dexp(-eincure))
c         else
c          ein=dexp(eincure)/(1.0d0+dexp(eincure))
c         endif
c         ecr(i1,k)=ecr(i1,k)+ein
c         if ((i1 .eq. 888) .or. (i1 .eq. n)) then
c          write(*,*) 'i=',i
c          write(*,*) 'i1=',i1,' k=',k,' cr=',ein
c          write(*,*) 'ecr=',ecr(i1,k)
c         endif
c        enddo
c211    continue

!------------------------------
!c.......computing DEVsurv
!------------------------------

        sumNC(1)=0.0d0 !loops the sum from i=1,n 
        sumNC(2)=0.0d0
        sumNC(3)=0.0d0
        sumNC(4)=0.0d0
        sumNC(5)=0.0d0
        do 511 i1=1,n
!         xx(1)=1.0d0
!         do j=2,9
!          xx(j)=xcure(j-1,i1)
!         enddo
         do 673 k=1,5 
!          einC = 0.0d0
!          do j=1,9
!           einC=einC+betaC(k,j)*xx(j) !posterior estimate ebetaC
!          enddo
!          if(einC .ge. 0.0d0) then
!           fbetaC(k)=1.0d0/(1.0d0+dexp(-einC)) !creates the pi function
!          else
!           fbetaC(k)=dexp(einC)/(1.0d0+dexp(einC)) !creates the pi function           
!          endif
          einL(k)=0.0d0
          do j=1,8
           einL(k)=einL(k)+xs(j,i1)*beta(k,j)
          enddo
          if (ts(i1) .le. s(k,1)) jg=1
          do 142 j1=1,(ng(k)-1)
           if ( (ts(i1) .gt. s(k,j1)) .and.
     1      (ts(i1) .le. s(k,j1+1)) ) then
             jg=j1+1
           endif
142       continue
          if (jg .eq. 1) then
           ein1=alam(k,1)*ts(i1)
          else
           ein1=alam(k,1)*s(k,1)
           do j2=2,jg-1
            ein1=ein1+alam(k,j2)*(s(k,j2)-s(k,j2-1))
           enddo
           ein1=ein1+alam(k,jg)*(ts(i1)-s(k,jg-1)) !cumulative haz
          endif
          ss(k)=-ein1*dexp(einL(k)) !log of the survival function
          h(k)=dlog(alam(k,jg))+einL(k) !log of hazard function
          prodint(k)=dexp(ss(k))
673      enddo
         !had to change the logic on calculating cpr here.
         cprnum2(1)=prodint(2)*prodint(3)*prodint(4)*prodint(5)
         cprnum2(2)=prodint(1)*prodint(3)*prodint(4)*prodint(5)
         cprnum2(3)=prodint(1)*prodint(2)*prodint(4)*prodint(5)
         cprnum2(4)=prodint(1)*prodint(2)*prodint(3)*prodint(5)
         cprnum2(5)=prodint(1)*prodint(2)*prodint(3)*prodint(4)   

         do 143 k=1,5
c          kp=3-k !k prime
          cprnum1=dexp(h(k))*dexp(ss(k))
c          cprnum2=fbetaC(kp)+(1.0d0-fbetaC(kp))
c     1         *dexp(ss(kp))
          cprnum=cprnum1*cprnum2(k)
!          cprnum=(1.0d0-fbetaC(k))*dexp(h(k))*dexp(ss(k)) !numerator
!     1     * fbetaC(kp)**cureD(i1,kp)
!     1     * ((1.0d0-fbetaC(kp))*dexp(ss(kp)))
!     1     **(1.0d0-cureD(i1,kp))
          cprdenom=0.0d0
c          cprdenom(2)=0.0d0

c          do kin=1,5
c           kinp=3-kin
c           cprdenom(kin)=cprdenom(kin)+(1.0d0-fbetaC(kin))*dexp(h(kin))
c     1      *dexp(ss(kin))*((1.0d0-fbetaC(kinp))*dexp(ss(kinp))
c     1      +fbetaC(kinp))
c          enddo
          do kin=1,5
           cprdenom1=dexp(h(kin))*dexp(ss(kin))
           cprdenom=cprdenom + cprdenom1*cprnum2(kin)
          enddo
c         cpr(k)=cprnum/(cprdenom(1)+cprdenom(2)) !conditional prob
         cpr(k)=cprnum/cprdenom !conditional prob
c         cpr(k)=cprnum - dlog(cprdenom) !log of conditional prob
          if (dold(i1) .lt. 6.0d0) then ! the case where u_i = 0 (known)
!         if (d(i1) .eq. ak) then 
           if (int(d(i1)) .eq. 0) then
            sumNC(k)= sumNC(k)+ss(k)
           else 
            if (int(d(i1)) .eq. k) then
             sumNC(k)=sumNC(k)+h(k)+ss(k)
            else 
              sumNC(k)=sumNC(k)+ss(k)
            endif 
           endif
          endif
c143      continue    
         if (dold(i1) .eq. 6.0d0) then !the case where u_i = 1 (no log)
c          sumNC(k)=sumNC(k)+(dlog(1.0d0-fbetaC(k))+h(k)+ss(k))*cpr(k)
c     1     + dlog(fbetaC(k)+(1.0d0-fbetaC(k))*dexp(ss(k)))*cpr(k) !kprime goes to other four k's
          sumNC(k)=sumNC(k)+(h(k)+ss(k))*cpr(k)
     1     +4.0d0*ss(k)*cpr(k) !kprime goes to other four k's
         endif
!........WAIC
!         sumw(i1)=sumw(i1)+dexp(sum0)
!         write(*,*)'sumw(i1)',sumw(i1)
143      continue
511     continue
        write(*,*) 'DEVsurv1',-2.0d0*sumNC(1)
        write(*,*) 'DEVsurv2',-2.0d0*sumNC(2)
        write(*,*) 'DEVsurv3',-2.0d0*sumNC(3)
        write(*,*) 'DEVsurv4',-2.0d0*sumNC(4)
        write(*,*) 'DEVsurv5',-2.0d0*sumNC(5)
        write(*,*) 'DEVsurvTot', -2.0d0*(sumNC(1)+sumNC(2)+sumNC(3)
     1       +sumNC(4)+sumNC(5))

        sDEVsurv1=sDEVsurv1-2.0d0*sumNC(1)
        sDEVsurv2=sDEVsurv2-2.0d0*sumNC(2)
        sDEVsurv3=sDEVsurv3-2.0d0*sumNC(3)
        sDEVsurv4=sDEVsurv4-2.0d0*sumNC(4)
        sDEVsurv5=sDEVsurv5-2.0d0*sumNC(5)
        sDEVsurvTot=sDEVsurvTot-2.0d0*(sumNC(1)+sumNC(2)
     1       +sumNC(3)+sumNC(4)+sumNC(5))
111    continue
       call rnget( iseed )

c......obtain posterior estimates
   
       do k=1,5
        do j=1,8
         ebeta(k,j)=sbeta(k,j)/real(nrep)
         sebeta(k,j)=dsqrt( (s2beta(k,j)-real(nrep)*ebeta(k,j)**2)
     1       / real(nrep-1) )
        enddo
       enddo

!       do k=1,5
!        do j=1,9
!         ebetaC(k,j)=sbetaC(k,j)/real(nrep)
!         sebetaC(k,j)=dsqrt( (s2betaC(k,j)-real(nrep)*ebetaC(k,j)**2)
!     1       / real(nrep-1) )
!        enddo
!       enddo

!       do 339 i1=1,n
!         do k=1,5
!          ecureD(i1,k)=scureD(i1,k)/real(nrep)
!        enddo
!339    continue 

       do k=1,5
        do j=1,ng(k)
         ealam(k,j)=salam(k,j)/real(nrep)
         sealam(k,j)=dsqrt( (s2alam(k,j)-real(nrep)*ealam(k,j)**2)
     1        / real(nrep-1) )
        enddo
       enddo

       !for c-index
       do k=1,5
        ecindex(k)=scindex(k)/real(nrep)
        secindex(k)=dsqrt( (s2cindex(k)-real(nrep)*ecindex(k)**2)
     1         / real(nrep-1) )
c.......weighted
c        ewcindex(k)=swcindex(k)/real(nrep)
c        sewcindex(k)=dsqrt( (s2wcindex(k)-real(nrep)*ewcindex(k)**2)
c     1         / real(nrep-1) )
       enddo

c......posterior ecr
c       open(unit=31,file='JMecr.txt',
c     1       access='sequential',status='unknown')   
c       do 79 i1=1,n
c         if (m(i1) .eq. 1) goto 79
c         do k=1,2
c          ecr(i1,k)=ecr(i1,k)/real(nrep)
c         enddo
c        write(31,751) i1,ecr(i1,1),ecr(i1,2)
c751     format(1x,I6,2f18.12)
c79     continue 
c       close(31)

!       do 119 i1=1,n
!         if(m(i1) .eq. 1) goto 119
!         do j=1,2
!          ethetaRi(i1,j)=sthetaRi(i1,j)/real(nrep)
!        enddo
!        write(222,753) idc(i1), (ethetaRi(i1,jmac),jmac=1,2)
!753     format(1x,I5,2f20.8)
!119    enddo
!       close(222)

c......output  MCMC chains   
       open(unit=27,file='JMlam.txt',
     1       access='sequential',status='unknown')
       do i=1,nrep
        write(27,656) i, (seqlam(jmac,ng(jmac),i),jmac=1,5)
656     format(1x,I5,5f20.8)
       enddo
       close(27)
      
       open(unit=291,file='JMbeta1.txt',
     1       access='sequential',status='unknown')
       do i=1,nrep
        write(291,6581) i,(seqbeta(1,jmac,i),jmac=1,8)
6581     format(1x,I5,8f20.8)
       enddo
       close(291)

       open(unit=292,file='JMbeta2.txt',
     1       access='sequential',status='unknown')
       do i=1,nrep
        write(292,6582) i,(seqbeta(2,jmac,i),jmac=1,8)
6582     format(1x,I5,8f20.8)
       enddo
       close(292)

       open(unit=293,file='JMbeta3.txt',
     1       access='sequential',status='unknown')
       do i=1,nrep
        write(293,6583) i,(seqbeta(3,jmac,i),jmac=1,8)
6583     format(1x,I5,8f20.8)
       enddo
       close(293)
       
       open(unit=294,file='JMbeta4.txt',
     1       access='sequential',status='unknown')
       do i=1,nrep
        write(294,6584) i,(seqbeta(4,jmac,i),jmac=1,8)
6584     format(1x,I5,8f20.8)
       enddo
       close(294)
       
       open(unit=295,file='JMbeta5.txt',
     1       access='sequential',status='unknown')
       do i=1,nrep
        write(295,6585) i,(seqbeta(5,jmac,i),jmac=1,8)
6585     format(1x,I5,8f20.8)
       enddo
       close(295)

       ! for c-index
	open(unit=31,file='JMcindex.txt',
     1       access='sequential',status='unknown')
       do i=1,nrep
        write(31,6591) i, (seqcindex(jmac,i),jmac=1,5)
6591	  format(1x,I5,5f20.8)
       enddo
       close(31)

!       open(unit=296,file='JMbetaC1.txt',
!     1       access='sequential',status='unknown')
!       do i=1,nrep
!        write(296,6586) i,(seqbetaC(1,jmac,i),jmac=1,9)
!6586     format(1x,I5,9f20.8)
!       enddo
!       close(296)

!       open(unit=297,file='JMbetaC2.txt',
!     1       access='sequential',status='unknown')
!       do i=1,nrep
!        write(297,6587) i,(seqbetaC(2,jmac,i),jmac=1,9)
!6587     format(1x,I5,9f20.8)
!       enddo
!       close(297)

!       open(unit=298,file='JMbetaC3.txt',
!     1       access='sequential',status='unknown')
!       do i=1,nrep
!        write(298,6588) i,(seqbetaC(3,jmac,i),jmac=1,9)
!6588     format(1x,I5,9f20.8)
!       enddo
!       close(298)

!       open(unit=299,file='JMbetaC4.txt',
!     1       access='sequential',status='unknown')
!       do i=1,nrep
!        write(299,6589) i,(seqbetaC(4,jmac,i),jmac=1,9)
!6589     format(1x,I5,9f20.8)
!       enddo
!       close(299)

!       open(unit=300,file='JMbetaC5.txt',
!     1       access='sequential',status='unknown')
!       do i=1,nrep
!        write(300,6590) i,(seqbetaC(5,jmac,i),jmac=1,9)
!6590     format(1x,I5,9f20.8)
!       enddo
!       close(300)  

c......obtain hpd intervals

c......95% HPD intervals for beta
       imean=1
       iseopt=1
       alphahpd=0.05d0
       do jmac1=1,5
        do jmac2=1,8
         write(*,*) 'hpd for beta j=',jmac1,jmac2
         do i=1,nrep
          ahpd(i)=seqbeta(jmac1,jmac2,i)
         enddo
         write(*,*) 'autocorrelation for beta',jmac1,jmac2
         call dacf(nrep,ahpd,iprint,iseopt,imean,xacmean,maxlag,
     1           acv,ac,seac)
         call hpd(nrep,alphahpd,ahpd,alow,aupp)
         betaupp(jmac1,jmac2)=aupp(1)
         betalow(jmac1,jmac2)=alow(1)
        enddo
       enddo

c......95% HPD intervals for betaC
!       imean=1
!       iseopt=1
!       alphahpd=0.05d0
!       do jmac1=1,5
!        do jmac2=1,9
!         write(*,*) 'hpd for betaC j=',jmac1,jmac2
!         do i=1,nrep
!          ahpd(i)=seqbetaC(jmac1,jmac2,i)
!         enddo
!         write(*,*) 'autocorrelation for betaC',jmac1,jmac2
!!         call dacf(nrep,ahpd,iprint,iseopt,imean,xacmean,maxlag,
!!     1           acv,ac,seac)
!         call hpd(nrep,alphahpd,ahpd,alow,aupp)
!         betaCupp(jmac1,jmac2)=aupp(1)
!         betaClow(jmac1,jmac2)=alow(1)
!        enddo
!       enddo

c......95% hpd for alam
       imean=1
       iseopt=1
       alphahpd=0.05d0
       do jmac1=1,5
        do jmac2=1,ng(jmac1)
         write(*,*) 'hpd for alam j=',jmac1,jmac2
         do i=1,nrep
          ahpd(i)=seqlam(jmac1,jmac2,i)
         enddo
         write(*,*) 'autocorrelation for alam',jmac1,jmac2
         call dacf(nrep,ahpd,iprint,iseopt,imean,xacmean,maxlag,
     1           acv,ac,seac)
         call hpd(nrep,alphahpd,ahpd,alow,aupp)
         alamupp(jmac1,jmac2)=aupp(1)
         alamlow(jmac1,jmac2)=alow(1)
        enddo
       enddo

       ! for c-index
c......95% HPD intervals for c-index
       imean=1
       iseopt=1
       alphahpd=0.05d0
       do jmac=1,5
        write(*,*) 'hpd for cindex j=',jmac
        do i=1,nrep
         ahpd(i)=seqcindex(jmac,i)
        enddo
        write(*,*) 'autocorrelation for cindex',jmac
c        call dacf(nrep,ahpd,iprint,iseopt,imean,xacmean,maxlag,
c     1            acv,ac,seac)
        call hpd(nrep,alphahpd,ahpd,alow,aupp)
        cindexupp(jmac)=aupp(1)
        cindexlow(jmac)=alow(1)
       enddo

c......compute DICsurv (replace ebeta, ealam, etc...)
         barDEVsurv1=sDEVsurv1/real(nrep)
         barDEVsurv2=sDEVsurv2/real(nrep)
         barDEVsurv3=sDEVsurv3/real(nrep)
         barDEVsurv4=sDEVsurv4/real(nrep)
         barDEVsurv5=sDEVsurv5/real(nrep)
         barDEVsurvTot=sDEVsurvTot/real(nrep)


         sumNC(1)=0.0d0 !loops the sum from i=1,n
         sumNC(2)=0.0d0
         sumNC(3)=0.0d0
         sumNC(4)=0.0d0
         sumNC(5)=0.0d0
         do 522 i1=1,n
c          einL(1)=0.0d0
c          einL(2)=0.0d0
c          einL(3)=0.0d0
c          einL(4)=0.0d0
c          einL(5)=0.0d0
!          xx(1)=1.0d0
!           do j=2,9
!            xx(j)=xcure(j-1,i1)
!           enddo
           do 253 k = 1,5
!            einC = 0.0d0
!            do j=1,9
!             einC=einC+ebetaC(k,j)*xx(j) !posterior estimate ebetaC
!            enddo
!            if(einC .ge. 0.0d0) then
!             fbetaC(k)=1.0d0/(1.0d0+dexp(-einC)) !creates the pi function
!            else
!             fbetaC(k)=dexp(einC)/(1.0d0+dexp(einC)) !creates the pi function           
!            endif
            einL(k)=0.0d0 
            do j=1,8
             einL(k)=einL(k)+xs(j,i1)*ebeta(k,j) !posterior est ebeta
            enddo
            if (ts(i1) .le. s(k,1)) jg=1
             do 147 j1=1,(ng(k)-1)
              if ( (ts(i1) .gt. s(k,j1)) .and.
     1         (ts(i1) .le. s(k,j1+1)) ) then
               jg=j1+1
              endif
 147         continue
            if (jg .eq. 1) then
             ein1=ealam(k,1)*ts(i1) !posterior est ealam
            else
             ein1=ealam(k,1)*s(k,1) !posterior est ealam
             do j2=2,jg-1
              ein1=ein1+ealam(k,j2)*(s(k,j2)-s(k,j2-1)) !post est ealam
             enddo
             ein1=ein1+ealam(k,jg)*(ts(i1)-s(k,jg-1)) !post est ealam
            endif
            ss(k)=-ein1*dexp(einL(k)) !log of the survival function
            h(k)=dlog(ealam(k,jg))+einL(k) !posterior est ealam
            prodint(k)=dexp(ss(k)) 
 253       enddo
          !had to change the logic on calculating cpr here.
          cprnum2(1)=prodint(2)*prodint(3)*prodint(4)*prodint(5)
          cprnum2(2)=prodint(1)*prodint(3)*prodint(4)*prodint(5)
          cprnum2(3)=prodint(1)*prodint(2)*prodint(4)*prodint(5)
          cprnum2(4)=prodint(1)*prodint(2)*prodint(3)*prodint(5)
          cprnum2(5)=prodint(1)*prodint(2)*prodint(3)*prodint(4) 
          
          do 157 k=1,5
c           kp=3-k !k prime
           cprnum1=dexp(h(k))*dexp(ss(k))
c           cprnum2=(1.0d0-fbetaC(kp))
c     1         *dexp(ss(kp))+fbetaC(kp)
           cprnum=cprnum1*cprnum2(k)
!           cprnum=(1.0d0-fbetaC(k))*dexp(h(k))*dexp(ss(k)) !num
!     1      * fbetaC(kp)**ecureD(i1,kp)
!     1      * ((1.0d0-fbetaC(kp))*dexp(ss(kp)))
!     1      **(1.0d0-ecureD(i1,kp))
           cprdenom=0.0d0
c           cprdenom(2)=0.0d0
c           do kin=1,2
c            kinp=3-kin
c            cprdenom(kin)=cprdenom(kin)+(1.0d0-fbetaC(kin))
c     1       *dexp(h(kin))*dexp(ss(kin))*((1.0d0-fbetaC(kinp))
c     1       *dexp(ss(kinp))+fbetaC(kinp))
c           enddo
           do kin=1,5
            cprdenom1=dexp(h(kin))*dexp(ss(kin))
            cprdenom=cprdenom + cprdenom1*cprnum2(kin)
           enddo
c           cpr(k)=cprnum/(cprdenom(1) + cprdenom(2)) !conditional probability
           cpr(k)=cprnum/cprdenom !conditional prob
c           cpr(k)=cprnum - dlog(cprdenom) !log of conditional prob
           if (dold(i1) .lt. 6.0d0) then ! the case where u_i = 0
            if (int(d(i1)) .eq. 0) then
             sumNC(k)= sumNC(k)+ss(k)
            else 
             if (int(d(i1)) .eq. k) then
              sumNC(k)=sumNC(k)+h(k)+ss(k)
             else 
               sumNC(k)=sumNC(k)+ss(k)
             endif 
            endif
           endif
           if (dold(i1) .eq. 6.0d0) then !the case where u_i = 1 (no log)
            sumNC(k)=sumNC(k)+(h(k)+ss(k))
     1       *cpr(k)+4.0d0*ss(k)*cpr(k) !special b/c separate
c     1      + (1.0d0-fbetaC(k))*cpr(k)*ss(k) !kprime goes to k
           endif
157       continue
522      continue
         
         DEVsurv1bar=-2.0d0*sumNC(1)
         DEVsurv2bar=-2.0d0*sumNC(2)
         DEVsurv3bar=-2.0d0*sumNC(3)
         DEVsurv4bar=-2.0d0*sumNC(4)
         DEVsurv5bar=-2.0d0*sumNC(5)
         DEVsurvTotbar=-2.0d0*(sumNC(1)+sumNC(2)+sumNC(3)+sumNC(4)
     1      +sumNC(5))

         write(*,*) 'DEVsurv1bar',DEVsurv1bar
         write(*,*) 'DEVsurv2bar',DEVsurv2bar
         write(*,*) 'DEVsurv3bar',DEVsurv3bar
         write(*,*) 'DEVsurv4bar',DEVsurv4bar
         write(*,*) 'DEVsurv5bar',DEVsurv5bar
         write(*,*) 'DEVsurvTotbar',DEVsurvTotbar

         pDsurv1=barDEVsurv1-DEVsurv1bar
         pDsurv2=barDEVsurv2-DEVsurv2bar
         pDsurv3=barDEVsurv3-DEVsurv3bar
         pDsurv4=barDEVsurv4-DEVsurv4bar
         pDsurv5=barDEVsurv5-DEVsurv5bar
         pDsurvTot=barDEVsurvTot-DEVsurvTotbar

         write(*,*) 'pDsurv1',pDsurv1
         write(*,*) 'pDsurv2',pDsurv2
         write(*,*) 'pDsurv3',pDsurv3
         write(*,*) 'pDsurv4',pDsurv4
         write(*,*) 'pDsurv5',pDsurv5
         write(*,*) 'pDsurvTot',pDsurvTot

         DICsurv1=DEVsurv1bar+2.0d0*pDsurv1
         DICsurv2=DEVsurv2bar+2.0d0*pDsurv2
         DICsurv3=DEVsurv3bar+2.0d0*pDsurv3
         DICsurv4=DEVsurv4bar+2.0d0*pDsurv4
         DICsurv5=DEVsurv5bar+2.0d0*pDsurv5
         DICsurvTot=DEVsurvTotbar+2.0d0*pDsurvTot
         
         write(*,*) 'DICsurv1',DICsurv1
         write(*,*) 'DICsurv2',DICsurv2
         write(*,*) 'DICsurv3',DICsurv3
         write(*,*) 'DICsurv4',DICsurv4
         write(*,*) 'DICsurv5',DICsurv5
         write(*,*) 'DICsurvTot',DICsurvTot


      open(unit=12,file='BSOoutput.out',!"Bayes Surv Only (No Cure)"
     1       access='sequential',status='unknown')
       write(12,*) 'Bayes Surv Only'
       write(12,*) 'sample size n=',n
c       write(12,*) 'actual sample size nstar=',nstar
       write(12,*) 'number of Gibbs iterations nrep=',nrep
       write(12,*) 'Number of burn-in iterations,'
       write(12,*) 'bin=',nbin
       write(12,*) 'Number of thinning=',nthin
       write(12,*) 'J=',ng 
       write(12,*) '-------------------------------------'
       write(12,*) 'MLE for lambdas'
       write(12,*) '-------------------------------------'
       write(12,*) 'j    s(j)'
       do j=1,ng(1)
        write(12,1005) j,s(1,j)
1005    format(1x,I5,1x,1f16.7)
       enddo
       do j=1,ng(2)
        write(12,1005) j,s(2,j)
       enddo
       do j=1,ng(3)
        write(12,1005) j,s(3,j)
       enddo
       do j=1,ng(4)
        write(12,1005) j,s(4,j)
       enddo
       do j=1,ng(5)
        write(12,1005) j,s(5,j)
       enddo
       do j=1,ng6
        write(12,1005) j,s6(j)
       enddo
       do k=1,5
        write(12,*) 'ntt=',ntt(k)
       enddo 

       write(12,*) '-------------------------------------'
       write(12,*) 'j1  j2    ebeta       sebeta        95% HPD Int.'

       do j1=1,5
        do j2=1,8
        write(12,1001) j1,j2,ebeta(j1,j2),sebeta(j1,j2),
     1     betalow(j1,j2),betaupp(j1,j2)
1001     format(1x,2I4,2f16.8,'(',f12.6,',',f12.6,')')
c 901   continue
        enddo
       enddo

!       write(12,*) '-------------------------------------'
!       write(12,*) 'j1  j2    ebetaC       sebetaC        95% HPD Int.'
!       do j1=1,5
!        do j2=1,9
!        write(12,1002) j1,j2,ebetaC(j1,j2),sebetaC(j1,j2),
!     1     betaClow(j1,j2),betaCupp(j1,j2)
!1002     format(1x,2I4,2f16.8,'(',f12.6,',',f12.6,')')
!        enddo
!       enddo
   
       write(12,*) '------------------------------------------'
       write(12,*) 'j      ealam       sealam         95% HPD Int.'
        do j1=1,5
         do j2=1,ng(j1)
         write(12,1003) j1,j2,ealam(j1,j2),sealam(j1,j2),
     1    alamlow(j1,j2),alamupp(j1,j2)
1003     format(1x,2I5,1x,2f16.8,'(',f12.6,',',f12.6,')')
         enddo
        enddo

c       if(compcind .eq. 1.0d0) then
        write(12,*) '------------------------------------------'
        write(12,*) 'j      ecindex       secindex        95% HPD Int.'
        do j1=1,5
         write(12,1004) j1,ecindex(j1),secindex(j1),
     1    cindexlow(j1),cindexupp(j1)
1004     format(1x,I5,5f16.8,'(',f12.6,',',f12.6,')')
        enddo

c        write(12,*) '------------------------------------------'
c        write(12,*) 'j      ewcindex      sewcindex'     
c        do j1=1,5
c         write(12,1009) j1,ewcindex(j1),secindex(j1)!,
c!     1    cindexlow(j1),cindexupp(j1)
c1009     format(1x,I5,5f16.8)!,'(',f12.6,',',f12.6,')')
c        enddo
c       endif

c       if(compDICsurv .eq. 1.0d0) then
       write(12,*) '------------------------------------------'
       write(12,*) 'barDEVsurv1=',barDEVsurv1
       write(12,*) 'barDEVsurv2=',barDEVsurv2
       write(12,*) 'barDEVsurv3=',barDEVsurv3
       write(12,*) 'barDEVsurv4=',barDEVsurv4
       write(12,*) 'barDEVsurv5=',barDEVsurv5
       write(12,*) 'barDEVsurvTot=',barDEVsurvTot
       write(12,*) 'DEVsurv1bar=',DEVsurv1bar
       write(12,*) 'DEVsurv2bar=',DEVsurv2bar
       write(12,*) 'DEVsurv3bar=',DEVsurv3bar
       write(12,*) 'DEVsurv4bar=',DEVsurv4bar
       write(12,*) 'DEVsurv5bar=',DEVsurv5bar
       write(12,*) 'DEVsurvTotbar=',DEVsurvTotbar
       write(12,*) 'pDsurv1=',pDsurv1
       write(12,*) 'pDsurv2=',pDsurv2
       write(12,*) 'pDsurv3=',pDsurv3
       write(12,*) 'pDsurv4=',pDsurv4
       write(12,*) 'pDsurv5=',pDsurv5
       write(12,*) 'pDsurvTot=',pDsurvTot
       write(12,*) 'DICsurv1=',DICsurv1
       write(12,*) 'DICsurv2=',DICsurv2
       write(12,*) 'DICsurv3=',DICsurv3
       write(12,*) 'DICsurv4=',DICsurv4
       write(12,*) 'DICsurv5=',DICsurv5
       write(12,*) 'DICsurvTot=',DICsurvTot
       write(12,*) '------------------------------------------'
c       write(12,*) 'lpml',alpml
c       endif

       write(12,*) '------------------------------------------'
       call idate(ntoday)
       call itime(now)
       write(12,*)
     1  'date: month/day/year and time: hour, minute, and second'
       write(12,*)
     1  'Begining at date=',ntodayb(2),'/',ntodayb(1),'/',ntodayb(3)
       write(12,*) 'time at',nowb(1),':',nowb(2),':',nowb(3)
       write(12,*)
     1   'Ending at date=',ntoday(2),'/',ntoday(1),'/',ntoday(3)
       write(12,*) 'time at',now(1),':',now(2),':',now(3)
       total = etime(elapsed)
       write(12, *) 'elapsed time in minutes'
       write(12,2216) total/60.0d0,elapsed(1)/60.0d0,elapsed(2)/60.0d0
2216   format('end: total=',f12.4,' user=',f12.4,
     1         ' system=', f12.4)
       write(12,*) '--------------------------------------------------'
       close(12)

       stop
       end

       include 'optim1.f'
       include 'gilks2.f'
       include 'hpd.f'
       include 'tnorm.f'
c      include 'GaussianHermite.f'
       include 'gibbsBSO_nocure.f' ! newest gibbs code (was 'gibbsJointNew_mc_newalgo.f')
