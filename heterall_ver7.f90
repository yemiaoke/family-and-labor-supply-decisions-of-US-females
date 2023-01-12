include 'mkl_vsl.f90'

program main
use mpi; USE MKL_VSL_TYPE; USE MKL_VSL
use state_space; use idxes; use params; use birth_prob; use errors; use inputs; use logl_array; use bhhh; use parallel

implicit none
TYPE (VSL_STREAM_STATE) :: stream
real ( kind = 8 ) bivnor
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
call CPU_TIME(start)
open(1, file = 'iters86.txt', status = 'old')
open(2, file = 'std.txt', status = 'old')
open (10, file='data4.csv') ! read data
do i=1,26500
    read(10,*), data_array_t(i,:)
enddo
close(10)
data_array = transpose(data_array_t)

! 63 parameters
theta0 = (/betae, betam, betak, betab, betaf, betap, muf, mup, betat, sig_labf, sig_labp, am, ak, aw, ae, etah, lamdam/)
theta_new = theta0; theta_old = theta0
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
! mpi parallelization begins
call MPI_Init(ierr)
call MPI_Comm_rank(MPI_COMM_WORLD, myid, ierr)
call MPI_Comm_size(MPI_COMM_WORLD, nprocs, ierr)

call RANDOM_SEED()
allocate(flags(nprocs))
! scatter sample(28,26500) to each process, i.e. into nprocs pieces of sample_sect
! Partition, i.e. get local number of agents and columns in the sample section
nagents = 1060/nprocs
if (mod(1060,nprocs)>myid) then; nagents = nagents + 1; endif
ncols = nagents*25
! Compute partitioned array
allocate(data_sect(28,ncols)); allocate(data_sect_t(ncols,28))

! Build arrays for mpi_scatterv
! rcounts containes all nsect's
! displs containes displacements of partitions in terms of columns
allocate(rcounts(nprocs),displs(nprocs))
displs(1) = 0
do i = 1, nprocs
    rcounts(i) = 1060/nprocs
    if (mod(1060,nprocs).gt.(i-1)) then; rcounts(i) = rcounts(i)+1; endif
    rcounts(i) = rcounts(i)*25
    if ((i-1).ne.0) then; displs(i) = displs(i-1)+rcounts(i-1); endif
enddo

! Convert from number of columns to number of elements
ncols = 28 * ncols
rcounts = 28 * rcounts
displs = 28 * displs

! Scatter array from root
call MPI_SCATTERV(data_array, rcounts, displs, MPI_REAL, data_sect, ncols, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

! reshape data sections
data_sect_t=transpose(data_sect)
allocate(sample(nagents,25,28))
do i=1,nagents
    sample(i,:,:)=data_sect_t(((i-1)*25+1):25*i,:)
enddo
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
! alternative set
i = 0
do e = 0, 1
    do m = 0, 1
        do f = 0, 1
            do p = 0, 1 
                do b = 0, 1    
                    if (f==1 .and. p==1) then; cycle; endif
                    i = i+1
                    alt_set(:,i) = (/e,m,f,p,b/)
                enddo
            enddo
        enddo
    enddo
enddo
j = 1
do i = 1,24
    m = alt_set(2,i)
    if (m==1) then; cycle; endif
    alt_set12(j) = i; j = j+1
enddo

! arrays for error terms
allocate (norm_lab(nagents,10,25,2), err_lab(nagents,10,25,2), uni_lab(nagents,10,25,2))
n_lab = nagents*10*25*2
allocate (norm_all(n_lab), uni_all(n_lab))

! initialize mkl rng
brng = VSL_BRNG_MT19937
errcode = vslnewstream(stream, brng, iseed+myid)

! call rng
method = VSL_RNG_METHOD_UNIFORM_STD_ACCURATE
errcode = vsRngUniform(method, stream, n_lab, uni_all, 0., 1.)
call vscdfnorminv(n_lab, uni_all, norm_all)

! random normal numbers
uni_lab = reshape(uni_all, (/nagents,10,25,2/))
norm_lab = reshape(norm_all, (/nagents,10,25,2/))

! arrays for storing log likelihood computed
allocate (logli_sect(64,nagents))
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
nf_errors = 100.; inewton = 0; ibhhh = 0; mome1 = 0.; mome2 = 0.; adam1 = 0.9; adam2 = 0.999; step_size = 0.
if (myid==0) then
    ifder = 1; linesearch = 0; stop_flag = 0; goback = 0; ll_old = -huge(1.)
endif

do while (maxval(abs(nf_errors))>=1e-6)
    call MPI_BCAST(theta_new, 63, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(theta_old, 63, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ifder, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(linesearch, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(goback, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(ibhhh, 1, MPI_INT, 0, MPI_COMM_WORLD, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    paraloop: do ibeta = 1, 64
        !if (myid==0) then; print *, 'ibeta', ibeta-1; endif
        theta_grad = theta_new
        if (ibeta>1) then
            !if (abs(theta_grad(ibeta-1))>0) then
                !delta(ibeta-1) = abs(theta_grad(ibeta-1))/100
            !else
                !delta(ibeta-1) = 0.01
            !endif
            delta(ibeta-1) = -max(abs(theta_grad(ibeta-1))/100, 0.005)
            theta_grad(ibeta-1) = theta_grad(ibeta-1) + delta(ibeta-1)
        end if
        betae = theta_grad(1:6); betam = theta_grad(7:14); betak = theta_grad(15:23); betab = theta_grad(24); betaf = theta_grad(25:29)
        betap = theta_grad(30:34); muf = theta_grad(35:42); mup = theta_grad(43:50); betat = theta_grad(51:53)
        sig_labf = theta_grad(54); sig_labp = theta_grad(55)
        am = theta_grad(56); ak = theta_grad(57); aw = theta_grad(58); ae = theta_grad(59); etah = theta_grad(60)
        lamdam = theta_grad(61:63)
        
        betam(2) = abs(betam(2)); betak(4) = abs(betak(4)); betak(7) = abs(betak(7)); betak(8:9) = -abs(betak(8:9))
        betab = -abs(betab); betaf(4:5) = abs(betaf(4:5)); betap(4:5) = abs(betap(4:5))
        muf(2) = abs(muf(2)); muf(4:6) = abs(muf(4:6)); muf(7:8) = -abs(muf(7:8))
        mup(2) = abs(mup(2)); mup(4:6) = abs(mup(4:6)); mup(7:8) = -abs(mup(7:8))
        betat = abs(betat); sig_labf = abs(sig_labf); sig_labp = abs(sig_labp); etah = abs(etah)
        
        penalty1 = 0; penalty2 = 0
        if (etah>1) then; penalty1 = (etah-1); etah = 1.; endif
        ! if (abs(betab)>0.) then; penalty2 = abs(betab)/10; endif
        err_lab(:,:,:,1) = norm_lab(:,:,:,1)*sig_labf; err_lab(:,:,:,2) = norm_lab(:,:,:,2)*sig_labp
!************************************************************************************************************************************
!************************************************************************************************************************************ 
!************************************************************************************************************************************
!************************************************************************************************************************************
        do iagent = 1,nagents ! iterate through agents
            ! characteristics
            blk = sample(iagent,1,9); hsp = sample(iagent,1,10); medu = sample(iagent,1,3); nsib = sample(iagent,1,6); mful = sample(iagent,1,11)
            mpar = sample(iagent,1,12); finc = sample(iagent,1,2)/10000; pqual = sample(iagent,1,4)/24; asv = sample(iagent,1,5)/100
            ! choice
            et = sample(iagent,:,21); mt = sample(iagent,:,23); ft = sample(iagent,:,16); pt = sample(iagent,:,17); bt = sample(iagent,:,27)
            ! state
            edut = sample(iagent,:,20); mdurt = sample(iagent,:,24); wexpt = sample(iagent,:,28); preg = sample(iagent,:,14)
            knumt = sample(iagent,:,13); kaget = sample(iagent,:,15)
            if (nsib>6) then; nsib = 6; endif
            do t = 1,25
                if (mdurt(t) > 4) then; mdurt(t) = 4; endif
                if (knumt(t) < 1) then; kaget(t) = 0; endif
                if (kaget(t) > 18) then; kaget(t) = 18; endif
                if (edut(t) > 18) then; edut(t) = 18; endif
            enddo
            ! wages
            wft = sample(iagent,:,18); wpt = sample(iagent,:,19); wht = sample(iagent,:,25)                                           
!************************************************************************************************************************************
!************************************************************************************************************************************
            ! terminal value   
            do itype = 1,2 ! iterate through types
                do iedu = 1, 13
                    do imdur = 1, 5
                        do iwexp = 1, 31
                            emax_all(iedu,imdur,iwexp,:,:) = dot_product(betat, (/iedu-1., imdur-1., (iwexp-1)*.5/)) 
                        enddo
                    enddo
                enddo
                do t = 25,1,-1 ! backward induction
                    if (mt(t)<0) then ! if1
                        logl_it(t) = 0.
                    else
                        ! fitted wages
                        wt_fit=wage_fit(wexpt(t), edut(t))
                        if (wht(t)>0) then
                            wh = wht(t)
                        else
                            wh = hwage_fit(edut(t))
                        endif
!************************************************************************************************************************************
!************************************************************************************************************************************
                        do i=1,24 ! iterate through discrete choices
                            e = alt_set(1,i); m = alt_set(2,i); f = alt_set(3,i); p = alt_set(4,i); b = alt_set(5,i)
                            if (e==et(t) .and. m==mt(t) .and. f==ft(t) .and. p==pt(t) .and. b==bt(t)) then
                                alt_t = i
                            endif
                            ! birth probability
                            probi = birthprob(b, t+15)
                            ! married women or single women who received a marriage offer, m=1 is for sure in their choice set
                            call find_emax(e, m, f, p, edut(t), mdurt(t), wexpt(t), knumt(t), kaget(t), emax1, emax0)
                            ! calculate utility
                            call u_l(f, p, e, ul); call u_b(b, ub); call u_m(m, f, p, mdurt(t), edut(t), um)
                            call u_k(knumt(t), e, m, f, p, uk); call u_e(e, edut(t), m, ue)
                            ! consumption
                            if (ft(t)>=1) then ! working full-time
                                do r = 1, 10
                                    call cspt_u(m, f, p, wft(t), wt_fit(2)*exp(err_lab(iagent,r,t,2)), wh, knumt(t), kaget(t), cspt_ps(r))
                                    call cspt_u(m, f, p, wft(t), wt_fit(2)*exp(-err_lab(iagent,r,t,2)), wh, knumt(t), kaget(t), cspt_ng(r))
                                enddo
                            elseif (pt(t)>=1) then ! working part-time
                                do r = 1, 10
                                    call cspt_u(m, f, p, wt_fit(1)*exp(err_lab(iagent,r,t,1)), wpt(t), wh, knumt(t), kaget(t), cspt_ps(r))
                                    call cspt_u(m, f, p, wt_fit(1)*exp(-err_lab(iagent,r,t,1)), wpt(t), wh, knumt(t), kaget(t), cspt_ng(r))
                                enddo
                            else
                                do r = 1, 10                                   
                                    call cspt_u(m, f, p, wt_fit(1)*exp(err_lab(iagent,r,t,1)), wt_fit(2)*exp(err_lab(iagent,r,t,2)),&
                                        & wh, knumt(t), kaget(t), cspt_ps(r))
                                    call cspt_u(m, f, p, wt_fit(1)*exp(-err_lab(iagent,r,t,1)), wt_fit(2)*exp(-err_lab(iagent,r,t,2)),&
                                        & wh, knumt(t), kaget(t), cspt_ng(r))
                                enddo
                            endif
                            do r = 1, 10 ! values
                                vri_ps(r,i) = ue+um+uk+ul+ub+cspt_ps(r)+disc*(probi*emax1+(1-probi)*emax0)
                                vri_ng(r,i) = ue+um+uk+ul+ub+cspt_ng(r)+disc*(probi*emax1+(1-probi)*emax0)
                                if (vri_ps(r,i)/=vri_ps(r,i)) then; print *, 'vri_ps nan'; stop; endif
                            enddo
                        enddo ! end loop alternatives
!************************************************************************************************************************************
!************************************************************************************************************************************
                        ! logit prob if receive a marriage offer    
                        do r=1,10
                            vmax_ps(r) = maxval(vri_ps(r,:)); vmax_ng(r) = maxval(vri_ng(r,:))
                            sumexp_ps(r) = sum(exp(vri_ps(r,:)-vmax_ps(r))); sumexp_ng(r) = sum(exp(vri_ng(r,:)-vmax_ng(r)))
                            logpr_ps1(r) = vri_ps(r,alt_t)-vmax_ps(r)-log(sumexp_ps(r))
                            logpr_ng1(r) = vri_ng(r,alt_t)-vmax_ng(r)-log(sumexp_ng(r))
                            logpr(r) = log(exp(logpr_ps1(r)-max(logpr_ps1(r),logpr_ng1(r)))+exp(logpr_ng1(r)-&
                                & max(logpr_ps1(r),logpr_ng1(r))))+max(logpr_ps1(r),logpr_ng1(r))
                        enddo
                        logit_prob1 = log(sum(exp(logpr-maxval(logpr))))+maxval(logpr)-log(20.)
                        
                        if ((t>1 .and. mt(t-1)==1) .or. mt(t)==1) then ! previously married or receive a marriage offer    
                            logit_prob = logit_prob1
                        else ! logit prob if recieve no marriage offer, m=1 is not in choice set
                            j = 1
                            do i = 1,24
                                if (alt_set(2,i)==0) then
                                    vri_single_ps(:,j) = vri_ps(:,i); vri_single_ng(:,j) = vri_ng(:,i)
                                    if (i==alt_t) then; alt_j = j; endif
                                    j = j+1
                                else
                                    cycle
                                endif
                            enddo
                            do r=1,10
                                vmax_ps(r) = maxval(vri_single_ps(r,:)); vmax_ng(r) = maxval(vri_single_ng(r,:))
                                sumexp_ps(r) = sum(exp(vri_single_ps(r,:)-vmax_ps(r))); sumexp_ng(r) = sum(exp(vri_single_ng(r,:)-vmax_ng(r)))
                                logpr_ps1(r) = vri_single_ps(r,alt_j)-vmax_ps(r)-log(sumexp_ps(r))
                                logpr_ng1(r) = vri_single_ng(r,alt_j)-vmax_ng(r)-log(sumexp_ng(r))
                                logpr(r) = log(exp(logpr_ps1(r)-max(logpr_ps1(r),logpr_ng1(r)))+exp(logpr_ng1(r)-&
                                    & max(logpr_ps1(r),logpr_ng1(r))))+max(logpr_ps1(r),logpr_ng1(r))
                                if (logpr_ps1(r)>0) then
                                    print *, logpr_ps1(r), vri_single_ps(r,alt_j), vmax_ps(r), log(sumexp_ps(r)), j, alt_t
                                    print *, vri_single_ps(r,:)
                                    stop
                                endif
                            enddo
                            logit_prob2 = log(sum(exp(logpr-maxval(logpr))))+maxval(logpr)-log(20.)    
                            mprob = marprob(t+15.)
                            loglit_m = max(logit_prob1,logit_prob2)
                            logit_prob = loglit_m + log(mprob*exp(logit_prob1-loglit_m) + (1-mprob)*exp(logit_prob2-loglit_m))
                        endif
                        
                        if (logit_prob>huge(1.)) then; print *, 'logit_prob inf', mt(t); stop; endif
                        if (logit_prob>0) then; print *, 'logit_prob positive', j, logit_prob1, logit_prob2; stop; endif
                        if (logit_prob/=logit_prob) then; print *, 'logit_prob error', mt(t); stop;endif
!************************************************************************************************************************************
!************************************************************************************************************************************
                    ! log likelihood contribution of agent i in period t
                    ! add wage shocks probability
                    probtf = 0; probtp = 0; lprobt = 0
                    if (ft(t)>=1) then
                        probtf = log_normal_pdf(log(wft(t))-log(wt_fit(1)), sig_labf)
                        lprobt = probtf
                    elseif (pt(t)>=1) then
                        probtp = log_normal_pdf(log(wpt(t))-log(wt_fit(2)), sig_labp)
                        lprobt = probtp
                    endif
                        
                    ! probability of giving birth
                    logprobt = 0    
                    if (t>=1) then
                        if (preg(t)>=1) then
                            gbprobt = birthprob(bt(t), t+15); logprobt = log(gbprobt)
                        endif
                    endif
                    if (logprobt/=logprobt .or. abs(logprobt)>huge(1.)) then; print *, 'logprobt error'; endif
                        
                    ! probability of a single woman receiving a marriage offer
                    logmprob = 0
                    if ((t==1 .and. mt(t)==1) .or. (t>1 .and. mt(t-1)==0 .and. mt(t)==1)) then
                        mprob = marprob(t+15.); logmprob = log(mprob)
                    endif
                    logl_it(t) = logit_prob + 10*logprobt + 10*lprobt + 10*logmprob
                endif ! end if1
                if (logl_it(t) /= logl_it(t)) then
                    print *, 'logit nan error'; print *, logprobt, lprobt
                endif
                if (abs(logl_it(t)) > huge(1.)) then
                    print *,  'logit inf error'; print *, logprobt, lprobt
                endif
!************************************************************************************************************************************
!************************************************************************************************************************************                
!************************************************************************************************************************************
                    ! update emax for grid state (edu, mdur, wexp, knum, kage)
                    do iedu=1,6
                        do iwexp=1,5
                            wt_fit_emax(iedu,iwexp,:) = wage_fit(wexp(iwexp),edu(iedu))
                        enddo
                    enddo
                    do iedu=1,6
                        do iwexp=1,5
                            do r=1,10
                                wf_emax_ps(iedu,iwexp,r) = wt_fit_emax(iedu,iwexp,1)*exp(err_lab(iagent,r,t,1))
                                wp_emax_ps(iedu,iwexp,r) = wt_fit_emax(iedu,iwexp,2)*exp(err_lab(iagent,r,t,2))
                                wf_emax_ng(iedu,iwexp,r) = wt_fit_emax(iedu,iwexp,1)*exp(-err_lab(iagent,r,t,1))
                                wp_emax_ng(iedu,iwexp,r) = wt_fit_emax(iedu,iwexp,2)*exp(-err_lab(iagent,r,t,2))
                                if (abs(wf_emax_ps(iedu,iwexp,r))>huge(1.)) then; print *, 'wf_emax_ps inf'; stop; endif
                            enddo
                        enddo
                    enddo                
                    do iedu=1,6 ! husband wage
                        wh_emax(iedu) = hwage_fit(edu(iedu))
                    enddo
                    do e = 0,1 ! ue
                        do m = 0,1
                            do iedu=1,6
                                call u_e(e, edu(iedu), m, ue_emax(e+1, m+1, iedu))
                            enddo
                        enddo
                    enddo
                    do i = 1,24 ! consumption
                        m = alt_set(2,i); f = alt_set(3,i); p = alt_set(4,i)
                        do iedu=1,6
                            do iwexp=1,5
                                do iknum=1,4
                                    do ikage = 1,4
                                        do r=1,10
                                            call cspt_u(m, f, p, wf_emax_ps(iedu,iwexp,r), wp_emax_ps(iedu,iwexp,r), wh_emax(iedu),&
                                                &knum(iknum), kage(ikage), c_emax_ps(iedu,iwexp,iknum,ikage,r,i))
                                            call cspt_u(m, f, p, wf_emax_ng(iedu,iwexp,r), wp_emax_ng(iedu,iwexp,r), wh_emax(iedu),&
                                                &knum(iknum), kage(ikage), c_emax_ng(iedu,iwexp,iknum,ikage,r,i))
                                            if (abs(c_emax_ps(iedu,iwexp,iknum,ikage,r,i))>huge(1.)) then
                                                print *, 'c_emax_ps inf'; stop
                                            endif
                                        enddo
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
                    do i=1,24 ! uk
                        e = alt_set(1,i); m = alt_set(2,i); f = alt_set(3,i); p = alt_set(4,i)
                        do iknum=1,4
                            call u_k(knum(iknum), e, m, f, p, uk_emax(iknum, i))
                        enddo
                    enddo
                    do i = 1,24 ! um
                        m = alt_set(2,i); f = alt_set(3,i); p = alt_set(4,i)
                        do imdur = 1,3
                            do iedu = 1,6
                                call u_m(m, f, p, mdur(imdur), edu(iedu), um_emax(imdur,iedu,i))
                            enddo
                        enddo
                    enddo
                    do i=1,24 ! ul
                        e = alt_set(1,i); f = alt_set(3,i); p = alt_set(4,i)
                        call u_l(f, p, e, ul_emax(i))
                        if (ul_emax(i)/=ul_emax(i)) then
                            print *, 'ul_emax nan'; stop
                        endif
                    enddo
!************************************************************************************************************************************                                              
                    do iknum=1,4
                        do ikage=1,4
                            do imdur=1,3
                                do iedu=1,6
                                    do iwexp=1,5
                                        do i=1,24
                                            e = alt_set(1,i); m = alt_set(2,i); f = alt_set(3,i); p = alt_set(4,i); b = alt_set(5,i)
                                            call u_b(b, ub)
                                            probi = birthprob(b, t+15)
                                            utl = ue_emax(e+1,m+1,iedu)+um_emax(imdur,iedu,i)+ul_emax(i)+ub+uk_emax(iknum,i)
                                            call find_emax(e,m,f,p,edu(iedu),mdur(imdur),wexp(iwexp),knum(iknum),kage(ikage),emax1,emax0)
                                            do r=1,10
                                                v_emax_ps(r,i) = utl+c_emax_ps(iedu,iwexp,iknum,ikage,r,i)+&
                                                    & disc*(probi*emax1+(1-probi)*emax0)
                                                v_emax_ng(r,i) = utl+c_emax_ng(iedu,iwexp,iknum,ikage,r,i)+&
                                                    & disc*(probi*emax1+(1-probi)*emax0)
                                                if (v_emax_ps(r,i) /= v_emax_ps(r,i) .or. v_emax_ng(r,i) /= v_emax_ng(r,i)) then
                                                    print *, 'v_emax_ps error'
                                                    stop
                                                endif
                                                if (abs(v_emax_ps(r,i))>huge(1.)) then
                                                    print *, 'v_emax_ps inf'
                                                endif
                                            enddo
                                        enddo
                                        do r=1,10
                                            vmax_ps(r)=maxval(v_emax_ps(r,:)); vmax_ng(r)=maxval(v_emax_ng(r,:))
                                            sumexp_ps(r) = sum(exp(v_emax_ps(r,:)-vmax_ps(r)))
                                            sumexp_ng(r) = sum(exp(v_emax_ng(r,:)-vmax_ng(r)))
                                            emaxr_ps(r)=0.5772+log(sumexp_ps(r))+vmax_ps(r)
                                            emaxr_ng(r)=0.5772+log(sumexp_ng(r))+vmax_ng(r)
                                            if (emaxr_ps(r) /= emaxr_ps(r)) then; print *, 'emaxps nan'; stop; endif
                                            if (emaxr_ng(r) /= emaxr_ng(r)) then; print *, 'emaxng nan'; stop; endif
                                        enddo
                                        emax_m1 = sum(emaxr_ps+emaxr_ng)/20 
                                        
                                        j = 1
                                        do i = 1,24
                                            if (alt_set(2,i)==0) then
                                                vri_single_ps(:,j) = v_emax_ps(:,i); vri_single_ng(:,j) = v_emax_ng(:,i)
                                                j = j+1
                                            else
                                                cycle
                                            endif
                                        enddo
                                        do r=1,10
                                            vmax_ps(r)=maxval(vri_single_ps(r,:)); vmax_ng(r)=maxval(vri_single_ng(r,:))
                                            sumexp_ps(r) = sum(exp(vri_single_ps(r,:)-vmax_ps(r)))
                                            sumexp_ng(r) = sum(exp(vri_single_ng(r,:)-vmax_ng(r)))
                                            emaxr_ps(r)=0.5772+log(sumexp_ps(r))+vmax_ps(r)
                                            emaxr_ng(r)=0.5772+log(sumexp_ng(r))+vmax_ng(r)
                                            if (emaxr_ps(r) /= emaxr_ps(r)) then; print *, 'emaxps nan'; stop; endif
                                            if (emaxr_ng(r) /= emaxr_ng(r)) then; print *, 'emaxng nan'; stop; endif
                                        enddo
                                        emax_m2 = sum(emaxr_ps+emaxr_ng)/20 
                                        
                                        if (imdur>1) then
                                            emax(iedu,imdur,iwexp,iknum,ikage) = emax_m1
                                        else
                                            mprob = marprob(t+15.)
                                            emax(iedu,imdur,iwexp,iknum,ikage) = emax_m1*mprob + emax_m2*(1-mprob)
                                        endif
                                        if (emax(iedu,imdur,iwexp,iknum,ikage) /= emax(iedu,imdur,iwexp,iknum,ikage)) then
                                            print *, 'emax nan 1:', iagent,t; stop
                                        endif
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
!************************************************************************************************************************************
!************************************************************************************************************************************
                    ! interpolate other state
                    do iedu = 6, 18
                        call loc_search(edu, 6, real(iedu), idx_edu); edul = edu(idx_edu(1)); eduh = edu(idx_edu(2))
                        do imdur = 0, 4
                            call loc_search(mdur, 3, real(imdur), idx_mdur); mdurl = mdur(idx_mdur(1)); mdurh = mdur(idx_mdur(2))
                            do iwexp = 0, 30
                                call loc_search(wexp, 5, real(iwexp*0.5), idx_wexp); wexpl = wexp(idx_wexp(1)); wexph = wexp(idx_wexp(2))
                                do iknum = 0, 7
                                    call loc_search(knum, 4, real(iknum), idx_knum); knuml = knum(idx_knum(1)); knumh = knum(idx_knum(2))
                                    do ikage = 0, 18
                                        call loc_search(kage, 4, real(ikage), idx_kage); kagel = kage(idx_kage(1)); kageh = kage(idx_kage(2))
                                        do i = 1, 2
                                            do j = 1, 2
                                                do k = 1, 2
                                                    do l = 1, 2
                                                        do r = 1, 2
                                                            v_grid(i,j,k,l,r) = emax(idx_edu(i), idx_mdur(l), idx_wexp(j), idx_knum(r),&
                                                                & idx_kage(k))
                                                            if (edul /= eduh) then
                                                                z_edu = iedu-(edul+eduh-edu(idx_edu(i)))
                                                            else
                                                                z_edu = 1
                                                            endif
                                                            if (wexpl /= wexph) then
                                                                z_wexp = iwexp*0.5-(wexpl+wexph-wexp(idx_wexp(j)))
                                                            else
                                                                z_wexp = 1
                                                            endif
                                                            if (kageh /= kagel) then
                                                                z_kage = ikage-(kageh+kagel-kage(idx_kage(k)))
                                                            else
                                                                z_kage = 1
                                                            end if
                                                            if (mdurh /= mdurl) then
                                                                z_mdur = imdur-(mdurh+mdurl-mdur(idx_mdur(l)))
                                                            else
                                                                z_mdur = 1
                                                            end if
                                                            if (knumh /= knuml) then
                                                                z_knum = iknum-(knumh+knuml-knum(idx_knum(r)))
                                                            else
                                                                z_knum = 1
                                                            end if
                                                            weight(i,j,k,l,r) = (z_edu*z_edu)*(z_wexp*z_wexp)*(z_kage*z_kage)*(z_mdur*z_mdur)*&
                                                                & (z_knum*z_knum)
                                                        enddo
                                                    enddo
                                                enddo
                                            enddo
                                        enddo
                                        sum_weight = sum(weight)
                                        weight = weight/sum_weight; aux_emax = v_grid*weight
                                        emax_all(iedu-5, imdur+1, iwexp+1, iknum+1, ikage+1) = sum(aux_emax)
                                    enddo
                                enddo
                            enddo
                        enddo
                    enddo
!************************************************************************************************************************************
!************************************************************************************************************************************
                enddo ! end loop periods
                ! loglikelihood contribution of agent i conditional on type
                loglit_type(itype) = sum(logl_it)
            enddo ! end loop types        
            loglit_m = maxval(loglit_type)
            logli_sect(ibeta, iagent) = loglit_m + log(etah*exp(loglit_type(1)-loglit_m) + (1-etah)*exp(loglit_type(2)-loglit_m))
            logli_sect(ibeta, iagent) = logli_sect(ibeta, iagent)-penalty1-penalty2
        enddo ! end loop agents
        if (ifder<=0) then; exit; endif ! just calculate loglikelihood, not grads
    enddo paraloop! end loop parameters
    !print *, 'id', myid, 'paraloop finished'    
    call mpi_barrier(MPI_COMM_WORLD, ierr)
!************************************************************************************************************************************
!************************************************************************************************************************************
    ! gather loglikelihood from each slave process to the master process(0)
    displs(1) = 0
    do i = 1, nprocs
        rcounts(i) = 1060/nprocs
        if (mod(1060,nprocs).gt.(i-1)) then; rcounts(i) = rcounts(i)+1; endif
        if ((i-1).ne.0) then; displs(i) = displs(i-1)+rcounts(i-1); endif
    enddo
    ncols = 64 * nagents
    rcounts = 64 * rcounts
    displs = 64 * displs

    call mpi_barrier(MPI_COMM_WORLD, ierr)
    call mpi_gatherv(logli_sect, ncols, MPI_REAL, logli, rcounts, displs, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call mpi_barrier(MPI_COMM_WORLD, ierr)
    !print *, 'id', myid, 'gatherv finished'
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
    ! bhhh updates on master process only
    if (myid==0) then
        ibhhh = ibhhh+1
        logli_t = transpose(logli)
        ll_theta(1) = sum(logli_t(:,1)); ll_new = ll_theta(1)

        print *, 'bhhh iterations:', ibhhh, inewton
        print *, 'ifder', ifder
        print *, 'logL_old:'
        print *, ll_old
        print *, 'logL_new:'
        print *, ll_new
        print *, 'parameters:'
        print '(5f20.5)', theta_new
    
        ! calculate numerical derivatives if flag 'ifder=1'
        if (ifder>0) then
            inewton = inewton+1; linesearch = 0; goback = 0; grad_indiv = 0.; grad= 0.
            do i = 2,64
                ll_theta(i) = sum(logli_t(:,i))
            enddo
            do iagent=1,1060
                do igrad = 1,63
                    grad_indiv(iagent,igrad) = (logli_t(iagent,igrad+1)-logli_t(iagent,1))/delta(igrad)
                enddo
            enddo
            grad = sum(grad_indiv, dim=1)/1060
            print *, 'grad:'
            print '(5f20.5)', grad
            
            ! approximate hessian
            do iagent=1,1060
                hess_indiv(iagent,:,:) = outer_product(grad_indiv(iagent,:))
            enddo
            hess = sum(hess_indiv, dim=1)/1060
            !print *, 'hess:'
            !print '(10f20.5)', hess
            ! inverse hessian
            hinv = hess
            call SGETRF(63,63,hinv,63,ipiv,info)
            if (info.ne.0) then; print *, 'Matrix is numerically singular'; endif
            call SGETRI(63,hinv,63,ipiv,work,63,info)
            if (info.ne.0) then; print *, 'Matrix inversion failed'; endif
            do i = 1,63; beta_sig(i) = sqrt(abs(hinv(i,i))); enddo
            write(2, '(5f20.5)') beta_sig
            write(2,*) '**********************************************************'
            write(2,*) '**********************************************************'
            write(2,*) '**********************************************************'
            
            mome1 = adam1*mome1+(1-adam1)*grad; mome2 = adam2*mome2+(1-adam2)*grad*grad
            m1hat = mome1/(1-adam1**inewton); m2hat = mome2/(1-adam2**inewton)
            ! step length
            
            if (step_size <0.01) then; step_size = step_size+0.002; endif
            !if (inewton>10) then; step_size = 0.01; endif
            step = step_size*m1hat/(sqrt(m2hat)+1e-8)
            !step = step_size*grad
            !step = matmul(hinv, grad)
            if (norm2(step)>1) then; step = step/norm2(step); endif
        endif
        
        write(1,*) 'bhhh iterations:', ibhhh, inewton
        write(1,*) 'logL:', ll_new
        write(1,*) 'gradient:'
        write(1,'(5f)') grad
        write(1,*) 'parameters:'
        write(1,'(5f20.5)') theta_new
        write(1,*) '**********************************************************'
    
        if (ll_new > ll_old) then
            stop_flag = 0
            if (goback==0 .and. ifder==0) then; linesearch = linesearch + 1; endif
            if (goback==0 .and. linesearch<3) then
                nf_errors = step
                theta_old = theta_new; ll_old = ll_new; theta_new = theta_new + step; ifder = 0
            else
                ifder = 1
                theta_old = theta_new; ll_old = -huge(1.); nf_errors = 100
            endif
        else
            if (linesearch>0) then
                theta_new = theta_old; ifder = 1; ll_old = -huge(1.)
                nf_errors = 100
            else
                if (goback<7) then
                    goback = goback + 1
                    print *, 'number of gobacks:', goback
                    step = step/2
                    theta_new = theta_old+step
                    linesearch = 0; ifder = 0
                    nf_errors = 100
                else
                    ll_theta(1) = ll_old
                    logmax_loc = maxloc(ll_theta, dim=1)
                    if (logmax_loc==1) then
                        if (stop_flag==1) then
                            print *, 'find max 1'
                            call CPU_TIME(finish)
                            open(3, file = 'end86.txt', status = 'old')
                            write(3,*) 'time used:', finish-start
                            write(3,*) 'logL:'
                            write(3,*) ll_old
                            write(3,*) 'parameters:'
                            write(3,'(5f20.5)') theta_old
                            write(3,*) 'newton updates:', inewton
                            write(3,*) 'bhhh iterations:', ibhhh
                            stop
                        endif
                        stop_flag = 1
                        theta_new = theta_old - step
                    else
                        theta_new = theta_old; theta_new(logmax_loc-1) = theta_old(logmax_loc-1) + delta(logmax_loc-1)
                    endif
                    ifder = 1
                    if (inewton<=20) then
                        nf_errors = 100
                    else
                        nf_errors = theta_old - theta_new
                    endif
                endif
            endif                
        endif
        print *, 'line search:', linesearch
        print *, 'step:', norm2(step)
        print '(5f20.5)', step
        print *, 'new parameters:'
        print '(5f20.5)', theta_new
        
        
        
    endif  
    ! broadcast differences from master process to slave process
    call MPI_BCAST(nf_errors, 63, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
enddo
if (myid==0) then
    print *, 'find max 2'
    call CPU_TIME(finish)
    open(3, file = 'end86.txt', status = 'old')
    write(3,*) 'time used:', finish-start
    write(3,*) 'logL:'
    write(3,*) ll_old
    write(3,*) 'parameters:'
    write(3,'(5f20.5)') theta_old
    write(3,*) 'newton updates:', inewton
    write(3,*) 'bhhh iterations:', ibhhh
endif
call MPI_Finalize(ierr) 
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
!************************************************************************************************************************************
contains
! utility
subroutine u_e(e,edu,m,ue)
implicit none
integer :: e, ifhg, m
real :: edu, age
real, intent(out) :: ue

ifhg = 0
if (edu>=12) then; ifhg = 1; endif
age = (t+15)/10.

if (e>=1) then
    ue = dot_product(betae, (/1., real(medu), finc, real(ifhg), age, real(m)/)) + (itype-1)*ae
else
    ue = 0
end if    
end subroutine u_e
!****************************************************   
subroutine u_m(m, f, p, mdur, edu, um)
implicit none
integer :: m, f, p, ifcg
real :: mdur, edu
real, intent(out):: um

ifcg = 0
if (edu>=16) then; ifcg = 1; endif

if (m>=1) then
    um = dot_product(betam, (/1., mdur, real(f), real(p), pqual, real(blk), real(hsp), real(ifcg)/)) + (itype-1)*am
else
    um = 0
endif
end subroutine u_m
!*******************************************************
function marprob(age)
implicit none
real:: age, z, marprob

z = dot_product(lamdam, (/1., age/10., age*age/100./))
marprob = 1/(1+exp(z))
end function marprob
!*******************************************************
subroutine u_k(knum, e, m, f, p, uk)
implicit none
real :: knum, part1, age
integer:: e, m, f, p
real, intent(out) :: uk

age = (t+15.)/10.
if (knum>=1) then
    part1 = dot_product(betak(1:8), (/1., real(nsib), real(e), real(m), real(f), real(p), age, age*age/)) + (itype-1)*ak
    uk = part1*knum + betak(9)*knum*knum
else
    uk = 0
end if
end subroutine u_k
!********************************************************
function birthprob(b, age)
implicit none
integer, intent(in):: b, age
real:: birthprob, z

if (b==1) then
    birthprob = 0
else
    z = dot_product(lamda, (/ 1., age/10., real(blk), real(hsp), blk*age/10., hsp*age/10. /))
    birthprob = 1/(1+exp(z))
endif
end function birthprob
!********************************************************************************************************************************** 
subroutine u_l(f, p, e, ul)
implicit none
integer f, p, e
real uf, up
real, intent(out):: ul

if (f>=1) then
    uf = dot_product(betaf, (/1., real(mful), real(e), real(blk), real(hsp)/))
else
    uf = 0
endif

if (p>=1) then
    up = dot_product(betap, (/1., real(mpar), real(e), real(blk), real(hsp)/))
else
    up = 0
end if
ul = uf+up
end subroutine u_l
!**********************************************************
subroutine u_b(b, ub)
implicit none
integer b
real, intent(out) :: ub

if (b>0) then
    ub = betab
else
    ub = 0
endif
end subroutine u_b
!**********************************************************
subroutine cspt_u(m, f, p, wf, wp, wh, knum, kage, cspt)
implicit none
integer m, f, p; real wf, wp, wh, knum, kage, inc, kcost, dispinc
real , intent(out) :: cspt

if (m>=1) then
    inc = wh + f*wf + p*wp
else
    inc = f*wf + p*wp
end if

kcost = 0
if (knum>0 .and. kage<18) then
    if (kage<5) then
        kcost = knum*0.101*inc
    else
        kcost = knum*0.067*inc
    endif
endif
if (f==0 .and. p==0) then
    kcost = 0.
elseif (p==1) then
    kcost = 0.5*kcost
endif

if (m>=1) then
    dispinc = (inc-kcost)*0.5
else
    dispinc = inc-kcost
end if

if (dispinc>1) then; cspt = log(dispinc); else; cspt = 0; endif
end subroutine cspt_u
!*************************************************************************************************************************************
! fitted wages
function wage_fit(wexp, edu)
implicit none
real wage_fit(2), wexp, edu
integer ifhsg, ifcg

ifhsg = 0; ifcg = 0
if (edu>=12 .and. edu<16) then; ifhsg = 1; endif
if (edu>=16) then; ifcg = 1; endif    

wage_fit(1) = exp(dot_product(muf, (/1., wexp/10., (edu-6.)/10., real(ifhsg), real(ifcg), asv, real(blk), real(hsp)/))&
    & + (itype-1)*aw)
wage_fit(2) = exp(dot_product(mup, (/1., wexp/10., (edu-6.)/10., real(ifhsg), real(ifcg), asv, real(blk), real(hsp)/))&
    & + (itype-1)*aw) 
end function wage_fit
!***************************************************************************************************************************
function hwage_fit(grade)
implicit none
real hwage_fit, grade, age, edu

edu = grade-6
age = (t+15)/10.

hwage_fit = exp(dot_product(muh, (/1., age, asv, edu, real(blk), real(hsp)/)))
end function hwage_fit
!***************************************************************************************************************************
subroutine find_emax(e, m, f, p, edux, mdurx, wexpx, knumx, kagex, emax1, emax0)
implicit none
integer, intent(in) :: e, m, f, p
real, intent(in) :: edux, mdurx, wexpx, knumx, kagex
real, intent(out) :: emax1, emax0
integer iedu1, imdur1, iwexp1, iknum1, iknum0, ikage0
real edu1, mdur1, wexp1, knum1, kage0, emax10, emax11, emax12, emax00, emax01, emax02

edu1 = min(edux+e, 18.); iedu1=int(edu1-5)

if (m==0) then; mdur1 = 0; else; mdur1 = min(mdurx+m, 4.); endif; imdur1=int(mdur1+1)

wexp1 = min(wexpx+f+0.5*p, 15.); iwexp1=int(2*wexp1+1)

iknum0 = int(knumx+1); knum1 = min(knumx+1, 7.); iknum1 = int(knum1+1)
if (knumx>=1) then; kage0 = min(kagex+1, 18.); else; kage0 = 0; endif; ikage0 = int(kage0+1)

emax0 = emax_all(iedu1, imdur1, iwexp1, iknum0, ikage0); emax1 = emax_all(iedu1, imdur1, iwexp1, iknum1, 1)
end subroutine find_emax
!********************************************************************************************************************************** 
function log_normal_pdf(x, sig)
implicit none
real, intent(in) :: x, sig
real :: log_normal_pdf, pi
pi = 3.141593
log_normal_pdf = -0.5*(x/sig)*(x/sig) - log(sig*sqrt(2*pi))
end function log_normal_pdf

function norm_pdf(x)
implicit none
real norm_pdf, x, pi
pi = 3.141593
norm_pdf = exp(-x*x/2)/sqrt(2*pi)
end function norm_pdf

function outer_product(a)
implicit none
real, dimension(:) :: a
real, allocatable :: outer_product(:,:)
integer i, j, n

n = size(a)
allocate (outer_product(n,n))
do i = 1, n
    do j = 1, n
        outer_product(i,j) = a(i)*a(j)
    enddo
enddo
end function outer_product
!**********************************************************************************************************************************
pure function linspace(from, to1, n)
real, intent(in):: from, to1
integer, intent(in):: n
real, allocatable:: linspace(:)
real :: range
integer :: i
    
range = to1 - from
allocate(linspace(n))
    
if (n == 1) then; linspace(1) = from; endif

do i = 1, n
    linspace(i) = from + range * (i - 1) / (n - 1)
enddo
end function
!**********************************************************************************************************************************
subroutine loc_search(array, n, value, idx)
implicit none
integer, intent(in) :: n
real, intent(in) :: array(n), value
integer, intent(out) :: idx(2)
integer :: i

i = 1
do while (i < n+1)
    if (value==array(i)) then
        idx(1) = i; idx(2) = i
        exit
    else if (value<array(i)) then
        idx(1) = i-1; idx(2)=i
        exit
    else
        i = i+1
    end if
end do
end subroutine loc_search
!**********************************************************************************************************************************
!**********************************************************************************************************************************
!**********************************************************************************************************************************
end program main