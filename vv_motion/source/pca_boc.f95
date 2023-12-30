module bocMod
  use matvec
  implicit none
  logical cflag ! flag to process cmh
  double precision,allocatable, dimension(:,:):: boc,boc0
  double precision,allocatable, dimension(:,:,:):: tva,tvb,cm
  double precision,dimension(2,4,3):: triad_avg, triad_ref
  double precision,dimension(2,2,3):: cm_avg, cm_ref
  ! frm0,frm1: first/last frame for analysis. From command line input
  !   Default: frm0=1, frm1=-1 (last frame)
  ! nframe, dfrm: determined in read_boc_data()
  ! pca_frm: number of frames to generate to show PCA motion (default=0)
  integer nframe, frm0, frm1, dfrm, pca_frm
  !integer runType;
  !! runType=1: perform PCA,
  !!        =2: read existing PC from files and project BOC data
  character (len=128):: idir,iprefx,osuffix,pdir
  ! pdir: directory where PCA data are to be read

  ! nboc number of boc data per frame (12 lines * 3 coord per line)
  integer, parameter:: nboc=12*3 
  integer, parameter:: nb=9 ! (3 arms)*(3 coords) (centroids are unused)
  integer, parameter:: nb2=2*nb
  ! tlen: length of triad's arm (proc_triad.py:write_pdb())
  double precision, parameter:: tlen=15. 

  ! PCA-related
  double precision, dimension(:), allocatable:: W, Sigma
  double precision, dimension(:), allocatable:: pca_std, pca_std_cm
  double precision, dimension(:,:), allocatable:: U, VT
  double precision pca_a(nb2,3,3), pca_b(nb2,3,3) 
  double precision pca_a_ref(nb2,3,3), pca_b_ref(nb2,3,3),ptheta
 ! For PCA of cm data
  double precision, dimension(:), allocatable:: W_cm, Sigma_cm
  double precision, dimension(:,:), allocatable:: U_cm, VT_cm
  ! pca_cm*(nb2,:,3): middle index for const-hinge-var
  double precision pca_cma(nb2,3,3), pca_cmb(nb2,3,3)
  double precision pca_cma_ref(nb2,3,3), pca_cmb_ref(nb2,3,3)

contains
  !-------------------------------------------------------
  subroutine read_boc_data(ifname)
    !use bocMod
    character(len=128),intent(in):: ifname
    integer i,j,idum
    character(len=32) :: sdum,sdum1
    double precision rdum
    
    open (unit=11, file=trim(ifname),status='old',action='read')
    read (11,*) sdum,sdum1,nframe ! read 1st line
    if (frm1 .eq. -1 ) frm1=nframe
    if (frm1 .le. frm0 ) then
       print *,'ERROR: frm1 .le. frm0'
       stop
    endif

    dfrm=frm1-frm0 +1
    
    allocate(boc(dfrm,nboc))
    allocate(boc0(dfrm,nboc))
    allocate(tva(dfrm,4,3))
    allocate(tvb(dfrm,4,3))
    allocate(cm(dfrm,4,3)) ! 4: ca,ha,cb,hb - 3: xyz-coord
    allocate(pca_std(nb2))
    allocate(pca_std_cm(nb2))
    
    do i = 1,3 ! read next 3 lines
       read (11,'(A)')
    enddo
    do i = 1, frm1 ! nframe
       read (11,*) sdum, idum
       !read (11,*)(boc(i,j),j=1,nboc)
       if (i .ge. frm0) then
          read (11,*)(boc(i-frm0+1,j),j=1,nboc)
       else
          read(11,*) (rdum, j=1,nboc)
       endif
    enddo
    boc0 = boc ! keep original data

  end subroutine read_boc_data
  
  !-------------------------------------------------------
  subroutine read_cm_data() ! call this after read_boc_data()
    !use bocMod
    character(len=2), dimension(4):: cfname
    integer i,j,idum,k
    character(len=32) :: sdum,sdum1
    double precision rdum

    cfname=(/ 'ca', 'ha', 'cb', 'hb' /)
    do k = 1, 4 ! read in 4 cm files
       sdum1=trim(idir)//"/cm_"//cfname(k)//".dat"
       print *,'Reading: ',sdum1
       open (unit=11, file=trim(sdum1),status='old',action='read')
       read (11,'(A)') sdum ! read 1st line
       do i = 1, frm1 ! nframe
          read (11,*) idum
          if (i .ne. idum) then
             print *,'ERROR: frame number mismatch.'
             stop
          endif
          if (i .ge. frm0) then
             read (11,*)(cm(i-frm0+1,k,j),j=1,3)
          else
             read(11,*) (rdum, j=1,3)
          endif
       enddo
       close(11)
    enddo
    !boc0 = boc ! keep original data

  end subroutine read_cm_data

  !-------------------------------------------------------
  subroutine read_std(ofn0,tmp_std)
    ! call this when post-processing w/o running pca
    double precision, dimension(:), intent(inout):: tmp_std
    integer i,j,idum,k,iost
    character(len=32) :: sdum
    character(len=256) :: sdum1
    character(len=128),intent(in):: ofn0    
    double precision rdum

    sdum=trim(idir)//"/pca"//trim(osuffix)//"/std"//trim(ofn0)//".dat"
    open (unit=12, file=sdum,iostat=iost,action='read')
    print *,'Reading: ',trim(sdum)

    tmp_std(k)=-1. ! initialize
    do
       read(12,'(A)',iostat=iost) sdum1
       if (iost > 0)  then
          print *,'error in reading ',trim(sdum)
          call exit(-1)
       else if (iost < 0) then ! eof
          exit
       else
          if (index(sdum1,'#')>0) cycle ! comment
          read(sdum1,*) i,tmp_std(i),rdum
       end if
    end do
    print *,'Number of data read: ',i
    
  end subroutine read_std
  
  !-------------------------------------------------------
  subroutine get_triad(flag)
    ! get triads from boc data & average
    character flag
    integer i,j
    logical Ldum,lsame
    Ldum=lsame(flag,'a')
    !do i=1,nframe ! assign to tva & tvb
    do i=1,dfrm ! assign to tva & tvb
       if (Ldum) then
          do j=1,4
             tva(i,j,:)=boc(i,(j-1)*3+1:j*3)
          enddo
       else
          do j=7,10
             tvb(i,j-6,:)=boc(i,(j-1)*3+1:j*3)
          enddo
       endif
    enddo

    do i=1,4
       do j=1,3
          if (Ldum) then
             triad_avg(1,i,j)=sum(tva(:,i,j)) / dfrm
          else
             triad_avg(2,i,j)=sum(tvb(:,i,j)) / dfrm
          endif
       enddo
    enddo
    
  end subroutine get_triad
  
  !-------------------------------------------------------
  subroutine triad_project(triad1,triad2)
    ! project triad1 to triad2. triad1(:,2:4,:) is overwritten
    integer i,j,k
    double precision,intent(inout), dimension(:,:,:):: triad1
    double precision,intent(in), dimension(:,:,:):: triad2    
    double precision, dimension(3):: vdum
    double precision rdum
    
    do i=1,dfrm
       rdum=0.D0
       do j=2,4
          do k=1,3
             vdum(k)=dot_product(triad1(i,j,:),triad2(i,k+1,:))
          enddo
          triad1(i,j,:)=vdum
          rdum = rdum + dot_product(vdum,vdum)
       enddo
       rdum = sqrt(rdum)
       do j=2,4
          triad1(i,j,:) = triad1(i,j,:) / rdum
       enddo
    enddo
    
  end subroutine triad_project
  
  !-------------------------------------------------------
  subroutine pca_cmh()
    ! perform PCA on cm data
    integer i,j,k, INFO, LWORK
    double precision f1(2,3)
    double precision, dimension(:,:),allocatable:: pca_tmp
    logical, dimension(2,3):: m1
    character(len=128):: sdum

    ! nb2 can be used here: (3 cm)*(xyz)*(tcra,tcrb)
    LWORK=5*dfrm
    allocate(W_cm(5*dfrm)); allocate(U_cm(dfrm,nb2))
    allocate(VT_cm(nb2,nb2)); allocate(Sigma_cm(nb2))
    
    allocate(pca_tmp(dfrm,nb2))

    do i=1,2 ! cm of const & hinge
       do j=1,3 ! xyz
          cm_avg(1,i,j)=sum(cm(:,i,j)) / dfrm ! ca,ha
          cm_avg(2,i,j)=sum(cm(:,i+2,j)) / dfrm ! cb,hb
       end do
    end do
    call write_cm_avg_pdb()

    !!!!!!!!!!!!
    ! prepare for PCA
    do i=1,dfrm ! pack arrays
       pca_tmp(i, 1:6 )=pack(cm(i,1:2,:)-cm_avg(1,1:2,:),.true.) ! TCRa
       pca_tmp(i, 7:9 )=pack(tva(i,1,:)-triad_avg(1,1,:),.true.)
       pca_tmp(i, 10:15)=pack(cm(i,3:4,:)-cm_avg(2,1:2,:),.true.) ! TCRb
       pca_tmp(i, 16:nb2 )=pack(tvb(i,1,:)-triad_avg(2,1,:),.true.)
    enddo

    ! perform SVD
    ! Arguments explained in: ~/computer/ftns/lapack/lapack-3.8.0/SRC/dgesvd.f
    CALL DGESVD( 'S', 'A', dfrm,nb2,pca_tmp,dfrm,Sigma_cm,U_cm, dfrm, VT_cm, &
         nb2, W_cm, LWORK, INFO )

    ! unpack into respective principal components
    m1= .true.; f1=0.;
    do i=1,nb2 ! i: PCA index
       ! VT_cm indices follow the order in pack() above
       pca_cma(i,1:2,:) = unpack(VT_cm(i,1:6),mask=m1,field=f1)        
       pca_cma(i,3,:) = VT_cm(i,7:9)
       pca_cmb(i,1:2,:) = unpack(VT_cm(i,10:15),mask=m1,field=f1)
       pca_cmb(i,3,:) = VT_cm(i,16:nb2)
    enddo

    ! print v-c pca
    do i=1,nb2 ! i: PCA index
       print *, "pc ", i
       print *, pca_cma(i,:,:) ! pca(pc1,triad arm,coors)
       print *, pca_cmb(i,:,:)              
    enddo
    
    ! Write data
    call write_pca_dir_cm_vmd()
    !call write_pca_dir_cm_dat()

    sdum="cm"
    call write_pca_dir(pca_cma,pca_cmb,Sigma_cm,sdum)
    sdum="_cm"
    call write_pca_traj_std(U_cm,Sigma_cm,10,sdum,pca_std_cm)
  end subroutine pca_cmh

  !-------------------------------------------------------
  subroutine pca_triad()
    ! perform PCA on triads
    integer i,j,k, INFO, LWORK
    double precision f1(3,3)
    double precision, dimension(:,:),allocatable:: pca_tmp
    logical, dimension(3,3):: m1
    character(len=128):: sdum
    
    LWORK=5*dfrm
    allocate(W(5*dfrm)); allocate(U(dfrm,nb2))
    allocate(VT(nb2,nb2)); allocate(Sigma(nb2))
    
    allocate(pca_tmp(dfrm,nb2))
    pca_tmp=0.

    call write_triad_avg_pdb()

    !!!!!!!!!!!!
    ! prepare for PCA
    !do i=frm0,nframe ! pack arrays
    m1=.true.
    do i=1,dfrm ! pack arrays
       pca_tmp(i,   1:nb  )=pack(tva(i,2:4,:)-triad_avg(1,2:4,:),m1) !.true.)
       pca_tmp(i,nb+1:nb2 )=pack(tvb(i,2:4,:)-triad_avg(2,2:4,:),m1) !.true.)

    enddo

    ! perform SVD
    ! Arguments explained in: ~/computer/ftns/lapack/lapack-3.8.0/SRC/dgesvd.f
    CALL DGESVD( 'S', 'A', dfrm,nb2,pca_tmp,dfrm,Sigma,U, dfrm, VT, &
         nb2, W, LWORK, INFO )

    ! unpack into respective principal components
    m1= .true.; f1=0.
    do i=1,nb2
       pca_a(i,:,:) = unpack(VT(i,1:nb),mask=m1,field=f1)        
       pca_b(i,:,:) = unpack(VT(i,nb+1:nb2),mask=m1,field=f1)
       !print *,norm2(VT(i,:)) ! check norms are 1
       print *, "pc ", i
       print *, pca_a(i,:,:) ! pca(pc1,triad arm,coors)
       print *, pca_b(i,:,:)       
    enddo

    ! Write data
    call write_pca_dir_vmd()
    sdum="triad"
    call write_pca_dir(pca_a,pca_b,Sigma,sdum)
    sdum=""
    call write_pca_traj_std(U,Sigma,10,sdum,pca_std)
  end subroutine pca_triad

  !-------------------------------------------------------
  subroutine project_boc(nb0)
    integer, intent(in):: nb0 ! number of modes to project
    integer i,j,k,p
    character(len=128):: sdum    
    double precision, dimension(:,:,:),allocatable:: tva0,tvb0,cma0,cmb0
    double precision, dimension(nb0)::tProj,cmProj

    allocate(tva0(dfrm,3,3))
    allocate(tvb0(dfrm,3,3))
    allocate(cma0(dfrm,3,3)) ! ca,ha,va - 3: xyz-coord
    allocate(cmb0(dfrm,3,3)) ! cb,hb,vb - 3: xyz-coord
    
    do i=1,dfrm ! subtract avg coords of ref set
       do j=1,3 ! triad arms
          do k=1,3 
             tva0(i,j,k)=tva(i,j+1,k)-triad_ref(1,j+1,k)
             tvb0(i,j,k)=tvb(i,j+1,k)-triad_ref(2,j+1,k)
          enddo
       enddo
       ! for const/hinge/var, subtract ref coords in order to see
       ! spatial displacement
       do j=1,2 ! const/hinge
          do k=1,3
             cma0(i,j,k)=cm(i,j,  k)-cm_ref(1,j,k)
             cmb0(i,j,k)=cm(i,j+2,k)-cm_ref(2,j,k)
          enddo
       enddo
       do k=1,3 ! var
          cma0(i,3,k)=tva(i,1,k)-triad_ref(1,1,k)
          cmb0(i,3,k)=tvb(i,1,k)-triad_ref(2,1,k)
       enddo
    enddo

    sdum=trim(pdir)//"/proj_triad.dat"
    print *,"Writing: ",trim(sdum)
    open (unit=11, file=trim(sdum),status='replace',action='write')
    sdum=trim(pdir)//"/proj_cm.dat"
    print *,"Writing: ",trim(sdum)
    open (unit=12, file=trim(sdum),status='replace',action='write')
    write (11,'("# frm0,dfrm,frm1,nframe ",4I6)') frm0,dfrm,frm1,nframe
    write (12,'("# frm0,dfrm,frm1,nframe ",4I6)') frm0,dfrm,frm1,nframe
    

    do i=1,dfrm
       do p=1,nb0 ! project
          tProj(p)=0.; cmProj(p)=0.
          do j=1,3
             do k=1,3
                tProj(p)=tProj(p)+tva0(i,j,k)*pca_a_ref(p,j,k)
                tProj(p)=tProj(p)+tvb0(i,j,k)*pca_b_ref(p,j,k)
                cmProj(p)=cmProj(p)+cma0(i,j,k)*pca_cma_ref(p,j,k)
                cmProj(p)=cmProj(p)+cmb0(i,j,k)*pca_cmb_ref(p,j,k)
             enddo
          enddo
       enddo ! do p=1,nb0 ! project
       do p=1,nb0 ! write
          write (11,'(F8.4)',advance='no') tProj(p)
          write (12,'(F8.4)',advance='no') cmProj(p)
       enddo
       write (11,'("")'); write (12,'("")')
    enddo ! do i=1,dfrm
    close(11); close(12)
  end subroutine project_boc
  !-------------------------------------------------------
  subroutine write_pca_dir_vmd()
    ! Write principal dir into vmd file
    integer i,j,k,idum
    double precision tdum(2,3,3)
    character(len=128):: sdum,file_id
    character(len=32):: cname(5) ! color names
    character(len=:),allocatable:: str ! for cname
        
    str="blue red orange cyan2 ochre"
    read(unit=str,fmt=*) cname ! set color names

    do i=1,nb2
       write(file_id,'(I2.2)') i ! add leading zero
       sdum=trim(idir)//"/pca"//trim(osuffix)//"/dir"//trim(file_id)//".vmd"
       open (unit=12, file=sdum,status='replace',action='write')
       write (12,'("# Rank ",I2," singular_val= ",E10.3)') i,Sigma(i)
       write (12,'("mol new")')
       write (12,'("graphics top color",A12)') trim(cname(1+modulo(i-1,5)))
       do k=1,2
          if (k==1) then
             tdum(1,:,:)=tlen*(triad_avg(k,2:4,:)+1.0*pca_a(i,:,:))
             tdum(2,:,:)=tlen*(triad_avg(k,2:4,:)-1.0*pca_a(i,:,:))
          else
             tdum(1,:,:)=tlen*(triad_avg(k,2:4,:)+1.0*pca_b(i,:,:))
             tdum(2,:,:)=tlen*(triad_avg(k,2:4,:)-1.0*pca_b(i,:,:))
          endif
          do j=1,3
             tdum(1,j,:) = tdum(1,j,:) + triad_avg(k,1,:)
             tdum(2,j,:) = tdum(2,j,:) + triad_avg(k,1,:)
          enddo
          
          do j=1,3
             !write (12,'("draw cylinder { ")',advance="no")
             write (12,'("draw arrow { ")',advance="no")
             write (12,'(3F7.2)',advance="no") (tdum(2,j,idum),idum=1,3)
             write (12,'(" } { ")',advance="no")
             write (12,'(3F7.2," }")') (tdum(1,j,idum),idum=1,3)
             !write (12,'(" } radius 0.3")')
          enddo
       enddo   ! do k=1,2
       write(sdum,'(I2.2)') i
       write (12,'("mol rename top pca",A2)') adjustl(sdum)
       close(12)
    enddo
  end subroutine write_pca_dir_vmd

  !-------------------------------------------------------
  subroutine read_pca_dir(p_a,p_b,ifn)
    ! Read principal dir of p_a,p_b from a text file
    double precision,dimension(nb2,3,3),intent(out):: p_a,p_b
    
    integer i,j,k
    character(len=128):: sdum
    character(len=128),intent(in)::ifn

    print *,"Reading: ",trim(ifn)
    open (unit=12, file=trim(ifn),status='old',action='read')
    read(12,'(A)') sdum ! read 1st line
    do i=1,nb2
       read (12,'(A)') sdum
       if (index(sdum,"Rank")==0) then ! assert
          print *, 'ERROR: Rank information missing'
          stop
       end if
       do j=1,3
          read (12,'(3E20.10)') (p_a(i,j,k),k=1,3)
       enddo
       do j=1,3
          read (12,'(3E20.10)') (p_b(i,j,k),k=1,3)
       enddo
    enddo
    close(12)
  end subroutine read_pca_dir

  !-------------------------------------------------------
  subroutine write_pca_dir(p_a,p_b,Sigma0,ofn)
    ! Write principal dir of p_a,p_b into a text file
    double precision,dimension(nb2,3,3),intent(in):: p_a,p_b
    double precision, intent(in):: Sigma0(nb2)
    
    integer i,j,k
    character(len=128):: sdum,sdum1
    character(len=16),intent(in)::ofn
        
    sdum=trim(idir)//"/pca"//trim(osuffix)//"/dir_"//trim(ofn)//".dat"
    print *,"Writing: ",trim(sdum)
    if (trim(ofn) .eq. "triad") then
       sdum1="# row 1-3: Va, row 4-6: Vb"
    else 
       sdum1="# row 1-3: TCRa [const, hinge, var], row 4-6: ditto for TCRb"
    end if
    open (unit=12, file=sdum,status='replace',action='write')
    write(12,'(72A)') trim(sdum1)
    do i=1,nb2
       write (12,'("# Rank ",I2," singular_val= ",E10.3)') i,Sigma0(i)
       do j=1,3
          write (12,'(3E20.10)',advance="no") (p_a(i,j,k),k=1,3)
          write (12,'("")')
       enddo
       do j=1,3
          write (12,'(3E20.10)',advance="no") (p_b(i,j,k),k=1,3)
          write (12,'("")')
       enddo
    enddo
    close(12)
  end subroutine write_pca_dir

  !-------------------------------------------------------
  subroutine write_pca_dir_cm_vmd()
    ! Write principal dir of cm data into vmd file
    integer i,j,k,idum
    double precision tdum(2,3,3)
    character(len=128):: sdum,file_id
    character(len=32):: cname(5) ! color names
    character(len=:),allocatable:: str ! for cname
        
    str="blue red orange cyan2 ochre"
    read(unit=str,fmt=*) cname ! set color names

    do i=1,nb2
       write(file_id,'(I2.2)') i ! add leading zero
       sdum=trim(idir)//"/pca"//trim(osuffix)//"/dir"//trim(file_id)//"_cm.vmd"
       open (unit=12, file=sdum,status='replace',action='write')
       write (12,'("# Rank ",I2," singular_val= ",E10.3)') i,Sigma_cm(i)
       write (12,'("mol new")')
       write (12,'("graphics top color",A12)') trim(cname(1+modulo(i-1,5)))
       do k=1,2
          if (k==1) then
             ! constant & hinge 
             tdum(1,1:2,:)=cm_avg(k,1:2,:)+tlen*pca_cma(i,1:2,:)
             tdum(2,1:2,:)=cm_avg(k,1:2,:)-tlen*pca_cma(i,1:2,:)
             ! variable domain
             tdum(1,3,:)=triad_avg(k,1,:)+tlen*pca_cma(i,3,:)
             tdum(2,3,:)=triad_avg(k,1,:)-tlen*pca_cma(i,3,:)
          else
             ! constant & hinge 
             tdum(1,1:2,:)=cm_avg(k,1:2,:)+tlen*pca_cmb(i,1:2,:)
             tdum(2,1:2,:)=cm_avg(k,1:2,:)-tlen*pca_cmb(i,1:2,:)
             ! variable domain
             tdum(1,3,:)=triad_avg(k,1,:)+tlen*pca_cmb(i,3,:)
             tdum(2,3,:)=triad_avg(k,1,:)-tlen*pca_cmb(i,3,:)
          endif
          
          do j=1,3
             !write (12,'("draw cylinder { ")',advance="no")
             write (12,'("draw arrow { ")',advance="no")
             write (12,'(3F7.2)',advance="no") (tdum(2,j,idum),idum=1,3)
             write (12,'(" } { ")',advance="no")
             write (12,'(3F7.2," }")') (tdum(1,j,idum),idum=1,3)
          enddo
       enddo   ! do k=1,2
       write(sdum,'(I2.2)') i
       write (12,'("mol rename top pca_cm",A2)') adjustl(sdum)
       close(12)
    enddo
  end subroutine write_pca_dir_cm_vmd

  !-------------------------------------------------------
  subroutine write_triad_avg_vmd()
    ! write average of triad into a vmd file
    integer idum,j,k
    character(len=128):: sdum
    double precision vdum(3)
    
    !sdum="./data_"//trim(idir)//"/avg_triad.vmd"
    sdum=trim(idir)//"/pca"//trim(osuffix)//"/avg_triad.vmd"
    print *,"Writing: ",trim(sdum)
    open (unit=12, file=sdum,status='replace',action='write')
    write(12,'("mol new")')
    write(12,'("graphics top color ochre")')
    do k=1,2
       do j=2,4
          vdum=triad_avg(k,1,:)+triad_avg(k,j,:)*tlen
          write (12,'("draw sphere { ")',advance="no") 
          write (12,'(3F8.2)',advance="no") (vdum(idum),idum=1,3)
          write(12,'(" } radius 1.0")')
       enddo
    enddo
    write (12,'("mol rename top avg")')
    close(12)
  end subroutine write_triad_avg_vmd
  
  !-------------------------------------------------------
  subroutine read_boc_avg(idir0,flag)
    ! Read average of cm & triad
    integer, intent(in):: flag
    integer:: i,j,k
    character(len=128),intent(in):: idir0
    character(len=128):: sdum

    sdum=trim(idir0)//"/boc_avg.dat"

    print *,"Reading: ",trim(sdum)
    open (unit=12, file=trim(sdum),status='old',action='read')
    read (12,'(A)') sdum
    do i=1,2
       if (cflag) then
          do j=1,2
             if (flag==0) then
                read (12,'(3E20.10)') (cm_ref(i,j,k),k=1,3)
             else
                read (12,'(3E20.10)') (cm_avg(i,j,k),k=1,3)
             endif
          enddo
       endif
       do j=1,4
          if (flag==0) then
             read (12,'(3E20.10)') (triad_ref(i,j,k),k=1,3)
          else
             read (12,'(3E20.10)') (triad_avg(i,j,k),k=1,3)
          endif
       enddo
    enddo
    close(12)
  end subroutine read_boc_avg
  
  !-------------------------------------------------------
  subroutine write_boc_avg()
    ! write average of cm & triad into a text file
    integer:: i,j,k
    character(len=128):: sdum

    sdum=trim(idir)//"/pca"//trim(osuffix)//"/boc_avg.dat"
    open (unit=12, file=sdum,status='replace',action='write')
    if (cflag) then
       write (12,'(A)') "# row 1-6: cm/ch/triad of TCRa, row: 7-12: same for TCRb"
    else
       write (12,'(A)') "# row 1-8: triad of TCRa/b"
    endif
    print *,"Writing: ",trim(sdum)
          
    do i=1,2
       if (cflag) then 
          do j=1,2
             write (12,'(3E20.10)',advance="no") (cm_avg(i,j,k),k=1,3)
             write (12,'("")')
          enddo
       endif
       do j=1,4
          write (12,'(3E20.10)',advance="no") (triad_avg(i,j,k),k=1,3)
          write (12,'("")')
       enddo
       !do k=1,3
       !print *,'e3: ', (triad_avg(i,4,k),k=1,3)  
       !enddo
    enddo
    close(12)
  end subroutine write_boc_avg
  
  !-------------------------------------------------------
  subroutine write_cm_avg_pdb()
    ! write average of cm into psf/pdb file
    integer:: idum,i,j,k,p
    integer:: natom=6 ! 2*(c-h-v)
    double precision:: occcupancy=0.5
    character(len=128):: sdum,segi(2),sdum1,chainID(2)
    double precision vdum(3)

    segi(1)='TCRA  '; segi(2)='TCRB  '
    chainID(1)='A'; chainID(2)='B'
    sdum=trim(idir)//"/pca"//trim(osuffix)//"/avg_cm.psf"
    open (unit=11, file=sdum,status='replace',action='write')

    write(11,'("PSF")')
    write(11,*)
    write(11,'(" !NTITLE")')
    write(11,'("*")')
    write(11,'(I6," !NATOM")') natom
    
    sdum=trim(idir)//"/pca"//trim(osuffix)//"/avg_cm.pdb"
    open (unit=12, file=sdum,status='replace',action='write')
    write(12,'("MODEL")')

100 FORMAT(A6,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3,2F6.2,6X,A4)
110 FORMAT(I8,1X,A4,1X,I4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8) !fm01 in psfres.src:psfwr2()
    do k=1,2
       do j=1,3
          if (j==3) then
             vdum=triad_avg(k,1,:)
          else
             vdum=cm_avg(k,j,:)
          endif
          p=j+3*(k-1)
          write(12,100) 'ATOM  ',p,' Z  ',' CM',trim(chainID(k)),j, &
               (vdum(idum),idum=1,3),1.0,0.5,segi(k)

          write(11,110) p,segi(k),j,' CM','Z',0,0.0,0.0,0
       enddo
    enddo

    write(11,*)
    write(11,'("       4 !NBOND")')
    i=0
    do k=1,2
       p=1+3*(k-1)
       do j=1,3
          if (j .ne. 3) then
             write(11,'(2I8)',advance="no") p+j-1,p+j
             i=i+1
          endif
          if (modulo(i,4)==0) write(11,*) ! newline
       end do
    end do
    
    write(12,'("ENDMDL")')
    close(11);    close(12)
  end subroutine write_cm_avg_pdb
  
  !-------------------------------------------------------
  subroutine write_triad_avg_pdb()
    ! write average of triad into psf/pdb file
    integer:: idum,i,j,k,p
    integer:: natom=8 ! 2*(centroid+3 arms)
    double precision:: occcupancy=0.5
    character(len=128):: sdum,segi(2),sdum1,chainID(2)
    double precision vdum(3)

    segi(1)='VA  '; segi(2)='VB  '
    chainID(1)='A'; chainID(2)='B'
    sdum=trim(idir)//"/pca"//trim(osuffix)//"/avg_triad.psf"
    open (unit=11, file=sdum,status='replace',action='write')

    write(11,'("PSF")')
    write(11,*)
    write(11,'(" !NTITLE")')
    write(11,'("*")')
    write(11,'(I6," !NATOM")') natom
    
    sdum=trim(idir)//"/pca"//trim(osuffix)//"/avg_triad.pdb"
    open (unit=12, file=sdum,status='replace',action='write')
    write(12,'("MODEL")')

100 FORMAT(A6,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3,2F6.2,6X,A4)
110 FORMAT(I8,1X,A4,1X,I4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8) !fm01 in psfres.src:psfwr2()
    do k=1,2
       do j=1,4
          if (j==1) then
             vdum=triad_avg(k,1,:)
          else
             vdum=triad_avg(k,1,:)+triad_avg(k,j,:)*tlen
          endif
          p=j+4*(k-1)
          write(12,100) 'ATOM  ',p,' S  ','TRI',trim(chainID(k)),j, &
               (vdum(idum),idum=1,3),1.0,0.5,segi(k)

          write(11,110) p,segi(k),j,'TRI','S',0,0.0,0.0,0
       enddo
    enddo

    write(11,*)
    write(11,'("       6 !NBOND")')
    i=0
    do k=1,2
       p=1+4*(k-1)
       do j=1,3
          write(11,'(2I8)',advance="no") p,p+j
          i=i+1
          if (modulo(i,4)==0) write(11,*) ! newline
       end do
    end do
    
    write(12,'("ENDMDL")')
    close(11);    close(12)
  end subroutine write_triad_avg_pdb
  
  !-------------------------------------------------------
  subroutine write_pca_traj_std(U0,Sig0,icut,ofn0,tmp_std)
    ! write principal component trajectories and their std
    ! icut: max number of principal components to write.
    ! icut<=0: write all components
    integer i,j,k,imax
    integer, intent(in):: icut
    double precision, intent(in), dimension(:):: Sig0
    double precision, intent(in), dimension(:,:):: U0
    character(len=128),intent(in):: ofn0
    double precision, intent(inout), dimension(:):: tmp_std
    
    double precision, dimension(:),allocatable:: tmp0,tmp_mean
    character(len=128):: sdum,file_id
    double precision rdum

    imax=icut;  if (icut .lt. 1) imax=nb2
    allocate(tmp_mean(imax))

    do i=1,imax
       write(file_id,'(I2.2)') i ! add leading zero
       sdum=trim(idir)//"/pca"//trim(osuffix)//"/traj"//trim(file_id) &
            //trim(ofn0)//".dat"
       open (unit=12, file=sdum,status='replace',action='write')
       write (12,'("# Rank ",I2," singular_val= ",E10.3)') i,Sig0(i)
       write (12,'("# frm0,dfrm,frm1,nframe ",4I6)') frm0,dfrm,frm1,nframe
       tmp_mean(i)=0.; tmp_std(i)=0.
       do k=1,dfrm
          rdum=U0(k,i)*Sig0(i)
          tmp_mean(i)=tmp_mean(i)+ rdum
          tmp_std(i)=tmp_std(i)+ rdum*rdum
          write(12,'(F8.4)') rdum
       enddo
       close(12)
    enddo

    do i=1,imax
       tmp_mean(i)=tmp_mean(i)/dfrm
       tmp_std(i)=tmp_std(i)/dfrm
       rdum=tmp_mean(i)
       tmp_std(i)=tmp_std(i)-rdum*rdum
       if (tmp_std(i)<0.) then
          print *,'ERROR: negative variance.'
          stop
       endif
       tmp_std(i)=sqrt(tmp_std(i))
    enddo

    sdum=trim(idir)//"/pca"//trim(osuffix)//"/std"//trim(ofn0)//".dat"
    open (unit=12, file=sdum,status='replace',action='write')
    sdum="# Input file: boc_"//trim(iprefx)//".dat"
    write (12,'(A)') trim(sdum)
    write (12,'("# frame_ini, frm_fin, nframe= ",3I6)') frm0, frm1,nframe
    !write (12,'("# nframe= ",I5)') nframe
    write (12,'("# STD and singular val of each principal component:")')
    do i=1,imax
       write (12,'(I2,"  ",E10.3,F10.4)') i, tmp_std(i),Sig0(i)
    enddo
    close(12)

    deallocate(tmp_mean);

  end subroutine write_pca_traj_std
  !-------------------------------------------------------
  
  subroutine write_pca_traj(icut)
    ! write principal component trajectories
    ! icut: max number of principal components to write.
    ! icut<=0: write all components
    integer i,j,k,imax
    integer, intent(in):: icut
    character(len=128):: sdum,file_id
    double precision vdum(3)

    imax=icut;  if (icut .lt. 1) imax=nb2
    do i=1,imax
       write(file_id,'(I2.2)') i ! add leading zero
       !sdum="./data_"//trim(idir)//"/pca/traj"//trim(file_id)//".dat"
       sdum=trim(idir)//"/pca"//trim(osuffix)//"/traj"//trim(file_id)//".dat"
       open (unit=12, file=sdum,status='replace',action='write')
       write (12,'("# Rank ",I2," singular_val= ",F6.3)') i,Sigma(i)
       write (12,'("# frm0,dfrm,frm1,nframe ",4I6)') frm0,dfrm,frm1,nframe
       do k=1,dfrm
          write(12,'(F8.4)') U(k,i)*Sigma(i)
       enddo
       close(12)
    enddo

  end subroutine write_pca_traj
  
  !-------------------------------------------------------
  subroutine write_triad_vmd(vmd_name)
    ! write triad trajectory as a vmd file
    integer i,j,idum
    double precision vdum(3)
    character(len=128),intent(in)::vmd_name
    character(len=16),allocatable:: sdum(:)
    character(len=16):: s1="red",s2="cyan2",s3="yellow2"
    allocate(sdum(3))
    sdum=[character(len=16):: s1,s2,s3]
    open (unit=12, file=trim(vmd_name),status='replace',action='write')
    write(12,'("mol new")')
    do j=2,4
       write(12,'("# Va triad arm ",I3)') j-1
       write(12,'("draw color ",A16)') sdum(j-1)
       do i=1,dfrm
          vdum=tva(i,1,:)+tva(i,j,:)*tlen 
          write (12,'("draw point { ")',advance="no") 
          write (12,'(3F8.2)',advance="no") (vdum(idum),idum=1,3)
          write(12,'(" }")')
       enddo
       write(12,*)
    enddo
    do j=2,4
       write(12,'("# Vb triad arm ",I3)') j-1
       write(12,'("draw color ",A16)') sdum(j-1)
       do i=1,dfrm
          vdum=tvb(i,1,:)+tvb(i,j,:)*15. ! 15: from proc_triad.py:write_pdb()
          write (12,'("draw point { ")',advance="no") 
          write (12,'(3F8.2)',advance="no") (vdum(idum),idum=1,3)
          write(12,'(" }")')
       enddo
       write(12,*)
    enddo

  end subroutine write_triad_vmd
  
  !-------------------------------------------------------
  subroutine write_cm_vmd(vmd_name)
    ! write triad trajectory as a vmd file
    integer i,j,idum
    double precision vdum(3)
    character(len=128),intent(in)::vmd_name
    character(len=16),allocatable:: sdum(:),sdum1(:)
    character(len=16):: s1="orange",s2="yellow2"
    character(len=16):: s11="ca",s12="ha",s13="cb",s14="hb"
    allocate(sdum(2))
    sdum=[character(len=16):: s1,s2]
    allocate(sdum1(4))
    sdum1=[character(len=16):: s11,s12,s13,s14]
    
    open (unit=12, file=trim(vmd_name),status='replace',action='write')
    write(12,'("mol new")')
    do j=1,4
       write(12,'("# subdomain ",A16)') sdum1(j)
       write(12,'("draw color ",A16)') sdum((j+1)/2)
       do i=1,dfrm
          write (12,'("draw point { ")',advance="no") 
          write (12,'(3F8.2)',advance="no") (cm(i,j,idum),idum=1,3)
          write(12,'(" }")')
       enddo
       write(12,*)
    enddo

    deallocate(sdum)
  end subroutine write_cm_vmd
  
  !-------------------------------------------------------
  subroutine write_pca_motion()
    ! write PCA motion into a pdb file.
    
    integer:: idum,j,k,p,q, ipca,ifrm
    integer:: natom=8 ! 2*(centroid+3 arms)
    ! std_max from: ~/Charmm/tcr/1oga/vab/wat0/anal/boc/data_vab/pca1/std
    double precision, parameter:: PI=3.14159265358979,std_max=0.219
    double precision:: occcupancy=0.5, thet0

    character(len=128):: sdum,segi(2),sdum1,chainID(2),file_id
    double precision vdum(3),pdum(3),thet1,rdum,rmax(3)
    ! rotation angle,increment, and axis:
    double precision pca_angle(2,3),pca_dthet(2,3),pca_axis(2,3,3) 

    sdum=''
    call read_std(sdum,pca_std) ! read std values
    
    print *,"Writing: PCA motion PDBs in ",trim(idir)//"/pca1/"

    thet0=ptheta*PI/180.
    segi(1)='VA  '; segi(2)='VB  '
    chainID(1)='A'; chainID(2)='B'
100 FORMAT(A6,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3,2F6.2,6X,A4)
110 FORMAT(I8,1X,A4,1X,I4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8) !fm01 in psfres.src:psfwr2()

    do ipca=1,3 ! create motion for PC1-PC3
       write(file_id,'(I2.2)') ipca ! add leading zero
       sdum=trim(idir)//"/pca"//trim(osuffix)//"/motion_triad"// &
            trim(file_id)//".pdb"
       open (unit=12, file=sdum,status='replace',action='write')
       rmax(ipca)=-1.0
       do j=1,3 ! find max component of ipca (j: triad arm)
          vdum=pca_a(ipca,j,:)
          rdum=dot_product(vdum,vdum)
          if (rdum > rmax(ipca)) rmax(ipca)=rdum
          pca_angle(1,j)=sqrt(rdum)
          vdum=pca_b(ipca,j,:)
          rdum=dot_product(vdum,vdum)
          if (rdum > rmax(ipca)) rmax(ipca)=rdum
          pca_angle(2,j)=sqrt(rdum)
       enddo ! find max component of ipca scaled by std
       rmax(ipca)=pca_std(ipca)/(std_max*sqrt(rmax(ipca)))
          
       do k=1,2
          do j=1,3
             pca_angle(k,j)=pca_angle(k,j)*thet0*rmax(ipca)
             pca_dthet(k,j)=pca_angle(k,j)/pca_frm
             if (k==1) then
                pca_axis(k,j,:)=cross(triad_avg(k,j+1,:),pca_a(ipca,j,:))
             else
                pca_axis(k,j,:)=cross(triad_avg(k,j+1,:),pca_b(ipca,j,:))
             endif
          enddo
       enddo
       
       do ifrm=-pca_frm,pca_frm ! total # of frames: 2*pca_frm+1
          write(12,'("MODEL")')
          do k=1,2
             vdum=triad_avg(k,1,:) ! centroid
             p=1+4*(k-1)
             write(12,100) 'ATOM  ',p,' S  ','TRI',trim(chainID(k)),1, &
                  (vdum(idum),idum=1,3),1.0,0.5,segi(k)

             do j=1,3 ! three arms
                rdum=ifrm*pca_dthet(k,j) ! actual rotation angle
                pdum=rotate_perp(pca_axis(k,j,:),triad_avg(k,j+1,:),rdum)
                vdum=triad_avg(k,1,:)+pdum*tlen
                p=j+1+4*(k-1)
                write(12,100) 'ATOM  ',p,' S  ','TRI',trim(chainID(k)),j+1, &
                     (vdum(idum),idum=1,3),1.0,0.5,segi(k)
             enddo
          enddo
          write(12,'("ENDMDL")')
       enddo ! ifrm=-pca_frm,pca_frm
       close(12)
    enddo ! ipca=1,3 ! create for PC1-PC3
    
  end subroutine write_pca_motion

  !-------------------------------------------------------
  subroutine write_pca_motion_cm()
    ! write PCA motion for cmh into a pdb file. This includes the triad
    
    integer:: idum,j,k,p,q, ipca,ifrm
    integer:: natom=8 ! 2*(centroid+3 arms)
    double precision, parameter:: PI=3.14159265358979
    double precision:: occcupancy=0.5, thet0, std_max(2)

    character(len=128):: sdum,segi(2),sdum1,chainID(2),file_id
    double precision vdum(3),vdum1(3),pdum(3),thet1,rdum,rmax(2)
    double precision vpos(2,3) ! position of variable domains
    ! rotation angle & axis
    double precision angle_triad(2,3),axis_triad(2,3,3)
    double precision angle_cm(2,3),axis_cm(2,2,3),dir_cm(2,2,3)

    ! std_max(1) from: ~/Charmm/tcr/1oga/tab/wat0/anal/boc/data_cab/pca1/std_cm
    ! std_max(2) from: ~/Charmm/tcr/1oga/tab/wat0/anal/boc/data_vab/pca1/std
    std_max=(/ 3.15,0.281 /)
    
    sdum=''
    call read_std(sdum,pca_std) ! read std values
    sdum="_cm"
    call read_std(sdum,pca_std_cm) ! read std values
    
    print *,"Writing: cmh PCA motion PDBs in ",trim(idir)//"/pca1/"

    thet0=ptheta*PI/180.
    segi(1)='VA  '; segi(2)='VB  '
    chainID(1)='A'; chainID(2)='B'
100 FORMAT(A6,I5,1X,A4,1X,A3,1X,A1,I4,4X,3F8.3,2F6.2,6X,A4)
110 FORMAT(I8,1X,A4,1X,I4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8) !fm01 in psfres.src:psfwr2()

    do ipca=1,3 ! create motion for PC1-PC3
       write(file_id,'(I2.2)') ipca ! add leading zero

       sdum=trim(idir)//"/pca"//trim(osuffix)//"/motion_triad"// &
            trim(file_id)//".pdb"
       open (unit=12, file=sdum,status='replace',action='write')

       sdum=trim(idir)//"/pca"//trim(osuffix)//"/motion_cm"// &
            trim(file_id)//".pdb"
       open (unit=14, file=sdum,status='replace',action='write')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! find max component for cmh & triad
       rmax=-1.0 ! rmax(1): max PCA for cmh, rmax(2): max PCA for triad
       do j=1,3 
          vdum=pca_cma(ipca,j,:)
          rdum=dot_product(vdum,vdum)
          if (rdum > rmax(1)) rmax(1)=rdum
          angle_cm(1,j)=sqrt(rdum)
          vdum=pca_cmb(ipca,j,:)
          rdum=dot_product(vdum,vdum)
          if (rdum > rmax(1)) rmax(1)=rdum
          angle_cm(2,j)=sqrt(rdum)
       enddo 
       rmax(1)=pca_std_cm(ipca)/(std_max(1)*sqrt(rmax(1)))

       do j=1,3 ! find max component for triad
          vdum=pca_a(ipca,j,:)
          rdum=dot_product(vdum,vdum)
          if (rdum > rmax(2)) rmax(2)=rdum
          angle_triad(1,j)=sqrt(rdum)
          vdum=pca_b(ipca,j,:)
          rdum=dot_product(vdum,vdum)
          if (rdum > rmax(2)) rmax(2)=rdum
          angle_triad(2,j)=sqrt(rdum)
       enddo 
       rmax(2)=pca_std(ipca)/(std_max(2)*sqrt(rmax(2)))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Determine cmh displacement amount and triad rotation angle/axis
       ! Axis for cm: get a perp vector from cm of C-domain to the pca_cm dir
       ! and also set center of rotation as this vector from cm_avg. 
       ! Note that the center of rotation for cm is not cm of C-domain.
       do k=1,2
          do j=1,3
             angle_cm(k,j)=angle_cm(k,j)*thet0*rmax(1)/pca_frm
             angle_triad(k,j)=angle_triad(k,j)*thet0*rmax(2)/pca_frm
             if (k==1) then
                axis_triad(k,j,:)=cross(triad_avg(k,j+1,:),pca_a(ipca,j,:))
             else
                axis_triad(k,j,:)=cross(triad_avg(k,j+1,:),pca_b(ipca,j,:))
             endif
             if (j>1) then ! find center of rotation and rot axis for cm
                if (k==1) then ! TCRa
                   vdum=pca_cma(ipca,j,:); call normalize(vdum)
                else  ! TCRb
                   vdum=pca_cmb(ipca,j,:); call normalize(vdum)
                endif                
                if (j==2) then
                   pdum=cm_avg(k,j,:)
                else 
                   pdum=triad_avg(k,1,:)
                endif
                vdum1=pdum-cm_avg(k,1,:) ! subtract C-domain pos
                ! dir_cm: vector from cm of C-domain to avg hinge/var cm.
                !         perp to pca_cm[ab]
                vdum1=vdum1-dot_product(vdum1,vdum)*vdum
                dir_cm(k,j-1,:)=vdum1
                ! center of rot: cm_ref_k,1/2,:): hinge/var
                cm_ref(k,j-1,:)=pdum-vdum1
                if (k==1) then ! TCRa
                   axis_cm(k,j-1,:)=cross(vdum1,pca_cma(ipca,j,:))
                else !TCRb
                   axis_cm(k,j-1,:)=cross(vdum1,pca_cmb(ipca,j,:))
                endif
             endif ! if (j>1) then
          enddo ! do j=1,3
       enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Write the motion
       do ifrm=-pca_frm,pca_frm ! total # of frames: 2*pca_frm+1
          write(12,'("MODEL")')
          write(14,'("MODEL")')
          do k=1,2
             p=1+3*(k-1)
             ! constant domain
             write(14,100) 'ATOM  ',p,' S  ',' CM',trim(chainID(k)),1, &
                  (cm_avg(k,1,idum),idum=1,3),1.0,0.5,segi(k)
             do j=1,2 ! rotate hinge & var
                rdum=ifrm*angle_cm(k,j+1) ! rotation angle
                pdum=rotate_perp(axis_cm(k,j,:),dir_cm(k,j,:),rdum) &
                     +cm_ref(k,j,:)
                p=j+3*(k-1)
                write(14,100) 'ATOM  ',p,' S  ',' CM',trim(chainID(k)),j, &
                     (pdum(idum),idum=1,3),1.0,0.5,segi(k)
             enddo

             ! rotate triad
             vdum=pdum ! store var domain cm position
             ! var domain cm
             p=1+4*(k-1)
             write(12,100) 'ATOM  ',p,' S  ','TRI',trim(chainID(k)),1, &
                  (pdum(idum),idum=1,3),1.0,0.5,segi(k)
             do j=1,3 ! 3 arms of triad
                rdum=ifrm*angle_triad(k,j) ! rotation angle
                vdum1=triad_avg(k,j+1,:)
                pdum=rotate_perp(axis_triad(k,j,:),vdum1,rdum)*tlen+vdum
                p=j+1+4*(k-1)
                write(12,100) 'ATOM  ',p,' S  ','TRI',trim(chainID(k)),j+1, &
                     (pdum(idum),idum=1,3),1.0,0.5,segi(k)
             enddo ! do j=1,3 ! 3 arms of triad
          enddo ! do k=1,2
          write(12,'("ENDMDL")')
          write(14,'("ENDMDL")')
       enddo ! ifrm=-pca_frm,pca_frm
       close(14)
       close(12)
    enddo ! ipca=1,3 ! create for PC1-PC3
    
  end subroutine write_pca_motion_cm

end module bocMod

  !************************************************************************
program pca_boc
  ! Principal component analysis of BOC
  !
  ! f95 pca_boc.f95 /home/hwm/computer/ftns/lapack/lapack-3.8.0/lib*.a matvec.f95 -o te
  ! 
  ! Usage: ./te -d dir_suffix -i input_prefx -o output_suffix -f0 frm0 -f1 frm1
  !
  use bocMod ! boc module defined after this program

  !tva,tvb,cha,chb,cma,cmb,

  integer idum,i,j
  character(len=128):: ifname,ofname,sdum,sdum1
  double precision rdum

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read command line argument
  call read_arg() !idir,iprefx)
  ifname = trim(idir)//"/boc_"//trim(iprefx)//".dat"
  print *,"# Input file: ", trim(ifname)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! read boc data and allocate memory
  dfrm=0 ! initialize
  call read_boc_data(ifname)
  print *,"# frm0, frm1, nframe= ", frm0, frm1, nframe
  call get_triad('a');  call get_triad('b')
  ofname=trim(idir)//"/pca"//trim(osuffix)//"/view_triad_traj.vmd"
  call write_triad_vmd(ofname)

  if (cflag) then
     call read_cm_data() ! call this after read_boc_data
     ofname=trim(idir)//"/pca"//trim(osuffix)//"/view_cm_traj.vmd"
     call write_cm_vmd(ofname)
  endif
  
  if (pdir=="") then 
     call pca_triad()
     if (cflag) call pca_cmh()
     call write_boc_avg()
     if (pca_frm > 0) then
        sdum1=trim(idir)//"/pca1/"
        sdum=trim(sdum1)//"dir_triad.dat"
        call read_pca_dir(pca_a,pca_b,sdum)
        sdum=trim(sdum1)//"dir_cm.dat"
        if (cflag) call read_pca_dir(pca_cma,pca_cmb,sdum)
        call read_boc_avg(sdum1,1) ! read to avg coord
        if (trim(idir) =="data_vab") then
           call write_pca_motion()
        else if (trim(idir)=="data_cab") then 
           call write_pca_motion_cm()
        else
           print *, "ERROR: Unknown idir: ",trim(idir)
           call exit(-1)
        endif
     endif
  else ! read existing PCA data to project OR generate movie 
     if (pca_frm == 0) then
        sdum=trim(pdir)//"/dir_triad.dat"
        call read_pca_dir(pca_a_ref,pca_b_ref,sdum)
        sdum=trim(pdir)//"/dir_cm.dat"
        call read_pca_dir(pca_cma_ref,pca_cmb_ref,sdum)
        call read_boc_avg(pdir,0) ! read to ref coord
        call project_boc(4)
     else ! generate movie of PCA motion
        sdum1=trim(idir)//"/pca1/"
        sdum=trim(sdum1)//"dir_triad.dat"
        call read_pca_dir(pca_a,pca_b,sdum)
        sdum=trim(sdum1)//"dir_cm.dat"
        call read_pca_dir(pca_cma,pca_cmb,sdum)
        call read_boc_avg(sdum1,1) ! read to avg coord
        if (trim(idir) =="data_vab") then
           call write_pca_motion()
        else if (trim(idir)=="data_cab") then 
           call write_pca_motion_cm()
        else
           print *, "ERROR: Unknown idir: ",trim(idir)
           call exit(-1)
        endif
     endif
  end if

end program pca_boc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine read_arg()
  use bocMod
  character(len=128) :: sdum, sdum1
  integer i, stat
  i=0
  idir=""; iprefx=""; osuffix=""; frm0=1; frm1=-1 ! default
  pdir=""
  ptheta=30. ! rot angle
  pca_frm=0 ! do not create PCA motion
  cflag=.true.
  
  do
     call get_command_argument(i, sdum)
     if (len_trim(sdum) == 0) exit
     if (trim(sdum) == '-d') then  ! input data directory
        i = i+1
        call get_command_argument(i, idir) 
     else if (trim(sdum) == '-i') then  ! input file prefix
        i = i+1
        call get_command_argument(i, iprefx) 
     else if (trim(sdum) == '-o') then  ! pca data output suffix
        i = i+1
        call get_command_argument(i, osuffix) ! initial frame to analyze
     else if (trim(sdum) == '-f0') then
        i = i+1
        call get_command_argument(i, sdum1) 
        read(sdum1,*,iostat=stat) frm0
     else if (trim(sdum) == '-f1') then
        i = i+1
        call get_command_argument(i, sdum1) 
        read(sdum1,*,iostat=stat) frm1
     else if (trim(sdum) == '-p') then
        i = i+1
        call get_command_argument(i, pdir) 
     else if (trim(sdum) == '-theta') then
        i = i+1
        call get_command_argument(i, sdum1)
        read(sdum1,*,iostat=stat) ptheta
     else if (trim(sdum) == '-fp') then
        i = i+1
        call get_command_argument(i, sdum1) 
        read(sdum1,*,iostat=stat) pca_frm
     else if (trim(sdum) == '-c') then
        i = i+1
        cflag=.false. ! do not process cmh

     endif
     i = i+1
     if (i>100) exit ! avoid infinite loop
  end do
  if ( idir == '') then
     print *,'ERROR: -d argument missing. Use "." for present dir.'
     stop
  endif
  if ( iprefx == '') then
     print *,'ERROR: -i argument missing'
     stop
  endif
  if ( osuffix == '') then
     print *,'ERROR: -o argument missing'
     stop
  endif
end subroutine read_arg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
