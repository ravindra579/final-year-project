! c ***********************************
module global
implicit none

integer::i_restart,i_step,Nmax,Nlg,Ng,lforce,iu1,id1,jb,jt,N1,&
N2,M1,M2,N,M,lspace,l_ABCN,l_CN
double precision::hc,ul,dl,radius,tol,betaxu,betaxd,betay,fr,&
rad_factor,pi,startingtime,Amp,St0,yfac,time,y_s,v_s,c_p,&
pseudo_fac,angle,dt,Re,Sts,rot,volume,coeffdm,pseudo_alpha,&
Fcx1,Fcy1,prx1,pry1,Str,Ampr,Ampr_dis,angler,r_dis,velocity_theta,&
pseudo_facx,v_sx,y_sx,pseudo_alphax,c_px,frx,Ampx,yfacx,anglex,Stsx
parameter(i_restart=1)
parameter(hc=23.d0,ul=8.d0,dl=25.d0,radius=0.5d0,tol=0.d0)
parameter(betaxu=1.06d0,betaxd=1.06d0,betay=1.06d0,fr=0.d0)
parameter(Re=200.d0,rad_factor=1.8d0,pi=4.d0*datan(1.d0))
parameter(N1=31,N2=194,M1=44,M2=207,startingtime=0.d0)
parameter(N=300,M=M1+M2-1,Ng=3,lspace=10)
parameter(Amp=0.d0,St0=0.d0)
parameter(Sts=fr*St0,dt=0.005d0,rot=0.d0)!rot=pi*fr*St0*0.4d0)
parameter(Str=0.d0,Ampr_dis=0.d0,Ampr=Ampr_dis*pi*Str)
parameter(l_ABCN=0,l_CN=1,yfac=Amp*(2.d0*pi*Sts)**2.d0,&
pseudo_alpha=2.d0*pi*Sts*Amp)
parameter(Ampx=0.d0,frx=0.d0,Stsx=frx*St0)
parameter(yfacx=Ampx*(2.d0*pi*Stsx)**2.d0,&
pseudo_alphax=2.d0*pi*Stsx*Ampx)
double precision,dimension(2,N+2,M+2)::x
double precision,dimension(2,N+1,M+1)::xcell,u,u_old
double precision,dimension(2,4,N,M)::s,diff
double precision,dimension(N+1,M+1)::vol_sol
double precision,dimension(N+1,M+1)::p,vol,ALW,ALS,ALP,AUE,AUN   
double precision,dimension(N+1,M+1)::app,awp,aep,asp,anp,apv,awv,aev&
,asv,anv
double precision,dimension(4,N+1,M+1)::F,F_old
double precision,dimension(2,N,M)::cvnp1,Fpr,del 
double precision,dimension(Ng,N,M)::w
double precision,dimension(N,M)::velx,vely,q
double precision,dimension(360)::anglefc,snor,dangle
integer,dimension(2,360)::ifl     
integer,dimension(2,Ng,N,M)::ingbr
integer,dimension(N+1,M+1)::icomput

      end module global
! c**************
program main
use global 
implicit none
integer::l,i,j
double precision::angle_p,angle_px,xb,yb
      
print*,'rotation',rot,Re
print*,'oscillation',Amp,fr
      
      call initialize
      
      do while(i_step.lt.100000000)
      
      time=time+dt
      i_step=i_step+1
      angle_p=angle      
      angle_px=anglex
      
      angle=2.d0*pi*Sts*time
      angler=2.d0*pi*Str*time
      anglex=2.d0*pi*Stsx*time
      y_s=Amp*dsin(angle)
      
      r_dis=-Ampr_dis*dcos(angler)
      velocity_theta=Ampr*dsin(angler)
      
      y_sx=Ampx*dsin(anglex)
      v_s=Amp*2.d0*pi*Sts*dcos(angle) 
      
      v_sx=Ampx*2.d0*pi*Stsx*dcos(anglex) 
      c_p=(pi/2.d0)*yfac*dsin(angle)
      
      c_px=(pi/2.d0)*yfacx*dsin(anglex)
      pseudo_fac=-yfac*(dsin(angle)+dsin(angle_p))
      
      pseudo_facx=-yfacx*(dsin(anglex)+dsin(angle_px))
      
      l=0
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.2)then
      l=l+1          
      xb=xcell(1,i,j)-del(1,i,j)
      yb=xcell(2,i,j)-del(2,i,j)      
      velx(i,j)=-2.d0*rot*yb-velocity_theta*dsin(anglefc(l))
	  vely(i,j)=2.d0*rot*xb+velocity_theta*dcos(anglefc(l))      
      endif
      enddo
      enddo   
      
        call momentum_init
      call previous

      if(l_CN.eq.1)then
      call CNLM
      elseif(l_CN.eq.2)then
      call CNNM
      endif
      if(l_ABCN.eq.1)call ABCN
      call unstd
      call intmd_result               
      
      enddo
      call output
       close(40);close(41);close(42);close(43)
      stop
      end
! c**************
subroutine initialize
      use global 
      implicit none
      open(40,file='cdcld.dat')
      open(41,file='cdclm.dat')
      open(42,file='velocityrms.dat')
      open(43,file='cfl.dat')
      call grid_unbounded
      call geo_para
      call cell_ident
      call mirror_ident
      call neighbors_new
      call forcecoeff_para
      call grid_plot
      call volume_mom
      call diagonals_SIP

      i_step=0
      if(i_restart.eq.1)then
      call init
      call vel_bc
      call p_bc
      call IBM_ubc
      call IBM_pbc
!       call flux_linear
      call mass_flux
      call startup_calc
      time=0.d0
      else
      call restart
      call vel_bc
      call p_bc
      call IBM_ubc
      call IBM_pbc
      call startup_calc
      time=startingtime      
      endif 
       
       angle=2.d0*pi*Sts*time
    
       anglex=2.d0*pi*Stsx*time
      return
      end
! c**********************
subroutine intmd_result
      use global 
      implicit none
        if(i_step.eq.i_step/1000*1000)then
        call save
        call final_data
        call position
        endif
 !       if(i_step.eq.110000)then
 !       call save1
  !      endif
      call force_coeff
      call momentum
      call cfl

      return
      end
! c**************
subroutine geo_para
      use global 
      implicit none
      integer::i,j
      double precision::xP,yP,xE,yE,xW,yW,xN,yN,xS,yS,&
     xne,yne,xnw,ynw,xse,yse,xsw,&
     ysw,xnt,ynt,xst,yst,xnb,ynb,xsb,ysb,&
     dx1,dx2,det,det1,det2,&
     a11,a12,a21,a22,a31,a32,alpha1,alpha2,&
     r1x,r1y,r2x,r2y,&
     uv11,uv12,uv21,uv22,dismin1,dis1,dismin2,dis2       
      
      do i=1,N+1
      do j=1,M+1
      xcell(1,i,j)=(x(1,i,j)+x(1,i+1,j)+x(1,i+1,j+1)+x(1,i,j+1))/4.d0
      xcell(2,i,j)=(x(2,i,j)+x(2,i+1,j)+x(2,i+1,j+1)+x(2,i,j+1))/4.d0
      enddo
      enddo

      do i=2,N
      do j=2,M

      xP=xcell(1,i,j)
      yP=xcell(2,i,j)

      xE=xcell(1,i+1,j)
      yE=xcell(2,i+1,j)
      
      xN=xcell(1,i,j+1)
      yN=xcell(2,i,j+1)

      xW=xcell(1,i-1,j)
      yW=xcell(2,i-1,j)

      xS=xcell(1,i,j-1)
      yS=xcell(2,i,j-1)

      xne=x(1,i+1,j+1)
      yne=x(2,i+1,j+1)
      xnw=x(1,i,j+1)
      ynw=x(2,i,j+1)
      xsw=x(1,i,j)
      ysw=x(2,i,j)
      xse=x(1,i+1,j)
      yse=x(2,i+1,j)
! c**east face
      dx1=dsqrt((xE-xP)**2.d0+(yE-yP)**2.d0)
      dx2=dsqrt((xne-xse)**2.d0+(yne-yse)**2.d0)
      uv11=(xE-xP)/dx1
      uv12=(yE-yP)/dx1
      uv21=(xne-xse)/dx2
      uv22=(yne-yse)/dx2      
      s(1,3,i,j)=yne-yse
      s(2,3,i,j)=xse-xne

      det=uv11*uv22-uv21*uv12
      det1=s(1,3,i,j)*uv22-s(2,3,i,j)*uv21
      det2=s(2,3,i,j)*uv11-s(1,3,i,j)*uv12

      alpha1=det1/det
      alpha2=det2/det

      diff(1,3,i,j)=alpha1/dx1
      diff(2,3,i,j)=alpha2/dx2
! c****north face
      dx1=dsqrt((xN-xP)**2.d0+(yN-yP)**2.d0)
      dx2=dsqrt((xne-xnw)**2.d0+(yne-ynw)**2.d0)
      uv11=(xN-xP)/dx1
      uv12=(yN-yP)/dx1
      uv21=(xne-xnw)/dx2
      uv22=(yne-ynw)/dx2      
      s(1,4,i,j)=ynw-yne
      s(2,4,i,j)=xne-xnw

      det=uv11*uv22-uv21*uv12
      det1=s(1,4,i,j)*uv22-s(2,4,i,j)*uv21
      det2=s(2,4,i,j)*uv11-s(1,4,i,j)*uv12

      alpha1=det1/det
      alpha2=det2/det

      diff(1,4,i,j)=alpha1/dx1
      diff(2,4,i,j)=alpha2/dx2
! c***west face
      dx1=dsqrt((xP-xW)**2.d0+(yP-yW)**2.d0)
      dx2=dsqrt((xnw-xsw)**2.d0+(ynw-ysw)**2.d0)
      uv11=(xP-xW)/dx1
      uv12=(yP-yW)/dx1
      uv21=(xnw-xsw)/dx2
      uv22=(ynw-ysw)/dx2      
      s(1,1,i,j)=ysw-ynw
      s(2,1,i,j)=xnw-xsw

      det=uv11*uv22-uv21*uv12
      det1=s(1,1,i,j)*uv22-s(2,1,i,j)*uv21
      det2=s(2,1,i,j)*uv11-s(1,1,i,j)*uv12

      alpha1=det1/det
      alpha2=det2/det

      diff(1,1,i,j)=alpha1/dx1
      diff(2,1,i,j)=alpha2/dx2
! c*****south face
      dx1=dsqrt((xP-xS)**2.d0+(yP-yS)**2.d0)
      dx2=dsqrt((xse-xsw)**2.d0+(yse-ysw)**2.d0)
      uv11=(xP-xS)/dx1
      uv12=(yP-yS)/dx1
      uv21=(xse-xsw)/dx2
      uv22=(yse-ysw)/dx2      
      s(1,2,i,j)=yse-ysw
      s(2,2,i,j)=xsw-xse

      det=uv11*uv22-uv21*uv12
      det1=s(1,2,i,j)*uv22-s(2,2,i,j)*uv21
      det2=s(2,2,i,j)*uv11-s(1,2,i,j)*uv12

      alpha1=det1/det
      alpha2=det2/det

      diff(1,2,i,j)=alpha1/dx1
      diff(2,2,i,j)=alpha2/dx2
      enddo
      enddo

      do i=1,N+1
      do j=1,M+1
      r1x=x(1,i+1,j)-x(1,i,j+1)
      r1y=x(2,i+1,j)-x(2,i,j+1)
      r2x=x(1,i+1,j+1)-x(1,i,j)
      r2y=x(2,i+1,j+1)-x(2,i,j)

      vol(i,j)=0.5d0*dabs(r1x*r2y-r1y*r2x)
      enddo
      enddo

! c****************************
      coeffdm=1.d0/Re

      dismin1=100.d0
      dismin2=100.d0
      do i=2,N
      dis1=dabs(xcell(1,i,2)-(-0.61d0))
      dis2=dabs(xcell(1,i,2)-(0.61d0))
      if(dis1.lt.dismin1)then
      dismin1=dis1
      iu1=i
      endif
      if(dis2.lt.dismin2)then
      dismin2=dis2
      id1=i
      endif
      enddo
      write(*,*)xcell(1,iu1,2),xcell(1,id1,2)

      dismin1=100.d0
      dismin2=100.d0
      do j=2,M
      dis1=dabs(xcell(2,2,j)-(-0.61d0))
      dis2=dabs(xcell(2,2,j)-(0.61d0))
      if(dis1.lt.dismin1)then
      dismin1=dis1
      jb=j
      endif
      if(dis2.lt.dismin2)then
      dismin2=dis2
      jt=j
      endif
      enddo
      write(*,*)xcell(2,2,jb),xcell(2,2,jt)

      return
      end
! c*********
subroutine grid_unbounded
      use global 
      implicit none
      integer::i,j,l
      double precision::upx,dnx,xu,yu,dx1,dx2,bx,by,&
     upy,dny,ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,&
     ay1,ay2,ay3,ay4,ay5,ay6,ay7,ay8,ay9      

      ax1=ul-rad_factor*radius
      ay1=(hc/2.d0)-rad_factor*radius
      ax2=ax1
      ay2=2.d0*rad_factor*radius
      ax3=ax1
      ay3=(hc/2.d0)-rad_factor*radius
      ax4=2.d0*rad_factor*radius
      ay4=ay1
      ax5=ax4
      ay5=ay2
      ax6=ax5
      ay6=ay3
      ax7=dl-rad_factor*radius
      ay7=ay4
      ax8=ax7
      ay8=ay5
      ax9=ax8
      ay9=ay6

! c******block 1**********
      
      dx1=1.d0/dfloat(N1-1)
      dx2=1.d0/dfloat(M1-1)
	  do i=1,N1
      do j=1,M1
      xu=dx1*dfloat(i-1)
      yu=dx2*dfloat(j-1)
      bx=(betaxu+1.d0)/(betaxu-1.d0)
      by=(betay+1.d0)/(betay-1.d0)
      upx=betaxu*(bx**xu-1.d0)
      dnx=1.d0+bx**xu
      upy=betay*(by**yu-1.d0)
      dny=1.d0+by**yu
      x(1,i+1,j+1)=-ul+ax1*(upx/dnx)
      x(2,i+1,j+1)=-(hc/2.d0)+ay1*(upy/dny)
      enddo
      enddo

! c*********block 2********

      
      dx1=1.d0/dfloat(N1-1)
      dx2=1.d0/dfloat(M2-M1)
	  do i=1,N1
      do j=M1,M2
      xu=dx1*dfloat(i-1)
      yu=dx2*dfloat(j-M1)
      bx=(betaxu+1.d0)/(betaxu-1.d0)
      upx=betaxu*(bx**xu-1.d0)
      dnx=1.d0+bx**xu
      x(1,i+1,j+1)=-ul+ax2*(upx/dnx)
      x(2,i+1,j+1)=x(2,2,M1+1)+ay2*yu
      enddo
      enddo

! c******block 3***************

      
      dx1=1.d0/dfloat(N1-1)
      dx2=1.d0/dfloat(M-M2)
	  do i=1,N1
      do j=M2,M
      xu=0.d0+dx1*dfloat(i-1)
      yu=0.d0+dx2*dfloat(j-M2)
      bx=(betaxu+1.d0)/(betaxu-1.d0)
      by=(betay+1.d0)/(betay-1.d0)
      upy=(betay+1.d0)-(betay-1.d0)*(by**(1.d0-yu))
      dny=1.d0+by**(1.d0-yu)
      upx=betaxu*(bx**xu-1.d0)
      dnx=1.d0+bx**xu
      x(1,i+1,j+1)=-ul+ax3*(upx/dnx)
      x(2,i+1,j+1)=x(2,2,M2+1)+ay3*(upy/dny)
      enddo
      enddo

! c******block 4**********
      
      
      dx1=1.d0/dfloat(N2-N1)
      dx2=1.d0/dfloat(M1-1)
	  do i=N1,N2
      do j=1,M1
      xu=dx1*dfloat(i-N1)
      yu=dx2*dfloat(j-1)
      by=(betay+1.d0)/(betay-1.d0)
      upy=betay*(by**yu-1.d0)
      dny=1.d0+by**yu
      x(1,i+1,j+1)=x(1,N1+1,2)+ax4*xu
      x(2,i+1,j+1)=-(hc/2.d0)+ay4*(upy/dny)
      enddo
      enddo

! c*********block 5********

      
      dx1=1.d0/dfloat(N2-N1)
      dx2=1.d0/dfloat(M2-M1)
	  do i=N1,N2
      do j=M1,M2
      xu=dx1*dfloat(i-N1)
      yu=dx2*dfloat(j-M1)
      x(1,i+1,j+1)=x(1,N1+1,2)+ax5*xu
      x(2,i+1,j+1)=x(2,2,M1+1)+ay5*yu
      enddo
      enddo

! c*****block 6****************

      
      dx1=1.d0/dfloat(N2-N1)
      dx2=1.d0/dfloat(M-M2)
	  do i=N1,N2
      do j=M2,M
      xu=dx1*dfloat(i-N1)
      yu=dx2*dfloat(j-M2)
      by=(betay+1.d0)/(betay-1.d0)
      upy=(betay+1.d0)-(betay-1.d0)*(by**(1.d0-yu))
      dny=1.d0+by**(1.d0-yu)
      x(1,i+1,j+1)=x(1,N1+1,2)+ax6*xu
      x(2,i+1,j+1)=x(2,2,M2+1)+ay6*(upy/dny)
      enddo
      enddo

! c*****block 7****************

      
      dx1=1.d0/dfloat(N-N2)
      dx2=1.d0/dfloat(M1-1)
	  do i=N2,N
      do j=1,M1
      xu=0.d0+dx1*dfloat(i-N2)
      yu=0.d0+dx2*dfloat(j-1)
      bx=(betaxd+1.d0)/(betaxd-1.d0)
      by=(betay+1.d0)/(betay-1.d0)
      upx=(betaxd+1.d0)-(betaxd-1.d0)*(bx**(1.d0-xu))
      dnx=1.d0+bx**(1.d0-xu)
      upy=betay*(by**yu-1.d0)
      dny=1.d0+by**yu
      x(1,i+1,j+1)=x(1,N2+1,2)+ax7*(upx/dnx)
      x(2,i+1,j+1)=-(hc/2.d0)+ay7*(upy/dny)
      enddo
      enddo

! c*****block 8******************

      
      dx1=1.d0/dfloat(N-N2)
      dx2=1.d0/dfloat(M2-M1)
	  do i=N2,N
      do j=M1,M2
      xu=0.d0+dx1*dfloat(i-N2)
      yu=0.d0+dx2*dfloat(j-M1)
      bx=(betaxd+1.d0)/(betaxd-1.d0)
      upx=(betaxd+1.d0)-(betaxd-1.d0)*(bx**(1.d0-xu))
      dnx=1.d0+bx**(1.d0-xu)
      x(1,i+1,j+1)=x(1,N2+1,2)+ax8*(upx/dnx)
      x(2,i+1,j+1)=x(2,N2+1,M1+1)+ay8*yu
      enddo
      enddo

! c*****block 9****************

      
      dx1=1.d0/dfloat(N-N2)
      dx2=1.d0/dfloat(M-M2)
	  do i=N2,N
      do j=M2,M
      xu=dx1*dfloat(i-N2)
      yu=dx2*dfloat(j-M2)
      bx=(betaxd+1.d0)/(betaxd-1.d0)
      upx=(betaxd+1.d0)-(betaxd-1.d0)*(bx**(1.d0-xu))
      dnx=1.d0+bx**(1.d0-xu)
      by=(betay+1.d0)/(betay-1.d0)
      upy=(betay+1.d0)-(betay-1.d0)*(by**(1.d0-yu))
      dny=1.d0+by**(1.d0-yu)
      x(1,i+1,j+1)=x(1,N2+1,2)+ax9*(upx/dnx)
      x(2,i+1,j+1)=x(2,N2+1,M2+1)+ay9*(upy/dny)
      enddo
      enddo

! c******************

      do j=2,M+1
      do l=1,2
      x(l,1,j)=x(l,2,j)
      x(l,N+2,j)=x(l,N+1,j)
      enddo
      enddo

      do i=1,N+2
      do l=1,2
      x(l,i,1)=x(l,i,2)
      x(l,i,M+2)=x(l,i,M+1)
      enddo
      enddo

      return
      end
! c**************
subroutine cell_ident
      use global 
      implicit none
      integer::i,j,l,lE,lW,lN,lS
      double precision::dis,dis1,delta
      
      do i=1,N+1
      do j=1,M+1
      dis=dsqrt(xcell(1,i,j)**2.d0+xcell(2,i,j)**2.d0)
      dis1=dabs(dis-radius)
      delta=dsqrt((x(1,i+1,j)-x(1,i,j))*(x(2,i,j+1)-x(2,i,j)))
      if((dis.gt.radius).and.(dis1.gt.tol*delta))then
!       if(dis.gt.radius)then
      icomput(i,j)=1
      else
      icomput(i,j)=0
      endif
      enddo
      enddo

! c*****************
      
      l=0
      do i=2,N
      do j=2,M
      if(icomput(i,j).ne.1)then
      lE=icomput(i+1,j)
      lW=icomput(i-1,j)
      lN=icomput(i,j+1)
      lS=icomput(i,j-1)
      if((lE.eq.1).or.(lW.eq.1).or.(lN.eq.1).or.(lS.eq.1))then
      icomput(i,j)=2
      l=l+1
      endif
      endif
      enddo
      enddo      
        Nlg=l
        write(*,*)"ghost cells",Nlg
      l=0
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then
      l=l+1
      endif
      enddo
      enddo
      Nmax=l
        write(*,*)"fluid cells",Nmax
      return
      end
! c***************
subroutine mirror_ident
      use global 
      implicit none
      integer::i,j
      double precision::xc,yc,xb,yb,dis
      
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.2)then
      xc=xcell(1,i,j)
      yc=xcell(2,i,j)
      
      dis=dsqrt(xc**2.d0+yc**2.d0)
      xb=xc*radius/dis
      yb=yc*radius/dis
      
      del(1,i,j)=xcell(1,i,j)-xb
      del(2,i,j)=xcell(2,i,j)-yb
 !     velx(i,j)=-2.d0*rot*yb
!	  vely(i,j)=2.d0*rot*xb      
      endif
      enddo
      enddo

      return
      end
! c***************
subroutine forcecoeff_para
      use global 
      implicit none
      integer::i,j,l
      double precision::xc,yc,ang,ang1,ang2,angmax,angmin,dismin1,dismin2,dis1,delta,dis 

      l=0
      angmax=-100.d0
      angmin=100.d0
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.2)then

      xc=xcell(1,i,j)
      yc=xcell(2,i,j)
      l=l+1
      snor(l)=dabs(dsqrt(xc**2.d0+yc**2.d0)-radius)
      ifl(1,l)=i
      ifl(2,l)=j
      
      if((xc.gt.0.d0).and.(yc.gt.0.d0))then
      anglefc(l)=datan(dabs(yc)/dabs(xc))
      elseif((xc.lt.0.d0).and.(yc.gt.0.d0))then
      anglefc(l)=pi-datan(dabs(yc)/dabs(xc))
      elseif((xc.lt.0.d0).and.(yc.lt.0.d0))then
      anglefc(l)=pi+datan(dabs(yc)/dabs(xc))
      elseif((xc.gt.0.d0).and.(yc.lt.0.d0))then
      anglefc(l)=2.d0*pi-datan(dabs(yc)/dabs(xc))
      endif
      
      if(anglefc(l).gt.angmax)angmax=anglefc(l)
      if(anglefc(l).lt.angmin)angmin=anglefc(l)
      
      endif
      enddo
      enddo
      lforce=l
     
      do l=1,lforce
      ang=anglefc(l)

      dismin1=100.d0
      dismin2=100.d0
      do i=1,lforce
      if(i.ne.l)then
      dis=dabs(ang-anglefc(i))

      if(ang.gt.anglefc(i))then 
      if(dis.lt.dismin1)then
      dismin1=dis
      ang1=anglefc(i)
      endif
      endif

      if(ang.lt.anglefc(i))then 
      if(dis.lt.dismin2)then
      dismin2=dis
      ang2=anglefc(i)
      endif
      endif
  
      endif
      enddo

      if(ang.gt.angmin)then
      ang1=ang1
      else
      ang1=angmax
      endif

      if(ang.lt.angmax)then
      ang2=ang2
      else
      ang2=angmin
      endif

      if((ang.gt.angmin).and.(ang.lt.angmax))then
      dangle(l)=(ang2-ang1)/2.d0
      else
      dangle(l)=pi+(ang2-ang1)/2.d0
      endif

      enddo      

      return
      end
! c***************
subroutine face_intpo(phi1,phi2,v1,v2,phi_int)
      use global
      implicit none
      double precision::phi1,phi2,v1,v2,phi_int
      
      phi_int=(phi1*v2+phi2*v1)/(v1+v2)

      return
      end
! c***********
subroutine conv_flux_HYB(phi,Fc,idd,jd,ld)
      use global 
      implicit none
      integer::idd,jd,ld
      double precision::Fc,cdse,cdsn,cdsw,cdss
      double precision,dimension(2,N+1,M+1)::phi

      call face_intpo(phi(ld,idd+1,jd),phi(ld,idd,jd),vol(idd+1,jd),vol(idd,jd),cdse)
      call face_intpo(phi(ld,idd,jd+1),phi(ld,idd,jd),vol(idd,jd+1),vol(idd,jd),cdsn)
      call face_intpo(phi(ld,idd-1,jd),phi(ld,idd,jd),vol(idd-1,jd),vol(idd,jd),cdsw)
      call face_intpo(phi(ld,idd,jd-1),phi(ld,idd,jd),vol(idd,jd-1),vol(idd,jd),cdss)

      Fc=cdse*F(3,idd,jd)+cdsn*F(4,idd,jd)+cdsw*F(1,idd,jd)+cdss*F(2,idd,jd)

      return
      end
! c***************
subroutine diff_flux(phi,Fd,idd,jd,ld)
      use global 
      implicit none
      integer::idd,jd,ld
      double precision::Fd
      double precision,dimension(2,N+1,M+1)::phi

      Fd=diff(1,3,idd,jd)*(phi(ld,idd+1,jd)-phi(ld,idd,jd))+&
     diff(1,4,idd,jd)*(phi(ld,idd,jd+1)-phi(ld,idd,jd))+&
     diff(1,1,idd,jd)*(phi(ld,idd,jd)-phi(ld,idd-1,jd))+&
     diff(1,2,idd,jd)*(phi(ld,idd,jd)-phi(ld,idd,jd-1))

      return
      end   
! c**********************
subroutine CNLM
      use global 
      implicit none
      integer::i,j,l
      double precision::rms,factor,fww,fss,fnn,aw,as,an,Fd,cv,df,pw,pe,pn,ps
      double precision,dimension(2,N,M)::ru,ru1
                            
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then
      do l=1,2
      call conv_flux_HYB(u,cv,i,j,l)
      call diff_flux(u,df,i,j,l)
      ru(l,i,j)=u_old(l,i,j)+(dt/(2.d0*vol(i,j)))*&
     (-cv+coeffdm*df)
     
      if(l.eq.1)then
	ru(l,i,j)=ru(l,i,j)-pseudo_facx*dt*0.5d0
	endif
	
        if(l.eq.2)then
	ru(l,i,j)=ru(l,i,j)-pseudo_fac*dt*0.5d0
	endif
              
              ru1(l,i,j)=ru(l,i,j)
              
      enddo
      endif
      enddo
      enddo
                     
      rms=1.d0
      do while(rms.gt.1.e-9)

      do i=2,N
      do j=2,M
      do l=1,4
      if(icomput(i,j).eq.1)then
      F_old(l,i,j)=F(l,i,j)
      endif
      enddo
      enddo
      enddo
      
! c**************************************
        do j=2,M
      factor=dt/(2.d0*vol(2,j))	
      fww=vol(2,j)/(vol(1,j)+vol(2,j))
      aw=factor*(F(1,2,j)*fww+coeffdm*diff(1,1,2,j))          
      ru1(1,2,j)=ru(1,2,j)-aw*(1.d0-v_sx)  
      ru1(2,2,j)=ru(2,2,j)+aw*v_s    
      enddo      
! c**************************************	
	do i=2,N
      factor=dt/(2.d0*vol(i,2))	
      fss=vol(i,2)/(vol(i,1)+vol(i,2))
      as=factor*(F(2,i,2)*fss+coeffdm*diff(1,2,i,2))
      	ru1(2,i,2)=ru(2,i,2)+as*v_s            
	  enddo      
! c***************************************	
	do i=2,N      
      factor=dt/(2.d0*vol(i,M))
      fnn=vol(i,M)/(vol(i,M+1)+vol(i,M))
      an=factor*(F(4,i,M)*fnn-coeffdm*diff(1,4,i,M))
      	ru1(2,i,M)=ru(2,i,M)+an*v_s            
	  enddo      
! c***************************************	        
!       do l=1,2
!       if((i.eq.2).and.(l.eq.1))then
!       factor=dt/(2.d0*vol(i,j))
! 	  fww=vol(i,j)/(vol(i-1,j)+vol(i,j))
!       aw=factor*(F(1,i,j)*fww+coeffdm*diff(1,1,i,j))
!       ru1(l,i,j)=ru(l,i,j)-aw
!       else
!       ru1(l,i,j)=ru(l,i,j)
!       endif
!       enddo
!       endif
!       enddo
!       enddo
      call BiCGSTAB(ru1)
      call flux_linear
!       call SIP
!       call BiCGSTABp
      call SIP_BiCGSTAB

      call mass_flux
      call flux_comp(rms)
!       write(*,*)"flux error",rms              
        print*,i_step,rms 
      enddo      
! c**************************************
        do j=2,M
      factor=dt/(2.d0*vol(2,j))	
      fww=vol(2,j)/(vol(1,j)+vol(2,j))
      aw=factor*(F(1,2,j)*fww+coeffdm*diff(1,1,2,j))          
      ru1(1,2,j)=ru(1,2,j)-aw*(1.d0-v_sx)      
      ru1(2,2,j)=ru(2,2,j)+aw*v_s    
      enddo      
! c**************************************	
	do i=2,N
      factor=dt/(2.d0*vol(i,2))	
      fss=vol(i,2)/(vol(i,1)+vol(i,2))
      as=factor*(F(2,i,2)*fss+coeffdm*diff(1,2,i,2))
      	ru1(2,i,2)=ru(2,i,2)+as*v_s            
	  enddo      
! c***************************************	
	do i=2,N      
      factor=dt/(2.d0*vol(i,M))
      fnn=vol(i,M)/(vol(i,M+1)+vol(i,M))
      an=factor*(F(4,i,M)*fnn-coeffdm*diff(1,4,i,M))
      	ru1(2,i,M)=ru(2,i,M)+an*v_s            
	  enddo      
! c***************************************	    

      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then
      do l=1,2

      if(l.eq.1)then
      call face_intpo(p(i-1,j),p(i,j),vol(i-1,j),vol(i,j),pw)
      call face_intpo(p(i+1,j),p(i,j),vol(i+1,j),vol(i,j),pe)
      Fd=pe*s(1,3,i,j)+pw*s(1,1,i,j)
     
      else
      call face_intpo(p(i,j-1),p(i,j),vol(i,j-1),vol(i,j),ps)
      call face_intpo(p(i,j+1),p(i,j),vol(i,j+1),vol(i,j),pn)
      Fd=pn*s(2,4,i,j)+ps*s(2,2,i,j)
      endif

      ru1(l,i,j)=ru1(l,i,j)-(dt/(2.d0*vol(i,j)))*(Fd+Fpr(l,i,j))
      Fpr(l,i,j)=Fd
      enddo
      endif
      enddo
      enddo

      call BiCGSTAB(ru1)
!       call IBM_pbc
!       call IBM_ubc

      return
      end
! c******************************
subroutine ABCN
      use global 
      implicit none
      integer::i,j,l
      double precision::cv,df,Fd,rms,pw,pe,pn,ps
      double precision,dimension(2,N,M)::ru
      
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then
      do l=1,2
      call conv_flux_HYB(u,cv,i,j,l)
      call diff_flux(u,df,i,j,l)
      ru(l,i,j)=u_old(l,i,j)+(dt/(2.d0*vol(i,j)))*&
     (-3.d0*cv+cvnp1(l,i,j)+coeffdm*df)
      cvnp1(l,i,j)=cv
      enddo
      endif
      enddo
      enddo
      
      call rhs_update(ru)
      call BiCGSTAB(ru)
      call flux_linear
!       call SIP
!       call BiCGSTABp
      call SIP_BiCGSTAB
      call mass_flux
      
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then
      do l=1,2

      if(l.eq.1)then
      call face_intpo(p(i-1,j),p(i,j),vol(i-1,j),vol(i,j),pw)
      call face_intpo(p(i+1,j),p(i,j),vol(i+1,j),vol(i,j),pe)
      Fd=pe*s(1,3,i,j)+pw*s(1,1,i,j)
      
      else
      call face_intpo(p(i,j-1),p(i,j),vol(i,j-1),vol(i,j),ps)
      call face_intpo(p(i,j+1),p(i,j),vol(i,j+1),vol(i,j),pn)
      Fd=pn*s(2,4,i,j)+ps*s(2,2,i,j)
      endif

      ru(l,i,j)=ru(l,i,j)-(dt/(2.d0*vol(i,j)))*(Fd+Fpr(l,i,j))
      Fpr(l,i,j)=Fd
      enddo
      endif
      enddo
      enddo

      call BiCGSTAB(ru)
!       call IBM_pbc           
      return
      end
! c***************
subroutine flux_linear
      use global 
      implicit none
      integer::i,j
      double precision::ue,vn,uw,vs
      
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then

      call face_intpo(u(1,i+1,j),u(1,i,j),vol(i+1,j),vol(i,j),ue)
      call face_intpo(u(2,i,j+1),u(2,i,j),vol(i,j+1),vol(i,j),vn)
      call face_intpo(u(1,i-1,j),u(1,i,j),vol(i-1,j),vol(i,j),uw)
      call face_intpo(u(2,i,j-1),u(2,i,j),vol(i,j-1),vol(i,j),vs)

      F(3,i,j)=ue*s(1,3,i,j)
      F(4,i,j)=vn*s(2,4,i,j)
      F(1,i,j)=uw*s(1,1,i,j)
      F(2,i,j)=vs*s(2,2,i,j)

      endif
      enddo
      enddo

      return
      end
! c**********************
subroutine mass_flux
      use global 
      implicit none
      integer::i,j
                  
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then

      F(3,i,j)=F(3,i,j)-dt*diff(1,3,i,j)*(p(i+1,j)-p(i,j))
      F(4,i,j)=F(4,i,j)-dt*diff(1,4,i,j)*(p(i,j+1)-p(i,j))
      F(1,i,j)=F(1,i,j)-dt*diff(1,1,i,j)*(p(i,j)-p(i-1,j))
      F(2,i,j)=F(2,i,j)-dt*diff(1,2,i,j)*(p(i,j)-p(i,j-1))

      endif
      enddo
      enddo

      return
      end
! c*****************      
subroutine BiCGSTAB(b)
      use global 
      implicit none
      integer::i,j,l,l1
      double precision::u_flux,apd,resid,Apc,Atc,rmsmax
      double precision,dimension(2,N,M)::b,rhs,r0,Ap1,At1
      double precision,dimension(2,N+1,M+1)::t1,p1
      double precision,dimension(2)::sum,sum1,sum2,beta,zi,alpha1,rms
      
      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      do l=1,2
      call A_operateu(u,u_flux,i,j,apd,l)
      if(icomput(i,j).eq.1)then
      rhs(l,i,j)=b(l,i,j)-u_flux
      elseif(icomput(i,j).eq.2)then
!       rhs(l,i,j)=-u_flux
        if(l.eq.1)rhs(l,i,j)=-velx(i,j)-u_flux
	if(l.eq.2)rhs(l,i,j)=-vely(i,j)-u_flux
      endif
      r0(l,i,j)=rhs(l,i,j)
      p1(l,i,j)=0.d0
      t1(l,i,j)=0.d0
      enddo
      endif
      enddo
      enddo

      l1=0
      do l=1,2
      rms(l)=1.d0
      beta(l)=0.d0
      zi(l)=0.d0
      enddo
      rmsmax=1.d0
      do while(rmsmax.gt.1.e-7)

      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      do l=1,2
      if(rms(l).gt.1.e-7)then
      p1(l,i,j)=rhs(l,i,j)+beta(l)*(p1(l,i,j)-zi(l)*Ap1(l,i,j))
      endif
      enddo
      endif
      enddo
      enddo

      do l=1,2
      sum(l)=0.d0
      sum1(l)=0.d0
      enddo
      do j=2,M
      do i=2,N
      if(icomput(i,j).ge.1)then
      do l=1,2
      if(rms(l).gt.1.e-7)then
      call A_operateu(p1,Apc,i,j,apd,l)
      sum(l)=sum(l)+rhs(l,i,j)*r0(l,i,j)
      sum1(l)=sum1(l)+Apc*r0(l,i,j)
      Ap1(l,i,j)=Apc
      endif
      enddo
      endif
      enddo
      enddo

      do l=1,2
      if(dabs(sum1(l)).gt.0.d0)then
      rms(l)=rms(l)
      else
      rms(l)=0.d0
      endif
      if(rms(l).gt.1.e-7)then
      alpha1(l)=sum(l)/sum1(l)
      endif
      enddo

      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      do l=1,2
      if(rms(l).gt.1.e-7)then
      t1(l,i,j)=rhs(l,i,j)-alpha1(l)*Ap1(l,i,j)
      endif
      enddo
      endif
      enddo
      enddo
      
      do l=1,2
      sum1(l)=0.d0
      sum2(l)=0.d0
      enddo
      do j=2,M
      do i=2,N
      if(icomput(i,j).ge.1)then
      do l=1,2
      if(rms(l).gt.1.e-7)then
      call A_operateu(t1,Atc,i,j,apd,l)
      sum1(l)=sum1(l)+Atc*t1(l,i,j)
      sum2(l)=sum2(l)+Atc*Atc
      At1(l,i,j)=Atc
      endif
      enddo
      endif
      enddo 
      enddo

      do l=1,2
      if(rms(l).gt.1.e-7)then
      zi(l)=sum1(l)/sum2(l)
      endif
      enddo

      do l=1,2
      sum1(l)=0.d0
      enddo
      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      do l=1,2
      if(rms(l).gt.1.e-7)then
      u(l,i,j)=u(l,i,j)+alpha1(l)*p1(l,i,j)+zi(l)*t1(l,i,j)
      rhs(l,i,j)=t1(l,i,j)-zi(l)*At1(l,i,j)
      sum1(l)=sum1(l)+r0(l,i,j)*rhs(l,i,j)
      endif
      enddo
      endif
      enddo
      enddo

      do l=1,2
      if(rms(l).gt.1.e-7)then
      beta(l)=(alpha1(l)*sum1(l))/(zi(l)*sum(l))
      endif
      enddo

      do l=1,2
      sum(l)=0.d0
      enddo
      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      do l=1,2
      if(rms(l).gt.1.e-7)then
      call A_operateu(u,u_flux,i,j,apd,l)
      if(icomput(i,j).eq.1)then
      resid=b(l,i,j)-u_flux
      elseif(icomput(i,j).eq.2)then
!       resid=-u_flux
        if(l.eq.1)resid=-velx(i,j)-u_flux
	if(l.eq.2)resid=-vely(i,j)-u_flux
      endif
      sum(l)=sum(l)+resid**2.d0
      endif
      enddo
      endif
      enddo
      enddo
 
      do l=1,2
      if(rms(l).gt.1.e-7)then
      rms(l)=dsqrt(sum(l)/dfloat(Nmax+Nlg))
      endif
      enddo
      rmsmax=dmax1(rms(1),rms(2))
        l1=l1+1
      if(l1.eq.l1/200*200)write(*,*)l1,rmsmax      

      enddo
!       write(*,*)l1,rmsmax

      call vel_bc
      call IBM_ubc

      return
      end
! c********************     
subroutine SIP_BiCGSTAB
      use global 
      implicit none
      integer::i,j,l1
      double precision::u_flux,apd,resid,Apc,Atc,p_flux,beta,rho1,rho1_old,sum1,alpha1,zi,rms,Ay,Az,residual,sum2
      double precision,dimension(N,M)::rhs,r0
      double precision,dimension(N+1,M+1)::t1,s1,p1,y1,z1,v1
            
      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      call A_operatep(p,p_flux,apd,i,j)
      if(icomput(i,j).eq.1)then
      rhs(i,j)=(F(1,i,j)+F(2,i,j)+F(3,i,j)+F(4,i,j))/dt-p_flux
      elseif(icomput(i,j).eq.2)then
      rhs(i,j)=-p_flux
      endif
      r0(i,j)=rhs(i,j)
      v1(i,j)=0.d0
      p1(i,j)=0.d0
      endif
      enddo
      enddo

      l1=0
      rms=1.d0
      rho1_old=1.d0
      alpha1=1.d0
      zi=1.d0
      do while(rms.gt.1.e-7)

      sum1=0.d0
      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      sum1=sum1+rhs(i,j)*r0(i,j)
      endif
      enddo
      enddo
      rho1=sum1

      beta=(rho1*alpha1)/(rho1_old*zi)

      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      p1(i,j)=rhs(i,j)+beta*(p1(i,j)-zi*v1(i,j))
      y1(i,j)=(p1(i,j)-ALS(i,j)*y1(i,j-1)-ALW(i,j)*y1(i-1,j))/ALP(i,j)
      endif
      enddo
      enddo

      do j=M,2,-1
      do i=N,2,-1
      if(icomput(i,j).ge.1)then
      y1(i,j)=y1(i,j)-AUN(i,j)*y1(i,j+1)-AUE(i,j)*y1(i+1,j)
      endif
      enddo
      enddo

      sum1=0.d0
      do j=2,M
      do i=2,N
      if(icomput(i,j).ge.1)then
      call A_operatep(y1,Ay,apd,i,j)
      sum1=sum1+Ay*r0(i,j)
      v1(i,j)=Ay
      endif
      enddo
      enddo
      alpha1=rho1/sum1

      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      s1(i,j)=rhs(i,j)-alpha1*v1(i,j)
      z1(i,j)=(s1(i,j)-ALS(i,j)*z1(i,j-1)-ALW(i,j)*z1(i-1,j))/ALP(i,j)
      endif
      enddo
      enddo
   
      do j=M,2,-1
      do i=N,2,-1
      if(icomput(i,j).ge.1)then
      z1(i,j)=z1(i,j)-AUN(i,j)*z1(i,j+1)-AUE(i,j)*z1(i+1,j)
      endif
      enddo
      enddo

      sum1=0.d0
      sum2=0.d0
      do j=2,M
      do i=2,N
      if(icomput(i,j).ge.1)then
      call A_operatep(z1,Az,apd,i,j)
      sum1=sum1+Az*s1(i,j)
      sum2=sum2+Az*Az
      t1(i,j)=Az
      endif
      enddo
      enddo
      zi=sum1/sum2

      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      p(i,j)=p(i,j)+alpha1*y1(i,j)+zi*z1(i,j)
      rhs(i,j)=s1(i,j)-zi*t1(i,j)
      endif
      enddo
      enddo
      rho1_old=rho1

      residual=0.d0
      do i=2,N
      do j=2,M
      if(icomput(i,j).ge.1)then
      call A_operatep(p,p_flux,apd,i,j)
      if(icomput(i,j).eq.1)then
      resid=(F(1,i,j)+F(2,i,j)+F(3,i,j)+F(4,i,j))/dt-p_flux
      elseif(icomput(i,j).eq.2)then
      resid=-p_flux
      endif
      residual=residual+resid**2.d0
      endif
      enddo
      enddo
      rms=dsqrt(residual/dfloat(Nmax+Nlg))
      l1=l1+1
!       if(l1.eq.l1/100*100)write(*,*)l1,rms      
!        write(*,*)l1,rms 
      enddo
!       write(*,*)l1
        
      call p_bc
      call IBM_pbc      
!         stop 
      return
      end      
! c***************
subroutine vel_bc
      use global 
      implicit none
      integer::i,j,l
      do i=2,N
      u(1,i,1)=u(1,i,2)
      u(2,i,1)=-v_s
      u(1,i,M+1)=u(1,i,M)
      u(2,i,M+1)=-v_s
      enddo

      do j=2,M
      do l=1,2
      if(l.eq.1)then
      u(l,1,j)=1.d0-v_sx
      else
      u(l,1,j)=-v_s
      endif
      u(l,N+1,j)=u(l,N,j)
      enddo
      enddo

      return
      end
! c*****************
subroutine p_bc
      use global 
      implicit none
      integer::i,j
      do i=2,N
      p(i,1)=p(i,2)
      p(i,M+1)=p(i,M)
      enddo

      do j=2,M
      p(1,j)=p(2,j)
      p(N+1,j)=0.d0
      enddo

      return
      end
! c****************      
subroutine previous
      use global 
      implicit none
      integer::i,j,l
      do i=1,N+1
      do j=1,M+1
      do l=1,2
      u_old(l,i,j)=u(l,i,j)
      enddo
      enddo
      enddo

      return
      end
! c**************
subroutine unstd
      use global 
      implicit none
      integer::i,j
      double precision::sumu,sumv,rmsu,rmsv,rms 
      
      sumu=0.d0
      sumv=0.d0
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then
      sumu=sumu+(u(1,i,j)-u_old(1,i,j))**2.d0
      sumv=sumv+(u(2,i,j)-u_old(2,i,j))**2.d0
      endif
      enddo
      enddo
      rmsu=dsqrt(sumu/dfloat(Nmax))
      rmsv=dsqrt(sumv/dfloat(Nmax))
      rms=dmax1(rmsu,rmsv)
      write(42,*)time,rmsu,rmsv
      if(dmax1(rmsu,rmsv).gt.1.d0)then
      print*,'problem1'
      stop
      endif     

      return
      end
! c**************
subroutine flux_comp(rmsd)
      use global 
      implicit none
      integer::i,j 
      double precision::rmsw,rmss,rmse,rmsn,rmsd,sumw,sume,sumn,sums
      sume=0.d0
      sumw=0.d0
      sumn=0.d0
      sums=0.d0
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then
      sumw=sumw+(F(1,i,j)-F_old(1,i,j))**2.d0
      sums=sums+(F(2,i,j)-F_old(2,i,j))**2.d0
      sume=sume+(F(3,i,j)-F_old(3,i,j))**2.d0
      sumn=sumn+(F(4,i,j)-F_old(4,i,j))**2.d0
      endif
      enddo
      enddo
      rmsw=dsqrt(sumw/dfloat(Nmax))
      rmss=dsqrt(sums/dfloat(Nmax))
      rmse=dsqrt(sume/dfloat(Nmax))
      rmsn=dsqrt(sumn/dfloat(Nmax))
      rmsd=dmax1(rmsw,rmss,rmse,rmsn)

      return
      end
! c**************
subroutine init
      use global 
      implicit none
      integer::i,j
      do i=1,N+1
      do j=1,M+1
      u(1,i,j)=0.d0
      u(2,i,j)=0.d0
      p(i,j)=0.d0 
!       do l=1,4
!       F(l,i,j)=0.d0
!       enddo
      enddo
      enddo

      return
      end
! c***************
subroutine save
      use global 
      implicit none
        integer::i,j,l
      open(21,file='save_var.dat') 
      do i=1,N+1
      do j=1,M+1
      write(21,100)u(1,i,j),u(2,i,j),p(i,j)
      enddo
      enddo
      close(21)

      open(22,file='save_fl.dat') 
      do i=2,N
      do j=2,M
      write(22,200)(F(l,i,j),l=1,4)
      enddo
      enddo
      close(22)

100   format(3(f23.16)) 
200   format(4(f23.16))             

      return
      end
! c****************
subroutine save1
      use global 
      implicit none
        integer::i,j,l
      open(25,file='save_var_300.dat') 
      do i=1,N+1
      do j=1,M+1
      write(25,100)u(1,i,j),u(2,i,j),p(i,j)
      enddo
      enddo
      close(25)

      open(26,file='save_fl_300.dat') 
      do i=2,N
      do j=2,M
      write(26,200)(F(l,i,j),l=1,4)
      enddo
      enddo
      close(26)

100   format(3(f23.16)) 
200   format(4(f23.16))             

      return
      end
! c****************
subroutine restart
      use global 
      implicit none
        integer::i,j,l
      open(23,file='save_var.dat') 
      do i=1,N+1
      do j=1,M+1
      read(23,100)u(1,i,j),u(2,i,j),p(i,j)
      enddo
      enddo
      close(23)

      open(24,file='save_fl.dat') 
      do i=2,N
      do j=2,M
      read(24,200)(F(l,i,j),l=1,4)
      enddo
      enddo
      close(24)

100   format(3(f23.16)) 
200   format(4(f23.16))             

      return
      end     
! c****************
subroutine output
      use global 
      implicit none
        integer::i,j,l
        double precision::si 
      open(1,file='res.dat')
      write(1,100)N-1,M-1
      do j=2,M
      do i=2,N
      write(1,200)(xcell(l,i,j),l=1,2),(u(l,i,j),l=1,2),p(i,j)
      enddo
      enddo
      close(1)

      open(2,file='si.dat')
      write(1,100)N-1,M-1
      do i=2,N
      si=0.d0
      do j=2,M
      if(icomput(i,j).eq.1)then
      si=si+(u(1,i,j)+u(1,i,j-1))*(xcell(2,i,j)-xcell(2,i,j-1))/2.d0
      write(2,*)(xcell(l,i,j),l=1,2),si
      else
      write(2,*)(xcell(l,i,j),l=1,2),100.d0
      endif
      enddo
      enddo
      close(2)

100   format("ZONE",1x,"I=",i3,1x,"J=",i3)
200   format(5(f16.8,1x))
      
      return
      end
! c*******************
subroutine grid_plot
      use global 
      implicit none
        integer::i,j,l
      open(1,file='grid.dat')
      write(1,100)N+2,M+2
      do j=1,M+2
      do i=1,N+2
      write(1,*)(x(l,i,j),l=1,2)
      enddo
      enddo
      close(1)

      open(11,file='fluid.dat')
      do i=1,N+1
      do j=1,M+1
      if(icomput(i,j).eq.1)then
      write(11,*)xcell(1,i,j),xcell(2,i,j)
      endif
      enddo
      enddo
      close(11)

      open(12,file='ghost.dat')
      do i=1,N+1
      do j=1,M+1
      if(icomput(i,j).eq.2)then
      write(12,*)xcell(1,i,j),xcell(2,i,j)
      endif
      enddo
      enddo
      close(12)

!       open(13,file='bndry.dat')
!       do i=1,N+1
!       do j=1,M+1
!       if(icomput(i,j).eq.2)then
!       write(13,*)ghdist(1,Ng,i,j),ghdist(2,Ng,i,j)
!       endif
!       enddo
!       enddo
!       close(13)
      
      open(14,file='solid.dat')
      do i=1,N+1
      do j=1,M+1
      if(icomput(i,j).eq.0)then
      write(14,*)xcell(1,i,j),xcell(2,i,j)
      endif
      enddo
      enddo
      close(14)

100   format("ZONE",1x,"I=",i3,1x,"J=",i3)

      return
      end
! c******************
subroutine diagonals_SIP
      use global 
      implicit none
      integer::i,j
      double precision::alpha,ap,aw,ae,as,an,factor
      alpha=0.9d0

      do i=1,N+1
      do j=1,M+1
      ALP(i,j)=0.d0
      ALW(i,j)=0.d0
      ALS(i,j)=0.d0
      AUE(i,j)=0.d0
      AUN(i,j)=0.d0
      enddo
      enddo

      do j=2,M
      do i=2,N
      if(icomput(i,j).ge.1)then

      ae=diff(1,3,i,j)
      aw=-diff(1,1,i,j)
      an=diff(1,4,i,j)
      as=-diff(1,2,i,j)
      ap=-(ae+aw+an+as)
      
      if(i.eq.2)then
      ap=ap+aw
      aw=0.d0
      endif

      if(j.eq.2)then
      ap=ap+as
      as=0.d0
      endif

      if(i.eq.N)then
      ae=0.d0
      endif

      if(j.eq.M)then
      ap=ap+an
      an=0.d0
      endif

      ALS(i,j)=as/(1.d0+alpha*AUE(i,j-1))

      ALW(i,j)=aw/(1.d0+alpha*AUN(i-1,j))

      ALP(i,j)=ap+alpha*(ALW(i,j)*AUN(i-1,j)+ALS(i,j)*AUE(i,j-1))-&
     ALW(i,j)*AUE(i-1,j)-ALS(i,j)*AUN(i,j-1)   

      AUE(i,j)=(ae-alpha*ALS(i,j)*AUE(i,j-1))/ALP(i,j)

      AUN(i,j)=(an-alpha*ALW(i,j)*AUN(i-1,j))/ALP(i,j)

      awp(i,j)=aw
      aep(i,j)=ae
      app(i,j)=ap
      anp(i,j)=an
      asp(i,j)=as

      factor=dt/(2.d0*vol(i,j))
      awv(i,j)=factor*coeffdm*diff(1,1,i,j)
	  aev(i,j)=-factor*coeffdm*diff(1,3,i,j)
	  anv(i,j)=-factor*coeffdm*diff(1,4,i,j)
	  asv(i,j)=factor*coeffdm*diff(1,2,i,j)
	  apv(i,j)=1.d0-factor*coeffdm*(diff(1,1,i,j)+diff(1,2,i,j)-&
     diff(1,3,i,j)-diff(1,4,i,j))

      endif
      enddo
      enddo

      return
      end
! c*****************
subroutine startup_calc
      use global 
      implicit none
      integer::i,j
      double precision::pw,ps,pe,pn
      
      do i=2,N
      do j=2,M
      if(icomput(i,j).eq.1)then
      call face_intpo(p(i-1,j),p(i,j),vol(i-1,j),vol(i,j),pw)
      call face_intpo(p(i,j-1),p(i,j),vol(i,j-1),vol(i,j),ps)
      call face_intpo(p(i+1,j),p(i,j),vol(i+1,j),vol(i,j),pe)
      call face_intpo(p(i,j+1),p(i,j),vol(i,j+1),vol(i,j),pn)        

      Fpr(1,i,j)=pe*s(1,3,i,j)+pw*s(1,1,i,j)
      Fpr(2,i,j)=pn*s(2,4,i,j)+ps*s(2,2,i,j)

      endif
      enddo
      enddo

      return
      end
! c****************************
subroutine rhs_update(b)
      use global 
      implicit none
      integer::i,j
      double precision::aw,factor,fww
      double precision,dimension(2,N,M)::b

      do j=2,M
      factor=dt/(2.d0*vol(2,j))
      if(l_CN.eq.2)then
	  fww=vol(2,j)/(vol(1,j)+vol(2,j))
      aw=factor*(F(1,2,j)*fww+coeffdm*diff(1,1,2,j))
      elseif(l_ABCN.eq.1)then
      aw=factor*coeffdm*diff(1,1,2,j)
      endif
      b(1,2,j)=b(1,2,j)-aw
      enddo
      
      return
      end
! c**********************      
subroutine force_coeff
      use global 
      implicit none
      integer::l
      double precision::sum1,sum2,um,vm,prf,cd,cl,surf,shf,c_v,e_v,pseudo,gfunc,gp,gfuncx,gpx,pseudox
      
      sum1=0.d0
      sum2=0.d0
      do l=1,lforce
      surf=radius*dangle(l)
      prf=p(ifl(1,l),ifl(2,l))*surf
!       prf=pbnd(l)*surf
      um=-u(1,ifl(1,l),ifl(2,l))
      vm=-u(2,ifl(1,l),ifl(2,l))

      shf=surf*(um*dsin(anglefc(l))-vm*dcos(anglefc(l)))/snor(l)

      sum1=sum1+(-prf*dcos(anglefc(l))+(1.d0/Re)*shf*dsin(anglefc(l)))
      sum2=sum2+(-prf*dsin(anglefc(l))-(1.d0/Re)*shf*dcos(anglefc(l)))

      enddo
      cd=2.d0*sum1
      cl=2.d0*sum2
!       pseudo=-2.d0*yfac*dsin(angle)
      gp=pseudo_alpha*dcos(angle)
      gfunc=gp**2.d0/(1.d0+gp**2.d0)
      
      gpx=pseudo_alphax*dcos(anglex)
      gfuncx=gpx**2.d0/(1.d0+gpx**2.d0)
!       print*,i_step,gp,gfunc
        pseudo=-c_p*gfunc
        
        pseudox=-c_px*gfuncx
        cl=cl-pseudo
        
        cd=cd-pseudox
      c_v=cl-c_p
      e_v=c_v*v_s
      write(40,200)time,cd,cl,c_p,c_v,e_v,y_s,r_dis
200   format(8(f14.8,1x))     

      return
      end
! c******************************************
subroutine momentum
      use global 
      implicit none
      integer::i,j,l
      double precision::unstdyx,unstdyy,xmomd,xmomu,xmomb,xmomt,&
      ymomd,ymomu,ymomb,ymomt,ue,uw,ve,vw,pfd,pfu,pft,pfb,prx,pry,&
      Fcx,Fcy,cd,cl,c_v,e_v,un,vn,us,vs,pe,pw,pn,ps
      
      unstdyx=0.d0
      unstdyy=0.d0
      do i=iu1,id1
      do j=jb,jt
      if(icomput(i,j).eq.1)then
      unstdyx=unstdyx+vol(i,j)*(u(1,i,j)-u_old(1,i,j))/dt
      unstdyy=unstdyy+vol(i,j)*(u(2,i,j)-u_old(2,i,j))/dt
      endif
      enddo
      enddo

      xmomt=0.d0
      xmomb=0.d0
      xmomu=0.d0
      xmomd=0.d0
      ymomt=0.d0
      ymomb=0.d0
      ymomu=0.d0
      ymomd=0.d0
      pfd=0.d0
      pfu=0.d0
      pft=0.d0
      pfb=0.d0
      cl=0.d0;cd=0.d0
      do j=jb,jt
      call face_intpo(u(1,id1+1,j),u(1,id1,j),vol(id1+1,j),&
      vol(id1,j),ue)
      call face_intpo(u(2,id1+1,j),u(2,id1,j),vol(id1+1,j),&
      vol(id1,j),ve)
      call face_intpo(u(1,iu1-1,j),u(1,iu1,j),vol(iu1-1,j),&
      vol(iu1,j),uw)
      call face_intpo(u(2,iu1-1,j),u(2,iu1,j),vol(iu1-1,j),&
      vol(iu1,j),vw)
      call face_intpo(p(id1+1,j),p(id1,j),vol(id1+1,j),&
      vol(id1,j),pe)
      call face_intpo(p(iu1-1,j),p(iu1,j),vol(iu1-1,j),&
      vol(iu1,j),pw)

      xmomd=xmomd+F(3,id1,j)*ue
      xmomu=xmomu+F(1,iu1,j)*uw
      ymomd=ymomd+F(3,id1,j)*ve
      ymomu=ymomu+F(1,iu1,j)*vw
      pfd=pfd+pe*s(1,3,id1,j)      
      pfu=pfu+pw*s(1,1,iu1,j)
      enddo

      do i=iu1,id1
      call face_intpo(u(1,i,jt+1),u(1,i,jt),vol(i,jt+1),vol(i,jt),un)
      call face_intpo(u(2,i,jt+1),u(2,i,jt),vol(i,jt+1),vol(i,jt),vn)
      call face_intpo(u(1,i,jb-1),u(1,i,jb),vol(i,jb-1),vol(i,jb),us)
      call face_intpo(u(2,i,jb-1),u(2,i,jb),vol(i,jb-1),vol(i,jb),vs)
      call face_intpo(p(i,jt+1),p(i,jt),vol(i,jt+1),vol(i,jt),pn)
      call face_intpo(p(i,jb-1),p(i,jb),vol(i,jb-1),vol(i,jb),ps)

      xmomt=xmomt+F(4,i,jt)*un
      xmomb=xmomb+F(2,i,jb)*us
      ymomt=ymomt+F(4,i,jt)*vn
      ymomb=ymomb+F(2,i,jb)*vs
      pft=pft+pn*s(2,4,i,jt)
      pfb=pfb+ps*s(2,2,i,jb)
      enddo

      Fcx=xmomd+xmomu+xmomt+xmomb
      Fcy=ymomt+ymomb+ymomd+ymomu
      prx=pfu+pfd
      pry=pfb+pft
      cd=(prx+prx1)+(Fcx+Fcx1)+2.d0*unstdyx
      cl=(pry+pry1)+(Fcy+Fcy1)+2.d0*unstdyy
      Fcx1=Fcx
      Fcy1=Fcy
      prx1=prx
      pry1=pry     
      
      cl=cl+pseudo_fac*volume
      
      cd=cd+pseudo_facx*volume
      c_v=-cl-c_p
      e_v=c_v*v_s

      write(41,200)time,-cd,-cl,c_p,c_v,e_v,y_s,r_dis
200   format(8(f14.8,1x))      
      
      return
      end
! c******************************************          
subroutine fdf(idd,jd,phid,ld,dxp,dxn,dyp,dyn,fx,fy)
        use global   
        implicit none
        integer::idd,jd,ld
        double precision::dxp,dxn,dyp,dyn,fx,fy
        double precision,dimension(2,N+1,M+1)::phid
            
      fx=(phid(ld,idd+1,jd)*dxn**2.d0-phid(ld,idd-1,jd)*dxp**2.d0-&
     phid(ld,idd,jd)*(dxn**2.d0-dxp**2.d0))/(dxp*dxn*(dxp+dxn))
      
      fy=(phid(ld,idd,jd+1)*dyn**2.d0-phid(ld,idd,jd-1)*dyp**2.d0-&
     phid(ld,idd,jd)*(dyn**2.d0-dyp**2.d0))/(dyp*dyn*(dyp+dyn))
                  
      return
      end
! ***************************
subroutine fds(idd,jd,phid,ld,dxp,dxn,dyp,dyn,sx,sy)
        use global 
        implicit none
      integer::idd,jd,ld
        double precision::dxp,dxn,dyp,dyn,sx,sy
        double precision,dimension(2,N+1,M+1)::phid
     
      sx=(phid(ld,idd+1,jd)*dxn+phid(ld,idd-1,jd)*dxp-&
     phid(ld,idd,jd)*(dxp+dxn))/(0.5d0*dxp*dxn*(dxp+dxn))

      sy=(phid(ld,idd,jd+1)*dyn+phid(ld,idd,jd-1)*dyp-&
     phid(ld,idd,jd)*(dyp+dyn))/(0.5d0*dyp*dyn*(dyp+dyn))
            
      return
      end
! ***************************
subroutine fdm(idd,jd,phid,ld,dxp,dxn,dyp,dyn,amxy)
        use global 
        implicit none
        integer::idd,jd,ld
        double precision::dxp,dxn,dyp,dyn,amxy
        double precision,dimension(2,N+1,M+1)::phid
            
      amxy=(phid(ld,idd+1,jd+1)-phid(ld,idd+1,jd-1)-phid(ld,idd-1,jd+1)+&
     phid(ld,idd-1,jd-1))/(dxp*dyp+dxp*dyn+dxn*dyp+dxn*dyn)
      
      return
      end
! ***************************
subroutine coefficients(idd,jd,dxp,dxn,dyp,dyn)
        use global 
        implicit none
        integer::idd,jd
        double precision::dxp,dxn,dyp,dyn
              
      dxp=xcell(1,idd+1,jd)-xcell(1,idd,jd)
      dxn=xcell(1,idd,jd)-xcell(1,idd-1,jd)

      dyp=xcell(2,idd,jd+1)-xcell(2,idd,jd)
      dyn=xcell(2,idd,jd)-xcell(2,idd,jd-1)

      return
      end
! ***************************************
subroutine coeff_p(idd,jd,dxp,dxn)
        use global 
        implicit none
        integer::idd,jd
        double precision::dxp,dxn        

      dxp=xcell(1,idd+1,jd)-xcell(1,idd,jd)
      dxn=xcell(1,idd,jd)-xcell(1,idd-1,jd)      
      
      return
      end
! ***************************************
subroutine fdf_p(idd,jd,phid,dxp,dxn,fx)
        use global      
        implicit none
        integer::idd,jd
        double precision::dxp,dxn,fx        
       double precision,dimension(N+1,M+1)::phid
            
      fx=(phid(idd+1,jd)*dxn**2.d0-phid(idd-1,jd)*dxp**2.d0-&
     phid(idd,jd)*(dxn**2.d0-dxp**2.d0))/(dxp*dxn*(dxp+dxn))          
                 
      return
      end
! ***************************
subroutine fds_p(idd,jd,phid,dxp,dxn,sx)
        use global   
        implicit none
       integer::idd,jd
        double precision::dxp,dxn,sx        
       double precision,dimension(N+1,M+1)::phid
     
      sx=(phid(idd+1,jd)*dxn+phid(idd-1,jd)*dxp-&
     phid(idd,jd)*(dxp+dxn))/(0.5d0*dxp*dxn*(dxp+dxn))
      
      return
      end
! ***************************           
subroutine IBM_pbc
        use global 
        implicit none
      integer::i,j,l,i1,j1,i2,j2,i3,j3
      double precision::Aphi,ap,first,second,amixed,&
      dxp1,dxn1,fx1,fy1,sx1,sy1,amxy1,&
      dxp2,dxn2,fx2,fy2,sx2,sy2,amxy2,&
      dxp3,dxn3,fx3,fy3,sx3,sy3,amxy3,xb,yb
      
      do i=2,N
      do j=2,M

      if(icomput(i,j).eq.2)then
      
      i1=ingbr(1,1,i,j);i2=ingbr(1,2,i,j);i3=ingbr(1,3,i,j)
!       i4=ingbr(1,4,i,j)
      j1=ingbr(2,1,i,j);j2=ingbr(2,2,i,j);j3=ingbr(2,3,i,j)
!       j4=ingbr(2,4,i,j)      
      
      xb=xcell(1,i,j)-del(1,i,j)
      yb=xcell(2,i,j)-del(2,i,j)
      
      call coeff_p(i1,j1,dxp1,dxn1)
      call coeff_p(i2,j2,dxp2,dxn2)
      call coeff_p(i3,j3,dxp3,dxn3)
!       call coeff_p(i4,j4,dxp4,dxn4)
    
      call fdf_p(i1,j1,p,dxp1,dxn1,fx1)
      call fdf_p(i2,j2,p,dxp2,dxn2,fx2)
      call fdf_p(i3,j3,p,dxp3,dxn3,fx3)
!       call fdf_p(i4,j4,p,dxp4,dxn4,fx4)

      fy1=-(fx1*xb)/yb
      fy2=-(fx2*xb)/yb
      fy3=-(fx3*xb)/yb
!       fy4=-(fx4*xb)/yb

!       first=(1.d0/q(i,j))*(w(1,i,j)*(fx1*del(1,i,j)+fy1*del(2,i,j))+w(2,i,j)*(fx2*del(1,i,j)+fy2&
!       *del(2,i,j))+w(3,i,j)*(fx3*del(1,i,j)+fy3*del(2,i,j))+w(4,i,j)*(fx4*del(1,i,j)+fy4*del(2,i,j)))
!         first=(1.d0/q(i,j))*(w(1,i,j)*(fx1*del(1,i,j)+fy1*del(2,i,j))+w(2,i,j)*(fx2*del(1,i,j)+fy2&
!       *del(2,i,j))+w(3,i,j)*(fx3*del(1,i,j)+fy3*del(2,i,j)))

        first=(1.d0/q(i,j))*(w(1,i,j)*(fx1*del(1,i,j)+fy1*del(2,i,j))+&
     w(2,i,j)*(fx2*del(1,i,j)+fy2*del(2,i,j))+&
     w(3,i,j)*(fx3*del(1,i,j)+fy3*del(2,i,j)))
      
      call fds_p(i1,j1,p,dxp1,dxn1,sx1)
      call fds_p(i2,j2,p,dxp2,dxn2,sx2)
      call fds_p(i3,j3,p,dxp3,dxn3,sx3)
!       call fds_p(i4,j4,p,dxp4,dxn4,sx4)
      
	amxy1=-(sx1*xb+fx1)/yb
	amxy2=-(sx2*xb+fx2)/yb
	amxy3=-(sx3*xb+fx3)/yb
! 	amxy4=-(sx4*xb+fx4)/yb	
	
	sy1=-(amxy1-fx1/yb)*(xb/yb)
	sy2=-(amxy2-fx2/yb)*(xb/yb)
	sy3=-(amxy3-fx3/yb)*(xb/yb)
! 	sy4=-(amxy4-fx4/yb)*(xb/yb)

!       second=(0.5d0/q(i,j))*(w(1,i,j)*(sx1*del(1,i,j)**2.d0+sy1*del(2,i,j)**2.d0)+w(2,i,j)*(sx2*&
!       del(1,i,j)**2.d0+sy2*del(2,i,j)**2.d0)+w(3,i,j)*(sx3*del(1,i,j)**2.d0+sy3*del(2,i,j)**2.d0&
!       +w(4,i,j)*(sx4*del(1,i,j)**2.d0+sy4*del(2,i,j)**2.d0)))        

!       amixed=(1.d0/q(i,j))*(w(1,i,j)*(amxy1*del(1,i,j)*del(2,i,j))+w(2,i,j)*(amxy2*del(1,i,j)*del(2,i,j))&
!       +w(3,i,j)*(amxy3*del(1,i,j)*del(2,i,j))+w(4,i,j)*(amxy4*del(1,i,j)*del(2,i,j)))

!         second=(0.5d0/q(i,j))*(w(1,i,j)*(sx1*del(1,i,j)**2.d0+sy1*del(2,i,j)**2.d0)+w(2,i,j)*(sx2*&
!       del(1,i,j)**2.d0+sy2*del(2,i,j)**2.d0)+w(3,i,j)*(sx3*del(1,i,j)**2.d0+sy3*del(2,i,j)**2.d0))
! 
!         amixed=(1.d0/q(i,j))*(w(1,i,j)*(amxy1*del(1,i,j)*del(2,i,j))+w(2,i,j)*(amxy2*del(1,i,j)*del(2,i,j))&
!       +w(3,i,j)*(amxy3*del(1,i,j)*del(2,i,j)))
      
      second=(0.5d0/q(i,j))*(w(1,i,j)*(sx1*del(1,i,j)**2.d0+&
     sy1*del(2,i,j)**2.d0)+w(2,i,j)*(sx2*del(1,i,j)**2.d0+&
     sy2*del(2,i,j)**2.d0)+w(3,i,j)*(sx3*del(1,i,j)**2.d0+&
     sy3*del(2,i,j)**2.d0))

      amixed=(1.d0/q(i,j))*(w(1,i,j)*(amxy1*del(1,i,j)*del(2,i,j))+&
     w(2,i,j)*(amxy2*del(1,i,j)*del(2,i,j))+&
     w(3,i,j)*(amxy3*del(1,i,j)*del(2,i,j)))
      
      p(i,j)=(1.d0/q(i,j))*(w(1,i,j)*p(i1,j1)+w(2,i,j)*p(i2,j2)+&
     w(3,i,j)*p(i3,j3))+first+second+amixed     
            
!       p(i,j)=(1.d0/q(i,j))*(w(1,i,j)*p(i1,j1)+w(2,i,j)*p(i2,j2)+w(3,i,j)*p(i3,j3)+w(4,i,j)*p(i4,j4))+first+second+amixed
         
      endif
      enddo
      enddo

      return
      end
! c****************
subroutine IBM_ubc
        use global  
        implicit none
      integer::i,j,l,i1,j1,i2,j2,i3,j3
      double precision::Aphi,VP,VE,VN,VS,fwp,fww,fep,fee,fnp,fnn,fsp,fss,aw,ae,as,an,ap,first,second,amixed,&
      dxp1,dxn1,dyp1,dyn1,fx1,fy1,sx1,sy1,amxy1,&
      dxp2,dxn2,dyp2,dyn2,fx2,fy2,sx2,sy2,amxy2,&
      dxp3,dxn3,dyp3,dyn3,fx3,fy3,sx3,sy3,amxy3
      
      do i=2,N
      do j=2,M
      
      if(icomput(i,j).eq.2)then

      i1=ingbr(1,1,i,j);i2=ingbr(1,2,i,j);i3=ingbr(1,3,i,j)
!       i4=ingbr(1,4,i,j)
      j1=ingbr(2,1,i,j);j2=ingbr(2,2,i,j);j3=ingbr(2,3,i,j)
!       j4=ingbr(2,4,i,j)      

      call coefficients(i1,j1,dxp1,dxn1,dyp1,dyn1)
      call coefficients(i2,j2,dxp2,dxn2,dyp2,dyn2)
      call coefficients(i3,j3,dxp3,dxn3,dyp3,dyn3)
!       call coefficients(i4,j4,dxp4,dxn4,dyp4,dyn4)
      
      do l=1,2
      call fdf(i1,j1,u,l,dxp1,dxn1,dyp1,dyn1,fx1,fy1)
      call fdf(i2,j2,u,l,dxp2,dxn2,dyp2,dyn2,fx2,fy2)
      call fdf(i3,j3,u,l,dxp3,dxn3,dyp3,dyn3,fx3,fy3)
!       call fdf(i4,j4,u,l,dxp4,dxn4,dyp4,dyn4,fx4,fy4)

      call fds(i1,j1,u,l,dxp1,dxn1,dyp1,dyn1,sx1,sy1)
      call fds(i2,j2,u,l,dxp2,dxn2,dyp2,dyn2,sx1,sy2)
      call fds(i3,j3,u,l,dxp3,dxn3,dyp3,dyn3,sx3,sy3)
!       call fds(i4,j4,u,l,dxp4,dxn4,dyp4,dyn4,sx4,sy4)

      call fdm(i1,j1,u,l,dxp1,dxn1,dyp1,dyn1,amxy1)
      call fdm(i2,j2,u,l,dxp2,dxn2,dyp2,dyn2,amxy2)
      call fdm(i3,j3,u,l,dxp3,dxn3,dyp3,dyn3,amxy3)
!       call fdm(i4,j4,u,l,dxp4,dxn4,dyp4,dyn4,amxy4)

!       first=(1.d0/q(i,j))*(w(1,i,j)*(fx1*del(1,i,j)+fy1*del(2,i,j))+w(2,i,j)*(fx2*del(1,i,j)+fy2*del(2,i,j))&
!       +w(3,i,j)*(fx3*del(1,i,j)+fy3*del(2,i,j))+w(4,i,j)*(fx4*del(1,i,j)+fy4*del(2,i,j)))
! 
!       second=(0.5d0/q(i,j))*(w(1,i,j)*(sx1*del(1,i,j)**2.d0+sy1*del(2,i,j)**2.d0)+w(2,i,j)*(sx2*del(1,i,j)**2.d0&
!       +sy2*del(2,i,j)**2.d0)+w(3,i,j)*(sx3*del(1,i,j)**2.d0+sy3*del(2,i,j)**2.d0)+w(4,i,j)*(sx4*del(1,i,j)**2.d0+sy4*del(2,i,j)**2.d0))
! 
!       amixed=(1.d0/q(i,j))*(w(1,i,j)*(amxy1*del(1,i,j)*del(2,i,j))+w(2,i,j)*(amxy2*del(1,i,j)*del(2,i,j))&
!       +w(3,i,j)*(amxy3*del(1,i,j)*del(2,i,j))+w(4,i,j)*(amxy4*del(1,i,j)*del(2,i,j)))

      first=(1.d0/q(i,j))*(w(1,i,j)*(fx1*del(1,i,j)+fy1*del(2,i,j))+&
     w(2,i,j)*(fx2*del(1,i,j)+fy2*del(2,i,j))+&
     w(3,i,j)*(fx3*del(1,i,j)+fy3*del(2,i,j)))

      second=(0.5d0/q(i,j))*(w(1,i,j)*(sx1*del(1,i,j)**2.d0+&
     sy1*del(2,i,j)**2.d0)+w(2,i,j)*(sx2*del(1,i,j)**2.d0+&
     sy2*del(2,i,j)**2.d0)+w(3,i,j)*(sx3*del(1,i,j)**2.d0+&
     sy3*del(2,i,j)**2.d0))

      amixed=(1.d0/q(i,j))*(w(1,i,j)*(amxy1*del(1,i,j)*del(2,i,j))+&
     w(2,i,j)*(amxy2*del(1,i,j)*del(2,i,j))+&
     w(3,i,j)*(amxy3*del(1,i,j)*del(2,i,j)))
                
!       u(l,i,j)=first+second+amixed
        if(l.eq.1)u(l,i,j)=velx(i,j)+first+second+amixed
      	if(l.eq.2)u(l,i,j)=vely(i,j)+first+second+amixed
      
      enddo
      
       endif
      enddo
      enddo

      return
      end
! c****************
subroutine A_operatep(phi,Aphi,ap,idd,jd)
      use global 
      implicit none
      integer::idd,jd,l,i1,j1,i2,j2,i3,j3
      double precision::Aphi,ap,first,second,amixed,&
      dxp1,dxn1,fx1,fy1,sx1,sy1,amxy1,&
      dxp2,dxn2,fx2,fy2,sx2,sy2,amxy2,&
      dxp3,dxn3,fx3,fy3,sx3,sy3,amxy3,xb,yb
      double precision,dimension(N+1,M+1)::phi
          
      if(icomput(idd,jd).eq.1)then

      Aphi=app(idd,jd)*phi(idd,jd)+aep(idd,jd)*phi(idd+1,jd)+anp(idd,jd)*phi(idd,jd+1)+&
     awp(idd,jd)*phi(idd-1,jd)+asp(idd,jd)*phi(idd,jd-1)
      ap=app(idd,jd)

      elseif(icomput(idd,jd).eq.2)then
      
      i1=ingbr(1,1,idd,jd);i2=ingbr(1,2,idd,jd);i3=ingbr(1,3,idd,jd)
!       i4=ingbr(1,4,i,j)
      j1=ingbr(2,1,idd,jd);j2=ingbr(2,2,idd,jd);j3=ingbr(2,3,idd,jd)
!       j4=ingbr(2,4,i,j)      
      
      xb=xcell(1,idd,jd)-del(1,idd,jd)
      yb=xcell(2,idd,jd)-del(2,idd,jd)
      
      call coeff_p(i1,j1,dxp1,dxn1)
      call coeff_p(i2,j2,dxp2,dxn2)
      call coeff_p(i3,j3,dxp3,dxn3)
!       call coeff_p(i4,j4,dxp4,dxn4)
    
      call fdf_p(i1,j1,phi,dxp1,dxn1,fx1)
      call fdf_p(i2,j2,phi,dxp2,dxn2,fx2)
      call fdf_p(i3,j3,phi,dxp3,dxn3,fx3)
!       call fdf_p(i4,j4,phi,dxp4,dxn4,fx4)
      
      fy1=-(fx1*xb)/yb
      fy2=-(fx2*xb)/yb
      fy3=-(fx3*xb)/yb
!       fy4=-(fx4*xb)/yb
        
!       first=(1.d0/q(i,j))*(w(1,i,j)*(fx1*del(1,i,j)+fy1*del(2,i,j))+w(2,i,j)*(fx2*del(1,i,j)+fy2&
!       *del(2,i,j))+w(3,i,j)*(fx3*del(1,i,j)+fy3*del(2,i,j))+w(4,i,j)*(fx4*del(1,i,j)+fy4*del(2,i,j)))

      first=(1.d0/q(idd,jd))*(w(1,idd,jd)*(fx1*del(1,idd,jd)+fy1*del(2,idd,jd))+&
     w(2,idd,jd)*(fx2*del(1,idd,jd)+fy2*del(2,idd,jd))+&
     w(3,idd,jd)*(fx3*del(1,idd,jd)+fy3*del(2,idd,jd)))
        
      call fds_p(i1,j1,phi,dxp1,dxn1,sx1)
      call fds_p(i2,j2,phi,dxp2,dxn2,sx2)
      call fds_p(i3,j3,phi,dxp3,dxn3,sx3)
!       call fds_p(i4,j4,phi,dxp4,dxn4,sx4)
      
	amxy1=-(sx1*xb+fx1)/yb
	amxy2=-(sx2*xb+fx2)/yb
	amxy3=-(sx3*xb+fx3)/yb
! 	amxy4=-(sx4*xb+fx4)/yb	
	
	sy1=-(amxy1-fx1/yb)*(xb/yb)
	sy2=-(amxy2-fx2/yb)*(xb/yb)
	sy3=-(amxy3-fx3/yb)*(xb/yb)
! 	sy4=-(amxy4-fx4/yb)*(xb/yb)

!       second=(0.5d0/q(i,j))*(w(1,i,j)*(sx1*del(1,i,j)**2.d0+sy1*del(2,i,j)**2.d0)+w(2,i,j)*(sx2*&
!       del(1,i,j)**2.d0+sy2*del(2,i,j)**2.d0)+w(3,i,j)*(sx3*del(1,i,j)**2.d0+sy3*del(2,i,j)**2.d0&
!       +w(4,i,j)*(sx4*del(1,i,j)**2.d0+sy4*del(2,i,j)**2.d0)))
! 
!       amixed=(1.d0/q(i,j))*(w(1,i,j)*(amxy1*del(1,i,j)*del(2,i,j))+w(2,i,j)*(amxy2*del(1,i,j)*del(2,i,j))&
!       +w(3,i,j)*(amxy3*del(1,i,j)*del(2,i,j))+w(4,i,j)*(amxy4*del(1,i,j)*del(2,i,j)))
!             
!       Aphi=(1.d0/q(i,j))*(w(1,i,j)*phi(i1,j1)+w(2,i,j)*phi(i2,j2)+w(3,i,j)*phi(i3,j3)+w(4,i,j)*phi(i4,j4))+first+second+amixed-phi(i,j)
      
      second=(0.5d0/q(idd,jd))*(w(1,idd,jd)*(sx1*del(1,idd,jd)**2.d0+&
     sy1*del(2,idd,jd)**2.d0)+w(2,idd,jd)*(sx2*del(1,idd,jd)**2.d0+&
     sy2*del(2,idd,jd)**2.d0)+w(3,idd,jd)*(sx3*del(1,idd,jd)**2.d0+&
     sy3*del(2,idd,jd)**2.d0))

      amixed=(1.d0/q(idd,jd))*(w(1,idd,jd)*(amxy1*del(1,idd,jd)*del(2,idd,jd))+&
     w(2,idd,jd)*(amxy2*del(1,idd,jd)*del(2,idd,jd))+&
     w(3,idd,jd)*(amxy3*del(1,idd,jd)*del(2,idd,jd)))
      
      Aphi=(1.d0/q(idd,jd))*(w(1,idd,jd)*phi(i1,j1)+w(2,idd,jd)*phi(i2,j2)+&
     w(3,idd,jd)*phi(i3,j3))+first+second+amixed-phi(idd,jd)
      
      
      ap=-1.d0

      endif

      return
      end 
! c**********************
subroutine A_operateu(phi,Aphi,idd,jd,ap,ld)
      use global 
      implicit none
      integer::idd,jd,ld,i1,j1,i2,j2,i3,j3
      double precision::Aphi,VP,VW,VE,VN,VS,fwp,fww,fep,fee,fnp,fnn,fsp,fss,factor,aw,ae,as,an,ap,first,second,amixed,&
      dxp1,dxn1,dyp1,dyn1,fx1,fy1,sx1,sy1,amxy1,&
      dxp2,dxn2,dyp2,dyn2,fx2,fy2,sx2,sy2,amxy2,&
      dxp3,dxn3,dyp3,dyn3,fx3,fy3,sx3,sy3,amxy3
      double precision,dimension(2,N+1,M+1)::phi
      
      if(icomput(idd,jd).eq.1)then

      VP=vol(idd,jd)
      VE=vol(idd+1,jd)
      VW=vol(idd-1,jd)
      VN=vol(idd,jd+1)
      VS=vol(idd,jd-1)
      
      factor=dt/(2.d0*VP)

      if(l_CN.ge.1)then
	  fwp=VW/(VW+VP)
	  fww=VP/(VW+VP)
	  
      fep=VE/(VE+VP)
	  fee=VP/(VE+VP)
	  
      fsp=VS/(VS+VP)
	  fss=VP/(VS+VP)
	  
      fnp=VN/(VN+VP)
      fnn=VP/(VN+VP)
	  
      aw=factor*F(1,idd,jd)*fww+awv(idd,jd)
      ae=factor*F(3,idd,jd)*fee+aev(idd,jd)
      an=factor*F(4,idd,jd)*fnn+anv(idd,jd)
      as=factor*F(2,idd,jd)*fss+asv(idd,jd)
      ap=factor*(F(1,idd,jd)*fwp+F(2,idd,jd)*fsp+F(3,idd,jd)*fep+F(4,idd,jd)*fnp)+apv(idd,jd)

      elseif(l_ABCN.eq.1)then
      aw=awv(idd,jd)
      ae=aev(idd,jd)
      an=anv(idd,jd)
      as=asv(idd,jd)
      ap=apv(idd,jd)
      endif
	 
	  if(idd.eq.2)then
      aw=0.d0  
      endif
	  
      if(jd.eq.2)then
      if(ld.eq.1)ap=ap+as
      as=0.d0  
      endif

      if(jd.eq.M)then
      if(ld.eq.1)ap=ap+an
      an=0.d0  
	  endif
    
	  if(idd.eq.N)then
      ap=ap+ae
      ae=0.d0  
      endif
	  
      Aphi=ap*phi(ld,idd,jd)+ae*phi(ld,idd+1,jd)+an*phi(ld,idd,jd+1)+&
     aw*phi(ld,idd-1,jd)+as*phi(ld,idd,jd-1)  
      
      elseif(icomput(idd,jd).eq.2)then

      i1=ingbr(1,1,idd,jd);i2=ingbr(1,2,idd,jd);i3=ingbr(1,3,idd,jd)
!       i4=ingbr(1,4,i,j)
      j1=ingbr(2,1,idd,jd);j2=ingbr(2,2,idd,jd);j3=ingbr(2,3,idd,jd)
!       j4=ingbr(2,4,i,j)      

      call coefficients(i1,j1,dxp1,dxn1,dyp1,dyn1)
      call coefficients(i2,j2,dxp2,dxn2,dyp2,dyn2)
      call coefficients(i3,j3,dxp3,dxn3,dyp3,dyn3)
!       call coefficients(i4,j4,dxp4,dxn4,dyp4,dyn4)         
      
      call fdf(i1,j1,phi,ld,dxp1,dxn1,dyp1,dyn1,fx1,fy1)
      call fdf(i2,j2,phi,ld,dxp2,dxn2,dyp2,dyn2,fx2,fy2)
      call fdf(i3,j3,phi,ld,dxp3,dxn3,dyp3,dyn3,fx3,fy3)
!       call fdf(i4,j4,phi,l,dxp4,dxn4,dyp4,dyn4,fx4,fy4)
      
      call fds(i1,j1,phi,ld,dxp1,dxn1,dyp1,dyn1,sx1,sy1)
      call fds(i2,j2,phi,ld,dxp2,dxn2,dyp2,dyn2,sx1,sy2)
      call fds(i3,j3,phi,ld,dxp3,dxn3,dyp3,dyn3,sx3,sy3)
!       call fds(i4,j4,phi,l,dxp4,dxn4,dyp4,dyn4,sx4,sy4)

      call fdm(i1,j1,phi,ld,dxp1,dxn1,dyp1,dyn1,amxy1)
      call fdm(i2,j2,phi,ld,dxp2,dxn2,dyp2,dyn2,amxy2)
      call fdm(i3,j3,phi,ld,dxp3,dxn3,dyp3,dyn3,amxy3)
!       call fdm(i4,j4,phi,l,dxp4,dxn4,dyp4,dyn4,amxy4)

!       first=(1.d0/q(i,j))*(w(1,i,j)*(fx1*del(1,i,j)+fy1*del(2,i,j))+w(2,i,j)*(fx2*del(1,i,j)+fy2*del(2,i,j))&
!       +w(3,i,j)*(fx3*del(1,i,j)+fy3*del(2,i,j))+w(4,i,j)*(fx4*del(1,i,j)+fy4*del(2,i,j)))
! 
!       second=(0.5d0/q(i,j))*(w(1,i,j)*(sx1*del(1,i,j)**2.d0+sy1*del(2,i,j)**2.d0)+w(2,i,j)*(sx2*del(1,i,j)**2.d0&
!       +sy2*del(2,i,j)**2.d0)+w(3,i,j)*(sx3*del(1,i,j)**2.d0+sy3*del(2,i,j)**2.d0)+w(4,i,j)*(sx4*del(1,i,j)**2.d0+sy4*del(2,i,j)**2.d0))
! 
!       amixed=(1.d0/q(i,j))*(w(1,i,j)*(amxy1*del(1,i,j)*del(2,i,j))+w(2,i,j)*(amxy2*del(1,i,j)*del(2,i,j))&
!       +w(3,i,j)*(amxy3*del(1,i,j)*del(2,i,j))+w(4,i,j)*(amxy4*del(1,i,j)*del(2,i,j)))

!         first=(1.d0/q(i,j))*(w(1,i,j)*(fx1*del(1,i,j)+fy1*del(2,i,j))+w(2,i,j)*(fx2*del(1,i,j)+fy2*del(2,i,j))&
!       +w(3,i,j)*(fx3*del(1,i,j)+fy3*del(2,i,j)))
! 
!       second=(0.5d0/q(i,j))*(w(1,i,j)*(sx1*del(1,i,j)**2.d0+sy1*del(2,i,j)**2.d0)+w(2,i,j)*(sx2*del(1,i,j)**2.d0&
!       +sy2*del(2,i,j)**2.d0)+w(3,i,j)*(sx3*del(1,i,j)**2.d0+sy3*del(2,i,j)**2.d0))
! 
!       amixed=(1.d0/q(i,j))*(w(1,i,j)*(amxy1*del(1,i,j)*del(2,i,j))+w(2,i,j)*(amxy2*del(1,i,j)*del(2,i,j))&
!       +w(3,i,j)*(amxy3*del(1,i,j)*del(2,i,j)))
        first=(1.d0/q(idd,jd))*(w(1,idd,jd)*(fx1*del(1,idd,jd)+fy1*del(2,idd,jd))+&
     w(2,idd,jd)*(fx2*del(1,idd,jd)+fy2*del(2,idd,jd))+&
     w(3,idd,jd)*(fx3*del(1,idd,jd)+fy3*del(2,idd,jd)))

      second=(0.5d0/q(idd,jd))*(w(1,idd,jd)*(sx1*del(1,idd,jd)**2.d0+&
     sy1*del(2,idd,jd)**2.d0)+w(2,idd,jd)*(sx2*del(1,idd,jd)**2.d0+&
     sy2*del(2,idd,jd)**2.d0)+w(3,idd,jd)*(sx3*del(1,idd,jd)**2.d0+&
     sy3*del(2,idd,jd)**2.d0))

      amixed=(1.d0/q(idd,jd))*(w(1,idd,jd)*(amxy1*del(1,idd,jd)*del(2,idd,jd))+&
     w(2,idd,jd)*(amxy2*del(1,idd,jd)*del(2,idd,jd))+&
     w(3,idd,jd)*(amxy3*del(1,idd,jd)*del(2,idd,jd)))
                
!       u(l,i,j)=first+second+amixed
!         if(l.eq.1)u(l,i,j)=velx(i,j)+first+second+amixed
!       	if(l.eq.2)u(l,i,j)=vely(i,j)+first+second+amixed
      Aphi=first+second+amixed-phi(ld,idd,jd)
      ap=-1.d0

      endif

      return
      end 
! c***************
subroutine final_data
        use global  
        implicit none
        integer::i,j,l
      double precision::dxp,dxn,dyp,dyn,fx,fy,ux,uy,vx,vy,sumu,sumv,flx,fly,si
      double precision,dimension(:,:,:),allocatable::u_i
      double precision,dimension(:,:),allocatable::qst,vorz,sif                
      
      allocate(u_i(2,N+1,M+1))          
      allocate(qst(N+1,M+1),vorz(N+1,M+1),sif(N+1,M+1))  

      do i=1,N+1
      do j=1,M+1      
      u_i(1,i,j)=u(1,i,j)+v_sx
      u_i(2,i,j)=u(2,i,j)+v_s    
      qst(i,j)=0.d0;vorz(i,j)=0.d0
      enddo
      enddo
      
      do i=2,N
      do j=2,M

      if(icomput(i,j).eq.1)then
      call coefficients(i,j,dxp,dxn,dyp,dyn)

      do l=1,2
      call fdf(i,j,u_i,l,dxp,dxn,dyp,dyn,fx,fy)

      if(l.eq.1)then
      ux=fx;uy=fy
      endif
      if(l.eq.2)then
      vx=fx;vy=fy
      endif
      
      enddo

      qst(i,j)=-(ux**2.d0+vy**2.d0+2.d0*uy*vx)/2.d0           
      vorz(i,j)=vx-uy        
      
      flx=x(1,i+1,j)-x(1,i,j)
      fly=x(2,i,j+1)-x(2,i,j)      
      sumu=sumu+u_i(1,i,j)*(dt/flx)
      sumv=sumv+u_i(2,i,j)*(dt/fly)
      
      endif
    
      enddo
      enddo
      
      do i=2,N
      si=0.d0
      do j=2,M
      if(icomput(i,j).eq.1)then
      si=si+(u_i(1,i,j)+u_i(1,i,j-1))*(xcell(2,i,j)-xcell(2,i,j-1))/2.d0
      sif(i,j)=si
      else
      si=1000.d0
      sif(i,j)=si
      endif
      enddo
      enddo
!             print*,y_s;stop 
      open(33,file='result.dat')
      write(33,100)N-1,M-1            
      do j=2,M
      do i=2,N
      write(33,201)xcell(1,i,j)+y_sx,xcell(2,i,j)+y_s,&
     u_i(1,i,j),u_i(2,i,j),p(i,j),qst(i,j),vorz(i,j),sif(i,j)
      enddo
      enddo
      close(33)     
      
      deallocate(u_i)
      deallocate(qst,vorz,sif)
100   format("ZONE",1x,"I=",i3,1x,"J=",i3)
201   format(8(f15.8,1x))

      return
      end 
! c***************
subroutine position
      use global  
      implicit none
      integer::j,Nth,Nr
      double precision::dtheta,xsph,ysph,theta
      
      Nr=1;Nth=360                  
      open(27,file='position.dat')     
      dtheta=2.d0*pi/dfloat(Nth-1)      
      write(27,100)1,Nth
      do j=1,Nth	
      theta=dtheta*dfloat(j-1)  
      xsph=radius*dcos(theta)
      ysph=radius*dsin(theta)
      write(27,*)xsph+y_sx,ysph+y_s
      enddo	
      close(27)      

100   format("ZONE",1x,"I=",i3,1x,"J=",i3)

      return
      end
! c*****************
subroutine volume_mom
      use global  
      implicit none
      integer::i,j 
      volume=0.d0
      do i=iu1,id1
      do j=jb,jt      
      if(icomput(i,j).eq.1)then      
	volume=volume+vol(i,j)
      endif
      enddo
      enddo
      	
      return
      end
! c******************
subroutine neighbors_new
      use global  
      implicit none
      integer::i,j,lg,l1,ll,i1,j1,itemp,jtemp 
      double precision::h1,h2,h3,Rdis,xb,yb,x1,y1,dis,temp         
      double precision,dimension(:),allocatable::xg
      integer,dimension(:),allocatable::in,jn
      
      open(300,file='weights1.dat')
      lg=0
      do i=2,N
      do j=2,M
      
      if(icomput(i,j).eq.2)then
      xb=xcell(1,i,j)-del(1,i,j)
      yb=xcell(2,i,j)-del(2,i,j)
      
      lg=lg+1

      l1=0
      do i1=max(i-lspace,2),min(i+lspace,N)
      do j1=max(j-lspace,2),min(j+lspace,M)
            
      if(icomput(i1,j1).eq.1)l1=l1+1     
      
      enddo
      enddo
      
      allocate(xg(l1),in(l1),jn(l1))
      l1=0
      do i1=max(i-lspace,2),min(i+lspace,N)
      do j1=max(j-lspace,2),min(j+lspace,M)
      if(icomput(i1,j1).eq.1)then
      x1=xcell(1,i1,j1)
      y1=xcell(2,i1,j1)      
      l1=l1+1
      dis=dsqrt((x1-xb)**2.d0+(y1-yb)**2.d0)
      xg(l1)=dis;in(l1)=i1;jn(l1)=j1
      endif
      enddo      
      enddo

      do ll=1,l1
      do i1=1,l1
      if(xg(i1).ge.xg(ll))then
      temp=xg(ll)
      xg(ll)=xg(i1)
      xg(i1)=temp
      itemp=in(ll);jtemp=jn(ll)
      in(ll)=in(i1);jn(ll)=jn(i1)
      in(i1)=itemp;jn(i1)=jtemp
      endif
      enddo
      enddo
  
      h1=xg(1);ingbr(1,1,i,j)=in(1)
      ingbr(2,1,i,j)=jn(1)

      h2=xg(2);ingbr(1,2,i,j)=in(2)
      ingbr(2,2,i,j)=jn(2)

      h3=xg(3);ingbr(1,3,i,j)=in(3)
      ingbr(2,3,i,j)=jn(3)
      
      Rdis=xg(4)
      if(xg(4).eq.xg(3))Rdis=xg(5)
!       i_flag=0;i1=0      
!       do while(i_flag.eq.0)
!       i1=i1+1
!       if(xg(i1).gt.h3)then
!       i_flag=1
!       Rdis=xg(i1)
!       endif
!       enddo
! 
!       h4=xg(4);ingbr(1,4,i,j)=in(4)
!       ingbr(2,4,i,j)=jn(4)
!       Rdis=xg(5)
!       Rdis=xg(4)
      
      w(1,i,j)=((Rdis-h1)/(Rdis*h1))**2.d0
      w(2,i,j)=((Rdis-h2)/(Rdis*h2))**2.d0
      w(3,i,j)=((Rdis-h3)/(Rdis*h3))**2.d0 
!       w(4,i,j)=((Rdis-h4)/(Rdis*h4))**2.d0

!       q(i,j)=w(1,i,j)+w(2,i,j)+w(3,i,j)+w(4,i,j)
        q(i,j)=w(1,i,j)+w(2,i,j)+w(3,i,j)
      write(300,100)lg,h1,h2,h3,Rdis
!         write(300,100)(w(l,i,j),l=1,3),q(i,j)
      deallocate(xg,in,jn)

      endif
      enddo
      enddo      
100   format(i3,1x,4(f20.13,1x))
      close(300)
      print*,'weights1'
 
      return
      end
! c***************
subroutine cfl
        use global  
        implicit none
        integer::i,j
        double precision::flx,fly,sumu,sumv,cflu,cflv
      
      sumu=0.d0
      sumv=0.d0
      do i=2,N
      do j=2,M      
      if(icomput(i,j).eq.1)then
      flx=x(1,i+1,j)-x(1,i,j)
      fly=x(2,i,j+1)-x(2,i,j)      
      sumu=sumu+u(1,i,j)*(dt/flx)
      sumv=sumv+u(2,i,j)*(dt/fly)      
      endif
      enddo
      enddo
      
      cflu=sumu/dfloat(Nmax)
      cflv=sumv/dfloat(Nmax)      
      write(43,*)time,cflu,cflv  

      return
      end 
! c***************
subroutine momentum_init
      use global 
      implicit none
      integer::i,j,l
      double precision::unstdyx,unstdyy,xmomd,xmomu,xmomb,xmomt,&
      ymomd,ymomu,ymomb,ymomt,ue,uw,ve,vw,pfd,pfu,pft,pfb,prx,pry,&
      Fcx,Fcy,un,vn,us,vs,pe,pw,pn,ps
      
!       unstdyx=0.d0
!       unstdyy=0.d0
!       do i=iu1,id1
!       do j=jb,jt
!       if(icomput(i,j).eq.1)then
!       unstdyx=unstdyx+vol(i,j)*(u(1,i,j)-u_old(1,i,j))/dt
!       unstdyy=unstdyy+vol(i,j)*(u(2,i,j)-u_old(2,i,j))/dt
!       endif
!       enddo
!       enddo

      xmomt=0.d0
      xmomb=0.d0
      xmomu=0.d0
      xmomd=0.d0
      ymomt=0.d0
      ymomb=0.d0
      ymomu=0.d0
      ymomd=0.d0
      pfd=0.d0
      pfu=0.d0
      pft=0.d0
      pfb=0.d0
      Fcx1=0.d0
      Fcy1=0.d0
      prx1=0.d0
      pry1=0.d0
      
      do j=jb,jt
      call face_intpo(u(1,id1+1,j),u(1,id1,j),vol(id1+1,j),&
      vol(id1,j),ue)
      call face_intpo(u(2,id1+1,j),u(2,id1,j),vol(id1+1,j),&
      vol(id1,j),ve)
      call face_intpo(u(1,iu1-1,j),u(1,iu1,j),vol(iu1-1,j),&
      vol(iu1,j),uw)
      call face_intpo(u(2,iu1-1,j),u(2,iu1,j),vol(iu1-1,j),&
      vol(iu1,j),vw)
      call face_intpo(p(id1+1,j),p(id1,j),vol(id1+1,j),&
      vol(id1,j),pe)
      call face_intpo(p(iu1-1,j),p(iu1,j),vol(iu1-1,j),&
      vol(iu1,j),pw)

      xmomd=xmomd+F(3,id1,j)*ue
      xmomu=xmomu+F(1,iu1,j)*uw
      ymomd=ymomd+F(3,id1,j)*ve
      ymomu=ymomu+F(1,iu1,j)*vw
      pfd=pfd+pe*s(1,3,id1,j)      
      pfu=pfu+pw*s(1,1,iu1,j)
      enddo

      do i=iu1,id1
      call face_intpo(u(1,i,jt+1),u(1,i,jt),vol(i,jt+1),vol(i,jt),un)
      call face_intpo(u(2,i,jt+1),u(2,i,jt),vol(i,jt+1),vol(i,jt),vn)
      call face_intpo(u(1,i,jb-1),u(1,i,jb),vol(i,jb-1),vol(i,jb),us)
      call face_intpo(u(2,i,jb-1),u(2,i,jb),vol(i,jb-1),vol(i,jb),vs)
      call face_intpo(p(i,jt+1),p(i,jt),vol(i,jt+1),vol(i,jt),pn)
      call face_intpo(p(i,jb-1),p(i,jb),vol(i,jb-1),vol(i,jb),ps)

      xmomt=xmomt+F(4,i,jt)*un
      xmomb=xmomb+F(2,i,jb)*us
      ymomt=ymomt+F(4,i,jt)*vn
      ymomb=ymomb+F(2,i,jb)*vs
      pft=pft+pn*s(2,4,i,jt)
      pfb=pfb+ps*s(2,2,i,jb)
      enddo

      Fcx1=xmomd+xmomu+xmomt+xmomb
      Fcy1=ymomt+ymomb+ymomd+ymomu
      prx1=pfu+pfd
      pry1=pfb+pft
!       cd=(prx+prx1)+(Fcx+Fcx1)+2.d0*unstdyx
!       cl=(pry+pry1)+(Fcy+Fcy1)+2.d0*unstdyy
!       Fcx1=Fcx
!       Fcy1=Fcy
!       prx1=prx
!       pry1=pry     
!       
!       cl=cl+pseudo_fac*volume
!       c_v=-cl-c_p
!       e_v=c_v*v_s
! 
!       write(41,200)time,-cd,-cl,c_p,c_v,e_v,y_s
! 200   format(7(f16.9,1x))      
      
      return
      end
! c******************************************
subroutine volume_fraction
      use global 
      implicit none
      integer::i,j
      double precision::x_max,y_max,x_min,y_min,dis_max,dis_min,x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,arear,area5,area6,area
      print*,'1'
      do i=1,N+1
      do j=1,M+1
      x_max=x(1,i,j)
      y_max=x(2,i,j)
      x_min=x(1,i+1,j+1)
      y_min=x(2,i+1,j+1)
      dis_max=dsqrt((y_max)**2.d0+(x_max)**2.d0)
      dis_min=dsqrt((y_min)**2.d0+(x_min)**2.d0)
      if(radius.ge.dis_max)then
      vol_sol(i,j)=0
      elseif(dis_min.ge.radius)then
      vol_sol(i,j)=1
      else
      x1=x(1,i,j)
      y1=x(2,i,j)
      x2=x(1,i+1,j)
      y2=x(2,i+1,j)
      x3=x(1,i+1,j+1)
      y3=x(2,i+1,j+1)
      x4=x(1,i,j+1)
      y4=x(2,i,j+1)
      arear=(x1-x2)*(y1-y4)
      x5=dsqrt(radius**2.d0-y3**2.d0)
      y5=dsqrt(radius**2.d0-x3**2.d0)
      x6=dsqrt(radius**2.d0-y1**2.d0)
      y6=dsqrt(radius**2.d0-x1**2.d0)
      area5=0.5*(x5-x3)*(y5-y3)
      area6=0.5*(x1-x6)*(y1-y6)
      if (((x3.ge.x5).and.(x5.ge.x4)).and.((y2.ge.y5).and.(y5.ge.y3)))then
      area=area5
      elseif(((x3.ge.x6).and.(x6.ge.x4)).and.((y2.ge.y6).and.(y6.ge.y3)))then
      area=arear-area6	
      else
      if ((x3.ge.x5).and.(x5.ge.x4))then
      area=0.5*(((x6-x3)*(y4-y5))-((x4-x5)*(y6-y3)))
      else
      area=0.5*(((x6-x3)*(y4-y5))-((x4-x5)*(y6-y3)))
      endif
      endif
      vol_sol(i,j)=area
      endif
      enddo
      enddo
      return
      end
!**************************************************** 
