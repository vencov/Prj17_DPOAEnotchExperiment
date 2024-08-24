function [G,Gs,x,W]= cgsgreenf(N, L, H);
%
%		[G,Gs,W,x]= CGSGREENF(N, L, H)
%
%
%		Input data:
%	      N = number of points on the basilar membrane (BM) [def. N=500]
%	      L = BM length [def. = 3.5 cm] 
%         H = scalae height [def. = 0.3 cm]
%
%		Output data:
%	      G = BM-BM Green's function matrix 
%         Gs=Stapes-BM Green's function
%         W = Sectional shear viscosity operator
%	      x = BM coordinate array [0 <= x < L]
plotflag = 1;

if nargin <3,
	H = 0.1; % spiral-canal-height  = 1/10 of BM length
end

if nargin <2,
    L = 3.5; % BM length in cm
end    
%H = H^2;
if nargin <1,
   N=800;
end   

x1=(0:N-1)./N;
x2=(1:N)./N;
x=L*((x1+x2)./2);

x=x(:)'; 

rho = 1;  
Ko=rho/H;
Co=2*rho;
dx=L/N;

[Ipo,Ibo]=G_peak_E(N,H,L);

G = zeros(N,N);
Go = zeros(N,N);
P=zeros(N,N);
B=zeros(N,N);
Go(1, :) = Ko*(2*L-x(1)-x-abs(x-x(1)));
P(1,:)=Ipo;
B(1,:)=Ibo(1:N);


for n=1:N,
    xo=x(n);
	Go(n, :) = Ko*(2*L-xo-x-abs(x-xo)); 
    P(n,1:n)=Ipo(n:-1:1);
    P(n,n+1:N)=Ipo(2:(N-n+1));
    B(n,:)=Ibo(n:N+n-1);
end

Go=Go*dx;
G=Go+P+B;
% G=Go;%+P+B;
P=P+B;
Gs=Co*(L-x); 
Gs=Gs';  
W = shearvisc(x,dx,N);

      

%------------------------------SubFunctions--------------------------------
function [Ipo,Ibo]=G_peak_E(N,H,L);
 if nargin<3
     L=3.5;
 end
 if nargin<2
   H=0.15;
 end
 if nargin<1
    N=600;
 end
  %H = 0.05;
   A=pi/H;
 %----------------Grid for P+Ipo------
   N1=N+1;
   N2=2*N;
   dx=L/N;
   Fac=((2*H)/(pi^2));
   xp(1)=0; 
   xp(2:N1)=((0:N1-2)+(1:N1-1))/2;
   xp=xp*dx;
   In=DilogRA(A*xp,400);
   In1=Fac*In;
   Ipo(1)=2*(In1(2)-In1(1));%symmetrical interval %-dx/2,dx2% around the peak
   Ipo(2:N)=In1(3:N1)-In1(2:N);
  %----------------Grid for B+Ibo------
   xb=0.5+(0:N2);%2*N length
   xb=xb*dx;
   In=DilogRA(A*xb,400);
   In1=((2*H)/(pi^2))*In;
   Ibo=In1(2:N2+1)-In1(1:N2);
   
   
function y=DilogRA(x,N)

 if nargin<2
    N=400;
 end
 if nargin<1
    x=(0:99)/99;
    x=x*10;
 end
 x=abs(x);
 y=zeros(size(x));
 for k=1:N
    y=y+(-exp(-k*x))/(k*k);
 end  


function W = shearvisc(x,dx,N)
  invdx2=ones(N,1)./(dx*dx);	
  W=diag(2*invdx2,0)-diag(invdx2(1:N-1),-1)-diag(invdx2(1:N-1),1);
  W(1,:)=0;%BC problem
  W(N,:)=0;%BC problem

