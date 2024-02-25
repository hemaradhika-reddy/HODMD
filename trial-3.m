clc;clear all;close all

f1 = 50;
fs = 3200;
T  = 1;
n = fs*T;
t = (0:n-1)/fs;
m = length(t);
dt = 1/((fs/f1)*f1); %1/fs;

%x = cos(2*pi*10*t) + 0.05*cos(2*3*pi*10*t) + 0.02*cos(2*pi*10*t*5);



%x = cos(2*pi*10*t) + 0.05*cos(2*0.6*pi*10*t) + 0.02*cos(2*3*pi*10*t);

%x = cos(2*pi*10*t) + 0.05*cos(2*pi*10*t*3) + 0.02*cos(2*pi*10*t*5) + 0.04*cos(2*pi*10*t*7) + 0.12*cos(2*pi*10*t*9) + 0.45*cos(2*pi*10*t*0.52) + 0.25*cos(2*pi*10*t*3.6) + 0.5*cos(2*pi*10*t*4.6);

%x = cos(2*pi*3.2*t) + 0.05*cos(2*pi*3.2*t*0.8) + 0.02*cos(2*pi*3.2*t*0.36);

x = cos(2*pi*t) + 0.05*cos(2*pi*3*t) + 0.02*cos(2*pi*t*5) +0.05*cos(2*pi*6.3*t);
x=   awgn(x,50);

%y = hankel(x);

%x=y;

plot(t,x)

%% Number of snapshots
nsnap=length(x);
V=x(1:nsnap);
Time = t;
%% DMD-d
d=1000;
varepsilon1=1e-10 ;%SVD tolerance
varepsilon=1e-3 ;%DMD tolerance

[M N]=size(V);

if d>1
    [Vreconst,deltas,omegas,amplitude,modes] =DMDd_SIADS(d,V,Time,varepsilon1,varepsilon);
else
    [Vreconst,deltas,omegas,amplitude] =DMD1_SIADS(V,Time,varepsilon1,varepsilon);
end 

figure; plot(real(modes));
figure;plot(V);hold on; 
figure;plot(real(Vreconst));

fest = omegas/(2*pi)
amps = amplitude*2

Est_out = [fest;amps]; % final estimated outputs

                            
%%%%%%%%%%%%%

function ContReconst=ContReconst_SIADS(t,t0,u,deltas,omegas)
[N,M]=size(u);
vv=zeros(M,1);
for m=1:M
 vv(m)=exp((deltas(m)+i*omegas(m))*(t-t0));   
end
ContReconst=u*vv;
end

%%%%%%%%%%%%%%%

function  [Vreconst,deltas,omegas,amplitude,modes] =DMDd_SIADS(d,V,Time,varepsilon1,varepsilon)
[J,K]=size(V);

%% STEP 1: SVD of the original data

[U,Sigma,T]=svd(V,'econ');
sigmas=diag(Sigma);
n=length(sigmas);

NormS=norm(sigmas,2);
p=0;
for k=1:n
    if norm(sigmas(k:n),2)/NormS>varepsilon1
        p=p+1;
    end
end

U=U(:,1:p);

%% Spatial complexity: kk
('Spatial complexity');
p;

%% Create reduced snapshots matrix
hatT=Sigma(1:p,1:p)*T(:,1:p)';
[N,~]=size(hatT);

%% Create the modified snapshot matrix
tildeT=zeros(d*N,K-d+1);
for ppp=1:d
    tildeT((ppp-1)*N+1:ppp*N,:)=hatT(:,ppp:ppp+K-d);
end

%% Dimension reduction
[U1,Sigma1,T1]=svd(tildeT,'econ');
sigmas1=diag(Sigma1);

Deltat=Time(2)-Time(1);
n=length(sigmas1);

NormS=norm(sigmas1,2);
p1=0;
for k=1:n
    RRMSEE(k)=norm(sigmas1(k:n),2)/NormS;
    if RRMSEE(k)>varepsilon1
        p1=p1+1;
    end
end

('Spatial dimension reduction');
p1;

U1=U1(:,1:p1);
hatT1=Sigma1(1:p1,1:p1)*T1(:,1:p1)';

%% Reduced modified snapshot matrix
[~,K1]=size(hatT1);
[tildeU1,tildeSigma,tildeU2]=svd(hatT1(:,1:K1-1),'econ');

%% Reduced modified Koopman matrix
tildeR=hatT1(:,2:K1)*tildeU2*inv(tildeSigma)*tildeU1';
[tildeQ,tildeMM]=eig(tildeR);
eigenvalues=diag(tildeMM);

M=length(eigenvalues);
qq=log(eigenvalues);
deltas=real(qq)/Deltat;
omegas=imag(qq)/Deltat;

Q=U1*tildeQ;
Q=Q((d-1)*N+1:d*N,:);
[NN,MMM]=size(Q);

for m=1:MMM
    NormQ=Q(:,m);
    Q(:,m)= Q(:,m)/norm(NormQ(:),2);
end

%% Calculate amplitudes
Mm=zeros(NN*K,M);
Bb=zeros(NN*K,1);
aa=eye(MMM);
for k=1:K
    Mm(1+(k-1)*NN:k*NN,:)=Q*aa;
    aa=aa*tildeMM;
    Bb(1+(k-1)*NN:k*NN,1)=hatT(:,k);
end

[Ur,Sigmar,Vr]=svd(Mm,'econ');
a=Vr*(Sigmar\(Ur'*Bb));

u=zeros(NN,M);
for m=1:M
    u(:,m)=a(m)*Q(:,m);
end
amplitude=zeros(M,1);

for m=1:M
    aca=U*u(:,m);
    amplitude(m)=norm(aca(:),2)/sqrt(J);
end

UU=[u;deltas';omegas';amplitude']';
UU1=sortrows(UU,-(NN+3));

UU=UU1';
u=UU(1:NN,:);
deltas=UU(NN+1,:);
omegas=UU(NN+2,:);
amplitude=UU(NN+3,:);
kk3=0;

for m=1:M
    if amplitude(m)/amplitude(1)>varepsilon
        kk3=kk3+1;
    else
    end
end

%% Spectral complexity: number of DMD modes.
('Spectral complexity');
kk3;
u=u(:,1:kk3);
deltas=deltas(1:kk3)
omegas=omegas(1:kk3)
amplitude=amplitude(1:kk3)

('Mode number, delta, omega, amplitude')
%DeltasOmegAmpl = [(1:kk3)',deltas',omegas',amplitude']

%% Reconstruction of the original snapshot matrix
hatTreconst=zeros(N,K);
for k=1:K
    hatTreconst(:,k)= ContReconst_SIADS(Time(k),Time(1),u,deltas,omegas);
end
    
Vreconst=U*hatTreconst;

%% Calculation of DMD modes
modes=zeros(J,kk3);
amplitude0=zeros(kk3,1);
for m=1:kk3
    NormMode=norm(U*u(:,m),2)/sqrt(J);
    amplitude0(m)=NormMode;
    modes(:,m)=U*u(:,m)/NormMode;
end

%If the calculation of the amplitudes is correct, ErrAmpl=0
%ErrAmpl=norm(amplitude(:,1:kk3)-amplitude0',2)
end

%%%%%%%%%%%%%%%%%%

function  [Vreconst,deltas,omegas,amplitude] =DMD1_SIADS(V,Time,varepsilon1,varepsilon)

[J,K]=size(V);
[U,Sigma,T]=svd(V,'econ');
sigmas=diag(Sigma);
Deltat=Time(2)-Time(1);
n=length(sigmas);

NormS=norm(sigmas,2);
kk=0;
for k=1:n
    if norm(sigmas(k:n),2)/NormS>varepsilon1
        kk=kk+1;
    end
end
%% Spatial complexity: kk
('Spatial complexity');
kk;

U=U(:,1:kk);
%% Create reduced snapshots matrix
hatT=Sigma(1:kk,1:kk)*T(:,1:kk)';
[N,K]=size(hatT);
[hatU1,hatSigma,hatU2]=svd(hatT(:,1:K-1),'econ');

%% Calculate Koopman operator
hatR=hatT(:,2:K)*hatU2*inv(hatSigma)*hatU1';
[Q,MM]=eig(hatR);

eigenvalues=diag(MM);

M=length(eigenvalues);
qq=log(eigenvalues);
deltas=real(qq)/Deltat;
omegas=imag(qq)/Deltat;

%% Calculate amplitudes
Mm=zeros(M*K,M);
Bb=zeros(M*K,1);
for k=1:K
    Mm(1+(k-1)*M:k*M,:)=Q*(MM^(k-1));
    Bb(1+(k-1)*M:k*M,1)=hatT(:,k);
end

[Ur,Sigmar,Vr]=svd(Mm,'econ');
a=Vr*(Sigmar\(Ur'*Bb));

u=zeros(M,M);
for m=1:M
    u(:,m)=a(m)*Q(:,m);
end

amplitude=zeros(M,1);
for m=1:M
    aca=U*u(:,m);
    amplitude(m)=norm(aca(:),2)/sqrt(J);
end
 
UU=[u;deltas';omegas';amplitude']';
UU1=sortrows(UU,-(M+3));

UU=UU1';
u=UU(1:M,:);
deltas=UU(M+1,:);
omegas=UU(M+2,:);
amplitude=UU(M+3,:);
kk2=0;
 
for m=1:M
    if amplitude(m)/amplitude(1)>varepsilon
        kk2=kk2+1;
    else
    end
end
%% Spectral complexity: number of DMD modes.
('Spectral complexity');
kk2;
u=u(:,1:kk2);
deltas=deltas(1:kk2);
omegas=omegas(1:kk2);
amplitude=amplitude(1:kk2);
('Mode number, delta, omega, amplitude')
DeltasOmegAmpl=[deltas',omegas',amplitude']

hatTreconst=zeros(N,K);
for k=1:K
    hatTreconst(:,k)= ContReconst_SIADS(Time(k),Time(1),u,deltas,omegas);
end

Vreconst=U*hatTreconst;

end