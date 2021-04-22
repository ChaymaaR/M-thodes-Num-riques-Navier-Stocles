clear all
clc
L=input('Longueur: ');
H=input('Largeur: ');
T=input('Temps: ');
I=input('Valeur de I: ');
J=input('Valeur de J: ');
dt=input('Pas de temps: ');
Re=input('Reynolds: ');
Up=input('Vitesse paroie supérieur: ');
dx=L./(I+1);
dy=L./(J+1);
%% DD2x
v=ones(I-1,1);
Id=eye(I,I);
D2x=(-2.*Id+diag(v,1)+diag(v,-1));
DD2x=zeros(I*J,I*J);
for j=1:J
Z=(j-1)*I+1:j*I;
DD2x(Z,Z)=D2x;
end
DD2x=(1./dx.^2).*DD2x
%% DD2y
Idd=eye(J,J);
w=ones(I*J-J,1);
D2y=-2.*Idd;
DD2y=zeros(I*J,I*J);
for i=1:I
ZZ=(i-1)*J+1:J*i;
DD2y(ZZ,ZZ)=D2y;
end
DD2y=DD2y+diag(w,J)+diag(w,-J);
DD2y=(1./dy.^2).*DD2y
%%A
A=(I./dt)-(1./Re).*(DD2x+DD2y)
%%Gxx
Gx=zeros(I,I+1);
Gxx=zeros(I*J,(J+1)*(I+1));
for i=1:I
Gx(i,i)=-1;
Gx(i,i+1)=1;
end
Gx;
for j=1:J
Z1=(j-1)*I+1:j*I;
Z2=(j-1)*(I+1)+1:j*(I+1);
Z3=j*(I+1)+1:(j+1)*(I+1);
Gxx(Z1,Z2)=Gx;
G2xx(Z1,Z3)=Gx;
end
Gxx=(1/(2*dx))*(Gxx+G2xx)
%%Gy
Gy=zeros(I,I+1);
Gyy=zeros(I*J,(J+1)*(I+1));
for i=1:I
Gy(i,i)=-1;
Gy(i,i+1)=1;
end
Gy;
for j=1:J
Z1=(j-1)*I+1:j*I;
Z2=(j-1)*(I+1)+1:j*(I+1);
Z3=j*(I+1)+1:(j+1)*(I+1);
Gyy(Z1,Z2)=Gx;
G2yy(Z1,Z3)=Gx;
end
Gyy=(1/(2*dy))*(Gyy+G2yy)
%%Dx
for i=1:I
Dx1(i,i)=1;
Dx1(i+1,i)=-1;
end
for j=1:J
Z1=(j-1)*(I+1);
Z2=(j-1)*I;
ZZ1=Z1+1:Z1+I+1;
ZZ2=Z2+1:Z2+I;
ZZ3=Z1+I+2:Z1+2*I+2;
ZZ4=Z2+1:Z2+I;
Dx(ZZ1,ZZ2)=Dx1;
Dx(ZZ3,ZZ4)=Dx1;
end
Dx=Dx./(2*dx)
%Dy
for i=1:I
Dy1(i,i)=1;
Dy1(i+1,i)=1;
end
for j=1:J
Z1=(j-1)*(I+1);
Z2=(j-1)*I;
ZZ1=Z1+1:Z1+I+1;
ZZ2=Z2+1:Z2+I;
ZZ3=Z1+I+2:Z1+2*I+2;
ZZ4=Z2+1:Z2+I;
Dy(ZZ1,ZZ2)=Dy1;
Dy(ZZ3,ZZ4)=-Dy1;
end
Dy=Dy./(2*dy)
%Bxx
U=zeros(I,J);
V=zeros(I,J);
Bx=zeros(I,J);
for i=2:I-1
j=1;
Bx(i,j)=(U(i,j)./dt)-(U(i,j)*(U(i+1,j)-U(i-1,j)))./(2*dx)-(U(i,j)*(U(i,j+1)))./(2*dy);
end
for i=2:I-1
j=J;
Bx(i,j)=(U(i,j)./dt)-(U(i,j)*(U(i+1,j)-U(i-1,j)))./(2*dx)-(U(i,j)*(Up-U(i,j-1)))./(2*dy);
end
for j=2:J-1
i=1;
Bx(i,j)=(U(i,j)./dt)-(U(i,j)*(U(i+1,j)))./(2*dx)-(U(i,j)*(U(i,j+1)-U(i,j-1)))./(2*dy);
end
for j=2:J-1
i=I;
Bx(i,j)=(U(i,j)./dt)-(U(i,j)*(-U(i-1,j)))./(2*dx)-(U(i,j)*(U(i,j+1)-U(i,j-1)))./(2*dy);
end
for j=1;
i=1;
Bx(i,j)=(U(i,j)./dt)-(U(i,j)*(U(i+1,j)))./(2*dx)-(V(i,j)*(U(i,j+1)))./(2*dy);
end
for j=J;
i=I;
Bx(i,j)=(U(i,j)./dt)-(U(i,j)*(-U(i-1,j)))./(2*dx)-(U(i,j)*(Up-U(i,j-1)))./(2*dy);
end
for j=2:J-1
for i=2:I-1
Bx(i,j)=(U(i,j)./dt)-(U(i,j)*(U(i+1,j)-U(i-1,j)))./(2*dx)-(U(i,j)*(U(i,j+1)-U(i,j-1)))./(2*dy);
end
end
Bxx=zeros(I*J,1);
for j=1:J
for i=1:I
Bxx(i*i,1)=Bx(i,j)
end
end
Bx
%Byy
U=zeros(I,J);
V=zeros(I,J);
By=zeros(I,J);
for i=2:I-1
j=1;
By(i,j)=(V(i,j)./dt)-(V(i,j)*(V(i+1,j)-V(i-1,j)))./(2*dx)-(V(i,j)*(V(i,j+1)))./(2*dy);
end
for i=2:I-1
j=J;
By(i,j)=(V(i,j)./dt)-(V(i,j)*(V(i+1,j)-V(i-1,j)))./(2*dx)-(V(i,j)*(-V(i,j-1)))./(2*dy);
end
for j=2:J-1
i=1;
By(i,j)=(V(i,j)./dt)-(V(i,j)*(V(i+1,j)))./(2*dx)-(V(i,j)*(V(i,j+1)-V(i,j-1)))./(2*dy);
end
for j=2:J-1
i=I;
By(i,j)=(V(i,j)./dt)-(V(i,j)*(-V(i-1,j)))./(2*dx)-(V(i,j)*(V(i,j+1)-V(i,j-1)))./(2*dy);
end
for j=1;
i=1;
By(i,j)=(V(i,j)./dt)-(V(i,j)*(V(i+1,j)))./(2*dx)-(V(i,j)*(V(i,j+1)))./(2*dy);
end
for j=J;
i=I;
By(i,j)=(V(i,j)./dt)-(V(i,j)*(-V(i-1,j)))./(2*dx)-(V(i,j)*(-V(i,j-1)))./(2*dy);
end
for j=2:J-1
for i=2:I-1
By(i,j)=(V(i,j)./dt)-(V(i,j)*(V(i+1,j)-V(i-1,j)))./(2*dx)-(V(i,j)*(V(i,j+1)-V(i,j-1)))./(2*dy);
end
end
Byy=zeros(I*J,1);
for j=1:J
for i=1:I
Byy(i*i,1)=By(i,j)
end
end
Byy
%Matrice final
Matrice_final=[A,0*eye(I*J,I*J),Gxx;0*eye(I*J,I*J),A,Gyy;Dx,Dy,(-0.0000001)*eye((J+1)*(I+1),(J+1)*(I+1))]
%C
C=0;
%P
U=zeros((I+1)*(J+1),1);
P=(Dx.*inv(A).*Bxx+Dy.*inv(A).*Byy-C).*inv(Dx.*inv(A).*Gyy+Dy.*inv(A).*Gyy)
%U
U=zeros(I*J,1);
U=A^(-1).*(Bxx-Gxx.*P)
%V
V=zeros(I*J,1);
V=A^(-1).*(Byy-Gyy.*P)