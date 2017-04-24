
clear
Nt=1;
ground=zeros(1,Nt);
eigvec=zeros(500,5);
for k=1:Nt
%% claim the variables
ymax=100;               %maximum of y
N=1000;                %number of points of y
N1=100;                 %number of points of phi
N2=50;                 %number of points of finer grid
A=zeros(N,N);           %Matrix for eigenvalue problem
a=zeros(1,N1);          %1:N
b=zeros(1,N1);          %discretization of phi
y=zeros(1,N);           %discretization of y
d=zeros(1,N);           %matrix of Coulumb potential
dN=ymax/N;
w=zeros(1,N1);          %weight
for i=1:N
    y(i)=dN*(i-1/2);
end
y1=y;
a=1:N1;

b=pi/(N1+1)*(a-(N1+1)/(2*pi) *sin((2*pi*a)/(N1+1)));
w=pi/(N1+1)*(1-cos(2*pi*a/(N1+1)));

d=zeros(N,N);
 for j=1:N
     %% old wrong
%  %      y2=y1(j)-dN/2+(1:N1)*dN/N1;
%     for i=1:N1
%         d(j,:)=d(j,:)+(1./(sqrt(y(j)^2+y1.^2-2*y1*y(j)*b(i))))*pi/N1*2.*y1;
% %         d(j,j)=d(j,j)+(1/(sqrt(y(j)^2+y2(i)^2-2*y2(i)*y(j)*b(i))))*pi/N1^2*2*y2(i);
%     end
%     
%%  finer grid for removing singularity
    y2=y1(j)-dN/2+((1:N2)-1/2)*dN/N2;
%%  calculate the nondiagonal terms
    for i=1:N1   %calculate the matrix as if singularity is not removed, to get nondiagonal terms
        d(j,:)=d(j,:)+(1./(sqrt(y(j)^2+y1.^2-2*y1*y(j)*cos(b(i)))) )*2.*y1*w(i);  %  integral before singularity removal
    end
    d(j,j)=0; %remove diagonal terms
%%  calculate diagonal terms
    for i1=1:N2
        d(j,j)=d(j,j)+sum((1./(sqrt(y(j)^2+y2(i1)^2-2*y2(i1)*y(j)*cos(b))))*2*y2(i1).*w)/N2;  
    end
 end
A=d;
A=-1*A*dN/pi;
for i=1:N
    A(i,i)=A(i,i)+y(i)^2;
end
eigen=eig(A);
result=sort(eigen(find(eigen<0))); %#ok<*FNDSB>
ground(k)=result(1);
[V,D]=eig(A);

% Vcompose=V(:,find(eigen<0));
% eigvec(1:N,k)=abs(Vcompose(:,1));
% hold on 
% plot(y,abs(Vcompose(:,1)))

end
% legend('100','200','300','400','500')
xlabel('y')
figure()
plot(100*(1:Nt),ground)
xlabel('Number of points')
ylabel('Ground states')


% gs(k)=result(1);
% plot(gs)

% st=-4*ones(1,k);
% figure(1);
% plot(st);
% hold on
% plot(ground)
% set(gca, 'xticklabel',N*(1:k)*dN)

% compose=sqrt(1./((-1)*result))-0.5;
[V,D]=eig(A);
Vcompose=V(:,find(eigen<0));
% S=size(result);
% x=1:S;
% figure(1);
% hold on
% te=['o','+','*','x','s','d'];
% t=3;
% scatter(x-1, log(-result),te(t))
% figure(2);
% hold on
% scatter(x-1, result, te(t));
% err=zeros(1,100);
% for i=1:100
%     err(i)=(A(i,i)-arr(i,i))/A(i,i)*100;
% end