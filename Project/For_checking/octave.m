
N = 11;


a = 0;
b = 2*pi;
x = linspace(a,b,N);

h = (b-a)/(length(x)-1);


epsilon = h/2;
T_final = 2;
dt = 0.01*h;


t_steps =round(T_final/dt);




function [ans] = CG0(A,b,tol = 0.0001)
    n=size(b,1);
    x = zeros(n,1);    
    r = b - A*x;
    p = r; 
    
    for k = 1:n      
      w = A*p;
   
      alpha =(r'*r)/(p'*w);
      
      
      x = x + alpha*p; 
           
      r_old=r;
     
     
      r = r - alpha*w; 
     
      if( norm(r) < tol )
        break;
      endif      
      B = (r'*r)/(r_old'*r_old);
      
      p = r + B*p; 
     
      
    endfor
    ans=x;
    itr=k;
endfunction






function [M] = Mass(x,h)
  n = length(x) - 1;
  M = zeros(n+1,n+1);
  for i = 1:n
    M(i,i) = M(i,i) + h/3;
    M(i,i+1) = M(i,i+1) + h/6;
    M(i+1,i) = M(i+1,i) + h/6;
    M(i+1,i+1) = M(i+1,i+1) + h/3;
    
  endfor
  M(1,1) = M(1,1) + h/3;
  M(n+1,n+1) = M(n+1,n+1) +h/3;
endfunction


function [S] = Diffusion(x,h)
  n = length(x) - 1;
  S = zeros(n+1,n+1);
  for i = 1:n
    S(i,i) = S(i,i) + 1/h;
    S(i,i+1) = S(i,i+1) - 1/h;
    S(i+1,i) = S(i+1,i) - 1/h;
    S(i+1,i+1) = S(i+1,i+1) + 1/h;
  endfor
  S(1,1) = S(1,1) + 1/h;
  S(n+1,n+1) = S(n+1,n+1) +1/h;
endfunction



function [A] = Advection(x,h)
  n = length(x) - 1;
  A = zeros(n+1,n+1);
  for i = 1:n
    A(i,i) = A(i,i) + 0;
    A(i,i+1) = A(i,i+1) + 1/2;
    A(i+1,i) = A(i+1,i) -1/2 ;
    A(i+1,i+1) = A(i+1,i+1) +0;
    
  endfor
  
endfunction


#initial_data = sin(x);
function [u0] = ID(x)
  u0 = sin(x);
endfunction




function F = f1(u,A,S,epsilon,t,h)

  epsilon;
  

  beta = 1/2 * u.^2;


  F = -A*beta' - epsilon*S*u' ;
endfunction
initial_data = ID(x);
U=zeros(2,length(x));
U(1,:) = u0 = initial_data;
M = Mass(x(2:N-1),h);
xx = x(2:N-1);
nn =length(xx); 
A = Advection(xx,h);
S = Diffusion(xx,h);

for k = 2:t_steps+1 
  t = dt*(k-1)
  display(k);
  

  
  k1 =   CG0(M, f1(u0(2:N-1)  ,A, S ,epsilon,t,h))*dt ;
  FF = f1(u0(2:N-1)  ,A, S ,epsilon,t,h);
  k2 =   CG0(M, f1(u0(2:N-1)+k1'/2  ,A ,S ,epsilon,t,h))*dt ;
  k3 =   CG0(M, f1(u0(2:N-1)+k2'/2  ,A ,S ,epsilon,t,h))*dt ; 
    
  k4 =   CG0(M, f1(u0(2:N-1)+k3'    ,A ,S ,epsilon,t,h))*dt ; 
  
  
  U(k,2:N-1)     = u0(2:N-1)+(k1'+2*k2'+2*k3'+k4')/6;
   
  U(k,1)         = 0;
  U(k,length(x)) = 0;
  u0             = U(k,:);
  
  
 
endfor

U;

%{
figure(1)
plot(x,U(t_steps+1,:),'c','LineWidth',1.5) 
hold on

xlabel('x')
ylabel('u(x,t=2)')
xlim([0 2*pi])
legend('N = 41','N = 81','N = 161','N = 321','N = 641');
#title('(a)')
set(gca,'FontSize',15)
%}
#0 0.001
#b r