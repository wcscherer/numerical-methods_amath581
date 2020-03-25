function [SolMat] = Crk_Nln(OpMat,InitVec,vspan,k,dt)
%Implements the backward euler method using the operator matrix OpMat on the
%input initial condition vector InitVec over the span vspan

N = size(OpMat,2);
Ut = [InitVec'];
B = (eye(N,N)-dt*k*OpMat/2);
A = (eye(N,N)+dt*k*OpMat/2);
t = 1;

while t <= size(vspan,2)-1
   
    Umid = (B\(A*Ut(:,end)));
     Ut = horzcat(Ut, Umid);
     t = t + 1;
    
end

SolMat = Ut;
end

