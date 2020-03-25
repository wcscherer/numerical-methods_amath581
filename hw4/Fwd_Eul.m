function [SolMat] = Fwd_Eul(OpMat,InitVec,vspan,k,dt)
%Implements the forward euler method using the operator matrix OpMat on the
%input initial condition vector InitVec over the span vspan
%   Detailed explanation goes here

Ut = [InitVec'];
t = 1;

while t <= size(vspan,2)-1

     Umid = Ut(:,end)+(k*dt)*OpMat*Ut(:,end);
     Ut = horzcat(Ut, Umid);
     t = t + 1;
    
    
end

SolMat = Ut;
end

