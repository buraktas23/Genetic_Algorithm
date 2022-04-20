function [kr] = kromozoms()
%UNTÝTLED2 Summary of this function goes here
%   Detailed explanation goes here


populations = unifrnd(1,100,1,10);


kr = zeros(6,10);

for i=1:6
   temp = randperm(10);
   
    for j=1:10
        kromozom(j)= populations(temp(j));        
    end
    kr(i,:) = kromozom;
end
    
end

