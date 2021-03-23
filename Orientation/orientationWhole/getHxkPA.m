function H = getHxkPA(x)

xa = x(3);  ya = x(6); za = x(9); 
q0 = x(10); q1 = x(11); q2 = x(12); q3 = x(13); 
g = 9.7932;

H = NaN(3,1);

% Q = [q0*q0 + q1*q1 - q2*q2 - q3*q3, 2*(q1*q2 + q0*q3), 2*(q1*q3 - q0*q2);
%      2*(q1*q2 - q0*q3), q0*q0 - q1*q1 + q2*q2 - q3*q3, 2*(q2*q3 + q0*q1);
%      2*(q1*q3 + q0*q2), 2*(q2*q3 - q0*q1), q0*q0 - q1*q1 - q2*q2 + q3*q3]
% 
% Q*[xa, ya, za + g]'

HN = NaN(3,1);

HN(1:3) = Q * [xa, ya, za + g]';

H(1,1) = HN(1);
H(2,1) = HN(2);
H(3,1) = HN(3);
% = (2*q1*q3-2*q0*q2)*g;
%  = (2*q0*q1+2*q2*q3)*g;
%  = (2*q0*q0+2*q3*q3-1)*g;

    
end