function steering_vector = steeringVector(M,theta,d,lambda)

m   = (0:M-1).';
tau = d*sin(theta)/lambda;
Tau = repmat(tau',M,1);

steering_vector = exp(-2*pi*1i*Tau.*repmat(m,1,size(theta,1)));

end