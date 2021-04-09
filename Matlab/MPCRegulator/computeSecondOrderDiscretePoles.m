function p = computeSecondOrderDiscretePoles(zeta,wn,h)

%Input:
%wn:         natural frequency
%alpha:      the factor with which you want to multiply the natural frequency
%h:          sampling period
%zeta:       damping

%Ouput:
%p:          2x1 column vector with corresponding discrete poles in the
%            unit dics

a1 =-2*exp(-zeta*wn*h)*cos(wn*h*sqrt(1-zeta^2));
a2 = exp(-2*zeta*wn*h);
p =roots([1 a1 a2]);

end

