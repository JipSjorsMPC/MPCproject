function p = computeRealDiscretePole(h,Re)
% Input: 
%Re:     [ vector containing the real poles]
%h:     sampling period 

%Output:
%y:      discrete pole location in the unit disc
p = zeros(length(Re),1);

    for i = 1:length(Re)
    p(i)= exp(Re(i)*h);
    end

end