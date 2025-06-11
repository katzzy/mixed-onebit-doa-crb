function [rx]=getrx(y,M1,M2,N,D,non_quan_bits,type,is_nested)

% Calculate the covariance matrix Ry
Ry = 1/N * (y * y');

if type == "mix quan"

    % Restore the pre-quantization covariance matrix
    % Quantization part
    Rx = sin(pi/2.*real(Ry)) + 1i * sin(pi/2.*imag(Ry));

    % Non-quantitative part
    for non_quan_bit=non_quan_bits
        for i = 1:M1+M2
            if ismember(i, non_quan_bits)
                Rx(i, non_quan_bit) = Ry(i,non_quan_bit);
                Rx(non_quan_bit, i) = Ry(non_quan_bit,i);
            else
                Rx(i, non_quan_bit) = sqrt(2.*pi)/2*real(Ry(i,non_quan_bit)) + 1i*sqrt(2.*pi)/2*imag(Ry(i,non_quan_bit));
                Rx(non_quan_bit, i) = sqrt(2.*pi)/2*real(Ry(non_quan_bit,i)) + 1i*sqrt(2.*pi)/2*imag(Ry(non_quan_bit,i));
            end
        end
    end
elseif type == "all quan"
    Rx = sin(pi/2.*real(Ry)) + 1i * sin(pi/2.*imag(Ry));
elseif type == "no quan"
    Rx = Ry;
else
    error("错误的量化类型");
end
if is_nested
    % Calculate the rx of the virtual array element
    rx = calc_virtual_rx(Rx, M1, M2, D);
else
    rx = Rx;
end
end



% Calculate the rx of the virtual array element
function [rx]=calc_virtual_rx(Rx, M1, M2, D)

% Calculate array position
d = zeros(1, M1 + M2);
for i = 1:M1+M2
    if i <= M1
        d(i) = i - 1;
    else
        d(i) = M1 + (i - M1 - 1)*(M1 + 1);
    end
end

N = zeros(D, 1);
rx = zeros(D,1);
for i = 1:M1+M2
    for j = 1:M1+M2
        temp = d(i)-d(j) + d(end) + 1;
        N(temp) = N(temp) + 1;
        rx(temp) = Rx(i,j) + rx(temp);
    end
end
rx = rx ./ N;

end