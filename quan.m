function [y]=quan(x, N, non_quan_bits,type)


if type == "mix quan"
    y = 1/sqrt(2)*(sign(real(x))+1i*sign(imag(x)));
    for non_quan_bit=non_quan_bits
        for j = 1:N
            y(non_quan_bit,j) = x(non_quan_bit,j);
        end
    end

elseif type == "all quan"
    y = 1/sqrt(2)*(sign(real(x))+1i*sign(imag(x)));


elseif type == "no quan"
    y = x;

else
    error("Incorrect quantization type!")
end
end

