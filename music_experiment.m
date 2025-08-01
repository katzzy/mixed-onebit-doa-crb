function [angle] = music_experiment(M1, M2, k, N, snr, angle_set, D, type, non_quan_bits)

angle_pi = angle_set./180.*pi;

x = xcreate(M1,M2,k,N,snr,angle_pi);

if type == "mix quan"
    y_mix_quan = quan(x, N, non_quan_bits,"mix quan");
    rx = getrx(y_mix_quan,M1,M2,N,D,non_quan_bits,"mix quan",true);

elseif type == "all quan"
    y_all_quan = quan(x, N, non_quan_bits,"all quan");
    rx = getrx(y_all_quan,M1,M2,N,D,non_quan_bits,"all quan",true);

elseif type == "no quan"
    y = quan(x, N, non_quan_bits,"no quan");
    rx = getrx(y,M1,M2,N,D,non_quan_bits,"no quan",true);
end


rx  = flip(rx);
angle = music(rx,k);


end