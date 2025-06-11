function [angle] = music(rx,k)

% figure
% pmusic(rx, k)
[SP1,f]=pmusic(rx, k);

[pks, locs1] = findpeaks(SP1, f);

if length(locs1) >= k
    [~, sortIdx] = sort(pks, 'descend');
    locs1 = locs1(sortIdx(1:k));
else
    if isempty(locs1)
        locs1 = 0;
    else
        figure
        pmusic(rx, k)
        error('Unable to find %d valid peaks, only found %d', k, length(locs1));
    end
end

for i=1:length(locs1)
    if locs1(i) < pi
        locs1(i) = locs1(i);
    else
        locs1(i) = locs1(i) - 2*pi;
    end
end

rad_angle = locs1/pi;

angle = real(asin(rad_angle)) * 180 / pi;

angle = sort(angle, "ascend");

angle = angle.';

end