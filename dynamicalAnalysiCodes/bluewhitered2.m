function c = bluewhitered2 
r = [1 0 0];       %# start
w = [1 1 1];    %# middle
b = [0 0 1];       %# end

%# colormap of size 64-by-3, ranging from red -> white -> blue
c1 = zeros(32,3); c2 = zeros(32,3);
for i=1:3
    c1(:,i) = linspace(b(i), w(i), 32);
    c2(:,i) = linspace(w(i), r(i), 32);
end
c = [c1(1:end-1,:);c2];

% figure
% surf(peaks), shading interp
% caxis([-8 8]), colormap(c), colorbar


% r = [1 0 0];       %# start
% w = [.9 .9 .9];    %# middle
% b = [0 0 1];       %# end
% 
% %# colormap of size 64-by-3, ranging from red -> white -> blue
% c1 = zeros(32,3); c2 = zeros(32,3);
% for i=1:3
%     c1(:,i) = linspace(r(i), w(i), 32);
%     c2(:,i) = linspace(w(i), b(i), 32);
% end
% c = [c1(1:end-1,:);c2];
% 
% surf(peaks), shading interp
% caxis([-8 8]), colormap(c), colorbar