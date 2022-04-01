addpath(genpath(pwd));
data_name = 'MouseKidney_green';
data_dir = ['C:\Users\freds\Desktop\Fourier Ptychography MATLAB\Data\' data_name '.mat']; 
load(data_dir);
for k = 1:size(imlow_HDR,3)
        figure(1);
        set(gcf,'outerposition',get(0,'ScreenSize'))
        imshow(imlow_HDR(:,:,k),[]);
        title(['raw image ' num2str(k)]); pause(0.1);
end