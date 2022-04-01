close all
clear
clc
addpath('Results')
%% blood
edge = 60;
figure;
try
    load('bloodsmear_red_result.mat');
    imBlood = zeros(size(him,1),size(him,1),3);
    subplot(221),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[-1 -1+pi]);title('red');
    imBlood(:,:,1)=abs(him)/max(abs(him(:)))*1.1;
catch
%
end
try
    load('bloodsmear_green_result.mat');
    subplot(222),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('green')
    imBlood(:,:,2)=abs(him)/max(abs(him(:)))*1.1;
catch
%
end
try
    load('bloodsmear_blue_result.mat');
    subplot(223),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('blue');
    imBlood(:,:,3)=abs(him)/max(abs(him(:)))*1.05;
    subplot(224),imshow(imBlood(edge+1:end-edge,edge+1:end-edge,:).^(1/0.75),[]);title('color');
catch
%
end

%% HE
edge = 60;
figure;
try
    load('HE_red_result.mat');
    imHE = zeros(size(him,1),size(him,1),3);
    subplot(221),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[-1 -1+pi]);title('red');
    imHE(:,:,1)=abs(him)/max(abs(him(:)))*1.1;
catch
%
end
try
    load('HE_green_result.mat');
    subplot(222),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('green')
    imHE(:,:,2)=abs(him)/max(abs(him(:)))*1.1;
catch
%
end
try
    load('HE_blue_result.mat');
    subplot(223),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('blue');
    imHE(:,:,3)=abs(him)/max(abs(him(:)))*1.05;
    subplot(224),imshow(imHE(edge+1:end-edge,edge+1:end-edge,:).^(1/0.75),[]);title('color');
catch
%
end

%% IHC
edge = 60;
figure;
try
    load('IHC_red_result.mat');
    imIHC = zeros(size(him,1),size(him,1),3);
    subplot(221),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[-1 -1+pi]);title('red');
    imIHC(:,:,1)=abs(him)/max(abs(him(:)))*1.1;
catch
%
end
try
    load('IHC_green_result.mat');
    subplot(222),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('green')
    imIHC(:,:,2)=abs(him)/max(abs(him(:)))*1.1;
catch
%
end
try
    load('IHC_blue_result.mat');
    subplot(223),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('blue');
    imIHC(:,:,3)=abs(him)/max(abs(him(:)))*1.05;
    subplot(224),imshow(imIHC(edge+1:end-edge,edge+1:end-edge,:).^(1/0.75),[]);title('color');
catch
%
end

%% Misc
edge = 60;
figure;
try
    load('blood_aberration_result.mat');
    imMisc = zeros(size(him,1),size(him,1),3);
    subplot(221),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[-1 -1+pi]);title('red');
    imMisc(:,:,1)=abs(him)/max(abs(him(:)))*1.1;
catch
%
end
try
    load('MouseKidney_green_result.mat');
    subplot(222),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('green')
    imMisc(:,:,2)=abs(him)/max(abs(him(:)))*1.1;
catch
%
end
try
    load('USAF_red_result.mat');
    subplot(223),imshow(angle(him(edge+1:end-edge,edge+1:end-edge)),[]);title('blue');
    imMisc(:,:,3)=abs(him)/max(abs(him(:)))*1.05;

catch
%
end