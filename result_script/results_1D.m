%output_2d_s is strain image for 1d calcualted with 73 window.

%SNR
window=output_2d_s(565:658,151:193);
spatial_avg=mean(window(:));
spatial_std=std(window(:));
snr=spatial_avg/spatial_std


%CNR
window_bg=window;
window_tg=output_2d_s(853:928,185:215);
mean_tg=mean(window_tg(:));
std_tg=std(window_tg(:));
mean_bg=spatial_avg;
std_bg=spatial_std;
cnr=sqrt((2*(mean_bg-mean_tg)^2)/(std_tg^2+std_bg^2))


%strain ratio
strain_ratio=mean_tg/mean_bg

%strain plot
%column
% for i=200:230
% figure,plot(output_2d_s(:,i)), title(i)
% end
figure,plot(output_2d_s(:,211))

% for i=970:1000
% figure,plot(output_2d_s(i,:)), title(i)
% end
