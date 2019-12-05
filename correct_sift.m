%% PROJECT 1 Analysis and Search of Visual Data - SIFT
%% Federico Favia, Mayank Gulati, September 2019

%clc
%clear all
%close all

%% reading image and converting to sing
I = imread('data1/obj1_5.jpg'); %reference image
G = rgb2gray(I); %grayscale
ref = single(G); %the reference image need to be normalized (single format)

figure2=figure()
imshow(uint8(ref));

%% Applying SIFT
%thresholds for sift: we computed them with trial-and-error method to
%achieve few hundreds of keypoints (250), best in term of representation
peakthresh = 8; %8
edgethresh = 1.8; %1.8

[sift_ref,sift_ref_desc] = vl_sift(ref,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);

h1 = vl_plotframe(sift_ref);
set(h1,'color','r','linewidth',3);

n_det_ref = size(sift_ref,2); %number of detected keypoints on reference image

title('SIFT best 250 keypoints, Peak Thresh=8, Edge Thresh=1.8')
saveas(figure2,'siftfkeypoints.png');

% %% scaling repeatibility
% m = 1.2; %scale factor
% counter = 0; %numerator of repeatibility (number of matches)
% 
% for i = 1 : 9 %nine scale factors
%     
%       scale(i) = m.^(i-1);
% 
%       sift_ref_scaled = scale(i) .* sift_ref(1:2,:); %scale the SIFT keypoints of the ref to predict 
%  
%       target = imresize(ref,m.^(i-1)); %scale the reference image
%       %imshow(uint8(target));     
%       [sift_tar,sift_tar_desc] = vl_sift(target,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);
%       %h2 = vl_plotframe(sift_tar);
%       n_det_tar = size(sift_tar,2); %number of detected keypoints on target
%       %hold on
%       %plot(sift_ref_scaled(1,:),sift_ref_scaled(2,:),'O');
%       
%       for j = 1 : n_det_ref
%          for k = 1 : n_det_tar
%                  if (abs(sift_ref_scaled(1,j) - sift_tar(1,k)) <= 2) && (abs(sift_ref_scaled(2,j) - sift_tar(2,k)) <= 2)
%                      counter = counter + 1; 
%                      break;
%                  end
%          end
%      end 
%       
%      repeatibility_scale(i) = counter ./ n_det_ref; %repeatibility formula
%      counter = 0; %reset counter of matches    
% end
% 
% % plot repeatibility over scaling
% figure1 = figure()
% plot(scale,repeatibility_scale,'-*b');
% hold on
% plot(scale, repeatibility_scale_surf,'-Or'); %taken from the other script
% title('SIFT vs SURF Repeatibility wrt scales')
% xlabel('Scale') 
% ylabel('Repeatibility(scale)') 
% legend('SIFT','SURF')
% saveas(figure1,'siftvssurf_repeat_scale.png');


% %% rotation repeatibility
% angle = 15; %rotation angle step (degree)
% counter = 0; %numerator of repeatibility (number of matches) 
% 
% for i = 1 : 25 %number of rotations
%     theta(i) = (i-1).*angle;
%     
%     rotation_matrix = [cosd(theta(i)) -sind(theta(i)); sind(theta(i)) cosd(theta(i))]; %rotation matrix
%     
%     %shift
%     sift_ref_shift(1,:) = sift_ref(1,:) - size(G,2)/2;
%     sift_ref_shift(2,:) = -( sift_ref(2,:) - size(G,1)/2 );
%     
%     sift_ref_rotated = rotation_matrix * sift_ref_shift(1:2,:); %rotate the SIFT keypoints of the ref to predict 
% 
%     target = imrotate(ref,theta(i),'bicubic'); %rotate the ref image
%     
%     %shifting back
%     sift_ref_rotated(1,:) = sift_ref_rotated(1,:) + size(target,2)/2;
%     sift_ref_rotated(2,:) = - (sift_ref_rotated(2,:) - size(target,1)/2 );
%     
%     %figure()
%     %imshow(uint8(target));     
%     [sift_tar,sift_tar_desc] = vl_sift(target,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);
%     n_det_tar = size(sift_tar,2); %number of detected keypoints on target
%     %h2 = vl_plotframe(sift_tar);
%     %hold on
%     %plot(sift_ref_rotated(1,:),sift_ref_rotated(2,:),'O')
%     
%     for j = 1 : n_det_ref
%         for k = 1 : n_det_tar
%                 if (abs(sift_ref_rotated(1,j) - sift_tar(1,k)) <= 2) && (abs(sift_ref_rotated(2,j) - sift_tar(2,k)) <= 2)
%                     counter = counter + 1; 
%                     break;
%                 end
%         end
%     end 
%      
%     repeatibility_rotation(i) = counter ./ n_det_ref; %repeatibility formula
%     counter = 0; %reset counter of matches
% end
% 
% % figure()
% % plot(sift_ref_rotated(1,:),sift_ref_rotated(2,:),'or'); %plot rotated sift of ref
% % hold on
% % plot(sift_tar(1,:),sift_tar(2,:),'*k'); %plot sift on rotated
% 
% %plot repeatibility over rotation angles
% figure2 = figure()
% plot(theta,repeatibility_rotation,'-*b');
% hold on;
% plot(theta,repeatibility_rotation_surf,'-Or'); %taken from surf script
% title('SIFT vs SURF Repeatibility wrt rotation')
% xlabel('Theta (deg)') 
% ylabel('Repeatibility(theta)') 
% legend('SIFT','SURF')
% saveas(figure2,'siftvssurf_repeat_rotation.png')
