%% PROJECT 1 -- SIFT

%% reading image and converting to sing
I = imread('data1/obj1_5.jpg');
G= rgb2gray(I); %grayscale
ref=single(G); %thre reference image need to be normalized (single format)
figure()
imshow(uint8(ref));

%% Applying SIFT
%thresholds for sift
peakthresh=8; %14.5
edgethresh=1.8; %2

[sift_ref,sift_ref_des] = vl_sift(ref,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);

h1 = vl_plotframe(sift_ref);
set(h1,'color','r','linewidth',3);

n_detected_original=size(sift_ref,2);

%% scaling repeatability
% m=1.2; %scale factor
%  for i=1:9
%      figure()
%      target = imresize(ref,m.^(i-1));
%      imshow(uint8(target));     
%      [sift_tar,sift_tar_des] = vl_sift(target,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);
%      h2 = vl_plotframe(sift_tar);
%      [matches,scores]= vl_ubcmatch(sift_ref_des,sift_tar_des,1.7);
%      %plot(matches(1,:),matches(2,:),'--');
%      repeatiblity_scale(i) = (size(matches,2))./(n_detected_original);
%  end

%% rotation repeatability
angle=15; %rotation step

for i=1:2
     %figure()
     x(i) = (i-1).*angle;
     target = imrotate(ref,x(i));
     imshow(uint8(target));     
     [sift_tar,sift_tar_des] = vl_sift(target,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);
     h2 = vl_plotframe(sift_tar);
     [matches,scores]= vl_ubcmatch(sift_ref_des,sift_tar_des,1.7);
     %plot(matches(1,:),matches(2,:),'--');
     repeatibility_rotation(i) = (size(matches,2))./(n_detected_original);
end
% figure()
% plot(x,repeatibility_rotation,'*');
% title('SIFT Repeatibility rotation')
% xlabel('Angle (degrees)') 
% ylabel('Repeatibility') 
