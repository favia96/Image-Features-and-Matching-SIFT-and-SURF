

%% PROJECT 1 -- SURF

%% reading image and converting 
I = imread('data1/obj1_5.jpg');
G= rgb2gray(I); %grayscale

figure()
imshow(uint8(ref));

%% Applying SURF

surf_ref = detectSURFFeatures(G);
surf_ref_best = surf_ref.selectStrongest(250);
%plot(surf.selectStrongest(250));

n_detected_original=133;

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

for i=1:25
     x(i) = (i-1).*angle;
     target = imrotate(G,x(i));

     surf_tar = detectSURFFeatures(target);
     surf_tar_best = surf_tar.selectStrongest(250);
     
     [f1,vpts1] = extractFeatures(G,surf_ref_best);
     [f2,vpts2] = extractFeatures(target,surf_tar_best);
     indexPairs = matchFeatures(f1,f2);
     
     repeatibility_rotation(i) = (size(indexPairs,1))./(n_detected_original);
end

%% plot repeatibility
figure()
plot(x,repeatibility_rotation);
title('SURF Repeatibility rotation')
xlabel('Angle (degrees)') 
ylabel('Repeatibility') 







