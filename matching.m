%% PROJECT 1 -- MATCHING
clc
clear all
close all

%% reading image and converting to sing
I = imread('data1/obj1_5.jpg'); %reference image
I2 = imread('data1/obj1_t1.jpg'); %target image
G = rgb2gray(I); %grayscale ref
G2 = rgb2gray(I2); %grayscale tar
ref = single(G); %the reference image need to be normalized (single format)
target = single(G2); %the target as well

figure()
%montage([G, G2]);

imshow([G, G2]);

%% Applying SIFT
%thresholds for sift: we computed them with trial-and-error method to
%achieve few hundreds of keypoints (250), best in term of representation
peakthresh = 8; %8
edgethresh = 1.8; %1.8

[sift_ref,sift_ref_desc] = vl_sift(ref,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);

h1 = vl_plotframe(sift_ref);
set(h1,'color','r','linewidth',3);

[sift_tar,sift_tar_desc] = vl_sift(target,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);


h2 = vl_plotframe(sift_tar);
set(h2,'color','y','linewidth',3);

%% 'fixed threshold' method
% for 
%     for
%         euclideanDist = sqrt(sum((sift_ref_desc(:,1) - sift_ref_desc(:,2)).^2));
%         
%         euclideanDist2 = sqrt(sum((sift_ref_desc(:,1) - sift_ref_desc(:,1)).^2));
%     end
% end



