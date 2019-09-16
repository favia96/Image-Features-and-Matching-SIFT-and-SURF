%% PROJECT 1 -- MATCHING
clc
clear all
close all

%% reading image and converting to sing
I = imread('data1/obj1_5.jpg'); %reference image
I2 = imread('data1/obj1_t1.jpg'); %target image
%I2 = imresize(I2,[size(I,1) size(I,2)]); %only if 2nd image has different dimension
G = rgb2gray(I); %grayscale ref
G2 = rgb2gray(I2); %grayscale tar
ref = single(G); %the reference image need to be normalized (single format)
target = single(G2); %the target as well

figure1 = figure()

imshow([G, G2]);

%% PART 1 - SIFT
%% Applying SIFT on reference image
%thresholds for sift: we computed them with trial-and-error method to
%achieve few hundreds of keypoints (250), best in term of representation
peakthresh = 8; %8
edgethresh = 1.8; %1.8

[sift_ref,sift_ref_desc] = vl_sift(ref,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);

n_det_ref = size(sift_ref,2); %number of detected keypoints on reference image

%plotting sift keypoints of reference image
h1 = vl_plotframe(sift_ref);
set(h1,'color','r','linewidth',3);

%% Applying SIFT on target image 
[sift_tar,sift_tar_desc] = vl_sift(target,'PeakThresh',peakthresh,'EdgeThresh',edgethresh);

n_det_tar = size(sift_tar,2); %number of detected keypoints on target image

%plotting the shifted sift keypoints of the target image (because it's a montage of himage)
sift_tar_plot = sift_tar;
sift_tar_plot(1,:) = sift_tar(1,:) + size(target,2); %with the shift
h2 = vl_plotframe(sift_tar_plot);
set(h2,'color','y','linewidth',3);
hold on;

%% 'fixed threshold' matching algorithm
%euclidean distance on feature (descriptor) space below a fixed thresh
%fixed_thresh_euc = 69; %trial-error (min=50.52, mean=99.04, median=98.39, max=141.57) if euclidean
%fixed_thresh_chi = 35;   %trial-error (min=17.5, mean=128.4, median=123.1, max=354) if chi-square
%counter_matches = 0;

 for i = 1 : n_det_ref
        for j = 1 : n_det_tar
           euclid_dist_sift(i,j) = sqrt(sum((sift_ref_desc(:,i) - sift_tar_desc(:,j)).^2));
%           %chi_square_dist_sift(i,j) = sum(((sift_ref_desc(:,i) - sift_tar_desc(:,j)).^2 ) ./ (sift_ref_desc(:,i) + sift_tar_desc(:,j) )) ./ 2;
%           %if(euclid_dist_sift(i,j) <= fixed_thresh_euc) 
%           %if(chi_square_dist_sift(i,j) <= fixed_thresh_chi) %in case of chi-square
%                 counter_matches = counter_matches + 1; 
%                 match_pairs_indexes(counter_matches,1) = i; %index of ref
%                 match_pairs_indexes(counter_matches,2) = j; %index of target                
%           end          
       end
 end

% for k = 1 : counter_matches
%     index_ref = match_pairs_indexes(k,1);
%     index_tar = match_pairs_indexes(k,2);
%     plot([sift_ref(1,index_ref)  sift_tar_plot(1,index_tar)],[sift_ref(2,index_ref) sift_tar_plot(2,index_tar)],'-b');
% end
% 
% title('SIFT Match: Fixed Threshold');
% saveas(figure1,'match_sift_fixed_thresh.png');

% %% 'nearest neighbour' matching algorithm
% %matching with the nearest neighbour in terms of euclidean distance
%  
% for l = 1 : n_det_ref
%     [m,n] = min(euclid_dist_sift(l,:));  
%     match_pairs_indexes(l,1) = l; %index of ref
%     match_pairs_indexes(l,2) = n; %index of target
% end
% 
% for k = 1 : size(match_pairs_indexes,1)
%    index_ref = match_pairs_indexes(k,1);
%    index_tar = match_pairs_indexes(k,2);
%    plot([sift_ref(1,index_ref)  sift_tar_plot(1,index_tar)],[sift_ref(2,index_ref) sift_tar_plot(2,index_tar)],'-b');
% end
% 
% title('SIFT Match: Nearest Neighbour');
% saveas(figure1,'match_sift_nearest_neighbour.png');

% %% 'nearest neighbour distance ratio' matching algorithm
% %matching with the nearest neighbour in terms of ratio of euclidean
%distance with 2nd nearest under a thresh
ratio_thresh = 0.89; %in literature the best is 0.8

for i = 1 : n_det_ref
    [m,n] = mink(euclid_dist_sift(i,:),2); %smallest 2
    if(m(1) / m(2)) <= ratio_thresh
        match_pairs_indexes(i,1) = i; %index of ref
        match_pairs_indexes(i,2) = n(1); %index of target
    end       
end

for k = 1 : size(match_pairs_indexes,1)
        index_ref = match_pairs_indexes(k,1);
        index_tar = match_pairs_indexes(k,2);
        if (index_ref~=0)
            plot([sift_ref(1,index_ref)  sift_tar_plot(1,index_tar)],[sift_ref(2,index_ref) sift_tar_plot(2,index_tar)],'-b');
        end  
end

title('SIFT Match: Nearest Neighbour Distance Ratio');
%saveas(figure1,'match_sift_nearest_neighbour_dist_ratio.png');


% %% PART 2 - SURF
% %% Applying SURF on reference image
% ref = G; %surf need only grayscale not single 
% target = G2;
% 
% surf_ref = detectSURFFeatures(ref);
% surf_ref = surf_ref.selectStrongest(250); %only 250 (few hundreds) strongest keypoints on ref
% [surf_ref_desc,vpts1] = extractFeatures(G,surf_ref);
% n_det_ref = surf_ref.Count; %number of detected keypoints on reference image
% surf_ref_desc = surf_ref_desc'; %transposed matrix
% 
% hold on
% surf_ref.plot; %plot surf keypoints on ref
% 
% % %% Applying SURF on target image 
% surf_tar = detectSURFFeatures(target);
% surf_tar = surf_tar.selectStrongest(250); %only 250 (few hundreds) strongest keypoints on target
% [surf_tar_desc,vpts2] = extractFeatures(G2,surf_tar);
% n_det_tar = surf_tar.Count; %number of detected keypoints on target image
% surf_tar_desc = surf_tar_desc'; %transposed matrix
% 
% % %plotting the shifted sift keypoints of the target image (because it's a montage of himage)
% surf_tar_plot = surf_tar;
% surf_tar_plot.Location(:,1) = surf_tar.Location(:,1) + size(target,2); %with the shift
% surf_tar_plot.plot; %plot surf keypoints on target
% 
% % hold on;
% 
% % %% 'nearest neighbour distance ratio' matching algorithm
% % %matching with the nearest neighbour in terms of ratio of euclidean
% % %distance with 2nd nearest under a thresh
% 
% for i = 1 : n_det_ref
%        for j = 1 : n_det_tar
%           euclid_dist_surf(i,j) = sqrt(sum((surf_ref_desc(:,i) - surf_tar_desc(:,j)).^2));
%           %chi_square_dist_surf(i,j) = sum(((surf_ref_desc(:,i) - surf_tar_desc(:,j)).^2 ) ./ (sift_ref_desc(:,i) + sift_tar_desc(:,j) )) ./ 2;     
%       end
% end
% 
% ratio_thresh = 0.8; %in literature the best is 0.8
%  
% for i = 1 : n_det_ref
%      [m,n] = mink(euclid_dist_surf(i,:),2); %smallest 2
%      if(m(1) / m(2)) <= ratio_thresh
%          match_pairs_indexes(i,1) = i; %index of ref
%          match_pairs_indexes(i,2) = n(1); %index of target
%      end       
% end
% 
% hold on;
% 
% for k = 1 : size(match_pairs_indexes,1)
%          index_ref = match_pairs_indexes(k,1);
%          index_tar = match_pairs_indexes(k,2);
%          if (index_ref~=0)
%              plot([surf_ref.Location(index_ref,1)  surf_tar_plot.Location(index_tar,1)],[surf_ref.Location(index_ref,2) surf_tar_plot.Location(index_tar,2)],'-b');
%          end  
% end
% title('SURF Match: Nearest Neighbour Distance Ratio');
% saveas(figure1,'match_surf_nearest_neighbour_dist_ratio.png');
