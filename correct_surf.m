%% PROJECT 1 -- SURF CORRECTED

%% reading image and converting to grayscale
I = imread('data1/obj1_5.jpg'); 
ref = rgb2gray(I); %grayscale

figure()
imshow(ref);

%% Applying SURF
surf_ref = detectSURFFeatures(ref);
surf_ref = surf_ref.selectStrongest(250); %only 250 (few hundreds) strongest keypoints

hold on
plot(surf_ref);

n_det_ref = 250; %number of detected keypoints on reference image

% %% scaling repeatibility
% m = 1.2; %scale factor
% counter = 0; %numerator of repeatibility (number of matches)
% 
% for i = 1:9 %nine scale factors
%     
%       scale(i) = m.^(i-1);
% 
%       surf_ref_scaled = scale(i) .* (surf_ref.Location)'; %scale the SURF keypoints of the ref to predict 
%  
%       target = imresize(ref,m.^(i-1)); %scale the reference image
%       
%       figure()
%       imshow(target);     
%       
%       surf_tar = detectSURFFeatures(target);
%       surf_tar = surf_tar.selectStrongest(250); %only 250 (few hundreds) strongest keypoints
%       
%       hold on
%       plot(surf_tar);
%       
%       n_det_tar = 250; %number of detected keypoints on target
%       hold on
%       plot(surf_ref_scaled(1,:),surf_ref_scaled(2,:),'O');
%       
%       for j = 1 : n_det_ref
%          for k = 1 : n_det_tar
%                  if (abs(surf_ref_scaled(1,j) - surf_tar.Location(k,1)) <= 2) && (abs(surf_ref_scaled(2,j) - surf_tar.Location(k,2)) <= 2)
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
% %% plot repeatibility over scaling
% figure1 = figure()
% plot(scale,repeatibility_scale,'-*g');
% title('SURF Repeatibility wrt scales')
% xlabel('Scale') 
% ylabel('Repeatibility(scale)') 
% saveas(figure1,'surf_repeat_scale.png')


%% rotation repeatibility
angle = 15; %rotation angle step (degree)
counter = 0; %numerator of repeatibility (number of matches) 

for i = 1:25 %number of rotations
    theta(i) = (i-1).*angle;
     
    rotation_matrix = [cosd(theta(i)) -sind(theta(i)); sind(theta(i)) cosd(theta(i))]; %rotation matrix
    
    %shift
    surf_ref_shift(1,:) = (surf_ref.Location(:,1))' - size(ref,2)/2;
    surf_ref_shift(2,:) = -( (surf_ref.Location(:,2))' - size(ref,1)/2 );

    surf_ref_rotated = rotation_matrix * surf_ref_shift(1:2,:); %rotate the SURF keypoints of the ref to predict 
    
    target = imrotate(ref,theta(i),'bicubic'); %rotate the ref image
     
    %shifting back
    surf_ref_rotated(1,:) = surf_ref_rotated(1,:) + size(target,2)/2;
    surf_ref_rotated(2,:) = - (surf_ref_rotated(2,:) - size(target,1)/2 );
    
    %figure()
    %imshow(target);     

    surf_tar = detectSURFFeatures(target);
    surf_tar = surf_tar.selectStrongest(250); %only 250 (few hundreds) strongest keypoints
       
    %hold on
    %plot(surf_tar);
       
    n_det_tar = 250; %number of detected keypoints on target
    %hold on
    %plot(surf_ref_rotated(1,:),surf_ref_rotated(2,:),'O');
       
    for j = 1 : n_det_ref
        for k = 1 : n_det_tar
                if (abs(surf_ref_rotated(1,j) - surf_tar.Location(k,1)) <= 2) && (abs(surf_ref_rotated(2,j) - surf_tar.Location(k,2)) <= 2)
                    counter = counter + 1; 
                    break;
                end
        end
    end 
     
    repeatibility_rotation(i) = counter ./ n_det_tar; %repeatibility formula
    counter = 0; %reset counter of matches
end

%% plot repeatibility over rotation angles
figure2 = figure()
plot(theta,repeatibility_rotation,'-Oy');
title('SURF Repeatibility wrt rotation')
xlabel('Theta (deg)') 
ylabel('Repeatibility(theta)') 
saveas(figure2,'surf_repeat_rotation.png')