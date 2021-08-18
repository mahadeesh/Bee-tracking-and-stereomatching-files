% This is a semi-automatic bee tracking program to track the freely flying
% bees 

% please see the readme file for complete usage information :
% https://github.com/mahadeesh/Bee-tracking-and-stereomatching-files 

% Inputs variable description:
%-----------------------------
% 1) (v)  - raw video.
% 2) (ov) - output video from image differencing step. See the help section for the detailed steps.
% 3) (frame number and frame count) - frame number from which you would like to start the tracking process.
% 4) (last_frame_number) - last frame number of the video(v). 
% 5) (CAM_INFO) - left or right camera information.
% 6) (bee_numbers) - unique ID number for each bee.  
% 7) (DATE) - date of the experiment or date of the data captured.
% 8) (bee_numbers_string) - variable which is used to indicate if the bee being is tracked is from left or right camera view

% Output variable description: 
%------------------------------
% 1) (J_BC_HT_CORDS_DATASET_#_CAM_##_REFER_BEE_###.mat)  - The final tracked 2D head and tail coordinates of a bee 
% 2) (OCCUL_FRAMES_DATASET_#_CAM_##_REFER_BEE_###.mat) -   Occlusion frames. See the help section for the detailed steps.
% 3) (RATIO_A_BY_B_DATASET_#_CAM_##_REFER_BEE_###.mat) - frame in which a bee image appear like a dot.
% 4) (ALL_FRAMES_DATASET_#_CAM_##_REFER_BEE_###.mat) - frame numbers of a track


% Sample inputs are provided at https://figshare.com/s/cb7edd04e5818607ea9a

% If there are questions in relation to the operation of this software,
% please contact: Mandiyam Mahadeeswara at m.mandiyam@uq.edu.au or mdevaraj@outlook.com

clear all;
clc;
close all;

% path of the raw and background subtracted video should be provided here

v = VideoReader('SAMPLE_INPUTS_FOR_BEE_TRACKING_PROGRAM/BC_EVENT_2_VIDEOS/MS_BC_2_CAM1_VIEW/SCI_DATA_BACKGROUND_SUB_RESULTS_CAM1_14682-16964_BC_EVENT_2.mp4'); % output video from imageJ
ov = VideoReader('SAMPLE_INPUTS_FOR_BEE_TRACKING_PROGRAM/BC_EVENT_2_VIDEOS/MS_BC_2_CAM1_VIEW/SCI_DATA_RAW_VIDEO_CAM1_frames 14682-16964_BC_EVENT_2.mp4'); % raw input video


frame_no = 0; frame_count = 0;   % enter the frame number from which tracking has to start 

last_frame_number = 1399;        % enter the last frame number of the raw video
CAM_INFO = 'DATASET_2_CAM_1';    % datset and camera view details
bee_numbers =  [1];              % enter the bee number
DATE = '16 FEB 2017';            % date of the experiment
bee_numbers_string = 'REFER_BEE';% variable which is used to indicate if the bee being is tracked is from left or right camera view

field_name_1 = 'Head_X_Cord'; field_name_2 = 'Head_Y_Cord'; field_name_3 = 'Tail_X_Cord'; field_name_4 = 'Tail_Y_Cord'; field_name_5 ='Bee_number'; 
field_name_6 ='Frame_number';

JUMP_HEAD_TAIL_CORD = struct(field_name_1,[],field_name_2,[],field_name_3,[],field_name_4,[],field_name_5,[],field_name_6,[]);

ov.CurrentTime = frame_no / ov.FrameRate;
v.CurrentTime = frame_no / v.FrameRate;
 
occluding_frames = []; 

JUMP_HT_CORDS_FILE_NAME = [sprintf('J_BC_HT_CORDS_%s_%s_%d',CAM_INFO,bee_numbers_string,bee_numbers) '.mat']; % output variable 1
% load(JUMP_HT_CORDS_FILE_NAME); 
 
OCCUL_INFO_FILE_NAME = [sprintf('OCCUL_FRAMES_%s_%s_%d',CAM_INFO,bee_numbers_string,bee_numbers) '.mat']; % output variable 2
% load(OCCUL_INFO_FILE_NAME); 

ALL_FRAMES_INFO_FILE_NAME = [sprintf('ALL_FRAMES_%s_%s_%d',CAM_INFO,bee_numbers_string,bee_numbers) '.mat']; % output variable 3
% load(ALL_FRAMES_INFO_FILE_NAME); 

CIRCLE_INFO_FILE_NAME = [sprintf('RATIO_A_BY_B_%s_%s_%d',CAM_INFO,bee_numbers_string,bee_numbers) '.mat']; % output variable 4
% load(CIRCLE_INFO_FILE_NAME); 


f1 = figure;
f1.Position = [-1790 -3 1540 804];

while hasFrame(ov)
    
    frame_count = frame_count +1
    
    img = readFrame(v);
    oImg = readFrame(ov);
    
    mask = double(imread('SAMPLE_INPUTS_FOR_BEE_TRACKING_PROGRAM/sci_data_mask1.png'))/255; % for CAM_1 view videos
  % mask = double(imread('SAMPLE_INPUTS_FOR_BEE_TRACKING_PROGRAM/sci_data_mask2.png'))/255; % for CAM_2 view videos
    
      
    img = rgb2gray(img);
    img = double(img).*mask;
    img = uint8(img);
         
    binimg = imbinarize(img,0.1);

    [B,L] = bwboundaries(binimg,'noholes');
    a = gca;
 
  ellipses = cell(size(B,1),1);

outer1_coords = repmat(struct('x',[],'y',[]),[size(B,1),1]);
outer2_coords = repmat(struct('x',[],'y',[]),[size(B,1),1]);

bee_count =0;

for k=1:length(B)
    
    bound = B{k};
    
    save_B{frame_count,k} = B{k};
 
    y = bound(:,1);
    x = bound(:,2);
    
    if (numel(x) >= 5) && (numel(x) <= 100)
        
       [ellipses{k} , elli_null] = fit_ellipse(x,y,a);
        
       elli_emp_count(frame_count,k) = elli_null;
       
       save_ellipses{frame_count,k} = ellipses{k};
        
        
        if (isempty(ellipses{k})) ||  elli_emp_count(frame_count,k) >= 1
                                   
            ellipses{k} = nan;
                       
            outer1_coords(k).x = nan;
            outer1_coords(k).y = nan;
        
            outer2_coords(k).x = nan;
            outer2_coords(k).y = nan;
            
            disp('Im empty: fit_elipse function is not returning head-tail coordinates');
                 
        else
        
        if ellipses{k}.phi > 0 && ellipses{k}.a > ellipses{k}.b
            
            phi = ellipses{k}.phi + pi/2;
            
        elseif ellipses{k}.phi < 0 && ellipses{k}.a < ellipses{k}.b
            
            phi = ellipses{k}.phi;
            
        elseif ellipses{k}.phi < 0 && ellipses{k}.a > ellipses{k}.b
            
            phi = ellipses{k}.phi + pi/2;
            
        elseif ellipses{k}.phi > 0 && ellipses{k}.a < ellipses{k}.b
            
            phi = ellipses{k}.phi;
            
        else
            
            phi = 0;
            
        end
    
        bee_count = bee_count +1;
      
       
        outer1_coords(k).x = ellipses{k}.X0_in + ellipses{k}.long_axis/2*sin(phi);
        outer1_coords(k).y = ellipses{k}.Y0_in + ellipses{k}.long_axis/2*cos(phi);
        
        outer2_coords(k).x = ellipses{k}.X0_in - ellipses{k}.long_axis/2*sin(phi);
        outer2_coords(k).y = ellipses{k}.Y0_in - ellipses{k}.long_axis/2*cos(phi);
        
        if outer1_coords(k).x < 2500
        
        one_prob_cord1_X(bee_count,frame_count) = outer1_coords(k).x ; 
        one_prob_cord1_Y(bee_count,frame_count) = outer1_coords(k).y ;
        
        two_prob_cord2_X(bee_count,frame_count) = outer2_coords(k).x ; 
        two_prob_cord2_Y(bee_count,frame_count) = outer2_coords(k).y ;
        
        each_elli_A_value(bee_count,frame_count) = ellipses{k}.a;
        each_elli_B_value(bee_count,frame_count) = ellipses{k}.b;
        
        CG_X_cord(bee_count,frame_count) = ellipses{k}.X0_in ;
        CG_Y_cord(bee_count,frame_count) = ellipses{k}.Y0_in ;
                
        else
            
            
        one_prob_cord1_X(bee_count,frame_count) = 0 ; 
        one_prob_cord1_Y(bee_count,frame_count) = 0 ;
        
        two_prob_cord2_X(bee_count,frame_count) = 0 ; 
        two_prob_cord2_Y(bee_count,frame_count) = 0 ;
        
        each_elli_A_value(bee_count,frame_count) = 0;
        each_elli_B_value(bee_count,frame_count) = 0;
        
        CG_X_cord(bee_count,frame_count) = 0 ;
        CG_Y_cord(bee_count,frame_count) = 0;
            
        end
        
     end
    
  end

end

num_bee_present_frame(frame_count)= bee_count;

for qq = 1:length(bee_numbers) 
    
    gg = bee_numbers(qq);
    
       if (frame_count - frame_no) == 1
                  
          imshow(oImg);
          fig = gcf;
          fig.Position = [-1790 -3 1540 804];
          hold on
          disp('First iteration: Click on the head and tail position of the bee which you would like to track');
          
          [c_first,r_first,p] = impixel;
          
          selected_bee_cord = [c_first(1) r_first(1)];
          
          % computing the selected bee's distance with all other head-tail points from the ellipse fit program, in
          % order to pick the correct coordinates
                    
       for mm = 1:bee_count  
    
       first_frame_dist_1(frame_count,mm) = sqrt((selected_bee_cord(1) - one_prob_cord1_X(mm,frame_count)).^2 + (selected_bee_cord(2) - one_prob_cord1_Y(mm,frame_count)).^2) ;
       first_frame_dist_2(frame_count,mm) = sqrt((selected_bee_cord(1) - two_prob_cord2_X(mm,frame_count)).^2 + (selected_bee_cord(2) - two_prob_cord2_Y(mm,frame_count)).^2) ;
         
       end
       
       
     first_frame_temp_dist_1 =  first_frame_dist_1(frame_count,:);
     cropped_first_frame_temp_dist_1 =  first_frame_temp_dist_1(first_frame_temp_dist_1 ~=0); 
     sort_first_frame_dist_1 = sort(cropped_first_frame_temp_dist_1,'ascend');
     
     first_frame_temp_dist_2 =  first_frame_dist_2(frame_count,:);
     cropped_first_frame_temp_dist_2 =  first_frame_temp_dist_2(first_frame_temp_dist_2 ~=0); 
     sort_first_frame_dist_2 = sort(cropped_first_frame_temp_dist_2,'ascend');
     
     % finding the selected bee's coordinates from the list of head-tail
     % coordinates
     
     first_frame_point_1_ind = find(first_frame_dist_1(frame_count,:) == sort_first_frame_dist_1(1));
     first_frame_point_2_ind = find(first_frame_dist_2(frame_count,:) == sort_first_frame_dist_2(1));
     
     % assigning selected bee's head and tail arbitrarily
          
      JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = one_prob_cord1_X(first_frame_point_1_ind,frame_count);
      JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = one_prob_cord1_Y(first_frame_point_2_ind,frame_count);

      JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = two_prob_cord2_X(first_frame_point_2_ind,frame_count);
      JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = two_prob_cord2_Y(first_frame_point_2_ind,frame_count);

      JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
      JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

      fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
      fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
      fin_position(gg,3) = 50 ;
      fin_label(gg,1) = gg; 
     
       % ellipse area calculation
       
       occu_area_ellipse(frame_count) = 3.14 * each_elli_A_value(first_frame_point_1_ind,frame_count) * each_elli_B_value(first_frame_point_1_ind,frame_count);
       RC_occu_area_ellipse(frame_count) = 0;
      
         occu_dist(frame_count,gg) = 50;
         FRAME_NUM(frame_count,gg) = frame_count;
         ratio_A_by_B(frame_count) = each_elli_A_value(first_frame_point_1_ind,frame_count)/each_elli_B_value(first_frame_point_1_ind,frame_count);
         
          
          continue;


      elseif (frame_count - frame_no) == 2
          

        prev_frame_tail_X_cord  = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg); 
        prev_frame_tail_Y_cord  = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg); 

        prev_frame_head_X_cord  = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg); 
        prev_frame_head_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg); 
        
        
        possib_1_ref_bee_HT_X_comp = prev_frame_head_X_cord - prev_frame_tail_X_cord;
        possib_1_ref_bee_HT_Y_comp = prev_frame_head_Y_cord - prev_frame_tail_Y_cord;
        
        possib_2_ref_bee_TH_X_comp = prev_frame_tail_X_cord - prev_frame_head_X_cord;
        possib_2_ref_bee_TH_Y_comp = prev_frame_tail_Y_cord - prev_frame_head_Y_cord;
        
         HT_dist(frame_count) = sqrt((possib_1_ref_bee_HT_X_comp^2) + (possib_1_ref_bee_HT_Y_comp^2));
        
        % two possible head-tail combinations and make them a unit vector
        
        possib_1_ref_vect = [possib_1_ref_bee_HT_X_comp possib_1_ref_bee_HT_Y_comp];
        possib_1_mag_ref_vect = rssq(possib_1_ref_vect);
        possib_1_unit_ref_vect = [ possib_1_ref_vect(1)/possib_1_mag_ref_vect  possib_1_ref_vect(2)/possib_1_mag_ref_vect]; 
        
        possib_2_ref_vect = [possib_2_ref_bee_TH_X_comp possib_2_ref_bee_TH_Y_comp];
        possib_2_mag_ref_vect = rssq(possib_2_ref_vect);
        possib_2_unit_ref_vect = [ possib_2_ref_vect(1)/possib_2_mag_ref_vect  possib_2_ref_vect(2)/possib_2_mag_ref_vect]; 
        
       % extracting the correct coordinates of the selected bee for frame 2 based on the nearest distance concept      
        
        for mm = 1:bee_count  
    
            second_frame_dist_1(frame_count,mm) = sqrt((prev_frame_head_X_cord - one_prob_cord1_X(mm,frame_count)).^2 + (prev_frame_head_Y_cord - one_prob_cord1_Y(mm,frame_count)).^2) ;
            second_frame_dist_2(frame_count,mm) = sqrt((prev_frame_head_X_cord - two_prob_cord2_X(mm,frame_count)).^2 + (prev_frame_head_Y_cord - two_prob_cord2_Y(mm,frame_count)).^2) ;
        end
        
        
     second_frame_temp_dist_1 =  second_frame_dist_1(frame_count,:);
     cropped_second_frame_temp_dist_1 =  second_frame_temp_dist_1(second_frame_temp_dist_1 ~=0); 
     sort_second_frame_dist_1 = sort(cropped_second_frame_temp_dist_1,'ascend');
     
     second_frame_temp_dist_2 =  second_frame_dist_2(frame_count,:);
     cropped_second_frame_temp_dist_2 =  second_frame_temp_dist_2(second_frame_temp_dist_2 ~=0); 
     sort_second_frame_dist_2 = sort(cropped_second_frame_temp_dist_2,'ascend');
     
        
      if sort_second_frame_dist_1(1) <= 50 || sort_second_frame_dist_2(1) <= 50 
            
         second_frame_point_1_ind = find(second_frame_dist_1(frame_count,:) == sort_second_frame_dist_1(1));
         second_frame_point_2_ind = find(second_frame_dist_2(frame_count,:) == sort_second_frame_dist_2(1));
                 
          assm_second_frame_head_X = one_prob_cord1_X(second_frame_point_1_ind,frame_count);
          assm_second_frame_head_Y = one_prob_cord1_Y(second_frame_point_1_ind,frame_count);
          
          assm_second_frame_tail_X = two_prob_cord2_X(second_frame_point_2_ind,frame_count);
          assm_second_frame_tail_Y = two_prob_cord2_Y(second_frame_point_2_ind,frame_count);
          
          head_velo_X_1 = (assm_second_frame_head_X - prev_frame_head_X_cord)/0.0050;
          head_velo_Y_1 = (assm_second_frame_head_Y - prev_frame_head_Y_cord)/0.0050;
          
          head_1 = sqrt((assm_second_frame_head_X - prev_frame_head_X_cord)^2 + (assm_second_frame_head_Y - prev_frame_head_Y_cord)^2);
          
          head_velo_X_2 = (assm_second_frame_head_X - prev_frame_tail_X_cord)/0.0050;
          head_velo_Y_2 = (assm_second_frame_head_Y - prev_frame_tail_Y_cord)/0.0050;
          
          head_2 = sqrt((assm_second_frame_head_X - prev_frame_tail_X_cord)^2 + (assm_second_frame_head_Y - prev_frame_tail_Y_cord)^2);
          
          tail_velo_X_1 = (assm_second_frame_tail_X - prev_frame_head_X_cord)/0.0050;
          tail_velo_Y_1 = (assm_second_frame_tail_Y - prev_frame_head_Y_cord)/0.0050;
          
          tail_1 = sqrt((assm_second_frame_tail_X - prev_frame_head_X_cord)^2 + (assm_second_frame_tail_Y - prev_frame_head_Y_cord)^2);
          
          tail_velo_X_2 = (assm_second_frame_tail_X - prev_frame_tail_X_cord)/0.0050;
          tail_velo_Y_2 = (assm_second_frame_tail_Y - prev_frame_tail_Y_cord)/0.0050;
          
          tail_2 = sqrt((assm_second_frame_tail_X - prev_frame_tail_X_cord)^2 + (assm_second_frame_tail_Y - prev_frame_tail_Y_cord)^2);
          
          vel_1(:,1) = [head_velo_X_1; head_velo_Y_1];
          vel_1(:,2) = [head_velo_X_2; head_velo_Y_2];
          vel_1(:,3) = [tail_velo_X_1; tail_velo_Y_1];
          vel_1(:,4) = [tail_velo_X_2; tail_velo_Y_2];
          
          % computing the displacements
          
          dist_HT_4_combi = [head_1 head_2 tail_1 tail_2]; 
                    
          long_dist_ind = find(dist_HT_4_combi == max(dist_HT_4_combi));
         
          % computing the reference motion vector
          
          dir_vect = [vel_1(1,long_dist_ind) vel_1(2,long_dist_ind)];
          mag_dir_vect = rssq(dir_vect);
          unit_dir_vect = [ dir_vect(1)/mag_dir_vect  dir_vect(2)/mag_dir_vect]; 
          
        % for frame 2 ... 
        
        possib_1_second_frame_body_X_comp = assm_second_frame_head_X - assm_second_frame_tail_X;
        possib_1_second_frame_body_Y_comp = assm_second_frame_head_Y - assm_second_frame_tail_Y;
          
        possib_1_2nd_frame_vect = [possib_1_second_frame_body_X_comp possib_1_second_frame_body_Y_comp];
        possib_1_mag_2nd_frame_vect = rssq(possib_1_2nd_frame_vect);
        possib_1_2nd_frame_vect_unit = [ possib_1_2nd_frame_vect(1)/possib_1_mag_2nd_frame_vect  possib_1_2nd_frame_vect(2)/possib_1_mag_2nd_frame_vect]; 
        
        possib_2_second_frame_body_X_comp =  assm_second_frame_tail_X - assm_second_frame_head_X;
        possib_2_second_frame_body_Y_comp =  assm_second_frame_tail_Y - assm_second_frame_head_Y;
          
        possib_2_2nd_frame_vect = [possib_2_second_frame_body_X_comp possib_2_second_frame_body_Y_comp];
        possib_2_mag_2nd_frame_vect = rssq(possib_2_2nd_frame_vect);
        possib_2_2nd_frame_vect_unit = [ possib_2_2nd_frame_vect(1)/possib_2_mag_2nd_frame_vect  possib_2_2nd_frame_vect(2)/possib_2_mag_2nd_frame_vect]; 
        
        % computing the angle between the reference motion vector and the
        % two possible head tail combinations for frame #2
        
        possib_1_dot_2nd_frame_vect = dot(unit_dir_vect,possib_1_2nd_frame_vect_unit);
        possib_2_dot_2nd_frame_vect = dot(unit_dir_vect,possib_2_2nd_frame_vect_unit);  
        
        possib_1_theta = acosd(possib_1_dot_2nd_frame_vect);
        possib_2_theta = acosd(possib_2_dot_2nd_frame_vect);
        
        
        if possib_2_theta <= 30
            
          JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = two_prob_cord2_X(second_frame_point_2_ind,frame_count);
          JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = two_prob_cord2_Y(second_frame_point_2_ind,frame_count);

          JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = one_prob_cord1_X(second_frame_point_1_ind,frame_count);
          JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = one_prob_cord1_Y(second_frame_point_1_ind,frame_count);

          JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
          JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

          fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
          fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
          fin_position(gg,3) = 1.5 ;
          fin_label(gg,1) = gg; 
          
          occu_area_ellipse(frame_count) = 3.14 * each_elli_A_value(second_frame_point_1_ind,frame_count) * each_elli_B_value(second_frame_point_1_ind,frame_count);
          RC_occu_area_ellipse(frame_count) = occu_area_ellipse(frame_count) - occu_area_ellipse(frame_count-1);
          
                    
        elseif possib_1_theta <= 30
            
          JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = one_prob_cord1_X(second_frame_point_1_ind,frame_count);
          JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = one_prob_cord1_Y(second_frame_point_1_ind,frame_count);

          JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = two_prob_cord2_X(second_frame_point_2_ind,frame_count);
          JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = two_prob_cord2_Y(second_frame_point_2_ind,frame_count);

          JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
          JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

          fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
          fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
          fin_position(gg,3) = 1.5 ;
          fin_label(gg,1) = gg; 
          
          occu_area_ellipse(frame_count) = 3.14 * each_elli_A_value(second_frame_point_1_ind,frame_count) * each_elli_B_value(second_frame_point_1_ind,frame_count);
          RC_occu_area_ellipse(frame_count) = occu_area_ellipse(frame_count) - occu_area_ellipse(frame_count-1);
            
        else
           
            disp('Change the thresholds');
            
        end
        
        % computing the angle between the reference motion vector and the
        % two possible head tail combination for frame # 1
        
        possib_1_dot_1st_frame_vect = dot(unit_dir_vect, possib_1_unit_ref_vect);
        possib_2_dot_1st_frame_vect = dot(unit_dir_vect, possib_2_unit_ref_vect); 
        
        possib_1_theta_1st_frame = acosd(possib_1_dot_1st_frame_vect);
        possib_2_theta_1st_frame = acosd(possib_2_dot_1st_frame_vect);
        
        
        if possib_1_theta_1st_frame <= 16
            
          JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg) = prev_frame_head_X_cord; % two_prob_cord2_X(first_frame_point_1_ind,frame_count-1);
          JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg) = prev_frame_head_Y_cord; % two_prob_cord2_Y(first_frame_point_1_ind,frame_count);

          JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg) = prev_frame_tail_X_cord ;% one_prob_cord1_X(first_frame_point_1_ind,frame_count);
          JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg) = prev_frame_tail_Y_cord ;% one_prob_cord1_Y(first_frame_point_1_ind,frame_count);

          JUMP_HEAD_TAIL_CORD(frame_count-1).Frame_number(gg) = frame_count-1;
          JUMP_HEAD_TAIL_CORD(frame_count-1).Bee_number(gg) = gg;

          fin_position_frame_1(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg); 
          fin_position_frame_1(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg);
          fin_position_frame_1(gg,3) = 1.5 ;
          fin_label(gg,1) = gg; 
          
          
        elseif possib_2_theta_1st_frame <= 16
            
          JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg) = prev_frame_tail_X_cord ; % one_prob_cord1_X(first_frame_point_1_ind,frame_count);
          JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg) = prev_frame_tail_Y_cord ; % one_prob_cord1_Y(first_frame_point_1_ind,frame_count);

          JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg) = prev_frame_head_X_cord ; % two_prob_cord2_X(first_frame_point_1_ind,frame_count);
          JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg) = prev_frame_head_Y_cord ; % two_prob_cord2_Y(first_frame_point_1_ind,frame_count);

          JUMP_HEAD_TAIL_CORD(frame_count-1).Frame_number(gg) = frame_count-1;
          JUMP_HEAD_TAIL_CORD(frame_count-1).Bee_number(gg) = gg;

          fin_position_frame_1(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg); 
          fin_position_frame_1(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg);
          fin_position_frame_1(gg,3) = 1.5 ;
          fin_label_frame_1(gg,1) = gg; 
            
          
        else
           
            disp('Change the thresholds');
            
        end
        
        
       imshow(oImg);
       fig = gcf;
       fig.Position = [-1790 -3 1540 804];
       hold on
       plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg),'o','Color','w','MarkerFaceColor','w','MarkerSize',2.5);
       hold on 
       plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg),'o','Color','b','MarkerFaceColor','b','MarkerSize',2.5);
       hold on
       plot(JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg),'o','Color','w','MarkerFaceColor','w','MarkerSize',2.5);
       hold on 
       plot(JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg),'o','Color','b','MarkerFaceColor','b','MarkerSize',2.5);
       hold on
       
       quest = sprintf('Check if the head and tail is reversed: If yes, enter "Y" or else "N"'); 
       str = input(quest,'s');
       
       if str == 'y'
           
         H_1F_X = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg);  H_1F_Y = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg); 
         T_1F_X = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg);  T_1F_Y = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg); 
         
         H_2F_X = JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg);  H_2F_Y = JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg); 
         T_2F_X = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg);  T_2F_Y = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg); 
           
           
          JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) =  H_2F_X; 
          JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) =  H_2F_Y; 

          JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = T_2F_X;
          JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = T_2F_Y;
          
          JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg) =  H_1F_X; 
          JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg) =  H_1F_Y; 

          JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg) = T_1F_X;
          JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg) = T_1F_Y;
          
           imshow(oImg);
           fig = gcf;
           fig.Position = [-1790 -3 1540 804];
           hold on
           plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg),'o','Color','m','MarkerFaceColor','m','MarkerSize',2.5);
           hold on 
           plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg),'o','Color','g','MarkerFaceColor','g','MarkerSize',2.5);
           hold on
           plot(JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg),'o','Color','m','MarkerFaceColor','m','MarkerSize',2.5);
           hold on 
           plot(JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg),'o','Color','g','MarkerFaceColor','g','MarkerSize',2.5);
          

       end

        
      end  
     
         RC_occu_area_ellipse(frame_count) = occu_area_ellipse(frame_count)- occu_area_ellipse(frame_count-1);
         occu_dist(frame_count,gg) = 50;
         FRAME_NUM(frame_count,gg) = frame_count;
          
          
      else
           
          
        if occu_dist(frame_count-1,gg) <= 10
            
          % when there is an occlusion event, code is not getting into this part of the if statement
            
           occu_second_frame_head_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg);
           occu_second_frame_head_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg);
           
           occu_first_frame_head_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Head_X_Cord(gg);
           occu_first_frame_head_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Head_Y_Cord(gg);
           
           occu_speed_12_X_cord =  (occu_second_frame_head_X_cord - occu_first_frame_head_X_cord)/0.0050;
           occu_speed_12_Y_cord =  (occu_second_frame_head_Y_cord - occu_first_frame_head_Y_cord)/0.0050;
                    
           occu_third_frame_head_X_cord = ((occu_speed_12_X_cord)*0.0050) + occu_second_frame_head_X_cord;
           occu_third_frame_head_Y_cord = ((occu_speed_12_Y_cord)*0.0050) + occu_second_frame_head_Y_cord;
          
           %%% for tail%%%%           
           
           occu_second_frame_tail_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg);
           occu_second_frame_tail_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg);
           
           occu_first_frame_tail_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Tail_X_Cord(gg);
           occu_first_frame_tail_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Tail_Y_Cord(gg);
           
           occu_speed_12_X_cord_tail =  (occu_second_frame_tail_X_cord - occu_first_frame_tail_X_cord)/0.0050;
           occu_speed_12_Y_cord_tail =  (occu_second_frame_tail_Y_cord - occu_first_frame_tail_Y_cord)/0.0050;
           
                     
           occu_third_frame_tail_X_cord = ((occu_speed_12_X_cord_tail)*0.0050) + occu_second_frame_tail_X_cord;
           occu_third_frame_tail_Y_cord = ((occu_speed_12_Y_cord_tail)*0.0050) + occu_second_frame_tail_Y_cord;
                      
           plot(occu_third_frame_head_X_cord,occu_third_frame_head_Y_cord,'o','Color','m','MarkerFaceColor','m','MarkerSize',4.5);
           hold on
           plot(occu_third_frame_tail_X_cord,occu_third_frame_tail_Y_cord,'o','Color','g','MarkerFaceColor','g','MarkerSize',4.5);
            
           JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = occu_third_frame_head_X_cord;
           JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = occu_third_frame_head_Y_cord;
           JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = occu_third_frame_tail_X_cord;
           JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = occu_third_frame_tail_Y_cord;
           JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
           JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;
           
           for mm = 1:bee_count  

               occu_third_frame_dist_1(frame_count,mm) = sqrt((occu_third_frame_head_X_cord - CG_X_cord(mm,frame_count)).^2 + (occu_third_frame_head_Y_cord - CG_Y_cord(mm,frame_count)).^2) ;
           
           end
           
           occu_third_frame_temp_dist_1 =  occu_third_frame_dist_1(frame_count,:);
           occu_cropped_third_frame_temp_dist_1 =  occu_third_frame_temp_dist_1(occu_third_frame_temp_dist_1 ~=0); 
           occu_sort_third_frame_dist_1 = sort(occu_cropped_third_frame_temp_dist_1,'ascend');

           occu_third_min_sort = min(occu_sort_third_frame_dist_1);
           
           occu_dist(frame_count,gg) = occu_third_min_sort;
           
           combined_occu_third_fram_dist = [occu_third_frame_dist_1(frame_count,:)];
           [occ_CP_bee_ind1,occ_CP_bee_ind2] = find( combined_occu_third_fram_dist <= occu_third_min_sort);
           
           occu_area_ellipse(frame_count) = 3.14 * each_elli_A_value(occ_CP_bee_ind2,frame_count) * each_elli_B_value(occ_CP_bee_ind2,frame_count);
           
           RC_occu_area_ellipse(frame_count) = occu_area_ellipse(frame_count) - occu_area_ellipse(frame_count-1);
           
           if RC_occu_area_ellipse(frame_count) <= 8 
               
               occu_dist(frame_count,gg) = 16;
               
           elseif RC_occu_area_ellipse(frame_count) >= 10 
               
               base_ellipse = occu_area_ellipse(frame_count-1);
               
               occu_area_ellipse(frame_count) = base_ellipse;
               
               occu_dist(frame_count,gg) = 5;
               
           end
           
                       
        else
            
           % tracking from 3rd frame onwards without any occlusion case
           % will begin here
           
           second_frame_head_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg);
           second_frame_head_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg);
           
           first_frame_head_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Head_X_Cord(gg);
           first_frame_head_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Head_Y_Cord(gg);
           
           speed_12_X_cord =  (second_frame_head_X_cord - first_frame_head_X_cord)/0.0050;
           speed_12_Y_cord =  (second_frame_head_Y_cord - first_frame_head_Y_cord)/0.0050;
                    
           third_frame_head_X_cord = ((speed_12_X_cord)*0.0050) + second_frame_head_X_cord;
           third_frame_head_Y_cord = ((speed_12_Y_cord)*0.0050) + second_frame_head_Y_cord;
          
           % for tail           
           
           second_frame_tail_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg);
           second_frame_tail_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg);
           
           first_frame_tail_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Tail_X_Cord(gg);
           first_frame_tail_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Tail_Y_Cord(gg);
           
           speed_12_X_cord_tail =  (second_frame_tail_X_cord - first_frame_tail_X_cord)/0.0050;
           speed_12_Y_cord_tail =  (second_frame_tail_Y_cord - first_frame_tail_Y_cord)/0.0050;
                                
           third_frame_tail_X_cord = ((speed_12_X_cord_tail)*0.0050) + second_frame_tail_X_cord;
           third_frame_tail_Y_cord = ((speed_12_Y_cord_tail)*0.0050) + second_frame_tail_Y_cord;
           
          imshow(oImg);
          fig = gcf;
          fig.Position = [-1790 -3 1540 804];
          hold on
                      
           plot(third_frame_head_X_cord,third_frame_head_Y_cord,'o','Color','w','MarkerFaceColor','w','MarkerSize',2.5);
           hold on
           plot(third_frame_tail_X_cord,third_frame_tail_Y_cord,'o','Color','b','MarkerFaceColor','b','MarkerSize',2.5);
           hold on   
           
           % eliminating the incorrect extrapolated points to pick the
           % correct bee
           
           EXTR_3_FRAME_X = third_frame_head_X_cord - third_frame_tail_X_cord;
           EXTR_3_FRAME_Y = third_frame_head_Y_cord - third_frame_tail_Y_cord; 
           
           EXTR_3_FRAME_VECT = [EXTR_3_FRAME_X  EXTR_3_FRAME_Y];
           MAG_EXTR_3_FRAME_VECT =  rssq(EXTR_3_FRAME_VECT);
           EXTR_3_FRAME_UNIT_VECT = [EXTR_3_FRAME_VECT(1)/MAG_EXTR_3_FRAME_VECT  EXTR_3_FRAME_VECT(2)/MAG_EXTR_3_FRAME_VECT];
           
           PREV_3_FRAME_X = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg); 
           PREV_3_FRAME_Y = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg);
           
           PREV_3_FRAME_VECT = [PREV_3_FRAME_X PREV_3_FRAME_Y];
           MAG_PREV_3_FRAME_VECT =  rssq(PREV_3_FRAME_VECT);
           PREV_3_FRAME_UNIT_VECT = [PREV_3_FRAME_VECT(1)/MAG_PREV_3_FRAME_VECT  PREV_3_FRAME_VECT(2)/MAG_PREV_3_FRAME_VECT];
           
           
           EXTR_3_FRAME_DOT_PROD = dot(PREV_3_FRAME_UNIT_VECT, EXTR_3_FRAME_UNIT_VECT);
           EXTR_3_FRAME_THETA(frame_count) = acosd(EXTR_3_FRAME_DOT_PROD);
           
           if EXTR_3_FRAME_THETA(frame_count) >= 30
               
                  plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg),'o','Color','m','MarkerFaceColor','m','MarkerSize',3.5);
                  hold on
                  plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg),'o','Color','g','MarkerFaceColor','g','MarkerSize',3.5)
                                    
                  disp('Click on the center of gravity(CG) of the selected bee (the extrapolated head-tail values might give incorrect coordinates)');
                  
                  [c_3_FRA,r_3_FRA,p] = impixel;
                  selected_bee_cord_EXTR_3_FRAME = [c_3_FRA(1) r_3_FRA(1)];
                  
                  third_frame_head_X_cord = selected_bee_cord_EXTR_3_FRAME(1) ; third_frame_head_Y_cord = selected_bee_cord_EXTR_3_FRAME(2);
              
           end
           
                      
           for mm = 1:bee_count  

               third_frame_dist_1(frame_count,mm) = sqrt((third_frame_head_X_cord - CG_X_cord(mm,frame_count)).^2 + (third_frame_head_Y_cord - CG_Y_cord(mm,frame_count)).^2) ;
               
           end

        
                third_frame_temp_dist_1 =  third_frame_dist_1(frame_count,:);
                cropped_third_frame_dist_1 =  third_frame_temp_dist_1(third_frame_temp_dist_1 ~=0);  
                sort_third_frame_min_dist_1 = sort(cropped_third_frame_dist_1,'ascend');
                
                dist_based_ind = find(third_frame_dist_1(frame_count,1:mm) == sort_third_frame_min_dist_1(1));
                
                combined_third_fram_dist = third_frame_dist_1(frame_count,:);
                [non_occ_CP_bee_ind1,non_occ_CP_bee_ind22] = find( combined_third_fram_dist <= 22);
                
        if ~isempty(non_occ_CP_bee_ind22)         
           
                % select the correct coordinates using the angle
                % information as well
                
               
                BODY_frame_X_1_ST_TIME = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg); 
                BODY_frame_Y_1_ST_TIME = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg); 
                                

                BODY_INFO_1_ST_TIME = [BODY_frame_X_1_ST_TIME  BODY_frame_Y_1_ST_TIME];
                MAG_BODY_INFO_1_ST_TIME =  rssq(BODY_INFO_1_ST_TIME);
                BODY_INFO_1_ST_TIME_UNIT = [BODY_INFO_1_ST_TIME(1)/MAG_BODY_INFO_1_ST_TIME  BODY_INFO_1_ST_TIME(2)/MAG_BODY_INFO_1_ST_TIME];
                                
                
                for hh = 1:length(non_occ_CP_bee_ind22) 
                    
                  non_occu_CP_bee_ind_array(1,hh) = combined_third_fram_dist(non_occ_CP_bee_ind22(hh));
                
                  opt_2_X_cord = two_prob_cord2_X(non_occ_CP_bee_ind22(hh),frame_count) - one_prob_cord1_X(non_occ_CP_bee_ind22(hh),frame_count);
                  opt_2_Y_cord = two_prob_cord2_Y(non_occ_CP_bee_ind22(hh),frame_count) - one_prob_cord1_Y(non_occ_CP_bee_ind22(hh),frame_count);

                  opt_2_vect = [opt_2_X_cord opt_2_Y_cord];
                  mag_opt_2_vect = rssq(opt_2_vect);
                  unit_opt_2_vect = [opt_2_vect(1)/mag_opt_2_vect opt_2_vect(2)/mag_opt_2_vect];

                  opt_1_X_cord = one_prob_cord1_X(non_occ_CP_bee_ind22(hh),frame_count) - two_prob_cord2_X(non_occ_CP_bee_ind22(hh),frame_count);
                  opt_1_Y_cord = one_prob_cord1_Y(non_occ_CP_bee_ind22(hh),frame_count) - two_prob_cord2_Y(non_occ_CP_bee_ind22(hh),frame_count);

                  opt_1_vect = [opt_1_X_cord opt_1_Y_cord];
                  mag_opt_1_vect = rssq(opt_1_vect);
                  unit_opt_1_vect = [opt_1_vect(1)/mag_opt_1_vect opt_1_vect(2)/mag_opt_1_vect]; 
                  
                    possib_1_dot_BODY_INFO_1_ST_TIME = dot(BODY_INFO_1_ST_TIME_UNIT, unit_opt_1_vect);
                    possib_2_dot_BODY_INFO_1_ST_TIME = dot(BODY_INFO_1_ST_TIME_UNIT, unit_opt_2_vect); 

                    possib_1_theta_BODY_INFO_1_ST_TIME(1,hh) = acosd(possib_1_dot_BODY_INFO_1_ST_TIME);
                    possib_2_theta_BODY_INFO_1_ST_TIME(1,hh) = acosd(possib_2_dot_BODY_INFO_1_ST_TIME);
                    
                
                end
                               
                possible_theta_1_ST_TIME = [possib_1_theta_BODY_INFO_1_ST_TIME;possib_2_theta_BODY_INFO_1_ST_TIME];
                [lw_x_pos,lw_y_pos] = find(possible_theta_1_ST_TIME == min(min(possible_theta_1_ST_TIME)));
                CORRECT_1_ST_TIME_IND = non_occ_CP_bee_ind22(lw_y_pos);
                IND_PICK_1_ST_TIME(frame_count) = CORRECT_1_ST_TIME_IND;
                
                % For every frame new ellipse area and rate of change is computed 
                                
                if CORRECT_1_ST_TIME_IND == dist_based_ind
                    
                    occu_area_ellipse(frame_count) = 3.14 * each_elli_A_value(CORRECT_1_ST_TIME_IND,frame_count) * each_elli_B_value(CORRECT_1_ST_TIME_IND,frame_count);
                    RC_occu_area_ellipse(frame_count) = occu_area_ellipse(frame_count) - occu_area_ellipse(frame_count-1);
                    ratio_A_by_B(frame_count) = each_elli_A_value(CORRECT_1_ST_TIME_IND,frame_count) / each_elli_B_value(CORRECT_1_ST_TIME_IND,frame_count); 
                    
                    
                elseif CORRECT_1_ST_TIME_IND ~= dist_based_ind
                    
                    
                          disp('Click on the CG of the selected bee (angle and distance indices mismatch)');
                          [c_ang,r_ang,p] = impixel;
                          ang_dist_mismatch_CG_X_cord = c_ang(1) ; ang_dist_mismatch_CG_Y_cord = r_ang(1);
                  
                              for nn = 1:bee_count  

                                   ang_dist_mismatch_cond(frame_count,nn) = sqrt((ang_dist_mismatch_CG_X_cord - CG_X_cord(nn,frame_count)).^2 + (ang_dist_mismatch_CG_Y_cord - CG_Y_cord(nn,frame_count)).^2) ;

                              end
                    
                            ang_dist_mismatch_dist =  ang_dist_mismatch_cond(frame_count,:);
                            cropped_ang_dist_mismatch_dist =  ang_dist_mismatch_dist(ang_dist_mismatch_dist ~=0);  
                            sort_ang_dist_mismatch_dist = sort(cropped_ang_dist_mismatch_dist,'ascend');

                            CORRECT_1_ST_TIME_IND = find(cropped_ang_dist_mismatch_dist == sort_ang_dist_mismatch_dist(1));

                            occu_area_ellipse(frame_count) = 3.14 * each_elli_A_value(CORRECT_1_ST_TIME_IND,frame_count) * each_elli_B_value(CORRECT_1_ST_TIME_IND,frame_count);
                            RC_occu_area_ellipse(frame_count) = occu_area_ellipse(frame_count) - occu_area_ellipse(frame_count-1);
                            ratio_A_by_B(frame_count) = each_elli_A_value(CORRECT_1_ST_TIME_IND,frame_count) / each_elli_B_value(CORRECT_1_ST_TIME_IND,frame_count); 
                    
                end
                
                
                % When the angle difference is less than 1 deg, swapping
                % happens leading to wrong bee selection, so such scenarios
                % are handled here
                
                COND_1 = length(min(possible_theta_1_ST_TIME));
                COND_2 = abs(diff(min(possible_theta_1_ST_TIME)));
                COND_3 = min(min(possible_theta_1_ST_TIME));
                
                
                if COND_1 == 2 && COND_2 <= 1 && COND_3 < 10
                    
                  imshow(oImg);
                  fig = gcf;
                  fig.Position = [-1790 -3 1540 804];
                  hold on
                  plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg),'o','Color','m','MarkerFaceColor','m','MarkerSize',3.5);
                  hold on
                  plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg),'o','Color','g','MarkerFaceColor','g','MarkerSize',3.5)
                                    
                  disp('Click on the CG of the selected bee (angle diffence error)');
                  
                  [c_tiny_ang,r_tiny_ang,p] = impixel;
                  selected_bee_cord = [c_tiny_ang(1) r_tiny_ang(1)];
                            
                   for mm = 1:length(non_occ_CP_bee_ind22)  

                   duo_small_angle_dist(1,mm) = sqrt((selected_bee_cord(1) - CG_X_cord(non_occ_CP_bee_ind22(1,mm),frame_count)).^2 + (selected_bee_cord(2) - CG_Y_cord(non_occ_CP_bee_ind22(1,mm),frame_count)).^2) ;
                   
                   end
                   
                   [small_duo_val,small_duo_ind] = find(duo_small_angle_dist == min(duo_small_angle_dist(1,1:mm)));
                   
                CORRECT_1_ST_TIME_IND = non_occ_CP_bee_ind22(small_duo_ind);
                occu_area_ellipse(frame_count) = 3.14 * each_elli_A_value(CORRECT_1_ST_TIME_IND,frame_count) * each_elli_B_value(CORRECT_1_ST_TIME_IND,frame_count);
                RC_occu_area_ellipse(frame_count) = occu_area_ellipse(frame_count) - occu_area_ellipse(frame_count-1);
                ratio_A_by_B(frame_count) = each_elli_A_value(CORRECT_1_ST_TIME_IND,frame_count) / each_elli_B_value(CORRECT_1_ST_TIME_IND,frame_count); 
                   
                end
                
                clear possib_1_theta_BODY_INFO_1_ST_TIME possib_2_theta_BODY_INFO_1_ST_TIME possible_theta_1_ST_TIME;
                   
           % no occlusion scenario or normal scenario will start from here
                 
           if abs(RC_occu_area_ellipse(frame_count)) < 20 && ( ratio_A_by_B(frame_count) <= 0.9 || ratio_A_by_B(frame_count) >= 1.3)
 
                       
             for uu = 1: length(CORRECT_1_ST_TIME_IND)
             
                                       
              opt_2_X_cord = two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count) - one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count);
              opt_2_Y_cord = two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count) - one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count);

              opt_2_vect = [opt_2_X_cord opt_2_Y_cord];
              mag_opt_2_vect = rssq(opt_2_vect);
              unit_opt_2_vect = [opt_2_vect(1)/mag_opt_2_vect opt_2_vect(2)/mag_opt_2_vect];

              opt_1_X_cord = one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count) - two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count);
              opt_1_Y_cord = one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count) - two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count);

              opt_1_vect = [opt_1_X_cord opt_1_Y_cord];
              mag_opt_1_vect = rssq(opt_1_vect);
              unit_opt_1_vect = [opt_1_vect(1)/mag_opt_1_vect opt_1_vect(2)/mag_opt_1_vect];   
                    
             
                confirm_prev_frame_head_X_cord  = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg); 
                confirm_prev_frame_head_Y_cord  = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg); 

                confirm_prev_frame_tail_X_cord  = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg); 
                confirm_prev_frame_tail_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg); 

                third_non_occu_BODY_frame_X = confirm_prev_frame_head_X_cord - confirm_prev_frame_tail_X_cord; 
                third_non_occu_BODY_frame_Y = confirm_prev_frame_head_Y_cord - confirm_prev_frame_tail_Y_cord; 

                third_non_occu_BODY_frame = [third_non_occu_BODY_frame_X  third_non_occu_BODY_frame_Y];
                mag_third_non_occu_BODY_frame =  rssq(third_non_occu_BODY_frame);
                third_non_occu_BODY_frame_UNIT = [third_non_occu_BODY_frame(1)/mag_third_non_occu_BODY_frame  third_non_occu_BODY_frame(2)/mag_third_non_occu_BODY_frame];

                possib_1_dot_3rd_frame_vect_non_occu(1,uu) = dot(third_non_occu_BODY_frame_UNIT, unit_opt_1_vect);
                possib_2_dot_3rd_frame_vect_non_occu(1,uu) = dot(third_non_occu_BODY_frame_UNIT, unit_opt_2_vect); 

                possib_1_theta_3rd_frame_non_occu(1,uu) = acosd(possib_1_dot_3rd_frame_vect_non_occu(1,uu));
                possib_2_theta_3rd_frame_non_occu(1,uu) = acosd(possib_2_dot_3rd_frame_vect_non_occu(1,uu));
              

             end
             
             possib_1_theta_MIN = min(possib_1_theta_3rd_frame_non_occu(1,1:uu));
             possib_2_theta_MIN = min(possib_2_theta_3rd_frame_non_occu(1,1:uu));
             
             possib_1_DOT_MAG_MIN = min(possib_1_dot_3rd_frame_vect_non_occu(1,1:uu));
             possib_2_DOT_MAG_MIN = min(possib_2_dot_3rd_frame_vect_non_occu(1,1:uu));
             
             
             if possib_1_theta_MIN < 25
                 
                 possib_MIN_theta_POS = find(possib_1_theta_3rd_frame_non_occu == possib_1_theta_MIN);
                 possib_MIN_theta_COUNT = 1;
             
             elseif possib_2_theta_MIN < 25
               
                 possib_MIN_theta_POS = find(possib_2_theta_3rd_frame_non_occu == possib_2_theta_MIN);
                 possib_MIN_theta_COUNT = 2;
                 
             else
                 
                 
                 possib_MIN_theta_POS =0;
                 possib_MIN_theta_COUNT = 0;
             end
             
        if  (possib_MIN_theta_COUNT) == 1 ||  (possib_MIN_theta_COUNT) == 2 || (possib_MIN_theta_COUNT) == 0 % We are checking both distance and angle to pick the right coordinates
                       
             if possib_MIN_theta_COUNT ==  1
                                                                  
              JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count);
              JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count);

              JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count);
              JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count);

              JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
              JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

              fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
              fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
              fin_position(gg,3) = 1.5 ;
              fin_label(gg,1) = gg; 
              
              occu_dist(frame_count,gg) = 16;
              
              elseif possib_MIN_theta_COUNT == 2
                        
              JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count);
              JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count);

              JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count);
              JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count);

              JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
              JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

              fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
              fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
              fin_position(gg,3) = 1.5 ;
              fin_label(gg,1) = gg; 
              
              occu_dist(frame_count,gg) = 16;
              
             else
                 
      % when angle difference is more than 25 degs, just use the magnitude to solve this problem

              opt_2_X_cord = two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count) - one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count);
              opt_2_Y_cord = two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count) - one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count);

              opt_2_vect = [opt_2_X_cord opt_2_Y_cord];
              mag_opt_2_vect = rssq(opt_2_vect);
              unit_opt_2_vect = [opt_2_vect(1)/mag_opt_2_vect opt_2_vect(2)/mag_opt_2_vect];

              opt_1_X_cord = one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count) - two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count);
              opt_1_Y_cord = one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count) - two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count);

              opt_1_vect = [opt_1_X_cord opt_1_Y_cord];
              mag_opt_1_vect = rssq(opt_1_vect);
              unit_opt_1_vect = [opt_1_vect(1)/mag_opt_1_vect opt_1_vect(2)/mag_opt_1_vect];   
                    
             
                confirm_prev_frame_head_X_cord  = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg); 
                confirm_prev_frame_head_Y_cord  = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg); 

                confirm_prev_frame_tail_X_cord  = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg); 
                confirm_prev_frame_tail_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg); 

                third_non_occu_BODY_frame_X = confirm_prev_frame_head_X_cord - confirm_prev_frame_tail_X_cord; 
                third_non_occu_BODY_frame_Y = confirm_prev_frame_head_Y_cord - confirm_prev_frame_tail_Y_cord; 

                third_non_occu_BODY_frame = [third_non_occu_BODY_frame_X  third_non_occu_BODY_frame_Y];
                mag_third_non_occu_BODY_frame =  rssq(third_non_occu_BODY_frame);
                third_non_occu_BODY_frame_UNIT = [third_non_occu_BODY_frame(1)/mag_third_non_occu_BODY_frame  third_non_occu_BODY_frame(2)/mag_third_non_occu_BODY_frame];

                possib_1_dot_3rd_frame_vect_non_occu = dot(third_non_occu_BODY_frame_UNIT, unit_opt_1_vect);
                possib_2_dot_3rd_frame_vect_non_occu = dot(third_non_occu_BODY_frame_UNIT, unit_opt_2_vect); 

                possib_1_theta_3rd_frame_non_occu = acosd(possib_1_dot_3rd_frame_vect_non_occu);
                possib_2_theta_3rd_frame_non_occu = acosd(possib_2_dot_3rd_frame_vect_non_occu);
                
                if  possib_1_theta_3rd_frame_non_occu  <=60
                    
                  JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count);
                  JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count);

                  JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count);
                  JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count);

                  JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
                  JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

                  fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
                  fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
                  fin_position(gg,3) = 1.5 ;
                  fin_label(gg,1) = gg; 
                        
                  occu_dist(frame_count,gg) = 16;
                    
                elseif possib_2_theta_3rd_frame_non_occu  <=60
                    
                  JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count);
                  JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count);

                  JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count);
                  JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count);

                  JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
                  JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

                  fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
                  fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
                  fin_position(gg,3) = 1.5 ;
                  fin_label(gg,1) = gg; 
              
                  occu_dist(frame_count,gg) = 16;
                   
                else
                    
   % when the two angles are greater than 60 degs, manually feed the input
                    
                  imshow(oImg);
                  fig = gcf;
                  fig.Position = [-1790 -3 1540 804];
                  hold on
                  plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg),'o','Color','m','MarkerFaceColor','m','MarkerSize',3.5);
                  hold on
                  plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg),'o','Color','g','MarkerFaceColor','g','MarkerSize',3.5)
                  
                  [my,pr,p] = impixel;

                  JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = my(1);
                  JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = pr(1);

                  JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = my(2);
                  JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = pr(2);

                  JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
                  JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

                  fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
                  fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
                  fin_position(gg,3) = 1.5 ;
                  fin_label(gg,1) = gg; 

                  occu_dist(frame_count,gg) = 17;
                   
                   
                end
                
               
              disp('Check the thresholds');
                  
              end
        
          
             if uu < 2 
             
                 occu_dist(frame_count,gg) = 17;
                 
             elseif uu >= 2 && abs(RC_occu_area_ellipse(frame_count)) <= 10
             
                 occu_dist(frame_count,gg) = 17;
                 
             elseif uu >= 2 && abs(RC_occu_area_ellipse(frame_count)) >= 10
                 
                 occu_dist(frame_count,gg) = 17;
                 
             end
             
             
     else
             
          disp('Check the thresholds');
          
     end
             
    elseif RC_occu_area_ellipse(frame_count) <= -20 
        
        % This condition arises when the area of the fitted ellipse 
        % suddenly goes down. Ex: When a bee is leaving the arena, half of
        % its body might still be visible. So this option checks and exits
        % the program if the bee is on the boundary lines.
        
        
       % bee area is getting smaller or black slit condition
       disp('Check, if the selected bee is on the arena boundary line?')
       disp('If yes, click on the head and tail location from the visible part of the selected bee image.')
              
          imshow(oImg);
          fig = gcf;
          fig.Position = [-1790 -3 1540 804];
          hold on
                 
          plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg),'o','Color','m','MarkerFaceColor','m','MarkerSize',3.5);
          hold on
          plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg),'o','Color','g','MarkerFaceColor','g','MarkerSize',3.5)
                    
          [c_slit,r_slit,p] = impixel;
         
          JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = c_slit(1);
          JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = r_slit(1);

          JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = c_slit(2);
          JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = r_slit(2);

          JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
          JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

          fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
          fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
          fin_position(gg,3) = 1.5 ;
          fin_label(gg,1) = gg; 
          
          
          estimated_b_axis_length = (sqrt((JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg))^2 + (JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg))^2))/2;
          estimated_a_axis_length = estimated_b_axis_length/2;
          occu_area_ellipse(frame_count) = 3.14 * estimated_a_axis_length * estimated_b_axis_length;
                    
          occu_dist(frame_count,gg) = 17;
       
        
           elseif abs(RC_occu_area_ellipse(frame_count)) > 20 || (RC_occu_area_ellipse(frame_count) <= -20)
          
        %  All occlusion cases are handled here as the difference in occulsion area from the current frame to previous frame has suddenly changed
             
           occu_second_frame_head_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg);
           occu_second_frame_head_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg);
           
           occu_first_frame_head_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Head_X_Cord(gg);
           occu_first_frame_head_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Head_Y_Cord(gg);
           
           occu_speed_12_X_cord =  (occu_second_frame_head_X_cord - occu_first_frame_head_X_cord)/0.0050;
           occu_speed_12_Y_cord =  (occu_second_frame_head_Y_cord - occu_first_frame_head_Y_cord)/0.0050;
                    
           occu_third_frame_head_X_cord = ((occu_speed_12_X_cord)*0.0050) + occu_second_frame_head_X_cord;
           occu_third_frame_head_Y_cord = ((occu_speed_12_Y_cord)*0.0050) + occu_second_frame_head_Y_cord;
          
           %% for tail           
           
           occu_second_frame_tail_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg);
           occu_second_frame_tail_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg);
           
           occu_first_frame_tail_X_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Tail_X_Cord(gg);
           occu_first_frame_tail_Y_cord = JUMP_HEAD_TAIL_CORD(frame_count-2).Tail_Y_Cord(gg);
           
           occu_speed_12_X_cord_tail =  (occu_second_frame_tail_X_cord - occu_first_frame_tail_X_cord)/0.0050;
           occu_speed_12_Y_cord_tail =  (occu_second_frame_tail_Y_cord - occu_first_frame_tail_Y_cord)/0.0050;
           
                     
           occu_third_frame_tail_X_cord = ((occu_speed_12_X_cord_tail)*0.0050) + occu_second_frame_tail_X_cord;
           occu_third_frame_tail_Y_cord = ((occu_speed_12_Y_cord_tail)*0.0050) + occu_second_frame_tail_Y_cord;
           

          imshow(oImg);
          fig = gcf;
          fig.Position = [-1790 -3 1540 804];
          hold on
%           plot(occu_third_frame_head_X_cord,occu_third_frame_head_Y_cord,'o','Color','c','MarkerFaceColor','c','MarkerSize',4.5);
%           hold on
%           plot(occu_third_frame_tail_X_cord,occu_third_frame_tail_Y_cord,'o','Color','y','MarkerFaceColor','y','MarkerSize',4.5);
%           hold on
          
          for pp = 1:4
    
          last_4_frames_occu_area_ellipse(1,pp) = occu_area_ellipse(frame_count-pp);

          end

          diff_occ_area = (diff(last_4_frames_occu_area_ellipse(1,1:pp)));
          count_of_diff_occ_area = nnz(~diff_occ_area);
          
          EXTRA_OCCL_COND = sqrt((occu_third_frame_head_X_cord - occu_third_frame_tail_X_cord)^2 + (occu_third_frame_head_Y_cord - occu_third_frame_tail_Y_cord)^2);      
          
           
          if EXTRA_OCCL_COND >=25 || count_of_diff_occ_area  == 3
              
              disp('Click on the head and tail location of the selected bee (incorrect extrapolated coordinates due to high number of occlusion events)');
              
              plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg),JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg),'-o','Color','m','MarkerFaceColor','m','MarkerSize',3.5);
              hold on
              plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg),JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg),'-o','Color','g','MarkerFaceColor','g','MarkerSize',3.5);
             
              [c_extr,r_extr,p] = impixel;
              
              JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = c_extr(1);
              JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = r_extr(1);
              JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = c_extr(2);
              JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = r_extr(2);
              JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
              JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;
              
          estimated_b_axis_length = (sqrt((JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg))^2 + (JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg))^2))/2;
          estimated_a_axis_length = estimated_b_axis_length/2;
          occu_area_ellipse(frame_count) = 3.14 * estimated_a_axis_length * estimated_b_axis_length;
               
              % occu_area_ellipse(frame_count) = occu_area_ellipse(frame_count);           
              % RC_occu_area_ellipse(frame_count) = RC_occu_area_ellipse(frame_count); 
              
              occu_dist(frame_count,gg) = 17;
              occluding_frames(frame_count,gg) = frame_count;
              save( OCCUL_INFO_FILE_NAME, 'occluding_frames');
            
          elseif EXTRA_OCCL_COND <=25
               
               JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = occu_third_frame_head_X_cord;
               JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = occu_third_frame_head_Y_cord;
               JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = occu_third_frame_tail_X_cord;
               JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = occu_third_frame_tail_Y_cord;
               JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
               JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

               current_occu_area_ellipse = occu_area_ellipse(frame_count-1);
               occu_area_ellipse(frame_count) = occu_area_ellipse(frame_count-1);           
               RC_occu_area_ellipse(frame_count) = RC_occu_area_ellipse(frame_count-1); 
               occu_dist(frame_count,gg) = 17;
               occluding_frames(frame_count,gg) = frame_count;
               save( OCCUL_INFO_FILE_NAME, 'occluding_frames');
               
          else
              
              disp('Check the thresholds');
             
          end
          
          
          else
               
          % when the bee image is appearing circular in shape
          
          disp('The selected bee is appearing circular in shape');
          
              circle_opt_2_X_cord = two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count) - one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count);
              circle_opt_2_Y_cord = two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count) - one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count);

              circle_opt_2_vect = [circle_opt_2_X_cord circle_opt_2_Y_cord];
              mag_opt_2_vect = rssq(circle_opt_2_vect);
              circle_unit_opt_2_vect = [circle_opt_2_vect(1)/mag_opt_2_vect circle_opt_2_vect(2)/mag_opt_2_vect];

              circle_opt_1_X_cord = one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count) - two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count);
              circle_opt_1_Y_cord = one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count) - two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count);


              circle_opt_1_vect = [circle_opt_1_X_cord circle_opt_1_Y_cord];
              mag_opt_1_vect = rssq(circle_opt_1_vect);
              circle_unit_opt_1_vect = [circle_opt_1_vect(1)/mag_opt_1_vect circle_opt_1_vect(2)/mag_opt_1_vect];


            circle_BODY_frame_X_1_ST_TIME = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg); 
            circle_BODY_frame_Y_1_ST_TIME = JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg); 


            circle_BODY_INFO_1_ST_TIME = [circle_BODY_frame_X_1_ST_TIME  circle_BODY_frame_Y_1_ST_TIME];
            circle_MAG_BODY_INFO_1_ST_TIME =  rssq(circle_BODY_INFO_1_ST_TIME);
            circle_BODY_INFO_1_ST_TIME_UNIT = [circle_BODY_INFO_1_ST_TIME(1)/circle_MAG_BODY_INFO_1_ST_TIME  circle_BODY_INFO_1_ST_TIME(2)/circle_MAG_BODY_INFO_1_ST_TIME];


                circle_possib_1_dot_BODY_INFO_1_ST_TIME = dot(circle_BODY_INFO_1_ST_TIME_UNIT, circle_unit_opt_1_vect);
                circle_possib_2_dot_BODY_INFO_1_ST_TIME = dot(circle_BODY_INFO_1_ST_TIME_UNIT, circle_unit_opt_2_vect); 

                circle_possib_1_theta_BODY_INFO_1_ST_TIME = acosd(circle_possib_1_dot_BODY_INFO_1_ST_TIME);
                circle_possib_2_theta_BODY_INFO_1_ST_TIME = acosd(circle_possib_2_dot_BODY_INFO_1_ST_TIME);
                    
                    
                    if circle_possib_1_theta_BODY_INFO_1_ST_TIME <= 40
                        
                        
                      JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count);
                      JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count);

                      JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count);
                      JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count);

                      JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
                      JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

                      fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
                      fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
                      fin_position(gg,3) = 1.5 ;
                      fin_label(gg,1) = gg; 

                      occu_dist(frame_count,gg) = 16;
                        
                        
                        
                    elseif circle_possib_2_theta_BODY_INFO_1_ST_TIME <=40
                        
                        
                      JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = two_prob_cord2_X(CORRECT_1_ST_TIME_IND,frame_count);
                      JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = two_prob_cord2_Y(CORRECT_1_ST_TIME_IND,frame_count);

                      JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = one_prob_cord1_X(CORRECT_1_ST_TIME_IND,frame_count);
                      JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = one_prob_cord1_Y(CORRECT_1_ST_TIME_IND,frame_count);

                      JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
                      JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

                      fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
                      fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
                      fin_position(gg,3) = 1.5 ;
                      fin_label(gg,1) = gg; 

                      occu_dist(frame_count,gg) = 16;
                                                
                        
                    else
                        
                      imshow(oImg);
                      fig = gcf;
                      fig.Position = [-1790 -3 1540 804];
                      hold on
                      plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg),'o','Color','m','MarkerFaceColor','m','MarkerSize',3.5);
                      hold on
                      plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg),'o','Color','g','MarkerFaceColor','g','MarkerSize',3.5)
                      hold on

                      [c_cir,r_cir,p] = impixel;

                      JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = c_cir(1);
                      JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = r_cir(1);

                      JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = c_cir(2);
                      JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = r_cir(2);

                      JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
                      JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;

                      fin_position(gg,1) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg); 
                      fin_position(gg,2) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
                      fin_position(gg,3) = 1.5 ;
                      fin_label(gg,1) = gg; 

                      occu_dist(frame_count,gg) = 17;
    
                    end
                    
          

           end  
        
       
        else
            
   
          imshow(oImg);
          fig = gcf;
          fig.Position = [-1790 -3 1540 804];
          hold on
          plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Head_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Head_Y_Cord(gg),'o','Color','m','MarkerFaceColor','m','MarkerSize',3.5);
          hold on
          plot(JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_X_Cord(gg), JUMP_HEAD_TAIL_CORD(frame_count-1).Tail_Y_Cord(gg),'o','Color','g','MarkerFaceColor','g','MarkerSize',3.5)
          hold on   
                        
        quest = sprintf('Check, if the bee has truely left the arena. If "yes", enter "y" or else "n"'); 
        str = input(quest,'s');
        
        % This condition helps to check if the bee in the arena or out of the
        % arena. This condition will be active when 'non_occ_CP_bee_ind22' is
        % empty
        
         if str == 'n'
                   
          [c_noise,r_noise,p] = impixel;
         
          JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) = c_noise(1);
          JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) = r_noise(1);

          JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg) = c_noise(2);
          JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg) = r_noise(2);

          JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg) = frame_count;
          JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg) = gg;
      
          estimated_b_axis_length = (sqrt((JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg))^2 + (JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg) - JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg))^2))/2;
          estimated_a_axis_length = estimated_b_axis_length/2;
          occu_area_ellipse(frame_count) = 3.14 * estimated_a_axis_length * estimated_b_axis_length;
          occu_dist(frame_count,gg) = 17;
         
         else
           
                    occu_dist(frame_count,gg) = -1;
                    break;
        end
       

       
  end
           
       end
        
       
   end
    
       
HEAD_X(frame_count,gg) = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg);
HEAD_Y(frame_count,gg) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
TAIL_X(frame_count,gg) = JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg);
TAIL_Y(frame_count,gg) = JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg);
FRAME_NUM(frame_count,gg) = JUMP_HEAD_TAIL_CORD(frame_count).Frame_number(gg);
BEE_NUM(frame_count,gg) = JUMP_HEAD_TAIL_CORD(frame_count).Bee_number(gg);   

pres_frame_HEAD_X_comp(frame_count,gg)  = JUMP_HEAD_TAIL_CORD(frame_count).Head_X_Cord(gg);
pres_frame_HEAD_Y_comp(frame_count,gg) = JUMP_HEAD_TAIL_CORD(frame_count).Head_Y_Cord(gg);
pres_frame_TAIL_X_comp(frame_count,gg) = JUMP_HEAD_TAIL_CORD(frame_count).Tail_X_Cord(gg);
pres_frame_TAIL_Y_comp(frame_count,gg) = JUMP_HEAD_TAIL_CORD(frame_count).Tail_Y_Cord(gg);

num_bee_prev_frame(frame_count) = num_bee_present_frame(frame_count);

save( JUMP_HT_CORDS_FILE_NAME, 'JUMP_HEAD_TAIL_CORD');
save( ALL_FRAMES_INFO_FILE_NAME, 'FRAME_NUM');
save( CIRCLE_INFO_FILE_NAME, 'ratio_A_by_B');

dist_1 = [];
dist_2 = [];
temp_theta = [];
third_frame_dist_1 = [];
third_frame_dist_2 = [];
occu_third_frame_dist_1 = [];
occu_third_frame_dist_2 = [];
combined_occu_third_fram_dist = [];
combined_third_fram_dist = [];

% Condition to check if the bee is in the arena or not
            
      if  HEAD_Y(frame_count,gg) >= 2  && HEAD_Y(frame_count,gg) <= 1080 && HEAD_X(frame_count,gg) <= 2 % left wall

           disp('The selected bee has left the arena');
           occu_dist(frame_count,gg) = -1;
           break;

      elseif HEAD_Y(frame_count,gg) >= 1078  && HEAD_Y(frame_count,gg) <= 1080 && HEAD_X(frame_count,gg) > 0 && HEAD_X(frame_count,gg) <= 2048 % straight bottom wall 

           disp('The selected bee has left the arena');
           occu_dist(frame_count,gg) = -1;
           break;

      elseif HEAD_Y(frame_count,gg) >= 0 && HEAD_Y(frame_count,gg) <= 1080 && HEAD_X(frame_count,gg) > 2046 && HEAD_X(frame_count,gg) <= 2048 % right wall

          disp('The selected bee has left the arena');
          occu_dist(frame_count,gg) = -1;
          break;
          
      elseif HEAD_Y(frame_count,gg) > 0  && HEAD_Y(frame_count,gg) <= 1.5 && HEAD_X(frame_count,gg) > 0 && HEAD_X(frame_count,gg) <= 2048 % straight top wall
          
          disp('The selected bee has left the arena');
          occu_dist(frame_count,gg) = -1;
          break;   

              
      end
 
      
              
end
fin_position = [];
fin_label = [];

frame_count

if  occu_dist(frame_count,gg) == -1 || frame_count == last_frame_number 
    
 break;   
   
end


end

%%


% do the interpolation operation

if ~isempty(occluding_frames)
    
to_be_interpo_occcu_frames = occluding_frames(:,gg)'; 
cropped_to_be_interpo_occcu_frames =  to_be_interpo_occcu_frames(to_be_interpo_occcu_frames ~=0); 

all_frame_nos = FRAME_NUM(:,gg)';
cropped_all_frame_nos = all_frame_nos(all_frame_nos ~=0);

uninterp_extract_cropped_frames = setdiff(cropped_all_frame_nos,cropped_to_be_interpo_occcu_frames);


for dd = 1:length(uninterp_extract_cropped_frames)
    
    frame_no_temp = uninterp_extract_cropped_frames(dd);
    
    new_cropped_head_X(1,dd) = JUMP_HEAD_TAIL_CORD(frame_no_temp).Head_X_Cord(gg);
    new_cropped_head_Y(1,dd) = JUMP_HEAD_TAIL_CORD(frame_no_temp).Head_Y_Cord(gg);
    new_cropped_tail_X(1,dd) = JUMP_HEAD_TAIL_CORD(frame_no_temp).Tail_X_Cord(gg);
    new_cropped_tail_Y(1,dd) = JUMP_HEAD_TAIL_CORD(frame_no_temp).Tail_Y_Cord(gg);
    
end
 
hx = uninterp_extract_cropped_frames;
hv = new_cropped_head_X;
hxq = cropped_to_be_interpo_occcu_frames;
hvq = interp1(hx,hv,hxq,'cubic'); 

hy = uninterp_extract_cropped_frames;
hu = new_cropped_head_Y;
hyp = cropped_to_be_interpo_occcu_frames;
huq = interp1(hy,hu,hyp,'cubic');

tail_x = uninterp_extract_cropped_frames;
tail_v = new_cropped_tail_X;
tail_xq = cropped_to_be_interpo_occcu_frames;
tail_vq = interp1(tail_x,tail_v,tail_xq,'cubic');

tail_y = uninterp_extract_cropped_frames;
tail_u = new_cropped_tail_Y;
tail_yq = cropped_to_be_interpo_occcu_frames;
tail_uq = interp1(tail_y,tail_u,tail_yq,'cubic');
 

f11 = figure;
f11.Position = [418 284 1228 675];
figure(f11),
imshow(oImg);
hold on
plot(new_cropped_head_X,new_cropped_head_Y,'-o','Color','w','MarkerSize',6.5);
hold on
plot(hvq,huq,'o','Color','g','LineWidth',1.5,'MarkerFaceColor','g','MarkerSize',5.5);
hold on
plot(new_cropped_tail_X,new_cropped_tail_Y,'Color','r')
hold on
plot(tail_vq,tail_uq,'o','Color','b','LineWidth',1.5,'MarkerFaceColor','b','MarkerSize',5.5);
hold off
 
for tt = 1:length(cropped_to_be_interpo_occcu_frames)
    
updated_cord_frame = cropped_to_be_interpo_occcu_frames(tt);
    
JUMP_HEAD_TAIL_CORD(updated_cord_frame).Head_X_Cord(gg) = hvq(tt);
JUMP_HEAD_TAIL_CORD(updated_cord_frame).Head_Y_Cord(gg) = huq(tt);
JUMP_HEAD_TAIL_CORD(updated_cord_frame).Tail_X_Cord(gg) = tail_vq(tt);
JUMP_HEAD_TAIL_CORD(updated_cord_frame).Tail_Y_Cord(gg) = tail_uq(tt);

    
end


save( JUMP_HT_CORDS_FILE_NAME, 'JUMP_HEAD_TAIL_CORD');

% start frame and end frame identification to plot the head and tail of the
% final plot

starting_frame_of_bee = min(cropped_all_frame_nos);
last_frame_of_bee = max(cropped_all_frame_nos);

for yy = starting_frame_of_bee:last_frame_of_bee
    
  
    head_X(yy) = JUMP_HEAD_TAIL_CORD(yy).Head_X_Cord(gg); tail_X(yy) =   JUMP_HEAD_TAIL_CORD(yy).Tail_X_Cord(gg);
    head_Y(yy) = JUMP_HEAD_TAIL_CORD(yy).Head_Y_Cord(gg); tail_Y(yy) =   JUMP_HEAD_TAIL_CORD(yy).Tail_Y_Cord(gg);
   
    body_len(yy) = sqrt((head_X(yy) - tail_X(yy))^2+ (head_Y(yy) - tail_Y(yy))^2);

    unit_body_X_cord(yy) = (tail_X(yy) - head_X(yy)) /  body_len(yy) ;
    unit_body_Y_cord(yy) = (tail_Y(yy) - head_Y(yy)) /  body_len(yy);
   
end

f12 = figure;
f12.Position = [418 284 1228 675];
figure(f12),
imshow(oImg);
hold on
plot(head_X,head_Y,'o','Color','g','MarkerSize',5.5,'LineWidth',1.5,'MarkerFaceColor','g');
hold on   
plot(tail_X,tail_Y,'o','Color','b','MarkerSize',2,'LineWidth',1.5,'MarkerFaceColor','b');
hold on
scale =0.2;
quiver(head_X,head_Y,unit_body_X_cord,unit_body_Y_cord,scale,'color','w','LineStyle','-','LineWidth',1.85,'ShowArrowHead','off');
% hold on

else
    
            all_frame_nos = FRAME_NUM(:,gg)';
            cropped_all_frame_nos = all_frame_nos(all_frame_nos ~=0);

            starting_frame_of_bee = min(cropped_all_frame_nos);
            last_frame_of_bee = max(cropped_all_frame_nos);

            for yy = starting_frame_of_bee:last_frame_of_bee


                head_X(yy) = JUMP_HEAD_TAIL_CORD(yy).Head_X_Cord(gg); tail_X(yy) =   JUMP_HEAD_TAIL_CORD(yy).Tail_X_Cord(gg);
                head_Y(yy) = JUMP_HEAD_TAIL_CORD(yy).Head_Y_Cord(gg); tail_Y(yy) =   JUMP_HEAD_TAIL_CORD(yy).Tail_Y_Cord(gg);

                body_len(yy) = sqrt((head_X(yy) - tail_X(yy))^2+ (head_Y(yy) - tail_Y(yy))^2);

                unit_body_X_cord(yy) = (tail_X(yy) - head_X(yy)) /  body_len(yy) ;
                unit_body_Y_cord(yy) = (tail_Y(yy) - head_Y(yy)) /  body_len(yy);

            end

            f12 = figure;
            f12.Position = [418 284 1228 675];
            figure(f12),
            imshow(oImg);
            hold on
            plot(head_X,head_Y,'o','Color','g','MarkerSize',5.5,'LineWidth',1.5,'MarkerFaceColor','g');
            hold on   
            plot(tail_X,tail_Y,'o','Color','b','MarkerSize',2,'LineWidth',1.5,'MarkerFaceColor','b');
            hold on
            scale =0.2;
            quiver(head_X,head_Y,unit_body_X_cord,unit_body_Y_cord,scale,'color','w','LineStyle','-','LineWidth',1.85,'ShowArrowHead','off');
            % hold on

end


% for points where the bee appears like a dot


starting_frame_of_bee_for_ratios = min(cropped_all_frame_nos);
last_frame_of_bee_for_ratios = max(cropped_all_frame_nos);

circle_count = 0;

for xx = starting_frame_of_bee_for_ratios:last_frame_of_bee_for_ratios
    
    if ( ratio_A_by_B(xx) >= 0.9 && ratio_A_by_B(xx) <= 1.3)

circle_count = circle_count+1;        
circular_bee_frames(circle_count) = xx;

JUMP_HEAD_TAIL_CORD(xx).Tail_X_Cord(gg) = JUMP_HEAD_TAIL_CORD(xx).Head_X_Cord(gg);
JUMP_HEAD_TAIL_CORD(xx).Tail_Y_Cord(gg) = JUMP_HEAD_TAIL_CORD(xx).Head_Y_Cord(gg) ;
        
plot(JUMP_HEAD_TAIL_CORD(xx).Tail_X_Cord(gg) ,JUMP_HEAD_TAIL_CORD(xx).Tail_Y_Cord(gg) ,'o','Color','m','MarkerSize',5.5,'LineWidth',1.5,'MarkerFaceColor','m');
hold on
    
    else
        
        
    end
   
end

save( JUMP_HT_CORDS_FILE_NAME, 'JUMP_HEAD_TAIL_CORD');


%%



