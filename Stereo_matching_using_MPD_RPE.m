% This is a program for stereomatching bees in the two camera views. 

% please see the readme file for complete usage information :
% https://github.com/mahadeesh/Bee-tracking-and-stereomatching-files 
	
% Inputs variable description: 
%------------------------------
% 1) (Ref_bee_no)  - unique bee IDs of bees flying in the left camera view
% 2) (opt_bee_numbers) - unique bee IDs of bees flying in the right camera view.
% 3) (part and part2) - body part of a bee. part is for head and part2 is for tail
% 4) (NO_REF_PTS) - predefined reference points for our experimental condition
% 5) (DATA_INFO) - Date of the data captured
% 6) (stereoParams_16_FEB_2017.mat) - stereo camera calibration file generated using Matlab's stereo calibration app.
% 7) (BC_stereo_cam_calib_trail_3_.mat) - stereo camera calibration file generated using Yves bouguet stereo calibration app. please see the readme file
% 8) myfilename4 - path of the 2D trajectory of the bee from right camera view(reference view), for which you would like to pick a matching pair from left camera view. 
% 9) myfilename3 - path of the 2D tracjectories of all bees flying the left camera view


% Output variable description: 
%-------------------------------
% 1) (MIN_DIST_REPROJ_ERR_INFO.mat)  - The final variable which will contain all minimum perpendicular distances and reprojection errors information
% 2) (MATCHING_BEE_MIN_DIST_REPROJ_ERR) - The correct matching bee's ID, minimum perpendicular distance value and reprojection error value

% See read the associated readme file once before running this program

% Sample inputs are also provided at https://github.com/mahadeesh/Bee-tracking-and-stereomatching-files 

% If there are questions in relation to the operation of this software,
% please contact: Mandiyam Mahadeeswara at m.mandiyam@uq.edu.au or mdevaraj@outlook.com

clear all;
close all;
clc;

Ref_bee_no =  [10 11 12 17 18 19];  
opt_bee_numbers = [1:6 8:16 18:99 101:103]; 

obj_no1 = 1; obj_no2=1;
part = 'HEAD'; 
part2 = 'TAIL';
NO_REF_PTS = 8;
DATA_INFO = '16 FEB 2017'; 
load('FEB_16_2017_STEREO_PARAMS.mat');
load('BC_stereo_cam_calib_trail_3_.mat');


field_name_1 = 'INDI_MIN_DIST'; field_name_2 = 'OPT_BEE_NUM'; field_name_3 = 'STD_VALUES';  
field_name_4 = 'INDI_RPE_ERR'; field_name_5 = 'RPE_OPT_BEE_NUM'; field_name_6 = 'RPE_STD_VALUES'; field_name_7 = 'FRA_LEN';

MIN_DIST_REPROJ_ERR_INFO = struct(field_name_1,[],field_name_2,[],field_name_3,[],field_name_4,[],field_name_5,[],field_name_6,[],field_name_7,[]);

MIN_DIST_REPROJ_ERR_INFO_FILE_NAME = [sprintf('MIN_DIST_REPROJ_ERR_INFO_%s',DATA_INFO) '.mat'];
% load(MIN_DIST_REPROJ_ERR_INFO_FILE_NAME);


for mm = 1:length(Ref_bee_no)
    
beestart_no = Ref_bee_no(mm); 
   
% a standard reference points for in our experiments and will not affect
% new set of inputs from the 

myfilename1 = sprintf('BC_16_FEB_LEFT_REF_POINTS_CAM1_BKP.txt');

left_view_ref_data = importdata(myfilename1); 
  
left_ref_value_x = left_view_ref_data(:,1);
left_ref_value_y = left_view_ref_data(:,2);


myfilename2 = sprintf('BC_16_FEB_RIGHT_REF_POINTS_CAM2_BKP.txt');

right_view_ref_data = importdata(myfilename2); 

right_ref_value_x = right_view_ref_data(:,1);
right_ref_value_y = right_view_ref_data(:,2);

left_ref_len = length(left_view_ref_data(:,1));
right_ref_len = length(right_view_ref_data(:,1));

accounted_bees =0;
 

for tt = 1:length(opt_bee_numbers)

 left_ref_len = 8;
 
 opt_act_num = opt_bee_numbers(tt);
 
% importing 2D position data of bees flying in the LEFT camera view
% provide the path of the sample inputs downloaded in your local foler

% myfilename3 = sprintf('ANALYSIS/DATA FROM SHARED DRIVE FROM SRINI 4 JUNE 2021/2D TRACKS RECONSTRUCTION/DS1_OTHER_BEES_2D_TRACK_CAM3/J_BC_HT_CORDS_CAM_3_%s_BEE_NUMBER_%d_f100.mat',DATA_INFO,opt_act_num);
myfilename3 = sprintf('ANALYSIS/WRITING/MANUSCRIPT 2/SCI DATA PROGRAMS FOR UPLOAD/SAMPLE_INPUTS_FOR_STEREOMATCHING_PROGRAM/DS1_OTHER_BEES_2D_TRACK_LEFT_CAM_VIEW/J_BC_HT_CORDS_CAM_3_%s_BEE_NUMBER_%d_f100.mat',DATA_INFO,opt_act_num);


LEFT_view_TRACK_data = load(myfilename3); 

LEFT_VIEW_DATA_LEN = length(LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD);

 for aa = 1:LEFT_VIEW_DATA_LEN
     
     if ~isempty(LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(aa).Head_X_Cord)
        
      opt_start_pos = (LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(aa).Frame_number(opt_act_num));
         
      break;
         
     end
    
 end

opt_end_pos = (LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(end).Frame_number(opt_act_num));

opt_trk_ct =0;
 
 for zz = opt_start_pos:opt_end_pos
  
     opt_trk_ct = opt_trk_ct+1;
     
    left_first_frame = LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(1).Frame_number(opt_act_num);
 
    left_new_yy = zz - left_first_frame + 1 ;   
     
    LEFT_TRAC_value_x(tt,opt_trk_ct) = LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(left_new_yy).Head_X_Cord(opt_act_num);
    LEFT_TRAC_value_y(tt,opt_trk_ct) = LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(left_new_yy).Head_Y_Cord(opt_act_num);
     
 end

OPT_VIEW_FIRST_FRAME = opt_start_pos;
OPT_VIEW_LAST_FRAME = opt_end_pos;


% importing 2D position data of bees flying in the RIGHT camera view 

% provide the path of the sample inputs downloaded in your local folder

% myfilename4 = sprintf('ANALYSIS/DATA FROM SHARED DRIVE FROM SRINI 4 JUNE 2021/2D TRACKS RECONSTRUCTION/DS1_REFERENCE_BEE_2D_TRACK_CAM2/J_BC_HT_CORDS_CAM_2_%s_BEE_NUMBER_BEE_%d_f1.mat',DATA_INFO,beestart_no);
myfilename4 = sprintf('ANALYSIS/WRITING/MANUSCRIPT 2/SCI DATA PROGRAMS FOR UPLOAD/SAMPLE_INPUTS_FOR_STEREOMATCHING_PROGRAM/DS1_REFERENCE_BEE_2D_TRACK_RIGHT_CAM_VIEW/J_BC_HT_CORDS_CAM_2_%s_BEE_NUMBER_BEE_%d_f1.mat',DATA_INFO,beestart_no);


RIGHT_view_TRACK_data = load(myfilename4);
RIGHT_VIEW_DATA_LEN = length(RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD);

 for bb = 1:RIGHT_VIEW_DATA_LEN
     
     if ~isempty(RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(bb).Head_X_Cord)
         
       
      ref_start_pos = (RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(bb).Frame_number(beestart_no));
         
      break;
         
     end
     
     
 end
 
 ref_end_pos = (RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(end).Frame_number(beestart_no));
 
 rt_trk_ct =0;
 
 for yy = ref_start_pos:ref_end_pos
  
     rt_trk_ct = rt_trk_ct+1;
     
    first_frame = RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(1).Frame_number(beestart_no);
 
    new_yy = yy - first_frame + 1 ;   
     
    RIGHT_TRAC_value_x(tt,rt_trk_ct) = RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(new_yy).Head_X_Cord(beestart_no);
    RIGHT_TRAC_value_y(tt,rt_trk_ct) = RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(new_yy).Head_Y_Cord(beestart_no);
    
    REF_frameslot(tt,rt_trk_ct) = yy;
     
 end

RIGHT_VIEW_FIRST_FRAME = ref_start_pos;
RIGHT_VIEW_LAST_FRAME = ref_end_pos;

REF_pos_len_info(tt) = rt_trk_ct;


% picking frame starting and ending point of each bee

frame_start_end_info(1,tt) = opt_end_pos;
frame_start_end_info(2,tt) = ref_end_pos;
frame_start_end_info(3,tt) = opt_act_num;

close_end_frame_info(1,tt) = abs(ref_end_pos - opt_end_pos);
close_end_frame_info(2,tt) = opt_act_num;

start_frame_info(1,tt) = abs(ref_start_pos - opt_start_pos);
start_frame_info(2,tt) = opt_act_num;

% reference data formation

  for mm = 1: left_ref_len

new_left_track_dat(mm,obj_no1,1) = left_ref_value_x(mm);
new_left_track_dat(mm,obj_no1,2) = left_ref_value_y(mm);
new_left_track_dat(mm,obj_no1,3) = mm;
new_right_track_dat(mm,obj_no2,1) = right_ref_value_x(mm);
new_right_track_dat(mm,obj_no2,2) = right_ref_value_y(mm);
new_right_track_dat(mm,obj_no2,3) = mm;


  end
  
  
 first_frame_info = [OPT_VIEW_FIRST_FRAME RIGHT_VIEW_FIRST_FRAME]; 
 FINAL_STARTING_FIRST_FRAME  = max(first_frame_info); 
 
 last_frame_info = [OPT_VIEW_LAST_FRAME RIGHT_VIEW_LAST_FRAME];
 FINAL_ENDING_LAST_FRAME  = min(last_frame_info); 
 
 ABS_FRAME_DIFF = FINAL_ENDING_LAST_FRAME - FINAL_STARTING_FIRST_FRAME;
 

 if (FINAL_ENDING_LAST_FRAME > FINAL_STARTING_FIRST_FRAME) && (ABS_FRAME_DIFF > 1)
 
accounted_bees = accounted_bees+1;
HP = left_ref_len;
head_count = 0;

for fframe =  FINAL_STARTING_FIRST_FRAME:FINAL_ENDING_LAST_FRAME
            
    left_ref_len = left_ref_len+1;
    
    HP = HP+1;
    
    head_count = head_count+1;
    
    lef_first_frame = LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(1).Frame_number(opt_act_num);
    lef_new_yy = fframe - lef_first_frame + 1 ; 
    
    new_left_track_dat(HP,obj_no1,1) = LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(lef_new_yy).Head_X_Cord(opt_act_num); 
    new_left_track_dat(HP,obj_no1,2) = LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(lef_new_yy).Head_Y_Cord(opt_act_num); 
    new_left_track_dat(HP,obj_no1,3) = fframe;

    right_first_frame = RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(1).Frame_number(beestart_no);
    rig_new_yy = fframe - right_first_frame + 1 ; 
    
    new_right_track_dat(HP,obj_no1,1) = RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(rig_new_yy).Head_X_Cord(beestart_no);
    new_right_track_dat(HP,obj_no1,2) = RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(rig_new_yy).Head_Y_Cord(beestart_no);
    new_right_track_dat(HP,obj_no1,3) = fframe;
        
    UNS_OPT_HEAD_X(tt,head_count) = new_left_track_dat(HP,obj_no1,1);
    UNS_OPT_HEAD_Y(tt,head_count) = new_left_track_dat(HP,obj_no1,2);
    
    UNS_RIGHT_TRAC_value_x(tt,head_count) = new_right_track_dat(HP,obj_no1,1);
    UNS_RIGHT_TRAC_value_y(tt,head_count) = new_right_track_dat(HP,obj_no1,2);
    
    frameslot(tt,head_count) = fframe;
            
    TP = HP+1;
    
    new_left_track_dat(TP,obj_no1,1) = LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(lef_new_yy).Tail_X_Cord(opt_act_num);
    new_left_track_dat(TP,obj_no1,2) = LEFT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(lef_new_yy).Tail_Y_Cord(opt_act_num);
    new_left_track_dat(TP,obj_no1,3) = fframe;
    
    new_right_track_dat(TP,obj_no1,1) = RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(rig_new_yy).Tail_X_Cord(beestart_no);
    new_right_track_dat(TP,obj_no1,2) = RIGHT_view_TRACK_data.JUMP_HEAD_TAIL_CORD(rig_new_yy).Tail_Y_Cord(beestart_no);
    new_right_track_dat(TP,obj_no1,3) = fframe;
    
    HP = TP;
    
end
  
pos_len_info(tt) = head_count;

OPT_HEAD_X(tt,1:pos_len_info(tt)) = smooth(UNS_OPT_HEAD_X(tt,1:head_count),29);
OPT_HEAD_Y(tt,1:pos_len_info(tt)) = smooth(UNS_OPT_HEAD_Y(tt,1:head_count),29);

SMO_RIGHT_TRAC_value_x(tt,1:pos_len_info(tt)) = smooth(UNS_RIGHT_TRAC_value_x(tt,1:head_count),29);
SMO_RIGHT_TRAC_value_y(tt,1:pos_len_info(tt)) = smooth(UNS_RIGHT_TRAC_value_y(tt,1:head_count),29);

for hh = 1:pos_len_info(tt)
    
    LOrig = [OPT_HEAD_X(tt,hh);OPT_HEAD_Y(tt,hh)];
    ROrig = [SMO_RIGHT_TRAC_value_x(tt,hh);SMO_RIGHT_TRAC_value_y(tt,hh)];
        
% computing the positions of the two cameras    

[LOrigMM, ROrigMM,C2_POS_WOR] = stereo_triangulation_camera_positions(LOrig, ROrig, rot_est, trans_est, fc_l, cc_l, kc_l, ac_l, fc_r, cc_r, kc_r, ac_r);

CAM1_VECT = [LOrigMM(1) LOrigMM(2) LOrigMM(3)]; 
CAM2_VECT =  [ROrigMM(1) ROrigMM(2) ROrigMM(3)];
C2_POS = [C2_POS_WOR(1) C2_POS_WOR(2) C2_POS_WOR(3)];
CAM1 = [0 0 0];

scale = 1.5;scale2 =1.5;
SCA_CAM1_VECT = scale*CAM1_VECT;
SCA_CAM2_VECT = scale2*CAM2_VECT;

NEW_SCALED_CAM2_VECT = [(C2_POS(1) + SCA_CAM2_VECT(1))  (C2_POS(2) + SCA_CAM2_VECT(2)) (C2_POS(3) + SCA_CAM2_VECT(3))]; % wrt to cam 2
REAL_CAM2_VECT = [(C2_POS(1) + CAM2_VECT(1))  (C2_POS(2) + CAM2_VECT(2)) (C2_POS(3) + CAM2_VECT(3))];


DIST_XL_XR(tt,hh) = sqrt((CAM1_VECT(1)-REAL_CAM2_VECT(1))^2 + (CAM1_VECT(2)-REAL_CAM2_VECT(2))^2 + (CAM1_VECT(3)-REAL_CAM2_VECT(3))^2);

% minimum perpendicular distance calculation

% line 1 details

a1 =  CAM1;  lamda = 1;
b1 = (SCA_CAM1_VECT - CAM1);
line1_eqn = a1 + lamda*b1;

% line 2 details

a2 = C2_POS; mew = 1;
b2 = (NEW_SCALED_CAM2_VECT - C2_POS); % 
line2_eqn =  a2 + mew*b2;

% minimum perpendicular distance estimation

perp_unit_vect_n = cross(b1,b2) ./ norm(cross(b1,b2)) ;
any_two_pt_vect = (NEW_SCALED_CAM2_VECT - SCA_CAM1_VECT);
short_dist_mpd = norm(dot(any_two_pt_vect,perp_unit_vect_n));

short_dist_using_mpd_method(tt,hh) = short_dist_mpd;

% reprojection error estimation

    MP1 = [OPT_HEAD_X(tt,hh) OPT_HEAD_Y(tt,hh)];
    MP2 = [SMO_RIGHT_TRAC_value_x(tt,hh) SMO_RIGHT_TRAC_value_y(tt,hh)];
    
[wrdpt_aft,reproj_aft(hh)] = triangulate(MP1,MP2,stereoParams_16_FEB_2017);

end

% each bee's mean MPD 

mean_short_dist_using_mpd(1,tt) = mean(short_dist_using_mpd_method(tt,1:hh)); % mean shortest distance
mean_short_dist_using_mpd(2,tt) = opt_act_num; % bee num
mean_short_dist_using_mpd(3,tt) = std(short_dist_using_mpd_method(tt,1:hh),0,2); % std value

F_start = frameslot(tt,1);
F_end =  frameslot(tt,head_count);

FRAME_POS = F_start:F_end;

FRAME_START_END(1,tt) = F_start;
FRAME_START_END(2,tt) = F_end;

% each bee's mean RPE

mean_reproj_error(1,tt) = mean(reproj_aft(1:hh));
mean_reproj_error(2,tt) = opt_act_num;
mean_reproj_error(3,tt) = std(reproj_aft(1:hh),0,2);

% all MPD and RPE values are stored here

MIN_DIST_REPROJ_ERR_INFO(beestart_no).INDI_MIN_DIST(opt_act_num,1:hh) = short_dist_using_mpd_method(tt,1:hh);
MIN_DIST_REPROJ_ERR_INFO(beestart_no).INDI_RPE_ERR(opt_act_num,1:hh) = reproj_aft(1:hh);

MIN_DIST_REPROJ_ERR_INFO(beestart_no).OPT_BEE_NUM(opt_act_num,1) = opt_act_num;
MIN_DIST_REPROJ_ERR_INFO(beestart_no).RPE_OPT_BEE_NUM(opt_act_num,1) = opt_act_num;

MIN_DIST_REPROJ_ERR_INFO(beestart_no).STD_VALUES(opt_act_num,1) = mean_short_dist_using_mpd(3,tt);
MIN_DIST_REPROJ_ERR_INFO(beestart_no).RPE_STD_VALUES(opt_act_num,1) = mean_reproj_error(3,tt);

MIN_DIST_REPROJ_ERR_INFO(beestart_no).FRA_LEN(opt_act_num,1) = pos_len_info(tt);

%****************************************************************************************
clear new_left_track_dat new_right_track_dat reproj_aft; 

 else

clear new_left_track_dat new_right_track_dat ;

end

end
%%
if accounted_bees >= 1
    
    % Picking the bee with minimum perpendicular distance value
 
    [A,B] = sort(mean_short_dist_using_mpd(1,:));
    
    for hh = 1:length(B)
        
    sorted_corres_BEE(1,hh) = A(hh);
    sorted_corres_BEE(2,hh) = opt_bee_numbers(1,B(hh));
    sorted_corres_BEE(3,hh) = mean_short_dist_using_mpd(3,B(hh));
    
    end
    sorted_corres_BEE;
    
sorted_bee_ind = find(sorted_corres_BEE(1,:) > 0);

VALID_MIN_DIST_INFO(beestart_no).MIN_DIST = sorted_corres_BEE(1,sorted_bee_ind(1):sorted_bee_ind(end));
VALID_MIN_DIST_INFO(beestart_no).OPT_BEE_NUM = sorted_corres_BEE(2,sorted_bee_ind(1):sorted_bee_ind(end));
VALID_MIN_DIST_INFO(beestart_no).STD_VALUES = sorted_corres_BEE(3,sorted_bee_ind(1):sorted_bee_ind(end)); 

MIN_DST_BEE(1,1) = sorted_corres_BEE(1,sorted_bee_ind(1));
MIN_DST_BEE(2,1) =  sorted_corres_BEE(2,sorted_bee_ind(1));
MIN_DST_BEE(1,2) = sorted_corres_BEE(1,sorted_bee_ind(2));
MIN_DST_BEE(2,2) =  sorted_corres_BEE(2,sorted_bee_ind(2));

% Picking the bee with minimum RPE value

    [A_reproj,B_reproj] = sort(mean_reproj_error(1,:));
    
    for hh = 1:length(B_reproj)
        
    sorted_corres_BEE_reproj(1,hh) = A_reproj(hh);
    sorted_corres_BEE_reproj(2,hh) = opt_bee_numbers(1,B_reproj(hh));
    sorted_corres_BEE_reproj(3,hh) = mean_reproj_error(3,B_reproj(hh));
    
    end
     
sorted_bee_ind_RPE = find(sorted_corres_BEE_reproj(1,:) > 0);

VALID_RPE_INFO(beestart_no).RPE_VALUE = sorted_corres_BEE_reproj(1,sorted_bee_ind_RPE(1):sorted_bee_ind_RPE(end));
VALID_RPE_INFO(beestart_no).OPT_BEE_NUM = sorted_corres_BEE_reproj(2,sorted_bee_ind_RPE(1):sorted_bee_ind_RPE(end));
VALID_RPE_INFO(beestart_no).STD_VALUES = sorted_corres_BEE_reproj(3,sorted_bee_ind_RPE(1):sorted_bee_ind_RPE(end)); 

RPE_ERR_VAL_BEE(1,1) = sorted_corres_BEE_reproj(1,sorted_bee_ind_RPE(1));
RPE_ERR_VAL_BEE(2,1) =  sorted_corres_BEE_reproj(2,sorted_bee_ind_RPE(1));
RPE_ERR_VAL_BEE(1,2) = sorted_corres_BEE_reproj(1,sorted_bee_ind_RPE(2));
RPE_ERR_VAL_BEE(2,2) =  sorted_corres_BEE_reproj(2,sorted_bee_ind_RPE(2));

% matching bee's ID, minimum distance value and reprojection error value will be stored 

MATCHING_BEE_MIN_DIST_REPROJ_ERR(1,beestart_no) = MIN_DST_BEE(1,1);
MATCHING_BEE_MIN_DIST_REPROJ_ERR(2,beestart_no) = MIN_DST_BEE(2,1);
MATCHING_BEE_MIN_DIST_REPROJ_ERR(3,beestart_no) = RPE_ERR_VAL_BEE(1,1);
MATCHING_BEE_MIN_DIST_REPROJ_ERR(4,beestart_no) = RPE_ERR_VAL_BEE(2,1);

save('MATCHING_BEE_MIN_DIST_REPROJ_ERR.mat');


disp('--------------DONE--------------------------------');

clear UNS_OPT_HEAD_X UNS_OPT_HEAD_Y UNS_RIGHT_TRAC_value_x UNS_RIGHT_TRAC_value_y OPT_HEAD_X OPT_HEAD_Y SMO_RIGHT_TRAC_value_x SMO_RIGHT_TRAC_value_y mean_reproj_error reproj_aft short_dist_from_mmy_pgm mean_short_dist_from_mmy_pgm sorted_corres_BEE;
    
end


end







 


