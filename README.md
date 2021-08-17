# Honeybee tracking and stereomatching program for reconstructing 3D trajectories of freely flying honeybees


## Introduction

Here are the files which I developed to generate time resolved or instantaneous 3D positions of bees flying in a bee cloud. The bee tracking program is a semi-automatic tracker of freely flying honeybees, which was developed using simple techniques such background subtraction, image differencing, interpolation, and extrapolation, to name a few. Furthermore, the problem of correspondence matching was solved in a simple way by using two different techniques namely minimum perpendicular distance method and minimum reprojection error method. Using the developed tracker and correspondence matching program, two large bee cloud datasets were digitised and could be accessed from here.

To know more about a bee cloud (which is different from a bee swarm), please see our previous work<sup>[1](https://www.nature.com/articles/s41598-018-35307-5)</sup>.  


## Highlights of the developed bee tracking and stereomatching program

1.	Despite being a semi-automatic bee tracking program, the final generated 3D position data contains head as well as tail locations of freely flying bees which makes it a unique and valuable dataset.

2.	The developed stereomatching code solves the correspondence matching problem just by using two cameras. 

3. The developed tracking program can be applied to track bees flying in a linear and curved tunnel like experimental setup.


## Prerequisites

All you need is MATLAB from MathWorks<sup>2</sup>. I have tested the code using MATLAB R2020 update 5 running on a macOS Catalina. All function files as given in the branch must be present.

## Contents of the files

1)	*Sci_data_bee_tracking_pgm_3.m* - The main bee tracking program which tracks individual bees flying in the bee cloud arena.

2)	*Stereo_matching_using_MPD_RPE.m* – a stereomatching code which matches every individual bee image in one camera view with its corresponding bee image in the other camera view. 

3)	*stereo_triangulation_camera_positions.m* – is used to calculate the positions of the two cameras with respect to the origin.

4)	*fit_ellipse.m* – is an open-source ellipse fitting program.


## How to perform image differencing operation using *ImageJ*

1.	Drag and drop the raw video in ImageJ (variable name – ‘v’ in the code)

2.	Go to Image option in the main menu >> stacks >> Z project >> select projection type as ‘median’ and click ‘OK’ >> image with only a background scene will be obtained (variable name – ‘ov’ in the code).

3.	Go to process in the main menu >> image calculator >> select operation as ‘difference’.

4.	From Image1 drop down list box, select the file name of the raw input video. Similarly, for Image2, select the name of the file obtained in step 2 and click ‘OK’.

5.	The resulting final video will eventually appear which will have to be saved using the ‘save as’ option from the main menu. Then, select ‘avi’ from the list of options which will save the final video in the destination folder as chosen by the user. 


## How to run the bee tracking program

1. Download the repository into a local folder

2. To run the bee tracking program, simply open the file *Sci_data_bee_tracking_pgm_3.m* in Matlab and click on the run button. Before clicking on the run command, make sure all necessary input videos/files are properly placed in the local repository folder. 

You can download all input video/files necessary to run the tracking program from [here](https://figshare.com/articles/media/Multi-Object_Tracking_in_Heterogeneous_environments_MOTHe_for_animal_video_recordings/11980356/3). Place all the downloaded videos/files in the folder where you have downloaded the repository.  

As this is a semi-automatic tracking program, the program prompts the user for inputs under certain miscellaneous conditions, please key in necessary data to smoothly continue the tracking process. When the program asks the user to input the head and tail location of a bee in a particular frame, first use the Zoom-in tool on the top right corner of the matlab figure window to zoom in to see the bee’s image closely, then click on the head and tail location of that bee. If you are running the tracking program on a laptop, sometimes there might be a delay of 1 to 2 seconds while clicking on the bee image. So, please wait before clicking on the bee image again.

Once the tracking program has finished tracking a bee, 4 different output files (.mat) will be saved in the downloaded repository folder. The details of the 4 output files are:

a) The tracked 2D coordinates of a bee can be obtained from ‘J\_BC\_HT\_CORDS\_DATASET\_#\_CAM\_##\_REFER\_BEE\_###.mat’.

b) The frames numbers at which occlusions happened will be stored in
 ‘OCCUL\_FRAMES\_DATASET\_#\_CAM\_##\_REFER\_BEE\_###.mat’.

c) The frames where a bee image would appear circular in shape will be saved in ‘RATIO\_A\_BY\_B\_DATASET\_#\_CAM\_##\_REFER\_BEE\_###.mat’. 

d) The list of frame numbers in which a bee was successfully tracked will be recorded in ALL\_FRAMES\_DATASET\_#\_CAM\_##\_REFER\_BEE\_###.mat.

\# - dataset number, \## - camera number (left or right camera), \### - bee number

#### Note: Using the developed tracking program, I was able to track bees flying in a bee cloud successfully. I have also tested this code by tracking individual bees flying in a curved tunnel. Please refer to the tracking results section to see the tracking output.

## How to track a bee from any frame number

To track bees from any point in the video, just key in the respective frame number from which you would like to begin the tracking process. 

Say for instance you are tracking a bee which is present in the arena from frame 1 to frame 800. Assume for some reason that matlab is not responding after 500 frames. In that case, you might need to close and reopen matlab again. Now, instead of starting from frame 1 again, you could continue the tracking process from frame 499. This could be done by simply updating the ‘frame_number’ and ‘frame_count’ variable with the frame number from which you would like to continue the tracking process. Also, uncomment the ‘load’ command for the 4 output variables so that the tracked information from frame 1 to frame 499 will be loaded into the program when executing the program from frame 500 onwards.


## Tracking results

### Bee cloud tracking results

Tracked head and tail positions of bees appearing in the left camera view

![](GIPHY_CAM_1_1.gif)

Tracked head and tail positions of bees appearing in the right camera view 

![](GIPHY_CAM_2_2.gif)

Tracked head and tail positions of a bee flying in a curved tunnel

![](GIPHY_CAM_5_5.gif)

## How to run the stereo-matching program 

Download the code for the stereomatching process from [here](https://github.com/mahadeesh/Bee-tracking-and-stereomatching-files). The stereomatching program along with the other related files are in ‘CODE_FOR_STEREOMATCHING’ folder. 

The sample inputs for executing the stereomatching program can be downloaded from ‘SAMPLE_INPUTS_FOR_STEREOMATCHING_PROGRAM’ folder given 
[here](https://figshare.com/s/cb7edd04e5818607ea9a ).

Please ensure that the ‘SAMPLE_INPUTS_FOR_STEREOMATCHING_PROGRAM’ folder is placed in the ‘CODE_FOR_STEREOMATCHING’ folder, to ensure that the Matlab programs have access to the sample input data via the specified directory paths.

From ‘CODE_FOR_STEREOMATCHING’ folder, open the file *Stereo_matching_using_MPD_RPE.m* in Matlab and click on the run button to execute the program. The inputs to this program are the 2D head and tail locations (in pixel coordinates) of bees in each camera view. The output will be a set of minimum perpendicular distances and reprojection error values. The correct matching bee for each reference bee will be the one with the lowest mean minimum perpendicular distance and mean reprojection error value. 

To use this stereomatching program, you must load the stereo camera’s internal and external parameters. Before executing this program, the user will have to generate two (stereo camera) calibration files which contain the internal and external parameters of the cameras. One stereo camera calibration file will be used to compute the minimum perpendicular distances, while the other file will be used to calculate the reprojection errors. The two calibration files are obtained by performing a calibration process, as described below.

The first calibration file is generated by using Yves Bouguet’s stereo calibration toolbox<sup>3</sup>. The second stereo calibration file is generated using Matlab’s stereo calibration application<sup>4</sup>. You can find the two calibration files by name *BC_stereo_cam_calib_trail_3_.mat* and *FEB_16_2017_STEREO_PARAMS.mat* in the SCI_DATA_CALIB_FILES folder, which were generated using our bee cloud videos. As these two files are already loaded in the stereomatching program, the user does not need to generate these calibration files while executing the sample data. However, for your data, you must generate your own stereo calibration files using Yves Bouguet and Matlab stereo calibration toolbox. The procedures for generating the stereo calibration files are explained in the respective websites<sup>3,4.
 
Please note that *BC_16_FEB_LEFT_REF_POINTS_CAM1_BKP.txt* and *BC_16_FEB_RIGHT_REF_POINTS_CAM2_BKP.txt* are the two files which were used to store calibration reference points while filming our data. These reference points are not included in the MPD and RPE calculations. When using this code for your data, you could simply retain these and run the code as this will not affect the output in any way.

Most importantly, the location or path of the sample input files in the local folder must be correctly supplied in ‘myfilename3’ and ‘*myfilename4*’ variables. In ‘*myfilename3*’ variable, the path of the folder which contains bee trajectories from the left camera view must be provided. Similarly, the location of the trajectories of bees flying in the right camera view must be provided in the ‘myfilename4’ variable. ‘DS1_REFER_BEES_2D_TRACK_RIGHT_CAM_VIEW’ is the folder where one can find the trajectories of bees flying in the right camera view. Similarly, ‘DS1_OTHER_BEES_2D_TRACK_LEFT_CAM_VIEW’ provides the trajectories of bees flying in the left camera view. The bee IDs of bees in these folders must be entered in ‘*Ref_bee_no*’ and ‘*opt_bee_numbers*’. The bee IDs for the sample data can be obtained from the name of the 2D trajectory file. 
 
For the sample input data, trajectories of 6 different bees flying in the right camera view (reference view) along with 100 trajectories of bees flying in the left camera view are provided.
 
Basically, the stereomatching program will pick each individual bee trajectory from the right camera view and will look for a corresponding match from 100 different bees in the left camera view. The correct matching bee is obtained by computing the minimum perpendicular distances and the reprojection errors. For each reference bee, a set of MPD and RPE values will be computed by the program which will be stored in '*MIN_DIST_REPROJ_ERR_INFO.mat*'. The correct corresponding pair is deduced by picking the pair with the lowest mean MPD and the lowest mean RPE. These values are available in '*MATCHING_BEE_MIN_DIST_REPROJ_ERR.mat*’.

 In the variable '*MATCHING_BEE_MIN_DIST_REPROJ_ERR.mat*', the column index pertains to the ID of the reference bee. The sample input data uses 6 reference bees (bee numbers 10,11,12, 17,18 and 19; hence, these are the only columns with non-zero entries.The values in the first row and third row indicate the obtained lowest mean minimum perpendicular distance and lowest mean reprojection error value, respectively. The numbers in the second and fourth row provide the IDs of the corresponding bee numbers for which lowest mean MPD and lowest mean RPE were obtained. 
 
If there are questions in relation to the operation of this software, contact Mandiyam Mahadeeswara at m.mandiyam@uq.edu.au or mdevaraj@outlook.com.


## References

1. Mahadeeswara, M. Y. & Srinivasan, M. V. Coordinated Turning Behaviour of Loitering Honeybees. Sci Rep 8, 16942, doi:10.1038/s41598-018-35307-5 (2018).
  
2. The Mathworks, Inc., Natick, Massachusetts, <https://www.mathworks.com/> (2021).

3. Bouguet, J.-Y. Complete Camera Calibration Toolbox for Matlab. Jean-Yves Bouguet’s Homepage <http://www.vision.caltech.edu/bouguetj/calib_doc/> (1999).

4. Matlab’s Stereo Camera Calibrator App <https://au.mathworks.com/help/vision/ug/stereo-camera-calibrator-app.html>.




