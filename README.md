# Bee tracking and stereomatching programs for 3D resconstruction of trajectories of freely flying honeybees

## Introduction

Here are the files which I developed to generate time resolved 3D position information of bees flying in a bee cloud. The bee tracking program is a semi-automatic tracker of freely flying honeybees, which was developed using simple techniques such background subtraction, image differencing, interpolation, and extrapolation, to name a few. Furthermore, the problem of correspondence matching is solved in a simple way by using two different techniques namely minimum perpendicular distance method and minimum reprojection error method. 

## Highlights of the developed bee tracking and stereomatching program

1. Despite being a semi-automatic bee tracking program, the final generated 3D position data contains head as well as tail locations of freely flying bees. Therefore, making it a unique and valuable dataset.

2. Solves the stereo-matching problem using two cameras. 

## Prerequisites

All you need is MATLAB from MathWorks<sup>1</sup>. I have tested the code using MATLAB R2020 update 5 running on a macOS Catalina. All function files as given in the branch must be present.

## Contents of the files

1)	*Sci_data_bee_tracking_pgm_3.m* - The main bee tracking program which tracks individual bees flying in the bee cloud arena.

2)	*Stereo_matching_using_MPD_RPE.m* – a stereomatching program which matches each bee in one camera view with its corresponding bee image in the other camera view. 

3)	*stereo_triangulation_camera_positions.m* – is used to calculate the positions of the two cameras with respect to the origin.

4)	*fit_ellipse.m* – is an open-source ellipse fitting program.

## How to run the bee tracking program

1. Download the repository into a folder

2. To run the bee tracking program, simply open the file *Sci_data_bee_tracking_pgm_3.m* in Matlab and click on the run button. Before clicking on the run command, make sure all necessary input videos/files are properly placed in the local repository folder. 

You can download all input video/files necessary to run the tracking program from [here](https://figshare.com/articles/media/Multi-Object_Tracking_in_Heterogeneous_environments_MOTHe_for_animal_video_recordings/11980356/3). Place all the downloaded videos/files in the folder where you have downloaded the repository.  

As this is a semi-automatic tracking program, the program prompts for the user inputs under certain miscellaneous conditions, please provide the necessary data to smoothly continue the tracking process. When the program asks the user to input the head and tail location a bee in a particular frame, first use the Zoom-in tool on the right top corner of the matlab figure window to zoom in to see the bee’s image closely, then click on the head and tail location of that bee. If you are running the tracking program on a laptop, sometimes there might be a delay of 1 to 2 seconds while clicking on the bee image during the zooming process. So, please wait before clicking on the bee image again. 

Once the tracking program has finished tracking a bee, 4 different output files (.mat) will be saved in the downloaded repository folder. The details of the 4 output files are:

The tracked 2D coordinates of a bee can be obtained from ‘J\_BC\_HT\_CORDS\_DATASET\_#\_CAM\_##\_REFER\_BEE\_###.mat’.

The frames numbers at which occlusion events occurred are available in ‘OCCUL\_FRAMES\_DATASET\_#\_CAM\_##\_REFER\_BEE\_###.mat’.

The frames where a bee image would appear circular in shape will be saved in ‘RATIO\_A\_BY\_B\_DATASET\_#\_CAM\_##\_REFER\_BEE\_###.mat’. 

The frame numbers for which the program executed without any errors are recorded in this variable ALL\_FRAMES\_DATASET\_#\_CAM\_##\_REFER\_BEE\_###.mat.

The list of frame numbers in which a bee was successfully tracked will be recorded in ALL\_FRAMES\_DATASET\_#\_CAM\_##\_REFER\_BEE\_###.mat.

\# - dataset number, \## - camera number (left or right camera), \### - bee number

## How to run the stereo-matching program 

Similarly, open the file *Stereo_matching_using_MPD_RPE.m* in Matlab and click run to execute the program. The input to this program is 2D head and tail locations (in pixel coordinates) of bees in each camera view. The output will be a set of minimum perpendicular distance and reprojection error values. The correct corresponding matching bee for each reference bee will be one with lowest minimum perpendicular distance and reprojection error value.

To use this stereo-matching program, you must load the stereo camera’s internal and external parameters. Before executing this program, the user will have to generate two (stereo camera) calibration files which consists of camera internal and external parameters. One stereo camera calibration file will be used in computing the minimum perpendicular distances, while the other one will be used in calculating the reprojection errors. One could obtain these two calibration files by performing a calibration process.

So, I generated the first calibration file using Yves Bouget’s stereo calibration toolbox<sup>2</sup>. The second stereo calibration file was generated using Matlab’s stereo calibration toolbox<sup>3</sup>. You can find the two calibration files by name *BC_stereo_cam_calib_trail_3_.mat* and *FEB_16_2017_STEREO_PARAMS.mat* in the repository, which were generated using our bee cloud videos. These two files are initially loaded in the stereomatching program. 

For your data, you must generate your own stereo calibration files using Yves Bouget and Matlab stereo calibration toolbox. The procedure to generate the stereo calibration files are explained in their websites<sup>2,3.

The variable in which all minimum perpendicular distances and reprojection error values are stored in *MIN_DIST_REPROJ_ERR_INFO.mat*.

## References

1. The Mathworks, Inc., Natick, Massachusetts, <https://www.mathworks.com/> (2021).

2. Bouguet, J.-Y. Complete Camera Calibration Toolbox for Matlab. Jean-Yves Bouguet’s Homepage <http://www.vision.caltech.edu/bouguetj/calib_doc/> (1999).

3. Matlab’s Stereo Camera Calibrator App <https://au.mathworks.com/help/vision/ug/stereo-camera-calibrator-app.html>.




