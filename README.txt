

================================================================================
________________________________________________________________________________
SOME DOCUMENTATION WITH SCHNITZCELLS (Tans lab version)
Initial version: 2014/04/16 Martijn Wehrens
________________________________________________________________________________

================================================================================

Contents:

Part 1a.    How to get started using main scripts.
Part 1b.    (Legacy) Documentation on how to use Excel file to analyze a 
            position.

Part 2.     Documentation about specific properties of the code.

Part X.     Some conventions used for coding.

================================================================================
________________________________________________________________________________
PART 1a
HOW TO ANALYSE A POSITION USING MATLAB MAIN SCRIPTS
________________________________________________________________________________

================================================================================

================================================================================
General information
================================================================================

(Note 2017: See the file 2017-04-21_schnitzcells_short_guide.docx for a quick guide on using the GUI.)

In end 2015, beginning of 2016, I made an effort to update the code in such a way that you can run your whole analysis from one script. Moreover, this script can also be controlled from a GUI (graphical user interface).

The main file you currently need to perform an analysis is 
>> Schnitzcells_masterscript
I've added quite some comments into this file. The final goal is to make this file self-explanatory. Also, the GUI should (in the future) contain enough comments & instructions such that it can be ran without needing external help files.
The main files of the GUI are
>> MW_GUI_schnitzcells 
Which can be edited easily with the command:
>> guide MW_GUI_schnitzcells 
The files linked to this script are MW_GUI_schnitzcells.m and MW_GUI_schnitzcells.fig.

In principle, a load of what is explained down here has now become legacy (Excel is not used to execute commands any more, only as a parameter file). However, the functions that are used are still the same, and a load of text below explains how these functions work. In the future, most of this text should be found in abovementioned script files. For now it remains here for reference.

Convenient scripts for further analysis are savingFluorDynamicsData.m and plottingGeneralDynamicData.m. These scripts can be found in the https://Leeuwenhoek@bitbucket.org/microscopeguerrillas/plotting_fluo_dyn_data.git repository, spefically the branch martijn_develop.

-MW 2016/01

================================================================================
________________________________________________________________________________
PART 1b
HOW TO USE EXCEL FILE AND MATLAB TO ANALYZE A POSITION
________________________________________________________________________________

================================================================================

================================================================================
General info on the Excel sheet
================================================================================

Generally, there are different matlab functions that perform the analysis. These functions add or edit files and folders in the data directory that contain analyzed data. A few files are plain text files, but the files containing the main body of data calculated are .mat files. The matlab functions also output images (to check the analysis intermediately) and plots (which are often the "end output" of the analysis). 

The Matlab functions take certain parameters as input. One of the most important ones is the "p" struct. This holds the "overhead"/"administrative" information, like what is the current analysis directory, date of the experiment, which microscope was used, which fluor etc. Most main functions take p as input, they need it to know which directory to edit and how.

Aside from parameters that are sometimes set (either constants at the top or sometimes hard-coded numbers!) in the function, there are also parameters that influence the segmentation and other parts of the analysis that can be given as input paramters. Important ones can be changed at the top of the Excel file, and are correspondingly updated in the rest of the sheet such that function are called with these parameter settings.

The Excel sheet contains essentially a list of matlab functions that need to be executed consecutively to perform the analysis. Since every dataset is analyzed with slightly different parameters (e.g. different range of analyzed frames, different fluor colors, ..) at the end of the analysis, the Excel sheet contains important experimental parameters that need to be stored. Therefore, a copy of the Excel sheet used for a specific dataset is typically stored in the directory of the dataset. (In fact, most dataset contain multiple "positions" - i.e. movies of growing colonies, respectively - and thus also multiple Excel sheets.)

Furthermore, because the analysis takes quite some effort to perform, a "preliminary" analysis is mostly performed. In such an analysis, not all frames are considered, and this analysis will only tell the user whether cells show normal growth behavior, and whether the dataset contains analyzable data. If that's the case, one can proceed with the full analysis. The matlab functions required for such a "preliminary" analysis are listed in a separate Excel file which contains the keyword "preliminary" (as opposed to the keyword "full_analysis").

A short outline of the general analysis:
- Files can be cropped to throw away area in the field of view that's not of interest.
- Movie is segmented.
- Segmentation is checked manually to correct for mistakes that the algorithm makes.
- Tracking of the cell lineages is performed. I.e. 'individuals' are identified in each movie frames. 
- Tracking is corrected manually. This might involve iteration of segmentation/tracking steps.
- Where applicable, fluorescent images are corrected.
- 'Schnitzcells' file is constructed. This contains a list of 'individual' bacteria (also known as 'schnitzes'), and the currently known paramters for that individual (e.g. cell size for each frame the bacteria is seen in). 
- Additional parameters are calculated based upon this "schnitzcells" struct. E.g. growth speed at each frame, enzyme concentration, enzyme production rate, etc.
- An additional normalization can be performed - which is referred to as "cyccor" - that normalizes for cell cycle effect. Here, parameters are normalized per schnitz, by dividing them by the individual-average behavior of that paramter. 
- One of the most important analysis which are typically performed involve making correlation functions for two measured parameters. This can however not be done straightforwardly on lineage data. Thus "branches" - each giving a full parameter trace of one lineage - are generated, which are weighed in order to calculate the correlation functions. (And can also be used for other analyses.)

A synopsis of the output generated:
- segmentation files (matrix, stored as .mat file)
- tracking files (text files, show ancestry relation frame to frame)
- schnitzcells .mat file (hold all params per individual over time)
- plots
(branches generated are not saved currently)

Where this output can be found:
Everything is stored per "position". Per default each position is stored in posX. However, the cropping function creates a folder "posXcrop", where data generated is stored in subdirectories:
*
./images 
    > directory contains the raw image files
./segmentation contains segmentation in a matrix format, stored in .mat file
    > depending on settings, subdir contains images that show intermediate steps of segmentation
    > also different parameter settings which lead to different segmentation can be found in subdirs. 
./data contains 
    > The tracking of cells from frame to frame
    > the file in which the central "schnitzcells" struct is stored (posXcrop-Schnitz.mat)
./analysis
    > contains output plots, sorted per topic.
    > also movies, part of the manual tracking correction, are found here


================================================================================
More specific info on each of the steps.
================================================================================

The Schnitzcells analysis constists of different steps. 

(Step 0: Cropping the images)

================================================================================
Step 1: Segmentation of the cells
================================================================================
The segmentation needs to be checked manually; There is a matlab function available to correct the segmentation (.._manualcheckseg()) and a function which tries to detect mistakes. It is convenient to combine the use of these two functions. The advised procedure is as follows:
    o First, run .._manualcheckseg() once to correct obvious errors. Cells that look incorrect to the algorithm are colored white by this function.
    o Second, run DJK_analyzeSeg(), this will identify cases which the algorithm thinks are incorrect. These can be corrected by running .._manualcheckseg() again.
    DJK_analyzeSeg() produces the following output which gives information about incorrect segmentation:
        >> analysis\segmentation\fitRange1_199
    This is also printed to screen, and gives a list of possible problems for tracking.
*
ADDITIONAL OPTIONS/PROBLEM SOLVING:
Before segmentation, the algorithm selects a piece it deems relevant to segment. This is done by max. contrast [I think; check]. Sometimes it crops out a piece of dirt instead of actual bacteria. In that case one can use either can set p.useFullImage = 1 with PN_segmoviephase_3colors.

*
NOTES: There are some "hidden" settings to the algorithm that determines whether cells might be incorrectly segmented. 
Settings in PN_imshowlabel:
- fractionbelowwhite = 0.2; 
    cells which are a fraction of <fractionbelowwhite> smaller than the median are marked white.



================================================================================
Step 2: Tracking, correcting again
================================================================================
- Now tracking can be performed. (Use .._tracker_centroid_vs_area().)
- Also here, problems can be detected by the algorithm. Execute ..._analyzeTracking() to run the detection algorithm and ..makeMovie() three times (with different settings) to produce images which highlight the problems. This produces the following output.
        >> ../analysis/tracking/manualRange...
        File which lists all problems.
        >> ../analysis/movies/
        Movies which highlight suspected errors.
        By default, three different movies get made. One problem movie enumerates schnitzcells w. issues. The second enumerates cells per frame movie (cell number [NOT schnitzcells]). The third 
% schnitzcells movie - numbers Schnitzcells (only for later)

Often problems will be caused by incorrect segmentation, and need to be resolved by running .._manualcheckseg() again. In some other cases, problems will need to be resolved by manually editing the files with the lineage links:
        << /data/
	The posXcrop-djk-output-YYY-to-ZZZ.txt files link schnitzes from frame YYY to frame ZZZ.  
        These files are structured as follows:
        pn p1 p2 ch
        Where pn is a parent of a nondividing cell, and ch its child in the next frame. If a cell divides, two lines are used. One links one child to the parent, and the other the other child to the parent. Somewhat unintuitively, the parent of child one is listed at p1, and the parent of child 2 is listed at p2 (but for one dividing cell, p1=p2). The child is always listed in the same place (but the ch is different for one dividing cell). 
        Note that there are two kinds of numbers: numbers for the Schnitzes (the "individuals") and a numbering per frame which simple numbers the areas in this frame identified as cells. Numbers in the /data/ folder files are necessarily these area numbers. Based upon the lineage structure defined by the posXcrop-djk-output-YYY-to-ZZZ.txt files, schnitzes are identified.

About schnitzcells structure
---
A minimal schnitz structure is described by P (the schnitz number of it's parent), children schnitz numbers D and E, sister schnitz number S, frame_nrs (an array) listing each frame this schnitz occurs in, cellno (an array) listing each cell number of the tracked cell within each of the frames, and N, the number of frames that distinct cell appears.        

================================================================================
Calculating parameters
================================================================================
The growth rate is calculated by the function DJK_addToSchnitzes_mu. Note that this function behaves somewhat inconvenient when multiple fluors are used. In general, this function creates the parameter muPXX_YY and muPXX_YY_all, where XX is the width of the fit window (# frames), and YY is the length parameter on which the calculation was based. muPXX_YY _only contains_ calculated mu values that match with the points at which fluor images where taken. If there were multiple fluor colors used, the muPXX_YY parameter added to the schnitzcells corresponds to the fluor listed first ("fluor1").


================================================================================
Fluor corrections
================================================================================

<Note that the word "flatfield" is sometimes used were actually "background" is the correct term. "Flatfield" is, in literature, used as synonym for "shading", but here only the term shading is used to designate shading.>

1. The function DJK_correctFluorImage_anycolor corrects for artefacts that arise when taking fluorescence images with the microscope. It substracts the camera noise, by substracting an image taken w/o lighting (we refer to this as flatfield, but this termonology is a bit undefined). 
2. Uneven illumination is compensated for by dividing by an image of only a fluorescent dye (normalized). 
3. The image is deconvolved to 'computate out' scattering effects.

Most important lines of this script, concerning steps 1-2, are:
fluor2image = fluor2image-flatfield_crop;
fluor2image = shading_mean.*fluor2image./shading_crop;


Different fluor parameters in schnitzcells struct
---
See also notes in DJK_addToSchnitzes_fluor_anycolor.
Excerpt:
% NOTES ON PROCEDURE AND OUTPUT
% ===
% In principle, pixels are selected from the fluor image using the 
% segmentation file, which contains the regions of the detected cells
% encoded as indices in a matrix. 
% There are however different corrections performed on the fluor images,
% and one can determine either the sum or the average of the fluor
% intensity within a cell. Furthermore, an additional method is introduced
% here to take only the fluor values from a subset of the detected cell
% area, namely only its central area, to prevent artifacts introduced by
% (incorrect) border detection.
% All these options lead to fluor values being calculated differently, and
% these different values are all stored.
% The parameters <fluor>1_suffix to <fluor>5_suffix (where <fluor> is a 
% capital letter encoding for the fluor, e.g. G, C, Y, ..) hold the 
% different consecutive corrections that were performed on the fluor image. 
% fluor6 holds the fluor value determined from the central area. The
% suffixes of these parameters tell you by which method the different pixel
% intensities are summarized (e.g. sum, mean, ..).

[TODO MW: see also notes 27-11-2014 for more info]

================================================================================

Creating branches
================================================================================
In order to be able to calculate cross correlations so-called branchgroups are created first. Key steps of that procedure are outlined below.
See also Daan Kiviet's thesis: "The lac Operon: Fluctuations, Growth and Evolution", p. 106.

To compute e.g. correlation functions, one should have a continuous function of two parameters for two parameters. E.g. growth rate mu and concentration of fluorescent reporter. This kind of data is however not available for lineages from the schnitzcells struct. Therefore, "branches" are created. These contains parameters X(t) per lineage. This necessitates that data is redundant. E.g. the trace X(t) for a 1st generation cell that has 100 offspring cells in the final colony, shuold be copied a 100 times because it exists in all lineages. 
|      ------           |      ------ branch1
| -----          =>     | =====  
|      ------           |      ------ branch2        
|______________         |______________         

There are however some complicating -but important - caveats/operations to this process, which will be outlined below.

In general, the function DJK_getBranches 'copies' X_{n<N}(t) data to X_{n<N,i}(t) i={1..F} such that X(t) data is available for each lineage. (With points early in time being shared among lineages.)

In principal, correlation functions can be determined from these "branches". However, there are additional important processes. 
1. Weighing of X(t) data to construct the branches
2. Substracting colony average at t
3. Creating branch groups
Typically in that order, although (3) and (2) can also be done the other way around.

Default order:
(1) Weighing.
In DJK_getBranches data is weighed such that the redundancy in the data is compensated for. This is done by the "3/4" weighing procedure (See eq. 6.6 Daan Kiviet's thesis). Different weighing schemes are however possible. 
(2) Substracting colony average 
If - and this behavior is indeed qualitatively observed - there are whole colony fluctuations (in growth), this means that correlations not corrected for this are not single cell correlations, but correlations due to fluctuations in the whole colony. To filter out that effect, for each point X(t) the time average is substracted. (Note that this may be a non-intuitive procedure.) I.e. for a branch j, X_j(t)-<X(t)>_j. The function performing this correction is called DJK_addToBranches_noise. The prefix "noise_" is added, because technically this parameter will represent the noise of X_j(t) around <X(t)>_j. 
This correction is quite important for how the final correlation will look. (One could plot the correlation function of the colony mean of the observables of interest to investigate the behavior that is filtered out.) 
(3) Creating branch groups
To get an estimate of the error in the correlation function, the colony branches are subdivided in groups, and error bars are calculated based on the correlations found for each subgroup. (This subgrouped branch data is referred to as branchgroups.) This is done by the function  DJK_divide_branch_data.
[Comment MW: This procedure seems rather ad hoc, but seems to relate to bootstrapping procedures often performed in statistics. TODO: find out how well it relates, and whether an improvement might be to perform something more resembling bootstrapping.] 
*
Functions involved in this process are DJK_getBranches, DJK_addToBranches_noise, DJK_trim_branch_data, DJK_divide_branch_data (alternatively: NW_divide_branch_data_inittime), DJK_plot_crosscorrelation_standard_error_store.


Alternative order:
In principle, the different subgroup branches should follow similar colony average trends. However (especially working with small datasets) this might not always be the case. To prevent artefacts from being introduced here again, steps (2) and (3) of the previous order can be inverted. 
The following functions should be used to follow that procedure:
DJK_getBranches (standard), DJK_trim_branch_data (standard), DJK_divide_branch_data (standard), NW_addToBranchGroups_noise (non standard),  DJK_plot_crosscorrelation_standard_error_store (standard).


================================================================================
________________________________________________________________________________
PART 2
SOME NOTES ON HOW MATLAB CODE WORKS
________________________________________________________________________________

================================================================================


See the "Schnitzcells_structure.xlsx" file for an overview of the matlab functions and a description.

SOME NOTES ON SEGMENTATION
===

THE L STRUCTURE

The "L" structure contains the segmented bacteria. This is simply a N * M matrix with zeros with each area that comprises a bacteria marked by a unique number that is not zero. N * M equals the pixel size of the picture.

L matrices are stored in \posXcrop\segmentation\posYcropsegZZZ.mat. These files contain a few matrices. 
- phsub     This contains the orginal phase ("ph") contrast image that was taken.
- LNsub     This contains the uncorrected segmented image.
- Lc        This contains the corrected segmented image.
Note that if Lc exists in the .m file, this means that the image has been corrected. 
See MW_illustrative_plotting_L_matrices.m for an example where these matrices are plotted. 

____________________________________________________________
PART X
SOME CONVENTIONS USED FOR CODING
____________________________________________________________

GENERAL STYLE
===
Please follow this styleguide:
http://www.mathworks.com/matlabcentral/fileexchange/2529-matlab-programming-style-guidelines
*
For constants that are set in functions, please use cap locks, and place them at the start of a function, so that users can easily find them. Preferably, also mention constants in the function description (i.e. first chunk of comments).


KEYWORDS AND STYLE IN THE CODE
===
TODO or "blubb" (depricated) indicates codes is "shaky", i.e. there is no bug observed, but there might be a better/neater way to program this. It might also indicate that there is a problem which needs to be solved.
*
There is no convention for sectioning yet.


CHANGING THE CODE
===
When you change a line in the code, please add your initials and the current date behind that line in comments, so people importing the changes know when something changed. (For clearity, the reason of the change can sometimes also be added.)




















