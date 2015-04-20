

____________________________________________________________
SOME DOCUMENTATION WITH SCHNITZCELLS (Tans lab version)
2014/04/16 Martijn Wehrens
____________________________________________________________



Contents:

Part 1.     Documentation on how to use Excel file to analyze a position.

Part 2.     Documentation about specific properties of the code.

Part X.     Some conventions used for coding.


____________________________________________________________
PART 1
HOW TO USE EXCEL FILE AND MATLAB TO ANALYZE A POSITION
____________________________________________________________

The Schnitzcells analysis constists of different steps. 

(Step 0: Cropping the images)

Step 1: Segmentation of the cells
===
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



Step 2: Tracking, correcting again
===
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
        These files are structured as follows:
        pn p1 p2 ch
        Where pn is a parent of a nondividing cell, and ch its child in the next frame. If a cell divides, two lines are used. One links one child to the parent, and the other the other child to the parent. Somewhat unintuitively, the parent of child one is listed at p1, and the parent of child 2 is listed at p2 (but for one dividing cell, p1=p2). The child is always listed in the same place (but the ch is different for one dividing cell). 
        Note that there are two kinds of numbers: numbers for the Schnitzes (the "individuals") and a numbering per frame which simple numbers the areas in this frame identified as cells. Numbers in the /data/ folder files are necessarily these area numbers.

About schnitzcells structure
---
A minimal schnitz structure is described by P (the schnitz number of it's parent), children schnitz numbers D and E, sister schnitz number S, frame_nrs (an array) listing each frame this schnitz occurs in, cellno (an array) listing each cell number of the tracked cell within each of the frames, and N, the number of frames that distinct cell appears.        



Fluor corrections
---

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

Creating branches
===
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





____________________________________________________________
PART 2
SOME NOTES ON HOW MATLAB CODE WORKS
____________________________________________________________


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




















