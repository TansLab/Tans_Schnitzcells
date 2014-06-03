

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
Before segmentation, the algorithm selects a piece it deems relevant to segment. This is done by max. contrast [I think; check]. Sometimes it crops out a piece of dirt instead of actual bacteria. In that case one can use either of following functions to instead analyze the whole image:
%NW_nocut_segmoviephase_3colors_diffStrains(...);
%NW_segmoviephase_3colors_diffStrains(...);
(Instead of PN_segmoviephase_3colors(...);.)

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

        

____________________________________________________________
PART 2
SOME NOTES ON HOW MATLAB CODE WORKS
____________________________________________________________


See the "Schnitzcells_structure.xlsx" file for an overview of the matlab functions and a description.

SOME NOTES ON SEGMENTATION
===

THE L STRUCTURE

The "L" structure contains the segmented bacteria. This is simply a N * M matrix with zeros with each area that comprises a bacteria marked by a unique number that is not zero. N * M equals the pixel size of the picture.



____________________________________________________________
PART X
SOME CONVENTIONS USED FOR CODING
____________________________________________________________


KEYWORDS AND STYLE IN THE CODE
===
TODO or "blubb" (depricated) indicates codes is "shaky", i.e. there is no bug
observed, but there might be a better/neater way to program this. It might also
indicate that there is a problem which needs to be solved.
*
There is no convention for sectioning yet.


CHANGING THE CODE
===
When you change a line in the code, please add your initials and the current 
date behind that line in comments, so people importing the changes know 
when something changed. (For clearity, the reason of the change can sometimes
also be added.)




















