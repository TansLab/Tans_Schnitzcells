clear all
close all

%pos2
p = DJK_initschnitz('pos2crop','2011-10-02','e.coli.AMOLF','rootDir','X:\colonies\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040]);
%4 - analyze segmentation and correct
DJK_analyzeSeg(p,'manualRange',[7:54:601]);
%5 - perform tracking, singleCell
DJK_trackcomplete(p,'trackRange',[7:54:601],'trackMethod','singleCell');
%6 - correcting fluor 
% fisrt color
optimalShift = DJK_getFluorShift(p,'manualRange', [7:54:601]);
load 'X:\schnitzcells\Daan_additions\shading\Correction_10Mhz_110801_200ms.mat' flatfield shading replace
load 'X:\schnitzcells\Daan_additions\psf\PSF_090402_centered.mat' PSF
DJK_correctFluorImage(p, flatfield, shading, replace,'manualRange', [7:54:601],  'fluorShift', optimalShift, 'deconv_func', @(im) deconvlucy(im, PSF));
DJK_analyzeFluorBackground(p,'manualRange', [7:54:601]);
% second color
optimalShift2 = DJK_getFluorShift_red(p,'manualRange', [7:54:601]);
load 'X:\schnitzcells\Daan_additions\shading\Correction_10Mhz_110801_200ms.mat' flatfield shading replace
load 'X:\schnitzcells\Daan_additions\psf\PSF_090402_centered.mat' PSF
DJK_correctFluorImage_red(p, flatfield, shading, replace,'manualRange', [7:54:601],  'fluor2Shift', optimalShift2, 'deconv_func', @(im) deconvlucy(im, PSF));
DJK_analyzeFluorBackground_red(p,'manualRange', [7:54:601]);
%7 - make & add to schnitzcells
DJK_compileSchnitzImproved_2colors(p,'quickMode',0);
DJK_addToSchnitzes_length(p);
DJK_addToSchnitzes_length(p, 'onScreen', 0);
DJK_addToSchnitzes_fluor(p, 'onScreen', 0);
DJK_addToSchnitzes_fluor_red(p, 'onScreen', 0);
[p,schnitzcells] = DJK_compileSchnitzImproved_2colors(p,'quickMode',1);
% Overview all data
fitTime = DJK_analyzeMu(p, schnitzcells, 'xlim', [0 1500], 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'av_mu_fitNew', 'xlim', [0 1500], 'ylim', [0 1.5], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'Y3_mean_all', 'xlim', [0 1500], 'ylim', [0 50], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'Y6_mean_all', 'xlim', [0 1500], 'ylim', [0 70], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'R3_mean_all', 'xlim', [0 1500], 'ylim', [0 50], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'R6_mean_all', 'xlim', [0 1500], 'ylim', [0 100], 'fitTime', fitTime, 'onScreen', 0);



%pos3
p = DJK_initschnitz('pos3crop','2011-10-02','e.coli.AMOLF','rootDir','X:\colonies\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040]);
%4 - analyze segmentation and correct
DJK_analyzeSeg(p,'manualRange',[15:54:555]);
%5 - perform tracking, singleCell
DJK_trackcomplete(p,'trackRange',[15:54:555],'trackMethod','singleCell');
%6 - correcting fluor 
% fisrt color
optimalShift = DJK_getFluorShift(p,'manualRange', [15:54:555]);
load 'X:\schnitzcells\Daan_additions\shading\Correction_10Mhz_110801_200ms.mat' flatfield shading replace
load 'X:\schnitzcells\Daan_additions\psf\PSF_090402_centered.mat' PSF
DJK_correctFluorImage(p, flatfield, shading, replace,'manualRange', [15:54:555],  'fluorShift', optimalShift, 'deconv_func', @(im) deconvlucy(im, PSF));
DJK_analyzeFluorBackground(p,'manualRange', [15:54:555]);
% second color
optimalShift2 = DJK_getFluorShift_red(p,'manualRange', [15:54:555]);
load 'X:\schnitzcells\Daan_additions\shading\Correction_10Mhz_110801_200ms.mat' flatfield shading replace
load 'X:\schnitzcells\Daan_additions\psf\PSF_090402_centered.mat' PSF
DJK_correctFluorImage_red(p, flatfield, shading, replace,'manualRange', [15:54:555],  'fluor2Shift', optimalShift2, 'deconv_func', @(im) deconvlucy(im, PSF));
DJK_analyzeFluorBackground_red(p,'manualRange', [15:54:555]);
%7 - make & add to schnitzcells
DJK_compileSchnitzImproved_2colors(p,'quickMode',0);
DJK_addToSchnitzes_length(p);
DJK_addToSchnitzes_length(p, 'onScreen', 0);
DJK_addToSchnitzes_fluor(p, 'onScreen', 0);
DJK_addToSchnitzes_fluor_red(p, 'onScreen', 0);
[p,schnitzcells] = DJK_compileSchnitzImproved_2colors(p,'quickMode',1);
% Overview all data
fitTime = DJK_analyzeMu(p, schnitzcells, 'xlim', [0 1500], 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'av_mu_fitNew', 'xlim', [0 1500], 'ylim', [0 1.5], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'Y3_mean_all', 'xlim', [0 1500], 'ylim', [0 50], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'Y6_mean_all', 'xlim', [0 1500], 'ylim', [0 70], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'R3_mean_all', 'xlim', [0 1500], 'ylim', [0 50], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'R6_mean_all', 'xlim', [0 1500], 'ylim', [0 100], 'fitTime', fitTime, 'onScreen', 0);



%pos5
p = DJK_initschnitz('pos5crop','2011-10-02','e.coli.AMOLF','rootDir','X:\colonies\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040]);
%4 - analyze segmentation and correct
DJK_analyzeSeg(p,'manualRange',[13:54:553]);
%5 - perform tracking, singleCell
DJK_trackcomplete(p,'trackRange',[13:54:553],'trackMethod','singleCell');
%6 - correcting fluor 
% fisrt color
optimalShift = DJK_getFluorShift(p,'manualRange', [13:54:553]);
load 'X:\schnitzcells\Daan_additions\shading\Correction_10Mhz_110801_200ms.mat' flatfield shading replace
load 'X:\schnitzcells\Daan_additions\psf\PSF_090402_centered.mat' PSF
DJK_correctFluorImage(p, flatfield, shading, replace,'manualRange', [13:54:553],  'fluorShift', optimalShift, 'deconv_func', @(im) deconvlucy(im, PSF));
DJK_analyzeFluorBackground(p,'manualRange', [13:54:553]);
% second color
optimalShift2 = DJK_getFluorShift_red(p,'manualRange', [13:54:553]);
load 'X:\schnitzcells\Daan_additions\shading\Correction_10Mhz_110801_200ms.mat' flatfield shading replace
load 'X:\schnitzcells\Daan_additions\psf\PSF_090402_centered.mat' PSF
DJK_correctFluorImage_red(p, flatfield, shading, replace,'manualRange', [13:54:553],  'fluor2Shift', optimalShift2, 'deconv_func', @(im) deconvlucy(im, PSF));
DJK_analyzeFluorBackground_red(p,'manualRange', [13:54:553]);
%7 - make & add to schnitzcells
DJK_compileSchnitzImproved_2colors(p,'quickMode',0);
DJK_addToSchnitzes_length(p);
DJK_addToSchnitzes_length(p, 'onScreen', 0);
DJK_addToSchnitzes_fluor(p, 'onScreen', 0);
DJK_addToSchnitzes_fluor_red(p, 'onScreen', 0);
[p,schnitzcells] = DJK_compileSchnitzImproved_2colors(p,'quickMode',1);
% Overview all data
fitTime = DJK_analyzeMu(p, schnitzcells, 'xlim', [0 1500], 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'av_mu_fitNew', 'xlim', [0 1500], 'ylim', [0 1.5], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'Y3_mean_all', 'xlim', [0 1500], 'ylim', [0 50], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'Y6_mean_all', 'xlim', [0 1500], 'ylim', [0 70], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'R3_mean_all', 'xlim', [0 1500], 'ylim', [0 50], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'R6_mean_all', 'xlim', [0 1500], 'ylim', [0 100], 'fitTime', fitTime, 'onScreen', 0);


%pos 6
p = DJK_initschnitz('pos6crop','2011-10-02','e.coli.AMOLF','rootDir','X:\colonies\', 'cropLeftTop', [1,1], 'cropRightBottom', [1392,1040]);
%4 - analyze segmentation and correct
DJK_analyzeSeg(p,'manualRange',[12:54:606]);
%5 - perform tracking, singleCell
DJK_trackcomplete(p,'trackRange',[12:54:606],'trackMethod','singleCell');
%6 - correcting fluor 
% fisrt color
optimalShift = DJK_getFluorShift(p,'manualRange', [12:54:606]);
load 'X:\schnitzcells\Daan_additions\shading\Correction_10Mhz_110801_200ms.mat' flatfield shading replace
load 'X:\schnitzcells\Daan_additions\psf\PSF_090402_centered.mat' PSF
DJK_correctFluorImage(p, flatfield, shading, replace,'manualRange',[12:54:606],  'fluorShift', optimalShift, 'deconv_func', @(im) deconvlucy(im, PSF));
DJK_analyzeFluorBackground(p,'manualRange', [12:54:606]);
% second color
optimalShift2 = DJK_getFluorShift_red(p,'manualRange',[12:54:606]);
load 'X:\schnitzcells\Daan_additions\shading\Correction_10Mhz_110801_200ms.mat' flatfield shading replace
load 'X:\schnitzcells\Daan_additions\psf\PSF_090402_centered.mat' PSF
DJK_correctFluorImage_red(p, flatfield, shading, replace,'manualRange',[12:54:606],  'fluor2Shift', optimalShift2, 'deconv_func', @(im) deconvlucy(im, PSF));
DJK_analyzeFluorBackground_red(p,'manualRange',[12:54:606]);
%7 - make & add to schnitzcells
DJK_compileSchnitzImproved_2colors(p,'quickMode',0);
DJK_addToSchnitzes_length(p);
DJK_addToSchnitzes_length(p, 'onScreen', 0);
DJK_addToSchnitzes_fluor(p, 'onScreen', 0);
DJK_addToSchnitzes_fluor_red(p, 'onScreen', 0);
[p,schnitzcells] = DJK_compileSchnitzImproved_2colors(p,'quickMode',1);
% Overview all data
fitTime = DJK_analyzeMu(p, schnitzcells, 'xlim', [0 1500], 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'av_mu_fitNew', 'xlim', [0 1500], 'ylim', [0 1.5], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'Y3_mean_all', 'xlim', [0 1500], 'ylim', [0 50], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'Y6_mean_all', 'xlim', [0 1500], 'ylim', [0 70], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'R3_mean_all', 'xlim', [0 1500], 'ylim', [0 50], 'fitTime', fitTime, 'onScreen', 0);
DJK_plot_avColonyOverTime(p, schnitzcells, 'R6_mean_all', 'xlim', [0 1500], 'ylim', [0 100], 'fitTime', fitTime, 'onScreen', 0);

