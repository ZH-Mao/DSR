# Direct-Simultaneous-Registration
This repository is an implementation of paper "DSR: Direct Simultaneous Registration for Multiple 3D Images".

1. The folder 'dataset' includes eight volumes for testing the code.
2. 'DSR_Main.m' is the main function of the algorithm.
3. 'Pairwise_Main.m' is used to get the relative poses between two images. The results are used as initial guess for 'DSR_Main.m'.
4. 'Registration.m' is used to align two images for evaluate the registration acuray visually.
5. 'SliceBrowser.m' is used to view the volumes via three orthogonal slices.
