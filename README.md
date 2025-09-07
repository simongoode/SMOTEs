# A Machine Learning Pipeline for Discovering Sub-Minute Optical Transients (SMOTEs)
This repository contains the code for the machine learning-accelerated search pipeline developed for the **Deeper, Wider, Faster (DWF) programme**. The primary goal of this project is to identify high-quality, single-detection transient candidates with durations on sub-minute timescales, a previously unexplored regime in optical astronomy.
This work is based on the paper: "_A Machine Learning empowered search for Sub-Minute Optical Transient Events with the Deeper, Wider, Faster programme_".

## Project Overview
Optical transient surveys are generating vast datasets, making manual inspection for rare events infeasible. This project introduces a novel pipeline that leverages machine learning to efficiently search for sub-minute optical transients. By analyzing data from the DWF programme, we have developed a method to filter through hundreds of thousands of light curves to find promising candidates.
The pipeline processes light curves, generates subtraction images to highlight transient events, and uses a Convolutional Neural Network (CNN) from the ROBOT (Removal of Bogus Transients) pipeline to classify candidates as either real or bogus. This automated approach allowed us to systematically explore sub-minute transients for the first time, leading to the discovery of two high-quality candidates.

## Key Features
 * Automated Candidate Selection: Identifies light curves with single detections that fit the "out-of-nowhere" criteria for sub-minute transients.
 * Image Processing: Generates template and subtraction images for robust candidate verification.
 * Machine Learning Classification: Utilizes a pre-trained CNN to assign a "real/bogus" score to each candidate, dramatically reducing the number of false positives.
 * Artifact Filtering: Implements several checks to rule out common false positives such as electronic artifacts, cosmic rays, and satellite glints.

## The Discovery Pipeline
The pipeline follows a multi-stage process to identify genuine transient candidates:
1. Light Curve Selection: From a dataset of 671,761 light curves, we select 385,775 candidates that have a single detection and are not at the beginning or end of an observation sequence.
2. Image Subtraction: For each candidate, a template image is created from exposures taken on the same night. This template is then subtracted from the candidate image to produce a difference image, highlighting the transient.
3. ML-Powered Filtering: The template, science, and subtraction images are passed to the ROBOT CNN classifier. We use a conservative decision boundary of 0.06 to minimize the false-negative rate, resulting in 5,477 filtered candidates.
4. Manual Inspection & Vetting: The final candidates are manually inspected to identify the most promising ones and to perform further analysis, including multi-wavelength counterpart searches and rejection of satellite debris glints.

## License
This project is licensed under the MIT License.

## Acknowledgments
This research was conducted as part of the Deeper, Wider, Faster (DWF) programme and was funded by the Australian Research Council Centre of Excellence for Gravitational Wave Discovery (OzGrav). We also acknowledge the use of data from the Dark Energy Camera (DECam).

