# Declining food availability and habitat shifts drive community responses to marine hypoxia

The physoxiaSSM folder contains all code necessary to run the model found in my paper.  Briefly, the code alters the baseline structure provided by the R package mizer (https://github.com/sizespectrum/mizer) by introducing oxygen dependence of several biological rates.  It allows species exposed directly to hypoxia to prioritize among physiological and behavioral rates, and compares results to the original model -- (Duskey, Casini, Limburg, and Gårdmark. In review. Declining food availability and habitat shifts drive community responses to marine hypoxia.) -- in which species were not allowed to prioritize.

## Abstract

Marine hypoxia has had major consequences for both economically and ecologically critical fish species around the world.  As hypoxic regions continue to grow in severity and extent, we must deepen our understanding of mechanisms driving population and community responses to major stressors.  It has been shown that food availability and habitat use are the most critical components of impacts on individual fish leading to observed outcomes at higher levels of organization.  However, differences within and among species in partitioning available energy for metabolic demands -- or metabolic prioritization -- in response to stressors are often ignored.  Here, I use both a multispecies size spectrum model and a meta-analysis to explore evidence in favor of metabolic prioritization in a community of commercially important fish species in the Baltic Sea.  Modeling results suggest that metabolic prioritization is an important component of the individual response to hypoxia, that it interacts with other components to produce realistic community dynamics, and that different species may prioritize differently.  It is thus suggested that declines in feeding activity, assimilation efficiency, and successful reproduction -- in addition to low food availability and changing habitat use -- are all important drivers of the community response to hypoxia.  Meta-analysis results also provide evidence that the dominant predator in the study system prioritizes among metabolic demands, and that these priorities may change as oxygen declines.  Going forward, experiments and models should explore how differences in priorities within and among communities drive responses to environmental degradation.  This will help management efforts to tailor recovery programs to the physiological needs of species within a given system.

## Data

All data used to perform the analysis are contained within the Data folder in this repository.  Their origin is briefly described in the .keep text file in the Data folder.

## Requirements

The size spectrum models contained within the code are run in the R package mizer.  Description, installation, updates, and helpful tips are provided at the following website: https://sizespectrum.org/mizer/

## Usage

Follow these steps prior to first time usage to ensure the code runs properly on your machine:

1. Download the "physoxiaSSM-main" zipped folder from Github in its entirety
2. Unzip and place "physoxiaSSM-main" in your preferred directory
3. Navigate to the Code folder and open the file "depend.R" in R or RStudio
4. Confirm that all listed packages are installed on your machine; uncomment and run lines corresponding to packages you have not yet installed, then re-comment
5. Change the "mypath" variable to the file path which contains the hypoxiaSSM-main folder
6. Save and close "depend.R"

Upon all subsequent uses, we recommend running depend.R first, or pasting your own working directory over the code in each relevant file you use.

NOTE: Nearly every folder contains a text file called ".keep" which describes the files and subfolders within.  Read these for a more detailed understanding of the repository and how to use it.

## Contact me

If you have questions or concerns regarding this code, or would like help in re-formatting it for your own use, please do not hesitate to contant the corresponding author at:

elizabeth.duskey@slu.se

## License

Copyright 2022 Elizabeth Duskey

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
