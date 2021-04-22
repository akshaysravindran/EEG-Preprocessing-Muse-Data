# EEG-Preprocessing-Muse-Data
Preprocessing steps to handle artifacts in EEG collected using Muse headset


This is a script based preprocessing flowchart to handle artifacts in EEG collected using the Muse EEG system

**The most important thing to note is that these steps are not universal. It really depends on the study, what specifically you are looking for and how your data looks like. Then step zero is always looking at the data in raw form. This is a general framework which might work for most studies, but some of these steps might not be ideal for certain studies


The overall flowchart I typically prefer for the muse system is given below <br/><br/><img src='/images/Muse_flowchart.png'> <br/><br/>


Here, in addition to regular filtering, we could use Artifact Subspace Reconstruction to remove a large portion of the burst artifacts. In studies involving minimal movements, a carefully selected threshold can serve as a good means to remove eye blnks. Typically, people using Muse might either neglect delta band to minimize ocular contribution or reject epochs/ segments with eyeblinks. However, this is extremely laborious and time consuming. It also removes a lot of useful data. If we select the threshold for ASR well, we can clean eye blinks while not affecting the cleaner section of the data


The overall flowchart I typically prefer for brainvision 32-64 channel system is given below <br/><br/><img src='/images/flowchart.png'> <br/><br/>
<br/> The figure below shows how a carefully selected threshold can remove blinks while not rejecting cleaner part. The figure contains EEG during an eyes closed and eyes open condition.
As can be seen, the blinks are removed whereas the section during eyes closed where there are no blinks, the data is preserved

<br/><img src='/images/ASR_cleaning.png'>


<br/> Here it is very critical that we select the threshold well. For e.g. a conservative threshod of say 5 standard deviation distorts/ reconstructs even the cleaner part of the data


<br/><img src='/images/ASR_cleaned_5.png'>
<br/> On the other hand, a lax threshold ends up not cleaning the data well. Many artifacts are retained.
<br/><img src='/images/ASR_cleaned_50.png'>

<br/> The default of value is never a good value. Make sure the value is ideal for your data, do empirical visual testing to select an ideal threshold for your data
<br/><img src='/images/ASR_cleaned_15.png'>

