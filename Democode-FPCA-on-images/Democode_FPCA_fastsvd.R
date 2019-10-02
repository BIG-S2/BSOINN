#This demo code is for reading imaging data and performing FPCA on the images. Please refer to Section 2.1 of the manuscript for the description of FPCA.
#Please adjust the code before use.

##load library
library(AnalyzeFMRI)  # this package (package version 1.1-20) contains useful functions that read imaging data. Note that the analyzed RAVENS image consists of two files, ".hdr" and ".img"; the two files should have the same name to represent one image.
library(corpcor)      # this package (package version 1.6.9) contains the function "fast.svd" that performs SVD on a matrix

##read imaging data
##suppose all image files (in formats of ".hdr" and ".img") are stored in the folder "D:\\images"
image_length = 300^2        #the dimensionality of an image
image_file_names = list.files("D:\\images") #read the file names of the images
image_file_names = image_file_names[2*(1:(length(image_file_names)/2))]  #get rid of the rundundancy of file names 

X = matrix(0,nrow = length(image_file_names), ncol = image_length)    #store images

for(i in 1:length(names)){
temp_image = f.read.nifti.volume(paste("D:\\images\\",names[i],sep = "")) # 'f.read.nifti.volume' reads images into memory; the function is provided in the package "AnalyzeFMRI".
X[i,] = as.vector(temp_image)  # store images
}

X = scale(X,center = T, scale=F)  # centralizing images

##FPCA on X using SVD
eigen_svd = fast.svd(t(X),tol=0.0001)  #  'fast.svd' performs SVD on X; the function is provided in the package "corpcor".
eigenimage_est = t(eigen_svd$u)   # Eigenimages are obtained. 
XI = t(t(eigen_svd$v)*eigen_svd$d)  # Eigenscores are obtained.

XI = t(t(XI)/apply(XI,2,sd)) # an optional procedure: the eigenscores can be standardized before the analysis. The first K columns of XI are useful eigenscores in our model.


