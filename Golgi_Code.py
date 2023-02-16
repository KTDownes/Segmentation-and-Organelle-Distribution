import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib import path
import scipy.interpolate as scinterp
from scipy import fftpack
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os #basic
import PIL #python imaging library
from PIL import Image 
import pandas as pd 
from skimage import data, img_as_float #converts image to floating point - will alter intensity values.
from skimage import exposure
import seaborn as sns; sns.set()
from skimage import io 
from skimage import transform, measure
from skimage.measure import regionprops, regionprops_table, label
from skimage import exposure 
from skimage.filters import try_all_threshold, threshold_otsu, threshold_yen, threshold_triangle
from skimage import morphology 
from scipy import ndimage as ndi
from skimage.segmentation import chan_vese, felzenszwalb, slic, quickshift, watershed,mark_boundaries
from matplotlib import cm
from matplotlib import interactive
import xlsxwriter 
from natsort import natsorted
import math 
from matplotlib.font_manager import FontProperties 

###############################################
# #function for finding contours and give back the correct ones
def getContours(gfp_outline):
   #get contours
    ret,thresh = cv2.threshold(gfp_outline,127,255,cv2.THRESH_OTSU)             #parameters can be changed to suit the dataset
    closing = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE,np.ones((10,10),np.uint8))
    contours, hierarchy = cv2.findContours(closing,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    outcon = []
    #find all contours from outline image and filter for the proper ones
    for i in range(0,len(contours)):
        area = cv2.contourArea(contours[i])
        if hierarchy[0,i,3] == -1 and area > 10000:
            outcon.append(contours[i])
    return outcon

# ###############################################
# #function for getting center of contour
def getCenter(contour):
    M = cv2.moments(contour)                       #cv2.moments returns key pixel "moments" = m values. 
    cX = int(M["m10"] / M["m00"])
    cY = int(M["m01"] / M["m00"])                  #centriod calc using m values             
    return [cX, cY]

################################################function for detecting particles
#function to return contours from the image
def detectParticles(img, img_crop, nperc):
    perc = 0.01              #perc variable, n perc is number of sampling i.e. 50% = 0.5. perc is how smooth
    imr,imc = img.shape      #find image shape
    #interpolate 3d image
    zim = img.copy()
    r,c = zim.shape
    x = np.zeros_like(zim)
    y = np.zeros_like(zim)
    for i in range(0,r):
        y[i,:] = i
    for i in range(0,c):
        x[:,i] = i
    nskip = int(imr*imc*nperc)      #use percentage of image
    xsparse = x.flatten()[::nskip]  #flatten to an array so only include sampled pixels
    ysparse = y.flatten()[::nskip]
    print('length of sparse:',len(xsparse))
    xnew = x
    ynew = y

    zsparse = zim.flatten()[::nskip]
    ######Radial basis function interpolation
    rbf = scinterp.Rbf(xsparse,ysparse,zsparse,epsilon=(perc*imr+perc*imc)/2,smooth=(perc*imr+perc*imc)/2)  #scipi package interpolates 2d surface
    zinterp = rbf(xnew,ynew)    #interpolated z values i.e. intensities   

    #process image for blob detection
    img2 = zim-zinterp  #this should remove background depending on the rbf variables
    img2[img2<0] = 0    #make negative values = 0 
    # img2 += np.abs(np.amin(img2))
    # img2 *= 255./np.amax(img2)

    #image fft filtering
    imfft = fftpack.fft2(img2)  #high frequency band pass filter ... remove high freq i.e. remove high repeating patterns i.e. backgroundy things
    
    #keep only fraction of coefficients
    keep_frac = 0.2
    imfft2 = imfft.copy()
    [row,col] = imfft2.shape

    #high frequency cut out
    imfft2[int(row*keep_frac):int(row*(1-keep_frac)),:] = 0
    imfft2[:,int(col*keep_frac):int(col*(1-keep_frac))] = 0

    #inverse fft to remake image
    im_ifft = fftpack.ifft2(imfft2).real
    im_ifft[im_ifft<0] = 0

    #contour detection on fft of image
    imfftth = im_ifft.copy()
    imfftth = img_crop*imfftth*(255.0/np.amax(img_crop*imfftth))
    contoursfft = measure.find_contours(imfftth,np.mean(imfftth)+np.std(imfftth))   #use inverse fft to find contours
    img2 = np.array(img2)
    return img2,rbf,contoursfft

###############################################MAIN
#upload file paths 
path = './'                     #File path where your images are
outpath = path+'Output/'
isExist = os.path.exists(outpath)                                   #File path where you want the output to go 
if not isExist:
    os.mkdir(outpath)
files = os.listdir(path)

###############################################################################################################################
#  Upload images for analysis. File naming is key. Create excel spreadsheets for output. 
gfp_paths = []
gol_paths = []
myc_paths = []
combined_paths = []
for names in files:
    if names.endswith("_GFP100ms.tif"):
        gfp_paths.append(names)
    if names.endswith("GM130.tif"):
        gol_paths.append(names)
    if names.endswith("_Mask100ms.tif"):
        myc_paths.append(names)
    if names.startswith("Composite_"):
        combined_paths.append(names)

gfp_paths = natsorted(gfp_paths)
gol_paths = natsorted(gol_paths)
myc_paths = natsorted(myc_paths)
combined_paths = natsorted(combined_paths)

#create a result file
with open(path+'data_golgi.csv','w') as file:
    file.write('label, intens, y1, x1, Ai, r, cellnum,exptype\n')
file.close()

###############################################################################################################################
#import images using file paths 
cellcount = 0
for pcount in tqdm(range(0,len(gfp_paths))):                  #upload images to code 
    gfp = cv2.imread('./'+gfp_paths[pcount],-1)
    gfp_outline = cv2.imread('./'+myc_paths[pcount],0)
    gol = cv2.imread('./'+gol_paths[pcount],-1)
    combine_original = cv2.imread('./'+combined_paths[pcount],-1)

###############################################################################################################################
    outcon =getContours(gfp_outline)                               #call cell segmentation function and retrieve cell outlines

###############################################################################################################################
    #cycle through each contour
    for i in range(0,len(outcon)):
        oneoutline = np.zeros_like(gfp_outline)                     #return array with same shape as gfp_outline, this is empty
        cv2.fillPoly(oneoutline, pts =[outcon[i]], color=(1))       #fills the array with the cell shape.
        [cc,cr] = getCenter(outcon[i])                              #defined function above returns centriod 
        (x,y),maxradius = cv2.minEnclosingCircle(outcon[i])         #makes a circle which contains the cell mask
        c,r,w,h = cv2.boundingRect(outcon[i])                       #return x,y = top left coordinate and w,h = width and height 
        ccreal = cc
        crreal = cr
        cc = cc - c
        cr = cr - r
        combine = np.copy(combine_original)                          #create a copy of image so original is not destroyed.
        imagetoshow = cv2.drawContours(combine, [outcon[i]], 0, (0,255,0), 3)

##########################################################################################################################################
        # #1st check point
        # plt.figure(figsize=(15,8))
        # plt.title(pcount)
        # plt.imshow(imagetoshow)
        # plt.show()

    #once happy with cell segmentation, continue to the following stage for intracellular structure segmentation
###############################################################################################################################
        gol_temp = gol[r:r+h,c:c+w]
        gfp_temp = gfp[r:r+h,c:c+w]
        oneoutline = oneoutline[r:r+h,c:c+w]                        #create golgi copy and crop to cell outline 
        seggolpic = oneoutline*gol_temp*(255.0/np.amax(oneoutline*gol_temp))

##########################################################################################################################################
# #the following section detects and segments structures within the cell.       
        img2,rbf,golcontours = detectParticles(gol_temp, oneoutline, 0.01) #optimised

        # plt.figure(figsize=(15,8))
        # plt.imshow(gol_temp), plt.grid(False)
        # for b in golcontours:
        #     plt.plot(b[:,1],b[:,0],'b', alpha = 0.5, linewidth=1)
        # plt.show()
        
##########################################################################################################################################
        #take golgi contours and calculate desired features
        mask_golgi = np.zeros_like(seggolpic)
        for cnt in golcontours: 
            c2in = []
            for icn in cnt:
                c2in.append([int(icn[1]),int(icn[0])])
            c2in = [np.array(c2in,dtype='int')]
            cv2.fillPoly(mask_golgi,c2in,color=(1))

        label_img = label(mask_golgi)
        regions = regionprops(label_img, intensity_image = seggolpic)
        centroids_x = []
        centroids_y = []
        x_val_area = []
        y_val_area = []
        area = []
        intensity = []
        blob = []

        fig, axs = plt.subplots(1, 3, figsize=(15,8))
        axs[0].imshow(imagetoshow), axs[0].grid(False), axs[0].set_title("Image_"+str(pcount)+"_Cell_"+str(i)), axs[0].axis('off')
        axs[1].imshow(gol_temp), axs[1].grid(False), axs[1].set_title('Golgi'), axs[1].axis('off')
        axs[1].plot(golcontours[0][:,1],golcontours[0][:,0],'b', linewidth=1, label= ('Golgi Particles'))
        for b in golcontours:
            axs[1].plot(b[:,1],b[:,0],'b',linewidth=1)
        fontP = FontProperties()
        fontP.set_size('x-small')
        axs[1].legend(bbox_to_anchor=(0.5,-0.01), loc='upper center',prop=fontP)
        axs[2].grid(False), axs[2].imshow(mask_golgi), axs[2].set_title('Golgi Mask with Contours + Centroids'), axs[2].axis('off')
        
        for props in regions:
            blob.append(props.label)
            y1, x1 = props.weighted_centroid
            centroids_x.append(x1)
            centroids_y.append(y1)
            intens = props.mean_intensity
            intensity.append(intens)
            Ai = props.area
            area.append(Ai)
            #calc distance
            x_val_area.append(x1*Ai)
            y_val_area.append(y1*Ai)
        overall_xval_area = np.sum(np.array(x_val_area))
        overall_yval_area = np.sum(np.array(y_val_area))
        overall_area = np.sum(np.array(area))
        
        weighted_x = overall_xval_area/overall_area
        weighted_y = overall_yval_area/overall_area

        with open("data_golgi.csv",'a') as f:
            for a in range(len(centroids_x)):
                intense = int(intensity[a])
                y1 = int(centroids_y[a])
                x1 = int(centroids_x[a])
                Ai = int(area[a])
                name = int(blob[a])
                p = math.sqrt((weighted_x - centroids_x[a])**2+(weighted_y - centroids_y[a])**2)
                r = p/(2*maxradius)
                f.write(str(name)+","+str(intense)+","+str(y1)+","+str(x1)+","+str(Ai)+","+str(r)+","+str(cellcount)+","+"\n")
        cellcount+=1
        
        axs[2].plot(np.asarray(centroids_x), np.asarray(centroids_y), '.g', markersize= 4, label = ('Particle, Intensity-Weighted Centroids'))
        axs[2].plot(weighted_x, weighted_y, '*r', markersize = 10, label = 'Area Weighted Centroid') #area weighted centroid, might want to also weight by intensity 
        axs[2].legend(bbox_to_anchor=(0.5,-0.01), loc='upper center',prop=fontP)
        plt.savefig(outpath+"I092_Golgi_Image"+str(pcount)+"_"+str(i)+".png",dpi=300)
        plt.close()
        # plt.show()
        # print(pcount)
