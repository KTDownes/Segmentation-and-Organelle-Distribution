############################################################################################################################### 
# import packages required
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
import xlsxwriter 
from tqdm import tqdm


##############################################################################################
# #function for finding contours and give back the correct ones
def getContours(gfp_outline):
   #get contours
    ret,thresh = cv2.threshold(gfp_outline,127,255,cv2.THRESH_OTSU)     #parameters can be changed to suit the dataset
    closing = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE,np.ones((10,10),np.uint8))
    contours, hierarchy = cv2.findContours(closing,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    outcon = []
    #find all contours from outline image and filter for the proper ones
    for i in range(0,len(contours)):
        area = cv2.contourArea(contours[i])
        if hierarchy[0,i,3] == -1 and area > 10000:
            # if area/cv2.arcLength(contours[i],True)>50:
            outcon.append(contours[i])
    return outcon

##############################################################################################
#file paths 
raw = []
path = './'                     #File path where your images are
files = os.listdir(path)

##############################################################################################
#  Upload images for analysis. File naming is key. 

gfp_paths = []
bio_paths = []
myc_paths = []
mask_paths = []
combined_paths = []

#change to the identifiers of your data set
for names in files:
    if names.endswith("_GFP100ms.tif"):             #channel 1
        gfp_paths.append(names)
    if names.endswith("_Biotin100ms.tif"):          #channel 2
        bio_paths.append(names)
    if names.endswith("_MYC100ms.tif"):             #channel 3
        myc_paths.append(names)
    if names.endswith("_Mask100ms.tif"):            #channel with best signal for identifying cells
        mask_paths.append(names)
    if names.startswith("Composite_"):              #RGB Image 
        combined_paths.append(names)

gfp_paths = natsorted(gfp_paths)
bio_paths = natsorted(bio_paths)
myc_paths = natsorted(myc_paths)
mask_paths = natsorted(mask_paths)
combined_paths = natsorted(combined_paths)

#import images using file paths 
cellcount = 0
for pcount in tqdm(range(0,len(gfp_paths))): 
    gfp = cv2.imread(path+gfp_paths[pcount],-1)
    gfp_outline = cv2.imread(path+mask_paths[pcount],0)
    myc = cv2.imread(path+myc_paths[pcount],-1)
    bio = cv2.imread(path+bio_paths[pcount],-1)
    combine_original = cv2.imread(path+combined_paths[pcount],-1)
   
##############################################################################################
    outcon = getContours(gfp_outline)           #call cell segmentation function and retrieve cell outlines

##############################################################################################
    #cycle through each contour
    for i in range(0,len(outcon)):
        oneoutline = np.zeros_like(gfp_outline)                     #return array with same shape as gfp_outline, this is empty
        cv2.fillPoly(oneoutline, pts =[outcon[i]], color=(1))       #fills the array with the cell shape.
        combine = np.copy(combine_original)
        imagetoshow = cv2.drawContours(combine, [outcon[i]], 0, (0,255,0), 3)
##############################################################################################
        # #1st check point
        # plt.figure(figsize=(15,8))        #show me the outlined cell. If not perfect tweak segmentation parameters. 
        # plt.title(pcount)                 #if cells are closely attached to each other, this will not be able to separate the two very successfully, when this occurs manual editting via FIJI can overcome this
        # plt.imshow(imagetoshow)
        # plt.show()

    #once happy with cell segmentation, continue to the following stage 
 
##############################################################################################
    #calculate area, average intensiy etc for each cell     
        mask = np.zeros_like(gfp_outline)
        cv2.fillPoly(mask,[outcon[i]],color=(1))
        label_img = label(mask)

        cells = pd.DataFrame(measure.regionprops_table(label_img, myc, cache = True, properties=('area','mean_intensity')))
        gfp_cells = pd.DataFrame(measure.regionprops_table(label_img, gfp, cache = True, properties=('area','mean_intensity')))
        bio_cells = pd.DataFrame(measure.regionprops_table(label_img, bio, cache = True, properties=('area','mean_intensity')))
        
        cells.insert(2, 'GFP_Intensity', gfp_cells['mean_intensity'])
        cells.insert(3, 'Bio_Intensity', bio_cells['mean_intensity'])
        cells.insert(0, 'Image', pcount)
        cells.insert(1,'Label', i)
        raw.append(cells)
        
#combine cell data for entire dataset  
data = pd.concat(raw)

workbook = xlsxwriter.Workbook('./Intensity.xlsx')
worksheet = workbook.add_worksheet()
data.to_excel('Intensity.xlsx')

