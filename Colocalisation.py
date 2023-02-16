############################################################################################################################### 
# import packages required
import numpy as np
import cv2
import matplotlib.pyplot as plt
from matplotlib import path
import scipy.interpolate as scinterp
from scipy import fftpack
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
from tqdm import tqdm

############################################################################################################################### 
#function for finding contours and give back the correct ones. Finds outline of the cell based on cytoplasmic or strong ER stain
def getContours(mask):
    #get contours
    ret,thresh = cv2.threshold(mask,0,255,cv2.THRESH_BINARY)        #0 can be adapted for further thresholding depending on stain
    closing = cv2.morphologyEx(thresh, cv2.MORPH_CLOSE,np.ones((10,10),np.uint8))               
    contours, hierarchy = cv2.findContours(closing,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
    outcon = []
    #find all contours from outline image and filter for the appropriate ones. 
    for i in range(0,len(contours)):
        area = cv2.contourArea(contours[i])
        if hierarchy[0,i,3] == -1 and area > 10000:
            outcon.append(contours[i])
    return outcon

###############################################################################################################################
#function for getting center of contour
def getCenter(contour):
    M = cv2.moments(contour)          #cv2.moments returns key pixel "moments" = m values. 
    cX = int(M["m10"] / M["m00"])
    cY = int(M["m01"] / M["m00"])       #centriod calc using m values             
    return [cX, cY]

####################################################################################################################################
def detectParticles(img, img_crop, nperc):
    perc = 0.01                                   #perc variable, n perc is number of sampling i.e. 50% = 0.5. perc is how smooth
    imr,imc = img.shape                           #find image shape
    #interpolate 3d image
    zim = img.copy()
    r,c = zim.shape
    x = np.zeros_like(zim)
    y = np.zeros_like(zim)
    for i in range(0,r):
        y[i,:] = i
    for i in range(0,c):
        x[:,i] = i
    nskip = int(imr*imc*nperc)                  #use percentage of image
    xsparse = x.flatten()[::nskip]              #flatten to an array so only include sampled pixels
    ysparse = y.flatten()[::nskip]
    xnew = x
    ynew = y
    zsparse = zim.flatten()[::nskip]

    ######Radial basis function interpolation
    rbf = scinterp.Rbf(xsparse,ysparse,zsparse,epsilon=(perc*imr+perc*imc)/2,smooth=(perc*imr+perc*imc)/2)  #scipi package interpolates 2d surface
    zinterp = rbf(xnew,ynew)                    #interpolated z values i.e. intensities 

    ##process image for blob detection
    imgx = zim-zinterp                          #this should remove background depending on the rbf variables
    imgx[imgx<0] = 0                            #make negative values = 0 
    # img2 += np.abs(np.amin(img2))
    # img2 *= 255./np.amax(img2)

    #image fft filtering
    imfft = fftpack.fft2(imgx)                  #high frequency band pass filter ... remove high freq i.e. remove high repeating patterns i.e. backgroundy things
    
    ##keep only fraction of coefficients
    keep_frac = 0.2
    imfft2 = imfft.copy()
    [row,col] = imfft2.shape

    ##high frequency cut out
    imfft2[int(row*keep_frac):int(row*(1-keep_frac)),:] = 0
    imfft2[:,int(col*keep_frac):int(col*(1-keep_frac))] = 0

    ##inverse fft to remake image
    im_ifft = fftpack.ifft2(imfft2).real
    im_ifft[im_ifft<0] = 0

    ##contour detection on fft of image
    imfftth = im_ifft.copy()
    imfftth = img_crop*imfftth*(255.0/np.amax(img_crop*imfftth))
    contoursfft = measure.find_contours(imfftth,np.mean(imfftth)+np.std(imfftth))   #use inverse fft to find contours
    imgx = np.array(imgx)
    return imgx,rbf,contoursfft

###############################################################################################################################
# #file paths 
path = './ '                                                        #File path where your images are
outpath = path+'Output/'
isExist = os.path.exists(outpath)                                   #File path where you want the output to go 
if not isExist:
    os.mkdir(outpath)
files = os.listdir(path)

###############################################################################################################################
#  Upload images for analysis. File naming is key. Create excel spreadsheets for output. 
mask_paths = []
img1_paths = []
img2_paths = []
img3_paths = []

for names in files:
    if names.startswith("1"):                   #Image from Channel 1 i.e. GFP signal
        img1_paths.append(names) 
for names in files:
    if names.startswith("2"):                   #Image from Channel 2 i.e. antibody stain
        img2_paths.append(names) 
### repeat this motif to analyse more channels 
for names in files:
    if names.startswith("Mask"):                #Image from Channel which best identifies cells. 
        mask_paths.append(names) 

img1_paths = natsorted(img1_paths)
img2_paths = natsorted(img2_paths)
mask_paths = natsorted(mask_paths)

with open(path+'data.csv','w') as file:
    file.write('image, cell, img1_img2, img2_img1, img1_puncta, img2_puncta, blob_1, blob_2, blob_ab\n') #create csv excel sheel with these columns. 
    # file.write('image, cell, img1_img2, img2_img1, img1_img3, img3_img1, img2_img3, img3_img2, img1_puncta, img2_puncta, img3_puncta, blob_1, blob_2, blob_3, blob_ab, blob_ac, blob_bc\n') #or these columns for triple comparison
file.close()
writer = pd.ExcelWriter(path+'area_pandas.xlsx', engine='xlsxwriter') #create an xslx excel spreadsheet to hold the detected particle areas
result_area = pd.DataFrame()

###############################################################################################################################
# # Main processing code. 
cellcount = 0                                               
for pcount in tqdm(range(0,len(img1_paths))):                  #upload images to code 
    img1 = cv2.imread(path+img1_paths[pcount],-1) 
    img2 = cv2.imread(path+img2_paths[pcount],-1) 
    # img3 = cv2.imread(path+img3_paths[pcount],-1) 
    mask = cv2.imread(path+mask_paths[pcount],0) 

    img1 = (img1/256).astype('uint8')                       #convert to 8 bit to reduce processing time. 
    img2 = (img2/256).astype('uint8')
    # img3 = (img3/256).astype('uint8')

    ###############################################################################################################################
    outcon =getContours(mask)                               #call cell segmentation function and retrieve cell outlines

    ###############################################################################################################################
    #cycle through each contour
    for i in range(0,len(outcon)):
        oneoutline = np.zeros_like(mask)                            #return array with same shape as image, this is empty
        cv2.fillPoly(oneoutline, pts =[outcon[i]], color=(1))       #fills the array with the cell shape.
        [cc,cr] = getCenter(outcon[i])                              #defined function above returns centriod          
        c,r,w,h = cv2.boundingRect(outcon[i])                       #return x,y = top left coordinate and w,h = width and height 
        ccreal = cc
        crreal = cr
        cc = cc - c
        cr = cr - r
        combine = np.copy(img1)                                     #create a copy of image so original is not destroyed. 
        imagetoshow = cv2.drawContours(combine, [outcon[i]], 0, (255,255,0), 3)            #draw cell outline onto image

        #reveal me for for 1st sanity check 
        # 1st check point
        # plt.figure(figsize=(15,8))                                #show me the outlined cell. If not perfect tweak segmentation parameters. 
        # plt.title(pcount)                                         #if cells are closely attached to each other, this will not be able to separate the two very successfully, when this occurs manual editting via FIJI can overcome this
        # plt.imshow(imagetoshow)
        # plt.show()

    #once happy with cell segmentation, continue to the following stage for intracellular structure segmentation
    ###############################################################################################################################
    #crop images to a box around the cell 
        img1_temp = img1[r:r+h,c:c+w]
        img2_temp = img2[r:r+h,c:c+w]
        # img3_temp = img3[r:r+h,c:c+w]
        oneoutline = oneoutline[r:r+h,c:c+w]                        

    ############################################################################################################################### 
    ##depending on the intensity and contrast of signal, structure detection can be alterred by classical filtering to aid further detection

        # otsu = threshold_otsu(img3_temp)          #otsu tends to be a good filter option, however others can be trialled here. 
        # img3_temp = img3_temp > otsu

        # otsu = threshold_otsu(img2_temp)
        # img2_temp = img2_temp > otsu

        # otsu = threshold_otsu(img1_temp)
        # img1_temp = img1_temp > otsu

    ###############################################################################################################################
    #the following section detects and segments structures within the cell. 

    #img1 
        segimg1 = oneoutline*img1_temp*(255.0/np.amax(oneoutline*img1_temp))
        imgx,rbf,img1contours = detectParticles(img1_temp, oneoutline, 0.01)       #0.01 = sampling density. varying this value will alter sensativity of detection. 
        
        # plt.figure(figsize=(15,8))                                               #visual check of structure detection. low throughput, keep silent for high throughput and check later. 
        # plt.imshow(img1_temp), plt.grid(False)
        # for b in img1contours:
        #     plt.plot(b[:,1],b[:,0],'b', alpha = 0.5, linewidth=1)
        # plt.show()
        
        mask_img1 = np.zeros_like(segimg1)                                          #cycle through structure contours are filter for characteristics
        for cnt in img1contours: 
            c2in = []
            for icn in cnt:
                c2in.append([int(icn[1]),int(icn[0])])
            c2in = [np.array(c2in,dtype='int')]
            area = cv2.contourArea(c2in[0])
            if area > 10:                                                           #i.e. structures with an area greater than 10 pixels to limit noise detection. 
                cv2.fillPoly(mask_img1,c2in,color=(1))
        a = mask_img1.astype('float')
       
        label_img_1 = label(mask_img1)                                              #calculate the size of each structure and number of structures per cell. 
        regions_1 = regionprops(label_img_1, intensity_image = mask_img1)
        blob_1 = []
        area_1 = []
        for props in regions_1:
            blob_1.append(props.label)
            Ai = props.area
            area_1.append(Ai)
        blob_1 = len(blob_1)
      
    ###############################################################################################################################
    #as above but for second channel
  #img2
        segimg2 = oneoutline*img2_temp*(255.0/np.amax(oneoutline*img2_temp))
        imgx,rbf,img2contours = detectParticles(img2_temp, oneoutline, 0.01) 

        # plt.figure(figsize=(15,8))
        # plt.imshow(img2_temp), plt.grid(False)
        # for b in img2contours:
        #     plt.plot(b[:,1],b[:,0],'b', alpha = 0.5, linewidth=1)
        # plt.show()
        
        mask_img2 = np.zeros_like(segimg2)
        for cnt in img2contours: 
            c2in = []
            for icn in cnt:
                c2in.append([int(icn[1]),int(icn[0])])
            c2in = [np.array(c2in,dtype='int')]
            area = cv2.contourArea(c2in[0])
            if area > 10:
                cv2.fillPoly(mask_img2,c2in,color=(1))
        
        b = mask_img2.astype('float')

        label_img_2 = label(mask_img2)
        regions_2 = regionprops(label_img_2, intensity_image = mask_img2)
        blob_2 = []
        area_2 = []
        
        for props in regions_2:
            blob_2.append(props.label)
            Ai = props.area
            area_2.append(Ai)
        blob_2 = len(blob_2)
        
    ###############################################################################################################################
    #as above but for third channel 
    # # #img3
    #     segimg3 = oneoutline*img3_temp*(255.0/np.amax(oneoutline*img3_temp))
    #     imgx,rbf,img3contours = detectParticles(img3_temp, oneoutline, 0.003) #optimised


    # #     # plt.figure(figsize=(15,8))
    # #     # plt.imshow(img3_temp), plt.grid(False)
    # #     # for b in img3contours:
    # #     #     plt.plot(b[:,1],b[:,0],'b', alpha = 0.5, linewidth=1)

    # #     # plt.show()
    
    #     mask_img3 = np.zeros_like(segimg3)
    #     for cnt in img3contours: 
    #         c2in = []
    #         for icn in cnt:
    #             c2in.append([int(icn[1]),int(icn[0])])
    #         c2in = [np.array(c2in,dtype='int')]
    #         area = cv2.contourArea(c2in[0])
    #         if area > 10:
    #             cv2.fillPoly(mask_img3,c2in,color=(1))
    #     c = mask_img3.astype('int')

    #     label_img_3 = label(mask_img3)
    #     regions_3 = regionprops(label_img_3, intensity_image = segimg3)
    #     blob_3 = []
    #     area_3 = []
        
    #     for props in regions_3:
    #         blob_3.append(props.label)
    #         Ai = props.area
    #         area_3.append(Ai)
    #     blob_3 = len(blob_3)
        
    ###############################################################################################################################
    #detect and segment all of these structures and then show them to me. Highthroughput = save images and check afterwards. 
        fig, axs = plt.subplots(2, 3, figsize=(15,8))
        axs[0,0].imshow(imagetoshow), axs[0,0].grid(False), axs[0,0].set_title("Image_"+str(pcount)+"_Cell_"+str(i)), axs[0,0].axis('off')
        axs[0,1].imshow(mask_img1), axs[0,1].grid(False), axs[0,1].set_title("Img1"), axs[0,1].axis('off')
        axs[0,2].imshow(img2_temp), axs[0,2].grid(False), axs[0,2].set_title("Img2"), axs[0,2].axis('off')
        axs[1,0].imshow(mask_img2), axs[1,0].grid(False), axs[1,0].set_title("Img2"), axs[1,0].axis('off')
        # axs[1,1].imshow(img3_temp), axs[1,1].grid(False), axs[1,1].set_title("Img3"), axs[1,1].axis('off')
        # axs[1,2].imshow(mask_img3), axs[1,2].grid(False), axs[1,2].set_title("Img3"), axs[1,2].axis('off')
        plt.savefig(outpath+'Output_Image_'+str(pcount)+'_Cell_'+str(i)+'.tif')
        plt.close()
        # plt.show()
    ###############################################################################################################################
    #happy with the segmented intracellular structures? Time to analyse them! 
    
        a_1 = np.count_nonzero(mask_img1) #img1 puncta                          #how many pixels in channel 1 are positive?
        b_1 = np.count_nonzero(mask_img2) #img2 puncta
        # c_1 = np.count_nonzero(mask_img3) #img3 puncta

        bothab = np.multiply(b,a)                                               #how many pixels are positive in both channels? 
        # bothac = np.multiply(a,c)
        # bothbc = np.multiply(c,b)       

        label_bothab = label(bothab)
        regions_ab = regionprops(label_bothab, intensity_image = bothab)        #calculate number of overlapping structures and their area. 
        blob_ab = []
        area_ab = []
        for props in regions_ab:
            blob_ab.append(props.label)
            Ai = props.area
            area_ab.append(Ai)
        blob_ab = len(blob_ab)
        # print(blob_ab)
        
        # label_bothac = label(bothac)                                          #repeat for all comparisons desired. 
        # regions_ac = regionprops(label_bothac, intensity_image = bothac)
        # blob_ac = []
        # area_ac = []
        # for props in regions_ac:
        #     blob_ac.append(props.label)
        #     Ai = props.area
        #     area_ac.append(Ai)
        # blob_ac = len(blob_ac)

        # label_bothbc = label(bothbc)
        # regions_bc = regionprops(label_bothbc, intensity_image = bothbc)
        # blob_bc = []
        # area_bc = []
        # for props in regions_bc:
        #     blob_bc.append(props.label)
        #     Ai = props.area
        #     area_bc.append(Ai)
        # blob_bc = len(blob_bc)

        overlap_ab = np.count_nonzero(bothab) #overlap a b                     #calculate percentage of channel 1 signal which is overlapping
        img12_coloc = overlap_ab / a_1 *100
        img21_coloc = overlap_ab / b_1 *100                                    #calculate percentage of channel 2 signal which is overlapping 

        # overlap_ac = np.count_nonzero(bothac)
        # img13_coloc = overlap_ac / a_1 *100
        # img31_coloc = overlap_ac / c_1 *100

        # overlap_bc = np.count_nonzero(bothbc)
        # img23_coloc = overlap_bc / b_1 *100
        # img32_coloc = overlap_bc / c_1 *100

        with open(path+"data.csv",'a') as f:                                    #input intracellular structure data into excel sheet
            lbl = pcount                                                        #image label
            cell = i                                                            #cell number
            img1_img2 = img12_coloc                                             #percentage of channel 1 colocalised to channel 2
            img1_puncta = a_1                                                   #number of structures in channel 1
            img2_img1 = img21_coloc                                             #percentage of channel 2 colocalised to channel 1 
            img2_puncta = b_1                                                   #number of structures in channel 2
            # img1_img3 = img13_coloc
            # img3_puncta = c_1
            # img3_img1 = img31_coloc
            # img2_img3 = img23_coloc
            # img3_img2 = img32_coloc
            blob1 = blob_1                                                      #number of distinct structures in channel 1
            blob2 = blob_2                                                      #number of distinct structures in channel 2
            # blob_3 = blob_3
            blobab = blob_ab                                                    #number of overlapping structures
            # blob_ac = blob_ac
            # blob_bc = blob_bc
             

            f.write(str(lbl)+","+str(cell)+","+str(img1_img2)+","+str(img2_img1)+","+str(img1_puncta)+","+str(img2_puncta)+","+str(blob1)+","+str(blob2)+","+str(blobab)+"\n")
            # f.write(str(lbl)+","+str(cell)+","+str(img1_img2)+","+str(img2_img1)+","+str(img1_img3)+","+str(img3_img1)+","+str(img2_img3)+","+str(img3_img2)+","+str(img1_puncta)+","+str(img2_puncta)+","+str(img3_puncta)+","+str(blob_1)+","+str(blob_2)+","+str(blob_3)+","+str(blob_ab)+","+str(blob_ac)+","+str(blob_bc)+"\n")
        
        #### data entry for the area of the segmented structures. 
        area_1 = pd.DataFrame(area_1, columns = ['Area_1'])                     #list of structures and their area per cell
        area_2 = pd.DataFrame(area_2, columns = ['Area_2'])
        # area_3 = pd.DataFrame(area_3, columns = ['Area_3'])
        area_ab = pd.DataFrame(area_ab, columns = ['Area_ab'])
        # area_ac = pd.DataFrame(area_ac, columns = ['Area_ac'])
        # area_bc = pd.DataFrame(area_bc, columns = ['Area_bc'])
        image = pd.Series([pcount])
        cell = pd.Series([i])
        areas = [image, cell, area_1, area_2, area_ab]
        # areas = [image, cell, area_1, area_2, area_3, area_ab, area_ac, area_bc]

        result = pd.concat(areas, axis=1)
        result_area = result_area.append(result)
    cellcount+=1                                                                #cycle through all cells in image
    print("Finished image_"+str(pcount))                                        #cycle through all images in database
result_area.to_excel(writer)
writer.save()
