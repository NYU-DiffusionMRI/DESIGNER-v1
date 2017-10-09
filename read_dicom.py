import os, sys
import dicom
import numpy as np

dpath = '/mccoy/input/test-file-array'
newDcmBase='/mccoy/output/hello'
dcmlist = os.listdir(dpath)

info = dicom.read_file(os.path.join(dpath,dcmlist[0]))
if info.ImageType.index("DIFFUSION") or info.ImageType.index("FMRI"):
    istimeseries = True
if info.ImageType.index("MOSAIC"):
    ismosaic = True
datasize = [info.Rows, info.Columns]
matrixsize = list(info.AcquisitionMatrix[i] for i in [0, 3])
NumberOfTiles = [int(x)/int(y) for x,y in zip(datasize, matrixsize)]
SeriesDescription = info.SeriesDescription
SeriesID = info.SeriesInstanceUID
SeriesNumber = info.SeriesNumber

newSeriesID = SeriesID[:-16] + str(np.int(np.floor(1000000000 + 8999999999*np.random.random()))) + SeriesID[-6:]
newSeriesDescription = 'Denoised_' + SeriesDescription
newSeriesNumber = SeriesNumber + 200

#DcmBaseName = os.path.basename(os.path.normpath(dpath))
#DcmBasePath = os.path.dirname(os.path.abspath(DcmBaseName))
#newDcmBase = os.path.join(DcmBasePath,'DN_' + DcmBaseName)
#if not os.path.exists(newDcmBase):
#    os.makedirs(newDcmBase)

if istimeseries:
    newSOPUID = []
    I = np.zeros((matrixsize[0], matrixsize[1], np.prod(NumberOfTiles), len(dcmlist)))
    
    for i in range(0, len(dcmlist)):
        info = dicom.read_file(os.path.join(dpath,dcmlist[i]))
        SOPUID = info.SOPInstanceUID
        newSOPUID.append(SOPUID[:-16] + str(np.int(np.floor(1000000000 + 8999999999*np.random.random()))) + SOPUID[-6:])
        
        data = info.pixel_array
        img = np.reshape(data, (datasize[0], datasize[1]), order='F')
        
        c = 0
        if ismosaic:
            for j in range(1, NumberOfTiles[0]+1):
                for k in range(1, NumberOfTiles[1]+1):
                    I[:,:,c,i] = img[(j-1)*matrixsize[0]:(j)*matrixsize[0], (k-1)*matrixsize[1]:(k)*matrixsize[1]]
                    c = c+1

#   do the denoising
#   Idn = MP(I,'full',kernel)
    Idn = I

    for i in range(0, len(dcmlist)):
        info = dicom.read_file(os.path.join(dpath,dcmlist[i]))
        info.SeriesDescription = newSeriesDescription
        info.SeriesNumber = newSeriesNumber
        info.SeriesID = newSeriesID
        info.SOPInstanceUID = newSOPUID[i]
        c = 0
        if ismosaic:
            for j in range(1, NumberOfTiles[0]+1):
                for k in range(1, NumberOfTiles[1]+1):
                    img[(j-1)*matrixsize[0]:(j)*matrixsize[0], (k-1)*matrixsize[1]:(k)*matrixsize[1]] = Idn[:,:,c,i]
                    c = c+1
        info.pixel_array = img
        info.PixelData = info.pixel_array.tostring()
        info.save_as(os.path.join(newDcmBase,'IM_00'+str(i)+'.dcm'))

sys.exit(0)





