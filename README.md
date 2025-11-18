# Adaptive Compression (ACA) Workflow for Large-scale Microscopy Data

This repository outlines a **two-stage process** for **rapid feature analysis (ACA)** on large-scale microscopy data, such as diSPIM data, LLSM data, confocal microscopy data, *et al*.

---

## 1. Sample Data

An example of a subregion intensity image from a whole stomach diSPIM dataset is provided in the `Example` folder. This image has already undergone **deconvolution** and **fusion** processing.

---

## 2. Analysis Workflow: Rapid Feature Analysis (ACA)

The rapid feature analysis using the ACA mapping process involves **image segmentation via U-Net** followed by **structural feature quantification in MATLAB**.

### Stage I: Image Segmentation (U-Net)

#### **Step 1: Segmentation Environment Setup (U-Net)**

Before segmentation, ensure your U-Net environment is properly configured. A high-performing, **pre-trained model** specifically optimized for diSPIM images is included in the `\U-Net_Pytorch\checkpoints` folder.

Follow the instructions below to set up the U-Net segmentation environment:

* **Requirements:**
    * Install **Python (3.6+)**, **CUDA**, and **PyTorch (1.13+)**.
    * Install dependencies:
        ```bash
        pip install -r requirements.txt
        ```

#### **Step 2: Image Segmentation**

1.  **Move the image** from the `Example` folder into the following directory:
    ```
    ACA mapping\U-Net_Pytorch\predict\
    ```
2.  Use the pre-trained model to segment the input image.
    **Execute the following command** in the terminal:
    ```bash
    python predict.py --model=checkpoints/LS_model.pth
    ```

The segmented image will be saved in the directory: `ACA mapping\U-Net_Pytorch\output\`

---

### 🛠️ Training

For specific types of microscopic images, the model can be trained using the following script.

**Usage:**

```bash
> python train.py -h
usage: train.py [-h] [--epochs E] [--batch-size B] [--learning-rate LR]
                [--load LOAD] [--scale SCALE] [--validation VAL] [--amp]

Train the UNet on images and target masks

optional arguments:
  -h, --help            show this help message and exit
  --epochs E, -e E      Number of epochs
  --batch-size B, -b B  Batch size
  --learning-rate LR, -l LR
                        Learning rate
  --load LOAD, -f LOAD  Load model from a .pth file
  --scale SCALE, -s SCALE
                        Downscaling factor of the images
  --validation VAL, -v VAL
                        Percent of the data that is used as validation (0-100)
  --amp                 Use mixed precision

```
> **Note on Training Parameters:**
>
> 1.  By default, the scale is **0.5**. To obtain better results (at the cost of higher memory usage), set it to **1**.
> 2.  **Automatic mixed precision** is recommended and available with the `--amp` flag.

### Data Structure Requirements

* **Data Folder Structure:** Input images and target masks should be placed in the `data/imgs` and `data/masks` folders, respectively. **Note** that these folders should **not** contain any sub-folders.
* **Image Format:** Images and masks must be **black and white**.
* **Custom Dataset:** You may use your own dataset if you ensure proper loading in `utils/data_loading.py`.

The **U-Net architecture** is based on the original paper by Olaf Ronneberger, Philipp Fischer, Thomas Brox:

> *U-Net: Convolutional Networks for Biomedical Image Segmentation* 


---

### Stage II: Feature Quantification (MATLAB)

#### **Step 3: Characteristics Quantification (ACA Mapping)**

After segmentation, the segmented output is processed using the **MATLAB main program 'main.m'** for feature calculation and mapping.

1.  **Run the MATLAB Script:**
    * Execute the main MATLAB program `'main.m'`.

2.  **Adjust Parameters:**
    Adjust the following parameters within the script's comments:

| Parameter | Description | Recommendation/Constraint |
| :--- | :--- | :--- |
| `'armdx'` & `'armdz'` | Define the size of the parameter calculation window. | **2 to 3** times the diameter of structure so as to provide optimal accuracy |
| `'ske_layer'` | Defines the layering method for rapid skeletonization. | Recommended that this value **does not exceed 32**. |
| `'loc'` | Determines the percentage of image information loss. | E.g., `0.1` corresponds to no more than 10% information loss. |
| `'win_th'` | The calculation window threshold used to determine the spatial/frequency domain calculation for **Local Coverage (LC)**. | N/A |

3.  **Output Values:**
    The average computed features are represented by the following variables:

| Feature | Variable |
| :--- | :--- |
| **Average Variance** | `'Vvaluelocal'` / `'Vvalueentire'` |
| **Average Waviness** | `'Meanwav'` |
| **Average Local Coverage** | `'Meanlc'` |

---

## 3. 结果 (Results)

The resulting **pseudo-color encoded maps** for the five calculated features are saved in the corresponding folder within: `ACA mapping\Results\`
