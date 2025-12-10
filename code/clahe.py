import cv2
from matplotlib import pyplot as plt

img_path = "/content/drive/MyDrive/raw_images/20251208_165259.bmp"

img = cv2.imread(img_path, cv2.IMREAD_GRAYSCALE)
if img is None:
    print("Cannot open the image")
else:
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8,8))
    clahe_img = clahe.apply(img)

    plt.figure(figsize=(12,6))

    plt.subplot(1,2,1)
    plt.title("Original")
    plt.imshow(img, cmap='gray')
    plt.axis('off')

    plt.subplot(1,2,2)
    plt.title("CLAHE")
    plt.imshow(clahe_img, cmap='gray')
    plt.axis('off')

    plt.tight_layout()
    plt.show()
