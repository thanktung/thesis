from pypylon import pylon
import cv2

def basler_capture_one():
    try:
        cam = pylon.InstantCamera(pylon.TlFactory.GetInstance().CreateFirstDevice())
        cam.Open()

        converter = pylon.ImageFormatConverter()
        converter.OutputPixelFormat = pylon.PixelType_BGR8packed
        converter.OutputBitAlignment = pylon.OutputBitAlignment_MsbAligned

        cam.StartGrabbingMax(1)
        grab = cam.RetrieveResult(3000, pylon.TimeoutHandling_ThrowException)

        if grab.GrabSucceeded():
            img = converter.Convert(grab).GetArray()
            grab.Release()
            cam.Close()
            return img

        grab.Release()
        cam.Close()
        return None

    except Exception as e:
        print("Basler capture error:", e)
        return None
