import os
import traceback
from totalsegmentator.python_api import totalsegmentator

# Python >= 3.10
# pip install TotalSegmentator
# https://github.com/wasserth/totalsegmentator
# A Nifti file or a folder (or zip file) with all DICOM slices of one patient is allowed as input.
# Output folder All vertebrae will be segmented with the prefix vertebrae_x in the output folder.

class SpineSegmentator:

    def __init__(self,
                 model_dir):

        self.model_dir = model_dir

        # Model path
        os.environ["TOTALSEG_HOME_DIR"] = self.model_dir
        # Temp files folder
        os.environ["TMPDIR"] = "D:\\Spine_Seg_Results\\temp_files"
        os.environ["TEMP"] = "D:\\Spine_Seg_Results\\temp_files"
        os.environ["TMP"] = "D:\\Spine_Seg_Results\\temp_files"

    def segment_spine(self,
                input_ct_path,
                output_dir):

        try:

            print("Model dir:", self.model_dir)
            print("Input CT:", input_ct_path)
            print("Output dir:", output_dir)

            os.makedirs(output_dir, exist_ok=True)

            totalsegmentator(
                input=input_ct_path,
                output=output_dir,
                task="total",
                fast=True
            )

            print("Segmentation finished")

        except Exception as e:

            print("ERROR:")
            print(e)

            traceback.print_exc()


segmentator = SpineSegmentator(
    model_dir="D:\\TotalSegmentator"
)

segmentator.segment_spine(
    input_ct_path="C:\\Users\Miguel\\Downloads\\Spine_CT.nii.gz",
    output_dir="D:\\Spine_Seg_Results"
)