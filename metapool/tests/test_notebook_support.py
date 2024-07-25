import sys
import pandas as pd
from pandas.testing import assert_frame_equal
import os

from unittest import TestCase, main
from metapool.notebook_support import _generate_sample_context


class NotebookSupportTests(TestCase):
    def test__generate_sample_context(self):
        # note that test must be run from the root of the metapool module
        # for this path to resolve
        path = os.path.dirname(__file__)
        plate_df_path = os.path.join(path,
            '../../notebooks/test_output/QC/YYYY_MM_DD_Celeste_Adaptation_df_A.txt')
            # 'notebooks/test_output/QC/YYYY_MM_DD_Celeste_Adaptation_df_A.txt'

        plate_df = pd.read_csv(plate_df_path, dtype="str", sep='\t')

        generated_context = _generate_sample_context(plate_df)
        self.assertTrue(False)
