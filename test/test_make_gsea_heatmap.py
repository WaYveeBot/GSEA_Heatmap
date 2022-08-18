import unittest
import logging
import fhtbioinfpy.setup_logger as setup_logger
import fhtbioinfpy.make_gsea_heatmap.make_gsea_heatmap as mgh
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.write_gct as write_gct
import pandas
import tempfile
import os

temp_wkdir_prefix = "TestMakeGSEAHeatmap"
logger = logging.getLogger(setup_logger.LOGGER_NAME)

class testMakeGSEAHeatmap(unittest.TestCase):
    def test_main(self):
        logger.debug("\n \n \n test_main \n \n ")
        
        heatmap_dir = "assets/orig_heatmaps"
        
        metadata_path = "assets/NS-22.0003_differential_metadata_r40x15.txt"
        with tempfile.TemporaryDirectory(prefix=temp_wkdir_prefix) as wkdir:
            logger.debug("wkdir: {}".format(wkdir))
            args = mgh.build_parser().parse_args(["--heatmap_dir", heatmap_dir, "--metadata_path", metadata_path, "--output_dir", wkdir, "--verbose"])

            mgh.main(args)

            #Retrieve files that the program outputted:
            outputted_file_names = [
                "NS-22.0003_heatmap_GSEA_fdr_r21633x15.gct",
                "NS-22.0003_heatmap_GSEA_matched_size_r21633x15.gct",
                "NS-22.0003_heatmap_GSEA_nes_r21633x15.gct",
                "NS-22.0003_heatmap_GSEA_pval_r21633x15.gct"
            ]
            outputted_file_paths = [os.path.join(wkdir, file) for file in outputted_file_names]

            #Check that files were ouputted sucessfully
            for file in outputted_file_paths:
                self.assertTrue(os.path.exists(file))

            #Retrieve sample files to compare against:
            expected_files = [os.path.join("assets", "example_output_heatmaps", file) for file in outputted_file_names]
            
            #Compare the two lists to make the contents of both files match:
            for out_file, exp_file in zip(outputted_file_paths, expected_files):
                logger.debug("Output file: {}".format(out_file))
                logger.debug("Expected file: {}".format(exp_file))
                self.assertEqual(parse.parse(out_file).data_df.shape, parse.parse(exp_file).data_df.shape)

    def test_heatmap_file_search(self):
        logger.debug("\n \n \n Testing heatmap_file_search: \n \n ")

        heatmap_dir = "assets/orig_heatmaps"
        heatmap_files = mgh.heatmap_file_search(heatmap_dir)
        self.assertEqual(len(heatmap_files), 4)
        
    def test_metadata_file_search(self):
        logger.debug("\n \n \n Testing metadata_file_search: \n \n ")

        metadata_path = "assets/NS-22.0003_differential_metadata_r40x15.txt"
        metadata_dataframe = mgh.load_differential_metadata(metadata_path)
        self.assertEqual(metadata_dataframe.shape, (40, 15))

    def test_parse_heatmaps(self):
        logger.debug("\n \n \n Testing Heatmap Parse: \n \n  ")
        heatmap_dir = "assets/orig_heatmaps"
        heatmap_files = mgh.heatmap_file_search(heatmap_dir)
        heatmap_GCT_list = mgh.parse_heatmaps(heatmap_files)    
        self.assertEqual(len(heatmap_GCT_list), 4)

    def test_make_new_heatmaps(self):

        test_files = [
                os.path.join("assets", "orig_heatmaps", "NS-22.0003_heatmap_GSEA_fdr_r21633x15.gct"),
                os.path.join("assets", "orig_heatmaps", "NS-22.0003_heatmap_GSEA_matched_size_r21633x15.gct"),
                os.path.join("assets", "orig_heatmaps", "NS-22.0003_heatmap_GSEA_nes_r21633x15.gct"),
                os.path.join("assets", "orig_heatmaps", "NS-22.0003_heatmap_GSEA_pval_r21633x15.gct")
            ]


        logger.debug("\n \n \n Testing make_new_heatmaps: \n \n ")
        test_heatmap_list = mgh.parse_heatmaps(test_files)
        metadata_path = "assets/NS-22.0003_differential_metadata_r40x15.txt"
        test_meta_file = pandas.read_csv(metadata_path, sep="\t", index_col=0)
        logger.debug("Heatmap shapes: \n {}".format([GCT.data_df.shape for _,GCT in test_heatmap_list]))
        test_heatmap_list_out = mgh.make_new_heatmaps(test_heatmap_list, test_meta_file)
        logger.debug("Output Heatmap shapes: \n {}".format([GCT.data_df.shape for _,GCT in test_heatmap_list_out]))
        #logger.debug("".format())
        #logger.debug("".format())
        #logger.debug("".format())
        #logger.debug("".format())


    def test_rename_columns(self):
        logger.debug("\n \n \n Testing rename columns \n \n ")
        expected_files = [
                os.path.join("assets", "example_output_heatmaps", "NS-22.0003_heatmap_GSEA_fdr_r21633x15.gct"),
                os.path.join("assets", "example_output_heatmaps", "NS-22.0003_heatmap_GSEA_matched_size_r21633x15.gct"),
                os.path.join("assets", "example_output_heatmaps", "NS-22.0003_heatmap_GSEA_nes_r21633x15.gct"),
                os.path.join("assets", "example_output_heatmaps", "NS-22.0003_heatmap_GSEA_pval_r21633x15.gct")
            ]

        test_files = [
                os.path.join("assets", "orig_heatmaps", "NS-22.0003_heatmap_GSEA_fdr_r21633x15.gct"),
                os.path.join("assets", "orig_heatmaps", "NS-22.0003_heatmap_GSEA_matched_size_r21633x15.gct"),
                os.path.join("assets", "orig_heatmaps", "NS-22.0003_heatmap_GSEA_nes_r21633x15.gct"),
                os.path.join("assets", "orig_heatmaps", "NS-22.0003_heatmap_GSEA_pval_r21633x15.gct")
            ]    
        
        metadata_path = "assets/NS-22.0003_differential_metadata_r40x15.txt"
        test_heatmap_list = mgh.parse_heatmaps(test_files)
        metadata_dataframe = mgh.load_differential_metadata(metadata_path)
        test_heatmap_list = mgh.make_new_heatmaps(test_heatmap_list, metadata_dataframe)
        logger.debug("Test Heatmap Shape: \n{}".format([test_heatmap.data_df.shape for _,test_heatmap in test_heatmap_list]))
        logger.debug("Test Heatmap Head: \n{}".format([test_heatmap.data_df.head() for _,test_heatmap in test_heatmap_list]))
        test_heatmap_list = mgh.rename_columns(test_heatmap_list)
        logger.debug("Test Heatmap Output Shape: \n{}".format([test_heatmap.data_df.shape for _,test_heatmap in test_heatmap_list]))
        logger.debug("Test Heatmap Output Head: \n{}".format([test_heatmap.data_df.head() for _,test_heatmap in test_heatmap_list]))
        logger.debug("Testing if columns are renamed properly: ")

        #Test to make sure columns are renamed correctly
        for i in range(0, len(test_heatmap_list)):
            output = test_heatmap_list[i][1]
            logger.debug("Comparing {} to expected".format(test_heatmap_list[i][0]))
            expected = parse.parse(expected_files[i])
            self.assertTrue(output.data_df.columns.equals(expected.data_df.columns))

    def test_write_output(self):
        logger.debug("\n \n \n Testing writing output \n \n ")
       
        heatmap_dir = "assets/orig_heatmaps"
        
        metadata_path = "assets/NS-22.0003_differential_metadata_r40x15.txt"
        with tempfile.TemporaryDirectory(prefix=temp_wkdir_prefix) as wkdir:
            heatmap_files = mgh.heatmap_file_search(heatmap_dir)
            metadata_dataframe = mgh.load_differential_metadata(metadata_path)
            heatmap_GCT_list = mgh.parse_heatmaps(heatmap_files)
            new_heatmap_list = mgh.make_new_heatmaps(heatmap_GCT_list, metadata_dataframe)
            new_heatmap_list = mgh.rename_columns(new_heatmap_list) 
            logger.debug("Heatmap list: {}".format(new_heatmap_list[0][1]))
            mgh.write_output(wkdir, new_heatmap_list)
        
        




if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()