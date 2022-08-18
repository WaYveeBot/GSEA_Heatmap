import logging
import fhtbioinfpy
import fhtbioinfpy.setup_logger as setup_logger
import argparse
import sys
import glob
import os.path
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.write_gct as write_gct
import pandas
import numpy
import datetime



logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--verbose", "-v", help="Whether to print a bunch of output.", action="store_true", default=False)
    inputArgGroup = parser.add_mutually_exclusive_group(required = True)
    inputArgGroup.add_argument("--heatmap_dir", help = "Directory where heatmap input files can be found.", type = str, default = None)
    inputArgGroup.add_argument("--input_heatmap_files", help = "Sepecify heatmap files directly for input.", type = str, nargs = "+", default = None)
    parser.add_argument("--metadata_path", help = "File path to differential metadata file.", type = str, required = True)
    parser.add_argument("--output_dir", help = "Directory where output heatmaps should be stored.", type = str, default = "./gsea_heatmaps") 
    #Could remove requirement and use a default output directory
    
    return parser
    
#Searches given directory for any files that match heatmap naming    
def heatmap_file_search(heatmap_dir):
    logger.debug("Heatmap direcory input: {}".format(heatmap_dir))
    input_search = os.path.join(heatmap_dir, "*_heatmap_GSEA_*_r*.gct")        
    logger.debug("Input search: {}".format(input_search))
    heatmap_files = glob.glob(input_search)
    heatmap_files.sort()
    logger.debug("len(heatmap_files): {}".format(len(heatmap_files)))
    logger.debug("\n heatmap_files[:3]: {}".format(heatmap_files[:3]))
    return heatmap_files

#Retrieves metadata csv file from given file path and converts to dataframe
def load_differential_metadata(metadata_path):
    logger.debug("Metadata_path: {}".format(metadata_path))
    differential_metadata_dataframe = pandas.read_csv(metadata_path, sep="\t", index_col="chd")
    logger.debug("Metadata Dataframe Shape: {}".format(differential_metadata_dataframe.shape))
    return differential_metadata_dataframe

#convert heatmap files to (file name, GCToo object) tuples
def parse_heatmaps(heatmap_files):
    heatmap_GCT_list = [(file_name, parse.parse(file_name)) for file_name in heatmap_files]
    logger.debug("Heatmap List length: {}".format(len(heatmap_GCT_list)))
    logger.debug("Heatmap Data Shape: {}".format([GCT.data_df.shape for _,GCT in heatmap_GCT_list]))
    return heatmap_GCT_list

def make_new_heatmaps(heatmap_GCT_list, differential_metadata_dataframe):                                    
    new_heatmap_list = []
    for file_name, GCT in heatmap_GCT_list[:]:
        file_name = os.path.basename(file_name)
        logger.debug("File Name: {}".format(file_name))
        transposed_meta_dataframe = differential_metadata_dataframe.transpose()
        transposed_meta_dataframe["new_index"] = None

        for cid in GCT.col_metadata_df.index[:]:
            dge_stat = GCT.col_metadata_df.dge_statistic[cid]
            logger.debug("DGE Statistic: {}".format(dge_stat))
            
            logger.debug("cid: {}".format(cid))

            match_name = "_".join(cid[len(dge_stat)+1:].split("_")[:-2])

            logger.debug("Match Name: {}".format(match_name))

            locations = transposed_meta_dataframe.dge_output_filename.str.startswith(match_name)
            sum_locations = sum(locations)

            if sum_locations != 1:
                msg = "Unable to find match between column of GSEA heatmap and differentual metadata. \nFile Name: {} \nColumn Identifier: {} \nMatch Name: {} \nSum Locations: {} \nPossible Matches: {}".format(file_name, cid, match_name, sum_locations, locations)
                logger.exception(msg)
                raise fhtbioinfpyMakeGSEAHeatmapCannotMatchToMetadataException(msg)
                
            
            transposed_meta_dataframe.loc[locations, "new_index"] = cid

        transposed_meta_dataframe.set_index("new_index", inplace = True)
        transposed_meta_dataframe.index_name = "cid"

        assert set(GCT.data_df.columns) == set(transposed_meta_dataframe.index)

        new_GCT = GCToo.GCToo(GCT.data_df, row_metadata_df = GCT.row_metadata_df, col_metadata_df = transposed_meta_dataframe.loc[GCT.data_df.columns])

        logger.debug("Original Heatmap: {}".format(GCT))
        logger.debug("New Heatmap: {}".format(new_GCT))

        new_heatmap_list.append((file_name, new_GCT))

    logger.debug("New Heatmap List Length: {}".format(len(new_heatmap_list)))
    return new_heatmap_list


def rename_columns(new_heatmap_list):
    for file_name, new_GCT in new_heatmap_list[:]:
        logger.debug("File name: {}".format(file_name))
        new_GCT.col_metadata_df.set_index("friendly_name", inplace = True)

        logger.debug("Column metadata head: \n {}".format(new_GCT.col_metadata_df.head()))
        new_GCT.data_df.columns = new_GCT.col_metadata_df.index

    return new_heatmap_list

def write_output(output_dir, new_heatmap_list):
    logger.debug("Entered/default output directory: {}".format(output_dir))
    logger.debug("Final output directory: {}".format(output_dir))

    if (os.path.exists(output_dir) == False):
        os.mkdir(output_dir)

    for file_name, new_GCT in new_heatmap_list[:]:
        logger.debug("Saving file: {}".format(file_name))
        write_gct.write(new_GCT, os.path.join(output_dir, file_name))
    
    return output_dir
    
class fhtbioinfpyMakeGSEAHeatmapCannotMatchToMetadataException(Exception):
    pass

def main(args):

    #Retrieve heatmap files
    heatmap_files = args.input_heatmap_files if args.input_heatmap_files is not None else heatmap_file_search(args.heatmap_dir)
    logger.debug("Heatmap files: {}".format(heatmap_files))
    #Retrieve metadata CSV file and convert to dataframe
    metadata_dataframe = load_differential_metadata(args.metadata_path)

    #Convert heatmap files to list of tuples
    heatmap_GCT_list = parse_heatmaps(heatmap_files)

    #
    new_heatmap_list = make_new_heatmaps(heatmap_GCT_list, metadata_dataframe)

    #
    new_heatmap_list = rename_columns(new_heatmap_list)
    
    #Write output to default or passed-in directory
    write_output(args.output_dir, new_heatmap_list)

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])

    setup_logger.setup(verbose=args.verbose)

    logger.debug("args:  {}".format(args))

    main(args)