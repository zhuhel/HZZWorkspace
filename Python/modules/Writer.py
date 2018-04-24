import Utilities, logging


def NP_latex_column(NP_data, filename, header, output_dir=".", prune=0.01):
    """
    A function to write out the systematic variations to a LaTeX formatted column.
    :NP_data: the input dictionary
    :filename: the output filename
    """
    logging.info("Writing LaTeX table of NP systematic variations to %s/%s", output_dir, filename)
    Utilities.check_and_mkdir(output_dir)
    entries = Utilities.data_to_plotpoints(NP_data, prune)

    
