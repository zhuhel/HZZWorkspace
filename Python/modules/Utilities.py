from ConfigParser import SafeConfigParser
import logging
import os


def config_dict(filename, tokenize=None):
    """
    A simple way to get a dictionary from a config file.
    If a list of strings is passed in tokenize, then those items will be split
    by commas and stripped of whitespace.
    """
    config = SafeConfigParser()
    config.optionxform = str
    return_dict = None

    logging.info("Attempting to read file %s...", filename)
    try:
        with open(filename) as f:
            config.readfp(f)
            return_dict = config._sections
            logging.info("Success!")
    except IOError:
        logging.error("Could not find file %s. Consider using absolute paths in your config file.", filename)
        raise

    if return_dict:
        for s, section in return_dict.iteritems():
            del section['__name__']
            if tokenize:
                assert isinstance(tokenize, list)
                for v, value in section.iteritems():
                    for tok in tokenize:
                        if v in tok[0]:
                            section[v] = [i.strip() for i in value.split(tok[1])]
        return return_dict
    else:
        logging.error("Config empty, something has gone wrong.")


def config_display(cfg_dict):
    """
    A nice way to display the given config dictionary.
    :param cfg_dict: the input config dictionary, as in the config_dict function above
    :return: a formatted string with the sections and values
    """
    return_str = "Listing config..."
    for s, section in cfg_dict.iteritems():
        return_str += "\n" + s + ":"
        for v, value in section.iteritems():
            return_str += "\n-- " + v + " = " + str(value)
    return return_str


def check_config(cfg_dict, keys, cfg_section='main'):
    """
    Check if the config dictionary contains the given keys.
    :param cfg_dict: the input config dictionary, as in the config_dict function above
    :param keys: the keys in cfg_section that need to be present
    :param cfg_section: the section of the config file to search
    :return: True if the config is fine, False otherwise
    """
    accept = True
    for k in keys:
        if k not in cfg_dict[cfg_section]:
            logging.error("Could not find '%s' option in '%s' section of config file", k, cfg_section)
            accept = False
    return accept


def check_np(np_dict):
    """
    Checks to see if the NP config file is appropriately designed.
    :param np_dict: The config dict coming from the config_dict function above
    :return: a boolean, True if the dict is fine, False if not
    """
    for s, section in np_dict.iteritems():
        var_sum = 0
        for v, value in section.iteritems():
            if "JET_JER_SINGLE_NP" in v:
                continue
            if "1up" in v:
                var_sum += 1
            elif "1down" in v:
                var_sum -= 1

        if var_sum == 0:
            return True
        else:
            return False


def check_and_mkdir(directory):
    """
    Asserts whether a given directory exists. If not, it is created.
    :directory: the directory to check.
    """
    if not os.path.exists(directory):
        logging.info("Creating directory '%s'", directory)
        os.makedirs(directory)


def data_to_plotpoints(NP_data, prune=0.01, skip="Nominal", replace_str="ATLAS_"):
    """
    Converts an NP dictionary to a list of tuples in the form 
    (name, percentage down variation, percentage up variation)
    :NP_data: the dictionary that has the NP data
    :prune: reject variations less than this value
    :skip_nominal: skips a nuisance parameter called this value
    :replace_str: replaces this text with blanks, to neaten the output
    :return: a list of tuples with the percentage variations.
    """
    plot_points = []
    for np, var in NP_data.iteritems():
        if np == skip:
            continue
        # Store the percentage deviation
        up_variation = 100*(max([var['up']['norm'] - 1.0, var['down']['norm'] - 1.0, 0.0]))
        down_variation = 100*(min([var['up']['norm'] - 1.0, var['down']['norm'] - 1.0, 0.0]))
        if up_variation - down_variation < prune:
            continue
        plot_points.append((np.replace("ATLAS_", ""), down_variation, up_variation))
    return plot_points
